import torch
import torch.nn as nn
import numpy as np
import pytorch_lightning as pl
from collections import defaultdict

from Models.Structure.structure import Graph
from Models.common import MLP

NONLINEARITIES = {
    "tanh": nn.Tanh(),
    "relu": nn.ReLU(),
    "softplus": nn.Softplus(),
    "elu": nn.ELU()
}

class Aggregator(nn.Module):
    def __init__(self, input_dim, output_dim, hidden_dim, num_layer=1, norm='None', act_last=False, act_fn='relu'):
        super().__init__()
        self.mlp = MLP(in_dim=input_dim, out_dim=output_dim, hidden_dim=hidden_dim, num_layer=num_layer, norm=norm, act_last=act_last, act_fn=act_fn)
        # self.norm = nn.BatchNorm1d(input_dim)

    def forward(self, x):
        x = self.mlp(torch.cat(x, dim=-1))
        return x

class SS(pl.LightningModule):
    def __init__(self, config):
        super().__init__()

        # Build the network
        self.config = config
        self.hidden_dim = config.model.hidden_dim

        cross_attention_config = config.model.cross_attention
        graph_config = config.model.graph
        transform_config = config.transform
        header_config = config.model.header

        # component flag
        self.use_gnn = config.model.use_gnn
        self.use_attention = True
        try:
            self.use_attention = cross_attention_config.use_attention
        except:
            pass

        # sequence_encoder
        ## environment_net
        self.graph_net = Graph(graph_config, transform_config)
        # output: [n_atom * embeding_length (10*128)]

        ## enzyme_net
        self.protein_mlp = nn.Linear(1280, self.hidden_dim)
        # output: [n_amino_acid * embeding_length (300*128)]

        # substarte atom based feature (sequence-based atom features, structure-based atom features, external atom based embedding (grover))
        self.atom_feature_mlp_dict = nn.ModuleDict(
            {
                self.config.data.atom_features[i-1]: nn.Sequential(
                    nn.Linear(int(self.config.data.atom_features[i]), self.hidden_dim),
                    nn.LayerNorm(self.hidden_dim) if self.config.model.feature_norm else nn.Identity()
                ) for i in range(1, len(self.config.data.atom_features), 2)
            }
        )
        self.atom_feature_aggregator = Aggregator(input_dim=self.hidden_dim*(int(self.config.model.use_gnn) + len(self.config.data.atom_features)//2), output_dim=self.hidden_dim, hidden_dim=self.hidden_dim, num_layer=1)
        # output: [n_atom * embeding_length (10*128)]

        # cross attention (reaction atom based feature, enzyme atom based feature)
        if self.use_attention:
            self.enzyme_attention = nn.MultiheadAttention(embed_dim=self.hidden_dim, num_heads=cross_attention_config. n_head, batch_first=True, dropout=cross_attention_config.dropout)
            self.reaction_attention = nn.MultiheadAttention(embed_dim=self.hidden_dim, num_heads=cross_attention_config.n_head, batch_first=True, dropout=cross_attention_config.dropout)

        # specificity prediction
        # (weighted sum substrate sequence embedding, weighted sum enzyme sequence embedding, avg substrate sequence embedding, avg enzyme sequence embedding, structure embedding, grover embedding (substrate-based), morgan embedding)

        self.feature_mlp_dict = nn.ModuleDict(
            {
                self.config.data.features[i-1]: nn.Sequential(
                    nn.Linear(int(self.config.data.features[i]), self.hidden_dim),
                    nn.LayerNorm(self.hidden_dim) if self.config.model.feature_norm else nn.Identity()
                ) for i in range(1, len(self.config.data.features), 2)
            }
        )
        self.specificity_header = Aggregator(input_dim=self.hidden_dim*(2 + self.use_attention * 2 + self.config.model.use_gnn + len(self.config.data.features)//2), output_dim=1, hidden_dim=self.hidden_dim, num_layer=header_config.num_layers, norm=header_config.norm, act_fn=header_config.act_fn)
        # PL 2.x: manual step output caching
        self.validation_step_outputs = []
        self.test_step_outputs = []

    def _get_graph_output(self, G, B):
        ligand_embedding = torch.zeros(B, 281, self.hidden_dim).to(self.device)
        output, (h, batch, index) = self.graph_net(G)
        ligand_embedding[batch, index, :] = h
        return output, ligand_embedding[:, :280, :]

    def forward(self, G):

        # TODO: change to max enzyme length and max reaction length
        x_pro = self.protein_mlp(G.embedding)
        x_pro = x_pro.view(-1, 1450, self.hidden_dim)

        # get structure embedding (mean and atom embedding)

        # get atom embedding (sequence (x_reaction), structure (structure_x_reaction), grover (grover))
        atom_features = []
        if self.config.model.use_gnn:
            str_mean, x_str = self._get_graph_output(G, x_pro.shape[0])
            atom_features.append(x_str)

        if "grover" in self.config.data.atom_features:
            grover = G.grover.view(-1, 280, int(self.config.data.atom_features[1]))
            grover = self.atom_feature_mlp_dict["grover"](grover)
            atom_features.append(grover)

        x_reaction = self.atom_feature_aggregator(atom_features)
    
        # start attention
        if self.use_attention:
            weighted_sum_pro, _ = self.enzyme_attention(x_pro, x_reaction, 
                                                        x_reaction, need_weights=True, key_padding_mask=G.reaction_padding_mask)
            weighted_sum_reaction, _ = self.reaction_attention(x_reaction, x_pro, 
                                                        x_pro, need_weights=True, key_padding_mask=G.enzyme_padding_mask)
       
        ## calculate embedding after attention
        ### shape:
        ### x.enzyme_padding_mask: [B, len_pro]
        ### x.reaction_padding_mask: [B, len_reaction]
        ### x_pro: [B, len_pro, len_embed]
        ### x_reaction: [B, len_reaction, len_embed]
        reaction_mask = (~G.reaction_padding_mask).unsqueeze(-1).float()
        enzyme_mask = (~G.enzyme_padding_mask).unsqueeze(-1).float()

        ### reaction_mask: [B, len_reaction, 1]
        ### enzyme_mask: [B, len_pro, 1]
        x_reaction = (x_reaction * reaction_mask).sum(dim=1) / reaction_mask.sum(dim=(1,2)).unsqueeze(-1)
        x_pro = (x_pro * enzyme_mask).sum(dim=1) / enzyme_mask.sum(dim=(1,2)).unsqueeze(-1)
        
        if self.use_attention:
            weighted_sum_reaction = (weighted_sum_reaction * reaction_mask).sum(dim=1) / reaction_mask.sum(dim=(1,2)).unsqueeze(-1)
            weighted_sum_pro = (weighted_sum_pro * enzyme_mask).sum(dim=1) / enzyme_mask.sum(dim=(1,2)).unsqueeze(-1)

        # specificity output
        embeddings = [x_pro, x_reaction]
        if self.use_attention:
            embeddings.extend([weighted_sum_pro, weighted_sum_reaction])

        if self.config.model.use_gnn:
            embeddings.extend([str_mean])
            # embeddings.append(str_mean)
            
        if "grover_mean" in self.config.data.features:
            embeddings.append(self.feature_mlp_dict["grover_mean"](G.grover_mean))
        if "morgan" in self.config.data.features:
            embeddings.append(self.feature_mlp_dict["morgan"](G.morgan))
        
        specificity_output = self.specificity_header(embeddings)
        return specificity_output, [G.tag, G.str_tag]

    def get_loss(self, x, logits, stage):
        logits = logits.squeeze(-1)

        loss = (torch.nn.BCEWithLogitsLoss(reduction='none')(logits.ravel(), x.label.ravel()) * x.sample_weight).mean()
        self.log(f"sp_loss/{stage}", loss)
        return loss

    def training_step(self, batch, batch_idx):
        logits, tag = self(batch)
        loss = self.get_loss(batch, logits, "train")
        return loss
        
    def evaluate(self, x, stage=None):
        sp_logits, tag = self(x)
        self.get_loss(x, sp_logits, stage)
        return sp_logits.detach().cpu(), x.label.detach().cpu(), tag
    
    def test_step(self, batch, batch_idx):
        output = self.evaluate(batch, 'test')
        self.test_step_outputs.append(output)
        return output

    def validation_step(self, batch, batch_idx):
        output = self.evaluate(batch, 'val')
        self.validation_step_outputs.append(output)
        return output

    def _print(self, a):
        a = a.tolist()
        a = [str(x) for x in a]
        print(",".join(a))

    def get_auc_aupr(self, logits, label, name, stage):
        # logits: [B]
        # label: [B]
        from sklearn import metrics
        logits = logits.ravel()
        label = label.ravel()
        
        fpr, tpr, thresholds = metrics.roc_curve(label, logits, pos_label=1)
        # precision, recall, thresholds = metrics.precision_recall_curve(label, logits)
        if name is not None:
            self.log(f"auc/{name}/{stage}", metrics.auc(fpr, tpr))
            self.log(f"aupr/{name}/{stage}", metrics.average_precision_score(label, logits))
            print(f"auc/{name}/{stage}", metrics.auc(fpr, tpr))
        else:
            self.log(f"auc/{stage}", metrics.auc(fpr, tpr))
            self.log(f"aupr/{stage}", metrics.average_precision_score(label, logits))
            print(f"auc/{stage}", metrics.auc(fpr, tpr))

        return metrics.auc(fpr, tpr), metrics.average_precision_score(label, logits)

    def cutomized_epoch_end(self, outputs, stage):
        sp_logits = np.concatenate([a[0] for a in outputs], axis=0)
        sp_labels = np.concatenate([a[1] for a in outputs], axis=0)
        
        self.get_auc_aupr(sp_logits, sp_labels, None, stage)
    
        logits_dict = defaultdict(list)
        labels_dict = defaultdict(list)
 
        for index, output in enumerate(outputs):
            for logits, label, tag_0, tag_1 in zip(output[0], output[1], output[2][0], output[2][1]):
                tag = (tag_0, tag_1)
                logits_dict[tag[0]].append(logits)
                labels_dict[tag[0]].append(label)
                logits_dict[tag[0] + "-" + tag[1]].append(logits)
                labels_dict[tag[0] + "-" + tag[1]].append(label)
                
                logits_dict[tag[1]].append(logits)
                labels_dict[tag[1]].append(label)

        for key, item in logits_dict.items():
            logits = np.array([t.item() for t in item])
            labels = np.array(labels_dict[key])
            self.get_auc_aupr(logits, labels, key, stage)

        self.logits = sp_logits
        self.labels = sp_labels
        self.logits_dict = logits_dict
        self.labels_dict = labels_dict

    def on_validation_epoch_start(self):
        self.validation_step_outputs = []

    def on_test_epoch_start(self):
        self.test_step_outputs = []

    def on_validation_epoch_end(self):
        self.cutomized_epoch_end(self.validation_step_outputs, 'val')
        self.validation_step_outputs = []

    def on_test_epoch_end(self):
        self.cutomized_epoch_end(self.test_step_outputs, 'test')
        self.test_step_outputs = []

    def configure_optimizers(self):
        from warmup_scheduler import GradualWarmupScheduler

        if self.config.training.optimizer == "adam":
            self.optimizer = optimizer = torch.optim.Adam(self.parameters(), lr=self.config.training.lr)
        elif self.config.training.optimizer == "sgd":
            self.optimizer = optimizer = torch.optim.SGD(self.parameters(), lr=self.config.training.lr)
        elif self.config.training.optimizer == "rmsport":
            self.optimizer = optimizer = torch.optim.RMSprop(self.parameters(), lr=self.config.training.lr)
        elif self.config.training.optimizer == "adamW":
            self.optimizer = optimizer = torch.optim.AdamW(self.parameters(), lr=self.config.training.lr)

        second_scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer,
                                                               factor=self.config.training.sched_factor,
                                                               patience=self.config.training.sched_patience,
                                                               min_lr=self.config.training.min_lr, mode='max'
                                                               )
        # first_scheduler = GradualWarmupScheduler(optimizer, multiplier=1, total_epoch=10, after_scheduler=second_scheduler)
        # print(first_scheduler.get_lr())
        lr_scheduler = {
          "scheduler": second_scheduler, 
          'interval': 'epoch',
          'frequency': self.config.training.val_frequency,
          'strict': True,
          "monitor": "aupr/val"
        }                
        return [optimizer], [lr_scheduler]
    
    def optimizer_step(self, epoch, batch_idx, optimizer, optimizer_closure):
        if epoch < self.config.training.warmup_epochs:
            lr = 1. * epoch / self.config.training.warmup_epochs * (self.config.training.lr - self.config.training.min_lr) + self.config.training.min_lr
            for pg in optimizer.param_groups:
                pg['lr'] = lr
        elif epoch == self.config.training.warmup_epochs:
            for pg in optimizer.param_groups:
                pg['lr'] = self.config.training.lr
        optimizer.step(closure=optimizer_closure)