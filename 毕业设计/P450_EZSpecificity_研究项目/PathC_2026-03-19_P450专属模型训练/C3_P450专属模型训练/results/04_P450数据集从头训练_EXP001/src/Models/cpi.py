import logging

import torch
import numpy as np
import torch.nn as nn

import pytorch_lightning as pl

class CPI(pl.LightningModule):
    """FFN."""

    def __init__(
        self,
        config,
        **kwargs,
    ):
        """__init__.
        Args:
        """
        super(CPI, self).__init__(**kwargs) 
        
        self.config = config

        self.best_val_auc = 0
        self.best_val_aupr = 0
        self.logits = None
        self.labels = None

        layers = config.model.layers
        hidden_size = config.model.hidden_size
        model_dropout = config.model.model_dropout
        num_tasks = config.model.num_tasks
        prot_emb_size = config.model.prot_emb_size
        chem_emb_size = config.model.chem_emb_size
        interaction_layer = config.model.interaction_layer

        self.interaction_layer = interaction_layer
        if interaction_layer == "concat":
            in_features = prot_emb_size + chem_emb_size
            self.input_layer = nn.Linear(
                in_features=in_features, out_features=hidden_size, bias=True
            )
        elif interaction_layer == "dot":
            self.input_layer_prot = nn.Linear(
                in_features=prot_emb_size, out_features=hidden_size, bias=True
            )
            self.input_layer_chem = nn.Linear(
                in_features=chem_emb_size, out_features=hidden_size, bias=True
            )
        else:
            raise ValueError(f"Unexpected value {interaction_layer}")

        ## Middle layers
        sequential_layers = []
        for i in range(layers - 1):
            # Add activation before each middle layer s.t. we can use this for
            # standard logistic regression without any activation
            sequential_layers.append(torch.nn.ReLU())
            new_layer = nn.Linear(
                in_features=hidden_size, out_features=hidden_size, bias=True
            )
            sequential_layers.append(new_layer)
            sequential_layers.append(nn.Dropout(p=model_dropout))

        self.inner_layers = nn.Sequential(*sequential_layers)

        # Have multiple outputs
        self.num_tasks = num_tasks
        self.out_features = self.num_tasks
        self.output_layer = nn.Sequential(
            torch.nn.ReLU(),
            nn.Linear(
                in_features=hidden_size, out_features=self.out_features, bias=True
            ),
        )

    def forward(self, batch, **kwargs):
        """Forward pass, return logits"""

        rxn_features, prot_features = batch
        
        feature_ar = []
        if self.interaction_layer == "concat":
            feature_ar.append(rxn_features)
            feature_ar.append(prot_features)

            # Concatenate the features then run module
            in_features = torch.cat(feature_ar, 1)

            output = self.input_layer(in_features)

        elif self.interaction_layer == "dot":
            feature_ar.append(prot_features)

            if len(feature_ar) <= 0:
                raise RuntimeError("Expected prot or ec features for ffn network")

            # Concatenate the features then run module
            rxn_repr = rxn_features

            protein_repr = torch.cat(feature_ar, 1)

            protein_repr = self.input_layer_prot(protein_repr)
            rxn_repr = self.input_layer_chem(rxn_repr)
            # Take element wise product between representations!
            output = torch.einsum("bd,bd->bd", protein_repr, rxn_repr)
        else:
            raise ValueError(f"Bad interaction layer value: {self.interaction_layer}")

        output = self.inner_layers(output)

        output = self.output_layer(output)

        return output
    
    def get_loss(self, label, logits):
        loss = nn.BCEWithLogitsLoss()(logits.ravel(), label.ravel())
        return loss

    def training_step(self, batch, batch_idx):
        reaction = batch.x
        enzyme = batch.enzyme
        label = batch.label
        loss = self.get_loss(label, self((reaction, enzyme)))
        return loss
        
    def evaluate(self, batch, stage=None):
        reaction = batch.x
        enzyme = batch.enzyme
        label = batch.label
        tag = batch.tag
        prediction = self((reaction, enzyme))
        return prediction.detach().cpu(), label.detach().cpu(), tag.detach().cpu()

    def test_step(self, batch, batch_idx):
        return self.evaluate(batch, 'test')

    def validation_step(self, batch, batch_idx):
        return self.evaluate(batch, 'val')
    
    def get_auc_aupr(self, logits, label, name, stage):
        # logits: [B]
        # label: [B]
        from sklearn import metrics
        logits = logits.ravel()
        label = label.ravel()
        
        fpr, tpr, thresholds = metrics.roc_curve(label, logits, pos_label=1)
        # precision, recall, thresholds = metrics.precision_recall_curve(label, logits)
        if name is not None:
            self.log(f"{stage}_{name}_auc", metrics.auc(fpr, tpr))
            self.log(f"{stage}_{name}_aupr", metrics.average_precision_score(label, logits))
        else:
            self.log(f"{stage}_auc", metrics.auc(fpr, tpr))
            self.log(f"{stage}_aupr", metrics.average_precision_score(label, logits))

        auc = metrics.auc(fpr, tpr)
        aupr = metrics.average_precision_score(label, logits)

        if stage == 'val' and name is None:
            if auc > self.best_val_auc:
                self.best_val_auc = auc
                self.best_val_aupr = aupr
            self.log(f"best_{stage}_auc", self.best_val_auc)
            self.log(f"best_{stage}_aupr", self.best_val_aupr)

        return auc, aupr

    def cutomized_epoch_end(self, outputs, stage):
        from collections import defaultdict

        logits = np.concatenate([a[0] for a in outputs], axis=0)
        labels = np.concatenate([a[1] for a in outputs], axis=0)
        tags = np.concatenate([a[2] for a in outputs], axis=0)

        self.labels = labels
        self.logits = logits

        self.get_auc_aupr(logits, labels, None, stage)

        logits_dict = defaultdict(list)
        labels_dict = defaultdict(list)

        # for i in range(len(outputs[0][2])):
        
        for index, (label, tag) in enumerate(zip(labels, tags)):
            # if label == 1 and tag != 5:
            #     for i in range(5):
            #         logits_dict[i].append(self.logits[index])
            #         labels_dict[i].append(self.labels[index])
            # else:
            logits_dict[tag].append(self.logits[index])
            labels_dict[tag].append(self.labels[index])
                
        for key, item in logits_dict.items():
            logits = np.array(item)
            labels = np.array(labels_dict[key])
            self.get_auc_aupr(logits, labels, key, stage)

        self.logits_dict = logits_dict
        self.labels_dict = labels_dict

    def validation_epoch_end(self, outputs):
        self.cutomized_epoch_end(outputs, 'val')
    
    def test_epoch_end(self, outputs):
        self.cutomized_epoch_end(outputs, 'test')

    def configure_optimizers(self):
        from warmup_scheduler import GradualWarmupScheduler

        if self.config.training.optimizer == "adam":
            self.optimizer = optimizer = torch.optim.Adam(self.parameters(), lr=self.config.training.lr)
        elif self.config.training.optimizer == "sgd":
            self.optimizer = optimizer = torch.optim.SGD(self.parameters(), lr=self.config.training.lr)
        elif self.config.training.optimizer == "rmsport":
            self.optimizer = optimizer = torch.optim.RMSprop(self.parameters(), lr=self.config.training.lr)

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
          "monitor": "val_auc"
        }                
        return [optimizer], [lr_scheduler]