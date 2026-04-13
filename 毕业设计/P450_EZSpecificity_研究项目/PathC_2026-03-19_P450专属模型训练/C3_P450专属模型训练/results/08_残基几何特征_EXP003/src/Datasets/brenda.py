import pytorch_lightning as pl
from torch_geometric.data import DataLoader
from torch_geometric.transforms import Compose
from rdkit import RDLogger

from Datasets.data_representer import get_representer
import Datasets.Structure.transforms as utils_trans
from Datasets.utils import read_datasets

RDLogger.DisableLog('rdApp.*')  

class Singledataset(pl.LightningDataModule):
    def __init__(self, config):
        super().__init__()
        
        self.config = config

        self.batch_size = config.data.batch_size
        self.n_worker = config.num_cpus
        self.train_df = read_datasets(config.data.train_data_path)
        self.val_df = read_datasets(config.data.val_data_path)
        self.test_df = read_datasets(config.data.test_data_path)
        
        self.train_transform = self.get_transform()
        self.val_transform = self.get_transform(is_train=False)

    def get_transform(self, is_train=True):
        if self.config.data.representer == 'cpi':
            transform = None
        else:
            transform = Compose([
                utils_trans.FeaturizeProteinAtom(),
                utils_trans.FeaturizeLigandAtom(),
                utils_trans.EdgeConnection(dist_noise=self.config.transform.dist_noise if is_train else False, cutoff=self.config.transform.cutoff, num_r_gaussian=self.config.transform.num_r_gaussian, k=self.config.transform.k)
            ])
        return transform
        
    def train_dataloader(self):
        data = get_representer(df=self.train_df, config=self.config, transform=self.train_transform, is_train=True)
        return DataLoader(data, batch_size=self.batch_size, shuffle=True, num_workers=self.n_worker, follow_batch=['ligand_index'])

    def val_dataloader(self):
        data = get_representer(df=self.val_df, config=self.config, transform=self.val_transform, is_train=False)
        return DataLoader(data, batch_size=self.batch_size, shuffle=False, num_workers=self.n_worker, follow_batch=['ligand_index'])

    def test_dataloader(self):
        data = get_representer(df=self.test_df, config=self.config, transform=self.val_transform, is_train=False)
        dataloader = DataLoader(data, batch_size=self.batch_size, shuffle=False, num_workers=self.n_worker, follow_batch=['ligand_index'])
        self.test_prediction_df = self.test_df.iloc[data.valid_idx].reset_index(drop=True)
        return dataloader

if __name__ == "__main__":
    pass
    # from easydict import EasyDict
    # import yaml
    # import warnings

    # import warnings

    # warnings.filterwarnings("ignore")

    # df = pd.read_csv("/work/yufeng/2022/enzyme_specificity/data/data.csv", sep=',')
    # config = load_config("/work/yufeng/2022/enzyme_specificity/src/Configs/specificity.yml")
    # reaction = Reaction(df=df, st=0, ed=len(df.index), config=config)
    # dataloader = DataLoader(reaction, batch_size=2, shuffle=True, num_workers=20)
    # for _ in dataloader:
    #     print(_)
    #     break
    # print(len(reaction))
    # print(reaction[0])