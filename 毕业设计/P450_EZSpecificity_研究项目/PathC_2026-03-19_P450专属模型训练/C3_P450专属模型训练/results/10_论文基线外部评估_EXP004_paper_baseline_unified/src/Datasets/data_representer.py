import lmdb
import pickle
from torch_geometric.data import Data
import torch.utils.data as data
import os
import torch
import numpy as np

from Datasets.utils import preprocess_enzyme_feature, preprocess_reaction_feature, torchify_dict, load_tensor, load_pickle, get_paths, check_paths_exist

from Datasets.Structure import StructureDataset

def get_representer(config, df, transform=None, is_train=True):
    if config.data.representer == 'structure_sequence':
        return StructureSequence(config=config, df=df, transform=transform, is_train=is_train)
    else:
        raise ValueError("The representer is not supported.")


class Reaction(data.Dataset):
    def __init__(self, df, config=None) -> None:
        super().__init__()

        self.df = df
        self.config = config
        self.enzyme_dbs = None
        self.reaction_dbs = None
        self.grover_dbs = None
        self.morgan_dbs = None
        self.db_key_dict = None

        self.valid_idx = None
        self._get_valid_idx()
        
        if "grover" in self.config.data.features or "grover_mean" in self.config.data.features:
            assert check_paths_exist(self.config.data.grover_path)

        if "morgan" in self.config.data.features:
            assert check_paths_exist(self.config.data.morgan_path)

        assert check_paths_exist(self.config.data.enzyme_lmdb_path)
        assert check_paths_exist(self.config.data.reaction_lmdb_path)

    def _check_valid_idx(self, idx):
        reation_idx = self.df.loc[idx, 'Substrate Index']
        enzyme_idx = self.df.loc[idx, 'Enzyme Index']
        dataset_id = self.df.loc[idx, 'dataset_id']

        key = str(reation_idx).encode()
        if key not in self.db_key_dict[f"reaction_{dataset_id}"]:
            print(f"reaction_{dataset_id} {reation_idx} not in db")
            return False
        
        key = str(enzyme_idx).encode()
        if key not in self.db_key_dict[f"enzyme_{dataset_id}"]:
            print(f"enzyme_{dataset_id} {enzyme_idx} not in db")
            return False            
            
        if "grover" in self.config.data.atom_features or "grover_mean" in self.config.data.features:
            key = str(reation_idx).encode()
            if key not in self.db_key_dict[f"grover_{dataset_id}"]:
                print(f"grover_{dataset_id} {reation_idx} not in db")
                return False
        
        return True

    def _get_valid_idx(self):
        if self.enzyme_dbs is None or self.reaction_dbs is None:
            self._connect_db()
        cnt = 0
        self.valid_idx = []
        for idx in range(self.df.index.shape[0]):
            if self._check_valid_idx(idx):
                self.valid_idx.append(idx)
            else:
                cnt += 1
        print(f"Invalid reaction data: {cnt}")
    
    def _add_key_dict(self, name, db):
        if db is None:
            return

        with db.begin() as txn:
            keys = list(txn.cursor().iternext(values=False))
            self.db_key_dict[name] = {k: i for i, k in enumerate(keys)}
        return 

    def _build_key_dict(self):
        if ("grover" in self.config.data.atom_features or "grover_mean" in self.config.data.features):
            for index, db in enumerate(self.grover_dbs):
                self._add_key_dict(f"grover_{index}", db)
        for index, db in enumerate(self.enzyme_dbs):
            self._add_key_dict(f"enzyme_{index}", db)
        for index, db in enumerate(self.reaction_dbs):
            self._add_key_dict(f"reaction_{index}", db)
    
    def _connect_db(self):
        # assert self.db is None, 'A connection has already been opened.'
        if self.grover_dbs is None and ("grover" in self.config.data.atom_features or "grover_mean" in self.config.data.features):
            self.grover_dbs = [lmdb.open(
                path,
                map_size=10*(1024*1024*1024),   # 10GB (reduced from 600GB for Windows compatibility)
                create=False,
                subdir=False,
                readonly=True,
                lock=False,
                readahead=False,
                meminit=False,
            ) for path in get_paths(self.config.data.grover_path)]

        if self.morgan_dbs is None and "morgan" in self.config.data.features:
            self.morgan_dbs = [np.load(path).astype(np.float32) for path in get_paths(self.config.data.morgan_path)]

        if self.enzyme_dbs is None:
            self.enzyme_dbs = [lmdb.open(
                path,
                map_size=10*(1024*1024*1024),   # 10GB (reduced from 600GB for Windows compatibility)
                create=False,
                subdir=False,
                readonly=True,
                lock=False,
                readahead=False,
                meminit=False,
            ) for path in get_paths(self.config.data.enzyme_lmdb_path)]

        if self.reaction_dbs is None:
            self.reaction_dbs = [lmdb.open(
                path,
                map_size=10*(1024*1024*1024),   # 10GB (reduced from 600GB for Windows compatibility)
                create=False,
                subdir=False,
                readonly=True,
                lock=False,
                readahead=False,
                meminit=False,
            ) for path in get_paths(self.config.data.reaction_lmdb_path)]
        
        self.db_key_dict = {}
        self._build_key_dict()


    def __len__(self):
        return self.df.index.shape[0]

    def getitem_with_real_idx(self, idx):
        reation_idx = self.df.loc[idx, 'Substrate Index']
        enzyme_idx = self.df.loc[idx, 'Enzyme Index']
        dataset_id = self.df.loc[idx, 'dataset_id']
        tag = self.df.loc[idx, 'tag']

        # if int(tag) == 0:
        #     ecnumbers = self.df.loc[idx, 'ecnumber'].split(".")
        #     fake_ecnumber = self.df.loc[idx, 'fake_ecnumber'].split(".")
        #     id = 4
        #     for i in range(len(ecnumbers)):
        #         if ecnumbers[i] != fake_ecnumber[i]:
        #             id = i
        #             break
        #     tag = str(id)
        # else:
        #     tag = str(5)

        tag = str(tag)

        substrate_features = {}
        
        if self.enzyme_dbs is None or self.reaction_dbs is None:
            self._connect_db()

        with self.reaction_dbs[dataset_id].begin(write=False) as txn:
            key = str(reation_idx).encode()
            value = txn.get(key)
            reaction_data = pickle.loads(value)
        
        with self.enzyme_dbs[dataset_id].begin(write=False) as txn:
            key = str(enzyme_idx).encode()
            value = txn.get(key)
            enzyme_data = pickle.loads(value)

        if "grover" in self.config.data.atom_features or "grover_mean" in self.config.data.features:
            with self.grover_dbs[dataset_id].begin(write=False) as txn:
                key = str(reation_idx).encode()
                value = txn.get(key)
                grover_data = pickle.loads(value)
                substrate_features["grover_mean"] = torch.from_numpy(grover_data['total_embedding'][np.newaxis, :])
                reaction_data["grover"] = torch.from_numpy(grover_data['embedding'])

        if "morgan" in self.config.data.features:
            substrate_features["morgan"] = torch.from_numpy(self.morgan_dbs[dataset_id][reation_idx][np.newaxis, :])
        
        # padding
        reaction_data = preprocess_reaction_feature(reaction_data, self.config.data.max_substrate_length)
        enzyme_data = preprocess_enzyme_feature(enzyme_data, self.config.data.max_enzyme_length)

        data = ReactionData.from_protein_ligand_dicts(torchify_dict(reaction_data), torchify_dict(enzyme_data), label=torch.from_numpy(np.array(self.df.loc[idx, 'Label'], dtype=float)), tag=tag, **substrate_features)

        return data

    def __getitem__(self, idx):
        if self.valid_idx is None:
            self._get_valid_idx()
        return self.getitem_with_real_idx(self.valid_idx[idx])

class StructureSequence(data.Dataset):
    def __init__(self, df, config, transform, is_train) -> None:
        super().__init__()

        self.df = df
        self.config = config
        
        transform = transform

        self.sequence_db = Reaction(config=config, df=df)
        self.structure_db = StructureDataset(config=config, df=df, transform=transform, is_train=is_train)

        self.valid_idx = None
        self._get_valid_idx()

    def _get_valid_idx(self):
        valid_idx_sequence = self.sequence_db.valid_idx
        valid_idx_structure = self.structure_db.valid_idx
        print(len(valid_idx_sequence))
        print(len(valid_idx_structure))
        self.valid_idx = list(set(valid_idx_sequence).intersection(set(valid_idx_structure)))
        print(len(self.valid_idx))

    def __len__(self):
        return len(self.valid_idx)

    def update_valid_idx(self, valid_idx):
        self.valid_idx = valid_idx

    def __getitem__(self, idx):
        if self.valid_idx is None:
            self._get_valid_idx()

        idx = self.valid_idx[idx]

        sequence_data = self.sequence_db.getitem_with_real_idx(idx)
        structure_data = self.structure_db.getitem_with_real_idx(idx)

        data = StructureSequenceData.from_sequence_structure(sequence_data, structure_data)
        
        if data.y is not None and data.label is not None:
            assert data.y == data.label
        
        data.active_site = None
        data.y = None
        return data
    
class ReactionData(Data):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @staticmethod
    def from_protein_ligand_dicts(protein_dict=None, ligand_dict=None, **kwargs):
        instance = ReactionData(**kwargs)

        if protein_dict is not None:
            for key, item in protein_dict.items():
                instance[key] = item

        if ligand_dict is not None:
            for key, item in ligand_dict.items():
                instance[key] = item

        for key, item in kwargs.items():
            instance[key] = item
        return instance

    def __inc__(self, key, value, *args, **kwargs):
        if key == 'edge_index':
            return self['element'].size(0)
        else:
            return super().__inc__(key, value)

class StructureSequenceData(Data):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @staticmethod
    def from_sequence_structure(sequence=None, structure=None, **kwargs):
        instance = StructureSequenceData(**kwargs)

        if sequence is not None:
            for key, item in sequence.to_dict().items():
                instance[key] = item

        if structure is not None:
            for key, item in structure.to_dict().items():
                instance[key] = item

        for key, item in kwargs.items():
            instance[key] = item
        return instance

    def __inc__(self, key, value, *args, **kwargs):
        if key == 'edge_index':
            return self['element'].size(0)
        elif key == 'ligand_index':
            return 0
        elif key == 'complex_edge_index':
            return self['ligand_x'].size(0)
        else:
            return super().__inc__(key, value, *args, **kwargs)