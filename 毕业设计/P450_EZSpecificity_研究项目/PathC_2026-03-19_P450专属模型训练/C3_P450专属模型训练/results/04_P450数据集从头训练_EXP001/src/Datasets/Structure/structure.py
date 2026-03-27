import os
import pickle
import lmdb
import torch
from torch.utils.data import Dataset
from tqdm.auto import tqdm
from easydict import EasyDict
import yaml
import sys
import ray

root_dir = "/projects/bbto/suyufeng/enzyme_specificity"
sys.path.append(f"{root_dir}/src")

from Datasets.utils import get_paths, check_paths_exist


from Datasets.Structure.protein_ligand import PDBProtein, parse_sdf_file_mol, get_zero_protein_feature
from Datasets.Structure.utils import StructureComplexData, torchify_dict

# @ray.remote
def _single_str_process(item):
    (index, config, reaction, protein, label, structure_index) = item
    print("Start processing {}".format(structure_index))
    if not os.path.exists(os.path.join(config.data.pdb_dir, f'{structure_index}.pdb')) or not os.path.exists(os.path.join(config.data.ligand_dir, f'{structure_index}.sdf')):
        return (1, None)
    try:
        pocket_fn = os.path.join(config.data.pdb_dir, f'{structure_index}.pdb')
        ligand_fn = os.path.join(config.data.ligand_dir, f'{structure_index}.sdf')
        pocket_dict = PDBProtein(pocket_fn).to_dict_atom()
        ligand_dict = parse_sdf_file_mol(ligand_fn, heavy_only=True)
        data = StructureComplexData.from_protein_ligand_dicts(
            protein_dict=torchify_dict(pocket_dict),
            ligand_dict=torchify_dict(ligand_dict),
        )
        assert data.protein_pos.size(0) > 0
        data.protein_filename = protein
        data.ligand_filename = reaction
        data.y = torch.tensor(float(label))
        return (0, (index, structure_index, data))
    except:
        return (2, None)

def str_process(config, df):
    processed_path = config.data.structure_processed_path
    db = lmdb.open(
        processed_path,
        map_size=10*(1024*1024*1024),   # 10GB
        create=True,
        subdir=False,
        readonly=False,  # Writable
    )

    # index = parse_pdbbind_index_file(self.index_path)

    parameters = []

    num_skipped_exist = 0
    num_skipped_parse = 0

    for index, (reaction, protein, label, structure_index) in tqdm(enumerate(zip(df['Substrate Index'].values, df['Enzyme Index'].values, df["Label"].values, df['Dock Index'].values))):
        structure_index = int(structure_index)
        if not os.path.exists(os.path.join(config.data.pdb_dir, f'{structure_index}.pdb')) or not os.path.exists(os.path.join(config.data.ligand_dir, f'{structure_index}.sdf')):
            continue
        parameters.append((index, config, reaction, protein, label, structure_index))


    from multiprocessing.pool import Pool


    # with db.begin(write=True, buffers=True) as txn:
    #     with Pool(60) as pool:
    #         for index, (cnt, item) in enumerate(pool.imap_unordered(_single_str_process, parameters)):
    #             print(f"Processing str {index}")
    #             if cnt == 0:
    #                 (index, structure_index, data) = item
    #                 txn.put(
    #                     key=str(structure_index).encode(),
    #                     value=pickle.dumps(data)
    #                 )
    #             else:
    #                 if cnt == 1:
    #                     num_skipped_exist += 1
    #                     print(f"non-existed skipped #{num_skipped_exist}: {structure_index}")
    #                 else:
    #                     num_skipped_parse += 1
    #                     print(f"parse failed skipped #{num_skipped_parse}: {structure_index}")

    for index, parameter in enumerate(parameters):
        print(f"Processing str {index}")
        (cnt, item) = _single_str_process(parameter)
        if cnt == 0:
            (index, structure_index, data) = item
            with db.begin(write=True, buffers=True) as txn:
                txn.put(
                    key=str(structure_index).encode(),
                    value=pickle.dumps(data)
                )
        else:
            if cnt == 1:
                num_skipped_exist += 1
                print(f"non-existed skipped #{num_skipped_exist}: {structure_index}")
            else:
                num_skipped_parse += 1
                print(f"parse failed skipped #{num_skipped_parse}: {structure_index}")

    print('num_skipped: ', num_skipped_exist, num_skipped_parse)

# @ray.remote(num_cpus=1)
def _single_seq_process(item):
    from rdkit import Chem
    from rdkit.Chem import AllChem

    try:
        index, smile = item
        mol = Chem.MolFromSmiles(smile)
        mol = Chem.AddHs(mol)

        AllChem.EmbedMolecule(mol,AllChem.ETKDG())
        mol = Chem.RemoveHs(mol)
        
        pocket_dict = get_zero_protein_feature()
        ligand_dict = parse_sdf_file_mol(None, mol, heavy_only=True)
        data = StructureComplexData.from_protein_ligand_dicts(
            protein_dict=torchify_dict(pocket_dict),
            ligand_dict=torchify_dict(ligand_dict),
        )
        data.protein_filename = None
        data.ligand_filename = index
        return index, data
    except:
        return index, None

def seq_process(config, df):
    processed_path = config.data.sequence_processed_path
    db = lmdb.open(
        processed_path,
        map_size=10*(1024*1024*1024),   # 10GB
        create=True,
        subdir=False,
        readonly=False,  # Writable
    )
    parameters = []
    for index, ligand in tqdm(enumerate(df['Substrate_SMILES'].values)):
        parameters.append((index, ligand))

    from multiprocessing.pool import Pool


    # with db.begin(write=True, buffers=True) as txn:
    #     with Pool(60) as pool:
    #         for index, (id, data) in enumerate(pool.imap_unordered(_single_seq_process, parameters)):
    #             print(f"Processing seq {index}")
    #             if data is not None:
    #                 txn.put(
    #                     key=str(id).encode(),
    #                     value=pickle.dumps(data)
    #                 )

    with db.begin(write=True, buffers=True) as txn:
        for index,  parameter in enumerate(parameters):
            print(f"Processing seq {index}")
            (id, data) = _single_seq_process(parameter)
            if data is not None:
                txn.put(
                    key=str(id).encode(),
                    value=pickle.dumps(data)
                )


class StructureDataset(Dataset):

    def __init__(self, config, df, transform=None, is_train=False):
        super().__init__()
        self.dbs = None
        self.df = df
        self.config = config

        self.complex_dbs = None
        self.ligand_dbs = None
        self.valid_idx = None
        self.high_quality_id_dicts = None
        self.db_key_dict = None
        self.idx_valid_dict = {}

        if is_train and self.config.data.full_data:
            self.fake_sequence_ratio = self.config.data.fake_sequence_ratio
        else:
            self.fake_sequence_ratio = 0.0

        self.transform = transform
        if check_paths_exist(config.data.structure_processed_path):
            assert "processed_path is not initialized"

        self._get_valid_idx()

    def _get_high_quality_id_dicts(self):
        self.high_quality_id_dicts = []
        
        for path in get_paths(self.config.data.high_quality_id_path):
            try:
                fin = open(path, "r")
                self.high_quality_id_dicts.append({int(line.strip()):True for line in fin})
                fin.close()
            except:
                self.high_quality_id_dicts.append(None)
    
    def check_high_quality(self, dataset_id, index):
        if self.high_quality_id_dicts is None:
            self._get_high_quality_id_dicts()
        return ((self.high_quality_id_dicts[dataset_id] is None) or (index in self.high_quality_id_dicts[dataset_id]))

    def _check_valid_idx(self, idx):
        structure_index = int(self.df['Dock Index'].values[idx])
        dataset_id = self.df['dataset_id'].values[idx]
        reaction = self.df['Substrate Index'].values[idx]

        # if self.check_high_quality(dataset_id, structure_index):
            # return 
        flag_complex = str(structure_index).encode() in self.db_key_dict[f"complex_{dataset_id}"]

        if not self.config.data.full_data:
            self.idx_valid_dict[idx] = (flag_complex, False)
            return flag_complex & self.check_high_quality(dataset_id, structure_index)
        
        flag_ligand = str(reaction).encode() in self.db_key_dict[f"ligand_{dataset_id}"]
        self.idx_valid_dict[idx] = (flag_complex, flag_ligand)

        return flag_complex | flag_ligand

    def _get_valid_idx(self):
        import numpy as np
        if self.complex_dbs is None or self.ligand_dbs is None:
            self._connect_db()

        cnt = 0
        self.valid_idx = []

        # unavailable_structure_indexs = []
        for idx in range(self.df.index.shape[0]):
            if self._check_valid_idx(idx):
                self.valid_idx.append(idx)
            else:
                cnt += 1
                # unavailable_structure_indexs.append(int(self.df['structure_index'].values[idx]))   
        
        # unavailable_structure_indexs = list(set(unavailable_structure_indexs))
        # np.savetxt("unavailable_structure_indexs.txt", unavailable_structure_indexs)
        print(f"#Invalid structures: {cnt}")

    def _add_key_dict(self, name, db):
        if db is None:
            return

        with db.begin() as txn:
            keys = list(txn.cursor().iternext(values=False))
            self.db_key_dict[name] = {k: i for i, k in enumerate(keys)}
        return 

    def _connect_db(self):
        assert self.complex_dbs is None, 'A connection has already been opened.'
        self.db_key_dict = {}
        self.complex_dbs = [lmdb.open(
            path,
            map_size=10*(1024*1024*1024),   # 10GB
            create=False,
            subdir=False,
            readonly=True,
            lock=False,
            readahead=False,
            meminit=False,
        ) for path in get_paths(self.config.data.structure_processed_path)]
        self.ligand_dbs = []
        for path in get_paths(self.config.data.sequence_processed_path):
            if os.path.exists(path):
                self.ligand_dbs.append(lmdb.open(
                    path,
                    map_size=10*(1024*1024*1024),   # 10GB
                    create=False,
                    subdir=False,
                    readonly=True,
                    lock=False,
                    readahead=False,
                    meminit=False,
                ))
            else:
                self.ligand_dbs.append(None)
        for i, db in enumerate(self.complex_dbs):
            self._add_key_dict(f"complex_{i}", db)
        for i, db in enumerate(self.ligand_dbs):
            self._add_key_dict(f"ligand_{i}", db)

    def __len__(self):
        if self.valid_idx is None:
            self._get_valid_idx()
        return len(self.valid_idx)
    
    def getitem_with_real_idx(self, idx):
        import random

        if self.complex_dbs is None or self.ligand_dbs is None:
            self._connect_db()

        structure_index = int(self.df['Dock Index'].values[idx])
        dataset_id = self.df['dataset_id'].values[idx]
        reaction = self.df['Substrate Index'].values[idx]
        
        flag_no_structure = False
        try:
            if self.config.data.no_structure == True:
                flag_no_structure = True
        except Exception as e:
            pass

        complex_random_value = random.random() if flag_no_structure else 1.0

        if ((self.check_high_quality(dataset_id, structure_index) and complex_random_value > self.fake_sequence_ratio) or (not self.idx_valid_dict[idx][1])) and self.idx_valid_dict[idx][0]:
            with self.complex_dbs[dataset_id].begin(write=False) as txn:
                data = txn.get(str(structure_index).encode())
                if data is None:
                    print(f"None: {structure_index}")
                    print(self.idx_valid_dict[idx])
                data = pickle.loads(data)
                data.str_tag = 'complex'
                data.sample_weight = torch.tensor(float(self.config.data.sample_weight[0]))
                data.mask_use_complex = torch.tensor([1.0] * (data.ligand_pos.size(0) + data.protein_pos.size(0)))
                data.mask_use_ligand = torch.tensor([0.0] * (data.ligand_pos.size(0) + data.protein_pos.size(0)))
        else:
            with self.ligand_dbs[dataset_id].begin(write=False) as txn:
                data = txn.get(str(reaction).encode())
                data = pickle.loads(data)
                data.str_tag = 'ligand'
                data.sample_weight = torch.tensor(float(self.config.data.sample_weight[1]))
                data.mask_use_complex = torch.tensor([0.0] * (data.ligand_pos.size(0) + data.protein_pos.size(0)))
                data.mask_use_ligand = torch.tensor([1.0] * (data.ligand_pos.size(0) + data.protein_pos.size(0)))
        
        data.id = idx
        if self.transform is not None:
            data = self.transform(data)
        return data

    def __getitem__(self, idx):
        if self.valid_idx is None:
            self._get_valid_idx()
        idx = self.valid_idx[idx]
        return self.getitem_with_real_idx(idx)
        
    # def __cat_dim__(self, key, value, *args, **kwargs):
    #     super().__cat_dim__(key, value, *args, **kwargs)
  
    # def __inc__(self, key, value, *args, **kwargs):
    #     super().__inc__(key, value, *args, **kwargs)
