import numpy as np
import torch
from tqdm import tqdm
from torch_scatter import scatter
import rdkit
import random
from rdkit import Chem
from rdkit.Chem.rdchem import BondType, HybridizationType
import pandas as pd
import pickle, os

from Datasets.const import restype_3to1, restype_name_to_atom14_names, letter_to_num

def get_paths(data):
    if type(data) == str:
        return [data]
    else:
        return data

def check_paths_exist(data):
    paths = get_paths(data)
    for path in paths:
        if not os.path.exists(path):
            return False
    return True

def read_datasets(data):
    result = []
    dataset_id = []

    for index, path in enumerate(get_paths(data)):
        df = pd.read_csv(path, sep=',')
        result.append(df)
        dataset_id.extend([index] * len(result[-1].index))
    
    df = pd.concat(result, ignore_index=True).reset_index(drop=True)
    df['dataset_id'] = dataset_id
    df['tag'] = dataset_id
    return df

def get_atom_coord(residue):
    atoms = restype_name_to_atom14_names[residue.get_resname()]
    coords = []
    mask = []
    for atom in atoms:
        if len(atom) == 0 or (not residue.has_id(atom)):
            coords.append(np.array([0, 0, 0]))
            mask.append(0)
        else:
            coords.append(residue[atom].get_coord())
            mask.append(1)
    return np.array(coords), np.array(mask)

def get_res_coord(residues):
    coords = []
    masks = []
    for residue in residues:
        coord, mask = get_atom_coord(residue)
        coords.append(coord)
        masks.append(mask)
    return np.array(coords)


def _normalize(tensor, dim=-1):
    '''
    Normalizes a `torch.Tensor` along dimension `dim` without `nan`s.
    '''
    return torch.nan_to_num(
        torch.div(tensor, torch.norm(tensor, dim=dim, keepdim=True)))


def _rbf(D, D_min=0., D_max=20., D_count=16, device='cpu'):
    '''
    From https://github.com/jingraham/neurips19-graph-protein-design
    
    Returns an RBF embedding of `torch.Tensor` `D` along a new axis=-1.
    That is, if `D` has shape [...dims], then the returned tensor will have
    shape [...dims, D_count].
    '''
    D_mu = torch.linspace(D_min, D_max, D_count, device=device)
    D_mu = D_mu.view([1, -1])
    D_sigma = (D_max - D_min) / D_count
    D_expand = torch.unsqueeze(D, -1)

    RBF = torch.exp(-((D_expand - D_mu) / D_sigma) ** 2)
    return RBF

def check_residue(residue):
    if residue.get_id()[0] != ' ':
        return False
    if not residue.get_resname() in restype_3to1:
        return False
    return True

def get_residue(chain):
    residues = []
    seqs = ""
    for residue in chain.get_residues():
        if check_residue(residue):
            residues.append(residue)
            seqs+=restype_3to1[residue.get_resname()]
    return residues, seqs

def change_inchi_to_smile(inchi):
    mol = rdkit.Chem.inchi.MolFromInchi(inchi)
    if mol is None:
        return "-"
    else:
        return Chem.MolToSmiles(mol)

def add_dataset_id(df, id):
    df['dataset_id'] = np.array([id] * len(df.index))
    return df

def generate_negative_sample(df, config=None, same_digits=None, num_negative_enzyme=None, random_negative_digit=None, has_positive_sample=True):
    if config is None:
        assert same_digits is not None and num_negative_enzyme is not None
    else:
        same_digits = config.data.sampling.same_digits
        num_negative_enzyme = config.data.sampling.num_negative_enzyme

    print("Start generating negative samples!")
    # calculate dict of (reaction -> [enzymes])
    reaction_dict = {}
    # enzyme_dict of (enzyme -> [ecnumbers])
    enzyme_dict = {}
    # ecnumber_dict of (ecnumber -> [enzymes])
    ecnumber_dict = {}
    # negative_dict of ([reaction, enzyme] -> true,false)
    negative_dict = {}

    for reaction, enzyme, ecnumber in zip(df['reaction'], df['enzyme'], df['ecnumber']):
        if reaction not in reaction_dict:
            reaction_dict[reaction] = []
        reaction_dict[reaction].append((enzyme, ecnumber))
        
        if enzyme in enzyme_dict:
            enzyme_dict[enzyme].append(ecnumber)
        else:
            enzyme_dict[enzyme] = [ecnumber]
        
        for i in range(0, 5):
            chuncked_ecnumber = ".".join(ecnumber.split(".")[:i])
            if chuncked_ecnumber not in ecnumber_dict:
                ecnumber_dict[chuncked_ecnumber] = []
            ecnumber_dict[chuncked_ecnumber].append(enzyme)

    for key in reaction_dict:
        reaction_dict[key] = {enzyme: True for enzyme, ecnumber in list(set(reaction_dict[key]))}

    # create negative data and positive data
    data = {
        'reaction': [],
        'enzyme': [],
        'ecnumber': [],
        'fake_ecnumber': [],
        'label': []
    }
    for reaction, enzyme, ecnumber in tqdm(zip(df['reaction'], df['enzyme'], df['ecnumber']), total=df['enzyme'].shape[0]):
        # positive sample
        if has_positive_sample:
            data['reaction'].append(reaction)
            data['enzyme'].append(enzyme)
            data['ecnumber'].append(ecnumber)
            data['fake_ecnumber'].append(ecnumber)
            data['label'].append(1)
        
        # negative sample
        
        # determine the difficulty of negative sample
        if random_negative_digit is None:
            negative_same_digit = same_digits
        else:
            negative_same_digit = random.randint(0, same_digits+1)

        current_n_negative_sample = 0
        sample_enzymes = ecnumber_dict[".".join(ecnumber.split(".")[:negative_same_digit])]
        if len(sample_enzymes) > num_negative_enzyme * 10:
            sample_enzymes = random.sample(sample_enzymes, num_negative_enzyme * 10)
        random.shuffle(sample_enzymes)

        for negative_enzyme in sample_enzymes:
            if current_n_negative_sample < num_negative_enzyme:
                if negative_enzyme not in reaction_dict[reaction] and (negative_enzyme, reaction) not in negative_dict:
                    data['reaction'].append(reaction)
                    data['enzyme'].append(negative_enzyme)
                    data['fake_ecnumber'].append(enzyme_dict[negative_enzyme][0])
                    data['ecnumber'].append(ecnumber)
                    data['label'].append(0)
                    negative_dict[(negative_enzyme, reaction)] = True
                    current_n_negative_sample += 1
    data = pd.DataFrame(data)
    return data

def download_uniprot_file():
    import re
    import requests
    from requests.adapters import HTTPAdapter, Retry

    re_next_link = re.compile(r'<(.+)>; rel="next"')
    retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))

    def get_next_link(headers):
        if "Link" in headers:
            match = re_next_link.match(headers["Link"])
            if match:
                return match.group(1)

    def get_batch(batch_url):
        while batch_url:
            response = session.get(batch_url)
            response.raise_for_status()
            total = response.headers["x-total-results"]
            yield response, total
            batch_url = get_next_link(response.headers)


    url = 'https://rest.uniprot.org/uniprotkb/search?compressed=false&fields=accession%2Creviewed%2Cid%2Clength%2Cft_act_site%2Cft_binding&format=tsv&query=%28%28ftlen_act_site%3A%5B%2A%20TO%2010000%5D%29%29&size=500'
    progress = 0
    with open(f"{root_dir}/data/brenda/raw_data/uniprot/importantsites.tsv", 'w') as f:
        for batch, total in get_batch(url):
            lines = batch.text.splitlines()
            if not progress:
                print(lines[0], file=f)
            for line in lines[1:]:
                print(line, file=f)
            progress += len(lines[1:])
            print(f'{progress} / {total}', end='\r')

def get_ligand_atom_features(rdmol):
    num_atoms = rdmol.GetNumAtoms()
    atomic_number = []
    aromatic = []
    sp, sp2, sp3 = [], [], []
    degree = []
    index_to_newindex_dict = {}
    for index, atom in enumerate(rdmol.GetAtoms()):
        atom_idx = atom.GetIdx()
        index_to_newindex_dict[atom_idx] = index
        atom = rdmol.GetAtomWithIdx(atom_idx)
        atomic_number.append(atom.GetAtomicNum())
        aromatic.append(1 if atom.GetIsAromatic() else 0)
        hybridization = atom.GetHybridization()
        sp.append(1 if hybridization == HybridizationType.SP else 0)
        sp2.append(1 if hybridization == HybridizationType.SP2 else 0)
        sp3.append(1 if hybridization == HybridizationType.SP3 else 0)
        degree.append(atom.GetDegree())
    node_type = torch.tensor(atomic_number, dtype=torch.long)

    row, col = [], []
    for bond in rdmol.GetBonds():
        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        start = index_to_newindex_dict[start]
        end = index_to_newindex_dict[end]
        row += [start, end]
        col += [end, start]
    row = torch.tensor(row, dtype=torch.long)
    col = torch.tensor(col, dtype=torch.long)
    hs = (node_type == 1).to(torch.float)
    num_hs = scatter(hs[row], col, dim_size=num_atoms).numpy()
    feat_mat = np.array([atomic_number, aromatic, degree, num_hs, sp, sp2, sp3], dtype=np.compat.long).transpose()
    return feat_mat

def parse_smile(string):
    mol = Chem.MolFromSmiles(string)

    tmp_mol = Chem.MolFromSmiles(string)
    Chem.SanitizeMol(tmp_mol)
    for atom in tmp_mol.GetAtoms():
        atom.SetAtomMapNum(0)
    canonical_mol = Chem.MolFromSmiles(Chem.MolToSmiles(tmp_mol))
    mapping_idx = [int(x) for x in tmp_mol.GetProp("_smilesAtomOutputOrder")[1:-1].split(",") if x.strip()]

    for atom in canonical_mol.GetAtoms():
        atom.SetAtomMapNum(mol.GetAtomWithIdx(mapping_idx[atom.GetIdx()]).GetAtomMapNum())
    mol = canonical_mol

    feat_mat = get_ligand_atom_features(mol)

    ptable = Chem.GetPeriodicTable()

    num_atoms = mol.GetNumAtoms()
    num_bonds = mol.GetNumBonds()

    element = []
    index_to_newindex_dict = {}
    newindex_to_index_dict = {}
    for index in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(index)
        index_to_newindex_dict[atom.GetAtomMapNum()] = atom.GetIdx()
        newindex_to_index_dict[atom.GetIdx()] = atom.GetAtomMapNum()
        atomic_number = atom.GetAtomicNum()
        element.append(atomic_number)
    element = np.array(element, dtype=np.compat.long)
    row, col, edge_type = [], [], []
    BOND_TYPES = {t: i for i, t in enumerate(BondType.names.values())}
    for bond in mol.GetBonds():
        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        row += [start, end]
        col += [end, start]
        edge_type += 2 * [BOND_TYPES[bond.GetBondType()]]
    edge_index = np.array([row, col], dtype=np.compat.long)
    edge_type = np.array(edge_type, dtype=np.compat.long)
    perm = (edge_index[0] * num_atoms + edge_index[1]).argsort()
    edge_index = edge_index[:, perm]
    edge_type = edge_type[perm]

    data = {
        'element': element,
        'edge_index': edge_index,
        'edge_type': edge_type,
        'atom_feature': feat_mat
    }
    return data, index_to_newindex_dict, newindex_to_index_dict, mol

def check_smile_equal(smilea, smileb):
    # smilea will contain index and smileb is clean smile
    a = Chem.MolFromSmiles(smilea)
    for atom in a.GetAtoms():
        atom.SetAtomMapNum(0)
    a = Chem.CanonSmiles(Chem.MolToSmiles(a))
    b = Chem.CanonSmiles(smileb)
    return a == b

def load_pickle(file_name):
    with open(file_name, 'rb') as f:
        return pickle.load(f)

def load_tensor(file_name, dtype):
    return [dtype(d) for d in np.load(file_name + '.npy', allow_pickle=True)]

def get_neighbor_list(data, idx, dict):
    neighbor_list = []
    for i in range(data["edge_index"].shape[1]):
        if data["edge_index"][0][i] == idx and data["edge_index"][1][i] in dict:
            neighbor_list.append((data["edge_index"][1][i], data["edge_type"][i]))
    return neighbor_list

def convert_protein_sequence_to_number(protein_sequence):
    protein_number = []
    for i in protein_sequence:
        protein_number.append(letter_to_num[i])
    return np.array(protein_number)

def torchify_dict(data):
    output = {}
    for k, v in data.items():
        if isinstance(v, np.ndarray):
            output[k] = torch.from_numpy(v)
        else:
            output[k] = v
    return output

def deredundant(values):
    # avoid using set(list) because it will change the order
    results = []
    value_dict = {}
    for value in values:
        if value not in value_dict:
            value_dict[value] = 1
            results.append(value)
    return np.array(results)


def preprocess_reaction_feature(data, max_reaction_length):
    data['reaction_padding_mask'] = np.zeros((max_reaction_length), dtype=bool)[None, ...]
    data['reaction_padding_mask'][0, data['element'].shape[0]:] = True
    if "grover" in data:
        data['grover'] = np.pad(data['grover'], ((0, max_reaction_length - data['grover'].shape[0]), (0, 0)), 'constant', constant_values=0)
    return data

def preprocess_enzyme_feature(data, max_enzyme_length):
    original_length = data['embedding'].shape[0]

    data['enzyme_padding_mask'] = np.zeros((max_enzyme_length), dtype=bool)[None, ...]
    data['enzyme_padding_mask'][0, data['embedding'].shape[0]:] = True
    data['embedding'] = np.pad(data['embedding'], ((0, max_enzyme_length - original_length), (0, 0)), 'constant', constant_values=0)
    return data
    
if __name__ == "__main__":
    pass