import sys
data_root_dir = "/scratch/bbto/suyufeng/tmp/enzyme_specificity"
src_root_dir = "/projects/bbto/suyufeng/enzyme_specificity"
sys.path.append(f"{src_root_dir}/src")
import numpy as np
import torch
from tqdm import tqdm
from rdkit import Chem, RDLogger
import ray
import pandas as pd
import lmdb
import pickle
import esm

from Datasets.utils import parse_smile, get_neighbor_list, check_smile_equal, convert_protein_sequence_to_number

# @ray.remote(num_cpus=1, max_calls=1)
def get_reaction_feature_single(item):
    # disable warning
    from transformers import logging
    logging.set_verbosity_error()
    RDLogger.DisableLog('rdApp.*')  

    # parameters
    # print(item)
    full_reaction, substrate, right, substrate_only = item
    substrate_only = True
    smiles = right.split(".")
    ratio_threshold=0.8
    miss_atom_threshold=2
    confident_threshold=0
    
    def generate_substrate_only_data(substrate_data):
        data = {
            'element': substrate_data['element'],
            'edge_index': substrate_data['edge_index'],
            'edge_type': substrate_data['edge_type'],
            'atom_feature': substrate_data['atom_feature'],
            'reaction_attention_label': np.zeros((substrate_data['element'].shape[0])),
            'num_nodes': np.array(substrate_data['element'].shape[0])
        }
        return data

    # delete molecular that doesn't contain any atom except h  (rxnmapper can't handle this)
    smile_temp = []
    for smile in smiles:
        mol = Chem.MolFromSmiles(smile)
        if mol is not None:
            mol = Chem.rdmolops.RemoveAllHs(mol)
            if mol.GetNumAtoms() > 0:
                smile_temp.append(smile)
    smiles = smile_temp

    # check valid of reaction
    # 1. substrate doesn't contain any atom except h
    mol = Chem.MolFromSmiles(substrate)
    if mol is not None:
        mol = Chem.rdmolops.RemoveAllHs(mol)
        if mol.GetNumAtoms() == 0:
            return full_reaction, substrate, "Substrate is H+", None
    else:
        return full_reaction, substrate, f"Substrate {substrate} is not parseable", None
    
    # product unknwon
    if substrate_only == True:
        return full_reaction, substrate, f"Product unknown", generate_substrate_only_data(parse_smile(substrate)[0])

    # # 2. product set empty
    # if len(smiles) == 0:
    #     return full_reaction, substrate, 'No valid product', generate_substrate_only_data(parse_smile(substrate)[0])
    
    # # get atom mapping via rxnmapper and create reaction attention label
    # try:
    #     from rxnmapper import RXNMapper
    #     rxn_mapper = RXNMapper()
    #     results = rxn_mapper.get_attention_guided_atom_maps([full_reaction])[0]
    # except Exception as e: 
    #     return full_reaction, substrate, 'Error in rxnmapper: %s' % e, generate_substrate_only_data(parse_smile(substrate)[0])

    # if results['confidence'] < confident_threshold:
    #     return full_reaction, substrate, 'Atom mapping confidence is too low: %f' % results['confidence'], generate_substrate_only_data(parse_smile(substrate)[0])
    
    # left, right = results['mapped_rxn'].split('>>')
    # reaction_attention_label = None

    # for substrate_smile in left.split('.'):
    #     # Check reactant equals substrate
    #     if check_smile_equal(substrate_smile, substrate) == False:
    #         continue
    #     # Two dicts: m2s and s2m
    #     # m2s: rxn_id -> canonicalized id
    #     # s2m: canonicalized id -> rxn_id

    #     substrate_data, substrate_m2s, substrate_s2m, _ = parse_smile(substrate_smile)
    #     substrate_matched_index_dict = {}
    #     if reaction_attention_label is None:
    #         reaction_attention_label = np.zeros((substrate_data['element'].shape[0]))

    #     for product_smile in right.split("."):
    #         product_data, product_m2s, product_s2m, product_mol = parse_smile(product_smile)
    #         for atom in product_mol.GetAtoms():
    #             product_m_id = atom.GetAtomMapNum()
    #             if product_m_id in substrate_m2s and product_m_id != 0:
                    
    #                 substrate_s_id = substrate_m2s[product_m_id]
    #                 substrate_matched_index_dict[substrate_s_id] = True
    #                 product_s_id = product_m2s[product_m_id]
                    
    #                 if not np.equal(substrate_data['atom_feature'][substrate_s_id], product_data['atom_feature'][product_s_id]).any():
    #                     reaction_attention_label[substrate_s_id] = 1
    #                 else:
    #                     substarte_neighbor = get_neighbor_list(substrate_data, substrate_s_id, substrate_s2m)
    #                     product_neighbor = get_neighbor_list(product_data, product_s_id, product_s2m)
    #                     if set(substarte_neighbor) != set(product_neighbor):
    #                         reaction_attention_label[substrate_s_id] = 1

    # if reaction_attention_label is None:
    #     return full_reaction, substrate, 'No matched substrate', generate_substrate_only_data(parse_smile(substrate)[0])

    # for i in range(len(reaction_attention_label)):
    #     if substrate_data['element'][i] == 1:
    #         reaction_attention_label[i] = 0

    # if len(substrate_matched_index_dict) / substrate_data['element'].shape[0] < ratio_threshold and substrate_data['element'].shape[0] - len(substrate_matched_index_dict) > miss_atom_threshold:
    #     return full_reaction, substrate, "Atom mapping ratio is too low: %f" % (len(substrate_matched_index_dict) / substrate_data['element'].shape[0]), generate_substrate_only_data(parse_smile(substrate)[0])

    # # if np.sum(reaction_attention_label) == 0:
    # #     return full_reaction, substrate, "No atom changed in reaction", generate_substrate_only_data(parse_smile(substrate)[0])

    # data = {
    #     'element': substrate_data['element'],
    #     'edge_index': substrate_data['edge_index'],
    #     'edge_type': substrate_data['edge_type'],
    #     'atom_feature': substrate_data['atom_feature'],
    #     'reaction_attention_label': reaction_attention_label,
    #     'num_nodes': np.array(substrate_data['element'].shape[0])
    # }

    # return full_reaction, substrate, None, data

def to_iterator(obj_ids):
    while obj_ids:
        done, obj_ids = ray.wait(obj_ids)
        yield ray.get(done[0])

def get_reaction_parameters_new_brenda(df_path, allow_no_match=False):
    df = pd.read_csv(df_path, sep=',').dropna(subset=['substrate', 'reaction'])
    parameters_dict = {}
    for substrate, reaction in zip(df['substrate'], df['reaction']):
        left = reaction.split(">>")[0]
        right = reaction.split(">>")[1]
        if reaction not in parameters_dict and len(substrate) < 275:
            parameters_dict[reaction] = (len(parameters_dict) - 1, (reaction, substrate, right, allow_no_match))
    return parameters_dict

def get_reaction_parameters_brenda(df_path, allow_no_match=False):
    df = pd.read_csv(df_path, sep=',').dropna(subset=['substrate', 'left', 'right'])
    parameters_dict = {}
    for left, substrate, right in zip(df['left'], df['substrate'], df['right']):
        reaction = left + '>>' + right
        if reaction not in parameters_dict and len(substrate) < 275:
            parameters_dict[reaction] = (len(parameters_dict) - 1, (reaction, substrate, right, allow_no_match))
    return parameters_dict

def get_reaction_parameters_halogenase(df_path):
    df = pd.read_csv(df_path, sep=',').dropna(subset=['Substrate0', 'ReactionInLine'])
    parameters_dict = {}
    for substrate, reaction in zip(df['Substrate0'], df['ReactionInLine']):
        if reaction not in parameters_dict:
            parameters_dict[reaction] = (len(parameters_dict) - 1, (reaction, substrate, reaction.strip().split(">>")[1], False))
    return parameters_dict

def get_reaction_parameters_new_brenda(df_path):
    df = pd.read_csv(df_path, sep=',').dropna(subset=['Substrate SMILES'])
    parameters_dict = {}
    for substrate in df['Substrate SMILES']:
        if substrate not in parameters_dict:
            parameters_dict[substrate] = (len(parameters_dict) - 1, (substrate + ">>" + substrate, substrate, substrate, False))
    return parameters_dict

def get_reaction_parameters_resume(df_path):
    df = pd.read_csv(df_path, sep=',')
    parameters_dict = {}
    for index, (substrate, reaction) in enumerate(zip(df['substrates'], df['reactions'])):
        parameters_dict[reaction] = (index, (reaction, substrate, reaction.strip().split(">>")[1], False))
    return parameters_dict

def get_reaction_parameters_halogenase_evaluation(df_path):
    df = pd.read_csv(df_path, sep=',').dropna(subset=['Canonical_SMILES'])
    parameters_dict = {}
    for substrate in df['Canonical_SMILES']:
        if substrate not in parameters_dict:
            parameters_dict[substrate] = (len(parameters_dict) - 1, (substrate + ">>" + substrate, substrate, substrate, True))
    return parameters_dict

def get_reaction_parameters_small_family(df_path):
    df = pd.read_csv(df_path, sep=',').dropna(subset=['SUBSTRATES'])
    parameters_dict = {}
    for substrate in df['SUBSTRATES']:
        if substrate not in parameters_dict:
            parameters_dict[substrate] = (len(parameters_dict) - 1, (substrate + ">>" + substrate, substrate, substrate, True))
    return parameters_dict

def get_reaction_parameters_df(df_path, substrate_only=False):
    df = pd.read_csv(df_path, sep=',')
    parameters_dict = {}
    for substrate in df['substrates']:
        if substrate not in parameters_dict:
            parameters_dict[substrate] = (len(parameters_dict), (substrate + ">>" + substrate, substrate, substrate, substrate_only))
    return parameters_dict

def get_reaction_parameters_resume_evaluation(df_path):
    df = pd.read_csv(df_path, sep=',')
    parameters_dict = {}
    for index, (substrate, reaction) in enumerate(zip(df['substrates'], df['reactions'])):
        parameters_dict[reaction] = (index, (reaction, substrate, reaction.strip().split(">>")[1], True))
    return parameters_dict

def get_reaction_feature(parameters_dict, save_lmdb_path, save_df_path, allow_no_match=False, resume=False):
    # ray.init(num_cpus=40, num_gpus=1)
    reactions = []
    substrates = []
    from multiprocessing.pool import Pool

    num_skipped = 0

    if resume:
        reaction_id_dict = {
            reaction: index for index, (reaction, substrate,_, _) in parameters_dict.values()
        }
    parameters = [parameter[1] for parameter in parameters_dict.values()]
    db = lmdb.open(
        save_lmdb_path,
        map_size=600*(1024*1024*1024),   # 600GB
        create=True,
        subdir=False,
        readonly=False, # Writable
    )
    with Pool(60) as pool:
        for full_reaction, substrate, tag, data in pool.imap_unordered(get_reaction_feature_single, parameters):
            if (data is not None and allow_no_match) or (tag is None):
                with db.begin(write=True, buffers=True) as txn:
                    if resume:
                        txn.put(
                            key = str(reaction_id_dict[full_reaction]).encode(),
                            value = pickle.dumps(data)
                        )
                    else:
                        txn.put(
                            key = str(len(reactions)).encode(),
                            value = pickle.dumps(data)
                        )
                        reactions.append(full_reaction)
                        substrates.append(substrate)
            else:
                num_skipped += 1
                print('Skipped (%d) %s %s' % (num_skipped, full_reaction, tag))
    db.close()

    if not resume:
        data = {
            "reactions": reactions,
            "substrates": substrates
        }
        df = pd.DataFrame(data)
        df.to_csv(save_df_path, sep=',', index=False)
        print("Save reaction list")
    ray.shutdown()

def get_enzyme_feature_small_famlity(df_path, save_df_path, save_lmdb_path):
    df = pd.read_csv(df_path, sep=',')
    env = lmdb.open(
        save_lmdb_path,
        map_size=600*(1024*1024*1024),   # 600GB
        create=True,
        subdir=False,
        readonly=False,  # Writable
    )
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    batch_converter = alphabet.get_batch_converter()
    model = model.to(torch.device("cuda:0"))

    data = []
    uniprot_dict = {}
    sequences = []
    uniprots = []

    if 'active_site' not in df.columns:
        print("No active site information! Fill with -1")
        df['active_site'] = np.ones_like(df['Enzyme_ID']).astype(int)

    # TODO: download sequence from uniprot
    for index, (sequence, uniprot, active_site) in enumerate(zip(df['SEQ'].values, df['Enzyme_ID'].values, df['active_site'].values)):
        if uniprot in uniprot_dict:
            continue
        if len(sequence) > 1000:
            continue
        sequences.append(sequence)
        uniprots.append(uniprot)
        uniprot_dict[uniprot] = (len(uniprots) - 1, active_site)
        data.append((uniprot, sequence))
        
    generate_esm_embedding(env, data, uniprot_dict)
    data = {
        "sequences": sequences,
        "uniprots": uniprots
    }
    df = pd.DataFrame(data)
    df.to_csv(save_df_path, sep=',', index=False)
    print("Save enzyme df")
        # break
    env.close()

def generate_esm_embedding(env, data, uniprot_dict, save_enzyme_df=None, mean_only=False):
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model = model.to(torch.device("cuda:0"))
    batch_converter = alphabet.get_batch_converter()
    for small_data in tqdm(data, total=len(data)):
        small_data = [small_data]
        atch_labels, batch_strs, batch_tokens = batch_converter(small_data)
        batch_tokens = batch_tokens.to(torch.device("cuda:0"))
        # Extract per-residue representations (on CPU)
        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[33], return_contacts=True)
        token_representations = results["representations"][33].cpu()

        # Generate per-sequence representations via averaging
        # NOTE: token 0 is always a beginning-of-sequence token, so the first residue is token 1.
        for i, (uniprot, seq) in enumerate(small_data):
            with env.begin(write=True, buffers=True) as txn:
                if mean_only:
                    embedding = token_representations[i, 1:1+len(seq)].detach().numpy().mean(axis=0)
                else:
                    embedding = token_representations[i, 1:1+len(seq)].detach().numpy()
                data = {
                    'embedding': embedding,
                    'active_site': uniprot_dict[uniprot][1],
                    'sequence': convert_protein_sequence_to_number(seq)
                }
                txn.put(key=str(uniprot_dict[uniprot][0]).encode(), value=pickle.dumps(data))
                
        batch_tokens = batch_tokens.to(torch.device("cpu"))

def get_enzyme_feature_from_df(enzyme_df_path, save_lmdb_path, duplicate=False, save_enzyme_df=None):
    df = pd.read_csv(enzyme_df_path, sep=',')
    env = lmdb.open(
        save_lmdb_path,
        map_size=600*(1024*1024*1024),   # 600GB
        create=True,
        subdir=False,
        readonly=False,  # Writable
    )
    data = []
    sequences = []
    uniprots = []
    uniprot_dict = {}
    for index, (sequence, uniprot) in enumerate(zip(df["sequences"], df["uniprots"])):
        if uniprot in uniprot_dict:
            # print(uniprot, "already in lmdb")
            # continue
            if duplicate:
                uniprot = str(index)
            else:
                print(uniprot, "already in lmdb")
                continue
        if len(sequence) > 1000:
            print(uniprot, "sequence too long")
            continue
        try:
            convert_protein_sequence_to_number(sequence)
        except:
            print(uniprot, "sequence contain non-standard amino acid")
            continue
        sequences.append(sequence)
        uniprots.append(uniprot)
        uniprot_dict[uniprot] =  (len(uniprot_dict), 1)
        data.append((uniprot, sequence))
    generate_esm_embedding(env, data, uniprot_dict)
    env.close()
    if save_enzyme_df is not None:
        data = {
            "sequences": sequences,
            "uniprots": uniprots
        }
        df = pd.DataFrame(data)
        df.to_csv(save_enzyme_df, sep=',', index=False)
        print("Save enzyme df")

def get_enzyme_feature_for_new_brenda(df_path, save_df_path, save_lmdb_path):
    df = pd.read_csv(df_path, sep=',')
    env = lmdb.open(
        save_lmdb_path,
        map_size=600*(1024*1024*1024),   # 600GB
        create=True,
        subdir=False,
        readonly=False,  # Writable
    )

    data = []
    uniprot_dict = {}
    sequences = []
    uniprots = []

    if 'active_site' not in df.columns:
        print("No active site information! Fill with -1")
        df['active_site'] = np.ones_like(df['Entry']).astype(int)
    
    for index, (sequence, uniprot, active_site) in enumerate(zip(df['Sequence'].values, df['Entry'].values, df['active_site'].values)):
        if uniprot in uniprot_dict:
            continue
        if len(sequence) > 1000:
            continue
        sequences.append(sequence)
        uniprots.append(uniprot)
        uniprot_dict[uniprot] = (len(uniprots) - 1, active_site)
        data.append((uniprot, sequence))
    generate_esm_embedding(env, data, uniprot_dict)    
    
    data = {
        "sequences": sequences,
        "uniprots": uniprots
    }
    df = pd.DataFrame(data)
    df.to_csv(save_df_path, sep=',', index=False)
    print("Save enzyme df")
        # break
    env.close()

def get_enzyme_feature(df_path, save_df_path, save_lmdb_path):
    df = pd.read_csv(df_path, sep=',').dropna(subset=['substrate', 'left', 'right'])
    env = lmdb.open(
        save_lmdb_path,
        map_size=600*(1024*1024*1024),   # 600GB
        create=True,
        subdir=False,
        readonly=False,  # Writable
    )

    df = pd.read_csv(df_path, sep=',')

    data = []
    uniprot_dict = {}
    sequences = []
    uniprots = []

    if 'active_site' not in df.columns:
        print("No active site information! Fill with -1")
        df['active_site'] = np.ones_like(df['uniprot']).astype(int)

    # TODO: download sequence from uniprot
    for index, (sequence, uniprot, active_site) in enumerate(zip(df['sequence'].values, df['uniprot'].values, df['active_site'].values)):
        if uniprot in uniprot_dict:
            continue
        if len(sequence) > 1000:
            continue
        sequences.append(sequence)
        uniprots.append(uniprot)
        uniprot_dict[uniprot] = (len(uniprots) - 1, active_site)
        data.append((uniprot, sequence))

    generate_esm_embedding(env, data, uniprot_dict)    
    
    data = {
        "sequences": sequences,
        "uniprots": uniprots
    }
    df = pd.DataFrame(data)
    df.to_csv(save_df_path, sep=',', index=False)
    print("Save enzyme df")
        # break
    env.close()

def get_enzyme_fature_tmap(fasta_path, save_lmdb_path):
    # disable warning
    from transformers import logging
    from Bio import SeqIO
    logging.set_verbosity_error()
    RDLogger.DisableLog('rdApp.*')  

    # parameters)
    env = lmdb.open(
        save_lmdb_path,
        map_size=600*(1024*1024*1024),   # 600GB
        create=True,
        subdir=False,
        readonly=False,  # Writable
    )
    data = []
    uniprot_dict = {}
    records = list(SeqIO.parse(fasta_path, "fasta"))
    
    for record in records[:200000]:
        uniprot = record.id.split("|")[1]
        sequence = str(record.seq)
        if uniprot in uniprot_dict:
            continue
        if len(sequence) > 1000:
            continue
        data.append((uniprot, sequence))
        uniprot_dict[uniprot] = len(uniprot_dict)
    generate_esm_embedding(env, data, uniprot_dict, mean_only=True)
    env.close()

if __name__ == "__main__":
    # Check get_reaction_feature_single function
    # object = get_reaction_feature_single.remote(('CC1(O)CC(=O)OC(O)C1.N=C(O)C1=CN([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)C=CC1>>N=C(O)c1ccc[n+]([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)c1',
    #                                      'CC1(O)CC(=O)OC(O)C1', 
    #                                      'N=C(O)c1ccc[n+]([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](OP(=O)(O)O)[C@@H]3O)[C@@H](O)[C@H]2O)c1'))
    # print(ray.get(object))
    
    # Create reaction features

    # specificity (full)
    # save_df_path = f"{root_dir}/data/new_brenda/reaction.csv"
    # save_lmdb_path = f"{root_dir}/data/new_brenda/reaction_features.lmdb"
    # input_path = f"{root_dir}/data/new_brenda/data.csv"
    # parameters_dict = get_reaction_parameters_new_brenda(df_path=input_path)
    # get_reaction_feature(parameters_dict=parameters_dict, save_lmdb_path=save_lmdb_path, save_df_path=save_df_path, allow_no_match=True)

    # selecitivity (available reaction)
    # save_df_path = f"{root_dir}/data/reaction.csv"
    # save_lmdb_path = f"{root_dir}/data/reaction_features.lmdb"
    # input_path = f"{root_dir}/data/brenda/reduced_data.csv"
    # parameters_dict = get_reaction_parameters_brenda(df_path=input_path)
    # get_reaction_feature(parameters_dict=parameters_dict, save_lmdb_path=save_lmdb_path, save_df_path=save_df_path, allow_no_match=True)
    
    # Check generated reaction features
    # save_lmdb_path = f"{root_dir}/data/halogenase/reaction_features.lmdb"
    # db = lmdb.open(
    #     save_lmdb_path,
    #     map_size=10*(1024*1024*1024),   # 10GB
    #     create=False,
    #     subdir=False,
    #     readonly=True,
    #     lock=False,
    #     readahead=False,
    #     meminit=False,
    # )
    # with db.begin() as txn:
    #     keys = list(txn.cursor().iternext(values=False))
    # for key in keys:
    #      with db.begin(write=False, buffers=True) as txn:
    #         # key = str(key).encode()
    #         value = txn.get(key)
    #         if value is None:
    #             raise KeyError
    #         print(key.decode())
    #         print(pickle.loads(value))
    #         break

    # Create enzyme features
    # save_df_path = f"{root_dir}/data/brenda/enzymes.csv"
    # save_lmdb_path = f"{root_dir}/data/brenda/enzyme_features.lmdb"
    # input_path = f"{root_dir}/data/brenda/data_cofactor.csv"
    # get_enzyme_feature(df_path=input_path, save_lmdb_path=save_lmdb_path, save_df_path=save_df_path)

    # Check generated enzyme features
    # save_lmdb_path = f"{root_dir}/data/enzyme_features.lmdb"
    # db = lmdb.open(
    #     save_lmdb_path,
    #     map_size=10*(1024*1024*1024),   # 10GB
    #     create=False,
    #     subdir=False,
    #     readonly=True,
    #     lock=False,
    #     readahead=False,
    #     meminit=False,
    # )
    # with db.begin() as txn:
    #     keys = list(txn.cursor().iternext(values=False))
    # print(keys)
    
    # Create enzyme features for small family dataset
    # for enzyme_tag in ["Duf", "Esterase", "Phosphatase", "Gt_acceptor", "Nitrilase", "Thiolase"]:
    # for enzyme_tag in ["halogenase"]:
    #     save_df_path = f"{root_dir}/data/small_family/{enzyme_tag}/enzymes.csv"
    #     save_lmdb_path = f"{root_dir}/data/small_family/{enzyme_tag}/enzyme_features.lmdb"
    #     input_path =  f"{root_dir}/data/small_family/{enzyme_tag}/data.csv"
    #     get_enzyme_feature_small_famlity(df_path=input_path, save_lmdb_path=save_lmdb_path, save_df_path=save_df_path)

    # Create enzyme feature for new Brenda dataset
    # save_df_path = f"{root_dir}/data/small_family/2023_brenda/enzymes.csv"
    # save_lmdb_path = f"{root_dir}/data/small_family/2023_brenda/enzyme_features.lmdb"
    # input_path = f"{root_dir}/data/small_family/2023_brenda/data.csv"
    # get_enzyme_feature_for_new_brenda(df_path=input_path, save_lmdb_path=save_lmdb_path, save_df_path=save_df_path)

    # Create enzyme feature from prepared enzyme df
    # save_enzyme_df = f"{root_dir}/data/small_family/Methyltransferse/enzymes_new.csv"
    # get_enzyme_feature_from_df(enzyme_df_path=f"{root_dir}/data/small_family/Methyltransferse/enzymes.csv", save_lmdb_path=f"{root_dir}/data/small_family/Methyltransferse/enzyme_features.lmdb", duplicate=True, save_enzyme_df=save_enzyme_df)
    
    # Create esm embedding for universal protein
    get_enzyme_fature_tmap(fasta_path=f"/scratch/bcao/suyufeng/uniprot_sprot.fasta", save_lmdb_path=f"/scratch/bcao/suyufeng/uniprot_features.lmdb")

    # Create reaction features for prepared reaction df
    # save_lmdb_path = f"{root_dir}/data/small_family/Methyltransferse/reaction_features.lmdb"
    # data_path = f"{root_dir}/data/small_family/Methyltransferse/reactions.csv"
    # parameters = get_reaction_parameters_df(data_path, substrate_only=True)
    # get_reaction_feature(parameters_dict=parameters, save_lmdb_path=save_lmdb_path, save_df_path=None, allow_no_match=True, resume=True)

    # Create reaction features for halogenase
    # df_path = f"{root_dir}/data/small_family/halogenase/data.csv"
    # save_lmdb_path = f"{root_dir}/data/small_family/halogenase/reaction_features.lmdb"
    # save_df_path = f"{root_dir}/data/small_family/halogenase/reactions.csv"
    # parameters = get_reaction_parameters_halogenase(df_path)
    # get_reaction_feature(parameters_dict=parameters, save_lmdb_path=save_lmdb_path, save_df_path=save_df_path, allow_no_match=False)

    # Create reaction features for halogenase evaluation dataset
    # save_lmdb_path = f"{root_dir}/data/halogenase/evaluation/reaction_features.lmdb"
    # save_df_path = f"{root_dir}/data/halogenase/evaluation/reactions.csv"
    # # parameters = get_reaction_parameters_halogenase_evaluation(f"{root_dir}/data/halogenase/evaluation/datas.csv")
    # parameters = get_reaction_parameters_resume_evaluation(f"{root_dir}/data/halogenase/evaluation/reactions.csv")
    # get_reaction_feature(parameters_dict=parameters, save_lmdb_path=save_lmdb_path, save_df_path=save_df_path, allow_no_match=True, resume=True)

    # Create reaction features for small family dataset
    # for enzyme_tag in ["Duf", "Esterase", "Phosphatase", "Gt_acceptor", "Nitrilase", "Thiolase"]:
    #     df_path = f"{root_dir}/data/small_family/{enzyme_tag}/data.csv"
    #     save_lmdb_path = f"{root_dir}/data/small_family/{enzyme_tag}/reaction_features.lmdb"
    #     save_df_path = f"{root_dir}/data/small_family/{enzyme_tag}/reactions.csv"
    #     parameters = get_reaction_parameters_small_family(df_path)
    #     get_reaction_feature(parameters_dict=parameters, save_lmdb_path=save_lmdb_path, save_df_path=save_df_path, allow_no_match=True, resume=False)

    # Create reaction features for new Brenda
    # df_path = f"{data_root_dir}/data/small_family/2023_brenda/data.csv"
    # save_lmdb_path = f"{data_root_dir}/data/small_family/2023_brenda/reaction_features.lmdb"
    # save_df_path = f"{data_root_dir}/data/small_family/2023_brenda/reactions.csv"
    # parameters = get_reaction_parameters_new_brenda(df_path)
    # get_reaction_feature(parameters_dict=parameters, save_lmdb_path=save_lmdb_path, save_df_path=save_df_path, allow_no_match=True, resume=False)

    pass