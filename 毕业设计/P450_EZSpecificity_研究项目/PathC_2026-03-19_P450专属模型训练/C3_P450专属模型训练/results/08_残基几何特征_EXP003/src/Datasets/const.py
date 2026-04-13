import numpy as np

global restype_3to1, restype_1to3, restype_name_to_atom14_names, restypes, res_angle_alt, n_res_chi, atom_type_masks, weights, letter_to_num, num_to_letter

letter_to_num = {'C': 4, 'D': 3, 'S': 15, 'Q': 5, 'K': 11, 'I': 9,
                       'P': 14, 'T': 16, 'F': 13, 'A': 0, 'G': 7, 'H': 8,
                       'E': 6, 'L': 10, 'R': 1, 'W': 17, 'V': 19, 
                       'N': 2, 'Y': 18, 'M': 12, 'X': 20, 'Z': 21, 'U': 22, 'B': 23}

num_to_letter = {v:k for k, v in letter_to_num.items()}

restypes = [
    'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P',
    'S', 'T', 'W', 'Y', 'V'
]

weights = {
    'A': 1,
    'R': 1,
    'N': 1,
    'D': 1,
    'C': 5,
    'Q': 1,
    'E': 1,
    'G': 1,
    'H': 5,
    'I': 1,
    'L': 1,
    'K': 1,
    'M': 5,
    'F': 1,
    'P': 1,
    'S': 1,
    'T': 1,
    'W': 5,
    'Y': 1,
    'V': 1,
}

restype_1to3 = {
    'A': 'ALA',
    'R': 'ARG',
    'N': 'ASN',
    'D': 'ASP',
    'C': 'CYS',
    'Q': 'GLN',
    'E': 'GLU',
    'G': 'GLY',
    'H': 'HIS',
    'I': 'ILE',
    'L': 'LEU',
    'K': 'LYS',
    'M': 'MET',
    'F': 'PHE',
    'P': 'PRO',
    'S': 'SER',
    'T': 'THR',
    'W': 'TRP',
    'Y': 'TYR',
    'V': 'VAL',
}

restype_name_to_atom14_names = {
    'ALA': ['N', 'CA', 'C', 'O', 'CB', '', '', '', '', '', '', '', '', ''],
    'ARG': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2',
            '', '', ''],
    'ASN': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'ND2', '', '', '', '', '',
            ''],
    'ASP': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'OD2', '', '', '', '', '',
            ''],
    'CYS': ['N', 'CA', 'C', 'O', 'CB', 'SG', '', '', '', '', '', '', '', ''],
    'GLN': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'NE2', '', '', '',
            '', ''],
    'GLU': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'OE2', '', '', '',
            '', ''],
    'GLY': ['N', 'CA', 'C', 'O', '', '', '', '', '', '', '', '', '', ''],
    'HIS': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2', '',
            '', '', ''],
    'ILE': ['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'CD1', '', '', '', '', '',
            ''],
    'LEU': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', '', '', '', '', '',
            ''],
    'LYS': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ', '', '', '', '',
            ''],
    'MET': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'SD', 'CE', '', '', '', '', '',
            ''],
    'PHE': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ',
            '', '', ''],
    'PRO': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', '', '', '', '', '', '', ''],
    'SER': ['N', 'CA', 'C', 'O', 'CB', 'OG', '', '', '', '', '', '', '', ''],
    'THR': ['N', 'CA', 'C', 'O', 'CB', 'OG1', 'CG2', '', '', '', '', '', '',
            ''],
    'TRP': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3',
            'CZ2', 'CZ3', 'CH2'],
    'TYR': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ',
            'OH', '', ''],
    'VAL': ['N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', '', '', '', '', '', '',
            ''],
    'UNK': ['', '', '', '', '', '', '', '', '', '', '', '', '', ''],
}

n_res_chi = {
    'A': 0,
    'R': 4,
    'N': 2,
    'D': 2,
    'C': 1,
    'Q': 3,
    'E': 3,
    'G': 0,
    'H': 2,
    'I': 2,
    'L': 2,
    'K': 4,
    'M': 3,
    'F': 2,
    'P': 2,
    'S': 1,
    'T': 1,
    'W': 2,
    'Y': 2,
    'V': 1,
}

restype_3to1 = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H',
    'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q',
    'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
}

res_angle_alt = np.ones((20, 19)).astype(np.int16)
res_angle_alt[3, 2+11] = res_angle_alt[3, 3+11] = -1
res_angle_alt[6, 4+11] = res_angle_alt[6, 5+11] = -1
res_angle_alt[13, 2+11] = res_angle_alt[13, 3+11] = -1
res_angle_alt[18, 2+11] = res_angle_alt[18, 3+11] = -1

atom_type_masks = np.zeros((20, 4, 14), dtype=int)
atom_type_one_hot = np.zeros((14, 4), dtype=int)

for key in restype_name_to_atom14_names.keys():
        if key == 'UNK':
                continue
        keyid = restypes.index(restype_3to1[key])
        atoms = restype_name_to_atom14_names[key]
        masks = []
        for c in ['C', 'N', 'O', 'S']:
                mask = np.zeros(14, dtype=int)
                for index, atom in enumerate(atoms):
                        if atom == 'C' or atom == 'CA':
                                continue
                        if len(atom) > 0 and atom[0] == c:
                                mask[index] = 1
                        else:
                                mask[index] = 0
                masks.append(mask)
        atom_type_masks[keyid] = np.array(masks)


        

