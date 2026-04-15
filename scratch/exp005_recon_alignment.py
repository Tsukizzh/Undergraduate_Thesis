"""
穿透验证：pocket PDB 的 author residue ID 能否直接索引 Enzymes.csv 对应酶的序列字符

拿 sample_000000 (enzyme_id=93, dock index=0)：
  1. 读 pocket PDB (structure/str_tmp_data/pocket/0.pdb)
  2. 解析出 (chain, resid, resname) 列表
  3. 去 Enzymes.csv 第 93 行拿序列
  4. 对每个口袋残基，检查 seq[resid-1] 的氨基酸是否和 resname 一致

如果 100% 一致，pocket_residue_idx 就是 `resid - 1`，任务可行。
如果失配，需要另外建立 UniProt residue -> seq position 映射。

多抽几个样本做这个检验。
"""
import csv
from pathlib import Path

# sample 0: enzyme_id=93, structure_index=0
# 我们需要知道 dock index 怎么对应到 pocket PDB 文件名
# 根据 sample 字段 dataset_id/enzyme_id/substrate_id，我们用 enzyme_id=93 对应 Enzymes.csv 第 93 行

ENZ_CSV = Path("/root/autodl-tmp/EZSpecificity/PathC/P450/data/Enzymes.csv")
POCKET_DIR = Path("/root/autodl-tmp/EZSpecificity/PathC/P450/data/structure/str_tmp_data/pocket")
DATA_CSV = Path("/root/autodl-tmp/EZSpecificity/PathC/P450/data/splits/random/testing_datas_0_pt.csv")

# 3 字母 -> 1 字母
AA3to1 = {
    'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H',
    'ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q',
    'ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y',
}


def parse_pocket_pdb(pdb_path):
    """返回 [(chain, resid, resname), ...] 按首次出现顺序排好序"""
    seen = set()
    residues = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            resname = line[17:20].strip()
            chain = line[21].strip()
            resid = int(line[22:26])
            key = (chain, resid)
            if key not in seen:
                seen.add(key)
                residues.append((chain, resid, resname))
    return residues


def load_csv_rows(path, key_col):
    with open(path, encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))


# 先载入 data.csv 查看哪些 dock_index 对应哪个 enzyme_index
data_rows = load_csv_rows(DATA_CSV, 'Dock Index')
enz_rows = load_csv_rows(ENZ_CSV, None)

print(f"data.csv rows: {len(data_rows)}")
print(f"data.csv columns: {list(data_rows[0].keys())}")
print(f"Enzymes.csv rows: {len(enz_rows)}")
print(f"Enzymes.csv columns: {list(enz_rows[0].keys())}")
print()

# 找 5 个不同酶的样本做穿透检验
# 简单策略：data.csv 前 100 行里挑 5 个不同 Enzyme Index
seen_enz = set()
samples_to_check = []
for row in data_rows:
    eid = int(row['Enzyme Index'])
    did = int(row['Dock Index'])
    if did < 0:
        continue
    if eid in seen_enz:
        continue
    seen_enz.add(eid)
    samples_to_check.append((did, eid))
    if len(samples_to_check) >= 5:
        break

print("=== 穿透验证 ===")
print(f"{'dock':>5} {'enz':>5} {'pdb_resid':>10}  {'pdb_aa':>7}  seq_aa_at_(resid-1)")
print("-" * 60)

for dock_index, enz_index in samples_to_check:
    pdb_path = POCKET_DIR / f"{dock_index}.pdb"
    if not pdb_path.exists():
        print(f"dock={dock_index}: pocket PDB 不存在")
        continue

    residues = parse_pocket_pdb(pdb_path)
    if not residues:
        print(f"dock={dock_index}: 空 pocket")
        continue

    # 酶序列
    seq = enz_rows[enz_index]['Protein sequence'].strip()
    uniprot = enz_rows[enz_index]['uniprots']

    print(f"\n--- dock={dock_index}, enz={enz_index}, uniprot={uniprot}, seq_len={len(seq)} ---")
    print(f"pocket 残基数: {len(residues)}")

    # 随机抽 10 个口袋残基做检查
    import random
    random.seed(dock_index)
    sample_res = random.sample(residues, min(10, len(residues)))

    n_match = 0
    for chain, resid, resname in sample_res:
        expected = AA3to1.get(resname, '?')
        # 试 resid-1
        if 0 <= resid - 1 < len(seq):
            seq_aa = seq[resid - 1]
        else:
            seq_aa = '(越界)'
        match = '✓' if seq_aa == expected else '✗'
        if match == '✓':
            n_match += 1
        print(f"  {dock_index:>5} {enz_index:>5} {resid:>10}  {resname} ({expected})   seq[{resid-1}]={seq_aa}  {match}")

    print(f"  匹配率: {n_match}/{len(sample_res)} ({100*n_match/len(sample_res):.0f}%)")
