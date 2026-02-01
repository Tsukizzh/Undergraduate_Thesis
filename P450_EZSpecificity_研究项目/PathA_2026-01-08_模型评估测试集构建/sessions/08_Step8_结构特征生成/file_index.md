# Step 8: 结构特征生成 - 文件索引

> **最后更新**: 2026-01-31

## 脚本文件

| 文件 | 路径 | 描述 |
|------|------|------|
| step8_extract_pocket_ligand.py | scripts/08_Step8_结构特征生成/ | 口袋/配体提取脚本（BioPython NeighborSearch） |
| step8_align_ligand.py | scripts/08_Step8_结构特征生成/ | 配体SMILES对齐脚本（RDKit模板匹配） |
| step8_generate_structure_lmdb.py | scripts/08_Step8_结构特征生成/ | LMDB生成脚本（StructureComplexData序列化） |
| test_extraction.py | scripts/08_Step8_结构特征生成/ | 小规模测试脚本 |
| verify_downstream.py | scripts/08_Step8_结构特征生成/ | 下游兼容性验证脚本 |

## 数据文件

### 输入

| 文件 | 路径 | 描述 |
|------|------|------|
| *.pdb / *.cif | data/01_Step1_PDB文件/ | 627个PDB结构文件（619个.pdb + 8个.cif） |
| dock_index_mapping.csv | data/03_Step3_配体处理与数据集构建/06_数据集切分/ | Dock_Index → PDB_ID, 配体CCD映射（539条） |
| Substrates.csv | data/04_Step4_格式修正后数据/ | 底物SMILES（436条） |

### 输出 - 主要

| 文件 | 路径 | 大小 | 描述 |
|------|------|------|------|
| structure_features.lmdb | data/08_Step8_结构特征/ | ~10 GB | StructureComplexData（511条记录） |
| high_quality_id.txt | data/08_Step8_结构特征/ | 2.5 KB | 有效Dock_Index列表（511个） |

### 输出 - 中间文件

| 文件/目录 | 路径 | 数量/大小 | 描述 |
|-----------|------|----------|------|
| pocket/*.pdb | data/08_Step8_结构特征/str_tmp_data/ | 539个 | 提取的口袋PDB（包括失败的空文件） |
| raw_ligand/*.sdf | data/08_Step8_结构特征/str_tmp_data/ | 532个 | 原始提取的配体SDF（v3.0修复altLoc后） |
| ligand/*.sdf | data/08_Step8_结构特征/str_tmp_data/ | 511个 | 对齐后的配体SDF（含AtomMapNum） |
| extraction_summary.csv | data/08_Step8_结构特征/str_tmp_data/ | - | 口袋/配体提取统计 |
| extraction_failures.csv | data/08_Step8_结构特征/str_tmp_data/ | - | 提取失败记录 |
| alignment_summary.csv | data/08_Step8_结构特征/str_tmp_data/ | - | 配体对齐统计 |
| lmdb_failures.csv | data/08_Step8_结构特征/str_tmp_data/ | - | LMDB生成失败记录 |

### 修改的源文件

| 文件 | 路径 | 修改内容 |
|------|------|----------|
| step8_align_ligand.py:52 | scripts/08_Step8_结构特征生成/ | encoding='utf-8-sig'（处理BOM） |
| step8_align_ligand.py:111-114 | scripts/08_Step8_结构特征生成/ | try-except包裹RemoveHs（处理无效valence） |
| step8_generate_structure_lmdb.py:68 | scripts/08_Step8_结构特征生成/ | StructureComplexData.__module__ = 'Datasets.Structure.utils' |
| step8_generate_structure_lmdb.py:71-79 | scripts/08_Step8_结构特征生成/ | 创建fake module hierarchy（Pickle兼容） |

## StructureComplexData 字段说明

| 字段 | 形状 | 类型 | 描述 |
|------|------|------|------|
| protein_element | (N_pocket,) | int | 口袋原子元素（原子序数） |
| protein_pos | (N_pocket, 3) | float32 | 口袋原子3D坐标 |
| protein_is_backbone | (N_pocket,) | bool | 是否骨架原子 |
| protein_atom_to_aa_type | (N_pocket,) | int | 氨基酸类型索引 |
| ligand_element | (N_ligand,) | int | 配体原子元素 |
| ligand_pos | (N_ligand, 3) | float32 | 配体原子3D坐标 |
| ligand_bond_index | (2, N_bonds) | int | 配体键连接 |
| ligand_bond_type | (N_bonds,) | int | 键类型 |
| ligand_atom_feature | (N_ligand, 5) | int | 原子特征矩阵 |
| ligand_index | (N_ligand,) | int | **SMILES原子索引映射**（关键字段） |
| ligand_center_of_mass | (3,) | float32 | 配体质心 |

## 会话文档

| 文件 | 路径 | 描述 |
|------|------|------|
| session_log.md | sessions/08_Step8_结构特征生成/ | 详细会话日志（~650行） |
| file_index.md | sessions/08_Step8_结构特征生成/ | 本文件索引 |

## 下游使用注意

### 1. 反序列化需要fake module注入

```python
import types, sys

if 'Datasets' not in sys.modules:
    sys.modules['Datasets'] = types.ModuleType('Datasets')
if 'Datasets.Structure' not in sys.modules:
    sys.modules['Datasets.Structure'] = types.ModuleType('Datasets.Structure')
if 'Datasets.Structure.utils' not in sys.modules:
    from torch_geometric.data import Data
    class StructureComplexData(Data):
        pass
    fake_utils = types.ModuleType('Datasets.Structure.utils')
    fake_utils.StructureComplexData = StructureComplexData
    sys.modules['Datasets.Structure.utils'] = fake_utils

# 然后才能 pickle.loads
```

### 2. 过滤到有效ID

```python
with open('high_quality_id.txt') as f:
    valid_ids = set(int(line.strip()) for line in f)

df = df[df['Dock Index'].isin(valid_ids)]  # 511条
```
