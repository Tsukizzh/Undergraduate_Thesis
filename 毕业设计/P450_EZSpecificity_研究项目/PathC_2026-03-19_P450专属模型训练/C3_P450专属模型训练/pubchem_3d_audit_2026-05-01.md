# PubChem CID 与 3D 构象覆盖率审计

**日期**：2026-05-01
**目的**：回答王老师 4/22 组会问题"对接出来的小分子奇怪，建议改用 PubChem 现成的三维结构"——量化 PubChem 在我们 P450 数据集上的覆盖率，决定能否替换现有 RDKit ETKDG 管线。

**审计方法**：claude + codex 多轮交叉核对（codex session `019dcfdf-2bd8-7143-941b-f8b12abbbbbb`），所有数字经 14 项独立闭合检查。

---

## 一句话结论

**87% 的训练用底物可以用 PubChem 现成 3D 结构，13% 必须保留 RDKit ETKDG 兜底。**

采用混合管线：能换就换，换不了用原来的。

---

## 关键数字

### PubChem CID 覆盖率（Stage 1+2）

数据源：`Substrates.csv`，2,125 条 SMILES，全部唯一。

| 类别 | 数量 | 占比 |
|---|---|---|
| 真 PubChem CID 命中 | 1,969 | 92.7% |
| ⤷ 其中互变异构碰撞（2 SMILES → 同 CID 160478） | 2 | — |
| ⤷ 实际不同 CID 数 | 1,968 | — |
| 通配符 SMILES（PubChem 拒绝） | 15 | 0.7% |
| PubChem 没收录（cid0_no_match） | 139 | 6.5% |
| 网络错误彻底失败 | 2 | 0.1% |
| **需 fallback 合成 ID** | **156** | **7.3%** |

### PubChem 3D 构象覆盖率（Stage A）

查询 1,968 个唯一 CID 的 `ConformerCount3D, IsomericSMILES, InChIKey, HeavyAtomCount, RotatableBondCount`，批量端点 100 CID/请求 × 20 批，**16 秒跑完**。

| 维度 | 有 3D | 无 3D | 覆盖率 |
|---|---|---|---|
| 全部 1,968 唯一 CID | 1,854 | 114 | 94.2% |
| 全部 1,969 底物条目 | 1,855 | 114 | 94.2% |
| **实际训练用底物（in_any_split=1，2,111 条）** | **1,849** | 113 + 149 = 262 | **87.6%** |

### InChIKey 一致性核对（PubChem vs 本地 RDKit）

15 条不一致——分两类：

| 类型 | 数量 | 处理 |
|---|---|---|
| **STEREO_ONLY**（同 connectivity 不同立体）| 9 | 立体异构差异，P450 对此敏感，**用本地 ETKDG** |
| **STRUCT_DIFF**（不同 connectivity）| 6 | PubChem 误匹配（互变异构归一化 / 怪异 valence SMILES），**必须用本地 ETKDG** |

实际训练集中 14 条 mismatch（idx=473 不在训练集）。

### 三档严格度下的覆盖率（实际训练用 2,111 条底物）

```
A. 原始可用（Raw availability）         : 1,849 / 2,111 = 87.59%
B. 立体异构容忍（Connectivity match）   : 1,846 / 2,111 = 87.45%
C. 严格身份验证（InChIKey exact match） : 1,837 / 2,111 = 87.02%  ← 论文报告用此
```

**论文叙事统一用 C 档：87.0%**——保守、可复现、对 P450 特异性研究最严谨。

---

## STRUCT_DIFF 红旗清单（6 条）

PubChem 把我们的 SMILES 匹配到了一个**化学式相同但 connectivity 不同**的 CID。这是因为 PubChem 在 SMILES → CID 查询时做了 tautomer / 位置异构归一化。

| substrate_idx | CID | 原始 SMILES | 性质 |
|---|---|---|---|
| 239 | 6755 | `O=C1C=C(O)C(=O)c2ccccc21` | 羟基萘醌位置异构（lawsone vs 2-OH-1,4-NQ） |
| 473 | 53356674 | `C=CC1=C(C)C2=Cc3...[Fe-]35...`（卟啉 Fe 复合物） | 不在训练集 |
| 905 | 20590621 | `[C]#CC#CC#CC#CC#CC#CC#[C]` | 多炔片段，怪异 valence |
| 1354 | 878 | `C[S]` | 硫的怪异 valence 表示 |
| 1423 | 160478 | `O=C1C(O)=CC(=O)c2c(O)cc(O)cc21` | 蒽醌 tautomer 漂移（同时是 STEREO 系统第二个互变异构对） |
| 1938 | 440559 | `CC(C)=CCC[C@@H](C)[C@H]1...`（甾醇衍生物） | tautomer / 立体差异 |

---

## STEREO_ONLY mismatch（9 条）

同连接性不同立体，用 ETKDG 兜底（保留原始 SMILES 的立体）。

idx 列表：142, 729, 730, 1206, 1207, 1355, 1466, 1593, 1594。

---

## 决策与执行计划

### 混合管线方案

```
对每条底物：
  if pubchem_3d_count > 0 and inchikey_match:
      用 PubChem 3D SDF (87.0% 走这条)
  else:
      用 RDKit ETKDG（管线原方法，13.0% 走这条）
```

### 执行步骤

1. ✅ **Manifest 已就绪**：`substrates_manifest.csv`（18 列 × 2125 行）+ `pubchem_3d_coverage.csv`（9 列 × 2125 行）已落盘 C3 目录
2. ⏳ 修对接管线 `postprocess_to_complex.py` 的元素错位 bug（独立任务，不在此审计内）
3. ⏳ 重做对接：用 manifest 决定每条底物的 3D 来源
4. ⏳ 重命名复合物 `{UniProt}_{final_substrate_id}.pdb`

---

## 产出清单（C3 目录下）

| 文件 | 内容 | 大小 |
|---|---|---|
| `substrates_manifest.csv` | 主 manifest，2125 行 × 18 列 | 488 KB |
| `pubchem_3d_coverage.csv` | 每条底物的 3D 信息（ConformerCount3D / PubChem InChIKey / 重原子数 / 旋转键数 / InChIKey 一致性 1/0） | 108 KB |
| `pubchem_3d_raw.json` | PubChem 批量 API 原始响应 | 508 KB |

---

## Manifest 字段说明（18 列）

```
1.  substrate_idx              0..2124 行号
2.  raw_smiles                 原始 SMILES
3.  rdkit_parse_ok             True/False
4.  rdkit_canonical_smiles     RDKit 非立体 canonical
5.  rdkit_isomeric_smiles      RDKit 立体 canonical
6.  rdkit_inchikey             本地算的 InChIKey
7.  cid                        PubChem CID（空 = 没拿到）
8.  cid_status                 hit / cid0_no_match / wildcard / rescued_retry / rescued_inchikey / still_error
9.  final_substrate_id         文件命名用：CID12345 / CID12345_<sha12>（碰撞）/ NO_CID_<sha12>
10. fallback_id                NO_CID_<sha12>（无论是否命中都生成）
11. cid_collision_group        若有 CID 碰撞，列出共享该 CID 的其他 substrate_idx
12. pubchem_query_method       hit / hit_with_collision / cid0_no_match / wildcard / ...
13. pubchem_return_count       PubChem 返回的 CID 数量
14. notes                      失败原因等
15. in_train                   1/0 是否在 train split
16. in_val                     1/0
17. in_test                    1/0
18. in_any_split               1/0 = 实际是否被训练用过
```

`pubchem_3d_coverage.csv` 增补列：
```
pubchem_conformer_count_3d     0 = 没 3D, ≥1 = 有 3D
pubchem_isomeric_smiles        PubChem 返回的立体 SMILES
pubchem_inchikey               PubChem 的 InChIKey
pubchem_heavy_atoms            重原子数（OMEGA 限 ≤50）
pubchem_rotatable_bonds        旋转键数（OMEGA 限 ≤15）
inchikey_match                 1 = 我们和 PubChem 算的 InChIKey 一致, 0 = 不一致
```

---

## 复盘（codex 多轮审核要点）

**Round 1**: codex 推荐批量端点 + 100 CID/批 + 同批次顺便查 InChIKey 做一致性核对——加速 84×（14 min → 16 sec）。

**Round 2**: codex 抓出"软标准化"陷阱——PubChem SMILES → CID 查询会做 tautomer 归一化，可能返回 connectivity 不同的 CID。要求做 InChIKey 一致性核对，分 STEREO_ONLY / STRUCT_DIFF 两类。

**Round 3**（待跑）: Stage B 下载所有 1,854 个有 3D 的 SDF，验证下载 → RDKit 解析 → 含 3D 坐标的实际操作链路。codex 强烈建议跑（10 分钟），但用户决定到 Stage A 已经足够回答老师的问题。

**审计闭合性**:
- Stage 1 桶互斥：15 + 29 + 2081 = 2125 ✓
- Stage 2 跨桶无泄漏：cid0 ∩ err = ∅ ✓
- 14 unused 行级闭合：7 真 CID + 7 fallback = 14 ✓
- final_substrate_id 全部唯一：2125 不重复 ✓
- 0-based 索引验证：服务器 ID 范围 0..2118 ⊂ CSV 0..2124 ✓
- 14 项独立断言全部通过

---

**审计完成**：2026-05-01
**结论可发布给王老师**：是
**下一步**：进入对接管线修复 + 重做对接（本审计已为重做提供了 manifest）
