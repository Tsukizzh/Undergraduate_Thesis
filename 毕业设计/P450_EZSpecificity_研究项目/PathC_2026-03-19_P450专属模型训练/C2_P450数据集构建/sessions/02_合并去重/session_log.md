# Phase 4: 合并去重 Session Log

**日期**: 2026-03-24
**目标**: 将S1-S9五个数据源合并去重，生成全局唯一的酶/化合物/交互表

---

## 1. 输入数据

| 数据源 | 酶 | 化合物 | 交互对 | 质量等级 |
|--------|-----|--------|--------|---------|
| S1_RCSB | 103 | 220 | 272 | Tier A |
| S2_ESIBank | 338 | 390 | 806 | Tier B |
| S3_P450Rdb | 857 | 1,492 | 2,798 | Tier B |
| S8_PlantP450DB | 578 | 295 | 979 | Tier B |
| S9_PCPD | 818 | 570 | 1,209 | Tier B |
| **合计** | **2,694** | **2,967** | **6,064** | |

## 2. 去重策略设计（Codex 3轮讨论）

### Round 1: 策略设计
- **酶去重主键**: UniProt ID（非空时）
- **酶去重fallback**: 序列SHA256 hash（PCPD有447个酶无UniProt）
- **化合物去重主键**: RDKit标准化 → Standard InChIKey
- **化合物去重fallback**: Canonical SMILES
- **交互去重**: (global_enzyme_id, global_compound_id) 对
- **S6_Figshare排除**: 仅6个酶+3,258化合物，辅助数据不进主benchmark

### Round 2: 边界情况审查
- **PCPD序列桥接**: 无UniProt酶如果序列完全匹配有UniProt的酶 → 合并并继承UniProt ID
- **同UniProt不同序列**: 合并，选最长序列为canonical，标记sequence_conflict
- **化合物标准化**: 保留立体异构（isomericSmiles=True），salt stripping默认关闭
- **InChIKey生成**: MolToInchi → InchiToInchiKey（修正Codex原型中MolToInchiKey的bug）
- **Label检查**: 验证所有label均为1（正样本），无冲突

### Round 3: 代码实现review
- Codex给出完整实现原型（~600行）
- 我修正了RDKit InChIKey API bug并简化代码
- 最终脚本: `scripts/02_合并去重/merge_dedup.py`

## 3. 去重结果

### 3.1 酶去重: 2,694 → 1,622

| 合并方式 | 数量 | 说明 |
|---------|------|------|
| uniprot_id | 1,457 | UniProt ID精确匹配 |
| uniprot_id + seq_bridge | — | 282个PCPD酶通过序列桥接继承UniProt |
| sequence_hash | 165 | 无UniProt，按序列hash合并 |
| **总计** | **1,622** | |

- 584个UniProt ID出现在>1个来源
- 序列冲突（同UniProt不同序列）: 11个，均为长度差<5%的minor冲突

### 3.2 化合物去重: 2,967 → 2,125

| 合并方式 | 数量 | 说明 |
|---------|------|------|
| inchikey | 2,110 | InChIKey精确匹配 |
| canonical_smiles | 15 | InChIKey生成失败，用canonical SMILES |
| **总计** | **2,125** | |

- 449个SMILES在多个来源中重复（raw match），经RDKit标准化后进一步合并
- 15个化合物含通配符（如`*`），InChIKey无法生成，用canonical SMILES fallback

### 3.3 交互去重: 6,064 → 4,751

- 1,313条重复交互被合并（21.7%去重率）
- 0条交互因缺失enzyme/compound映射被跳过
- 0条label冲突
- 所有交互label=1（纯正样本集）

## 4. 输出文件

```
data/combined/
├── global_enzymes.csv       (1,622行)   ← 主表
├── global_compounds.csv     (2,125行)   ← 主表
├── global_interactions.csv  (4,751行)   ← 主表
├── enzyme_xref.csv          (2,694行)   ← 源→全局映射
├── compound_xref.csv        (2,967行)   ← 源→全局映射
├── interaction_evidence.csv (6,064行)   ← 每条原始证据
└── merge_audit.csv          (343行)     ← 冲突与桥接记录
```

### 字段说明

**global_enzymes.csv**: global_enzyme_id, canonical_uniprot_id, canonical_p450_symbol, canonical_species, canonical_species_class, canonical_ec_number, canonical_sequence, sequence_hash, sequence_length, sources, source_count, member_count, merge_basis, sequence_conflict, species_conflict, ec_conflict, symbol_aliases, species_aliases, ec_number_all

**global_compounds.csv**: global_compound_id, canonical_name, canonical_smiles, standard_inchi, standard_inchikey, sources, source_count, member_count, merge_basis, structure_parse_status, name_aliases, original_smiles_count

**global_interactions.csv**: global_interaction_id, global_enzyme_id, global_compound_id, label, sources, supporting_source_count, evidence_count, best_quality_tier, has_pmid_any, has_products_any, has_pdb_any, max_num_reactions, conflicting_labels

## 5. 目录结构清理

同时进行了C2目录整理：
- 删除空占位目录（旧编号02_格式标准化、03_合并去重等）
- 统一编号: scripts/01_数据下载 → 02_合并去重 → 03_负样本与Split
- 删除空的skipped数据源目录（BM3Variants, BRENDA, CYPED等）
- data/P450_Family/ 空子目录删除，待Split阶段重建

## 6. 待验证事项

- [ ] Codex Round 4: 验证去重结果正确性（抽样检查）
- [ ] 11个同序列不同UniProt的酶：split时需按sequence_hash分组防泄露

## 7. 下一步

- Phase 5: 负样本生成（双向随机采样）+ 4种Split × 4折CV
- Path A的245条抑制剂负样本 → `data/combined/negatives/biological.csv`
