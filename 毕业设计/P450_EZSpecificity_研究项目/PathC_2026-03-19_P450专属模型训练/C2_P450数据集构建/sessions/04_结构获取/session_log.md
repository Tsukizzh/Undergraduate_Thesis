# Phase 6: 结构获取与对接 Session Log

**日期**: 2026-03-24
**状态**: 进行中

---

## 1. 目标

为1,622个P450酶获取含血红素的3D结构，为2,125个底物生成3D坐标，最终完成52,254个酶-底物对的分子对接。

## 2. 已完成

### 2.1 酶结构来源覆盖率统计

1,622个酶按结构来源分为四类：

| 类别 | 酶数量 | data.csv行数 | 占比 | 说明 |
|------|--------|-------------|------|------|
| S1实验PDB | 103 | 6,802 | 13.0% | RCSB晶体结构,已有文件 |
| PCPD预测PDB | 246 | 7,435 | 14.2% | Maestro生成,全部含血红素,已有文件 |
| AlphaFill下载成功 | 602 | 19,481 | 37.3% | 从AlphaFill网站批量下载,含血红素 |
| **小计(已有结构)** | **951** | **33,718** | **64.5%** | |
| AlphaFill失败(真404) | 525 | — | — | AlphaFill数据库中不存在 |
| AlphaFill无血红素 | 15 | — | — | AlphaFill有结构但无HEM模板 |
| 无UniProt(需ColabFold) | 131 | 2,306 | 4.4% | 无真实UniProt ID |
| **小计(待处理)** | **671** | **18,536** | **35.5%** | |

### 2.2 AlphaFill批量下载

- 脚本: `scripts/04_结构获取/download_alphafill.py`
- 输出: `data/structures/alphafill/cif/` (602个CIF文件) + `data/structures/alphafill/json/` (元数据)
- Manifest: `data/structures/manifests/receptor_manifest.csv`

### 2.3 发现并修复的Bug

**SEQHASH假UniProt Bug**: 24个酶的`canonical_uniprot_id`字段存储了内部序列哈希(`SEQHASH_...`)或GenBank编号(`KAJ5738141`)，不是真正的UniProt ID。

- 原因: Phase 4 merge_dedup.py缺少UniProt格式验证
- 影响: 这24个被错误分类为"AlphaFill候选"，去查询API当然返回404
- 修复: 在manifest中将这24个重新分类为`needs_colabfold`
- 正确数字: 真正有UniProt的酶=1,433个(不是之前的1,457个)

### 2.4 PCPD PDB文件质量检查

- 423个PCPD PDB文件，**全部含血红素**(HEM残基)
- 290个标注为Maestro/Schrodinger生成，133个无头信息但也是预测结构
- **0个实验结构**，全部是计算预测的

### 2.5 PyMOL手动血红素移植教程

用户亲手操作了一次完整的血红素移植流程：
1. 加载目标结构(Q6LEN2, AlphaFold裸结构)
2. 加载模板结构(P05181/CYP2E1, AlphaFill含血红素)
3. `align`对齐两个结构(RMSD=8.5Å, 因为选的模板不够近缘)
4. `create`提取血红素(43个原子)
5. `save`合并保存为target_with_heme.pdb
6. 学会了Wizard → Measurement测量原子距离
7. 学会了label命令显示残基名

## 3. 当前困难

### 3.1 540个酶需要补血红素

525个AlphaFill 404 + 15个无血红素 = **540个酶**有UniProt但没有含血红素的结构。

**解决方案**: 从AlphaFold DB下载裸结构 → 脚本自动移植血红素(和用户在PyMOL手动做的一样,但自动选最近缘模板)。

**关键**: 模板选择很重要。RMSD < 3Å才算好的对齐。需要为每个目标酶找到序列最相似的、已有血红素的P450结构作为模板。我们已经有951个含血红素的结构可以作为模板库。

### 3.2 131个酶需要ColabFold预测

没有真正的UniProt ID → 不能从AlphaFold DB下载 → 需要用ColabFold从序列预测3D结构 → 再补血红素。

**解决方案**: Google Colab上运行ColabFold,输入序列,等待预测。

### 3.3 底物3D坐标还未生成

2,125个底物的SMILES需要转成3D坐标(SDF/PDBQT格式),还没开始。

## 4. 下一步计划

| 优先级 | 任务 | 方式 | 预计耗时 |
|--------|------|------|---------|
| 1 | 批量下载540个酶的AlphaFold裸结构 | 脚本,本地 | ~15分钟 |
| 1 | 生成2,125个底物的3D坐标 | 脚本,本地 | ~30分钟 |
| 2 | 批量移植血红素到540个裸结构 | 脚本(自动align+copy) | ~1-2小时 |
| 3 | ColabFold预测131个酶的结构 | Google Colab | 数小时 |
| 4 | 清理349个现成PDB(去杂质) | 脚本 | ~30分钟 |
| 5 | 受体准备(加氢,PDBQT,定义盒子) | 脚本 | 数小时 |
| 6 | Pilot对接(100-500对,CPU Vina) | 本地 | 数小时 |
| 7 | 全量对接(52,254对,AutoDock-GPU) | Cloud-2服务器 | 数小时-1天 |

优先级1的两个任务可以并行。
