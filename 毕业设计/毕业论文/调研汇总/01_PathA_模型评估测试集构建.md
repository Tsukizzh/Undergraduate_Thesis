# Path A 深度调研报告：P450酶底物特异性预测模型评估

调研日期：2026-04-23
调研来源：Path A 完整 session logs 与进度日志

---

## 1. 数据来源与规模

**独立测试集构成**：
- **原始数据**：RCSB PDB 中的 P450-配体复合物晶体结构
- **初始记录**：**740条** 酶-配体对
- **Photosystem I 污染排除**：**58条** 非 P450 污染（来自 PDB 6ZZX 的 bug）→ **682条** 纯 P450
- **PDB 文件总数**：**627 个** 唯一结构
- **唯一酶序列数**：**292 条**（远大于 153 个 UniProt ID，因包含同一蛋白的突变体）

**为什么不用 ESIBank 训练集中的 P450**：
- ESIBank P450 来自 BRENDA（文献报道的催化关系），我们来自 PDB（晶体共结构）
- 为避免数据泄漏明确排除 ESIBank 条目
- 酶重叠验证：我们 148 个 UniProt ID 与 ESIBank 的 389 个 P450 **0% 重叠**

---

## 2. 配体分类与处理流程

**四类分类标准（682 条最终分布）**：

| 类别 | 数量 | 占比 | 定义 |
|------|------|------|------|
| SUBSTRATE 底物 | 271 | 39.7% | 酶催化反应的反应物（文献+结构或 Kcat/Km） |
| INHIBITOR 抑制剂 | 245 | 35.9% | 能结合但不能被催化（Type-I/II 光谱或 Ki/IC50） |
| EXCLUDE 排除 | 143 | 21.0% | 辅因子（HEM）、溶剂、缓冲液 |
| PRODUCT 产物 | 23 | 3.4% | 酶催化反应的产物 |

**三方深度验证流程（Sub-task 04）**：
1. **并行验证阶段**：6 窗口并行处理 682 条，标记 39 条 `needs_human_review=true`
2. **深度复审阶段**：Claude 主导 + Codex 结构分析 + Gemini 文献搜索 → 三方辩论共识
3. **审核结果**：维持原分类 35 条，修改 4 条（如 REC_0352: SUBSTRATE→INHIBITOR，基于 Fe 距离 7.55Å；REC_0434 基于 PMID:19801656）

**Photosystem I 污染详解**：
- **发现**：UniProt W8SY74、W8SUA3 各 29 条（PDB 6ZZX）
- **根因**：RCSB 下载脚本 bug（遍历所有 polymer entities 而非仅 P450 匹配链）
- **处理**：58 条全部标记非 P450 删除，740→682
- **质量影响**：Enzymes.csv 294→292

---

## 3. 特征生成链路

| Step | 任务 | 输入 | 输出 | 完成度 |
|------|------|------|------|--------|
| Step 5 | ESM 酶嵌入 | 292 条酶序列 | enzyme_features.lmdb (~700MB) | 100% |
| Step 6 | 反应图特征 | 436 条 SMILES | reaction_features.lmdb | 100% |
| Step 7.1 | Morgan 指纹 | 436 条 SMILES | morgan_fingerprint.npy (436×1024) | 100% |
| Step 7.2 | GROVER 指纹 | 436 条 SMILES | grover_fingerprint.lmdb (~10GB) | 100% |
| Step 8 | 结构特征 | 539 条 PDB | structure_features.lmdb (517 条) | 96%（22 条对齐失败） |

**关键卡点及解决**：
- **Step 7.2 GROVER 生成**：内存爆炸（98% RAM）→ 自定义 `build_vocab_low_memory(num_workers=1)`；LMDB 中文路径报错→使用临时 ASCII 路径
- **Step 8 配体对齐**：altLoc 过滤过严（v3.0）→ 添加 `_normalize_altloc` 处理 `.` 和 `?`；RDKit 配体解析失败（v3.1）→ `MolFromPDBBlock` 添加 element 从列 77-78

---

## 4. 模型推理与结果

**5 个 Run 的 Checkpoint 配置**：
- Run01：`best-checkpoint`（原始 EZSpecificity 模型）
- Run02-05：`best-checkpoint-v1` 至 `best-checkpoint-v4`（4 个其他预训练 checkpoint）

**最终评估结果（Run01，基准）**：
- **样本统计**：517 条预测（B3 完整数据集：271 正 + 268 负（含 23 产物 + 245 抑制剂））
- **核心指标**：**AUC-ROC = 0.6636**（AUPR = 0.6360）
- **阈值问题**：默认 0.5 不匹配 → 最优阈值仅 **0.0379**，模型对 P450 数据严重保守（logit 均值 -4.01，概率均值 14.2%）
- **准确率细节**：Accuracy=0.5164, Precision=0.6341, Recall=0.0996（极低）, F1=0.1722
- **正负样本比**：271:268 ≈ 1:0.99，相对平衡

---

## 5. 关键结论

### 为什么结果比论文低（0.6636 vs 0.7198）？

**论文基准对标（Nature Figure 3a）**：
- Random（都见过）：0.8927
- Unknown enzyme：0.7976
- **Unknown enzyme + substrate**：**0.7198** ← 最接近我们的场景

**低 5.6% 原因（三方共识）**：

| 原因 | 重要程度 | 说明 |
|------|:------:|------|
| 任务不匹配 | ★★★ | 论文负样本=随机配对（通常"无关"），我们=抑制剂（能结合但不催化）→ 区分难度大幅提升 |
| 酶完全未知 | ★★ | 148 个测试酶与 ESIBank 0% 重叠 |
| 域偏移 | ★★ | PDB 配体（抑制剂为主）vs BRENDA 底物的化学空间差异 |
| 概率校准失效 | ★ | 最优阈值 0.04 远低于 0.5 |

### 这个数字合理吗？

**合理，完全在预期范围内。** 我们的场景介于论文的 "Unknown enzyme"(0.7976) 和 "Unknown enzyme+substrate"(0.7198) 之间，0.6636 比 0.7198 低 5.6% 反映了**负样本定义差异**导致的难度提升。

---

## 6. 关键文件路径

| 文档 | 路径 |
|------|------|
| 总进度 | `PathA_.../进度日志.md` |
| 特征生成详解 | `PathA特征生成流程详解_从数据到模型输入.md` |
| Step 3 配体分类 | `sessions/03_Step3.../01_配体分类审核/session_log.md` |
| Step 9 推理 | `sessions/09_Step9.../Run01_best-checkpoint/session_log.md` |
| Step 10 结果分析 | `sessions/10_Step10.../Run01_best-checkpoint/session_log.md` |
| v1 废弃说明 | `_v1_按UniProt去重_Step1至3_158PDB/README.md` |

**核心数据集位置**：
- 酶序列：`data/04_Step4_格式修正后数据/Enzymes.csv`（292 条）
- 底物 SMILES：`data/04_Step4_格式修正后数据/Substrates.csv`（436 条）
- B3 完整测试集：`data/03_Step3.../06_数据集切分/B3_完整数据集_271pos268neg/data.csv`（539 条）
- 最终预测：`sessions/09_Step9.../Run01_best-checkpoint/predictions.csv`（517 条, AUC=0.6636）
