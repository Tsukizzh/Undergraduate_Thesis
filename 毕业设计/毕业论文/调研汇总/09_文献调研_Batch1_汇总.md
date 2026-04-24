# 文献调研 Batch 1 — 收敛汇总

**用途**：为第 3 章（P450 专属基准数据集构建）提供数据库 / 对接工具 / 指纹库 / 特征生成模型 / 化学标识符的权威引文
**来源**：用户提供 4 版 AI 检索结果，本文档按以下规则收敛
**收敛规则**：
1. 多版本冲突时，优先最新 NAR database issue（2025 > 2024 > 2023 > 2021）
2. DOI 必须可验证；不可验证则标注并保留备用引用
3. 作者名以 Nature 原文 PDF 为准（EZSpecificity 等关键论文）
4. 无 formal paper 的资源（RDKit、ESIBank）按软件引用或通过母论文引用
5. 收敛后共 **23** 条（**含 2 条可选历史引用 Berman 2000 + Nelson 2004**）

---

## 一、最终清单（按章节定位排序）

| # | BibTeX key | 论文简写 | 用途 | 冲突说明 |
|---|------------|----------|------|----------|
| 1 | `UniProt2025` | UniProt 2025 NAR | §3.1 酶序列来源 | v1/v3 分别为 2025/2023 版，取 **2025** |
| 2 | `Burley2025` | RCSB PDB 2025 NAR | §3.1 晶体结构来源 | v1 取 2025 版（2023 版亦可用） |
| 3 | `Berman2000` | 原始 PDB 2000 NAR | §3.1 历史基础引用（可选）| 各版本一致 |
| 4 | `Blum2025` | InterPro 2025 NAR | §3.2 P450 功能域验证 | v1 取 2025 版（v2/v3 取 2023 Paysan-Lafosse 版，我们用最新） |
| 5 | `Chang2021` | BRENDA 2021 NAR | §3.1 EZSpecificity 训练集来源讨论 | 各版本一致 |
| 6 | `Nelson2009` | P450 Homepage | §3.1 P450 命名规范 | 各版本一致 |
| 7 | `Wang2021PCPD` | PCPD 数据库 | §3.1 S9 源（中科院天津）| v1/v3 一致，v2 误给 Nelson 2004 |
| 8 | `Nelson2004` | 植物 P450 比较基因组学 | §3.1 S8 PlantP450DB 的命名背书（可选）| v2 版独有，作为 PlantP450DB 源的背景文献 |
| 9 | `Varadi2024` | AlphaFoldDB 2024 NAR | §3.5 AI 预测结构（可选）| v3 取 2024 版（2022 版亦可） |
| 10 | `EZSpec2025` | **Cui 2025 Nature (SOTA)** | **§3.3 + §3.6 + 贯穿全论文** | 各版本 DOI 一致；v1 作者列表与 PDF 一致 ✅ |
| 11 | `UniDock2023` | Uni-Dock 2023 JCTC | §3.5 GPU 对接 | 各版本一致 |
| 12 | `Vina2010` | AutoDock Vina 原始 | §3.5 对接方法溯源 | 各版本一致 |
| 13 | `Vina2021` | AutoDock Vina 1.2 | §3.5 对接方法溯源 | 各版本一致 |
| 14 | `Meeko2025` | Meeko 2025 JCIM | §3.5 配体预处理 | 各版本一致 |
| 15 | `AlphaFill2023` | AlphaFill 辅因子回填 | §3.5 血红素回填（可选）| 各版本一致 |
| 16 | `ECFP2010` | Morgan / ECFP | §3.6 Morgan 指纹 1024-bit | 各版本一致 |
| 17 | `RDKit` | RDKit 软件引用 | §3.6 分子处理 | 各版本一致（无正式论文，用 Zenodo DOI） |
| 18 | `Biopython2009` | BioPython | §3.5 PDB 解析 | 各版本一致 |
| 19 | `ESM2023` | **ESM-2 Lin 2023 Science** | §3.6 酶序列嵌入 | 各版本一致 |
| 20 | `ESM2021` | ESM Rives 2021 PNAS | §3.6 ESM 源头（可选）| 各版本一致 |
| 21 | `GROVER2020` | GROVER Rong 2020 NeurIPS | §3.6 分子嵌入 | 各版本一致（无 DOI，NeurIPS 出版） |
| 22 | `InChI2015` | InChI 2015 J Cheminform | §3.2 底物去重 | 各版本一致 |
| 23 | `SMILES1988` | SMILES 原始 | §3.1 底物表示 | 各版本一致 |

---

## 二、关键决策说明

### 2.1 EZSpecificity 作者列表（以 Nature 原文 PDF 为准）

组会 PDF 明确记载 Nature 2025 论文的完整作者列表：

> **Haiyang Cui, Yufeng Su, Tanner J. Dean, Tianhao Yu, Zhengyi Zhang, Jian Peng, Diwakar Shukla, Huimin Zhao**

v3 版本给出的 "Su, Yue" 与 "Dean, Thomas J." 均为 OCR 错误，本 bib 以 PDF 为准。

### 2.2 ESIBank（A8）

所有版本均确认 ESIBank 无独立发表论文，通过 EZSpecificity 母论文（`EZSpec2025`）引用。我们论文中的措辞建议：

> ……EZSpecificity 论文构建了 ESIBank 酶-底物复合物数据库\cite{EZSpec2025}……

### 2.3 AlphaFoldDB 版本选择

v1 给出 2022 版（50(D1)），v3 给出 2024 版（52(D1)）。本 bib 取 **2024** 版，因为我们实际使用 AlphaFoldDB 时确实使用的是 2023 年以后的结构覆盖规模。不强制；如未使用 AlphaFoldDB 提供的 P450 结构，可删除。

### 2.4 InterPro 版本选择

v1 取 Blum 2025 (53(D1) D444–D456)，v2/v3 取 Paysan-Lafosse 2023 (51(D1) D418–D427)。两者都是 InterPro 官方数据库论文，取最新版 Blum 2025。

### 2.5 Meeko

v1/v3 一致给出 Santos-Martins 2025 JCIM 65(24) 13045–13050。v3 附备注指出若工作完成于 2025 年前，可改引 GitHub 仓库。本项目 2026 年完成，直接引 Meeko 2025 论文。

### 2.6 RDKit 引用方式

RDKit 官方无发表论文，统一采用 software 引用 + Zenodo DOI：

```bibtex
@software{RDKit,
  author = {Landrum, Greg and others},
  title = {RDKit: Open-source cheminformatics},
  url = {https://www.rdkit.org},
  doi = {10.5281/zenodo.591637}
}
```

### 2.7 PlantP50DB（哥本哈根大学 S8 源）

无独立发表论文。本 bib 以 `Nelson2004` 作为背景引用（Nelson 等 2004 年对拟南芥和水稻 727 个 P450 的比较基因组学分析为植物 P450 系统命名奠基）。正文描述为：

> ……PlantP450DB 整合了以 Nelson 等系统化命名为基础的植物 P450 数据\cite{Nelson2004}……

---

## 三、已放入 `reference.bib` 的 BibTeX 条目

以上 23 条已全部追加写入 `毕业设计/毕业论文/reference.bib` 的 "% ===== Batch 1 (P450 论文 Ch3 引文) =====" 区块下。章节正文可直接使用 `\cite{BibTeX key}` 调用。

---

## 四、未来 Batch 2 / Batch 3 规划

- **Batch 2（Ch4/Ch5 用）**：
  - EGNN (Satorras et al. 2021 ICML) — Ch2 §2.3 + Ch3 §3.6 结构通道
  - GVP (Jing et al. 2021 ICLR) — Ch5 §5.3 残基级几何向量感知
  - Transformer / 多头注意力 (Vaswani et al. 2017) — Ch2 §2.3 + Ch3 §3.6
  - AdamW / Loshchilov 2019 — Ch4 §4.1 优化器
  - PyTorch Lightning — Ch4 §4.1 训练框架
  - AUC-ROC / AUPR 在不平衡数据上的语义引用 — Ch4 §4.1

- **Batch 3（Ch1/Ch2 用）**：
  - P450 生物学与药物代谢综述（Guengerich / Ortiz de Montellano）
  - 植物次生代谢 P450 综述（Bak / Nelson / Werck-Reichhart）
  - CLEAN (Yu et al. 2023 Science) — Ch1 §1.2 + Ch2 §2.1
  - ESP (Kroll et al. 2023 Nat Commun) — Ch1 §1.2 + Ch2 §2.1
  - CPI-prediction / DeepEC / GraphBind 等既往基线
  - 中文综述（加分项，Ch1 §1.2 引 1-2 篇）

---

**收敛日期**：2026-04-25
**责任人**：Claude Code + 用户四版 AI 检索协作
**版本**：v1.0
