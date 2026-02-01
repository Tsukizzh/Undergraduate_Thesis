# P450-EZSpecificity 研究项目（初步探索阶段仓库）

基于 EZSpecificity 模型的 P450 酶特异性预测研究。

## 项目背景

EZSpecificity 是一个基于交叉注意力机制的 SE(3)-等变图神经网络，用于预测酶-底物特异性（Nature 2025）。本项目专注于 P450 细胞色素酶家族的特异性研究。

## 目录结构

```
EZSpecificity_Project/
├── src/                                    # EZSpecificity源代码
│   ├── Models/                             # 模型实现
│   ├── Datasets/                           # 数据处理
│   ├── example.ipynb                       # 推理示例
│   └── environment.yml                     # Conda环境
│
├── P450_EZSpecificity_研究项目/             # 核心研究目录（四路径计划）
│   ├── PathA_2026-01-08_模型评估测试集构建/
│   │   ├── data/                           # 数据文件（Step 1-5产出）
│   │   │   ├── 01_Step1_PDB文件/           # 627个PDB结构
│   │   │   ├── 02_Step2_酶序列/            # 292条唯一酶序列
│   │   │   ├── 03_Step3_配体处理与数据集/   # 6套数据集方案
│   │   │   ├── 04_Step4_CSV格式修正/       # 格式化后CSV
│   │   │   └── 05_Step5_ESM酶嵌入/         # enzyme_features.lmdb
│   │   ├── scripts/                        # 执行脚本
│   │   ├── sessions/                       # 执行日志
│   │   ├── reports/                        # 分析报告
│   │   ├── source_data/                    # 源数据（682条）
│   │   └── 进度日志.md
│   ├── P450研究四路径计划.md
│   └── 全局进度日志.md
│
├── 提取P450过程日志/                        # 研究过程记录（按时间排序）
│   ├── 2025-12-31_01-30_仓库结构扫描/
│   ├── 2025-12-31_02-00_数据来源分析/
│   ├── 2025-12-31_02-10_数据结构详解/
│   ├── 2025-12-31_02-20_P450识别结果/
│   ├── 2026-01-02_01-46_P450精确验证/
│   ├── 2026-01-02_23-00_底物数据整合/
│   ├── 2026-01-03_03-10_P450物种分类/
│   ├── 2026-01-03_05-30_InterPro验证补充/
│   ├── 2026-01-04_01-00_P450_PDB实验结构_ML数据集构建/
│   ├── 2026-01-21_实验二v2_数据重构/       # 按(UniProt+配体)去重
│   └── 项目进度日志.md
│
├── P450_EZSpecificity完整研究手册_终极整合版.md
├── 路径选择与核心概念FAQ.md
├── CLAUDE.md                               # Claude Code工作指南
└── README.md                               # 本文件
```

## 关键数据

| 数据集 | 数量 | 说明 |
|--------|------|------|
| ESIBank 中的 P450 | 389 个酶 | 训练集中已有的 P450 酶 |
| 独立测试集（v2） | 682 条酶-配体对 | 153个唯一P450，627个PDB |
| 配体分类完成 | 682 条 | SUBSTRATE/INHIBITOR/PRODUCT/EXCLUDE |
| 非P450污染排除 | 58 条 | Photosystem I（740→682） |

## 当前进度

- [x] P450 酶识别与验证
- [x] 底物数据整合
- [x] PDB 实验结构收集
- [x] 四路径研究计划制定
- [x] 路径A Step 1-2：数据收集与序列提取（v2完成，682条）
- [x] 路径A Step 3：配体处理与数据集构建（6套方案完成）
- [x] 路径A Step 4：CSV格式修正（列名/列顺序/编码统一）
- [x] 路径A Step 5：ESM酶嵌入生成（292条→enzyme_features.lmdb）
- [ ] 路径A Step 6：底物特征生成（反应图+Morgan+GROVER）🚀 **待执行**
- [ ] 路径A Step 7-10：结构特征+模型推理+结果分析
- [ ] 路径B：P450 专属数据集构建
- [ ] 路径C：P450 专属模型训练
- [ ] 路径D：区域选择性预测

## 参考文献

Cui, Y., et al. (2025). Enzyme specificity prediction using cross-attention graph neural networks. *Nature*.
