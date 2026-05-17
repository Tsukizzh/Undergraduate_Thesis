# Q9 — 底物分类 + 下游预测（Stage 1 模型）

## 老师原话

> 已经把底物进行大致分类，后续尝试酶能否催化哪类底物

## 状态

🟡 **部分完成**：分类已做（C3 v6 FINAL）；下游模型未做

## 实施位置

> ⚠️ "分类部分"在 `PathC_.../C3_P450专属模型训练/C3 底物分类` 已完成。
> 本目录跟踪"下游模型"部分。

## 已完成（C3 v6 FINAL，2026-03-27 ~ 31）

### 数据规模
- 2,125 化合物 → **1,870 confirmed (88.0%) + 255 other (12.0%) + 63 multi-label**

### 8 类多标签

| 类别 | 数量 |
|---|---:|
| Terpenoid 萜类 | 484 |
| Amino_acid 氨基酸 | 388 |
| Fatty_acid 脂肪酸 | 278 |
| Alkaloid 生物碱 | 251 |
| Steroid 甾体 | 211 |
| Phenylpropanoid 苯丙素 | 176 |
| Polyketide 聚酮 | 137 |
| Other 其他 | 255 |

### 流水线
P1 Gold → P2 NPC Superclass → P3 Pathway → P4 AA SMARTS → P5 Other + 9 corrections + 352 agent review

### 输出
- `substrate_multilabel_FINAL.csv`
- 分类器：`classify_multilabel.py`

## 待做（Stage 1 下游模型）

### 任务定义
- **输入**：酶序列（ESM 嵌入）
- **输出**：8 维多标签向量（sigmoid + BCEWithLogitsLoss）
- **任务**：判断"该酶能催化哪类底物"

### 与原 Stage 2 模型（EXP001~005）的关系
- **Stage 2**：酶 + 底物 → 二分类（特异性预测）
- **Stage 1**：酶 → 底物类别（粗粒度筛选）
- **联合用法**：Stage 1 先粗筛酶可能催化的类别，Stage 2 在该类别内做精细打分

## 架构候选

### Plan A：纯序列 MLP
- ESM-2 (1280) → MLP → 8 dim sigmoid
- 简单基线

### Plan B：复用 EXP001 编码器
- 复用已训练的 P450 EXP001 ckpt 的酶 encoder（128 dim）
- → 新加 8 dim 多标签头
- 参数高效，可单 GPU 快速训练

### Plan C：联合训练
- Stage 1 + Stage 2 共享酶 encoder，多任务损失
- 可能互相增强

## 训练数据准备

### 多标签构造
- 对每个酶，统计其在正样本对中出现过的所有底物，回查每个底物的 8 类标签
- 取 OR 聚合 → 该酶的多标签向量
- 例：酶 X 催化 4 个 Terpenoid + 2 个 Alkaloid → 标签 = [1, 0, 0, 1, 0, 0, 0, 0]

### 切分
- 沿用 EXP006 cluster split（避免家族泄露）

## 工作量估计

- Plan A: 1–2 天
- Plan B: 2–3 天
- Plan C: 5–7 天

## 待办

- [ ] 确认 Stage 1 是否纳入毕业论文核心贡献
- [ ] 选定架构（A/B/C）
- [ ] 准备酶-类别多标签训练数据
- [ ] 训练 + 评估
- [ ] 设计与 Stage 2 的联合分析（Recall@K / 类别先验对 Stage 2 的提升）

## 与其他问题的关联

- **Q6 + Q12**：Stage 1 的训练标签依赖清洗后的正样本
- **Q2（EXP006）**：cluster split 是 Stage 1 评估的合适基准

## 变更日志

| 日期 | 变更 |
|---|---|
| 2026-05-08 | session 创建（指向 PathC/C3 分类已完成部分） |

