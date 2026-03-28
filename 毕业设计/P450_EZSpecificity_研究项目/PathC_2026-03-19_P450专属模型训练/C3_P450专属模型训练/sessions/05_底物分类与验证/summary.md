# Step 5: 底物分类与多轮验证

> **日期**: 2026-03-27 ~ 2026-03-29
> **状态**: 多轮验证完成，审计准确率 88.5%，待合并 8 类+修正后重新审计

## 做了什么

对 P450 数据集中 2,125 个底物进行化学类别分类，经历 4 大阶段验证。

## 四个阶段

### 阶段 1: NPClassifier + Codex 修正（3-27）
- NPClassifier API 自动分 15 类
- Codex 4 轮讨论设计修正规则（结构否决+重推导+救回）
- 修正 275 个化合物

### 阶段 2: 多源验证管线（3-28）
- 结构 SMARTS 验证器（15 类）→ Tier 1-4 分层
- ClassyFire 批量查询 → 855/2,110 命中
- 共识引擎：Tier 1(541) + Tier 2(794) + Tier 3(512) + Tier 4(278)

### 阶段 3: Agent 文献验证（3-28~29）
- **Round 1**: 15 Agent 并行验证 450 个高风险化合物，~480 次 web 搜索 → 190 个重分类
- **Round 2**: 10 Agent 并行验证 341 个 Unclassified，~340 次 web 搜索 → 50 个救回
- **二次验证**: 4 Agent 复查 106 个 medium confidence → 22 个分歧由 Codex 仲裁

### 阶段 4: 准确率审计（3-29）
- 5 Agent 并行审计 200 个分层随机抽样
- **结果：177/200 = 88.5%**
- 23 个错误主要集中在 Alkaloid/Amino_acid 边界（7个）和 Unclassified 过于保守（9个）

## 总工作量

| 项目 | 数量 |
|------|------|
| Opus Agent | ~34 个 |
| Web 搜索 | ~1,820 次 |
| Codex 讨论 | ~8 轮 |
| 覆盖率 | 2,125/2,125 = 100% |

## 当前 15 类分布

| 类别 | 数量 | | 类别 | 数量 |
|------|------|-|------|------|
| Unclassified | 398 | | Sesquiterpenoid | 98 |
| Alkaloid | 295 | | Triterpenoid | 97 |
| Amino_acid | 261 | | Polyketide | 95 |
| Steroid | 211 | | Monoterpenoid | 79 |
| Fatty_acid | 203 | | Phenylpropanoid | 79 |
| Diterpenoid | 137 | | Terpenoid_other | 66 |
| | | | Flavonoid | 50 |
| | | | Macrolide | 33 |
| | | | Coumarin | 23 |

## 审计发现的核心问题

**问题不是"分类太细"，而是"边界规则不统一"**——从未正式定义每个类别的精确边界。

## 待做

1. 合并到 8 类（甾体/萜类/生物碱/氨基酸/脂肪酸/苯丙素广义/聚酮广义/其他）
2. 按正式分类手册修正 + 重新审计
3. 进入阶段 1 模型训练（酶序列→底物类别预测）

## 文件

- 最终分类: `data/05_底物分类/substrate_classifications_FINAL.csv`
- 分类手册: Codex 制定，见 session_log 第十九节
- 创新性分析: `C3_创新性分析.md`
