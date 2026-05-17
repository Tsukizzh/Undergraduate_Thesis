# Q2 — 序列相似度划分

## 老师原话

> 按照序列相似度进行划分

## 状态

🟢 **进行中（PathC/EXP006）**

## 实施位置

> ⚠️ 本问题的实际工作目录在
> `PathC_2026-03-19_P450专属模型训练/C3_.../experiments/EXP006_cluster_split_allfix_unified/`
> 本目录仅作为 PathD 索引，**不重复执行**。

## 当前进度（截至 2026-05-08）

| 阶段 | 状态 | 备注 |
|---|---|---|
| 骨架搭建 | ✅ | 从 EXP001_allfix_unified 模板复制 |
| FASTA 导出 | ✅ | 1622 酶 → 1609 唯一序列 / 11 重复组 / 24 酶 |
| MMseqs2 安装 | ✅ | `/root/miniconda3/bin/mmseqs` |
| 4 阈值聚类 | ✅ | id80 / id60 / id40 / id30 |
| TSV 解析为标准 CSV | ✅ | `data/clusters/mmseqs/*/clusters_*.csv` |
| 7 张审计图 | ✅ | fig0–fig6 |
| Cluster split 生成 | ⏳ | 待选定阈值后启动 |
| pt_cache 重建 | ⏳ | run_train.sh 中 CACHE 为 TODO 占位 |
| 训练 | ⏳ | 配方沿用 EXP003_fixed |
| 评估 | ⏳ | — |

## 聚类摘要（1622 酶）

| 阈值 | clusters | 最大簇 | singletons |
|---|---:|---:|---:|
| id80 | 1139 | 19 | 900 |
| id60 | 758 | 38 | 496 |
| id40 | 328 | 148 | 162 |
| id30 | 140 | 481 | 65 |

## 待澄清点

1. **选哪个阈值做主实验？**
   - id40 平衡选项（328 簇 / 162 单簇）
   - id30 严格选项（140 簇 / 最大 481）
2. **样本基础**：复用 EXP001_allfix_unified 的 unified intersection（22083/11008/11000）还是按 1622 酶 × 全量 substrate 重新生成？
3. **是否一次性跑 4 个阈值**（id80/60/40/30）做完整泛化曲线？

## 待办

- [ ] 选定主阈值
- [ ] 编写 cluster split 生成脚本
- [ ] 重建 pt_cache
- [ ] DDP 训练
- [ ] 与 random_split / all_split 结果对比 → 输出"序列相似度 vs AUC"曲线

## 与其他问题的关联

- **Q12（负样本杂泛性）**：cluster split 把同一簇的酶统一分到 train/val/test 一侧，能否也统一处理同簇的"假负样本"？

## 变更日志

| 日期 | 变更 |
|---|---|
| 2026-05-08 | session 创建（指向 PathC/EXP006 实施位置） |

