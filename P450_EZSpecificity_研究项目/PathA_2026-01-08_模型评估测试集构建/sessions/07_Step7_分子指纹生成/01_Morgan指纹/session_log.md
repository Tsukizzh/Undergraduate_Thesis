# Step 7.1: Morgan指纹生成 - 详细会话日志

> **执行日期**: 2026-01-30
> **版本**: v1.0
> **执行者**: Claude Code + Codex + Gemini (多模型协作)

---

## 概述

**为全部436个底物生成1024位Morgan指纹（radius=2），产出numpy数组供模型分子特征流水线使用。**

---

## 背景知识

### 什么是Morgan指纹？

Morgan指纹（也称为ECFP - 扩展连接指纹）是一种将分子结构编码为固定长度二进制向量的方法。

**类比**：想象用放大镜观察一个分子。对于每个原子，你检查它在特定半径范围内的"邻域"。原子及其周围环境形成一个独特的模式。Morgan指纹收集分子中所有这样的模式，并将它们哈希到固定长度的位向量中。

### 关键参数

| 参数 | 值 | 含义 |
|------|----|----|
| radius | 2 | 查看每个原子周围2个键的范围（相当于ECFP4）|
| nBits | 1024 | 输出向量长度（1024个二进制位）|
| dtype | int8 | 每个位存储为0或1 |

### 为什么是1024位？

EZSpecificity模型配置中指定了 `features: ["morgan", 1024]`。这意味着模型期望输入恰好是1024维的Morgan指纹作为其输入特征之一。

---

## 输入

### 文件: `Substrates.csv`

| 字段 | 描述 | 示例 |
|------|------|------|
| Substrate_Index | 底物ID（从0开始）| 0, 1, 2, ..., 435 |
| Substrate_SMILES | 分子SMILES表示 | COc1ccccc1O |

**总计: 436个底物**

---

## 输出

### 文件: `morgan_fingerprint.npy`

| 属性 | 值 |
|------|-----|
| 形状 | (436, 1024) |
| 数据类型 | int8 |
| 文件大小 | 436.1 KB |
| 索引对齐 | 第i行 = Substrate_Index i |

### 抽样检查（前3个底物）

| Substrate_Index | SMILES（截断）| 置位数 |
|-----------------|--------------|--------|
| 0 | COc1ccccc1O | 16 |
| 1 | COc1cccc(OC)c1O | 19 |
| 2 | CN[C@H](C(=O)N[C@H](CO)Cc1c[nH]c2ccccc12)C(C)C | 56 |

### 统计摘要

| 指标 | 值 |
|------|-----|
| 全零行数 | 0（所有底物均已处理）|
| 平均置位数 | 37.2 |
| 最小置位数 | 3 |
| 最大置位数 | 105 |

---

## 模型如何加载Morgan指纹

来自 `src/Datasets/data_representer.py`:

```python
# 加载（第112行）
self.morgan_dbs = [np.load(path).astype(np.float32) for path in ...]

# 访问（第189行）
substrate_features["morgan"] = torch.from_numpy(
    self.morgan_dbs[dataset_id][substrate_index][np.newaxis, :]
)
# 结果: 形状为(1, 1024)的tensor，数据类型float32
```

**关键点**：行索引必须等于Substrate_Index。我们的脚本在生成指纹前按Substrate_Index排序，确保正确对齐。

---

## 执行

### 命令
```bash
D:\anaconda3\envs\torch\python.exe step7_1_generate_morgan.py
```

### 结果
```
Morgan Fingerprint Generation Complete
============================================================
Input:         .../04_Step4_格式修正后数据/Substrates.csv
Output:        .../07_Step7_分子指纹/morgan_fingerprint.npy
Substrates:    436
Fingerprint:   radius=2, nBits=1024
Shape:         (436, 1024)
Dtype:         int8
File size:     436.1 KB
All-zero rows: 0
Avg bits set:  37.2
Min bits set:  3
Max bits set:  105
============================================================

Elapsed: 0.2s
```

### 硬件
- 仅使用CPU（无需GPU）
- 处理时间：0.2秒

---

## 多模型审查

| 轮次 | 审查者 | 发现 | 结果 |
|------|--------|------|------|
| 1 | Codex | 形状/数据类型/索引对齐验证通过 | 批准 |
| 2 | Codex | 确认 np.load().astype(float32) 兼容性 | 批准 |

---

## 验证清单

| 检查项 | 结果 |
|--------|------|
| 形状匹配 (436, 1024) | 通过 |
| 数据类型为 int8 | 通过 |
| 无全零行 | 通过 |
| 连续索引 0..435 | 通过 |
| 文件可被 np.load() 读取 | 通过 |
| 与模型配置一致 | 通过 |

---

## 工具与环境

| 工具 | 版本 | 用途 |
|------|------|------|
| Python | 3.12.5 (torch环境) | 运行时 |
| RDKit | 2024.03+ | SMILES解析、指纹生成 |
| NumPy | 1.26.4 | 数组存储 |
| Pandas | 2.2.2 | CSV读取 |

---

## 文件位置

```
PathA_2026-01-08_模型评估测试集构建/
├── data/
│   ├── 04_Step4_格式修正后数据/
│   │   └── Substrates.csv                    ← 输入（436个底物）
│   └── 07_Step7_分子指纹/
│       └── morgan_fingerprint.npy            ← 输出 (436, 1024) int8
└── scripts/
    └── 07_Step7_分子指纹生成/
        └── step7_1_generate_morgan.py        ← 脚本
```

---

**文档版本**: v1.0
**最后更新**: 2026-01-30
