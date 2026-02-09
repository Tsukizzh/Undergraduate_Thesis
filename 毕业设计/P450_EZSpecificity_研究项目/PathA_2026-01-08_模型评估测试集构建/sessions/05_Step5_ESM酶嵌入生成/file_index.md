# Step 5 文件索引

> **更新时间**: 2026-01-30
> **状态**: ✅ 已完成

---

## 执行摘要

| 项目 | 数值 |
|------|------|
| 输入酶数量 | 292 |
| 成功生成 | 292 (100%) |
| 输出文件 | enzyme_features.lmdb |
| 文件大小 | ~700 MB |
| 处理时间 | ~2-3分钟 |
| GPU使用 | RTX 4070 SUPER |

---

## 脚本文件

| 文件 | 用途 |
|------|------|
| [step5_generate_esm_embeddings.py](../../scripts/05_Step5_ESM酶嵌入生成/step5_generate_esm_embeddings.py) | ESM酶嵌入生成主脚本 |

---

## 输入数据

| 文件 | 来源 | 说明 |
|------|------|------|
| [Enzymes.csv](../../data/04_Step4_格式修正后数据/Enzymes.csv) | Step 4 | 292条酶序列 |

---

## 输出数据

| 文件 | 记录数 | 说明 |
|------|--------|------|
| [enzyme_features.lmdb](../../data/05_Step5_ESM酶嵌入/enzyme_features.lmdb) | 292 | ESM-2蛋白质嵌入 |

---

## LMDB格式说明

### Key-Value结构

| 字段 | 类型 | Shape | 说明 |
|------|------|-------|------|
| key | bytes | - | `str(Enzyme_Index).encode()` |
| embedding | np.ndarray | (L, 1280) | per-residue ESM-2嵌入 |
| active_site | int | - | 活性位点标记（默认1） |
| sequence | np.ndarray | (L,) | 数字编码序列 |

### 数字编码映射

与 `src/Datasets/const.py:letter_to_num` 一致：

```python
{'C': 4, 'D': 3, 'S': 15, 'Q': 5, 'K': 11, 'I': 9,
 'P': 14, 'T': 16, 'F': 13, 'A': 0, 'G': 7, 'H': 8,
 'E': 6, 'L': 10, 'R': 1, 'W': 17, 'V': 19,
 'N': 2, 'Y': 18, 'M': 12, 'X': 20, 'Z': 21, 'U': 22, 'B': 23}
```

---

## 与其他步骤的关联

```
Enzymes.csv (Enzyme_Index)
    ↓ Step 5
enzyme_features.lmdb (key = Enzyme_Index)
    ↓ 模型推理
data.csv (Enzyme Index) → 查找对应嵌入
```

---

## 技术规格

| 项目 | 值 |
|------|-----|
| ESM模型 | esm2_t33_650M_UR50D |
| 参数量 | 650M |
| 输出维度 | 1280 |
| 提取层 | 第33层 |
| 序列长度限制 | < 1000 aa |

---

**文档版本**: v1.0
