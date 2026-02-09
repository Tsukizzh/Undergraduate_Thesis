# Step 9: 模型推理 Session Log

## 基本信息

| 项目 | 内容 |
|------|------|
| **任务** | 使用预训练 EZSpecificity 模型对 P450 测试集进行推理 |
| **日期** | 2026-01-31 |
| **版本** | v1.0 |
| **状态** | ✅ 完成 |

## 执行概要

### 输入数据
- **测试集**: 539 条酶-底物对 (data.csv)
- **有效样本**: 517 条 (经 high_quality_id.txt 过滤)
- **特征来源**:
  - 酶嵌入: enzyme_features.lmdb (Step 5)
  - 反应图: reaction_features.lmdb (Step 6)
  - GROVER指纹: grover_fingerprint.lmdb (Step 7)
  - Morgan指纹: morgan_fingerprint.npy (Step 7)
  - 结构特征: structure_features.lmdb (Step 8)

### 输出结果
- **预测文件**: `data/09_Step9_模型推理/predictions.csv`
- **记录数**: 517 条
- **输出列**: Dock Index, Enzyme Index, Substrate Index, Label, score, logit, prob

## 技术细节

### Windows 兼容性修复

原始代码存在以下 Windows 不兼容问题：

| 问题 | 解决方案 |
|------|---------|
| `set_sharing_strategy('file_descriptor')` - Linux only | 直接移除 |
| DDP 策略在 Windows 崩溃 | 使用纯 forward pass 替代 trainer.test() |
| `num_workers=16` 在 Windows 卡死 | 设置 `config.num_cpus = 0` |
| `full_data=True` 需要不存在的 str_features.lmdb | 设置 `full_data = False` |
| numpy int32 → torch.one_hot 失败 | 添加 `.long()` 类型转换 |

### 源代码修改

修改 `src/Datasets/Structure/transforms.py`:

```python
# Line 32: 添加 .long() 转换
amino_acid = F.one_hot(data.protein_atom_to_aa_type.long(), num_classes=self.max_num_aa)

# Line 135: 添加 .long() 转换
ligand_bond_feature_real = F.one_hot((data.ligand_bond_type - 1).long(), num_classes=6)
```

**原因**: Windows 上 `numpy.dtype=int` 是 32 位，而 `torch.one_hot()` 要求 64 位整数。

### 依赖安装

```bash
pip install pytorch-lightning==1.9.0
pip install warmup_scheduler
pip install ray
pip install cython
```

## 推理结果统计

### 基本统计

| 指标 | 值 |
|------|---|
| 总预测数 | 517 |
| 正样本 (Label=1) | 261 |
| 负样本 (Label=0) | 256 |
| 处理时间 | ~4 分钟 |

### Logit 分布

| 百分位 | 值 |
|--------|------|
| Min | -21.60 |
| P10 | -8.77 |
| P25 | -6.27 |
| P50 (中位数) | -3.57 |
| P75 | -1.46 |
| P90 | -0.16 |
| Max | 2.82 |

### Probability 分布

| 指标 | 值 |
|------|------|
| Min | 0.0000 |
| Max | 0.9435 |

## 多模型协作

### 第一轮: Codex 技术分析
- SESSION_ID: `019c10eb-d5d4-7a02-9227-1e377e3057d7`
- 发现: DDP 兼容性、num_workers 问题、full_data 需求

### 第二轮: Gemini 方案审核
- SESSION_ID: `9b866390-e260-4482-aad6-3afdabac8605`
- 确认: 纯 forward pass 方案可行，`shuffle=False` 保证对齐

### 第三轮: Codex 代码审核
- 确认脚本结构正确
- 建议添加空数据集边界处理

## 运行环境

| 组件 | 版本/规格 |
|------|----------|
| GPU | NVIDIA GeForce RTX 4070 SUPER (12.9 GB) |
| CUDA | 12.1 |
| Python | 3.12 (D:/anaconda3/envs/torch) |
| PyTorch | 2.3.0 |
| pytorch-lightning | 1.9.0 |

## 文件清单

### 生成文件
- `scripts/09_Step9_模型推理/step9_inference.py` - 推理脚本
- `data/09_Step9_模型推理/predictions.csv` - 预测结果

### 修改文件
- `src/Datasets/Structure/transforms.py` - 添加 Windows int64 兼容修复
