# GROVER LMDB 对齐 Bug 发现记录

**日期**: 2026-04-14
**发现节点**: EXP002a_fixed 训练中（尚未完成），用户追问 "EXP003 的角度是否是对的酶"，深挖后意外查出 GROVER 通道的第二类对齐 bug
**严重程度**: 高 —— 影响 EXP001/EXP002a/EXP002b/EXP003/EXP003_fixed/EXP002a_fixed 的全部底物 GROVER 嵌入

---

## 一、数据链路完整梳理

本次排查倒逼我把整条特征生成链梳理清楚，作为长期参考。

### 1.1 三份底层词汇表

所有特征生成和训练都基于三份独立的 CSV 词汇表，每行一个唯一实体，**行号本身即为全局 ID**：

| 表 | 行数 | 列 | 唯一 ID 含义 |
|---|---|---|---|
| `Enzymes.csv` | 1,622 | Protein sequence | 行号 = Enzyme Index |
| `Substrates.csv` | 2,125 | Substrate_SMILES | 行号 = Substrate Index |
| `data.csv` | 52,254 | Dock Index, Enzyme Index, Substrate Index, Label | 行号 = Dock Index = 训练样本 ID |

data.csv 每行就是一个训练样本，通过三个 index 引用三份词汇表和对应的结构文件：

```
data.csv 第 i 行
  ├─ Enzyme Index    → Enzymes.csv 某行        → 酶序列
  ├─ Substrate Index → Substrates.csv 某行     → 底物 SMILES
  └─ Dock Index      → complex/{Dock Index}.pdb → 对接复合物
```

### 1.2 五条特征生成管线（Phase 7）

训练需要的特征分 5 类，由 5 个独立脚本生成，各自选择存储后端与 key 规则：

| 通道 | 脚本 | 存储 | 查询索引 | 应有 key | 实际 key | 状态 |
|---|---|---|---|---|---|---|
| **ESM 酶嵌入** | phase7_step2_esm.py | `enzymes.lmdb` | Enzyme Index | Enzyme Index | 压缩计数器 | ❌ **Bug 1** |
| **GROVER 底物嵌入** | phase7_step5_grover.sh | `grover_fingerprint.lmdb` | Substrate Index | Substrate Index | gcsv 行号 | ❌ **Bug 2** |
| **Morgan 指纹** | phase7_step34_graph_morgan.py | `morgan_fingerprint.npy` | Substrate Index | 数组第 N 行 | 数组第 N 行 | ✅ |
| **底物图特征** | phase7_step34_graph_morgan.py | `reaction_features.lmdb` | Substrate Index（校验用）| Substrate Index | 压缩计数器 | ⚠️ **Bug 3**（仅筛选，不污染特征） |
| **结构特征** | phase7_step6_structure.py | `structure_features.lmdb` | Dock Index | Dock Index | Dock Index | ✅ |

### 1.3 训练样本装配（build_pt_cache.py）

构建每个 .pt 样本时按三个 index 并行查五处特征源：

```python
enz_idx  = data_csv[i]["Enzyme Index"]
sub_idx  = data_csv[i]["Substrate Index"]
dock_idx = data_csv[i]["Dock Index"]

esm_emb   = enzymes.lmdb.get(str(enz_idx))            # ❌ Bug 1 → 错酶
grover    = grover.lmdb.get(str(sub_idx))             # ❌ Bug 2 → 错底物
morgan    = morgan.npy[sub_idx]                       # ✅
ligand    = structure.lmdb.get(str(dock_idx))["ligand"]  # ✅
pocket    = structure.lmdb.get(str(dock_idx))["pocket"]  # ✅（含残基角度）
```

打包为 .pt 文件后写入 `pt_cache/`，训练 DataLoader 直接读 .pt —— **所有污染都发生在这一步，后续训练完全无法察觉。**

---

## 二、两次 Bug 的共同根源

Bug 1（ESM）和 Bug 2（GROVER）本质是同一类错误：

> **写 LMDB 时用"顺序计数器"作 key，而非"原始行号"。**
> **只要在写入前跳过过任何一行，计数器就和行号错位。**

| Bug | 代码片段 | 跳行原因 | 跳行数 | 错位起点 |
|---|---|---|---|---|
| ESM Bug 1 | `uniprot_dict[str(idx)] = (len(uniprot_dict), 1)` → `txn.put(str(uniprot_dict[uniprot][0]))` | 44 非标准氨基酸 + 1 长度超限 | 45 | 第一个被跳酶之后 |
| GROVER Bug 2 | `grover_substrates.csv` 手动删掉 `*[H]` 行后喂给 GROVER，GROVER 按 CSV 行顺序写连续 key | `*[H]` 触发 GROVER index-out-of-bounds，处理方式是删行 | 1 | Substrate Index 8 之后 |

正确模式（Morgan 就是反例）：**失败填零不跳行**，或 **用原始行号作 key**。

---

## 三、GROVER Bug 的硬证据

### 3.1 CSV 行数对比

```
Substrates.csv:          2125 行（0 NaN，0 重复）
grover_substrates.csv:   2124 行
```

### 3.2 CSV 逐行对齐

| 行 | Substrates.csv | grover_substrates.csv |
|---|---|---|
| 0-7 | 完全一致 | 完全一致 |
| 8 | `*[H]` | `*c1c(*)c(*)c(C2Oc3...`（= Substrates.csv[9]）|
| 9+ | — | 从 Substrates.csv[10] 开始继续下移 |

结论：`grover_substrates.csv` 就是 `Substrates.csv` 删掉行 8 (`*[H]`) 后的版本。

### 3.3 原子数实证（决定性证据）

对 `grover_fingerprint.lmdb` 全部 2124 个 key，读 `data["embedding"].shape[0]`（chemprop 实际处理的分子原子数），比对三个候选 SMILES 的 `GetNumAtoms()`：

| 对齐假设 | 匹配率 |
|---|---|
| `lmdb[N] == Substrates.csv[N]`（"无 bug"）| 352/2124 = **16.57%**（巧合） |
| `lmdb[N] == Substrates.csv[N+1]`（shift）| 2117/2124 = **99.67%** |
| `lmdb[N] == grover_substrates.csv[N]` | **2124/2124 = 100.00%** |

首 20 key 过渡点：

```
k=0..7:  匹配 orig[N]   ← 正常
k=8..19: 匹配 orig[N+1] ← 从 k=8 起全部错位 1 格
```

断点精确落在 k=8（即 `*[H]` 被删的位置）。

### 3.4 Codex 独立复核

> "Effectively airtight. The 100% match to gcsv[N] is the strongest evidence — `data["embedding"].shape[0]` is produced from the molecule GROVER actually processed. The exact transition at k=8 is the clincher. A chemprop atom count could coincidentally match the wrong SMILES, but that cannot explain all three signals simultaneously."

替代解释已逐一排除：
- `grover_vocab` 不在训练管线中消费（只用于 GROVER 内部 token 化）
- DataLoader `shuffle=False`，顺序与 CSV 一致
- `a_scope` / `id+st` 在 batch 内保持 embedding 与 key 一致，无法恢复原始 Substrate Index

---

## 四、错位结构与影响范围

```
LMDB key    →  grover_substrates.csv row  →  真实 Substrates.csv 行
key "0"     →  row 0                       →  Substrate Index 0   ✓
...
key "7"     →  row 7                       →  Substrate Index 7   ✓
key "8"     →  row 8                       →  Substrate Index 9   ✗
key "9"     →  row 9                       →  Substrate Index 10  ✗
...
key "2123"  →  row 2123                    →  Substrate Index 2124 ✗
(key "2124" 不存在)
```

**训练时的实际后果** —— 代码 `txn.get(str(sub_idx).encode())`：

| Substrate Index | 查到的 GROVER | 数量 | 占比 |
|---|---|---|---|
| 0-7 | 自己的（正确）| 8 | 0.38% |
| 8 | Substrate Index 9 的 | 1 | 0.05% |
| 9-2123 | Substrate Index N+1 的 | 2115 | 99.53% |
| 2124 | LMDB 缺失 → 样本过滤 | 1 | 0.05% |

**99.6% 的底物加载了错误的 GROVER 嵌入。**

---

## 五、Codex 对影响程度的判断

> "Fixing GROVER could help, but I would expect a smaller gain than the ESM fix. This corruption is severe in coverage but milder in semantics than random noise, because adjacent substrates in your dataset may be chemically similar. So the model may have learned to partially tolerate it, especially with Morgan and structure channels still correct. I would expect degradation, not total destruction."

直觉印证：Substrates.csv 前几行肉眼可见是同一甾体骨架的衍生物，邻近底物往往结构相似，shifted-by-1 GROVER 比随机乱配伤害小，但仍是**结构性污染**而非无害噪声。

---

## 六、受影响实验清单

所有实验都经由 `grover_fingerprint.lmdb` 构建 pt_cache，因此全部中招：

| 实验 | ESM 状态 | GROVER 状态 | 备注 |
|---|---|---|---|
| EXP001 | ❌ Bug 1 | ❌ Bug 2 | 旧 baseline |
| EXP002a (Fe/HEM) | ❌ Bug 1 | ❌ Bug 2 | Test 0.7816 |
| EXP002b (LR tuning) | ❌ Bug 1 | ❌ Bug 2 | Test 0.7889 |
| EXP003 (residue geom) | ❌ Bug 1 | ❌ Bug 2 | Test 0.7914 |
| EXP003_fixed | ✅ 已修 | ❌ **Bug 2 未修** | Test 0.8943，**GROVER 仍污染** |
| EXP002a_fixed（训练中）| ✅ 已修 | ❌ **Bug 2 未修** | 本次发现时正在跑 |
| EXP001_fixed（待跑）| — | — | 尚未启动 |

**EXP003_fixed Test AUC=0.8943 的绝对数值仍然是带 GROVER 污染的。但 EXP003→EXP003_fixed 的 +0.1029 增量确实是 ESM 修复单独贡献的。**

---

## 七、修复方案

两条路：

### 方案 A：完整重跑 GROVER（Codex 首选）
- 在完整 `Substrates.csv` 上跑 GROVER
- 显式处理 `*[H]` 崩溃（跳过但写占位 + 映射表）
- 产物 `grover_fingerprint_fixed.lmdb`，key 直接等于 Substrate Index
- 干净、不脆弱，但要加载 GROVER checkpoint、再次面对 `*[H]` 崩溃

### 方案 B：直接 rekey 现有 LMDB（首推）
- 现有 LMDB 每个 key 的 **内容是对的**，只是 **key 贴错标签**
- 新建 `grover_fingerprint_fixed.lmdb`，遍历旧 LMDB：
  - key "k" 的 SMILES = `grover_substrates.csv[k]`
  - 对应真实 `Substrate Index` = k（当 k<8）或 k+1（当 k≥8）
  - 把旧 LMDB key "k" 的 payload 原封不动写入新 LMDB 的真实 Substrate Index key
- 新 LMDB 共 2124 个 key：{0..7, 9..2124}，key "8" 缺失
- 然后用 fix_cache_overlay 模式重建 `pt_cache_heme_fixed` / `pt_cache_geom_fixed`，validity check 会过滤 Substrate Index 8 相关样本（已知解析不了）
- **纯文件层面 rekey，几十秒完成，无需加载 GROVER 模型**

**推荐方案 B**：已经 empirical 证明内容正确，仅 key 映射错，无重跑必要。

---

## 八、pt_cache 层的连锁影响

修好 GROVER LMDB 还不够，下游 cache 需要同步重建：

```
grover_fingerprint.lmdb (错 key)
    ↓
pt_cache/ (flatbin + per-sample .pt)       ← 所有实验共享
pt_cache_heme/ (EXP002a/b 的 Fe/HEM overlay)
pt_cache_geom/ (EXP003 的残基几何 overlay)

↑ 三个 cache 都以旧 grover LMDB 为 upstream，必须在修好 LMDB 后重建对应的 fixed 版本
```

**需要重建**:
- `pt_cache_fixed/` （shared enzyme + substrate GROVER + substrate Morgan + complex graph）
- `pt_cache_heme_fixed/` （overlay 上 Fe/HEM 后的版本）
- `pt_cache_geom_fixed/` （overlay 上残基几何后的版本）

沿用 EXP003_fixed 的 fix_cache_overlay 非破坏式模式：旧 cache 保留作对照，新 cache 通过 symlink + filter 生成，成本较低。

---

## 九、重跑计划

修复后必须重跑 fixed baseline 系列才能得到可发表的数据：

1. **修 GROVER**（方案 B，~30 秒）
2. **rebuild pt_cache_fixed / pt_cache_heme_fixed / pt_cache_geom_fixed**（fix_cache_overlay 改版，~30 分钟）
3. **重跑 EXP001_fixed / EXP002a_fixed / EXP003_fixed**（训练配方仍沿用 EXP003_fixed：lr=4e-4/warmup=12/wd=1e-5/bs=88，4×RTX5090，每个 ~1-1.5 小时）
4. 消融链应为：EXP001_fixed（feature_dim=28）→ EXP002a_fixed（31，加 Fe/HEM）→ EXP003_fixed（37，加残基角度）

EXP003_fixed Test AUC=0.8943 会被一个"全部 fixed"版本替代，可能更高也可能差不多——取决于 GROVER 在 P450 数据里的实际贡献。

---

## 十、决策点（待用户裁决）

1. **是否立即停 EXP002a_fixed？** 它正在 GROVER 污染的数据上训练，继续跑出来的权重仍需废弃
2. **修复用方案 A 还是 B？** 推荐 B（纯 rekey）
3. **重跑顺序？** 建议 EXP001_fixed → EXP002a_fixed → EXP003_fixed 依次

---

## 十一、教训与规则

1. **写 LMDB 必须用原始行号作 key，禁止用压缩计数器。**
2. **跳过失败行时必须记录映射或填充占位，不得静默删行。**
3. **所有 feature LMDB 在使用前必须做 key 与词汇表的 1:1 对齐校验。**
4. **发现一类 bug 后必须全链路审计同类型写入逻辑**——Bug 1 和 Bug 2 结构完全一致却间隔发现，说明当时修 ESM 时未做系统化审计。
5. **Morgan 的"失败填零"模式是值得全链路推广的正确模式。**
