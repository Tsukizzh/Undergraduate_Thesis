# Q7 — 原论文是否剔除 Heme 排查

## 老师原话

> 查验原论文是否把 Heme 给剔除掉

## 状态

✅ **已完成（2026-05-16）** — Claude + Codex 双方独立审计交叉验证一致

## 最终结论（三方判定）

| 问题 | 答案 | 关键证据 |
|---|---|---|
| **Q-A** 原模型蛋白侧 GNN 能否看见 HEM/HETATM？ | ❌ **看不见** | `protein_ligand.py:67` 严格 `line[0:6].strip() == 'ATOM'` 过滤，HETATM 在解析时被静默丢弃 |
| **Q-B** 蛋白侧原子词汇表是否含 Fe（Z=26）？ | ❌ **不含** | `transforms.py:23` `[1, 6, 7, 8, 16, 34]` 仅 H/C/N/O/S/Se，6 元素 |
| **Q-C** 是否有任何代码路径把 HEM-Fe 信息送入蛋白/结构通道？ | ❌ **完全没有** | 全代码库审计 brenda.py / structure.py / create_features.py / example.ipynb / main_testing.py 均无 |

**配体侧词汇表含 Fe（Z=26）**，但 P450 底物绝大多数无铁，所以这条通道在 ESIBank 训练数据中几乎不被激活。

## 端到端数据流（修正后）

```
对接 PDB（含 HEM 的 HETATM 行）
   ↓
preprocess.ipynb / example.ipynb cell 14   ── 口袋提取
   - ATOM 行：距离配体 < 10Å + 非 H 原子 → 写入 pocket.pdb
   - HETATM 行：⭐ 无条件写入 pocket.pdb（含 HEM-Fe 等全部辅因子）
   - 文档明确："Any cofactor-related atoms must begin with HETATM" ← 原作者主观上想保留
   ↓
pocket.pdb（HEM 还在里面！）
   ↓
PDBProtein._enum_formatted_atom_lines  (protein_ligand.py:67)
   - 严格 line[0:6].strip() == 'ATOM'
   - ⛔ HETATM 在这里被静默丢弃
   ↓
protein_dict（HEM 已消失）
   ↓
FeaturizeProteinAtom  (transforms.py:23)
   - 词汇表 [H, C, N, O, S, Se] —— 无 Fe
   - 即使 HEM-Fe 侥幸进来也会变成全零 one-hot
   ↓
GNN（永远看不到 HEM 和 Fe）
```

## 关键证据片段

### 证据 ① — 蛋白原子词汇表（`src/Datasets/Structure/transforms.py:19-35`）

```python
class FeaturizeProteinAtom(object):
    def __init__(self):
        self.atomic_numbers = torch.LongTensor([1, 6, 7, 8, 16, 34])    # H, C, N, O, S, Se
        self.max_num_aa = 21
```

**Fe (Z=26) 不在内**。

### 证据 ② — 配体原子词汇表（`src/Datasets/Structure/transforms.py:38-42`）

```python
class FeaturizeLigandAtom(object):
    def __init__(self):
        self.atomic_numbers = torch.LongTensor([1, 6, 7, 8, 9, 15, 16, 17, 26, 35, 53])
        # H, C, N, O, F, P, S, Cl, Fe, Br, I
```

**配体侧有 Fe**，但 P450 底物几乎无铁。

### 证据 ③ — 氨基酸词汇表（`src/Datasets/Structure/protein_ligand.py:20-29`）

```python
AA_NAME_SYM = {
    'ALA': 'A', 'CYS': 'C', ..., 'TRP': 'W', 'TYR': 'Y',
    'UNK': 'X'
}
AA_NAME_NUMBER = {k: i for i, (k, _) in enumerate(AA_NAME_SYM.items())}
```

仅 20 标准 AA + UNK = 21 类，**HEM 不在字典里**。即使 HEM 残基侥幸进入解析层也会被打成 UNK：

```python
# protein_ligand.py:111-114
if atom['res_name'] not in self.AA_NAME_NUMBER:
    self.atom_to_aa_type.append(self.AA_NAME_NUMBER['UNK'])
```

### 证据 ④ — PDB 解析器只接受 ATOM 记录（`src/Datasets/Structure/protein_ligand.py:65-94`）

```python
def _enum_formatted_atom_lines(self):
    for line in self.block.splitlines():
        if line[0:6].strip() == 'ATOM':   # ⛔ HETATM 不匹配
            ...
            yield {...}
        elif line[0:6].strip() == 'HEADER':
            yield {...}
        elif line[0:6].strip() == 'ENDMDL':
            break
```

严格 `'ATOM'` 匹配，HETATM 被静默丢弃。

### 证据 ⑤ — 口袋提取的 bug-or-intent 双重性（`src/Datasets/Structure/preprocess.ipynb` cell 2 + `src/example.ipynb` cell 14）

```python
for line in protein_lines:
    if "ATOM" in line:
        if 'H' not in line[12:16].strip() and distance(...) < 10:
            fin.write(line)
    elif "HETATM" in line or "ENDMDL" in line:
        fin.write(line)   # ⭐ 无条件写入，HEM 全部进 pocket.pdb
```

**关键事实**："HETATM" 不含 "ATOM" 子串（HET-A-T-M 第三位是 T 不是 O），所以 HETATM 走 elif 分支，**被无条件写入 pocket.pdb**（没有距离过滤）。

### 证据 ⑥ — 数据加载层 transforms 流水线（`src/Datasets/brenda.py:31-34`）

```python
transform = Compose([
    utils_trans.FeaturizeProteinAtom(),
    utils_trans.FeaturizeLigandAtom(),
    utils_trans.EdgeConnection(...)
])
```

数据加载层无 cofactor 特殊处理，全部依赖上面的 transforms。

### 证据 ⑦ — 序列通道（ESM）也不携带辅因子信息（`src/Datasets/create_features.py:314-339`）

```python
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
...
data = {
    'embedding': embedding,
    'active_site': uniprot_dict[uniprot][1],
    'sequence': convert_protein_sequence_to_number(seq)
}
```

ESM 通道只看原始氨基酸序列，不含 HEM/Fe 上下文。

## 意外发现 —— 原论文的实现层 bug

原作者**主观上想保留辅因子**（example.ipynb 文档明确写道："Any cofactor-related atoms must begin with the standard PDB identifier 'HETATM'. You don't have to contain cofactor in the docking file."），并且 pocket 提取脚本**确实把 HEM 写入了 pocket.pdb**。

**但下游 `PDBProtein` 解析器用了严格的 `'ATOM'` 匹配，把 HETATM 全部丢弃**。

这是**实现层面的 intent vs implementation gap**，不是有意的设计选择。
对论文叙事价值：揭示原模型一个隐性缺陷 —— **P450 子集的 Fe 催化中心信号在不知不觉中被丢失**。

## 对 Q1 实验设计的直接影响

✅ **Q1 是个清晰、有意义的实验**。在 ESIBank 全集上训练"加 Fe/HEM 嵌入"版本对比原版，预期：
- 加 Fe 后跨家族区分能力↑（因为含铁酶 vs 不含铁酶有了显式信号）
- P450 子集内 A ≈ B（因为 P450 内部 Fe 是常量、无信息量）

**改造点（最小化，与 EXP002a 在 P450 子集上做过的改动一致）**：
1. `protein_ligand.py:67` 改 `'ATOM'` → 允许 `'HETATM'`
2. `protein_ligand.py:20` `AA_NAME_SYM` 加 HEM 及其变体（HEM/HEC/HEB/HEA 等）
3. `transforms.py:23` 蛋白原子词汇表 `[1,6,7,8,16,34]` → `[1,6,7,8,16,26,34]`（加 Fe）
4. ESIBank 全集重训 + 评估

**复用度提示**：EXP002a 在 P450 干净数据上 Test AUC -0.005（feature_dim 28→31）。该负面结果恰好印证老师推论的方向："对 P450 而言 Fe 是常量、无信息量"。但 P450 子集内无对照，必须用 ESIBank 全集才能形成"含铁酶 vs 不含铁酶"的天然对比。

## 排查路径

### 1. 论文 + supplementary
- 论文：Cui et al. *Nature* **647**, 639–647 (2025) [doi:10.1038/s41586-025-09697-2]
- supplementary 中关于 ligand parsing / pocket extraction / atom vocabulary 的章节

### 2. 代码

| 文件 | 排查重点 |
|---|---|
| `src/Datasets/Structure/protein_ligand.py` | 口袋提取逻辑：是否过滤 HETATM？是否保留 HEM 残基？|
| `src/Datasets/Structure/transforms.py` | 原子词汇表：原始版本是 H/C/N/O/S/Se（Fe 是 EXP002a 加的）|
| `src/Datasets/brenda.py` | 数据集加载：是否在数据层过滤辅因子？|
| `src/Datasets/create_features.py` | ESM 特征 & 反应特征：是否处理 HEM 上下文？|

### 3. 实证测试

直接读论文 ckpt 的训练数据 LMDB 中某条 P450 样本：
- 检查口袋原子里是否含 HEM 残基的 Fe 原子
- 如有：未剔除；如无：剔除

## 待办

- [x] 逐文件审查：`protein_ligand.py` / `transforms.py` / `brenda.py` / `create_features.py` / `structure.py` / `example.ipynb` / `main_testing.py`
- [x] Claude 主审 + Codex 独立交叉验证（SESSION_ID: 019e2f0d-cbb6-7f40-9c11-6e3318d3066e）
- [x] 输出三方判定（Q-A / Q-B / Q-C）
- [ ] （可选）实证测试：拿 P450 样本读 LMDB 验证蛋白原子里确实无 Fe
- [ ] （可选）读论文 supplementary 中关于 ligand parsing 章节，看作者是否承认这个 gap

## 输出物

- [x] 本 session_log 已是完整 audit 报告（含代码引用 + 证据片段 + 实现层 bug 分析 + 对 Q1 影响）

## 与其他问题的关联

- ✅ **解锁 Q1**：Q1 实验设计依赖的本问题结论已敲定 → 改造点最小化
- ✅ **重新解读 EXP002a**：EXP002a 在 P450 子集上 -0.005 正是老师推论方向的预期（"P450 内 Fe 是常量"），但要支持论文主张必须做 ESIBank 全集对照

## 变更日志

| 日期 | 变更 |
|---|---|
| 2026-05-08 | session 创建 |
| 2026-05-16 | Claude + Codex 双方审计完成，结论锁定，本 session_log 升级为完整 audit 报告 |

