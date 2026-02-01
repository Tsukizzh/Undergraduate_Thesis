# Session记录 - Step 2: 提取酶序列

## 基本信息
- 时间：2026-01-09
- 步骤：Step 2
- 任务：提取酶序列（生成Enzymes.csv）

---

## 一、这一步要做什么？为什么要做？

### 1.1 目标
EZSpecificity模型需要三个输入文件：
1. `Enzymes.csv` - **酶的氨基酸序列** ← 本步骤
2. `Substrates.csv` - 底物的SMILES表示 (Step 3)
3. `data.csv` - 酶-底物对的索引关系 (Step 4)

### 1.2 什么是氨基酸序列？
蛋白质（酶）由氨基酸链组成。氨基酸序列就是用20个字母表示蛋白质的"一维结构"。

例如：`MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV...`

每个字母代表一种氨基酸：
- M = 甲硫氨酸 (Methionine)
- K = 赖氨酸 (Lysine)
- T = 苏氨酸 (Threonine)
- ...共20种标准氨基酸

### 1.3 为什么需要序列？
EZSpecificity模型使用ESM-2（一个蛋白质语言模型）将序列转换为向量表示，用于预测酶-底物特异性。

---

## 二、两种获取序列的方式

### 2.1 方式A：从PDB文件提取（初始方案）

**PDB文件是什么？**
PDB (Protein Data Bank) 文件包含蛋白质的3D原子坐标，是X射线晶体学或冷冻电镜实验的结果。

**提取逻辑**：
```
PDB文件 → 读取ATOM记录 → 提取残基名称 → 转换为单字母代码 → 序列
```

**问题**：PDB文件中的序列是"observed sequence"（实际解析出结构的部分），可能缺少：
- N端（蛋白质开头）的柔性区域
- C端（蛋白质结尾）的柔性区域
- 中间未解析的loop区域

### 2.2 方式B：使用源数据中的序列（最终方案）

**源数据是什么？**
`source_data/01_核心数据/ML训练数据集_186个.csv` 是之前从RCSB API获取的数据，包含一列 `sequence`。

**这个序列是什么？**
是"deposited sequence"（作者提交的完整序列），比PDB文件中的observed sequence更完整。

---

## 三、初始方案的问题发现

### 3.1 我最初做了什么？
1. 编写了 `step2_extract_sequence.py` 脚本
2. 使用BioPython库解析158个PDB/CIF文件
3. 对每个文件，选择最长的蛋白质链
4. 将三字母氨基酸代码转换为单字母代码
5. 生成了Enzymes.csv（158条序列）

### 3.2 用户发现了什么问题？
用户问："你确定Enzymes.csv里的序列和PDB ID对应上吗？"

我进行了验证，比较了：
- 源数据中的序列长度（来自RCSB API）
- 从PDB文件提取的序列长度

### 3.3 验证结果
```
PDB      源数据长度    提取长度    差异
1F4T        368         367       -1
1N97        389         385       -4
1NR6        473         462      -11
1ODO        408         396      -12
...
```

**结论**：150/158个序列长度不一致！提取的序列普遍更短。

### 3.4 具体例子：1F4T
```
源数据序列末尾：...RLVVRLKSNE  (368个氨基酸)
提取序列末尾：  ...RLVVRLKSN   (367个氨基酸)
```
差异：PDB文件中缺少C端最后一个残基 "E"（谷氨酸）。

**原因**：晶体结构中，C端的最后一个残基可能因为柔性太大，无法解析出电子密度。

---

## 四、方案修正

### 4.1 决策
放弃从PDB文件提取，直接使用源数据中的 `sequence` 列。

### 4.2 为什么这样做？
1. 源数据的序列是RCSB API返回的"deposited sequence"，是完整的
2. 保证序列和PDB ID的对应关系（源数据中有明确的pdb_id列）
3. 避免PDB文件解析可能带来的其他问题

### 4.3 实施步骤
1. 读取源数据 `ML训练数据集_186个.csv`
2. 筛选出158个有效PDB（排除与ESIBank重叠的26个）
3. 按PDB文件名排序
4. 提取 `sequence` 列，写入 `Enzymes.csv`
5. 生成审计文件 `step2_pdb_sequence_mapping.csv`，记录PDB ID和序列的对应关系

---

## 五、数据清洗（之前Step 1遗留的问题）

### 5.1 发现非P450数据
在审查数据时，Codex发现有20个PDB的title不含"P450"或"cytochrome"关键词。

### 5.2 验证结果
- 18个是P450（只是命名特殊，如"fatty acid decarboxylase"）
- 2个不是P450：
  - **6ZZX**：Photosystem I（光合系统I，完全不是P450！）
  - 8D8P：fatty acid decarboxylase（经查证，这是P450家族成员）

### 5.3 6ZZX为什么会混进来？
回溯到之前的RCSB搜索任务，发现是RCSB数据库的注释错误或搜索API的误匹配。

### 5.4 处理
- 从源数据中删除6ZZX（186条→184条）
- 从PDB文件目录中删除6ZZX.pdb（159个→158个）

---

## 六、最终产出

### 6.1 生成的文件

**Enzymes.csv**（158行）
- 只包含一列：`Protein sequence`
- 每行是一个P450酶的完整氨基酸序列
- 按PDB文件名字母顺序排列

**step2_pdb_sequence_mapping.csv**（158行）
- 审计文件，记录每行序列对应哪个PDB
- 列：row_index, file_name, pdb_id, sequence_length
- 用于追溯和验证

### 6.2 序列统计
- 数量：158条
- 长度范围：368-786 氨基酸
- 平均长度：437 氨基酸

### 6.3 数据流总结
```
源数据 ML训练数据集_186个.csv
    │
    ├─ 移除6ZZX（非P450）→ 184条
    │
    ├─ 排除26个PDB（与ESIBank重叠）→ 158条
    │
    └─ 提取sequence列 → Enzymes.csv (158条)
                │
                └─ 审计文件记录PDB ID对应关系
```

---

## 七、与Codex/Gemini的讨论要点

1. **序列来源**：最终决定用源数据而非PDB文件提取
2. **非P450识别**：Codex帮助验证了20个可疑记录
3. **数据一致性**：用户发现对应问题，促使方案修正

---

## 八、经验教训

1. **不要假设数据正确**：即使是从权威来源（RCSB）获取的数据，也可能有错误（如6ZZX）
2. **验证对应关系**：生成数据后，要验证索引/ID的对应关系是否正确
3. **保留审计文件**：方便后续追溯和验证

---

## 产生的文件
见 file_index.md
