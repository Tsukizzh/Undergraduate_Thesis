# 并行验证操作完整指南

## 总览

本方案通过6个并行Claude Code窗口同时验证P450配体分类，大幅提升处理速度。

### 文件结构

```
scripts/03_Step3_配体处理与数据集构建/02_配体分类并行验证/
├── 00_一键准备.py                 # 一键准备：分块 + 生成SOP
├── 03_合并结果.py                 # 合并所有窗口的结果
├── README.md                      # 本说明文件
└── sop_prompt_chunk_01.md ~ 06.md # 自动生成的SOP Prompt

data/03_Step3_配体处理与数据集构建/02_配体分类并行验证/
├── 01_输入数据_chunks/            # 分块后的输入数据
│   ├── chunk_01.csv               # 窗口1的输入
│   ├── chunk_02.csv               # 窗口2的输入
│   └── ...
├── 02_输出数据_results/           # 各窗口的输出结果
│   ├── chunk_01_results/
│   │   └── verified_results.jsonl # 窗口1的验证结果
│   ├── chunk_02_results/
│   │   └── verified_results.jsonl # 窗口2的验证结果
│   └── ...

sessions/03_Step3_配体处理与数据集构建/02_配体分类并行验证/
└── session_log.md                 # 操作日志和文件索引

reports/03_Step3_配体处理与数据集构建/
├── 配体分类验证结果_并行处理_merged.jsonl  # 合并后的完整结果
├── 配体分类验证结果_并行处理_merged.csv    # CSV格式
└── 配体分类验证统计_并行处理.json          # 统计报告
```

---

## 完整流程

### 步骤1：一键准备

运行一键准备脚本，自动完成分块和SOP生成：

```bash
cd "C:\Users\Administrator\Desktop\EZSpecificity_Project\毕业设计\P450_EZSpecificity_研究项目\PathA_2026-01-08_模型评估测试集构建\scripts\03_Step3_配体处理与数据集构建\02_配体分类并行验证"

python 00_一键准备.py
```

**输出：**
- `data/.../01_输入数据_chunks/chunk_01.csv ~ 06.csv` - 6个分块CSV文件
- `data/.../02_输出数据_results/chunk_01_results ~ 06_results/` - 6个结果输出文件夹
- `scripts/.../sop_prompt_chunk_01.md ~ 06.md` - 6个SOP Prompt文件

---

### 步骤2：打开6个并行窗口

**窗口分配：**

| 窗口 | Chunk | 记录数（约） | 操作 |
|------|-------|--------------|------|
| 窗口1 | chunk_01 | ~100条 | 复制`sop_prompt_chunk_01.md` |
| 窗口2 | chunk_02 | ~100条 | 复制`sop_prompt_chunk_02.md` |
| 窗口3 | chunk_03 | ~100条 | 复制`sop_prompt_chunk_03.md` |
| 窗口4 | chunk_04 | ~100条 | 复制`sop_prompt_chunk_04.md` |
| 窗口5 | chunk_05 | ~100条 | 复制`sop_prompt_chunk_05.md` |
| 窗口6 | chunk_06 | ~100条 | 复制`sop_prompt_chunk_06.md` |

**操作步骤：**

1. 在VS Code中打开 `sop_prompt_chunk_01.md`
2. 复制全部内容 (Ctrl+A → Ctrl+C)
3. 打开新的Claude Code窗口
4. 粘贴内容作为第一条消息
5. 等待Claude开始处理
6. 重复步骤1-5，打开其他5个窗口

---

### 步骤3：监控进度

每个窗口会：
1. 读取对应的chunk CSV文件
2. 逐条验证记录
3. 实时追加写入 `verified_results.jsonl`
4. 每完成10条报告一次进度

**预计时间：**
- 每条记录约2-5分钟
- 每个窗口约114条 × 3分钟 = 6小时
- 并行处理总时间：约6-7小时（取决于最慢的窗口）

**监控技巧：**
- 可以随时打开 `verified_results.jsonl` 查看进度
- 文件每完成一条就追加一行，便于实时查看
- 如果某个窗口卡住，可以关闭后重新打开并继续（文件已追加的不会重复）

---

### 步骤4：合并结果

当所有6个窗口完成后：

```bash
python 03_合并结果.py
```

**输出：**
- `reports/03_Step3_配体处理与数据集构建/配体分类验证结果_并行处理_merged.jsonl`
- `reports/03_Step3_配体处理与数据集构建/配体分类验证结果_并行处理_merged.csv`
- `reports/03_Step3_配体处理与数据集构建/配体分类验证统计_并行处理.json`

**统计信息包括：**
- 总记录数
- 分类分布（SUBSTRATE/INHIBITOR/PRODUCT/EXCLUDE）
- 置信度分布（HIGH/MEDIUM/LOW）
- 证据等级分布（A/B/C/D）
- 陷阱检测统计
- 需要人工审核的记录数

---

## 测试验证结果

**在5条测试记录上验证了流程：**

| # | 酶-配体 | 验证结果 | 人工审核答案 | 匹配? | 关键点 |
|---|---------|----------|--------------|-------|--------|
| 1 | CYP2B4 + Paroxetine | INHIBITOR | INHIBITOR | ✅ | - |
| 2 | CYP101A1 + TCZ | EXCLUDE | EXCLUDE | ✅ | 检测到野生型vs突变体陷阱 |
| 3 | CYP2D6 + Ajmalicine | INHIBITOR | INHIBITOR | ✅ | - |
| 4 | P450TOL + TOLUENE | SUBSTRATE | SUBSTRATE | ✅ | - |
| 5 | TleB + B9R | SUBSTRATE | SUBSTRATE | ✅ | 检测到底物vs产物混淆陷阱 |

**成功率：5/5 (100%)** ✅

---

## 常见问题

### Q1: 窗口可以少于6个吗？
A: 可以。修改 `00_一键准备.py` 中的 `NUM_CHUNKS = 6` 改为你想要的数量。

### Q2: 如果某个窗口卡住或崩溃怎么办？
A: JSONL格式支持断点续传。重新打开窗口，粘贴相同的SOP，Claude会检测到已处理的记录并跳过。

### Q3: 如何检查某个chunk的进度？
A: 打开 `data/.../02_输出数据_results/chunk_XX_results/verified_results.jsonl` 文件，行数 = 已完成的记录数。

### Q4: 能否更改方法论？
A: 可以。修改 `00_一键准备.py` 中的 `generate_sop()` 函数。

### Q5: 如果两个窗口意外处理了同一条记录？
A: 合并脚本会检测重复的 `record_id` 并报警。需要手动清理重复项。

---

## 下一步

完成验证后：

1. **质量检查**：抽查10-20条"需人工审核"的记录
2. **与原始分类对比**：生成变更报告
3. **导师汇报准备**：使用统计数据和典型案例

---

**创建日期**：2026-01-26
**测试状态**：✅ 5条记录测试通过，100%准确率
**准备状态**：✅ 可以开始处理完整682条数据
