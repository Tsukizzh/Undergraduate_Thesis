# Step 1 v2 执行日志：下载PDB文件

> **执行时间**: 2026-01-21
> **版本**: v2 (基于实验二v2的740条记录)

---

## 一、执行结果摘要

| 项目 | 数值 |
|------|------|
| 需要的唯一PDB数 | 627 |
| 从v1归档复制 | 145 |
| 新下载 | 482 |
| 下载失败 | 0 |
| **最终PDB文件数** | **627** |

### 文件格式分布

| 格式 | 数量 | 说明 |
|------|------|------|
| PDB格式 (.pdb) | 619 | 标准格式 |
| CIF格式 (.cif) | 8 | 大型结构，RCSB只提供CIF |

---

## 二、详细执行过程

### 2.1 输入分析

**输入文件**: `source_data/01_核心数据/独立测试集_740条.csv`

```python
# 分析输入
df = pd.read_csv("独立测试集_740条.csv")
print(f"总记录数: {len(df)}")  # 740

unique_pdbs = df['pdb_id'].str.upper().unique()
print(f"唯一PDB数: {len(unique_pdbs)}")  # 627
```

**为什么740条记录只需要627个PDB？**
```
原因：一个PDB文件可以包含多个配体

示例：
PDB 1AKD 包含2个配体：
  - camphor (底物)
  - 某抑制剂
在740条记录中算作2条，但只需下载1个PDB文件

计算：740条记录 / 627个PDB ≈ 1.18 配体/PDB
```

### 2.2 检查已有资源

**v1归档位置**: `_archived_v1_按UniProt去重_Step1至3/data_v1/01_Step1_PDB原始文件/PDB结构文件_158个/`

```python
# 检查可复用的PDB
existing_files = list(OLD_PDB_DIR.glob("*.*"))
existing_pdb_ids = set(f.stem.upper() for f in existing_files)
print(f"v1归档中的PDB数: {len(existing_pdb_ids)}")  # 158

overlap = required_pdb_ids & existing_pdb_ids
missing = required_pdb_ids - existing_pdb_ids
print(f"可复用: {len(overlap)}")  # 145
print(f"需下载: {len(missing)}")  # 482
```

**复用分析**:
| 项目 | 数值 | 说明 |
|------|------|------|
| v1归档中的PDB | 158 | v1版本的全部PDB |
| v2需要且v1有 | 145 | 可复用 |
| v1有但v2不需要 | 13 | 被ESIBank排除规则过滤掉的 |
| v2需要但v1没有 | 482 | 需要新下载 |

### 2.3 复制已有文件

```python
# 复制可复用的PDB
copied = 0
for f in existing_files:
    pdb_id = f.stem.upper()
    if pdb_id in required_pdb_ids:
        shutil.copy2(f, NEW_PDB_DIR / f.name)
        copied += 1
print(f"复制完成: {copied}个")  # 145
```

### 2.4 下载缺失文件

```python
# 下载逻辑
for pdb_id in sorted(missing):
    pdb_lower = pdb_id.lower()

    # 优先尝试PDB格式
    url_pdb = f"https://files.rcsb.org/download/{pdb_lower}.pdb"
    response = requests.get(url_pdb, timeout=30)

    if response.status_code == 200:
        save(pdb_id + ".pdb")
    else:
        # 大结构只有CIF格式
        url_cif = f"https://files.rcsb.org/download/{pdb_lower}.cif"
        response = requests.get(url_cif, timeout=30)
        if response.status_code == 200:
            save(pdb_id + ".cif")
        else:
            failed.append(pdb_id)
```

**下载统计**:
| 阶段 | 数量 |
|------|------|
| 需要下载 | 482 |
| PDB格式成功 | 474 |
| CIF格式成功 | 8 |
| 失败 | 0 |

---

## 三、输出文件详情

### 3.1 输出位置

`data/01_Step1_PDB文件/` (627个文件)

### 3.2 文件命名规范

| 格式 | 命名 | 示例 |
|------|------|------|
| PDB | `{PDB_ID}.pdb` | `1AKD.pdb` |
| CIF | `{PDB_ID}.cif` | `7QAN.cif` |

### 3.3 CIF格式文件列表

以下8个结构较大，RCSB只提供CIF格式：

```
7QAN.cif
8VK6.cif
9CV8.cif
9KPU.cif
9MS2.cif
（以及新下载的3个）
```

---

## 四、数据来源说明

```
627个PDB文件的来源：

┌─────────────────────────────────────────────────────────────┐
│ 从v1归档复制: 145个                                          │
│   - 来源: _archived_v1_按UniProt去重_Step1至3/.../PDB结构文件_158个/        │
│   - 说明: v1的158个PDB中，有145个在v2的627个需求中           │
│   - 另外13个v1的PDB不在v2需求中（可能被ESIBank排除规则过滤）   │
├─────────────────────────────────────────────────────────────┤
│ 新下载: 482个                                                │
│   - 来源: RCSB PDB API                                       │
│   - URL格式: https://files.rcsb.org/download/{pdb_id}.pdb    │
│   - 如果PDB格式不可用，则下载CIF格式                          │
└─────────────────────────────────────────────────────────────┘
```

---

## 五、v1 vs v2 对比

| 项目 | v1 | v2 | 变化 |
|------|----|----| -----|
| 输入记录数 | 158 | 740 | +582 (+368%) |
| PDB文件数 | 158 | 627 | +469 (+297%) |
| 数据来源 | ML训练数据集_186个 | 独立测试集_740条 | 不同 |
| 去重规则 | 按UniProt | 按(UniProt + Ligand) | 不同 |

**为什么PDB数量增长不如记录数增长快？**
- 记录数增长368%（158→740）
- PDB数增长297%（158→627）
- 原因：v2中同一个PDB可能对应多条记录（多配体）

---

## 六、验证检查

### 6.1 文件完整性

| 检查项 | 预期 | 实际 | 状态 |
|--------|------|------|------|
| 总文件数 | 627 | 627 | ✅ |
| 复制文件数 | 145 | 145 | ✅ |
| 下载文件数 | 482 | 482 | ✅ |
| 下载失败数 | 0 | 0 | ✅ |

### 6.2 格式检查

| 检查项 | 预期 | 实际 | 状态 |
|--------|------|------|------|
| PDB格式数 | ~620 | 619 | ✅ |
| CIF格式数 | ~8 | 8 | ✅ |
| 文件名规范 | {PDB_ID}.{ext} | 符合 | ✅ |

---

## 七、脚本位置

`scripts/step1_download_pdb_v2.py`

---

## 八、下一步

Step 2: 从CSV的sequence字段生成Enzymes.csv

---

**执行者**: Claude Code + Codex + Gemini
**日期**: 2026-01-21
**文档版本**: v2.0
