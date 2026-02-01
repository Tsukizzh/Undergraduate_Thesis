# Step 7.1: Morgan指纹 - 文件索引

> **最后更新**: 2026-01-30

## 脚本文件

| 文件 | 路径 | 描述 |
|------|------|------|
| step7_1_generate_morgan.py | scripts/07_Step7_分子指纹生成/ | Morgan指纹生成脚本 |

## 数据文件

### 输入
| 文件 | 路径 | 描述 |
|------|------|------|
| Substrates.csv | data/04_Step4_格式修正后数据/ | 436个底物及其SMILES |

### 输出
| 文件 | 路径 | 大小 | 描述 |
|------|------|------|------|
| morgan_fingerprint.npy | data/07_Step7_分子指纹/ | 436.1 KB | Morgan指纹 (436, 1024) int8格式 |

## 会话文档
| 文件 | 路径 | 描述 |
|------|------|------|
| session_log.md | sessions/07_Step7_分子指纹生成/01_Morgan指纹/ | 本次会话日志 |
| file_index.md | sessions/07_Step7_分子指纹生成/01_Morgan指纹/ | 本文件索引 |
