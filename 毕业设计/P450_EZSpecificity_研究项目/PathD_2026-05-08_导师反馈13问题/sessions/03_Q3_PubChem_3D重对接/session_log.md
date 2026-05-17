# Q3 — PubChem 3D 配体重对接

## 老师原话

> 优化分子对接中的配体 3D 结构来源（换成 PubChem 中的 3D sdf），并且进行重新对接。

## 状态

🔴 **待启动**（已有部分审计：`pubchem_3d_audit_2026-05-01.md`）

## 背景

### 当前配体 3D 来源
- RDKit 生成的 3D 构象（ETKDG + UFF/MMFF 优化）
- 部分通过 Meeko 修复（处理 GROVER 失败 + RDKit 解析失败）
- **问题**：RDKit 生成的构象只是经验力场最低能量构象，不一定接近实验态

### PubChem 3D 优势
- PubChem 提供"3D conformer"——基于 OMEGA / 实验晶体的优化构象
- 多构象选项（多 pose 重对接可提高发现真实结合模式概率）

## 已有材料

- ✅ `PathC/.../C3_P450专属模型训练/pubchem_3d_audit_2026-05-01.md`：PubChem 3D 覆盖率审计
- ✅ `pubchem_3d_coverage.csv`：每个底物的 PubChem CID + 是否有 3D
- ✅ `pubchem_3d_raw.json`：API 原始返回数据

## 待澄清点

1. **PubChem 3D 覆盖率多少？**（需要先读 audit 报告）
2. **缺失的怎么 fallback**？继续用 RDKit？还是用别的工具（OMEGA / UCSF Chimera ETKDG-Boltzmann）？
3. **重对接范围**：只重做正样本（4751 对）还是全部 50180 对（含负样本）？
4. **重对接对模型实验的影响**：触发 Q1 / Q2 / Q4 / 后续所有结构通道实验的全链路重训

## 工作量估计

| 阶段 | 工作量 |
|---|---|
| 读取已有 audit + 补全缺失底物的 PubChem 3D | 0.5–1 天 |
| 重对接 50,180 对（Uni-Dock fast）| 12–24 小时 GPU |
| 重建 complex PDB | 4–6 小时 |
| 重建 LMDB + pt_cache | 6–8 小时 |
| 重训对比实验 | 1–2 天 |

## 待办

- [ ] 读取 `pubchem_3d_audit_2026-05-01.md` 确认覆盖率
- [ ] 决定 fallback 策略
- [ ] 跑 PubChem 3D 批量下载补全
- [ ] 重对接（先选小批量验证流程）
- [ ] 重建 LMDB + pt_cache
- [ ] 重训 EXP001_allfix_unified（基线对比）→ 评估

## 与其他问题的关联

- ⚠️ **优先级敏感**：Q3 应早于 Q1 / Q4 / Q13 完成，否则后者基于"旧对接"的实验结果会作废
- **Q13** Fe-催化原子距离筛选 必须用准确的对接构象（取决于 Q3 是否完成）

## 变更日志

| 日期 | 变更 |
|---|---|
| 2026-05-08 | session 创建 |

