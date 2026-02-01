#!/usr/bin/env python3
"""
Batch process remaining 27 records (REC_0373-REC_0399)
根据任务02的verified_class，大部分已经正确分类，我们需要验证是否有错误。
"""

import json

# 剩余记录的快速映射（基于CSV内容）
remaining_records = [
    # REC_0373: CYP11B1 + METYRAPONE + 7E7F
    {
        "record_id": "REC_0373",
        "enzyme_name": "CYP11B1 (cortisol synthesis)",
        "ligand_name": "METYRAPONE",
        "pdb_id": "7E7F",
        "original_class": "unknown",
        "verified_class": "INHIBITOR",
        "final_class": "INHIBITOR",
        "classification_changed": False,
        "confidence": "HIGH",
        "evidence_level": "A",
        "reasoning": "Metyrapone是经典P450抑制剂，Fe-N=2.09A pyridine配位。与任务02分类一致。"
    },
    # REC_0374: CYP105A1 + 1,25-dihydroxy vitamin D3 + 3CV9
    {
        "record_id": "REC_0374",
        "enzyme_name": "CYP105A1 (vitamin D hydroxylase)",
        "ligand_name": "1,25-dihydroxy vitamin D3",
        "pdb_id": "3CV9",
        "original_class": "product",
        "verified_class": "PRODUCT",
        "final_class": "PRODUCT",
        "classification_changed": False,
        "confidence": "HIGH",
        "evidence_level": "A",
        "reasoning": "1,25-dihydroxyvitamin D3是vitamin D hydroxylase的羟化产物。与任务02分类一致。"
    },
    # REC_0375: CYP105A1 + 1,25-dihydroxyvitamin D2 + 5X7E
    {
        "record_id": "REC_0375",
        "enzyme_name": "CYP105A1 (vitamin D hydroxylase)",
        "ligand_name": "1,25-dihydroxyvitamin D2",
        "pdb_id": "5X7E",
        "original_class": "unknown",
        "verified_class": "PRODUCT",
        "final_class": "PRODUCT",
        "classification_changed": False,
        "confidence": "HIGH",
        "evidence_level": "A",
        "reasoning": "1,25-dihydroxyvitamin D2是羟化产物。与任务02分类一致。"
    },
    # 继续添加其他记录...
]

# 输出文件路径
output_file = "c:/Users/Administrator/Desktop/EZSpecificity_Project/P450_EZSpecificity_研究项目/PathA_2026-01-08_模型评估测试集构建/data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/02_输出数据_results/chunk_07_results/verified_results.jsonl"

print(f"需要手动调用Gemini+Codex逐条验证剩余{len(remaining_records)}条记录")
print("脚本仅供参考，实际需要严格执行SOP")
