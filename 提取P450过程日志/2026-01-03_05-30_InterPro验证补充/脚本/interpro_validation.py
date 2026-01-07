#!/usr/bin/env python3
"""
InterPro验证补充脚本
验证P450识别的完整性，检查是否有遗漏

任务1: 随机抽样被排除酶，检查是否含有P450 InterPro域
任务2: 检查额外InterPro子类（IPR002397, IPR002403等）
任务3: 用Pfam PF00067交叉验证

2026-01-03
"""

import csv
import json
import os
import random
import sqlite3
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

# 路径配置
BASE_DIR = Path(r"C:\Users\Administrator\Desktop\EZSpecificity_Project")
LOG_DIR = BASE_DIR / "提取P450过程日志"
VERIFY_DIR = LOG_DIR / "2026-01-03_05-30_InterPro验证补充"
DATA_DIR = VERIFY_DIR / "数据文件"

# 输入文件
NON_P450_CSV = LOG_DIR / "2026-01-02_01-46_P450精确验证" / "非P450酶列表_24815个.csv"
P450_CSV = LOG_DIR / "2026-01-02_01-46_P450精确验证" / "P450酶列表_最终版389个.csv"
CACHE_DB = LOG_DIR / "2026-01-02_01-46_P450精确验证" / "UniProt查询缓存.db"

# P450 InterPro域（原始筛选使用的）
ORIGINAL_P450_INTERPRO = {"IPR001128", "IPR036396", "IPR002401"}

# 额外的P450 InterPro子类
ADDITIONAL_P450_INTERPRO = {
    "IPR002397",  # Cytochrome P450, E-class, group II
    "IPR002403",  # Cytochrome P450, E-class, group IV
    "IPR002394",  # Cytochrome P450, B-class (细菌)
    "IPR008067",  # Cytochrome P450, CYP4 family
    "IPR008068",  # Cytochrome P450, CYP6 family
}

# 所有P450 InterPro域（原始 + 额外）
ALL_P450_INTERPRO = ORIGINAL_P450_INTERPRO | ADDITIONAL_P450_INTERPRO

# Pfam P450域
P450_PFAM = "PF00067"


def load_non_p450_enzymes():
    """加载非P450酶列表"""
    enzymes = []
    with open(NON_P450_CSV, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            enzymes.append(row)
    return enzymes


def load_p450_enzymes():
    """加载P450酶列表"""
    enzymes = []
    with open(P450_CSV, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            enzymes.append(row)
    return enzymes


def check_interpro_domains(interpro_str, domains_to_check):
    """检查InterPro字符串是否包含指定的域"""
    if not interpro_str:
        return set()
    found = set()
    for domain in domains_to_check:
        if domain in interpro_str.upper():
            found.add(domain)
    return found


def query_uniprot_pfam(uniprot_ids, batch_size=100):
    """查询UniProt获取Pfam信息"""
    results = {}

    for i in range(0, len(uniprot_ids), batch_size):
        batch = uniprot_ids[i:i + batch_size]
        print(f"  查询Pfam信息: {i+1}-{min(i+batch_size, len(uniprot_ids))}/{len(uniprot_ids)}")

        # 构建查询
        id_query = " OR ".join([f"accession:{acc}" for acc in batch])
        query = f"({id_query})"

        params = {
            "query": query,
            "fields": "accession,xref_pfam",
            "format": "tsv",
            "size": str(len(batch)),
        }

        url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(params)

        try:
            req = urllib.request.Request(url, headers={
                "User-Agent": "InterProValidation/1.0",
                "Accept": "text/plain"
            })
            with urllib.request.urlopen(req, timeout=60) as resp:
                tsv = resp.read().decode("utf-8")

            lines = [ln for ln in tsv.strip().split("\n") if ln]
            if len(lines) > 1:
                reader = csv.DictReader(lines, delimiter="\t")
                for row in reader:
                    acc = row.get("Entry", "").strip().upper()
                    pfam = row.get("Pfam", "") or ""
                    if acc:
                        results[acc] = pfam

            time.sleep(0.3)  # 礼貌性延时

        except Exception as e:
            print(f"    警告: 查询失败 - {e}")
            for acc in batch:
                results[acc] = ""

    return results


def task1_sample_excluded_enzymes(non_p450_list, sample_size=100):
    """任务1: 随机抽样被排除酶，检查P450域"""
    print("\n" + "="*60)
    print("任务1: 随机抽样被排除酶，检查P450 InterPro域")
    print("="*60)

    # 随机抽样
    random.seed(42)  # 固定种子以便复现
    sample = random.sample(non_p450_list, min(sample_size, len(non_p450_list)))

    # 检查每个样本
    results = []
    with_any_p450_domain = 0
    with_original_domains = 0
    with_additional_domains = 0

    for enzyme in sample:
        uniprot_id = enzyme["uniprot_id"]
        interpro = enzyme.get("xref_interpro", "") or ""

        # 检查原始P450域
        original_found = check_interpro_domains(interpro, ORIGINAL_P450_INTERPRO)
        # 检查额外P450域
        additional_found = check_interpro_domains(interpro, ADDITIONAL_P450_INTERPRO)
        # 所有找到的
        all_found = original_found | additional_found

        has_any = len(all_found) > 0

        result = {
            "uniprot_id": uniprot_id,
            "protein_families": enzyme.get("protein_families", ""),
            "xref_interpro": interpro,
            "original_p450_domains": ";".join(sorted(original_found)) if original_found else "",
            "additional_p450_domains": ";".join(sorted(additional_found)) if additional_found else "",
            "has_any_p450_domain": has_any,
        }
        results.append(result)

        if has_any:
            with_any_p450_domain += 1
        if original_found:
            with_original_domains += 1
        if additional_found:
            with_additional_domains += 1

    # 输出结果
    print(f"\n抽样数量: {len(sample)}")
    print(f"含有任何P450 InterPro域: {with_any_p450_domain} ({100*with_any_p450_domain/len(sample):.1f}%)")
    print(f"  - 含原始筛选域 (IPR001128/036396/002401): {with_original_domains}")
    print(f"  - 含额外子类域 (IPR002397/002403等): {with_additional_domains}")

    # 保存详细结果
    output_path = DATA_DIR / "任务1_被排除酶抽样检查_100个.csv"
    with open(output_path, "w", newline="", encoding="utf-8-sig") as f:
        fieldnames = ["uniprot_id", "protein_families", "xref_interpro",
                      "original_p450_domains", "additional_p450_domains", "has_any_p450_domain"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print(f"\n详细结果已保存: {output_path}")

    return {
        "sample_size": len(sample),
        "with_any_p450_domain": with_any_p450_domain,
        "with_original_domains": with_original_domains,
        "with_additional_domains": with_additional_domains,
        "percentage": 100 * with_any_p450_domain / len(sample),
    }


def task2_check_additional_interpro(non_p450_list):
    """任务2: 检查额外InterPro子类是否有遗漏"""
    print("\n" + "="*60)
    print("任务2: 检查额外InterPro子类 (IPR002397, IPR002403等)")
    print("="*60)

    # 统计每个额外域的出现情况
    domain_counts = defaultdict(list)

    for enzyme in non_p450_list:
        uniprot_id = enzyme["uniprot_id"]
        interpro = enzyme.get("xref_interpro", "") or ""

        for domain in ADDITIONAL_P450_INTERPRO:
            if domain in interpro.upper():
                domain_counts[domain].append({
                    "uniprot_id": uniprot_id,
                    "protein_families": enzyme.get("protein_families", ""),
                    "xref_interpro": interpro,
                })

    # 输出结果
    print(f"\n在{len(non_p450_list)}个被排除酶中检查额外InterPro域:")
    print("-" * 40)

    total_with_additional = 0
    all_potential_missed = []

    for domain in sorted(ADDITIONAL_P450_INTERPRO):
        count = len(domain_counts[domain])
        print(f"  {domain}: {count}个")
        total_with_additional += count
        all_potential_missed.extend(domain_counts[domain])

    print("-" * 40)
    print(f"  总计: {total_with_additional}个")

    # 去重（一个酶可能有多个额外域）
    unique_enzymes = {e["uniprot_id"]: e for e in all_potential_missed}

    print(f"\n唯一酶数量（可能被遗漏）: {len(unique_enzymes)}")

    # 保存结果
    if unique_enzymes:
        output_path = DATA_DIR / "任务2_额外InterPro域检查结果.csv"
        with open(output_path, "w", newline="", encoding="utf-8-sig") as f:
            fieldnames = ["uniprot_id", "protein_families", "xref_interpro"]
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(unique_enzymes.values())
        print(f"详细结果已保存: {output_path}")

    return {
        "domain_counts": {k: len(v) for k, v in domain_counts.items()},
        "total_with_additional": total_with_additional,
        "unique_enzymes": len(unique_enzymes),
    }


def task3_pfam_crossvalidation(non_p450_list, p450_list, sample_size=200):
    """任务3: Pfam PF00067交叉验证"""
    print("\n" + "="*60)
    print("任务3: Pfam PF00067交叉验证")
    print("="*60)

    # 从被排除酶中抽样
    random.seed(43)
    sample = random.sample(non_p450_list, min(sample_size, len(non_p450_list)))
    sample_ids = [e["uniprot_id"] for e in sample]

    # 查询Pfam信息
    print(f"\n查询{len(sample_ids)}个被排除酶的Pfam信息...")
    pfam_results = query_uniprot_pfam(sample_ids)

    # 检查有多少含有PF00067
    with_pf00067 = []
    for enzyme in sample:
        uniprot_id = enzyme["uniprot_id"]
        pfam = pfam_results.get(uniprot_id, "")
        if P450_PFAM in pfam.upper():
            with_pf00067.append({
                "uniprot_id": uniprot_id,
                "protein_families": enzyme.get("protein_families", ""),
                "xref_interpro": enzyme.get("xref_interpro", ""),
                "xref_pfam": pfam,
            })

    print(f"\n抽样数量: {len(sample)}")
    print(f"含有PF00067 (P450 Pfam域): {len(with_pf00067)} ({100*len(with_pf00067)/len(sample):.1f}%)")

    # 同时验证已识别的P450酶是否都有PF00067
    print(f"\n验证已识别的389个P450酶的Pfam...")
    p450_ids = [e["uniprot_id"] for e in p450_list]
    p450_pfam_results = query_uniprot_pfam(p450_ids)

    p450_with_pf00067 = 0
    p450_without_pf00067 = []
    for enzyme in p450_list:
        uniprot_id = enzyme["uniprot_id"]
        pfam = p450_pfam_results.get(uniprot_id, "")
        if P450_PFAM in pfam.upper():
            p450_with_pf00067 += 1
        else:
            p450_without_pf00067.append({
                "uniprot_id": uniprot_id,
                "protein_families": enzyme.get("protein_families", ""),
                "xref_pfam": pfam,
            })

    print(f"已识别P450中含有PF00067: {p450_with_pf00067}/{len(p450_list)} ({100*p450_with_pf00067/len(p450_list):.1f}%)")
    print(f"已识别P450中缺少PF00067: {len(p450_without_pf00067)}个")

    # 保存结果
    if with_pf00067:
        output_path = DATA_DIR / "任务3_Pfam交叉验证_被排除酶中含PF00067.csv"
        with open(output_path, "w", newline="", encoding="utf-8-sig") as f:
            fieldnames = ["uniprot_id", "protein_families", "xref_interpro", "xref_pfam"]
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(with_pf00067)
        print(f"\n被排除酶中含PF00067详情已保存: {output_path}")

    return {
        "sample_size": len(sample),
        "excluded_with_pf00067": len(with_pf00067),
        "excluded_pf00067_percentage": 100 * len(with_pf00067) / len(sample),
        "p450_with_pf00067": p450_with_pf00067,
        "p450_without_pf00067": len(p450_without_pf00067),
        "p450_pf00067_percentage": 100 * p450_with_pf00067 / len(p450_list),
    }


def main():
    print("="*60)
    print("InterPro验证补充 - 检查P450识别的完整性")
    print("="*60)

    # 确保输出目录存在
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    # 加载数据
    print("\n加载数据...")
    non_p450_list = load_non_p450_enzymes()
    p450_list = load_p450_enzymes()
    print(f"  非P450酶: {len(non_p450_list)}个")
    print(f"  P450酶: {len(p450_list)}个")

    # 执行三个任务
    results = {}

    results["task1"] = task1_sample_excluded_enzymes(non_p450_list, sample_size=100)
    results["task2"] = task2_check_additional_interpro(non_p450_list)
    results["task3"] = task3_pfam_crossvalidation(non_p450_list, p450_list, sample_size=200)

    # 生成总结
    print("\n" + "="*60)
    print("验证总结")
    print("="*60)

    summary = {
        "验证时间": time.strftime("%Y-%m-%d %H:%M:%S"),
        "task1_抽样检查": {
            "说明": "随机抽样100个被排除酶，检查是否含有P450 InterPro域",
            "结果": f"{results['task1']['with_any_p450_domain']}个含有P450域 ({results['task1']['percentage']:.1f}%)",
            "结论": "假阴性风险极低" if results['task1']['with_any_p450_domain'] == 0 else f"发现{results['task1']['with_any_p450_domain']}个可能遗漏",
        },
        "task2_额外InterPro": {
            "说明": "检查被排除酶中是否含有额外P450 InterPro子类",
            "结果": results['task2']['domain_counts'],
            "唯一酶数量": results['task2']['unique_enzymes'],
            "结论": "无遗漏" if results['task2']['unique_enzymes'] == 0 else f"发现{results['task2']['unique_enzymes']}个可能遗漏",
        },
        "task3_Pfam交叉验证": {
            "说明": "用Pfam PF00067进行交叉验证",
            "被排除酶抽样": f"{results['task3']['sample_size']}个",
            "含PF00067": f"{results['task3']['excluded_with_pf00067']}个 ({results['task3']['excluded_pf00067_percentage']:.1f}%)",
            "已识别P450含PF00067": f"{results['task3']['p450_with_pf00067']}/{len(p450_list)} ({results['task3']['p450_pf00067_percentage']:.1f}%)",
        },
    }

    # 保存JSON结果
    json_path = DATA_DIR / "验证结果汇总.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"\n结果已保存: {json_path}")

    return results


if __name__ == "__main__":
    main()
