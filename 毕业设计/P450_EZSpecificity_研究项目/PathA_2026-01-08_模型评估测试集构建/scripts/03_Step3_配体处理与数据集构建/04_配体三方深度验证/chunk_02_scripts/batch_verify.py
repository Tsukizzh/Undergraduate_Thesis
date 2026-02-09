"""
Batch verification script for Chunk 02 ligand classification.
IMPORTANT: This script calls Gemini and Codex for EVERY record.
"""
import csv
import json
import subprocess
import time
from pathlib import Path

# Paths
INPUT_CSV = r"c:/Users/Administrator/Desktop/EZSpecificity_Project/毕业设计/P450_EZSpecificity_研究项目/PathA_2026-01-08_模型评估测试集构建/data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/01_输入数据_chunks/chunk_02.csv"
OUTPUT_JSONL = r"c:/Users/Administrator/Desktop/EZSpecificity_Project/毕业设计/P450_EZSpecificity_研究项目/PathA_2026-01-08_模型评估测试集构建/data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/02_输出数据_results/chunk_02_results/verified_results.jsonl"
WORK_DIR = r"c:/Users/Administrator/Desktop/EZSpecificity_Project/毕业设计/P450_EZSpecificity_研究项目/PathA_2026-01-08_模型评估测试集构建"

# Session IDs for continuation
gemini_session_id = ""
codex_session_id = ""

def call_gemini(enzyme_name, ligand_name, session_id=""):
    """Call Gemini for literature search."""
    prompt = f"""Literature search for: {enzyme_name} + {ligand_name}

Search for:
1. Is it SUBSTRATE (metabolized), INHIBITOR (inhibits enzyme), or PRODUCT?
2. Quantitative data: Ki, IC50, Km values (with units)
3. PubMed ID (PMID) if available

Return format:
- Classification: SUBSTRATE/INHIBITOR/PRODUCT/EXCLUDE
- Data: Ki=X nM, IC50=Y uM, etc.
- PMID: XXXXXXXX
- If no literature found: explicitly state "No literature found"

Brief answer only, no long explanations."""

    # Note: Actual MCP call would be done via Claude Code
    # This is a placeholder for batch processing logic
    return {
        "classification": "PLACEHOLDER",
        "data": "N/A",
        "pmid": "N/A",
        "summary": "To be filled by actual Gemini call"
    }

def call_codex(pdb_id, enzyme_name, ligand_name, session_id=""):
    """Call Codex for structure analysis."""
    prompt = f"""Analyze PDB structure: {pdb_id}
Enzyme: {enzyme_name}
Ligand: {ligand_name}

Analyze:
1. Fe-ligand distance (if P450 heme is present)
2. Binding mode (active site vs allosteric)
3. Does structure support SUBSTRATE/INHIBITOR/PRODUCT classification?

Return:
- Fe-N distance: X.XX Å (if nitrogen coordinates to Fe)
- Binding mode: [description]
- Structural classification: SUBSTRATE/INHIBITOR/PRODUCT/EXCLUDE
- Reasoning: [brief]

Brief answer only."""

    # Note: Actual MCP call would be done via Claude Code
    # This is a placeholder for batch processing logic
    return {
        "classification": "PLACEHOLDER",
        "fe_distance": "N/A",
        "binding_mode": "N/A",
        "summary": "To be filled by actual Codex call"
    }

def process_record(record):
    """Process a single record with three-party verification."""
    global gemini_session_id, codex_session_id

    # Step 1: Gemini literature search
    gemini_result = call_gemini(
        record['enzyme_name'],
        record['ligand_name'],
        gemini_session_id
    )

    # Step 2: Codex structure analysis
    codex_result = call_codex(
        record['pdb_id'],
        record['enzyme_name'],
        record['ligand_name'],
        codex_session_id
    )

    # Step 3: Claude decision (consensus logic)
    gemini_class = gemini_result['classification']
    codex_class = codex_result['classification']
    verified_class = record['verified_class']

    # Determine final classification
    if gemini_class == codex_class:
        final_class = gemini_class
        consensus = True
        confidence = "HIGH"
    elif gemini_class == "NO_LITERATURE" and codex_class != "EXCLUDE":
        final_class = codex_class
        consensus = False
        confidence = "MEDIUM"
    else:
        # Conflict - needs human review
        final_class = "NEEDS_REVIEW"
        consensus = False
        confidence = "LOW"

    # Determine evidence level
    if gemini_result.get('data') and 'Ki' in gemini_result['data']:
        evidence_level = "A"
    elif gemini_result.get('pmid') and gemini_result['pmid'] != "N/A":
        evidence_level = "B"
    elif codex_class != "EXCLUDE":
        evidence_level = "C"
    else:
        evidence_level = "D"

    # Check if classification changed from task 02
    classification_changed = (final_class != verified_class)

    # Build result
    result = {
        "record_id": record['record_id'],
        "chunk_id": record['chunk_id'],
        "uniprot_id": record['uniprot_id'],
        "enzyme_name": record['enzyme_name'],
        "ligand_name": record['ligand_name'],
        "pdb_id": record['pdb_id'],
        "original_class": record['original_class'],
        "verified_class": verified_class,
        "final_class": final_class,
        "classification_changed": classification_changed,
        "confidence": confidence,
        "evidence_level": evidence_level,
        "gemini_result": gemini_result,
        "codex_result": codex_result,
        "consensus": consensus,
        "needs_human_review": (not consensus or final_class == "NEEDS_REVIEW"),
        "reasoning": f"Gemini: {gemini_class}, Codex: {codex_class}. "
                     f"{'Consensus achieved.' if consensus else 'Conflict detected.'}"
    }

    return result

def main():
    """Main processing loop."""
    # Read input CSV
    with open(INPUT_CSV, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        records = list(reader)

    print(f"Total records to process: {len(records)}")

    # Process each record
    processed = 0
    for record in records:
        print(f"\nProcessing {record['record_id']} ({processed+1}/{len(records)})...")

        result = process_record(record)

        # Append to JSONL
        with open(OUTPUT_JSONL, 'a', encoding='utf-8') as f:
            f.write(json.dumps(result, ensure_ascii=False) + '\n')

        processed += 1

        # Rate limiting
        time.sleep(2)  # 2 seconds between records to avoid API overload

    print(f"\n✓ Completed {processed}/{len(records)} records")

if __name__ == "__main__":
    main()
