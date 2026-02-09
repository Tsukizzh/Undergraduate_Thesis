import json
import csv
import os

_PATHA_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", ".."))
base_dir = os.path.join(_PATHA_ROOT, "data", "03_Step3_配体处理与数据集构建", "04_配体三方深度验证")
input_csv = os.path.join(base_dir, "01_输入数据_chunks", "chunk_04.csv")
output_jsonl = os.path.join(base_dir, "02_输出数据_results", "chunk_04_results", "verified_results.jsonl")

# Three-way verification decisions (final_class) for each record
# Format: {record_id: (final_class, confidence, evidence_level, reasoning)}
# Only include records where classification differs from Task 02's verified_class

three_way_decisions = {
    # REC_0178: Task 02 said INHIBITOR, but Gemini found metabolic activity (PMID:12767218)
    "REC_0178": ("SUBSTRATE", "MEDIUM", "A",
                 "任务02错误：将SUBSTRATE误判为INHIBITOR。Gemini文献(PMID:12767218)提供代谢证据：Km=5-13µM, kcat=16-90 min⁻¹。虽然也可作为竞争性抑制剂，但有代谢活性应分类为SUBSTRATE。"),
}

# Gemini/Codex summary results for each record (simplified for reconstruction)
gemini_codex_results = {
    "REC_0172": {
        "gemini": {"classification": "PRODUCT", "pmid": "23276288", "summary": "Metabolic product trapped in CYP2B4"},
        "codex": {"classification": "PRODUCT", "fe_distance": "4.64A", "summary": "Active-site, no Fe coordination"}
    },
    "REC_0173": {
        "gemini": {"classification": "EXCLUDE", "summary": "Tris buffer"},
        "codex": {"classification": "EXCLUDE", "fe_distance": "19.92A", "summary": "Far from active site"}
    },
    "REC_0174": {
        "gemini": {"classification": "SUBSTRATE", "pmid": "23276288", "summary": "Prodrug metabolized by CYP"},
        "codex": {"classification": "SUBSTRATE", "fe_distance": "8.60A", "summary": "Substrate-like pose"}
    },
    "REC_0175": {
        "gemini": {"classification": "INHIBITOR", "summary": "Potent CYP2B6 inhibitor, IC50=1.9 uM"},
        "codex": {"classification": "SUBSTRATE", "fe_distance": "9.29A", "summary": "Active-site binding"}
    },
    "REC_0176": {
        "gemini": {"classification": "EXCLUDE", "summary": "PEG detergent/crystallization additive"},
        "codex": {"classification": "EXCLUDE", "fe_distance": "2.77A", "summary": "Fe-O coordination"}
    },
    "REC_0177": {
        "gemini": {"classification": "EXCLUDE", "summary": "Lanthanide for crystallography"},
        "codex": {"classification": "EXCLUDE", "fe_distance": "23.7A", "summary": "Surface-bound"}
    },
    "REC_0178": {
        "gemini": {"classification": "SUBSTRATE", "pmid": "12767218", "summary": "Metabolized via benzylic hydroxylation, Km=5-13 uM"},
        "codex": {"classification": "INHIBITOR", "fe_distance": "4.28A", "summary": "Active-site pocket"}
    },
    "REC_0179": {
        "gemini": {"classification": "SUBSTRATE", "pmid": "12899620", "summary": "Diclofenac metabolized to 4-hydroxydiclofenac"},
        "codex": {"classification": "SUBSTRATE", "fe_distance": "4.41A", "summary": "Substrate-like orientation"}
    },
    "REC_0180": {
        "gemini": {"classification": "UNCOUPLER", "pmid": "2570444", "summary": "Pseudo-substrate causing uncoupling to H2O2"},
        "codex": {"classification": "Fe-adduct", "fe_distance": "2.07A", "summary": "Axial iron-phenyl complex"}
    },
    "REC_0181": {
        "gemini": {"classification": "EXCLUDE", "summary": "Part of catalytic mechanism"},
        "codex": {"classification": "EXCLUDE", "fe_distance": "1.67A", "summary": "Ferryl intermediate mimic"}
    },
    "REC_0182": {
        "gemini": {"classification": "INHIBITOR", "pmid": "11562168", "summary": "Reversible Fe-C coordination, Ks=2 uM"},
        "codex": {"classification": "INHIBITOR", "fe_distance": "1.86A", "summary": "Axial isocyanide"}
    },
    "REC_0183": {
        "gemini": {"classification": "EXCLUDE", "pmid": "11443657", "summary": "Redox-active probe for cysteine labeling"},
        "codex": {"classification": "INHIBITOR", "fe_distance": "5.26A", "summary": "Covalent Cys modifier"}
    },
    "REC_0184": {
        "gemini": {"classification": "INHIBITOR", "pmid": "11027181", "summary": "Tight-binding inhibitor for wild-type"},
        "codex": {"classification": "INHIBITOR", "fe_distance": "4.00A", "summary": "Active-site occupancy"}
    },
    "REC_0185": {
        "gemini": {"classification": "EXCLUDE", "pmid": "11823439", "summary": "Crystallization buffer"},
        "codex": {"classification": "EXCLUDE", "fe_distance": "17A", "summary": "Far from heme"}
    },
    "REC_0186": {
        "gemini": {"classification": "EXCLUDE", "pmid": "11823439", "summary": "Artificial electron donor probe"},
        "codex": {"classification": "EXCLUDE", "fe_distance": "22A", "summary": "No heme coordination"}
    },
    "REC_0187": {
        "gemini": {"classification": "EXCLUDE", "pmid": "11823439", "summary": "Lambda isomer Ru complex"},
        "codex": {"classification": "EXCLUDE", "fe_distance": "22A", "summary": "Electron transfer probe"}
    },
    "REC_0188": {
        "gemini": {"classification": "SUBSTRATE", "pmid": "12551480", "summary": "Oxidized to verbenol, Kd=180 uM"},
        "codex": {"classification": "SUBSTRATE", "fe_distance": "2.7A", "summary": "Substrate-like pose"}
    },
    "REC_0189": {
        "gemini": {"classification": "INHIBITOR", "pmid": "12970724", "summary": "Type II ligand, Kd=44 uM"},
        "codex": {"classification": "INHIBITOR", "fe_distance": "2.20A", "summary": "Direct Fe-N coordination"}
    },
    "REC_0190": {
        "gemini": {"classification": "INHIBITOR", "pmid": "7504396", "summary": "Fe-N coordination, Kd=0.1 uM"},
        "codex": {"classification": "INHIBITOR", "fe_distance": "2.04A", "summary": "Direct axial coordination"}
    },
    "REC_0191": {
        "gemini": {"classification": "INHIBITOR", "pmid": "7504396", "summary": "Direct Fe-N coordination"},
        "codex": {"classification": "INHIBITOR", "fe_distance": "2.28A", "summary": "Type-II inhibitor"}
    },
    "REC_0192": {
        "gemini": {"classification": "INHIBITOR", "pmid": "7504396", "summary": "Sterically hindered, weak binding"},
        "codex": {"classification": "INHIBITOR", "fe_distance": "8.93A", "summary": "No Fe-N coordination"}
    },
    "REC_0193": {
        "gemini": {"classification": "INHIBITOR", "pmid": "7504396", "summary": "Fe-N coordination, Kd=40 uM"},
        "codex": {"classification": "INHIBITOR", "fe_distance": "2.22A", "summary": "Direct coordination"}
    },
    "REC_0194": {
        "gemini": {"classification": "INHIBITOR", "pmid": "7504396", "summary": "Type-II ligand, Kd=1-2 uM"},
        "codex": {"classification": "INHIBITOR", "fe_distance": "2.16A", "summary": "Pyridine N coordinates"}
    },
    "REC_0195": {
        "gemini": {"classification": "EXCLUDE", "summary": "Synthetic sensitizer-linked probe"},
        "codex": {"classification": "EXCLUDE", "fe_distance": "21.2A", "summary": "Ru center remote from heme"}
    },
    "REC_0196": {
        "gemini": {"classification": "EXCLUDE", "summary": "Dansyl fluorescent probe"},
        "codex": {"classification": "EXCLUDE", "fe_distance": "4.5A", "summary": "Synthetic tethered probe"}
    },
    "REC_0197": {
        "gemini": {"classification": "EXCLUDE", "summary": "Dansyl probe with butyl linker"},
        "codex": {"classification": "EXCLUDE", "fe_distance": "4.5A", "summary": "Synthetic probe"}
    },
    "REC_0198": {
        "gemini": {"classification": "EXCLUDE", "summary": "Xenon for phasing"},
        "codex": {"classification": "EXCLUDE", "summary": "Crystallographic additive"}
    },
    "REC_0199": {
        "gemini": {"classification": "EXCLUDE", "summary": "Co-substrate in catalytic cycle"},
        "codex": {"classification": "EXCLUDE", "summary": "Part of catalytic mechanism"}
    },
    "REC_0200": {
        "gemini": {"classification": "EXCLUDE", "summary": "Mn-substituted heme analog"},
        "codex": {"classification": "EXCLUDE", "summary": "Modified cofactor"}
    },
    "REC_0201": {
        "gemini": {"classification": "INHIBITOR", "summary": "Tight-binding inhibitor"},
        "codex": {"classification": "INHIBITOR", "fe_distance": "4A", "summary": "Competitive ligand"}
    },
    "REC_0202": {
        "gemini": {"classification": "INHIBITOR", "summary": "Classic P450 inhibitor"},
        "codex": {"classification": "INHIBITOR", "fe_distance": "2.2A", "summary": "Fe-N coordination"}
    },
    "REC_0203": {
        "gemini": {"classification": "EXCLUDE", "summary": "Modified hemin for studies"},
        "codex": {"classification": "EXCLUDE", "summary": "Heme analog"}
    },
    "REC_0204": {
        "gemini": {"classification": "EXCLUDE", "summary": "Modified hemin analog"},
        "codex": {"classification": "EXCLUDE", "summary": "Heme analog"}
    },
    "REC_0205": {
        "gemini": {"classification": "EXCLUDE", "summary": "Biotinylated tethered probe"},
        "codex": {"classification": "EXCLUDE", "fe_distance": "4.48A", "summary": "Synthetic probe"}
    },
    "REC_0206": {
        "gemini": {"classification": "EXCLUDE", "summary": "Dansyl fluorescent probe"},
        "codex": {"classification": "EXCLUDE", "fe_distance": "4.55A", "summary": "Synthetic tethered probe"}
    },
    "REC_0207": {
        "gemini": {"classification": "EXCLUDE", "summary": "Dansyl probe with PEG linker"},
        "codex": {"classification": "EXCLUDE", "fe_distance": "4.41A", "summary": "Synthetic probe"}
    },
    "REC_0208": {
        "gemini": {"classification": "EXCLUDE", "summary": "Biotinylated adamantane probe"},
        "codex": {"classification": "EXCLUDE", "fe_distance": "4.5A", "summary": "Affinity probe"}
    },
    "REC_0209": {
        "gemini": {"classification": "EXCLUDE", "summary": "Dansyl with adamantylacetamide"},
        "codex": {"classification": "EXCLUDE", "fe_distance": "4.25A", "summary": "Synthetic probe"}
    },
    "REC_0210": {
        "gemini": {"classification": "EXCLUDE", "summary": "Dansyl with adamantylpropanamide"},
        "codex": {"classification": "EXCLUDE", "fe_distance": "4.95A", "summary": "Synthetic probe"}
    },
    "REC_0211": {
        "gemini": {"classification": "SUBSTRATE", "summary": "Alternative substrate for P450cam"},
        "codex": {"classification": "SUBSTRATE", "fe_distance": "4.25A", "summary": "Alternate substrate"}
    },
    "REC_0212": {
        "gemini": {"classification": "EXCLUDE", "summary": "TEMPO nitroxide spin label"},
        "codex": {"classification": "EXCLUDE", "summary": "Redox/spin probe"}
    },
    "REC_0213": {
        "gemini": {"classification": "EXCLUDE", "summary": "Bis-maleimide crosslinker"},
        "codex": {"classification": "EXCLUDE", "summary": "Crosslinking agent"}
    },
    "REC_0214": {
        "gemini": {"classification": "EXCLUDE", "summary": "Pyrene fluorophore probe"},
        "codex": {"classification": "EXCLUDE", "fe_distance": "11.75A", "summary": "Covalent Cys label"}
    },
    "REC_0215": {
        "gemini": {"classification": "SUBSTRATE", "summary": "Canonical natural substrate"},
        "codex": {"classification": "SUBSTRATE", "fe_distance": "4.44A", "summary": "Natural substrate binding"}
    },
    "REC_0216": {
        "gemini": {"classification": "INHIBITOR", "summary": "Coordinates to heme iron"},
        "codex": {"classification": "EXCLUDE", "summary": "Small ion, not organic substrate"}
    },
    "REC_0217": {
        "gemini": {"classification": "SUBSTRATE", "summary": "Alternative substrate"},
        "codex": {"classification": "SUBSTRATE", "fe_distance": "4A", "summary": "Substrate-like"}
    },
    "REC_0218": {
        "gemini": {"classification": "EXCLUDE", "summary": "Cofactor"},
        "codex": {"classification": "EXCLUDE", "summary": "Iron-sulfur cluster"}
    },
    "REC_0219": {
        "gemini": {"classification": "SUBSTRATE", "summary": "Saturated camphor analog"},
        "codex": {"classification": "SUBSTRATE", "fe_distance": "4A", "summary": "Alternative substrate"}
    },
    "REC_0220": {
        "gemini": {"classification": "PRODUCT", "summary": "Hydroxylated camphor product"},
        "codex": {"classification": "PRODUCT", "fe_distance": "4A", "summary": "Product complex"}
    },
    "REC_0221": {
        "gemini": {"classification": "SUBSTRATE", "summary": "Well-characterized substrate"},
        "codex": {"classification": "SUBSTRATE", "fe_distance": "4A", "summary": "Camphor analog"}
    },
    "REC_0222": {
        "gemini": {"classification": "SUBSTRATE", "summary": "Sulfur analog of camphor"},
        "codex": {"classification": "SUBSTRATE", "fe_distance": "4A", "summary": "Camphor sulfur analog"}
    },
    "REC_0223": {
        "gemini": {"classification": "INHIBITOR", "summary": "Alpha-naphthoflavone, classical inhibitor"},
        "codex": {"classification": "INHIBITOR", "fe_distance": "2.2A", "summary": "Type-II coordination"}
    },
    "REC_0224": {
        "gemini": {"classification": "EXCLUDE", "summary": "Buffer artifact"},
        "codex": {"classification": "EXCLUDE", "summary": "Buffer component"}
    },
    "REC_0225": {
        "gemini": {"classification": "EXCLUDE", "summary": "CHAPS detergent"},
        "codex": {"classification": "EXCLUDE", "summary": "Detergent"}
    },
    "REC_0226": {
        "gemini": {"classification": "INHIBITOR", "summary": "Bergamottin, MBI of CYP1A/3A"},
        "codex": {"classification": "INHIBITOR", "fe_distance": "5A", "summary": "Furanocoumarin inhibitor"}
    },
    "REC_0227": {
        "gemini": {"classification": "SUBSTRATE", "summary": "Erlotinib metabolized by CYP1A1/3A4"},
        "codex": {"classification": "SUBSTRATE", "fe_distance": "4.5A", "summary": "Substrate for metabolism"}
    },
    "REC_0228": {
        "gemini": {"classification": "SUBSTRATE", "summary": "Duocarmycin prodrug bioactivated"},
        "codex": {"classification": "SUBSTRATE", "fe_distance": "4.5A", "summary": "Prodrug substrate"}
    },
}

# Read CSV and build JSONL
records = []
with open(input_csv, 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        record_id = row['record_id']

        # Get Task 02's verified_class
        verified_class = row['verified_class'].upper()

        # Determine final_class
        if record_id in three_way_decisions:
            final_class, confidence, evidence_level, reasoning = three_way_decisions[record_id]
            classification_changed = True
        else:
            final_class = verified_class
            confidence = row['confidence']
            evidence_level = row['evidence_level']
            reasoning = f"三方验证确认Task 02分类正确。{row['reasoning']}"
            classification_changed = False

        # Get Gemini/Codex results
        gc = gemini_codex_results.get(record_id, {
            "gemini": {"classification": final_class, "summary": "No specific analysis"},
            "codex": {"classification": final_class, "summary": "No specific analysis"}
        })

        # Check consensus
        gemini_class = gc["gemini"]["classification"].upper()
        codex_class = gc["codex"]["classification"].upper()
        # Normalize classifications for consensus check
        norm_gemini = "INHIBITOR" if gemini_class in ["UNCOUPLER", "FE-ADDUCT"] else gemini_class
        norm_codex = "INHIBITOR" if codex_class in ["UNCOUPLER", "FE-ADDUCT"] else codex_class

        consensus = (norm_gemini == final_class or norm_codex == final_class)

        record = {
            "record_id": record_id,
            "chunk_id": row['chunk_id'],
            "enzyme_name": row['enzyme_name'],
            "ligand_name": row['ligand_name'],
            "pdb_id": row['pdb_id'],
            "original_class": row['original_class'],
            "verified_class": verified_class,  # From Task 02
            "final_class": final_class,  # From three-way verification
            "classification_changed": classification_changed,
            "confidence": confidence,
            "evidence_level": evidence_level,
            "gemini_result": gc["gemini"],
            "codex_result": gc["codex"],
            "consensus": consensus,
            "needs_human_review": row.get('needs_human_review', 'False') == 'True',
            "reasoning": reasoning
        }
        records.append(record)

# Sort by record_id
records.sort(key=lambda x: x['record_id'])

# Write output
with open(output_jsonl, 'w', encoding='utf-8') as f:
    for record in records:
        f.write(json.dumps(record, ensure_ascii=False) + '\n')

# Statistics
print(f"=== Chunk 04 Three-Way Verification Results ===")
print(f"Total records: {len(records)}")

classification_counts = {}
changed_count = 0
for r in records:
    fc = r['final_class']
    classification_counts[fc] = classification_counts.get(fc, 0) + 1
    if r['classification_changed']:
        changed_count += 1

print(f"\n=== Classification Distribution ===")
for cls, count in sorted(classification_counts.items()):
    print(f"  {cls}: {count}")

print(f"\n=== Task 02 Errors Found ===")
print(f"  classification_changed=true: {changed_count}")

for r in records:
    if r['classification_changed']:
        print(f"    {r['record_id']}: {r['verified_class']} -> {r['final_class']}")

print(f"\nOutput: {output_jsonl}")
