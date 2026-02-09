import json

# Define records REC_0094 to REC_0103
records = [
    {
        "record_id": "REC_0094",
        "chunk_id": "02",
        "uniprot_id": "A6ZSR0",
        "enzyme_name": "CYP51 (S. cerevisiae)",
        "ligand_name": "(R)-2-(2,4-Difluorophenyl)-1,1-difluoro-3-(1H-tetrazol-1-yl)-1-(5-(4-(2,2,2-trifluoroethoxy)phenyl)pyridin-2-yl)propan-2-ol",
        "pdb_id": "5UL0",
        "original_class": "product",
        "verified_class": "INHIBITOR",
        "final_class": "INHIBITOR",
        "classification_changed": True,
        "confidence": "HIGH",
        "evidence_level": "A",
        "gemini_result": {
            "classification": "INHIBITOR",
            "data": "MIC=0.002-0.03ug/mL",
            "pmid": "28869766",
            "summary": "VT-1161 tetrazole antifungal drug candidate targeting CYP51 (Warrilow et al., Antimicrob. Agents Chemother., 2017)"
        },
        "codex_result": {
            "classification": "INHIBITOR",
            "fe_distance": "2.09A (Fe-N)",
            "binding_mode": "Tetrazole N coordinates Fe, active-site",
            "summary": "Type-II inhibitor"
        },
        "consensus": True,
        "needs_human_review": False,
        "reasoning": "任务02错误：将INHIBITOR误判为PRODUCT。Gemini (PMID:28869766) confirms VT-1161 is potent antifungal inhibitor. Codex (Fe-N=2.09A) shows type-II tetrazole coordination. Correcting to INHIBITOR."
    },
    {
        "record_id": "REC_0095",
        "chunk_id": "02",
        "uniprot_id": "A7U483",
        "enzyme_name": "CYP17A2 (zebrafish steroid 17a-hydroxylase)",
        "ligand_name": "Abiraterone",
        "pdb_id": "4R20",
        "original_class": "unknown",
        "verified_class": "INHIBITOR",
        "final_class": "INHIBITOR",
        "classification_changed": False,
        "confidence": "HIGH",
        "evidence_level": "A",
        "gemini_result": {
            "classification": "INHIBITOR",
            "data": "IC50=5-10nM",
            "pmid": "21935924",
            "summary": "Abiraterone inhibitor for CYP17A2 (Ohe et al., Mol. Cell. Endocrinol., 2011)"
        },
        "codex_result": {
            "classification": "INHIBITOR",
            "fe_distance": "1.92A (Fe-N)",
            "binding_mode": "Steroid nitrogen coordinates Fe, active-site",
            "summary": "Type-II steroidal inhibitor"
        },
        "consensus": True,
        "needs_human_review": False,
        "reasoning": "Gemini (IC50=5-10nM, PMID:21935924) and Codex (Fe-N=1.92A) both confirm INHIBITOR. Abiraterone is steroidal CYP17 inhibitor."
    },
    {
        "record_id": "REC_0096",
        "chunk_id": "02",
        "uniprot_id": "A7U483",
        "enzyme_name": "CYP17A2 (zebrafish steroid 17a-hydroxylase)",
        "ligand_name": "MERCURY (II) ION",
        "pdb_id": "4R20",
        "original_class": "unknown",
        "verified_class": "EXCLUDE",
        "final_class": "EXCLUDE",
        "classification_changed": False,
        "confidence": "HIGH",
        "evidence_level": "A",
        "gemini_result": {
            "classification": "INHIBITOR_GENERAL",
            "data": "General heavy metal toxicity",
            "pmid": "12529988",
            "summary": "Mercury inhibits P450s generally (general toxicology)"
        },
        "codex_result": {
            "classification": "EXCLUDE",
            "fe_distance": "19.36A (Fe-Hg)",
            "binding_mode": "Solvent, not in active site",
            "summary": "Crystallization artifact, far from heme"
        },
        "consensus": False,
        "needs_human_review": False,
        "reasoning": "Codex more accurate: Mercury is crystallization phasing artifact (Fe-Hg=19.36A, solvent position). Adopting EXCLUDE over Gemini general toxicology classification."
    },
    {
        "record_id": "REC_0097",
        "chunk_id": "02",
        "uniprot_id": "A7U483",
        "enzyme_name": "CYP17A2 (zebrafish steroid 17a-hydroxylase)",
        "ligand_name": "PROGESTERONE",
        "pdb_id": "4R21",
        "original_class": "substrate",
        "verified_class": "SUBSTRATE",
        "final_class": "SUBSTRATE",
        "classification_changed": False,
        "confidence": "HIGH",
        "evidence_level": "A",
        "gemini_result": {
            "classification": "SUBSTRATE",
            "data": "Km=0.2-1uM",
            "pmid": "15601824",
            "summary": "Progesterone substrate for CYP17A2 17-hydroxylase (Zhou et al., J. Endocrinol., 2004)"
        },
        "codex_result": {
            "classification": "SUBSTRATE",
            "fe_distance": "4.70A (no Fe-N)",
            "binding_mode": "Active-site",
            "summary": "Substrate positioning"
        },
        "consensus": True,
        "needs_human_review": False,
        "reasoning": "Gemini (Km=0.2-1uM, PMID:15601824) and Codex (Fe=4.70A) both confirm SUBSTRATE. Major physiological substrate."
    },
    {
        "record_id": "REC_0098",
        "chunk_id": "02",
        "uniprot_id": "A9ERX9",
        "enzyme_name": "P450_A9ERX9",
        "ligand_name": "MYRISTIC ACID",
        "pdb_id": "6GK6",
        "original_class": "unknown",
        "verified_class": "SUBSTRATE",
        "final_class": "SUBSTRATE",
        "classification_changed": False,
        "confidence": "HIGH",
        "evidence_level": "A",
        "gemini_result": {
            "classification": "SUBSTRATE",
            "data": "Kd=1-5uM",
            "pmid": "21357618",
            "summary": "Myristic acid substrate for CYP153A-like fatty acid hydroxylase (Scheps et al., Org. Biomol. Chem., 2011)"
        },
        "codex_result": {
            "classification": "SUBSTRATE",
            "fe_distance": "3.92A (no Fe-N)",
            "binding_mode": "Active-site",
            "summary": "Substrate positioning"
        },
        "consensus": True,
        "needs_human_review": False,
        "reasoning": "Gemini (Kd=1-5uM, PMID:21357618) and Codex (Fe=3.92A) both confirm SUBSTRATE."
    },
    {
        "record_id": "REC_0099",
        "chunk_id": "02",
        "uniprot_id": "A9FDB7",
        "enzyme_name": "P450_A9FDB7",
        "ligand_name": "PROGESTERONE",
        "pdb_id": "6F88",
        "original_class": "substrate",
        "verified_class": "SUBSTRATE",
        "final_class": "SUBSTRATE",
        "classification_changed": False,
        "confidence": "HIGH",
        "evidence_level": "C",
        "gemini_result": {
            "classification": "ERROR_MISIDENTIFIED",
            "data": "Gemini误认为是Progesterone Receptor",
            "pmid": "21131705",
            "summary": "Gemini error: confused with nuclear receptor"
        },
        "codex_result": {
            "classification": "SUBSTRATE",
            "fe_distance": "4.19A (no Fe-N)",
            "binding_mode": "Active-site",
            "summary": "Substrate positioning"
        },
        "consensus": False,
        "needs_human_review": False,
        "reasoning": "Gemini misidentified enzyme as nuclear receptor. Codex correct: P450_A9FDB7 is P450 enzyme, progesterone is SUBSTRATE (Fe=4.19A, active-site binding). Adopting Codex classification."
    },
    {
        "record_id": "REC_0100",
        "chunk_id": "02",
        "uniprot_id": "B2HH81",
        "enzyme_name": "P450_B2HH81",
        "ligand_name": "2-[BIS-(2-HYDROXY-ETHYL)-AMINO]-2-HYDROXYMETHYL-PROPANE-1,3-DIOL",
        "pdb_id": "7TEF",
        "original_class": "product",
        "verified_class": "EXCLUDE",
        "final_class": "EXCLUDE",
        "classification_changed": True,
        "confidence": "HIGH",
        "evidence_level": "A",
        "gemini_result": {
            "classification": "INHIBITOR/CRYSTAL_ADDITIVE",
            "data": "Weak interaction, buffer molecule",
            "pmid": "20400951",
            "summary": "Bis-Tris buffer often binds in P450 crystallization"
        },
        "codex_result": {
            "classification": "EXCLUDE",
            "fe_distance": "4.49A (N at 7.86A)",
            "binding_mode": "Active-site proximity but buffer molecule",
            "summary": "Buffer molecule, not substrate"
        },
        "consensus": True,
        "needs_human_review": False,
        "reasoning": "任务02错误：将buffer molecule Bis-Tris误判为PRODUCT。Gemini and Codex both confirm this is crystallization buffer, not product. Correcting to EXCLUDE."
    },
    {
        "record_id": "REC_0101",
        "chunk_id": "02",
        "uniprot_id": "B2HHT9",
        "enzyme_name": "CYP124A1 (M. marinum sterol/terpenoid)",
        "ligand_name": "(2E,6E)-3,7,11-trimethyldodeca-2,6,10-trien-1-yl acetate",
        "pdb_id": "8FJO",
        "original_class": "unknown",
        "verified_class": "SUBSTRATE",
        "final_class": "SUBSTRATE",
        "classification_changed": False,
        "confidence": "HIGH",
        "evidence_level": "A",
        "gemini_result": {
            "classification": "SUBSTRATE",
            "data": "Kd=2.4±0.3uM",
            "pmid": "25143376",
            "summary": "Farnesyl acetate substrate, terminal hydroxylation (Lazos et al., J. Biol. Chem., 2014)"
        },
        "codex_result": {
            "classification": "SUBSTRATE",
            "fe_distance": "5.10A (no Fe-N)",
            "binding_mode": "Active-site",
            "summary": "Substrate positioning"
        },
        "consensus": True,
        "needs_human_review": False,
        "reasoning": "Gemini (Kd=2.4uM, PMID:25143376) and Codex (Fe=5.10A) both confirm SUBSTRATE."
    },
    {
        "record_id": "REC_0102",
        "chunk_id": "02",
        "uniprot_id": "B2HHT9",
        "enzyme_name": "CYP124A1 (M. marinum sterol/terpenoid)",
        "ligand_name": "(2E,6E)-3,7,11-trimethyldodeca-2,6,10-trien-1-ol",
        "pdb_id": "8FKB",
        "original_class": "product",
        "verified_class": "SUBSTRATE",
        "final_class": "SUBSTRATE",
        "classification_changed": True,
        "confidence": "HIGH",
        "evidence_level": "A",
        "gemini_result": {
            "classification": "SUBSTRATE",
            "data": "Kd=1.5-3uM, terminal hydroxylation to ω-hydroxyfarnesol",
            "pmid": "25143376",
            "summary": "Farnesol is substrate for further hydroxylation (Lazos et al., J. Biol. Chem., 2014)"
        },
        "codex_result": {
            "classification": "SUBSTRATE",
            "fe_distance": "4.35A (no Fe-N)",
            "binding_mode": "Active-site",
            "summary": "Substrate positioning"
        },
        "consensus": True,
        "needs_human_review": False,
        "reasoning": "任务02错误：将SUBSTRATE误判为PRODUCT。Gemini (PMID:25143376) confirms farnesol is substrate for terminal hydroxylation to ω-hydroxyfarnesol. Codex (Fe=4.35A) supports substrate classification. Correcting to SUBSTRATE."
    },
    {
        "record_id": "REC_0103",
        "chunk_id": "02",
        "uniprot_id": "B2HHT9",
        "enzyme_name": "CYP124A1 (M. marinum sterol/terpenoid)",
        "ligand_name": "(3beta,8alpha,9beta)-3-hydroxycholest-5-en-7-one",
        "pdb_id": "8GDI",
        "original_class": "substrate",
        "verified_class": "SUBSTRATE",
        "final_class": "SUBSTRATE",
        "classification_changed": False,
        "confidence": "HIGH",
        "evidence_level": "A",
        "gemini_result": {
            "classification": "SUBSTRATE",
            "data": "Kd=1.1uM",
            "pmid": "25143376",
            "summary": "7-ketocholesterol substrate for CYP124A1 (Lazos et al., J. Biol. Chem., 2014)"
        },
        "codex_result": {
            "classification": "SUBSTRATE",
            "fe_distance": "4.22A (no Fe-N)",
            "binding_mode": "Active-site",
            "summary": "Substrate positioning"
        },
        "consensus": True,
        "needs_human_review": False,
        "reasoning": "Gemini (Kd=1.1uM, PMID:25143376) and Codex (Fe=4.22A) both confirm SUBSTRATE."
    }
]

# Append to JSONL file
output_file = r"c:/Users/Administrator/Desktop/EZSpecificity_Project/毕业设计/P450_EZSpecificity_研究项目/PathA_2026-01-08_模型评估测试集构建/data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/02_输出数据_results/chunk_02_results/verified_results.jsonl"

with open(output_file, 'a', encoding='utf-8') as f:
    for record in records:
        f.write(json.dumps(record, ensure_ascii=False) + '\n')

print(f"Successfully appended {len(records)} records (REC_0094 to REC_0103)")
