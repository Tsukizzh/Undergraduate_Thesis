"""
Extract substrate names from PCPD Function text.
Stage 1: Rule-based extraction (regex patterns).
Output: substrate_extraction.jsonl with status for each entry.
"""
import json, re, os

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "downloads", "PCPD")

# Load JSON data
with open(os.path.join(BASE, "resource_new.json"), encoding="utf-8") as f:
    json_data = {e["ID"]: e for e in json.load(f)}

# Load parsed FASTA data
with open(os.path.join(BASE, "fasta_parsed.json"), encoding="utf-8") as f:
    fasta_data = json.load(f)

print(f"JSON entries: {len(json_data)}")
print(f"FASTA entries: {len(fasta_data)}")

# ── Filters ──
UNIPROT_RE = re.compile(r'^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$|^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-Z0-9]{10}$')
CYP_RE = re.compile(r'^CYP\d+[A-Z0-9]*$', re.IGNORECASE)
GENE_RE = re.compile(r'^[a-z]{2,5}[A-Z0-9\-]*$')

def is_identifier(token):
    """Check if a token is a UniProt ID, CYP name, or gene name."""
    t = token.strip().rstrip(';').strip()
    if UNIPROT_RE.match(t):
        return True
    if CYP_RE.match(t):
        return True
    if GENE_RE.match(t) and len(t) < 12:
        return True
    if t.endswith('-like') or t.endswith(' like'):
        return True
    if re.match(r'^\d+_\d+_CYP', t):
        return True
    return False

def normalize_text(text):
    """Fix common issues in function text."""
    text = re.sub(r'andcis-', 'and cis-', text)
    text = re.sub(r'andbeta-', 'and beta-', text)
    text = re.sub(r'andalpha-', 'and alpha-', text)
    text = re.sub(r'and(\S)', r'and \1', text)
    return text

def extract_substrates(func_text):
    """Extract substrate names from function text. Returns (specific, classes, status)."""
    if not func_text or func_text.lower() == 'none':
        return [], [], 'none'

    text = normalize_text(func_text)

    # Filter out identifier-only entries
    parts = re.split(r'[;]', text)
    clean_parts = [p.strip() for p in parts if p.strip() and not is_identifier(p.strip())]
    if not clean_parts:
        return [], [], 'gene_only'

    text = '; '.join(clean_parts)
    specific = []
    classes = []

    # Pattern 1: "X monooxygenase/hydroxylase/oxidase/..." and "X side-chain cleavage enzyme"
    enzyme_suffix = r'(?:\d+[a-z]*-)?(?:mono(?:oxygenase)?|oxygenase|hydroxylase|oxidase|epoxidase|reductase|synthase|demethylase|dehydrogenase|desaturase|lyase|isomerase|peroxidase|dioxygenase)\b'
    m = re.match(rf'^(.+?)\s+{enzyme_suffix}', text, re.IGNORECASE)
    if m:
        substrate = m.group(1).strip()
        substrate = re.sub(r'^(P450\S*|cytochrome)\s*;?\s*', '', substrate, flags=re.IGNORECASE).strip()
        # Handle "14-alpha" prefix in "Sterol 14-alpha demethylase" - substrate is before the position
        substrate = re.sub(r'\s+\d+[-]?(?:alpha|beta)\s*$', '', substrate).strip()
        if substrate and len(substrate) > 2 and not is_identifier(substrate):
            specific.append(substrate)

    # Pattern 1b: "X side-chain cleavage enzyme"
    m = re.match(r'^(.+?)\s+side[- ]chain\s+cleavage', text, re.IGNORECASE)
    if m:
        substrate = m.group(1).strip()
        if substrate and not is_identifier(substrate):
            specific.append(substrate)

    # Pattern 2: "hydroxylation/oxidation/epoxidation/reduction of X"
    for m in re.finditer(r'(?:hydroxylation|oxidation|epoxidation|reduction|demethylation|dealkylation|cleavage|metabolism|conversion|transformation|O-dealkylation|N-oxidation|dehydrogenation)\s+of\s+(.+?)(?:\s+(?:to|by|at|in|from|with|via|using|resulting)\b|$|[,;])', text, re.IGNORECASE):
        substrate = m.group(1).strip()
        # Clean up
        substrate = re.sub(r'\s+(?:at the|at|the)\s+.*', '', substrate)
        if substrate and not is_identifier(substrate):
            specific.append(substrate)

    # Pattern 3: "catalyzes the X of Y" or "catalyzes X"
    for m in re.finditer(r'catalyz\w+\s+(?:the\s+)?(?:\w+\s+of\s+)?(.+?)(?:\s+(?:to|by|at|in|from|with|resulting)\b|$|[,;])', text, re.IGNORECASE):
        substrate = m.group(1).strip()
        substrate = re.sub(r'^(?:hydroxylation|oxidation|epoxidation|H2O2-dependent epoxidation)\s+of\s+', '', substrate, flags=re.IGNORECASE)
        if substrate and len(substrate) > 2 and not is_identifier(substrate):
            specific.append(substrate)

    # Pattern 4: "from X to Y" -> substrate is X
    for m in re.finditer(r'from\s+(.+?)\s+to\s+', text, re.IGNORECASE):
        substrate = m.group(1).strip()
        if substrate and not is_identifier(substrate):
            specific.append(substrate)

    # Pattern 5: Pipe-delimited substrates
    if '|' in text:
        for part in text.split('|'):
            part = part.strip()
            if part and not is_identifier(part):
                specific.append(part)

    # Pattern 6: "including A, B, and C" or "including A and B"
    for m in re.finditer(r'including\s+(.+?)(?:\.|$)', text, re.IGNORECASE):
        items = re.split(r',\s*|\s+and\s+', m.group(1))
        for item in items:
            item = item.strip()
            if item and len(item) > 2 and not is_identifier(item):
                specific.append(item)

    # Pattern 7: "activity on/toward/against X"
    for m in re.finditer(r'activity\s+(?:on|toward|against|with)\s+(.+?)(?:\s+(?:and|to|by|at|in)\b|$|[,;])', text, re.IGNORECASE):
        substrate = m.group(1).strip()
        if substrate and not is_identifier(substrate):
            specific.append(substrate)

    # Expand coordinated forms: "alpha- and beta-ionone" -> ["alpha-ionone", "beta-ionone"]
    expanded = []
    for s in specific:
        m = re.match(r'(alpha|beta|gamma|delta|cis|trans|ent|epi)-?\s+and\s+(alpha|beta|gamma|delta|cis|trans|ent|epi)-(\S+)', s)
        if m:
            expanded.append(f"{m.group(1)}-{m.group(3)}")
            expanded.append(f"{m.group(2)}-{m.group(3)}")
        else:
            expanded.append(s)
    specific = expanded

    # Classify: specific vs class
    CLASS_KEYWORDS = ['fatty acid', 'terpene', 'terpenoid', 'alkaloid', 'flavonoid', 'steroid',
                      'phenylpropanoid', 'coumarin', 'macrolide', 'herbicide', 'xenobiotic',
                      'drug', 'several', 'various', 'number of', 'long-chain', 'short-chain',
                      'medium-chain', 'polyunsaturated', 'derivatives', 'analogs', 'substrates']

    final_specific = []
    for s in specific:
        s = s.strip().rstrip('.').strip()
        if not s or len(s) < 3:
            continue
        if any(kw in s.lower() for kw in CLASS_KEYWORDS):
            classes.append(s)
        else:
            final_specific.append(s)

    # Deduplicate
    final_specific = list(dict.fromkeys(final_specific))
    classes = list(dict.fromkeys(classes))

    if final_specific:
        return final_specific, classes, 'specific'
    elif classes:
        return [], classes, 'class_only'
    else:
        return [], [], 'unclear'

# ── Process all entries ──
results = []
stats = {'none': 0, 'gene_only': 0, 'specific': 0, 'class_only': 0, 'unclear': 0}

for cyp_id, jdata in json_data.items():
    func_json = jdata.get("Function", "")
    func_fasta = fasta_data.get(cyp_id, {}).get("func_from_fasta", "")
    species = fasta_data.get(cyp_id, {}).get("species", "")
    kingdom = jdata.get("Kindom", "")

    # Try JSON first
    specific, classes, status = extract_substrates(func_json)

    # If unclear/gene_only, try FASTA
    if status in ('unclear', 'gene_only', 'none') and func_fasta:
        spec2, cls2, stat2 = extract_substrates(func_fasta)
        if stat2 in ('specific', 'class_only'):
            specific = spec2
            classes = cls2
            status = stat2

    stats[status] += 1
    results.append({
        "cyp_id": cyp_id,
        "kingdom": kingdom,
        "species": species,
        "function_json": func_json,
        "function_fasta": func_fasta,
        "substrates_specific": specific,
        "substrates_class": classes,
        "status": status
    })

# Write results
output = os.path.join(BASE, "substrate_extraction.jsonl")
with open(output, "w", encoding="utf-8") as f:
    for r in results:
        f.write(json.dumps(r, ensure_ascii=False) + "\n")

print(f"\n=== EXTRACTION RESULTS ===")
for k, v in sorted(stats.items(), key=lambda x: -x[1]):
    print(f"  {k}: {v}")
print(f"  Total: {sum(stats.values())}")

# Count total substrate names
all_specific = set()
for r in results:
    all_specific.update(r["substrates_specific"])
print(f"\nUnique specific substrate names: {len(all_specific)}")

# Show samples
print("\nSamples:")
for r in results[:5]:
    if r["status"] == "specific":
        print(f"  {r['cyp_id']}: {r['substrates_specific']}")
        break
for r in results:
    if r["status"] == "class_only":
        print(f"  {r['cyp_id']}: class={r['substrates_class']}")
        break
for r in results:
    if r["status"] == "unclear":
        print(f"  {r['cyp_id']}: func={r['function_json'][:80]}")
        break
for r in results:
    if r["status"] == "gene_only":
        print(f"  {r['cyp_id']}: func={r['function_json'][:80]}")
        break
