"""
Post-processing cleanup for PCPD substrate extraction results.
Fixes issues found by Codex audit.
"""
import json, re, os

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "downloads", "PCPD")
INPUT = os.path.join(BASE, "substrate_extraction_final.jsonl")
OUTPUT = os.path.join(BASE, "substrate_extraction_cleaned.jsonl")

# ── Blacklist: things that are NOT substrates ──
BLACKLIST = {
    "An Organic Molecule", "NADPH", "NADH", "NAD+", "NADP+",
    "O2", "H2O2", "H2O", "oxygen", "molecular oxygen",
    "Bifunctional cytochrome P450/NADPH--P450",
    "Bifunctional cytochrome P450",
    "cytochrome P450", "P450",
}

BLACKLIST_PATTERNS = [
    r'^P450[\-\s]',              # P450-terp, P450cam etc (enzyme names)
    r'^CYP\d',                    # CYP IDs
    r'^[A-Z][0-9][A-Z0-9]{3}[0-9]$',  # UniProt IDs
    r'^[A-Z0-9]{10}$',           # Long accessions
    r'^pdb[\s_]',                 # PDB references
    r'^\d+[A-Z]{3}\d*$',         # PDB IDs like 8EUL
    r'hydroxylation',             # Reaction descriptions, not substrates
    r'epoxidation',
    r'oxidation',
    r'reduction',
    r'demethylation',
    r'dealkylation',
    r'regio-',
    r'stereoselective',
    r'^both\s',
    r'^many\s',
    r'^several\s',
    r'^various\s',
    r'^a\s+nitrogen',             # "a nitrogen atom"
    r'^the\s+(?!vitamin)',        # "the natural product X" but keep "the" in some cases
    r'reactions?\s+involved',
    r'participates',
    r'^\d+alpha-hydroxylation',
    r'^\d+beta-hydroxylation',
]

def is_blacklisted(name):
    name_stripped = name.strip()
    if name_stripped in BLACKLIST:
        return True
    for p in BLACKLIST_PATTERNS:
        if re.search(p, name_stripped, re.IGNORECASE):
            return True
    return False

def clean_substrate(name):
    """Clean a substrate name."""
    s = name.strip()
    # Remove leading determiners
    s = re.sub(r'^(?:the|a|an)\s+(?:natural\s+product\s+)?', '', s, flags=re.IGNORECASE)
    # Remove trailing position descriptors
    s = re.sub(r'\s+C\d+\s+methyl$', '', s)
    s = re.sub(r'\s+\d+alpha$|\s+\d+beta$', '', s)
    # Remove "hydroxylation" suffix if attached
    s = re.sub(r'\s+hydroxylation$', '', s, flags=re.IGNORECASE)
    # Remove duplicate names (e.g., "Mevastatin|Mevastatin")
    if '|' in s:
        parts = [p.strip() for p in s.split('|')]
        parts = list(dict.fromkeys(parts))
        s = parts[0] if len(set(p.lower() for p in parts)) == 1 else '|'.join(parts)
    # Remove semicolon-contaminated suffixes
    if ';' in s:
        parts = s.split(';')
        # Keep only parts that look like chemical names
        clean_parts = []
        for p in parts:
            p = p.strip()
            if p and not is_blacklisted(p) and len(p) > 2:
                clean_parts.append(p)
        s = clean_parts[0] if clean_parts else s.split(';')[0].strip()
    # Fix tokenization issues like "Prostagland in" -> "Prostaglandin"
    s = re.sub(r'Prostagland\s+in\b', 'Prostaglandin', s)
    s = re.sub(r'Dehydroepiand\s+rosterone', 'Dehydroepiandrosterone', s)
    s = s.strip().rstrip('.').strip()
    return s

def expand_coordinated(substrates):
    """Expand coordinated forms like 'both alpha- and beta-ionone'."""
    expanded = []
    for s in substrates:
        # "both alpha- and beta-ionone"
        m = re.match(r'(?:both\s+)?(alpha|beta|gamma|cis|trans|ent|epi)[-\s]+and\s+(alpha|beta|gamma|cis|trans|ent|epi)[-\s]*(\S+)', s, re.IGNORECASE)
        if m:
            expanded.append(f"{m.group(1)}-{m.group(3)}")
            expanded.append(f"{m.group(2)}-{m.group(3)}")
        else:
            expanded.append(s)
    return expanded

# ── Process ──
entries = []
with open(INPUT, encoding="utf-8") as f:
    for line in f:
        entries.append(json.loads(line))

fixes = {"blacklisted": 0, "cleaned": 0, "expanded": 0, "status_fixed": 0}

for e in entries:
    subs = e.get("substrates_specific", [])

    # Step 1: Clean each substrate
    cleaned = []
    for s in subs:
        s_clean = clean_substrate(s)
        if is_blacklisted(s_clean):
            fixes["blacklisted"] += 1
            continue
        if s_clean != s:
            fixes["cleaned"] += 1
        cleaned.append(s_clean)

    # Step 2: Expand coordinated forms
    expanded = expand_coordinated(cleaned)
    if len(expanded) != len(cleaned):
        fixes["expanded"] += 1

    # Step 3: Filter empty/short
    final = [s for s in expanded if s and len(s) > 2]

    # Step 4: Deduplicate (case-insensitive)
    seen = set()
    deduped = []
    for s in final:
        key = s.lower()
        if key not in seen:
            seen.add(key)
            deduped.append(s)

    e["substrates_specific"] = deduped

    # Step 5: Fix status
    if e["status"] == "specific" and not deduped:
        if e.get("substrates_class"):
            e["status"] = "class_only"
        else:
            e["status"] = "no_substrate"
        fixes["status_fixed"] += 1
    elif e["status"] == "class_only" and not e.get("substrates_class") and deduped:
        e["status"] = "specific"
        fixes["status_fixed"] += 1

# ── Stats ──
final_stats = {}
for e in entries:
    s = e["status"]
    final_stats[s] = final_stats.get(s, 0) + 1

all_subs = set()
for e in entries:
    all_subs.update(e.get("substrates_specific", []))

print(f"Fixes applied:")
print(f"  Blacklisted removed: {fixes['blacklisted']}")
print(f"  Names cleaned: {fixes['cleaned']}")
print(f"  Coordinated expanded: {fixes['expanded']}")
print(f"  Status fixed: {fixes['status_fixed']}")

print(f"\n=== FINAL CLEANED STATS ===")
for k, v in sorted(final_stats.items(), key=lambda x: -x[1]):
    print(f"  {k}: {v}")
print(f"  Total: {sum(final_stats.values())}")
print(f"\nUnique specific substrate names: {len(all_subs)}")

# Save
with open(OUTPUT, "w", encoding="utf-8") as f:
    for e in entries:
        f.write(json.dumps(e, ensure_ascii=False) + "\n")
print(f"Saved to {OUTPUT}")
