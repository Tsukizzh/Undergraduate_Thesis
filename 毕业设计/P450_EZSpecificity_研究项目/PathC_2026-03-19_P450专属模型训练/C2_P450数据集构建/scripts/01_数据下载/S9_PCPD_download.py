"""
Download all PCPD data: FASTA sequences, PDB structures, reaction images.
Supports resume (skips existing files).
"""
import json, os, urllib.request, time, sys

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "downloads", "PCPD")
JSON_FILE = os.path.join(BASE, "resource_new.json")
CDN = "https://d1en57qlwrlmqu.cloudfront.net/P450/20240612"

FASTA_DIR = os.path.join(BASE, "FASTA")
PDB_DIR = os.path.join(BASE, "PDB")
IMG_DIR = os.path.join(BASE, "img_reaction")

os.makedirs(FASTA_DIR, exist_ok=True)
os.makedirs(PDB_DIR, exist_ok=True)
os.makedirs(IMG_DIR, exist_ok=True)

with open(JSON_FILE, encoding="utf-8") as f:
    entries = json.load(f)

print(f"Total entries: {len(entries)}")

def download(url, dest):
    """Download a file. Returns True if successful."""
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "PCPD-Download/1.0"})
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = resp.read()
            with open(dest, "wb") as f:
                f.write(data)
            return True
    except Exception as e:
        return False

# Download all three types for each entry
fasta_ok, fasta_skip, fasta_fail = 0, 0, 0
pdb_ok, pdb_skip, pdb_fail = 0, 0, 0
img_ok, img_skip, img_fail = 0, 0, 0

for i, entry in enumerate(entries):
    cyp_id = entry["ID"]

    # FASTA
    fasta_path = os.path.join(FASTA_DIR, f"{cyp_id}.fasta")
    if os.path.exists(fasta_path):
        fasta_skip += 1
    else:
        if download(f"{CDN}/FASTA/{cyp_id}.fasta", fasta_path):
            fasta_ok += 1
        else:
            fasta_fail += 1
        time.sleep(0.05)

    # PDB
    pdb_path = os.path.join(PDB_DIR, f"{cyp_id}.pdb")
    if os.path.exists(pdb_path):
        pdb_skip += 1
    else:
        if download(f"{CDN}/PDB/{cyp_id}.pdb", pdb_path):
            pdb_ok += 1
        else:
            pdb_fail += 1
        time.sleep(0.05)

    # Reaction image
    img_path = os.path.join(IMG_DIR, f"{cyp_id}.png")
    if os.path.exists(img_path):
        img_skip += 1
    else:
        if download(f"{CDN}/img_reaction/{cyp_id}.png", img_path):
            img_ok += 1
        else:
            img_fail += 1
        time.sleep(0.05)

    if (i + 1) % 100 == 0:
        print(f"  [{i+1}/{len(entries)}] FASTA:{fasta_ok}ok/{fasta_fail}fail  PDB:{pdb_ok}ok/{pdb_fail}fail  IMG:{img_ok}ok/{img_fail}fail")

print(f"\n=== DONE ===")
print(f"FASTA: {fasta_ok} downloaded, {fasta_skip} skipped, {fasta_fail} failed")
print(f"PDB:   {pdb_ok} downloaded, {pdb_skip} skipped, {pdb_fail} failed")
print(f"IMG:   {img_ok} downloaded, {img_skip} skipped, {img_fail} failed")
