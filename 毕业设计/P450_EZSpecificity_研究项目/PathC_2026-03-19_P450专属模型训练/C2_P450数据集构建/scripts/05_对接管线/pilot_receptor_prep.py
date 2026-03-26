"""
Phase 6 Stage C PILOT: Receptor PDB Cleanup + PDBQT
Tests 1 enzyme from each of 4 source types to validate the full pipeline.
Codex-reviewed: 2026-03-25
"""

import pandas as pd
import numpy as np
import subprocess
import shutil
import tempfile
import warnings
from pathlib import Path
from Bio.PDB import MMCIFParser, PDBParser, PDBIO, Select

warnings.filterwarnings('ignore')

BASE = Path(__file__).parent.parent.parent
DATA = BASE / "data"
STRUCTS = DATA / "structures"
MANIFESTS = STRUCTS / "manifests"
REGISTRIES = DATA / "registries"
OUT_CLEAN = STRUCTS / "receptors_clean"
OUT_PDBQT = STRUCTS / "receptors_pdbqt"
OUT_CLEAN.mkdir(parents=True, exist_ok=True)
OUT_PDBQT.mkdir(parents=True, exist_ok=True)

MEEKO = r"D:\anaconda3\envs\torch\Scripts\mk_prepare_receptor.exe"


def normalize_pdb(input_pdb, output_pdb):
    """Post-process cleaned PDB to fix common format issues for Meeko compatibility.

    Fixes:
    1. MSE (selenomethionine) → MET (methionine): Se → S
    2. HEM element column: fix X/blank/F → correct element (FE, C, N, O, S)
    3. Blank chain IDs → "A"

    Codex-reviewed: uses constrained element mapping, not generic atom_name[0].
    """
    def infer_hem_element(atom_name):
        """Infer element for HEM atom from atom name (constrained mapping)."""
        name = atom_name.strip().upper()
        if name == "FE":
            return "FE"
        if name.startswith("C"):
            return " C"
        if name.startswith("N"):
            return " N"
        if name.startswith("O"):
            return " O"
        if name.startswith("S"):
            return " S"
        if name.startswith("H"):
            return " H"
        return None  # unknown, leave as-is

    with open(input_pdb) as f:
        lines = f.readlines()

    fixed_lines = []
    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            fixed_lines.append(line)
            continue

        # Ensure line is at least 80 chars (pad if needed)
        line = line.rstrip("\n").ljust(80) + "\n"

        resname = line[17:20].strip()
        atom_name = line[12:16].strip()
        element = line[76:78].strip()
        chain_id = line[21:22]

        # Fix 1: MSE → MET
        if resname == "MSE":
            line = line[:17] + "MET" + line[20:]
            if atom_name == "SE":
                line = line[:12] + " SD " + line[16:]
                line = line[:76] + " S" + line[78:]
            if line.startswith("HETATM"):
                line = "ATOM  " + line[6:]

        # Fix 2: HEM element column
        if line[17:20].strip() == "HEM":
            inferred = infer_hem_element(atom_name)
            if inferred is not None:
                if element in ("X", "", "F") or (atom_name == "FE" and element != "FE"):
                    line = line[:76] + inferred.rjust(2) + line[78:]

        # Fix 3: Blank chain ID → "A"
        if chain_id == " ":
            line = line[:21] + "A" + line[22:]

        # Fix 5: Fix protein atom element column for standard amino acids
        # BioPython CIF→PDB conversion often writes atom name as element (CA→"CA", CE→"CE", etc.)
        # Standard amino acids only contain elements: C, N, O, S, H
        # Any 2-char "element" like CA, CB, CD, CE, CG, CZ, NE, NH, NZ, OG, SD etc. is wrong
        if resname != "HEM":
            name = atom_name.strip().upper()
            correct_element = None
            if name.startswith("C"):
                correct_element = " C"
            elif name.startswith("N"):
                correct_element = " N"
            elif name.startswith("O"):
                correct_element = " O"
            elif name.startswith("S"):
                correct_element = " S"
            elif name.startswith("H"):
                correct_element = " H"
            if correct_element and line[76:78] != correct_element:
                line = line[:76] + correct_element + line[78:]

        fixed_lines.append(line)

    # Fix 4: Reorder HEM atoms to canonical v2 order (Meeko requires specific order)
    # Canonical order from a working v2 PDBQT-compatible PDB
    CANONICAL_HEM_ORDER = [
        " O2D", " CGD", " O1D", " CBD", " CAD", " C3D", " C2D", " CMD",
        " C4D", " CHA", " ND ", " C1D", " CHD", "FE  ", " NB ", " C4B",
        " C3B", " CAB", " CBB", " C2B", " CMB", " C1B", " CHB", " NC ",
        " C4C", " C3C", " CAC", " CBC", " C2C", " CMC", " C1C", " CHC",
        " NA ", " C1A", " C4A", " C3A", " CMA", " C2A", " CAA", " CBA",
        " CGA", " O1A", " O2A",
    ]

    # Separate HEM and non-HEM lines
    hem_lines_map = {}  # atom_name_4char -> line
    non_hem_lines = []
    hem_insert_pos = None

    for i, line in enumerate(fixed_lines):
        if line.startswith(("ATOM", "HETATM")) and line[17:20].strip() == "HEM":
            atom_name_4 = line[12:16]
            hem_lines_map[atom_name_4] = line
            if hem_insert_pos is None:
                hem_insert_pos = len(non_hem_lines)
        else:
            non_hem_lines.append(line)

    if hem_lines_map and hem_insert_pos is not None:
        # Reorder HEM atoms according to canonical order
        reordered_hem = []
        for canon_name in CANONICAL_HEM_ORDER:
            if canon_name in hem_lines_map:
                reordered_hem.append(hem_lines_map[canon_name])
        # Add any remaining HEM atoms not in canonical order (e.g. hydrogens from PCPD)
        canonical_set = set(CANONICAL_HEM_ORDER)
        for name, line in hem_lines_map.items():
            if name not in canonical_set:
                reordered_hem.append(line)

        # Insert reordered HEM at the original position
        fixed_lines = non_hem_lines[:hem_insert_pos] + reordered_hem + non_hem_lines[hem_insert_pos:]

    with open(output_pdb, "w") as f:
        f.writelines(fixed_lines)


# Residues to remove
WATER = {"HOH", "WAT", "H2O", "DOD"}
IONS_BUFFERS = {"NA", "CL", "MG", "CA", "ZN", "K", "MN", "FE2", "FE3", "CO",
                "NI", "CU", "SO4", "PO4", "GOL", "EDO", "ACT", "DMS", "BME",
                "CIT", "TRS", "MPD", "PEG", "EPE", "MES", "IMD", "FMT"}


class ProteinHemSelect(Select):
    """BioPython selector: keep protein residues + one specific HEM residue."""
    def __init__(self, keep_chain_ids, hem_chain_id, hem_resseq):
        self.keep_chain_ids = keep_chain_ids
        self.hem_chain_id = hem_chain_id
        self.hem_resseq = hem_resseq

    def accept_chain(self, chain):
        return chain.id in self.keep_chain_ids

    def accept_residue(self, residue):
        resname = residue.get_resname().strip()
        chain_id = residue.get_parent().id
        hetflag = residue.id[0]

        # Keep the selected HEM
        if resname == "HEM" and chain_id == self.hem_chain_id and residue.id[1] == self.hem_resseq:
            return True
        # Keep standard amino acids (not HETATM)
        if hetflag == " " or hetflag == "H_MSE":  # MSE = selenomethionine
            if resname not in WATER and resname not in IONS_BUFFERS:
                return True
        # Remove everything else (water, ions, other ligands)
        return False


def find_hem_candidates(structure):
    """Find all HEM residues with Fe atoms in the structure."""
    candidates = []
    for model in structure:
        for chain in model:
            for res in chain:
                if res.get_resname().strip() == "HEM":
                    # Check both element field AND atom name for Fe
                    # Some PDB files have element="F" instead of "FE" (column alignment issue)
                    fe_atoms = [a for a in res if a.element.strip().upper() == "FE"
                                or a.get_name().strip().upper() == "FE"]
                    if fe_atoms:
                        fe = fe_atoms[0]
                        candidates.append({
                            "chain_id": chain.id,
                            "resseq": res.id[1],
                            "fe_coord": fe.get_coord(),
                            "n_atoms": len(list(res.get_atoms())),
                        })
    return candidates


def select_best_hem(candidates, manifest_fe=None):
    """Select the best HEM from candidates, preferring manifest Fe coordinates."""
    if not candidates:
        return None
    if len(candidates) == 1:
        return candidates[0]

    if manifest_fe is not None and not np.any(np.isnan(manifest_fe)):
        # Pick the HEM closest to manifest Fe coordinates
        dists = [np.linalg.norm(c["fe_coord"] - manifest_fe) for c in candidates]
        best_idx = np.argmin(dists)
        best = candidates[best_idx]
        best["manifest_dist"] = dists[best_idx]
        return best

    # Fallback: pick the one with most atoms (most complete HEM)
    return max(candidates, key=lambda c: c["n_atoms"])


def find_main_protein_chain(structure):
    """Find the chain with the most amino acid residues."""
    best_chain = None
    best_count = 0
    for chain in structure[0]:
        aa_count = sum(1 for r in chain if r.id[0] == " ")
        if aa_count > best_count:
            best_count = aa_count
            best_chain = chain.id
    return best_chain, best_count


def cleanup_structure(source_path, source_type, manifest_row, enzyme_index):
    """Clean a structure file: keep protein + 1 HEM, remove everything else."""
    result = {
        "enzyme_index": enzyme_index,
        "source_type": source_type,
        "source_path": str(source_path),
    }

    # Parse
    if source_type == "alphafill_cif":
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    try:
        structure = parser.get_structure("receptor", str(source_path))
    except Exception as e:
        result["status"] = f"parse_failed: {e}"
        return result

    # Find main protein chain
    main_chain, n_residues = find_main_protein_chain(structure)
    result["main_chain"] = main_chain
    result["n_protein_residues"] = n_residues

    # Find HEM candidates
    hem_candidates = find_hem_candidates(structure)
    result["n_hem_candidates"] = len(hem_candidates)

    if not hem_candidates:
        result["status"] = "no_hem_found"
        return result

    # Select best HEM
    manifest_fe = None
    if manifest_row is not None:
        try:
            fx, fy, fz = float(manifest_row.get("heme_fe_x", np.nan)), \
                         float(manifest_row.get("heme_fe_y", np.nan)), \
                         float(manifest_row.get("heme_fe_z", np.nan))
            manifest_fe = np.array([fx, fy, fz])
        except (ValueError, TypeError):
            pass

    best_hem = select_best_hem(hem_candidates, manifest_fe)
    result["selected_hem_chain"] = best_hem["chain_id"]
    result["selected_hem_resseq"] = best_hem["resseq"]
    result["fe_x"], result["fe_y"], result["fe_z"] = best_hem["fe_coord"]

    # For AlphaFill CIF: build new structure with single-char chain IDs
    # to avoid PDB format limitations (max 1-char chain ID)
    if source_type == "alphafill_cif":
        from Bio.PDB import StructureBuilder
        sb = StructureBuilder.StructureBuilder()
        sb.init_structure("clean")
        sb.init_model(0)

        # Copy protein chain → "A"
        sb.init_chain("A")
        src_chain = structure[0][main_chain]
        for res in src_chain:
            rn = res.get_resname().strip()
            if res.id[0] == " " and rn not in WATER and rn not in IONS_BUFFERS:
                sb.init_seg(" ")
                sb.init_residue(rn, res.id[0], res.id[1], res.id[2])
                for atom in res:
                    sb.init_atom(atom.get_name(), atom.get_coord(), atom.get_bfactor(),
                                 atom.get_occupancy(), atom.get_altloc(),
                                 atom.get_fullname(), atom.element)

        # Copy selected HEM → also in chain "A" (same chain as protein, like v2 format)
        hem_chain = structure[0][best_hem["chain_id"]]
        for res in hem_chain:
            if res.get_resname().strip() == "HEM" and res.id[1] == best_hem["resseq"]:
                sb.init_seg(" ")
                sb.init_residue("HEM", "H_HEM", best_hem["resseq"], " ")
                for atom in res:
                    sb.init_atom(atom.get_name(), atom.get_coord(), atom.get_bfactor(),
                                 atom.get_occupancy(), atom.get_altloc(),
                                 atom.get_fullname(), atom.element)

        new_structure = sb.get_structure()
        keep_chains = {"A"}
        best_hem["chain_id"] = "A"

        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False, prefix="clean_") as f:
            tmp_pdb = f.name
        io = PDBIO()
        io.set_structure(new_structure)
        io.save(tmp_pdb)
    else:
        # For PDB sources: use selector approach
        keep_chains = {main_chain, best_hem["chain_id"]}
        selector = ProteinHemSelect(keep_chains, best_hem["chain_id"], best_hem["resseq"])

        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False, prefix="clean_") as f:
            tmp_pdb = f.name
        io = PDBIO()
        io.set_structure(structure)
        io.save(tmp_pdb, selector)

    # Normalize PDB format (fix element columns, MSE→MET, blank chains)
    normalize_pdb(tmp_pdb, tmp_pdb)

    # Copy to final location
    out_pdb = OUT_CLEAN / f"{enzyme_index}.pdb"
    shutil.copy2(tmp_pdb, str(out_pdb))
    Path(tmp_pdb).unlink()

    # Validate
    parser2 = PDBParser(QUIET=True)
    s2 = parser2.get_structure("check", str(out_pdb))
    n_hem_final = sum(1 for r in s2.get_residues() if r.get_resname().strip() == "HEM")
    n_fe_final = sum(1 for a in s2.get_atoms() if a.element.strip().upper() == "FE"
                     or a.get_name().strip().upper() == "FE")
    n_protein_final = sum(1 for r in s2.get_residues() if r.id[0] == " ")

    result["clean_pdb"] = str(out_pdb)
    result["n_hem_final"] = n_hem_final
    result["n_fe_final"] = n_fe_final
    result["n_protein_final"] = n_protein_final

    if n_hem_final != 1 or n_fe_final != 1:
        result["status"] = f"validation_failed: HEM={n_hem_final}, Fe={n_fe_final}"
        return result

    result["status"] = "ok"
    return result


# HEM PDBQT template: atom_name → (charge, atom_type)
# Extracted from a working v2 Meeko output
HEM_PDBQT_TEMPLATE = {
    "CBC": ("0.009", "C"), "CAC": ("0.003", "C"), "C3C": ("0.003", "A"),
    "C2C": ("-0.018", "A"), "CMC": ("0.047", "C"), "C1C": ("0.049", "A"),
    "CHC": ("0.106", "C"), "C4B": ("0.213", "C"), "NB": ("-0.204", "N"),
    "C1B": ("0.209", "C"), "CHB": ("0.105", "C"), "C4A": ("0.048", "A"),
    "C3A": ("-0.022", "A"), "CMA": ("0.047", "C"), "C2A": ("-0.019", "A"),
    "CAA": ("0.044", "C"), "CBA": ("0.057", "C"), "CGA": ("0.042", "C"),
    "O1A": ("-0.550", "OA"), "O2A": ("-0.550", "OA"), "C1A": ("0.049", "A"),
    "CHA": ("0.105", "C"), "C4D": ("0.209", "C"), "C3D": ("0.034", "C"),
    "CAD": ("0.049", "C"), "CBD": ("0.058", "C"), "CGD": ("0.042", "C"),
    "O1D": ("-0.550", "OA"), "O2D": ("-0.550", "OA"), "C2D": ("0.030", "C"),
    "CMD": ("0.052", "C"), "C1D": ("0.209", "C"), "ND": ("-0.205", "N"),
    "FE": ("1.241", "Fe"), "NC": ("-0.354", "N"), "C4C": ("0.053", "A"),
    "CHD": ("0.106", "C"), "NA": ("-0.355", "N"), "C2B": ("0.034", "C"),
    "CMB": ("0.052", "C"), "C3B": ("0.055", "C"), "CAB": ("0.008", "C"),
    "CBB": ("0.009", "C"),
}


def test_meeko(clean_pdb_path, enzyme_index):
    """Prepare receptor PDBQT using split-and-merge strategy:
    1. Remove HEM from PDB → protein-only PDB
    2. Run Meeko on protein-only → protein.pdbqt
    3. Build HEM PDBQT block from clean PDB coordinates + v2 charge template
    4. Merge protein.pdbqt + HEM PDBQT → final receptor.pdbqt
    """
    result = {"enzyme_index": enzyme_index}

    with tempfile.TemporaryDirectory(prefix="meeko_rec_") as tmpdir:
        # Step 1: Create protein-only PDB (no HEM)
        with open(str(clean_pdb_path)) as f:
            all_lines = f.readlines()

        protein_lines = [l for l in all_lines if not (l.startswith(("ATOM", "HETATM")) and l[17:20].strip() == "HEM")]
        hem_lines = [l for l in all_lines if l.startswith(("ATOM", "HETATM")) and l[17:20].strip() == "HEM"]

        tmp_protein_pdb = Path(tmpdir) / "protein_only.pdb"
        with open(str(tmp_protein_pdb), "w") as f:
            f.writelines(protein_lines)

        # Step 2: Run Meeko on protein-only
        # -a: allow bad residues (skip incomplete amino acids, safe for protein-only)
        cmd = [MEEKO, "--read_pdb", str(tmp_protein_pdb), "-o", str(Path(tmpdir) / "protein"), "-p", "-a"]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=120)

        protein_pdbqt = Path(tmpdir) / "protein.pdbqt"
        if proc.returncode != 0 or not protein_pdbqt.exists():
            result["meeko_status"] = f"protein_meeko_failed"
            result["meeko_stderr"] = proc.stderr[:500]
            return result

        # Step 3: Build HEM PDBQT block from coordinates + template charges
        hem_pdbqt_lines = []
        # Get max serial from protein PDBQT
        with open(str(protein_pdbqt)) as f:
            prot_lines = f.readlines()
        max_serial = 0
        for l in prot_lines:
            if l.startswith(("ATOM", "HETATM")):
                try:
                    max_serial = max(max_serial, int(l[6:11].strip()))
                except ValueError:
                    pass

        serial = max_serial + 1
        for hem_line in hem_lines:
            atom_name = hem_line[12:16].strip()
            if atom_name in HEM_PDBQT_TEMPLATE:
                charge, atype = HEM_PDBQT_TEMPLATE[atom_name]
                # Build PDBQT line: same as PDB but with charge and type columns
                resname = "HEM"
                chain = hem_line[21:22]
                resseq = hem_line[22:26]
                x = hem_line[30:38]
                y = hem_line[38:46]
                z = hem_line[46:54]
                occ = hem_line[54:60] if len(hem_line) > 59 else " 1.00 "
                bfac = hem_line[60:66] if len(hem_line) > 65 else " 0.00 "

                pdbqt_line = f"ATOM  {serial:5d} {atom_name:4s} {resname:3s} {chain}{resseq}    {x}{y}{z}{occ}{bfac}    {float(charge):+6.3f} {atype:2s}\n"
                hem_pdbqt_lines.append(pdbqt_line)
                serial += 1

        # Step 4: Merge protein PDBQT + HEM PDBQT
        # Insert HEM before END
        merged_lines = []
        for l in prot_lines:
            if l.strip() == "END":
                merged_lines.extend(hem_pdbqt_lines)
            merged_lines.append(l)
        if not any(l.strip() == "END" for l in prot_lines):
            merged_lines.extend(hem_pdbqt_lines)

        # Write final PDBQT
        out_pdbqt = OUT_PDBQT / f"{enzyme_index}.pdbqt"
        with tempfile.NamedTemporaryFile(suffix=".pdbqt", delete=False, prefix="merged_") as f:
            tmp_merged = f.name
        with open(tmp_merged, "w") as f:
            f.writelines(merged_lines)
        shutil.copy2(tmp_merged, str(out_pdbqt))
        Path(tmp_merged).unlink()

        # Validate
        with open(str(out_pdbqt)) as f:
            final_content = f.read()
        fe_lines = [l for l in final_content.split("\n") if "FE" in l[12:16].upper() and "HEM" in l]
        hem_atom_count = sum(1 for l in final_content.split("\n") if l.startswith("ATOM") and "HEM" in l[17:20])

        result["meeko_status"] = "ok"
        result["pdbqt_size"] = out_pdbqt.stat().st_size
        result["fe_in_pdbqt"] = len(fe_lines) > 0
        result["hem_atoms_in_pdbqt"] = hem_atom_count
        result["fe_lines"] = fe_lines[:3]
        result["pdbqt_path"] = str(out_pdbqt)

    return result


def main():
    print("=" * 60)
    print("Phase 6 Stage C PILOT: Receptor Prep")
    print("=" * 60)

    manifest = pd.read_csv(MANIFESTS / "receptor_manifest.csv")

    # Select 1 pilot from each source
    pilots = []

    # 1. v2 transplant (simplest)
    v2 = manifest[manifest["status"] == "heme_transplanted"].iloc[0]
    pilots.append(("v2_transplant", v2, STRUCTS / "heme_transplant" / "pdb" /
                   Path(v2["heme_transplant_pdb_path"]).name))

    # 2. AlphaFill CIF
    af = manifest[manifest["status"] == "downloaded_with_heme"].iloc[0]
    pilots.append(("alphafill_cif", af, STRUCTS / "alphafill" / "cif" /
                   Path(str(af["alphafill_cif_path"])).name))

    # 3. PCPD PDB
    pcpd = manifest[manifest["status"] == "pcpd_predicted"].iloc[0]
    pcpd_name = Path(str(pcpd["existing_pdb_path"])).name
    pcpd_path = BASE / "downloads" / "PCPD" / "PDB" / pcpd_name
    pilots.append(("pcpd_pdb", pcpd, pcpd_path))

    # 4. Experimental PDB — search in multiple locations
    exp = manifest[manifest["status"] == "experimental_pdb"].iloc[0]
    exp_pdb_id = str(exp["existing_pdb_path"]).replace("s1_rcsb/", "")
    exp_path = None
    for candidate in [
        BASE / "downloads" / "RCSB" / "PDB" / f"{exp_pdb_id}.pdb",
        BASE.parent.parent / "PathA_2026-01-08_模型评估测试集构建" / "data" / "01_Step1_PDB文件" / f"{exp_pdb_id}.pdb",
    ]:
        if candidate.exists():
            exp_path = candidate
            break
    if exp_path is None:
        exp_path = BASE / "downloads" / "RCSB" / "PDB" / f"{exp_pdb_id}.pdb"  # will fail gracefully
    pilots.append(("experimental_pdb", exp, exp_path))

    for source_type, row, src_path in pilots:
        eidx = row["global_enzyme_id"]
        print(f"\n{'='*40}")
        print(f"Source: {source_type} | Enzyme: {eidx}")
        print(f"Path: {src_path}")
        print(f"Exists: {src_path.exists()}")

        if not src_path.exists():
            print(f"  SKIP: file not found")
            continue

        # Step 1: Cleanup
        print(f"\n  [Cleanup]")
        result = cleanup_structure(src_path, source_type, row, eidx)
        for k, v in result.items():
            if k != "source_path":
                print(f"    {k}: {v}")

        if result["status"] != "ok":
            print(f"  FAIL at cleanup: {result['status']}")
            continue

        # Step 2: Meeko
        print(f"\n  [Meeko]")
        meeko_result = test_meeko(result["clean_pdb"], eidx)
        for k, v in meeko_result.items():
            print(f"    {k}: {v}")


if __name__ == "__main__":
    main()
