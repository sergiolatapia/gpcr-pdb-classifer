"""
Step 2 -- which chain of the entry is the GPCR.

Every chain of >= 80 CA is CE-aligned to a 7TM reference panel; among chains with
a valid 7TM alignment (>= 130 residues, RMSD <= 5.0 A) the longest alignment wins,
then the lowest RMSD, then the most resolved residues.

The algorithm is copied unchanged from pdb_classification_test/
2_gpcr_chain_classification.py. Benchmarked on 1,712 GPCRdb entries it never
picked a non-GPCR chain: every disagreement with GPCRdb's preferred chain was an
identical copy of the receptor, a chain-ID naming difference (AAA vs A), or the
other GPCR subunit of a class C heterodimer. Only the reference file lookup
differs: files are pinned below instead of searched for in several folders.
"""

import os

from pymol import cmd

# Comprehensive GPCR Reference Panel across classes
REFS_CONFIG = {
    "ClassA_ADRB2": ("2rh1", "chain A"),
    "ClassA_RHO":   ("1f88", "chain A"),
    "ClassA_OPRM1": ("4dkl", "chain A"),
    "ClassA_S1PR1": ("5zkp", "chain A"),
    "ClassB1_GCGR": ("4l6r", "chain A"),
    "ClassB1_GLP1R": ("6x18", "chain R"),
    "ClassB2_AGRD1":("7wu2", "chain R"),
    "ClassC_GRM1":  ("4or2", "chain A"),
    "ClassD1_STE2": ("7ad3", "chain A"),
    "ClassF_SMO":   ("4jkv", "chain A"),
    "ClassO1_OR51E2":("8f76", "chain A"),
    "ClassO2_consOR1":("8uxy", "chain A"),
    "ClassT2_TAS2R46":("7xp5", "chain R")
}

# Single-chain reference files are loaded whole; the full 1F88 entry is cut to chain A,
# exactly as the benchmark script did.
REF_FILES = {
    "2rh1": "ClassA_2rh1_A.pdb", "1f88": "1f88.cif", "4dkl": "ClassA_4dkl_A.pdb",
    "5zkp": "ClassA_5zkp_A.pdb", "4l6r": "ClassB1_GCGR_4l6r_A.pdb", "6x18": "ClassB1_GLP_6X18_R.pdb",
    "7wu2": "ClassB2_7wu2_R.pdb", "4or2": "ClassC_mGlu_4OR2_A.pdb", "7ad3": "ClassD1_7ad3_A.pdb",
    "4jkv": "ClassF_4jkv_A.pdb", "8f76": "ClassO1_8f76_A.pdb", "8uxy": "ClassO2_8uxy_A.pdb",
    "7xp5": "ClassT2_TAS2R46_7xp5_R.pdb",
}


def load_refs(ref_dir):
    for ref_name, (pdb_code, sel) in REFS_CONFIG.items():
        obj_name = f"ref_{ref_name}"
        fname = REF_FILES[pdb_code]
        cmd.load(os.path.join(ref_dir, fname), f"raw_{obj_name}")
        if fname.endswith(".pdb"):      # pre-cut single-chain reference
            cmd.create(obj_name, f"raw_{obj_name} and polymer.protein and name CA+C+N+O")
        else:
            cmd.create(obj_name, f"raw_{obj_name} and ({sel}) and polymer.protein and name CA+C+N+O")
        cmd.delete(f"raw_{obj_name}")


def process_single_cif(cif_path):
    """
    Evaluates every chain in target CIF using CEAlign against the reference panel.
    Returns: (pdb_id, best_chain, best_rmsd, best_aln_len, best_ref)
    """
    filename = os.path.basename(cif_path)
    pdb_id = filename[:4].upper()

    try:
        cmd.load(cif_path, "target_struct")
    except Exception as e:
        return (pdb_id, None, 999.0, 0, f"Load Error: {e}")

    chains = cmd.get_chains("target_struct")
    if not chains:
        cmd.delete("target_struct")
        return (pdb_id, None, 999.0, 0, "No chains found")

    candidate_chains = []

    for ch in chains:
        ca_count = cmd.count_atoms(f"target_struct and chain {ch} and name CA")
        if ca_count < 80:  # Skip small peptides, tags, nanobody fragments
            continue

        best_ch_rmsd = float('inf')
        best_ch_aln_len = 0
        best_ch_ref = ""

        for ref_name in REFS_CONFIG.keys():
            ref_obj = f"ref_{ref_name}"
            try:
                aln = cmd.cealign(ref_obj, f"target_struct and chain {ch}")
                rmsd = aln['RMSD']
                aln_len = aln['alignment_length']

                # Check if this reference gives better alignment for this chain
                # Priority: longer alignment length with reasonable RMSD (< 5.0 A)
                if (aln_len > best_ch_aln_len and rmsd <= min(best_ch_rmsd + 0.5, 5.0)) or (aln_len == best_ch_aln_len and rmsd < best_ch_rmsd):
                    best_ch_rmsd = rmsd
                    best_ch_aln_len = aln_len
                    best_ch_ref = ref_name
                elif best_ch_aln_len < 130 and rmsd < best_ch_rmsd:
                    best_ch_rmsd = rmsd
                    best_ch_aln_len = aln_len
                    best_ch_ref = ref_name
            except Exception:
                pass

        candidate_chains.append({
            "chain": ch,
            "ca_count": ca_count,
            "rmsd": best_ch_rmsd,
            "aln_len": best_ch_aln_len,
            "ref": best_ch_ref
        })

    cmd.delete("target_struct")

    if not candidate_chains:
        return (pdb_id, None, 999.0, 0, "No chains passed residue filter")

    # Filter 1: Valid 7TM alignments (aln_len >= 130 and rmsd <= 5.0 A)
    valid_gpcr_chains = [c for c in candidate_chains if c["aln_len"] >= 130 and c["rmsd"] <= 5.0]

    if valid_gpcr_chains:
        # If multiple chains qualify (e.g. homodimers), rank by:
        # 1. Alignment length (higher is better)
        # 2. RMSD (lower is better)
        # 3. CA count (more resolved residues is better)
        valid_gpcr_chains.sort(key=lambda x: (-x["aln_len"], x["rmsd"], -x["ca_count"], x["chain"]))
        best = valid_gpcr_chains[0]
        return (pdb_id, best["chain"], best["rmsd"], best["aln_len"], best["ref"])
    else:
        # Fallback: Pick chain with highest aln_len or lowest RMSD
        candidate_chains.sort(key=lambda x: (-x["aln_len"], x["rmsd"], -x["ca_count"]))
        best = candidate_chains[0]
        return (pdb_id, best["chain"], best["rmsd"], best["aln_len"], f"Fallback: {best['ref']}")
