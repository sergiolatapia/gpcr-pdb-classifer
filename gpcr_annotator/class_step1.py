"""
Step 1 -- GPCR class (A, B1, B2, C, D1, F, O1, O2, T2).

Every protein chain of >= 100 CA is CE-aligned to a panel of class reference
structures; the class of the best valid 7TM alignment (>= 140 aligned residues,
longest alignment preferred, RMSD within 0.5 A) wins.

The algorithm is copied unchanged from pdb_classification_test/
1_gpcr_class_classification.py, which was benchmarked on 1,712 GPCRdb entries:
98.0% correct if class O counts as class A, 77.6% strictly -- step 1 calls many
class A receptors "ClassO1/ClassO2" (356 such calls on that set, 8 true odorant
receptors). Only the reference file lookup differs: files are pinned below
instead of found by globbing.
"""

import os

from pymol import cmd

refs_config = {
    "ClassA": [
        ("2rh1", "A", "Amine_ADRB2"),
        ("1u19", "A", "Opsin_RHO"),
        ("4dkl", "A", "Peptide_OPRM1"),
        ("5zkp", "A", "Lipid_S1PR1")
    ],
    "ClassB1": [
        ("4l6r", "A", "Secretin_GCGR"),
        ("6x18", "R", "Secretin_GLP1R")
    ],
    "ClassB2": [("7wu2", "R", "Adhesion_AGRD1")],
    "ClassC": [
        ("4or2", "A", "Glutamate_GRM1"),
        ("7dd5", "A", "Calcium_CaSR")
    ],
    "ClassD1": [("7ad3", "A", "Fungal_STE2")],
    "ClassF":  [("4jkv", "A", "Frizzled_SMO")],
    "ClassO1": [("8f76", "A", "Odorant_OR51E2")],
    "ClassO2": [("8uxy", "A", "Odorant_consOR1")],
    "ClassT2": [
        ("7xp5", "R", "Taste_TAS2R46"),
        ("9iiw", "R", "Taste_TAS2R14")
    ]
}

REF_FILES = {
    "2rh1": "ClassA_2rh1_A.pdb", "1u19": "ClassA_RHO_1U19_A.pdb", "4dkl": "ClassA_4dkl_A.pdb",
    "5zkp": "ClassA_5zkp_A.pdb", "4l6r": "ClassB1_GCGR_4l6r_A.pdb", "6x18": "ClassB1_GLP_6X18_R.pdb",
    "7wu2": "ClassB2_7wu2_R.pdb", "4or2": "ClassC_mGlu_4OR2_A.pdb", "7dd5": "ClassC_CaSR_7DD5_A.pdb",
    "7ad3": "ClassD1_7ad3_A.pdb", "4jkv": "ClassF_4jkv_A.pdb", "8f76": "ClassO1_8f76_A.pdb",
    "8uxy": "ClassO2_8uxy_A.pdb", "7xp5": "ClassT2_TAS2R46_7xp5_R.pdb", "9iiw": "ClassT2_TAS2R14_9IIW_R.cif",
}


def load_refs(ref_dir):
    for cls, ref_list in refs_config.items():
        for pdb, chain, name in ref_list:
            cmd.load(os.path.join(ref_dir, REF_FILES[pdb]), f"ref_{cls}_{name}")


def process_cif_structural(cif_path):
    filename = os.path.basename(cif_path)
    pdb_id = filename[:4].upper()

    try:
        cmd.load(cif_path, "target")
    except Exception:
        return None

    chains = cmd.get_chains("target")
    if not chains:
        cmd.delete("target")
        return None

    best_cls = None
    best_chain = None
    best_ref_name = ""
    best_rmsd = float('inf')
    best_aln_len = 0
    scores = {}

    # Evaluate every protein chain in the target structure
    for ch in chains:
        ca_count = cmd.count_atoms(f"target and chain {ch} and name CA")
        if ca_count < 100:  # Skip small peptides/nanobodies
            continue

        for cls, ref_list in refs_config.items():
            for pdb, chain, name in ref_list:
                ref_obj = f"ref_{cls}_{name}"
                try:
                    aln = cmd.cealign(ref_obj, f"target and chain {ch}")
                    rmsd = aln['RMSD']
                    aln_len = aln['alignment_length']

                    # Track best RMSD per class
                    if cls not in scores or rmsd < scores[cls]:
                        scores[cls] = rmsd

                    # For overall best class, require valid 7-TM alignment (>= 140 residues)
                    if aln_len >= 140:
                        if (aln_len > best_aln_len and rmsd <= best_rmsd + 0.5) or (aln_len == best_aln_len and rmsd < best_rmsd):
                            best_rmsd = rmsd
                            best_aln_len = aln_len
                            best_cls = cls
                            best_chain = ch
                            best_ref_name = f"{name}({pdb}_{chain})"
                except Exception:
                    pass

    cmd.delete("target")

    # Fallback if no chain passed >= 140 threshold
    fallback = False
    if best_cls is None and scores:
        best_cls = min(scores, key=scores.get)
        best_rmsd = scores[best_cls]
        fallback = True

    res_dict = {
        "pdb_id": pdb_id,
        "best_class": best_cls,
        "best_chain": best_chain if best_chain else "NA",
        "best_ref": best_ref_name,
        "best_rmsd": best_rmsd if best_rmsd != float('inf') else -1.0,
        "best_aln_len": best_aln_len,
        "no_valid_7tm_alignment": fallback,
    }
    for cls in refs_config.keys():
        res_dict[f"rmsd_{cls}"] = scores.get(cls, -1.0)

    return res_dict
