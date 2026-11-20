"""
Run the three steps: class -> GPCR chain -> state.
Only experimental PDB entries are accepted; other files are rejected before step 1.

Each step gets a fresh PyMOL session loaded with only its own
references, so it runs exactly as it did in the benchmark scripts.
"""

import os
import re
import sys

from . import class_step1 as class_step, chain_step2 as chain_step
from . import state_step3 as SE
from .transducer import entry_info

PKG_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(os.path.dirname(PKG_DIR), "data")
STATE_MODEL = os.path.join(DATA_DIR, "state_model.pkl")
DEFAULT_CACHE = os.path.join(os.path.dirname(PKG_DIR), "pdb_cache")
PDB_ID = re.compile(r"^[0-9][A-Za-z0-9]{3}$")


def _start_pymol():
    import pymol
    pymol.pymol_argv = ["pymol", "-qc"]
    pymol.finish_launching()
    from pymol import cmd
    return cmd


def _fresh_session(cmd, cache_dir):
    cmd.reinitialize()
    cmd.feedback("disable", "all", "everything")
    # PyMOL downloads definitions of residues it has no template for (e.g. UNK,
    # UNL) into fetch_path, which defaults to the current directory. Bond
    # assignment only; CA-based alignments are unaffected.
    cmd.set("fetch_path", cache_dir)


def fetch_cif(pdb_id, cache_dir):
    """Download <ID>.cif from RCSB into cache_dir (reused if already there)."""
    import requests
    if not PDB_ID.match(pdb_id):
        raise ValueError(f"not a PDB ID: {pdb_id!r}")
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.join(cache_dir, f"{pdb_id.lower()}.cif")
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path
    r = requests.get(f"https://files.rcsb.org/download/{pdb_id.upper()}.cif", timeout=60)
    r.raise_for_status()
    with open(path + ".part", "w") as f:
        f.write(r.text)
    os.replace(path + ".part", path)
    return path


def load_state_model(path=STATE_MODEL):
    """The trained state calibration: per-group activation cuts and the class A
    gradient-boosting model (see README)."""
    import pickle
    import sklearn
    with open(path, "rb") as f:
        params = pickle.load(f)
    trained = params.get("meta", {}).get("sklearn_version")
    if trained and trained != sklearn.__version__:
        print(f"warning: the state model was saved with scikit-learn {trained}, this is "
              f"{sklearn.__version__}; install scikit-learn=={trained} if results look wrong",
              file=sys.stderr)
    return params


def annotate(cif_paths, cache_dir=DEFAULT_CACHE, log=lambda msg: print(msg, file=sys.stderr)):
    """Returns one result dict per input path, in input order."""
    os.makedirs(cache_dir, exist_ok=True)
    cmd = _start_pymol()
    res = {}
    for p in cif_paths:
        entry_id, method, why = entry_info(p)
        res[p] = {"input": p, "pdb_id": entry_id or os.path.basename(p)[:4].upper(), "warnings": []}
        if not method:
            res[p]["error"] = why
    todo = [p for p in cif_paths if not res[p].get("error")]

    log("step 1/3: class")
    _fresh_session(cmd, cache_dir)
    class_step.load_refs(os.path.join(DATA_DIR, "ref_class"))
    for p in todo:
        r1 = class_step.process_cif_structural(p)
        if r1 is None:
            res[p]["error"] = "could not read the file"
            continue
        if r1["best_class"] is None:
            res[p]["error"] = "no protein chain of >= 100 residues could be aligned; not a GPCR entry"
            continue
        res[p]["gpcr_class"] = r1["best_class"]
        res[p]["class_ref"] = r1["best_ref"]
        if r1["no_valid_7tm_alignment"]:
            res[p]["warnings"].append("no chain gave a valid 7TM alignment - probably not a GPCR")
        if r1["best_class"] in ("ClassO1", "ClassO2"):
            res[p]["warnings"].append("class O calls are unreliable: most are class A receptors")

    log("step 2/3: GPCR chain")
    _fresh_session(cmd, cache_dir)
    chain_step.load_refs(os.path.join(DATA_DIR, "ref_chain"))
    for p in todo:
        if res[p].get("error"):
            continue
        _, chain, _, _, ref = chain_step.process_single_cif(p)
        res[p]["gpcr_chain"] = chain
        if chain is None:
            res[p]["error"] = ref
        elif str(ref).startswith("Fallback"):
            res[p]["warnings"].append("chain chosen by fallback (no chain with a valid 7TM alignment)")

    log("step 3/3: active / inactive")
    _fresh_session(cmd, cache_dir)
    SE.load_refs(os.path.join(DATA_DIR, "ref_state"))
    params = load_state_model()
    rows = {}
    for p in todo:
        if res[p].get("error"):
            continue
        rows[p] = SE.process(p, res[p]["gpcr_class"], res[p]["gpcr_chain"])
    SE.attach_probs(list(rows.values()), params)
    for p, r in rows.items():
        state, source, conf = SE.decide(r, params)
        res[p].update({"state": state, "evidence": source, "confidence": conf,
                       "p_active": r.get("p_active"), "partner": r.get("transducer_kind", ""),
                       "experimental_method": r.get("exptl_method", "")})
        if source.startswith("unmeasured"):
            res[p]["warnings"].append("state could not be measured; reported as Intermediate (low confidence)")
    return [res[p] for p in cif_paths]
