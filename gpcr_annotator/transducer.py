"""
transducer.py -- is an active-state signalling partner bound in this PDB entry?
(Used by step 3; also checks that the input is an experimental PDB entry.)

A G protein, an arrestin or a GRK can only couple to an ACTIVE receptor, so their
presence in a deposited structure is a can be a clue of the receptor state.
This is checked first; the structural TM6 measurement is only needed for entries
without such a partner (apo, ligand-only, or antibody-bound structures).

Only partners that REQUIRE an active receptor are counted:
  G_alpha          any G alpha (UniProt GNA*), mini-G, engineered/chimeric G alpha,
                   and G alpha C-terminal peptides
  Nb35/scFv16      the two G-protein-binding antibody fragments of cryo-EM complexes
  arrestin         beta-arrestin 1/2, visual (S-) and cone arrestin
  GRK              G protein-coupled receptor kinases
  active_nanobody  intracellular nanobodies published as active-state stabilisers
                   (Nb80, Nb6B9, Nb9-8, Nb39, Nb.AT110i1)

Deliberately NOT counted, because they can be found in inactive structures too:
  G beta / G gamma on their own (GNB5 + RGS7 is bound to INACTIVE GPR158/GPR179),
  generic "nanobody"/"scFv"/Fab names (Nb6, Nb60, Nb2, Nb52 and most Fabs stabilise
  INACTIVE receptors; "Nanobody 6" and "Nanobody6" are different proteins with
  opposite effects), RAMPs, DARPins, BRIL-binding Fabs.

Measured on the 1,692 two-state entries of the GPCRdb benchmark: fires on 1,022/1,232
Active and 0/460 Inactive entries. The whitelist was assembled while looking at that
set, but every member is there for a mechanistic reason, not a fitted one.

Only experimental PDB entries are annotated (mmCIF _exptl.method present).
entry_info() is checked before any step runs; files without an experimental
method (predicted models, coordinates-only files) are rejected.

Dependencies: biopython.
"""

import re
import warnings

UNP_GALPHA = re.compile(r"^GNA(S|I|O|Q|L|Z|T|1[0-6])", re.I)
UNP_ARRESTIN = re.compile(r"^ARR(B1|B2|S|C)_", re.I)
UNP_GRK = re.compile(r"^(GRK\d|ARBK\d|RK)_", re.I)

DESC_NB35_SCFV16 = re.compile(r"\bnb ?-?35\b|nanobody[- ]?35|scfv[- ]?16|svfv[- ]?16", re.I)
DESC_ACTIVE_NB = re.compile(
    r"\bnb ?-?80\b|nanobody[- ]?80|nb6b9|\bnb ?9-8\b|nanobody 9-8|\bnb ?-?39\b|nanobody[- ]?39|nb\.?at110", re.I)
DESC_ARRESTIN = re.compile(r"arrestin", re.I)
DESC_GRK = re.compile(r"receptor kinase|rhodopsin kinase|adrenergic receptor kinase", re.I)
# "guanine nucleotide-binding protein" unless the same clause says beta/gamma subunit,
# plus the many names engineered G alphas are deposited under.
DESC_GALPHA = re.compile(
    r"guanine nucleotide[- ,]*binding protein(?![^,]*subunit (beta|gamma))(?![^,]*\b(beta|gamma)[- ]?\d)"
    r"|mini[- ]?g|\bminig|\bmg(s|i|o|q)\b|gs-mini|g-?alpha|galpha|alpha subunit of gs|gs protein alpha"
    r"|g protein subunit (alpha|q|s|13)|g subunit q|g protein alpha|\bg protein$|g protein ct2|gnas complex",
    re.I)
RECEPTOR_NAME = re.compile(r"receptor(?! kinase)|smoothened|frizzled|opsin|gpr\d", re.I)

# Crystallisation / expression fusion partners: never the receptor's own UniProt entry.
FUSION_UNP = re.compile(
    r"^(LYS_BPT4|ENLYS_BPT4|C562_ECOLX|RUBR_CLOPA|MALE_ECOLI|MALE_ECO57|GLGA_PYRAB|FLAV_DESVH"
    r"|FLAV_DESDE|TRXA_ECOLI|CYB5_|GFP_|B2MG_|SUMO|PSBQ|DAPB)", re.I)


def _as_list(d, key):
    v = d.get(key, [])
    return v if isinstance(v, list) else [v]


def _read_cif(cif_path):
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return MMCIF2Dict(cif_path)


def experimental_method(d):
    """The deposited experimental method ('X-RAY DIFFRACTION', 'ELECTRON MICROSCOPY', ...),
    or "" for a computational model (AlphaFold, ColabFold, ModelCIF, coordinates-only file)."""
    if any(k.startswith("_ma_") for k in d):          # ModelCIF: AlphaFold DB / AF3 server
        return ""
    methods = [m for m in _as_list(d, "_exptl.method") if m not in ("?", ".", "")]
    if any("predict" in m.lower() or "model" in m.lower() for m in methods):
        return ""
    return ";".join(methods)


def entry_info(cif_path):
    """(PDB ID, experimental method, "") for an experimental PDB entry;
    (PDB ID or "", "", reason) otherwise. The ID is _entry.id from the file."""
    try:
        d = _read_cif(cif_path)
    except Exception as exc:
        return "", "", f"could not read the mmCIF header ({exc})"
    ids = [x for x in _as_list(d, "_entry.id") if x not in ("?", ".", "")]
    entry_id = ids[0].upper() if ids else ""
    method = experimental_method(d)
    if not method:
        return entry_id, "", ("not an experimental PDB entry (no experimental method in the file); "
                              "predicted models are not supported")
    return entry_id, method, ""


def read_entities(cif_path, d=None):
    """Polymer entities of an mmCIF file: {entity_id: {desc, chains, refs}}."""
    if d is None:
        d = _read_cif(cif_path)
    ents = {}
    for i, typ, desc in zip(_as_list(d, "_entity.id"), _as_list(d, "_entity.type"),
                            _as_list(d, "_entity.pdbx_description")):
        if typ == "polymer":
            ents[i] = {"desc": desc, "chains": [], "refs": []}
    for i, strands in zip(_as_list(d, "_entity_poly.entity_id"), _as_list(d, "_entity_poly.pdbx_strand_id")):
        if i in ents:
            ents[i]["chains"] = strands.split(",")
    for i, db, code in zip(_as_list(d, "_struct_ref.entity_id"), _as_list(d, "_struct_ref.db_name"),
                           _as_list(d, "_struct_ref.db_code")):
        if i in ents and db == "UNP":
            ents[i]["refs"].append(code)
    return ents


def classify_entity(ent):
    """Active-state partner category of one polymer entity, or None."""
    desc, refs = ent["desc"], ent["refs"]
    for r in refs:
        if UNP_GALPHA.match(r):
            return "G_alpha"
        if UNP_ARRESTIN.match(r):
            return "arrestin"
        if UNP_GRK.match(r):
            return "GRK"
    if DESC_NB35_SCFV16.search(desc):
        return "Nb35/scFv16"
    if DESC_ACTIVE_NB.search(desc):
        return "active_nanobody"
    if DESC_ARRESTIN.search(desc):
        return "arrestin"
    if DESC_GRK.search(desc):
        return "GRK"
    if DESC_GALPHA.search(desc):
        # a receptor entity whose name merely contains "G protein" is not a G protein
        if RECEPTOR_NAME.search(desc) and not re.search(r"guanine nucleotide", desc, re.I):
            return None
        return "G_alpha"
    return None


def receptor_uniprot(ents, chain):
    """UniProt entry name of the receptor chain (fusion partners skipped), or the entity name."""
    for ent in ents.values():
        if chain in ent["chains"]:
            for r in ent["refs"]:
                if not FUSION_UNP.match(r):
                    return r.upper()
            return re.sub(r"\W+", "_", ent["desc"].lower())[:40]
    return ""


def detect(cif_path, receptor_chain=None):
    """
    Returns {"experimental": bool, "exptl_method": str,
             "transducer": bool, "transducer_kind": "G_alpha;Nb35/scFv16"|"",
             "transducer_desc": str, "receptor_uniprot": str}.

    The partner check is only meaningful for experimental entries, where a missing
    partner reflects what was actually deposited; without an experimental method
    "transducer" is always False. Never raises.
    """
    out = {"experimental": False, "exptl_method": "", "transducer": False,
           "transducer_kind": "", "transducer_desc": "", "receptor_uniprot": ""}
    try:
        d = _read_cif(cif_path)
        ents = read_entities(cif_path, d)
    except Exception:
        return out
    out["exptl_method"] = experimental_method(d)
    out["experimental"] = bool(out["exptl_method"])
    if receptor_chain:
        out["receptor_uniprot"] = receptor_uniprot(ents, receptor_chain)
    if not out["experimental"]:
        return out
    kinds, descs = set(), []
    for ent in ents.values():
        k = classify_entity(ent)
        if k:
            kinds.add(k)
            descs.append(ent["desc"])
    out["transducer"] = bool(kinds)
    out["transducer_kind"] = ";".join(sorted(kinds))
    out["transducer_desc"] = "; ".join(descs[:3])
    return out


if __name__ == "__main__":
    import sys
    for p in sys.argv[1:]:
        print(p, detect(p))
