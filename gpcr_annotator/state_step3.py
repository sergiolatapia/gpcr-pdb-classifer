"""
Step 3 -- active / inactive state of the receptor chain (method: README.md).

Measurement code from pdb_classification_test/A3_state_refpair_eval.py; the
calibration (cuts and class A model) comes from data/state_model.pkl.
"""

import itertools
import os
import sys
from collections import defaultdict

import numpy as np
from pymol import cmd

from .transducer import detect

# key: (inactive file, chain, active file, chain)
REF_PAIRS = {
    "A_RHO": ("Inactive_ClassA_RHO_1U19.cif", "A", "Active_ClassA_RHO_3PQR.cif", "A"),
    "A_Gi":  ("Inactive_ClassA_Gi_4DKL.cif", "A", "Active_ClassA_Gi_6DDE.cif", "R"),
    "A_Gs":  ("Inactive_ClassA_Gs_2VT4.cif", "B", "Active_ClassA_Gs_7JJO.cif", "R"),
    "B1":    ("Inactive_ClassB1_5VEW.cif", "A", "Active_ClassB1_6X18.cif", "R"),
    "B2":    ("Inactive_ClassB2_8IKJ.cif", "R", "Active_ClassB2_8IKL.cif", "R"),
    "D":     ("Inactive_ClassD_7QA8.cif", "A", "Active_ClassD_7AD3.cif", "A"),
    "F_SMO": ("Inactive_ClassF_SMO_5L7D.cif", "A", "Active_ClassF_SMO_6OT0.cif", "R"),
    "F_FZD": ("Inactive_ClassF_FZD_6BD4.cif", "A", "Active_ClassF_FZD_8WMA.cif", "A"),
}
CLASS_C_PAIR = ("Inactive_ClassC_7MTQ.cif", "A", "Active_ClassC_7MTS.cif", "A")
TMD_TEMPLATE = ("Inactive_ClassA_RHO_1U19.cif", "A")   # rhodopsin 7TM, used to bracket the TM domain

CLASS_TO_REFS = {
    "ClassA":  ["A_RHO", "A_Gi", "A_Gs"],
    "ClassO1": ["A_RHO", "A_Gi", "A_Gs"],
    "ClassO2": ["A_RHO", "A_Gi", "A_Gs"],
    "ClassB1": ["B1", "B2"],
    "ClassB2": ["B1", "B2"],
    "ClassD1": ["D"],
    "ClassD":  ["D"],
    "ClassF":  ["F_SMO", "F_FZD"],
}   # ClassC: dimer measurement; ClassT2: no inactive reference, not measured

CORE_CUT = 2.0         # A; reference-pair residues closer than this form the static core
MOBILE_CUT = 4.0       # A; residues displaced more than this form the mobile set
MIN_BLOCK = 5          # shortest contiguous mobile stretch kept
MIN_CORE_PAIRS = 40
MAX_EXTRAP = 20        # residues a mobile probe may be walked from its core anchor
GAP_TOL = 10           # unmodelled residues the walk may cross
MIN_MOBILE_PAIRS = 5
MIN_RECEPTOR_CA = 150

QC_MAX_PERP = 0.60     # quality gate for a structural reading
QC_MIN_MOBILE = 8
QC_MAX_CORE_RMSD = 3.5

RHO_TM4, RHO_TM6 = 161, 267   # class C probe points: rhodopsin W4.50 and P6.50 equivalents
CLASSC_WIN = 2
CLASSC_MAX_SEP = 45.0  # A; farther chains are not a dimer partner
CLASSC_MAX_RMS = 8.0
CLASSC_MIN_ALIGN = 120

MEAS_GROUP = {"ClassA": "A", "ClassO1": "A", "ClassO2": "A", "ClassB1": "B", "ClassB2": "B",
              "ClassC": "C", "ClassD1": "D", "ClassD": "D", "ClassF": "F"}
DEFAULT_THRESHOLD = 0.5        # groups without a fitted cut
AMBIGUOUS_P_BAND = 0.10        # class A: |P(active) - 0.5| below this -> Intermediate
UNMEASURABLE_RECEPTORS = {"GP156_HUMAN": "GPR156: no inactive structure to calibrate against"}

AXES = {}


# ---------------------------------------------------------------- geometry

def kabsch(P, Q):
    """Rigid transform taking P onto Q: (R, t, rmsd)."""
    cp, cq = P.mean(0), Q.mean(0)
    H = (P - cp).T @ (Q - cq)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    t = cq - R @ cp
    fit = (R @ P.T).T + t
    return R, t, float(np.sqrt(((fit - Q) ** 2).sum(1).mean()))


def iterative_core(P, Q, cut=CORE_CUT, iters=12, floor=20):
    """Superpose P onto Q on the pairs that fit, so moving parts don't bias the fit."""
    mask = np.ones(len(P), bool)
    for _ in range(iters):
        R, t, _ = kabsch(P[mask], Q[mask])
        new = np.linalg.norm((R @ P.T).T + t - Q, axis=1) < cut
        if new.sum() < floor or (new == mask).all():
            break
        mask = new
    R, t, _ = kabsch(P[mask], Q[mask])
    return mask, R, t


def contiguous_blocks(nums, gap=4):
    blocks, cur = [], [nums[0]]
    for x in nums[1:]:
        if x - cur[-1] <= gap:
            cur.append(x)
        else:
            blocks.append(cur)
            cur = [x]
    blocks.append(cur)
    return blocks


def median_offset(pairs):
    offs = sorted(q - r for r, q in pairs)
    return offs[len(offs) // 2]


# ---------------------------------------------------------------- PyMOL helpers

def _ca_object(path, chain, obj):
    """CA-only object for one chain. Drops coordinate-less placeholders (occupancy 0),
    and residues >= 1000 (fusion partners) when the receptor itself is numbered below 1000."""
    cmd.delete("_raw")
    cmd.load(path, "_raw")
    base = f"_raw and chain {chain} and polymer.protein and name CA and alt ''+A"
    sel = f"{base} and q > 0.01"
    if cmd.count_atoms(sel) == 0:
        sel = base
    if cmd.count_atoms(f"({sel}) and not resi 1000-") >= MIN_RECEPTOR_CA:
        sel = f"({sel}) and not resi 1000-"
    cmd.delete(obj)
    cmd.create(obj, sel)
    cmd.delete("_raw")
    return cmd.count_atoms(obj)


def _residues(obj):
    out = []
    cmd.iterate_state(1, obj, "out.append(resv)", space={"out": out})
    return out


def _coords(obj):
    return np.asarray(cmd.get_coords(obj, 1), dtype=float)


def _cealign_pairs(ref_obj, mob_obj):
    """(pairs, rmsd, n_aligned); pairs index into each object's coordinate array."""
    cmd.delete("_aln")
    r = cmd.cealign(ref_obj, mob_obj, object="_aln")
    raw = cmd.get_raw_alignment("_aln")

    def idxmap(obj):
        ids = []
        cmd.iterate_state(1, obj, "ids.append(index)", space={"ids": ids})
        return {ix: i for i, ix in enumerate(ids)}

    mr, mm = idxmap(ref_obj), idxmap(mob_obj)
    pairs = []
    for col in raw:
        d = dict(col)
        if ref_obj in d and mob_obj in d:
            a, b = mr.get(d[ref_obj]), mm.get(d[mob_obj])
            if a is not None and b is not None:
                pairs.append((a, b))
    cmd.delete("_aln")
    return pairs, float(r["RMSD"]), int(r["alignment_length"])


# ---------------------------------------------------------------- reference axes

def build_axis(key, f_inactive, ch_inactive, f_active, ch_active, refs_dir):
    oi, oa = f"{key}_I", f"{key}_A"
    _ca_object(os.path.join(refs_dir, f_inactive), ch_inactive, oi)
    _ca_object(os.path.join(refs_dir, f_active), ch_active, oa)

    # 7TM window: union of both members' alignment to the rhodopsin 7TM
    lo, hi = None, None
    for obj in (oi, oa):
        pairs, _, _ = _cealign_pairs("TMD_TEMPLATE", obj)
        rl = _residues(obj)
        hit = sorted(rl[b] for _, b in pairs)
        if len(hit) < 100:
            raise RuntimeError(f"{key}: only {len(hit)} residues bracket the 7TM")
        lo = hit[0] if lo is None else min(lo, hit[0])
        hi = hit[-1] if hi is None else max(hi, hit[-1])
    for obj in (oi, oa):
        cmd.create(obj, f"{obj} and resi {lo}-{hi}")

    ri, ra = _residues(oi), _residues(oa)
    ci, ca = _coords(oi), _coords(oa)
    pos_a = {n: k for k, n in enumerate(ra)}          # same receptor: author numbering matches
    idx = [(k, pos_a[n]) for k, n in enumerate(ri) if n in pos_a]
    if len(idx) < 80:
        raise RuntimeError(f"{key}: only {len(idx)} residues shared between the pair")
    P = ca[[j for _, j in idx]]
    Q = ci[[k for k, _ in idx]]
    resnums = np.array([ri[k] for k, _ in idx])

    mask, R, t = iterative_core(P, Q)
    A_fit = (R @ P.T).T + t
    mob = np.linalg.norm(A_fit - Q, axis=1) > MOBILE_CUT
    if mob.sum() < MIN_BLOCK:
        raise RuntimeError(f"{key}: no mobile region above {MOBILE_CUT} A")
    blocks = [b for b in contiguous_blocks(sorted(resnums[mob].tolist())) if len(b) >= MIN_BLOCK]
    if not blocks:
        raise RuntimeError(f"{key}: mobile residues too scattered")
    best = max(blocks, key=len)                        # the cytoplasmic end of TM6
    mob_set = set(best)
    return {
        "key": key,
        "obj_inactive": oi,
        "core_resnums": set(resnums[mask].tolist()),
        "mobile_resnums": list(best),
        "mob_xyz_I": {int(n): Q[i] for i, n in enumerate(resnums) if n in mob_set},
        "mob_xyz_A": {int(n): A_fit[i] for i, n in enumerate(resnums) if n in mob_set},
    }


def classc_protomer(path, chain, tag):
    """TM4 and TM6 probe points of one protomer, in the file's own frame."""
    obj = f"C_{tag}"
    try:
        n = _ca_object(path, chain, obj)
    except Exception:
        return None
    if n < MIN_RECEPTOR_CA:
        cmd.delete(obj)
        return None

    # align a copy: cealign moves the mobile object and would erase the dimer geometry
    # (no leading underscore: PyMOL hides such objects from get_raw_alignment)
    cmd.delete("PROBE")
    cmd.create("PROBE", obj)
    try:
        pairs, rms, al = _cealign_pairs("TMD_TEMPLATE", "PROBE")
    except Exception:
        cmd.delete("PROBE")
        cmd.delete(obj)
        return None
    tr, rp = _residues("TMD_TEMPLATE"), _residues("PROBE")
    aligned = sorted((tr[a], rp[b]) for a, b in pairs)
    cmd.delete("PROBE")
    if al < CLASSC_MIN_ALIGN or rms > CLASSC_MAX_RMS:
        cmd.delete(obj)
        return None

    rq, cq = _residues(obj), _coords(obj)
    qpos = {v: i for i, v in enumerate(rq)}
    out = {"obj": obj}
    for name, anchor in (("TM4", RHO_TM4), ("TM6", RHO_TM6)):
        near = [(t, q) for t, q in aligned if abs(t - anchor) <= 8]
        if len(near) < 3:
            near = sorted(aligned, key=lambda pr: abs(pr[0] - anchor))[:5]
            if not near or min(abs(pr[0] - anchor) for pr in near) > 25:
                cmd.delete(obj)
                return None
        off = median_offset(near)
        got = []
        for win in (CLASSC_WIN, CLASSC_WIN + 2, CLASSC_WIN + 4):
            got = [qpos[anchor + d + off] for d in range(-win, win + 1) if anchor + d + off in qpos]
            if len(got) >= 3:
                break
        if len(got) < 3:
            cmd.delete(obj)
            return None
        out[name] = cq[got].mean(0)
    return out


def _closest_pair(me, others):
    """(sep, chain, d4, d6) of the protomer closest to `me`, or None."""
    best = None
    for c, p in others:
        d4 = float(np.linalg.norm(me["TM4"] - p["TM4"]))
        d6 = float(np.linalg.norm(me["TM6"] - p["TM6"]))
        if best is None or min(d4, d6) < best[0]:
            best = (min(d4, d6), c, d4, d6)
    return best


def build_classc_axis(refs_dir):
    """delta = d(TM6-TM6) - d(TM4-TM4) of the inactive and active class C references."""
    f_i, ch_i, f_a, ch_a = CLASS_C_PAIR
    out = {}
    for tag, fname in (("I", f_i), ("A", f_a)):
        path = os.path.join(refs_dir, fname)
        cmd.delete("_raw")
        cmd.load(path, "_raw")
        chains = list(cmd.get_chains("_raw"))
        cmd.delete("_raw")
        prot = {}
        for ch in chains:
            p = classc_protomer(path, ch, f"REF{tag}{ch}")
            if p:
                prot[ch] = p
        if len(prot) < 2:
            raise RuntimeError(f"class C reference {fname}: fewer than 2 readable protomers")
        best = None
        for a, b in itertools.combinations(sorted(prot), 2):
            cand = _closest_pair(prot[a], [(b, prot[b])])
            if best is None or cand[0] < best[0]:
                best = cand
        for p in prot.values():
            cmd.delete(p["obj"])
        out[tag] = best[3] - best[2]
    return {"delta_I": out["I"], "delta_A": out["A"]}


def load_refs(refs_dir):
    _ca_object(os.path.join(refs_dir, TMD_TEMPLATE[0]), TMD_TEMPLATE[1], "TMD_TEMPLATE")
    AXES.clear()
    for key, files in REF_PAIRS.items():
        try:
            AXES[key] = build_axis(key, *files, refs_dir)
        except Exception as exc:
            print(f"[axis {key}] FAILED: {exc}", file=sys.stderr)
    try:
        AXES["C_dimer"] = build_classc_axis(refs_dir)
    except Exception as exc:
        print(f"[axis C_dimer] FAILED: {exc}", file=sys.stderr)


# ---------------------------------------------------------------- measuring a query

def score_class_c(cif_path, chain, axis):
    """Activation from the dimer interface: (result, error)."""
    cmd.delete("_raw")
    try:
        cmd.load(cif_path, "_raw")
        chains = list(cmd.get_chains("_raw"))
    except Exception as exc:
        return None, f"load: {exc}"
    finally:
        cmd.delete("_raw")

    me = classc_protomer(cif_path, chain, "SELF")
    if me is None:
        return None, f"chain {chain} does not present a readable 7TM"
    best = None
    for c in chains:
        if c == chain:
            continue
        p = classc_protomer(cif_path, c, "PART")
        if p is None:
            continue
        cand = _closest_pair(me, [(c, p)])
        cmd.delete(p["obj"])
        if cand[0] <= CLASSC_MAX_SEP and (best is None or cand[0] < best[0]):
            best = cand
    cmd.delete(me["obj"])
    if best is None:
        return None, "no partner protomer within range (monomer, ECD-only or TMD-only entry)"
    delta = best[3] - best[2]
    activation = float((delta - axis["delta_I"]) / (axis["delta_A"] - axis["delta_I"]))
    return {"ref": "C_dimer", "activation": activation, "qc_pass": True}, None


def transfer_mobile(axis, ri, rq, core_pairs):
    """Map the reference's mobile residues onto the query by walking along the chain
    from nearby core anchors (CEalign is unreliable inside the moving helix)."""
    ipos = {n: i for i, n in enumerate(ri)}
    qpos = {n: i for i, n in enumerate(rq)}
    anchors = sorted((ri[a], rq[b]) for a, b in core_pairs)
    if not anchors:
        return {}
    out = {}
    for m in axis["mobile_resnums"]:
        near = [(rn, qn) for rn, qn in anchors if abs(rn - m) <= MAX_EXTRAP]
        if len(near) < 3:
            near = sorted(anchors, key=lambda p: abs(p[0] - m))[:3]
        off = median_offset(near)
        rn, qn = min(near, key=lambda p: abs(p[0] - m))
        if abs(m - rn) > 2 * MAX_EXTRAP:
            continue
        qm = m + off
        if qm not in qpos:
            continue
        # numbering must be (nearly) unbroken on both sides over the whole walk
        lo, hi = sorted((rn, m))
        if sum(1 for x in range(lo, hi + 1) if x not in ipos) > GAP_TOL:
            continue
        qlo, qhi = sorted((qn, qm))
        if sum(1 for x in range(qlo, qhi + 1) if x not in qpos) > GAP_TOL:
            continue
        out[m] = qpos[qm]
    return out


def score_against(axis, qobj):
    pairs, _, _ = _cealign_pairs(axis["obj_inactive"], qobj)
    if not pairs:
        return None
    ri, ci = _residues(axis["obj_inactive"]), _coords(axis["obj_inactive"])
    rq, cq = _residues(qobj), _coords(qobj)

    core_pairs = [(a, b) for a, b in pairs if ri[a] in axis["core_resnums"]]
    if len(core_pairs) < MIN_CORE_PAIRS:
        return None
    R, t, core_rmsd = kabsch(cq[[b for _, b in core_pairs]], ci[[a for a, _ in core_pairs]])

    mmap = transfer_mobile(axis, ri, rq, core_pairs)
    if len(mmap) < MIN_MOBILE_PAIRS:
        return None
    keys = sorted(mmap)
    Qm = (R @ cq[[mmap[k] for k in keys]].T).T + t
    ref_i = np.array([axis["mob_xyz_I"][k] for k in keys])
    ref_a = np.array([axis["mob_xyz_A"][k] for k in keys])

    c_I, c_A = ref_i.mean(0), ref_a.mean(0)
    span = float(np.linalg.norm(c_A - c_I))
    if span < 1e-6:
        return None
    u = (c_A - c_I) / span
    v = Qm.mean(0) - c_I
    d_i = float(np.sqrt(((Qm - ref_i) ** 2).sum(1).mean()))
    d_a = float(np.sqrt(((Qm - ref_a) ** 2).sum(1).mean()))
    return {
        "ref": axis["key"],
        "activation": float(v @ u / span),                       # 0 = inactive ref, 1 = active ref
        "perp": float(np.linalg.norm(v - (v @ u) * u) / span),   # distance off the axis
        "rmsd_score": d_i / (d_i + d_a) if (d_i + d_a) > 1e-6 else 0.5,
        "core_rmsd": core_rmsd,
        "n_core": len(core_pairs),
        "n_mobile": len(keys),
    }


def _qc(r):
    return r["perp"] <= QC_MAX_PERP and r["n_mobile"] >= QC_MIN_MOBILE and r["core_rmsd"] <= QC_MAX_CORE_RMSD


def process(cif_path, cls, chain):
    """Partner check + structural reading for the receptor chain."""
    base = {"class": cls, "chain": chain, **detect(cif_path, chain)}
    if cls == "ClassC":
        ax = AXES.get("C_dimer")
        if ax is None:
            return {**base, "error": "class C dimer axis failed to build"}
        res, err = score_class_c(cif_path, chain, ax)
        return {**base, "error": err} if err else {**base, **res, "error": None}

    ref_keys = [k for k in CLASS_TO_REFS.get(cls, []) if k in AXES]
    if not ref_keys:
        return {**base, "error": f"no reference pair for {cls}"}
    qobj = "QRY"
    try:
        n = _ca_object(cif_path, chain, qobj)
        if n < 100:
            return {**base, "error": f"chain {chain} yields only {n} usable CA atoms"}
    except Exception as exc:
        return {**base, "error": f"load: {exc}"}

    cands = []
    for k in ref_keys:
        try:
            r = score_against(AXES[k], qobj)
        except Exception:
            r = None
        if r:
            cands.append(r)
    cmd.delete(qobj)
    if not cands:
        return {**base, "error": "no usable alignment to any reference"}

    # report the reference giving a valid reading, then the best core fit
    best = min(cands, key=lambda r: (0, r["core_rmsd"], -r["n_core"]) if _qc(r) else (1, r["perp"], -r["n_core"]))
    acts = sorted(r["activation"] for r in cands)
    mid = len(acts) // 2
    consensus = acts[mid] if len(acts) % 2 else 0.5 * (acts[mid - 1] + acts[mid])
    return {**base, **best,
            "qc_pass": bool(_qc(best)),
            "activation_consensus": float(consensus),
            "activation_spread": float(acts[-1] - acts[0]),
            **{f"act_{r['ref']}": r["activation"] for r in cands},
            "error": None}


# ---------------------------------------------------------------- decision

def _model_X(rows, feats):
    return np.array([[np.nan if r.get(f) is None else float(r.get(f)) for f in feats] for r in rows], float)


def _measured(r):
    return not (r.get("transducer") or r.get("error") or r.get("activation") is None
                or r.get("receptor_uniprot") in UNMEASURABLE_RECEPTORS)


def attach_probs(rows, params):
    """Set r["p_active"] where the class A model decides (batched); None elsewhere."""
    by = defaultdict(list)
    for r in rows:
        r["p_active"] = None
        g = MEAS_GROUP.get(r.get("class"), r.get("class"))
        if _measured(r) and g in params["models"]:
            by[g].append(r)
    for g, rs in by.items():
        m = params["models"][g]
        for r, p in zip(rs, m["estimator"].predict_proba(_model_X(rs, m["features"]))[:, 1]):
            r["p_active"] = float(p)


def decide(r, params):
    """(state, evidence, confidence) for a row from process(); call attach_probs first."""
    if r.get("transducer"):
        return "Active", "partner:" + r.get("transducer_kind", ""), "high"
    if not _measured(r):
        why = UNMEASURABLE_RECEPTORS.get(r.get("receptor_uniprot")) or r.get("error") or "no structural reading"
        return "Intermediate", "unmeasured: " + why, "low"
    conf = "high" if r.get("qc_pass") else "low"
    p = r.get("p_active")
    if p is not None:
        if abs(p - 0.5) < AMBIGUOUS_P_BAND:
            return "Intermediate", f"structure (no partner), borderline p_active={p:.2f}", conf
        return ("Active" if p > 0.5 else "Inactive"), f"structure (no partner), model p_active={p:.2f}", conf
    group = MEAS_GROUP.get(r["class"], r["class"])
    t = params["cuts"].get(group, DEFAULT_THRESHOLD)
    state = "Active" if r["activation"] > t else "Inactive"
    return state, f"structure (no partner), activation {r['activation']:.2f} vs cut {t:.2f}", conf
