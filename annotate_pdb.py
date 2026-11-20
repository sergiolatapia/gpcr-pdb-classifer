#!/usr/bin/env python3
"""
Command-line entry point: collects the inputs (PDB IDs, mmCIF files, folders), runs
the three steps in gpcr_annotator/pipeline.py on them, and writes the output.
Experimental PDB entries only; predicted models (e.g. AlphaFold) are rejected.

    python annotate_pdb.py 7MTS                      # download from RCSB and annotate
    python annotate_pdb.py 2RH1 3SN6 path/to/x.cif   # IDs and local mmCIF files can be mixed
    python annotate_pdb.py path/to/folder/           # every .cif in a folder
    python annotate_pdb.py 2RH1 7MTS -o my_results   # choose the output folder

Output folder (default ./gpcr_annotation_output):
    summary.tsv              PDB, class, chain, state, confidence, warning
    <PDB>_<state>.cif        a copy of each input structure, e.g. 9IIW_inactive.cif

State is Active, Inactive or Intermediate (low confidence or not measurable).
See README.md for how each is decided and accuracy expectations.
"""

import argparse
import os
import sys

from gpcr_annotator.pipeline import annotate, fetch_cif, PDB_ID

COLS = ["pdb_id", "gpcr_class", "gpcr_chain", "state", "confidence", "warning"]

WARNING = ("This structure might not be a GPCR, or is a GPCR of an underrepresented "
           "subfamily. State prediction is unreliable.")


def needs_warning(r):
    """No valid 7TM alignment (class or chain fell back), or the state could not be measured."""
    return (any("7TM" in w for w in r.get("warnings", []))
            or str(r.get("evidence", "")).startswith("unmeasured"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("entries", nargs="+", help="PDB IDs (fetched from RCSB), mmCIF files of PDB entries and/or folders of them")
    ap.add_argument("-o", "--outdir", default="gpcr_annotation_output",
                    help="output folder (default: ./gpcr_annotation_output)")
    ap.add_argument("--cache", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "pdb_cache"),
                    help="where downloaded mmCIF files are kept (default: ./pdb_cache)")
    args = ap.parse_args()

    paths = []
    for e in args.entries:
        if os.path.isdir(e):
            found = sorted(os.path.join(e, f) for f in os.listdir(e) if f.lower().endswith(".cif"))
            if not found:
                print(f"{e}: folder has no .cif files, skipped", file=sys.stderr)
            paths.extend(os.path.abspath(f) for f in found)
        elif os.path.isfile(e):
            paths.append(os.path.abspath(e))
        elif PDB_ID.match(e):
            try:
                paths.append(fetch_cif(e, args.cache))
            except Exception as exc:
                print(f"{e}: download failed ({exc})", file=sys.stderr)
        else:
            print(f"{e}: neither a file nor a PDB ID, skipped", file=sys.stderr)
    paths = list(dict.fromkeys(paths))          # the same file given twice is annotated once
    if not paths:
        sys.exit("nothing to annotate")

    results = annotate(paths, cache_dir=args.cache)

    print(f"\n{'PDB':<6}{'class':<9}{'chain':<7}{'state':<14}confidence")
    for r in results:
        if r.get("error"):
            print(f"{r['pdb_id']:<6}ERROR: {r['error']}")
            continue
        print(f"{r['pdb_id']:<6}{r['gpcr_class']:<9}{r['gpcr_chain']:<7}{r['state']:<14}"
              f"{r.get('confidence') or '-'}")
        if needs_warning(r):
            print(f"{'':<6}WARNING: {WARNING}")

    import shutil
    os.makedirs(args.outdir, exist_ok=True)
    with open(os.path.join(args.outdir, "summary.tsv"), "w") as f:
        f.write("\t".join(COLS) + "\n")
        for r in results:
            if r.get("error"):
                f.write(f"{r['pdb_id']}\t\t\tERROR\t\t{r['error']}\n")
                continue
            row = dict(r, warning=WARNING if needs_warning(r) else "")
            f.write("\t".join("" if row.get(c) is None else str(row.get(c)) for c in COLS) + "\n")
            shutil.copy2(r["input"], os.path.join(args.outdir, f"{r['pdb_id']}_{r['state'].lower()}.cif"))
    print(f"\nresults in: {os.path.abspath(args.outdir)}/ (summary.tsv + one renamed .cif per structure)")


if __name__ == "__main__":
    main()
