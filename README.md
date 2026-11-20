# gpcr_pdb_annotator

Given a GPCR structure from the PDB, report:

1. its **GPCR class** (A, B1, B2, C, D1, F, O1, O2, T2),
2. which **chain** is the receptor, and
3. whether the receptor is **Active**, **Inactive**, or **Intermediate** (borderline).

It is meant for annotating new PDB releases as they appear, so the AF MultiState
template database can be updated without waiting for GPCRdb.

Only **experimental PDB entries** are accepted. A file without an experimental method
(for example an AlphaFold or other predicted model) is reported as an error and not
annotated.

## Install

```bash
conda install -c conda-forge pymol-open-source
pip install -r requirements.txt
```

`requirements.txt` pins scikit-learn to the version the state model was saved with
(see [Calibration](#calibration)).

## Use

```bash
python annotate_pdb.py 7MTS                      # PDB ID: downloaded from RCSB into ./pdb_cache
python annotate_pdb.py 2RH1 3SN6 path/to/x.cif   # IDs and local mmCIF files can be mixed
python annotate_pdb.py path/to/folder/           # every .cif in a folder
python annotate_pdb.py 2RH1 7MTS -o my_results   # choose the output folder
```

```
PDB   class    chain  state         confidence
2RH1  ClassA   A      Inactive      high
3SN6  ClassA   R      Active        high
7MTQ  ClassC   B      Inactive      high
```

Results go to an output folder (default `./gpcr_annotation_output`):

- `summary.tsv`: PDB, class, chain, state, confidence, warning;
- `<PDB>_<state>.cif`: a copy of each input structure, e.g. `2RH1_inactive.cif`.

Runtime is about 5 s per entry on one CPU core.

## How each answer is reached

**Class.** Every protein chain is structurally aligned (PyMOL CEalign) to 15
class reference structures. The class with the best valid 7TM alignment wins.

**Chain.** Every chain is aligned to a 13-structure 7TM panel. Among chains with a
valid 7TM alignment, the longest and tightest wins.

**State.** Decided in this order:

1. **Partner.** If the entry contains a partner that only binds an active receptor
   (G protein including mini-G, Nb35/scFv16, arrestin, GRK, or an active-state nanobody
   such as Nb80/Nb6B9), the state is **Active**. On the benchmark this rule never
   fired on an inactive structure (0 of 460).
2. **Structure.** Otherwise the cytoplasmic end of TM6 is measured against matched
   active/inactive reference pairs for the receptor's class:
   - **Class A + odorant:** a gradient-boosted model combines the readings against
     three class A reference pairs (rhodopsin, Gi-coupled, Gs-coupled) with quality
     measures. When its probability of Active is between 0.40 and 0.60 the call is
     **Intermediate**, meaning "too close to call", not a conformational label.
     The model is trained on the 535 benchmark class A + odorant entries that have
     no partner (157 Active, 378 Inactive, 129 receptors), the same population it
     is applied to.
   - **Class C:** read from the dimer interface (TM6–TM6 vs TM4–TM4 distance between
     protomers), Active below −36 Å.
   - **Class F:** a fitted cut on the TM6 activation axis.
   - **Classes B1, B2, D1:** Active past the midpoint between the inactive and active
     reference (too few structures without a partner to fit a cut).
3. **Intermediate, low confidence** when the state cannot be measured (the output
   carries a warning). That happens for:
   - taste receptors (T2) without a partner, since no inactive reference exists;
   - GPR156 without a partner, since no inactive GPR156 structure exists;
   - class C entries with a single protomer;
   - the few structures that cannot be aligned to any reference.

`confidence` is `low` when the structural measurement failed its quality check.

### The TM6 measurement

Each class has a reference pair: the same receptor solved inactive and active.

| Class | Inactive | Active | Receptor |
|---|---|---|---|
| A | 1U19 | 3PQR | rhodopsin |
| A | 4DKL | 6DDE | μ-opioid (Gi-coupled) |
| A | 2VT4 | 7JJO | β1-adrenergic (Gs-coupled) |
| B1 | 5VEW | 6X18 | GLP-1R |
| B2 | 8IKJ | 8IKL | CD97 |
| C | 7MTQ | 7MTS | mGlu2 |
| D1 | 7QA8 | 7AD3 | Ste2 |
| F | 5L7D | 6OT0 | Smoothened |
| F | 6BD4 | 8WMA | Frizzled-4 |

Odorant receptors use the class A pairs. B1 and B2 entries are read against both
B pairs, because step 1 sometimes swaps the two.

For each pair, no generic residue numbering needed:

1. Both members are trimmed to the 7TM domain (by alignment to rhodopsin), dropping
   fusion partners such as T4L or BRIL.
2. They are superposed on their **static core**, found iteratively (pairs within
   2 Å), so the part that moves does not bias the fit.
3. Residues still more than 4 Å apart form the **mobile set**. In practice this is
   the cytoplasmic end of TM6.
4. The centres of the mobile set in the inactive and active member define the
   **activation axis**.

The query chain is aligned to the inactive member and superposed on core residues
only. Its mobile residues are located by walking along the chain from nearby core
residues, because the alignment itself is unreliable inside the moving helix. The
centre of those residues is projected onto the axis: **activation** 0 means where the
inactive reference sits, 1 where the active one does.

A reading passes the quality check (`confidence` high) when the query's TM6 lies close
to the axis, at least 8 mobile residues could be placed, and the core fits within
3.5 Å. For class A all three pairs are read; every reading goes into the model.

**Class C** is read differently, because TM6 barely moves within one protomer;
activation rotates the dimer from a TM4–TM4 to a TM6–TM6 interface. The tool measures
the TM6–TM6 minus TM4–TM4 distance between the two protomers, at the positions
equivalent to rhodopsin W4.50 and P6.50, and scales it between 7MTQ (inactive) and
7MTS (active).

## Accuracy

Benchmark: 1,692 GPCRdb entries with an Active/Inactive label (labels as of
2026-11-19). The figures are from cross-validation that holds out whole receptors,
so no receptor is ever used both to fit and to score. They are the share of
Active/Inactive calls that are correct.

"Borderline" is the model's Intermediate band; "Not measurable" entries are also
reported as Intermediate (low confidence). Neither counts towards accuracy.

| Class | Entries | Accuracy | Borderline | Not measurable |
|---|---|---|---|---|
| A + odorant | 1,353 | 98.1% | 5 | 6 |
| B1 | 147 | 96.6% | 0 | 0 |
| B2 | 34 | 100% | 0 | 2 |
| C | 103 | 96.9% | 0 | 7 |
| D1 | 1 | 100% | 0 | 0 |
| F | 33 | 87.9% | 0 | 0 |
| T2 | 18 | 100% | 0 | 2 |
| **All** | **1,692** | **97.7%** | **5** | **17** |

- **Class step:** right for 98.0% of entries if class O counts as class A.
- **Chain step:** picked a receptor chain in every one of 1,712 entries.

## Known limitations

- **Class O calls are unreliable.** Step 1 labels many class A receptors as O1/O2
  (356 calls on the benchmark, 8 true odorant receptors). The output flags this with a
  warning. It does not affect the state call, because odorant receptors are measured
  as class A.
- **Fusion constructs.** Class A structures with a crystallisation fusion (T4L/BRIL,
  usually in ICL3) cause 24 of the 26 class A errors: 7.6% error rate vs 1.0% without.
- **CaSR with agonist and inhibitor both bound** (e.g. 7M3E, 7DD5) has the active dimer
  geometry and is called Active, although it is labelled Inactive.
## Calibration

The activation cuts and the class A model ship pre-trained in `data/state_model.pkl`,
saved with scikit-learn 1.9.0. With a different scikit-learn version the tool prints a
warning; install the pinned version if that happens.

## Layout

```
annotate_pdb.py                command-line entry point
gpcr_annotator/
  pipeline.py                  runs the three steps
  class_step1.py               step 1, same algorithm as the benchmark script
  chain_step2.py               step 2, same algorithm as the benchmark script
  state_step3.py               step 3 measurement and decision rule
  transducer.py                active-state partner detection (used by step 3) and
                               the experimental-entry check
data/
  ref_class/  ref_chain/  ref_state/   reference structures for each step
  state_model.pkl              trained state calibration (cuts + class A model)
```

## License

MIT; see [LICENSE](LICENSE).
