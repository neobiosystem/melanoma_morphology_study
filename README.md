# Closed-loop melanoma morphology study — reproduction guide

A four-step closed loop that links clinical dermoscopy images to a mechanistic
cell-dynamics model and back again:

```
  clinical photos ──> morphology            (Part 1, complete)
                          │
                          v
               cell-dynamics CA            (Part 2, complete — negative result)
                          │
                          v
          mechanism ↔ morphology alignment (Part 3, complete)
                          │
                          v
              emergent mechanism            (Part 4, planned)
```

Parts 1–3 are complete. This folder is a **self-contained mirror** of the source
tree with exactly two things removed: the 21 GB raw image/mask archive (replaced
by download links below) and large derived graphics you can regenerate. It is
organized so the three parts keep their sibling layout — Part 2 and Part 3
`addpath` Part 1/Part 2 as siblings, so the relative paths below must be
preserved for the pipeline to run.

## Repository layout

```
for_git/
├── 黑色素瘤研究计划.md         # the four-step study plan (Chinese)
├── Part1/                     # clinical image → morphology
│   ├── 01-clinical-image-extraction-en.md   # Part 1 spec
│   ├── config/default_config.m              # every Part 1 path/threshold/hyperparameter
│   ├── +morph/ +seg/ +stats/ +agree/ +io/   # measurement & statistics packages
│   ├── scripts/  tests/
│   ├── data/                                # manifest.csv + isic_metadata.csv (NO raw)
│   └── results/                             # expected outputs (see below)
├── Part2/                     # cellular automaton
│   ├── 02-cell-dynamics-simulation-en.md
│   ├── config/sim_config.m
│   ├── +ca/ +bridge/ +comp/ +phen/
│   ├── scripts/  tests/
│   └── results/
└── Part3/                     # mechanism → morphology alignment
    ├── 03-simulation-image-alignment-en.md
    ├── config/part3_config.m
    ├── +ca3/ +mech/ +align/
    ├── scripts/  tests/
    └── results/
```

## Requirements

- **MATLAB R2023a** (Image Processing Toolbox required; Statistics and Machine
  Learning Toolbox for the tests/analyses that call `runtests`).
- A parallel pool for the sweeps. This machine is limited to **8 GB RAM**, so
  the configs use **4 workers** (`sim_config.m`), not 6 — raise it only if your
  machine has more memory.
- All paths resolve from `fileparts(mfilename('fullpath'))`, and every entry
  script below adds its own paths, so you can run them from any working
  directory as long as the `Part1/Part2/Part3` sibling structure is intact.

## Getting the data (Part 1 only)

Part 1 measures real lesions, so it needs the ISIC 2018 challenge archive. The
URLs are the single source of truth in `Part1/config/default_config.m`:

| file | URL |
|---|---|
| images | `https://isic-challenge-data.s3.amazonaws.com/2018/ISIC2018_Task1-2_Training_Input.zip` |
| masks | `https://isic-challenge-data.s3.amazonaws.com/2018/ISIC2018_Task1_Training_GroundTruth.zip` |

Unzip them so the raw folders land exactly where the config expects:

```
Part1/data/raw/ISIC2018_Task1-2_Training_Input/     # *.jpg images
Part1/data/raw/ISIC2018_Task1_Training_GroundTruth/ # *_segmentation.png masks
```

The three-tier grade (`benign` / `in_situ` / `invasive`) is not in the challenge
archive; it is joined back from the ISIC Archive API (`collectionId = 63`) into
`Part1/data/isic_metadata.csv`, which is already included here. To re-fetch it:
`run Part1/scripts/fetch_cohort_metadata.m`. Part 2 and Part 3 are pure
simulation and run **without** the raw archive.

## Reproduction order

**Part 1 — morphology.** Measure the cohort on both segmentation routes
(expert ground-truth masks and automatic Otsu), then derive the between-grade
statistics and figures:

```matlab
run('Part1/scripts/run_analysis.m');      % writes descriptors.csv + qc_report.csv
run('Part1/scripts/export_posthoc.m');    % writes posthoc_all.csv
run('Part1/scripts/make_figures.m');      % fig1–fig6
run('Part1/scripts/robustness_dual_route.m');   % route sign-reversal figure + table
```

`run_analysis` accepts `'Limit', N, 'Routes', {'expert'|'otsu'}` for a fast
smoke pass before the full cohort.

**Part 2 — cell dynamics (negative result).** Sweep four candidate mechanisms
and score their morphology through the Part 1 bridge (no images needed):

```matlab
mechanism_comparison(1024, false);        % full 4-arm sweep → mechanism_trajectories.csv + part2_handoff.json
```

For chunked/background runs, use the per-arm form
`mechanism_comparison(1024, false, 'nutrient_fingering')` and then
`run('Part2/scripts/assemble_results.m')`.

**Part 3 — alignment.** Phase A compares the two new mechanisms (push
displacement vs. growth-field heterogeneity); Phase B sweeps and sensitivity-analyses
the winner:

```matlab
run('Part3/scripts/run_discovery.m');        % Phase A → discovery_trajectories.csv + part3_handoff.json
run('Part3/scripts/run_phaseB.m');           % Phase B → phase_*.csv + sobol_*.json (requires oat_*.csv first)
```

Each part has a `scripts/report.m` that renders the paper-style prose report
(R2023a has no plain-text Live Code, so use `publish`):

```matlab
publish('Part3/scripts/report.m', struct('format','html','outputDir','Part3/results/report'))
```

## Verification

The test suites use synthetic shapes only (no raw data), so they run anywhere:

```matlab
run('Part1/scripts/run_tests.m');   % box-counting vs. filled disk, descriptor checks
run('Part2/scripts/run_tests.m');   % lattice isotropy, measurement-bridge regression
run('Part3/scripts/run_tests.m');   % homogeneous-limit anisotropy, distance hand-check
```

All suites must be green before trusting any result. The load-bearing checks are
the box-counting/disk test (Part 1), the isotropy test (Part 2), and the
homogeneous-limit + distance hand-computation tests (Part 3).

## What the included `results/` folders are for

Each `PartN/results/` holds the **expected outputs** — the exact CSV/JSON/handoff
tables and figures this pipeline produced — so you can (a) publish the `report.m`
prose immediately, and (b) diff your own reruns against the reported numbers.
They are small (~3 MB total). Not included (regenerable): the per-lesion QC
montages (`make_qc_montages.m`), the rendered share page (`build_share_page.py`),
and the curated exemplar images (`make_classic_examples.m`).
