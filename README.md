# KS-FIM

Rank-correct Fisher Information Matrix observability comparison, Cartesian vs. Kustaanheimo-Stiefel coordinates

Copyright (c) 2026 Harishkumar Sellamuthu — MIT License (see `LICENSE`)

## 1. Overview

KS-FIM checks a specific claim from a 2021 COSPAR presentation
("Regularized Orbit Observability for Resident Space Objects,"
PEDAS.1-0023-21): that satellite-tracking observability, measured as the
determinant of a Fisher Information Matrix (FIM) built from range
measurements, is "stronger" when computed in Kustaanheimo-Stiefel (KS)
regularized coordinates than in Cartesian coordinates.

The claim is methodologically suspect on its face: KS coordinates are
4-dimensional (`u₁..u₄`) while physical position is only 3-dimensional —
the extra dimension is a well-known redundant "gauge" direction (the KS
map's Hopf-fibration structure). A 4×4 FIM built from range measurements,
which only ever depend on the true 3-DOF position, should therefore be
structurally rank-deficient — and comparing its raw determinant directly
against a 3×3 Cartesian determinant, without correcting for that
dimensional mismatch, is not a like-for-like comparison. This repo builds
both FIMs independently (using KSROP's own tested KS-transformation code
as ground truth) and checks empirically whether that's actually what's
going on.

**Status: Phase 1 complete.** Result: the rank-deficiency hypothesis is
confirmed in every tested case (28/28 tests passing, including a direct
reproduction of the presentation's own 4 case-study objects), and the
rank-corrected KS/Cartesian ratio behaves exactly as a coordinate-scaling
artifact would — constant across very different station geometries for a
fixed satellite position — not as evidence of genuine new information from
using KS coordinates. See `ALGORITHM.md` §8 for full findings, or issue #6
for the same writeup as the durable GitHub record of the finding.

Independent of every other repo under `GitHub\` except KSROP (reused for
the KS transformation, same pattern as KS-Pc/OREM).

## 2. Project Structure

```
KS-FIM/
├── src/
│   ├── linalg.F       Jacobi eigenvalue solver + small determinant
│   │                    helpers for real symmetric 3x3/4x4 matrices
│   └── fim.F           FIM builders (Cartesian analytical gradient,
│                        KS finite-difference gradient), rank-corrected
│                        reduced-determinant helpers
├── app/
│   └── ksfim_case_study.F   reproduces the COSPAR presentation's own
│                              4 HEO case-study objects
├── test/
│   ├── test_fim_cartesian.F            hand-computable sanity check
│   ├── test_fim_ks_rank.F              core rank-deficiency test
│   └── test_fim_reduced_consistency.F  Jacobian-scaling-artifact test
├── input/
│   └── const_new.dat   physical constants (from KSROP)
├── fpm.toml
└── ALGORITHM.md         full methodology + findings writeup
```

## 3. Quick Start

```bash
fpm build --compiler ifx
fpm test --compiler ifx
fpm run ksfim_case_study --compiler ifx
```

Expect 28/28 tests passing and, for each of the 4 case-study objects, a
printed comparison of the raw (presentation-style) and rank-corrected
KS/Cartesian determinant ratios.

## 4. Building / Setup

Fortran, `fpm`-based, fixed-form source with implicit typing (matches
KSROP/OREM/KS-Pc convention). Verified with Intel `ifx` 2025.0 on Windows
via the Visual Studio + Intel oneAPI environment (`vcvars64.bat` +
Intel Fortran `vars.bat`, both required in the same shell session before
invoking `fpm`). Not yet verified with `gfortran`.

## 5. Running

- `fpm run ksfim_case_study --compiler ifx` — the only executable; prints
  the 4-object case-study comparison described in Quick Start.
- No CLI arguments or config files — all inputs (orbital elements, station
  coordinates) are hardcoded from the source presentation's own published
  tables, directly in `app/ksfim_case_study.F`.

## 6. Testing

`fpm test --compiler ifx` — **28/28 tests passing** as of the last run
documented here (2026-08-05):
- `test_fim_cartesian`: 3/3 — Cartesian FIM builder matches a
  hand-computable orthogonal-line-of-sight geometry exactly.
- `test_fim_ks_rank`: 24/24 — `rank(F_ks) = min(nstat, 3)` across 3
  orbital regimes (LEO/GTO/high-e) × 4 station counts (2-5), never rank 4.
- `test_fim_reduced_consistency`: 1/1 — the rank-corrected KS/Cartesian
  ratio is constant to 7-8 significant figures across 5 very different
  station configurations at a fixed satellite position.

## 7. Inputs & Outputs

**Inputs** (see `ALGORITHM.md` §3): Keplerian elements, station geodetic
coordinates, per-station range-measurement σ — all hardcoded per case
study, transcribed directly from the source presentation's own tables (not
re-derived).

**Outputs**: console output only (no file I/O) — per-object KS
eigenvalues, empirical rank, and both raw and rank-corrected det(F)
comparisons. `output/` exists per repo convention but is currently unused.

## 8. Known Issues / Limitations

See `ALGORITHM.md` §9 for the full technical detail. Tracked as issues:
#1 GEO case study (object 28868) not yet reproduced, #2 multi-revolution
observability time series not yet implemented, #3 GMST-accurate station
placement (deliberate Phase 1 simplification — doesn't affect the
structural findings), #4 numerical (not yet cross-checked analytical)
KS-space gradient, #5 open question on whether any correctly-normalized
representation-dependent observability effect exists beyond what this
repo falsified.

## 9. Version History

- **2026-08-05** — Initial implementation. `src/linalg.F` (Jacobi
  eigenvalue solver), `src/fim.F` (Cartesian + KS FIM builders,
  rank-corrected reduced-determinant helpers), full test suite (28/28
  passing), and `app/ksfim_case_study.F` reproducing the source
  presentation's 4 HEO case-study objects. Core finding: the presentation's
  raw 4×4-vs-3×3 determinant comparison is not methodologically sound —
  the KS-space FIM is structurally rank-deficient (confirmed empirically
  in every tested case), and the rank-corrected KS/Cartesian ratio behaves
  as a predictable coordinate-scaling artifact, not genuine new
  observability information. See `ALGORITHM.md` §8.

## 10. Dependencies / References

- **KSROP** (`git+https://github.com/hari251086/KSROP`, tag `v2.2.0`) —
  KS transformation (`car2ks`/`ks2car`), Keplerian-Cartesian conversion
  (`oe2car`), physical constants.
- **Source material under evaluation**: "Regularized Orbit Observability
  for Resident Space Objects," H. Sellamuthu, 43rd COSPAR Scientific
  Assembly, PEDAS.1-0023-21 (2021). Local copy at `E:\Research\
  Conferences\c2021 - COSPAR\Paper and Presentation\DROPBOX-31122020\`
  (not part of this repo).
- Related, same-author KS-regularization lineage cited by the source
  presentation: Sellamuthu & Sharma, "Orbit theory with lunar perturbation
  in terms of Kustaanheimo–Stiefel regular elements," J. Guid. Control
  Dyn. 40 (2017) 1272-1277; "Hybrid orbit propagator for small spacecraft
  using Kustaanheimo–Stiefel elements," J. Spacecraft Rockets 55 (2018)
  1282-1288. See also this workspace's `reference_wider_reentry_materials
  _survey` memory for the fuller personal-research lineage this work sits
  in.
