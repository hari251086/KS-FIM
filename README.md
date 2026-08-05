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

**Status: Phase 1-4 complete (65/65 tests).** The rank-deficiency
hypothesis is confirmed in every tested case, including a direct
reproduction of the presentation's own 4 HEO case-study objects plus its
GEO case, and the rank-corrected KS/Cartesian ratio behaves exactly as a
coordinate-scaling artifact would — constant across very different
station geometries for a fixed satellite position, not evidence of
genuine new information from using KS coordinates. Phase 2 (issues #3,
#4) added an independent analytical-gradient cross-check (agrees with the
numerical one to full displayed precision) and GMST-accurate per-epoch
station placement — which left the rank-corrected ratio bit-for-bit
unchanged while the *raw* ratio changed substantially, even flipping sign
for one object. Phase 3 (issue #1) added the GEO case and an extreme
1-station edge test — confirming `rank(F_ks)=rank(F_cartesian)` down to
rank 1, and generalizing a caveat the source presentation only noted for
near-degenerate station pairs (raw determinant → 0) into "true at any low
station count." Phase 4 (issue #2) derived and verified an exact closed
form, `reduced_det_ks/det_cartesian = 64·r³`, and confirmed it holds to
`~1e-9` relative precision at every point across 2 full orbital
revolutions of a real case-study object — the ratio isn't just constant
across station geometry, its exact value is predictable from the
satellite's instantaneous radius alone. See `ALGORITHM.md` §8 for full
findings, or issue #6 for the same writeup as the durable GitHub record of
the finding.

Independent of every other repo under `GitHub\` except KSROP (reused for
the KS transformation, same pattern as KS-Pc/OREM).

## 2. Project Structure

```
KS-FIM/
├── src/
│   ├── linalg.F       Jacobi eigenvalue solver + small determinant
│   │                    helpers for real symmetric 3x3/4x4 matrices
│   ├── fim.F           FIM builders (Cartesian analytical gradient, KS
│   │                    finite-difference AND analytical gradient),
│   │                    rank-corrected reduced-determinant helpers
│   └── timeconv.F       GMST (IAU-1982) + geodetic-to-ECI conversion
├── app/
│   ├── ksfim_case_study.F   reproduces the COSPAR presentation's own
│   │                          4 HEO case-study objects (GMST-accurate)
│   │                          plus its GEO case (issue #1)
│   └── ksfim_timeseries.F   multi-revolution closed-form-ratio check,
│                              144 samples over 2 revs of object 35497
│                              (issue #2)
├── test/
│   ├── test_fim_cartesian.F            hand-computable sanity check
│   ├── test_fim_ks_rank.F              core rank-deficiency test
│   ├── test_fim_reduced_consistency.F  Jacobian-scaling-artifact test
│   ├── test_fim_analytical_vs_fd.F     analytical-vs-numerical gradient
│   │                                    cross-check (issue #4)
│   ├── test_timeconv.F                 GMST/geodetic-to-ECI checks
│   │                                    (issue #3)
│   ├── test_fim_geo_edge.F             1-/2-station rank-tracking edge
│   │                                    case (issue #1)
│   └── test_fim_ratio_formula.F        exact ratio=64*r^3 closed-form
│                                        check (issue #2)
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
fpm run ksfim_timeseries --compiler ifx
```

Expect 65/65 tests passing and, for each of the 4 HEO case-study objects
plus the GEO case, a printed comparison of the raw (presentation-style)
and rank-corrected KS/Cartesian determinant ratios (both finite-difference
and analytical). `ksfim_timeseries` prints the max relative deviation from
the exact closed-form ratio `64·r³` across 2 full revolutions (expect
`~1e-9`).

## 4. Building / Setup

Fortran, `fpm`-based, fixed-form source with implicit typing (matches
KSROP/OREM/KS-Pc convention). Verified with Intel `ifx` 2025.0 on Windows
via the Visual Studio + Intel oneAPI environment (`vcvars64.bat` +
Intel Fortran `vars.bat`, both required in the same shell session before
invoking `fpm`). Not yet verified with `gfortran`.

## 5. Running

- `fpm run ksfim_case_study --compiler ifx` — prints the 4-object
  case-study comparison described in Quick Start.
- `fpm run ksfim_timeseries --compiler ifx` — prints and writes
  `output/ksfim_timeseries_35497.csv`, the 144-sample multi-revolution
  closed-form-ratio check described in Quick Start.
- No CLI arguments or config files — all inputs (orbital elements, station
  coordinates) are hardcoded from the source presentation's own published
  tables, directly in `app/ksfim_case_study.F` and `app/ksfim_timeseries.F`.

## 6. Testing

`fpm test --compiler ifx` — **65/65 tests passing** as of the last run
documented here (2026-08-05):
- `test_fim_cartesian`: 3/3 — Cartesian FIM builder matches a
  hand-computable orthogonal-line-of-sight geometry exactly.
- `test_fim_ks_rank`: 24/24 — `rank(F_ks) = min(nstat, 3)` across 3
  orbital regimes (LEO/GTO/high-e) × 4 station counts (2-5), never rank 4.
- `test_fim_reduced_consistency`: 1/1 — the rank-corrected KS/Cartesian
  ratio is constant to 7-8 significant figures across 5 very different
  station configurations at a fixed satellite position.
- `test_fim_analytical_vs_fd`: 9/9 (issue #4) — the analytical KS-space
  gradient agrees with the finite-difference one to ~`1e-10` relative,
  including matching rank and reduced determinant.
- `test_timeconv`: 4/4 (issue #3) — `gmst_deg` matches the published
  J2000.0 reference value to `1e-4` deg; `geodetic_to_eci` sanity checks.
- `test_fim_geo_edge`: 4/4 (issue #1) — `rank(F_ks)=rank(F_cartesian)`
  holds down to the extreme edge case of a single observing station.
- `test_fim_ratio_formula`: 20/20 (issue #2) — the exact closed form
  `reduced_det_ks/det_cartesian = 64·r³` holds to `~1e-8`–`1e-10` relative
  precision across 10 varied orbits/station geometries.

## 7. Inputs & Outputs

**Inputs** (see `ALGORITHM.md` §3): Keplerian elements, station geodetic
coordinates, per-station range-measurement σ — all hardcoded per case
study, transcribed directly from the source presentation's own tables (not
re-derived).

**Outputs**: `ksfim_case_study` is console output only. `ksfim_timeseries`
prints a summary and writes `output/ksfim_timeseries_35497.csv` (per-step
`r`, rank, raw/reduced determinants, actual vs. predicted ratio, reldiff).

## 8. Known Issues / Limitations

See `ALGORITHM.md` §9 for the full technical detail. Tracked as issues:
#5 open question on whether any correctly-normalized representation-dependent
observability effect exists beyond what this repo falsified — the only
issue still open. #1 (GEO case) resolved in Phase 3; #2 (multi-revolution
time series) resolved in Phase 4; #3 (GMST-accurate station placement) and
#4 (analytical-gradient cross-check) resolved in Phase 2.

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
- **2026-08-05 (Phase 2)** — Issues #3 and #4 resolved same day.
  `src/timeconv.F` (`gmst_deg`, `geodetic_to_eci`) adds real per-object
  GMST-accurate station placement to the case study; `range_grad_ks_
  analytical`/`fim_build_ks_analytical` (`src/fim.F`) add an independent
  analytical gradient (textbook `dr=2L(u)du` identity) alongside the
  original finite-difference one. Both strengthen the finding rather than
  just extending it: the analytical gradient agrees with the numerical one
  to full displayed precision (not just a numerical-differencing
  artifact), and the GMST fix left the rank-corrected KS/Cartesian ratio
  bit-for-bit unchanged (depends only on satellite position, exactly as
  predicted) while the *raw* ratio changed substantially under the same
  station-geometry change — even flipping sign for object 42928, which is
  physically nonsensical for a real information metric. 41/41 tests
  passing. See `ALGORITHM.md` §8.
- **2026-08-05 (Phase 3)** — Issue #1 resolved. Added the source
  presentation's GEO case (object 28868, ANIK F1R — a representative
  near-GEO state is used since the source doesn't publish exact orbital
  elements for this case) with its own 1-/2-station sub-cases, plus a
  dedicated regression test (`test_fim_geo_edge`) at the same regime.
  Confirms `rank(F_ks)=rank(F_cartesian)` down to the extreme edge case of
  a single station (rank 1), and generalizes a caveat the source
  presentation only noted for near-degenerate close station pairs (raw
  determinant → 0) into "true at any low station count." Also caught and
  fixed a units mistake found while adding this case: station σ had been
  set to the station *altitude* values instead of the slide's published
  `σ=10 km` — a magnitude-only correction, confirmed not to affect any
  rank or ratio-based finding. 45/45 tests passing. See `ALGORITHM.md` §8.
- **2026-08-05 (Phase 4)** — Issue #2 resolved, plus a bonus exact
  derivation. Derived and verified the closed form
  `reduced_det_ks/det_cartesian = 64·r³` (`test_fim_ratio_formula`, 20/20,
  10 varied orbits/station geometries), then added
  `app/ksfim_timeseries.F` tracking object 35497 across 2 full revolutions
  (144 samples) — the closed-form prediction holds at every sampled point,
  max relative deviation `1.474×10⁻⁹`. This subsumes and strengthens the
  earlier "constant across station geometry" finding (Phase 1) into an
  exact, position-only-dependent formula, and shows the observability
  variation the source presentation's own plots show over a revolution is
  fully explained by `r`'s own variation along the orbit. 65/65 tests
  passing. See `ALGORITHM.md` §5, §8.

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
