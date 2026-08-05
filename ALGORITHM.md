# ALGORITHM.md

## 1. Overview

KS-FIM investigates a specific claim: that observability of a resident
space object (RSO), measured via the determinant of a Fisher Information
Matrix (FIM) built from range measurements, is "stronger" when computed in
Kustaanheimo-Stiefel (KS) regularized coordinates than in Cartesian
coordinates. The claim originates from a 2021 COSPAR presentation,
"Regularized Orbit Observability for Resident Space Objects" (PEDAS.1-0023
-21, H. Sellamuthu, `E:\Research\Conferences\c2021 - COSPAR\`), which
reports raw 4×4 KS-space FIM determinants many orders of magnitude larger
than the corresponding 3×3 Cartesian determinants for the same tracking
scenario. This repo builds both FIMs independently, using KSROP's own
tested `car2ks`/`ks2car` transformation as ground truth, to check whether
that comparison is methodologically sound.

Independent of every other repo under `GitHub\` except for its one real
dependency, KSROP (reused for the KS transformation and orbital-element
conversion — same reuse pattern as KS-Pc and OREM).

## 2. Problem Statement

A Fisher Information Matrix for `n` range measurements of a physical
position `r` is `F = Σᵢ (1/σᵢ²) (∂hᵢ/∂r)(∂hᵢ/∂r)ᵀ`, with `hᵢ(r) = |r - Rᵢ|`
the range to station `i`. In Cartesian coordinates, `r` is a genuine
3-dimensional state and `F` is 3×3.

The KS transformation maps the same physical position to a 4-dimensional
state `u = (u₁,u₂,u₃,u₄)` via `r = (x,y,z,0) = L(u)u`, where `L(u)` is the
4×4 Levi-Civita/KS matrix. This is the well-known **Hopf-fibration
redundancy** of KS space: `r` depends on only 3 true degrees of freedom,
and `u` carries one extra "gauge" (phase) direction that leaves `r`
unchanged. The presentation's own methodology works around this
redundancy at the *transformation* level (a `u₃=0` or `u₄=0` gauge-fixing
condition, per its "Alternate KS transformation" slide) but does not
address it at the *FIM comparison* level: it builds a genuine 4×4 `F` in
`u`-space and compares its raw determinant directly against the 3×3
Cartesian determinant for the same scenario.

**The question this repo answers**: is that comparison valid, or does it
conflate a coordinate-representation artifact (the redundant 4th
dimension, and/or unnormalized magnitude scaling between differently-
dimensioned/differently-scaled coordinate systems) with genuine gain in
observability?

## 3. Inputs

- Satellite osculating Keplerian elements `(a, e, i, RAAN, AOP, M)` — km,
  dimensionless, degrees — converted to Cartesian state via KSROP's
  `oe2car`.
- Observing-station positions, either directly in the same inertial frame
  as the satellite state (test suite) or geodetic lat/lon/altitude on a
  spherical Earth (`app/ksfim_case_study.F`, transcribed directly from the
  COSPAR presentation's own station table — Kwajalein, Ascension Island,
  SHAR, Tenerife).
- Per-station range-measurement standard deviation `σᵢ` (km).

## 4. Core Algorithm

1. **State conversion** — `oe2car` (KSROP) for the satellite's Cartesian
   state; `car2ks` (KSROP) for its KS state `u`, at the oscillator
   frequency `w = sqrt(-0.5·μ/a)`.
2. **Cartesian FIM** — `fim_build_cartesian` (`src/fim.F`): for each
   station, the range gradient is the exact analytical unit line-of-sight
   vector `(r - Rᵢ)/|r - Rᵢ|`. `F = Σᵢ (1/σᵢ²) grad·gradᵀ`, 3×3.
3. **KS-space FIM** — `fim_build_ks` (`src/fim.F`): for each station, the
   range gradient w.r.t. `u` is computed by **central finite difference**
   through KSROP's own `ks2car` forward map — `h(u) = |r(u) - Rᵢ|`,
   perturbing each of the 4 components of `u` independently. Deliberately
   numerical rather than the presentation's analytical `∂u/∂r = L(u)⁻¹/2`
   chain rule, so this doesn't silently inherit any error in that
   derivation — it re-derives the same gradient independently from
   already-tested transformation code. `F = Σᵢ (1/σᵢ²) grad·gradᵀ`, 4×4.
4. **Eigendecomposition** — `jacobi_eig` (`src/linalg.F`, classic cyclic
   Jacobi rotation method, since `F` is always real symmetric
   positive-semidefinite by construction). Used instead of a raw
   determinant because a determinant can be a deceptively small-but-
   nonzero float even when an eigenvalue is genuinely zero to numerical
   precision — the actual question ("is this matrix rank-deficient") is
   an eigenvalue question, not a determinant question.
5. **Rank-corrected comparison** — `fim_reduced_det_ks` /
   `fim_reduced_det_cartesian` (`src/fim.F`): count eigenvalues above a
   relative threshold (`10⁻⁸` × largest) to get an empirical rank, then
   take the product of only those eigenvalues (not always top-3/top-3 —
   with fewer than 3 independent range directions, e.g. 2 stations, even
   the Cartesian FIM is rank-deficient, and comparing at mismatched ranks
   produces meaningless results, confirmed as a real failure mode during
   this repo's own test development — see `test/test_fim_reduced_
   consistency.F`'s header).

## 5. Key Equations / Physics

**Why the 4×4 KS-space FIM should be rank ≤ 3, not 4**: since
`grad_ks = Jᵀ · grad_cartesian` where `J = ∂r/∂u` (3×4, because `r` has
only 3 independent components), the full KS-space FIM is
`F_ks = Jᵀ F_cartesian J`. A matrix formed this way can never have rank
exceeding `rank(F_cartesian) ≤ 3`, regardless of how many stations
contribute or their geometry — this is a structural property of the
transformation, not an empirical coincidence.

**Why the "advantage" should be a predictable scale factor, not new
information**: the KS matrix satisfies the standard identity
`L(u)ᵀL(u) = |u|² I₄` — the KS map is *conformal* (a uniformly-scaled
orthogonal transformation) at any fixed `u`. Consequently, for a fixed
satellite position, the ratio `(reduced det F_ks) / (reduced det
F_cartesian)`, when both are reduced to their *matching* true rank, should
be a constant depending only on the satellite's position (via `|u|²=r`),
independent of which stations are used or how many — **not** a genuine
improvement in how well the position is determined.

## 6. Outputs

`app/ksfim_case_study.F` prints, per case-study object: the 4 KS
eigenvalues (descending), the empirical rank, `det(F_cartesian)`, the raw
4×4 `det(F_ks)` (the quantity the original presentation reported), the
rank-corrected reduced `det(F_ks)`, and both the raw and corrected
KS/Cartesian ratios — so the inflated raw number and the corrected one are
visible side by side rather than the corrected analysis replacing the
original silently.

## 7. Complexity & Performance

Trivial — small dense linear algebra (3×3/4×4 matrices), no propagation,
no optimization loop. Single-threaded; the 4-core cap (`GitHub\CLAUDE.md`
§1) is not a binding constraint here.

## 8. Validation & Accuracy

3/3 (`test_fim_cartesian`, hand-computable orthogonal-LOS geometry),
24/24 (`test_fim_ks_rank`, 3 orbital regimes × 4 station counts each,
checking both rank and eigenvalue-gap magnitude), 1/1
(`test_fim_reduced_consistency`) — **28/28 passing**.

**Empirical findings** (this repo's own runs, not carried over from the
presentation):
- `test_fim_ks_rank`: across LEO/GTO/high-eccentricity test cases and 2-5
  station counts, `rank(F_ks)` is *exactly* `min(nstat, 3)` in every case
  — never 4. Where a 4th eigenvalue exists numerically, it is 8-18 orders
  of magnitude smaller than the largest (e.g. `3.7×10⁸` vs `-7.1×10⁻⁹`) —
  consistent with floating-point noise around a true zero, not a genuine
  information direction.
- `test_fim_reduced_consistency`: for one fixed satellite position, the
  rank-corrected ratio `reduced_det_ks / reduced_det_cartesian` across 5
  station configurations spanning 3-4 stations, very different geometries
  and σ values, was **`3.7257368×10¹⁴` in every case, agreeing to 7-8
  significant figures** (max relative deviation `5.2×10⁻⁸`). This is a
  quantitative confirmation that the KS/Cartesian ratio is a fixed
  function of satellite position, not of the tracking geometry — exactly
  what a coordinate-scaling artifact predicts, and not what a genuine
  observability improvement from station geometry would look like.
- `ksfim_case_study`: reproducing the presentation's own 4 HEO objects
  (35497, 37151, 39615, 42928) with its own published orbital elements and
  station table (simplified: no GMST epoch alignment — see the app's
  header comment for why this doesn't affect the structural question),
  `rank_est = 3` in all 4 cases, with the 4th eigenvalue `10⁻⁸` to
  `10⁻¹⁰` relative to the largest. The rank-corrected KS/Cartesian ratio
  differs *between* objects (`9.6×10¹³` to `2.4×10¹³`, roughly a 4×
  spread) — consistent with the ratio depending on each object's own
  position/`|u|`, exactly as predicted, rather than being universal.

**Conclusion**: the original presentation's raw-4×4-determinant comparison
is not methodologically sound. The 4×4 KS-space FIM is empirically
confirmed rank-deficient (rank 3, structurally, in every tested case), and
once correctly rank-matched against the Cartesian FIM, the resulting ratio
behaves exactly as a coordinate/Jacobian scaling artifact would — constant
across station geometry, varying only with satellite position. This
doesn't necessarily mean representation-dependent observability is never a
real effect for nonlinear systems in general (a legitimate research
question), but this specific presentation's evidence for a genuine "KS
gives better observability" effect does not survive a rank-correct,
same-basis comparison.

## 9. Known Limitations

- **No GMST/epoch-accurate station placement.** Stations are placed at
  their published geodetic coordinates in the same inertial frame as the
  satellite, without rotating for Earth's sidereal orientation at the
  actual observation epoch. This was a deliberate Phase 1 simplification
  (see `app/ksfim_case_study.F` header) — it does not affect the
  structural rank-deficiency/Jacobian-scaling findings above, which are
  properties of the coordinate transformation itself, not of tracking
  geometry — but it does mean the absolute magnitudes in this repo's case
  study are not a bit-for-bit reproduction of the presentation's own
  numbers.
- **No time-series / multi-revolution observability tracking.** The
  presentation's own plots show observability evolving over ~2 orbital
  revolutions; this repo checks single-instant snapshots only, sufficient
  for the structural question under test but not a full reproduction of
  the original study's scope.
- **Numerical (not analytical) KS-space gradient.** `range_grad_ks` uses
  central finite differences rather than an analytical Jacobian. This was
  a deliberate choice (independence from the presentation's own possibly-
  flawed analytical derivation) but means gradient accuracy is bounded by
  finite-difference truncation/roundoff (~`10⁻⁶` relative step), not
  exact to machine precision.
- **Not yet extended to the presentation's GEO case (object 28868).** Only
  the 4 HEO objects are reproduced; the GEO single/two-station case from
  the same presentation is not yet checked.
- **Doesn't (yet) address whether representation-dependent local
  linearization is a real effect elsewhere.** This repo falsifies one
  specific presentation's specific comparison; it does not itself explore
  whether there are other, correctly-normalized senses in which KS
  coordinates could offer a genuine local-conditioning advantage for
  nonlinear filtering (a legitimate, different research question, out of
  scope for Phase 1).

## 10. Dependencies

- **KSROP** (`git+https://github.com/hari251086/KSROP`, tag `v2.2.0`) —
  `car2ks`/`ks2car` (the KS transformation itself, `src/Subrouts.F`),
  `oe2car` (Keplerian-to-Cartesian conversion), `init_constants`
  (physical constants incl. `R_Earth`, `amue`, read from
  `input/const_new.dat`, copied from KSROP's own `input/` at repo setup
  time — same pattern as KS-Pc/OREM).
- **Source material**: "Regularized Orbit Observability for Resident
  Space Objects," H. Sellamuthu, 43rd COSPAR Scientific Assembly,
  PEDAS.1-0023-21 (2021) — the presentation this repo evaluates. Located
  at `E:\Research\Conferences\c2021 - COSPAR\Paper and Presentation\
  DROPBOX-31122020\` (not part of this repo; a local research-archive
  file, referenced for orbital elements/station data only).
