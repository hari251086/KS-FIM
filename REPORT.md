# KS-FIM Investigation Report: Methodology, Reasoning, and Analysis

**Repo:** `hari251086/KS-FIM` | **Author:** Harishkumar Sellamuthu
**Investigation period:** 2026-08-05 (single session, Phases 1-6)
**Status:** Core question closed and falsified (exact, reproducible).
Extension question (issue #5) left open, inconclusive.

---

## 1. Executive Summary

A 2021 COSPAR presentation ("Regularized Orbit Observability for Resident
Space Objects," PEDAS.1-0023-21) claimed that satellite-tracking
observability — measured as the determinant of a Fisher Information
Matrix (FIM) built from range measurements — is quantitatively "stronger"
when the FIM is computed in Kustaanheimo-Stiefel (KS) regularized
coordinates than in Cartesian coordinates, for the same physical tracking
scenario.

This investigation built an independent, from-scratch implementation of
both FIMs against KSROP's own tested KS-transformation code, and tested
the claim empirically and analytically across six escalating phases. The
result:

- **The claim, as the source measured it (raw 4×4 vs. 3×3 determinant
  comparison), is false.** The 4×4 KS-space FIM is structurally
  rank-deficient (rank ≤ 3) in every tested case — the "advantage" the
  presentation reports is not new physical information, it is an
  unnormalized units/dimension artifact.
- **The artifact has an exact closed form.** Once correctly rank-matched,
  the KS/Cartesian information ratio equals exactly `64·r³` (`r` = the
  satellite's instantaneous orbital radius), independent of station
  geometry, count, or measurement noise. This was verified to
  10⁻⁸–10⁻¹¹ relative precision across dozens of orbits and a full
  144-point, 2-revolution time series.
- **The artifact extends to conditioning, not just magnitude.** The FIM's
  condition number — a dimensionless, rank-matched metric relevant to
  numerical filter behavior — is *also* exactly representation-invariant.
- **A deeper, genuinely different question remains open.** Whether actual
  *sequential nonlinear filtering* (an EKF tracking a satellite over
  time, not a single-instant FIM snapshot) shows any representation-
  dependent effect was tested directly with a full covariance Kalman
  filter. The result is **inconclusive** — a numerically jagged
  divergence pattern, most consistent with finite-difference numerical
  noise rather than a real physical effect. This is reported honestly as
  unresolved, not forced into a false conclusion either way.

92/92 automated tests pass. All code, tests, and documentation are
committed to `HS-dev` and `main` (fast-forward merged after every phase).

---

## 2. Origin and Motivating Claim

The source material is a COSPAR conference presentation reporting FIM
observability metrics for four real HEO objects (35497, 37151, 39615,
42928) and one GEO object (28868), tracked by named ground stations
(Kwajalein, Ascension, SHAR, Tenerife, and others for the GEO case). Its
core empirical claim: the raw determinant of a 4×4 KS-space FIM is
routinely 10⁴–10¹⁸ times larger than the corresponding 3×3 Cartesian
FIM's determinant for the identical tracking geometry, interpreted as KS
coordinates providing genuinely superior orbit-determination
observability.

**Two structural red flags motivated building this repo instead of
accepting the claim at face value:**

1. **Dimensional mismatch.** Physical position has exactly 3 true degrees
   of freedom. The KS map `r = (x, y, z, 0) = L(u)u` embeds this 3-DOF
   position into a 4-dimensional `u`-space via the Hopf fibration — the
   4th dimension is a well-known redundant "gauge" direction with no
   independent physical content. A range measurement `h(r) = |r − R_station|`
   depends only on physical position, so a FIM built from such
   measurements, however it's parameterized, can never carry information
   about a direction the measurement doesn't see. A 4×4 matrix built this
   way should therefore be structurally rank-deficient (rank ≤ 3), not
   full rank as raw nonzero 4×4 determinants would imply.
2. **Unnormalized magnitude comparison.** Comparing a raw 4×4 determinant
   (units: km⁸/σ⁸, roughly) against a raw 3×3 determinant (km⁶/σ⁶)
   directly is not a like-for-like comparison — it mixes physical units
   across a coordinate change with its own Jacobian scale factor. A
   magnitude difference under such a comparison is the *expected*
   signature of a units artifact, not evidence of new information.

Both suspicions turned out to be correct and became the spine of Phases
1 and 4 respectively.

---

## 3. Methodology

### 3.1 General approach

The guiding principle throughout: **never take either the source
presentation's numbers or this repo's own derivations on faith — build
independent ground truth and check everything against it, twice where
possible.**

- **Ground truth for the KS transformation**: KSROP's own `car2ks`/
  `ks2car` routines (`src/Subrouts.F`), already independently tested
  across hundreds of cases in that repo, not a re-derivation of the KS
  map. This repo's `fim.F` builds FIMs *on top of* KSROP's transformation
  rather than reimplementing it.
- **Two independent gradient derivations.** The KS-space range gradient
  needed for each FIM was computed two ways that don't share an
  assumption: (a) central finite difference through the real
  `ks2car` map (Phase 1, deliberately *not* the source presentation's own
  analytical `du/dr` chain-rule slide, so this doesn't inherit any error
  in that derivation), and (b) an independently-derived analytical
  gradient from the textbook KS differential identity `dr = 2·L(u)·du`
  (Phase 2). Agreement between the two to ~10⁻¹⁰ relative precision
  (Phase 2) is what let later phases trust the analytical form for
  performance and cleaner derivations without re-litigating Phase 1.
- **Rank via eigenvalues, never via raw determinant.** A cyclic Jacobi
  eigenvalue solver (`src/linalg.F`, `jacobi_eig`) was built from scratch
  (no existing eigensolver anywhere in KSROP) specifically because a raw
  determinant can be a small but nonzero float even when an eigenvalue is
  genuinely zero to numerical precision — the actual rank-deficiency
  question is an eigenvalue question. All comparisons in this repo use a
  "reduced determinant" (product of only the eigenvalues above a
  `1e-8`-relative threshold) rather than the raw determinant of the full
  matrix.
- **Rank-matched, not fixed-size, comparisons.** An early implementation
  bug (caught during Phase 1 test development) hardcoded "always compare
  the top 3 eigenvalues." With fewer than 3 independent range directions
  (e.g. 2 stations), the true rank is 2 on *both* sides, and forcing a
  3rd eigenvalue into the product corrupted the comparison, including
  producing nonsense negative ratios. The fix — always use
  `min(true rank on each side)` — became a standing rule for every
  subsequent phase.
- **Real case-study reproduction before synthetic generalization.** Each
  phase's synthetic test cases (varied orbits, station counts) were
  paired with reproduction of the source presentation's own named objects
  and stations wherever the source published enough data to do so
  (`ksfim_case_study.F`), so findings are anchored to the actual claim
  being tested, not only to convenient synthetic geometry.

### 3.2 Toolchain

- Fortran, `fpm` package manager, fixed-form source with implicit typing
  (`implicit double precision (a-h,o-z)`), matching the KSROP/OREM/KS-Pc
  house style.
- Compiled and tested with Intel `ifx` 2025.0 on Windows (Visual Studio +
  Intel oneAPI environment variables required in the same shell session).
- KSROP consumed as a real `fpm` git dependency (tag `v2.2.0`), not
  hand-copied files.

### 3.3 Escalation structure

Each phase was scoped to *strengthen or extend* the prior phase's finding
rather than simply add more of the same, and each phase's own results
were used to sharpen the next question:

| Phase | Question | Method | Outcome |
|---|---|---|---|
| 1 | Is the 4×4 KS FIM really rank-deficient, and is the ratio a units artifact? | Finite-difference gradient, from-scratch eigensolver | Yes to both — falsified |
| 2 | Is Phase 1 a finite-differencing artifact? Does real station geometry (GMST) change anything? | Independent analytical gradient; real per-epoch station rotation | No and no — strengthened |
| 3 | Does the finding hold at the extreme edge (GEO, 1 station)? | Extended case study + dedicated edge test | Yes, generalizes a caveat the source itself only partially noted |
| 4 | Is the "constant ratio" actually a *predictable* number? | Analytical derivation from KS conformal-map identities | Yes — exact closed form `64r³`, verified analytically and empirically |
| 5 | Is there *any* correctly-normalized comparison where KS wins? | Condition number (dimensionless, rank-matched) | No — also exactly invariant, a direct corollary of Phase 4 |
| 6 | Does *sequential nonlinear filtering* (not a single snapshot) show an effect the FIM can't see? | Full covariance EKF over a real multi-step timeline | **Inconclusive** — first non-clean result, reported as such |

---

## 4. Phase-by-Phase Analysis

### 4.1 Phase 1 — Core Finding: Rank Deficiency and the Units Artifact

**Built**: `src/linalg.F` (Jacobi eigensolver), `src/fim.F`
(`fim_build_cartesian`, `fim_build_ks` via finite-difference gradient,
`fim_reduced_det_ks`/`fim_reduced_det_cartesian`), `app/ksfim_case_study.F`
(reproduces the source's 4 real HEO objects with their published elements
and station table).

**Findings**:
- `rank(F_ks) = min(nstat, 3)` exactly, in every tested case (3 orbital
  regimes × 4 station counts synthetically, plus all 4 real case-study
  objects) — never rank 4. Where a 4th eigenvalue exists numerically it
  is 8–18 orders of magnitude smaller than the largest, consistent with
  floating-point noise around a true zero.
- For a *fixed* satellite position, the rank-corrected ratio
  `reduced_det_ks / reduced_det_cartesian` across 5 wildly different
  station configurations (different counts, geometries, σ) was
  `3.7257368×10¹⁴` in *every* case, agreeing to 7–8 significant figures
  — as clean a confirmation as an empirical test can produce that the
  "advantage" is a predictable, position-dependent Jacobian/conformal-
  scaling artifact (`L(u)ᵀL(u) = |u|²I₄`), not new physical information.
- Across the 4 different real case-study objects, the corrected ratio
  *did* vary (9.6×10¹³ to 2.4×10¹³) — consistent with depending on each
  object's own position, not universal across objects, which is the
  correct signature (constant only for a *fixed* position across
  *different* station geometries, the actual claim under test).

**Verdict at this stage**: not a rejection of KS regularization itself —
its propagator-accuracy benefits are real, well-established, and
unrelated to this finding — but a specific falsification of "KS gives
better FIM-based observability" as the source measured it.

### 4.2 Phase 2 — Independent Analytical Cross-Check + Real Station Geometry

**Built**: `src/timeconv.F` (`gmst_deg`, IAU-1982 formula, verified
against the published J2000.0 reference value to `1e-4` deg;
`geodetic_to_eci`), `range_grad_ks_analytical`/`fim_build_ks_analytical`
(`src/fim.F`, via the textbook `dr = 2·L(u)·du` differential identity).

**Reasoning**: Phase 1's gradient was computed by finite-differencing
`ks2car`. Two questions this doesn't answer on its own: (a) is the result
a numerical-differencing artifact, and (b) does the earlier
fixed-offset station placement (not real per-object GMST rotation) matter?

**Findings**:
- The analytical gradient reproduces the finite-difference reduced
  determinant to full displayed precision (~10⁻¹⁰ relative, identical
  digits) for both synthetic cases and all 4 real objects — Phase 1's
  findings are not a numerical-differencing artifact.
- Adding real per-object epoch/GMST station placement left the
  **rank-corrected** ratio bit-for-bit *unchanged* for every object,
  exactly as predicted (it depends only on satellite position). The
  **raw** ratio, by contrast, changed substantially under the identical
  geometry change — object 42928's raw ratio flipped sign entirely
  (`+1.50×10³ → −2.05×10⁴`). A negative "observability" is physically
  nonsensical for any legitimate information metric — a second,
  independent line of evidence (beyond the magnitude argument) that the
  raw comparison is numerically unsound, not just unnormalized.

41/41 tests passing.

### 4.3 Phase 3 — GEO Case and the Low-Station-Count Edge

**Built**: extended `ksfim_case_study.F` with the source's GEO object
(28868, ANIK F1R — no orbital elements published for this case in the
source, so a representative near-GEO state is used, documented in the
app's header) and its 1-/2-station sub-cases (Kermit; Kermit+Fairbanks);
`test_fim_geo_edge.F`.

**Findings**:
- `rank(F_ks) = rank(F_cartesian) = nstat` holds exactly down to a
  *single* observing station (rank 1) — the raw 4×4 determinant is
  numerically zero (~10⁻⁵⁷ to ~10⁻²⁷) at both 1 and 2 stations.
- This generalizes a caveat the source presentation itself only noted for
  near-degenerate *close station pairs* ("if two stations are very near,
  observability can go to zero, misconstrued as a single station") into a
  structural fact: the raw metric collapses to zero at *any* low station
  count, not just near-coincident geometry.
- A units bug was caught and fixed here: `ksfim_case_study.F`'s station σ
  had been mistakenly set to the station *altitude* column values
  (0.01–2.39 km) instead of the source's actual published `σ = 10 km`,
  present since Phase 1. Confirmed magnitude-only (σ enters both FIMs
  identically per station, so it cannot change rank or the rank-corrected
  ratio) — fixed anyway for numerical fidelity to the published numbers.

45/45 tests passing.

### 4.4 Phase 4 — The Exact Closed Form

**Reasoning**: Phase 1 established the ratio is *constant* across station
geometry for a fixed position. That raises an obvious next question: is
the constant itself *predictable*, or just empirically stable? This phase
derives it analytically first, then verifies.

**Derivation**: `L(u)` satisfies the standard KS identity
`L(u)ᵀL(u) = |u|²I₄` — `L` is a uniformly-scaled orthogonal (conformal)
map at any fixed `u`. Its 3-row physical submatrix `L₃(u)` therefore
satisfies `L₃(u)L₃(u)ᵀ = |u|²I₃`, meaning all 3 singular values of
`L₃(u)` equal `|u|` exactly. Since `F_ks = JᵀF_cartesianJ` with
`J = ∂r/∂u = 2·L₃(u)` (the Phase 2 differential identity), `F_ks`'s
nonzero eigenvalues equal exactly `4|u|² ·` the corresponding
`F_cartesian` eigenvalues. Combined with the KS regularizing identity
`|u|² = r`:

```
reduced_det_ks / det_cartesian = (4|u|²)³ = 64·r³
```

— exactly, independent of station geometry, count, or σ; a single number
determined entirely by the satellite's instantaneous distance from the
primary.

**Built**: `test_fim_ratio_formula.F` (10 varied orbits × station
geometries), `app/ksfim_timeseries.F` (object 35497 tracked across 2 full
revolutions, 144 samples, mean-anomaly-stepped via `oe2car`).

**Findings**:
- `test_fim_ratio_formula`: 20/20 passing, reldiff `~1e-8`–`1e-10` in
  every case, including a direct check of the `|u|² = r` identity itself.
- `ksfim_timeseries`: the closed-form prediction held at *every one* of
  144 sampled points across 2 full revolutions (`r` ranging 6635.62–
  26894.78 km, perigee to apogee) — max relative deviation
  `1.474×10⁻⁹`. This directly explains why the source presentation's own
  plots show observability "evolving" over a revolution: it is fully
  explained by `r`'s own variation along the orbit, not an independent
  KS-space time-dependent effect.

65/65 tests passing.

### 4.5 Phase 5 — Condition Number Invariance

**Reasoning**: Issue #5 (filed during initial scoping) asked whether
*any* correctly-normalized, rank-matched comparison could reveal a real
KS advantage. Phase 4's result already implies an answer for the next
natural dimensionless candidate: since *every* nonzero `F_ks` eigenvalue
scales by the *same* factor (`4|u|²`) relative to its `F_cartesian`
counterpart, their *ratios to each other* — the condition number,
`λ_max/λ_min` among the rank-matched eigenvalues — must be exactly
invariant too, not just approximately.

**Built**: `test_fim_condition_number.F` (same 10 orbits/station
geometries as `test_fim_ratio_formula`).

**Findings**: `condition_number(F_ks) = condition_number(F_cartesian)`
exactly, confirmed to `~1e-8`–`1e-11` relative precision (limited only by
eigensolver precision) across all 10 cases. KS coordinates offer no
conditioning advantage in the linear/FIM sense either — closing off
essentially every reading of "is there a real KS advantage" that a
single-instant FIM-based comparison could support.

75/75 tests passing.

### 4.6 Phase 6 — Sequential EKF Study (Inconclusive)

**Reasoning**: Phases 1–5 exhaustively address the *single-instant, linear
FIM* question. A genuinely different question remained: does actual
sequential nonlinear filtering — an EKF tracking a satellite over a real
multi-step timeline, with repeated re-linearization — show a
representation-dependent effect that a snapshot analysis structurally
cannot see? This was explicitly scoped by the user as worth a real
investment ("Full EKF simulation study") rather than deferred.

**Design principles**:
- **Isolate linearization behavior from propagator accuracy.** KS's
  propagator-accuracy advantage (smoother numerical integration near
  periapsis, etc.) is real, well-established, and *not* what's under
  test — conflating the two would make any result uninterpretable. The
  reference trajectory therefore uses an **exact** two-body Kepler flow
  map (`kepler_flow_cartesian`: `car2oe` → mean-anomaly advance →
  `oe2car`), never numerically integrated, so any effect found is
  attributable purely to linearization/filtering, not propagation error.
- **Reuse validated building blocks wherever possible.** The measurement
  update reuses the exact `range_grad_cartesian`/`range_grad_ks_analytical`
  rows already validated in Phases 1–4. State-transition matrices are
  computed by central finite difference of the exact flow map
  (`stm_cartesian`, `stm_ks`), following the same finite-differencing
  philosophy established in Phase 1 for `range_grad_ks`, rather than a
  newly hand-derived Jacobian.
- **A physically fair, matched initial condition.** `P0_cart` is a
  diffuse 1 km/1 m/s prior; `P0_ks` is *not* an independently chosen
  prior but `P0_cart` mapped through `car2ks_jacobian` (a finite-
  difference Jacobian of `car2ks` itself) — an honest test needs the same
  physical uncertainty represented in both coordinate systems, not two
  separately-guessed starting points that could bias the comparison
  either way.
- **Standard covariance Kalman recursion**, `Q = 0` (no artificial
  process noise), so the comparison isolates pure linearized-uncertainty
  transport: predict `P ← ΦPΦᵀ`; update via sequential scalar
  measurements (`src/kalman.F`, one station at a time — mathematically
  equivalent to a batch update for independent measurements, but needs
  only a scalar division, never an N×N inverse).

**Bug found #1 (real, previously latent)**: `w = dsqrt(-0.5d0*(amue/
a_kep(1)))`, used in *every* prior phase (11 occurrences across 10
files), computes `sqrt` of a negative number — `NaN`. The correct
formula, matching KSROP's own `driver_KS.F`/`test_subrouts.F` convention,
is `w = sqrt(0.25·μ/a) = sqrt(−E/2)` (both the sign *and* the coefficient
were wrong). This never affected any Phase 1–5 finding: `w` only enters
the KS *velocity* (`us`) via `ks2car`/`car2ks`; every prior computation
depends only on *position* (`u`, via `ks2car`'s `x` output, which is
`w`-independent by construction). Phase 6 is the first computation in
this repo to actually use KS velocity — which is what exposed a bug that
had been silently producing `NaN` in an unused variable for five phases.
Fixed repo-wide; the full 92-test suite (up from 75) still passes 100%
after the fix, empirically confirming the bug was genuinely inert
everywhere it had previously appeared.

**Bug found #2 (metric design flaw, caught before publishing a result)**:
An initial naive comparison of the *full* 6×6 (Cartesian) vs. 8×8 (KS)
state covariance's condition number produced an enormous, implausible
divergence (up to `8.4×10⁵`× relative difference). Diagnosis: this
comparison conflates position with velocity (different physical units in
each representation — km/s vs. the KS `us` scale) and ignores that KS's
`u`-block carries a structural gauge redundancy Cartesian's `r`-block
does not — the two full-state matrices are not dimensionally the same
kind of object even when representing identical physical uncertainty.
This is the same category of mistake Phase 1 originally diagnosed in the
*source presentation's* raw-determinant comparison, now recurring in this
repo's own draft methodology. **Fix**: Schur-complement position
marginalization — `P_pos = P_pp − P_pv·P_vv⁻¹·P_vp` — computed via a new
`mat_inverse` (Gauss-Jordan with partial pivoting, `src/linalg.F`),
yielding a position-only covariance in each representation (3×3
Cartesian, 4×4 KS, rank-matched via the same eigenvalue-threshold logic
used throughout the repo). This is the direct sequential analog of Phase
5's already-proven, dimensionally-fair metric.

**Result, even after both fixes**: **inconclusive.**
- `rank(P_pos_ks) ≠ rank(P_pos_cart)` at 71 of 143 steps (KS drops to
  rank 2 while Cartesian holds steady at rank 3).
- The condition-number ratio ranges from `~0.1×` to `~7.9×10⁴×` across
  the run — not converging toward a constant the way every Phase 1–5
  result did.
- **Decisive diagnostic**: the divergence pattern is *jagged*, not
  smooth. `cond_cart` varies gently step to step, tracking the smoothly-
  varying orbital geometry, while `cond_ks` spikes by 2–3 orders of
  magnitude between *adjacent* steps (e.g. steps 57→58→59:
  `1.6×10³ → 5.0×10⁷ → 1.9×10³`). A real physical effect tracking a
  smooth orbital trajectory should itself vary smoothly; this pattern is
  instead the classic signature of numerical sensitivity compounding
  through 144 sequential finite-difference-STM and Schur-complement-
  inversion operations.

**Disposition**: this result does **not** confirm a representation-
dependent nonlinear-filtering effect, and does **not** confirm its
absence either — unlike Phase 5's exact, reproducible-to-10⁻¹¹ result,
the finite-difference-STM approach used here is not precise enough to
distinguish a genuine effect from numerical noise at the observed scale.
A decisive answer would need an **analytically-derived** (not finite-
differenced) KS state transition matrix — a materially larger derivation
effort, left as explicitly unstarted future work rather than attempted
under this session's scope. This is reported honestly as inconclusive,
consistent with the standing principle applied throughout this
investigation: only publish findings that are clean and reproducible.

92/92 tests passing (12 new: `test_dynamics` — exact-flow round-trip,
energy conservation, and the `2·L₃(u)·(∂u/∂r) = I₃` right-inverse
identity; `test_kalman` — scalar update vs. a hand-computable 1D case;
`test_mat_inverse` — Gauss-Jordan vs. `A·A⁻¹ = I`).

---

## 5. Bugs Found and Fixed — Full Ledger

A chronological record of every real defect caught during this
investigation, in the order found:

| # | Phase | Defect | Root cause | Impact | Fix |
|---|---|---|---|---|---|
| 1 | 1 | Rank-mismatch in reduced-determinant comparison | Hardcoded "always top-3 eigenvalues" regardless of true rank | With <3 stations, corrupted comparisons, including negative ratios | Use `min(true rank)` on each side |
| 2 | 3 | Station σ set to altitude values (0.01–2.39 km) instead of published `σ=10 km` | Copy/reuse error from the settings box shared with the HEO cases | Magnitude-only (σ enters both FIMs identically) | Corrected to published value |
| 3 | 6 | `w = sqrt(-0.5·μ/a)` computes `NaN` — wrong sign *and* coefficient | Formula never cross-checked against KSROP's own convention; silently inert since nothing used KS velocity | None on Phases 1–5 (position-only); would have corrupted any future KS-velocity use | Fixed repo-wide to `sqrt(0.25·μ/a)`, all 11 occurrences |
| 4 | 6 | Full 6×6/8×8 state condition-number comparison isn't dimensionally fair | Conflates position/velocity units, ignores KS's gauge redundancy | Produced a spurious `~10⁵`×–`10⁸`× "finding" before the fix | Schur-complement position marginalization |

Notably, bugs #1 and #4 are the *same category of error* — an
unnormalized or rank-mismatched comparison producing an artificially
inflated "difference" — recurring first in the source presentation's own
methodology (which this repo exists to critique) and then, independently,
in this repo's own draft Phase 6 metric. Catching #4 before publishing a
result is itself evidence the investigation's own standards were applied
consistently, not just to the object of study.

---

## 6. Key Mathematical Results

**Rank-deficiency (Phase 1, proven analytically in `ALGORITHM.md` §5)**:
since `F_ks = JᵀF_cartesianJ` with `J = ∂r/∂u` (3×4, because physical
position `r` has only 3 independent components), the full KS-space FIM
can never exceed `rank(F_cartesian) ≤ 3`, regardless of station count or
geometry — a structural property of the transformation, not an empirical
coincidence.

**Exact closed-form ratio (Phase 4)**:
```
reduced_det_ks / det_cartesian = 64·r³
```
where `r` is the satellite's instantaneous orbital radius. Derived from
`L(u)ᵀL(u) = |u|²I₄` (KS conformality) and the regularizing identity
`|u|² = r`.

**Condition-number invariance (Phase 5, a direct corollary of the
above)**:
```
condition_number(F_ks) = condition_number(F_cartesian)   [exactly]
```
because every nonzero `F_ks` eigenvalue scales by the identical factor
`4|u|²` relative to its `F_cartesian` counterpart, so their ratios to each
other cancel that factor.

**Right-inverse identity (Phase 6, used to validate `car2ks_jacobian`)**:
since `r(u(r)) ≡ r` identically along whatever specific branch `car2ks`'s
algebraic construction picks of the KS map's gauge circle,
```
2·L₃(u)·(∂u/∂r) = I₃   [exactly, regardless of gauge branch]
```
— note this is *not* the same as matching the Moore-Penrose pseudo-
inverse `L₃(u)ᵀ/(2|u|²)`, which assumes a particular (minimum-norm)
branch that `car2ks`'s actual construction has no reason to follow; an
early draft of this test incorrectly assumed the pseudo-inverse form and
had to be corrected (see `test_dynamics.F` for the full reasoning).

---

## 7. Test Coverage Summary

92/92 passing as of the final commit (`c5107e6`):

| Test | Count | Phase | What it checks |
|---|---|---|---|
| `test_fim_cartesian` | 3 | 1 | Cartesian FIM vs. hand-computable orthogonal-LOS geometry |
| `test_fim_ks_rank` | 24 | 1 | `rank(F_ks)=min(nstat,3)` across 3 orbital regimes × 4 station counts |
| `test_fim_reduced_consistency` | 1 | 1 | Rank-corrected ratio constant across 5 station configurations |
| `test_fim_analytical_vs_fd` | 9 | 2 | Analytical vs. finite-difference gradient agreement |
| `test_timeconv` | 4 | 2 | GMST vs. published J2000.0 reference; geodetic-to-ECI sanity |
| `test_fim_geo_edge` | 4 | 3 | Rank tracking down to 1-station edge case |
| `test_fim_ratio_formula` | 20 | 4 | Exact `64r³` closed form across 10 orbits/geometries |
| `test_fim_condition_number` | 10 | 5 | Exact condition-number invariance, same 10 cases |
| `test_dynamics` | 12 | 6 | Exact-flow round-trip, energy conservation, right-inverse identity |
| `test_kalman` | 3 | 6 | Scalar Kalman update vs. hand-computable case, trace/symmetry |
| `test_mat_inverse` | 2 | 6 | Gauss-Jordan inverse vs. `A·A⁻¹=I` |

---

## 8. Final Verdict

**On the original claim**: falsified, cleanly and on multiple independent
lines of evidence. The COSPAR presentation's raw 4×4-vs-3×3 FIM
determinant comparison is not methodologically sound. The 4×4 KS-space
FIM is empirically and analytically confirmed rank-deficient (structural
rank ≤ 3) in every tested case spanning HEO and GEO regimes, 1–5 station
counts, two independent gradient derivations, with and without GMST-
accurate geometry, and 10+ varied orbital configurations. Once correctly
rank-matched, the resulting KS/Cartesian ratio behaves exactly as a
coordinate/Jacobian scaling artifact would: an exact, closed-form function
of satellite position alone (`64r³`), constant across station geometry,
and with condition number exactly invariant too. This does not mean
representation-dependent observability is never a real effect for
nonlinear systems in general — a legitimate, different research question
— but this specific presentation's evidence for a genuine "KS gives
better observability" effect does not survive a rank-correct,
same-basis, dimensionally-fair comparison, examined from every angle a
single-instant FIM analysis can support.

**On the extension question (issue #5)**: genuinely unresolved. Phase 5
closed the linear-FIM reading definitively (condition number is also
exactly invariant). Phase 6 built the actual nonlinear sequential-
filtering study needed to test the deeper question directly, but the
result was numerically inconclusive rather than a clean confirmation or
rejection. This is reported as an open problem with a clear, scoped path
to resolution (an analytical KS state transition matrix), not glossed
over or forced into a false conclusion.

---

## 9. Open Questions / Future Work

1. **Issue #5 (only remaining open issue)**: does sequential nonlinear
   filtering show a real representation-dependent effect? Resolving this
   needs an analytically-derived KS state transition matrix (replacing
   Phase 6's finite-difference `stm_ks`) precise enough to distinguish a
   genuine effect from the numerical noise that currently dominates the
   Phase 6 result. This is a nontrivial derivation (the KS EOM's
   variational/STM equations), not a quick follow-on — treat as a fresh,
   explicitly-scoped investigation if picked up again, not "one more
   phase in the same style" as Phases 1–5.
2. Not investigated and out of scope for this repo: whether *any*
   correctly-normalized sense of representation-dependent observability
   exists for nonlinear orbital filtering *in general* (beyond this one
   presentation's specific claim) — a legitimate, broader research
   question this repo's falsification does not itself resolve.

---

## 10. Reproducibility

```bash
git clone https://github.com/hari251086/KS-FIM
cd KS-FIM
fpm build --compiler ifx
fpm test --compiler ifx            # 92/92 passing
fpm run ksfim_case_study --compiler ifx
fpm run ksfim_timeseries --compiler ifx
fpm run ksfim_ekf_study --compiler ifx
```

Requires Intel `ifx` 2025.0 (Visual Studio + Intel oneAPI environment
variables in the same shell session on Windows); not yet verified with
`gfortran`. Dependency: KSROP (`git+https://github.com/hari251086/KSROP`,
tag `v2.2.0`).

Full technical writeup with equations and per-phase validation detail:
`ALGORITHM.md`. Project overview and version history: `README.md`. Durable
GitHub record of the findings: issues #6 (verdict, updated every phase)
and #5 (the still-open extension question).
