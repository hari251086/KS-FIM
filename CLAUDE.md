# CLAUDE.md

Cross-repo rules live in `GitHub\CLAUDE.md` — not repeated here, only
what's specific to KS-FIM. Fortran-specific gotchas are in the shared,
tree-walked `.claude/rules/fortran-ks-gotchas.md`.

## Project

Research repo testing the 2021 COSPAR claim that KS-space Fisher
Information Matrix observability beats Cartesian for orbit determination.
Result: **falsifies** the claim (see `ALGORITHM.md`). Phase 6 sequential
EKF-style study (issue #5) is inconclusive, not a confirmed finding — an
analytical KS state-transition matrix is still needed to settle it.

## Build / Test

```bash
fpm build --compiler ifx
fpm test  --compiler ifx
fpm run ksfim_case_study    --compiler ifx
fpm run ksfim_timeseries    --compiler ifx
fpm run ksfim_ekf_study     --compiler ifx
```
Only `ifx` is verified — **not yet verified with `gfortran`**, unlike
most other repos in this lineage. No CI configured.
92/92 tests passing.

## Key code

- `src/` / `app/` — case-study, timeseries, and EKF-study drivers
- `test/` — 9 suites (`test_fim_cartesian`, `test_fim_ks_rank`,
  `test_fim_analytical_vs_fd`, etc.)
- All inputs (orbital elements, station coordinates) are hardcoded from
  the source presentation's own published tables directly in `app/*.F`
  — no CLI args or config files.

## Always / never

- Don't report the Phase 6 EKF result as a confirmed finding — it's
  explicitly inconclusive (see `ALGORITHM.md` §5 for why) pending an
  analytical KS STM (issue #5, open).
- Don't assume `gfortran` builds cleanly here without checking first —
  it's never been verified on this repo, unlike KSROP/KS-Pc/OREM.
