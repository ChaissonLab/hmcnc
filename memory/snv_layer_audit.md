---
name: SNV Layer Audit
description: Findings from audit of the SNV/BAF layer in hmmcnc.cpp — what's working, dead code, and bugs
type: project
originSessionId: 24fd0286-cf0c-4690-a12c-65c01a50c958
---
Audit of the SNV layer in `src/hmmcnc.cpp` as of Phase 5.

## What works

- **`StoreSNVs`** (line ~1359) — extracts heterozygous positions from BAM pileup vs reference; only runs in BAM mode (`-a`), not bed-input mode (`-b`).
- **Post-hoc SNV scoring** (line ~3543–3570): after Viterbi calls, each called interval accumulates `binoP[curCN-1][totCov][alt]` vs diploid baseline `binoP[1][totCov][alt]`; sets `FAIL` if CN model fits worse than diploid.

## Dead code

- **`CombineEmissions`** (line 508) — defined, never called.
- **`CSEmisP`** (line 537) — defined, never called.
- **`covSnvTransP` / `snvSnvTransP`** — initialized in `InitParams` but never passed to any F-B or BW function. The BW E/M steps are purely coverage + clipping; SNVs have zero effect on training.

## Bugs

- **`binoP` is only 3 states** (`binoP.resize(3)`, indices 0/1/2 = CN=1/2/3). Post-hoc scorer does `binoP[curCN-1][totCov][alt]`, which is out-of-bounds for CN≥4 (silent undefined behavior).
- **`alt` vs `min(ref,alt)` inconsistency**: `binoP` was designed for `m = min(ref,alt)` as the symmetric minor-allele count (in `CSEmisP`), but the post-hoc scorer passes the raw `alt` count directly.
- **Hardcoded BAF p values**: CN=1→0.1, CN=2→0.5, CN=3→0.66. CN=3 at 0.66 only models the 2-alt/3-total orientation; the 1-alt/3-total orientation (p=0.33) is not covered.

## Possible directions (not yet decided)

1. Fix post-hoc filter only — correct 3-state bound, use `min(ref,alt)`, model both CN=3 orientations.
2. Remove SNV scoring entirely — it's not contributing reliably; simplify post-processing.
3. Integrate BAF into HMM emissions — add per-bin BAF alongside coverage; `CombineEmissions`/`CSEmisP` were designed for this but never wired in.
4. Use external SNV VCF instead of pileup-based `StoreSNVs` (better for long reads with higher base error rate).

**Why:** SNV layer was audited to assess whether it's contributing to calls or needs redesign before Phase 6 validation.
**How to apply:** Before touching the SNV layer, confirm which direction the user wants to pursue.
