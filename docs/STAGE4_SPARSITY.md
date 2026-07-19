# Stage 4, note: symbolic sparsity cannot be bit-identical (and why that's fine)

**Status of the code as of this note:** NEMO reproduces NMGC-2.0 bit-identically
(the stage-2 equivalence result, worst relative difference 0). Nothing here has
changed that. This note is about a change Task 4 *proposes*, not a defect in the
current code.

## What Task 4 asks for

Task 4 has two independent parts:

1. **Precompute the RHS index maps.** `set_dependant_rates` recovers a grain rank
   per call via ~98 string comparisons and 14 `read(c_i,'(I2)')` parses, even
   though `GRAIN_RANK(nb_reactions)` is already built once at init and already
   used elsewhere. Wire the RHS to the precomputed map instead of re-deriving it.
2. **Symbolic Jacobian sparsity.** Replace the numerical sparsity inference in
   `set_work_arrays`/`count_nonzeros` (evaluate the Jacobian, keep entries with
   `|value| > 1e-99`) with a pattern built once from the reaction list.

The spec asks to re-verify **bit-identical** after both. Part 1 can meet that.
**Part 2 cannot**, for a structural reason measured below.

## The measurement

Instrumented `set_work_arrays` on the reference 0D 2-grain case (`-O0`) to record
the pattern it hands DLSODES at every solver step, and built the symbolic pattern
from the reaction list three ways for comparison.

Per-step **numerical** pattern (what the current code actually uses), 9 steps:

| step | nnz (`>1e-99`) | nnz (`≠0` exactly) | pattern changed? |
|---|---|---|---|
| 1 | 27128 | 27160 | — |
| 2–3 | 27274 | 27288 | yes |
| 4–6 | 27276 | 27290 | yes |
| 7 | 27275 | 27290 | yes |
| 8 | 27275 | 27290 | yes (different again) |
| 9 | 27274 | 27290 | yes |

At least **5 distinct sparsity structures across 9 steps.** The pattern is
*value-dependent*: the `1e-99` cut drops structurally-nonzero-but-tiny entries,
and which entries are tiny depends on the abundances, which move along the run.
(Note `nnz(>1e-99) < nnz(≠0)` at every step by 14–32: the cut removes genuinely
nonzero entries, not just exact zeros, and the removed set shifts step to step.)

**Symbolic** pattern from the reaction list (built once):

| pattern | nnz |
|---|---|
| naive — mark every compound of every reaction using species *j* | 29179 |
| + within-reaction cancellation (a species that is both a reactant and a product gets `+tmp` and `−tmp`, exact bit-zero) | 29176 (−3) |
| + drop identically-zero-rate (disabled) reaction channels | 27296 (−1880) |
| numerical structural union (max `nnz(≠0)` over steps) | 27290 |

So a *rate-aware* symbolic build lands within **6 entries** of the numerical
structural pattern — getting symbolic sparsity *correct* is easy. The within-
reaction cancellation subtlety is real but tiny (3 entries); the large gap versus
the naive build (1880) is disabled/zero-rate channels.

## The conclusion

Symbolic sparsity is inherently a **bit-changing** optimization. Its purpose is to
replace a per-step-varying, threshold-pruned structure with *one fixed superset*
that DLSODES factorizes once and reuses — that reuse is the speedup, and it
necessarily changes the last bits of the trajectory. A fixed pattern cannot equal
a pattern that changes every step, so no amount of care makes Part 2 bit-identical
to the current code (and therefore to NMGC-2.0).

This means the literal Task-4 requirement — *symbolic sparsity **and** re-verify
bit-identical* — is internally inconsistent for Part 2. The two parts need
**different** acceptance tests:

* **Part 1 (index maps):** pure refactor → gate on `cmp` byte-equality against a
  baseline captured from the current code (equivalently, the equivalence test).
* **Part 2 (symbolic sparsity):** gate on (i) the symbolic pattern being a proven
  **superset** of the numerical pattern at every step (correctness: nothing real
  dropped), and (ii) the reference trajectory agreeing to a tight **relative
  tolerance** (suggest ~1e-10), *not* `cmp`. The derived round-trip
  (`tests/roundtrip_derived.sh`) stays byte-identical regardless, since it is a
  self-consistency check unaffected by the sparsity method.

## Decision needed from the design thread

Confirm the Part 2 acceptance criterion (superset + tight-rtol, in place of
bit-identical). Once confirmed, Part 1 can be implemented and proven bit-identical
immediately; Part 2 follows under the relaxed gate.

(Reproducing the measurement: the instrumentation was a throwaway block in
`set_work_arrays` plus a `sparsity_symbolic_diagnostic` routine; it was reverted
and is not in the tree. It can be reinstated from this note if the numbers need
re-checking on a different case.)
