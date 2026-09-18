# NEMO — Phase II Rung 5: GTODN site enumeration & classification

**Status:** decisions LOCKED by the design thread (all six flags resolved). This is
the checked artifact for the Rung-5 site classification and the record of decisions.
**Branch/HEAD at analysis:** `phase2-coagulation` @ `8f546bc`.
**Scope:** v1 = **2-phase** (`is_3_phase = 0`), and — now enforced by a hard scope
guard (Flag ③④⑤) — `modify_rate_flag = 0`, `is_er_cir = 0` whenever coagulation is on.
**Method:** every read of `GTODN` / `GTODN_FIXED` / the derived `ab_lay`, `SUMLAY`
in the RHS builders was traced to its reaction type and reduced to the **net power
of the grain abundance** `Y_tot(k) = Y(GRAIN_k⁰)+Y(GRAIN_k⁻)` in the resulting
`reaction_rates(J)` coefficient, after algebraic cancellations.

---

## 0. The form (adopted) and the corrected floor rule

For any site whose coefficient obeys `coeff ∝ Y_tot(k)^ν` with `ν` locally
constant, the missing column-`GRAIN_k` Jacobian entry is

```
∂F/∂Y(GRAIN_k⁰) = ∂F/∂Y(GRAIN_k⁻) = ν · F / Y_tot(k)
```

`F` = reaction flux (the `tmp_value` `get_jacobian` already forms); both charge
columns identical because the divisor is the total. `ν = +1` linear, `ν = −1`
reciprocal. This is adopted as the canonical form (it subsumes the brief's
`−rate/Y²`).

**Floor rule (corrected — load-bearing).** The floor indicator is **not** keyed on
the sign of `ν`. It is keyed on **whether the coefficient routes through a
`max(Y_tot, floor)` divisor**:

> **Apply the floor (entry = 0 below floor) wherever `max(Y_tot, floor)` appears in
> the coefficient; otherwise no floor.**

- Accretion: bare `×Y_tot`, divides by nothing → **no floor** (`ν=+1`).
- LH: `×GTODN = ×1/max(Y_tot,floor)` → **floored** (`ν=−1`).
- Photodesorption: routes through the floored `SUMLAY` → **floored** even though
  `ν=+1`. (My earlier "floor only for `ν<0`" was correct for the A/B pair but
  breaks here; the divisor-routing rule is exact.)

Below the floor, `max(Y_tot,floor)` is constant, so any coefficient routed through
it has zero grain-derivative there.

---

## 1. Full enumeration

Lines are HEAD `8f546bc`, `src/ode_solver.f90` unless noted.

| # | line(s) | rxn type | how `Y_tot(k)` enters | ν | floor? | class | v1 | Rung-5 action |
|---|---------|----------|-----------------------|---|--------|-------|----|---------------|
| A | 1859–1860 | 99 accretion | explicit `/GTODN` ⇒ `×Y_tot` | +1 | no | **linear** | ✔ | live read + analytic `+F/Y_tot`, both cols |
| B | 2049 | 14 LH surface | explicit `×GTODN` ⇒ `×1/Y_tot` | −1 | yes | **reciprocal (stiff)** | ✔ | floored live read + analytic `−F/Y_tot`, 0 below floor, both cols; gated `modify_rate_flag==0` |
| C | 1884 | 66 photodesorp (ext UV) | via `MLAY/SUMLAY(k)` (floored), cap active when `SUMLAY≥MLAY` | +1 | yes (via SUMLAY) | **in-scope (Flag ②)** | ✔ | RHS read already live (Rung 3); **add** analytic `+F/Y_tot`, 0 below floor, 0 when cap inactive; both cols |
| D | 1903 | 67 photodesorp (CR UV) | same as C | +1 | yes (via SUMLAY) | **in-scope (Flag ②)** | ✔ | same as C |
| E | 1534 | (derived) `ab_lay = nb_sites/GTODN`, **frozen** | feeds A caps, ER | +1 | — | **frozen partial-coupling (Flag ①)** | ✔ | none; documented + cap-binding diagnostic |
| F | 1552–1556 | (derived) `SUMLAY/sumlaysurf/sumlaymant`, live+floored | Rung 3 | — | — | **Rung 3 RHS; entries via C/D** | ✔ | do **not** re-touch RHS read; entries supplied through C/D |
| G | 2099 | 30 ER, `ab_surf≤ab_lay` | `TSQ/GTODN/ab_lay`, `ab_lay=nb_sites/GTODN` ⇒ **cancels** | 0 | — | **cancels** | ✖ guard | none (scope guard) |
| H | 2102 | 30 ER, `ab_surf>ab_lay` | `TSQ/GTODN/ab_surf` ⇒ `×Y_tot` | +1 | no | **deferred v2** | ✖ guard | none (scope guard) |
| I | 2434,2439 | 14 (modify test) | `×GTODN` with `Y·GTODN<1` threshold | piecewise | — | **deferred v2** | ✖ guard | none (scope guard) |
| J | 2445,2451 | 14 (modify rewrite) | `/GTODN` into `DIFF` ⇒ cancels outer `×GTODN` | 0 (PICK) | — | **deferred v2** | ✖ guard | none (scope guard) |
| K | 2203–2204 | (derived) 3-phase `sumlay*` frozen | 3-phase | — | — | **deferred v2** | ✖ guard | none (scope guard) |
| L | 2256 | 40 surf→mant | `alpha ∝ ab_surf·GTODN/nb_sites` then `/ab_surf`, tangled | mixed | — | **deferred v2** | ✖ guard | none (scope guard) |
| M | 988,992 | 10/11 H2 single-grain | scalar `GTODN_FIXED`, `IS_GRAIN_REACTIONS==0` | n/a | — | **out of scope** | ✖ | none |

**In-scope set = three reaction classes, one form:** accretion (99, `ν=+1`, no
floor), photodesorption (66/67, `ν=+1`, floored-via-SUMLAY), LH (14, `ν=−1`,
floored). Sub-finding (LH mantle branch): when a reaction carries a mantle (`'K'`)
compound and `sumlaymant>1`, the LH `×GTODN` and the `/sumlaymant` divide (both
routed through the same floored divisor) **cancel to `ν=0`**. Mantle (`'K'`)
species exist only in the full/3-phase network — `reduced_CHO` and the fiducial
fixture have **zero** — so this is unreachable in v1 and is additionally blocked by
the `is_3_phase=0` guard. **Recommendation:** when coding B, mirror the RHS
`sumlaymant` branch (zero the entry when the divide is active) so 021 stays
FD-exact if a mantle-bearing network is ever run; cheap insurance, not a v1
requirement.

---

## 2. Decisions recorded against each flag

**Unifying form — ADOPTED.** `ν·F/Y_tot`, both charge columns identical,
`−F/Y_tot` (F = flux = `tmp_value`) for the reciprocal. Floor rule restated as
"floor wherever `max(Y_tot,floor)` appears" (§0).

**② SUMLAY Jacobian entries — 021 OWES THEM (forced, not optional).** Rung 4's
GRAIN_RANK superset already registered the 66/67 `(compound, GRAIN_k)` slots; the
FD oracle (022) fills them; production 021 is MITER=1 all-analytic with no
per-entry hybrid, so every registered, FD-filled entry must be analytic in 021 or
gate 5a fails. Photodesorption is a `ν=+1` site (same `+F/Y_tot` as accretion) but
routed through the floored `SUMLAY`, so it carries the floor indicator (and the
`SUMLAY≥MLAY` cap gate). Consequences recorded:
- **Landmine 7 corrected (design thread owns the erratum):** "monolayer count is
  already done" referred to the *RHS read* (Rung 3), which was conflated with the
  Jacobian entry. The entry never existed. Corrected rule: **do not re-touch the
  SUMLAY RHS read; DO supply its Jacobian entries.** Adding
  `d(photodesorption)/dY(GRAIN_k)` is not re-touching the count and is not a
  double-count (each reaction gets its own entry).
- **Gate 5a must include an ice-rich state**, not only the ice-free dust-only
  fiducial (where `SUMLAY≈0` leaves C/D unexercised). The ice-rich state is added
  to the 5a state set.

**③④⑤ confounders — ONE SCOPE GUARD.** Under coagulation, **hard-assert (stop, not
warn)** `modify_rate_flag==0 .and. is_er_cir==0 .and. is_3_phase==0`, with a
message naming the unsupported combination. This makes v1 scope explicit and
enforced (fails loudly rather than silently mis-computing a Jacobian), and matches
the paper's stated v1 scope. The LH reciprocal entry (B) is additionally gated on
`modify_rate_flag==0` (locally safe/redundant). Modified-LH / ER / 3-phase
Jacobians are **deferred to v2** explicitly (rows G–L).

**① ab_lay — FROZEN (confirmed), with two conditions.**
- (a) Recorded here as a **deliberate partial-coupling**: a frozen surface-site
  cap, the analogue of Rung 3's frozen surface area, and the v2 hook. Not silent.
- (b) The "second-order" claim is to be **checked, not asserted**: add a diagnostic
  that the accretion monolayer cap (`ab_lay`) rarely binds in the fiducial runs.
  Rationale: coagulation reduces total grain surface (`∝a²` falls faster than
  number rises), so a frozen `ab_lay` over-estimates available surface — harmless
  only if the cap seldom binds. If the diagnostic shows frequent binding,
  **escalate**; otherwise frozen is validated.

**⑥ grain_abundance_floor — PROMOTE TO A `parameters.in` KNOB (default 2.99e-39).**
Reasons: (a) enables the floor-insensitivity study (0.1× / 10×) without
recompiling — the actual proof the floor is chemically negligible; (b) one knob
enforces that the Rung-3 divisor floor and the Rung-5 rate-site floors share a
single value; (c) the paper implies tunability. Default unchanged, so existing
runs are unaffected. Paper wording to be corrected by the author: "set empirically
from the coupled models" → "set low enough to be chemically negligible — verified
by insensitivity of the results over [range] — and high enough to keep the
reciprocal terms finite" (it is a physically-motivated ~10 grains/AU³ threshold,
not a data fit). The insensitivity run is a validation-phase deliverable, noted now.

---

## 3. Net scope (as decided)

Three reaction classes (99 no-floor, 66/67 floored-via-SUMLAY, 14 floored), one
derivative form (`ν·F/Y_tot`, floor where a `max(Y_tot,floor)` divisor appears),
one confounder scope-guard (hard-assert under coag), `ab_lay` frozen + one
cap-binding diagnostic, and `grain_abundance_floor` promoted to a shared knob.
Gate 5a state set extended to include an ice-rich state.