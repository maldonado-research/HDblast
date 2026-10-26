# Scope and defect review of the full23 response candidate

26 October 2026 UTC. Separate internal AI-assisted review. Scope: C1 regularity, causal support, applicability of the retained full-space contraction and angular evaluation constant, the two-term defect identity, and treatment of the finite Robin mismatch. This review read `compute/full23_response.py`, `compute/angular_ridge.py`, their normal/optimized reports and profiles, the initial profile definition, and the retained primary proofs. It did not rerun a frozen producer, re-prove upstream PDE existence, or audit the initial/coefficient interval machinery assigned to other reviewers.

**Conclusion:** no blocking scope or defect-logic issue found. The saved result is a conditional full23 C1 response approximation with error at most the exact rational number in `total_X_error_upper` (approximately `1.1787044290072725e28`). Its accuracy is established by a whole-cone error ledger; the six saved samples only evaluate that candidate. It does not certify pointwise second-derivative Einstein residuals or an exact nonlinear family.

## Domain and support

The domain is precisely the retained

`D+={0<=x<=3/2, max(0,x-5/4)<=r<=7/4-x}`.

For a backward characteristic predecessor `(s,v)` after half-line reflection, `s+|v|<=x+r` and `s-|v|<=x-r`; these inequalities keep it inside the same domain. The shell history ends at `x=5/4`, so no unprovided shell data at later times enter the candidate. Initial evaluation arguments have magnitude at most `7/4`.

The actual first/second data and exact sources have the retained zero wedge `x+r<1`. The nominal displacement and all three continuous P1 profiles are identically zero through `r=65/64`; both endpoints of the first large P1 cell are exactly zero. Their free wave and first feedback consequently vanish on the stronger wedge `x+r<65/64`. Thus the error with the actual response vanishes on `x+r<1` and is eligible for the retained zero-wedge estimates. No global exterior support assertion beyond the registered causal footprint is required or made.

The independent bounded checker `check_scope.py` verifies the zero profile rows, TT initial table zeros, domain endpoints and ordering. It also checks that the cutoff polynomial has value zero/one and matching derivatives through order four at its joins.

## C1 regularity and boundary compatibility

The nominal displacement is the expression-defined rational profile, not its P1 table. Its registered formulas involve at most three derivatives of a packet with C4 cutoff joins, so it is C1 across those joins. The candidate velocity is continuous P1 with shared node values. Both have a zero neighborhood of the origin, so their specified odd/even half-line reflections create no origin regularity defect. The d'Alembert free wave is therefore C1 on the physical domain.

The independent P1 displacement and radial-derivative tables used in feedback are continuous; the latter is not asserted to be the derivative of the former. These tables enter only the continuous **forcing** `F(r+x)+G(r-x)` and continuous shell lift datum. The reflected Duhamel operator maps such forcing to C1 on the physical half-domain, with the standard one-sided boundary interpretation; an odd reflected forcing may have a jump at the shell without invalidating this physical C1 conclusion. The boundary lift is C1 because its continuous datum vanishes near the initial corner. Thus `V=B0+B1` is C1, has the homogeneous Dirichlet rows exactly, and belongs to the zero-wedge contraction space. No C2 or globally smooth claim is justified merely by these P1 profiles.

The candidate need not satisfy the final Robin equation exactly. The Banach space underlying the retained integral equation imposes homogeneous Dirichlet and the zero wedge; Robin data are imposed by its boundary-lift map. This is why a finite feedback candidate is admissible for the error argument.

## Full-H contraction and two-term error

`T/analytic/LOCAL_TOTAL_REVIEW.md`, lines 43–80, explicitly proves the zero-wedge value and trace estimates in the **full retained angular Hilbert space**, not only its mean block:

`sup_D+ |W|_H <= (3/4)||W||_X`,
`sup_shell |W|_H <= (1/4)||W||_X`,
`||K f||_X <= (3/2)sup|f|_H`, and `||L g||_X <= sup|g|_H`.

With the retained operator bounds these give

`q=(3/2)(2/10000+3/(4*1000000))+1/(4*2000)=3409/8000000<1`.

Using `b=U_free+KF+LG`, `U=b+TU`, and `V=B0+B1`, exact algebra gives

`U-V=T^2 U+(I+T)(b-B0)+(TB0-B1)`.

Consequently `E=q^2*Ucap+(1+q)*e0+e1` is valid when the three component bounds are valid. The new `Ucap=(X0+e0)/(1-q)` follows directly from `||b||<=X0+e0`; it does not import a mean-only response cap as a full23 estimate.

The entire omitted exact forcing is charged to `e0` through the retained full-H bound

`S=(3/2)*17387513156573+3851400=52162547172519/2`.

This omission includes the actual TT forcing. Zero candidate TT slots therefore do not assert a zero actual TT response. The target can legitimately pass without resolving the much smaller forced tensor signal.

For `TB0-B1`, the exact displacement and its derivative are replaced only within feedback; the same P1 velocity is already used in B0 and B1. Thus no extra velocity interpolation term is missing from the feedback profile error. The coefficient errors multiply B0 with the correct zero-wedge factors; the independent profile errors enter as

`df_profile=(nA+nB)*epsUr+nD*epsU`,
`dg_profile=nM*epsU`.

Applying K and L gives `e1=(3/2)*df+dg`. This is the two-term identity route; no extra omitted-iteration residual or inverse factor should be appended to its final error. Initial/coefficient interval correctness and the ridge formulas themselves remain the separate numerical review obligations.

## Finite Robin mismatch

At the shell, the free Neumann part has zero normal derivative and the first lift supplies `barM*hatB0`. The candidate's actual Robin mismatch is therefore

`R_b(V)=PN V_r-M PN V-G = barM*hatB0-M*(B0+B1)-G`.

Its terms are not silently set to zero: first-step coefficient/profile errors enter e1, omitted exact forcing enters e0, and missing feedback enters the controlled tail. A useful additional **conditional consequence**, without directly evaluating this residual, is obtained from the exact target's Robin law:

`R_b(V)=-(PN (U-V)_r-M PN(U-V))`,
`sup_shell |R_b(V)|_H <= (1+||M||/4) E <= (1+1/8000) E`.

For the saved E this is approximately `1.1788517670608983e28`. The new `SCOPE_CHECKS.json` records the exact rational bound. It is a norm consequence of the response error certificate and retained exact physical boundary theorem; it is not an independently evaluated cancellation or a pointwise Einstein constraint check.

## Angular evaluation and physical scope

The retained N selected-space closure proves that the exact target belongs to the 23-feature space E_v. The full23 reconstruction and its x/r derivatives use precisely those same fixed angular templates. Thus the error also lies in E_v, and the probability-normalized positive Frobenius evaluation theorem `C*=sqrt(321/19)` applies. Its primary sources are N/evaluation-bound/EVALUATION_BOUND.md and N/angular-review/ANGULAR_EVALUATION_REVIEW.md. No dimension-only homogeneous-kernel argument is being used.

Combining it with the zero wedge gives `sup|U-V|_F <= (3/4) C* E` and `sup|(U-V)_x|_F, sup|(U-V)_r|_F <= C* E`. These are amplitude/frame-component scaled coordinate derivatives. Physical covariant or ordinary-time derivative bounds require the appropriate `ell`, `rho` and connection conversions. The saved per-amplitude sample enclosures use the weaker factor `3E/sqrt(g_i)` for values and `E/sqrt(g_i)` for derivatives; this is conservative and consistent.

The C1 certificate permits first-derivative gauge or boundary residual bounds only after including those conversions and quadratic-gauge source errors. It does not bound pointwise second derivatives, so it does not directly certify pointwise Einstein residuals. Quantitative C2 control or a declared weak norm would be required for that stronger claim. All inherited exact subsidiary and physical completion assumptions remain conditional.

## Derived integral residual and intrinsic angular derivative

The reviewed lemmas `INTEGRAL_RESIDUAL_LEMMA.md` and `ANGULAR_DERIVATIVE_LEMMA.md` give two further conditional consequences without a new evolution. For the same candidate, the complete integral residual satisfies `R=e0+(1+q)e1+q^2 X0`; hence its resolvent error is `R/(1-q)`. Exact arithmetic checks `R/(1-q)-E=2q e1/(1-q)>=0`, and both error bounds pass the `10^29` target. This residual is bounded through the certified defect ledger; it is not a direct pointwise Einstein residual evaluation. No extra tail should be added to its resolvent bound.

For the error's intrinsic angular derivative `D_kappa` on the unit S3, the SO(4)-invariant ambient section space has dimension 287 and rough spectral cap 24. The derivative image is invariant of dimension at most 287, so its homogeneous evaluation bound and the zero wedge give `(3/4)sqrt(287*24)E < (255/4)E`. This is a positive Frobenius norm summed over angular derivative directions and tensor slots. It is not the physical spacetime covariant derivative, which needs warp and connection conversions. The extended `SCOPE_CHECKS.json` records exact rational integral residual, resolvent, difference-identity and intrinsic-angular bounds.

## Executed checks and limits

`check_scope.py` uses only the standard library, imports no project module and invokes no old producer. It independently recomputes q, S, e0/e1/Ucap/total-error joins, the rational C* upper bound, both pointwise error factors, cutoff endpoint jets, zero initial profile support and matching normal/optimized response reports. `SCOPE_CHECKS.json` records the examined code/report/profile hashes and exact conclusions. These checks passed. They do not replace the other reviewers' initial, interval-coefficient, matrix-norm or exact-kernel audits.

The portable checker requires a fresh output outside the capsule. For example, from the capsule root:

```bash
python -I -S -B structure/response-review/check_scope.py \
  --capsule-root . \
  --normal-results compute/response-normal \
  --optimized-results compute/response-optimized \
  --output /tmp/fresh-response-scope.json
```

`--normal-results` defaults to `compute/response-normal` under the capsule root. Omit `--optimized-results` to skip only the normal/optimized byte comparison; a reproduction wrapper can supply freshly generated result directories. No absolute input or output path is embedded in the report. The portable I/O change preserves the earlier fields; the subsequently requested derived-bound extension adds the exact integral-residual and intrinsic-angular consequences above.
