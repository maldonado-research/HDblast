# Separating the forced increment in the retained formal second response

This conditional corollary uses the accepted causal construction, its
registered gauge and second initial prescription, and the actual bulk/shell
source bounds. It does not compute a waveform or assert an exact nonlinear
solution, a nonzero forced increment, or an invariant backreaction observable.

On `D={0<=x<=3, r>=0, x+r<=8}`, the complete 23-amplitude formal response
satisfies

    U=U_free+KF+LG+TU,       ||T||_X<=q=2109/1000000,
    ||V||_X=max(sup|V|_H/3,sup|V_x|_H,sup|V_r|_H).

H uses probability-normalized sphere measure and the positive metric
Frobenius norm plus the scalar square. The full operator is
`TV=K(A V_x+B V_r+D0 V)+L(M P_N V|r=0)`. K uses the registered odd/even
reflection `J=P_N-P_D`; `Lg(x,r)=-integral_0^max(x-r,0)g(s)ds`.
The retained estimates are `||Kf||_X<=3 sup|f|_H` and
`||Lg||_X<=sup|g|_H`. All physical and normal harmonic shell sources occur in

    G=-ell rho(-2QK,-Qp,GN),       GN=-2 tr_gamma QK-2QC_y.

The actual first master retains its dynamic shell construction and supplies
the actual F,G. No shell contribution or Robin term is omitted below.
The retained data and sources vanish on the flat initial shell wedge; in
particular G(0)=0. K and L then produce C1 fields with zero Cauchy data,
and T preserves this compatible class, so the decomposition is in X.

**Proposition.** Define the auxiliary reduced-operator decomposition

    S=KF+LG,
    U_hom=(I-T)^(-1)U_free,       W=(I-T)^(-1)S.

Then `U=U_hom+W`, W has zero initial data, and the accepted bounds
`F_H=17387513156573`, `G_H=3851400` imply

    ||W||_X <= (3F_H+G_H)/(1-q)
             =52162543321119000000/997891 <5.228e13,              (1)
    sup_D |mean_w(W)| <=156487629963357000000/997891
                       <1.5684e14.                              (2)

**Proof.** Since `q<1`, the same bounded inverse applies to both right-hand
sides and their sum. The complete source obeys
`||S||_X<=3F_H+G_H=52162543321119`. The resolvent norm is at most
`1/(1-q)`, proving (1). The sphere mean functional has H-norm one, and the
X norm assigns values weight three, proving (2). No angular pointwise
evaluation factor enters. The proof includes LG and the Robin part of T.

The accepted evolved scalar mean has opposite strict margins `1.85e30` on
these complete rectangles:

| Mean sign | Scaled time x | Scaled radius r |
|---|---|---|
| Positive | `[1/1000,1/200]` | `[73/64,585/512]` |
| Negative | `[1/1000,1/200]` | `[951/512,119/64]` |

The exact upper bound (2), divided by this certified margin, is less than
`8.477e-17`. Therefore these signs persist in U_hom. This is a comparison
of an upper bound with a certified margin, not a measured ratio of two
responses. The large mean signs do not resolve the forced increment W.

U_hom is not U_free: it retains all homogeneous propagation, couplings and
Robin terms. The registered second initial data can already encode quadratic
constraint and gauge-completion information from the first pair. Zeroing
F,G does not make U_hom another physical Einstein–scalar shell solution.
Neither U_hom nor W is claimed separately to satisfy the complete sourced
physical constraint system.

This global bound does not exclude local importance of W where U_hom is
small or vanishes, or within a weak modal sector. The bound alone does not
establish that W is nonzero or determine its sign or profile. All conclusions
remain conditional on the accepted smooth construction and estimates.

For a separate numerical computation of W, the exact two-term approximation
obeys

    ||W-(S+TS)||_X <=q^2(3F_H+G_H)/(1-q)<2.326e8,                (3)

because `W-(S+TS)=T^2W`. This tail excludes input, quadrature and
reconstruction errors. Those must be added, or included once in a full
integral residual. A total-response error certificate of `1e29` in X can
certify preservation of the large mean signs but cannot certify resolving W
at the scale (1). It does not prove that the actual numerical error is large
or that a more accurate computation is impossible.

The accompanying method note gives a continuous C0-source-to-X residual
route without derivatives of actual F, subject to complete source error and
corner compatibility checks. An actual W feature must exceed its certified
error before a nonzero or sign claim is justified. No such candidate or
feature is produced here.

`check_budget.py` authenticates local copies of the accepted reports and
checks the rational inequalities. It does not replay their expensive
scientific proofs. `INPUT_PINS.json` records their hashes and original
locations; normal and optimized isolated runs agree. These files are
separate from the immutable geometric-source capsule.
