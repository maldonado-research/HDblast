# Independent first-shell source and algebra review

30 October 2026. The reviewed first-shell implementation binding passes.
No physical sign, unit-normal, tension-normalization, ell-scaling, or
tangential-gauge cancellation defect was found in the final checker.

This review independently derived the direct physical scalar bridge,
reviewed the source/AST joins and checker algebra, and ran fresh normal and
optimized subprocesses with all controls. It does not import the checker
as a scientific library or claim to supply another trace solution.

## Direct derivation

Use the copy-even one-copy outward residuals
`R_ab=K_ab−sigma*gamma_ab/6` and
`P=n(phi)+sigma_prime/2`. The background has `sigma=6a`,
`sigma_prime=−2p`. For the longitudinal first field,
`h_yy=−4uY`, `h_ab=2uY*gamma_ab` and the scalar amplitude chi,
variation of the unit normal gives

```
R1_ab=(u_y+2a*u+p*chi/3)*Y*gamma_ab,
P1=chi_y+2p*u+(sigma_second/2)*chi.
```

The metric factor vanishes directly from
`chi=−3(rho^2*u)_y/(rho^2*p)`. The term `2p*u` in P1 is the normal
lapse contribution and is essential.

Set `h=V_prime/p`, `s=sigma_second/2`,
`k=rho*(5a/2−h)`, and `C=3alpha/(rho^3*p)`. The authenticated alpha and
background recurrences imply

```
alpha_y/alpha=h−5a/2,
C_y/C=−3a/2,
k_y+k^2/rho=rho*Q−2rho*p^2/3−5/(2rho),
Q=2h^2−10a*h+35a^2/4−V_second−p^2/6.
```

Applying these directly to
`chi=C*theta*(W_r/ell+k*W)` gives, before imposing any field equation,

```
P1 / [3alpha*theta/(ell^2*rho^4*p)]
 = −W_rr+ell*rho*Bcal*W_r
   +ell^2*(rho^2*Q+rho*Bcal*k−5/2)*W,
Bcal=h−4a+s.
```

The retained shell coefficients are
`mu^−1=rho*Bcal`, `gamma=rho*(h−5a/2)=−k`. Therefore the preceding
expression is exactly `Rmaster−Dtrace`, where

```
Rmaster=W_xx−W_rr+ell^2*(rho^2*Q−35Hubble^2/4)*W,
Dtrace=W_xx−(ell/mu)*(W_r−ell*gamma*W)
       +ell^2*(5/2−35Hubble^2/4)*W.
```

This derives the final checker identity independently of its polynomial
engine. The W_xx and Hubble terms cancel in the difference. No bulk PDE was
used to produce the direct physical normal condition. Both the actual
bulk master PDE and the actual dynamic trace law are needed to conclude
that its right hand side vanishes. The formal jet witness W=Wr=0,
Wxx=Wrr=1 correctly leaves `Rmaster=0`, `Dtrace=1`, and nonzero P1.

## Current source and gauge join

The final checker authenticates all listed bytes and exact tree membership
for the geometric-source (168 members), interior-master (40 members), and
regular-shell-spectra (129 members) closures, before AST use and after
readback. It rejects capsule symlinks and tree symlinks. It executes only
isolated authenticated function/expression ASTs in exact arithmetic.

The source `frame_amplitudes` scalar coefficient, background R/B/k and
alpha recurrence, actual Q expression, and registered tension polynomial
are correctly isolated. The fixed background certificate's exact
gamma/inverse-mu/Bcal definitions are explicitly checked and evaluated;
these are not independent choices of interval endpoints. The actual
affine detuning has zero second derivative, giving
`sigma_second=4+4eta`, consistently `s=2+2eta`.

For the present gauge pullback, the shell conditions are f=0 and
b_r=d_r=0, with tangential derivatives. The checker leaves the values b,d
and the normal derivative f_r free. In the scalar boundary row,
`(p*f)_y=p*f_y` cancels the additional normal-lapse contribution
`−p*(h_yy)_gauge/2`. The metric normal rows likewise cancel the normal
derivative of the `a*f` tangential perturbation against the lapse term.
The mixed shell metric values vanish, and the remaining offdiagonal and
anisotropic normal rows vanish. The code correctly tests that nonzero
coefficients of free b,d,f_r remain in the reconstructed fields before
forming the physical residuals. No gauge-wave PDE is used in these
boundary covariance identities.

## Independent replay and scope

The checker was frozen at SHA256
`e6a9f67038e241c1081e0ef916ebe3482b028b5ed4950b876577eb4a1c409623`
and remained byte-identical through this review's runs. Fresh commands
`python -I -S -B` and `python -I -S -B -O`, each with `--controls`,
passed all 33 identity groups and 12 corrupted subprocesses. I independently
confirmed the exact expected first failed identity and nonzero coefficient
remainder for every mutation. Each failed subprocess leaves no report.
The two full reviewer reports are byte-identical, 11,715 bytes, SHA256
`c0dc60926573d0181d03c4744c25d031dd9cd294da44ee4ca3a9fb5cf081dd59`.
The machine receipt is `SHELL_REVIEW_RECEIPT.json`; reviewer outputs are
`SHELL_FROZEN_REVIEW_NORMAL.json` and
`SHELL_FROZEN_REVIEW_OPTIMIZED.json`.

The older exact physical junction derivation and the smooth actual
master/dynamic-boundary/gauge solutions remain analytic inputs. This is a
finite implementation-to-proof coefficient join. It supplies no new
existence theorem, nonlinear family, second-order evolution result,
physical approximation bound, or cosmological-origin conclusion.
