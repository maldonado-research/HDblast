# Independent Sol acceptance of the quantitative response enclosure

15 October 2026 UTC. Read-only review of the primary result and arithmetic.
The earlier Sol analytic proof, audit and manifest have not been edited.

**Decision: accept the stated quantitative interval as a conditional analytic
result for the newly specified state.** At
\(\delta=1/1000\), \(\ell=10^{-10}\) and \(T_1=3\ell\), the reviewed argument
and exact arithmetic support

\[
\boxed{10^{15}<\mathcal C_{T_1}<4\times10^{17}.}
\]

This acceptance is conditional on the retained analytic regular-branch,
localized exact-data, completed-gauge/uniqueness and finite-time exact-family
hypotheses. It accepts a quantitative interval for the leading physical
coefficient, rather than merely a proposed route to such an interval. It does
not accept a computed physical waveform, a user-supplied dimensional scale, or
an effective nonzero range of the nonlinear parameter \(\varepsilon\).
Finite arithmetic programs do not formally certify the analytic hypotheses or
the PDE proof.

## Findings checked independently

1. **Physical output and normalization.** The inherited coupled reconstruction
   gives \(S=3\cosh^{-3/2}t\,D/(\rho_b^3p_b)\), with
   \(D=w_z-\gamma w_b\), \(\gamma=C_b/2\). Thus the proper-time/volume weight
   cancels exactly and
   \(\mathcal C=3\pi^2\|D\|_2^2/(\rho_b^2p_b^2)\). The approximation is to
   the coupled master with its retained physical gauge/embedding map. No test
   scalar replaces that map. The selected dimensionless \(t,z,\ell\), and
   normalized initial tangent in the model units, supply no observational
   physical length or time scale.

2. **Coherent shell coefficients and residual signs.** The branch uses
   \(\beta=\mu(35\operatorname{sech}^2t/4-25/4)+\gamma\),
   \(\mu=(\rho_bB_b)^{-1}>0\). For \(e=w-u\), the free odd-extension
   comparator has \(u_b=u_{tt,b}=0\), so its raw boundary residual is
   \(-\widehat D\). Accordingly the true error equation is
   \(-e_z+\mu e_{tt}+\beta e=G\), **\(G=+\widehat D\)**. The bulk source is
   \(F=-Wu\). This agrees with both the energy source \(+Gv\) and the
   characteristic ODE \(\mu v'+v+\beta q=J_++G\).

3. **Energy and derivative trace.** Integration by parts gives the stated
   nonpositive moving-edge flux
   \(-[(e_t-e_z)^2+\lambda e^2]/2\). Positive bulk and shell shifts bound
   both cross terms with the stated \(a\), and combined Cauchy–Schwarz bounds
   the forcing by \(\sqrt{f^2+g^2/m}\,X\). The additional characteristic
   identity, together with Minkowski and \(z=t-s\), controls the normal
   derivative. The proof does not assume that bulk energy alone bounds that
   derivative and does not silently omit a boundary acceleration residual.

4. **Comparator bounds.** The smooth odd extensions give the exact data after
   conjugation and \(\widehat D=q'+e^{-t}q\). The bound
   \(|q'|\le16/\ell\), support width \(\ell\), and rise from zero to one
   over \(\ell/2\) justify both norm bounds used in the primary note. The
   stronger full-bump lower bound is unnecessary for the claimed interval.
   The comparator is not asserted to solve the dynamic shell boundary law.

5. **Actual branch/collar input.** Starting from the pinned analytic theorem's
   endpoint inequalities, a separate Fraction-only calculation reconstructs
   adequate shell bounds and all four strict first-exit inclusions. It
   independently recovers \(\mu>1/1000\), \(\mu<1/400\), \(|\gamma|<250\),
   \(|\beta|<1000\), the squared prefactor bounds and
   \(|W|<10^6\) on \(0\le z\le8\ell\). It uses the near-\(\phi=+1\)
   analytic branch, its affine detuning, and \(\sigma''=4+4\eta\). No
   floating historical trajectory or S85 alternative supplies an enclosure.

6. **Whole correction support and strict time margin.** The energy-domain
   length \(L=8\ell\) is explicitly distinguished from the outer endpoint of
   the nonlinear correction annulus. The conditional geometry
   \(A_0=\ell/2\), \(L_{\rm corr}=3\ell\),
   \(T_{\rm evol}=7\ell/2\), \(c_*=5/4\), \(Z=8\ell\) obeys
   \(L_{\rm corr}+c_*T_{\rm evol}=59\ell/8<Z\) and
   \(T_1<T_{\rm evol}\). It requires the retained inverse/localization
   hypotheses for that full annulus and sufficiently small exact-family
   parameters to ensure the speed bound. No effective value of that smallness
   is inferred from the linear certificate.

## Independent coarse arithmetic

The new checker imports neither primary production program. It uses only
exact rational arithmetic and elementary analytic comparisons. A simpler
outward chain than the primary program is sufficient:

\[
f<29,\quad h<5060000000100,\quad
e^{aT/2}<1001/1000,\quad B_*<1520,
\]
\[
\sqrt{T(1000+250^2)}<7/1600,\qquad
\Delta<\frac{71060000087}{10000000000}<8.
\]

Machin's identity with eight arctangent terms suffices for
\(80000<K<125000\). The comparator/error triangle then gives
\(140000<\|D\|_2<1600010\). Consequently the separate calculation proves

\[
1568000000000000<\mathcal C_{T_1}<320004000012500000,
\]

which implies the primary note's displayed interval. This corroborates the
arithmetic using coarser constants rather than copying its large rational
endpoints. The primary normal and optimized receipts match exactly and pass
25 finite groups. The separate review passes 35 checks in both modes with
identical receipts.

## Accepted pins

All SHA-256 values below identify the reviewed bytes. Subsequent substantive
changes require a new review.

| Artifact | SHA-256 |
| --- | --- |
| `SHELL_RESPONSE_ERROR_CERTIFICATE.md` | `5f01be536c77123e9f667f4cd5514b52067641b2b389ba4e6eb655cb24654022` |
| `certify_response.py` | `f634507663644e204d65d12c996fce2b160fee36331d066f654ed5d44f37f386` |
| `check_response_certificate.py` | `953c8a51b21ed7c2cc3c51517a6df6d96ae0166c6f6a3990a6790fd68b7d3d36` |
| `archive/ACTUAL_BACKGROUND_AND_RIGOROUS_COLLAR_BOUNDS.md` | `0832a16c4bdf18cc188861a2c2860e7e6ff736dec39426d6ca7c589dc0434d5c` |
| `archive/ANALYTIC_SHELL_BOX_DELTA_1E3.json` | `8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8` |
| `archive/derive_shell_box.py` | `3a4d2a71847114b8dce9f910679a1f4ce14898212001b0355279ec633ebed110` |
| Retained Jacobi endpoint theorem | `f1e2571a551ec122dd7755a13145bae73bd92bd7f7426bdc3002e2f52f8dbfc8` |
| Retained exact rescaled certificate | `22c6cd4e505ee2c7e549420355b94a5b3c079ebb6945314f5b285c13bb67835d` |
| Original Sol analytic note, unchanged | `b8c77c810c70aa286227c5dabe010a53942af84067eb03961f0186419a5ea96a` |

Machine-readable pins and exact independent bounds are in
`INDEPENDENT_QUANTITATIVE_REVIEW.json`. The evaluator is
`independent_quantitative_review.py`; it is a new cross-review artifact and
does not alter the previously pinned Sol proof files.
