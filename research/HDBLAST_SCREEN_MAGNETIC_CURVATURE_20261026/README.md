# Conditional registered screen magnetic curvature

This checkpoint evaluates a selected angular projection of the second parameter derivative of five-dimensional Weyl curvature. It combines the previously bounded tensor response with the complete quadratic first-metric contribution, including the observer and screen corrections. Under the retained premises, the selected coefficient is bounded strictly away from zero on the exterior angular sphere labelled by `(x,r)=(1/4,5/4)`, before shell or reflected contact.

This is a formal curvature-response result for specified registered observers. The compatible full second-order Einstein–matter boundary problem remains unfinished. No finite-amplitude family of exact nonlinear solutions, gauge-invariant measurement, gravitational radiation, force, detector response or higher-dimensional origin of the Big Bang is established.

## Defined quantity and observer

For the future unit normal `U` to the registered `t=constant` slices, let `n` be the positive unit spatial normal to `y=constant` **within those slices**. The screen consists of their three angular tangent directions. Use the outward sphere orientation induced by ambient `EPS[0,1,2,3]=+1` and define

\[
H_{ab}=\tfrac12\epsilon_{(a}{}^{cd}C_{b)Ucd}.
\]

Parentheses denote normalized symmetrization. This three-screen symmetric tensor is a projection of five-dimensional Weyl curvature. It is not the full magnetic Weyl tensor on the observer's four-dimensional rest space. The normal is not the unnormalized coordinate vector `partial_y`, nor the five-dimensional gradient of `y` without spatial projection.

Use `g(lambda)=g0+lambda*q+lambda^2*k/2+...`; every displayed `H2` is the second parameter derivative, not the coefficient of `lambda^2`. The background Weyl vanishes. The entire first variation of this screen tensor vanishes for the registered scalar first field, without a matter-equation shortcut. Nevertheless the full first Weyl tensor need not vanish, so its coupling to the first observer/screen variations must be retained at second order.

No geometric clock, observer or point transport has been constructed. Consequently the vanishing of `H0,H1` is not used to claim second-order gauge invariance for this coordinate prescription. Different observer prescriptions can change `H2`.

## Complete selected coefficient

Let `R=rho*cosh(t)`, `theta=tanh(t)`, `t=ell*x`, and let `A,B,C,D,E,L,M` be the parent first-metric amplitudes in its registered frame. Subscripts `t,y` below are ordinary physical coordinate derivatives. For the frozen scalar seed `Y=n0^2-n1^2`, exact all-sphere identities give

\[
\operatorname{curl}T_+=3T_+,\qquad
\operatorname{curl}T_-=-3T_-,\qquad
\epsilon_{(a}{}^{cd}(D_cY)(D_{b)}D_dY)=2(T_+-T_-).
\]

With the equal TT response derivative `u_x` from the preceding capsule,

\[
\Pi_{\rm selected}H_2=h(T_+-T_-),\qquad
h=\frac{3u_x}{2\ell\rho^2\cosh t}+\frac{4F}{R},
\]

\[
F=\frac{(D+4M-L)M_t-2ML_t}{4\rho}
-\frac{E(D+L)}{2R}
+\left(\frac C R-\frac{M_y}{2}\right)
 \left[\frac{R(C_t-\theta C)}{2\rho}-\frac{RE_y}{2}+\frac B2\right].
\]

Here the selected projection is defined by amplitudes `6 mean(H:Tplus)` and `6 mean(H:Tminus)`, with normalized sphere measure `mean(1)=1`. It must not be identified with an arbitrary full TT harmonic space. The scalar quadratic tensor itself closes exactly on this pair. The background linear operator contributes no scalar/vector second-metric terms to this projection; the frozen angular curl map supplies its TT normalization and sign.

The derivation uses the exact ADM/Codazzi identity. Symmetric screen dualization cancels all Schouten/Ricci subtractions. The calculation retains the lapse and shift, radial shape, angular connection, and screen-volume/raising corrections. The radial amplitude `A` and background radial expansion cancel only after these terms are combined. The final map has 14 quadratic Laurent monomials; it is not obtained by dropping individually small terms.

## Certified interval and meaning

The endpoints below are outward decimal displays of exact rational enclosures, retaining all authenticated actual first-field errors:

| Quantity | Enclosure |
|---|---|
| Full `Tplus` coefficient `h` | `[-0.024929793753, -0.004128493946]` |
| Full `Tminus` coefficient `-h` | `[0.004128493946, 0.024929793753]` |
| Complete scalar quadratic contribution `4F/R` | `[-0.000000000000000000690462, 0.000000000000000001996501]` |

The quadratic contribution has absolute value below `2e-18`, and its ratio to the absolute linear term is below `5e-16`, in this registered normalization. It cannot cancel the signed linear contribution at this sphere. This conclusion concerns the **complete** projected coefficient, rather than the previously isolated linear term.

Exact angular normalization gives `mean |Pi_selected H2|^2=h^2/3>0`. Pointwise, its squared norm is `2*h^2*S*(1-S)`, where `S=n0^2+n1^2`. It vanishes on the two great circles `S=0` and `S=1`; nonzero angular projection does not mean nonzero at every angular location. Reversing orientation reverses the signed coefficients and preserves the nonvanishing conclusion.

The correction is a second parameter derivative about a model background. Its numerical value is not an observed field strength or energy. Any physical realization additionally requires the missing nonlinear/constraint/boundary and observer-anchoring work.

## Verification and reproduction

Read `theory/DERIVATION.md`, `angular/DERIVATION.md`, and `review/INDEPENDENT_SCREEN_AUDIT.md` for the mathematical arguments. Exact quotient-polynomial and separately written homogeneous-polynomial calculations verify the sphere identities at all angular points. A separate standard-library coordinate program computes all five-dimensional connection, Ricci, scalar and Weyl components for generic rational fixtures, including a combined `q+kTT` fixture. It checks the full second-variation normalization and detects omission of observer or epsilon corrections. Those finite coordinate fixtures support the analytic derivation; they are not an all-spacetime proof by sampling.

The numeric exporter imports the retained actual C3 error bounds and the previously certified `u_x` interval. It authenticates ten complete frozen capsules before importing their project code and after execution. No older source or response campaign is rerun. Normal and optimized isolated Python runs, including relocated source copies, give identical canonical results. `check_final_coefficient.py` also proves exact Laurent-polynomial equality between the geometric implementation and the encoded numerical map.

`MANIFEST.json` authenticates every source and reference file except its own bytes. The companion ZIP contains exactly this capsule, with identical members. Obtain its manifest and ZIP hashes from the external publication record. Then run from a checkout containing the ten upstream capsules:

```sh
python -I -S -B research/HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026/reproduce.py \
  --repo-root /path/to/HDblast \
  --manifest-sha256 MANIFEST_SHA256 \
  --archive /path/to/HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026.zip \
  --archive-sha256 ZIP_SHA256 \
  --output-dir /fresh/external/replay-directory
```

The parent directory of the fresh output must exist. `REPLAY_PLAN.json` lists every direct run, expected wrong-formula rejection, reference byte comparison and Python program role. The entry authenticates all complete input/source inventories and ZIP bytes before executing project code, and reauthenticates them afterward. Raw execution logs are generated only in external replay output; they are not scientific source files. Runtime requires only the Python standard library.

Internal AI-assisted reviews and computational certificates do not replace independent specialist review, proof-assistant formalization, or assessment of novelty and significance. See the bounded primary-source comparisons in `literature/`; downloaded paper bodies and private editorial history are excluded.

No new Zenodo version or journal submission accompanies this checkpoint. The next scientific gate is to specify and certify a genuinely geometric observer/clock/point prescription and its transport corrections, while continuing the full compatible second-order boundary problem.
