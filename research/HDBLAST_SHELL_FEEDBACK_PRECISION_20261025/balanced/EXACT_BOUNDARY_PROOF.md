# Fixed-budget exact-boundary repair experiment

This experiment evaluates a new continuous candidate while preserving the
frozen bulk feedback and velocity table. It establishes no sign for the exact
response or an intrinsic observable.

The new candidate is

\[
 V_c=B_0+K\widehat f_{\rm old}
       +L\{M_*P_N[u_0+P v_{\rm old}](x)\},
 \qquad B_0=\operatorname{free}(u_0,v_{\rm old}).
\]

Here \(u_0\) is the exact registered rational displacement, \(v_{\rm old}\)
is the frozen continuous P1 velocity, and \(P\) denotes integration from
zero. The frozen candidate instead uses its P1 displacement in the Robin
forcing. No change is made to the frozen capsules, their profile tables,
or their published error claims.

## Balanced shell evaluation

Put \(f=u_{0,J}\), \(m=(M_*)_{Jw}\), and

\[
 C(r)=\frac{3}{2\ell\rho_{\rm nom}(r)p_{\rm nom}(r)}.
\]

The registered displacement formula gives the exact nominal identity
\(u_{0,w}=C f_r\). It uses derivatives of the registered nominal
polynomials; no exact-background ODE identity is substituted. Both \(C\)
and \(m\) are negative. Since \(f\) vanishes on \([0,65/64]\),

\[
 f(x)-m\int_0^x u_{0,w}(r)\,dr
   =[1-mC(x)]f(x)+m\int_{65/64}^x C'(r)f(r)\,dr.
\]

The only other nonzero entry in the J row of \(M_*\) is \((M_*)_{JK}\),
and \(u_{0,K}=0\). The shell lift has the sign
\(Lg(x,0)=-\int_0^xg(s)\,ds\). Therefore the complete J component of
the new candidate is exactly

\[
 (V_c)_J(x,0)=[1-mC(x)]f(x)
 +m\int_{65/64}^x C'f
 +(Pv_{\rm old})_J(x)+(K\widehat f_{\rm old})_J(x,0)
 -(M_*)_{JK}(P^2v_{\rm old})_K(x)
 -m(P^2v_{\rm old})_w(x).
\]

The velocity primitives and bulk term are evaluated by the authenticated
frozen rational ridge helper. For the endpoint factor, the numerator of

\[
 1-mC=\frac{2\ell\rho_{\rm nom}p_{\rm nom}-3m}
             {2\ell\rho_{\rm nom}p_{\rm nom}}
\]

is formed as exact polynomial coefficients before interval evaluation.

## Fixed quadrature and outward arithmetic

The complete budget was declared before evaluation: 256 subintervals,
equivalently 128 Simpson panels. They are split into 16 subintervals on
\([65/64,66/64]\), 112 on \([66/64,9/8]\), and 128 on
\([9/8,5/4]\). This gives 257 distinct point evaluations and 128 whole-panel
derivative evaluations. There is no adaptive refinement or mesh escalation.

Let \(h(r)=mC'(r)f(r)\). On every Simpson panel \([a,b]\),

\[
 \left|\int_a^b h(r)\,dr
 -\frac{b-a}{6}\{h(a)+4h((a+b)/2)+h(b)\}\right|
 \le\frac{(b-a)^5}{2880}\sup_{[a,b]}|h^{(4)}|.
\]

The frozen outward interval implementation uses exact fractions on the
\(2^{-112}\) grid. The unchanged registered expressions for \(f\) are
extracted by AST from the authenticated producer, stopping before its unused
later calculations. Order-two input jets supply point values. Order-six
input jets supply the four derivatives of \(f\) required by the remainder.
The jets differentiate in \(z=\ell r\), so \(4!\ell^4c_4\) encloses the
ordinary fourth r derivative of \(h\). The polynomial cutoff at \(66/64\)
is a segment endpoint; each panel uses one analytic branch. One-sided
derivative enclosures on adjoining branches suffice for Simpson's theorem.

All seven frozen capsule closures are authenticated before project imports
and again after the computation. The report contains every panel's rational
Simpson interval, fourth derivative interval, and remainder bound.

## Results and limits

The following integer intervals round the reported rational candidate
intervals outward:

| Shell time x | New exact-boundary candidate J interval | Simpson remainder upper |
| --- | --- | --- |
| 9/8 | [−49087477769334561, −49087475281526126] | 1243904217 |
| 5/4 | [−197742099538614228, −197742091652732819] | 3942940704 |

Thus the new candidate is negative at these two points. The candidate-only
change from the frozen P1 boundary evaluation is approximately
\(-2.8450589277\times10^{18}\) at \(9/8\) and
\(-1.3087437399\times10^{18}\) at \(5/4\).

The separate component comparison bounds available at completion,
approximately \(1.033806\times10^{18}\) and \(1.411113\times10^{19}\),
exceed these candidate magnitudes. Adding those total-response errors and
the quadrature uncertainty leaves intervals containing zero. The repaired
candidate sign therefore establishes neither the exact J sign nor an
observable sign. No nonlinear family or physical-origin inference follows.

Reproduction:

```sh
python exact_boundary.py --repo-root ../repo --output NEW_REPORT.json
```

The output path must be new. This runs only the bounded new quadrature;
it imports authenticated helpers but does not execute the old campaigns.

Artifact pins at completion:

- `exact_boundary.py`: `c01a8f5f95737dcacebdc26379ae0ee5863af9c1977890b3901ad6113103e619`
- `EXACT_BOUNDARY_NORMAL.json`: `36cdd588603c2c15226a6793430d8ed3ff3f2101e94afca1f8744a73fce9e276`

An independent reviewer checked the identity, AST extraction, derivative
scaling, branch splitting and scope, and replayed every stored Simpson
remainder and cumulative interval assembly using exact rational arithmetic.
