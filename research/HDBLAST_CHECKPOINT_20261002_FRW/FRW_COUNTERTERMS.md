# Positive-reference curved-background subtraction

This is a candidate construction on a sufficiently smooth spatially flat FRW background, not an executed absolute source on the archived shell. It extends the previous flat matching without a denominator that vanishes at x=M²=0.

Let u=a chi and use conformal time. The exact physical mode equation is
\[
u_k''+[k^2+a^2x-a''/a]u_k=0,\qquad uu'^*-u'u^*=i.
\]
Choose a fixed positive reference r and define
\[
\omega^2=k^2+a^2r,\quad h=x-r,\quad
\Delta=a^2h-a''/a,\quad {\cal H}=a'/a.
\]
Count h as adiabatic order 2 and each derivative as +1. This is subtraction bookkeeping; it does not assume the physical crossing is slow.

The formal WKB frequency has
\[
W=\omega+W_2+W_4,\quad
W_2=\frac{\Delta}{2\omega}-\frac{\omega''}{4\omega^2}
+\frac{3\omega'^2}{8\omega^3},
\]
\[
W_4=-\frac{W_2^2}{2\omega}-\frac{W_2''}{4\omega^2}
+\frac{\omega''W_2}{4\omega^3}
+\frac{3\omega'W_2'}{4\omega^3}
-\frac{3\omega'^2W_2}{4\omega^4}.
\]
For this WKB mode the per-mode stress expressions are
\[
\rho_W=\frac1{4a^4}\left[
W+\frac{k^2+a^2x+{\cal H}^2}{W}
+\frac{{\cal H}W'}{W^2}+\frac{W'^2}{4W^3}\right],
\]
\[
p_W=\frac1{4a^4}\left[
W+\frac{{\cal H}^2-k^2/3-a^2x}{W}
+\frac{{\cal H}W'}{W^2}+\frac{W'^2}{4W^3}\right].
\]
Expand the complete expressions consistently through total order four. Do not simply substitute a truncated W into unexpanded rational functions. The field-square subtraction is
\[
Q_{\rm sub}=\frac1{2a^2}\left(\frac1\omega-\frac{W_2}{\omega^2}\right).
\]
The paired mass-squared source is S=Q/2. Exact physical modes, not their WKB approximation, supply the bare observables.

## Exchange identity

The complete formal WKB expressions satisfy the mode energy-exchange identity up to the mode-equation residual. Matching successive WKB orders removes that residual through the necessary orders. Sorting by grade leaves
\[
\rho_0'+3{\cal H}(\rho_0+p_0)=0,
\]
\[
\rho_2'+3{\cal H}(\rho_2+p_2)=h'Q_0/2,\qquad
\rho_4'+3{\cal H}(\rho_4+p_4)=h'Q_2/2.
\]
There are no higher-grade terms in the derivative of these truncated expressions. Thus their sum supplies the corresponding subtraction identity. code/verify_frw_counterterms.py evaluates these identities on 24 exact rational local Taylor jets, including x=0 cases, and checks the flat reduction on 12 further jets. Those selected exact checks are not a universal symbolic proof, field evolution or convergence test. Deliberate pressure/current omissions must be detected.

## Flat reduction and remaining matching

For a=1 and omega=sqrt(k²+r), the subtractions reduce to
\[
\rho_{\rm sub}=\omega/2+h/(4\omega)-h^2/(16\omega^3),
\]
\[
p_{\rm sub}=\frac{k^2}{6\omega}-\frac{hk^2}{12\omega^3}
+\frac{h^2k^2}{16\omega^5}
-\frac{h''(2k^2/3+r)}{16\omega^5},
\quad Q_{\rm sub}=\frac1{2\omega}-\frac h{4\omega^3}.
\]
Subtracting these reproduces the preceding direct flat prescription, including its curvature-pressure term.

Flat agreement fixes the displayed potential and derivative-dependent pressure convention, but cannot determine a constant F(r)R coefficient: its flat variation vanishes. The earlier declared F(r)=0 must remain an explicit curved-action matching condition, and the complete scheme mapping must be established. Likewise finite curvature-squared coefficients were not fixed by the flat calculation. A common covariant action must specify any finite deltaF R, alpha R² and beta C² terms before a physical curved-source claim. Ward consistency and a flat limit do not remove those choices.

Positive r makes subtraction denominators regular at x=0 for a>0. It does not guarantee that W0+W2+W4 is positive at low momentum. Do not use that truncated series blindly as the all-k initial vacuum. State preparation and a positive low-momentum completion require their own specification.

## A controlled next experiment

Before using a smoothed archive, a useful fully specified test has an exactly static past and future. Let f(s)=exp(-1/s) for s>0 and zero otherwise, and B(s)=f(s)/[f(s)+f(1-s)]. Define
\[
a(\eta)=\exp[A B(\eta/T)],\qquad
x(\eta)=r[2B(\eta/T)-1]^2.
\]
Both profiles are C-infinity and exactly constant outside [0,T], while x crosses zero. Exact static incoming data eliminate the finite-time vacuum ambiguity. This is a proposed new prescribed-background control, not an archived HDBLAST solution. Full stress/current integration, momentum-cutoff convergence and curved finite matching remain to be executed.

