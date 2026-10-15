# Validated quadratic material tension contribution at the shell

15 October 2026. Conditional mathematics for the registered five-dimensional
Einstein–scalar shell model. AI-assisted derivation and independent internal
checks; no external peer review or novelty claim.

## Result and relation to the second-order problem

For the previously certified wider quadrupole state, this checkpoint encloses
the **explicit tension nonlinearity** in the physical material scalar junction,
its angular components, its complete proper-time integrals, and its small
postpulse tail. It also encloses the corresponding constitutive tension term.
The independently derived fourth-moment integral is new numerical information
relative to the retained first-variation norm certificate.

The earlier full-contact checkpoint already constructs the complete second
Taylor coefficient conditionally, including scalar degrees 0, 2, 4 and the
degree-2 tensor sector, under its stated initial-data and evolution assumptions.
That construction is preserved. This calculation supplies quantitative input
for its boundary rows at the wider state; it does not evaluate its full `U2`.

In particular, the explicit material contribution below is **not the complete
right side of the equations for second-order bulk-coordinate unknowns**.
Normal normalization, metric, embedding and pullback contributions enter when
material jets are expressed using those unknowns. Section 3 gives the actual
coordinate counterpart, already derived in the retained boundary system.

## 1. Registered state, derivative convention and physical clock

Use the same near-\(\phi=+1\) background, \(\delta=1/1000\),
\(c=96714321425692/161839258930731\), and outward normals on both copies.
Write \(\eta=\phi-1\). The actual common tension is

\[
\sigma(\phi)=\frac23+2\eta^2+\frac23\eta^3
 +\delta(1+c+c\eta),\qquad
\sigma''_b=4+4\eta_b>0,\quad \sigma'''=4.
\tag{1}
\]

It is not the tension from the older alternative model. The authenticated
shell box gives
\(3.999651068<\sigma''_b<3.999676637\).
All second coefficients are derivatives:

\[
U_\varepsilon=U_0+\varepsilon U_1+\tfrac12\varepsilon^2 U_2+\cdots.
\tag{2}
\]

For \(\ell=10^{-6}\), set \(t=\tau/\rho_b\),
\(x=t/\ell\), and retain the full window \(0\le x\le3\).
The quadrupole is \(Y=n_1^2-n_2^2\) on the unit three-sphere. Its angular
matrices commute, \(A=B=\operatorname{diag}(1,-1,0,0)\).
The registered initial displacement is \(Q(z/\ell)\) and velocity is
\(e^{-z}Q(z/\ell)\) for the normalized master field, where
\(Q(x)=\exp[4-1/((x-1)(2-x))]\) on \((1,2)\), flat zero otherwise.
This is the coupled metric–scalar perturbation retained in the prior proof.

The physical shell scalar first variation is

\[
\psi_1(\tau,\Omega)=S(t)Y,\qquad
E(x):=\ell S(\ell x)=G\cosh^{-3/2}(\ell x)H(x),\quad
G=\frac3{\rho_b^3p_b}<0.
\tag{3}
\]

Here \(H(x)=\ell(w_z-\gamma w_b)(\ell x,0)\) belongs to the **actual
coupled first variation**, with its potential and time-dependent boundary
residuals enclosed. The prior file supplies 4353 continuous cells: one before
contact, 4096 over the pulse, and 256 over \([2,3]\). Every old cell is
authenticated and retained byte for byte. Numerical values below use model
units; no observed duration or amplitude calibration is assigned.

## 2. Material scalar row, metric row and tension trace

The exact physical junctions in the outward-normal SUM convention are

\[
\sum_iK_{i,ab}-\frac{\sigma(\psi)}3h_{ab}=0,\qquad
\sum_ip_i+\sigma'(\psi)=0.
\tag{4}
\]

Twice differentiating the scalar law gives

\[
\sum_ip_{i,2}+\sigma''_b\psi_2
 =Q_\phi,\qquad \boxed{Q_\phi=-4\psi_1^2=-4S(t)^2Y^2.}
\tag{5}
\]

For identical copies the entire row is divided by two, so the right side is
\(-2\psi_1^2\). The quantities \(p_{i,2}\) and \(\psi_2\) are full
material derivatives; they already contain second normal, embedding, clock
and pullback effects. Equation (5) isolates a constitutive nonlinearity rather
than eliminating these effects.

For completeness the full material metric row is

\[
\sum_iK_{i,2,ab}-\frac{\sigma_b}3h_{2,ab}
 -\frac{\sigma'_b}3\psi_2h_{0,ab}
 =\frac{2\sigma'_b}3\psi_1h_{1,ab}
  +\frac{\sigma''_b}3\psi_1^2h_{0,ab}.
\tag{6}
\]

For a pure tension stress \(\mathcal S_{ab}=-\sigma(\psi)h_{ab}\), the
physical mixed trace satisfies
\(\Theta=-\mathcal S^a{}_a/4=\sigma(\psi)\), hence

\[
\Theta_2=\sigma'_b\psi_2+T_{\rm quad},\qquad
T_{\rm quad}=\sigma''_b\psi_1^2.
\tag{7}
\]

The positive sign of \(T_{\rm quad}\) fixes neither the total \(\Theta_2\)
nor the sign of a second-order metric or cosmological backreaction. Lower-index
stress additionally contains \(-2\sigma'_b\psi_1h_1-\sigma_bh_2\).
Since \(\sigma''_b=4\phi_b\), these two explicit contributions obey the
pointwise identity \(T_{\rm quad}=-\phi_bQ_\phi\).
Integrals here use the **background** shell measure, so they are not complete
second variations of an integral with perturbed volume or moving endpoints.

Because \(\psi_0\) and \(p_{i,0}=p_b=-\sigma'_b/2\) are constant on the shell, a first clock/label relabeling
does not change \(\psi_1\) or \(Q_\phi\). It changes the second jets by
\(\psi_2\mapsto\psi_2+2\mathcal L_{\xi_1}\psi_1\) and
\(p_{i,2}\mapsto p_{i,2}+2\mathcal L_{\xi_1}p_{i,1}\).
The change of the left side of (5) is twice the Lie derivative of the
vanishing first junction row. This covariance does not make the individual
second jets clock independent.

## 3. Explicit distinction from the coordinate source

In the retained longitudinal first representative,
\(h_{1,yy}=-4uY\), \(h_{1,ya}=0\),
\(h_{1,ab}=2uYh_{0,ab}\), \(\chi_1=vY\), with zero first shell
displacement. Thus the shell trace \(v_b=S(t)\) in this representative and clock, but this
identification must not be silently transferred to another gauge.
Set \(a=\rho_y/\rho\), \(p=\phi_{0,y}\).
The one-copy first shell rows are

\[
u_y+2au+pv/3=0,\qquad
v_y+2pu+\sigma''_bv/2=0.
\tag{8}
\]

Direct differentiation of the unit normal gives
\(D^2p=(4uv_y+12pu^2)Y^2\) and
\(D^2K_{ab}=(4uu_y+20au^2)Y^2h_{0,ab}\).
After moving these terms to the right and using (8), the complete scalar
and isotropic metric source coefficients in this original representative are

\[
\boxed{S_p=-4pu^2+2\sigma''_buv-2v^2,\qquad
S_K=\frac{\sigma''_b}6v^2-12au^2.}
\tag{9}
\]

These equations, including their normal terms, are independently reproduced
by the archive checker. This checkpoint encloses the material contribution
\(-2v^2\), **not the whole numerical function \(S_p\)**. Evolving `U2`
also needs the complete bulk sources, all shell/constraint rows and specified
second initial data. Their prior conditional construction supplies a framework;
a numerical response with certified residuals remains to be evaluated.

## 4. Exact angular decomposition and power

Define
\(J=n_1^2+n_2^2\), \(Z_2=J-1/2\), and
\(H_4=Y^2-J/2+1/12\). Then

\[
Y^2=\frac16+\frac12Z_2+H_4,\qquad
\Delta Z_2=-8Z_2,\quad\Delta H_4=-24H_4.
\tag{10}
\]

The coefficient functions of \(Q_\phi\) in this precise basis are
\(-2S^2/3,-2S^2,-4S^2\). Squaring a degree-2 perturbation therefore forces
a nonzero scalar monopole as well as degrees 2 and 4. This angular closure
was already known in the retained construction; here it is evaluated with
the actual certified waveform.

With normalized spherical averages,
\(\langle Y^4\rangle=3/40\),
\(\langle Z_2^2\rangle=1/12\), and
\(\langle H_4^2\rangle=19/720\).
The mutually orthogonal squared angular powers of the three projections,
divided by their sum, are respectively

\[
\boxed{\frac{10}{27},\quad\frac5{18},\quad\frac{19}{54}.}
\tag{11}
\]

These fractions concern this source and this angular state; they are not
fractions of a nonlinear solution's energy or observable cosmic anisotropy.

## 5. Proper-time source integrals and fourth moment

The background spacetime measure on the shell is

\[
d\mu_0=\rho_b^4\cosh^3t\,dt\,d\Omega,
\qquad K=\frac{3\pi^2}{\rho_b^2p_b^2}.
\tag{12}
\]

The previous response certificate encloses
\(C=(K/\ell)\int_0^3H^2dx\), with
\(6.06\times10^{11}<C<7.44\times10^{11}\).
Exact sign and normalization then give

\[
\int Q_\phi\,d\mu_0=-4C,\qquad
\|Q_\phi\|_{L^1(d\mu_0)}=4C,\qquad
\int T_{\rm quad}\,d\mu_0=\sigma''_b C.
\tag{13}
\]

The source's squared \(L^2\) norm requires the independently computed
fourth moment, rather than just \(C\):

\[
I_4=\int_0^3\cosh^{-3}(\ell x)H(x)^4dx,
\qquad
\|Q_\phi\|_2^2=\frac{36KG^2}{5\ell^3}I_4,
\quad
\|T_{\rm quad}\|_2^2=\frac{9(\sigma''_b)^2KG^2}{20\ell^3}I_4.
\tag{14}
\]

For example, \(\int_{S^3}Y^4d\Omega=3\pi^2/20\) while
\(\int Y^2d\Omega=\pi^2/3\); their ratio is \(9/20\).
This and \(dt=\ell dx\) explain both the coefficient and inverse-cubic
width in (14). Exact outward rational arithmetic establishes the strict
rounded bounds:

| Quantity, full window \(0\le x\le3\) | Certified strict lower bound | Certified strict upper bound |
|---|---:|---:|
| \(I_4\) | 86.46 | 87.43 |
| \(\int Q_\phi d\mu_0\) | \(-2.972\times10^{12}\) | \(-2.427\times10^{12}\) |
| \(\|Q_\phi\|_1\) | \(2.427\times10^{12}\) | \(2.972\times10^{12}\) |
| \(\|Q_\phi\|_2^2\) | \(5.778\times10^{21}\) | \(8.644\times10^{21}\) |
| \(\int T_{\rm quad}d\mu_0\) | \(2.427\times10^{12}\) | \(2.971\times10^{12}\) |
| \(\|T_{\rm quad}\|_2^2\) | \(5.777\times10^{21}\) | \(8.642\times10^{21}\) |

Exact rational endpoints, rather than these rounded displays, are authoritative
in `computed/QUADRATIC_FORCING_ENCLOSURE.json`. Their apparent magnitude depends
on the registered narrow-pulse normalization and model units. It does not set
a finite perturbation amplitude or a physical energy scale.

## 6. Continuous tail and enclosure method

Throughout \(2\le x\le3\), the retained strict theorem gives
\(3.2\times10^{-8}<H(x)<5.3\times10^{-8}\).
Its sharper exact interval, the negative \(G\) interval and the positive
clock factor imply

\[
\boxed{-2.3\times10^{-19}<
 \ell^2\langle Q_\phi(\tau=\rho_b\ell x,\cdot)\rangle_{S^3}
 <-6.7\times10^{-20},\qquad 2\le x\le3.}
\tag{15}
\]

Thus the angular mean remains strictly negative over the entire tail. The
pointwise angular field is nonpositive and vanishes on the nodal set \(Y=0\);
there is no claim of strict negativity at every angle.

The new CSV carries all six scalar/tension coefficient intervals, multiplied
by \(\ell^2\), on all 4353 inherited continuous cells. Endpoint integers
are divided by \(2^{112}\); time endpoints by 4096. Powers of intervals
crossing zero have lower bound zero. The fourth-moment sums use exact rational
powers and the whole-cell factor
\([(1-(\ell x_{\rm right})^2)^2,1]\). No time samples stand in for
continuous enclosures.

For an elementary clock proof, the power series gives
\(\cosh t\le(1-t^2/2)^{-1}\) when \(0\le t^2<1\).
Convexity implies \((1-t^2/2)^{3/2}\ge1-3t^2/4\ge1-t^2>0\).
Consequently \(\cosh^{-3/2}t\ge1-t^2\), and squaring gives the factor
used above. Sol independently uses a sharper positive series bound and checks
that its resulting interval is contained in this enclosure.

## 7. Reproduction and next scientific step

The source includes authenticated minimal public inputs, independently derived
Astra/Sol/archive notes and checkers, an independent complete-cell reconstruction,
and a fresh normal/optimized replay with genuine wrong-formula controls. The
replay roster names each calculation and failure reason; its finite checks
validate arithmetic and identities, not the analytic theorem assumptions.

The next quantitative task is to reconstruct the remaining first-order metric
trace needed in (9), enclose the complete coordinate boundary and bulk source
functions, then evaluate the conditional second-order evolution with specified
second initial data and certify its residuals. A useful resulting observable is
the full physical monopole of \(\Theta_2\), including \(\sigma'_b\psi_2\)
and the proper clock/embedding transformations. An effective finite-amplitude
remainder would require a separate nonlinear estimate.

No exact nonlinear cosmological spacetime, higher-dimensional origin of the
Big Bang, observational evidence, priority claim or APS readiness follows from
this checkpoint. No new Zenodo version or APS submission accompanies it.
