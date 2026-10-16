# Independent review of the full second-response operator and norms

16 October 2026. This review accepts the exact reduced linear operator,
the finite angular Hilbert space, the frame conversion, and the quantitative
majorants in `evolution/FULL_SECOND_RESPONSE_CONSTRUCTION.md`, conditional on
its authenticated background, first solution, source and initial-data inputs.
It does not independently certify the qualitative regularity, subsidiary
boundary argument or existence proof; those have a separate fresh review.
This is an internal AI-assisted mathematical audit, not external peer review.

The accompanying standalone checker imports no producer code. Its direct
coordinate finite-jet computation is separate from the covariant recipe.
Its angular branching is computed from the diagonal-spin restriction rule;
its exact arithmetic joins read copied, SHA256-pinned actual reports. Finite
jet witnesses detect algebra mistakes but do not prove a differential identity
at every point. The following derivations supply the all-point arguments.

## 1. Curvature signs and the scalar derivative coupling

Use the convention

\[
 R^a{}_{bcd}=\partial_c\Gamma^a_{db}-\partial_d\Gamma^a_{cb}
 +\Gamma^a_{ce}\Gamma^e_{db}-\Gamma^a_{de}\Gamma^e_{cb},
 \qquad \mathrm{Ric}_{ab}=R^c{}_{acb}.
\]

For a covariant metric variation \(k\), differentiating the metric
compatibility identity gives

\[
 P^a{}_{bc}=\tfrac12g^{ad}
 (\nabla_bk_{cd}+\nabla_ck_{bd}-\nabla_dk_{bc}).
\]

Differentiating Ricci gives
\(\mathrm{Ric}'_{ab}=\nabla_cP^c{}_{ab}-\nabla_bP^c{}_{ac}\).
Commuting the two covariant derivatives on the two slots of \(k\) yields

\[
 \mathrm{Ric}'_{ab}=-\tfrac12\Box k_{ab}
 +\nabla_{(a}dC[k]_{b)}
 +\tfrac12(\mathrm{Ric}_a{}^ck_{cb}+\mathrm{Ric}_b{}^ck_{ca})
 -R_{acbd}k^{cd},
\qquad dC[k]_a=\nabla^bk_{ab}-\tfrac12\nabla_a\mathrm{tr}\,k.
\]

Parentheses contain one half. Subtracting the variation of the registered
relative-harmonic reduction removes exactly the displayed symmetric
derivative. There is no residual derivative of a background gauge vector,
which vanishes identically. The curvature operator is therefore exactly
\(L_R\) of construction (1), with its displayed minus Riemann term. As a
quick structural check \(L_R(g_0)=0\), as required by constant metric scaling.

Direct variation of the scalar wave operator is

\[
 D(\Box\phi)[k,w]=\Box w-k^{ab}\nabla_a\nabla_b\phi_0
                   -dC[k]^a\nabla_a\phi_0.
\]

This follows by differentiating both inverse metric and Christoffel symbol;
the contraction \(g^{bc}P^a{}_{bc}\) equals \(dC[k]^a\). Consequently
\(S'[k,w]=\Box w-k:\nabla^2\phi_0-dC[k]\cdot d\phi_0-V''w\).
Combining this with \(E'_{\rm red}[k,w]+Q_{\rm red}=0\) and
\(S'[k,w]+Q_\phi=0\) gives precisely both equations (1), including the
factor two on the metric source and the sign of the scalar source. A
second harmonic constraint is not permission to delete the \(dC[k]\)
coupling: it equals \(-Q_C\), rather than zero, after the subsidiary step.

The checker tests this derivation on two generic curved, nonnormal metric
two-jets in dimension five. Both first and second derivatives of the metric
are nonzero. Direct inverse-metric and connection coefficient differentiation
is compared against the covariant recipe with the derivative-slot connection
term included. There are 162 exact rational comparisons in total. The
curvature-sign, derivative-slot, symmetrization and scalar-\(dC\) mutations
fail these actual identities, rather than checking a metadata declaration.

## 2. The frame wave operator

Write \(D_a\) for the unit-sphere connection, including all angular tensor
slots. On frame components

\[
 \nabla_y=\partial_y,\qquad
 \nabla_0=\rho^{-1}\partial_t+\Omega_0,\qquad
 \nabla_a=(\rho\cosh t)^{-1}D_a+\Omega_a.
\]

For the rough wave operator one must also subtract the connection in the
derivative slot: \(\Box=\sum_A\epsilon_A
(\nabla_A\nabla_A-\nabla_{\nabla_Ae_A})\). Its radial and temporal terms are

\[
 \partial_y^2-\rho^{-2}\partial_t^2
 -2\rho^{-1}\Omega_0\partial_t-\Omega_0^2+a\partial_y.
\]

Here \(\partial_t\Omega_0=0\). Summing the angular terms gives

\[
 (\rho c)^{-2}\Delta_\kappa+rac2{\rho c}\sum_a\Omega_aD_a
 +\sum_a\Omega_a^2-3h(\rho^{-1}\partial_t+\Omega_0)
 +3a\partial_y.
\]

The angular structural maps \(\Omega_a\) are made from the invariant
angular metric and its identity tensor. They are sphere-parallel when the
angular label and all input/output tensor slots are included. This justifies
the absence of a hidden \(D_a\Omega_a\) term; it does not differentiate
normal-frame Christoffel symbols as constants.

Adding gives \(+4a\partial_y\). Since
\(\partial_y=-\rho^{-1}\partial_z\),
\(\partial_y^2+4a\partial_y=
\rho^{-2}\partial_z^2-3a\rho^{-1}\partial_z\).
Multiplication by \(\rho^2\) therefore gives exactly construction (2),
including \(Z=-\Omega_0^2-3h\Omega_0+\sum\Omega_a^2\).
Connection on a covariant rank-two tensor is the sum of the two one-slot
actions. Thus \(\|\Omega_0\|\le2|a|\) and
\(\|\Omega_a\|\le2(|a|+|h|)\) in the positive Frobenius norm.

For additional clarity, the only scalar gauge contraction appearing in (1)
is explicitly

\[
\begin{split}
 dC[k]_y={}&\tfrac12\partial_y(k_{yy}+k_{00}-\sum_a k_{aa})
 -\rho^{-1}\partial_tk_{y0}+(\rho c)^{-1}\sum_aD_ak_{ya}\\
 &+a(4k_{yy}+k_{00}-\sum_a k_{aa})-3hk_{y0}.
\end{split}
\]

This checks the scalar derivative and connection terms independently of a
choice of harmonic coordinates or tensor basis.

## 3. Finite angular space, spectrum and evaluation

The cover of the connected sphere isometry group is
\(SU(2)_L\times SU(2)_R\), with diagonal \(SU(2)\) as the isotropy group.
The real degree-two scalar representation is \((1,1)\), dimension nine.
The elementary symmetric-square identity and the spin-one tensor product give

\[
\mathrm{Sym}^2(1,1)=(0,0)\oplus(1,1)\oplus(2,2)
                    \oplus(2,0)\oplus(0,2).
\]

Indeed the symmetric square of spin one is spin zero plus spin two, while
its exterior square is spin one. The five representation dimensions add to
\(1+9+25+5+5=45\). This is a representation list, not the section-space
dimension. Every equivariant bilinear differential source of the first pair
factors through this symmetric square: polarization makes the two first
arguments symmetric even when the displayed local formula initially does not.

The fibre of the metric/scalar bundle, decomposed under isotropy, is
\(5[0]\oplus2[1]\oplus[2]\). The five spin-zero components are
\(k_{yy},k_{y0},k_{00},\mathrm{tr}_\kappa k,w\). The two spin-one copies are
\(k_{ya},k_{0a}\); the last copy is the angular trace-free symmetric tensor.
Restriction of \((j_L,j_R)\) to the diagonal group contains precisely
the spins \(|j_L-j_R|,\ldots,j_L+j_R\), once each. Sections of the associated
spin-\(s\) bundle have one copy of an irrep for each occurrence of spin \(s\)
in this restriction: an equivariant section is equivalently a function
\(f(gh)=\tau_s(h)^{-1}f(g)\), so its multiplicity is the dimension of the
intertwiner space from the restricted irrep to \(\tau_s\).

Applying that rule to the five allowed irreps gives

| Isotropy spin | Allowed section dimensions | Fibre copies | Contribution |
|---|---|---:|---:|
| 0 | 1, 9, 25 | 5 | 175 |
| 1 | 9, 25 | 2 | 68 |
| 2 | 9, 25, 5, 5 | 1 | 44 |

Thus the full dimension is 287. The one Dirichlet scalar and one Dirichlet
one-form give dimension \(35+34=69\), with complement 218. The scalar
gradient and trace-free Hessian identify the spin-one and diagonal spin-two
copies; these maps are nonzero for degrees two and four. The off-diagonal
irreps \((2,0),(0,2)\) contain spin two but no spin one, hence their divergence
vanishes; they supply the two five-dimensional TT sectors. No coexact-vector
irrep occurs. Counting off-diagonal metric entries twice changes the invariant
inner product, not these multiplicities.

For the unit sphere, the horizontal Casimir is
\((J_L-J_R)^2=2(C_L+C_R)-(J_L+J_R)^2\).
On isotropy-spin-\(s\) equivariant functions the last operator is \(s(s+1)\).
The canonical symmetric-space connection is the sphere Levi-Civita connection
and its induced tensor connections. Its rough-Laplacian eigenvalue is therefore

\[
 2[j_L(j_L+1)+j_R(j_R+1)]-s(s+1).
\]

The normalization is fixed by the scalar \((1/2,1/2)\) eigenvalue three on
the unit three-sphere. The permitted scalar values are \(0,8,24\), one-form
values \(6,22\), and trace-free tensor values \(2,18,6,6\). Thus
\(\|\Delta_\kappa\|\le24\). Integration by parts in the positive fibre
metric gives \(\|D_\kappa v\|^2=-\langle v,\Delta_\kappa v\rangle\), so
\(\|D_\kappa\|\le\sqrt{24}\) with its derivative slot summed.

Let \(s_j\) be any orthonormal basis of the 287-dimensional space for sphere
probability measure. The diagonal evaluation kernel is
\(K(p)=\sum_js_j(p)\otimes s_j(p)^*\). Its trace is constant by transitivity
and invariance of the positive fibre metric. Its integral is 287, so its
constant trace is 287. Since \(K(p)\) is positive,
\(\|\mathrm{ev}_p\|^2=\|K(p)\|\le\mathrm{tr}K(p)=287<17^2\).
This bound applies equally to \(U_x,U_r\), which remain in the same angular
space. It is not a scalar-component spectral cutoff.

## 4. Bounds on every operator block

The pinned background report controls the entire \(8\ell\) collar.
The checker verifies its domain and the coarse bounds
\(129\le\rho\le131\), \(|a|\le14/125\),
\(|a_y|\le10^{-4}\), \(|p|\le7/50000\).
The potential polynomial is differentiated independently with rational
coefficients. This yields \(|V|<1, |V'|<10^{-3}, |V''|<4\).
Using \(p_y=V'-4ap\) gives
\(\|\nabla^2\phi_0\|_F\le|p_y|+2|ap|<1\).

Warped sectional components are \(-(a_y+a^2)\) and \(\rho^{-2}-a^2\).
Every frame curvature entry is at most
\(R_*=a^2+|a_y|+129^{-2}<13/1000\). The warped Ricci endomorphism is
diagonal with norm at most \(4R_*\). The full Riemann contraction on a
25-entry tensor array has operator norm at most \(25R_*\), by its Frobenius
matrix norm. Hence \(\|L_R\|\le29R_*<1\). Restriction to symmetric tensors
cannot increase that norm.

For each of the scaled time and radial derivative blocks, the metric
\(dw\) coupling has norm at most \(4\ell\rho|p|\), and the scalar \(dC\)
coupling at most \(3\ell\rho|p|\), as follows from the explicit formula
above. Adding block norms is safe on the Hilbert direct sum. The angular
metric coupling costs at most \(20\ell^2\rho|p|\); the scalar divergence
costs at most \(10\ell^2\rho|p|\), using
\(\sqrt3\sqrt{24}<10\). Its connection part costs at most
\(\ell^2\rho^2|p|(8|a|+6|h|)\); summing the one time and three angular
connection actions suffices for this conservative bound.

The angular connection cross term obeys

\[
\Big\|\frac{2\rho}{c}\sum_a\Omega_aD_a\Big\|
 \le 2\rho\,[\sqrt3\,2(|a|+|h|)]\sqrt{24}
 <40\rho(|a|+|h|).
\]

The remaining zeroth blocks have bounds respectively
\(2\rho^2\) (curvature), \(4\rho^2/3\) (metric potential),
\(4\rho^2\) (scalar-to-metric potential, since
\((4/3)\|g_0\|_F=(4/3)\sqrt5<4\)), \(4\rho^2\) (scalar potential),
and \(\rho^2\) (scalar Hessian contraction). Together with \(24\), the
angular connection term, and \(\rho^2Z_*\), these are exactly the terms of
construction (4). No dimension factor of 287 is required for these bounded
geometric operations.

The frame Robin norm uses the norm two of the tangential metric injection
and trace. The separate terms sum to
\(\ell\rho(10|a|+2|\sigma'|+|p|/2+|\sigma''|/2)\).
The actual signed background tube is
\(-9\cdot10^{-5}\le\phi_0-1\le-8\cdot10^{-5}\), so \(0<\phi_0<1\).
It follows that \(|\sigma''|/2=2\phi_0<2\), validating the final constant
two in construction (9). An unsigned bound alone would instead require
2.00018; the signed tube is authenticated and checked explicitly here.

The independent rational block sum agrees exactly with the producer:

\[
 A_*=\frac{58816389}{10^{12}},\quad
 B_*=\frac{2207219}{5\cdot10^{10}},\quad
 D_*=\frac{299136118485884125939}{1386750000000000000000000000},\quad
 M_*=\frac{40880253}{10^{11}}.
\]

They are strictly below \(10^{-4},10^{-4},10^{-6},1/2000\).
The signed-tube clarification changes none of these numbers.

## 5. Initial/source joins and the quantitative conclusion

Unlike a prose-only cap lookup, this audit reads the pinned actual initial
spacetime frame report. Its component caps are
\(1.06\cdot10^{24},6.4\cdot10^{29},5.13\cdot10^{33},3.06\cdot10^{39}\)
for \(k,\nabla k,w,dw\). The metric Frobenius factor is at most five.
Radial parallelism gives \(U_{0r}=-\ell\rho\nabla_yU_0\); the time frame
conversion gives \(k_{0x}=\ell\rho(\nabla_0k_0-\Omega_0k_0)\).
In a norm bound the subtracted connection contributes its absolute size,
\(2\ell\rho|a|\|k_0\|\). These reproduce the producer initial joins and
\(\|U_{\rm free}\|_X<10^{36}\).

The copied bulk and shell reports are checked by SHA256 and read as exact
rationals. They imply the rounded construction bounds on
\(Q_{\rm red},Q_\phi,Q_R,Q_P,Q_C\), and
\(|G_N|\le8\max|Q_R|+2|Q_C|<1.3\cdot10^{10}\).
The positive source norm therefore gives exactly
\(F_H\le171610068644000<2\cdot10^{14}\) and
\(G_H\le3851400<4\cdot10^6\).

With the integral/trace estimates stated in the construction, the norm
contraction is

\[
 q=3(10^{-4}+10^{-4}+3\cdot10^{-6})+3/2000
   =2109/10^6,
 \qquad (1-q)^{-1}=10^6/997891.
\]

The value weight in the norm is essential: the boundary trace costs three,
not one. The resulting exact bound

\[
 \frac{10^{36}+3(2\cdot10^{14})+4\cdot10^6}{1-q}
 <1.02\cdot10^{36}
\]

agrees with the producer report. Multiplying by \(3\cdot17\) for values,
or by 17 for scaled first derivatives, proves the respective
\(6\cdot10^{37}\) and \(2\cdot10^{37}\) caps. These are conservative
upper bounds on the integral-defined formal coefficient, not waveform values,
nonlinear remainder estimates, or bounds for second derivatives.

## 6. Replay and scope

Run `python check_evolution_review.py --output NEW_PATH` and the same command
under `python -O`. Their reports are byte-identical. Nine direct mutations
are each rejected with the recorded specific reason and no output report:
curvature sign, derivative-slot connection, reduction symmetrization, scalar
gauge coupling, missing TT sectors, scalar-component cutoff, unsupported
evaluation factor, omitted trace weight, and wrong initial norm conversion.
`EVOLUTION_REVIEW_RECEIPT.json` pins the proof, checker, inputs, both reports
and the explicit control replay.

No outstanding blocker was found within this review's scope. The construction
still depends on its retained background/master and accepted physical initial
and boundary data. The finite-dimensional integral and subsidiary arguments
have their own proof and fresh review; this audit does not substitute algebraic
tests for those arguments. No frozen predecessor artifact was changed.
