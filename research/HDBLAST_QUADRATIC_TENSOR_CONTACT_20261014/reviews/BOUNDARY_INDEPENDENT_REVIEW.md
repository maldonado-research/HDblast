# Independent review of the forced tensor contact note

Date: 2026-10-14 UTC. This is an internal mathematical review, not external
peer review. The reviewer independently derived and directly checked the
bulk tensor operator and quadratic source in the adjacent `astra-bulk/`
directory before reading the boundary note.

Reviewed file:
`boundary/FORCED_TENSOR_CONTACT_IBVP.md`

Reviewed SHA-256:
`2ce6de9a3915466c060b3cd7d1b47646ff3fba318966f7f3c6b5cdfa9c3e19c3`.

**Disposition: no fatal gap found in the stated finite-mode longitudinal-gauge
forced tensor theorem.** The source, canonical conjugation, shell signs,
characteristic construction, and auxiliary energy are consistent. One
nonfatal clarification is recommended for bookkeeping when switching
between longitudinal and material Gaussian gauges; see paragraph 4 below.
This review does not certify the coupled second-order shell problem.

## 1. Operator and quadratic source normalization

The boundary note's equation (3) agrees with the independent direct
five-dimensional Christoffel/Ricci calculation in
`astra-bulk/verify_tensor_operator.py`, which checked all 15 independent
Einstein residual components using the actual Weyl-generated tensor
harmonic in Hopf coordinates. For the second derivative convention,

\[
 \delta E^{\rm TT}_{ij}
 =-\frac12\rho^2\cosh^2t\,\mathcal L_TH\,T_{ij},\qquad
 \mathcal L_TH=H_{yy}+4aH_y
 -\rho^{-2}(H_{tt}+3\tanh t H_t+\lambda\operatorname{sech}^2tH).
\]

The de Sitter tensor curvature shift is \(\lambda=k_T^2+2\), hence
\(\lambda=8\) for the specified \(l=2\) tensor. This is not \(k_T^2=6\).
The factor \(-1/2\) is fixed by the direct check, and the Taylor convention
is the same in both notes. Therefore

\[
 D^2E^{\rm TT}=CT\quad\Longrightarrow\quad
 \mathcal L_TH=\frac{2C}{\rho^2\cosh^2t}.
\]

The independently verified complete source is
\(C=16u^2+8v^2/3\) for \(\Psi=uY\), \(\chi=vY\). The displayed
time-evolution source in Section 2,
\(-(32u^2+16v^2/3)\operatorname{sech}^2t\), has the correct sign and
normalization. No matter-only substitution is being accepted.

## 2. Canonical transformation and shell flux

Using \(\partial_y=\rho^{-1}\partial_x\), the tensor equation before
conjugation is

\[
 H_{tt}+3\tanh t H_t-H_{xx}-3\rho' H_x
 +\lambda\operatorname{sech}^2tH=-\rho^2J.
\]

Set \(P=\rho^{3/2}\cosh^{3/2}t\), \(v=PH\),
\(S=(\log P)_x=3\rho'/2\). Direct differentiation gives

\[
 W_l=S^2+S_x-\frac94+
       \left(\lambda+\frac34\right)\operatorname{sech}^2t,
 \qquad F=-\rho^{7/2}\cosh^{3/2}t\,J
          =-2\rho^{3/2}\cosh^{-1/2}t\,C.
\]

Every power, the \(-9/4\), and the \(+3/4\) curvature correction agree.
For \(H_y|_b=g\),

\[
 v_x-S_bv=\rho_b^{5/2}\cosh^{3/2}t\,g.
\]

With inward \(z=-x\), this becomes \(-v_z+\alpha v=G\) with
\(\alpha=-S_b\). The energy boundary term is consequently
\(+\alpha|v(0)|^2/2\), usually negative here; changing its sign would
destroy the flux cancellation. The note keeps the correct sign.

## 3. Which shell law is proved

I compared Section 3 with the frozen shell-Noether convention and the
material supplement's exact trace-reversed Israel law,

\[
 K_{+,ab}+K_{-,ab}=\frac{\sigma}{3}h_{ab}.
\]

Both radial coordinates increase outward. The sum of normal derivatives
therefore has the sign used in equation (6). Under the stated first-order
even longitudinal representative, \(\zeta_1=0\), the induced metric and
first radial extrinsic-curvature terms are conformal to the angular
background. All quadratic tension and lapse terms have zero spatial TT
projection. The linear second-order shift contributes an angular symmetric
derivative; a second-order normal displacement contributes a Hessian and
trace; tangential reidentification contributes a symmetric derivative.
Their spatial TT projections vanish on the closed sphere.

For the TT metric derivative itself,

\[
 K_{2,ij}^{\rm TT}
 =\rho^2\cosh^2t\,(aH+H_y/2)T_{ij}.
\]

Matching makes the two traces equal. Using \(\sigma_b/3=2a_b\) cancels
the two warp terms, leaving \(H_{+,y}+H_{-,y}=0\). Thus the even
amplitude is Neumann and the odd amplitude is Dirichlet. Neither
\(\partial_y g_{2,ij}^{\rm TT}=0\) nor a generic Eulerian homogeneous
Neumann condition would be justified.

The note correctly avoids extending this homogeneous Eulerian law to
first-order odd graph states. A copy-odd first derivative can have an even
second derivative; its quadratic embedding terms still require treatment.
The homogeneous even result is sufficiently qualified by equation (5).

## 4. Gauge conversion: correct rule, one useful clarification

The material pullback formula

\[
 g_2^{\rm mat}=g_2^w+2\mathcal L_Vg_1^w
 +\mathcal L_V^2g_0+\mathcal L_Wg_0
\]

agrees with the material supplement. The angular TT projection of
\(\mathcal L_Wg_0\), and of its radial derivative, vanishes. Thus the
inhomogeneous boundary traces in equation (9) are determined entirely by
first-order gauge/embedding data. The longitudinal-to-Gaussian generator
has \(V^y=O(y-y_b)\), \(V^a=O((y-y_b)^2)\) and zero boundary value;
the quadratic TT correction consequently has \(B|_b=B_y|_b=0\), as stated.

The interior correction \(B\) generally does not vanish. Accordingly the
source in a Gaussian material gauge need not equal the explicit
longitudinal source from Section 2. The note's source-shift rule is correct,
but the distinction could be made explicit as

\[
 H^{\rm mat}=H^{\rm long}+B_{\rm long},\qquad
 J^{\rm mat}=J^{\rm long}+\mathcal L_TB_{\rm long},
\]
\[
 H^{\rm mat}=H^w+B_w,\qquad
 J^w=J^{\rm mat}-\mathcal L_TB_w.
\]

This prevents inserting the unchanged numerical longitudinal expression
for \(C\) after a nontrivial first-order Gaussian gauge change. It is a
clarification of a rule already present, not a failure of the stated
longitudinal theorem. No second-order gauge-invariance claim is needed.

## 5. Half-line construction, corner compatibility, and control energy

The characteristic equations, directions, and boundary ODE are correct:

\[
 p_t-p_z=F-Wv,\quad m_t+m_z=F-Wv,\quad
 q'=p(0)-\alpha q+G,\quad m(0)=p(0)-2\alpha q+2G.
\]

Here \(p\) leaves the physical half-line at the shell and \(m\) enters.
There is one prescribed incoming characteristic per parity. Substitution
of the outgoing integral into the incoming integral and the trace ODE
indeed gives an \(O(\delta)\) contraction in the joint \((v,q)\) uniform
norm; the unsmoothed \(\alpha q\) boundary term acquires its required
time factor in the final integration. The defects
\(v_z-(p-m)/2\) and \(v|_0-q\) have zero time derivative once the
first-order system holds. With their zero initial values, no hidden
constraint is lost.

For smooth compatible data, differentiated Volterra equations supply the
successive derivatives. Compatibility across the characteristic from the
corner is necessary; the note assumes it explicitly and provides the
vanishing-collar data for which all orders hold. The proof does not use a
bare Galerkin energy bound as an unsupported smoothness theorem. Its
classical existence argument is the finite-mode Volterra construction.

The Robin form is semibounded. For instance, when \(\alpha<0\), take
\(\epsilon=(2|\alpha|)^{-1}\) in the stated trace inequality; then a bulk
shift \(K>\|W_-\|_\infty+2\alpha^2+1\) gives a positive control norm.
The integration-by-parts boundary flux cancels exactly with the derivative
of \(\alpha q^2/2\). The remaining terms in equation (15) are bounded by
the control energy and \(\|F\|_2^2\). Smooth boundary lifts have the stated
sign and correctly require two time derivatives of the forcing. No
unsupported boundary-\(L^2\) estimate or spectral positivity claim is used.

## 6. Actual background domain and what crosses contact

Compact initial support \(z\le L\) and the stipulated source support give
\(z\le L+t\). Thus a finite slab stays in
\(r\ge b\exp[-L-t_*]>0\), uses only the actual positive-\(s\) background,
and never sees the arbitrary auxiliary extension. This is a controlled
collar result, not a global ambient-time moving-ball theorem or a claim
about a cone endpoint at finite global de Sitter time.

The contact formula
\(t_c=\log(b/R_0)\),
\(T_c=b\sinh t_c=(b^2-R_0^2)/(2R_0)\) agrees with ambient null propagation.
The inherited first-order construction can actually cover a postcontact
interval in this chart. If its hemisphere conformal time is \(\tau\),

\[
 \tau(x,t)=\arctan(e^{x+t})-\arctan(e^{x-t})
 \le\arctan(\sinh t)\quad (x\le0,t\ge0).
\]

Therefore choosing
\(t_*<\operatorname{arsinh}(\tan\tau_*)\) keeps the entire global-time
collar slab inside the inherited \(|\tau|\le\tau_*\) solution. Since
the inherited first contact occurs strictly before \(\tau_*\), such
\(t_*\) can be chosen beyond \(t_c\). Section 8's controlled-slab
qualification is consequently realizable, not merely formal.

The resulting statement is a smooth forced tensor **coefficient** through
contact, with specified compatible tensor initial data. It neither proves
a nonzero boundary trace exactly at the support edge nor supplies the
missing scalar, vector, gauge, and constraint boundary equations. It is
not an exact nonlinear shell-evolution theorem.

## Evidence used

- Independently derived `astra-bulk/FULL_QUADRATIC_TT_BULK_SOURCE.md` and
  the all-15-component direct tensor-operator checks, passed under Python
  3.12.14 and SymPy 1.14.0 in normal and optimized modes.
- The inherited `HDBLAST_SHELL_NOETHER_20261013/manuscript/SHELL_NOETHER.md`
  for outward-normal and Israel conventions.
- The readback `HDBLAST_MATERIAL_JUNCTIONS_20261014/MATERIAL_SECOND_ORDER_JUNCTIONS.md`
  for material derivatives, trace reversal, and the exact quadratic tension
  source.
- The inherited `HDBLAST_REGULAR_SHELL_CONTACT_20261013/existence/FINITE_SLAB_SHELL_CONTACT_EXISTENCE.md`
  for the actual first-order domain, compatible collar class, and controlled
  postcontact existence.

No edits were made to the boundary mathematical note during this review.

## Final readback acceptance

The final boundary note was read back on 2026-10-14 UTC and its SHA-256
was directly verified as

`e9651b98a5a0b50e2607a5d83e7dfc944c5528ae6de5461d67a384708039f2c5`.

**Accepted at this final hash.** The following additions were reviewed:

1. Section 4 now explicitly distinguishes the actual material source from
   the longitudinal source and gives both
   \(J_{\rm mat}=J_{\rm long}+\mathcal L_TB_{\rm long}\) and
   \(J_w=J_{\rm mat}-\mathcal L_TB_w\). This resolves the clarification
   recommended in paragraph 4 of this review.
2. Section 8 supplies the concrete zero-TT initial-data choice. On the
   background time-symmetric slice, an angular TT metric variation has
   zero full spatial trace/divergence and zero contraction with the radial
   scalar gradient and angular background Ricci tensor. It therefore lies
   in the linear Hamiltonian kernel. An angular TT momentum density has
   zero full spatial covariant divergence and lies in the linear momentum
   kernel. Because the background momenta vanish, these subtractions do
   not introduce cross-coupling into the other linearized constraint.
   They leave the already fixed quadratic constraint source unchanged.
   In the specified first-order longitudinal gauge, first-order metric,
   lapse, and kinetic products are angular trace, so zero TT momentum is
   equivalent to the stated zero TT metric velocity. The paragraph does
   not assert freedom to vary arbitrary other constraint data.
3. The explicit controlled interval
   \(t_c<t_*<\operatorname{arsinh}(\tan\tau_*)\) and its coordinate
   inequality agree with the domain calculation in paragraph 6 above.
   This makes the postcontact interval concrete while retaining the
   positive-\(s\), finite-mode collar scope.

No new test run was needed for this readback. The original independent
bulk and operator checks remain applicable. No fatal gap or unresolved
review comment remains for the stated forced tensor coefficient theorem;
the full coupled second-order shell problem remains outside its scope.
