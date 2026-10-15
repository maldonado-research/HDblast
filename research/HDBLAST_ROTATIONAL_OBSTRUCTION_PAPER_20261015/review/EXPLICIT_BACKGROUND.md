# Direct nonvacuity example for the revised shell obstruction draft

Derived independently during the fresh review on 15 October 2026. This note gives an elementary exact background for the paper's hypotheses. It does **not** identify this choice of potential/tension with the registered HDBLAST background, certify that earlier branch, prove shell-contact evolution, or claim a new exact-solution discovery. Length and scalar normalizations below are fixed dimensionless units, with the manuscript's convention \(8\pi G=1\).

## Bulk fields and potential

In five ambient dimensions with \(\eta=\operatorname{diag}(-1,1,1,1,1)\), set

\[
 s=X^AX_A=R^2-T^2,\qquad z=1+s>0,
\]

\[
 g_{AB}=z^{-4}\eta_{AB},\qquad
 \Omega=z^{-2},\qquad
 \Phi=\sqrt6\log z,\qquad
 U(\phi)=84e^{2\phi/\sqrt6}-36e^{3\phi/\sqrt6}.
\]

For a general conformal factor \(g=e^{2u(s)}\eta\) in five dimensions,

\[
 R_{AB}=12(u_s^2-u_{ss})X_AX_B
       +(-16u_s-4su_{ss}-12su_s^2)\eta_{AB}.
\]

Here \(u=-2\log z\), so \(u_s=-2/z\), \(u_{ss}=2/z^2\), and therefore

\[
 R_{AB}=\frac{24}{z^2}X_AX_B
       +\frac{32-24s}{z^2}\eta_{AB}.
\]

Since \(\Phi_s=\sqrt6/z\),

\[
 \Phi_A\Phi_B=\frac{24}{z^2}X_AX_B.
\]

Evaluating the potential on the background gives

\[
 U(\Phi)=84z^2-36z^3=(48-36s)z^2,
\]

and thus

\[
 \frac23U(\Phi)g_{AB}=\frac{32-24s}{z^2}\eta_{AB}.
\]

Consequently \(R_{AB}-\Phi_A\Phi_B-(2/3)U(\Phi)g_{AB}=0\) identically.

The scalar equation also checks directly. For any radial scalar \(v(s)\),

\[
 \Box_gv=\Omega^{-2}\{4sv_{ss}+10v_s+12su_sv_s\}.
\]

Substitution gives

\[
 \Box_g\Phi=\sqrt6(10-18s)z^2.
\]

On the other hand,

\[
 U'(\Phi)=\frac{168z^2-108z^3}{\sqrt6}
          =\sqrt6(10-18s)z^2,
\]

which proves the scalar equation without using a contracted Bianchi implication.

## Two identical interiors and the physical shell

Use two identical regions \(s\le b^2\), joined at the timelike hyperboloid \(W:s=b^2\), with \(b>0\). Their outward unit normals are both represented in their own interior coordinates by

\[
 n^A=\frac{X^A}{b\Omega_b},\qquad \Omega_b=(1+b^2)^{-2}.
\]

This vector has unit positive norm and points to increasing \(s\). The induced metric and scalar agree exactly between the copies. Let \(\phi_b=\sqrt6\log(1+b^2)\). The conformal transformation of the flat hyperboloid curvature gives

\[
 K_{ab}=a_bh_{ab},\qquad
 a_b=\frac1{\Omega_b}\left(\frac1b+2b\frac{\Omega_s}{\Omega}\Big|_b\right)
     =\frac{(1-3b^2)(1+b^2)}b.
\]

The shell has dimension four, hence \(K=4a_b\). The outward scalar derivative on either copy is

\[
 p_b=n\Phi=\frac{2b\Phi_s(b^2)}{\Omega_b}
     =2\sqrt6\,b(1+b^2).
\]

The manuscript's junction residuals consequently are

\[
 I_{ab}=(-6a_b+\sigma(\phi_b))h_{ab},\qquad
 P=4\sqrt6\,b(1+b^2)+\sigma'(\phi_b).
\]

Both vanish if the fixed smooth tension satisfies

\[
 \sigma(\phi_b)=6a_b,\qquad
 \sigma'(\phi_b)=-4\sqrt6\,b(1+b^2).
\]

No values of \(\sigma\) away from \(\phi_b\) are needed for this exact background. One may choose it affine globally.

At the convenient radius \(b=1/2\),

\[
 a_b=\frac58,\quad p_b=\frac{5\sqrt6}4,\quad
 \phi_b=\sqrt6\log\frac54,
\]

and a completely specified choice is

\[
 \sigma(\phi)=\frac{15}4-\frac{5\sqrt6}2
           \left(\phi-\sqrt6\log\frac54\right).
\]

Then \(\sigma(\phi_b)=15/4>0\), \(\sigma'(\phi_b)=-5\sqrt6/2\), and the two identical-copy metric/scalar junction equations vanish exactly.

## Full short slab, time symmetry, and witness annulus

At fixed \(T\), each interior section is the compact four-ball

\[
 0\le R\le\sqrt{b^2+T^2}.
\]

For \(|T|\le T_0<1\), its full physical range satisfies

\[
 -T_0^2\le s\le b^2,
\]

so \(1+s\) stays uniformly positive. The Cartesian metric, scalar, inverse metric and all derivatives are smooth on each closed interior piece of the doubled slab, including the center and \(s=0\). Across the material shell, the matched metric is piecewise smooth with the extrinsic-curvature jump fixed by the junction equations. The scalar derivative \(f=\sqrt6/(1+s)\) never vanishes within either interior. The shell is timelike, since its displayed normal is spacelike. Matching is exact and no external spatial boundary is introduced by the construction.

The \(T=0\) section is time symmetric: \(\partial_T\Omega=\partial_T\Phi=0\), the spatial shift vanishes, the canonical gravitational/scalar momenta vanish, and the shell is momentarily at rest. All canonical and cut data are invariant under the spatial rotations used in the paper. Thus the invariant-density hypotheses and the actual-background hypothesis are realized simultaneously.

For this example the coefficients in the manuscript simplify to

\[
 L=-\frac4{1+s},\qquad
 m_r=\frac38r^3(1+r^2)^4,\qquad
 n_r=6r^7(1+r^2)^2.
\]

The diagonal radial zeroth-order coefficient is

\[
 m_r(2+9r^2L)+n_r
  =m_r\frac{2-32r^2-18r^4}{(1+r^2)^2}.
\]

For \(0<r\le1/5\), its numerator is positive: it decreases with \(r\) and at \(r=1/5\) equals \(432/625\). Thus every nonzero \(Q\in C_c^\infty((1/10,1/5))\) gives

\[
 \mathcal B_{\rm rad}(Q,Q)
 =\int\left[m_r(rQ')^2+
       m_r\frac{2-32r^2-18r^4}{(1+r^2)^2}Q^2\right]\,dr>0
\]

directly, without needing the thin-annulus estimate. The paper's estimate remains valid and more general. If a particular smooth bump is desired, take

\[
 Q(r)=
 \begin{cases}
 \exp\!\left[-\dfrac1{(r-11/100)(19/100-r)}\right],&11/100<r<19/100,\\
 0,&\text{otherwise}.
 \end{cases}
\]

This function is smooth and flat at the endpoints. With \(b=1/2\), \(r_1=1/10\), \(r_2=1/5\), a precontact choice \(0<T_*<1/20\) satisfies both support bounds in the draft and remains inside the globally regular slab. Taking the paper's \(H=n_1^2-n_2^2\), \(J=x_1\partial_2-x_2\partial_1\), and any nonzero real \(\omega\) therefore produces its strictly nonzero rotational obstruction in this completely specified example.

This establishes nonemptiness of the paper's assumed background class and one explicit admissible precontact witness. It does not supply a nonlinear continuation through that forbidden tangent, a contact theorem, or a result for another fixed choice of \(U,\sigma\).
