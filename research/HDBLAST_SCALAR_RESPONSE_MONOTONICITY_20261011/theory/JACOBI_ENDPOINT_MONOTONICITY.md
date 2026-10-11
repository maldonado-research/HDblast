# Regular-cone Jacobi identities and a stronger uniform branch estimate

11 October 2026. Local mathematical continuation for review. This file does not
modify or replace the accepted project files. Its conclusions concern the static
Einstein–scalar boundary problem and its specified small Dirichlet tube. No
global stability, dynamical, cosmological, or external novelty claim is made.

## 1. Conclusions

The sufficient ``strict downward shell crossing'' inequality in Section 6 of
the 8 October uniform branch proof can be deleted. Under its other construction
and scalar-junction hypotheses, the metric residual along the scalar-junction
locus satisfies the exact sign estimate

\[
 \frac{d}{dT}\left[a(T;\beta(T))-\frac{\sigma(\beta(T))}{6}\right]
 \leq -\rho(T;\beta(T))^{-2}<0.                       \tag{1}
\]

The first inequality is strict for every nonconstant scalar interior. The
argument also applies to the derivative norm with an arbitrary positive scale
\(s\), with its correctly propagated constants given below. In particular,
the conservative former slope constant is unnecessary; it was an artifact of
estimating two terms separately that cancel at leading order.

A separate comparison estimate improves the linear Dirichlet-to-Neumann
remainder to

\[
 0\leq\lambda-\frac{F_\lambda'(T)}{F_\lambda(T)}
 \leq\frac{\lambda}{k(w-3k)}H_0(T)^2.               \tag{2}
\]

Consequently the constant called \(C_H\) in the baseline proof may be replaced by
\[
 C_H^{\rm new}=\frac{\lambda}{k(w-3k)}.             \tag{3}
\]

Neither improvement changes the boundary expansion coefficients. Their utility
is a larger interval and smaller *proved upper bounds* on the remainders, when
the remaining sufficient inequalities are evaluated.

## 2. Equations, variables, and normalization

Let \(W\in C^4([-r,r])\), \(f\in C^3([-r,r])\),
\[
 W(0)=3k,\quad W'(0)=0,\quad W''(0)=w>4k>0,
 \quad f(0)=f_0>0,
\]
and set
\[
 U=\frac12(W')^2-\frac23W^2,\quad
 m^2=U''(0)=w(w-4k),\quad\lambda=w-4k,\quad
 \sigma=2W+\delta f.
\]
Primes on \(U,W,f,\sigma\) mean scalar derivatives; primes on radial quantities
mean derivatives in \(y\). Put
\[
 p=\eta',\qquad a=\rho'/\rho,\qquad \mathcal H=\rho^{-2}.
\]
For \(y>0\), the regular interior obeys
\[
 p'=U'(\eta)-4ap,\qquad
 a'=-\mathcal H-\frac{p^2}{3},\qquad
 a^2=\mathcal H+\frac{p^2}{12}-\frac{U(\eta)}6.       \tag{4}
\]
Equivalently,
\[
 a'=-a^2-\frac{p^2}{4}-\frac{U(\eta)}6.             \tag{5}
\]
The cone data are \(\rho(0)=0,\rho'(0)=1,p(0)=0\). Its sole scalar datum is
\(c=\eta(0)\). Radius normalization is fixed: it is not an extra Jacobi parameter.
At the endpoint, \(\beta=\eta(T)\) and \(H^2=\mathcal H(T)\).

The weighted construction supplies a \(C^1\) family of regular interiors at
each finite \(T>0\), with \(\eta(T;T,\beta)=\beta\). All partial derivatives in
\(\beta\) below hold \(T\) fixed. The scalar junction is
\[
 g(T,\beta)=p(T;T,\beta)+\frac12\sigma'(\beta)=0,
 \quad
 g_\beta=p_\beta+\frac12\sigma''>0.                 \tag{6}
\]
In the baseline sufficient conditions, \(g_\beta\geq w/2\).

## 3. Cone Jacobi fields and their endpoint identities

Differentiate the cone solution with respect to \(c\), at fixed \(y\):
\[
 h=\eta_c,\quad q=p_c=h',\quad
 v=(\log\rho)_c,\quad b=a_c=v'.
\]
For finite regular segments \(\rho>0\), cone parameter dependence and ordinary
ODE parameter dependence justify these derivatives. Their equations are
\[
 q'=U''h-4aq-4pb,\qquad
 b'=2\mathcal H v-\frac23pq,\qquad v'=b.            \tag{7}
\]
The differentiated constraint is
\[
 2ab=-2\mathcal H v+\frac{pq-U'h}{6}.               \tag{8}
\]
Define the two Jacobi combinations
\[
 j=b+\frac{ph}{3},\qquad
 \mathcal D=h p'-p q.                              \tag{9}
\]
Direct substitution of (4), (7), and (8) gives
\[
 j'+4aj=-2\mathcal H v,                            \tag{10}
\]
\[
 \mathcal D'+4a\mathcal D
       =4p\bigl(\mathcal H h+p j\bigr).             \tag{11}
\]
For clarity, (8) implies
\(pq-p'h=12aj+12\mathcal H v\). Substituting this into
\(j'=2\mathcal H v+(p'h-pq)/3\) proves (10).
For (11), differentiate \(hp'-pq\), cancel \(qp'-p'q\), and use
\(p''=U''p-4a'p-4ap'\). The result is
\(-4a\mathcal D-4a'ph+4p^2b\), which is (11).

The regular cone expansions, including their parameter derivatives, are
\[
 \begin{aligned}
 \rho&=y-U(c)y^3/36+O(y^5),&
 p&=U'(c)y/5+O(y^3),\\
 h&=1+U''(c)y^2/10+O(y^4),&
 q&=U''(c)y/5+O(y^3),\\
 v&=-U'(c)y^2/36+O(y^4),&
 b&=-U'(c)y/18+O(y^3).
 \end{aligned}                                    \tag{12}
\]
Only the displayed leading terms and smaller remainders are needed. They follow
from the regular cone integral equations and their differentiated equations.
Thus
\[
 j=U'(c)y/90+O(y^3),\qquad
 \mathcal D=U'(c)/5+O(y^2).                         \tag{13}
\]
In particular, *both* \(\rho^4j\) and \(\rho^4\mathcal D\) tend to zero at
the cone; \(\mathcal D\) itself need not tend to zero. Integrating (10)–(11)
with their proper cone boundary terms yields
\[
 j(y)=-2\rho(y)^{-4}\int_0^y\rho(t)^2v(t)\,dt,      \tag{14}
\]
\[
 \mathcal D(y)=4\rho(y)^{-4}
   \int_0^y\rho(t)^4p(t)
        [\mathcal H(t)h(t)+p(t)j(t)]\,dt.           \tag{15}
\]
The second integrand is \(O(t^3)\) at zero, so no singular boundary term is
hidden in this formula.

Whenever \(h(T)\ne0\), the cone-to-endpoint map is locally invertible and
\[
 p_\beta(T)=q(T)/h(T),\quad
 a_\beta(T)=b(T)/h(T),\quad
 Q(T):=a_\beta(T)+p(T)/3=j(T)/h(T).                 \tag{16}
\]
The regular Dirichlet family already ensures this nondegeneracy: writing
\(c=c(T,\beta)\), differentiation of its endpoint condition gives
\(h(T)c_\beta=1\). The positivity proof below additionally gives \(h(T)\geq1\).

Changing \(T\) at fixed \(\beta\) gives \(c_T=-p/h(T)\), and therefore
\[
 \partial_Tp\vert_\beta=p'-p p_\beta=\mathcal D/h(T),
 \quad
 \partial_Ta\vert_\beta=-H^2-p^2/3-p a_\beta.       \tag{17}
\]
All statements here concern finite \(T\); no implicit-function theorem at
\(T=\infty\) is used.

Along (6), \(\beta_T=-(\mathcal D/h(T))/g_\beta\). For
\(\mathscr R(T)=a(T;T,\beta(T))-\sigma(\beta(T))/6\), use
\(\sigma'=-2p\) to obtain the exact cancellation
\[
 \mathscr R_T=-H^2+Q(T)(\beta_T-p).                 \tag{18}
\]
This identity holds everywhere on the scalar-junction locus, not merely at a
metric-junction root.

## 4. Convexity gives the signs needed in (18)

Assume the scalar stays in an interval containing zero on which
\(U''>0\), and assume \(\rho>0\) away from the cone. If \(c\ne0\), let
\(\epsilon=\operatorname{sign}c\). The scalar divergence equation
\[
 (\rho^4p)'=\rho^4U'(\eta)
\]
and \(U'(0)=0\) show that \(\epsilon p>0\),
\(\epsilon U'(\eta)>0\), and \(\epsilon\eta>0\) for \(y>0\) on the segment:
\(\eta\) starts with sign \(c\), moves away from zero, and a first failure
of that assertion contradicts the integral formula for \(p\).

Put \(z=-\epsilon b\). Differentiation of (5), rather than of the constraint
form of \(a'\), gives the particularly simple positive system
\[
 (\rho^4q)'=\rho^4(U''h+4|p|z),\qquad
 (\rho^2z)'=\rho^2(|p|q/2+|U'|h/6),\qquad h'=q.    \tag{19}
\]
The cone terms in the first two integrated equations vanish. The coefficients
coupling different unknowns are nonnegative. By (12), for sufficiently small
positive \(y\) all of \(h,q,z\) are positive. Start at any such point; the
ordinary cooperative linear system preserves nonnegativity, as can also be
seen by a first-zero argument after applying the displayed integrating factors.
Thus on the full segment
\[
 h\geq1,\qquad q>0,\qquad z>0.
\]
This argument avoids assigning arbitrary zero data to a singular Jacobi
equation. The strict center coefficients fix the regular solution. Since
\(v(0)=0\), \(v'=b=-\epsilon z\) implies \(\epsilon v<0\).
Equation (14) now gives \(\epsilon j>0\), hence \(p j>0\).
Equation (15) then gives \(\epsilon\mathcal D>0\).

It follows from (16)–(17) that \(Q\) has the sign of \(p\) and
\(\partial_Tp\vert_\beta\) also has the sign of \(p\). With \(g_\beta>0\),
\(\beta_T\) has the opposite sign. In particular
\[
 Q(\beta_T-p)<0\quad(c\ne0),
\]
which proves (1), with strict first inequality. This proof requires no sign
assumption on \(\sigma''\) beyond the actual positive scalar-junction
derivative \(g_\beta\). If \(c=0\), the interior is the exact scalar-free cone,
\(p=b=v=j=0\), and (18) gives equality \(\mathscr R_T=-H^2\).

## 5. Strict convexity is already implied by the construction

This step is essential: Section 4 adds no convexity hypothesis to the uniform
theorem's surviving sufficient conditions. Let
\[
 M_2=\sup_{[-r,r]}|U''|,\quad M_3=\sup_{[-r,r]}|U'''|,
 \quad L_1=M_3/2,
\]
and choose \(\nu=\lambda/2\), \(\mu_\nu=\nu(\nu+4k)\),
\(D=m^2-\mu_\nu>0\). With arbitrary derivative scale \(s>0\), the relevant
weighted norm and constants are
\[
 \|\eta\|_{T,s}=
 \max\{\sup|\eta|/\omega_T,\sup|\eta'|/(s\omega_T)\},
 \quad \omega_T=F_\nu/F_\nu(T),
\]
\[
 \begin{aligned}
 M&=\max(1,\lambda/s),&B&=2M,\\
 C_G&=\max\{D^{-1},(1+m^2/D)/(4ks)\},&
 C_F&=s^2/4+M_2/12,\\
 C_N&=L_1B^2+4C_F(s/k)B^3b,&
 L_N&=2L_1Bb+20C_F(s/k)B^2b^2.
 \end{aligned}                                    \tag{20}
\]
The tube conditions include \(Bb\leq r\),
\(\|\eta\|_{T,s}\leq B|\beta|\leq Bb\), and
\(q_{\rm contr}:=C_G L_N\leq1/2\). Since
\(C_G\geq D^{-1}\) and \(L_N\geq M_3Bb\),
\[
 M_3Bb\leq D/2,
 \quad
 U''(x)\geq m^2-M_3|x|
      \geq m^2-D/2=(m^2+\mu_\nu)/2>0
 \quad (|x|\leq Bb).                              \tag{21}
\]
Thus every interior in the tube satisfies the strict convexity used above.

There is also a useful metric consequence. Since \(U'(0)=0\), strict
convexity on this interval implies \(U(\eta)\geq U(0)=-6k^2\).
Thus the Volterra source
\(\mathcal F=p^2/4+[U(\eta)+6k^2]/6\) is nonnegative along the actual
solution. The already established positivity of \(R=\rho/r_0\), together
with \(R'=-\int_0^y\partial_yK(y,t)\mathcal F(t)R(t)\,dt\) and
\(\partial_yK\geq0\), implies \(R'\leq0\). Since \(R(0)=1\),
\[
 0<R(y)\leq1,\qquad H_0(T)^2=R(T)^2H^2\leq H^2.   \tag{21a}
\]
This sign improvement is only asserted on fixed-point solutions; the
contraction estimates on arbitrary candidate functions retain their original
absolute-value bounds.

The norm propagation in (20) follows by replacing \(|p|\leq kB|\beta|\omega_T\)
by \(|p|\leq sB|\beta|\omega_T\). The metric source has bound
\(C_FB^2\beta^2\omega_T^2\), its logarithmic derivative still has geometric
denominator \(k\), and the product \((a-a_0)p\) therefore acquires \(s/k\).
The physical derivative Green bound is unchanged, and division by \(s\)
gives exactly the displayed \(C_G\). In particular setting \(s=k\) recovers
the original proof.

Keep every sufficient condition used for the positive metric, Dirichlet
contraction, scalar-junction monotonicity, and the opposite signs of the metric
residual at \(T=1/k\) and \(T\to\infty\). Delete only the former condition
\(C_{\rm slope}\bar\delta\leq a_1/4\). The remaining construction proves a
scalar-junction root \(\beta(T,\delta)\) for all \(T\geq1/k\), and Section 4
proves that its metric residual is strictly decreasing. Its opposite endpoint
signs therefore give exactly one metric root in the same tube. The strictly
negative derivative supplies the ordinary \(C^1\) implicit-function theorem
at every positive \(\delta\), just as in the baseline proof.

## 6. A quantitative alternative that does not use sign arguments

This independent estimate isolates the curvature factor even without convexity.
It can be useful for other settings. Let \(R=\rho/r_0\),
\(r_0=\sinh(ky)/k\), and let \(J\) be the uniform bound
\(\int_0^T\omega_T^2\leq J\) used in the construction. Suppose its Volterra
smallness bound is \(\varepsilon_V=C_FB^2b^2J/(2k)\leq1/4\).
Then \(2/3\leq R\leq4/3\), and differentiation of
\(R=1-\mathcal V(\mathcal F R)\), where
\(\mathcal F=p^2/4+(U+6k^2)/6\), yields
\[
 \|R_\beta\|_\infty\leq
 \frac{16}{9}\frac{C_FB^2J}{k}|\beta|.
\]
Indeed \(|\mathcal F_\beta|\leq2C_FB^2|\beta|\omega_T^2\), the inverse
Volterra norm is at most \(4/3\), the metric factor at most \(4/3\), and
the kernel at most \(1/(2k)\). Hence, for the endpoint-normalized variation
\(v_\beta=\partial_\beta\log\rho\),
\[
 \|v_\beta\|_\infty\leq C_v|\beta|,
 \qquad C_v=\frac{8C_FB^2J}{3k}.                  \tag{22}
\]
If \(a\geq a_{\min}>0\), then
\(\int_0^T\rho^2\leq\rho(T)^2/(2a_{\min})\). Dividing (14) by \(h(T)\)
therefore proves the explicit curvature-factor bound
\[
 |a_\beta+p/3|\leq\frac{C_v}{a_{\min}}H^2|\beta|. \tag{23}
\]
The existing metric estimate permits
\[
 a_{\min}=k-C_FB^2b^2/k\geq k/2.
\]
With the baseline constants
\[
 \begin{aligned}
 C_P&=sC_GB[2L_1B+20C_F(s/k)B^2b],\\
 \bar a&=k\coth(1)+C_FB^2b^2/k,\\
 C_{Tp}&=M_2+4\bar a sB+sB(\lambda+C_Pb),\\
 C_{T\beta}&=2C_{Tp}/w,
 \end{aligned}
\]
one has \(|\beta_T-p|\leq(C_{T\beta}+sB)|\beta|\). For
\(|\beta|\leq L\delta\), (18) gives
\[
 \mathscr R_T\leq-H^2(1-\theta),\qquad
 \theta=\frac{C_v}{a_{\min}}(C_{T\beta}+sB)L^2\bar\delta^2. \tag{24}
\]
Thus \(\theta<1\) is a sufficient explicit substitute for the former slope
condition, independently of Sections 4–5. For rational computation,
\(\coth(1)<4/3\) gives a rational upper replacement for \(\bar a\).

## 7. Improved linear endpoint comparison

Let \(F_\lambda\) be the positive regular solution of
\[
 F_\lambda''+4a_0F_\lambda'-m^2F_\lambda=0,
 \quad F_\lambda(0)=1,\ F_\lambda'(0)=0,
 \quad a_0=k\coth(ky).
\]
Its logarithmic derivative \(z=F_\lambda'/F_\lambda\) obeys
\(z'=m^2-4a_0z-z^2\) and \(0\leq z\leq\lambda\), by the positive divergence
equation and the first-crossing test at \(z=\lambda\). Put
\[
 d=a_0-k>0,\quad e=\lambda-z\geq0,\quad E=e/d.
\]
The exact equations are
\[
 d'=-d(d+2k),\quad
 e'=-(w+z)e+4dz,
 \quad
 E'=4\lambda-2(w-3k)E+dE(E-3).                   \tag{25}
\]
At the regular cone \(d\sim1/y\), \(z=O(y)\), so \(E\to0\). Let
\(C=2\lambda/(w-3k)\). Since \(w>4k\), one has \(0<C<2\). At the
barrier \(E=C\), the first two terms of (25) cancel while
\(dC(C-3)<0\). Starting at a sufficiently small positive radius, the usual
first-crossing argument therefore proves \(0\leq E<C\) for every finite
positive radius. Since
\(H_0^2=a_0^2-k^2=d(a_0+k)\) and \(a_0+k\geq2k\), this yields (2).

At a nonlinear endpoint, (21a) gives \(H_0^2=R(T)^2H^2\leq H^2\).
The unchanged nonlinear endpoint bound
\[
 |p(T)-z(T)\beta|\leq C_D\beta^2,
 \qquad C_D=sC_GC_N,
\]
therefore gives
\[
 |p(T)-\lambda\beta|\leq
 C_H^{\rm new}H^2|\beta|+C_D\beta^2,              \tag{26}
\]
with (3). If the sign information (21a) is unavailable in a different problem,
the absolute-value metric bound \(R\leq4/3\) still gives the weaker constant
\(16\lambda/[9k(w-3k)]\). Unlike the former exponential estimate, (26) holds at every
\(T>0\); the branch construction may still use \(T\geq1/k\) for its other
conditions.

## 8. Remainder constants and scope

Choose any \(L\geq |f'(0)|/w\), and use the notation
\[
 a_1=kf_0/3,\quad d_0=2(w-2k),
 \quad\alpha=-f'(0)/(2d_0),
\]
The choice \(L=|f'(0)|/w\) is valid because \(g(T,0)=\delta f'(0)/2\)
and \(g_\beta\geq w/2\) imply
\(|\beta(T,\delta)|\leq |f'(0)|\delta/w\). The baseline extra \(1/w\)
in \(L\) is not needed. If \(f'(0)=0\), then \(\beta=0\) is exactly the
unique scalar-junction root, and any independently positive Dirichlet tube
radius \(b\) may be used; taking \(L=0\) is harmless in the remainder formulas.

Define the endpoint functions
\[
 A(x)=W(x)f(x)/9-W'(x)f'(x)/12,\quad
 B_0(x)=f(x)^2/36-f'(x)^2/48.
\]
With \(W_3=\sup|W'''|\), \(F_2=\sup|f''|\), all on \([-r,r]\), the
improved sufficient remainder constants are
\[
 K_\eta^{\rm new}=
 \frac{C_H^{\rm new}(3a_1/2)L+(C_D+W_3/2)L^2+F_2L/2}{d_0},
\]
\[
 K_H^{\rm new}=|A'(0)|K_\eta^{\rm new}
                 +\frac12\sup|A''|L^2+\sup|B_0'|L.
\]
They give the same uniform estimates
\[
 |\beta-\alpha\delta|\leq K_\eta^{\rm new}\delta^2,
\]
\[
 \left|H^2-a_1\delta-
 \left[\frac{f_0^2}{36}-\frac{k f'(0)^2}{24(w-2k)}\right]\delta^2\right|
 \leq K_H^{\rm new}\delta^3.
\]
For polynomial endpoint functions \(A(x)=\sum_i A_i x^i\) and
\(B_0(x)=\sum_i B_i x^i\), an additional exact refinement is available.
For a fixed interval endpoint \(\bar\delta\), put
\[
 E_* = \min\{L,|\alpha|+K_\eta^{\rm new}\bar\delta\},
\]
\[
 K_H^{\rm poly}=|A_1|K_\eta^{\rm new}
  +\sum_{i\geq2}|A_i|E_*^i\bar\delta^{i-2}
  +\sum_{i\geq1}|B_i|E_*^i\bar\delta^{i-1}.        \tag{27}
\]
Indeed \(|\beta|\leq E_*\delta\), and the exact shell constraint is
\(H^2=\delta A(\beta)+\delta^2 B_0(\beta)\). Subtract its constant and
linear-in-\(\beta\) terms, use
\(|\beta-\alpha\delta|\leq K_\eta^{\rm new}\delta^2\), and bound each
remaining polynomial term. This proves the same curvature remainder with
\(K=\min(K_H^{\rm new},K_H^{\rm poly})\). No series truncation is involved.

For the distinct comparator with scalar-independent detuning \(f\equiv f_0\), its
curvature is exactly \(a_1\delta+f_0^2\delta^2/36\). Strict lower curvature
therefore follows only when the additionally evaluated margin
\(k f'(0)^2/[24(w-2k)]-K\bar\delta>0\) is positive. This final
sign claim must not be inferred from existence or monotonicity alone.

## 9. Complete sufficient-condition checklist and construction summary

Here is a complete way to apply this continuation without having to infer which
old inequality survives. Choose \(r,b,s,\bar\delta>0\) and
\(L\geq|f'(0)|/w\), with the model assumptions of Section 2.
Use (20), the definition of \(C_P\) in Section 6, and
\[
 C_0=4k\coth(1)+\nu,\quad
 \zeta=\frac{\mu_\nu[1-\exp(-C_0/k)]}{C_0},\quad
 J=\frac2k+\frac1{2\zeta}.
\]
Define derivative bounds on \([-r,r]\):
\[
 W_2=\sup|W''|,\quad W_3=\sup|W'''|,\quad
 F_2=\sup|f''|,\quad F_{\rm abs}=\sup|f|,
 \quad C_K=L\sup|A'|+\sup|B_0|.
\]
The following eleven conditions suffice; the sixth is retained as a strict
first-stage feasibility condition, although the scalar-junction derivative bound
is what is directly used for each \(\delta\):

| Condition | Required inequality |
|---|---|
| Field tube | \(Bb\leq r\) |
| Metric Volterra inverse | \(C_FB^2b^2J/(2k)\leq1/4\) |
| Scalar contraction | \(C_GL_N\leq1/2\) |
| Scalar ball preservation | \(C_GC_Nb\leq M\) |
| Positive metric derivative | \(C_FB^2b^2/k\leq k/2\) |
| Scalar derivative feasibility | \((C_P+W_3)b<w/2\) |
| Scalar root inside tube | \(L\bar\delta\leq b\) |
| Scalar-junction derivative | \((C_P+W_3)b+\bar\delta F_2/2\leq w/2\) |
| Positive reduced curvature | \(C_K\bar\delta\leq a_1/2\) |
| Positive residual at \(Y=1/k\) | \((3a_1/2)\bar\delta<(9/16)H_0(Y)^2\) |
| Positive shell tension | \(W_2b^2/6+\bar\delta F_{\rm abs}/6\leq k/2\) |

An upper bound for any displayed nonnegative constant is sufficient. In
particular, wholly rational replacements are
\[
 \zeta\geq\frac{(15/16)\mu_\nu}{16k/3+\nu},\qquad
 J\leq\frac2k+\frac{16k/3+\nu}{(15/8)\mu_\nu},
 \qquad H_0(1/k)^2>k^2/2.
\]
They follow from \(\coth(1)<4/3\), \(C_0/k>4\), \(\exp(4)>16\),
and \(\sinh(1)^2<2\). Thus the residual condition can be checked by the
stronger rational inequality \((3a_1/2)\bar\delta<9k^2/32\).

For completeness, the construction behind this checklist has the following
steps. For \(\omega_T=F_\nu/F_\nu(T)\), positivity of
\((r_0^4F_\nu')'=\mu_\nu r_0^4F_\nu\) gives
\(0<\omega_T\leq1\), \(0\leq\omega_T'\leq\nu\omega_T\), and the
large-radius lower bound \(F_\nu'/F_\nu\geq\zeta\) for \(y\geq2/k\)
gives \(\int_0^T\omega_T^2\leq J\). The scalar Green inverse of
\(\mathcal L_m v=v''+4a_0v'-m^2v\), with regular center data and \(v(T)=0\),
has weighted norm at most \(C_G\). Its scalar bound follows by the maximum
principle from \(\mathcal L_m\omega_T=-D\omega_T\); integration of its
divergence equation gives the derivative bound in (20).

For any trial scalar in \(\|\eta\|_{T,s}\leq B|\beta|\), the metric is the
unique Volterra solution
\[
 R(y)=1-\int_0^y
 \frac{\sinh(k(y-t))\sinh(kt)}{k\sinh(ky)}
 \mathcal F(t)R(t)\,dt.
\]
The kernel is between zero and \(1/(2k)\), and its radial derivative is
\([r_0(t)/r_0(y)]^2\). The second checklist condition gives
\(2/3\leq R\leq4/3\) and
\(|a-a_0|\leq C_FB^2\beta^2\omega_T^2/k\).
For two trial scalars at norm distance \(d_\eta\), the corresponding estimate
is \(|a-\widetilde a|\leq4C_FB|\beta|d_\eta\omega_T^2/k\).

Writing the scalar equation as
\[
 \eta=\beta F_\lambda/F_\lambda(T)
 +G_T\{U'(\eta)-m^2\eta-4(a-a_0)\eta'\}
\]
gives exactly the source bound \(C_N\beta^2\), Lipschitz bound \(L_N\),
and norm bounds (20). The first five checklist conditions therefore construct
a unique regular scalar in this ball, with positive metric and \(a\geq k/2\),
uniformly in finite \(T\). The common \(Bb\) ball and the invertible
contraction linearization give \(C^1\) parameter dependence and
\(\|\eta_\beta\|_{T,s}\leq B\). Differentiating the scalar fixed-point
equation gives \(|p_\beta-z(T)|\leq C_P|\beta|\); its undifferentiated
remainder gives \(|p-z(T)\beta|\leq C_D\beta^2\).

The scalar derivative conditions then give
\(g_\beta\geq z(T)+w-(C_P+W_3)b-\bar\delta F_2/2\geq w/2\), so the
intermediate value theorem constructs its unique root within
\(|\beta|\leq L\delta\). On this root locus the exact constraint gives
\[
 (a-\sigma/6)(a+\sigma/6)
  =H^2-\delta A(\beta)-\delta^2B_0(\beta).
\]
The positive-curvature condition puts the last two terms together in
\([a_1\delta/2,3a_1\delta/2]\). The final two conditions and
\(2/3\leq R\leq4/3\) give positive residual at \(T=1/k\) and negative
residual as \(T\to\infty\), for each \(\delta>0\). Sections 3–5 give its
strict decrease, hence its unique zero in the tube. Section 7 and Taylor
expansion of the exact boundary constraint give Section 8's remainder bounds.

This proof summary uses no shooting data, continuation interpolants, or
floating-point residuals. It proves every \(0<\delta\leq\bar\delta\)
satisfying the sufficient inequalities, rather than just selected samples.

## 10. Registered exact-decimal model: arithmetic corollary

For the specified rational polynomial model, take
\[
 k=1/9,\quad W(\eta)=1/3+\eta^2+\eta^3/3,
 \quad f(\eta)=1+c+c\eta,
 \quad c=\frac{2}{1.0357712571566784}-\frac43,
\]
where the displayed decimal is interpreted exactly as a rational number.
Choose
\[
 r=1/100,\quad s=14/9,\quad B=2,\quad
 L=c/2=\frac{48357160712846}{161839258930731},\quad
 \bar\delta=1/500,\quad b=L/500.
\]
The separate self-contained exact-arithmetic producer
`../bounds/evaluate_rescaled_bounds.py` and its output
`../bounds/RESCALED_BOUND_CERTIFICATE.json` verify all eleven conditions.
They also directly verify polynomial convexity bounds, as redundant checks
of the general implication (21). The improved endpoint constant is exactly
\(C_H^{\rm new}=42/5\). The arithmetic gives the rational upper bounds
\[
 K_\eta^{\rm new}<3.196,\qquad
 K=K_H^{\rm poly}<0.254,
 \qquad \frac{c^2}{384}-K\bar\delta>0.000423.       \tag{28}
\]
Here the terminating decimals in (28) are also exact rational comparisons,
not rounded substitutions into an inequality. Approximate displays are
\(K_\eta^{\rm new}\simeq3.1951835148\),
\(K\simeq0.2533589973\); the actual exact fractions are in the certificate.

Consequently, for **every** \(0<\delta\leq0.002\), there is exactly one
regular static shell with \(T\geq1/k\) in the stated Dirichlet tube, and
\[
 \left|\beta+\frac{9c}{64}\delta\right|<3.196\delta^2,
\]
\[
 \left|H^2-\frac{1+c}{27}\delta
    -\left[\frac{(1+c)^2}{36}-\frac{c^2}{384}\right]\delta^2\right|
 <0.254\delta^3.                                  \tag{29}
\]
At the largest detuning these bounds are respectively
\(1.2784\times10^{-5}\) and \(2.032\times10^{-9}\).
Relative to the distinct model with scalar-independent detuning
\(f_{\rm metric}(\eta)=1+c\), its exactly known curvature is
\[
 H_{\rm metric}^2=\frac{1+c}{27}\delta
                   +\frac{(1+c)^2}{36}\delta^2,
\]
and (28) proves the strict statement
\[
 H^2-H_{\rm metric}^2<-0.000423\delta^2<0
 \qquad (0<\delta\leq0.002).                       \tag{30}
\]
This corollary covers the three previously investigated detunings 0.0005,
0.001, and 0.002, and all positive detunings below them in the displayed
interval. It certifies the mathematical branch and its analytic error bound;
it does not certify the floating-point shooting values individually or
identify their parameter-rounding conventions with the exact-decimal model.

All uniqueness conclusions remain restricted to the constructed tube. The
proof neither rules out other large-field interiors nor analyzes perturbative
stability. Exact arithmetic evaluations and independent internal review are
separate artifacts; they do not establish external peer review or novelty.
