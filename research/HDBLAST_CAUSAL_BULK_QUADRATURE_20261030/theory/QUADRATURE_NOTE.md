# Causal panel enclosure of the first actual-bulk correction

This calculation concerns the same conditional reduced system and the same
fixed actual coefficient and affine-source functions as the earlier
actual-bulk update. It evaluates an enclosure of the correction

\[
C=W_1-W_0=(I-\mathcal L M P_N)^{-1}K r_0,\qquad
r_0=(L_a-L_0)B_s+L_a(W_0-B_s)+F.
\]

The cancellation-free residual formula avoids subtracting two large copies of
the free wave. The interior evaluator bounds the complete modal
source on every integration panel. It retains the actual Robin matrix and
the actual affine sources. Their uniform magnitude bounds enclose the same
fixed functions; they are not replaced by zero functions or midpoint values.
The affine bulk source itself is not newly evaluated pointwise.

At the shell point \((a,0)=(9/8,0)\), set
\(\alpha=s+y\), \(\theta=(y-s)/(s+y)\). The active physical source triangle
becomes \(1\leq\alpha\leq a\), \(-1\leq\theta\leq1\), with positive measure
\(ds\,dy=(\alpha/2)d\alpha\,d\theta\). The inherited zero wedge is a separate
support premise. Its validity does not follow from a small source norm.

For every rectangular panel the checker multiplies a whole-panel source
interval by the exact positive mass of
\(\alpha(a-\alpha)^j/(2j!)\). Summing these products encloses

\[
S_j=\int_{s,y\geq0,\,1<s+y<a}
       \frac{(a-s-y)^j}{j!}\,r_0(s,y)\,ds\,dy.
\]

This is a range enclosure of an integral, including quadrature error. It does
not infer error from agreement between sampled meshes. The exact partition
cover and positive masses are checked using rational arithmetic. Splitting at
the initial-profile joins improves ranges while preserving the complete
physical domain.

Folding the Robin history gives the trace enclosure

\[
P_N C(a,0)=P_N\sum_{j=0}^{N}(-M)^jS_j+\text{tail}.
\]

The complete physical source is coupled before applying the output projector.
The reflection is not moved through a coupling matrix. Repeated interval
matrix products include the same fixed actual matrix; discarded correlations
can widen the result but do not freeze its entries.

Writing \(d=a-1\), the exact whole-triangle weight is

\[
A_j=\frac{d^{j+1}}{j!}
     \left(\frac{a}{j+1}-\frac{d}{j+2}\right).
\]

If \(R_{\rm local}\) is the inherited genuine pointwise source norm cap and
\(m\) bounds the actual Robin operator, the omitted trace has norm at most

\[
R_{\rm local}\,
\frac{m^{N+1}A_{N+1}}{1-md/(N+2)}.
\]

Indeed, the ratios of successive positive tail terms are bounded by
\(md/(j+1)\). This bound encloses the infinite integral-defined correction,
not merely the displayed finite polynomial.

The complete physical null-edge integral is independently enclosed as

\[
p(a)=\frac a2\int_{-1}^1r_0\!\left(
\frac{a(1-\theta)}2,\frac{a(1+\theta)}2\right)d\theta.
\]

The shell first jets then satisfy

\[
P_N C_x=P_Np-MP_NC,\quad P_NC_r=MP_NC,
\qquad P_D C=P_DC_x=0,\quad P_D C_r=P_Dp.
\]

The checker also intersects its numerical correction boxes with the prior
analytic correction bounds and reports both raw and intersected intervals.
An intersection is used only after checking overlap. The separate transfer
calculation can combine these correction intervals with the earlier
reconstruction and the authenticated global target-error bound.

For the published mesh, 16 alpha panels and 32 theta panels are supplemented
by the exact join splits, giving 544 physical source panels and 34 edge
panels. With 128 primitive panels per analytic branch, all 53 nontrivial
correction boxes are strictly narrower than the earlier symmetric analytic
ones. The finite meshes are used to obtain rigorous whole-panel bounds;
their agreement is not the justification for the error enclosure.

In the registered scalar-mean component, exact rational comparisons of the
reported endpoints give the readable strict bounds

\[
1.596\times10^{21}<C_{w:0}(9/8,0)<6.733\times10^{21},
\]
\[
8.355\times10^{22}<\partial_x C_{w:0}(9/8,0)<1.058\times10^{23},
\]
\[
4.147\times10^{17}<\partial_r C_{w:0}(9/8,0)<1.750\times10^{18}.
\]

All three correction intervals exclude zero. Across the complete table, 16
value intervals, 16 time-derivative intervals and 17 radial-derivative
intervals exclude zero. These are statements about this integral-defined
coordinate correction. In particular, the earlier negative scalar-mean
residual at the endpoint does not determine the sign of the causal correction:
the latter depends on the entire physical source triangle and Robin history.
The source value itself has not been inferred from these correction signs.

The coarse mesh is retained as a diagnostic. The final interval enclosure,
not a mesh-difference estimate, establishes the displayed signs. The
scalar-mean correction remains many orders of magnitude smaller than the
prior scalar-mean reconstruction at this coordinate sample.

These are conditional enclosures of one reduced-system integral update at
one shell point. They do not establish a self-consistent second-order solution,
a quantitative second-derivative error, an Einstein or constraint residual,
an exact nonlinear spacetime, a physical force, or a higher-dimensional origin
of the Big Bang.
