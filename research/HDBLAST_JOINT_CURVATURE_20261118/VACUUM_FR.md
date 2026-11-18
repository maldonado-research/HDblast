# Testing four-dimensional metric f(R) gravity with the shell curvature

The gravitational comparison uses the Ricci scalar of the physical four-dimensional metric. A minimally coupled detector is assumed to measure that same metric in both descriptions. On the registered shell the background scalar and intrinsic Ricci scalar are both constant. Their linear perturbations therefore have no first-order coordinate ambiguity. Let $r(x)Y$ denote the scalar-curvature perturbation on the proper observer clock, where $x=q/\ell$, $\ell=10^{-6}$ and $-\Delta Y=8Y$. The shell constraint gives $r(x)=c_b S(-\ell x)$ with $c_b=2p_b B_b\ne0$. This relation makes the previously enclosed time dependence available for a test of a gravitational observable. [Registered field and clock](https://github.com/maldonado-research/HDblast/blob/59651ab35b954ef268c81a771e9ec8064691f9ad/research/HDBLAST_DIMENSIONAL_COMPARISON_20261118/SINGLE_FIELD_TEST.md)

## What a vacuum metric theory predicts

Use a local metric action $\int \sqrt{-g}\, f(R)/(2\kappa_4^2)$ with a smooth $f$ near $R_0=12/\rho_b^2$. The same de Sitter background must solve the vacuum equations, so $R_0 f'_0=2f_0$, with the subscript denoting evaluation at $R_0$. In signature $(-+++)$, the field equation with minimally coupled matter has trace
$$
3\Box f'(R)+R f'(R)-2 f(R)=\kappa_4^2 T.
$$
At linear order around vacuum, $f'_0$ is constant and varying the wave operator on that constant contributes zero. Thus
$$
3 f''_0\Box\delta R+(R_0 f''_0-f'_0)\delta R=\kappa_4^2\delta T .
$$
If $f''_0=0$ and $f'_0\ne0$, the source-free equation already forces $\delta R=0$. If $f''_0\ne0$, the curvature on any source-free interval must solve the scalar wave equation with constant squared mass $m_s^2=(f'_0/f''_0-R_0)/3$. For example $f(R)=R-2\Lambda+\alpha R^2$ with $\alpha\ne0$ has $m_s^2=1/(6\alpha)$. The degenerate case $f'_0=f''_0=0$ has no linear trace constraint from this equation. [Radhakrishnan et al., January 2026 preprint, Eqs. (26)–(29), (34)–(37)](https://arxiv.org/html/2601.01028v1#S2.E37)

Put $C=\cosh^{3/2}(\ell x)$ and $v=Cr$. For a curvature mode free on both test intervals, $v''+[P+\ell^2 k(x)]v=0$, where $P=(\rho_b\ell)^2m_s^2$ and $k(x)=(35/4)\operatorname{sech}^2(\ell x)-9/4$. On the target, $v$ is a constant nonzero multiple of the earlier certified waveform. Multiplication by that constant cancels in each test functional
$$
K_i[v]=-\frac{\int \psi_i'' v}{\int\psi_i v}
       -\ell^2\frac{\int k\psi_i v}{\int\psi_i v}.
$$
The two compact quadratic spline weights have zero value and first derivative at their outer endpoints and nonzero denominators on the target. If the mode were free there, each $K_i$ would have to equal the same $P$. The certified intervals are $(73.97,84.65)$ and $(45.64,60.23)$. Hence **no constant squared mass, of either sign, reproduces this intrinsic curvature waveform on both supports**. This excludes the regular vacuum metric $f(R)$ comparisons above and includes histories with a drive outside both supports if the same linearized equation still holds on them. Tensor waves in that theory also satisfy the linear trace constraint. The argument makes no assumptions about the candidate's otherwise free initial data. [Exact spline calculation and inputs](https://github.com/maldonado-research/HDblast/blob/59651ab35b954ef268c81a771e9ec8064691f9ad/research/HDBLAST_DIMENSIONAL_COMPARISON_20261118/SHAPE.json)

## What a matter source would have to supply

For $f''_0\ne0$, retaining the projected trace of a first-order matter source gives
$$
v''+[P+\ell^2 k(x)]v=-\frac{(\rho_b\ell)^2\kappa_4^2}{3f''_0}\, C\,\delta T_Y .
$$
The source on the right is expressed on the registered clock. Two integrations by parts show that its weighted integral divided by $\int\psi_i v$ equals $P-K_i[v]$. The disjoint intervals therefore require
$$
\max_{i=A,B}\left|
\frac{(\rho_b\ell)^2\kappa_4^2}{3f''_0}
\frac{\int\psi_i C\,\delta T_Y}{\int\psi_i v}\right|>6.87 .
$$
This is a necessary dimensionless ratio for a candidate matching the target. If the uniformly normalized curvature trace $\ell C r/(G c_b)$ differs from the registered normalized waveform by at most $1/2000$, the corresponding bound remains above $6.18$. An upper cap on physical matter would require independent units, coupling and source information. An evolving background curvature, additional gravitating fields, nonminimal matter couplings or a response that mixes several geometric quantities would change the comparison equation. [Uncertainty model and weak-test proof](https://github.com/maldonado-research/HDblast/blob/59651ab35b954ef268c81a771e9ec8064691f9ad/research/HDBLAST_DIMENSIONAL_COMPARISON_20261118/SINGLE_FIELD_TEST.md)
