# A finite spectral cutoff with a uniform scalar-error bound

This quantitative illustration retains the packet’s projected preparation and the fixed readout from the four-dimensional tower. Write $r=z/\ell$ on $[0,6]$, $x=t/\ell$, $\Lambda=\mu/\ell$, $b_0=\ell\beta_0$, $V_0=\ell^2V$. In scaled variables the static positive operator $A_0$ has mass-space $\mathcal H_0=L^2(0,6)\oplus\mathbb R_\Lambda$, form $a_0=\int(|u_r|^2+V_0|u|^2)dr+b_0|u(0)|^2$ and eigenvalues $\ell^2\lambda_j$. The physical scalar is $G\mu\cosh^{-3/2}(\ell x)B_0((5/2)I-A_0/\ell^2)Z(x)$, where $B_0$ extracts the endpoint. Throughout the admitted slab $\|A_0 Z(x)\|_{\mathcal H_0}<417$. The earlier proof established this graph estimate on its extended half-line. The registered solution vanishes in a neighborhood of $r=6$ there, so restriction preserves both the operator expression and its norm. Also $V_0<103000\ell^2$, $b_0<300\ell$, and $A_0\ge5547\ell^2 I$. [Graph and coefficient bounds, §§1, 3, 4, 6](../HDBLAST_PAST_SHELL_CONTACT_20261104/primary/PROOF.md) · [Gain bound, §1](../HDBLAST_CHARACTERISTIC_ODD_TRACE_20261104/primary/PROOF.md)

Let $P_{>K}$ be the spectral projection onto scaled eigenvalues above $K$. Approximate the endpoint Riesz vector $p=(0,1/\Lambda)$ by the form-domain vector $v_h=((1-r/h)_+/\Lambda,1/\Lambda)$. Direct integration gives $\|p-v_h\|^2=h/(3\Lambda^2)$ and $a_0[v_h]\le(1/h+103000\ell^2 h/3+b_0)/\Lambda^2$. The spectral theorem bounds $\|P_{>K}v_h\|\le K^{-1/2}\sqrt{a_0[v_h]}$. Consequently

$$
 \|P_{>K}p\|\le\Lambda^{-1}\left[\sqrt{h/3}+K^{-1/2}\sqrt{h^{-1}+103000\ell^2h/3+300\ell}\right].
$$

Take $K=10^{16}$, $h=10^{-8}$ and $\ell=10^{-6}$. The first term in brackets is less than $0.6\times10^{-4}$ since $1/3<0.6^2$. The second is less than $1.01\times10^{-4}$ by squaring the positive quantities. In spectral variables $0<\lambda_{\rm scaled}/\ell^2-5/2\le\lambda_{\rm scaled}/\ell^2$, so the graph vector paired against $P_{>K}p$ has norm below $417/\ell^2$. Since $|G|<3/250$ and $\mu/\Lambda=\ell$, the fixed-preparation truncated output $S_K$ satisfies throughout the admitted clock interval

$$
 |S-S_K| < (3/250)(417/\ell)(1.61\times10^{-4})=805.644<806.
$$

An upper bound on the number of retained species follows from min–max. If the subspace with eigenvalues $\le K$ has dimension $d$, imposing $u(0)=0$ leaves at least $d-1$ dimensions. On that subspace the mass-space norm is the bulk $L^2$ norm, both endpoints vanish, and $\int|u_r|^2\le a_0[u]\le K\int|u|^2$. The ordinary Dirichlet eigenvalues on $[0,6]$ are $(n\pi/6)^2$, so $d\le1+\lfloor6\sqrt K/\pi\rfloor\le200,000,000$. The bound is deliberately loose and establishes sufficiency only. In particular the two window means, whose absolute values exceed 14,000, retain their signs with these fixed projected data. No numerical mass list, minimal species count or physical observation threshold follows from this budget.
