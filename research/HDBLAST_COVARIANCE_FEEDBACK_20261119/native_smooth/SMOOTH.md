# Quiet shells with the same initial sources and their first derivatives

Two Gaussian preparations can agree on the initial parent energy, pressure, scalar force and the first time derivatives of all three and still generate distinct solutions. On the fine B1 grid their later quiet-shell expansion rates differ by $-4.49560\times10^{-7}$ in units of $H_0$ at $s=2.10$. Their inferred parent energy densities differ by $-0.134698\,H_0^4$. Both effects have the same sign on the coarser grid. This follow-up uses the four-preparation study's closed Gaussian parent model and its common $T=1.11$ fork. [Fork and model](../native/NATIVE.md) · [Computed differences](SMOOTH.json)

## Two preparations sharing the initial derivative data

We keep $X=\langle\chi^2\rangle$ and $P=\langle\Pi^2\rangle$ at every sampled momentum in both preparations. Define $t=(q/8)^2$ and $h(t)=t^2-b_1t-b_0$, with coefficients making $h$ orthogonal to $1$ and $t$ on the Gauss nodes with positive weights $w_q C_{\rm old}(q)$. Using the maximum absolute $h$ on the continuous interval $0\le t\le1$, set

$$
C_\pm(q)=\pm0.9\,C_{\rm old}(q)\,h(t)/\max|h|.
$$

This ensures $|C_\pm|<|C_{\rm old}|$ and permits a valid mixed Gaussian state at every node. The preparations are fitted on each mesh to the first two covariance moments of that mesh's common fork. Their initial one-point sources agree with each other and with those of the retained and erased states. The initial source derivatives agree with those of the erased state. The preceding natural evolution approaches the fork with the retained state's derivatives, so both replacements share the same derivative change at the preparation event.

At fixed common $X,P,a,\phi$, the initial differences in $j'$ and $p'$ depend only on $\sum w\Delta C$ and $\sum w q^2\Delta C$. Both sums vanish to quadrature roundoff by construction. The unchanged energy source and its derivative then match as well. The three Neumann data and their three first coordinate-time derivatives coincide for the two preparations. Higher moments can distinguish them. For instance $\sum w q^4(C_+-C_-)=10.26050$ on the fine grid. These signed quadrature moments use the same 24 nodes and the archived dimensionless $q$. [Initial covariance data](FORK.json)

## Evolution

We independently advance the quiet native shell from each state on both meshes. At the same elapsed proper clock, the following differences are $C_+$ minus $C_-$.

| At $s=2.10$ | Coarse wall | Fine wall |
|---|---:|---:|
| $\Delta(H/H_0)$ | −4.52395 × 10⁻⁷ | −4.49560 × 10⁻⁷ |
| Shell scalar $\Delta\phi_b$ | +1.10831 × 10⁻⁷ | +1.10061 × 10⁻⁷ |
| $\Delta\rho/H_0^4$ | −0.135548 | −0.134698 |
| $\Delta p/H_0^4$ | +0.335528 | +0.332872 |
| $\Delta j$ in hatted units | −0.925154 | −0.918315 |

The two-grid differences estimate sensitivity to the discretization. The construction demonstrates nonuniqueness of the later stress and expansion given just these initial moments in the modeled, finite-band scalar population. Both preparations here are engineered interventions in a common quantum state. [Four native wall series](WALLS.zip) · [Construction and solver interface](smooth.py)

~~~sh
python -I -B reproduce.py /path/to/HDblast
OPENBLAS_NUM_THREADS=1 python -I -B -O reproduce.py /path/to/HDblast --full --jobs 2
~~~
