# A finite-precision test of bandwidth and preparation

A tenth time sample after nine silent samples places a lower bound on the prepared amplitudes of free four-dimensional fields with restricted masses. The test permits signed, nonreciprocal couplings, arbitrary initial displacements and velocities, and arbitrarily many species. It uses ten point readouts at $x=k/8$ for $k=0,\ldots,9$, including the preparation time $x=0$. All clocks, gains and angular normalizations are those of the published model comparison.

Let $A=(\ell/G)S(-\ell x)$ be the target in calibrated physical-amplitude units and $u=\cosh^{3/2}(\ell x)A$. The target is zero for $0\le x\le1$ by initial separation, finite propagation and continuity. At $T=9/8$ the two adjacent closed waveform-cell enclosures and the past-transfer bound give $u(T)>0.365949$. These are predictions of the model. [Support proof](https://github.com/maldonado-research/HDblast/blob/20f431f2f6fd8ae9d040d0e8b87bddb41a2f3fa2/research/HDBLAST_PAST_SHELL_CONTACT_20261104/primary/PROOF.md) · [Normalization and transfer](https://github.com/maldonado-research/HDblast/blob/20f431f2f6fd8ae9d040d0e8b87bddb41a2f3fa2/research/HDBLAST_DIMENSIONAL_COMPARISON_20261118/SINGLE_FIELD_TEST.md)

Write a free four-dimensional candidate in the physical units as $A_{\rm c}=\sum_j c_j a_j(x)$, with $P_j=(m_j\rho_b\ell)^2\ge0$ and $v_j=\cosh^{3/2}(\ell x)a_j$. Define its initial frequency $\omega_j=\sqrt{P_j+13\ell^2/2}$ and the total absolute oscillation amplitude
$$
 B=\sum_j b_j,\qquad b_j=|c_j|\sqrt{v_j(0)^2+v_j'(0)^2/\omega_j^2}.
$$
The derivative is in the scaled clock $x$. The couplings include the fixed gain. $B$ has the same units as $A$ and is assumed finite. Its numerical value is independent of the signs chosen for the individual fields. For a calibrated positive frequency floor $\omega_{\min}$, the initial energy of the previous test obeys $B\le D\sqrt{2E_0}/\omega_{\min}$ by Cauchy–Schwarz and $K_j(0)=P_j+8\ell^2\ge\omega_j^2$. This supplies one way to translate an independent energy and detector calibration into an upper bound on $B$. [Energy and detector conventions](https://github.com/maldonado-research/HDblast/issues/151)

Freeze the transformed modal equation at $x=0$. Its constant-frequency solution with the same initial data is $v_j^0=v_j(0)\cos(\omega_jx)+v_j'(0)\sin(\omega_jx)/\omega_j$. The exact coefficient of $v_j$ differs from $\omega_j^2$ by $-(35/4)\ell^2\tanh^2(\ell x)$. Its absolute value is at most $(35/4)\ell^4 T^2$ on $[0,T]$. In Duhamel's formula use $|\sin(\omega_j z)|/\omega_j\le z$ and $\int_0^T z\,dz=T^2/2$. With $\eta=(35/8)\ell^4 T^4$ and $\ell=10^{-6}$ this gives
$$
\|c_j(v_j-v_j^0)\|_\infty\le \frac{\eta}{1-\eta} b_j <10^{-22}b_j
$$
for $b_j>0$, with zero difference for $b_j=0$. This bound is uniform in mass. The readout error is at most $10^{-22}B$.

Use the ninth forward difference at zero with step $h=1/8$. Its coefficients $(-1)^{9-k}\binom9k$ have absolute sum $512$. On $\exp(i\omega x)$ the difference is $(e^{i\omega h}-1)^9$. Thus each frozen real modal contribution has difference at most $b_j(2|\sin(\omega_j h/2)|)^9$. If $|A_{\rm c}-A|\le\delta$ separately at each of the ten samples, the errors in $u$ are bounded by $\delta/f$, where $f=1-10^{-11}<\operatorname{sech}^{3/2}(\ell x)$ throughout. The nine early target values vanish, so its ninth difference equals its tenth sample $u(T)$. The resulting necessary inequality for any such candidate is
$$
 u(T)\le \sum_j b_j(2|\sin(\omega_j h/2)|)^9
              +512(10^{-22}B+\delta/f).
$$
For an independently imposed frequency ceiling $\omega_j\le\Omega$, the sum is bounded above by $B(\Omega h)^9$. Inserting the exact lower endpoint for $u(T)$ yields the following strict lower budgets, rounded down.

| Frequency ceiling in the $x$ clock | Exact target at the ten samples | Deterministic error $\delta=1/2000$ at each sample in $A$ units |
| --- | --- | --- |
| $\Omega=2$ | $B>95{,}931$ | $B>28{,}822$ |
| $\Omega=4$ | $B>187.36$ | $B>56.29$ |

This is an amplitude cost for reproducing ten registered model values with a specified bandwidth. The two ceilings illustrate the dependence on an independently chosen mass range. No observational bandwidth, source-amplitude budget or measurement uncertainty has been supplied for an actual cosmological source. The original exact tower has unbounded frequencies. The table assumes point samples with the stated timing and amplitude tolerance.

## An implementation with finite exposures

A second version uses the physical trace $A$ averaged uniformly over each of nine fixed windows $I_k=[9k/64,(9k+1)/64]$ for $k=0,\ldots,8$. The first eight windows lie within the silent interval. The last, $[9/8,73/64]$, has model mean greater than $0.51393$ from the integrals of the certified physical-amplitude cell bounds.

Take the eighth difference of the observed means with stride $9/64$. Its coefficients have absolute sum 256. Averaging each frozen mode over a common width amounts to multiplying its complex amplitude by the mean of $e^{i\omega s}$ over that exposure, an additional factor of modulus at most one. The bound on the exact-to-frozen transformed oscillator remains below $10^{-22} B$ up to $x=73/64$. In the physical readout, replacing the dilution factor by one costs at most a further $10^{-11} B$, since every frozen mode's amplitude is bounded by its $b_j$. This gives
$$
 \overline A_{I_8}\le B[(9\Omega/64)^8+256(2\times10^{-11})]+256\delta
$$
for an additional deterministic error of at most $\delta$ on each of the nine physical means. At $\delta=1/2000$ the necessary budget is **$B>9{,}856$ for $\Omega=2$** or **$B>38.50$ for $\Omega=4$**. This version needs uniform exposure weights and timing, and makes no assumption about pointwise measurement error within a window.

The numerical check reads the previous comparison's three authenticated input files and uses exact rational arithmetic. [Replay](bandwidth.py) · [Rational outputs](BANDWIDTH.json)
