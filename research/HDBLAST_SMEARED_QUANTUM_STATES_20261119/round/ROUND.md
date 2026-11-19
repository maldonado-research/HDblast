# Preparing the daughter field on the rolling sphere

This calculation asks how much the scalar daughter stress changes when an observer prepares the field by minimizing the energy measured over a slightly different time interval around the bounce. A shorter measurement makes a different choice of quantum fluctuations. Both signs of the resulting energy difference are possible. The background and the scalar action are held fixed within each comparison.

## Background, clock and density

We use the two unseeded rolling closed backgrounds of the preceding energy-exchange study, with the same fitted potential, launch displacements, first-minimum times and $H_0/M_{\rm Pl}=0.0088134414$ (reduced Planck mass). Write $x=H_0t_E$, $a=H_0a_{E,\rm phys}$, and $H=(\partial_x a)/a$. Velocity and extrinsic curvature vanish at launch. Time reversal gives $q(-x)=q(x)$ and $a(-x)=a(x)$ for the contracting side. The propagation through the bounce is regular. The archived fit is an interpolated numerical background. The outputs describe a test field on these fixed geometries.

The minimally coupled massless Jordan scalar has Einstein-frame action $-\frac12\int \sqrt{-g_E}Zg_E^{ab}\partial_a\chi\partial_b\chi$, where $Z=f_{\rm top}/f(q)>0$. A real orthonormal scalar harmonic on the unit three-sphere with $n\ge1$ has oscillator mass $M=a^3Z$, frequency $\omega=k/a$, $k=\sqrt{n(n+2)}$ and degeneracy $(n+1)^2$. The volume is $V=2\pi^2a^3$. Summing over a full shell gives a homogeneous density. The $n=0$ harmonic is a free particle. Its identical, uncorrelated state in all preparations cancels from every reported difference. We also run a constant-$Z$ control on each rolling metric.

The primary observer measures the Einstein-frame physical energy density. The even sampling probability in proper time has density $w_T(x)\propto\exp[-1/(1-(x/T)^2)]$ inside $|x|<T$ and vanishes outside. The half-widths are $T=0.10,0.20,0.30$. We minimize $\int w_T \rho_E\,dx$. Time reflection makes the integrated mixed position-momentum coefficient vanish. A mixed Gaussian or additional coherent excitation cannot lower the positive quadratic energy. Its minimum has a pure squeezed covariance at the bounce, with $2M\omega\langle Q^2\rangle=e^u$, $2\langle P^2\rangle/(M\omega)=e^{-u}$, and zero symmetrized position-momentum covariance. The covariance determinant is $1/4$ in units $\hbar=1$.

A second observer's clock follows minimally coupled matter. We choose its constant units to agree with the Einstein clock at launch, $g_{J,0}=Bg_E$ with $B=Z/Z(0)$. Then $d\tau_J=\sqrt B\,dx$ and $\rho_{J,0}=B^{-2}\rho_E$ for this convention and scalar normalization. The corresponding $T=0.20$ minimization uses the weight $w_{.20}(\tau_J)B^{-3/2}dx$ multiplying $\rho_E$. The JSON compares this preparation to the Einstein $T=0.20$ state on the same events, and reports density differences measured by each observer. The archived gravitational density fractions always use the Einstein-frame parent density as denominator.

## Calculation

Factor out the oscillatory phase in the instantaneous ladder variable as $b(x)=e^{-ik\eta}(A b_0+B_1 b_0^\dagger)$ with $d\eta/dx=1/a$. For $s=\partial_x\ln Z$ the interaction-picture equations are

$$
A'=(H+s/2)e^{2ik\eta} B_1^*,\qquad
B_1'=(H+s/2)e^{2ik\eta} A^*,
\qquad A(0)=1,\quad B_1(0)=0.
$$

The subscript on $B_1$ distinguishes the oscillator mixing amplitude from the metric ratio $B$. Set $C=1+2|B_1|^2$ and $D=2\operatorname{Re}(AB_1^*)$. The mean energy per real oscillator is $E=\omega(C\cosh u+D\sinh u)/2$. The positive smeared quadratic has its minimum at

$$
u_T=-\operatorname{artanh}\frac{\int w_T a^{-4} D\,dx}{\int w_T a^{-4} C\,dx}.
$$

The factor $a^{-4}$ includes the physical volume and the instantaneous frequency. Substituting the matter weight gives the matter-clock minimizer. Evenness allows both integrals to be computed on the positive half. The program directly computes differences of hyperbolic functions in the two preparations and avoids subtracting divergent bare vacuum sums.

For one oscillator let $I=V_{\rm osc}-K_{\rm osc}$ denote potential energy minus kinetic energy. Then $p=(E-2I)/(3V)$, and the work per unit dimensionless time is split as

$$
E'=sI-H(E-2I).
$$

The first term comes from the changing matter coefficient and the second from the geometry. The same identity holds for differences and sums over degeneracies. The reported work is integrated over the whole spatial sphere. A constant $Z$ gives zero direct coupling work at every time. The conversion to a density fraction uses $\Delta\rho/\rho_{\rm parent}=(\Delta\rho/H_0^4)(H_0/M_{\rm Pl})^2/[3(U+q'^2/6)]$, where $U$ is the dimensionless potential used by the background solver.

## Numerical checks and scope

The program compares the reconstructed solution to the archived crossing, tests the background interpolation at midpoints, checks $|A|^2-|B_1|^2=1$, verifies both signs of the smeared variational inequality, and records the two work integrals and their residual on two time intervals. In unit-curvature-scale exact closed de Sitter space it also compares the integrated result against direct smearing of the analytic modes
$\chi_n=(L\cos\eta+i\sin\eta)e^{-iL\eta}/\sqrt{2Lk^2}$ with $L=n+1$ and $\eta=\arctan(\sinh x)$. The different expressions for the bounce covariance agree to better than $10^{-12}$ for all 15 tested shell/window combinations.

The JSON records sharp cutoffs at several shell labels. Pressure weights oscillating kinetic and gradient energies with opposite signs and is substantially more cutoff-sensitive than density near passage. Exact low-energy states on smooth spacetimes motivate the construction. The finite numerics on this archived interpolant provide finite-cutoff stress differences. An absolute stress tensor would also require a specified local renormalization prescription. The differences alone do not give a particle yield, an absolute transferred energy, a temperature, or evidence for the astronomical origin of our universe.

## First minimum on the rolling backgrounds

All entries use complete shells through $n=2048$. Each cell subtracts the Einstein $T=0.20$ preparation on the same background. Density fractions use the parent density at first passage. The direct work integrates from the bounce to that event in units of $H_0$. A positive work difference means one preparation received more direct work than the other over this interval.

| $\phi_m$ | Preparation | $\Delta\rho_E/\rho_{\rm parent}$ | $\Delta\rho_{J,0}/H_0^4$ | Direct work difference |
|---|---|---:|---:|---:|
| +0.5 | Einstein 0.10 | $+2.428\times10^{-8}$ | $+1.190\times10^{-3}$ | $-0.01175$ |
| −0.5 | Einstein 0.10 | $+4.006\times10^{-8}$ | $+1.071\times10^{-3}$ | $+0.006295$ |
| +0.5 | Einstein 0.30 | $-1.089\times10^{-9}$ | $-5.339\times10^{-5}$ | $+0.006255$ |
| −0.5 | Einstein 0.30 | $-1.499\times10^{-8}$ | $-4.005\times10^{-4}$ | $-0.003641$ |
| +0.5 | Matter 0.20 | $+9.67\times10^{-12}$ | $+4.74\times10^{-7}$ | $+1.96\times10^{-5}$ |
| −0.5 | Matter 0.20 | $+7.42\times10^{-12}$ | $+1.98\times10^{-7}$ | $+6.08\times10^{-6}$ |

For the narrow Einstein preparation, constant-$Z$ controls give density fractions $+3.671\times10^{-8}$ and $+3.442\times10^{-8}$ in branch order. Both have zero direct work. Thus a difference of this size can also arise from preparation and expansion alone. Among these six inverse-$f$ comparisons, increasing the cutoff from 1536 to 2048 changes the first-passage fractions by at most $5.1\times10^{-14}$. This observed cutoff difference is not a bound on the infinite tail. In the narrow preparation the first-passage pressure differences in units of $H_0^4$ are approximately $-1.52\times10^{-4}$ and $+2.02\times10^{-4}$. The last cutoff increment changes these values by approximately $0.8\times10^{-6}$ and $1.7\times10^{-6}$. Pressure should be interpreted at the corresponding precision.

For a resolution check one can run `python round/sle.py --nmax 512 --order 320 --rtol 3e-13 --block 72 --output /tmp/round-refined.json` and compare at the same cutoff. In our check against the order-416, tolerance-$2\times10^{-11}$ integration, the maximum absolute change across both couplings, branches, sample times and preparations was $7.4\times10^{-12}$ in $\Delta\rho/H_0^4$ and $2.2\times10^{-8}$ in $\Delta p/H_0^4$. The JSON keeps extra digits for numerical replay.
