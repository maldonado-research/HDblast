# Independent prehistory enclosure and physics audit

7 October 2026. This is internal independent AI-assisted mathematics, not external peer review. No physical source callback or saved quantum array was evaluated in this preparation. The candidate encloses a uniquely prescribed analytic target at nine rational momenta. It does not enclose the error of the retained binary80 incoming states or a continuum momentum integral.

## Target and source references

The actual reference producer is `research/HDBLAST_CHECKPOINT_20261002_METRIC_RESPONSE_MEMORY_FOLLOWUP/independent/forced_metric.py` in `/workspace/HDblast`, SHA256 `4fa625f2b534d741e2a48d6f0bfe84a89ee17bf306805906f962a8552e005906`.

* Lines 26–30 fix represented epsilon, initial time −6, observations, cutoffs and sources.
* Lines 163–189 define B=exp[1−1/(1−z²)], z=s+4, the two sources h=B and h=zB, their exact compact support, and g=4L²h−2Lh′−h″, L=−1/s.
* Lines 192–198 define the stable entire drift kernel Phi=(exp(2ikτ)−1)/(2ik), including its removable k=0 limit.
* Lines 253–256 identify delta_v=v0*u, delta_v′=v0*(w−ik*u), so normalized U=u/epsilon and W=w/epsilon are the variables relevant here.
* Lines 359–368 construct GL8 local-time quadrature and exact zero u,w at −6. Lines 387–395 propagate those approximations. Lines 397–400 check the sampled Wronskian residual; this is not an incoming-state error estimate.

The model derivation is `research/HDBLAST_CHECKPOINT_20261002_METRIC_RESPONSE/theory/primary/METRIC_ACTION_AND_CONTACT_DERIVATION.md`, lines 18–32: at x=r=2, H=1, xi=0, v0=exp(−ikη)/sqrt(2k) is the fixed incoming BD mode, with zero forced variation before the compact metric perturbation. This is a transported BD prescription, not a finite-time WKB or instantaneous-vacuum initialization.

Consequently, for a=−9/2,

    W(a,k)=−∫[-5,a] exp(2ik(a−s)) g(s) ds,
    U(a,k)=−∫[-5,a] Phi(a−s,k) g(s) ds.

The original rounded modes approximate this target through source, phase, GL8 and arithmetic operations. Their values are not decoded in this task. A new target enclosure alone does not establish their error.

## Exact obstruction: normalization is blind to real prehistory error

Let mu be any finite real signed measure on a bounded prehistory interval. It may represent a continuous residual, an exact real quadrature, or their difference. Replacing g(s)ds by dmu gives

    Re U=−∫sin(2kτ)/(2k)dmu,
    W=−∫[cos(2kτ)+i sin(2kτ)]dmu,
    c=Re U−Im W/(2k)=0.

Thus any real quadrature, including an arbitrarily inaccurate one or omission of the entire prehistory, has exactly zero canonical defect. After the incoming time, each resulting homogeneous difference obeys the same direct Ward identity. Those checks cannot certify the prehistory quadrature. This is a necessary limitation in the actual fixed-BD calculation, not an assertion that its prescribed state is arbitrary.

The discrete algebra is stronger: writing arbitrary E=C+iS and drift=(E−1)/(2ik), the update Unew=U+drift*W, Wnew=E*W preserves c for every real C,S. The property does not require |E|=1. A paired, nonunitary wrong phase can pass this scalar normalization invariant. The exact checker includes E=0 as an adversarial control, rejected by the actual propagation equation.

Rounded E, drift and contractions need their own arithmetic enclosure. The target c=0 theorem is not permission to project retained rounded u,w onto c=0.

## Infrared behavior and direct operators

Write A=W/(2ik), p=Re A, q=Im A. Then

    p=−∫sin(2kτ)/(2k)dmu,
    q= ∫cos(2kτ)/(2k)dmu.

U and W are entire functions of complex k for a finite measure of bounded support, by the uniformly convergent exponential series on compact k sets. A generally has a simple pole: kA→i*mu(total)/2. Therefore a uniform sigma for A is not automatic, despite exact normalization. This is a representation pole, not a divergent finite-band stress. The October 7 theorem explicitly allows weighted infrared envelopes, including O(1/k).

The k=0 probe is this analytic continuation of the normalized response kernels. The plane-wave normalization v0=exp(−ikt)/sqrt(2k) is not evaluated at k=0, and the probe is not a physical zero-mode stress evaluation.

For k>0, let M0=∫|dmu| and M1(t)=∫(t−s)|dmu|, with t at or after the support. The exact correlation gives |p|≤M1 and |q|≤M0/(2k). The direct density and pressure are

    R=3L²p/(2k)+Lq,
    P=−(2k/3+L²/(2k))p+Lq.

Their measure-weighted single-source kernels, before division by 2*Pi², are

    k²R_kernel=−3L² sin(2kτ)/4 + Lk cos(2kτ)/2,
    k²P_kernel=(k²/3+L²/4) sin(2kτ) + Lk cos(2kτ)/2.

Both are entire and vanish at k=0. With the same fixed positive Pi as the inherited model,

    |R_K|≤K²(3L²M1+LM0)/(8Pi²),
    |P_K|≤M1*K⁴/(12Pi²)+K²(L²M1+LM0)/(8Pi²).

These are useful if an independently small residual measure norm is available. The total variation of a quadrature atomic measure minus g(s)ds generally does not decrease under refinement; using that loose norm does not certify convergence. The present candidate instead encloses analytic source polynomials and their full continuous residuals.

## Flat endpoint cap, independently integrated by parts

Use x=s+5, delta=1/128, D=delta(2−delta)=255/16384, Lc=1/(5−delta), H=(3/8)^63 and p_delta=2(1−delta)/D². The flat source at x=0 has h=h′=0. Its endpoint satisfies Bdelta≤H because its exponent is −16129/255<−63 and e>8/3. B increases on this cap. For rho=0 for B and rho=1 for zB,

    |h|≤H, |h′(capend)|≤H(p_delta+rho), L≤Lc.

For E=exp(2ik(a−s)), the exact identity is

    ∫E g = [E(−h′−2Lh−2ikh)] + ∫E(4k²−4ikL+6L²)h.

For Phi=(E−1)/(2ik),

    ∫Phi g = [Phi(−h′−2Lh)−E h]
                + ∫[6L²Phi−(2ik+2L)E]h.

The brackets are endpoint differences on the cap. Since |E|=1 and |Phi|≤a−s≤1/2 for real k,

    Wcap≤H[p_delta+rho+2Lc+2k
                    +delta(4k²+4kLc+6Lc²)],
    Ucap≤H[(p_delta+rho+2Lc)/2+1
                    +delta(3Lc²+2k+2Lc)].

These are positive analytic errors from omitted source, not numerical evaluations or an assumption of zero source on the cap. The independent source omission has c=0 in exact arithmetic, but the candidate exports independent U/W boxes and does not narrow them using that correlation. At k=0 only the exactly proved imaginary-zero property is used.

## Interior analytic source and directed-rational route

Start x_left=delta and repeatedly choose x_right=min(3*x_left/2,1/2). There are eleven parents. Their center is c=(x_left+x_right)/2−5, halfwidth Hpanel=(x_right−x_left)/2 and analytic radius R=x_left/2, hence Hpanel/R≤1/2. Put

    rho_disk=1−(c+5)+R<1,
    T=−c−R>0, Ddisk=1−rho_disk².

On each complex disk, |B|<2 follows from Re[1/(1−z²)]>1/2. With B′/B=−2z/(1−z²)² and B″/B=(6z⁴−2)/(1−z²)⁴, the conservative disk bounds used here are

    MB=2[4/T²+4rho_disk/(T Ddisk²)+2/Ddisk²
                       +8rho_disk²/Ddisk³+4rho_disk²/Ddisk⁴],
    MzB=rho_disk*MB+2[2/T+4rho_disk/Ddisk²].

Cauchy's tail of degree112 is at most M*2^−112. This new geometric proof is necessary: the frozen radius1/8, M64 central proof was only established on [−9/2,−7/2] and cannot be extended to the nonanalytic flat endpoint.

`incoming_dyadic.py` independently constructs all 113 forcing coefficients from exact reciprocal/exponential formal recurrences, including source coefficients through degree114. It encloses the scalar exponential by range reduction into [−1,0], the exact alternating series of even degree200, and outward dyadic squaring. Chosen polynomial coefficients use 512-bit dyadic points; all coefficient radii, scalar-exponential enclosure errors and source tails enter a uniform parent source bound. The two historical pure-rational helpers are copied unchanged and pinned as baseline files; no primary Arb helper is imported.

Each parent polynomial is exactly translated onto child time intervals of width≤1/64. The mode solver chooses a degree160 complex polynomial with 512-bit dyadic coefficients. Its independently recomputed defect is against W′=2ikW−g and U′=W, with all source coefficients present. The exact real-frequency propagator has modulus one, so integrating the absolute polynomial defect supplies W error; the additional (h−s) weight supplies U error. There is no division by k. Coefficient rounding is included by the true-equation defect. Endpoint point shifts and upward radius rounding are added explicitly. The actual full-width maximum phase is8; the solver does not use the primary route's centered-phase series or cap4 assumption.

For each whole state, model, ODE-defect, endpoint point-rounding and radius-rounding components telescope exactly. The omitted cap is added at the final time. If d=1 at k=0 and d=2 otherwise, exported rectangles satisfy the checked exact identity

    exported L1 radius = d*(cap+model+defect+point_round+radius_round)
                         + export_excess,
    0≤export_excess<d*2^−512.

The radius is the actual sum of real/imag halfwidths, not a requested tolerance. The candidate keeps all18 source×momentum cases and both U/W quantities. Root's prospective gate is 1e−20 for these36 exported radii. This is a new narrow target-enclosure gate and does not change the inherited2e−8 full-integral pressure/contact gate.

API: `run(auth)` returns `(data,budget)` and requests root authorization before each of22 actual source constructions. `run_fabricated()` constructs only explicit degree112 rational polynomial fixtures on the same schedule. Planned hard limits are900seconds,512MiB RSS and20MiB serialized output per route; two complete fabricated rehearsals must pass before physical evaluation. Actual source evaluation additionally requires prospective public freeze, exact readback and independent review through the root wrapper.

## Positivity, ultraviolet behavior, and scientific limits

Normalization, positivity and Hadamard admissibility would still not select a state by themselves. For example a smooth momentum bump beta(k) supported strictly between sampled momenta and away from k=0 can define v_epsilon=sqrt(1+epsilon²|beta|²)v0+epsilon*beta*v0*. Its Wronskian is exactly one. It gives a positive quasifree mode construction and changes the two-point function smoothly, preserving the Hadamard singularity if the reference state is Hadamard. It leaves all sampled modes untouched and yields arbitrary first-order beta. In finite spatial volume with compact momentum support only finitely many oscillator squeezes are involved, so they are unitarily implementable. In infinite homogeneous volume a finite particle density does not by itself establish a global Fock-unitary map; no such claim is needed here. The fixed BD pre-source prescription excludes this alternative, which is why enforcing that prescription analytically matters.

The present nine-probe target enclosure and continuous source-model error do not prove all-k accuracy of the mode approximation, an ultraviolet tail, a nonlinear positive state beyond the declared perturbation order, or the retained native implementation's error. At finite epsilon the exact normalized family also has quadratic terms that a first-order stress cannot silently omit. The full twelve-case certificate remains UNRESOLVED, original rounded incoming-state error NOT_ENCLOSED, metric calibration FAIL, higher-dimensional Big Bang origin NOT_ESTABLISHED and external novelty NOT_ASSESSED.
