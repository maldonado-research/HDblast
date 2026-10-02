# Archived v1.2.1-era equation register — historical source reconstruction

This companion preserves the larger **EA-001 through EA-097** register from the FTL-only September 8 handoff. It is a historical reconstruction and is **not the corrected current implementation**. These identifiers form a different namespace from the detailed `EA-SR-*` / `EA-PCI-*` register.

Read the [detailed account and corrections](RESEARCH_ACCOUNT_AND_EQUATIONS.md) first. In particular, keep occupancy separate from first-moment severity, qualify cosmology normalization, and do not interpret classifier cross-fit stability as proof of overlap accuracy. Any historical use of “current,” “primary,” or “strict” in the source describes its source-era workflow; the later reported audit defects remain applicable.

No executable kit or observational posterior bundle is included. Mathematical formulas are conditional on their accompanying definitions. Missing proof steps are not supplied by this publication. GitHub math delimiters are a formatting adaptation only.

---


The following registry preserves substantive equations available in the conversation and FTL files. “Status” is summarized in prose after each equation. Units are dimensionless unless noted.

### EA-001 — FTL speed definition

$$
u=Wc,\quad W>1,\quad \delta=W-1
$$

Proposed notation in this research; W dimensionless; u speed; c speed of light; delta dimensionless. Used throughout.

### EA-002 — Relative-speed domain

$$
\beta\in[0,1)
$$

Standard SR notation; beta=v/c dimensionless.

### EA-003 — SR paradox wedge boundary

$$
\beta_c(W)=\frac{2W}{1+W^2}
$$

Proposed/adapted antitelephone boundary used by the audit; paradox-capable region is beta > beta_c(W).

### EA-004 — Paradox strength

$$
P(\beta,W)=1-\frac{W^2(1-\beta^2)}{(W-\beta)^2}\quad(\beta>\beta_c(W));\qquad P=0\ \text{otherwise}.
$$

Proposed audit functional based on exact SR wedge geometry. Dimensionless. Valid for W>1 and beta in [0,1).

### EA-005 — Alternative paradox-strength form

$$
P(\beta,W)=\frac{\beta[(1+W^2)\beta-2W]}{(W-\beta)^2}
$$

Algebraically equivalent inside the wedge; useful for beta monotonicity.

### EA-006 — Paradox occupancy

**Historical labeling error:** the displayed integral is a first positive-strength moment under the clipped P convention in EA-004. It is not wedge occupancy. The source heading and equation are preserved; use Q_occ and M_1 as distinguished in the corrected guide.

$$
Q[\rho]=\int P(\beta,W)\rho(\beta,W)\,d\beta\,dW
$$

Proposed population functional. Dimensionless; rho is a probability density over configurations.

### EA-007 — SR severity functional

$$
J[\rho]=\int P(\beta,W)^2\rho(\beta,W)\,d\beta\,dW
$$

Central J-first SR control variable. Dimensionless; current preferred SR object.

### EA-008 — Small-FTL variance proxy

$$
\varepsilon_{\rm FTL}^2=\langle(W-1)^2\rangle=\langle\delta^2\rangle
$$

Proxy only, not primary outside small-FTL regime.

### EA-009 — Small-FTL Q approximation

$$
Q\simeq \frac12\varepsilon_{\rm FTL}^2
$$

Toy/proxy relation; source record gives finite validity bands.

### EA-010 — Small-FTL J approximation

$$
J\simeq \frac16\varepsilon_{\rm FTL}^2
$$

Proxy relation; superseded as primary by J-first formulation.

### EA-011 — Hidden SGWB SBPL spectrum

$$
\Omega_{\rm gw}(f)=\Omega_k\left[\frac{(f/f_k)^{\alpha_1\Delta}+(f/f_k)^{\alpha_2\Delta}}{2}\right]^{-1/\Delta}
$$

Modeling ansatz for hidden stochastic gravitational-wave background. Frequencies in Hz; Omega dimensionless energy density per log frequency.

### EA-012 — Extra-radiation mapping

$$
\Delta N_{\rm eff}=\mathcal K\int \Omega_{\rm gw}(f)\,d\ln f
$$

Standard-style cosmology compression used by the audit; integral over log-frequency support.

### EA-013 — Radiation conversion constant

$$
\mathcal K=\frac87\left(\frac{11}{4}\right)^{4/3}
$$

Dimensionless conversion factor used in the package.

### EA-014 — Energy normalization

$$
\varepsilon_{\rm energy}=\frac{\Delta N_{\rm eff}}{0.3}
$$

Guardrail normalization; 0.3 is a conservative approximate cap used in the project.

### EA-015 — PTA axis

$$
y=1-\mathrm{PCI}
$$

PCI is posterior consistency/overlap in [0,1]; y is disagreement axis.

### EA-016 — Ringdown axis

$$
z=\frac{\varepsilon_{\rm pooled}}{\sigma_\varepsilon}
$$

Dimensionless pooled ringdown deviation; assumes harmonized epsilon definition.

### EA-017 — Einstein vector

$$
X=(x,y,z)
$$

Core diagnostic vector.

### EA-018 — Einstein distance

$$
D_E=\|X\|=\sqrt{x^2+y^2+z^2}
$$

Dimensionless scalar diagnostic.

### EA-019 — Einstein allowed unit ball

$$
D_E\le 1
$$

Audit convention. Failure indicates the proposed configuration exceeds unit budget.

### EA-020 — Data floor

$$
D_{\min}=\sqrt{y^2+z^2}
$$

Observational contribution even when x=0.

### EA-021 — Canonical radicand

$$
R=1-y^2-z^2=1-(1-\mathrm{PCI})^2-z^2
$$

v1.2.1 canonical observational object.

### EA-022 — Slack/headroom

$$
x_{\max}=\sqrt{R}\quad(R\ge0)
$$

Valid only for nonnegative R; invalid mass reported separately.

### EA-023 — Causality-cosmology action

$$
S=\varepsilon_{\rm energy}J
$$

Intermediate action-like audit quantity.

### EA-024 — Action scale

$$
S_0=\frac16
$$

Chosen so small-FTL mapping reduces to previous proxy form.

### EA-025 — J-first x mapping

$$
x=\frac{S}{S_0}=6\varepsilon_{\rm energy}J=20\Delta N_{\rm eff}J
$$

Current preferred x mapping.

### EA-026 — Small-FTL x proxy

$$
x\simeq \varepsilon_{\rm energy}\varepsilon_{\rm FTL}^2
$$

Superseded as primary; only valid in small-FTL proxy regime.

### EA-027 — Envelope bound

$$
J\le J_{\rm bound}(\Delta N_{\rm eff})=\frac{x_{\max}}{20\Delta N_{\rm eff}}
$$

Current J-first no-free-FTL-lunch envelope.

### EA-028 — Original epsilon bound

$$
\varepsilon_{\rm FTL}^2\le \frac{\sqrt{1-(1-\mathrm{PCI})^2-(\varepsilon_{\rm pooled}/\sigma_\varepsilon)^2}}{\varepsilon_{\rm energy}}
$$

Historical proxy bound; superseded by J-first except where small-FTL proxy is valid.

### EA-029 — Large-W paradox limit

$$
P(\beta,W)\to \beta^2\quad(W\to\infty)
$$

SR limiting case supporting Testability Gate.

### EA-030 — Severity ceiling

$$
J_{\max}=E[\beta^4]
$$

Asymptotic maximum severity for a beta population.

### EA-031 — Existence gate condition

$$
J_{\max}>J_{\rm bound}(\Delta N_{\rm eff})
$$

Finite bound exists only under this inequality.

### EA-032 — Delta N_eff threshold

$$
\Delta N_{\rm eff}^{\rm thr}=\frac{x_{\max}}{20J_{\max}}
$$

Testability Gate threshold.

### EA-033 — Beta(a,b) fourth moment

$$
E[\beta^4]=\frac{a(a+1)(a+2)(a+3)}{(a+b)(a+b+1)(a+b+2)(a+b+3)}
$$

Standard beta-distribution moment used for kinematic models.

### EA-034 — Tail-fraction diagnostic

$$
f_{\rm rel}\,\beta_{\rm rel}^4\gtrsim \frac{x_{\max}}{20\Delta N_{\rm eff}}
$$

Approximate requirement for a relativistic subpopulation to reach testability.

### EA-035 — Einstein budget fractions

$$
f_{\rm FTL}=\frac{x^2}{D_E^2},\quad f_{\rm PTA}=\frac{y^2}{D_E^2},\quad f_{\rm RD}=\frac{z^2}{D_E^2}
$$

Diagnostic partition; fractions sum to one when D_E>0.

### EA-036 — Covariance-weighted tension

$$
\mathcal T_E=X^\top C^{-1}X,\quad D_E^{({\rm cov})}=\sqrt{\mathcal T_E}
$$

A2 statistical upgrade once joint covariance is available.

### EA-037 — Exact finite-W severity curve

$$
J(W)=\int_{\beta_c(W)}^1 P(\beta,W)^2 f_\beta(\beta)\,d\beta
$$

Finite-speed SR severity for a chosen beta distribution.

### EA-038 — W monotonicity of P

$$
\frac{\partial P}{\partial W}=\frac{2W\beta(1-\beta^2)}{(W-\beta)^3}>0
$$

Derived for 0<beta<1, W>1 inside wedge.

### EA-039 — W monotonicity of wedge boundary

$$
\frac{d\beta_c}{dW}=\frac{2(1-W^2)}{(1+W^2)^2}<0\quad(W>1)
$$

Wedge expands as W grows.

### EA-040 — Exact critical W

$$
J(W_*)=J_{\rm bound}(\Delta N_{\rm eff})
$$

Unique when 0<J_bound<J_max.

### EA-041 — Practical finite-W threshold

$$
\Delta N_{\rm eff}^{\rm practical,exact}(W_{\rm cap})=\frac{x_{\max}}{20J(W_{\rm cap})}
$$

Threshold for reachability below a chosen W cap.

### EA-042 — Beta moment general

$$
E[\beta^k]=\frac{(a)_k}{(a+b)_k}
$$

Standard rising-factorial moment.

### EA-043 — Large-W coefficient

$$
c_1=4(E[\beta^3]-E[\beta^5])
$$

Controls first approach to J_max.

### EA-044 — Beta c1 closed form

$$
c_1=\frac{4ab(a+1)(a+2)(2a+b+7)}{(a+b)(a+b+1)(a+b+2)(a+b+3)(a+b+4)}
$$

Beta-family expression.

### EA-045 — Stiffness parameter

$$
\Xi=\frac{c_1}{J_{\max}}=\frac{4b(2a+b+7)}{(a+3)(a+b+4)}
$$

Finite-W stiffness; larger means slower saturation.

### EA-046 — Practical/existence ratio

$$
\frac{\Delta N_{\rm eff}^{\rm practical}(W_{\rm cap})}{\Delta N_{\rm eff}^{\rm thr}}=1+\frac{\Xi}{W_{\rm cap}}+O(W_{\rm cap}^{-2})
$$

Large-W approximation.

### EA-047 — Finite-W asymptotic

$$
J(W)=J_{\max}-\frac{c_1}{W}+O(W^{-2})
$$

Large-W expansion.

### EA-048 — Asymptotic W star

$$
W_*\approx \frac{c_1}{J_{\max}-J_{\rm bound}}
$$

Interpretive approximation near the gate.

### EA-049 — Small-FTL wedge width

$$
1-\beta_c(1+\delta)=\frac{\delta^2}{2}+O(\delta^3)
$$

Boundary-layer result for delta<<1.

### EA-050 — Tail density

$$
f(\beta)\sim c(1-\beta)^{b-1}\quad(\beta\to1^-)
$$

Assumption for small-FTL Tail Theorem.

### EA-051 — Small-FTL Q tail scaling

$$
Q(\delta)\simeq \frac{c}{b}\left(\frac{\delta^2}{2}\right)^b
$$

Asymptotic conditional on tail model.

### EA-052 — Small-FTL J tail scaling

$$
J(\delta)\simeq c\left(\frac{\delta^2}{2}\right)^b\frac{2}{b(b+1)(b+2)}
$$

Asymptotic conditional on tail model.

### EA-053 — Tail ratio

$$
\frac{J}{Q}\simeq \frac{2}{(b+1)(b+2)}
$$

Universal ratio under tail model.

### EA-054 — Posterior phase thresholds

$$
t_{\rm exist}=20\Delta J_{\max},\quad t_{\rm practical}=20\Delta J(W_{\rm cap})
$$

Thresholds in x_max space.

### EA-055 — Phase probabilities in X

$$
p_{\rm closed}=P[X>t_{\rm exist}],\quad p_{\rm practical}=P[X\le t_{\rm practical}],\quad p_{\rm latent}=P[t_{\rm practical}<X\le t_{\rm exist}]
$$

Posterior phase theorem before radicand refinement.

### EA-056 — Quantile frontiers

$$
\Delta_{\rm exist}^{(p)}=\frac{Q_X(p)}{20J_{\max}},\quad \Delta_{\rm practical}^{(p)}=\frac{Q_X(p)}{20J(W_{\rm cap})}
$$

Probability-calibrated thresholds.

### EA-057 — Finite-W overhead

$$
\Lambda(W_{\rm cap})=\frac{J_{\max}}{J(W_{\rm cap})}\ge1
$$

Pure-SR multiplicative practicality overhead.

### EA-058 — Lambda factorization

$$
\Delta_{\rm practical}^{(p)}=\Lambda(W_{\rm cap})\Delta_{\rm exist}^{(p)}
$$

A2 data enters scale; SR finite-W enters overhead.

### EA-059 — Beta-class rectangle

$$
\mathcal R=[a_-,a_+]\times[b_-,b_+]
$$

Family of beta kinematic priors.

### EA-060 — Beta derivative of P

$$
\frac{\partial P}{\partial\beta}=\frac{2W^2(W\beta-1)}{(W-\beta)^3}>0
$$

Inside wedge; supports corner extrema.

### EA-061 — Kinematic envelope existence interval

$$
\Delta_{\rm exist}^{(p)}(\mathcal R)\in\left[\frac{Q_X(p)}{20J_{\max}(a_+,b_-)},\frac{Q_X(p)}{20J_{\max}(a_-,b_+)}\right]
$$

Beta-class envelope.

### EA-062 — Kinematic ambiguity factors

$$
A_{\rm exist}=\frac{J_{\max}(a_+,b_-)}{J_{\max}(a_-,b_+)},\quad A_{\rm practical}=\frac{J(W_{\rm cap};a_+,b_-)}{J(W_{\rm cap};a_-,b_+)}
$$

Separates data uncertainty from kinematic prior uncertainty.

### EA-063 — Moment-only lower envelope

$$
\underline J_W(m)=\operatorname{co}[g_W](m),\quad m=E[\beta^4]
$$

Family-free lower bound; co is lower convex envelope.

### EA-064 — Moment-only severity interval

$$
\underline J_W(m)\le J(W)\le m
$$

Sharp given only m=J_max.

### EA-065 — Zero lower-envelope condition

$$
\underline J_W(m)=0\quad\text{whenever}\quad m\le \beta_c(W)^4
$$

Finite-W practicality unresolved from m alone.

### EA-066 — Moment-only practical interval

$$
\Delta_{\rm practical}(X;m,W)\in\left[\frac{X}{20m},\frac{X}{20\underline J_W(m)}\right]
$$

Conservative family-free practical threshold interval.

### EA-067 — PCI overlap

$$
\mathrm{PCI}(p,q)=\int \min(p,q)\,d\theta
$$

Overlap coefficient.

### EA-068 — PCI-TV identity

$$
y=1-\mathrm{PCI}=\mathrm{TV}(p,q)
$$

Statistical interpretation of PTA disagreement.

### EA-069 — PCI ladder

$$
\mathrm{PCI}^{all}_{K,2D}\le\mathrm{PCI}^{min\ pair}_{2D}\le\mathrm{PCI}^{mean\ pair}_{2D}\le\mathrm{PCI}^{mean\ pair}_{1D}
$$

Primary/sensitivity hierarchy.

### EA-070 — Generic threshold map in PCI

$$
\Theta_K(\pi,z)=\frac{\sqrt{2\pi-\pi^2-z^2}}{20K}
$$

Used for PTA optimism and estimator taxes.

### EA-071 — PCI threshold derivative

$$
\frac{\partial\Theta_K}{\partial\pi}=\frac{1-\pi}{20K\sqrt{2\pi-\pi^2-z^2}}>0
$$

Threshold increases with PCI in valid domain.

### EA-072 — PCI threshold concavity

$$
\frac{\partial^2\Theta_K}{\partial\pi^2}=-\frac{1-z^2}{20K(2\pi-\pi^2-z^2)^{3/2}}<0
$$

Concavity in PCI.

### EA-073 — Universal PTA optimism tax

$$
\Omega_{\rm PTA}=\sqrt{\frac{2\pi_2-\pi_2^2-z^2}{2\pi_1-\pi_1^2-z^2}}
$$

K cancels; same tax for all SR thresholds.

### EA-074 — Threshold factorization

$$
\frac{\Theta}{\Theta_0}=\frac{X}{X_0}\frac{K_0}{K}
$$

Audit budget factorization.

### EA-075 — Log-budget share

$$
s_i=\frac{\ln\Omega_i}{\ln\Omega}
$$

Multiplicative uncertainty-budget partition.

### EA-076 — Posterior slack in Pi/Z

$$
X=\sqrt{2\Pi-\Pi^2-Z^2}
$$

A2 slack with PCI random variable Pi and ringdown random variable Z.

### EA-077 — Point-estimate optimism ordering

$$
E[X]\le X_{2\rm mom}\le X_{\rm mean}
$$

Jensen/curvature consequence in v3.0 note.

### EA-078 — Second-order propagation

$$
E[X]\approx X_\mu-\frac{\mathrm{tr}\Sigma}{2X_\mu}-\frac{v^\top\Sigma v}{2X_\mu^3},\quad v=(1-\mu_\Pi,-\mu_Z)^\top
$$

Local covariance propagation.

### EA-079 — Central leverage shares

$$
s_{\rm PTA}=\frac{y^2}{y^2+z^2},\quad s_{\rm RD}=\frac{z^2}{y^2+z^2}
$$

Which observational central value steers thresholds.

### EA-080 — Observational polar variables

$$
D=\sqrt{y^2+z^2},\quad \phi=\arctan(z/y),\quad x_{\max}=\sqrt{1-D^2}
$$

Polar observational geometry.

### EA-081 — Polar leverage

$$
s_{\rm PTA}=\cos^2\phi,\quad s_{\rm RD}=\sin^2\phi
$$

Angle controls PTA/RD burden.

### EA-082 — Polar posterior tax

$$
E[x_{\max}]\approx x_{\max}(\mu)-\frac{\sigma_r^2}{2x_{\max}(\mu)^3}-\frac{\sigma_t^2}{2x_{\max}(\mu)}
$$

Radial uncertainty is amplified near D -> 1.

### EA-083 — Radicand invalid mass

$$
p_{\rm invalid}=P(R<0)
$$

Canonical validity gate.

### EA-084 — Radicand practical phase

$$
p_{\rm practical}(\Delta)=P(0\le R\le(20\Delta K_p)^2)
$$

Validity-gated practical mass.

### EA-085 — Radicand latent phase

$$
p_{\rm latent}(\Delta)=P((20\Delta K_p)^2<R\le(20\Delta K_e)^2)
$$

Validity-gated latent mass.

### EA-086 — Radicand closed phase

$$
p_{\rm closed}(\Delta)=P(R>(20\Delta K_e)^2)
$$

Validity-gated closed mass.

### EA-087 — Classifier mixture

$$
m(\theta)=K^{-1}\sum_i p_i(\theta)
$$

Equal mixture of PTA posteriors.

### EA-088 — Classifier posterior ratio

$$
\eta_i(\theta)=P(i\mid\theta)=\frac{p_i(\theta)}{\sum_j p_j(\theta)}
$$

Density-ratio identity.

### EA-089 — K-way classifier PCI

$$
\mathrm{PCI}_{all}=K E_{\theta\sim m}[\min_i\eta_i(\theta)]
$$

Primary classifier PCI identity.

### EA-090 — Pairwise classifier PCI

$$
\mathrm{PCI}_{ij}=2E_{m_{ij}}[\min(\eta_i,\eta_j)]
$$

Pairwise overlap identity.

**Pairwise normalization clarification (2026-09-30 editorial review):** the probabilities in EA-090 must be recomputed for the two-class mixture $m_{ij}=(p_i+p_j)/2$: $\eta_i^{(ij)}=p_i/(p_i+p_j)$ and $\eta_j^{(ij)}=p_j/(p_i+p_j)$. The identity is $\mathrm{PCI}_{ij}=2\mathbb E_{m_{ij}}[\min(\eta_i^{(ij)},\eta_j^{(ij)})]$. For $K>2$, substituting the K-class probabilities defined in EA-088 does not in general recover pairwise overlap. The historical display above is preserved; this note states the required convention.

### EA-091 — Cross-fit PCI estimator

$$
\widehat{\mathrm{PCI}}^{CF}_{all}=K\frac1n\sum_r\min_i\widehat\eta_i^{(-s(r))}(\theta_r)
$$

Out-of-fold estimator to reduce leakage.

### EA-092 — Leakage diagnostic

$$
\Delta_{\rm leak}=\widehat{\mathrm{PCI}}^{IN}_{all}-\widehat{\mathrm{PCI}}^{CF}_{all}
$$

In-sample minus cross-fit diagnostic.

### EA-093 — Ringdown fixed-effect mean

$$
\mu_{\rm FE}=\frac{\sum_iw_i\mu_i}{\sum_iw_i}
$$

Inverse-variance pooling; assumes harmonized epsilon.

### EA-094 — Ringdown fixed-effect sigma

$$
\sigma_{\rm FE}=\left(\sum_iw_i\right)^{-1/2}
$$

Pooled standard error.

### EA-095 — Ringdown fixed-effect z

$$
z_{\rm FE}=\mu_{\rm FE}/\sigma_{\rm FE}
$$

Ringdown contribution to Einstein vector.

### EA-096 — Cochran heterogeneity Q

$$
Q=\sum_iw_i(\mu_i-\mu_{\rm FE})^2
$$

Ringdown heterogeneity diagnostic.

### EA-097 — I squared

$$
I^2=\max\left(0,\frac{Q-(N-1)}{Q}\right)
$$

Approximate heterogeneity fraction.


### A3.98 Corrected, superseded, or restricted formulas

- The small-FTL `epsilon_FTL^2` proxy is retained only as a restricted approximation. It is not the primary current SR severity measure.
- The original epsilon-bound form is retained historically but should not be used outside validated small-FTL regimes.
- Knee-frequency PTA overlap is retained as screening/triage only; it is not strict headline PCI.
- Histogram PCI is retained as stress-test/sensitivity only in the frozen strict path; calibrated cross-fit classifier PCI is the primary strict estimator.
- Clipped `sqrt(max(0,R))` is not a headline validity-gated phase object when `P(R<0)` is nontrivial.

### A3.99 Derivation links and limiting cases

Key limiting cases:

- `W -> 1+`: `J(W) -> 0`; no SR wedge severity in the luminal limit.
- `W -> infinity`: `P -> beta^2`, so `J -> E[beta^4]`.
- `R < 0`: observational compression invalidates the remaining slack; invalid mass must be reported.
- `Delta N_eff < Delta N_eff^thr`: the gate is closed under the stated beta model.
- `Delta N_eff >= Delta N_eff^practical(W_cap)`: the model is constrainable below the finite speed cap under the stated assumptions.
