# Core mathematics of the FTL Einstein Audit

This guide summarizes conditional mathematics preserved in the September 8, 2026 source snapshot. It is a navigation aid, not a new derivation claimed as a discovery. The full source register and its corrections are in [RESEARCH_ACCOUNT_AND_EQUATIONS.md](RESEARCH_ACCOUNT_AND_EQUATIONS.md). The `EA-*` identifiers are editorial references introduced in that handoff.

## 1. Geometry: when the reply precedes the message

Alice and Bob coincide at time zero. In Alice's frame, Bob moves away at speed $v$, and Alice sends at $T>0$. The outgoing hypothetical signal has speed $u>c$:

$$
x=u(t-T),\qquad x_B=vt.
$$

Their intersection gives $t_1=uT/(u-v)$ and $x_1=vt_1$. With $\gamma=(1-v^2/c^2)^{-1/2}$, Bob's reception time is $t'_1=t_1/\gamma$. He replies immediately at speed $u$ toward Alice in his own frame. Alice follows $x'_A=-vt'$, and the reply follows $x'=-u(t'-t'_1)$, giving $t'_2=ut'_1/(u-v)$ and $t_2=t'_2/\gamma$.

Writing $W=u/c>1$ and $\beta=v/c\in[0,1)$ therefore gives **EA-SR-003**:

$$
\frac{t_2}{T}=\frac{W^2(1-\beta^2)}{(W-\beta)^2}.
$$

For positive $\beta$, comparing $t_2$ with $T$ gives **EA-SR-004**:

$$
t_2<T
\iff\beta\big[(1+W^2)\beta-2W\big]>0
\iff\beta>\beta_c(W)=\frac{2W}{1+W^2}.
$$

This assumes flat spacetime, collinear inertial motion, immediate response, and equal hypothetical FTL speeds in the respective senders' frames. A preferred-frame rule, delayed reply, or different propagation law requires its own calculation. The construction does not establish a physically available communication channel.

## 2. Population severity and the Q notation correction

Define **EA-SR-005** and **EA-POP-001–004**:

$$
P=1-\frac{t_2}{T},\qquad P_+=\max(P,0),
$$

$$
Q_{\rm occ}=\mathbb E[\mathbf 1_{P>0}],\qquad
M_1=\mathbb E[P_+],\qquad
J=\mathbb E[P_+^2].
$$

The expectation uses a specified normalized population density $\rho(\beta,W)$. These quantities are dimensionless. Occupancy counts configurations; $M_1$ weights their positive strength; $J$ weights squared strength. For $Q_{\rm occ}>0$, $\sqrt{J/Q_{\rm occ}}$ is the conditional RMS positive strength.

Historical sources use **Q** for both occupancy and the first moment. This guide uses separate symbols as an editorial disambiguation, not an assertion that the original record used consistent notation.

For uniform $\beta$ and fixed $W=1+\delta$, the exact occupancy is

$$
Q_{\rm occ}=1-\beta_c(W)=\frac{\delta^2}{2+2\delta+\delta^2}.
$$

Its leading coefficient is $1/2$. The later report of a first-moment leading coefficient $1/4$ concerns a different quantity; it cannot be described as correcting this exact occupancy formula. The private release-specific derivations remain unavailable.

## 3. Severity ceiling and the strict Testability Gate

For a fixed observer-speed distribution, $P\to\beta^2$ as $W\to\infty$, so **EA-GATE-001** gives

$$
J_{\max}=\mathbb E[\beta^4].
$$

For a Beta distribution with $a,b>0$, the standard moment formula applied here yields

$$
J_{\max}=\frac{a(a+1)(a+2)(a+3)}{(a+b)(a+b+1)(a+b+2)(a+b+3)}.
$$

The uniform case $a=b=1$ gives $1/5$. This is a population assumption, not a measured distribution of a hidden sector.

If an allowance is $J_{\rm bound}$, a finite crossing in the stated increasing-severity model requires $J_{\rm bound}<J_{\max}$. Equality only reaches the infinite-speed asymptote. A finite practical cap $W_{\rm cap}$ creates a stronger requirement involving $J(W_{\rm cap})$. Neither boundary is a universal test of every possible FTL theory.

## 4. Proposed audit geometry

The author chooses an archived extra-radiation guardrail $\Delta N_{\rm eff}^{\max}=0.3$ and energy fraction $\varepsilon_{\rm energy}=\Delta N_{\rm eff}/0.3$. The proposed coupling is

$$
x=6\varepsilon_{\rm energy}J=20\Delta N_{\rm eff}J.
$$

With $y=1-\mathrm{PCI}$ and a harmonized ringdown deviation $z=\varepsilon_{\rm pooled}/\sigma_\varepsilon$, the audit adopts

$$
X=(x,y,z),\qquad D_E=\sqrt{x^2+y^2+z^2}\le1.
$$

The coordinates are dimensionless but do not share an automatically justified statistical scale. This is an author-defined construction. Its unit boundary is not a documented confidence level, and the physical proportionality in $x$ is not derived from an FTL-sector action.

For $R=1-y^2-z^2\ge0$ and $\Delta N_{\rm eff}>0$,

$$
x_{\max}=\sqrt R,\qquad J_{\rm bound}=\frac{\sqrt R}{20\Delta N_{\rm eff}}.
$$

Invalid $R<0$ draws must remain separately counted; $\Delta N_{\rm eff}=0$ is a separate boundary case. Dependence between uncertain PCI and ringdown values affects the distribution of $R$.

## 5. Posterior overlap: an exact identity versus an estimated answer

For compatible normalized densities $p_1,\ldots,p_K$ over a common parameter space,

$$
\mathrm{PCI}_{\rm all}=\int\min_i p_i(\theta)\,d\theta.
$$

Under the equal-prior mixture $m=K^{-1}\sum_i p_i$ and true Bayes class probabilities $\eta_i=p_i/\sum_jp_j$,

$$
\mathrm{PCI}_{\rm all}=K\,\mathbb E_{\theta\sim m}[\min_i\eta_i(\theta)].
$$

This is an exact probability identity. Replacing $\eta_i$ with fitted classifier probabilities introduces estimation error. Cross-fitting can prevent direct training/evaluation leakage when groups are handled correctly; it cannot by itself certify accuracy. The later reported synthetic counterexample and Brier-certificate gap are discussed in the evidence document.

## 6. Radiation normalization remains unresolved

The historical record displays a conversion proportional to an integral of $\Omega_{\rm gw}(f)$ over $d\ln f$. The standard radiation prefactor belongs to a gravitational-wave-to-photon density ratio. If $\Omega_{\rm gw}$ is the conventional critical-density fraction, the photon-density conversion must also appear.

The original convention, integration limits, spectral exterior tails, epoch, and transfer history must be specified before the calculation can support an empirical constraint. This publication preserves the issue rather than silently repairing the archived equation.

## Symbols and units

| Symbol | Meaning | Units or condition |
|---|---|---|
| $c,u,v$ | Light speed, hypothetical signal speed, relative observer speed | Speed |
| $T,t_1,t_2$ | Send, reception, return times in the stated frame | Time |
| $W,\beta,\delta$ | Signal factor, observer-speed fraction, $W-1$ | Dimensionless |
| $P,P_+,Q_{\rm occ},M_1,J$ | Defined causal-loop measures | Dimensionless; population/model dependent |
| $a,b$ | Beta-distribution shape parameters | Positive and dimensionless |
| $\Delta N_{\rm eff}$ | Extra-radiation quantity in a specified convention | Dimensionless |
| $\varepsilon_{\rm energy}$ | Archived radiation-budget fraction | Not the charged-lepton epsilon of other programs |
| PCI | Posterior-overlap index | Between zero and one; not an FTL probability |
| $\varepsilon_{\rm pooled},\sigma_\varepsilon$ | Harmonized ringdown deviation and assigned uncertainty | Same units; normally dimensionless fractional deviation |
| $x,y,z,R,D_E$ | Proposed diagnostic coordinates and functions | Dimensionless; acceptance/calibration assumptions required |
