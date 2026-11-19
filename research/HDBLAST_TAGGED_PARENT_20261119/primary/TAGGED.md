# A tagged scalar parent and a massless Dirac pair on the four recorded walls

A resting scalar prepared with one extra quantum produces a calculable excess of massless fermion–antifermion pairs when a weak Yukawa interaction is switched on for a finite time. Across the four recorded HDblast wall trajectories, this excess is lower for the nearby positive pulse, higher for the nearby negative pulse, and modestly higher for the distant positive pulse. The same signs occur for two interaction durations, two wall meshes, and an explicit daughter momentum screen. The magnitude depends on each choice.

## The experiment

Use the archived reciprocal $y=1,b=1$ wall histories as prescribed geometry. The archive's $y$ labels its historical simulation. Here the new probe coupling $y_{\rm p}$ is independent and infinitesimal. Evolve a free minimally coupled scalar at zero comoving momentum from its physical Hamiltonian ground preparation at $s_i\simeq0.769$ using $\mu=m/H_0=100(\phi_b-0.5)$ and $s=H_0\tau$. A one-particle excitation uses that same mode function. Each calculation compares initial occupation one with initial occupation zero on the identical wall, using a fermion conformal vacuum in both cases. All other scalar modes have the same initial state in this subtraction.

Normalize $v=a\chi$ by $v v^{*\prime}-v^* v'=i$ in conformal time, and by $1/\sqrt{V}$ in a large comoving normalization box. The initial conditions are $\chi_i=(2a_i^3|\mu_i|)^{-1/2}$ and $\dot\chi_i=-i|\mu_i|\chi_i$. The mode satisfies

$$
v''+\left(a^2\mu^2-a''/a\right)v=0 .
$$

A single massless Dirac species couples through $y_{\rm p}f(s)\chi\bar\psi\psi$. The interaction vanishes outside $[s_1,s_2]$ and has $f=\sin^4[\pi(s-s_1)/(s_2-s_1)]$ inside. Wide means $[1.15,2.05]$ and narrow means $[1.25,1.95]$. These are two physical proper-time switches. Both the coupling and its first three derivatives vanish at the endpoints. The pulse and quiet members of a comparison see the same switch as a function of proper time.

Define $g(\eta)=f(s(\eta))v(\eta)$ and $G(\omega)=\int g(\eta)e^{+i\omega\eta}d\eta$. The final pair has equal and opposite comoving momenta of magnitude $\omega/2$. The massless spin trace gives a factor two after fermion wave-function normalization. The angularly integrated single-particle momentum measure thus gives $\omega^2d\omega/(8\pi^2)$ in the continuum box limit. The two excess channels obey

$$
\frac{dN_{\rm ann,stim}}{y_{\rm p}^2}
=\frac{\omega^2 |G(\pm\omega)|^2}{8\pi^2}d\omega,\qquad
\frac{dE_{\rm pair}}{H_0 y_{\rm p}^2}
=\frac{\omega}{a(s_{\rm cut})}\frac{dN}{y_{\rm p}^2}.
$$

Annihilation takes the tagged parent from occupation one to zero. Creating another quantum has squared matrix element two from occupation one and one from the reference vacuum, leaving one unit of stimulated excess. Both channels enter inclusive pair number with the same sign. Their difference gives excess net loss of parent occupation in this perturbation calculation. The energy is per initial extra particle, evaluated at $s_{\rm cut}=2.10$ after the switch is off. It counts both daughters. No physical volume multiplier or population density enters these per-parent results.

Parseval yields an independent proper-time number check,

$$
\frac{N_{\rm ann}+N_{\rm stim}}{y_{\rm p}^2}
=\frac{1}{4\pi}\int |(fv)'|^2 d\eta
=\frac{1}{4\pi}\int\frac{|a\dot f\,v+fv'|^2}{a}ds .
$$

For a resting free parent of constant mass $m$ in flat spacetime it gives the exact switched result $[m\int f^2\,ds+m^{-1}\int\dot f^2\,ds]/(8\pi)$. The code checks this analytic normalization with $m=32$.

## Results at the common cut

All numbers in the table are for the fine wall mesh. The quiet columns are inclusive pair count divided by $y_{\rm p}^2$ and pair energy divided by $y_{\rm p}^2H_0$. The pulse columns are percent changes relative to quiet, again in number / energy order.

| Proper-time window | Total comoving pair energy | Quiet number / energy | Near positive (%) | Near negative (%) | Far positive (%) |
|---|---|---:|---:|---:|---:|
| Wide | Full spectrum | 0.494140 / 10.528530 | −0.06742 / −0.06723 | +0.47033 / +0.61017 | +0.06164 / +0.10570 |
| Wide | $\omega\le60$ | 0.491716 / 10.401708 | −0.07045 / −0.07472 | +0.47967 / +0.63478 | +0.06067 / +0.10378 |
| Narrow | Full spectrum | 0.408355 / 8.252938 | −0.07120 / −0.09463 | +0.47150 / +0.70501 | +0.05721 / +0.10018 |
| Narrow | $\omega\le60$ | 0.406972 / 8.180152 | −0.07345 / −0.10076 | +0.47828 / +0.72483 | +0.05643 / +0.09847 |

On the quiet wall the wide-window count splits into 0.339193 from annihilation and 0.154948 from stimulated creation. The narrow-window values are 0.277675 and 0.130680. The instantaneous flat-width estimates $\int f^2|\mu|\,ds/(8\pi)$ would give 0.154793 and 0.115062. The exact switched probe retains parent-mode evolution, finite-time energy resolution and the stimulated contribution. Its coefficients multiply $y_{\rm p}^2$.

The screen retains each daughter's comoving momentum up to 30. On the quiet fine mesh it discards 0.491% of the wide number and 1.205% of the wide energy. For the narrow switch those fractions are 0.339% and 0.882%. The largest physical momentum admitted by the screen during any switched history is $39.2282\,H_0$, below the archived reference value $53.6107\,H_0$. The unscreened results include the entire numerical spectrum. This test of a shared screen makes no assumption about an ultraviolet completion.

These results apply to one extra rest-parent state, a specified switching function and prescribed background trajectories. Gravity and the coupling switch can transfer energy. No coupled stress response, thermalization, finite-density population yield, or expansion history is inferred from a single tagged-parent result.

## Reproduction

From this directory, use a checkout of HDblast with the 2026-11-18 source-work wall ZIP and the fine/coarse 2026-11-18 pulse-response rolling summaries. The code verifies their SHA-256 digests. NumPy and SciPy are required.

```sh
export OPENBLAS_NUM_THREADS=1
python -I -B tagged.py "$REPO" > TAGGED.json
python -I -B tagged.py "$REPO" --samples 4096 --padding 64 > /tmp/hdb-tagged-refined.json
python -I -B check_tagged.py TAGGED.json /tmp/hdb-tagged-refined.json > CHECK.json
```

The proper-time solver evolves the field and its canonical momentum. The conformal-time solver independently evolves the rescaled field and its momentum, with proper time as a dynamical variable. Both use DOP853 and cubic Hermite wall interpolation with the stored slopes. Across all histories their relative endpoint mode discrepancy is below 1.84 × 10⁻⁹. The maximum absolute Wronskian error is 8.98 × 10⁻¹². The proper-time integral agrees with the independently sampled Fourier count to 8.53 × 10⁻¹⁰. Doubling both sampling and zero padding changes the inclusive number by at most 1.09 × 10⁻¹¹ and the final energy by 3.40 × 10⁻⁹. The largest change in a quoted pulse contrast between the two archived wall meshes is 0.00271 percentage point.
