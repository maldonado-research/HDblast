# How long can the finite photon share stay above one half after injection ends?

In the finite model with maximum photon shell $L_{\max}=48$ and total change $\delta=3$ in $\ln I$, radiation already has more than half the energy when the gauge coupling stops changing. Both fitted branches keep that majority through the entire expanding interval tested below. Raising the cutoff to 128 at the same $\delta$ puts the peak below one half on both branches. The predicted majority therefore depends on how many photon shells the model retains.

We begin at the first outgoing crossing of 80% of the directed field displacement to the minimum. This is the upper edge of the switch. All photon occupations are then constant while the field stays beyond the switch. The isotropic, vacuum-subtracted energy of the finite transverse shell ensemble obeys $r a^4=C$, $p_\gamma=r/3$, where $C$ comes from the coupled oscillator solve. Expansion continues according to

$$
q'=v,\quad v'=-3hv-3u_{,q},\quad a'=ha,\quad
h'=a^{-2}-v^2/2-2r,\quad r=C/a^4.
$$

We stop at the first radius maximum or the first return to the 80% surface in the opposite field direction. Re-entering the switch can resume energy exchange. No continuation across that surface is included here. The program also records a calibrated field range stop or a time limit if either is reached first. Neither limit was reached in the reported trajectories.

| Fitted branch | $L_{\max}$ | $\delta$ | Share at switch end | Highest share after switch | Share at stop | Stop | $\Delta\ln a$ with share ≥ 1/2 |
|---|---:|---:|---:|---:|---:|---|---:|
| $+0.5$ | 48 | 3 | 0.63614 | 0.64234 | 0.64111 | Radius maximum | 0.02170 |
| $-0.5$ | 48 | 3 | 0.51883 | 0.52609 | 0.51670 | Switch re-entry | 0.06168 |
| $+0.5$ | 128 | 3 | 0.36248 | 0.36922 | 0.35792 | Radius maximum | 0 |
| $-0.5$ | 128 | 3 | 0.24996 | 0.25607 | 0.24838 | Switch re-entry | 0 |
| $+0.5$ | 32 | 1.5 | 0.51254 | 0.52000 | 0.51506 | Radius maximum | 0.02749 |
| $-0.5$ | 32 | 1.5 | 0.34801 | 0.35433 | 0.34648 | Switch re-entry | 0 |
| $+0.5$ | 48 | 0.1 | 0.00575 | 0.00590 | 0.00511 | Switch re-entry | 0 |
| $-0.5$ | 48 | 0.1 | 0.00694 | 0.00711 | 0.00689 | Switch re-entry | 0 |

The durations in the last column begin at the switch end for all three above-half rows. They only count the source-free part of the expansion. At $L_{\max}=48$, $\delta=3$, the positive branch runs from $x=0.68477$ to $0.98086$, and the negative branch from $x=0.72193$ to $0.97584$.

A follow-up test deliberately searched for a threshold crossing *after* the switch. We chose $L_{\max}=32$, $\delta=1.535$ after inspecting nearby $\delta$ values on the positive branch and then evaluated both branches. The positive-branch fraction begins at 0.49648, crosses one half at $x=0.58286$, reaches 0.50391 at $x=0.67469$, then falls below one half at $x=0.86965$. Radiation holds a majority for $0.28679$ model-time units and $\Delta\ln a=0.01747$. At first turnaround the fraction is 0.49853. The negative branch at the same settings peaks at 0.33525. This post hoc example tests how the fraction behaves near equality in a specified finite model.

During this phase the varying fraction is a result of expansion alone. With $e_\phi=u+v^2/6$, $w_\phi=(v^2/6-u)/e_\phi$ and $F=r/(r+e_\phi)$, component conservation gives

$$
F'=3hF(1-F)(w_\phi-1/3).
$$

For an expanding interval with both densities positive, stationary fractions occur at $v^2/6=2u$. We locate these events, evaluate the extrema and endpoints, and find any interior equality $r=e_\phi$ between successive extrema. The duration counted between successive half crossings is $\Delta N=\ln(a_{\rm end}/a_{\rm start})$ whenever the fraction is above one half within that interval. The photon density itself decreases throughout expansion.

Two independent evolutions use $(\ln a,h,C)$ and $(a,a',r)$ with different integrators, step sizes and tolerances. Across all ten cases their largest difference in selected times, states, peak fractions and expansion durations is below $2.0\times10^{-10}$. All first minima are also present in the complete photon-oscillator integration. Matching those independently evolved events gives a maximum deviation below $8.8\times10^{-12}$ in the stated model units. A finite source-free radiation majority establishes neither a thermal photon distribution nor sustained dominance past the end of the chosen interval.

[Source-free solver](coast.py) · [Full oscillator data](FEEDBACK.json)
