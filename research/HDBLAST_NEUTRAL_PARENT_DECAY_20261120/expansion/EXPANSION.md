# How much expansion can follow the heating?

The coupled solution reaches more than 90% particle internal energy very close to maximum size. We can ask a sharper question. From that point, could a small permanent vacuum component keep a closed wall expanding if energy exchange with the bulk stopped? For this comparison we keep the total wall energy at the matching time fixed. A second calculation determines the radius needed at the same heating density for a radiation-filled wall to double in size. [Calculation](fate.py) · [Values](FATE.json)

## The late closed-wall comparison

Let $a_*,e_*,m_b$ be the radius, total energy density above the background tension, and bulk mass at a time when the wall is expanding. Energy density is in units of that tension. We choose the first archived post-exit expanding sample with at least 90% internal particle energy. For the scalar-threshold simulation it is at wall time 11.28, with $a_*=9.6884$, $e_*=0.00537510$, $m_b=0.44516$ and $H_*=0.01327$. The Friedmann equation reproduces the sampled $H_*^2$ within $3.0\times 10^{-12}$. The first crossing in the finer integration is at time 11.26695. [Coupled solution](../coupled/DECAY.json)

After this slice, assume a constant bulk mass, no further exchange, a permanent positive vacuum density $\lambda$, and a separately conserved material component. For constant pressure-to-density ratio $w$, put $n=3(1+w)$, $y=a_*/a$, and $e(y)=\lambda+(e_*-\lambda)y^n$. The exact Friedmann function in the study units is

$$
F(y)=H^2=e(y)[2+e(y)]+(m_b/a_*^4)y^4-y^2/a_*^2.
$$

All choices of $\lambda$ at fixed $a_*,e_*,m_b$ give the same initial expansion rate. The vacuum allocation is hypothetical. The computed rolling field has already reached zero potential on this slice. Ordinary pressureless matter gives $n=3$, and radiation gives $n=4$. The square in the Friedmann equation is retained.

For $n\ge3$ and nonnegative $m_b$, $F'(y)/y$ increases with $y$. Explicitly, it equals $2n(e_*-\lambda)[1+e(y)]y^{n-2}+4m_by^2/a_*^4-2/a_*^2$. Thus the minimum over the future expansion lies at a unique interior stationary point or an endpoint. At every $0<y<1$, increasing $\lambda$ increases $F$ at fixed initial total density. Solving $F=F'=0$ gives the boundary between a finite-time turnaround and continued expansion. At the boundary the solution approaches the double root asymptotically. The code independently recovers the threshold from the maximum of the vacuum fraction required at each radius.

| Matter in the continued solution | Critical vacuum fraction at the matching time | Additional radius growth before turnaround with 10% vacuum | With zero vacuum |
| --- | ---: | ---: | ---: |
| Pressureless matter | 25.61% | 2.36% | 1.63% |
| Radiation | 41.34% | 1.03% | 0.82% |

For any material with nonnegative pressure and no subsequent energy input, conservation gives $d(a^3e_{\rm mat})/da=-3a^2p_{\rm mat}\le0$. Its density is therefore bounded by the pressureless comparison at every larger radius, even when its equation of state varies. In this class of continuations the 25.61% threshold is a necessary lower bound on a permanent vacuum share to expand indefinitely from this slice. A measured internal share above 90% leaves less than 10% of the same energy budget for such a component. The coupled system, which still evolves its own fields and bulk, grows another 0.76% from this sample before contracting. The previous clock-activated trajectory gives thresholds of 24.40% and 39.94% at its own first sampled 90% crossing. The bound depends on the slice and assumptions. A later supply of energy or an evolving component with negative pressure requires its own dynamics.

## Target radius for a longer radiation era

A second comparison keeps the heating density and bulk mass fixed and varies the radius at which heating finishes. The starting expansion rate now follows from the new Friedmann constraint. Take zero vacuum and a desired radius growth factor $b$. At the endpoint the density is $e_f=e_*/b^n$. Write $S=e_f(2+e_f)$. The large-radius turning point satisfies $0=S-u+m_bu^2$ with $u=1/(b a_*)^2$. For the measured positive mass and the outer branch $a^2>2m_b$, the minimum starting radius to reach $b$ is

$$
a_{*,\min}=\frac{1}{b}\sqrt{\frac{1+\sqrt{1-4m_b S}}{2S}}.
$$

For $\lambda=0$ and $n\ge3$, $F(y)/y^2$ increases with $y$. Reaching the endpoint suffices to check the entire intervening expansion. At the matching density a radiation-filled wall needs an initial radius of **19.29** to double, **26.22** for one e-fold, and **193.72** for three e-folds. The measured starting radius is **9.69**. Equality puts the first turnaround exactly at the target. These are targets for a new coupled source and expansion calculation. Achieving the required radius while producing the same heating density remains an open problem.
