# Explicit formulas and their admissibility conditions

8 October 2026, America/Los_Angeles. The following elementary derivations are supplied for checking the next calculation. They are distinct from the primary authors' statements in `PRIMARY_LITERATURE_REVIEW_20261008.md`, and from an executed numerical certificate. All bounds used in a certificate must be evaluated with enclosing arithmetic on the actual target, input and domain.

## Residual propagation for a linear state

Let `y′=A(t)y+b(t)` on `[t0,T]`, with an absolutely continuous approximation `v`, and define its continuous residual

`r(t)=v′(t)−A(t)v(t)−b(t)`.

For `e=y−v` and exact fundamental matrix `Φ(t,s)`, variation of constants gives

`e(t)=Φ(t,t0)e(t0)−∫[t0,t] Φ(t,s)r(s) ds`.

If `||Φ(t,s)||≤K(t,s)`, `||e(t0)||≤δ0`, and `||r(s)||≤R(s)` throughout the domain, then

`||e(t)||≤K(t,t0)δ0+∫[t0,t]K(t,s)R(s)ds`.

An exact or validated propagator bound is preferable for an oscillatory system. A general fallback uses a matrix measure/logarithmic norm `μ(A)` in the **declared norm**, giving `K(t,s)≤exp(∫[s,t]μ(A(u))du)`. Replacing it by `||A||` remains valid but can be excessively pessimistic. Node residuals, an unbounded interpolant between nodes, and double-precision evaluation do not establish `R(s)`.

For a differentiable positive Hermitian weight `H(t)`, define `||e||H=(e†H e)^(1/2)`. If

`H′+A†H+HA ≤ 2a(t)H`

in Hermitian matrix order, then

`||e(t)||H(t) ≤ exp(∫[t0,t]a)δH0 + ∫[t0,t]exp(∫[s,t]a)||r(s)||H(s)ds`.

The inequality follows by differentiating `e†H e` and applying Cauchy–Schwarz, with the zero-norm case obtained by continuity. Its use requires enclosing positivity of `H`, the matrix-order inequality, and all change-of-norm factors. Complex states can equivalently be split into real and imaginary components.

## The registered transformed mode system

For the study's transformed system `U′=W`, `W′=2ikW−g`, define interpolant residuals `rU=vU′−vW` and `rW=vW′−2ikvW+g`. Its exact free propagation over `d≥0` is

`E(d)=exp(2ikd)`, `P(d)=d exp(ikd)sinc(kd)`,

`S(d)=[[1,P(d)],[0,E(d)]]`.

For real `k>0`, `|E|=1`, `|P|≤p(d)=min(d,1/k)`. The continuous bounds therefore imply

`|eW(t)|≤δW0+∫[t0,t]RW(s)ds`,

`|eU(t)|≤δU0+p(t−t0)δW0+∫[t0,t]RU(s)ds+∫[t0,t]p(t−s)RW(s)ds`.

At `k=0`, use the continuous limit `P(d)=d` and `p(d)=d`. In the Euclidean state norm the exact operator norm is `(sqrt(4+|P|²)+|P|)/2`. Thus the homogeneous stability factor does not grow as `exp(kT)`. These are the same elementary semigroup ingredients used by the study's mathematical agent; that agent's `mathematics/LATER_TRAJECTORY_THEOREM.md` declares the observable-specific result and archive obstruction. Rounded snapshots do not determine the differentiable path between them, so these formulas cannot be applied to an unspecified archived interpolant. Source enclosure errors must enter the continuous residual rather than being silently set to zero.

## Oscillatory second-order modes

For `q″+Ω(t)²q=f`, with `Ω(t)>0`, let `x=(q,q′)` and use `H=diag(Ω²,1)`. The homogeneous oscillator generator satisfies `A†H+HA=0`; the only growth term is `H′`. One may take `a=|Ω′|/Ω`. Thus a state residual `(r1,r2)` contributes `sqrt(Ω²|r1|²+|r2|²)` under the preceding integral bound. For a scalar interpolant with its exact derivative, `r1=0` and `r2=v″+Ω²v−f`.

For constant `Ω=ω>0`, the exact propagator is an isometry in this energy norm. It follows directly that

`sqrt(ω²|e_q(t)|²+|e_q′(t)|²) ≤ sqrt(ω²|e_q(t0)|²+|e_q′(t0)|²)+∫|r2(s)|ds`.

Consequently `|e_q|≤E/ω`, `|e_q′|≤E`. The free oscillation causes no exponential factor in this norm. This is a mathematical bound conditional on the complete continuous residual; it does not imply that any sampled trajectory is enclosed. It excludes `ω=0`, and a variable frequency approaching zero needs a different uniformly positive weight. If a perturbed oscillator has a negative or sign-changing squared frequency, this positive-energy argument cannot simply be assumed.

Given enclosing state errors `|δq|≤a`, `|δp|≤b` and bounded approximations, quadratic observable differences obey, for example,

`||q|²−|v|²|≤2|v|a+a²`,

`|Re(q p*)−Re(v w*)|≤|v|b+|w|a+ab`.

The full observable still requires its momentum weights, contacts, integration and ultraviolet errors; these inequalities alone do not certify an observable gate.

## Half-space scalar Dirichlet-to-Neumann and Robin resolvents

Declare flat signature `(+,−,−,−,−)`, `y≥0`, one-sided bulk, free real scalar `Φ`, `M>0`, boundary real scalar `q`, `m>0`, `c>0`, and coupling `g`. Take

`Sbulk=½∫dy d4x[(∂tΦ)²−|∇xΦ|²−(∂yΦ)²−M²Φ²]`,

`Sboundary=∫d4x{½[(∂tq)²−|∇xq|²−m²q²]−cΦ0²/2+gqΦ0+Jq}`.

Use Fourier convention `exp(−iωt+ik·x)`. At `y=0` the **outward** unit normal is `−∂y`. The bulk boundary variation is `+Φy(0)δΦ0`, so stationarity gives

`Φy(0)=cΦ0−gq`, `(∂t²−Δx+m²)q−gΦ0=J`.

The no-incoming/outgoing retarded solution is `Φ(y)=Φ0 exp(−sR y)` with

`sR(ω,k)=sqrt(k²+M²−(ω+i0)²)`.

The square root is analytic for `Imω>0`, positive on the Euclidean/Laplace axis and continued from there; for `ω>sqrt(k²+M²)`, `sR=−ip`, `p=sqrt(ω²−k²−M²)>0`. The outward Dirichlet-to-Neumann symbol is `sR`. Eliminating the boundary value gives

`GR=1/(c+sR)`, `Φ0=g GR q`,

`DR=k²+m²−(ω+i0)²−g² GR`.

The induced self-energy enters the inverse brane propagator with a minus sign in this convention. A formal single-field in-out on-shell action does not automatically select this retarded prescription or a noise state. Initial bulk excitations add a homogeneous boundary force and must be declared.

## Exact Robin spectrum and the four-dimensional replica

For `c>0`, the half-line Robin condition is `u′p(0)=c up(0)` and has continuum modes

`up(y)=sqrt(2/π)[p cos(py)+c sin(py)]/sqrt(p²+c²)`, `p≥0`.

They have `up(0)²=(2/π)p²/(p²+c²)` and no bound state: a decaying `e^(−ay)` would require `−a=c`, impossible for positive `a,c`. The boundary Green function is

`GR=∫[0,∞]dp up(0)²/[k²+M²+p²−(ω+i0)²]`.

Set `u=M²+p²`. The mass-squared spectral density is

`ρ(u)=sqrt(u−M²)/[π(c²+u−M²)]` for `u>M²`, and zero below threshold.

For Euclidean `a=sqrt(z+M²)>0`, direct integration verifies the normalization:

`∫[M²,∞]du ρ(u)/(u+z) = (2/π)∫[0,∞]dp p²/[(p²+c²)(p²+a²)] =1/(c+a)`.

The same equality extends to the analytic retarded domain. The four-dimensional fields `Xp(x)` with mass squared `M²+p²`, continuum measure `dp`, and coupling `g up(0) q Xp` reproduce the same resolvent exactly. Matching their Gaussian incoming state also reproduces the same retained noise/influence functional. This is an explicit change of mode representation, not evidence for a fifth dimension. A finite number of closed free quadratic fields gives a rational fixed-`k` resolvent, so it cannot reproduce the branch cut exactly on a domain. Approximate finite-band matching and interacting four-dimensional continua are separate admissible rivals.

## Stability, flux and thermal noise

For fields vanishing at large `y`,

`∫dy(Φy²+M²Φ²)=∫dy(Φy+MΦ)²+MΦ0²`.

The boundary potential is therefore bounded below by one half of `m²q²+(c+M)Φ0²−2gqΦ0`. The strict positive-energy domain is `g²<m²(c+M)`. Equality admits a zero-frequency threshold mode; it is excluded when claiming strict stability. The model's dimensional units follow its action: in natural units `[Φ]=mass^(3/2)`, `[q]=mass`, `[c]=mass`, `[g]=mass^(3/2)`, so both sides of the stability inequality have mass cubed.

For real classical fields define the bulk energy current toward positive `y` by `jy=−Φt Φy`. The wave equation gives local conservation `∂t ebulk+∇x·jx+∂y jy=0`. On a finite interval `[0,Y]`, integrating gives `dEbulk/dt=jy(0)−jy(Y)` after vanishing transverse flux. Including `cΦ0²/2−gqΦ0` with the boundary energy gives `dEboundary/dt=∫J qt−jy(0)`. Thus the full signed ledger is `dEtotal/dt=∫J qt−jy(Y)`. Outgoing radiation has positive `jy(Y)`; a positive damping term alone does not establish this complete ledger or a cosmological heat source.

For a stationary Gaussian bath in a KMS state of inverse temperature `β`, define the force symmetrized covariance by `NF=½⟨{F,F}⟩` and retarded bath response by `GR=iθ(t)⟨[Φ0(t),Φ0(0)]⟩`. With these conventions,

`Im GR(ω,k)=sign(ω) sqrt(ω²−k²−M²)/[c²+ω²−k²−M²]`

above the frequency threshold, and zero below it. The force noise spectrum is

`NF(ω,k)=g² coth(βω/2) Im GR(ω,k)`.

The product is nonnegative and even in `ω`. At zero temperature `coth→signω`; noise remains nonzero. This formula is conditional on the stationary bath and its stated response convention. It must not be imposed on an arbitrary nonequilibrium or correlated preparation without deriving that state's covariance and initial terms.

At fixed `k`, `ImGR~1/ω` as positive frequency tends to infinity, so raw equal-time vacuum noise has a logarithmic ultraviolet divergence. Use noise as a distribution and define smooth spacetime-smeared observables or a declared regulator. No finite raw variance, stress renormalization or full quantum thermal history is established by the toy formula. The recent exponential-decomposition theorem's compact/exponential spectral-tail premise is unmet without such additional work.
