# Independent analytic infrared coherence audit

This note performs no source, time, momentum, or mode numerical evaluation.
Let s=epsilon f be a real smooth compact source, F0=integral f(t)dt,
F1=integral t f(t)dt, and Fhat(2k)=integral f(t)exp(-2ikt)dt. After the pulse,
write the canonical mode as [alpha exp(-ik eta)+beta exp(ik eta)]/sqrt(2k),
with alpha=1+epsilon alpha1+epsilon^2 alpha2+... and
beta=epsilon beta1+epsilon^2 beta2+.... Direct first-order forcing gives

    alpha1=-iF0/(2k), beta1=iFhat(2k)/(2k),
    u1=alpha1+beta1 exp(2ik eta) -> F1-eta F0 as k->0.

For a nonzero nonnegative compact pulse, F0>0. The isolated occupation
contribution to physical variance at second order is

    Qocc^(2)=epsilon^2/(2pi^2 a^2) integral k |beta1|^2 dk,

whose infrared leading term is
epsilon^2 F0^2/(8pi^2 a^2) integral dk/k. It is logarithmically divergent.
The canonical excitation energy is instead

    Ecan^(2)=epsilon^2/(2pi^2) integral k^3 |beta1|^2 dk,

with infrared behavior epsilon^2 F0^2/(8pi^2) integral k dk, which is finite.
This canonical quantity is not the complete physical minimal stress.

The full variance must retain coherence. Using the Wronskian normalization,
its change is (2pi^2 a^2)^(-1) integral k[|beta|^2+
Re(alpha beta* exp(-2ik eta))]dk. At second order,

    |beta1|^2 ~ F0^2/(4k^2),
    alpha1 beta1* ~ -F0^2/(4k^2) + O(1/k).

These leading logarithmic terms cancel. The additional beta2* term must not
be omitted from a full second-order response. Because the full u1 is bounded
on each fixed finite interval, its second Born forcing gives beta2=O(1/k),
with a purely imaginary leading coefficient for real f. That coefficient
has bounded real contribution at fixed finite observation time.
Combine occupation and coherence with the same infrared regulator before
removing it; separately divergent integrated pieces cannot be added as numbers.

More directly, u2 obeys u2''-2ik u2'=-f u1 with zero initial perturbation.
The Volterra kernel sin(k Delta)/k has the regular limit Delta. Thus u1 and
u2 remain bounded as k->0 on the fixed compact time interval, and the full
second-order variance integrand is

    epsilon^2 k[|u1|^2+2Re(u2)]/(4pi^2 a^2),

which is infrared integrable. The positive-reference subtraction is also
regular at k=0. A late-time limit, a varied initial state, or an operation
discarding coherence requires a separate analysis; one cannot exchange those
limits with this fixed-time cancellation silently.

The first-order alpha/beta truncation is not a normalized finite-amplitude
state. The second-order Wronskian condition is

    2Re(alpha2)+|alpha1|^2-|beta1|^2=0.

For the registered odd signed pulse uB(u), F0=0, so its isolated quadratic
occupation variance is infrared finite at this order. That cancellation does
not generally force the full higher-order zero-frequency transfer to vanish.
No occupation-only variance or first-order physical stress result establishes
heating, particle yield, or thermalization.
