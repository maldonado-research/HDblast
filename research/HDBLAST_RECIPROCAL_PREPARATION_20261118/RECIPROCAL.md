# Two finite exposures distinguish an aligned displacement release

Can the shell’s quiet interval and subsequent pulse be produced by displacing any number of free four-dimensional scalar fields through the same channel that later measures them? For the preparation below the answer is no. In two finite exposures a reciprocal release obeys $|L|\le1.053|E|$. The registered shell has $E=0$ and $4.11449<L<4.12482$.

## Preparation and clock

Use the same quadrupole $Y$, de Sitter background and elapsed clock $x=q/\ell\in[0,X]$, with $\ell=10^{-6}$ and $X=25/16$. The physical target in the previous edition’s normalized amplitude units is $A=(\ell/G)S(-\ell x)$. For a real canonical free field of scaled mass squared $P_j=(m_j\rho_b\ell)^2\ge0$, write its quadrupole as $a_j(x)Y$. Its equation in this clock is

$$
 a_j''+3H a_j'+K_j a_j=0,\qquad H=\ell\tanh(\ell x),\qquad
 K_j=P_j+8\ell^2\operatorname{sech}^2(\ell x).
$$

Hold the clock and the readout $F=\sum_j c_j a_j$ fixed. Specify the release at $s\in[0,31/32]$. Each field starts from rest, $a_j'(s)=0$. The readout weights $w_j=c_j a_j(s)$ all have the same sign and $\sum_j|w_j|<\infty$. This includes a reciprocal preparation specified independently of the signal. Minimize the instantaneous quadratic energy at $s$ under a linear bias $-J\sum_j c_j a_j$, then remove the bias and evolve freely. With a common positive angular and volume normalization absorbed in $J$, this gives $a_j(s)=Jc_j/K_j(s)$ and $w_j=Jc_j^2/K_j(s)$. The choice of masses, couplings and overall nonzero gain is arbitrary. An impulsive force generically prepares velocities and is a different experiment.

For each mode of unit displacement and zero velocity let $a(s)=1$, $a'(s)=0$. Direct differentiation gives

$$
 \frac{d}{dx}\left(a^2+\frac{a'^2}{K}\right)
 =-2H\left(2+\frac{P}{K}\right)\frac{a'^2}{K}\le0 .
$$

Thus $|a(x)|\le1$ throughout the post-release interval. At a perfectly timed release sample this already gives $|F(x)|\le|F(s)|$. The shell is quiet through $x=1$ and its uniform mean on $[5/4,11/8]$ exceeds $4.01324$. Comparing that mean with the release value requires a maximum absolute mismatch greater than $2.00662$. [Prior waveform, proper-clock units and window means](https://doi.org/10.5281/zenodo.24351909) · [Propagation and quiet interval, §4](https://github.com/maldonado-research/HDblast/blob/20f431f2f6fd8ae9d040d0e8b87bddb41a2f3fa2/research/HDBLAST_PAST_SHELL_CONTACT_20261104/primary/PROOF.md)

## Two positive shutters

The early shutter begins at the release. Its raw shape for $s\le x\le s+h$ is $p(x)=3[1-(x-s)/h]^2/h$, with $h=1/32$. The late shutter $v$ is a symmetric triangle on $[5/4,11/8]$, of width $b=1/8$ and peak $2/b$. Both shapes have unit integral and are zero elsewhere. Put $C(x)=\cosh^{3/2}(\ell x)$. The two reported **physical** exposures of any signal $F$ are

$$
 E(F)=\frac{\int p C F}{\int p C},\qquad
 L(F)=\frac{\int v C F}{\int v C}.
$$

The factor $C$ prescribes a known positive adjustment to each shutter. It varies by less than $2\times10^{-11}$ from unity on the horizon. Both physical shutters are normalized and have finite duration. They can also be implemented by reweighting a resolved trace. The clock registration associates the early shutter with the actual release. The bounds cover every stated $s$ separately.

For a unit mode write $y=C(x)a(x)/C(s)$, $\omega^2=P+13\ell^2/2$ and $d(x)=(35/4)\ell^2\tanh^2(\ell x)$. The conjugated equation, data and bounds are

$$
 y''+\omega^2 y=d(x)y,\qquad y(s)=1,\qquad y'(s)=\beta_s=\tfrac32\ell\tanh(\ell s),
 \quad 0\le d(x)\le d_*=(35/4)\ell^4X^2 .
$$

From $|a|\le1$ and $1-10^{-11}<1/C(x)\le1$, we have $|y|\le\Gamma=(1-10^{-11})^{-1}$. For $e(x)=y(x)-\cos(\omega(x-s))$, the retarded sine formula and $|\sin(\omega t)|/\omega\le t$ give $|e|\le\eta=\beta_*h+\Gamma d_*h^2/2$ throughout the early shutter, where $\beta_*=3\ell^2(1-h)/2$. A second bound retains control for arbitrarily large mass. Integrate $e''+\omega^2 e=dy$ against $p$ twice by parts on $[s,s+h]$. The sole boundary contribution to $\omega^2\int pe$ is $p(s)e'(s)$, and $\int_s^{s+h}|p''|=6/h^2$. With $z=\omega h$ this gives

$$
 z^2\left|\int pe\right|\le \xi=3h\beta_*+6\eta+h^2\Gamma d_* .
$$

The exact flat cosine mean in that shutter is $f(z)=6(1-\sin z/z)/z^2$. For $0<z\le1$, integrating $\cos t\ge1-t^2/2$ gives $f(z)\ge1-z^2/20\ge19/20$. For $z\ge1$, $\sin z/z\le101/120$. Indeed sinc decreases on $[1,2]$, $\sin1\le1-1/6+1/120$, and $\sin z/z\le1/2$ for $z\ge2$. Hence $z^2f(z)\ge19/20$ for $z\ge1$. The early conjugated mean is strictly positive for every mass.

The late shutter gives $|\int vy|\le\Gamma$. Its continuous piecewise-linear shape has second distributional derivative of total variation $16/b^2$, including the two endpoint corners. Testing the equation against it gives $z^2|\int vy|\le\Gamma(16h^2/b^2+h^2d_*)$. Dividing the respective low- and high-frequency bounds yields

$$
 \frac{|\int vy|}{\int py}\le
 \max\left\{\frac{\Gamma}{19/20-\eta},
 \frac{\Gamma(16 h^2/b^2+h^2d_*)}{19/20-\xi}\right\}<1.053.
$$

In the ratio of normalized physical exposures the additional factor is $\int pC/\int vC\le1$, since $C$ increases and the early shutter precedes the late shutter. Each unit release therefore has $E(a)>0$ and $|L(a)|\le1.053 E(a)$. The same-sign weights give **$|L(F)|\le1.053|E(F)|$**. Absolute summability and $|a|\le1$ justify the sum and the integrations for countably many species. The estimate permits unbounded masses. An ordinary unweighted early box has a cosine transform that changes sign; this argument uses the specified positive-response shutter.

## Shell comparison and precision

Finite propagation makes $A$ identically zero on the entire early shutter for each admitted release time. The late numerator for $A$ is the triangular mean of the prior normalized past, dilution-removed waveform $u=C A$. The exact replay widens each future master interval by $3\times10^{-6}$ in $u$ units, integrates the triangle exactly over all 512 intersecting time cells and divides by the enclosure of $\int vC$. This gives

$$
 E(A)=0,\qquad 4.11449<L(A)<4.12482.
$$

If the model is required to match the two target means within errors $e_E,e_L$, its inequality requires $e_L+1.053e_E\ge L(A)$. Equal allowable discrepancies must exceed $2.00413$ in the fixed $A$ units. Expressed relative to the actual late magnitude, the necessary maximum discrepancy is at least $1/(1+1.053)>48.70\%$. All target bounds rescale with a common gain, and the relative separation is invariant. These are deterministic tolerances on the model-versus-target means. If actual data have their own errors and the template has separate errors, combine them in the corresponding exposure budgets.

For signed displacement weights with zero release velocities define $M$ as the smaller of the two sums of $|w_j|E(a_j^{\rm unit})$ over positive and negative weights. Then $|L(F)|\le1.053(|E(F)|+2M)$. Exact agreement requires $M>1.95369$ in the normalized early-exposure units. Even with deviations of $1/2000$ allowed in each exposure, $M>1.95321$. Preparing velocities or continuing to drive fields after release changes the model. The general four-dimensional realization previously constructed uses its own initial phase-space data. The inequality above neither determines those data nor provides independent physical evidence that a cosmological observer could prepare or measure this reciprocal channel. [Full-trace comparison construction](https://github.com/maldonado-research/HDblast/blob/20f431f2f6fd8ae9d040d0e8b87bddb41a2f3fa2/research/HDBLAST_DIMENSIONAL_COMPARISON_20261118/TOWER.md)
