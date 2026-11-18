# An in/out experiment with stationary past and future

A real crossing of the massless point is unnecessary for excitation of a time-dependent oscillator. The earlier D-Blast scalar-pulse calculation illustrates this with an exactly solvable mass profile. Here we test a different question. If an external apparatus applies a mass modulation measured from the B1 shell, then returns the mass to its initial value, does a strictly positive-frequency incoming mode have an outgoing negative-frequency part? In this explicitly specified stationary experiment the answer is yes for **every** comoving wave number in the finite band $5\le k\le18$. The following bound applies to the declared mass protocol, including its rounding and return switch. The expanding shell and its later evolution are tested separately. [Earlier project calculation](../literature/PORTFOLIO.md) · [Pulse coefficients](PULSES.json) · [Mode and detector studies](../README.md)

## Definition of the comparison

Use dimensionless time $s$ and frequencies in units of $H_0$. Set $m_0=100(0.5-\phi_b(0))$ and $\omega^2=k^2+m_0^2$. On each measured history integrate the shell proper clock from the lapse and interpolate $\phi_b$ linearly onto $s_j=jh$ with $h=1/400$ and $j=42,\ldots,58$. The coefficients $c_j=10^4(\phi_p-\phi_q)(\phi_p+\phi_q-1)$ are the differences of the squared masses between a pulse and its quiet control at equal proper time, rounded to 12 significant figures in the published table. Here $p,q$ label pulse and quiet. The frozen background and a smooth return are choices of the comparison.

Use the centered cardinal cubic B-spline $B_3$. It is the convolution of four unit-width centered top hats of unit area. Prescribe

$$V(s)=\sum_j c_j B_3((s-jh)/h),\qquad y''+(\omega^2+V)y=0.$$

The modulation and its first two derivatives vanish outside $[0.10,0.15]$. The spline uses the measured values as coefficients and smooths the profile over nearby samples. Since the splines are nonnegative and their full integer translates sum to one, $m_0^2+V\ge m_0^2+\min(0,c_j)>0$. Both asymptotic intervals have the same definite frequency. This permits an unambiguous in/out particle definition for this externally driven oscillator. [Authenticated extraction](make_pulse.py)

## A lower bound without a weak-drive assumption

Normalize the incoming solution to $e^{-i\omega s}$ before the pulse and write the late solution as $\alpha e^{-i\omega s}+\beta e^{i\omega s}$. The retarded sine Green function gives a first-order term and an exact remainder bound:

$$
\beta_1=\frac{i h}{2\omega}\operatorname{sinc}^4(\omega h)
               \sum_j c_j e^{-2i\omega jh},\qquad
|\beta-\beta_1|\le\tfrac12(e^L-1-L),\quad
L=\frac{h\sum_j |c_j|}{\omega}.
$$

To see the bound, the retarded equation has integral kernel $-\sin(\omega(s-t))/\omega$. Grönwall's inequality bounds the deviation from the incoming solution up to time $t$ by $e^{L(t)}-1$, with $L(t)=\int^t|V|/\omega$. Inserting this deviation in the exact expression for the outgoing coefficient and integrating gives the stated exponential remainder. The leading Fourier integral follows from the four top-hat convolutions.

The executable bound uses exact rational arithmetic on the published decimal coefficients. For the whole band, $50<\omega<54$. After rotating the first-order sum about $s=0.12$ and aligning its real part with the sign of $\sum c_j$, the largest angular distance is below 3 radians. For nonnegative aligned contributions, monotonicity of cosine on this range and its alternating Taylor lower bound through degree six give a lower projection at the largest allowed angle. Negative contributions use $\cos x\le1$. Also $\operatorname{sinc}(\omega h)\ge1-(54h)^2/6>0$. Finally $\tfrac12(e^L-1-L)\le L^2/[4(1-L/3)]$ for $L<3$. Subtracting this error from the strictly positive projection certifies the result. Rounded lower endpoints in the report are rounded down, upper endpoints up. [Rational certificate](bounds.py) · [All eight protocol bounds](BOUNDS.json)

## Results and energy source

For the finest positive pulse the exact occupation $|\beta|^2$ of each mode lies between **0.000002147118 and 0.000007681557** throughout the band. Direct evolution gives $2.82753\times10^{-6}$ at $k=10$. The opposite pulse gives $2.82041\times10^{-6}$ at that wave number. At half the positive initial scalar amplitude the result is $7.06441\times10^{-7}$. These are per-mode occupation numbers.

For one real scalar, integrating $k^2|\beta|^2/(2\pi^2)$ and its energy-weighted counterpart over the stated band gives $n/H_0^3=0.00025343$ and $\rho/H_0^4=0.0131354$ for the finest positive protocol. Numerical propagation uses exact oscillator transfers for midpoint-held frequencies, with 400 and 800 time steps. The largest sampled occupation difference between those resolutions is $2.95\times10^{-11}$ for this case. A square mass step checks the transfer against its analytic solution. The change in mode energy agrees with the work supplied at all changes of the prescribed squared mass, $\sum\Delta V |y|^2/(4\omega)$, including the final return. The apparatus supplies that work. These finite-band numbers cannot be booked as energy delivered by an HDblast shell. [Numerical check and source-work test](scatter.py) · [Results](SCATTER.json)

An independent adaptive Runge–Kutta solver integrates the smooth coefficient directly. At $k=5,10,18$ it agrees with the 800-step transfer occupation to within $1.1\times10^{-11}$ for the finest positive, negative and half-positive protocols. [Independent solver](independent.py) · [Computed occupations](scatter-independent.json)
