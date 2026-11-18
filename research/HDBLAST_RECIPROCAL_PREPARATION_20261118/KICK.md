# A collocated velocity impulse gives a detectable early exposure

Consider free, canonical four-dimensional scalar fields in the registered quadrupole sector on the same de Sitter background and clock. An instantaneous source at elapsed time zero changes their velocities in proportion to one fixed real vector of couplings. Initial displacements vanish. A detector measures the same linear combination of the fields with one arbitrary overall gain. There is no later forcing. This is a specific source and receiver model. An unrestricted preparation allows different coupling vectors and phases.

Use $x=q/\ell$, $\ell=10^{-6}$ and $Q=25/16$. For any physical quadrupole coefficient $f(x)$ define two unit-normalized finite exposures,

$$
T[f]=2\int_0^1(1-x)f(x)\,dx,\qquad
L[f]=8\int_{5/4}^{11/8}f(x)\,dx.
$$

**Every collocated velocity impulse of this kind satisfies $|L[f]|\leq9|T[f]|$.** This holds for every nonnegative squared mass and an arbitrary number of species. The registered linear target $A(x)=(\ell/G)S(-\ell x)$ has $T[A]=0$ by finite propagation and $L[A]>4.013$. Thus every such comparator has a maximum absolute error over the two exposures of at least $L[A]/10>0.4013$ in the registered normalization. That is at least 10% of the target late mean for any common amplitude calibration. Only averages enter the test. [Quiet interval, §4](https://github.com/maldonado-research/HDblast/blob/20f431f2f6fd8ae9d040d0e8b87bddb41a2f3fa2/research/HDBLAST_PAST_SHELL_CONTACT_20261104/primary/PROOF.md) · [Physical readout and clock](https://github.com/maldonado-research/HDblast/blob/20f431f2f6fd8ae9d040d0e8b87bddb41a2f3fa2/research/HDBLAST_DIMENSIONAL_COMPARISON_20261118/SINGLE_FIELD_TEST.md) · [Exact late mean](https://github.com/maldonado-research/HDblast/blob/20f431f2f6fd8ae9d040d0e8b87bddb41a2f3fa2/research/HDBLAST_DIMENSIONAL_COMPARISON_20261118/SHAPE.json)

## One unit impulse on the physical clock

Let $a_P(0)=0$, $a'_P(0)=1$, $P=(m\rho_b\ell)^2\geq0$. A prime denotes differentiation in $x$. With $v=\cosh^{3/2}(\ell x)a_P$, the free equation is

$$
v''+[\omega^2-b(x)]v=0,\quad
\omega^2=P+13\ell^2/2,\quad
b(x)=\tfrac{35}{4}\ell^2\tanh^2(\ell x).
$$

First take $s(x)=\sin(\omega x)/\omega$, the unit-velocity impulse when $b=0$. Write $B=\min(1,\omega^{-2})$. The exact early exposure is

$$
T[s]=\frac{2}{\omega^2}\left(1-\frac{\sin\omega}{\omega}\right)\geq\frac{19}{60}B.
$$

For $0<\omega\leq1$, the alternating Taylor bound on $\sin\omega/\omega$ gives $1-\sin\omega/\omega\geq\omega^2(1/6-1/120)$. For $1\leq\omega\leq2$, sinc decreases because $\sin\omega-\omega\cos\omega$ starts at zero and has positive derivative $\omega\sin\omega$ on $(0,2]$. Hence $\sin\omega/\omega\leq\sin1\leq101/120$. For $\omega\geq2$, use $\sin\omega/\omega\leq1/\omega\leq1/2$.

For the late mean, $|L[s]|\leq\min(11/8,1/\omega,16/\omega^2)$. The three terms follow from $|s(x)|\leq x$, $|\sin|\leq1$ and the cosine difference when integrating the sine over width $1/8$. On the bands cut at $1,2,4,10$, bounds for $|L[s]|/T[s]$ are respectively $165/38$, $120/19$, $4$, $20/3$, and $80/9$. In the first band divide $11/8$ by the early floor. In the second use $1/\omega$ and the sinc bound at one. In the next two bands use $1/\omega$ and $1-\sin\omega/\omega\geq1-1/\omega_{\min}$ with lower edges two and four. In the last band use the cosine difference and $1-1/10$. These estimates give $|L[s]|\leq(80/9)T[s]$ for the entire mass axis.

## Control of curvature and dilution at arbitrary mass

Put $\eta=(35/4)\ell^4Q^2$ and $M=\min(Q,1/\omega)$. Duhamel's identity gives $v=s+\int_0^x s(x-y)b(y)v(y)\,dy$ with $|b|\leq\eta$. Hence $\|v\|_\infty\leq M/(1-\eta QM)$ and

$$
\|v-s\|_\infty\leq\frac{\eta Q^3}{1-\eta Q^2}B.
$$

Here $Q\geq1$ and $M^2\leq Q^2B$. To return to the physical field, let $r=\operatorname{sech}^{3/2}(\ell x)$. On $[0,Q]$ one has $0<r\leq1$, $|1-r|\leq d=3\ell^2Q^2/4$ and $|r'|\leq e=3\ell^2Q/2$, by $\log\cosh t\leq t^2/2$ and $|\tanh t|\leq |t|$. The physical error decomposes as $a_P-s=r(v-s)+(r-1)s$.

Let $g$ be either normalized exposure weight, extended by zero. Its total variation is at most 16, including the endpoint jumps, and its integral is one. Extend $r$ constantly beyond the considered interval when calculating variation. For $\omega\leq1$, $|\int g(r-1)s|\leq Qd$. For $\omega\geq1$, one integration by parts against $\sin(\omega x)$ gives $|\int g(r-1)s|\leq(16d+e)/\omega^2$. This retains the two inverse powers needed for arbitrarily large masses. Thus both exposure errors are at most $\alpha B$ with

$$
\alpha=\frac{\eta Q^3}{1-\eta Q^2}+\max(Qd,16d+e)<10^{-9}.
$$

Consequently $T[a_P]>0$ and $|L[a_P]|\leq (80/9)T[s]+\alpha B\leq9T[a_P]$, because $(9-80/9)(19/60)>10\alpha$. The replay checks all constants and the signed target window using exact rational numbers. [Calculator](kick.py) · [Exact output](KICK.json)

## Source and receiver

For couplings $d_j$, collocation means initial velocities $k d_j$ and measured field $h\sum_j d_j a_j$ with constants $k,h$. Thus the coefficient $c_j=kh d_j^2$ multiplying each unit-impulse trace has one common sign. Positivity of $T[a_{P_j}]$ and the triangle inequality give the asserted exposure bound. Finite detector norm and finite kick energy, as when $d\in\ell^2$, suffice for convergence. More generally, for a fixed-sign spectral measure the inequality passes to any limit of finite modal exposures with finite early sum. Matching the two target means to accuracy $\epsilon$ would require $L[A]\leq|L[A]-L[f]|+9|T[f]|\leq10\epsilon$.

If a source and receiver couple with independently signed products $c_j$, define $N_\pm=\sum_{\pm c_j>0}|c_j|T[a_{P_j}]$, assuming both sums are finite. Then $T[f]=N_+-N_-$ and $|L[f]|\leq9(N_++N_-)$. Agreement to absolute exposure error $\epsilon$ requires

$$
\min(N_+,N_-)\geq \frac{L[A]-10\epsilon}{18}.
$$

At $\epsilon=1/2000$ in registered $A$ units, each sign must supply more than 0.2226 units of early integrated response. This quantifies the cancellation a four-dimensional kick would need to conceal. The exact field-tower replica uses correlated initial displacements and velocities taken from a radial packet and satisfies different preparation assumptions. Launch time, intrinsic source duration, source and receiver coupling vectors, and relative gain between exposures all belong to the physical specification. An initially silent source that acts later can evade this test. The comparison uses the existing calculated linear waveform and has no observed cosmological source input.

## A one-signed pulse with unknown launch time and duration

The exclusion survives a common source acting over an interval. Start the fields at rest and let one fixed source and receiver vector couple to a nonnegative temporal drive supported anywhere in $0\leq x\leq1/4$. A common sign reversal is allowed. The drive may be an integrable function or a finite positive measure, including discrete impulses. Assume retarded linear evolution with the source entering each field's acceleration in its stated coupling. For this source class,

$$
|L[f]|\leq12|T[f]|,\qquad
\max(|L[f]-L[A]|,|T[f]|)\geq L[A]/13>0.3087.
$$

To prove the bound, place a unit velocity jump at any launch time $s\in[0,1/4]$ and put $h=1-s$, so $3/4\leq h\leq1$. Its comparison in the flat temporal equation is zero before $s$ and is $s_\omega(x-s)=\sin(\omega(x-s))/\omega$ afterward. The early triangle now gives

$$
T_s=\frac{2h}{\omega^2}\left(1-\frac{\sin(\omega h)}{\omega h}\right)
    \geq h^3(19/60)B.
$$

Indeed $T_s$ is $h^3$ times the zero-launch triangle evaluated at $p=\omega h$, and $\min(1,p^{-2})\geq B$. The late mean obeys $|L_s|\leq\min(11/8,h/p,16h^2/p^2)$. Split the $p$ axis at the same band edges $1,2,4,10$. The five earlier ratio bounds are multiplied respectively by $h^{-3},h^{-2},h^{-2},h^{-2},h^{-1}$. At the smallest allowed $h=3/4$, their maximum is $320/27<12$.

The curved response for a unit velocity jump is $r(x)v_s(x)/r(s)$ for $x\geq s$, where $v_s$ solves the conjugated equation with data $(0,1)$ at $s$. The earlier integral bound on $v_s-s_\omega(x-s)$ holds on the shorter interval. As $r(s)>1/2$ and $r$ is decreasing for $x\geq0$, the bounds on $|r(x)/r(s)-1|$ and its derivative on $[s,Q]$ are $2d$ and $2e$. Cutting the two weights to $x\geq s$ preserves their integral bound of one and variation bound of 16. The two physical exposure errors are therefore bounded by $2\alpha B$. The inequality $(12-320/27)(3/4)^3(19/60)>26\alpha$ proves $T[a_{P,s}]>0$ and $|L[a_{P,s}]|\leq12T[a_{P,s}]$. Integrating over a one-signed temporal measure and summing squared couplings preserves the result. This proof addresses bounded uncertainty in the initial launch interval. Later emission, temporal sign changes, and a separate transfer function in the detector change the comparison.

For drives ending at any specified $1-\delta<1$ with $0<\delta\leq1$, a simpler uniform bound is $|L[f]|\leq(51/\delta^3)|T[f]|$. For an individual launch use $y=(x-s)/h$ on the early window, $0\leq y\leq1$. In this variable the conjugated unit-jump solution is $h$ times a solution with initial velocity one, squared frequency $p^2$ and residual potential $h^2 b(s+hy)$. The triangle acquires the factor $h^3$. The earlier Duhamel and dilution estimates apply to this unit window, with bounds at most $\eta$ on the residual potential and $2d,2e$ on the dilution difference and its $y$ derivative. Thus the positive physical early exposure is at least $h^3(19/60-2\alpha)\min(1,p^{-2})$. On the original clock the physical late exposure has absolute value at most $(16+2\alpha)B$. Since $h\geq\delta$, $\min(1,p^{-2})\geq B$, and $(16+2\alpha)/(19/60-2\alpha)<51$, the result follows. Its guarantee becomes weak for sources extending close to the front.

An inverse use of this last inequality places a bound on the last time at which the source acts. Consider the same nonzero single-sign drive and require both reported exposures to be within $1/2000$ of the model values in calibrated $A$ units. If the drive ends a distance $\delta$ before the front at $x=1$, then $\delta^3\le 51(1/2000)/(L[A]-1/2000)<(0.1853)^3$. Its support must therefore extend beyond $x=0.8147$. The clock and tolerances are in model units. This necessary condition alone cannot determine a launch history or show that an allowed later drive fits the waveform.
