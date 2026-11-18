# A shape test for one freely propagating four-dimensional scalar

The measured quantity here is the coefficient of a specified angular pattern on a three-dimensional spatial cut of the model shell. The earlier checkpoint identified its actual linear waveform on a common observer clock. This comparison asks which four-dimensional equations could produce that waveform. The scope is a canonical, minimally coupled scalar on the *same* de Sitter background, a fixed clock registration, a constant mass and a linear readout with common gain. [Field and clock identification](../HDBLAST_SHELL_RESPONSE_20261118/BRIDGE.md)

## A gain-independent ratio

Use $q=-\tau/\rho_b$, $x=q/\ell$, $\ell=10^{-6}$ and $0\le x\le25/16$. Let $E_\epsilon$ denote the normalized actual cut area and transport $Y=n_1^2-n_2^2$ with the original labels. Define the angular regression slope

$$ b_\epsilon=\operatorname{Cov}_\epsilon(Y,\Phi_\epsilon)/\operatorname{Var}_\epsilon(Y). $$

The denominator is positive on each complete cut since $Y$ is nonconstant and the area density is positive. For a readout $g\Phi_\epsilon+d(q)$ the slope is exactly $g b_\epsilon$ for any nonzero constant $g$ and any angularly uniform offset $d$. The background field is constant and $E_0[Y^2]=1/6$, so the derivative of the slope at zero amplitude is $S(-q)$. The established parameter regularity makes $b_\epsilon/\epsilon\to S(-q)$ uniformly. Hold the template and clock fixed throughout this comparison. [Cut regularity, symmetry and variance](../HDBLAST_SHELL_RESPONSE_20261118/VARIANCE.md)

Let $W_1=[5/4,11/8]$ and $W_2=[49/32,25/16]$ in $x$. Divide each window integral by its own width before taking $R=\overline{S}_{W_2}/\overline{S}_{W_1}$, where the integrand is $S(-\ell x)$. To preserve correlations in the old enclosures, write

$$ u(x)=\ell\cosh^{3/2}(\ell x)S(-\ell x)/G,\qquad G=3/(\rho_b^3p_b)<0. $$

The certified **future** master intervals contain $H(x)=\ell\cosh^{3/2}(\ell x)S(+\ell x)/G$. Their `H_lo,H_hi` columns already include the discrepancy between the comparator and the actual master. The past–future transfer theorem gives $|S(-q)-S(+q)|<1/40$ on this domain. Since $-G\ge2493872451/250000000000$ and $\log\cosh q\le q^2/2$, one has $1-10^{-11}<\operatorname{sech}^{3/2} q\le1$ and $|u-H|<3\times10^{-6}$. Thus each continuous time cell gives an enclosure for $u$. The same $G/\ell$ multiplies the two physical window means and cancels exactly. Multiplying each $u$ interval by the positive $\operatorname{sech}^{3/2}$ interval, averaging by overlap and dividing intervals yields

$$ -0.36754 < R < -0.36459. $$

The first window mean is nonzero, so the limiting ratio is well defined. Uniform convergence and the strict margins also imply inclusion in this open band for a sufficiently small nonzero finite amplitude of the registered family. A common gain cancels for each amplitude. The numerical size of that amplitude neighborhood is unknown. A change in relative gain or clock calibration between windows would change the test. [Physical gain and continuous cells, §§1,4](../HDBLAST_INERTIAL_SHELL_WAVEFORM_20261015/INERTIAL_SHELL_WAVEFORM.md) · [Past transfer, §§1,6](../HDBLAST_CHARACTERISTIC_ODD_TRACE_20261104/primary/PROOF.md) · [Exact ratio](SHAPE.json)

## Two integrals of the free field equation

For a four-dimensional field $a(x)Y$ on $h=\rho_b^2(-dq^2+\cosh^2q\,\kappa)$, put $v=\cosh^{3/2}(\ell x)a(x)$ and $P=(m\rho_b\ell)^2$. The angular eigenvalue is eight. Removing the $3\tanh q$ time derivative from the Klein–Gordon equation gives, wherever the projected source vanishes,

$$ v''+[P+\ell^2 k(x)]v=0,\qquad k(x)=(35/4)\operatorname{sech}^2(\ell x)-9/4. $$

A common nonzero gain between the trial field and the target slope is allowed. If their waveforms coincide where tested, $v$ is proportional to $u$ there. For any compactly supported $C^1$ test $\psi$ below, two integrations by parts identify the same constant $P$ through

$$ K_\psi[u]=\frac{-\int \psi'' u}{\int\psi u}-\ell^2\frac{\int k\psi u}{\int\psi u}=P. $$

Choose $h_*=1/32$. The weight $\psi$ is the convolution of three indicators of $[0,h_*]$, shifted to the start of the test. It is a nonnegative piecewise quadratic with integral $h_*^3$. Its second derivative on the three adjacent subintervals takes the values $1,-2,1$. Both the weight and its first derivative are zero at the outer endpoints. Internal derivative jumps occur only in the *second* derivative, so the integrations by parts introduce no point masses. All integrals use the waveform values directly. No differencing of noisy samples is needed.

| Test | Support in $x$ | Outward enclosure for the required $P$ |
| --- | --- | --- |
| A | $[5/4,43/32]$ | $(73.97,84.65)$ |
| B | $[11/8,47/32]$ | $(45.64,60.23)$ |

The normalized $u$ is uniformly positive on both supports. The replay multiplies the 384 actual cell intervals in each support by the exact integrals of $\psi$, or by the signed lengths for $-\psi''$. It uses $-9/4\le k\le13/2$ in the positive weighted average. Each ratio is independent of a common multiplicative factor, including its sign. For every single constant $P$ the larger of the two absolute residuals $|P-K_\psi[u]|$ exceeds **6.87**. Different initial position and velocity, or a drive between the supports, cannot change the free equation on the two supports. This exclusion concerns the shape of the leading waveform. The theory still permits a time-dependent mass, forcing on a test support, other readouts or additional dynamical fields. [Exact integral calculation](shape.py) · [Rational endpoints](SHAPE.json)

As a deterministic stability check, widen each normalized $u$ enclosure by $1/2000$ on both sides. The required mass intervals stay inside $(73.47,85.15)$ and $(44.76,61.11)$, giving a minimum maximum residual above **6.18**. This models a uniform error on the normalized waveform at fixed clock. It supplies neither an instrumental noise distribution nor an uncertainty on a cosmological measurement. Because the weights involve only bounded integrals and the denominators stay away from zero after normalization, the strict separation also persists for the normalized slope of the registered family at sufficiently small nonzero amplitude on the same comparison clock. This last consequence gives no evaluated amplitude radius.

## Several fields with restricted masses or cancellation

The two weights also give a necessary condition for many free fields. Suppose their dimensionless squared masses $P_j$ lie in one interval of width $\Delta$, and their readout is $v=\sum_j c_jv_j$ with constant coefficients. For each test $i=A,B$, put $d_{ij}=\int\psi_i c_jv_j$, $D_i=\sum_jd_{ij}\ne0$ and $C_i=\sum_j|d_{ij}|/|D_i|$. This is the amount of cancellation among the signed contributions to that window. A value of one means they have the same sign. Use finite sums or assume absolute convergence and sufficient convergence to apply the weak equation term by term. Each component is free on both test supports.

The common geometric term in the equation is linear in the field. Consequently $K_i[v]=\sum_j P_j d_{ij}/D_i$. Subtract the midpoint of the allowed mass-squared interval from each $P_j$ and use the triangle inequality in the two tests. This gives

$$|K_A[v]-K_B[v]|\le\tfrac12\Delta(C_A+C_B)\le\Delta\max(C_A,C_B).$$

For a mixture matching the registered trace, the exact numerical intervals require $\Delta\max(C_A,C_B)>13.74$. With the stated uniform waveform error the necessary bound is still $>12.37$. In particular, if the modal contributions have a common sign in each test, any number of free fields needs a squared-mass spread greater than 13.74 at the ideal registered waveform. A narrow allowed mass range would require compensating contributions of different signs. Applying this inequality to data requires an independently justified restriction on the mass spread or on the allowed cancellations. The unrestricted tower imposes neither restriction. [Exact interval endpoints](SHAPE.json) · [Four-dimensional field tower](TOWER.md)
