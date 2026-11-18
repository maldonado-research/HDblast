# A four-dimensional comparison with a controlled energy cost

One can represent the linear shell signal with four-dimensional canonical scalar fields whose total initial mode energy is within a factor $1+3/(2\cdot5547)$ of the collar's initial mode energy. This representation uses a local detector that responds to acceleration as well as position and velocity. Its output has unbounded sensitivity to high temporal frequencies at an instant. Each fixed nonzero exposure width has a finite readout bound. The distinction matters when assigning the source and detector budgets in a test of the shell.

## Preparation and energy

Use the finite collar and its real eigenmodes from the published comparison. Let $q=-t$ be elapsed de Sitter time, $T(q)=\tanh q$ and $f(q)=\cosh^{-3/2}q$. The Hilbert-space-valued displacement has spectral coefficients $y_j(q)=\langle E_j,W(-q)\rangle$. The quadrupole amplitudes of independent canonical four-dimensional fields are $a_j=f y_j$. They obey
$$
a_j''+3T a_j'+[m_j^2+8\operatorname{sech}^2 q]a_j=0,
\qquad m_j^2=\lambda_j-13/2 .
$$
Primes in this note denote $q$ derivatives. Use the angular normalization in which the common integral of $Y^2$ is divided out, and the dimensionless background metric and canonical coefficients of the compared equations. Denote collar Cauchy data by $(U,V)=(W(0),\partial_q W(-q)|_{q=0})$. Then
$$
2E_B=\|V\|^2+\|A^{1/2}U\|^2,\qquad
2E_4(0)=\sum_j(|a_j'(0)|^2+[m_j^2+8]|a_j(0)|^2)
       =2E_B+\tfrac32\|U\|^2 .
$$
The identity follows from $f(0)=1$, $f'(0)=0$ and Parseval. In particular $E_B\leq E_4(0)\leq[1+3/(2\cdot5547)]E_B$. This is uniform on the natural finite-energy data $U\in D(A^{1/2})$, $V\in\mathcal H$. Conversion from these two model quadratic forms to physical energies in common units would require both actions, their normalizations and the model length scales. [Radial operator, normalization and mass dictionary](https://github.com/maldonado-research/HDblast/blob/20f431f2f6fd8ae9d040d0e8b87bddb41a2f3fa2/research/HDBLAST_DIMENSIONAL_COMPARISON_20261118/TOWER.md)

## Detection

Let $h_j=\sqrt\mu\, e_j(0)$. Endpoint Parseval gives $\sum h_j^2=1$. A receiver can first form the collective coordinate $X(q)=\sum_j h_ja_j(q)$, using fixed square-summable couplings. Its local differential readout is
$$
S(-q)=G\sqrt\mu\,\mathcal L X(q),\qquad
\mathcal L X=X''+3T X'+(4-8T^2)X .
$$
Indeed the four-dimensional equation gives $\mathcal L a_j=(5/2-\lambda_j)a_j$, so this observable exactly reproduces the registered endpoint derivative and shell field on the admitted clock interval. The change from $t$ to $q=-t$ preserves the friction sign because $\tanh$ is odd. The smooth registered packet supports all the required derivatives. Smooth data in the collar's graph spaces also give pointwise convergence, while the identity for rough energy data can be read distributionally by taking finite spectral sums. This assertion concerns the finite-collar model. Its identification with the untruncated shell for a given physical solution requires the same domain-of-dependence conditions as the preceding edition.

For a smooth exposure weight $\eta$ compactly supported inside the clock domain, integration by parts yields
$$
\int \eta S=G\sqrt\mu\int X\,[\eta''-(3T\eta)'+(4-8T^2)\eta],
$$
where $S$ under either exposure integral means $S(-q)$. Thus the exposure can use the ordinary position readout $X$ with a fixed, signed temporal filter. Define $E_4(q)$ with the time-dependent mode stiffness $m_j^2+8\operatorname{sech}^2q$. Multiplying the field equations by $a'_j$ gives the sum of nonpositive terms $-3T|a'_j|^2-8\operatorname{sech}^2(q)T|a_j|^2$ in its derivative for $q\geq0$. Writing $m_*^2=5547-13/2$, one obtains $|X'|\leq\sqrt{2E_4(0)}$ and $|X|\leq\sqrt{2E_4(0)}/m_*$. These estimates and the integrated identity pass to energy data by spectral truncation. In particular a smooth exposure is bounded by $|G|\sqrt{2\mu E_4(0)}\,\|\eta''-(3T\eta)'+(4-8T^2)\eta\|_1/m_*$. For a normalized rectangular exposure of width $L>0$, integrating the acceleration once and bounding its boundary velocities yields
$$
\left|\frac1L\int_{q_1}^{q_1+L}S(-q)\,dq\right|
\leq |G|\sqrt{2\mu E_4(0)}
\left(\frac2L+3+\frac4{m_*}\right).
$$
The response to independent controls would have to use the same preparation and temporal filter. An experimental realization would require measured bandwidth, gain, timing and noise properties for that filter. The source couples to interior data and the receiver to the endpoint. Their couplings can therefore differ.

## Why instantaneous detection changes the comparison

Write $b_j=G\mu(5/2-\lambda_j)e_j(0)$ for the instantaneous displacement readout coefficients. The Riesz representative $R=(0,1/\mu)$ of endpoint evaluation belongs to $\mathcal H$ but lies outside $D(A^{1/2})$, whose endpoint component must equal the trace of the bulk function. Hence $\sum_j\lambda_j e_j(0)^2=\infty$ by the form-domain characterization. Since all $\lambda_j\geq5547$, it follows that $H_N=\sum_{j\leq N}|b_j|^2/\lambda_j$ tends to infinity. Set $U_N=\sum_{j\leq N}b_jE_j/(\lambda_j\sqrt{H_N})$ and $V=0$. These smooth finite sums have $\|A^{1/2}U_N\|=1$, yet their initial readout is $\sqrt{H_N}$. Thus a bounded position-only detector and a preparation map uniformly bounded in canonical energy cannot jointly reproduce the instantaneous derivative output for every collar datum. The example uses unrestricted boundary-compatible data at $q=0$. It imposes no lower bound for the single compactly supported interior packet. The previous construction rescales the canonical modes by $5/2-\lambda_j$ and keeps a position-only detector. Its cost is finite for that smooth packet. [Previous bounded-detector construction](https://doi.org/10.5281/zenodo.24351909)
