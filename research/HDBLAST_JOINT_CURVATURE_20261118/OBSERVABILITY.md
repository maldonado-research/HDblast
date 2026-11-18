# Gravitational cost during the initial quiet interval

The registered nonlinear shell agrees with its background for the entire observer interval $0\le x\le h=19/40$, where $q=\ell x$ and $\ell=10^{-6}$. Its intrinsic curvature also agrees there. A four-dimensional scalar can be invisible to a linear detector because different fields cancel. This calculation tests what those fields do to the intrinsic Einstein-frame scalar curvature while their detector output is quiet. [Exact nonlinear shadow, §§3–4](https://github.com/maldonado-research/HDblast/blob/59651ab35b954ef268c81a771e9ec8064691f9ad/research/HDBLAST_EFFECTIVE_OBSERVER_20261116/EFFECTIVE_OBSERVER_DOMAIN_AND_RESPONSE_BOUNDS.md) · [Late exposure and normalization](https://github.com/maldonado-research/HDblast/blob/59651ab35b954ef268c81a771e9ec8064691f9ad/research/HDBLAST_RECIPROCAL_PREPARATION_20261118/KICK.md)

## A single free mode cannot stay arbitrarily small in early mean square

Use the physical four-dimensional de Sitter metric with radius $\rho_b$ and a canonical quadrupole $a(x)Y$. For squared mass $P=(m\rho_b\ell)^2\ge0$, let $v=\cosh^{3/2}(\ell x)a$, $\omega^2=P+13\ell^2/2$ and $b(x)=(35/4)\ell^2\tanh^2(\ell x)$. The mode obeys $v''+\omega^2v=b(x)v$. Let $L$ denote the uniform mean on $[5/4,11/8]$ and $e(f)=(h^{-1}\int_0^h f^2 dx)^{1/2}$. For all frequencies, amplitudes and phases,
$$
|L(a)|\le 8 e(a).
$$
The number eight is a sufficient constant. Its derivation uses the finite exposure at both ends.

Set $x_c=h/2$ and $p=\omega h$. Let $f$ solve the constant-frequency equation sharing the two exact Cauchy data of $v$ at $x_c$. For $p\le1/2$ write $f=A\cos(pz)+B\sin(pz)/p$ with $z=(x-x_c)/h$. In the early interval $z$ is uniform on $[-1/2,1/2]$, so the two terms are orthogonal. The inequalities $\cos^2(pz)\ge1-(pz)^2$ and $\sin^2(pz)/p^2\ge z^2-p^2z^4/3$ give
$$
e(f)^2\ge (47/48) A^2+(79/960)B^2 .
$$
Here the even moments of $z$ are $1/12$ and $1/80$ at orders two and four. In the late exposure $z\ge0$ and its average is $43/19$. Therefore
$$
|L(f)|/e(f)\le\sqrt{48/47+(960/79)(43/19)^2}<7.954
$$
when $e(f)>0$. The same argument on the entire interval up through $11/8$ with $|z|\le91/38$ gives $\|f\|_\infty<9 e(f)$ for a nonzero solution. For $p\ge1/2$ write the coefficients of cosine and sine without the divisor $p$. The smaller eigenvalue of their early mean-square Gram matrix is $(1-|\sin p/p|)/2\ge79/3840$. Sinc decreases on $[1/2,2]$, the alternating Taylor upper bound gives $\sin(1/2)/(1/2)\le1841/1920$, and $|\sin p/p|\le1/2$ for $p\ge2$. The pointwise ratio to the early RMS is then at most $\sqrt{3840/79}<7$. An identically zero $f$ gives zero initial data and therefore zero $v$.

For $X=11/8$, one has $0\le b(x)\le(35/4)\ell^4X^2$. Apply Duhamel’s identity forward or backward from $x_c$, using $|\sin(\omega t)|/\omega\le |t|$. With $k=(35/8)\ell^4 X^4$, the supremum of $|v-f|$ is at most $k\|f\|_\infty/(1-k)$. Also $0\le1-\operatorname{sech}^{3/2}(\ell x)\le3\ell^2 X^2/4$. Consequently $\|a-f\|_\infty\le\alpha e(f)$ with $\alpha=9[k/(1-k)+3\ell^2X^2/4]<10^{-10}$. This implies $e(a)\ge(1-\alpha)e(f)$ and $|L(a)|\le(7.954+\alpha)e(f)$. Their ratio is below eight.

## The signed curvature measurement

Take a four-dimensional Einstein model with fixed cosmological constant $\Lambda$ and coupling $\kappa>0$ defined by $G_{ab}+\Lambda g_{ab}=\kappa T_{ab}$. The real canonical fields have constant stationary backgrounds. Absorb their constant vacuum potential into $\Lambda$. Their leading perturbations are $\epsilon a_j(x)Y(n)$ and their linear equations are free in the mass eigenbasis with nonnegative squared masses. Allow independently signed initial data. The scalar detector reads $A_{\rm c}=\sum_j c_j a_j$ in the registered $A=(\ell/G)S$ units. Write $D^2=\sum_j c_j^2$. Assume $0<D<\infty$ and sufficient energy and convergence to form the trace and perform the integrations; finite smooth towers satisfy these conditions.

Let $u=n_1^2+n_2^2$, $Y=n_1^2-n_2^2$ and $w=3-4u$. Denote normalized expectation over the actual observer cut by $E_\epsilon$. On background cuts, $E_0[wY^2]=0$, $E_0[wu]=1/6$ and $|\nabla Y|^2=4(u-Y^2)$. For the minimally coupled fields the trace law is $R-4\Lambda=\kappa[\sum (\partial\chi_j)^2+4(V-V_0)]$. At quadratic order all time-derivative and potential terms have angular factor $Y^2$. The signed cut expectation leaves exactly the positive expression
$$
Q_E=\lim_{\epsilon\to0}\frac{1}{h\epsilon^2}\int_0^h E_\epsilon[w(R-4\Lambda)]dx
=\frac{2\kappa}{3\rho_b^2 h}\sum_j\int_0^h \operatorname{sech}^2(\ell x)a_j(x)^2dx .
$$
Variation of the volume density drops out because the background and linear curvature contrasts vanish. The same cancellation makes this coefficient insensitive to a smooth first-order clock and label displacement. This is a measurement of intrinsic spacetime scalar curvature on the evolving cuts. It requires access to that curvature, not only to the cut area or an uncalibrated scalar trace.

By Cauchy–Schwarz and the single-mode estimate,
$$
|L(A_{\rm c})|^2\le D^2\sum_j|L(a_j)|^2\le64 D^2\sum_j e(a_j)^2 .
$$
The shell's late normalized scalar mean exceeds 4.01324558. If the four-dimensional fit misses that mean by at most 0.0005 in the stated units, then
$$
\frac{D^2Q_E}{\kappa}\ \ge\
\frac{2}{3\rho_b^2}
[1-(\ell h)^2]\frac{(4.01324558-0.0005)^2}{64}
> 9.9363\times10^{-6}.
$$
This bound applies in the specified background-length, field and detector units. The exact shell target is $Q_E=0$. A nonzero coupling and finite fixed $D$ give a strict mismatch, though its numerical size can approach zero over models with unbounded $D^2/\kappa$. To use the threshold with finite data needs an independent bound on $D^2/\kappa$ and errors on the quadratic curvature coefficient. An observation of linear curvature alone would already test the minimally coupled model when its background sources are stationary.

If the measured metric is a known conformal function of the measured scalar times the Einstein-frame metric, the bound directly concerns the curvature reconstructed in that Einstein frame. Exact agreement of the physical metric and scalar with their constant-scalar background on the entire quiet interval, through second order with the necessary derivative control, also implies quiet reconstructed Einstein-frame curvature. An error on scalar amplitude alone does not bound spacetime derivatives in a conformal curvature transformation. Independent hidden conformal fields, additional angular sectors or stresses, later forcing and arbitrary metric readouts change this test.
