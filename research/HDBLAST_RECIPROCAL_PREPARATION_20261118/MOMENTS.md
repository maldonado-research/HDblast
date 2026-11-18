# How the exact four-dimensional realization prepares a silent detector

The known four-dimensional realization uses a source supported in the interior of the radial collar and a readout at its endpoint. Their separation has a spectral signature. In that realization the signed displacement and velocity couplings each change sign infinitely often as the mass increases. For each sequence, summing its absolute weights after multiplying them by any growing exponential of mass gives an infinite result. These are statements about the ideal linear response, on the registered clock.

## The separated source and endpoint

Use the positive self-adjoint collar operator $A$, eigenpairs $(\lambda_j,E_j)$, finite endpoint mass $\mu$ and bounded endpoint evaluation $L(f,b)=b$ of the earlier construction. The real eigenfunctions have $LE_j=e_j(0)$. Order them by $\lambda_j$. The endpoint equation reads $-e'_j(0)+\beta_0 e_j(0)=\mu\lambda_j e_j(0)$. If $e_j(0)=0$, both Cauchy data for the second-order radial eigenvalue equation would vanish. Thus every endpoint value is nonzero. The same boundary equation fixes the radial solution up to a scalar for a given eigenvalue, so all eigenvalues are simple. They tend to infinity by compact resolvent. [Operator, initial data and readout](https://github.com/maldonado-research/HDblast/blob/20f431f2f6fd8ae9d040d0e8b87bddb41a2f3fa2/research/HDBLAST_DIMENSIONAL_COMPARISON_20261118/TOWER.md)

Let $F$ denote either the nonzero initial displacement $(B,0)$ or the initial velocity $(-e^{-z}B,0)$ in elapsed $q$. Up to common normalization, define
$$
w_j=e_j(0)(5/2-\lambda_j)\langle E_j,F\rangle.
$$
Multiplication by the clock and detector constants leaves the following zero and sign-change statements invariant. The factor $5/2-\lambda_j$ never vanishes. Completeness and the nonvanishing endpoint projections imply that each of the two choices gives a nonzero sequence.

For every nonnegative integer $k$ one has
$$
\sum_j |w_j|\lambda_j^k<\infty,\qquad
\sum_j w_j\lambda_j^k=L A^k((5/2)I-A)F=0 .
$$
The first bound follows by Cauchy–Schwarz from $\sum_j|e_j(0)|^2=1/\mu$ and $F\in D(A^{k+1})$. For the second, repeated applications of the local radial differential operator preserve the interior support of each smooth profile. The endpoint of the resulting pair is zero. Thus the detector is orthogonal to every polynomial spectral weighting of each prepared vector. This conclusion uses support information available before evaluation of the later waveform.

If a nonzero sequence had only finitely many sign changes after dropping its zero entries, choose one polynomial root strictly between the eigenvalues at each reversal. The overall sign of the resulting polynomial can be chosen so every nonzero $w_j p(\lambda_j)$ is positive. Its convergent sum would be positive, contradicting the vanishing moments. Therefore both signed sequences have infinitely many reversals. The product $e_j(0)\langle E_j,F\rangle$ is unchanged when the sign of an eigenfunction is reversed.

For any $a>0$ both choices obey $\sum_j |w_j|e^{a\sqrt{\lambda_j}}=\infty$. Otherwise $C(s)=\sum_j w_j\cos(s\sqrt{\lambda_j})$ would be analytic in $|\operatorname{Im}s|<a$. Its even derivatives at zero vanish by the moment identities and its odd derivatives vanish by evenness. It would therefore be identically zero on the real axis. But the long-time mean of $2 C(s)\cos(s\sqrt{\lambda_h})$ equals $w_h$ for every $h$, by distinct positive frequencies and absolute summability. This would force the sequence to vanish. The distinction between $\sqrt{\lambda_j}$ and the physical mass is an affine shift in squared mass and a fixed change of units. Existence of some exponential moment is invariant under this change.

## Other smooth four-dimensional preparations

For the transformed Klein–Gordon amplitude $v_j=\cosh^{3/2}(\ell x)a_j$, write $v_j''+(P_j+h(x))v_j=0$, where $h(x)=\ell^2[(35/4)\operatorname{sech}^2(\ell x)-9/4]$ is even and $P_j=(m_j\rho_b\ell)^2\ge0$. With fixed readout $c_j$, group the real weights $c_jv_j(0)$ and $c_jv'_j(0)$ at each distinct squared mass. Only these grouped weights enter the scalar signal. Assume their series permit differentiation term by term to every order near zero. Absolute polynomial moments of every order for the grouped absolute weights suffice, by the mode equation and bounded real-time oscillator propagation on the compact clock interval.

If the observed sum is silent on the initial open clock interval, all of its one-sided derivatives at zero vanish, also after multiplication by $\cosh^{3/2}(\ell x)$. The mode equation and evenness of $h$ give $v_j^{(2k)}(0)=p_k(P_j)v_j(0)$ and $v_j^{(2k+1)}(0)=r_k(P_j)v'_j(0)$ with degree $k$ and leading coefficient $(-1)^k$ in each polynomial. This follows by differentiating the equation successively. Triangular elimination yields vanishing polynomial moments separately for the two grouped measures. A later nonzero signal requires at least one nonzero grouped measure. For that measure, a finite partition into ordered mass intervals of alternating fixed sign with separating roots between sign blocks is impossible. For discrete masses ordered without finite accumulation, this means infinitely many sign changes. The cosine-transform argument also precludes an absolute exponential moment in $\sqrt P$ for the active measure. At $P=0$ the zero-frequency coefficient is extracted by the mean itself.

For example, this excludes an exact delayed positive signal from a collocated displacement-only preparation $v_j(0)=g c_j$, or from a velocity-only kick $v'_j(0)=g c_j$ with zero displacement and nonzero gain. For these two collocated examples, the initial value or the initial derivative already suffices when the indicated sum exists. The moment argument additionally covers smooth preparations with a finite number of reversals in modal sign or exponential damping of the signed modal weights. The sign restriction applies to mass-grouped preparation–receiver products. Independent dimensions, the signs of individual eigenfunctions, or arbitrary coefficient representations are not measured by it.

A finite-precision test needs further numerical restrictions. A finite subset of the actual tower can approximate the delayed signal uniformly on this compact clock interval even though each truncation has finitely many sign changes and finite exponential moments. Thus exact silence and infinite moment tests alone give no positive uniform error margin over models with arbitrary budgets. The result locates the difference between an independently prescribed collocated source and the exact separated-source mode construction. It supplies no calibration law for an astronomical source.

The square-root exponential condition also appears in the classical Stieltjes moment problem. Here the two positive parts of a signed mass measure have the same polynomial moments. The proof above works directly with the discrete shell spectrum. [Yarovaya, Stoyanov and Kostyashin, §5, equation (42)](https://arxiv.org/html/1912.00160v1#S5.E42)
