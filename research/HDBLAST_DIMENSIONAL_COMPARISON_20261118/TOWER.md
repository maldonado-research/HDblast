# A source-free four-dimensional realization of the linear shell trace

The registered first-order quadrupole can be reproduced on the complete admitted clock interval by a countable collection of free positive-mass scalar fields in four dimensions. The initial amplitudes encode the specified packet. For the two window means alone, at most two fields suffice if their preparations can be adjusted. These results concern the observed linear scalar in the selected angular sector. Physical preparation, nonlinear field dynamics and the probabilities of measured outcomes require further input.

## The radial problem and the observation

Use dimensionless coordinate time $t$, inward depth $z$, $\ell=10^{-6}$, and the background shell clock $q=-t\in[0,Q]$, $Q=25\ell/16$. At first order the full proper-cut scalar is $S(t)Y(n)$ with $Y=n_1^2-n_2^2$, $-\Delta_\kappa Y=8Y$. Write $a=35/4$ and $\beta_0=\gamma+5\mu/2$. The registered master and its boundary satisfy

$$
w_{tt}-w_{zz}+[V(z)-a\tanh^2t]w=0,\qquad
-w_z+\mu w_{tt}+[\beta_0-\mu a\tanh^2t]w=0\quad(z=0).
$$

The physical coefficient is $S(t)=G\cosh^{-3/2}t\,(w_z(t,0)-\gamma w(t,0))$ in the registered normalization. Both initial profiles are smooth and supported in $[\ell,2\ell]$. They are $w(0,z)=B(z)$ and $w_t(0,z)=e^{-z}B(z)$ for the archived bump. In elapsed past time the velocity is $-e^{-z}B$. Finite propagation gives $z\le 2\ell+|t|$. Take a fixed outer wall $Z=6\ell$, with $w(t,Z)=0$. On $|t|\le Q$ the actual master vanishes in an open neighborhood of this wall, so uniqueness gives exactly the same trace there. The artificial wall specifies the comparator spectrum. Its masses need not describe the global regular cone. [Registered equation and observable, §1](../HDBLAST_CHARACTERISTIC_ODD_TRACE_20261104/primary/PROOF.md) · [Support and energy, §§2–4](../HDBLAST_PAST_SHELL_CONTACT_20261104/primary/PROOF.md) · [Actual observer identification](../HDBLAST_SHELL_RESPONSE_20261118/BRIDGE.md)

## A positive endpoint-mass operator

In unscaled $z$, give $\mathcal H=L^2(0,Z)\oplus\mathbb R$ norm $\|(f,b)\|^2=\int_0^Z|f|^2dz+\mu|b|^2$. Use the form

$$
\mathfrak a[(f,b)]=\int_0^Z(|f'|^2+V|f|^2)dz+\beta_0|b|^2,
\qquad f\in H^1(0,Z),\quad f(Z)=0,\quad b=f(0).
$$

The retained collar bounds give $V=\rho^2 Q_{\rm bg}>129^2/3=5547$, $0<\mu<1/400$ and $\beta_0>180$. Hence $\beta_0/\mu>72000$ and $\mathfrak a\ge5547\|\cdot\|^2$. The form domain is dense. Smooth compactly supported bulk functions approximate arbitrary bulk $L^2$ data, and a smooth shrinking layer at zero can separately approximate the desired endpoint component in $\mathcal H$. The form norm is equivalent to the $H^1$ norm on the trace-linked domain by boundedness of $V$ and the trace theorem. Its embedding into $\mathcal H$ is compact by compactness in bulk $L^2$ and finite-dimensionality of the endpoint component. The closed form therefore defines a self-adjoint operator with compact resolvent. Integrating by parts identifies

$$
 A(f,b)=(-f''+Vf,\ [-f'(0)+\beta_0 b]/\mu),\qquad
 D(A)=\{(f,f(0)): f\in H^2(0,Z),\ f(Z)=0\}.
$$

The boundary derivative enters with a minus sign for inward $z$. The endpoint row belongs to the operator value. This is the finite-interval form of the project's earlier endpoint-mass construction. The earlier global cone spectrum and this auxiliary discrete spectrum use different outer domains. [Earlier operator and unscaled coefficient bounds, §§1–3](../HDBLAST_PAST_SHELL_CONTACT_20261104/primary/PROOF.md) · [Regular-cone domain and spectrum, §6](../HDBLAST_REGULAR_SHELL_SPECTRA_20261013/domain/SCALAR_SELF_ADJOINT_DOMAIN_AND_GAP.md)

## Local four-dimensional fields and their readout

Choose a real orthonormal eigenbasis $E_j=(e_j,e_j(0))$ of $A$, with eigenvalues $\lambda_j\ge5547$. For $W(t)=(w(t,\cdot),w(t,0))$, the dynamic boundary and the bulk equation have the same scalar temporal shift. Thus

$$
 W_{tt}+[A-a\tanh^2 t\, I]W=0,\qquad
y_j(t)=\langle E_j,W(t)\rangle,\quad
y_j''+[\lambda_j-a\tanh^2t]y_j=0.
$$

The amplitudes and velocities at zero are the projections of the two stated initial profiles with zero endpoint values. This defines one linear preparation map for the registered master data, with no time-dependent forcing. Introduce free scalar fields $\chi_j(t,n)=\cosh^{-3/2}t\,y_j(t)Y(n)$ on the dimensionless four-dimensional de Sitter metric $-dt^2+\cosh^2t\,\kappa$. For a field in this angular sector the minimally coupled Klein–Gordon equation is

$$
 \chi_{j,tt}+3\tanh(t)\chi_{j,t}+[8\operatorname{sech}^2t+m_j^2]\chi_j=0.
$$

After removing the dilution factor $\cosh^{-3/2}t$, its time coefficient is $m_j^2+13/2-(35/4)\tanh^2t$. Set $m_j^2=\lambda_j-13/2\ge11081/2>0$. The mass squared in the physical shell metric is $m_j^2/\rho_b^2$. Each species evolves locally from its Cauchy data. The species label indexes independent fields on the common four-dimensional spacetime.

The radial eigenvalue row gives $e_j'(0)-\gamma e_j(0)=\mu(5/2-\lambda_j)e_j(0)$. A fixed linear readout of the fields is consequently

$$
 \sum_{j\ge1}G\mu(5/2-\lambda_j)e_j(0)\,\chi_j(t,n)=S(t)Y(n).
$$

For convergence use the endpoint projection $B_0(f,b)=b$, bounded on $\mathcal H$ by $\mu^{-1/2}$. The full readout on $W$ is $G\mu\cosh^{-3/2}t\,B_0((5/2)I-A)W$ before multiplying by $Y$. It is continuous in the graph norm of $A$. Both smooth initial profiles belong to every $D(A^k)$ since repeated application of the static operator keeps their supports inside the same interior interval, with zero endpoint jets. The wave propagators of $A$ preserve each graph space. Duhamel’s formula there with the bounded scalar potential gives $W\in C([-Q,0];D(A^k))$ for every $k$. For instance the initial velocity term $A^{-1/2}\sin(t\sqrt A) W_t(0)$ is bounded in each such norm, and the Volterra iteration converges uniformly on this compact time interval. Spectral truncations converge to the identity in the $D(A)$ norm and have norm at most one there, hence converge uniformly on the compact image of $W$. Applying the bounded readout gives uniform convergence of the observed finite sums in $t,n$. This establishes the trace equality without interchanging an uncontrolled normal derivative and an infinite sum. [Related endpoint readout, §6, equation (15)](../HDBLAST_PAST_SHELL_CONTACT_20261104/primary/PROOF.md)

There is also a formulation with a bounded detector across species. Give the free canonical fields the rescaled initial amplitudes $\psi_j=(5/2-\lambda_j)\chi_j$. They satisfy the same independent Klein–Gordon equations. The detector coefficients are now $d_j=G\mu e_j(0)$. The Riesz representative of $B_0$ in $\mathcal H$ is $(0,1/\mu)$, so Parseval gives $\sum_j e_j(0)^2=1/\mu$ and $\sum_j |d_j|^2=G^2\mu$. The extra factor of $\lambda_j$ in the prepared data has a finite energy cost for the smooth packet. On the compact clock interval the dimensionless de Sitter scale factor and its reciprocal are bounded. Up to a common finite angular normalization the sum of the instantaneous free-field Cauchy energies is bounded there by a constant times $\sum_j (1+\lambda_j)^2(|y_j'|^2+(1+\lambda_j)|y_j|^2)$. This sum is uniformly finite by the graph-space solution for the displacement and velocity, using the smooth power-domain data. The detector and preparations thus have finite norms for this construction. No externally calibrated cap on either norm, nor an energy equivalence between the two physical models, is assumed.

The exact tower thus also matches the quadratic coefficient of the observed variance, $S^2/6$, and each signed linear window mean. A finite truncation with the original projected data approximates these quantities arbitrarily closely on the compact clock interval. The initial projection can be highly correlated across species. A claim about the probability of preparing these data would require a specified preparation law.

## How much information do two windows carry?

Let $L_1,L_2$ be the two normalized window means of the scalar quadrupole coefficient, continuous linear functionals on $C([0,Q])$. Let $r_j(q)$ be the signed coefficient of the $j$th output term evaluated at $t=-q$, and $v_j=(L_1r_j,L_2r_j)$. The target $v=(L_1 S(-\cdot),L_2 S(-\cdot))$ is the limit of finite sums of $v_j$. Their linear span in $\mathbb R^2$ is closed. Choose a basis from among these vectors. It has at most two members, and $v$ is a linear combination of those members with real coefficients. Multiplying the initial displacement and velocity of each selected four-dimensional species by its coefficient reproduces **both signed means exactly**, and hence their ratio when the first mean is nonzero. No numerical choice of masses, amplitudes or preparation-energy ceiling is supplied by this existence proof.

Fixing the two means leaves other observables free. More generally $k$ continuous linear summaries admit the same argument with at most $k$ freely reweighted species. Exact vanishing on an entire early time interval supplies infinitely many constraints. Each finite free sum here is real analytic in time because its differential equations have real analytic time coefficients. If that sum vanishes on an open interval, it vanishes identically on the connected clock domain. The registered linear scalar trace is exactly quiet before $q=\ell$ and nonzero later, so an exact full-waveform match in this free family requires the infinite limit. The projected truncations nevertheless approximate the complete compact-time trace to any prescribed positive continuous-norm tolerance. The gap in ideal analytic properties therefore supplies no uniform positive separation against unbounded species count at finite precision. [Quiet interval, §4](../HDBLAST_PAST_SHELL_CONTACT_20261104/primary/PROOF.md) · [Nonzero interval and cut variance](../HDBLAST_SHELL_RESPONSE_20261118/VARIANCE.md)

The equality here covers the registered classical linear scalar readout on the observer’s background clock. Agreement on its leading variance follows algebraically. Equality of nonlinear finite-amplitude spectra, induced metric observables, detector noise, physical source histories, or dynamics in other angular sectors would require separate constructions.
