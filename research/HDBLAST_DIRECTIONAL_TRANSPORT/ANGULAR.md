# Directional and repeated-state bounds for a sharp two-body source

The HDblast calculation compares a prescribed daughter birth rate with the total flux of states through a moving momentum shell. The same continuum approximation allows three further tests. Can uniform births fit at each angle? Can another angular distribution preserve the assumed birth momentum? How many crossing states have already been used? The derivation uses the specified field history, sharp birth momentum and collisionless translation in comoving momentum. Scattering and finite spectral width change the state map. [Published sweep calculation](https://github.com/maldonado-research/HDblast/blob/dc99ff5c5894abcff3e11e02e4b97e19da31de00/research/HDBLAST_COHERENT_DECAY_20261120/kinetics/SWEEP.md)

## Instantaneous directional bound

The shell in conserved momentum coordinates has center $b(t)e_z$, radius $k(t)>0$ and outward unit vector $n$. Its normal speed is $\dot k+\dot b u$, with $u=n_z$. Reflecting $u$ when necessary permits us to use $c=|\dot k|$, $d=|\dot b|$ and $g(u)=|c+du|$ for the angular integrals. Write $\langle f\rangle=\tfrac12\int_{-1}^1 f(u)du$. With the published spin, charge and flavor multiplicities, the state density is $C=12/(4\pi)=3/\pi$ per unit coordinate momentum volume. Allowing vacant states on either side gives instantaneous capacity $12k^2\langle g\rangle$. The total requested rate $P$ defines $h=P/(12k^2)$. For uniform emission the minimum shortage in daughters per proper time is

$$
D_{\rm iso}=12k^2\langle(h-g)_+\rangle
\ \geq\ (P-12k^2\langle g\rangle)_+ .
$$

The inequality follows from convexity of the positive part. For $d>0$, the angle-averaged shortage is $[G(c+d)-G(c-d)]/(2d)$, where $G(x)=[(x+h)_+^2-2x_+^2+(x-h)_+^2]/2$. Its derivative is the triangle $(h-|x|)_+$. Integrating the overlapping pieces of this triangle directly avoids cancellation at large $d/h$. For $d=0$ the average is $(h-c)_+$. At $d=c>0$ and $0<h\leq 2c$ the missing fraction is $h/(4c)$. Whenever $d>c$ and $h>0$, tangency gives a nonzero interval of angles with a deficit.

If $\langle g\rangle\geq h>0$, the smallest possible total-variation distance from uniform emission among angular densities $w$ satisfying $\langle w\rangle=1$ and $hw\leq g$ is precisely $\langle(1-g/h)_+\rangle$. Every feasible distribution must move at least this much probability out of crowded angles. Retaining $\min(1,g/h)$ and redistributing the missing probability in proportion to $(g/h-1)_+$ attains the bound whenever redistribution is required. This can change the directed momentum at birth.

## Which angular sources preserve the model's birth momentum?

The charged-stream closure injects zero directed momentum into each charge class. Back-to-back two-body kinematics guarantees zero momentum for the pair as a whole. A charge-dependent forward bias could still inject oppositely directed momenta into the individual streams. Preserving each stream's source moment requires $\langle uj\rangle=0$, where $0\leq j(u)\leq g(u)$. Azimuthal averaging preserves this constraint and the total rate. The maximum possible rate is $12k^2F_0$. When $cd=0$, the capacity is even and $F_0=\langle g\rangle$. Otherwise the unique $t\in(0,1)$ solving the first equation gives the maximum in the second,

$$
 \int_{-1}^{t}u\,|c+du|\,du=0,\qquad
 F_0=\tfrac12\int_{-1}^{t}|c+du|\,du .
$$

The first integral up to zero is negative. Its value at one is positive because $g(u)>g(-u)$ for $0<u\leq 1$. The integral increases strictly for positive upper limit. For each admissible $j$ the zero-momentum condition gives $\langle j\rangle=\langle(1-u/t)j\rangle\leq\langle(1-u/t)_+g\rangle=F_0$. Taking $j=g$ for $u<t$ and zero otherwise attains equality. Scaling this choice attains every smaller nonnegative rate. At $d=c>0$ one obtains $t=1/2$ and $F_0=9c/16$, whereas $\langle g\rangle=c$.

If each stream's distribution must also be even in $u$, its maximum density is $\min(g(u),g(-u))$. The rate ceiling becomes $12k^2F_{\rm even}$ with

$$
 F_{\rm even}=\int_0^1|c-du|du=
 \begin{cases}c-d/2,&0\leq d\leq c,\\d/2-c+c^2/d,&d>c.\end{cases}
$$

At $d=c$, $F_{\rm even}=c/2$. Rejecting only the locally overfilled part of a uniform request generally leaves an anisotropic subset with directed momentum. Its remaining count can therefore exceed the count attainable while preserving zero birth momentum.

## Exact collisionless counting across repeated crossings

For one charge sign, identify the conserved label $q=b(t)e_z+k(t)n$. Use signed $\dot b,\dot k$ and signed $u=n_z$ here. The map $(t,n)\mapsto q$ has unsigned volume Jacobian $k^2|\dot k+\dot b u|$ relative to $dt\,d\Omega$. On a fixed time window, assume a piecewise smooth history with finitely many regular preimages almost everywhere and no positive source measure mapped onto a set of zero volume. Changing variables with multiplicity gives the requested population per available label,

$$
 f_{\rm req}(q)=\sum_{t_i:\ |q-b(t_i)e_z|=k(t_i)}
 \frac{h(t_i)}{|\dot k(t_i)+\dot b(t_i)u_i|},\qquad
 u_i=\frac{q_z-b(t_i)}{k(t_i)} .
$$

An isolated tangency can produce an integrable density singularity. A stationary shell with positive births over a nonzero interval instead deposits singular source measure on a surface. All such requests exceed its zero three-dimensional continuum capacity and must be counted separately. Opposite charges follow reflected maps, giving identical availability for the two members of a back-to-back pair in the symmetric reference model. Multiplication by $C$ gives the published count of individual daughters.

With vacant states at the start, the amount requested in excess of Pauli capacity is exactly

$$
 D_{\rm mem}=C\int_{\mathbb R^3}[f_{\rm req}(q)-1]_+\,d^3q .
$$

This is the smallest cumulative count that must be rejected or reassigned when partial occupations may be accepted up to one per label. Accepting arrivals until each label fills attains it in this counting model. The formula omits quantum interference between passages and the feedback of rejected births. If initial occupations are known, replace one by $1-f_{\rm initial}(q)$. For nonnegative $x_i$, $[\sum_i x_i-1]_+\geq\sum_i[x_i-1]_+$. Applying the Jacobian formula yields $D_{\rm mem}\geq\int D_{\rm iso}\,dt$. Supplying all requested births on the original labels requires at least $D_{\rm mem}$ cumulative vacancy creation in overrequested labels, counting vacancies in particle units.

A cheaper bound uses the swept union $U$. Its capacity is $C\,\operatorname{Vol}(U)$, so $D_{\rm mem}\geq[\int Pdt-C\,\operatorname{Vol}(U)]_+$. To calculate the union, fix cylindrical radius $r$. On each connected time interval with $k(t)\geq r$, the two possible heights are $z_\pm(t)=b(t)\pm\sqrt{k(t)^2-r^2}$. Each continuous branch covers the interval between its least and greatest height. Merge these intervals across both branches and all connected time intervals. Calling the union length $\ell(r)$ gives $\operatorname{Vol}(U)=2\pi\int_0^\infty r\,\ell(r)dr$.

For fixed $k>0$ and a continuous center with total range $B=\max b-\min b$, the result is

$$
 \operatorname{Vol}(U)=
 \begin{cases}
 2\pi k^2B-\pi B^3/12,&0\leq B\leq2k,\\
 \pi k^2 B+4\pi k^3/3,&B\geq2k .
 \end{cases}
$$

The height intervals have length $B$ with centers separated by $2\sqrt{k^2-r^2}$, hence union length $B+\min(B,2\sqrt{k^2-r^2})$. Integration gives the expression including continuity at $B=2k$. Backtracking inside the same range adds zero new states. Even monotonic translation can encounter the same label twice on different sides of the sphere. With fixed center and strictly increasing radius, each label in the spherical layer has one crossing and the two deficit calculations coincide.

A two-pass control isolates the effect of memory. Hold the center fixed, expand from $k_1>0$ to $k_2>k_1$, then retrace the radii. Take $h=3|\dot k|/4$ on both passes, with zero source at the turnaround. Each instantaneous isotropic demand uses 75% of its local capacity, so $\int D_{\rm iso}\,dt=0$. Almost every label in the layer is requested twice, and its total requested occupation is $3/2$. A third of the requested births must then be rejected or supplied with vacancies. The union bound is exact for this control.

There is also a finite-width union bound without a specified line shape. For support $k_-(t)\leq |q-b(t)e_z|\leq k_+(t)$, use the union of these spherical bands. When $r\leq k_+$, let $s_+=\sqrt{k_+^2-r^2}$ and $s_-=\sqrt{(k_-^2-r^2)_+}$. At that instant the accessible heights comprise $b+[-s_+,-s_-]$ and $b+[s_-,s_+]$. Their merged union over time gives the same cylindrical volume formula. Any width leaving $C\,\operatorname{Vol}(U)<\int Pdt$ fails the total-count test. Sufficient capacity alone does not determine the production line shape or the occupied distribution.
