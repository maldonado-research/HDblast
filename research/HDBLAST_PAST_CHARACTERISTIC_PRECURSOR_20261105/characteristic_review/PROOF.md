# A local two-characteristic precursor of the registered shell packet

5 November 2026. Independent conditional derivation. This note concerns the
registered *first* coupled master and a short interval before its supplied
Cauchy section. It constructs linear characteristic traces; it does not select
a remote physical source, complete second-order constraints, or prove a new
nonlinear shell evolution theorem.

## 1. Actual equation, scale and the available past coefficients

The retained master reference, equations (10)--(13), fixes

```
w_tt-w_zz+mathcalW(t,z)w=0,
w(0,z)=q_ell(z), w_t(0,z)=exp(-z)q_ell(z),
mathcalW=35 sech(t)^2/4-25/4+C_z/2+C^2/4+2rho^2 p^2/3.
```

Here z increases away from the shell; ell=10^-6; the smooth bump is supported
on [ell,2ell] and has maximum one. The longitudinal harmonic coefficient is
`V=alpha(z)cosh(t)^(-3/2)w`, with alpha>0 on the retained regular collar.
The full metric/scalar first solution is reconstructed from V and its derivatives.

The potential report authenticates an actual static collar 0<=z<=8ell and
the rectangle 0<=t<=3ell there, with |mathcalW|<102000. Its source derivation
also gives the exact identity

```
mathcalW(t,z)=rho(z)^2 Qbg(z)-(35/4)tanh(t)^2.
```

Consequently mathcalW is exactly even in t. Thus the *master* bound extends
to |t|<=3ell on that same static collar by the explicit formula. This is
not a transfer of the different future five-mean matrix caps to the past.
The regular background/static-collar and analytic reconstruction premises
remain retained inputs. A separate fresh interval audit can re-establish the
wide collar; no narrow historical collar is silently enlarged here.

Set x=t/ell,r=z/ell and W(x,r)=w(ell x,ell r). Then

```
W_xx-W_rr+v(x,r)W=0,
v=ell^2 mathcalW(ell x,ell r), |v|<=mu=103/10^9,
W(0,r)=Q(r), W_x(0,r)=ell exp(-ell r)Q(r),
Q(r)=exp(4-1/((r-1)(2-r))) for 1<r<2, zero otherwise.
```

The shared simple bound 103000 is conservative relative to the retained
102000 bound. The potential is smooth and even in x on the actual collar.

## 2. An actual backward linear solution without shell contact

Fix any a<1. Extend the potential boundedly and smoothly as an auxiliary
whole-line coefficient, leaving it unchanged where
1-a<=r<=2+a, -a<=x<=0. Such an extension is made with cutoffs inside the
positive collar; it asserts no physical negative-depth background. Cauchy
finite propagation confines W to

```
1-|x|<=r<=2+|x|, -a<=x<=0.
```

Thus the auxiliary extension does not affect this solution. It stays a positive
distance from the shell r=0 and within r<3. All its spatial and temporal jets
vanish on an open shell collar. The retained first master reconstruction then
gives a linear Einstein--scalar perturbation whose first physical junction rows
are identically background there. The old straight-shell/gauge identification
is still a retained analytic premise; no second canonical completion follows.

For tau=-x the backward free solution is explicitly

```
B(-tau,r)=[Q(r-tau)+Q(r+tau)]/2
         -(ell/2) integral_(r-tau)^(r+tau) exp(-ell s)Q(s) ds.
```

Let K be the usual retarded Duhamel operator for the backward time tau,
with kernel one half and the actual potential v. Then

```
W=B-KW, ||K||_C0<=beta_a=mu a^2/2,
||B||_C0<=1+a ell,
||W||_C0<=B_a=(1+a ell)/(1-beta_a),
||W-B||_C0<=beta_a B_a.
```

All beta_a<1 here. The uniformly convergent Neumann series
`W=sum_(n>=0)(-K)^n B` constructs the actual backward linear solution.
Its tail after n=N is bounded by
`beta_a^(N+1)(1+a ell)/(1-beta_a)`. Smooth coefficients and compatible
smooth compact Cauchy data give smoothness by differentiated Volterra equations.
This is an analytic construction, without a new numerical field query.

## 3. The half-time characteristic diamond and its three decisive values

Take a=1/2 and define

```
u=x+r, vnull=x-r,
D={u>=1, vnull>=-2, u+vnull<=0}.
```

Its past apex is A=(-1/2,3/2), and its future base is exactly x=0,1<=r<=2,
the whole registered packet. Its two past faces are u=1 and vnull=-2. The
diamond stays in 1<=r<=2, while the enclosing backward solution has support
in 1/2<=r<=5/2. There is no shell contact.

Put

```
beta=mu/8, Bhalf=(1+ell/2)/(1-beta), epsilon=beta Bhalf.
```

At L=(-1/4,5/4), on u=1, the free endpoint arguments are 1 and 3/2.
At R=(-1/4,7/4), on vnull=-2, they are 3/2 and 2. Hence

```
B(L),B(R)>=1/2-ell/4,
W(L),W(R)>=1/2-ell/4-epsilon >1/2-10^-6=0.499999.       (1)
```

At A the displacement endpoints are 1 and 2, both flat zeros, so

```
B(A)=-(ell/2) integral_1^2 exp(-ell s)Q(s) ds.
```

For 5/4<=s<=7/4, `(s-1)(2-s)>=3/16`; therefore Q(s)>=exp(-4/3)>1/4.
The last strict inequality follows from the rational Taylor bound exp(4/3)<4.
Also exp(-ell s)>=1-2ell. Integrating over this interval of length1/2 gives

```
B(A)<=-ell(1-2ell)/16,
W(A)<=-ell(1-2ell)/16+epsilon <-1/25000000=-4*10^-8.     (2)
```

All displayed numerical comparisons use exact fractions in the companion
arithmetic checker. The negative apex is a property of the supplied nonzero
first velocity together with the small scaled-potential remainder. Zero
initial velocity would be a different state and would not support this proof.

## 4. Compatible Goursat traces and unavoidable variation on both faces

Define the *actual* characteristic precursors by restriction of that solution:

```
F(u)=W((u-2)/2,(u+2)/2), 1<=u<=2,       [vnull=-2],
G(vnull)=W((1+vnull)/2,(1-vnull)/2), -2<=vnull<=-1, [u=1],
c=W(A)=F(1)=G(-2).
```

Their equality at A is exact, not an independently chosen zero corner.
They are constructed by the convergent Cauchy Neumann series above; for
example their free approximations are

```
F_free(u)=Q(u)/2-(ell/2) integral_u^2 exp(-ell s)Q(s) ds,
G_free(vnull)=Q(-vnull)/2-(ell/2) integral_1^(-vnull) exp(-ell s)Q(s) ds,
||F-F_free||,||G-G_free||<=epsilon.
```

These are the values of the normalized master, not the pulled-back physical
shell scalar used in the earlier second-order affine-null calculation.

In null coordinates the equation is
`4W_(u vnull)+v((u+vnull)/2,(u-vnull)/2)W=0`. It gives the exact Goursat formula

```
W(u,vnull)=F(u)+G(vnull)-c
  -(1/4) integral_1^u integral_(-2)^vnull
     v((s+t)/2,(s-t)/2) W(s,t) dt ds.                    (3)
```

Every integration rectangle stays within D. Because
`(u-1)+(vnull+2)<=1`, its area is at most1/4; the integral operator norm
in (3) is therefore at most mu/16<1. This is a unique smooth local Goursat
solution for these compatible traces, and by uniqueness reproduces exactly
the registered Cauchy packet and velocity on the future base. It does not
require data at a physical timelike boundary.

Equations (1)--(2) show

```
F(1)=G(-2)<0, F(3/2)>0, G(-3/2)>0.
```

Thus neither trace is constant. Each has an interior zero before its midpoint.
The mean-value theorem even gives a point on each first half-face with

```
F'(u)>0.999998, G'(vnull)>0.999998,
Jplus=2F'(u)>1.999996 on vnull=-2,
Jminus=2G'(vnull)>1.999996 on u=1,
```

where Jplus=W_x+W_r and Jminus=W_x-W_r. These inequalities are pointwise
existence assertions on different points of the two past faces, not uniform
lower bounds on a field and not physical energy fluxes. A prescription that
sets the relevant whole-face tangential principal derivative to zero, or sets
either face trace equal to the common corner constant, cannot produce this
packet on this diamond. This is a local finite-past obstruction to those
precise data prescriptions, not to every meaning of an incoming state.

## 5. Independent persistence argument on past Cauchy sections

For completeness there is a different norm proof, independent of the primary
L2 Duhamel bound. Put DQ=||Q_r||_L2>0 and
`N(tau)^2=||W||_2^2+||W_tau||_2^2+||W_r||_2^2`.
Endpoint Cauchy gives ||Q||_2<=DQ/2, hence

```
N(0)^2< (9/4)DQ^2,
|d(N^2)/dtau|<=(1+mu)N^2.
```

For 0<=tau<=1/4, `(1+mu)tau/2<1/4` and exp(1/4)<4/3, so N(tau)<2DQ.
Characteristic Duhamel with translation preserving L2 gives

```
||Jplus(-tau)||_2,||Jminus(-tau)||_2
 >=(1-ell/2-mu/2)DQ > (1-10^-6)DQ >0.
```

The initial lower bound uses the actual velocity and ||Q||_2<=DQ/2.
This proves no current-type single principal channel occurs on any entire
Cauchy section of that last quarter-time slab. The primary L2 Volterra route
can extend the common 10^-6 guard to every fixed tau<1 while retaining the
strictly positive shell gap. No claim at shell contact tau=1 is made here.

## 6. Physical meaning and limits

The two-face zeros/sign changes survive an arbitrary positive smooth
conjugating factor, in particular the actual factor alpha cosh^(-3/2)t
relating W to the registered longitudinal master. They also survive an
orientation-preserving reparametrization of each face. The quantitative
characteristic weights, derivative sizes and L2 constants do not survive
arbitrary such conjugations unchanged. A master normalization is not a
general bulk gauge invariant, physical stress flux, particle count or
Noether energy.

The first reconstructed fields have no shell contact on this local past
slab. That verifies their first linear shell compatibility through the open
background collar; it does not complete second-order constraints, a nonlinear
family on a chosen past domain, a retarded source, or a matter-production
sector. Earlier scattering or a remote incoming history can still feed these
two compatible traces. Choosing another earlier cut can change what data
are called incoming. The result does not establish such a history or exclude
all of them, and supplies no higher-dimensional Big Bang mechanism.
