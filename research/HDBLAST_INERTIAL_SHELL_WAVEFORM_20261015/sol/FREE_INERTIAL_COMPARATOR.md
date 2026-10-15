# Exact free inertial comparison wave

15 October 2026 UTC. Conditional linear comparison calculation for the new
width ell=10^-6, T=3ell and the retained branch. This note and its checker do
not transfer the narrow-state response coefficient, establish a nonlinear
epsilon radius or give an observational time scale.

Let f(z)=Q(z/ell), g(z)=exp(-z)f(z), gamma=C_b/2,
beta0=gamma+5mu/2, and j0(t)=f'(t)+g(t). The free comparison wave w0 solves

    w0_tt-w0_zz=0,
    -w0_z+mu*b0''+beta0*b0=0, b0(t)=w0(t,0),
    w0(0,z)=f(z), w0_t(0,z)=g(z).

Its boundary characteristic is w0_t+w0_z=j0, hence

    mu*b0''+b0'+beta0*b0=j0,
    D0=j0-b0'-gamma*b0.

Let k be the impulse kernel, satisfying

    mu*k''+k'+beta0*k=0, k(0)=0, k'(0)=1/mu.

Because f is smooth and zero before ell, integration by parts is exact:

    b0(t)=integral_0^t [k'(t-s)+exp(-s)k(t-s)]f(s) ds,
    b0'(t)=f(t)/mu
      +integral_0^t [k''(t-s)+exp(-s)k'(t-s)]f(s) ds.

Write a=1/mu, c=beta0/mu=gamma*a+5/2 and

    P(r,s)=(a-gamma-exp(-s))*k'(r)
             +(c-gamma*exp(-s))*k(r).

The exact normalized trace is therefore

    H0(x)=Q'(x)+ell*(exp(-ell*x)-1/mu)*Q(x)
      +ell^2*integral_1^min(x,2) P(ell*(x-y),ell*y)*Q(y) dy.

This expression avoids quadrature of Q'. The velocity data are retained in
both the local exponential and the kernel. In particular the leading tail
factor is `(1-mu*(gamma+1))/mu^2`; omitting its `+1` would discard the
initial-velocity contribution.

For the authenticated shell intervals, a-gamma-1>0. The validated kernel
has k'>0 and k>=0 on [0,T]. Thus P>0 on the needed rectangle and the
comparison tail is positive for every x in [2,3]. The enclosure does not
assume a tail sign for the actual coupled solution; that requires its
separate residual estimate.

The standalone module uses 112-bit outward integer fixed-point arithmetic.
For y=(1000k,k'), the matrix infinity norm is at most 1000 on the whole
parameter box. Its degree-eight Taylor enclosure has remainder, in k', at
most

    2*a_upper*(1000*r_upper)^9/9!.

The factor two is rigorous because 1000T=3/1000<1/2 and
exp(1000T)<=1/(1-1000T)<2. A separate bound

    k'(r)>=a_lower-a_upper*(1000T)/(1-1000T)>0

holds throughout [0,T] and justifies k>=0 even when evaluating a subinterval
whose left endpoint is positive. Every arithmetic operation rounds outward;
the exponential routine retains a nonzero upper bound for very small tails.

The bump primitive is enclosed by 4096 continuous interval rectangles.
Since Q>=0, its primitive is monotone. On a time cell [l,r], the integral
term is enclosed by the positive kernel range on
`0<=delay<=ell*(r-1)`, `ell<=s<=ell*min(r,2)`, multiplied by the primitive
range from l to r. This avoids a quadratic quadrature loop. Preparation took
about 3.6 seconds and evaluating all 4096 time cells took about 13 seconds
in one initial benchmark; timings are not certificate inputs.

The reproduced continuous comparison-tail interval is approximately

    3.23575026180407e-8 <= H0(x) <= 5.27014430788258e-8,
    for every 2<=x<=3.

Only the rational endpoints in the JSON decide bounds. The displayed
decimals are summaries. The comparison center is negative, approximately
`[-0.000484475,-0.000439180]`, reflecting the retained local inertial term.

There is also a simpler independent analytic sign check. The retained shell
box implies `440<=a<=486`, `0<gamma<=248`. The same short-time bound gives
`k'>438`, so `P>(440-248-1)*438=83658`. On
`11/8<=y<=13/8`, the bump satisfies
`Q(y)>=exp(-4/15)>=1-4/15=11/15`. Therefore its integral is at least
`11/60`, and the free comparison tail is greater than
`ell^2*83658*11/60=1.53373e-8`. This weaker sign bound requires no
quadrature and retains the velocity contribution.

The module API is `load_shell(source_root)`, `BumpMass(count=4096)` and
`comparator_cell(left,right,mu,gamma,mass)`. Returned intervals are integer
pairs divided by `S=2**112`. `q_cell(left,right)` returns continuous Q,Q'
ranges. `kernel` and `positive_kernel` are also available for review.

Reproduce from any working directory, with an extracted prior capsule or a
minimal tree containing the one pinned shell JSON at its recorded path:

```sh
python -I -B inertial_comparator.py --source-root /path/to/capsule
python -I -O -B inertial_comparator.py --source-root /path/to/capsule
python -I -B check_comparator.py --source-root /path/to/capsule
```

The checker creates a relocated one-input fixture, compares fresh normal
and optimized outputs, and rejects four mathematical faults in both modes:
reversed impulse initial velocity, reversed tail kernel sign, omitted local
inertia and reversed bump derivative. Changed authenticated input is also
rejected. All checks use explicit exceptions and survive Python optimization.
No archived producer is imported or executed.
