# Independent piecewise-ridge audit of the candidate kink

This check uses the published `RESPONSE_PROFILES.json` (SHA-256
`00165ede35939b768150a9cb9b04cd0b01f0ec6be6d54c8dbcf1b37b90959fdc`) as
exact rational data and imports no project code. It reconstructs F and G,
including the single and double integrals of the P1 velocity, by direct
segment integration, independently of the frozen ridge implementation.
It tests the K:0 row at x=1/1000, r=1499/1000, where r+x=3/2 is a P1 knot.

At x<r the shell lift is zero, and the unreflected characteristic formula is

```
B1_x = 1/2 {x F(a) + [PG(a)-PG(q)]/2
             + [PF(a)-PF(q)]/2 + x G(q)},
a=r+x, q=r-x.
```

The primitive of F and G is continuous and their ridge profiles are
continuous; hence the only B1_xx jump when a crosses a single P1 knot while
q remains in one cell is `x Delta F'/2`. The exact free-wave K:0 displacement
is identically zero in `initial/full23_initial.py` (the first mean_u entry);
the exact velocity integral gives `Delta B0_xx = Delta v'/2`. The P1 vector
formula for F gives

```
Delta F' = [(A+B)(Delta du'+Delta v') + D Delta u']/2,
Delta V_xx = (Delta v'+x Delta F')/2.
```

In the 0<q<a region, each P1 u,du,v is linear, V=integral v is quadratic,
F and G are quadratic, and their primitives PF and PG are cubic. Thus B1_x
as a function of x at fixed r is cubic on each side. Exact four-point
one-sided formulas yield its derivative without a numerical limit, as long as
all stencil points remain in the adjacent cells; the companion check
`check_candidate_kink.py` verifies this explicitly. Here the left/right
B1_xx stencils come from the independently integrated segment expressions;
they give exactly `x Delta F'/2`. Both terms add to a nonzero exact rational,
approximately `3.6483656164e20` in registered scaled-coordinate units.
It is a kink in this nominal approximate candidate, not a conclusion about
the exact response or a physical force.

Run from any directory with Python 3 (stdlib only):

```
python audit_b1.py --repo-root /path/to/HDblast
python -O audit_b1.py --repo-root /path/to/HDblast
```

The check verifies the source hash, original profile label and formula,
nonzero total jump, and equality between an exact derivative of independent
segment-integrated expressions and the slope-jump formula. The underlying
half-line Green-kernel identities and numerical accuracy of the candidate
remain inherited premises. This does not recompute the actual source or
establish an exact nonlinear solution.
