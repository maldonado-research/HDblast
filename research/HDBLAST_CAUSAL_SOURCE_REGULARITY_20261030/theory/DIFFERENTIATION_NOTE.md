# Independent differentiation audit for prescribed wave forcing

For the full-line principal wave propagator

```
Kf(x,r)=(1/2) integral_0^x integral_(r-x+s)^(r+x-s) f(s,y) dy ds,
```

time differentiation satisfies

```
partial_x Kf = K(partial_x f) + Bfree[0,f(0)].
```

The second term is a free wave with zero displacement and velocity f(0,r).
It is essential when the initial forcing is nonzero. This is an identity
before applying any physical output reflection, coupling, or boundary
reconstruction; it does not require commuting those operations.

`check_time_differentiation.py` expands the original inner antiderivative at
its moving radial endpoints and integrates every resulting s power exactly.
It checks equality of the resulting univariate x polynomials for each
monomial x^a r^b with 0<=a,b<=3 at three rational radii, with separate
principal wave and zero Cauchy checks at three rational events. A mixed
polynomial check combines initial-nonzero and initial-zero forcing terms.
The case f=1 explicitly rejects omission of the free initial velocity term.

For a free wave Y with displacement v and velocity w and an H-isometric
reflection where applicable, suppose sup||v||<=V, sup||v_r||<=V_r and
sup||w||<=W. The d'Alembert formula and triangle inequality give, on a
retained time interval of length T=3/2,

```
sup||Y|| <= V+T W,
sup||Y_x||, sup||Y_r|| <= V_r+W,
X(Y) <= max((V+(3/2)W)/3,V_r+W).
```

A constant w with v=0 witnesses failure if W is omitted from the derivative
branch. On the retained initial interval [-T,T], the smooth displacement
v(r)=V_r*r^2/(2T), together with constant w=W, attains Y_x(T,0)=V_r+W.
A constant w also witnesses failure of using a unit-time value increment
on the length-3/2 interval. The report includes these manufactured checks;
it does not turn the triangle bound into a claim about a particular actual
HDBLAST field.

Normal and optimized reports match byte for byte. No scientific imports,
floating point decisions, or assert statements are used. This audit does
not address physical-side regularity across an output parity interface,
which requires the separate boundary-source argument.
