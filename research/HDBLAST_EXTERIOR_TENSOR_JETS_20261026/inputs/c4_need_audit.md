# Exact derivative requirement of the complete TT source

The entire authenticated map `metric_forcing/ab/5` has 78 terms and is exactly
equal to row 6. Its only second physical derivative symbols are `Myy`, `Myt`,
and `Mtt`. This statement is obtained by scanning every power in every term,
with no source term omitted. The source map SHA-256 is
`34b10b62d0b36a4506a35407a591a9fdff319971994f25e0b5f031f48eb4f615`.

In the registered reconstruction, `M=2d` contains no derivative of the actual
gauge amplitude `d`. The other frame amplitudes `A,B,C,D,E,L` contain at most
one derivative of the master or gauge fields. They occur in this TT row only
undifferentiated or with a single physical derivative `y` or `t`. The scalar
symbol `s` occurs undifferentiated and contains at most one master derivative.

Consequently a single source-coordinate derivative `partial_r` needs at most
third derivatives of the actual first master and gauge solutions. For the
second derivative symbols of M this is explicit:

* `(Myy)_r` uses `d_rrr` together with differentiated background factors.
* `(Myt)_r` uses `d_xrr` together with differentiated background factors.
* `(Mtt)_r` uses `d_xxr`.

All other metric symbols need at most two derivatives before `partial_r` and
at most three afterwards. One must differentiate the physical conversions
`partial_y=-(ell*rho)^(-1) partial_r` and `partial_t=ell^(-1) partial_x`,
including every occurrence of rho. The coefficient symbols `rho,a,ay,Vp` and
the time factors must also be differentiated, and the constant ell has zero
radial derivative. Background ODE derivatives through order three suffice;
the authenticated background routine supplies order four.

There is therefore no need to assert a new C4 actual-solution error bound.
The authenticated actual C3 errors are precisely sufficient for this selected
radial source derivative. Constructing the derivatives of the actual smooth
functions uses those C3 error enclosures coefficient by coefficient, rather
than differentiating independently chosen interval values. Any evaluator
which forms the full physical third jet indiscriminately must avoid silently
using its unavailable fourth master/gauge coefficients: only the symbols
listed by this exact dependency audit are authorized.

The premise remains the retained smooth actual background/master/common-gauge
solutions and their authenticated integral/PDE error certificates. This audit
does not independently prove existence and does not itself enclose the source
derivative. `c4_need_audit.py` authenticates the map bytes before parsing them,
imports no frozen Python, checks all 78 terms, and emits the dependency table.
