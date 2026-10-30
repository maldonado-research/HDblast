# Inclusion under local refinement

Fix the target, nominal initial functions and their three telescoping
difference channels `D_c,k,i(r)`. Suppose an old parent interval `P` has a valid
component cap `B_P,c,k,i`, and an exact finite cover `I_j` of `P` has newly
computed caps `B_j,c,k,i`. At any point of a leaf both bounds apply, so

```
|D_c,k,i(r)| <= min(B_P,c,k,i,B_j,c,k,i),  r in I_j.
```

Together with all untouched parent intervals, this is a continuous cover of
the old domain. For each component, taking the maximum over the complete
cover gives a valid global magnitude cap `C_c,k,i`. The positive diagonal
Gram norm then gives

```
||D_c,k(r)||_H <= sqrt(sum_i H_i C_c,k,i^2),
||sum_c D_c,k(r)||_H <= sum_c sqrt(sum_i H_i C_c,k,i^2).
```

Each square root is rounded upward with exact rational arithmetic. The
different channel maxima and norm radii need not be attained at one common
point; the triangle bound remains valid. Alternative cellwise norm bounds
could be valid too, but the present calculation keeps the original
max-before-Gram convention to make the change directly auditable.

The adaptive selection rule affects computation time, not this inclusion
argument. The published calculation stops at a finite cover and verifies its
rational endpoints, positive lengths, lack of gaps/overlaps and analytic
branch boundaries. The fixed branch endpoints are included. The inherited
cutoff has matching derivatives through order four there, covering the
initial quantities used by the paired evaluator. Its regularity and paired
algebra assumptions remain inherited premises.

For context, the registered cutoff

```
S(s)=126s^5-420s^6+540s^7-315s^8+70s^9
```

has the exact identities

```
S'(s)=630s^4(1-s)^4,
S(s)=s^5[1+5(1-s)+15(1-s)^2+35(1-s)^3+70(1-s)^4],
1-S(s)=(1-s)^5[1+5s+15s^2+35s^3+70s^4].
```

Expansion and differentiation verify these polynomial identities. They expose
cancellation in the alternating power representation near `s=1` and explain
why local interval subdivision can help. The new numerical certificate uses
the **unchanged original evaluator**, including its original cutoff expression
and exponential error estimates. Factored cutoff evaluations and sharper
signed kernel inequalities were considered separately and are not credited
with the published improvement. This is standard interval inclusion and
polynomial algebra, not a new mathematical method or an additional physical
assumption.
