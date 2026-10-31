# Only two position variations are needed for the coefficient

Along the background radial geodesic let `c=cosh(chi)`, `s=sinh(chi)`,
`V^t=c/rho`, `V^y=s`. Write its first coordinate position and velocity
variations as `xi^t Y,xi^y Y,delta V^t Y,delta V^y Y` in the t,y slots.
The angular trajectory variation exists but does not couple into these two
equations at first order: the background mixed base/angular Christoffels
with two radial velocity slots vanish. The corresponding first connection
forcing in the t,y equations uses only A,B,D and their first t,y derivatives.

Define `X=rho xi^t`, `Z=xi^y`, `P=rho delta V^t`, `Q=delta V^y`; the
symbol Z here is a variation coefficient, not the inward coordinate z.
With primes denoting actual proper-time derivatives,

```
X' = a s X+P,                     Z'=Q,
P' = -a s P-2a c Q-2a_y c s Z-S_P,
Q' = -2a c P-(2a²+a_y)c² Z-S_Q,

S_P = a B-(c² D_t-s² A_t)/(2rho)-c s D_y-s² B_y,
S_Q = c²[B_t/rho-a(A+D)-D_y/2]
      +c s[A_t/rho-2aB]+s² A_y/2.
```

Here `A_t,A_y` denote derivatives of the metric amplitude A, not the Weyl
invariants named `A_t,A_y` in the main formula. The API avoids this collision:
metric derivatives are `m['At'],m['Ay']`, while curvature combinations are
returned separately by `first_weyl_invariants`.

To verify signs directly, the nonzero base connection coefficients of the
background are `Gamma^t_ty=a` and `Gamma^y_tt=rho²a`. Their first variations
per unit harmonic Y are

```
delta Gamma^t_tt=-D_t/2+a rho B,
delta Gamma^t_ty=-D_y/2,
delta Gamma^t_yy=-(B_y+aB)/rho+A_t/(2rho²),
delta Gamma^y_tt=rho B_t-rho²[D_y/2+a(A+D)],
delta Gamma^y_ty=A_t/2-rho aB,
delta Gamma^y_yy=A_y/2.
```

Contracting them with `(c/rho,s)` gives S_P after multiplication by rho
and S_Q directly. The identity `c²-s²=1` produces the isolated `aB` term
in S_P. Differentiating the background connection and velocity gives the
drift terms above; differentiating the factor rho contributes `a s P` and
`a s X` where appropriate.

For fixed initial chart labels, xi starts at zero. In general the prescribed
initial unit normal gives a nonzero initial velocity variation. On the
specified launch tube near scaled radius 1/8, the retained first initial
seed is identically zero because its support is `[1,2]`; the selected exact
curve also has its correction annulus away from this tube. Thus all four
initial coefficients above are zero for the actual selected family.

An inexpensive uniform bound follows without solving the angular trajectory.
On a proper-time interval `[0,T]`, put

```
S=max(sup |S_P|,sup |S_Q|),
M=max(sup(|a s|+2|a c|),sup 2|a c|),
K=max(sup 2|a_y c s|,sup |(2a²+a_y)c²|),
b=sup |a s|.
```

For `W=max(sup|P|,sup|Q|)`, integration gives
`sup|Z|<=T W` and `W<=T S+M T W+K T² W/2`. Consequently, whenever the
displayed denominators are positive,

```
W <= T S/(1-M T-K T²/2),
sup|Z| <= T W,
sup|X| <= T W/(1-bT).
```

These are outward bounds for the actual first variation when the source,
background, and trajectory tube are enclosed. Using whole-box actual jet
ranges for S is valid; comparator jets alone would not be. T must bound
proper elapsed time, not the scaled coordinate endpoint 5/4. This lemma
does not establish the Gaussian congruence Jacobian: its separate background
certificate and a common small-family inverse-map margin are still required.
