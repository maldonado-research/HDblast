# Join the background overcover to the selected exact-family slab

This clarification changes neither frozen background report. Their SHA-256
remains `abd05a8a93fd2f543a112c344ca7f92cf95ce6d56e379cb80ba1e7d9dfb94882`.
The auxiliary background tube up to proper time `1/5000`, bounded by
coordinate `x<=8/5`, is a valid **background overcover**. The selected
exact family is only source-bound on `0<=t<=3ell/2`. Its Gaussian observer
construction is therefore restricted here to

```
0<=tau<=tau_h=1/5500,
ell/16<=zi<=3ell/16, whole S3.
```

This restricted tube contains the target sphere with strict margin and
lies strictly inside the retained exact-family coordinate-time slab. No
extension of that family beyond `t=3ell/2` is used or asserted.

The source is
`research/HDBLAST_EXACT_FAMILY_TT_TRANSFER_20261031/DERIVATION.md`, sections
3–4: its selected correction annulus is `ell/2<z<5ell/2`, and its prescribed
slab is `Tphys=3ell/2`, `Z=6ell`. The source has 35768 bytes and SHA-256
`69bfe846cc4fda8ec6a981969fbdd12fff7862b5ca7cf31a19d696b0a839165a`.
`FAMILY_DOMAIN_SOURCE_PINS.json` pins this source, its domain receipt/gate,
the selected-source pin files, and the original exact-family capsule
manifest. `check_family_domain_join.py` authenticates them by relative
repository paths, checks their manifest entries, and uses exact fractions
only. No repository Python is imported.

The frozen target proper-time upper endpoint obeys the simpler strict
bound

```
tau_dagger_hi < 131*(5ell/4) = 131/800000 < 1/5500,
1/5500 - 131/800000 = 159/8800000 > 0.
```

Using the actual receipt endpoint gives the larger exact headroom

```
1/5500 - tau_dagger_hi
 =138589358236169298911263847684649
  /7139408180485387989229432452677632000 > 0.
```

The auxiliary background bound `|chi|<=15H`, `H=8ell/5`, gives
`|chi|<=3/125000`. For `h^2/2<1`, the exponential series implies
`cosh h<=1/(1-h^2/2)`: each coefficient satisfies
`1/(2n)!<=1/2^n`. Therefore

```
cosh chi <= C_chi = 31250000000/31249999991.
```

Since `dt/dtau=cosh(chi)/rho`, `rho>=129`, and initial `t=0`, the whole
restricted background tube has

```
t <= tau_h*C_chi/129 = 62500000/44343749987229,
x <= 62500000000000/44343749987229 < 10/7 < 3/2.
```

The strict physical-time margin to the selected endpoint is

```
3ell/2 - tau_h*C_chi/129
 =8031249961687/88687499974458000000 > 0.
```

An enclosing box for the selected-family observer tube is thus
`0<=x<=10/7`, `1/32<=r<=1/4`, with coordinate-time margin `1/14` in
scaled units to the old slab endpoint. Its spatial box lies strictly in
`0<r<6`. The first-field certificate's domain is also respected.
The frozen background determinant, Gram, clock, and radial-margin bounds
remain valid on this smaller tube.

Under the retained exact C2 family/regularity theorem, ODE and initial-label
Jacobi dependence and the strict margins give a common sufficiently small
parameter neighborhood on this restricted tube. Initial coordinate time
is fixed at zero; positivity of `dt/dtau` keeps it nonnegative. The strict
upper-time and spatial margins keep the perturbed tube inside the already
supplied exact-family slab. Compact inverse-map and screen margins persist
there. No effective parameter radius is claimed.

The frozen report's `family_persistence` entry is a conditional ODE
stability statement. One of its explicit premises is **"whole family
uniformly defined on an open neighborhood of the compact tube"**. That
premise is not supplied for the entire auxiliary `x<=8/5` overcover by the
selected theorem. This note applies that stability statement only to the
restricted `tau<=1/5500` tube just proved to fit inside the old slab. It must
not be read as granting exact-family existence through `x=8/5`. At the
initial face `t=0`, a neighborhood means a neighborhood relative to the
forward spacetime/one-sided initial-value domain. Standard one-sided ODE
dependence suffices; no negative-time exact-family extension is required.

The initial observer variation is exactly zero on the whole label band.
The source curve is

```
x_e=x0+e*v_seed+A*u(e),
supp(v_seed) subset {ell<=z<=2ell},
supp(A*u(e)) subset closure{ell/2<z<5ell/2}.
```

Hence every parameter's canonical initial data equal the background on
an open neighborhood of all these initial labels. The geometric gap is

```
ell/2 - max(zi) = ell/2 - 3ell/16 = 5ell/16 = 1/3200000 > 0.
```

The source's prescribed initial lapse and shift are
`N_e=rho(1+e*u*Y)`, `beta_e=0`. Since `u` and all its local spatial
jets vanish below the seed support, these equal their background values
throughout that neighborhood. The local relative-harmonic completion
also preserves this exact open zero collar. Thus the future unit initial
normal is exactly `rho^-1 partial_t` for every parameter there, and the
fixed initial label embedding is independent of parameter. It follows
exactly that

```
xi(0)=0, v(0)=0,
partial_zi xi(0)=partial_zi v(0)=0,
partial_Omega xi(0)=partial_Omega v(0)=0.
```

In particular the first-position variation bound has `W0=0`. This fact
uses the selected curve's full correction support and prescribed initial
gauge, rather than merely the support of the first seed.

Normal and optimized receipts `FAMILY_DOMAIN_JOIN_NORMAL.json` and
`FAMILY_DOMAIN_JOIN_OPTIMIZED.json` are byte-identical, with SHA-256
`a2aae6d7bd1e522e1cf176349c675fc87b281bf0e8296b6b22e904fc3b49eeee`.
Reproduce them with:

```sh
python -I -B check_family_domain_join.py --repo-root /path/to/repository
python -I -O -B check_family_domain_join.py --repo-root /path/to/repository
```
