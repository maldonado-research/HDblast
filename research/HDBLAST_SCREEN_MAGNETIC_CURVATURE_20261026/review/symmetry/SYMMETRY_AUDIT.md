# Independent scalar-square magnetic TT selection audit

The registered scalar seed does **not** have a symmetry that forces its
quadratic magnetic screen TT projection to vanish. The frozen seed is
`Y=n0^2-n1^2`, whose trace-free ambient matrix has spectrum `(1,-1,0,0)`.
Its repeated zero eigenvalue gives an O(2) stabilizer; it is not the usual
SO(3)-axisymmetric quadrupole, whose spectrum is proportional to
`(3,-1,-1,-1)`. An equal standing scalar seed in both slots of a symmetric
quadratic map remains capable of generating the selected tensor sector.

Under `Spin(4)=SU(2)_L x SU(2)_R`, scalar degree two is `(1,1)`. Exchanging
the two factors of its tensor square acts on `(J_L,J_R)` by
`(-1)^(J_L+J_R)`. Consequently

```
Sym^2(1,1)=(0,0)+(1,1)+(2,2)+(2,0)+(0,2),
dim=1+9+25+5+5=45.
```

The last two summands are the two degree-two TT representations. This
argument incorporates the symmetry in the two inputs; it does not rely
on treating independent angular quadrupoles as available.

The ambient reflection `R=diag(1,1,1,-1)` fixes the seed. With covector
slots transformed as well as polynomial arguments, it exchanges the
frozen `Tplus,Tminus` templates. A polar metric response therefore obeys
equal selected TT coefficients. A magnetic tensor defined with one
screen epsilon is an axial tensor: the same reflection instead allows
opposite coefficients. Equality of the metric coefficients does not imply
that their symmetric curl is zero. Reversing the fixed orientation changes
all magnetic signs but not this conclusion.

There is a direct exact witness. On the unit sphere, use
`epsilon_ikl=sum_m epsilon_mikl n_m`, `epsilon_0123=+1`, and normalized
symmetrization. Define

```
S_ab=(D_aY)(D_bY),
M_ab=epsilon_(a^{cd}D_c S_{b)d}
    =epsilon_(a^{cd}(D_cD_{b)}Y)(D_dY).
```

In the actual frozen template normalization the all-point identities are

```
Pi_TT S=-(2/3)(Tplus+Tminus),
M=-2 Tplus+2 Tminus,
tr M=0,   div M=0,
<Tplus,Tplus>=<Tminus,Tminus>=1/6.
```

The second equality holds as the complete tensor, with no remaining
scalar or vector component. It proves that quadratic angular scalar jets
can generate magnetic TT. It does not compute the geometric coefficient
of this angular tensor in the physical `Q_H[q,q]`; cancellations among
that coefficient's metric, observer, screen and volume terms remain a
separate calculation.

For comparison, the same witness is identically zero for
`Y_axis=3n0^2-n1^2-n2^2-n3^2`. In that case the gradient supplies one
distinguished tangent direction and the Hessian is a combination of its
square and the sphere metric; the epsilon symmetrization vanishes. More
generally, any local natural axial symmetric two-tensor constructed from
this one axisymmetric scalar and invariant background data vanishes:
the O(2) isotropy on the transverse tangent plane admits only diagonal
SO(2)-invariant symmetric tensors, while its orientation reflection
requires an axial tensor to change sign. This stronger symmetry is absent
for the registered seed.

The linear second-metric isolation is also sound at the representation
level. Scalar sections have `(j,j)` labels, coexact-vector sections have
`|j_L-j_R|=1`, and TT sections have `|j_L-j_R|=2`. Every background-natural
linear operator assembled from the sphere metric, derivatives and epsilon
is an SO(4) intertwiner. The epsilon changes reflection parity, without
changing these SO(4) labels. Thus the TT projection of the linearized
screen magnetic operator cannot receive any scalar or vector second
metric amplitude. This argument includes all section multiplicities and
does not assume the two TT template labels are curl eigenlabels.

`check_scalar_magnetic_symmetry.py` uses its own standard-library rational
polynomial calculus, imports no frozen source, and checks the displayed
identities in the sphere polynomial quotient, including exact trace,
divergence, template norm and reflection-slot transformations. The
normal and optimized reports are identical. Six separate intended-error
subprocesses reject an incorrect magnetic sign, an incorrect electric
normalization, and omitted reflection signs on covector slots. See
`SYMMETRY_NORMAL.json`, `SYMMETRY_OPTIMIZED.json`, and
`SYMMETRY_CONTROLS.json`.

The checker validates its output before the scientific calculation. It
rejects existing destinations, symlink leaves and ancestors, writes inside
its own source directory, writes under any
`research/HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026` capsule, and missing
parents. It creates the report exclusively. `run_symmetry_controls.py`
reproduces both scientific reports, six intended-error rejections, and
twelve early output-guard rejections in a fresh external output directory.
The guard probes include a scientific mutation and require the output
guard to reject first. Source hashes and sentinel contents are checked
after the run. The scientific goldens and receipts contain no host paths.

The frozen definitions were read in
`research/HDBLAST_NUMERICAL_RESPONSE_INPUTS_20261016/angular/angular_exact.py`
and the frozen selected reflection and scalar/vector/TT block argument
were read in
`research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026/structure/exact_angular_structure.py`.
No frozen file, source calculation, response campaign, or remote state
was modified.
