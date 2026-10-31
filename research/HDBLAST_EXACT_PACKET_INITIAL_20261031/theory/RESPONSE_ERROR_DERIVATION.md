# Downstream audit for the frozen exact-packet displacement experiment

31 October 2026. Read-only review of the source tree
`/workspace/hdblast-reduced-publish-31a`; all new artifacts are confined to
`/workspace/hdblast-response-gate-audit-31b`. No Git or source mutation was
performed. This note audits transfer assumptions and arithmetic structure;
the producer and independent evaluator remain responsible for proving the
new full-cover initial and local interval gates.

## Decision and stopping boundary

No fatal regularity or transfer incompatibility was found in the specified
candidate: complete the registered nominal displacement, pack all23 slots,
apply the five mean overrides, then extend each completed component by its
own exact degree-two Taylor polynomial at R=7/4; omit displacement
convolution; retain the byte-defined old smoothed P1 velocity. This provides
the displacement C2 and velocity C1 required by the original physical-side
C2 response proof. It supplies no higher global displacement regularity.

The benchmark must stop before response evolution unless the independent
local and complete initial gates pass. Failure of an enclosure to prove a
threshold is inconclusive unless a certified lower bound violates it. The
audit does not replace these gates with the theorem's conditional validity.

The factor-two thresholds are experimental acceptance criteria frozen in
`NEXT_QUANTITATIVE_TARGET.md`, not factors produced by the wave theorem.
Use exact retained rational halves in `retained_constants.json`. Their
source reports are byte-identical to the 30k normal replay references.
Readable displays are in `RETAINED_CONSTANTS_NOTE.md`.

## Initial regularity and trace inventory

The raw functions are q=exp(4-1/((r-1)(2-r))) for r>1 and zero for r<=1,
and E=exp(-ell r), ell=10^-6. Nonzero nominal rho and p on the complete
footprint are essential, including the edge branch. With those denominators
separated from zero, exact q is smooth and flat at 1. Products, inverses of
the regular coefficient functions, and finitely many derivatives preserve
flatness. This includes all five mean overrides, scalar slots and structural
zeros. Every required displacement jet vanishes at 1. The two old artificial
cutoff radii are evaluation boundaries of the same analytic expression;
their exact left/right jet differences are zero, not merely overlapping
intervals.

The registered construction differentiates the packet before reconstruction.
The highest raw packet orders in completed displacement values are: zero
for N:0,2,4; two for J:0,2,4; three for w:0,2,4; K, L, vector, H and TT
displacements are structural zero. Consequently q through order five is
sufficient to compute the completed displacement through rr. E enters
f_s algebraically and is differentiated once in f_y, hence the scalar
displacement through rr needs E through radial order three. These counts
do not license a C2 raw-packet extension.

For each completed component f_i and h>=0 define

```
f_i,ext(R+h)=f_i(R-)+f_i,r(R-) h+f_i,rr(R-) h^2/2.
```

The interior and exterior traces at orders 0,1,2 are exactly the same
analytic quantities; their differences are identically zero. Structural
zero components stay zero. The exterior third derivative is zero and the
interior third derivative generally is not. Thus the global auxiliary
displacement is C2 and generally not C3. The outer extension is a definition
of the completed datum, never a packet substituted back into reconstruction.
An unbounded quadratic exterior is harmless here: free-wave initial
arguments on D+ lie in [-R,R], and the new displacement is not convolved.
Only the retained compact footprint is used in the value and slope caps.

The unchanged velocity is the complete stored rational all23 P1 vector,
extended at R with its final P1 slope and convolved with

```
k_delta(s)=35/(32 delta)*(1-(s/delta)^2)^3, |s|<=delta,
delta=1/2048.
```

It is not newly reconstructed from the exact packet. Its first two stored
node vectors are zero and its first nonzero segment begins at 65/64, so
the smoothed velocity is identically zero through
65/64-delta=2079/2048>1. Its reflected origin collar is flat. The kernel
is C2 and the P1 distributional second derivative is a finite sum of slope
jumps. Thus v_s'' is a finite sum of C2 translates and v_s is C4. Around R,
the final P1 cell and the exterior use the same affine function, so the
smoothed velocity equals that affine function: its 0/1/2 traces are the
last node, last secant, and zero, exactly on both sides. The last cell is
wider than twice the kernel radius. This also certifies C2 velocity if a
presentation states that all initial data are C2. The full per-component
velocity traces and unchanged error charges were independently inventoried
in the producer's `inventory-velocity.json`; its hash is recorded in the
source-pins artifact.

Both displacement and velocity therefore vanish for r<=1. The free wave,
local coefficient multiplication, physical reflected K, causal L, and
Volterra inverse P preserve the wedge x+r<=1. This re-establishes the
whole-domain value factor 3/4 and shell factor 1/4, with admitted shell
history ending at x=5/4 and active shell length 1/4. Displacement convolution
would invalidate this argument by spreading the exact support edge.

## What the original C2 theorem actually needs

The original-source theorem uses initial displacement through rr and
velocity through r. A time derivative of the prescribed forcing, rather
than a radial derivative of the initial curvature, produces the candidate
second jets. Hence no global completed-displacement third derivative or
new actual sixth derivative is needed for this displacement-only comparison.
No old nominal M3 convolution allowance is attached to the new displacement.
The old complete velocity nodal-radius, interpolation, slope-jump and
convolution charges remain.

The actual coefficients/source have the retained continuous physical-side
regularity, actual F_x, G_x, and compatible flat initial/shell corner. The
new initial C2/C1 data preserve that corner. Differentiation gives a
physical-side C1 problem for the time derivative; the principal equation
recovers rr continuously. The original review's differentiated Picard
iterates have a summable n q^n estimate, so this is also a regularity
argument, not an estimate on an assumed derivative. One-sided shell and
domain-edge traces are included. C2 across the artificial doubled shell
is neither assumed nor claimed: an odd reflected output source may jump
there. The fixed target's regularity and analytic physical/constraint
premises remain inherited; this auxiliary initial datum is not thereby
proved to satisfy the physical constraints exactly.

## Complete Cauchy and boundary cancellations

Use actual La=A partial_x+B partial_r+D, frozen constant L0, and
P=(I-L M PN)^-1. Define the new fields, with the original actual F,G,M,

```
f0=L0 Bs,
W0=P[Bs+K f0+LG],
W1=P[Bs+K(La W0+F)+LG],
Dfield=W0-Bs,  C=W1-W0,
r0=La Dfield+(La-L0)Bs+F,
C=P K r0,  r1=La C.
```

K and the compatible boundary correction have zero Cauchy displacement and
velocity. Hence W0, W1 and Bs share u_new, v_old and their tangential
radial first jet. Dfield and C have zero initial value, x derivative and
r derivative. They need not have zero initial acceleration:

```
Dfield_xx(0)=f0(0),
C_xx(0)=r0(0),
r1(0)=A C_x(0)+B C_r(0)+D C(0)=0.
```

This common first Cauchy jet proves the required exact r1 cancellation.
It does not set F(0) or either correction acceleration to zero.
Equivalently, partial_x Kf=K f_x+Bfree[0,f(0)]. The free initial-velocity
term is mandatory in both differentiated corrections.

Subtract the free-wave shell law before differentiating Dfield. Its
boundary source is M PN Bs+G, and Dfield_x therefore has source
M PN Bs_x+G_x. By contrast, C uses the same actual affine G and M in both
terms; its boundary source and differentiated boundary source are zero.
The target-error e=U-W1 likewise has zero boundary source and derivative.
Actual M_x=0. G_x is charged in Dfield_x, never silently omitted or
incorrectly reintroduced as a final boundary-error charge.

## Minimal valid full-response recalculation

All constants below refer to the same all23 positive H norm and D+.
Let T=3/2, a,b,d,m bound actual A,B,D,M, a0,b0,d0 the frozen operators,
delta_a,delta_b,delta_d their differences, and

```
c=a+b+3d/4,  c0=a0+b0+3d0/4,
comm=cap(A_x)+cap(B_x)+3cap(D_x)/4,
bdry=m/4,  Pinv=1/(1-bdry),  q=3409/8000000,
betaR=(q-bdry)/(1-bdry).
```

Authenticate/recheck the actual matrix/Gram/operator constants; they are
unchanged by the candidate. Supply fresh complete absolute candidate caps
U,U1,U2,V,V1 and distinct paired initial error caps eu,eur,eurr,ev,evr.
Use errors only in the target-error equation, and absolute caps in the
candidate forcing ledger. Recompute

```
Bs=max((U+T V)/3,U1+V),
I=max((eu+T ev)/3,eur+ev),
Qfree=U2+V1.
```

One may recompute the old sharp shell-history construction. A simpler
fully valid fresh bound, avoiding reuse of changed shell samples, is

```
D0X=Pinv*((T c0+bdry)*Bs+Gcap),
R=T*((delta_a+delta_b+3delta_d/4)*Bs+c*D0X)+KFcap,
CX=Pinv*R,
E=(I+betaR*R)/(1-q),  F0=c*CX.
```

Here D0X bounds X(W0-Bs), not X(W0). Gcap=487122 and
KFcap=52162547172519/2+Gcap are conservative retained source-response
caps. The latter is an X(KF) cap, not a pointwise F magnitude. The direct
global D0X formula follows from
Dfield=P[K L0 Bs+L(M PN Bs+G)]; it preserves the actual Robin inverse.
It is slightly weaker than the old specialized shell history but is
legitimate for new data. A sharper independently re-proved D0X or R may
replace it; old shell point values cannot simply be inherited.

For the differentiated correction ledger use the original pointwise
source caps FH=17387513156573, Fx=905405545954284 and Gx=58807538:

```
f00=a0 V+b0 U1+d0 U,
f0x=(a0+b0)*Qfree+d0*Bs,
D0xX=Pinv*(f00+T*f0x+m*Bs+Gx),
r00=delta_a V+delta_b U1+delta_d U+FH,
r0x=(a+b)*D0xX+(d+comm)*D0X
     +(delta_a+delta_b)*Qfree+(delta_d+comm)*Bs+Fx,
CxX=Pinv*(r00+T*r0x),
F1=(a+b)*CxX+(d+comm)*CX.
```

These are the original proof's cancellation-preserving expressions. They
avoid paying separately for large absolute W0 and Bs second jets. Retain
the initial terms f00 and r00, actual Gx and m*Bs in D0xX, coefficient
commutator, both inverse types, and the unchanged velocity caps.

Finally, because r1(0)=0 and the final boundary error is zero,

```
w0=eurr+a*ev+b*eur+d*eu,
J=max((ev+T*w0)/3,evr+w0),
Z=(J+T*(F1+comm*E))/(1-q),
Bxx=Bxr=Z,  Brr=Z+c*E+F0.
```

Compare E to the old exact C1 cap and each new Bab to its own exact old
half. A single old common readable bound is not the threshold. Divide by
the same global actual-target lower scale only after reporting absolutes.
The scale cannot replace the local denominator, normalize xx/xr separately,
or certify pointwise relative accuracy.

No causal-source-smoothed Z_epsilon bound is used. In particular its
filtered_F_x_H=66768050521240320 and initial forcing F(0) belong to another
candidate. The shared old absolute initial jets may authenticate old
data, but must be replaced by new absolute displacement jets here.

The standalone `response_ledger_oracle.py` implements these equations using
exact fractions and source pins, and independently replays the retained
exact C2 ledger before permitting new arithmetic. It is an arithmetic
oracle conditional on certified supplied caps and premises; it is not an
independent interval evaluator or gate authorization.

## Affected outputs and scope

Removing displacement convolution changes the initial function wherever
the old kernel averaged it, including regions where unsmoothed formulas
agree. Rebuild the finite response definitions and all affected point,
quadrature and constraint tables, including (9/8,0), after initial gates
pass. The old outer-collar causal equality argument does not protect these
points against this new displacement change. A source-adapter identity
or actual q3 bound is not a new response approximation.

Frame/natural-equation or observable claims require their own conversion,
connection, subsidiary, history, clock, quadrature and nonzero-reference
ledgers. A successful experimental certificate comparison proves neither
a factor-two reduction of unknown global true error nor a nonlinear or
cosmological conclusion. The local upper-versus-old-lower gate, if proved,
does establish the stated actual local improvement for this candidate.

## Executed response arithmetic for the passing initial campaign

The final stable producer report is
`/workspace/hdblast-exact-packet-31b/normal/INITIAL.json`, SHA256
`692ac13f72354c526a0108ddae32800a9f22efb1a8afe4e987cbf403b5a296c3`.
It reports a passing 128-bit, 49-leaf full-cover initial campaign. Its
paired velocity errors equal the retained rational values exactly. The
response oracle independently rechecks the initial acceptance inequalities
and all23 local positive-Gram aggregation before any response arithmetic;
the independent interval evaluator checks the expression/coverage itself.

Using precisely the fresh global Dfield route frozen above, the complete
new conditional response certificate passes every exact benchmark gate:

| Quantity | New readable strict upper | Retained acceptance maximum (display) |
|---|---:|---:|
| Complete C1 X error | 1.603e23 | 3.54033516534e25 |
| xx error in H | 2.135e25 | 3.32099897079e26 |
| xr error in H | 2.135e25 | 3.32099897079e26 |
| rr error in H | 2.144e25 | 3.32121561210e26 |

Every decision is an exact rational comparison; the decimal table is only
a display. The exact rr upper is approximately 2.143224747842756e25, and
xx=xr approximately 2.134073127028094e25. These are large absolute
scaled-coordinate bounds. The same global actual-target scale gives
rr/max-second-seminorm quotient below 6.597e-9; it retains the same limited
global meaning as the previous quotient. The old/new rr upper-certificate
ratio exceeds 30.99, and the C1 upper-certificate ratio exceeds 220.97.
Neither ratio compares the unknown global true errors.

The complete initial local H error upper at P* is approximately
3.032832000704e-13, divided by the certified actual target lower approximately
0.6942780997055390. Its ratio is below 4.369e-13 and the error is strictly
less than half the old complete-H true-error lower. This is a local
upper-versus-lower comparison, distinct from the global certificate ratios.
The exact target and error enclosures and all23 squared contributions are
retained in the producer input and independently aggregated in RESPONSE.

The new full second-error majorant is no longer mostly the initial
displacement rr error. Its unchanged velocity-gradient input is about
9.67170e24; the differentiated correction-forcing cap F1 is about
7.68077e24 and receives the factor T=3/2 in the final numerator. The new
initial displacement rr-error is about 1.38790e23. These are proof-majorant
channels, not a decomposition of unknown actual error. No velocity change
or additional bulk iteration was made.

The oracle first reproduces the exact retained C2 ledger at every reported
node. Its new `normal/RESPONSE.json` and `optimized/RESPONSE.json` are
byte-identical at SHA256
`991ee9e04d17a879e4e65f3f68ad8f9e76929b4c10d4698d2af06dd6588f3273`.
The executed oracle is SHA256
`c74b48b392d9f3741edc1c80e19aaa1eed3578775d93520c35c7e7fa1126f211`.
Canonical input pins use stable source-relative paths and roles; physical
runtime paths are checked internally and do not change report bytes.
An additional real replay relocates the executed oracle, candidate input,
and every scientific report to different paths and a sparse source root;
the complete normal and optimized report bytes remain identical. See
`portability/PORTABILITY.json`.
Eleven deliberately altered formula executions in each Python mode all
fail with an exact nonzero arithmetic witness and no success artifact;
see `controls/CONTROLS.json`. Some variants add a conservative charge or
lose a cancellation rather than falsifying an upper bound; those are
rejected as departures from this frozen canonical recipe. Omission variants
demonstrate the nonzero missing charges. Producer/independent interval
formula, coverage and source-role controls remain separate.

No response point integration, quadrature table or constraint/natural
residual was produced here. The output is a conditional response-bound
corollary for `HDBLAST_EXACT_PACKET_INITIAL_20261031` using the specified
new initial candidate and original actual target/source.
