# Independent bounded feasibility audit for width 10^-6

15 October 2026 UTC. Planning calculation only. No frozen scientific output
was edited and no historical or reference-only producer was executed.

**The same retained shell box and outer rectangle suffice to extend the
analytic collar from 8e-10 to 8e-6.** The uniform waveform majorant also
remains small enough to make a wider reconstruction a reasonable next
calculation. This is a bounded feasibility deduction for the specified next
width, not a completed wider waveform or a new publication checkpoint.

Set `ell=1/1000000`, `T=3ell`, and `L=8ell=1/125000`. The actual background
is the same retained near-phi=+1 analytic branch at exact `delta=1/1000` and
`c=96714321425692/161839258930731`. The exact shell intervals are taken
from the authenticated prior JSON; historical floating trajectories are not
inputs. The new checker independently evaluates the ODE and coefficients.

## Collar first-exit calculation

Use the same closed outer rectangle

| Variable | Outer interval |
| --- | --- |
| rho | [129,131] |
| eta | [-9/100000,-8/100000] |
| p | [-14/100000,-12/100000] |
| a | [111/1000,112/1000] |

Inward depth satisfies

\[
\rho_z=-\rho^2a,\quad\eta_z=-\rho p,\quad
p_z=-\rho(U_\eta-4ap),\quad a_z=\rho^{-1}+\rho p^2/3.
\]

The independently constructed potential polynomial is

\[
U=-2/27+(14/9)\eta^2+(50/27)\eta^3-\eta^4/6
 -(4/9)\eta^5-(2/27)\eta^6.
\]

It follows directly from `U=W_eta^2/2-2W^2/3`,
`W=1/3+eta^2+eta^3/3`. Horner interval evaluation of its exact derivative
gives the following displays. Decisions use exact rational endpoints.

| Component | RHS interval, approximately | Lower first-exit margin | Upper first-exit margin |
| --- | --- | --- | --- |
| rho | [-1922.032,-1847.151] | 0.909105510698679 | 1.07496116366244 |
| eta | [0.01548,0.01834] | 2.76721226168564e-6 | 6.94067738314359e-7 |
| p | [0.0240106266109535,0.0296950799497111] | 2.87649840336551e-6 | 4.10086663074235e-6 |
| a | [0.00763420698625954,0.00775279385116279] | 0.000377370423640991 | 0.000622566559246576 |

For every component, the complete prefix enclosure
`I_shell+[0,L]F(I_outer)` lies strictly inside the outer rectangle. If the
actual background first exited before depth `L`, integrating its bounded
derivative to that exit would instead put it strictly inside the rectangle.
The vector field is smooth near this compact positive-rho rectangle, so
ODE continuation covers the full depth. Moreover `p<0` throughout the
rectangle, keeping every coefficient division by p safe.

The eta upper margin is the smallest relevant face margin here; its exact
value is

\[
234015434201092441/337165122772356250000000>0.
\]

For this fixed rectangle and this one-step interval argument, the supremum
of strict certified depths is

\[
L_* = 5669686020285051/123672167032900272500
 \simeq4.58444786430948\,10^{-5}.
\]

The eta upper face is limiting. Every `0<=L<L_*` passes the strict
first-exit enclosure; equality does not pass strict inclusion. If requiring
`L=8ell`, the corresponding sufficient threshold is
`ell<L_*/8`, approximately `5.73055983038685e-6`. These thresholds concern
this fixed rectangle and interval proof. They do not locate an actual
background exit, singularity or maximal existence depth. The explicitly
tested overextension `L=5e-5` fails this rectangle's eta upper inclusion.

## Coefficient bounds on the wider collar

Using the exact derivative-free identity with `h=U_eta/p`,

\[
\mathcal W=-\frac{35}{4}\tanh^2t+
\rho^2[2h^2-10ah+35a^2/4-U_{\eta\eta}-p^2/6],
\]

the new interval calculation gives, for every real t throughout `0<=z<=8e-6`,
approximately

\[
11699.5678886766\le\mathcal W\le101462.879778632.
\]

The displayed values are not rounded decision endpoints. The unchanged
exact shell box also supplies `mu>0.00205977121593046`,
`202.086050569103<=gamma<=247.367804503084`, and
`202.071852682969<=beta(t)<=247.373483657537` for all real t.
Thus the old coarse bounds remain valid on the extended collar:

\[
m=1/1000,\ M=10^6,\ B=1000,\ \Gamma=250.
\]

The independently recomputed intervals additionally justify the convenient
stronger rounded bounds

\[
m=1/500,\ M=102000,\ B=250,\ \Gamma=250.                       \tag{1}
\]

The small difference from the archived p_z and potential intervals is
expected: this audit evaluates the exactly expanded potential polynomial
by Horner's rule, while the historical certificate used a natural
derivative-product form. Both interval evaluations are outward; the new
one is slightly tighter.

## Same waveform majorant at the proposed width

For the width-rescaled initial bump, the coupled linear data would again
be `w(0,z)=Q(z/ell)` and `w_t(0,z)=exp(-z)Q(z/ell)`. The independent
characteristic/kernel theorem applies to a compatible coupled linear
solution with those data on the new collar. Define

\[
\widehat H(x)=Q'(x)+\ell e^{-\ell x}Q(x),\quad H(x)=\ell D(\ell x),
\quad D=w_z-\gamma w_b.
\]

For `T=3ell`, put

\[
r=MT^2/2,\quad s=BT^2/(2m),\quad
d=(1-r)(1-s)-Tr/m,
\]
\[
A_*=(1+T)(1-s+T/m)/d,\quad V_*=(1+T)/(md),
\]
\[
\eta=\ell[MT A_*+(1+\Gamma T)V_*],
\]
\[
V_{\rm tail}=\ell/m^2+\ell/m+(r/m)A_*+sV_*,\quad
\eta_{\rm tail}=\ell[MT A_*+V_{\rm tail}+\Gamma TV_*].
\]

Both choices of coefficient bounds have positive diagonal and determinant.
The same theorem then gives a global error `|H-Hhat|<=eta` on `[0,3]`,
exact prefront zero on `[0,1]`, and `|H|<=eta_tail` on `[2,3]`.
The support-aware squared norm error is bounded by
`eta^2+eta_tail^2`. Exact evaluation gives:

| Bounds used | Determinant | Global error eta, approximately | Tail bound, approximately | Strict L2 error upper bound |
| --- | --- | --- | --- | --- |
| Retained coarse | 0.99999098652025 | 0.00100377104524695 | 4.76904524694947e-6 | 0.0010038 |
| Stronger rounded (1) | 0.9999989778117581875 | 0.000500683472665010 | 9.32472665010112e-7 | 0.000500685 |

For example, the retained-coarse exact global error is
`2007523995553919/1999981973040500000`. The new audit checks the coupled
majorant equations directly; omitting the determinant feedback term fails
the supersolution check. These estimates are error tolerances for a
proposed wider reconstruction, not sampled waveform values. A post-pulse
bound does not imply that the true tail vanishes.

The wider collar and small linear majorants resolve two concrete feasibility
questions. A full next reconstruction must still instantiate the requested
wider state and verify its retained localization/evolution conditions,
produce any desired continuous-cell waveform export, and recompute its
own response integral. **No C interval from the narrow state is transferred.**
No effective nonlinear epsilon, nonlinear remainder or observational time
scale follows from this audit. The wider first-exit argument is conditional
on the same actual-background theorem; finite arithmetic does not formally
verify that theorem or the analytic PDE/kernel proof.

## Reproduction and independent comparison

Planning scripts and outputs may reside outside the unchanged science
capsule. From any working directory, point `--source-root` to a separately
extracted `HDBLAST_UNIFORM_FIRST_VARIATION_WAVEFORM_20261015` capsule:

```sh
python -I -B /path/to/evidence/wider-state-sol/check_wider_feasibility.py --source-root /path/to/HDBLAST_UNIFORM_FIRST_VARIATION_WAVEFORM_20261015
python -I -O -B /path/to/evidence/wider-state-sol/check_wider_feasibility.py --source-root /path/to/HDBLAST_UNIFORM_FIRST_VARIATION_WAVEFORM_20261015 --output WIDER_FEASIBILITY_OPTIMIZED.json
```

Both runs pass 38 exact checks. Their mathematical outputs coincide after
the optimization-mode metadata field is removed. Exact intervals, prefix
images, positive face margins, depth thresholds, coefficient boxes and
majorants are retained in the JSON outputs. The checker authenticates its
mathematical inputs and imports no producer. An independent Astra
calculation subsequently agreed on the limiting eta-face threshold and
the exact retained-coarse global/tail errors. That agreement supplements
the written proof and exact checks; it is internal review, not external
peer review or formal proof certification.

`--source-root` is optional in the original task directory layout. Only
the three authenticated paths listed in `input_pins` are needed beneath
the supplied directory; the checker reads no other capsule file. Its
outputs always go beside the checker. A portability replay copies the
planning checker outside the task tree and supplies a separate minimal
input tree containing exactly those three files. The replay verifies the
same mathematical output and does not add planning inputs or outputs to
the frozen science capsule.
