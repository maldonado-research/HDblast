#!/usr/bin/env python3
"""Bind the displayed causal inequalities to the accepted rational report."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

REPORT_SHA = 'beff99f6d7c2792cdcfc1540c8949036ca454ae36418e5ea5a80f51ec453b969'

BODY = r'''Set $t_0=1/4$, $r_0=5/4$, and
\[
 C=\{0\le x\le t_0,\ |r-r_0|\le t_0-x\}.
\]
Every backward cone from a point of $C$ stays inside $C$. Moreover
$r-x\ge1$ and $x+r\le3/2<8$ there. Thus the restriction of the actual
global auxiliary $W$ satisfies
\[
 W|_C=\KK F|_C+\KK(AW_x+BW_r+D_0W)|_C.
\]
Both the affine shell lift and the Robin feedback vanish exactly by domain
of dependence. The material-shell first solution still determines $F$;
no boundary law or bulk coupling has been changed.

Use the same norm restricted to $C$. Endpoint differentiation and the
retarded integral give
\[
 \norm{\KK f}_{X_C}\le t_0\norm{f}_{C^0(C;\HH)},\qquad
 |\KK f(t_0,r_0)|_{\HH}\le\frac{t_0^2}{2}\norm{f}_{C^0(C;\HH)}.
\]
The retained individual full-state coefficient bounds imply
\[
 \sigma=10^{-4}+10^{-4}+3\times10^{-6}=\frac{203}{10^6},
 \quad q_C=t_0\sigma=\frac{203}{4000000}<1,
 \quad\norm{W}_{X_C}\le B_C:=\frac{t_0F_{\HH,C}}{1-q_C}.
\]
This follows directly from the restriction inequality for the already
specified $W$; a separate local inverse theorem is unnecessary. All 23
amplitudes enter the source norm and the coefficient majorants. At the apex,
the complete feedback allowance is
\[
 |\lambda(W-\KK F)(t_0,r_0)|\le
 E:=\frac{t_0^2}{2}\sigma B_C
   =\frac{203F_{\HH,C}}{127993504}.
\]
It bounds all feedback iterations together.

The calculation partitions $C$ into 20 exact-area intersections of a
uniform rectangular mesh with spacing $h=t_0/4=1/16$. In time row
$i=0,\ldots,3$, the radial indices are $k=i,\ldots,7-i$. The two end
intersections have area $h^2/2$ and every intermediate intersection has
area $h^2$. Their union is $C$ and their interiors are disjoint; the total
area is $t_0^2=1/16$. Source intervals on the rectangular supersets enclose
the actual source on each intersection, including the inherited first-field
$C^3$ uncertainty. These are range-times-area integral bounds, not sampled
quadrature errors. In particular,
\[
 \frac12\sum_j |C_j|\,l_j
 \le\lambda(\KK F)(t_0,r_0)
 \le\frac12\sum_j |C_j|\,u_j,
 \qquad l_j\le\lambda(F)|_{C_j}\le u_j.
\]
Negative cell lower endpoints are included in the signed sum.
The full-state source cap uses exact probability-sphere Gram matrices,
including both mixed-block multiplicities and both tensor orientations.
The final accepted cap is $F_{\HH,C}=17381837650838$.

\begin{center}
\begin{tabular}{ll}
\toprule
Quantity & Safe outward inequality\\
\midrule
$\lambda(\KK F)(t_0,r_0)$ & $5.4\times10^7<\lambda(\KK F)<6.4\times10^{10}$\\
Complete feedback allowance & $E<2.757\times10^7$\\
$\lambda(W)(t_0,r_0)$ & $2.64\times10^7<\lambda(W)<6.4\times10^{10}$\\
\bottomrule
\end{tabular}
\end{center}
The reproducibility reports retain the exact rational endpoints. The interval
is wide; its strict positive lower endpoint is the certified feature.

For a continuous-region consequence, set $\epsilon=10^{-6}$ and
\[
 R_\epsilon=\{|x-1/4|\le\epsilon,\ |r-5/4|\le\epsilon\}.
\]
This entire rectangle lies in the enlarged backward cone with apex
$(t_+,r_0)=(t_0+2\epsilon,r_0)$, whose minimum $r-x$ is
$1-2\epsilon>0$. The accepted global source cap gives
\[
 B_+=\frac{t_+\,17387513156573}{1-t_+\sigma},\qquad
 |\lambda(W)(x,r)-\lambda(W)(t_0,r_0)|\le2\epsilon B_+.
\]
Both mean derivatives are bounded by $B_+$ and the joining segment stays
inside the convex enlarged cone. The allowance is below $8.695\times10^6$.
Consequently, throughout the complete closed rectangle,
\begin{equation}
 \boxed{\quad1.7\times10^7<\lambda(W)(x,r)<6.5\times10^{10}.
 \quad}                                                    \label{eq:rectangle}
\end{equation}
No derivative bound for $F$, point-sampling sign inference, or shell-contact
approximation enters this extension. The result concerns the covariant $yy$
frame component's probability-sphere mean in the registered gauge. The
variables and field normalization are the retained scaled ones, with no SI
calibration or physical-force interpretation.
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    raw = args.report.read_bytes()
    if hashlib.sha256(raw).hexdigest() != REPORT_SHA:
        raise ValueError('unaccepted causal report bytes')
    data = json.loads(raw)
    lo, hi = map(Q, data['actual_W_yy_probability_mean_interval'])
    slo, shi = map(Q, data['evaluated_Duhamel_mean_interval'])
    rectangle = data['continuous_rectangle_corollary']
    rlo, rhi = map(Q, rectangle['actual_W_yy_probability_mean_interval'])
    checks = {
        'born_display_lower': slo > 54000000,
        'born_display_upper': shi < 64000000000,
        'feedback_display_upper': Q(data['coupling_point_value_error_upper']) < 27570000,
        'point_display_lower': lo > 26400000,
        'point_display_upper': hi < 64000000000,
        'rectangle_display_lower': rlo > 17000000,
        'rectangle_display_upper': rhi < 65000000000,
        'rectangle_variation_display_upper': Q(rectangle['mean_variation_upper']) < 8695000,
        'global_norm_bridge_value_weight_three': lo / 3 > 8800000,
        'exact_twenty_cells': data['partition']['cells'] == 20,
        'full_coupling_retained': data['full_coupling']['all23_bulk_sectors_retained'] is True,
    }
    if not all(checks.values()):
        raise ValueError('displayed inequality does not follow from accepted report')
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir/'CAUSAL_RESULT.tex').write_text(BODY)
    report = {
        'status': 'PASS_DISPLAYED_CAUSAL_AND_NORM_INEQUALITIES',
        'accepted_report_sha256': REPORT_SHA,
        'generator_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'causal_section_sha256': hashlib.sha256(BODY.encode()).hexdigest(),
        'checks': checks,
        'exact_global_X_norm_lower_from_value': str(lo/3),
        'scope': 'Display arithmetic only; accepted analytic and actual-source interval certificate is separate.'
    }
    (args.output_dir/'MANUSCRIPT_ARITHMETIC.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
