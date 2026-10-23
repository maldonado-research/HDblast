Portable rational branch-comparison replay
=========================================

Run with Python3 (standard library only):

  python replay_branch.py
  python -O replay_branch.py
  python -I -S -B verify_portable.py --output /tmp/PORTABLE_CHECKS.json

The first two commands replay all exact polynomial, residual, damping, barrier,
shell, and necessary-root-enclosure arithmetic and compare every exact numerical
field with NUMERICAL_REFERENCE.json. Optional --output FILE writes the complete
conditional report. Optional --compare-branch-report FILE additionally compares
the corresponding exact fields with a full original branch report. No absolute
path is encoded, and neither script imports project or third-party packages.

verify_portable.py copies only replay_branch.py and NUMERICAL_REFERENCE.json
into a fresh temporary directory. It checks normal and optimized execution,
byte-identical reports from an unrelated working directory, and wrong-coefficient,
excessive-damping, inward-rounding, and corrupted-reference controls. All eight
negative controls must fail at the intended explicit require() check. Each child
uses -I -S -B, plus -O in optimized mode, to isolate imports and suppress cache
writes. --output selects the receipt path; its default is PORTABLE_CHECKS.json
in the caller's working directory. The source directory may be read only.
No proof
check uses Python assert, so python -O cannot remove proof checks.

Upstream premises, not conclusions proved here
----------------------------------------------
The retained theorem supplies existence and uniqueness in its stated regular
Dirichlet tube; the beta0 +/- (799/250)*delta^2 enclosure; regular cone data and
expansions; scalar sign/monotonicity; and the positive-a constraint branch.
This portable replay neither authenticates their provenance nor reproves the
theorem. It does not prove an independently admissible cone or a globally
unique invariant manifold. The exact model is embedded in the script.

The normalized slope z=p/eta is defined on the retained nonzero negative scalar
branch. With x=eta, h=rho^-2, M=U'/x, and a>0,
 x'=xz, h'=-2ah, z'=M-4az-z^2.
The script checks 0<M<4 and U<0 for |x|<=9e-5. Given the upstream regular cone
and scalar monotonicity, z enters at0, stays in[0,2], a remains positive, and
|x| stays below the endpoint magnitude. Thus the path crosses h=1/100 before
the shell. On that smaller h domain the constraint supplies the proved a bounds.

Degree4 graph and exact defect
------------------------------
Z_ROWS and A_ROWS are the complete rational coefficient listings in the script;
each line is i j coefficient for x^i*h^j. They independently match all seven
exact polynomial fields U,M,Z,A,B,R,D of the parent graph certificate. P=xZ is
degree5. There is no unspecified truncated remainder: the replay constructs
the full residual polynomials
 B=2hZ_h-4Z,
 R=M-Z^2-xZ_x*Z+A*B,
 D=h+x^2*Z^2/12-U/6-A^2,
and checks exact vanishing of their coefficients through total degree4.

For e=z-Z, the exact secant identity is
 e'=R+B*D/(a+A)-K*e,
 K=4a+z+Z+xZ_x-B*x^2*(z+Z)/(12*(a+A)).
The last term is favorable because B<0 and z+Z>=0. Coefficientwise bounds give
K>=37/20 on the declared domain. The polynomial majorant of
|R|+|B||D|/(amin+Amin) bounds all x powers by eps=9e-5 and reduces h powers
above5 using h<=h0=1/100. If Cj are those coefficients, define
 Ej=Cj/((37/20)-2*j*(3/20)),
 E(h)=2*(h/h0)^6+sum_j Ej*h^j.
The denominators are positive. Moreover (37/20)-12*(3/20)=1/20>0, so
E'+K*E is at least the residual forcing along the actual characteristic.
At entry |e|<2, and scalar comparison proves |e|<=E(h).

Necessary enclosure of an already existing root
------------------------------------------------
The shell constraint defines h=H(beta) exactly. Let
g(beta)=beta*Z(beta,H(beta))+W'(beta)+delta*c/2.
At the retained existing scalar-junction root, |g|<=eps*E(Hmax). The replay
checks a positive lower bound for g' over the retained beta interval. Four
rational seed steps select a center; the mean-value theorem gives radius
(|g(seed)|+eps*E(Hmax))/min(g'). Outward terminating-rational endpoints enclose
this radius. The final radius is1901/200000000000000000000 =9.505e-18.

This is a necessary enclosure joined to retained existence. It is not a new
endpoint-map interval-Newton inclusion. It supplies no propagated initial-data
error, source/solution error, total response, or observable sign conclusion.
