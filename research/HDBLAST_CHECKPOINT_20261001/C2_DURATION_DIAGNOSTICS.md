# C2 — Radiation dominance without vacuum cancellation

**Status:** post hoc diagnostic and exact frozen-model continuation, with exploratory archived-data examples. The original registered B3 labels are preserved. No external mathematical novelty is claimed.

## Why the earlier fraction can mislead

The exact shell identity in the old normalization is

H² = Lrad + Qrad + V + Weyl + Dphi,

Lrad=sigma R/18, Qrad=R²/36,
V=sigma²/36-sigma_phi²/48+U/6,
Dphi=v²/12-sigma_phi Jphi/24-Jphi²/48.

Jphi=kappa5² j. In the friction model j=Yv. Negative V can make H² small while substantial competing terms remain. Therefore Omega_r=(Lrad+Qrad)/H² can exceed one; that fact is not evidence for a clean radiation era.

Define a conservative fraction for a conventional low-energy radiation regime:

F_abs=Lrad/[Lrad+Qrad+|V|+|Weyl|+v²/12+|sigma_phi Jphi|/24+Jphi²/48].

Proposed post hoc instantaneous thresholds are F_abs>=0.9 and |Weyl|/Lrad<=0.03, with H>0. Quadratic radiation is counted separately because the target is a conventional low-energy expansion law. These thresholds are diagnostic choices, not new observational constraints.

The scalar-free F0 obtained by omitting the final three nonnegative terms is an upper bound on F_abs. A failure of F0 at a saved point excludes F_abs there. It cannot exclude every earlier or later point.

## Saved plateau checks

Here q=Qrad/Lrad=R/(2sigma). The table uses linear radiation in the numerator.

| Saved run | F0 upper fraction | q | Original classification |
|---|---:|---:|---|
| Y=0.5 dc=0.01 | 0.720377358 | 6.891490e-4 | FAIL-Weyl |
| Y=0.7 dc=0.01 | 0.768840497 | 6.741209e-4 | FAIL-Weyl |
| Y=1.5 dc=0.01 | 0.791417212 | 4.969050e-4 | PASS-conservative |
| Y=2 dc=0.0001 | 0.765992242 | 3.974307e-4 | PASS-combined |
| Y=2 dc=0.01 | 0.765449166 | 3.961361e-4 | PASS-combined |
| Y=3 dc=0.01 | 0.705380431 | 2.769229e-4 | PASS-combined |
| Y=5 dc=0.01 | 0.588847024 | 1.626028e-4 | UNRELIABLE |

All seven saved points fail the proposed 0.9 threshold. Conditional or unreliable baseline cases retain their limitations. This table is not a whole-trajectory classification.

For the finest Y=2, dc=0.01 point:
Omega_r=1.3574158351536791, Omega_vac=-0.38632816250071517,
r=Weyl/(Lrad+Qrad)=0.021299555782589328.
|V|/(Lrad+Qrad)=0.28460561052536804 and |V|/Lrad=0.2847183530796413.
F0=0.7654491662636435. A total-radiation numerator instead gives 0.7657523883048801; that is a different diagnostic.

The full raw-array result and its execution status are in TRAJECTORY_AUDIT.md. Its new gate is defined in TRAJECTORY_PROTOCOL.md.

## Exact source-free frozen continuation

Assume flat homogeneous expansion, fixed positive tension and V=-L<0, vanishing matter transfer and scalar force, and conserved Weyl. Small velocity alone is insufficient: a nonzero scalar force or exchange can change these quantities.

Let N=ln(a/a_s), A=Lrad_s+Weyl_s, B=Qrad_s. Then

h²=-L+A exp(-4N)+B exp(-8N).

The old B4 differential equation already contains this quadratic-density term. The present closed-form duration check retains it; this is not a new physical term or an audit finding that B4 forgot it.

For x=exp(4N),

(dot x)²=16(-Lx²+Ax+B),  ddot x=8A-16Lx,

x_turn=[A+sqrt(A²+4LB)]/(2L),
N_turn=ln(x_turn)/4,
Delta tau_turn=acos[(2L-A)/sqrt(A²+4LB)]/(4 sqrt(L)),

on the initial expanding branch x=1. The initial condition requires A+B-L>0. If A>0, the quadratic correction relative to the linear-only continuation obeys

0 <= DeltaN <= LB/(4A²).

For Y2, using the **local saved plateau vacuum** gives exact N_turn=0.31936177170908386; treating all initial radiation as a^-4 gives 0.3194316926118252. The archived numerical maximum is 0.3195150947587777. These nearby values do not certify a new continuum solution. Using the different B4 static-vacuum value yields approximately 0.3193674073 and 0.3194373288 instead; do not mix the vacuum choices.

A deliberately large quadratic control L=0.2,A=1,B=10 gives N_turn=ln(10)/4=0.5756462732485115. A separately implemented RK4 second-order evolution agrees and conserves its invariant; see saved output.

## Frozen dominance window

Let q=B/Lrad_s, w=Weyl_s/Lrad_s and ell=|V|/Lrad_s. The scalar-free F0>=0.9 condition is

|w| + q exp(-4N) + ell exp(4N) <= 1/9.

Its endpoints follow from ell x²-(1/9-|w|)x+q<=0. There may be both an early high-density entry and a late vacuum exit. For a future continuation restrict to N>=0. Negative roots relative to a saved plateau are a mathematical extrapolation, not actual earlier data.

For the Y2 plateau the frozen scalar-free gate lies entirely before N=0; it never passes in that assumed frozen future. This still does not describe the earlier rolling trajectory.

## Energy injection with endpoint quantities

The radiation Ward identity during rolling is Rdot+4HR=S, so

a_f^4 R_f=a_s^4 R_s+J,  J=integral a^4 S dtau.

For a specified **linear-radiation** endpoint target |V_f|/(alpha_f R_f)<=epsilon, alpha_f=sigma(phi_f)/18>0, it requires

J>=a_f^4 |V_f|/(epsilon alpha_f)-a_s^4 R_s,

with a zero lower bound if the right side is negative and only nonnegative injection is allowed. Endpoint V and tension may change; they must be evaluated there. This is not sufficient to produce a sustained era or to close the Weyl/source equations.

For a target whose numerator is **total** radiation alpha_f R_f+beta R_f², use instead

R_required=2(|V_f|/epsilon)/[alpha_f+sqrt(alpha_f²+4beta|V_f|/epsilon)],
J>=a_f^4 R_required-a_s^4 R_s.

The larger linear requirement must not be called a necessary bound for the quadratic criterion.

Conditionally keeping the Y2 plateau tension/vacuum unchanged, one additional e-fold with the chosen linear ratio <=0.1 demands J/R_s>=154.45095358631988 (a_s=1). This is an endpoint budget target, not a realized source. No nonzero jv should be assumed for a strictly frozen scalar.

The next calculation must determine source, vacuum and Weyl together from consistent dynamics. The archived stand-in transfer law alone does not establish a physical reheating mechanism.
