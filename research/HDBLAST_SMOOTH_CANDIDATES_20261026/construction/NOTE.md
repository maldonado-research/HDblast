# C4 smoothing of the frozen full23 initial profiles

This is a replacement **nominal approximation**, not an exact solution. It
changes the displacement and velocity used by the prior free wave; any
comparison with the actual registered initial data requires a new error ledger.
The accompanying candidate omits the affine bulk/shell sources, just as the
prior finite candidate did. Its inhomogeneous boundary defect is not established to vanish.

## Construction

Take the published 23-vector continuous P1 tables p_u,p_v on [0,R], R=7/4,
with their Gram norm, node set and mixed odd/even radial reflection parity.
Both tables vanish identically on [0,65/64]. Continue the final P1 interval
affinely beyond R. Reflect each continuation with its declared parity to the
negative axis. Choose delta=1/2048 and let

    k_delta(s)=35/(32 delta) [1-(s/delta)^2]^3, |s|<=delta,
               0 otherwise.

The kernel is even, nonnegative, mass one and C2. Let u_s=k_delta*p_u,
v_s=k_delta*p_v and du_s=partial_r u_s. The kernel is C2 and p_u,p_v are
continuous P1; distributionally their second derivatives are finite sums of
point masses. Hence u_s'' and v_s'' are finite linear combinations of C2
kernel translates, so u_s,v_s are C4 globally. Their even/odd parity is
preserved, and du_s has opposite parity. They vanish identically for
|r|<65/64-delta=2079/2048>1. The TT entries remain zero.

The affine continuation is used only to define convolution near the outer
endpoint; the candidate's domain and comparisons remain the frozen causal
footprint |r|<=R. Do not silently treat this continuation as the actual field
outside R. Its global linear growth is harmless for convolution with compact
kernel.

Each P1 table can be written as a sum of hinge functions
`p(r)=sum_i Delta m_i (r-r_i)_+` on the positive half-line, since the first
segment is identically zero. Because the kernel width is smaller than the
zero collar, this formula produces the correct reflected convolution on
r>=0. Its smoothing is a sum of `Delta m_i h_delta(r-r_i)`, with

    h_delta(z)=integral k_delta(s) (z-s)_+ ds.

For z<=-delta, h=0; for z>=delta, h=z. On the inner interval h is an exact
rational polynomial of degree 8 in z/delta. More generally, for n>=0,

    H_n(z)= integral k_delta(s) (z-s)_+^n/n! ds.

The evaluator integrates the polynomial exactly and uses H_{k+1} for the kth
primitive of u_s,v_s; for du_s it uses H_k. The integral is normalized at r=0
because of the open zero collar. The checker proves the polynomial recurrence
`dH_n/dz=H_{n-1}`, the unit mass/even moments, and matching of h through
four derivatives at both z=+-delta by exact rational identities. All outputs
at rational arguments are exact rational numbers.

## New B0+B1 candidate and boundary

Build the new free term from the same smoothed displacement and velocity:

    B0(x,r) = [u_s(r+x)+u_s(r-x)
                    + integral_{r-x}^{r+x} v_s(s)ds]/2,

where p_u,p_v and their primitives are reflected with the original 23-slot
mixed parity. For the first feedback use the frozen constant A,B,D,M matrices,
these same u_s,v_s, and the derivative-consistent du_s. Construct F and G
as in the frozen ridge formula, applying the matrix action **before** the
mixed reflection. Evaluate the reflected bulk Duhamel term plus the homogeneous
Robin lift using the authenticated `angular_ridge.py` implementation. The
B0+B1 expression is the new candidate; it is not an iteration to convergence.

Since u_s,v_s are C4, the interior reflected wave/Duhamel formula is at least
C2 in x,r; the compact zero collar makes joins at the initial reflected corner
compatible. The zero input on Dirichlet rows and parity imply P_D B0=0 at the
shell. For the first feedback, the inherited half-line operator yields

    P_D B1=0,   P_N (B1)_r = M P_N B0       at r=0,

with natural one-sided derivatives. Thus the total candidate has exact
homogeneous Dirichlet rows, but its Robin residual against the full equation
includes `-M P_N B1 - G`; no exact inhomogeneous boundary solution is claimed.
The standard half-line smoothness claim assumes the inherited mixed reflection
Duhamel identities and enough one-sided data regularity. The fully coupled
reflected forcing need not be smooth across r=0 for odd rows; smoothness there
is asserted only from the physical side. Near x=r, the branch-specific terms
involve q=r-x and vanish with their normalized primitives for |q|<1, by the
zero collar, so the branch switch introduces no join kink. It does not prove
that the unknown exact solution or its full C2 error is bounded.

## Reproduction and scope

Run

    python smooth_candidate.py --repo-root /path/to/HDblast
    python -O smooth_candidate.py --repo-root /path/to/HDblast

The standard-library checker pins the frozen profile and ridge code before
executing the latter, verifies the exact polynomial identities and C4 joins,
zero collar, TT zeros, derivative consistency, and one active-shell mixed
boundary example. It computes B0+B1 at the earlier kink point. The general
smoothness claim follows from the kernel and Green-kernel argument, not from
this finite point. The approximate field depends on the prior P1 input; the
program does not recompute that input, the actual field, the bulk quadratic
source or a physical Einstein residual.
