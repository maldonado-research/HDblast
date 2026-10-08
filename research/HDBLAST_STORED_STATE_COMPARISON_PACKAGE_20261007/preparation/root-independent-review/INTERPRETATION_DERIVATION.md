# Interpretation of the registered stored-state comparison

This derivation applies to the complete registered serialized outputs. It
introduces no physical-source construction, original-array read, target
evaluation, or new physical observable. The output receipt separately records
which completed run was checked.

## An interval for the maximum L1 error

Fix a prefix and one complex coordinate X, either U or W. Let the exact saved
value be s_j, and let the certified target rectangle have midpoint t_j and
real/imaginary half-widths r_j,1 and r_j,2. Put r_j=r_j,1+r_j,2. The saved-minus-
target error rectangle has midpoint m_j=s_j-t_j and the same half-widths.
The stored rectangle upper bound is exactly

    B_j = max|delta X_j,real| + max|delta X_j,imag|
        = |m_j,real| + |m_j,imag| + r_j.

For the true target inside its certified rectangle, the reverse triangle
inequality and triangle inequality give

    B_j - 2r_j <= ||delta X_j||_1 <= B_j.

The registered summary stores S=max_j B_j. Let R be the reported maximum
complete exported L1 radius for the same X over all registered nodes. Then
R>=r_j for every node in every prefix, and selecting any node attaining S
proves

    max(0,S-2R) <= max_j ||delta X_j||_1 <= S.

All quantities can be evaluated as exact rational numbers from the serialized
summary. Using the global R is conservative. This is the Cartesian complex L1
norm |Re z|+|Im z|. It is not the complex modulus |z|; the corresponding modulus
lower bound is max(0,S-2R)/sqrt(2), while S remains a modulus upper bound.
No interpolation between retained momentum nodes is involved.

## Exact canonical residual and its perturbative scope

For real forcing g, U'=W and W'=2ikW-g imply

    c = Re U - Im W/(2k),
    c' = Re W - Im(2ikW-g)/(2k) = 0.

The prescribed target has exactly zero initial data, hence c_target=0. Therefore
the exact error invariant is the saved invariant itself:

    c_error = c_saved = Re(u_1/epsilon) - Im(w_1/epsilon)/(2k).

The implementation forms d=2k c_saved first, retains the saved invariant, and
does not project it to zero. The stable phase amplitude obeys
kA=-i delta W/2. Neither identity licenses dropping the constant imaginary
component of delta U from the incoming-state norm.

For f=exp(-ik eta)(1+epsilon U)/sqrt(2k), with the stated real epsilon,

    (f f*' - f' f*)/i
      = 1 + 2epsilon c
          + epsilon^2 [|U|^2 - Im(U* W)/k].

Thus the reported maximum of 2epsilon|c| is the exact absolute first-order
Wronskian defect. It does not report the full finite-epsilon normalization
defect, and no cancellation with the quadratic term is assumed.

## Units and finite-sum scope

The exact represented constants are

    epsilon = 3777893186295716171 / 37778931862957161709568,
    Pi = 14488038916154245685 / 4611686018427387904.

The normalized state coordinates are U=u_1/epsilon and W=w_1/epsilon.
Multiplying a normalized incoming norm by the positive exact epsilon gives
the corresponding raw u_1 or w_1 norm. The inherited positive weights enter
through mu_j=w_j k_j^2/(2Pi^2) and q_j=w_j/(4Pi^2). No weights are regenerated
or normalized to an assumed cutoff integral.

Transport quantities are the inherited scaled a_0^4/epsilon linear stress
convention. They are finite sums of the normalized linear response operators;
epsilon is not inserted again inside those operators. Their homogeneous
propagation bounds concern two exact linear evolutions from the compared
incoming states, with identical subsequent forcing and complete contacts.
They do not certify later saved numerical trajectories or the inhomogeneous
forcing/contact evaluation. Multiplication by epsilon/a_0^4 converts the
normalized stress difference to its raw linear stress convention at that time.

The finite source/grid/prefix universe consists of four distinct capsules and
twelve prefixes, with 49,152 capsule-node occurrences and 86,016 prefix node
appearances. Node coverage and exact finite sums do not imply a continuum
momentum error bound, quadrature accuracy, an ultraviolet completion, or a
full pressure/contact certificate. The separate historical calibration,
origin, and novelty conclusions are not changed by this state comparison.

This mathematical interpretation was independently checked against the frozen
transport and caller code. It is not an additional independent construction of
the physical target enclosures; those rely on the registered analytic proof,
source certificate, and pinned worker.
