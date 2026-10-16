# Quadratic gauge covector: universal identity and conditional reduction

This working derivation fixes the meaning of the gauge source used by the
geometric evaluator. It is an algebraic implementation clarification; no claim
of an externally new curvature or gauge-variation formula is made.

Let the reference connection be the Levi-Civita connection of a fixed smooth
nondegenerate metric g0, and put g(e)=g0+e q on an interval on which it remains
nondegenerate. Every derivative below without a superscript uses the reference
connection. Define

    P_bcd = (nabla_c q_db + nabla_d q_cb - nabla_b q_cd)/2,
    P^a_cd = (g0)^ab P_bcd,
    C^a(g) = g^cd (Gamma[g]-Gamma[g0])^a_cd,
    C_b(g) = g_ba C^a(g).

The exact difference-of-connections identity gives

    (Gamma[g(e)]-Gamma[g0])^a_cd = e g(e)^ab P_bcd.

Consequently lowering with the *varying* metric cancels the inverse metric in
this index, and yields the exact, off-shell identity

    C_b(g(e)) = e g(e)^cd P_bcd.

Writing q^cd=(g0)^ca(g0)^db q_ab, it follows that

    (C_b)'(0)  = (g0)^cd P_bcd = (C1)_b,
    (C_b)''(0) = -2 q^cd P_bcd.                         (1)

No harmonic condition or field equation is needed for (1). In contrast, if
QCup^a denotes the second variation of the *upper* gauge vector, then

    QCup^a = -2 q^a_d C1^d - 2 q^cd P^a_cd,
    (g0)_ba QCup^a = -2 q_bd C1^d - 2 q^cd P_bcd.

The varying-lowering contribution +2 q_ba C1^a is exactly what distinguishes
(C_b)'' from this background-lowered upper Hessian. The distinction must remain
visible when testing arbitrary metric jets.

For a reduction Ered_ab(g,phi)=E_ab(g,phi)-nabla^g_(a C_b)(g), with parentheses
including one half, the affine-path Hessian is universally

    Qred_ab = QE_ab - nabla^0_(a (C_b)''(0)
                       + 2 P^c_ab (C1)_c.              (2)

Here C(g0)=0, and the final term is the second variation of the connection in
the covariant derivative acting on a covector. Equation (2) is the appropriate
off-shell manufactured-test formula. A reduction defined instead with a fixed
reference derivative would have a different off-shell cross term; implementations
must state which derivative they use.

For the registered actual first field, C1=0 as a field. Thus (1) coincides with
the background-lowered upper QC, its reference derivatives coincide as well,
and (2) reduces to the retained on-solution convention

    QEred_ab = QE_ab - nabla^0_(a QC_b),
    QC_b = -2 q^cd P_bcd.                              (3)

The first harmonic solve supplies the field identity and hence its derivative
identity. Arbitrary independent metric jets do not supply that premise. A
rectangular interval extension of the polynomial in (3) can still enclose the
source of the actual harmonic field, provided the actual jets are enclosed;
the off-shell points in that rectangle need not themselves solve the gauge
equations. Such enclosure logic must not be confused with a universal identity
between (2) and (3) on arbitrary jets.

The expansions are in second parameter derivatives, g=g0+e q+e^2 k/2. The
coefficient of e^2 is one half of each displayed Hessian. Raising a lower time
index in the orthonormal Lorentz frame introduces a minus sign. All ambient
sphere tensors must be tangent before applying the frame contractions.
