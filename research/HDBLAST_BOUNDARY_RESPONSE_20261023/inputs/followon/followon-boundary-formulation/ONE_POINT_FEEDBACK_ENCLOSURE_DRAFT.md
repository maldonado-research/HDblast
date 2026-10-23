# DRAFT: a boundary-inclusive one-point enclosure before C1 reconstruction

PROPOSED_UNVERIFIED_FOLLOWON. Not part of the frozen checkpoint. No source value or response has been evaluated here. This conditional argument uses only the retained operator definitions, coefficient caps and exact zero wedge. It requires independent analytic review before use.

At the proposed apex (a,b)=(3/2,1/4), d=a-b=5/4. Its closed physical backward cone is backward-closed including shell histories. Write S=KF+LG and W=S+TW, with the SAME reflected K, negative L lift and Robin feedback as the retained construction. In X over this cone, ||K||<=a and ||L||<=max(1,d/3)=1. With sigma=203/10^6, ||M||<=m=1/2000 and boundary-value trace norm<=3,

    q_D=a sigma+3m=3609/2000000<1,
    S_D_cap=a F_H+G_H=52162547172519/2,
    ||W||_X,D <= B_D=S_D_cap/(1-q_D)=52162547172519000000/1996391.

The complete point-value feedback is bounded without numerically reconstructing W:

    |W(a,b)-S(a,b)|_H <= [(a^2/2)sigma+3dm] B_D
                        =877739181271977213/15971128 (approximately 5.495787031e+10).

The retained exact zero wedge x+r<1 improves this. W and its first derivatives vanish there; therefore its bulk feedback f=A W_x+B W_r+D0 W vanishes there, and its boundary feedback g=M P_N W is zero for s<1. The total unfused direct+reflected integration area is a^2. Each cone removes the same triangle s+u<1 of area1/2 (since d>1), so the remaining K norm weight is (a^2-1)/2=5/8. The boundary lift has effective length d-1=1/4. Thus

    |W(a,b)-S(a,b)|_H <= E_D=[(5/8)sigma+3(1/4)m] B_D
                         =209432626897663785/15971128 (approximately 1.311320195e+10).

The mean functional lambda=probability-S3 mean of covariant yy has H dual norm1, so any honest interval [S_minus,S_plus] obtained from actual source interval integration immediately gives

    lambda W(a,b) in [S_minus-E_D,S_plus+E_D].

This cheap first test needs continuous source enclosures and exact folded cell/boundary integration, but no numerical C1 reconstruction of W. It may be too broad to resolve a sign; no sign is asserted. Certified C1 reconstruction remains necessary for the proposed evaluated feedback iterates/full numerical waveform, or use a rigorously bounded full residual. Do not confuse a point-value allowance with the X norm or multiply it by another factor3. The factor3 was already used in the Robin trace bound; H-dual lambda adds no factor.

Premises and signs: FULL_SECOND_RESPONSE_CONSTRUCTION eqs10-13, exact zero-wedge restriction at lines219-223, full coefficient bounds and F_H/G_H from the frozen capsules. Lg=-integral g, M=-ell rho B_y. Triangle folding and zero-wedge removal change no boundary condition. This is a standard conditional norm application, not a novelty or physical-force claim.
