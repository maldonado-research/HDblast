# Direct, geometry-weighted Taylor integration: independent bound review

This is a mathematical design and fabricated-fixture review. No registered physical source, source arrays, response arrays, or scientific output were evaluated. It is not a certificate for the HDBLAST experiment.

## 1. Assumptions and target

On a real cell `t=t0+x`, `0<=x<=h`, suppose

`u'=w`, `w'=i Ω w−ε g`,

where Ω is real and constant on this cell; k and ε are real and nonzero. All numerical parameters must be exact or enclosed. The geometry varies:

`L(t)=−1/t`,

`A(t)=3 L(t)^3/(k ε)`, `B(t)=L(t)^2/(k ε)`, `C(t)=L(t)/ε`, `D(t)=B(t)−i C(t)`.

The modal action is

`J = ∫ Re[A(t) u(t)+D(t) w(t)] dt`.

The interval cannot contain `t=0`. Freezing A or D at a cell midpoint changes the target and is not allowed. The construction below directly integrates the actual geometry-weighted action. It does not use a response-density endpoint primitive.

## 2. Polynomial reconstruction with an explicit ODE defect

Let the real source be enclosed by `g(x)=p(x)+e(x)`, where `p=Σ a_j x^j` and

`|e(x)| <= E_g(x) = Σ e_j x^j`, `e_j>=0`.

This includes a uniform remainder as `e_0`, interval errors of polynomial coefficients as `e_j`, and errors inherited from nested source moments. Choose degree N at least the degree of p. Choose central initial modes `û0, ŵ0` with absolute complex errors at most `δu, δw`. Define

`W_N=Σ_{j=0}^N w_j x^j`, `w_0=ŵ0`,

`(j+1) w_{j+1} = i Ω w_j−ε a_j` for `0<=j<N`,

and `U_{N+1}=û0+Σ_{j=0}^N w_j x^(j+1)/(j+1)`.

Thus `U'=W` exactly. Compute the full defect polynomial

`R=W'−i Ω W+ε p`.

All its coefficients below degree N vanish algebraically. The degree-N coefficient is `−i Ω w_N+ε a_N`. If a polynomial coefficient is an interval, bound the resulting full defect instead of asserting cancellation. Put `|R(x)| <= E_R(x)=Σ r_j x^j`, `r_j>=0`, and `c_j=r_j+|ε| e_j`.

Real Ω is essential: the propagator has modulus one. Variation of constants gives

`|w−W|(x) <= δw + Σ c_j x^(j+1)/(j+1)`,

`|u−U|(x) <= δu + x δw + Σ c_j x^(j+2)/[(j+1)(j+2)]`.

There is no exponential Grönwall factor. These are genuine Taylor/reconstruction remainder bounds, not estimates from convergence of floating-point samples. At Ω=0 a sufficiently high-degree forced polynomial is exact. There is no division by Ω and no low-frequency singularity.

## 3. Geometry polynomials and their exact rational remainders

Let `c=t0+h/2`, `z=x−h/2`, `a=h/2`, `r=a/|c|<1`, and `m=P+1`. Expand each inverse power around the cell midpoint:

`t^(−p) = c^(−p) Σ_{n=0}^P (−1)^n binom(p+n−1,n) (z/c)^n + remainder`, `p=1,2,3`.

Its uniform absolute remainder is bounded by `|c|^(−p) T_p(r)`, where

`T_1 = r^m/(1−r)`,

`T_2 = r^m[(m+1)/(1−r)+r/(1−r)^2]`,

`T_3 = r^m[binom(m+2,2)/(1−r)+(2m+3)r/(2(1−r)^2)+r(1+r)/(2(1−r)^3)]`.

These follow by differentiating the geometric series; every term displayed is nonnegative. They are rational for dyadic or rational geometry and do not require a transcendental library. Convert the centered polynomials to powers of x by a finite binomial identity. Form `A_P,B_P,C_P,D_P` with the exact signed normalizations above. Denote their uniform remainders by `ρA,ρB,ρC`; a conservative complex remainder is `ρD=ρB+ρC`.

For `d=min(|t0|,|t0+h|)>0` on a cell of one sign,

`Abar=3/(|k ε| d^3)`,

`Dbar=1/(|k ε| d^2)+1/(|ε| d)`

bound the true coefficients. Hence `Apbar=Abar+ρA`, `Dpbar=Dbar+ρD` bound their polynomial approximations.

## 4. Direct integral and a complete local error budget

Compute exactly, or with outward interval arithmetic,

`J_P = Re ∫ [A_P U+D_P W] dx`.

If the product polynomial has coefficients `f_j`, its exact integral is `Σ f_j h^(j+1)/(j+1)`. This is the geometry-weighted nested polynomial moment; A(t) and D(t) have not been frozen. There is no numerical quadrature error for the polynomial itself. Rounding enclosures must be included explicitly as `B_arith` when rational arithmetic is not used.

Let `G_j=|a_j|+e_j`, and let `U0bar, W0bar` bound the true initial complex moduli. Then

`Wtruebar=W0bar+|ε|Σ G_j h^(j+1)/(j+1)`,

`Utruebar=U0bar+h W0bar+|ε|Σ G_j h^(j+2)/[(j+1)(j+2)]`.

Use the exact decomposition

`A u+D w−A_P U−D_P W`

`=A_P(u−U)+D_P(w−W)+(A−A_P)u+(D−D_P)w`.

A sufficient local enclosure is `J in J_P ± B_local`, with

`B_state = Apbar[h δu+h² δw/2+Σ c_j h^(j+3)/((j+1)(j+2)(j+3))]`

`          +Dpbar[h δw+Σ c_j h^(j+2)/((j+1)(j+2))]`,

`B_geometry=h(ρA Utruebar+ρD Wtruebar)`,

`B_local=B_state+B_geometry+B_arith`.

This decomposition avoids bounding a high-phase Taylor polynomial by the sum of its large coefficient magnitudes in the geometry error. It bounds the actual states by the unitary propagator instead. An alternative decomposition using true A,D and bounds for U,W is also valid; the minimum of independently valid bounds may be taken.

## 5. Inherited mode uncertainty and nested source moments

The displayed `δu,δw` must include the enclosure inherited from the preceding cell, not just local arithmetic. At a cell boundary the conservative carry is

`δw_next=δw+Σ c_j h^(j+1)/(j+1)+outward endpoint evaluation radius`,

`δu_next=δu+h δw+Σ c_j h^(j+2)/((j+1)(j+2))+outward endpoint evaluation radius`.

All local action enclosures can then be summed. Errors are deterministic bounds; correlations must not be discarded by statistical root-sum-of-squares rules.

If a nested source moment is `q(x)=q0+∫ f`, with input enclosure `|f−f_P|<=ρf` and inherited moment error `δq0`, then `|q−q_P|<=δq0+xρf`. It contributes a constant and a linear coefficient to `E_g` after including its exact source prefactor. A second nested integral contributes the corresponding quadratic envelope, including its inherited initial moment. For a product, use `|ab−ã b̃|<=|ã|δb+|b̃|δa+δaδb` with certified bounds, then propagate that envelope. A Cauchy truncation bound alone does not include initial nested-moment uncertainty or rounding of source coefficients.

## 6. High and low phase; genuine quadrature certification

For small `|Ω|h` this recurrence needs no phase quotient. For high phase, choose a predeclared sufficiently large N or subdivide the cell, retain all inherited mode enclosures, and verify the resulting `B_local`. If the remainder is too large, the result is inconclusive; a producer tolerance or matching two floating-point answers is not a substitute for the bound. Every operation used to obtain defects, parameter ranges and final sums must be exact or outward rounded.

If Ω itself is enclosed around a central value, add `δΩ |W_N(x)|` to the ODE-defect envelope. The same unit-modulus argument holds for any real, even varying, true Ω. If Ω may have an imaginary part, this proof must be replaced by a suitable growth bound.

The formulas for the reconstruction and the displayed geometry-tail prefactors assume exact k and ε. If either is uncertain, include its effect in the full ODE defect and in the coefficient remainders as well: for example, an ε radius contributes `δε |p|` to the state defect and the reciprocal prefactors in A,D require outward interval bounds. Such parameter uncertainty cannot be covered merely by rounding the final answer.

The exact polynomial integration plus `B_state+B_geometry` is already a certified local quadrature remainder for the true direct integrand. An optional derivative-based check has the following form: for a degree-q Taylor polynomial of the entire direct integrand centered in the cell and a certified `(q+1)`st derivative bound M, the integrated remainder is at most `M h^(q+2)/[2^(q+1)(q+1)!(q+2)]`. Merely differentiating a sampled interpolant does not certify M.

## 7. Analytic source disk and the proposed M=64

The following check concerns the source formulas supplied to this independent reviewer. It assumes `z=η+4`, or an equivalent unit derivative mapping (`dz/dη=1`), and the stated complex domains. Nonunit source scaling requires its derivative factors and a new bound.

Let `B(z)=exp(1−1/(1−z²))`. Real cell centers lie in `[-1/2,1/2]` and the complex analytic disk radius is `R=1/8`. Hence `|z|<=5/8`, `|1−z²|>=39/64`, and `|η|>=27/8`. The only bump singularities are `z=±1`, outside these disks; the geometry pole η=0 is also outside them.

Put `v=z²`, `s=25/64`. For `|v|<=s`,

`Re[1/(1−v)] >= 1/(1+s)=64/89`.

For example, multiplying the desired inequality by `|1−v|²(1+s)` leaves `s+(1−s)Re(v)−|v|²>=s−(1−s)s−s²=0`.

Consequently `|B|<=exp(25/89)<4/3`. A purely rational proof of the last step uses `n!>=2·3^(n−2)` for `n>=2`:

`exp(x)<=1+x+x²/[2(1−x/3)]` at `x=25/89`,

whose exact value is `57051/43076<4/3`.

Writing `d=1−z²`, the derivatives are

`B' = −2z B/d²`,

`B'' = B[4z²/d^4−2/d²−8z²/d³]`.

For `h=B` and `h=zB`, respectively, use `g=4L²h−2Lh'−h''`, `|L|<=8/27`, and the triangle inequality. Exact rational bounds are

`|g_B| <= 2737816576/62462907 < 44 < 64`,

`|g_zB| <= 2321190208/62462907 < 38 < 64`.

Therefore M=64 is a valid common Cauchy bound under the stated domains and unit source scaling. This proof is considerably sharper than using `|B|<=exp(1+64/39)` or even `|B|<=e`; those weaker choices do not justify M=64.

With panel halfwidth `a=1/128`, radius `R=1/8`, and a degree-24 source Taylor polynomial, the uniform analytic truncation remainder is

`ρg <= 64 (a/R)^25/(1−a/R)=1/(15·2^90)`.

Coefficient enclosures, inherited nested moments and all subsequent integration errors remain separate additions. This rational source-domain proof is not a numerical evaluation of the physical bump.

## 8. Fabricated checks and limitations

The accompanying `test_direct_taylor_bounds.py` uses only the Python standard library, exact fractions, explicitly fabricated rational polynomials and parameters. It tests full recurrence defects, variable inverse-geometry tail bounds, polynomial integration, low/zero/high phase cases, inherited initial-mode uncertainty, a linear nested-moment error envelope, and the rational analytic source inequalities. A much higher-order independently enclosed reconstruction checks the direct integral fixtures. There are no calls to the physical source or scientific input arrays.

Passing these checks validates this design's algebra and fixtures. It does not validate the registered HDBLAST calculation, choose a scientific acceptance target, establish novelty, or provide an observational Big Bang mechanism.
