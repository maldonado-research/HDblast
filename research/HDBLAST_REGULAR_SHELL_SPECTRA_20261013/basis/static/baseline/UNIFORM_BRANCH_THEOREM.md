# A uniform regular-cone branch theorem: proposed analytic proof

Research continuation, 8 October 2026. This is a mathematical proof candidate for adversarial review, not an APS acceptance claim or evidence for a higher-dimensional origin of the Big Bang. No uploaded numerical producer is used in this proof. Novelty relative to curved Einstein–scalar brane literature remains unassessed.

## Statement and conventions

Shift the scalar stationary point to zero. Suppose W is C^4 and f is C^3 on [-r,r],

    W(0)=3k, W'(0)=0, W''(0)=w>4k>0, f(0)=f0>0.

Put U=W'^2/2-2W^2/3, m²=w(w-4k), lambda=w-4k>0. The radial equations are

    eta'=p,
    p'=U'(eta)-4a p,
    rho''=-rho[p²/4+U(eta)/6],  a=rho'/rho,

with rho(0)=0, rho'(0)=1, p(0)=0. The shell at T satisfies

    a(T)=W(beta)/3+delta f(beta)/6,
    p(T)=-W'(beta)-delta f'(beta)/2,  beta=eta(T).

There is a positive delta_* determined by the sufficient inequalities below such that every 0<delta<=delta_* has a regular solution with T>=1/k and |beta|<=L delta. It is unique in the explicitly specified small weighted Dirichlet tube below. The shell moves to T=infinity as delta approaches zero. The branch is C^1 for positive delta. Uniform constants K_eta and K_H below imply

    |beta + f'(0) delta/[4(w-2k)]| <= K_eta delta²,

    |H² - k f0 delta/3
          - [f0²/36-k f'(0)²/(24(w-2k))]delta²| <= K_H delta³.

Here H²=1/rho(T)². This is a local static mathematical result, not a dynamical stability or cosmological result. Constants are deliberately conservative, and the theorem does not assert that delta_* covers the uploaded numerical sweep.

## 1. Uniform linear inverse on expanding intervals

Let r0(y)=sinh(ky)/k, a0=k coth(ky), and

    L_m v=v''+4a0 v'-m²v.

For any nu>0, let F_nu be the regular solution of

    F_nu''+4a0 F_nu'-mu_nu F_nu=0,
    F_nu(0)=1, F_nu'(0)=0, mu_nu=nu(nu+4k).

The divergence-form equation (r0^4 F_nu')'=mu_nu r0^4 F_nu proves F_nu>0 and F_nu'>=0. Its logarithmic derivative z_nu obeys

    z_nu'=mu_nu-4a0 z_nu-z_nu²,  0<=z_nu<=nu.

The upper bound follows by the first-crossing test at z_nu=nu because a0>=k. Choose nu=lambda/2, D=m²-mu_nu>0, and s_T(y)=F_nu(y)/F_nu(T). Thus 0<s_T<=1 and 0<=s_T'<=nu s_T.

Define

    C0=4k coth(1)+nu,
    zeta=mu_nu [1-exp(-C0/k)]/C0,
    J=2/k+1/(2zeta).

For y>=2/k, the divergence-form formula, restricted to [y-1/k,y], gives z_nu(y)>=zeta: on that interval r0(t)/r0(y)>=exp[-k coth(1)(y-t)] and F_nu(t)/F_nu(y)>=exp[-nu(y-t)]. Consequently

    integral_0^T s_T(y)^2 dy <= J

for every T. The bound is conservative but independent of T.

Let G_T be the inverse of L_m with v'(0)=0, v(T)=0. The maximum principle applies at the regular center by its even radial extension: L_m s_T=-D s_T. Therefore

    |G_T N| <= (||N/s_T||_infinity/D) s_T.

Integrating (r0^4 v')'=r0^4(m²v+N), and using r0(t)/r0(y)<=exp[-k(y-t)] for 0<t<=y, proves

    |(G_T N)'| <= [1+m²/D] ||N/s_T||_infinity s_T/(4k).

For the norm ||v||_T=max(sup|v|/s_T,sup|v'|/(k s_T)), the operator norm is bounded by

    C_G=max(1/D, [1+m²/D]/(4k²)).

There is no inverse constant growing with T. Existence of this regular Green inverse follows directly by the two fundamental solutions or the finite-interval radial boundary problem; uniqueness follows from the same maximum principle.

Let F_lambda denote the regular solution at mass m². Comparing their Riccati equations gives z_lambda>=z_nu: the larger positive mass gives the initially larger derivative, and a first later crossing would have positive derivative of the difference. Integrating these logarithmic derivatives backward from T yields phi_T=F_lambda/F_lambda(T)<=s_T. Equivalently, the maximum principle applied to L_m(phi_T-s_T)=D s_T gives the same conclusion. Thus 0<phi_T<=s_T and 0<=phi_T'<=lambda phi_T, hence ||phi_T||_T<=M=max(1,lambda/k).

## 2. Eliminate the metric by a Volterra equation

For a candidate eta with p=eta', put

    F(y)=p(y)²/4+[U(eta(y))+6k²]/6,
    R=rho/r0.

The metric equation, with its exact cone data, is equivalent to

    R(y)=1-integral_0^y K(y,t) F(t) R(t) dt,
    K(y,t)=sinh(k(y-t))sinh(kt)/(k sinh(ky)).

Two useful exact facts are

    0<=K<=1/(2k),
    partial_y K=[r0(t)/r0(y)]².

Use the finite local bounds

    M2=sup_|x|<=r |U''(x)|,
    M3=sup_|x|<=r |U'''(x)|,
    L1=M3/2, C_F=k²/4+M2/12.

Set B=2M. On the ball ||eta||_T<=B|beta|, |beta|<=b, Bb<=r, we have

    |F|<=C_F B² beta² s_T²,
    |F_eta-F_tilde_eta|<=2 C_F B|beta| d s_T²,

where d=||eta-tilde_eta||_T. Choose

    epsilon=C_F B² b² J/(2k)<=1/4.

The Volterra operator has sup norm <=epsilon. Its Neumann inverse gives

    2/3<=R<=4/3,
    |a-a0|=|R'/R|<=C_F B² beta² s_T²/k.

For two candidates in the same ball, the resolvent identity and the exact derivative of K give

    |a_eta-a_tilde_eta|<=4 C_F B |beta| d s_T²/k.

For clarity, this last bound uses ||Delta R||<= (16/9)(J/(2k)) 2 C_F B|beta|d, min R>=2/3, max R<=4/3; the two contributions Delta R'/R and R' Delta R/(R R_tilde) sum to at most 2(2 C_F B|beta|d)/k. No assumption that F has one sign is used.

## 3. Uniform nonlinear Dirichlet contraction

For each finite T, use the Banach space of real C^1 functions on [0,T] with eta'(0)=0, equipped with ||.||_T. The inverse and the homogeneous solution map into this regular space. The scalar equation is equivalent to

    eta=beta phi_T+G_T N(eta),
    N=U'(eta)-m²eta-4(a-a0)eta'.

Uniform estimates on the ball are

    ||N/s_T||<=C_N beta²,
    C_N=L1 B²+4 C_F B³ b,

    ||[N(eta)-N(tilde_eta)]/s_T||<=L_N d,
    L_N=2 L1 B b+20 C_F B² b².

In addition to the previous smallness conditions impose

    q=C_G L_N<=1/2,
    C_G C_N b<=M,
    C_F B² b²/k<=k/2.

The map preserves the ball and contracts it. It gives, for each T>0 and |beta|<=b, a unique solution in that weighted ball. Metric positivity R>=2/3 and a>=k/2 follow. The integral equations imply the regular center expansions rho=y+O(y³), eta=eta_h+O(y²); bootstrap supplies the differentiability required by the equations. The Hamiltonian constraint is exact: its derivative vanishes under the two radial equations and its cone value is zero.

The statement for beta=0 is the exact solution eta=0, rho=r0. For nonzero beta, the theorem only claims uniqueness within the stated weighted tube; it does not exclude large-field solutions or unrelated global branches.

There is also a common ball ||eta||_T<=Bb, independent of beta for |beta|<=b: the same estimates with |beta| replaced by b show that all maps preserve this ball and contract there with constant q. The solutions in the smaller B|beta| balls are consequently the same unique fixed points in the common ball. This supplies C^1 dependence on beta, including at beta=0, by the contraction/implicit-function theorem on a fixed Banach space. Strict versions of the ball-preservation inequalities can be chosen without affecting feasibility. Differentiating the map at its solution gives ||eta_beta||_T<=B and

    |p_beta(T)-D_T|<=C_P |beta|,
    D_T=F_lambda'(T)/F_lambda(T),
    C_P=k C_G B [2 L1 B+20 C_F B² b].

The endpoint scalar momentum and metric derivative obey

    |p(T)-D_T beta|<=C_D beta², C_D=k C_G C_N,
    |a_beta(T)|<=C_a |beta|, C_a=4 C_F B²/k.

For completeness the regular cone IVP has the local integral equations

    rho(y)=y-integral_0^y (y-t)rho(t)[p(t)²/4+U(eta(t))/6]dt,
    p(y)=rho(y)^(-4) integral_0^y rho(t)^4 U'(eta(t))dt,
    eta(y)=eta_h+integral_0^y p(t)dt.

The contraction at the center can be made explicit without assuming a small Lipschitz constant in unsuitable raw variables. Set r_c(y)=rho(y)/y=1+y²R_c(y), B_c(y)=p(y)/y, and express

    eta(y)=eta_h+y² integral_0^1 u B_c(yu)du,
    R_c(y)=-integral_0^1 (1-u)u r_c(yu)
                  [y²u² B_c(yu)²/4+U(eta(yu))/6]du,
    B_c(y)=integral_0^1 u^4[r_c(yu)/r_c(y)]^4 U'(eta(yu))du.

On bounded continuous R_c,B_c with r_c bounded away from zero, every functional variation in these two right sides carries y². For a sufficiently short interval their joint map is a contraction; its center values are R_c(0)=-U(eta_h)/36 and B_c(0)=U'(eta_h)/5. This proves the local regular-cone IVP and its C^1 dependence on eta_h. Beyond a positive coordinate the ordinary smooth IVP theorem applies. Locally in any finite T>0, rescaling the Dirichlet Volterra/Green equations to [0,1] and using their invertible contraction linearization proves C^1 dependence on T. No uniform norm for the differentiated rescaled operators is needed. The moving-endpoint identities below supply the uniform derivative bounds actually used.

## 4. Large-radius linear logarithmic derivative

For z=F_lambda'/F_lambda, 0<=z<=lambda. Writing e=lambda-z yields the exact equation

    e'=-(w+z)e+4(a0-k)z.

For T>=Y=1/k, a0-k<=2k exp(-2ky)/(1-exp(-2)) for y>=Y. Integrating the resulting differential inequality gives

    0<=lambda-D_T<=C_L exp(-2kT),
    C_L=lambda exp(2)+8k lambda/[(1-exp(-2))(w-2k)].

Since H0²=1/r0(T)²>=4k² exp(-2kT) and H0²=R(T)² H²<=16H²/9, the safe nonlinear-endpoint bound is

    |p(T)-lambda beta|<=C_H H² |beta|+C_D beta²,
    C_H=4 C_L/(9k²).

The factor 16/9 is retained because F need not be positive. The strict w>4k assumption gives a positive growing exponent and a denominator w-2k bounded away from zero. This estimate is uniform as T approaches infinity.

## 5. Scalar junction and existence of a shell

Put f1=f'(0), L=(|f1|+1)/w, and define local bounds

    W2=sup|W''|, W3=sup|W'''|,
    F1=sup|f'|, F2=sup|f''|, F_abs=sup|f|, W_abs=sup|W|,

all on [-r,r]. Choose a positive preliminary delta_bar such that

    L delta_bar<=b,
    (C_P+W3)b+delta_bar F2/2<=w/2.

For each T>=Y, the scalar-junction function

    g_T(beta,delta)=p(T;beta)+W'(beta)+delta f'(beta)/2

has derivative >=w/2. Its value at zero is delta f1/2. The intermediate value theorem, with the stated slope bound, gives a unique root beta(T,delta) satisfying |beta|<=|f1|delta/w<=L delta. These scalar roots stay in the uniform Dirichlet tube.

Define the exact reduced boundary functions

    A(x)=W(x)f(x)/9-W'(x)f'(x)/12,
    B_fun(x)=f(x)²/36-f'(x)²/48,
    a1=k f0/3,
    C_K=L sup|A'|+sup|B_fun|.

Shrink delta_bar until C_K delta_bar<=a1/2 and

    (3 a1/2)delta_bar < (9/16) H0(Y)².

Also require sigma/6=W/3+delta f/6>0 on |beta|<=L delta, for example W_abs is not enough but W(beta)>=3k-W2 b²/2 and

    k-W2 b²/6-delta_bar F_abs/6>=k/2.

On the scalar-junction locus the Hamiltonian constraint implies

    (a-sigma/6)(a+sigma/6)=H²-delta A(beta)-delta² B_fun(beta).

The exact reduced curvature K=delta A+delta²B_fun lies in [a1 delta/2,3a1 delta/2]. At T=Y, H²>=9H0(Y)²/16>K. As T tends to infinity, H²<=9H0(T)²/4 tends uniformly to zero while K>=a1 delta/2. Thus the metric residual has a root.

## 6. Local uniqueness, without a singular IFT at delta=0

For fixed beta, regard the regular solution as parametrized by its central scalar and endpoint T. Along the Dirichlet family, the chain rule gives the exact identities

    partial_T p|beta=U'(beta)-4ap-p p_beta,
    partial_T a|beta=-H²-p²/3-p a_beta.

These use the exact bulk relation a'=-H²-p²/3. The required parametrization is nondegenerate: eta_beta(T)=1 and regular solutions have one scalar cone parameter; the C^1 Dirichlet construction therefore has nonzero endpoint derivative with respect to that parameter.

The scalar-junction derivative is

    beta_T= -[U'(beta)-4ap-p p_beta]
             /[p_beta+W''(beta)+delta f''(beta)/2].

Let

    a_bar=k coth(1)+C_F B² b²/k,
    C_Tp=M2+4a_bar kB+kB(lambda+C_P b),
    C_Tbeta=2 C_Tp/w.

Then |beta_T|<=C_Tbeta |beta|. For the metric residual F=a-W/3-delta f/6,

    F_T=-H²-p²/3-p a_beta
          +(a_beta-W'/3-delta f'/6) beta_T.

Consequently at every possible root

    F_T<=-a1 delta/2+C_slope delta²,

where

    C_slope=kB C_a L²+[(C_a+W2/3)L+F1/6] C_Tbeta L.

Impose C_slope delta_bar<=a1/4. Every root is then a strict downward crossing. Since the continuous residual begins positive and ends negative, it has exactly one root: two strict downward crossings would require an intermediate zero that is not a strict downward crossing. The ordinary IFT at each positive delta proves C^1 branch dependence there. It is never applied at T=infinity.

## 7. Uniform scalar and curvature remainders

Let

    d0=lambda+w=2(w-2k), alpha=-f1/(2d0),
    K_eta=[C_H (3a1/2)L+(C_D+W3/2)L²+F2 L/2]/d0.

The scalar junction, Section 4, and the root bound H²<=3a1 delta/2 give

    |beta-alpha delta|<=K_eta delta².

Since

    A(0)=a1,
    A'(0)=f1(4k-w)/12,
    B_fun(0)=f0²/36-f1²/48,

Taylor's theorem in the exact boundary identity gives the completely explicit bound

    K_H=|A'(0)|K_eta + (sup|A''|/2)L² + sup|B_fun'| L,

    |H²-a1 delta-[A'(0)alpha+B_fun(0)]delta²|<=K_H delta³.

Algebra reduces A'(0)alpha+B_fun(0) to

    f0²/36-k f1²/[24(w-2k)].

Thus the previously conditional expansion follows on an actually constructed regular local branch, for every positive delta satisfying the explicit sufficient inequalities. The construction does not certify the much larger detunings of a numerical experiment unless these inequalities (or sharper validated estimates) are evaluated there.

## Review status and limitations

This proof is pending independent adversarial review. All constants above are formulas in bounded derivatives and elementary positive functions; no claimed numerical value of delta_* is supplied here. The inequalities are feasible: first choose b>0 sufficiently small (all b-dependent left sides tend to zero), then delta_bar>0 sufficiently small. To make the scalar derivative condition feasible, the additional choice (C_P+W3)b<w/2 must be included in the first-stage choice of b. Set delta_*=delta_bar after all strict inequalities hold.

This result neither establishes external novelty nor solves stability, finite detuning error certification at the uploaded sample values, energy transfer, reheating, observational discrimination, or an origin of the Big Bang. The only global-in-radius assertion is existence and uniqueness inside this specified small weighted neighborhood of the regular AdS cone.
