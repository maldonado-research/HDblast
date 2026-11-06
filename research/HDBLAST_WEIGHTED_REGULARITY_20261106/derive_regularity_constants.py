"""Exact low-order scalar coercivity constants and both high-order joins.

Continuous integration by parts, dotted closure and Fourier membership
arguments are in the companion proofs, not certified by finite arithmetic.
"""
from fractions import Fraction as Q
from math import factorial, comb


def need(condition, message):
    if not condition:
        raise ValueError(message)


def scalar_coercivity():
    checks=[]
    def check(condition,name):
        need(condition,name)
        need(name not in checks,'Duplicate scalar coercivity check')
        checks.append(name)
    ell=Q(1,10**6);kap=ell/100
    R0,R1=Q(129),Q(131)
    a0,a1=Q(111,1000),Q(112,1000)
    p0,p1=Q(3,25000),Q(7,50000)
    V1=Q(28,100000)
    sqrt_tail_lower=Q(433,500)
    check(sqrt_tail_lower**2<Q(3,4),'tail square-root lower bound')
    kinv=R1*ell/(2*sqrt_tail_lower*16)
    check(kinv<Q(473,10**8),'inverse tail-coordinate derivative')
    measure_coefficient=2*a1*R1*ell/sqrt_tail_lower
    check(measure_coefficient<8,'weighted radial density nonincreasing for t>=4')
    alpha=Q(1,2)+Q(1,24)+a1*kinv
    check(alpha<Q(55,100),'tail first-derivative coefficient')
    Cmax=Q(19)/R0**2+p1**2/3+3*a1**2
    check(Cmax*kinv**2<Q(1,100),'tail zeroth-order coefficient')
    check(1/(1-Q(55,100)-Q(1,100))<Q(5,2),'Hardy absorption constant')
    b_lambda=a1+Q(21)/(3*a0*R0**2)
    check(b_lambda<Q(116,1000),'radial reconstruction value coefficient all modes')
    check(2*kinv/(3*a0)<Q(285,10**7),'tail derivative adjoint coefficient')
    check(Q(116,1000)*kinv<Q(55,10**8),'tail derivative potential coefficient')
    tail_B=Q(5,2)*R1**2*ell**2/(2*Q(3,4))
    check(tail_B<Q(3,10**8),'tail angular-trace forcing coefficient')
    tail_A=Q(5,2)*4*Q(285,10**7)*625
    check(tail_A+Q(3,10**8)<Q(181,1000),'complete tail energy coefficient')
    tail_N=Q(5,2)*(4*Q(55,10**8)+6+Q(55,100)*2)
    check(tail_N<18,'complete tail interior coefficient')
    check(2*Q(181,1000)<Q(4,10),'two-tail assembly energy coefficient')
    check(1+2*18==37,'two-tail assembly interior coefficient')

    # Compact interior cutoff r in[1/6,1/5].
    check(45*2/(R0*ell)<700000,'interior cutoff derivative')
    check(36*2/(R0*ell)<560000,'interior logarithmic weight derivative')
    cmax=2*R1**2/3
    bmax=2*a1*R1**2/p0**2
    bpmax=2*R1**2/p0**2*(1/R0**2+p1**2/3+10*a1**2+2*a1*V1/p0)
    check(cmax<11441,'monopole c coefficient')
    check(bmax<Q(267)*10**9,'monopole b coefficient')
    check(bpmax<Q(1545)*10**9,'monopole b derivative coefficient')
    check(a1+1/(a0*R0**2)<Q(113,1000),'monopole radial reconstruction coefficient')
    check(Q(7,4)**2>3,'angular trace square-root ceiling')
    K=(Q(11441)+Q(1545)*10**9+2*(Q(267)*10**9)*560000
       +(Q(267)*10**9)*Q(113,1000)+2*(Q(267)*10**9)*700000
       +Q(7,4)*a1*(Q(267)*10**9))
    check(1296*K<9*10**20,'interior monopole mixed energy coefficient')
    check(Q(8,10)+37**2*9*10**20<Q(13)*10**23,'closed monopole coercivity coefficient')

    for lam,gain,total in [(8,6650,154000000),(24,1873,43400000)]:
        check(6*R1**4/(lam*(lam-3))<gain**2,f'angular tracefree gain lambda{lam}')
        check(Q(4,10)+37*625*gain<total,f'closed nonconstant coercivity lambda{lam}')
    check(kap**-4==10**32,'defining-function scale conversion')
    return {'status':'PASS','named_checks':checks,'named_check_count':len(checks),
            'exact':{'interior_K':str(K),'interior_1296K':str(1296*K),
                     'lambda0_coercivity':str(Q(13)*10**55),
                     'lambda8_coercivity':str(154000000*10**32),
                     'lambda24_coercivity':str(43400000*10**32)}}


def join(source):
    names=[]
    def check(condition,name):
        need(condition,name)
        need(name not in names,'Duplicate regularity join check')
        names.append(name)
    check(source['status']=='PASS','derived source status')
    check(source['orders']=={'background':16,'global_bump':16,'source_Dy':12,
                            'source_Ds':12,'momentum_C_Ds':14,
                            'known_ADM_conversion_profiles_Dy':13},
          'source and coefficient derivative margins')
    check(source['conventions']['Ds']=='q Dy'
          and source['conventions']['momentum_G']==
          'kappa r rho^(3/2) exp(1/r) (r_lambda,sqrt(lambda) g_lambda/rho)',
          'actual transformed source and density normalization')
    for key,values in source['source_Dy_derivative_upper'].items():
        check(len(values)==13 and all(Q(v)>=0 for v in values),
              'actual source proper derivative envelope through12 '+key)
    scalar=source['scalar_H12_rational_majorant']
    for key,limit in (('f0',10**24),('f2',10**25),('f4',10**25)):
        check(Q(scalar[key])<limit,'actual scalar transformed H12 bound '+key)
    for mode in ('0','2','4'):
        check(Q(source['momentum_H12_rational_majorant'][mode])<10**17,
              'actual momentum transformed H12 bound '+mode)
    c=list(map(Q,source['momentum_C_Ds_opnorm_upper']))
    check(len(c)>=15 and all(v<Q(1,10**10) for v in c[:15]),
          'whole annulus momentum matrix coefficient jets through14')
    # Low-order momentum gap, preserving the actual rho³dy measure.
    t=Q(4,3)
    exp4=sum((t**j/factorial(j) for j in range(6)),Q(0))
    exp4+=(t**6/factorial(6))/(1-t/7)
    low_checks=[
        ('momentum tangential algebraic coefficient36',4*Q(10,3)*36**2>131**2),
        ('momentum radial algebraic coefficient8',Q(1000,111)*(Q(1,3)+Q(5*36,3*129))<8),
        ('momentum joint algebraic coefficient37',8**2+36**2<37**2),
        ('momentum monopole algebraic coefficient4',12*Q(111,1000)**2*4**2>1),
        ('momentum derivative lower coefficient one fifth',2*(Q(112,1000)**2+Q(24,129**2))<Q(1,25)),
        ('momentum tail coordinate inverse coefficient',131**2<2*100**2),
        ('momentum tail square root envelope',2<Q(10,7)**2),
        ('momentum tail logmeasure monotonicity',2*Q(112,1000)*Q(1,10**6)*131*Q(10,7)<Q(42,10**6)),
        ('momentum tail absorption factor',Q(1,5)*Q(1,10**4)*Q(1,4)==Q(5,10**6)),
        ('momentum inverse tail absorption',1/(1-Q(5,10**6))<Q(100001,100000)),
        ('momentum global coercivity coefficient1000',333+Q(100001,100000)*(666+Q(1,10000))<1000),
        ('momentum exponential four thirds ceiling4',exp4<4),
        ('momentum selected natural energy below6e33',36*Q(222,100)*10**15*Q(16,9)*4*10**16<6*10**33),
        ('momentum selected weighted potential below6e52',1000*10**16*(6*10**33)==6*10**52),
    ]
    for name,condition in low_checks:
        check(condition,name)
    # Scalar natural energy is bounded by the same selected compact comparison.
    exp1=sum(Q(1,factorial(j)) for j in range(11))
    exp1+=Q(1,factorial(11))/(1-Q(1,12))
    expthird=sum(Q(1,3**j*factorial(j)) for j in range(9))
    expthird+=Q(1,3**9*factorial(9))/(1-Q(1,30))
    check(exp1<Q(272,100) and Q(4,3)**4*expthird<Q(442,100),
          'scalar natural energy exponential and contrast ceilings')
    kappa=Q(1,10**8)
    scalar_energy=kappa**-4*Q(272,100)*Q(442,100)*(558*10**37)
    scalar_potential=13*10**23*kappa**-4*scalar_energy
    check(scalar_potential<9*10**128,'selected scalar weighted potential below9e128')
    scalar_high=20*(3*10**25+9*10**128)
    check(scalar_high<2*10**130,'actual scalar transformed H16 below2e130')
    check(Q(445,100)*10**2000*scalar_high<10**2131,
          'selected scalar weighted H16 and physical metric H14 geometric ceiling')
    # Momentum whole-line normal coefficients derive from the computed C jets.
    cstar=max(c[:14]);k1=4*cstar;k0=2**12*cstar*cstar+2*cstar
    eps=max(k1,k0)
    check(eps<Q(1,10**6),'momentum normal lower coefficient absorption hypothesis')
    theta=Q(4096,10**6)
    check(Q(1,2)/(1-theta)<Q(503,1000)
          and theta/(1-theta)<Q(4114,10**6),
          'momentum Fourier highest derivative absorption')
    check(15*(1+Q(503,1000)+Q(4114,10**6))**2<10**2,
          'momentum H14 derivative square sum factor10')
    check(max(comb(14,j) for j in range(15))<64**2,
          'momentum Bessel convention conversion factor64')
    momentum_high=10*(6*10**52+3*10**17)
    check(momentum_high<10**54,'actual momentum transformed H14 below1e54')
    # Full physical/weighted tensor connection and angular reconstruction.
    physical=[
        ('momentum normalized degree2 polynomial budget',Q(5,4)**2*48<32**2),
        ('momentum normalized degree4 polynomial budget',Q(61,12)**2*Q(720,19)<32**2),
        ('momentum mixed angular tensor budget',2*Q(32*8)**2/8<1000**2),
        ('momentum tracefree angular tensor budget',Q(32*(64+24),5)<1000),
        ('momentum initial matrix Taylor jet budget',8*54+4<1000),
        ('momentum physical covariant recurrence budget',max(824+215*j for j in range(13))<4000),
        ('momentum physical output coefficient ceiling',10**6*4000**13<10**53),
        ('momentum physical exponential multiplier ceiling',24**24<10**34),
        ('momentum output Sobolev square sum factor',14<4**2),
        ('momentum potential Sobolev square sum factor',15<4**2),
        ('momentum weighted potential covariant recurrence budget',max(663+215*j for j in range(14))<4000),
        ('momentum weighted potential coefficient ceiling',10**4*4000**14<10**55),
    ]
    for name,condition in physical:
        check(condition,name)
    for key in ('q_over_rho','q_times_a','q_y','alpha_Ds_log_U','alpha_physical_Ds_log_m'):
        values=source['whole_annulus_geometric_coefficient_Ds_upper'][key]
        check(len(values)>=15 and all(Q(v)<18 for v in values[:15]),
              'momentum physical whole annulus coefficient jets '+key)
    check(all(v<18 for v in c[:15]),'momentum physical matrix coefficient jets')
    check(10**300*momentum_high<10**354 and 10**60*momentum_high<10**114,
          'selected physical momentum H13 and weighted potential H14 ceilings')
    check(10**354<10**355,'rounded public physical momentum H13 ceiling')
    return {'status':'PASS','named_checks':names,'named_check_count':len(names),
            'exact':{'scalar_natural_energy_upper':str(scalar_energy),
                     'scalar_weighted_potential_upper':str(scalar_potential),
                     'scalar_transformed_H16_upper':str(scalar_high),
                     'scalar_geometric_H14_public_upper':str(10**2131),
                     'momentum_C_Ds_max':str(cstar),'momentum_K1_jet_upper':str(k1),
                     'momentum_K0_jet_upper':str(k0),'momentum_normal_epsilon_upper':str(eps),
                     'momentum_transformed_H14_upper':str(momentum_high),
                     'momentum_geometric_H13_proved_upper':str(10**354),
                     'momentum_geometric_H13_public_upper':str(10**355),
                     'momentum_geometric_weighted_H14_upper':str(10**114)}}


def physical_cauchy(source):
    names=[]
    def check(condition,name):
        need(condition,name)
        need(name not in names,'Duplicate physical Cauchy check')
        names.append(name)
    profiles=source['quadratic_ADM_conversion_profile_Dy_upper']
    for key,ceiling in (('u_w_velocity_over_rho',10**73),('u_vt_over_rho',10**84)):
        values=list(map(Q,profiles[key]))
        check(len(values)==14 and all(value<=ceiling for value in values),
              'known ADM conversion radial profile through13 '+key)
    # Spatial dimension four: tr(h1)=(-4+3*2)uY and volume'/volume=tr/2.
    volume_first=Q(-4+3*2,2)
    check(volume_first==1,'canonical volume first derivative coefficient')
    # With t=-3w/(2rho), differentiate e A^(3/2)B^(-3/2)
    # and e A^(1/2)B^(-1/2), A'= -4uY, B'=2uY.
    radial=2*Q(4,3)*(Q(3,2)*(-4)-Q(3,2)*2)*Q(-3,2)
    angular=2*Q(-2,3)*(Q(1,2)*(-4)-Q(1,2)*2)*Q(-3,2)
    check(radial==36 and angular==-6,'physical ADM known quadratic metric momentum coefficients')
    check(-2*volume_first==-2,'physical normal scalar momentum quadratic coefficient')
    trace_eigenvalue=2*(1-Q(4,3))
    check(trace_eigenvalue==Q(-2,3) and trace_eigenvalue**2<2**2,
          'physical ADM linear tensor conversion spectral norm2')
    check(131**4*Q(1,10**6)<18**2 and 14<4**2
          and Q(445,100)<5 and 36**2+3*6**2<200**2,
          'known ADM full geometric tensor and support allowances')
    Bcross=5*200*4*18*14*1000**13*(10**80)**13*10**200*10**84
    check(Bcross<10**1600,'known ADM full geometric H13 cross budget')
    I2=10**2131+2*10**355+2*Bcross
    check(I2<10**2132,'full physical second Cauchy jet per copy bound')
    check(2*I2**2<(10**2133)**2,'full physical second Cauchy jet reflected union bound')
    return {'status':'PASS','named_checks':names,'named_check_count':len(names),
            'exact':{'known_K2_radial_profile_coefficient':str(radial),
                     'known_K2_angular_profile_coefficient':str(angular),
                     'known_Pi2_profile_coefficient':'-2',
                     'known_ADM_cross_H13_upper':str(Bcross),
                     'combined_per_copy_triangle_upper':str(I2),
                     'physical_second_Cauchy_per_copy_public_upper':str(10**2132),
                     'physical_second_Cauchy_reflected_union_public_upper':str(10**2133)}}
