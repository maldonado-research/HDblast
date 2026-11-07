"""Exact arithmetic supporting the conditional background linear estimate.

No giant exponential is evaluated as floating point. The logarithmic
bound is compared exactly after canceling its shared power of ten.
"""
from fractions import Fraction as F
from math import factorial


def evaluate(root,data,check):
    C=10**data['coefficient_decade']
    omega=F(data['bundle_connection_bound'])
    domega=F(data['bundle_connection_derivative_bound'])
    gamma=F(data['base_connection_entrysum_bound'])
    qentry=F(data['inverse_entrysum_bound'])
    fiber=F(data['fiber_norm_conversion'])
    frame=F(data['moving_frame_conversion'])
    lt=F(data['lower_t_multiplier'])*C
    lw=F(data['lower_w_multiplier'])*C
    l0=F(data['lower_0_multiplier'])*C
    # Temporal product connection and its base Christoffels vanish. The
    # spatial bound includes the whole map from four E fibers to E.
    check('temporal_lower_map',fiber*C<=lt)
    connection_spatial=frame*qentry*(2*omega+gamma)
    check('spatial_connection_reserve',connection_spatial<=F(data['spatial_connection_reserve']))
    check('spatial_lower_map',fiber*frame*C+connection_spatial<lw)
    zero_bound=fiber*data['state_dimension']*C+fiber*omega*C+qentry*(domega+omega*omega)
    check('zeroth_lower_map',zero_bound<l0)
    check('entrysum_inverse_bound',F(25,16000)<qentry)
    check('symmetric_tensor_fiber_conversion',2<=fiber* fiber)
    # Mixed coefficient derivatives after the actual first-jet substitution.
    connection_derivative=10**data['connection_derivative_decade']
    check('connection_Cauchy_ordinary_derivatives',10**6*factorial(16)*100**16<connection_derivative)
    inverse_derivative=10**data['inverse_ordinary_derivative_decade']
    check('inverse_ordinary_state_derivatives',16**2*factorial(16)*10**data['inverse_aggregate_decade']<inverse_derivative)
    check('five_factor_derivative_partitions',5**16<10**12)
    ordinary=10**20*(C*connection_derivative**2+inverse_derivative*connection_derivative**2)
    check('substituted_coefficient_ordinary_budget',ordinary<10**data['substituted_ordinary_decade'])
    check('required_mixed_order',12+2<=16)
    check('required_connection_order',12+1<=16)
    for stage in range(12):
        check('coefficient_covariant_slot_branch:'+str(stage),5*(stage+1+32*16)<10**4)
    check('coefficient_covariant_branch_uniform',5*(13+32*16)==2625 and 2625<10**4)
    fixed=10**6*10**data['substituted_ordinary_decade']*(10**4*connection_derivative)**12
    check('fixed_covariant_coefficient_budget',fixed<10**data['fixed_covariant_decade'])
    check('fixed_slot_norm_reserve',8*4<10**6)
    check('moving_frame_Cauchy_derivatives',8*factorial(19)*100**19<10**data['moving_frame_derivative_decade'])
    moving=3**12*fixed*10**(2*data['moving_frame_derivative_decade'])
    check('moving_covariant_coefficient_budget',moving<10**data['coefficient_H12_decade'])
    # Actual background g^-1=rho^-2 diag(-1,1,1,1,1) in the moving frame.
    beta=F(1,129);a=F(14,125)
    cd=2*beta*a
    check('actual_background_divergence',cd==F(data['background_divergence_bound']))
    t=F(data['candidate_time'])
    check('fixed_time_interval',t==F(13,8000000))
    check('tanh_upper_use',0<t<F(1,100))
    # At=Ct=Bd=0 at the background; |tanh t|<=t suffices for moving term.
    K=F(135)+18000*cd+36000*lt+18000*lw+135*l0+F(9,4)*t
    check('background_K_budget',K<10**data['K_decade'])
    check('K_coefficient_multiplier',36000*F(data['lower_t_multiplier'])+
          18000*F(data['lower_w_multiplier'])+135*F(data['lower_0_multiplier'])==F(239400))
    # With K<10^(Cdec+6), T/2*10^6=13/16. No conversion of the
    # approximately 150000-digit exponent to a decimal string is needed.
    check('K_relative_decade',data['K_decade']==data['coefficient_decade']+6)
    check('log_G_relative_decade',data['log_G_decade']==data['coefficient_decade'])
    log_ratio=t*F(10**(data['K_decade']-data['log_G_decade']),2)
    check('log_G_exact_ratio',log_ratio==F(13,16))
    check('log_G_budget',log_ratio<1)
    check('first_energy_norm_conversion',F(data['first_norm_prefactor'])**2>2*269**2)
    check('first_energy_input_conversion',F(data['first_norm_prefactor'])>=1)
    check('first_H1_velocity_envelope',F(data['first_norm_prefactor'])*10**data['initial_first_decade']<10**data['first_propagated_decade'])
    check('outer_wall_support_arithmetic',F(5,2)+2*F(13,8)<6)
    return {'lower_t_bound':str(data['lower_t_multiplier'])+' * 10^'+str(data['coefficient_decade']),
        'lower_w_bound':str(data['lower_w_multiplier'])+' * 10^'+str(data['coefficient_decade']),
        'lower_0_bound':str(data['lower_0_multiplier'])+' * 10^'+str(data['coefficient_decade']),
        'background_Cd_exact_upper':str(cd),'K_bound':'10^'+str(data['K_decade']),
        'log_G_bound':'10^'+str(data['log_G_decade']),
        'log_G_exact_ratio_to_bound':str(log_ratio),
        'first_H1_plus_velocity_bound':'10^'+str(data['first_propagated_decade'])+' * exp(10^'+str(data['log_G_decade'])+')',
        'coefficient_covariant_derivatives_through_12_bound':'10^'+str(data['coefficient_H12_decade']),
        'scope':'Conditional background first linear wave only; shell/outer flux premises and finite regularity remain analytic assumptions.'}
