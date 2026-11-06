"""Exact global first-packet majorants and arithmetic transport ledger.

The geometric norm and derivative recurrences are analytic premises from
the proof. This code recomputes their finite integer/rational inequalities.
"""
from fractions import Fraction as F
from math import factorial
import importlib.util


def evaluate(root,data,check):
    engine=root/'inputs/derive_source_jets_06e.py'
    spec=importlib.util.spec_from_file_location('transport_pinned_source_jets',engine)
    if spec is None or spec.loader is None:
        raise ValueError('Source Taylor engine unavailable')
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    check('source_sixteen_stored_orders',m.N>=16)
    check('fixed_ell',m.ELL==F(data['ell']))
    j=m.background()
    bump,_=m.bump_bounds()
    Q,_=m.compose_bump(bump,j['b'])
    src=m.source_jets(j,Q)
    rho,a,u=j['rho'],j['a'],src['u']
    w=src['w_velocity']
    profiles={
        'u_value':(u,14),
        'scalar_value':(src['v'],14),
        'rho2u_value':(rho**2*u,14),
        'w_velocity':(w,13),
        'scalar_velocity':(src['vt'],13),
        'rho2w_velocity':(rho**2*w,13),
        'mixed_radial_velocity':(rho**2*(u.dy()+4*a*u),13),
    }
    ceilings={}
    for name,(jet,order) in profiles.items():
        # source_jets consumes one derivative for scalar values/velocities;
        # the order sixteen store therefore covers every entry requested.
        bound=max(factorial(n)*jet.c[n]/m.ELL**n for n in range(order+1))
        check('global_profile:'+name,bound<10**data['first_profile_decade'])
        ceiling=m.decade_ceiling(bound)
        check('strict_profile_ceiling:'+name,bound<10**ceiling)
        ceilings[name]={'through_order':order,'decade_ceiling':ceiling}
    for stage in range(14):
        check('raw_derivative_slot_branch:'+str(stage),
              stage+1+1+3*(2+stage)<=60)
    raw_first=4*10**6*1000**14*10**(14*data['warped_background_y_decade']+
             data['first_profile_decade']+data['angular_decade'])
    check('raw_first_C_and_H_budget',raw_first<10**data['raw_first_decade'])
    check('geometric_support_volume',F(40*131**4,10**6)<10**6)
    check('volume_measure_conversion',129**4>1)
    check('rho_a_z_product_budget',2**14*10**(2*data['background_z_decade'])<10**data['conversion_atom_decade'])
    for stage in range(14):
        check('conversion_derivative_slot_branch:'+str(stage),
              stage+2+10*(2+stage)<data['conversion_branch_factor'])
    conversion=10**3*data['conversion_branch_factor']**14*10**(14*data['conversion_atom_decade'])
    check('warped_to_product_conversion',conversion<10**data['conversion_decade'])
    linear=10**6*2**14*10**(data['flow_decade']+data['normalization_decade'])
    hessian=10**6*3**14*10**(data['flow_decade']+2*data['normalization_decade'])
    check('linear_map_C14_budget',linear<10**data['linear_map_decade'])
    check('hessian_map_C14_budget',hessian<10**data['hessian_map_decade'])
    check('shared_map_budget',max(data['linear_map_decade'],data['hessian_map_decade'])<data['shared_map_decade'])
    check('first_converted_budget',data['raw_first_decade']+data['conversion_decade']<=data['converted_first_decade'])
    check('second_converted_budget',data['raw_second_decade']+data['conversion_decade']<=data['converted_second_decade'])
    check('linear_Leibniz_factor',2**28*15<10**(2*data['linear_Leibniz_decade']))
    check('bilinear_Leibniz_factor',3**28*15<10**(2*data['bilinear_Leibniz_decade']))
    first=10**(data['linear_Leibniz_decade']+data['shared_map_decade']+data['converted_first_decade'])
    check('first_one_or_two_copy_rounding',2*first<10**data['final_first_decade'])
    second_linear=10**(data['linear_Leibniz_decade']+data['shared_map_decade']+data['converted_second_decade'])
    second_quadratic=3*10**(data['bilinear_Leibniz_decade']+data['shared_map_decade']+2*data['converted_first_decade'])
    check('second_quadratic_budget',second_quadratic<10**data['quadratic_budget_decade'])
    check('second_combined_intermediate_budget',second_linear+second_quadratic<10**data['quadratic_budget_decade'])
    check('second_one_or_two_copy_rounding',2*(second_linear+second_quadratic)<10**data['final_second_decade'])
    # Value H14 needs z14/state2 (total16); velocity H13 needs z13/state2
    # (total15), after exact even-frame time cancellation at t=0.
    for label,order in [('second_value',16),('second_velocity',15)]:
        check('retained_flow_derivative_order:'+label,order<=data['flow_total_order'])
    check('support_inside_literal_theta_collar',F(5,2)*m.ELL<6*m.ELL)
    check('support_away_from_boundary',F(1,2)*m.ELL>0)
    for normal_order in range(data['flow_total_order']+1):
        for state_order in range(data['flow_total_order']-normal_order+1):
            check('flow_mixed_jet:'+str(normal_order)+':'+str(state_order),
                  8*factorial(normal_order)*factorial(state_order)*10**(6*normal_order)*16000**state_order
                  <10**data['flow_decade'])
    return {'first_profile_ceilings':ceilings,
        'raw_first_C14_C13_and_W1_bound':'10^'+str(data['raw_first_decade']),
        'retained_raw_W2_bound':'10^'+str(data['raw_second_decade']),
        'geometric_to_product_norm_factor_bound':'10^'+str(data['conversion_decade']),
        'first_straightened_Cauchy_norm_bound_one_or_two_copies':'10^'+str(data['final_first_decade']),
        'second_straightened_Cauchy_norm_bound_one_or_two_copies':'10^'+str(data['final_second_decade']),
        'continuous_premises':'Analytic chart, tensor norm conversion and product recurrences remain proof premises.'}
