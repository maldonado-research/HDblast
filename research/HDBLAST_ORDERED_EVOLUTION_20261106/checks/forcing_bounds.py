"""Finite arithmetic for the analytic H12 quadratic forcing input map.

The Sobolev embedding, tensor formula and product allocation arguments
are supplied in QUADRATIC_FORCING.md, not inferred from these checks.
"""
from fractions import Fraction as F


def evaluate(data,check):
    Z=F(data['collar_length'])
    check('actual_collar_length',Z==F(6,10**6))
    for rank in range(17):
        check('tensor_sphere_embedding_rank:'+str(rank),(1+3*rank*rank)+4*rank+2<1000)
    check('radial_Hilbert_trace_factor',1/Z+2+Z<10**6)
    check('rounded_H3_embedding',10**6<10**data['H3_embedding_decade'])
    check('bilinear_H12_product_factor',4*3**12*10**data['H3_embedding_decade']<10**data['bilinear_H12_decade'])
    check('linear_H12_multiplier_factor',2**24*13<10**(2*data['linear_H12_decade']))
    check('acceleration_multiplier_factor',2**48*13<10**18)
    check('covariant_coefficient_conversion',150200+64*12+6<151000<data['coefficient_decade'])
    check('acceleration_budget',2+2*10**219<10**data['acceleration_multiplier_decade'])
    check('full_quadratic_forcing_budget',10**data['bilinear_H12_decade']*(9+2*10**data['acceleration_multiplier_decade'])<10**data['forcing_multiplier_decade'])
    forcing_exponent=data['forcing_multiplier_decade']+2*data['coefficient_decade']
    check('forcing_coefficient_exponent',forcing_exponent==data['forcing_map_decade'])
    initial_exponent=forcing_exponent+2*data['initial_N14_decade']
    check('initial_forcing_exponent',initial_exponent==data['initial_exact_decade'])
    check('initial_forcing_rounded_budget',initial_exponent<data['initial_rounded_decade'])
    # Derivative allocation for U in H14 and W in H12. When all twelve
    # derivatives hit U, put W in L-infinity instead.
    for j in range(13):
        h14_h12=[]
        h13_h13=[]
        for left in range(j+1):
            for right in range(j-left+1):
                if left<=11:
                    valid=left+3<=14 and right<=12
                else:
                    valid=left==12 and right==0 and 3<=12
                h14_h12.append(valid)
                h13_h13.append(min(left,right)+3<=13 and max(left,right)<=13)
        check('H14_H12_all_allocations_order:'+str(j),all(h14_h12))
        check('H13_H13_all_allocations_order:'+str(j),all(h13_h13))
    return {'tensor_embedding_max_rank':16,'H3_Linfinity_bound':'10^'+str(data['H3_embedding_decade']),
            'coefficient_C':'10^'+str(data['coefficient_decade']),
            'acceleration_H12_multiplier':'10^'+str(data['acceleration_multiplier_decade'])+' C',
            'quadratic_forcing_H12_map':'10^'+str(forcing_exponent)+' N14(t)^2',
            'initial_forcing_H12_exact_decade':initial_exponent,
            'initial_forcing_H12_rounded_bound':'10^'+str(data['initial_rounded_decade']),
            'scope':'Evaluated initial force and conditional H12 forcing map; evolved N14(t) is not bounded by this replay.'}
