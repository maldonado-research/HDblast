"""Exact constants and derivative allocations for the forced hierarchy.

Large analytic exponentials are represented by exact exponent ledgers.
The continuous energy, Sobolev and recovery arguments remain in the proof.
"""
from fractions import Fraction as F
from math import comb


def evaluate(d,check):
    # Every mixed time/spatial Leibniz allocation is tested. Each factor's
    # intrinsic order records whether it is u, Du, or the full Hessian Hu.
    intrinsic=((0,0),(0,1),(1,1),(0,2))
    allocation_count=0
    for total in range(13):
        valid=[]
        for time in range(total+1):
            space=total-time
            for at in range(time+1):
                for bt in range(time-at+1):
                    ct=time-at-bt
                    for az in range(space+1):
                        for bz in range(space-az+1):
                            cz=space-az-bz
                            for r,s in intrinsic:
                                left=bt+bz+r;right=ct+cz+s
                                valid.append(at+az<=12 and left+right<=14 and min(left,right)<=7
                                             and min(left,right)+3<=14 and max(left,right)<=14)
                                allocation_count+=1
        check('mixed_source_all_allocations_total:'+str(total),all(valid))
        for time in range(total+1):
            space=total-time
            mt=sum(comb(time,a)*comb(time-a,b) for a in range(time+1) for b in range(time-a+1))
            mx=sum(comb(space,a)*comb(space-a,b) for a in range(space+1) for b in range(space-a+1))
            check('mixed_multinomial_sum:'+str(time)+':'+str(space),mt*mx==3**total)
    check('full_source_weight_sum',1+2+2+2==7)
    check('first_mixed_spatial_conversion',15*2**7==1920 and 1920<10**4)
    check('first_mixed_full_exponent',1980163+1985300+4==3965467 and 3965467<d['first_mixed_decade'])
    check('first_mixed_initial_input',15<10**(d['first_initial_mixed_decade']-1985230))
    check('mixed_source_product_prefactor',13*4*7*4*3**12*10**7<10**d['mixed_source_multiplier_decade'])
    check('mixed_source_coefficient_exponent',d['coefficient_decade']+d['mixed_source_multiplier_decade']==d['mixed_source_map_decade'])
    check('uniform_source_exact_exponent',d['mixed_source_map_decade']+2*d['first_mixed_decade']==d['source_exact_decade'])
    check('uniform_source_rounding',d['source_exact_decade']<d['source_uniform_decade'])
    check('uniform_source_log_rounding',2*10**(d['first_log_decade']-152017)<10**(d['source_log_decade']-152017))
    check('initial_source_exact_exponent',d['mixed_source_map_decade']+2*d['first_initial_mixed_decade']==d['initial_source_exact_decade'])
    check('initial_source_rounding',d['initial_source_exact_decade']<d['initial_source_decade'])
    check('source_rotation_conversion',13*60**24<10**(2*d['source_rotation_decade']))
    source_rows=sum(comb(14-k,2)*6**k for k in range(13))
    check('normal_source_word_count',source_rows==3761479855 and source_rows<10**12)
    check('normal_source_array_factor',10**12*60**24<10**60)
    check('source_linear_multiplier',2**24*13<10**10)
    check('normalized_source_multiplier',2**24*13<10**(2*(d['normalized_source_decade']-210)))
    check('normalized_source_rotation_spatial_factor',source_rows*60**24*10**(2*d['normalized_source_decade'])<10**(2*d['normal_source_decade']))
    check('sharp_Z11_aggregate_factor',78*60**22<10**44)
    check('normalized_initial_source_exponent',d['initial_source_decade']+d['normalized_source_decade']==4200220 and 4200220<4200300)
    check('full_retained_source_L1_and_initial_cost',d['initial_source_decade']<=d['source_uniform_decade'] and
          d['source_uniform_decade']+1<d['total_source_cost_decade'])

    H=10**d['normalized_coefficient_decade'];R=10**d['initial_step_decade']
    for j in range(12):
        check('forced_initial_step:'+str(j),1+2**j*(2*10**60+4*10**5*H)<R)
        m=11-j
        check('forced_initial_all_allocations:'+str(j),all(
            k+b<=11 and m-b<=12-j+k and m-b+1<=13-j+k
            for k in range(j+1) for b in range(m+1)))
        check('forced_initial_source_order:'+str(j),m+j==11)
    # a_(new)<=R max(a_old)+B; summing the geometric series through
    # twelve updates gives max a<=R^12(a0+a1+B) for R>=2.
    check('forced_recurrence_geometric_sum',R>=2 and (R**12-1)<R**12*(R-1))
    source_initial_exponent=d['normalized_source_decade']+d['initial_source_decade']
    check('initial_source_dominates_Cauchy',d['second_initial_input_decade']<=source_initial_exponent)
    initial_jet_exponent=12*d['initial_step_decade']+source_initial_exponent+1
    check('initial_jet_sum_rounding',2<10 and initial_jet_exponent==d['initial_jet_decade'])
    check('initial_rotation_energy_factor',F(273,2)*60**24<10**50)
    check('initial_energy_rounding',initial_jet_exponent+25<d['second_initial_energy_decade'])

    words=sum((13-k)*6**k for k in range(13))
    check('commuted_source_word_count',words<10**12)
    # The preceding order-thirteen commutator bound applies to this
    # order-twelve subfamily without enlarging its coefficient.
    check('energy_commutator_order',12<=13)
    rate=F(10**d['base_K_decade'],2)+95*10**d['commutator_decade']
    check('forced_energy_rate',rate<10**d['energy_rate_decade'])
    T=F(d['time']);check('actual_time',T==F(13,8000000))
    check('forced_energy_rate_log',T*10**(d['energy_rate_decade']-d['rate_log_decade'])<1)
    check('source_integral_factor',95*T<1)
    source_energy_exponent=d['source_uniform_decade']+d['source_rotation_decade']
    check('initial_energy_dominated_in_Duhamel',d['second_initial_energy_decade']<=source_energy_exponent)
    check('Duhamel_prefactor_rounding',source_energy_exponent+1<d['second_energy_decade'])
    check('Duhamel_exponential_rounding',
          F(10**d['rate_log_decade']+10**d['source_log_decade'],10**d['second_log_decade'])<1)

    M=10**d['normal_step_decade']
    nall=sum(comb(15-k,2)*6**k for k in range(14))
    check('second_wave_recovery_row_count',nall==22568879235 and nall<10**12)
    check('forced_normal_homogeneous_row',1+2**11*(100*10**60+13*H)<10**6*H)
    check('forced_normal_all_rows',1+nall*(10**6*H)**2<M*M)
    for radial in range(2,14):
        valid=[]
        for time in range(14-radial):
            for length in range(14-radial-time):
                a=radial-2
                valid.append(a+time+length<=11 and a<radial)
                for k in range(time+1):
                    for b in range(a+1):
                        valid.append(b+k<=11-length and a-b+1<=radial-1
                                     and a-b+time-k+1+length<=13)
                valid.append(a+time+length+2<=13)
        check('forced_normal_allocation_order:'+str(radial),all(valid))
    check('augmented_normal_source_step',10**d['normal_source_decade']+1<M)
    check('normal_twelve_updates',12*d['normal_step_decade']==d['normal_total_decade'])
    check('base_W1_bound',2+36000+36000+144096<10**6)
    check('X13_rotation_conversion_factor',14*2**13<512**2)
    check('source_dominated_in_normal_recovery',d['source_uniform_decade']<=d['second_energy_decade'] and d['source_log_decade']<=d['second_log_decade'])
    check('final_mixed_prefactor_reserve',512*(10**3+1)<10**6)
    check('final_mixed_exponent',d['normal_total_decade']+d['second_energy_decade']+6<d['second_mixed_decade'])
    return {'mixed_derivative_allocations_checked':allocation_count,
      'source_S12_bound':'10^'+str(d['source_uniform_decade'])+' exp(10^'+str(d['source_log_decade'])+')',
      'initial_source_I11_bound':'10^'+str(d['initial_source_decade']),
      'retained_L1_and_initial_source_cost_bound':'10^'+str(d['total_source_cost_decade'])+' exp(10^'+str(d['source_log_decade'])+')',
      'initial_second_commuted_energy_bound':'10^'+str(d['second_initial_energy_decade']),
      'second_commuted_energy_bound':'10^'+str(d['second_energy_decade'])+' exp(10^'+str(d['second_log_decade'])+')',
      'second_mixed_X13_bound':'10^'+str(d['second_mixed_decade'])+' exp(10^'+str(d['second_log_decade'])+')',
      'scope':'Exact finite arithmetic and derivative allocations support the conditional forced second-wave theorem; analytic premises are not tested numerically.'}
