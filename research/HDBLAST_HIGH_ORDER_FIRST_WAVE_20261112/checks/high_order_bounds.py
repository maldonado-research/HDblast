"""Exact finite constants and all derivative allocations; no PDE sampling."""
from fractions import Fraction as F
from math import comb,factorial


def evaluate(d,check):
    D=10**d['coefficient_decade'];H=10**d['normalized_coefficient_decade']
    check('raw_mixed_order_for_extension',13+2<=d['raw_coefficient_order'])
    check('connection_derivative_order',13+1<=16)
    check('order13_covariant_branch',5*(14+32*16)<10**4)
    check('order13_fixed_and_moving_conversion',150200+64*13+6+120+7<d['coefficient_decade'])
    check('time_frame_Cauchy_jets',8*factorial(19)*100**19<10**60)
    check('normalized_maps',2**13*10**210*D<H)
    check('squared_sech_time_jets',4*factorial(13)*100**13<10**60)
    check('linear_multiplier',2**24*13<10**10)
    R=10**d['initial_step_decade']
    for j in range(13):
        check('initial_step:'+str(j),1+2**j*(2*10**60+4*10**5*H)<R)
        m=12-j
        check('initial_all_derivative_allocations:'+str(j),all(
            k+b<=12 and m-b<=13-j+k and m-b+1<=14-j+k
            for k in range(j+1) for b in range(m+1)))
    check('initial_max_jet_exponent',13*d['initial_step_decade']+d['initial_input_decade']==d['initial_jet_decade'])
    check('reverse_rotation_rank_bound',2+12*16**2<60**2)
    check('initial_energy_finite_sum',315*60**26<10**60)
    check('sharp_initial_energy_multiplier',F(315,2)*60**26<10**50)
    check('sharp_initial_energy_exponent',d['initial_jet_decade']+25<1985300)
    check('initial_energy_rounding',d['initial_jet_decade']+30<d['initial_energy_decade'])
    nwords=sum((14-k)*6**k for k in range(14))
    check('commuted_word_count',nwords<10**13)
    check('commutator_one_word',2**13*(4000*10**60+538*D)<10**7*D)
    check('commutator_all_words',nwords*(10**7*D)**2<10**(2*d['commutator_decade']))
    for p in range(14):
        check('time_commutator_no_derivative_loss:'+str(p),all(p-k+i+1<=13
              for i in range(14-p) for k in range(1,p+1)))
    rate=F(10**d['base_K_decade'],2)+95*10**d['commutator_decade']
    check('energy_rate',rate<10**d['energy_rate_decade'])
    T=F(d['time']);check('actual_time',T==F(13,8000000))
    check('log_growth',T*10**(d['energy_rate_decade']-d['log_evolution_decade'])<1)
    check('angular_first_second_recovery',1+13**2+40**2<50**2)
    nall=sum(comb(16-k,2)*6**k for k in range(15))
    check('all_recovery_rows',nall==135413275530 and nall<10**12)
    check('normal_single_row',1+2**12*(100*10**60+13*H)<10**6*H)
    check('normal_all_rows',1+10**12*(10**6*H)**2<10**(2*d['normal_step_decade']))
    for r in range(2,15):
        allocations=[]
        for p in range(15-r):
            for length in range(15-r-p):
                a=r-2
                allocations.append(a+p+2+length<=14 and a<r)
                for k in range(p+1):
                    allocations.append(a+p-k+length+2<=14)
                    for b in range(a+1):
                        allocations.append(b+k<=12-length and a-b+1<=r-1
                            and a-b+p-k+1+length<=14)
        check('normal_allocation_order:'+str(r),all(allocations))
    check('base_W1',2+36000+36000+144096==216098 and 216098<10**6)
    for j in range(14):check('angular_covariant_slot:'+str(j),6*(j+2)**2<39**2)
    check('final_spatial_binomial_factor',4*sum(1601**m for m in range(15))<10**60)
    check('sharp_full_spatial_factor',2**13<128**2 and 128+128==256)
    check('sharp_recovery_exponent',13*d['normal_step_decade']+3==1980163)
    check('sharp_N14_prefactor',1985300+1980163+3==3965466
          and 3965466<d['final_N14_decade'])
    check('recovery_exponent',30+13*d['normal_step_decade']+3==d['normal_recovery_decade'])
    check('N14_rounding',d['normal_recovery_decade']+d['initial_energy_decade']<d['final_N14_decade'])
    check('forcing_prefactor_rounding',d['prior_forcing_decade']+2*d['final_N14_decade']<d['final_forcing_decade'])
    check('forcing_exponential_rounding',2*10**d['log_evolution_decade']<10**d['final_forcing_log_decade'])
    check('quadrupole_angular_H14',sum(8**j for j in range(15))<3000000**2)
    check('quadrupole_radial_to_full_H14',sum(9**j for j in range(15))==25736391511831
          and sum(9**j for j in range(15))<10**14)
    return {'initial_temporal_jet_bound':'10^'+str(d['initial_jet_decade']),
      'initial_commuted_energy_bound':'10^'+str(d['initial_energy_decade']),
      'normal_recovery_multiplier':'10^'+str(d['normal_recovery_decade']),
      'first_wave_N14_bound':'10^'+str(d['final_N14_decade'])+' exp(10^'+str(d['log_evolution_decade'])+')',
      'whole_slice_H12_forcing_bound':'10^'+str(d['final_forcing_decade'])+' exp(10^'+str(d['final_forcing_log_decade'])+')',
      'scope':'Finite arithmetic/allocation checks support the written conditional energy, equivariance and Sobolev proofs.'}
