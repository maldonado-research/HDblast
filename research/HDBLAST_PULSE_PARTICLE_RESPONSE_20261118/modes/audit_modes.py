"""Standard-library reconstruction of the Gaussian midpoint and numerical summary."""
import argparse
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent


def covariance(mode):
    v,p=complex(*mode[:2]),complex(*mode[2:])
    return abs(v)**2,(v.conjugate()*p).real,abs(p)**2


def det(m):
    return m[0]*m[2]-m[1]**2


def dot_symmetric(m,n):
    return m[0]*n[0]+2*m[1]*n[1]+m[2]*n[2]


def covariance_checks(results):
    discrepancy=purity=matched=0.
    count=0
    for grid in results.values():
        for cut in grid.values():
            for seed in ('instant','wkb'):
                for pair in cut[seed]['selected'].values():
                    cp=covariance(pair['pulse_mode'])
                    cq=covariance(pair['quiet_mode'])
                    total=[a+b for a,b in zip(cp,cq)]
                    purity=max(purity,abs(det(cp)-.25),abs(det(cq)-.25))
                    d=det(total)
                    discrepancy=max(discrepancy,abs((d-1)-pair['relative']))
                    common=[total[2]/math.sqrt(d),-total[1]/math.sqrt(d),total[0]/math.sqrt(d)]
                    for cc in (cp,cq):
                        occupation=(dot_symmetric(common,cc)-1)/2
                        matched=max(matched,abs(occupation-pair['minimax_common_occupation']))
                    count+=1
    if max(discrepancy,purity,matched)>1e-11:
        raise ValueError('Gaussian covariance benchmark failed')
    return dict(pairs=count,max_purity_defect=purity,max_mismatch_determinant_error=discrepancy,
                max_equal_midpoint_occupation_error=matched)


def summary(report):
    results=report['results']
    out={'covariance_checks':covariance_checks(results)}
    rows=[]
    for path in ('g1/plus','g2/plus','g4/plus','g8/plus','g8/minus',
                 'g8/half_plus','g8/half_minus','g8/without_remote_controls','g8/midpoint_plus'):
        cut=results[path]['0.15']
        instant=cut['instant']
        mode=instant['selected']['10.0']
        band=instant['band']['5-18/32']
        rows.append(dict(case=path,equal_proper_time_T=cut['T'],k10_mismatch=mode['relative'],
                      k10_common_floor=mode['minimax_common_occupation'],k10_instant_diff=mode['delta_n'],
                      k10_wkb_end_diff=mode['delta_wkb_n'],band_mismatch=band['relative'],
                      band_floor=band['minimax_common_occupation'],delta_band_variance=band['delta_field_variance'],
                      band_excitation_energy_difference=band['delta_energy_excitation'],
                      band_zero_point_difference=band['delta_energy_zero_point'],
                      band_total_energy_difference=band['delta_energy_total']))
    out['fixed_clock']=rows
    p,m,hp,hm=[results['g8/'+x]['0.15']['instant']['band']['5-18/32'] for x in ('plus','minus','half_plus','half_minus')]
    def ratio(field,sign):
        return (p[field]+sign*m[field])/(hp[field]+sign*hm[field])
    out['amplitude_control']={'even_relative_ratio':ratio('relative',1),
                              'odd_variance_ratio':ratio('delta_field_variance',-1),
                              'even_variance_ratio':ratio('delta_field_variance',1)}
    cut=results['g8/plus']['0.15']
    band=cut['instant']['band']['5-18/32']
    wkb=cut['wkb']['band']['5-18/32']
    out['comparison']={'prearrival_k10_mismatch':results['g8/plus']['0.09']['instant']['selected']['10.0']['relative'],
               'early_k10_mismatch':results['g8/plus']['0.12']['instant']['selected']['10.0']['relative'],
               'quiet_variance':band['field_variance_quiet'],
               'fractional_variance_change':band['delta_field_variance']/band['field_variance_quiet'],
               'alternative_field_momentum_map_rel':band['relative_field'],
               'wkb_initial_band_mismatch':wkb['relative'],'wkb_initial_variance_diff':wkb['delta_field_variance'],
               'max_interpolation_sample_screened_frequency':cut['max_screened_physical_frequency_at_sampled_nodes'],
               'max_wronskian_error':max(c[s]['max_wronskian_error'] for v in results.values()
                                          for c in v.values() for s in ('instant','wkb')),
               'fine_grid_relative_change_band_mismatch':
                    abs(band['relative']/results['g4/plus']['0.15']['instant']['band']['5-18/32']['relative']-1),
               'alternate_stepper_relative_change_band_mismatch':
                    abs(band['relative']/results['g8/midpoint_plus']['0.15']['instant']['band']['5-18/32']['relative']-1)}
    out['quadrature_comparison']={key:abs(value-cut['instant']['band']['5-18/16'][key])
                                  for key,value in band.items()}
    if out['comparison']['max_wronskian_error']>1e-10 or band['relative']<=0:
        raise ValueError('Mode or norm checks failed')
    return out


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report',nargs='?',type=Path,default=HERE/'MODES.json')
    print(json.dumps(summary(json.loads(parser.parse_args().report.read_text())),indent=2,allow_nan=False))
