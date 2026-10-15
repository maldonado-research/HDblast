#!/usr/bin/env python3
"""Pin-specific Astra arithmetic/consistency cross-review; no producer import.

The continuous cell-range formulas and characteristic proof were reviewed
analytically; this script independently checks export and norm arithmetic.
"""
import csv
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
EXPECTED={
    'UNIFORM_FIRST_VARIATION_WAVEFORM.md':'6eada42c742aec6141925c244df21dc263e3555ae2a064e1e5d3fb7dcc9a82ad',
    'waveform_intervals.py':'f17f3a54d7aa34a0012b52eaccefb7e1e7c319c32ef51ed1286e0c53954f5593',
    'compute_waveform.py':'0eaa55feaf321c6ad57ef3334c9b351a18daba25f4ec1b737d7aaa19b61fa1b4',
    'check_waveform.py':'78d0789a84da840d5a7352d96df5c4c6c61787eca060bd33ea707ad8b531138e',
}


def need(condition,message):
    if not condition:
        raise RuntimeError(message)


def ceildiv(n,d):
    return -((-n)//d)


def main():
    checks=[]
    def check(name,condition):
        need(condition,name)
        checks.append(name)
    for name,digest in EXPECTED.items():
        check('expected_pin:'+name,sha256((ROOT/name).read_bytes()).hexdigest()==digest)

    report=json.loads((ROOT/'computed/WAVEFORM_ENCLOSURE.json').read_text())
    rows=list(csv.DictReader((ROOT/'computed/WAVEFORM_CELLS.csv').open()))
    S=report['cell_data']['value_denominator']
    N=report['cell_data']['time_grid_denominator']
    check('fixed_arithmetic_and_grid',S==2**96 and N==4096)
    check('csv_rows_and_hash',len(rows)==4098 and
          sha256((ROOT/'computed/WAVEFORM_CELLS.csv').read_bytes()).hexdigest()
          ==report['cell_data']['sha256'])
    for name,digest in report['source_pins'].items():
        check('producer_source_pin:'+name,sha256((ROOT/name).read_bytes()).hexdigest()==digest)
    check('gain_input_pin',sha256((ROOT/report['gain_input']['path']).read_bytes()).hexdigest()
          ==report['gain_input']['sha256'])

    values={k:F(v) for k,v in report['error_constants'].items()}
    ell=F(1,10**10);T=3*ell;m=F(1,1000);M=10**6;B=1000;gamma=250
    U=1+T;r=M*T*T/2;s=B*T*T/(2*m);den=(1-r)*(1-s)-T*r/m
    A=U*(1-s+T/m)/den;V=U/(m*den)
    eta=ell*(M*T*A+(1+gamma*T)*V)
    Vtail=ell/m**2+ell/m+r*A/m+s*V
    tail=ell*(M*T*A+Vtail+gamma*T*V)
    advertised_eta=F(10000004,10**14);advertised_tail=F(4752,10**17)
    zeta=F(10000005,10**14)
    check('positive_majorant',r<1 and s<1 and den>0)
    check('majorant_exact_values',all(values[k]==v for k,v in {
        'ell':ell,'T':T,'r':r,'s':s,'determinant':den,'A_upper':A,'V_upper':V,
        'eta_upper':eta,'published_eta':advertised_eta,'eta_tail_upper':tail,
        'published_tail':advertised_tail,'norm_error_upper':zeta}.items()))
    check('strict_error_theorems',eta<advertised_eta and tail<advertised_tail
          and eta**2+tail**2<zeta**2)

    def dyadic(a,b=None):
        a=F(a);b=a if b is None else F(b)
        return (a.numerator*S//a.denominator,ceildiv(b.numerator*S,b.denominator))
    def times(a,b):
        ends=[x*y for x in a for y in b]
        return min(ends)//S,ceildiv(max(ends),S)
    def sqr(a):
        low=0 if a[0]<=0<=a[1] else min(a[0]**2,a[1]**2)
        return low//S,ceildiv(max(a[0]**2,a[1]**2),S)
    def root(a):
        k=isqrt(a.numerator*S*S//a.denominator)
        return F(k,S),F(k,S) if F(k,S)**2==a else F(k+1,S)

    gain=tuple(map(F,report['physical_gain_interval']))
    K=tuple(map(F,report['physical_norm_prefactor_interval']))
    check('signed_gain_and_positive_norm_prefactor',gain[0]<gain[1]<0<K[0]<K[1])
    err=dyadic(-advertised_eta,advertised_eta);tailiv=dyadic(-advertised_tail,advertised_tail)
    gainiv=dyadic(*gain)
    grid=[(0,N)]+[(N+i,N+i+1) for i in range(N)]+[(2*N,3*N)]
    sums=[0,0]
    for i,(row,pair) in enumerate(zip(rows,grid)):
        a,b=int(row['left_grid']),int(row['right_grid'])
        need((a,b)==pair,'cell grid or coverage mismatch at '+str(i))
        comp=(int(row['comparison_lo']),int(row['comparison_hi']))
        actual=(int(row['H_lo']),int(row['H_hi']))
        physical=(int(row['ell_S_lo']),int(row['ell_S_hi']))
        need(comp[0]<=comp[1] and actual[0]<=actual[1] and physical[0]<=physical[1],
             'reversed cell at '+str(i))
        if i==0:
            need(comp==(0,0) and actual==(0,0),'prefront cell is nonzero')
        elif i==len(rows)-1:
            need(comp==(0,0) and actual==tailiv,'tail cell mismatch')
        else:
            need(actual==(comp[0]+err[0],comp[1]+err[1]),'pointwise error addition mismatch')
            squared=sqr(comp);sums[0]+=squared[0];sums[1]+=squared[1]
        factor=dyadic(1-(ell*F(b,N))**2,1)
        need(physical==times(gainiv,times(factor,actual)),
             'signed physical reconstruction mismatch at '+str(i))
    checks.extend(['all_4098_cells_cover_domain','all_cell_error_additions',
                   'all_signed_physical_cell_multiplications','prefront_and_tail_cells'])
    comparison_integral=tuple(F(v,N*S) for v in sums)
    quadrature=report['quadrature']
    check('comparison_integral_from_exported_cells',
          comparison_integral==tuple(map(F,quadrature['comparison_integral_interval'])))
    norms=(max(F(0),root(comparison_integral[0])[0]-zeta),
           root(comparison_integral[1])[1]+zeta)
    actual_integral=tuple(v*v for v in norms)
    check('actual_norm_triangle_inequality_arithmetic',
          actual_integral==tuple(map(F,quadrature['actual_H_integral_interval'])))
    C=(K[0]*actual_integral[0]/ell,K[1]*actual_integral[1]/ell)
    check('proper_time_jacobian_and_prefactor',C==tuple(map(F,quadrature['C_interval'])))
    check('advertised_response_interval',F(606,100)*10**15<C[0]<C[1]<F(744,100)*10**15)
    rising=report['sign_witnesses']['rising'];falling=report['sign_witnesses']['falling']
    check('normalized_trace_sign_witnesses',F(rising['H_interval'][0])>0
          and F(falling['H_interval'][1])<0)
    check('physical_sign_reversal',F(rising['ell_S_interval'][1])<0
          and F(falling['ell_S_interval'][0])>0)

    replay_normal=(HERE/'ROOT_WAVEFORM_CHECK_NORMAL.json').read_bytes()
    replay_optimized=(HERE/'ROOT_WAVEFORM_CHECK_OPTIMIZED.json').read_bytes()
    replay=json.loads(replay_normal)
    check('root_checker_normal_optimized_identical',replay_normal==replay_optimized)
    check('root_checker_pass',replay['status']=='PASS' and replay['finite_check_groups']>=40)
    peak=json.loads((HERE/'peak_analysis/certificate.json').read_text())
    check('peak_actual_error_not_overclaimed',not peak['actual_waveform_corollary']['critical_point_uniqueness_certified']
          and not peak['actual_waveform_corollary']['actual_zero_uniqueness_certified']
          and not peak['actual_waveform_corollary']['absolute_maximum_lobe_identified'])
    peakmax=tuple(map(F,peak['enclosures']['Actual waveform maximum height, conditional']))
    peakmin=tuple(map(F,peak['enclosures']['Actual waveform minimum height, conditional']))
    check('main_rounded_peak_heights',F('4.2356403')<peakmax[0]<peakmax[1]<F('4.2356406')
          and F('-4.2356406')<peakmin[0]<peakmin[1]<F('-4.2356403'))

    files=list(EXPECTED)+[
        'plot_waveform.py','computed/WAVEFORM_CELLS.csv','computed/WAVEFORM_ENCLOSURE.json',
        'figures/SHELL_FIRST_VARIATION_WAVEFORM.png','figures/SHELL_FIRST_VARIATION_WAVEFORM.svg',
        'archive/PHYSICAL_GAIN_AUDIT.md','archive/PHYSICAL_GAIN_ENCLOSURE.json',
        'sol/UNIFORM_WAVEFORM_CERTIFICATE.md','astra/UNIFORM_POINTWISE_WAVEFORM.md',
        'astra/UNIFORM_POINTWISE_NORMAL.json','astra/UNIFORM_POINTWISE_OPTIMIZED.json',
        'astra/peak_analysis/peak_analysis.txt','astra/peak_analysis/verify_peaks.py',
        'astra/peak_analysis/certificate.json','astra/peak_analysis/certificate_optimized.json',
        'astra/ROOT_WAVEFORM_CHECK_NORMAL.json','astra/ROOT_WAVEFORM_CHECK_OPTIMIZED.json',
        'astra/review_root_waveform.py','astra/FINAL_ROOT_WAVEFORM_REVIEW.md',
    ]
    result={
        'status':'CONDITIONAL_ACCEPT_NO_BLOCKER',
        'reviewer':'Astra independent analytic derivation and final cross-review',
        'main_sha256':EXPECTED['UNIFORM_FIRST_VARIATION_WAVEFORM.md'],
        'checks':checks,
        'root_checker_finite_groups':replay['finite_check_groups'],
        'all_exported_cells_checked':len(rows),
        'visual_review':'PNG viewed: physical sign reversal, envelope, illustrative centerline, model-time axes, unresolved tail and retained-hypothesis qualifiers are consistent.',
        'analytic_review':[
            'positive boundary damping and integration-by-parts cancellation',
            'continuous displacement/velocity majorant and incoming derivative trace',
            'finite-propagation prefront and history-preserving tail bounds',
            'support-aware observation norm and exact K/ell scaling',
            'negative correlated physical gain and cosh conjugation',
            'continuous-cell interval formulas, support-edge limits, exponential bounds',
            'actual peak value-gap exclusion with no actual uniqueness or lobe-selection claim',
            'retained exact-family, gauge, nonlinear causal-margin, and scale limitations',
        ],
        'scope':'Acceptance of the specified conditional first-variation result and finite arithmetic; not formal verification of the PDE, nonlinear hypotheses, or independent external peer review.',
        'blocking_findings':[],
        'reviewed_files':{name:sha256((ROOT/name).read_bytes()).hexdigest() for name in files},
    }
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__=='__main__':
    main()
