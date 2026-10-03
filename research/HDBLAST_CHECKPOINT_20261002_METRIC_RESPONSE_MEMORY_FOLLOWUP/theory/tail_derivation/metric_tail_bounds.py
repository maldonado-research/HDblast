"""Deferred constructive omitted-UV-band bounds for homogeneous metric response.

Import performs no physical source, response, mode or interval evaluation.
Call only after all sources and inputs are publicly frozen. Global rational
derivative caps are conservative; they do not give precise continuum stress.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import inf,nextafter
from pathlib import Path
import json

INTERVAL_DECIMAL_PRECISION=60
OBSERVATIONS=tuple(F(x) for x in ('-5.5','-4.5','-4','-3.5','-2.5','-1.5'))

@lru_cache(maxsize=1)
def _inputs():
    root=Path(__file__).resolve().parent
    cap=json.loads((root/'DERIVATIVE_CERTIFICATE.json').read_text())
    local=json.loads((root/'METRIC_TAIL_INPUTS.json').read_text())
    if cap['status']!='PASS' or local['status']!='PASS':raise RuntimeError('Analytic input certificate failed')
    return cap,local

def _rat(iv,x):
    q=F(str(x));return iv.mpf(q.numerator)/iv.mpf(q.denominator)

def _upper(x):
    if x==0:return 0.0
    return nextafter(float(x.b),inf)

def _moment_contact_upper(iv,rows,L,h,K):
    ans=iv.mpf(0)
    for power,terms in rows.items():
        n=int(power)
        if n<5:raise RuntimeError('Omitted contact moment is not integrable')
        coefficient=sum((_rat(iv,abs(F(t['coefficient'])))*L**t['L_power']*h[t['h_index']] for t in terms),iv.mpf(0))
        ans+=coefficient*K**(3-n)/(n-3)
    return ans/(2*iv.pi**2)

def _ratio_upper(iv,spec,L,v,one_minus_v,h):
    numerator=iv.mpf(0)
    for term in spec['numerator']:
        value=_rat(iv,abs(F(term['coefficient'])))*L**term['L_power']*v**term['v_power']
        if term['h_index'] is not None:value*=h[term['h_index']]
        numerator+=value
    denominator=iv.mpf(0)
    for coefficient in spec['denominator_coefficients']:denominator=denominator*v+_rat(iv,coefficient)
    if not denominator.a>0:raise RuntimeError('Baseline denominator is not positive')
    return one_minus_v**spec['one_minus_v_power']*numerator/(denominator*iv.pi**2)

def metric_tail_bounds(source,eta,cutoff):
    """Upper bounds for the removed k>K band in the frozen metric calibration.

    q=a² deltaQ/epsilon, with conformal q derivatives; rho,p=a⁴ delta(rho,p)
    /epsilon; current=a² deltaj/epsilon=q for b=1,r=2,delta phi=0.
    Q0,rho0,p0 bound physical H=1 baseline differences. They are not divided
    by epsilon. Only the omitted UV band is enclosed; all finite integration
    and arithmetic errors require separate evidence.
    """
    if source not in ('positive_B','signed_uB'):raise ValueError('Unknown metric pulse')
    time=F(str(eta));k=F(str(cutoff))
    if time not in OBSERVATIONS:raise ValueError('Observation outside frozen six-point domain')
    if k<=0:raise ValueError('Cutoff must be positive')
    from mpmath.ctx_iv import MPIntervalContext
    iv=MPIntervalContext();iv.dps=INTERVAL_DECIMAL_PRECISION
    cert,local=_inputs();caps=cert['envelopes'][source]
    K=_rat(iv,k);L=-1/_rat(iv,time);a=L;M2=2*L*L
    radial=iv.sqrt(K*K+M2);v=K/radial
    om=M2/(radial*(radial+K))
    om3=om*(1+v+v*v)
    zero=iv.mpf(0)
    local_nonzero=(-5<time<-3)
    h=[_rat(iv,row['sup_upper_integer']) if local_nonzero else zero for row in caps['h_derivative_caps']]
    g=[_rat(iv,q) if local_nonzero else zero for q in caps['g_sup_upper_rational']]
    N=[_rat(iv,q) if time>-5 else zero for q in caps['g_N_upper_rational']]
    basic=[(3*g[j]*M2/2+N[j]/4)/(16*iv.pi**2*K*K) for j in range(3)]
    C=1/(8*iv.pi**2)
    qproxy=[basic[0],basic[1]+C*L*g[0]*om3,
            basic[2]+C*(2*L*g[1]*om3+L*L*g[0]*om*(3*v**4+3*v**3+v*v+v+1))]
    base={name:_ratio_upper(iv,spec,L,v,om,h) for name,spec in local['baseline_specs'].items()}
    density_contact=L*(3*L*g[0]+g[1])*om3/(96*iv.pi**2)
    J5=om3/6
    J7=om**2*(3*v**3+6*v*v+4*v+2)/60
    J9=om**3*(15*v**4+45*v**3+48*v*v+24*v+8)/840
    pressure_contact=(70*L*L*g[0]*J9+(30*L*L*g[0]+10*L*g[1])*J7+(g[2]+L*g[1]+9*L*L*g[0])*J5)/(48*iv.pi**2)
    rho_proxy=(3*L*L*qproxy[0]+L*qproxy[1])/2+a*a*g[0]*base['Q0']/2+density_contact
    p_proxy=(qproxy[2]+3*L*qproxy[1]+3*L*L*qproxy[0])/6+a*a*g[0]*base['Q0']/6+pressure_contact
    contact={name:_moment_contact_upper(iv,rows,L,h,K) for name,rows in local['contact_inventory'].items()}
    q_reference=[_ratio_upper(iv,spec,L,v,om,h) for spec in local['q_reference_derivative_specs']]
    qlocal=[contact[name]+q_reference[j] for j,name in enumerate(('q','q_prime','q_second'))]
    rho_local=contact['rho']+4*h[0]*L**4*base['rho0']+(h[2]+4*L*h[1])*L*L*base['Q0']/2
    p_local=contact['p']+4*h[0]*L**4*base['p0']+h[2]*L*L*base['Q0']/2
    bound={'q':qproxy[0]+qlocal[0],'q_prime':qproxy[1]+qlocal[1],'q_second':qproxy[2]+qlocal[2],
           'rho':rho_proxy+rho_local,'p':p_proxy+p_local,
           'current':qproxy[0]+qlocal[0],**base}
    result={name:_upper(value) for name,value in bound.items()}
    result.update({'source':source,'eta':str(eta),'K':str(cutoff),'interval_dps':INTERVAL_DECIMAL_PRECISION,
      'normalization':{'q':'a²deltaQ/epsilon; conformal derivatives for q_prime,q_second','rho,p':'a⁴delta(rho,p)/epsilon','current':'a²deltaj/epsilon=q; fixed b1,r2,delta_phi0','Q0,rho0,p0':'Physical H1 baseline difference, not divided by epsilon'},
      'canonical_proxy_tail_upper':{name:_upper(qproxy[j]) for j,name in enumerate(('q','q_prime','q_second'))},
      'proxy_stress_tail_upper':{'rho':_upper(rho_proxy),'p':_upper(p_proxy)},
      'metric_contact_tail_upper':{'q':_upper(qlocal[0]),'q_prime':_upper(qlocal[1]),'q_second':_upper(qlocal[2]),'rho':_upper(rho_local),'p':_upper(p_local)},
      'global_g_N_upper_rational':caps['g_N_upper_rational'],
      'global_g_sup_upper_rational':caps['g_sup_upper_rational'],
      'local_source_jets_exact_zero':not local_nonzero,
      'source_history_exact_zero':time<=-5,
      'method':'Exhaustive rational bump extrema certificates; exact Leibniz global derivative norms; 60-digit directed interval propagation; integrable rational metric-contact moments and factored baseline tails.',
      'uncertainty_scope':'Omitted k>K band only. Conservative global bounds are loose and do not certify precise continuum pressure or sign. Finite-mode/time/momentum quadrature and response arithmetic errors are separate.'})
    return result
