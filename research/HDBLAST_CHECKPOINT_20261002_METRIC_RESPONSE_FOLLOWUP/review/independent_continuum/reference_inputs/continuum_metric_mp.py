"""Genuine arbitrary-precision continuum route; import performs no evaluation.

The canonical source and all endpoint/constant terms are recomputed in each
registered mpmath context. No binary64 source jet enters these integrals.
"""
from __future__ import annotations
import math
import mpmath as mp
from source_jet import POLYNOMIALS

CONTRACT={'method':'mpmath_tanh_sinh_endpoint_subtraction_squared',
          'decimal_precisions':[50,70], 'maxdegree':10,
          'panels':['0','0.5','1'], 'reported_precision':70,
          'error_policy':'max_precision_gap_quad_estimates_and_float_conversion',
          'source_arithmetic':'analytic_integer_polynomials_mpmath_only'}


def require(ok,message):
    if not ok:raise RuntimeError(message)


def serialize_error_allowances(allowances):
    """Upward float serialization with the legacy validator's exact division.

    The inputs are empirical estimates, not certified integration errors.
    This helper evaluates no source or quadrature and accepts synthetic vectors.
    """
    denominator=8*math.pi**2
    raw=[];normalized=[]
    for allowance in allowances:
        require(mp.isfinite(allowance) and allowance>=0,'Invalid empirical error allowance')
        equivalent=float(allowance*mp.mpf(denominator))
        estimate=equivalent/denominator
        require(math.isfinite(equivalent) and math.isfinite(estimate),'Nonfinite serialized allowance')
        # At most a few adjacent floats are needed to avoid downward rounding
        # in either conversion or the final binary64 division.
        for _ in range(4):
            if mp.mpf(estimate)>=allowance:break
            equivalent=math.nextafter(equivalent,math.inf)
            estimate=equivalent/denominator
        require(mp.mpf(estimate)>=allowance,'Upward error serialization failed')
        raw.append(equivalent);normalized.append(estimate)
    return raw,normalized


def source_jet_mp(source,eta,guard):
    guard()
    require(source in ('positive_B','signed_uB'),'Unregistered high-precision source')
    u=eta+mp.mpf(4)
    if abs(u)>=1:return [mp.mpf(0)]*6
    denominator=(1-u)*(1+u)
    bump=mp.exp(-u*u/denominator)
    derivatives=[]
    for n,coefficients in enumerate(POLYNOMIALS):
        polynomial=mp.mpf(0)
        for coefficient in coefficients:polynomial=polynomial*u+coefficient
        derivatives.append(bump*polynomial/denominator**(2*n))
    if source=='positive_B':return derivatives
    return [u*derivatives[n]+(n*derivatives[n-1] if n else 0) for n in range(6)]


def canonical_forcing_jet_mp(jet,eta):
    L=-1/eta
    return [4*sum(math.comb(n,j)*math.factorial(n-j+1)*L**(n-j+2)*jet[j] for j in range(n+1))
            -2*sum(math.comb(n,j)*math.factorial(n-j)*L**(n-j+1)*jet[j+1] for j in range(n+1))-jet[n+2]
            for n in range(4)]


def forcing_jet_mp(source,eta,guard):
    return canonical_forcing_jet_mp(source_jet_mp(source,eta,guard),eta)


def log_history_mp(source,eta,order,endpoint_jet,config,guard):
    guard()
    require(order in (1,2,3),'Unregistered high-precision derivative order')
    if eta<=-5:return mp.mpf(0),mp.mpf(0)
    panels=[mp.mpf(x) for x in config['panels']]
    if eta<-3:
        ell=eta+5
        endpoint=endpoint_jet[order]
        logell=mp.log(ell)
        def integrand(z):
            guard()
            if z==0:return mp.mpf(0)
            difference=forcing_jet_mp(source,eta-ell*z*z,guard)[order]-endpoint
            return 2*ell*z*difference*(logell+2*mp.log(z))
        integral,error=mp.quad(integrand,panels,method='tanh-sinh',maxdegree=config['maxdegree'],error=True)
        integral+=ell*endpoint*(logell-1)
    else:
        def integrand(z):
            guard()
            # Flat compact-source endpoints, including eta=-3, are exactly0.
            if z==0 or z==1:return mp.mpf(0)
            t=-5+2*z
            return 2*forcing_jet_mp(source,t,guard)[order]*mp.log(eta-t)
        integral,error=mp.quad(integrand,panels,method='tanh-sinh',maxdegree=config['maxdegree'],error=True)
    constant=mp.log(mp.sqrt(2)/(-eta))+mp.euler+1
    guard()
    return integral+constant*endpoint_jet[order-1],error


def one_precision(source,eta_text,dps,config,guard):
    guard()
    with mp.workdps(dps):
        eta=mp.mpf(eta_text)
        endpoint=forcing_jet_mp(source,eta,guard)
        memories=[log_history_mp(source,eta,n,endpoint,config,guard) for n in (1,2,3)]
        L=-1/eta
        pref=1/(8*mp.pi**2)
        qjet=[-pref*memories[0][0],-pref*(memories[1][0]+L*endpoint[0]),-pref*(memories[2][0]+2*L*endpoint[1]+L*L*endpoint[0])]
        require(all(mp.isfinite(x) for x in qjet+endpoint+[entry for pair in memories for entry in pair]),'Nonfinite raw high-precision data')
        return {'dps':dps,'qjet':qjet,'quadrature_estimates':[pref*x[1] for x in memories],
                'memories':memories,'endpoint_forcing_jet':endpoint,
                'record':{'dps':dps,'qjet_decimal':[mp.nstr(x,dps) for x in qjet],
                          'quadrature_estimates_decimal':[mp.nstr(pref*x[1],dps) for x in memories],
                          'log_history_decimal':[mp.nstr(x[0],dps) for x in memories],
                          'raw_log_history_quad_errors_decimal':[mp.nstr(x[1],dps) for x in memories],
                          'endpoint_forcing_jet_decimal':[mp.nstr(x,dps) for x in endpoint]}}


def continuum_q_jet_mp(source,eta,config,guard):
    guard()
    require(config==CONTRACT,'High-precision continuum contract changed')
    runs=[one_precision(source,str(eta),dps,config,guard) for dps in config['decimal_precisions']]
    low,high=runs
    with mp.workdps(config['reported_precision']):
        values=[float(x) for x in high['qjet']]
        gaps=[abs(high['qjet'][n]-low['qjet'][n]) for n in range(3)]
        conversions=[abs(mp.mpf(values[n])-high['qjet'][n]) for n in range(3)]
        allowances=[max(gaps[n],conversions[n],low['quadrature_estimates'][n],high['quadrature_estimates'][n]) for n in range(3)]
        raw_allowances,errors=serialize_error_allowances(allowances)
        memories=[(float(x[0]),raw_allowances[n]) for n,x in enumerate(high['memories'])]
        record={'method':config['method'],'decimal_precisions':config['decimal_precisions'],
                'maxdegree':config['maxdegree'],'panels':config['panels'],
                'reported_precision':config['reported_precision'],'precision_runs':[x['record'] for x in runs],
                'qjet_precision_gaps_decimal':[mp.nstr(x,config['reported_precision']) for x in gaps],
                'qjet_float_conversion_errors_decimal':[mp.nstr(x,config['reported_precision']) for x in conversions],
                'effective_raw_memory_allowances_decimal':[mp.nstr(mp.mpf(x),config['reported_precision']) for x in raw_allowances],
                'error_policy':config['error_policy'],
                'saved_memory_error_scope':'Empirical raw-memory-equivalent allowance representing the complete qjet precision gap, both quadrature estimates and float conversion; true raw mp.quad estimates are preserved in precision_runs.',
                'uncertainty':'Precision differences, tanh-sinh quadrature estimates and final float conversion are empirical numerical diagnostics, not interval error certificates.'}
    guard()
    require(all(math.isfinite(x) for x in values+errors),'Nonfinite high-precision continuum output')
    return values,errors,memories,record
