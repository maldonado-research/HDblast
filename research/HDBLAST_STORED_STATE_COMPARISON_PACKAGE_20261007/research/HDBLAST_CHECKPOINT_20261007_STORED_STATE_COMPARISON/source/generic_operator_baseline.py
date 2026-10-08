#!/usr/bin/env python3
"""Exact-rational polynomial-defect operator; no physical source callback.

This file only executes fabricated fixtures. Its generic functions require a
separately proved source polynomial enclosure. It does not establish that
input enclosure by sampling or trust coefficient agreement as a proof.
"""
from fractions import Fraction as Q
from math import comb, factorial
import hashlib
import json
from pathlib import Path
import resource
import sys
import time


ZERO = (Q(0), Q(0))


def add(a, b):
    return (a[0]+b[0], a[1]+b[1])


def scale(a, q):
    return (a[0]*q, a[1]*q)


def iomega(a, omega):
    return (-omega*a[1], omega*a[0])


def norm_upper(a):
    return abs(a[0])+abs(a[1])


def evaluate(p, x):
    out = ZERO
    for value in reversed(p):
        out = add(scale(out, x), value)
    return out


def centered_to_left(coefficients, half_width):
    """Exact translation p(x-half_width); no source evaluation."""
    out = [Q(0)]*len(coefficients)
    for n, a in enumerate(coefficients):
        for j in range(n+1):
            out[j] += a*comb(n,j)*(-half_width)**(n-j)
    return tuple(out)


def exact_integral(coefficients, width):
    return sum((a*width**(n+1)/Q(n+1)
                for n,a in enumerate(coefficients)), Q(0))


def chosen_dyadic(value, bits):
    """Choose one exact dyadic point; its error is certified by the residual."""
    unit=1<<bits
    n=(value.numerator*unit)//value.denominator
    return Q(n,unit)


def chosen_complex(a,bits):
    return (chosen_dyadic(a[0],bits),chosen_dyadic(a[1],bits))


def outward_radius(value,bits=512):
    """Exact upward dyadic rounding prevents prefix-denominator growth."""
    if not isinstance(value,Q) or value<0:
        raise ValueError('Nonnegative exact radius required')
    unit=1<<bits
    return Q((value.numerator*unit+value.denominator-1)//value.denominator,unit)


def defect_certificate(g, source_error, width, omega, degree,
                       solver_omega=None, coefficient_bits=None,
                       incoming_w=ZERO, incoming_u=ZERO,
                       incoming_w_error=Q(0),incoming_u_error=Q(0)):
    """Certify zero-input source operators, conditional on |g_real-g|≤error.

    w'=i*omega*w-g_real, u'=w. At zero incoming modes -w(width)=Mexp,
    -u(width)=integral phi1(i*omega*(width-s))*(width-s)*g_real(s) ds.
    A nonzero incoming state instead encloses negative endpoint state; its
    source-operator interpretation needs a documented prefix that began at
    zero, as in the fabricated whole-interval construction below.
    No division by omega is used. source_error is a uniform absolute bound.
    """
    if not isinstance(source_error,Q) or source_error<0:
        raise ValueError('Nonnegative exact source enclosure required')
    if not isinstance(width,Q) or width<=0:
        raise ValueError('Positive exact width required')
    if not isinstance(omega,Q) or degree<0:
        raise ValueError('Exact real omega/nonnegative degree required')
    if any(not isinstance(a,Q) for a in g):
        raise ValueError('Only exact chosen real polynomial coefficients accepted')
    sw = omega if solver_omega is None else solver_omega
    if not isinstance(sw,Q):
        raise ValueError('Exact real solver omega required')
    if min(incoming_w_error,incoming_u_error)<0:
        raise ValueError('Nonnegative inherited error required')
    w = [incoming_w]
    for n in range(degree):
        gn = g[n] if n<len(g) else Q(0)
        wn=scale(add(iomega(w[n], sw),(-gn,Q(0))),Q(1,n+1))
        w.append(wn if coefficient_bits is None else chosen_complex(wn,coefficient_bits))
    # U'=W remains an exact identity; rational divisions are never rounded here.
    u = [incoming_u]+[scale(a,Q(1,n+1)) for n,a in enumerate(w)]
    # Authenticate the actual chosen polynomial against the true equation.
    residue=[]
    for n in range(max(len(w),len(g))):
        derivative=scale(w[n+1],Q(n+1)) if n+1<len(w) else ZERO
        frequency=iomega(w[n],omega) if n<len(w) else ZERO
        gn=g[n] if n<len(g) else Q(0)
        residue.append(add(add(derivative,scale(frequency,Q(-1))),(gn,Q(0))))
    bounds=[norm_upper(a) for a in residue]
    rw = incoming_w_error+source_error*width
    ru = incoming_u_error+width*incoming_w_error+source_error*width**2/2
    for n,qn in enumerate(bounds):
        rw += qn*width**(n+1)/Q(n+1)
        ru += qn*width**(n+2)/Q((n+1)*(n+2))
    me = scale(evaluate(w,width),Q(-1))
    md = scale(evaluate(u,width),Q(-1))
    m0 = exact_integral(g,width)
    return {
        'source_error':source_error,'width':width,'omega':omega,
        'degree':degree,'coefficient_bits':coefficient_bits,'M0':(m0,source_error*width),
        'incoming_state_scope':('zero_anchored_source_operator'
            if incoming_w==ZERO and incoming_u==ZERO and incoming_w_error==0 and incoming_u_error==0
            else 'propagated_prefix_or_general_state_must_declare_origin'),
        'Mexp':(me,rw),'Duhamel_u':(md,ru),
        'residual_coefficients':tuple(residue),
        'chosen_w_coefficients':tuple(w),
        'chosen_u_coefficients':tuple(u),
    }


def rational(value):
    return {'numerator':str(value.numerator),'denominator':str(value.denominator)}


def rectangle(center, radius):
    return {'real_lower':rational(center[0]-radius),
            'real_upper':rational(center[0]+radius),
            'imag_lower':rational(center[1]-radius),
            'imag_upper':rational(center[1]+radius),
            'complex_norm_error_upper':rational(radius)}


def output_moments(receipt):
    m0,r0=receipt['M0']
    return {'M0':{'lower':rational(m0-r0),'upper':rational(m0+r0)},
            'Mexp':rectangle(*receipt['Mexp']),
            'Duhamel_u':rectangle(*receipt['Duhamel_u'])}


def read_q(obj):
    if set(obj)!= {'numerator','denominator'}:
        raise ValueError('Exact endpoint encoding required')
    if not isinstance(obj['numerator'],str) or not isinstance(obj['denominator'],str):
        raise ValueError('Endpoint integer strings required')
    den=int(obj['denominator'])
    if den<=0:
        raise ValueError('Positive denominator required')
    return Q(int(obj['numerator']),den)


def validate_serialized(receipt, encoded):
    expected=output_moments(receipt)
    # Contract is exact here; a later outward decimal rendition needs its own guard.
    if encoded!=expected:
        raise ValueError('Altered or inward/incorrect serialized enclosure')
    for name in ('Mexp','Duhamel_u'):
        got=encoded[name]
        if read_q(got['real_lower'])>read_q(got['real_upper']):
            raise ValueError('Reversed interval')
    return True


def constant_source_truth(omega,width,terms=200):
    """Different exact entire-kernel series, for fabricated g=1 only.

    ME=sum (iω)^n H^(n+1)/(n+1)!;
    Du=sum (iω)^n H^(n+2)/(n+2)!.
    Real-phase integral Taylor remainders bound these directly, including ω=0.
    """
    power=(Q(1),Q(0)); me=ZERO; du=ZERO
    for n in range(terms+1):
        me=add(me,scale(power,width**(n+1)/factorial(n+1)))
        du=add(du,scale(power,width**(n+2)/factorial(n+2)))
        power=iomega(power,omega)
    x=abs(omega*width)
    re=width*x**(terms+1)/factorial(terms+2)
    ru=width**2*x**(terms+1)/factorial(terms+3)
    return (me,re),(du,ru)


def fabricated_tests():
    checks=[]; mutations=[]
    def require(ok,name):
        if not ok: raise RuntimeError(name)
        checks.append(name)
    h=Q(1,64)
    for omega in (Q(0),Q(1,2**39),Q(1,2**11),Q(1,2),Q(2),Q(32),Q(128),Q(256),Q(512)):
        r=defect_certificate((Q(1),),Q(0),h,omega,120)
        truth=constant_source_truth(omega,h)
        for field,(tc,tr) in zip(('Mexp','Duhamel_u'),truth):
            center,rr=r[field]
            require(max(abs(center[0]-tc[0]),abs(center[1]-tc[1]))+tr<=rr,
                    f'constant forcing independently enclosed {field}, omega={omega}')
        require(r['M0']==(h,Q(0)),f'exact M0, omega={omega}')
        enc=json.loads(json.dumps(output_moments(r)))
        require(validate_serialized(r,enc),f'exact serialization, omega={omega}')
    polynomial=(Q(2),Q(-3),Q(5,7),Q(-1,11))
    zero=defect_certificate(polynomial,Q(0),h,Q(0),120)
    require(zero['Mexp']==((exact_integral(polynomial,h),Q(0)),Q(0)),
            'zero phase polynomial ME exact')
    expected_u=sum((a*h**(n+2)/Q((n+1)*(n+2))
                    for n,a in enumerate(polynomial)),Q(0))
    require(zero['Duhamel_u']==((expected_u,Q(0)),Q(0)),
            'zero phase polynomial nested Duhamel exact')
    # Source uncertainty is never lost in a phase/moment reduction.
    uncertain=defect_certificate(polynomial,Q(1,10**9),h,Q(128),120)
    require(uncertain['Mexp'][1]>=h/Q(10**9),'source error propagated to ME')
    require(uncertain['Duhamel_u'][1]>=h*h/Q(2*10**9),'source error propagated to nested moment')
    # A wrong phase is tested against the declared equation, not its own recurrence.
    wrong=defect_certificate((Q(1),),Q(0),h,Q(512),120,solver_omega=Q(-512))
    require(wrong['Mexp'][1]>Q(1,10**8),'wrong phase yields wide true-equation defect')
    mutations.append('wrong phase sign cannot pass width gate')
    low=defect_certificate((Q(1),),Q(0),h,Q(512),2)
    require(low['Mexp'][1]>Q(1,10**8),'insufficient degree yields wide defect')
    mutations.append('insufficient degree cannot pass width gate')
    biased=defect_certificate((Q(1),),Q(1,100),h,Q(0),120)
    require(biased['M0'][1]>Q(1,10**8),'known biased source cannot pass width gate')
    mutations.append('declared biased source cannot pass width gate')
    encoded=output_moments(zero)
    encoded['Mexp']['real_lower']=rational(zero['Mexp'][0][0]+Q(1,10**30))
    try: validate_serialized(zero,encoded)
    except ValueError: mutations.append('inward serialization rejected')
    else: raise RuntimeError('inward serialization survived')
    require(centered_to_left((Q(1),Q(2),Q(3)),Q(1,2))==
            (Q(3,4),Q(-1),Q(3)),'exact centered-to-left translation')
    for radius in (Q(0),Q(1,7),Q(1,15*2**90),Q(123,101**200)):
        rounded=outward_radius(radius)
        require(radius<=rounded<radius+Q(1,2**512),
                'outward exact bounded dyadic radius '+str(radius.numerator))
    return checks,mutations


def fabricated_benchmark(models=None):
    """Fixed synthetic shape:2 sources x64panels x9 rational phase probes.

    Coefficients are explicit fabricated polynomials, not B or zB and never a
    registered source callback. This benchmarks the generic operator only.
    """
    start=time.monotonic(); count=0; checksum=hashlib.sha256(); max_me=Q(0)
    max_total_me=Q(0); max_total_mu=Q(0); whole_cases=0
    for source in range(2):
        for omega in (Q(0),Q(1,2**39),Q(1,2**11),Q(1,2),Q(2),Q(32),Q(128),Q(256),Q(512)):
            incoming_w=ZERO;incoming_u=ZERO;rw=Q(0);ru=Q(0)
            for panel in range(64):
                if models is None:
                    g=tuple(Q(((-1)**(n+source))*(panel+1),
                              (n+1)*(source+2)*2**(3*n)) for n in range(25))
                    error=Q(1,15*2**90)
                else:
                    g,error=models[source][panel]
                result=defect_certificate(g,error,Q(1,64),omega,96,
                                          coefficient_bits=512)
                payload=json.dumps(output_moments(result),sort_keys=True,separators=(',',':'))
                validate_serialized(result,json.loads(payload))
                checksum.update(payload.encode());count+=1
                max_me=max(max_me,result['Mexp'][1])
                # An independent whole-interval candidate carries point modes
                # and proved norm radii. State quantization is explicitly added.
                whole=defect_certificate(g,error,Q(1,64),omega,96,
                    coefficient_bits=512,incoming_w=incoming_w,incoming_u=incoming_u,
                    incoming_w_error=rw,incoming_u_error=ru)
                endw=scale(whole['Mexp'][0],Q(-1))
                endu=scale(whole['Duhamel_u'][0],Q(-1))
                incoming_w=chosen_complex(endw,512)
                incoming_u=chosen_complex(endu,512)
                rw=outward_radius(whole['Mexp'][1]+norm_upper(add(endw,scale(incoming_w,Q(-1)))))
                ru=outward_radius(whole['Duhamel_u'][1]+norm_upper(add(endu,scale(incoming_u,Q(-1)))))
            max_total_me=max(max_total_me,rw)
            max_total_mu=max(max_total_mu,ru)
            whole_payload=json.dumps({'Mexp':rectangle(scale(incoming_w,Q(-1)),rw),
                         'Duhamel_u':rectangle(scale(incoming_u,Q(-1)),ru)},
                         sort_keys=True,separators=(',',':'))
            checksum.update(whole_payload.encode());whole_cases+=1
    return {'scope':'fabricated generic operator only; no source-model construction or physical arrays',
            'cases':count,'sources':2,'panels_per_source':64,'phase_probes':9,
            'whole_interval_cases':whole_cases,'operator_certificates':2*count,
            'source_degree':24,'mode_degree':96,'coefficient_bits':512,
            'wall_seconds':time.monotonic()-start,
            'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'stable_serialized_moment_sha256':checksum.hexdigest(),
            'maximum_ME_radius':rational(max_me),
            'maximum_whole_interval_ME_radius':rational(max_total_me),
            'maximum_whole_interval_Mu_radius':rational(max_total_mu)}


def main():
    if len(sys.argv)!=3 or sys.argv[1]!='--fabricated-only':
        raise SystemExit('Usage: generic_operator.py --fabricated-only FRESH_RECEIPT.json')
    target=Path(sys.argv[2])
    if target.exists(): raise SystemExit('Fresh receipt required')
    checks,mutations=fabricated_tests()
    benchmark=fabricated_benchmark()
    receipt={'status':'PASS_FABRICATED_EXACT_OPERATOR_ONLY','python_optimization':sys.flags.optimize,
             'checks':checks,'mutation_controls':mutations,'benchmark':benchmark,
             'physical_source_evaluations':0,'retained_arrays_decoded':0,
             'primary_helper_imports':0,'source_model_certificate':'NOT_IMPLEMENTED_IN_THIS_GENERIC_EXECUTABLE',
             'full_twelve_case_direct_pressure_certificate':'UNRESOLVED',
             'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':receipt['status'],'checks':len(checks),'mutations':len(mutations),
                      'fabricated_cases':benchmark['cases'],'wall_seconds':benchmark['wall_seconds'],
                      'peak_rss_kib':benchmark['peak_rss_kib']}))


if __name__=='__main__':
    main()
