"""Independent incoming BD enclosure candidate; no source evaluation at import.

Physical entry is run(auth), under the parent's prospective registration guard.
run_fabricated() never constructs the registered bump or reads retained modes.
"""
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import hashlib
import json
import resource
import time
import generic_operator_baseline as op
import source_algebra_baseline as algebra

SOURCES=('positive_B','signed_uB')
MOMENTA=(Q(0),Q(1,2**40),Q(1,2**12),Q(1,4),Q(1),Q(16),Q(64),Q(128),Q(256))
SOURCE_DEGREE=112
MODE_DEGREE=160
BITS=512
EXP_DEGREE=200
DELTA=Q(1,128)
MAX_WIDTH=Q(1,64)
WALL_BUDGET=900
RSS_BUDGET_KIB=512*1024
OUTPUT_BUDGET_BYTES=20*1024*1024


def canonical(q): return str(q.numerator)+'/'+str(q.denominator)


def resource_guard(started):
    elapsed=time.monotonic()-started
    rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if elapsed>WALL_BUDGET or rss>RSS_BUDGET_KIB:
        raise RuntimeError('Registered independent resource budget exceeded')
    return {'wall_seconds':elapsed,'peak_rss_kib':rss}


def panels():
    left=DELTA
    out=[]
    while left<Q(1,2):
        right=min(3*left/2,Q(1,2))
        cx=(left+right)/2; half=(right-left)/2; radius=left/2
        rho=1-cx+radius; separation=5-cx-radius; den=1-rho*rho
        if not 0<half<=radius/2 or not 0<rho<1 or separation<=0:
            raise ValueError('Invalid geometric source disk')
        mb=2*(4/separation**2+4*rho/(separation*den**2)+2/den**2+8*rho**2/den**3+4*rho**2/den**4)
        mz=rho*mb+2*(2/separation+4*rho/den**2)
        pieces=(right-left)/MAX_WIDTH
        count=(pieces.numerator+pieces.denominator-1)//pieces.denominator
        width=(right-left)/count
        out.append({'left':left-5,'right':right-5,'center':cx-5,
                    'half':half,'radius':radius,'M':(mb,mz),'children':count,'width':width})
        left=right
    if len(out)!=11 or sum((r['right']-r['left'] for r in out),Q(0))!=Q(1,2)-DELTA:
        raise ValueError('Wrong registered source partition')
    return out


def exp_reduced_interval(q):
    """Directed dyadic enclosure: alternating exact series then squaring."""
    if not isinstance(q,Q) or q>0:
        raise ValueError('Exact nonpositive exponent required')
    reduced=q; squarings=0
    while reduced<-1:
        reduced/=2; squarings+=1
    lo,hi=algebra.exp_rational_interval(reduced,EXP_DEGREE)
    lo=op.chosen_dyadic(lo,BITS)
    hi=-op.chosen_dyadic(-hi,BITS)
    for _ in range(squarings):
        lo=op.chosen_dyadic(lo*lo,BITS)
        hi=-op.chosen_dyadic(-hi*hi,BITS)
    if not 0<=lo<=hi<=1+Q(1,2**BITS):
        raise ValueError('Invalid nonpositive exponential enclosure')
    return lo,hi,squarings


def physical_model(row,source,auth):
    # Parent authenticates the actual callbacks; no self-issued authorization.
    if not callable(auth):
        raise ValueError('Root authenticated source callback required')
    auth('real_source_taylor',{'source':source,'center':canonical(row['center']),
                             'scope':'UNIQUE_ANALYTIC_BD_PREHISTORY_ONLY'})
    if source not in SOURCES:
        raise ValueError('Unknown source')
    degree=SOURCE_DEGREE; center=row['center']; z=center+4
    inv=algebra.reciprocal_coefficients((1-z*z,-2*z,Q(-1)),degree+2)
    exponent=(1-inv[0],)+tuple(-v for v in inv[1:])
    lo,hi,squarings=exp_reduced_interval(exponent[0])
    beta=algebra.normalized_exp_coefficients(exponent,degree+2)
    h=beta if source==SOURCES[0] else tuple(z*beta[n]+(beta[n-1] if n else 0) for n in range(degree+3))
    ls=tuple(-(-Q(1))**n/center**(n+1) for n in range(degree+1))
    l2=tuple(sum((ls[j]*ls[n-j] for j in range(n+1)),Q(0)) for n in range(degree+1))
    gamma=tuple(4*sum((l2[j]*h[n-j] for j in range(n+1)),Q(0))
                -2*sum((ls[j]*(n-j+1)*h[n-j+1] for j in range(n+1)),Q(0))
                -(n+1)*(n+2)*h[n+2] for n in range(degree+1))
    chosen=[]; coefficient_error=Q(0)
    for n,g in enumerate(gamma):
        lower,upper=(g*lo,g*hi) if g>=0 else (g*hi,g*lo)
        point=op.chosen_dyadic((lower+upper)/2,BITS)
        chosen.append(point)
        coefficient_error+=max(point-lower,upper-point)*row['half']**n
    tail=row['M'][SOURCES.index(source)]/2**SOURCE_DEGREE
    error=op.outward_radius(tail+coefficient_error,BITS)
    return tuple(chosen),error,{'analytic_tail':canonical(tail),
        'coefficient_error':canonical(op.outward_radius(coefficient_error,BITS)),
        'uniform_source_error':canonical(error),'exp_squarings':squarings}


def fabricated_model(row,sid):
    """Exact finite polynomial with samedegree/conditioning; no real source."""
    # Same declared fabricated target as the independent Arb route, constructed
    # here directly as exact ordinary-power coefficients rather than a series.
    sign=1 if sid==0 else -1
    coefficients=tuple(Q(sign**n,n+1)/(2*row['half'])**n for n in range(SOURCE_DEGREE+1))
    return coefficients,Q(0),{'scope':'EXPLICIT_FINITE_POLYNOMIAL_FIXTURE','analytic_tail':'0/1','coefficient_error':'0/1','uniform_source_error':'0/1'}


def translate(coefficients,shift):
    """Horner translation p(x+shift), exact rational finite polynomial."""
    out=[]
    for value in reversed(coefficients):
        previous=out
        out=[Q(0)]*(len(previous)+1)
        for j,co in enumerate(previous):
            out[j]+=shift*co;out[j+1]+=co
        out[0]+=value
    return tuple(out)


def cap_bounds(source,k):
    """Exact IBP cap envelope; no bump/phase/source evaluation."""
    if source not in SOURCES or not isinstance(k,Q) or not 0<=k<=256:
        raise ValueError('Registered source and momentum required')
    D=DELTA*(2-DELTA); lc=1/(5-DELTA); H=Q(3,8)**63
    pd=2*(1-DELTA)/D**2; rho=0 if source==SOURCES[0] else 1
    w=H*(pd+rho+2*lc+2*k+DELTA*(4*k*k+4*k*lc+6*lc*lc))
    u=H*((pd+rho+2*lc)/2+1+DELTA*(3*lc*lc+2*k+2*lc))
    return u,w


def rectangle(center,radius,real=False):
    def interval(c):
        return {'lo':canonical(op.chosen_dyadic(c-radius,BITS)),
                'hi':canonical(-op.chosen_dyadic(-c-radius,BITS))}
    return {'real':interval(center[0]),'imag':{'lo':'0/1','hi':'0/1'} if real else interval(center[1])}


def exported_radius(rect):
    return sum(((Q(rect[key]['hi'])-Q(rect[key]['lo']))/2 for key in ('real','imag')),Q(0))


def execute(auth=None,fabricated=False):
    started=time.monotonic(); geom=panels(); models={}; model_evidence=[]
    for sid,source in enumerate(SOURCES):
        children=[]
        for parent_id,row in enumerate(geom):
            coeff,error,ev=fabricated_model(row,sid) if fabricated else physical_model(row,source,auth)
            for child_id in range(row['children']):
                left=row['left']+child_id*row['width']
                local=translate(coeff,left-row['center'])
                children.append((local,error,row['width'],parent_id,child_id))
            model_evidence.append({'source':source,'parent':parent_id,
                'center':canonical(row['center']),'interval':[canonical(row['left']),canonical(row['right'])],
                'half_width':canonical(row['half']),'analytic_radius':canonical(row['radius']),
                'disk_bound':canonical(row['M'][sid]),'children':row['children'],
                'chosen_coefficients_sha256':hashlib.sha256(json.dumps([canonical(c) for c in coeff],separators=(',',':')).encode()).hexdigest(),**ev})
            resource_guard(started)
        models[source]=children
    whole_rows=[];error_rows=[];whole_budgets=[]
    for source in SOURCES:
        for k in MOMENTA:
            iw=op.ZERO;iu=op.ZERO;rw=Q(0);ru=Q(0)
            categories={name:{'U':Q(0),'W':Q(0)} for name in ('model','defect','point','radius')}
            for coeff,error,width,parent_id,child_id in models[source]:
                result=op.defect_certificate(coeff,error,width,2*k,MODE_DEGREE,
                    coefficient_bits=BITS,incoming_w=iw,incoming_u=iu,
                    incoming_w_error=rw,incoming_u_error=ru)
                ew=op.scale(result['Mexp'][0],Q(-1));eu=op.scale(result['Duhamel_u'][0],Q(-1))
                nw=op.chosen_complex(ew,BITS);nu=op.chosen_complex(eu,BITS)
                sw=op.norm_upper(op.add(ew,op.scale(nw,Q(-1))))
                su=op.norm_upper(op.add(eu,op.scale(nu,Q(-1))))
                nrw=op.outward_radius(result['Mexp'][1]+sw,BITS)
                nru=op.outward_radius(result['Duhamel_u'][1]+su,BITS)
                increments={
                    'model':{'W':error*width,'U':error*width**2/2},
                    'defect':{'W':result['Mexp'][1]-rw-error*width,
                              'U':result['Duhamel_u'][1]-ru-width*rw-error*width**2/2},
                    'point':{'W':sw,'U':su},
                    'radius':{'W':nrw-result['Mexp'][1]-sw,'U':nru-result['Duhamel_u'][1]-su}}
                for name in categories:
                    categories[name]['U']+=width*categories[name]['W']+increments[name]['U']
                    categories[name]['W']+=increments[name]['W']
                if sum((v['W'] for v in categories.values()),Q(0))!=nrw or sum((v['U'] for v in categories.values()),Q(0))!=nru:
                    raise RuntimeError('Independent whole error budget does not telescope')
                if any(value<0 for row in increments.values() for value in row.values()):
                    raise RuntimeError('Negative independent error component')
                error_rows.append({'source':source,'momentum':canonical(k),'parent':parent_id,'child':child_id,
                    'width':canonical(width),'source_W':canonical(op.outward_radius(error*width,BITS)),
                    'source_U':canonical(op.outward_radius(error*width**2/2,BITS)),
                    'defect_W':canonical(op.outward_radius(result['Mexp'][1]-rw-error*width,BITS)),
                    'defect_U':canonical(op.outward_radius(result['Duhamel_u'][1]-ru-width*rw-error*width**2/2,BITS)),
                    'point_shift_W':canonical(op.outward_radius(sw,BITS)),
                    'point_shift_U':canonical(op.outward_radius(su,BITS)),
                    'radius_round_W':canonical(nrw-result['Mexp'][1]-sw),
                    'radius_round_U':canonical(nru-result['Duhamel_u'][1]-su)})
                iw,iu,rw,ru=nw,nu,nrw,nru
                resource_guard(started)
            cu,cw=(Q(0),Q(0)) if fabricated else cap_bounds(source,k)
            urect=rectangle(iu,ru+cu,k==0);wrect=rectangle(iw,rw+cw,k==0)
            whole_rows.append({'source':source,'momentum':canonical(k),'U':urect,'W':wrect,
                'total_absolute_radii':{'U':canonical(exported_radius(urect)),'W':canonical(exported_radius(wrect))},
                'cap_error':{'U':canonical(cu),'W':canonical(cw)},
                'interior_norm_radii':{'U':canonical(ru),'W':canonical(rw)}})
            dimensions=1 if k==0 else 2
            excess={'U':exported_radius(urect)-dimensions*(ru+cu),
                    'W':exported_radius(wrect)-dimensions*(rw+cw)}
            if any(value<0 or value>=Q(dimensions,2**BITS) for value in excess.values()):
                raise RuntimeError('Outward export bound or complete radius failed')
            budget={'source':source,'momentum':canonical(k),
                'cap_disk_radius':{'U':canonical(cu),'W':canonical(cw)},
                'complete_output_L1_radius':{'U':canonical(exported_radius(urect)),'W':canonical(exported_radius(wrect))},
                'export_excess_L1_radius':{key:canonical(value) for key,value in excess.items()},
                'output_real_dimensions':dimensions,
                'inflation_factor':canonical(Q(dimensions)),
                'budget_identity':'completeL1 = real_dimensions*(cap + model + ODEdefect + point_round + radius_round) + export_excess'}
            for category,field in (('model','source_model_disk_radius'),('defect','ode_defect_disk_radius'),('point','point_round_disk_radius'),('radius','radius_round_disk_radius')):
                budget[field]={key:canonical(value) for key,value in categories[category].items()}
            whole_budgets.append(budget)
    payload={'schema_version':1,'scope':'FABRICATED_ONLY' if fabricated else 'BD_PREHISTORY_TARGET_AT_FIXED_RATIONAL_PROBES',
        'target':'U=u/represented_epsilon, W=w/represented_epsilon; exact source target and exact zero before−5',
        'domain':['-5/1','-9/2'],'archive_arrays_decoded':0,'source_callbacks':0 if fabricated else 22,
        'whole_rows':whole_rows,'configuration':'INDEPENDENT_DYADIC512_SOURCE112_ODE160',
        'original_binary80_state_error':'NOT_ENCLOSED','full_twelve_case_certificate':'UNRESOLVED'}
    evidence={'method':'independent exact rational source recurrences and polynomial true-equation defect',
        'schema_version':1,'scope':payload['scope'],'configuration':payload['configuration'],
        'source_callbacks':payload['source_callbacks'],'archive_arrays_decoded':0,
        'original_binary80_state_error':'NOT_ENCLOSED','full_twelve_case_certificate':'UNRESOLVED',
        'fabricated_source_provider':fabricated,'physical_source_evaluations':0 if fabricated else 22,
        'source_degree':SOURCE_DEGREE,'mode_degree':MODE_DEGREE,'coefficient_bits':BITS,
        'exp_degree':EXP_DEGREE,'geometric_panels':len(geom),'ode_children_per_source':len(models[SOURCES[0]]),
        'source_model_rows':model_evidence,'mode_error_rows':error_rows,'whole_rows':whole_budgets,
        'retained_arrays_decoded':0,'primary_helper_imports':0,
        'global_state_origin':'zero polynomial state at capend, independent omitted real-source cap enclosed at finalendpoint',
        'resource':resource_guard(started),'outward_export_bits':BITS,
        'exported_radius_definition':'L1 sum of exported real/imag rectangle halfwidths',
        'zero_component_policy':'Exact rational sums of individually nonnegative checked increments; source-model zero only for the exact finite-polynomial fabricated fixture; no artificial padding.',
        'k0_imaginary_projection':'exact real source and zero initial data theorem only; never original rounded inputs'}
    if len(json.dumps([payload,evidence],separators=(',',':')).encode())>OUTPUT_BUDGET_BYTES:
        raise RuntimeError('Independent output budget exceeded')
    return payload,evidence


def run(auth): return execute(auth,False)


def run_fabricated(): return execute(None,True)
