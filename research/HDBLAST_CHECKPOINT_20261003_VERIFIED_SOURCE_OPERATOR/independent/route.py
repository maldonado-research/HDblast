"""Independent exact-rational source/operator route; no evaluation at import."""
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import resource
import time
import generic_operator as op
import source_models as sm

SOURCES=('positive_B','signed_uB')
MOMENTA=(Q(0),Q(1,2**40),Q(1,2**12),Q(1,4),Q(1),Q(16),Q(64),Q(128),Q(256))
WIDTH=Q(1,64)
HALF=Q(1,128)
DEGREE=96
BITS=512


def canonical(q):
    return str(q.numerator)+'/'+str(q.denominator)


def interval_rect(center,radius,real=False):
    if radius<0: raise ValueError('Nonnegative exact interval radius required')
    def bounds(c):
        return {'lo':canonical(op.chosen_dyadic(c-radius,BITS)),
                'hi':canonical(-op.chosen_dyadic(-c-radius,BITS))}
    return {'real':bounds(center[0]),
            'imag':{'lo':'0/1','hi':'0/1'} if real else
                   bounds(center[1])}


def exported_radius(rect):
    """Actual L1 half-width of the exported outward rectangle."""
    return sum(((Q(rect[key]['hi'])-Q(rect[key]['lo']))/2
                for key in ('real','imag')),Q(0))


def bound(q):
    """Evidence values are explicit upper bounds, not unbounded exact strings."""
    return canonical(op.outward_radius(q,BITS))


def digest_model(coeff,error):
    value={'coefficients':[canonical(v) for v in coeff],'error':canonical(error)}
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def model_rows(auth):
    out={s:[] for s in SOURCES};evidence=[]
    for source in SOURCES:
        for panel in range(64):
            center=Q(-9,2)+(2*panel+1)*HALF
            auth('real_source_taylor',{'source':source,'center':str(center)})
            model=sm.registered_source_model(center,source,auth.authorization)
            parts=[];ev={'source':source,'panel':panel,'center':canonical(center)}
            for key in ('forcing','Lg'):
                part=model[key]
                left=op.centered_to_left(part['center_coefficients'],HALF)
                error=part['uniform_source_error']
                parts.append((left,error))
                ev[key]={'analytic_tail':canonical(part['analytic_tail']),
                         'coefficient_error_uniform':canonical(error-part['analytic_tail']),
                         'uniform_error':canonical(error),
                         'left_chosen_model_sha256':digest_model(left,error)}
            out[source].append(tuple(parts));evidence.append(ev)
    return out,evidence


def fake_model_rows():
    # Explicit fake model provider; it never calls auth or registered_source_model.
    import fabricated_complete_route as fixture
    force,work,construction=fixture.fabricated_models()
    out={s:[] for s in SOURCES};evidence=[]
    for sid,source in enumerate(SOURCES):
        for panel in range(64):
            forcing=force[sid][panel];lg=work[sid][panel]
            out[source].append((forcing,lg))
            evidence.append({'source':source,'panel':panel,
                'scope':'EXPLICIT_FABRICATED_FINITE_POLYNOMIAL_PROVIDER',
                'forcing':{'uniform_error':canonical(forcing[1]),
                           'left_chosen_model_sha256':digest_model(*forcing)},
                'Lg':{'uniform_error':canonical(lg[1]),
                      'left_chosen_model_sha256':digest_model(*lg)}})
    return out,evidence,construction


def _execute(models,source_evidence,started,fabricated,construction=None):
    panels=[];wholes=[];work_panels=[];work_wholes=[];errors=[]
    for source in SOURCES:
        # Source-work is computed from its own source polynomial, not a mode primitive.
        work_total=Q(0);work_radius=Q(0)
        work_radius_round=Q(0)
        m0_total=Q(0);m0_radius=Q(0);m0_radius_round=Q(0)
        for panel,((g,gerror),(lg,lgerror)) in enumerate(models[source]):
            a=Q(-9,2)+panel*WIDTH;b=a+WIDTH
            center=op.exact_integral(lg,WIDTH);radius=lgerror*WIDTH
            work_rect=interval_rect((center,Q(0)),radius,True)
            work_panels.append({'source':source,'panel':panel,'interval':[canonical(a),canonical(b)],
                'moment':work_rect,'total_absolute_radius':canonical(exported_radius(work_rect))})
            work_total+=center
            next_radius=op.outward_radius(work_radius+radius,BITS)
            work_radius_round+=next_radius-(work_radius+radius)
            work_radius=next_radius
            m0_total+=op.exact_integral(g,WIDTH)
            next_radius=op.outward_radius(m0_radius+gerror*WIDTH,BITS)
            m0_radius_round+=next_radius-(m0_radius+gerror*WIDTH)
            m0_radius=next_radius
        work_rect=interval_rect((work_total,Q(0)),work_radius,True)
        work_wholes.append({'source':source,'interval':['-9/2','-7/2'],
            'moment':work_rect,'total_absolute_radius':canonical(exported_radius(work_rect))})
        whole_for_source=[]
        by_panel=[[] for _ in range(64)]
        for k in MOMENTA:
            incoming_w=op.ZERO;incoming_u=op.ZERO;rw=Q(0);ru=Q(0)
            radius_round_w=Q(0);radius_round_u=Q(0)
            point_round_w=Q(0);point_round_u=Q(0)
            for panel,((g,error),_) in enumerate(models[source]):
                local=op.defect_certificate(g,error,WIDTH,2*k,DEGREE,coefficient_bits=BITS)
                whole=op.defect_certificate(g,error,WIDTH,2*k,DEGREE,coefficient_bits=BITS,
                    incoming_w=incoming_w,incoming_u=incoming_u,
                    incoming_w_error=rw,incoming_u_error=ru)
                a=Q(-9,2)+panel*WIDTH;b=a+WIDTH
                m0,r0=local['M0'];me,re=local['Mexp'];mu,rm=local['Duhamel_u']
                moment_rects={'M0':interval_rect((m0,Q(0)),r0,True),
                               'Mexp':interval_rect(me,re,k==0),'Mu':interval_rect(mu,rm,k==0)}
                by_panel[panel].append({'source':source,'panel':panel,'momentum':canonical(k),
                    'interval':[canonical(a),canonical(b)],
                    'moments':moment_rects,
                    'total_absolute_radii':{name:canonical(exported_radius(rect))
                                            for name,rect in moment_rects.items()}})
                endw=op.scale(whole['Mexp'][0],Q(-1));endu=op.scale(whole['Duhamel_u'][0],Q(-1))
                nextw=op.chosen_complex(endw,BITS);nextu=op.chosen_complex(endu,BITS)
                shiftw=op.norm_upper(op.add(endw,op.scale(nextw,Q(-1))))
                shiftu=op.norm_upper(op.add(endu,op.scale(nextu,Q(-1))))
                nwr=op.outward_radius(whole['Mexp'][1]+shiftw,BITS)
                nur=op.outward_radius(whole['Duhamel_u'][1]+shiftu,BITS)
                extra_w=nwr-whole['Mexp'][1]-shiftw
                extra_u=nur-whole['Duhamel_u'][1]-shiftu
                errors.append({'source':source,'panel':panel,'momentum':canonical(k),
                    'source_ME_error':bound(error*WIDTH),
                    'source_Mu_error':bound(error*WIDTH**2/2),
                    'local_ME_defect_error':bound(re-error*WIDTH),
                    'local_Mu_defect_error':bound(rm-error*WIDTH**2/2),
                    'prefix_ME_defect_increment':bound(whole['Mexp'][1]-rw-error*WIDTH),
                    'prefix_Mu_defect_increment':bound(whole['Duhamel_u'][1]-ru-WIDTH*rw-error*WIDTH**2/2),
                    'inherited_w_radius':bound(rw),'inherited_u_radius':bound(ru),
                    'chosen_w_point_shift':bound(shiftw),'chosen_u_point_shift':bound(shiftu),
                    'w_radius_upward_rounding':bound(extra_w),
                    'u_radius_upward_rounding':bound(extra_u)})
                incoming_w=nextw;incoming_u=nextu;rw=nwr;ru=nur
                point_round_w+=shiftw;point_round_u+=shiftu
                radius_round_w+=extra_w;radius_round_u+=extra_u
            moment_rects={'M0':interval_rect((m0_total,Q(0)),m0_radius,True),
                           'Mexp':interval_rect(op.scale(incoming_w,Q(-1)),rw,k==0),
                           'Mu':interval_rect(op.scale(incoming_u,Q(-1)),ru,k==0)}
            whole_for_source.append({'source':source,'momentum':canonical(k),
                'interval':['-9/2','-7/2'],
                'moments':moment_rects,
                'total_absolute_radii':{name:canonical(exported_radius(rect))
                                       for name,rect in moment_rects.items()}})
        for rows in by_panel:panels.extend(rows)
        wholes.extend(whole_for_source)
    payload={'schema_version':1,'scope':'UNIFORM_ANALYTIC_MODEL_PLUS_FIXED_RATIONAL_MOMENT_PROBES',
             'configuration':'INDEPENDENT_DYADIC512_SOURCE24_ODE96','archive_arrays_decoded':0,
             'panel_rows':panels,'whole_rows':wholes,
             'source_work_panel_rows':work_panels,'source_work_whole_rows':work_wholes}
    evidence={'method':'independent exact polynomial ODE defect and rational source model',
              'source_degree':24,'mode_degree':96,'scalar_exp_degree':200,'dyadic_bits':512,
              'source_model_rows':source_evidence,'mode_error_rows':errors,
              'physical_source_evaluations':0 if fabricated else 128,
              'fabricated_source_provider':fabricated,'retained_arrays_decoded':0,'primary_helper_imports':0,
              'global_state_origin':'exact zero before the first source panel',
              'source_forcing_disk_bound':'64/1','Lg_disk_bound':'32/1',
              'evidence_error_values':'each error component is independently rounded upward on512bit grid',
              'exported_radius_definition':'L1 sum of actual real/imag rectangle halfwidths',
              'exported_endpoint_rounding_per_component_upper':'1/'+str(2**512),
              'resource':{'wall_seconds':time.monotonic()-started,
                          'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},
              'construction_benchmark':construction}
    # Root's bounded worker/schema validator remains authoritative.
    return payload,evidence


def run(auth):
    started=time.monotonic()
    models,evidence=model_rows(auth)
    return _execute(models,evidence,started,False)


def run_fabricated():
    started=time.monotonic()
    models,evidence,construction=fake_model_rows()
    return _execute(models,evidence,started,True,construction)
