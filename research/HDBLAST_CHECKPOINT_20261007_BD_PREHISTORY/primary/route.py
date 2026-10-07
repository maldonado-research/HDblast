"""Prepared primary route. CLI permits fabricated fixtures only before freeze."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import resource
import time
from flint import arb, acb, arb_series, ctx
if __package__:
    from . import verified_incoming as v
else:
    import verified_incoming as v

SOURCES = ('positive_B','signed_uB')
MOMENTA = (F(0),F(1,2**40),F(1,2**12),F(1,4),F(1),F(16),F(64),F(128),F(256))
CONFIGURATION = 'PRIMARY_ARB512_SOURCE112_PHASE128'
GATE = F(1,10**20)

def require(condition,message):
    if not condition:
        raise ValueError(message)

def source_disks(panels):
    return {'W':sum((2*p.half_width*p.remainder for p in panels),F(0)),
            'U':sum((2*p.half_width*(v.B-p.center)*p.remainder for p in panels),F(0))}

def phase_disks(panels,k,family):
    w,u = F(0),F(0)
    if not k:
        return {'U':u,'W':w}
    for p in panels:
        h = p.half_width
        distance = v.B-p.center-h
        c_l1 = sum((v.upper(x) for x in p.coefficients),F(0))
        w += h*c_l1*family.e_tail
        u += h*h*c_l1*family.q_tail
        if 2*k*distance<=1:
            from math import factorial
            u += distance*F(3,factorial(v.PHASE_DEGREE+2))*(2*h*c_l1)
    return {'U':u,'W':w}

def run_primary(authorize,provider=None,fabricated=False):
    require(callable(authorize),'Root authorization required')
    require(not fabricated or provider is not None,'Fabricated route requires explicit provider')
    require(fabricated or provider is None,'Physical provider must be fixed registered source')
    provider = provider or v.registered_parent
    data = {'schema_version':1,'configuration':CONFIGURATION,
            'scope':'FABRICATED_ONLY' if fabricated else 'BD_PREHISTORY_TARGET_AT_FIXED_RATIONAL_PROBES',
            'archive_arrays_decoded':0,'source_callbacks':0 if fabricated else 22,
            'original_binary80_state_error':'NOT_ENCLOSED','full_twelve_case_certificate':'UNRESOLVED',
            'whole_rows':[],'source_model_rows':[],
            'model_scope':'Continuous source and cap representation truncation for 0<=k<=256; full arithmetic only at nine rational probes.'}
    budget = {'schema_version':1,'configuration':CONFIGURATION,'whole_rows':[],
              'analytic_model':{'source_degree':v.DEGREE,'phase_degree':v.PHASE_DEGREE,
               'precision_bits':v.PRECISION_BITS,'parents':11,'cap_width':v.qstring(v.DELTA),
               'max_child_width':'1/64','local_phase_cap':'4/1','integrand_phase_cap':'8/1',
               'exp_majorant':'4096/1','canonical_c_projection':False},
              'gate':'1e-20 complete exported U/W complex L1 rectangle radius',
              'budget_semantics':'Exact exported radius relation: complete=baseline+factor*(cap+source_model+phase_model)+export_excess, with factor1 at k0 exact real or2 for the complex rectangle. Baseline encloses every source-coefficient, affine-recentering, finite-kernel and global-transcendental arithmetic operation. Model contributions are proved disk bounds and export_excess covers the final outward inflation/serialization. These bounds are conservative, not individually measured errors.',
              'physical_source_evaluations':0 if fabricated else 22,'archive_arrays_decoded':0}
    started = time.monotonic()
    child_count = None
    with ctx.workprec(v.PRECISION_BITS):
        family = v.Kernels()
        for source in SOURCES:
            parents = []
            children = []
            for index,geom in enumerate(v.PARENTS):
                parent = provider(index,source,authorize)
                majorant = geom['positive_majorant'] if source=='positive_B' else geom['signed_majorant']
                require(parent.parent_index==index and parent.center==geom['center'] and
                        parent.half_width==geom['half_width'] and
                        parent.remainder==majorant/2**v.DEGREE,'Unregistered source panel')
                parents.append(parent)
                new_children = v.child_panels(parent)
                children.extend(new_children)
                data['source_model_rows'].append({'source':source,'parent':index,
                     'center':v.qstring(parent.center),'half_width':v.qstring(parent.half_width),
                     'disk_radius':v.qstring(geom['disk_radius']),
                     'source_majorant':v.qstring(majorant),'source_uniform_tail':v.qstring(parent.remainder),
                     'coefficient_count':len(parent.coefficients),'children':len(new_children),
                     'coefficients':[v.rectangle(coefficient)['real'] for coefficient in parent.coefficients]})
            require(child_count is None or child_count==len(children),'Source-dependent child geometry')
            child_count = len(children)
            sd = source_disks(children)
            cap = v.cap_bounds(source)
            for k in MOMENTA:
                baseline = v.incoming(children,k,family,False)
                pd = phase_disks(children,k,family)
                row = {'source':source,'momentum':v.qstring(k),'total_absolute_radii':{},
                       'cap_error':{name:v.qstring(value) for name,value in cap.items()}}
                baseline_radii = {}
                inflation = {}
                excess = {}
                for name in ('U','W'):
                    error = cap[name]+sd[name]+pd[name]
                    value = baseline[name]
                    if not k:
                        value = acb(value.real,arb(0))
                    baseline_radius = v.radius(v.rectangle(value))
                    target = v.inflate_complex(value,error)
                    # Exact real-target identity at k=0 only; no mode/Wronskian projection.
                    if not k:
                        target = acb(target.real,arb(0))
                    row[name] = v.rectangle(target)
                    row['total_absolute_radii'][name] = v.qstring(v.radius(row[name]))
                    baseline_radii[name] = v.qstring(baseline_radius)
                    inflation[name] = v.qstring(error)
                    export_excess = F(row['total_absolute_radii'][name])-baseline_radius-(2 if k else 1)*error
                    require(export_excess>=0,'Outward export lost a positive model radius')
                    excess[name] = v.qstring(export_excess)
                data['whole_rows'].append(row)
                budget['whole_rows'].append({'source':source,'momentum':v.qstring(k),
                    'source_model_disk_radius':{name:v.qstring(value) for name,value in sd.items()},
                    'phase_model_disk_radius':{name:v.qstring(value) for name,value in pd.items()},
                    'cap_disk_radius':row['cap_error'],
                    'coefficient_and_arithmetic_baseline_L1_radius':baseline_radii,
                    'baseline_L1_radius':baseline_radii,'inflation_disk_radius':inflation,
                    'inflation_factor':'2/1' if k else '1/1','export_excess_L1_radius':excess,
                    'complete_output_L1_radius':row['total_absolute_radii']})
    require(len(data['whole_rows'])==18 and len(data['source_model_rows'])==22,'Incomplete fixed universe')
    max_radius = max(F(value) for row in data['whole_rows'] for value in row['total_absolute_radii'].values())
    budget['maximum_complete_L1_radius'] = v.qstring(max_radius)
    budget['status'] = 'PASS_NARROW_INCOMING_TARGET_WIDTH_GATE' if max_radius<=GATE else 'UNRESOLVED_CERTIFICATE'
    budget['child_panels_per_source'] = child_count
    budget['wall_seconds'] = time.monotonic()-started
    budget['peak_rss_kib'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return data,budget

def fabricated_parent(index,source,authorize):
    authorize('fabricated_source',{'source':source,'parent':index})
    p = v.PARENTS[index]
    sign = 1 if source=='positive_B' else -1
    coefficients = tuple(v.ball(F(sign**j,2**j*(j+1))) for j in range(v.DEGREE+1))
    majorant = p['positive_majorant'] if source=='positive_B' else p['signed_majorant']
    return v.Panel(p['center'],p['half_width'],coefficients,majorant/2**v.DEGREE,index)

def manufactured_operations():
    """Unrelated positive-center manufactured source jets only, no bump."""
    previous_cap = ctx.cap
    with ctx.workprec(v.PRECISION_BITS):
        try:
            ctx.cap=v.DEGREE+3
            for source in SOURCES:
                for p in v.PARENTS:
                    t=arb_series([v.ball(2+p['center']+5),v.ball(p['half_width'])],prec=v.DEGREE+3)
                    manufactured=(t*(2+t*t).inv()).exp()
                    if source=='signed_uB':
                        manufactured=t*manufactured
                    coefficients=v.forcing_coefficients(t,manufactured,p['half_width'])
                    require(len(coefficients)==v.DEGREE+1 and all(x.is_finite() for x in coefficients),
                            'Manufactured jet rehearsal incomplete')
        finally:
            ctx.cap=previous_cap
    return {'scope':'UNRELATED_POSITIVE_CENTER_MANUFACTURED_JETS','operations':22,
            'coefficient_count':v.DEGREE+1,'registered_source_calls':0}

def run_fabricated():
    rehearsal=manufactured_operations()
    def auth(name,detail):
        require(name=='fabricated_source','Physical source callback in fabricated route')
    data,budget=run_primary(auth,fabricated_parent,True)
    budget['manufactured_operation_rehearsal']=rehearsal
    return data,budget

def run(auth):
    return run_primary(auth)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--fabricated',action='store_true')
    parser.add_argument('--output-directory',type=Path,required=True)
    args=parser.parse_args()
    require(args.fabricated,'Only root authenticated wrapper may run physical callback')
    require(not args.output_directory.exists(),'Refusing to replace existing output')
    started=time.monotonic()
    data,budget=run_fabricated()
    require(budget['status']=='PASS_NARROW_INCOMING_TARGET_WIDTH_GATE','Fabricated width gate failed')
    args.output_directory.mkdir(parents=True)
    for filename,value in [('PRIMARY.json',data),('PRIMARY_ERROR_BUDGET.json',budget)]:
        (args.output_directory/filename).write_text(json.dumps(value,sort_keys=True)+'\n')
    receipt={'status':'PASS_COMPLETE_FABRICATED_PRIMARY_ROUTE','physical_source_evaluations':0,
        'archive_arrays_decoded':0,'rows':len(data['whole_rows']),'parent_source_rows':22,
        'wall_seconds':time.monotonic()-started,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'output_bytes':sum(x.stat().st_size for x in args.output_directory.iterdir()),
        'route_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'engine_sha256':hashlib.sha256(Path(v.__file__).read_bytes()).hexdigest(),
        'maximum_complete_L1_radius':budget['maximum_complete_L1_radius']}
    (args.output_directory/'RECEIPT.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps(receipt,sort_keys=True))

if __name__=='__main__':
    main()
