"""Meaningful fabricated target checks and negative controls, no bump callback."""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import hashlib
import json
import sys
from flint import arb, acb, ctx
import verified_incoming as v
import route

checks=[]
controls=[]

def require(value,label):
    if not value:
        raise RuntimeError(label)
    checks.append(label)

def rejects(function,label):
    try:
        function()
    except (ValueError,RuntimeError):
        controls.append(label)
        return
    raise RuntimeError('Mutation accepted: '+label)

def encoded_overlap(encoded,reference):
    other=v.rectangle(reference)
    return all(F(encoded[name]['lo'])<=F(other[name]['hi']) and
               F(other[name]['lo'])<=F(encoded[name]['hi']) for name in ('real','imag'))

def constant_parent(index,source,authorize):
    authorize('fabricated_source',{'source':source,'parent':index})
    p=v.PARENTS[index]
    majorant=p['positive_majorant'] if source=='positive_B' else p['signed_majorant']
    sign=1 if source=='positive_B' else -1
    coefficients=(v.ball(sign),)+tuple(v.ball(0) for _ in range(v.DEGREE))
    return v.Panel(p['center'],p['half_width'],coefficients,majorant/2**v.DEGREE,index)

def fabricated_auth(name,detail):
    if name!='fabricated_source':
        raise RuntimeError('Physical source attempt forbidden in fabricated tests')

with ctx.workprec(v.PRECISION_BITS):
    # A rational proof for every exponential tail constant, without exp calls.
    require(sum(F(1,factorial(j)) for j in range(4))==F(8,3),'exp positive sum to degree3')
    require(F(8,3)+F(5,96)==F(87,32)<F(11,4),'e upper bound from geometric tail')
    require(F(11,4)**8<F(4096),'exp8 majorant4096')
    require(F(87,32)<3,'stable entire Phi majorant3')
    table=v.mixed_table()
    for j in (0,1,2,7,24,56,112):
        for m in (0,1,2,32,64,128,129):
            require(table[j][m]==v.mixed_exact(j,m),f'exact mixed recurrence j{j}m{m}')
    require(len(table)==113 and all(len(row)>=130 for row in table),'all113 source powers supported')
    # Affine recentering on each child must preserve a chosen exact parent
    # polynomial; test a sparse degree112 term as well as a low-degree jet.
    p=v.PARENTS[-2]
    coeff=[v.ball(0)]*113
    for j in (0,1,3,112):
        coeff[j]=v.ball(F((-1)**j,j+1))
    parent=v.Panel(p['center'],p['half_width'],tuple(coeff),F(0),9)
    children=v.child_panels(parent)
    for ci,child in enumerate(children):
        for x in (F(-1),F(0),F(1)):
            eta=child.center+child.half_width*x
            y=(eta-parent.center)/parent.half_width
            expected=sum((F((-1)**j,j+1)*y**j for j in (0,1,3,112)),F(0))
            actual=sum((coefficient*v.ball(x**j) for j,coefficient in enumerate(child.coefficients)),arb(0))
            require(F(str(actual.lower().fmpq()))<=expected<=F(str(actual.upper().fmpq())),
                    f'affine sparse degree112 child{ci}x{x}')
    data,budget=route.run_primary(fabricated_auth,constant_parent,True)
    require(budget['status']=='PASS_NARROW_INCOMING_TARGET_WIDTH_GATE','complete constant fixture width gate')
    distance=v.B-v.A
    for row in data['whole_rows']:
        k=F(row['momentum'])
        sign=1 if row['source']=='positive_B' else -1
        if not k:
            expected_w=acb(v.ball(-sign*distance),arb(0))
            expected_u=acb(v.ball(-sign*distance*distance/2),arb(0))
        else:
            omega=acb(0,v.ball(2*k))
            phi=((omega*v.ball(distance)).exp()-1)/omega
            expected_w=-sign*phi
            expected_u=-sign*(phi-v.ball(distance))/omega
        require(encoded_overlap(row['W'],expected_w),'independent constant W '+row['source']+row['momentum'])
        require(encoded_overlap(row['U'],expected_u),'independent constant U '+row['source']+row['momentum'])
        for name in ('U','W'):
            require(v.radius(row[name])==F(row['total_absolute_radii'][name]),'actual export radius '+name+row['source']+row['momentum'])
        if not k:
            require(row['U']['imag']==row['W']['imag']=={'lo':'0/1','hi':'0/1'},'exact k0 reality '+row['source'])
    for row in budget['whole_rows']:
        for name in ('U','W'):
            inflation=sum(F(row[key][name]) for key in ('cap_disk_radius','source_model_disk_radius','phase_model_disk_radius'))
            require(inflation==F(row['inflation_disk_radius'][name]),'model budget sum '+name+row['source']+row['momentum'])
            expected=F(row['baseline_L1_radius'][name])+F(row['inflation_factor'])*inflation+F(row['export_excess_L1_radius'][name])
            require(expected==F(row['complete_output_L1_radius'][name]),'exact budget radius relation '+name+row['source']+row['momentum'])
    k0=next(row for row in data['whole_rows'] if row['source']=='positive_B' and row['momentum']=='0/1')
    k16=next(row for row in data['whole_rows'] if row['source']=='positive_B' and row['momentum']=='16/1')
    require(not encoded_overlap(k0['W'],acb(v.ball(distance),arb(0))),'reject source work/sign reversal at k0')
    controls.append('source sign reversal excluded by independent constant reference')
    require(not encoded_overlap(k0['W'],acb(v.ball(-distance-F(1,1000)),arb(0))),'reject biased source reference')
    controls.append('biased source excluded by independent constant reference')
    omega=acb(0,v.ball(32))
    true_w=-((omega*v.ball(distance)).exp()-1)/omega
    require(not encoded_overlap(k16['W'],true_w.conjugate()),'reject phase conjugation at k16')
    controls.append('phase conjugation excluded by independent constant reference')
    family=v.Kernels()
    rejects(lambda:family.at(F(257),F(1,128)),'out-of-band momentum')
    rejects(lambda:family.at(F(256),F(1,64)),'unproved local phase cap')
    rejects(lambda:v.Panel(F(-4),F(1,128),(arb(0),)*112,F(0),0),'truncated source coefficients')
    rejects(lambda:v.Panel(F(-4),F(1,128),(arb(0),)*113,F(-1),0),'negative source remainder')
    rejects(lambda:v.registered_parent(0,'positive_B',None),'real callback without authorization')
    rejects(lambda:v.registered_parent(0,'positive_B',fabricated_auth),'real callback through fabricated authorization')
    auth_events=[]
    def callback_shape_guard(name,detail):
        require(name=='real_source_taylor','real callback event label')
        require(set(detail)=={'source','center','scope'},'common root callback key shape')
        require(detail['scope']=='UNIQUE_ANALYTIC_BD_PREHISTORY_ONLY','common root callback scope')
        auth_events.append(detail)
        raise RuntimeError('Shape rehearsal stops before first real coefficient')
    rejects(lambda:v.registered_parent(0,'positive_B',callback_shape_guard),
            'root callback shape rehearsal rejects before real coefficients')
    require(len(auth_events)==1,'single pre-coefficient shape rehearsal event')
    def bad_parent(index,source,auth):
        parent=constant_parent(index,source,auth)
        return v.Panel(parent.center,parent.half_width,parent.coefficients,parent.remainder/2,index)
    rejects(lambda:route.run_primary(fabricated_auth,bad_parent,True),'incorrect source remainder registration')
    rejects(lambda:route.run_primary(fabricated_auth,constant_parent,False),'physical callback provider override')

result={'status':'PASS_FABRICATED_TARGET_AND_MUTATION_CHECKS','checks_passed':len(checks),
        'rejected_controls_count':len(controls),'checks':checks,'rejected_controls':controls,
        'physical_source_evaluations':0,'archive_arrays_decoded':0,'python_optimization':sys.flags.optimize,
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'route_sha256':hashlib.sha256(Path(route.__file__).read_bytes()).hexdigest(),
        'engine_sha256':hashlib.sha256(Path(v.__file__).read_bytes()).hexdigest()}
destination=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).with_name('FABRICATED_TESTS.json')
destination.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps({key:result[key] for key in ('status','checks_passed','rejected_controls_count','physical_source_evaluations','archive_arrays_decoded')}))
