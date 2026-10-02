"""Pre-freeze synthetic algebra/format checks. No registered source/mode evaluation."""
import ast
import math
from pathlib import Path
import tempfile
import numpy as np
import validate_stress as v
import raw_stress_audit as r

checks=[]
def check(name,condition):
    v.require(condition,name);checks.append(name)

def rejects(name,function):
    try:function()
    except RuntimeError:checks.append(name)
    else:raise RuntimeError('Mutation survived: '+name)

rejects('explicit guards remain enabled',lambda:v.require(False,'synthetic false'))
rejects('scalar NaN rejection',lambda:v.close(float('nan'),0,1,'nan'))
rejects('scalar infinity rejection',lambda:v.close(float('inf'),0,1,'infinite'))
rejects('changed scalar rejection',lambda:v.close(1,2,.1,'changed'))
rejects('changed raw shape',lambda:r.near_array(np.zeros(2),np.zeros(3),'shape'))
rejects('raw NaN rejection',lambda:r.near_array(np.array([float('nan')]),np.zeros(1),'nonfinite'))
rejects('raw alteration rejection',lambda:r.near_array(np.ones(1),np.zeros(1),'altered'))
rejects('odd Simpson endpoint rejection',lambda:r.simpson(np.zeros(4),3,.1))
# Synthetic polynomial ledger, unrelated to either registered pulse.
x=np.arange(9,dtype=np.longdouble)/8
actual=r.simpson(x**3,8,np.longdouble(1)/8)
check('Simpson cubic identity',abs(actual-np.longdouble(1)/4)<1e-18)
# Arbitrary synthetic variables at an unregistered time; no response/source eval.
eta=-2.0; f=[.7,-.2,.3,0.,0.,0.];values={'q':.11,'q_prime':-.07,'q_second':.09}
c=v.closed(eta,f,values,10);a=L=.5
trace=values['q_second']/2-L*values['q_prime']-3*L*L*values['q']-a*a*c['Q0']*f[0]+c['anomaly']
check('independent finite-contact trace identity',abs(-c['rho']+3*c['p']-trace)<1e-15)
check('quadratic current contact included',c['current']==values['q']+c['Q0']*f[0]/4)
records={};v.mutation(records,'synthetic tiny wrong contact',1e-14,2e-11,'algebraic')
check('small-contact sensitivity honestly recorded',records['synthetic tiny wrong contact']['maximum_absolute_residual']==1e-14 and not records['synthetic tiny wrong contact']['operational_gate_exceeded_somewhere'])
v.guard_control(records,'synthetic failed state',lambda:v.require(False,'state'))
check('format control explicitly classified',records['synthetic failed state']['rejected'] and 'no altered-state physical run' in records['synthetic failed state']['kind'])
for name in ('validate_stress.py','raw_stress_audit.py'):
    tree=ast.parse((Path(__file__).parent/name).read_text())
    check('no assertions '+name,not any(isinstance(n,ast.Assert) for n in ast.walk(tree)))
    check('no producer imports '+name,not any(isinstance(n,ast.ImportFrom) and n.module in ('stress_primary','forced_stress','source_jet') for n in ast.walk(tree)))
with tempfile.TemporaryDirectory(prefix='stress-guard-',dir='/tmp') as directory:
    existing=Path(directory)
    rejects('existing output directory rejected',lambda:v.validate_output_path(existing))
    file=existing/'existing.json';file.write_text('{}')
    rejects('existing output file rejected',lambda:v.validate_output_path(file))
    rejects('in-checkpoint output rejected',lambda:v.validate_output_path(v.ROOT/'not-created-validator-guard'))
    check('fresh external output accepted',v.validate_output_path(existing/'fresh-output')==existing/'fresh-output')
print('PASS',len(checks),'synthetic algebra/format/static checks; no registered numerical evaluations')
