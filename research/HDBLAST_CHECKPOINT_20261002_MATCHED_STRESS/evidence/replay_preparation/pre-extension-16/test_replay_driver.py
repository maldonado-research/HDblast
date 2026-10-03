"""Static/synthetic replay checks. Executes no child command or physical model."""
import ast
import copy
import importlib.util
from pathlib import Path
import tempfile
spec=importlib.util.spec_from_file_location('replay',Path(__file__).with_name('replay_stress.py'))
r=importlib.util.module_from_spec(spec);spec.loader.exec_module(r)
checks=[]
def check(name,condition):r.require(condition,name);checks.append(name)
def rejects(name,function):
    try:function()
    except RuntimeError:checks.append(name)
    else:raise RuntimeError('Mutation survived: '+name)
plan=r.command_plan(Path('/synthetic/checkpoint'),Path('/synthetic/fresh'),Path('/synthetic/mirror'))
check('sixteen exact commands',len(plan)==16 and len({n for n,c in plan})==16)
check('one primary invocation',sum(n=='primary' for n,c in plan)==1)
check('one independent invocation',sum(n=='independent_modes' for n,c in plan)==1)
check('seven optimized checks',sum(c[0]=='-O' for n,c in plan)==7)
check('validators require public freeze',all('--public-freeze-commit' in c and r.FREEZE in c for n,c in plan if n.startswith('validation')))
check('validator outputs are directories',all(Path(c[c.index('--output')+1]).suffix=='' for n,c in plan if n.startswith('validation')))
check('tail proof uses self-contained mirror',all(c[c.index('--repo-root')+1]=='/synthetic/mirror' for n,c in plan if n.startswith('stress_tail_algebra')))
controls=[{'rejected':True} for i in range(5)]+[{'maximum_absolute_residual':1,'nonzero_algebraic_sensitivity_witness':True} for i in range(15)]
normal={'status':'PASS','counts':r.EXPECTED_COUNTS,'controls':controls,'python_optimization':0}
optimized={**normal,'python_optimization':1};r.compare_checks(normal,optimized);checks.append('normal optimized comparison accepts expected evidence')
bad=copy.deepcopy(optimized);bad['counts']['core_comparisons']=179
rejects('changed comparison count rejected',lambda:r.compare_checks(normal,bad))
bad=copy.deepcopy(optimized);bad['python_optimization']=0
rejects('wrong optimization mode rejected',lambda:r.compare_checks(normal,bad))
bad=copy.deepcopy(optimized);bad['controls'][0]['rejected']=False
rejects('control result mismatch rejected',lambda:r.compare_checks(normal,bad))
with tempfile.TemporaryDirectory(dir='/tmp',prefix='stress-replay-guard-') as directory:
    root=Path(directory);(root/'regular').write_text('synthetic')
    check('regular input accepted',r.safe_path(root,'regular')==root/'regular')
    rejects('relative input escape rejected',lambda:r.safe_path(root,'../outside'))
    rejects('absolute input rejected',lambda:r.safe_path(root,'/tmp/elsewhere'))
    (root/'link').symlink_to(root/'regular')
    rejects('symlink input rejected',lambda:r.safe_path(root,'link'))
tree=ast.parse(Path(r.__file__).read_text())
check('explicit exceptions under optimization',not any(isinstance(n,ast.Assert) for n in ast.walk(tree)))
check('no numerical producer imports',not any(isinstance(n,ast.ImportFrom) and n.module in ('stress_primary','forced_stress') for n in ast.walk(tree)))
print('PASS',len(checks),'static/synthetic replay checks; zero child commands executed')
