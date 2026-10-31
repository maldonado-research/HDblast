#!/usr/bin/env python3
"""Bounded exact-packet initial experiment; no response evolution.

All decisions are rational; exponential intervals enclose the one exact
analytic function. Shared packet and E differences are identically zero.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
FREEZE_PIN = 'aa9a7eae96aae988d3d1a25d83dcd6b436893a9464f7e3813436472a8e137ec1'
CORE = 'research/HDBLAST_HIGHER_INITIAL_JETS_20261030/checks/initial_second_core.py'
PROFILES = 'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026/initial/results/INITIAL_PROFILES.json'
PROFILES_PIN = '555f911f6c336862dfacb6aa56ef5179ac06701812d593b5bfbbf6ea9bd4d5d1'

def need(ok, text):
    if not ok:
        raise RuntimeError(text)

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def encode(x):
    if isinstance(x, Q): return str(x)
    if isinstance(x, dict): return {k: encode(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [encode(v) for v in x]
    return x

def sqrt_bound(q, bits, upper=True):
    need(q >= 0, 'negative square')
    den = 1 << bits
    k = math.isqrt(q.numerator*den*den//q.denominator)
    if upper and Q(k*k, den*den) < q: k += 1
    return Q(k, den)

def cover(rows, lo, hi):
    return bool(rows) and rows[0][0] == lo and rows[-1][1] == hi and all(a < b for a,b in rows) and all(a[1] == b[0] for a,b in zip(rows, rows[1:]))

def authenticate(root, frozen):
    result = {}
    for rel, pin in frozen['input_sha256'].items():
        path = root/rel
        need(path.is_file() and not path.is_symlink() and digest(path) == pin, 'source hash '+rel)
        result[rel] = pin
    need(digest(root/PROFILES) == PROFILES_PIN, 'velocity profile hash')
    result[PROFILES] = PROFILES_PIN
    # Authenticate every member of every inherited manifest before imports.
    for rel in sorted(frozen['input_sha256']):
        if not rel.endswith('/MANIFEST.json'): continue
        path = root/rel
        manifest = json.loads(path.read_bytes())
        for row in manifest['files']:
            name = Path(row['path'])
            need(not name.is_absolute() and '..' not in name.parts, 'manifest member path')
            member = path.parent/name
            need(member.is_file() and not member.is_symlink(), 'manifest member type '+str(member))
            need(member.stat().st_size == row['bytes'] and digest(member) == row['sha256'], 'manifest member hash '+str(member))
    return result

class ExactPacket:
    def __init__(self, root, bits):
        spec = importlib.util.spec_from_file_location('exact_packet_authenticated_core', root/CORE)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self.h = module.HigherInitial(root)
        self.c = self.h.ctx
        self.I, self.J, self.P = self.h.I, self.h.J, self.h.P
        self.ell = self.h.ell
        self.bits = bits
        # Context's original convention check happens before this explicit
        # registered refinement of the outward grid. No input bytes change.
        self.c.init.GRID = 1 << bits
        self.c.init.expminus = self.expminus
        self.c.init.expbox = self.expbox
        self.weights = list(self.c.weights)

    def expminus(self, q):
        need(q >= 0, 'expminus nonnegative argument')
        if q == 0: return self.I(1)
        t, halvings = q, 0
        while t > 1: t /= 2; halvings += 1
        degree = 2
        while Q(1, math.factorial(degree+1)) >= Q(1, 1 << (self.bits+16)):
            degree += 2
        upper = sum(((-t)**k/math.factorial(k) for k in range(degree+1)), Q(0))
        lower = upper-t**(degree+1)/math.factorial(degree+1)
        need(0 < lower <= upper, 'alternating exponential enclosure')
        out = self.I(lower, upper)
        for _ in range(halvings): out = out*out
        return out

    def expbox(self, q):
        return self.I(self.expminus(q.hi).lo, self.expminus(q.lo).hi)

    def norm(self, caps):
        return sqrt_bound(sum(w*x*x for w,x in zip(self.weights,caps)), self.bits)

    def packet(self, lo, hi):
        if hi <= 1: return self.h.constant(0)
        need(lo >= 1 and hi <= Q(7,4), 'packet footprint')
        return self.h.bump(self.I(lo,hi), lo == 1)

    def evaluate(self, lo, hi):
        bg = self.h.background(self.I(lo,hi))
        q = self.packet(lo,hi)
        E0 = self.expbox(self.ell*self.I(lo,hi))
        E = self.J([(-1)**k*E0/math.factorial(k) for k in range(6)])
        zero = self.h.constant(0)
        shared_q = self.P(q,q,[zero,zero,zero])
        shared_E = self.P(E,E,[zero,zero,zero])
        u, unused_reconstructed_velocity = self.c._registered(bg,shared_E,shared_q)
        need(all(z.n >= 2 for z in u), 'completed displacement through rr')
        scales = [Q(1),self.ell,2*self.ell**2]
        actual = [[s*z.a.c[k] for z in u] for k,s in enumerate(scales)]
        nominal = [[s*z.b.c[k] for z in u] for k,s in enumerate(scales)]
        channels = [[[s*z.d[ch].c[k] for z in u] for k,s in enumerate(scales)] for ch in range(3)]
        need(all(v.lo == v.hi == 0 for row in channels[1] for v in row), 'q and E shared functional differences exactly zero')
        return actual,nominal,channels

    def denominators(self):
        bg = self.h.background(self.I(0,Q(7,4)))
        rows = {}
        for key in ('rho','p'):
            a,b = bg[key].a.val(),bg[key].b.val()
            need((a.lo > 0 or a.hi < 0) and (b.lo > 0 or b.hi < 0), 'uniform actual and nominal '+key+' denominator separation')
            rows[key] = {'actual':a.row(),'nominal':b.row(),'nominal_distance_from_zero':min(abs(b.lo),abs(b.hi))}
        return rows

def velocity(root, evaluator, old):
    p = json.loads((root/PROFILES).read_bytes())
    nodes = list(map(Q,p['nodes']))
    vt = [list(map(Q,row)) for row in p['profiles']['v']]
    need(p['labels'] == list(evaluator.c.op.LABELS) and list(map(Q,p['gram'])) == evaluator.weights, 'velocity full23 Gram')
    need(nodes[0] == 0 and nodes[-1] == Q(7,4) and all(a < b for a,b in zip(nodes,nodes[1:])), 'velocity complete nodes')
    need(all(z == 0 for row in vt[:2] for z in row) and nodes[1] == Q(65,64), 'velocity flat collar')
    slopes = [[(v1[i]-v0[i])/(b-a) for i in range(23)] for a,b,v0,v1 in zip(nodes,nodes[1:],vt,vt[1:])]
    nodecap = max(evaluator.norm(row) for row in vt)
    derivative = max(evaluator.norm(row) for row in slopes)
    width = Q(1,2048)
    need(Q(65,64)-width > 1, 'unchanged velocity vanishes on [0,1] after convolution')
    return {'absolute_H_caps':{'v':nodecap+derivative*width*Q(35,128),'v_r':derivative},
        'error_H_caps':{'v0':old['initial_error_H_caps']['v0'],'v0_r':old['initial_error_H_caps']['v0_r']},
        'profile_sha256':PROFILES_PIN,'delta_v':width,'P1_exterior':'last cell slope continued; H-isometric parity reflection',
        'kernel':'35/(32*delta_v)*(1-(s/delta_v)^2)^3 on |s|<=delta_v; zero outside',
        'kernel_mass':1,'kernel_absolute_first_moment_per_width':Q(35,128),
        'regularity':'P1 is continuous and Lipschitz; compact C2 kernel implies smoothed velocity C2, hence required C1; differentiation uses locally integrable P1 slopes',
        'error_ledger_retained_exactly':True,'support_zero_through':Q(65,64)-width}

def local(e, joins):
    a,n,d = e.evaluate(Q(65,64),Q(65,64))
    errors = [d[0][2][i]+d[2][2][i] for i in range(23)]
    rows = []
    lower2,upper2,error2 = Q(0),Q(0),Q(0)
    for label,w,actual,candidate,error in zip(e.c.op.LABELS,e.weights,a[2],n[2],errors):
        low = 0 if actual.lo <= 0 <= actual.hi else min(abs(actual.lo),abs(actual.hi))
        lower2 += w*low*low; upper2 += w*actual.mag()**2; error2 += w*error.mag()**2
        rows.append({'label':label,'Gram':w,'actual_u_rr':actual.row(),'candidate_u_rr':candidate.row(),
            'paired_actual_minus_candidate_u_rr':error.row(),'actual_H_squared_lower_contribution':w*low*low,
            'actual_H_squared_upper_contribution':w*actual.mag()**2,'error_H_squared_upper_contribution':w*error.mag()**2})
    denlo,denup,num = sqrt_bound(lower2,e.bits,False),sqrt_bound(upper2,e.bits),sqrt_bound(error2,e.bits)
    oldlower = Q(joins['joins'][0]['fixed_width_error_H_norm_strict_lower'])
    need(Q(joins['joins'][0]['r']) == Q(65,64), 'old local witness same point')
    need(denlo > 0, 'local actual H denominator positive')
    return {'point':[Q(0),Q(65,64)],'rows':rows,'actual_H_norm_interval':[denlo,denup],
        'new_error_H_upper':num,'old_error_H_strict_lower':oldlower,'old_error_relative_to_current_target_H_strict_lower':oldlower/denup,
        'new_relative_error_H_upper':num/denlo,'relative_gate':num <= denlo/10,'true_error_improvement_gate':num < oldlower/2,
        'pass':num <= denlo/10 and num < oldlower/2}

def campaign(root,bits,old,combined,joins):
    e = ExactPacket(root,bits)
    den = e.denominators()
    vel = velocity(root,e,old)
    loc = local(e,joins)
    limits = [Q(old['initial_error_H_caps']['e0']),Q(old['initial_error_H_caps']['e0_r']),Q(combined['refined_uniform_initial_u_rr_H_upper'])/2]
    seeds = [(Q(0),Q(1)),(Q(1),Q(65,64))]+[(Q(k,64),Q(k+1,64)) for k in range(65,112)]
    need(cover(seeds,Q(0),Q(7,4)), 'seed exact cover')
    leaves,splits = [],[]
    stack = [(lo,hi,0,i,None) for i,(lo,hi) in reversed(list(enumerate(seeds)))]
    while stack:
        lo,hi,depth,seed,parent = stack.pop()
        if hi <= 1:
            actual,nominal,channels = e.evaluate(lo,hi)
        else:
            actual,nominal,channels = e.evaluate(lo,hi)
        absolute = [[v.mag() for v in row] for row in nominal]
        raw = [[[v.mag() for v in row] for row in channel] for channel in channels]
        caps = raw if parent is None else [[[min(x,y) for x,y in zip(row,oldrow)] for row,oldrow in zip(ch,oldch)] for ch,oldch in zip(raw,parent)]
        norms = [[e.norm(row) for row in ch] for ch in caps]
        missed = any(norms[ch][k] > limits[k]/2 for ch in (0,2) for k in range(3))
        if missed and depth < 16 and len(seeds)+len(splits) < 4096:
            mid = (lo+hi)/2
            splits.append({'r':[lo,hi],'depth':depth,'seed':seed,'midpoint':mid,'component_caps':caps})
            stack.append((mid,hi,depth+1,seed,caps));stack.append((lo,mid,depth+1,seed,caps))
        else:
            leaves.append({'r':[lo,hi],'depth':depth,'seed':seed,'absolute_component_caps':absolute,
                'absolute_H_caps':[e.norm(row) for row in absolute],'paired_channel_component_caps':caps,
                'paired_channel_H_caps':norms,'threshold_satisfied':not missed,
                'stop_reason':'threshold' if not missed else ('depth_budget' if depth == 16 else 'leaf_budget')})
    need(cover([z['r'] for z in leaves],Q(0),Q(7,4)), 'complete exact leaf cover')
    need(len(leaves) == len(seeds)+len(splits) and len(leaves) <= 4096 and max(z['depth'] for z in leaves) <= 16, 'frozen budget respected')
    channel_H = [[max(z['paired_channel_H_caps'][ch][k] for z in leaves) for k in range(3)] for ch in range(3)]
    errors = [sum(channel_H[ch][k] for ch in range(3)) for k in range(3)]
    absolute = [max(z['absolute_H_caps'][k] for z in leaves) for k in range(3)]
    regions = []
    for name,lo,hi in [('exact_zero',Q(0),Q(1)),('flat_edge_to_old_cutoff',Q(1),Q(65,64)),('former_cutoff_band',Q(65,64),Q(66,64)),('active_interior',Q(66,64),Q(7,4))]:
        rows = [z for z in leaves if lo <= z['r'][0] and z['r'][1] <= hi]
        need(cover([z['r'] for z in rows],lo,hi), 'region exact cover '+name)
        regions.append({'name':name,'r':[lo,hi],'leaf_count':len(rows),'rr_absolute_error_H_upper':sum(max(z['paired_channel_H_caps'][ch][2] for z in rows) for ch in range(3)),
            'relative_ratio':'undefined at flat support edge; unresolved elsewhere without a uniform positive target denominator'})
    ea,en,ed = e.evaluate(Q(7,4),Q(7,4))
    endpoint = []
    for i,label in enumerate(e.c.op.LABELS):
        endpoint.append({'label':label,'interior_exact_jet_symbols':['f_'+str(i)+'(R)','f_'+str(i)+"_r(R)",'f_'+str(i)+'_rr(R)'],
            'candidate_trace_enclosures':[en[k][i].row() for k in range(3)],'exterior_trace_minus_interior_exact':[0,0,0],
            'paired_error_enclosures':[(ed[0][k][i]+ed[2][k][i]).row() for k in range(3)]})
    regions.append({'name':'outer_endpoint','r':[Q(7,4),Q(7,4)],'rr_absolute_error_H_upper':e.norm([(ed[0][2][i]+ed[2][2][i]).mag() for i in range(23)]),'relative_ratio':'unresolved; absolute bound reported'})
    initialpass = all(x <= y for x,y in zip(errors,limits))
    e.h.authenticate()
    return {'bits':bits,'status':'INITIAL_AND_LOCAL_GATES_PASS' if loc['pass'] and initialpass else 'INCONCLUSIVE_INITIAL_OR_LOCAL_GATES',
        'labels':list(e.c.op.LABELS),'Gram':e.weights,'denominators':den,'local':loc,'velocity':vel,
        'initial_error_H_caps':dict(zip(['e0','e0_r','e0_rr'],errors)) | vel['error_H_caps'],
        'candidate_initial_absolute_H_caps':dict(zip(['u','u_r','u_rr'],absolute)) | vel['absolute_H_caps'],
        'initial_gate_limits':dict(zip(['e0','e0_r','e0_rr'],limits)),'initial_gate_pass':initialpass,
        'paired_channel_names':['background_coefficient_difference','exact_shared_q_E_zero','actual_ODE_remainder'],
        'paired_channel_H_caps':channel_H,'regions':regions,'endpoint_componentwise_quadratic_matching':endpoint,
        'seeds':seeds,'leaf_count':len(leaves),'split_count':len(splits),'leaves':leaves,'splits':splits,
        'unresolved_cells':[z['r'] for z in leaves if not z['threshold_satisfied']],
        'budget':{'bits':bits,'max_leaves':4096,'max_depth':16},'no_response_evolution_performed':True}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo-root',type=Path,required=True)
    ap.add_argument('--output-dir',type=Path,required=True)
    args = ap.parse_args()
    need(not args.output_dir.exists(), 'fresh absent output path required')
    self_before = digest(Path(__file__))
    need(digest(HERE/'FREEZE.json') == FREEZE_PIN, 'frozen definition hash')
    frozen = json.loads((HERE/'FREEZE.json').read_bytes())
    before = authenticate(args.repo_root,frozen)
    base = args.repo_root/'research/HDBLAST_REFINED_ERROR_TRANSFER_20261030/reference'
    old = json.loads((base/'C2.json').read_bytes());combined = json.loads((base/'COMBINED.json').read_bytes())
    joins = json.loads((args.repo_root/'research/HDBLAST_INITIAL_CURVATURE_REFINEMENT_20261030/reference/JOINS.json').read_bytes())
    reports = []
    for bits in (128,192):
        report = campaign(args.repo_root,bits,old,combined,joins)
        reports.append(report)
        if report['status'] == 'INITIAL_AND_LOCAL_GATES_PASS': break
    need(before == authenticate(args.repo_root,frozen), 'all input files unchanged after execution')
    need(digest(Path(__file__)) == self_before and digest(HERE/'FREEZE.json') == FREEZE_PIN, 'producer and frozen map unchanged after execution')
    result = {'schema':'HDBLAST_EXACT_PACKET_INITIAL_EXPERIMENT_V1','candidate':'W1_exact_packet',
        'status':reports[-1]['status'],'freeze_sha256':FREEZE_PIN,'producer_sha256':self_before,
        'authentication_before_after':before,'campaigns':reports,'target_and_original_F_Fx_G_M_unchanged':True,
        'response_evolution_performed':False,'exact_shared_cancellation':'q_actual=q_candidate and E_actual=E_candidate as analytic functions; their paired difference jets are exactly zero',
        'regularity':'completed displacement C2 across componentwise quadratic outer extension; no C3 claim; flat at1. Unchanged smoothed P1 velocity C2, required C1.',
        'scope':frozen['scope']}
    args.output_dir.mkdir(parents=True)
    (args.output_dir/'INITIAL.json').write_text(json.dumps(encode(result),sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'campaigns':[r['bits'] for r in reports],'leaf_counts':[r['leaf_count'] for r in reports],
        'initial_error_H_caps':{k:float(Q(v)) for k,v in reports[-1]['initial_error_H_caps'].items()},
        'local_relative_upper':float(reports[-1]['local']['new_relative_error_H_upper'])},sort_keys=True))

if __name__ == '__main__': main()
