#!/usr/bin/env python3
"""Prospectively fixed saved-mode ledger diagnostic. Never evolves a source."""
from __future__ import annotations
import argparse
from decimal import Decimal
from fractions import Fraction
import hashlib
import importlib.util
import io
import json
import math
from pathlib import Path
import resource
import re
import sys
import time
import zipfile
import mpmath
import numpy as np

LD = np.longdouble
SOURCES = ('positive_B', 'signed_uB')
SETTINGS = ('coarse', 'fine')
CUTOFFS = (64, 128, 256)
LEVELS = ((80, 266), (100, 333))
EPS_RATIO = (3777893186295716171, 37778931862957161709568)
PI_RATIO = (14488038916154245685, 4611686018427387904)
MEMBERS = ('k', 'momentum_weights', 'u_4', 'w_4', 'u_5', 'w_5',
           'observation_eta', 'history_eta', 'history_values',
           'history_ledger_integrand', 'history_baseline_contact',
           'history_source_jet', 'history_forcing_jet', 'history_geometry',
           'history_cutoffs', 'history_quantity_names', 'history_geometry_names')
QUANTITIES = ('q','q_prime','q_second','rho','p','Q0','rho0','p0','current')
GEOMETRY = ('a0','L','L_prime','L_second','L_third','M','M_prime','M_second','M_third','M_fourth')
GATES = {'profile':'0.0000002', 'D_cont':'0.0000002', 'E_flow':'0.0000002',
         'decomposition':'0.000000000001', 'precision_gap':'0.000000000001',
         'serialization':'0.000000000001', 'attribution':'0.000002'}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(65536), b''):
            h.update(block)
    return h.hexdigest()


def authenticated_helper(root, registration_sha256, freeze_commit):
    """Pin the verifier BEFORE importing it; the outer driver pins this entry."""
    require(isinstance(registration_sha256,str) and re.fullmatch(r'[0-9a-f]{64}',registration_sha256),
            'Explicit canonical registration SHA256 required')
    require(isinstance(freeze_commit,str) and re.fullmatch(r'[0-9a-f]{40}',freeze_commit),
            'Explicit complete public freeze commit required')
    root=Path(root).absolute()
    for p in (root,*root.parents):
        require(not p.is_symlink(),'Symlink checkpoint ancestry rejected')
    registration_path=root/'FULL_REGISTRATION.json'
    helper_path=root/'code'/'ledger_integrity.py'
    source_path=root/'code'/'diagnostic_primary.py'
    for p in (registration_path,helper_path,source_path,root/'code'):
        require(not p.is_symlink(),'Symlink bootstrap source rejected')
    require(registration_path.is_file() and sha(registration_path)==registration_sha256,
            'Bootstrap registration byte hash mismatch')
    registration=json.loads(registration_path.read_text())
    require(registration.get('schema_version')==1 and isinstance(registration.get('files'),dict),
            'Bootstrap registration schema mismatch')
    for relative,path in (('code/ledger_integrity.py',helper_path),('code/diagnostic_primary.py',source_path)):
        require(path.is_file() and sha(path)==registration['files'].get(relative),
                'Bootstrap registered source hash mismatch '+relative)
    require(Path(__file__).absolute()==source_path,'Execute the registered checkpoint entry point')
    spec=importlib.util.spec_from_file_location('ledger_integrity',helper_path)
    require(spec is not None and spec.loader is not None,'Verifier loader unavailable')
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.verify_frozen


def fraction_of_mpf(x):
    sign, man, exp, _ = x._mpf_
    numerator = -man if sign else man
    return Fraction(numerator << exp, 1) if exp >= 0 else Fraction(numerator, 1 << -exp)


def represented(ctx, value):
    """Exact LD binary ratio, then one operation in the declared MP context."""
    require(isinstance(value,np.longdouble) and np.isfinite(value), 'Native finite binary80 scalar required before conversion')
    n, d = value.as_integer_ratio()
    require(d > 0 and d & (d-1) == 0, 'Binary80 ratio denominator is not a power of two')
    if n:
        magnitude=abs(n)
        trailing=(magnitude & -magnitude).bit_length()-1
        mantissa=magnitude >> trailing
        require(mantissa.bit_length() <= 64, 'Binary80 mantissa exceeds64 bits')
        restored=np.ldexp(LD(-mantissa if n<0 else mantissa),trailing-(d.bit_length()-1))
    else:
        restored=LD('-0') if np.signbit(value) else LD('0')
    require(restored == value and np.signbit(restored) == np.signbit(value), 'LD exact-ratio/sign roundtrip failed')
    return ctx.mpf(n)/ctx.mpf(d)


def complex_represented(ctx, value):
    return ctx.mpc(represented(ctx, value.real), represented(ctx, value.imag))


def round_shift(n, p):
    """Round signed integer / 2**p to nearest integer, ties to even."""
    q = n >> p
    remainder = n-(q << p)
    half = 1 << (p-1)
    return q + int(remainder > half or (remainder == half and q & 1))


def quantize(x, p):
    sign, man, exp, _ = x._mpf_
    n = -man if sign else man
    exponent = exp+p
    return n << exponent if exponent >= 0 else round_shift(n, -exponent)


def dyadic(ctx, n, p):
    return ctx.mpf(n)/ctx.mpf(1 << p)


def simpson_to(history, end, dt):
    require(end % 2 == 0 and end >= 0, 'Simpson endpoint is not an even global index')
    return dt/3*(history[0]+history[end]
                 +4*np.sum(history[1:end:2],axis=0,dtype=LD)
                 +2*np.sum(history[2:end:2],axis=0,dtype=LD))


def input_metadata(root):
    root = Path(root)
    manifest = json.loads((root/'INPUT_MANIFEST.json').read_text())
    require(manifest['schema_version'] == 1, 'Unknown input manifest schema')
    require(manifest['scientific_calibration_status'] == 'FAIL_UNCHANGED', 'Prior failure status lost')
    require(tuple(manifest['source_constants']['epsilon']['exact_represented_binary80_ratio']) == EPS_RATIO,
            'Wrong represented epsilon')
    require(tuple(manifest['source_constants']['pi']['exact_represented_binary80_ratio']) == PI_RATIO,
            'Wrong represented pi')
    require(len(manifest['capsules']) == 4, 'Exactly four immutable capsules required')
    caps = {}
    for source in SOURCES:
        for setting in SETTINGS:
            name = f'metric_modes_{source}_{setting}.npz'
            matches = [c for c in manifest['capsules'] if c['capsule_path'] == name]
            require(len(matches) == 1, 'Missing or duplicated capsule '+name)
            cap = matches[0]
            require(len(cap['members']) == len(MEMBERS) and
                    {m['name'] for m in cap['members']} == {n+'.npy' for n in MEMBERS},
                    'Incomplete or duplicated declared member manifest')
            path = root/name
            require(not path.is_symlink() and path.is_file(), 'Capsule path must be a regular file')
            require(path.stat().st_size == cap['capsule_bytes'] and sha(path) == cap['capsule_sha256'],
                    'Capsule hash mismatch '+name)
            with zipfile.ZipFile(path) as z:
                require(len(z.namelist()) == len(MEMBERS) and set(z.namelist()) == {n+'.npy' for n in MEMBERS}, 'Wrong capsule member set')
                require(z.testzip() is None, 'Capsule CRC failed')
                for m in cap['members']:
                    raw = z.read(m['name'])
                    require(len(raw) == m['bytes'] and hashlib.sha256(raw).hexdigest() == m['sha256'],
                            'Member byte mismatch '+m['name'])
                    stream = io.BytesIO(raw)
                    version = np.lib.format.read_magic(stream)
                    shape, order, dtype = np.lib.format._read_array_header(stream,version)
                    require(list(shape) == m['header']['shape'] and dtype.str == m['header']['dtype_descriptor']
                            and bool(order) == m['header']['fortran_order'], 'NPY header mismatch')
            caps[(source,setting)] = path
    require(np.finfo(LD).nmant == 63 and np.dtype(LD).itemsize == 16, 'Producer binary80/x86-extended platform required')
    require(LD('0.0001').as_integer_ratio() == EPS_RATIO, 'Runtime epsilon differs')
    require(np.arccos(LD(-1)).as_integer_ratio() == PI_RATIO, 'Runtime pi differs')
    return manifest, caps


def validate_arrays(a, setting):
    n = 8192 if setting == 'coarse' else 16384
    count = 577 if setting == 'coarse' else 1153
    samples = 129 if setting == 'coarse' else 257
    require(set(a) == set(MEMBERS), 'Loaded member set differs')
    for key in ('k','momentum_weights'):
        require(a[key].shape == (n,) and a[key].dtype == np.dtype(LD), 'Wrong node/weight shape or precision')
        require(np.all(np.isfinite(a[key])) and np.all(a[key] > 0), 'Nonpositive/nonfinite nodes/weights')
    require(np.all(np.diff(a['k']) > 0) and a['k'][-1] < 256, 'Momentum ordering/range failed')
    for key in ('u_4','w_4','u_5','w_5'):
        require(a[key].shape == (n,) and a[key].dtype == np.dtype(np.clongdouble)
                and np.all(np.isfinite(a[key])), 'Wrong retained complex mode schema')
    require(tuple(a['history_quantity_names']) == QUANTITIES, 'Quantity labels differ')
    require(tuple(a['history_geometry_names']) == GEOMETRY, 'Geometry labels differ')
    require(tuple(a['history_cutoffs']) == CUTOFFS, 'Cutoff labels differ')
    require(a['history_eta'].shape == (count,) and a['history_eta'].dtype == np.dtype(LD), 'Wrong history times')
    dt = LD(1)/(128 if setting == 'coarse' else 256)
    require(np.array_equal(a['history_eta'], LD(-6)+np.arange(count,dtype=LD)*dt), 'Wrong global history grid')
    require(tuple(a['observation_eta']) == (-5.5,-4.5,-4,-3.5,-2.5,-1.5), 'Wrong observations')
    shapes = {'history_values':(count,3,9), 'history_ledger_integrand':(count,3),
              'history_baseline_contact':(count,3),'history_source_jet':(count,6),
              'history_forcing_jet':(count,4),'history_geometry':(count,10)}
    for key, shape in shapes.items():
        require(a[key].shape == shape and a[key].dtype == np.dtype(LD)
                and np.all(np.isfinite(a[key])), 'Wrong/nonfinite history array '+key)
    indices = np.flatnonzero((a['history_eta'] >= LD('-2.5')) & (a['history_eta'] <= LD('-1.5')))
    require(len(indices) == samples and np.all(np.diff(indices) == 1), 'Wrong inclusive profile interval')
    ia, ib = int(indices[0]), int(indices[-1])
    for key in ('history_source_jet','history_forcing_jet','history_baseline_contact'):
        require(np.count_nonzero(a[key][ia:ib+1]) == 0, 'Source/contact does not vanish on selected interval')
    expected_nodes = (2048,4096,8192) if setting == 'coarse' else (4096,8192,16384)
    require(tuple(int(np.count_nonzero(a['k'] < cutoff)) for cutoff in CUTOFFS) == expected_nodes,
            'Fixed momentum prefix sizes differ')
    return indices, dt


def moment_kernel(a, setting, dps, p, started=None):
    """Quantized three-moment recurrence, exact integer accumulation in three bins."""
    ctx = mpmath.mp.clone()
    ctx.dps = dps
    eps = ctx.mpf(EPS_RATIO[0])/EPS_RATIO[1]
    pi = ctx.mpf(PI_RATIO[0])/PI_RATIO[1]
    indices, dt_ld = validate_arrays(a,setting)
    times = [represented(ctx,t) for t in a['history_eta'][indices]]
    dt = represented(ctx,dt_ld)
    La, Lb = -1/times[0], -1/times[-1]
    length = len(times)
    bins = [[[0]*length for _ in range(3)] for _ in range(3)]
    constants = [[ctx.zero,ctx.zero] for _ in range(3)]
    flows = [[ctx.zero,ctx.zero,ctx.zero] for _ in range(3)]
    # Upper bounds relative to exact phases of the rounded work-context seed/rotation.
    bounds = [ctx.zero]*3
    endpoint_checks = [ctx.zero]*3
    for i in range(len(a['k'])):
        if started is not None and i % 256 == 0:
            require(time.monotonic()-started <= 900, 'Shared twelve-case wall budget exceeded')
            require(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss <= 262144, 'RSS budget exceeded')
        k, weight = represented(ctx,a['k'][i]), represented(ctx,a['momentum_weights'][i])
        omega = 2*k
        mu = weight*k*k/(2*pi*pi)
        ua, wa = complex_represented(ctx,a['u_4'][i]), complex_represented(ctx,a['w_4'][i])
        ub, wb = complex_represented(ctx,a['u_5'][i]), complex_represented(ctx,a['w_5'][i])
        d = wa/(ctx.j*omega)
        c = ua-d  # No Wronskian/normalization projection.
        alpha = mu*d/(2*k*eps)
        b = 0 if k < 64 else (1 if k < 128 else 2)
        constants[b][0] += mu*c.real/(2*k*eps)
        constants[b][1] += mu*k*c.real/eps
        rotation = ctx.exp(ctx.j*omega*dt)
        er, ei = quantize(rotation.real,p), quantize(rotation.imag,p)
        zr, zi = quantize(alpha.real,p), quantize(alpha.imag,p)
        omega_i, kappa_i = quantize(omega,p), quantize(4*k*k/3,p)
        # Exact dyadic magnitude upper bounds, then a conservative MP expression.
        mag_seed_i = math.isqrt(zr*zr+zi*zi)+1
        mag_rotation_i = math.isqrt(er*er+ei*ei)+1
        seed_mag = dyadic(ctx,mag_seed_i,p)
        rho = max(ctx.one,dyadic(ctx,mag_rotation_i,p))
        gamma = ctx.mpf(2)/(1 << p) # exceeds sqrt(2)/(2Q), seed/step complex rounding.
        # Exact rational envelope against the work-context seed/rotation treated
        # as exact dyadics. No claim about their transcendental/context errors.
        Q=1 << p
        nsteps=length-1
        drift=max(Q,mag_rotation_i+2)-Q
        require(nsteps*drift<Q,'Quantized phase norm bound is outside the frozen regime')
        growth=Fraction(Q,Q-nsteps*drift) # (1+drift/Q)^n <= 1/(1-n*drift/Q).
        alpha_bound=Fraction(mag_seed_i+2,Q)
        e_bound=Fraction(2,Q)*(1+nsteps*(alpha_bound+1))*growth
        actual_bound=alpha_bound*growth+e_bound
        eb=abs(fraction_of_mpf(omega))*e_bound+Fraction(2,Q)*(actual_bound+1)
        eu=abs(fraction_of_mpf(4*k*k/3))*e_bound+Fraction(2,Q)*(actual_bound+1)
        rp_bound=4*fraction_of_mpf(Lb)**2*e_bound+2*abs(fraction_of_mpf(Lb))*eb+eu
        # Inflate once by16 context eps to preserve the upper-bound direction
        # through integer-to-MP conversion, division, multiplication and summing.
        bound_mp=ctx.mpf(rp_bound.numerator)/rp_bound.denominator
        bounds[b]=(bounds[b]+bound_mp)*(1+16*ctx.eps)
        for j in range(length):
            bins[b][0][j] += zr
            bins[b][1][j] += round_shift(-omega_i*zi,p)
            bins[b][2][j] += round_shift(kappa_i*zr,p)
            if j+1 != length:
                nr = round_shift(zr*er-zi*ei,p)
                ni = round_shift(zr*ei+zi*er,p)
                zr,zi = nr,ni
        E = ctx.exp(ctx.j*omega*(times[-1]-times[0]))
        du = ub-ua-ctx.expm1(ctx.j*omega*(times[-1]-times[0]))*d
        dw = wb-E*wa
        flow = mu*((2*k*k+3*Lb*Lb)*du.real/k-dw.imag-Lb*dw.real/k)/(2*eps)
        triangle = mu*((2*k*k+3*Lb*Lb)*abs(du)/k+(1+Lb/k)*abs(dw))/(2*eps)
        flows[b][0] += flow
        flows[b][1] += triangle
        # Direct per-mode primitive independently phase-evaluated, compared to recurrence aggregate.
        primitive = mu*(3*c.real*(Lb*Lb-La*La)
            +(d*(3*(Lb*Lb*E-La*La)-ctx.j*omega*(Lb*E-La))).real)/(2*k*eps)
        flows[b][2] += primitive
        endpoint_checks[b] += abs(dyadic(ctx,zr,p)- (alpha*E).real)
    rows = []
    cq=ce=flow=triangle=direct_i=endpoint_gap=ctx.zero
    sums = [[0]*length for _ in range(3)]
    cumulative_bound = ctx.zero
    for b, cutoff in enumerate(CUTOFFS):
        cq += constants[b][0]
        ce += constants[b][1]
        flow += flows[b][0]
        triangle += flows[b][1]
        direct_i += flows[b][2]
        endpoint_gap += endpoint_checks[b]
        cumulative_bound += bounds[b]
        for j in range(length):
            for m in range(3):
                sums[m][j] += bins[b][m][j]
        A,B,U = [[dyadic(ctx,v,p) for v in row] for row in sums]
        R=[]; P=[]; F=[]
        for j,t in enumerate(times):
            L=-1/t
            R.append(ce+3*L*L*(cq+A[j])-L*B[j])
            P.append(ce/3-L*L*(cq+A[j])-U[j]-L*B[j])
            F.append(6*L**3*(cq+A[j])+3*L*U[j]+2*L*L*B[j])
        I = 3*(Lb*Lb-La*La)*cq+3*(Lb*Lb*A[-1]-La*La*A[0])-(Lb*B[-1]-La*B[0])
        rows.append({'K':cutoff,'times':times,'R':R,'P':P,'F':F,'I_ab':direct_i,
                     'E_flow':flow,'triangle_bound':triangle,
                     'independent_phase_primitive':direct_i,
                     'recurrence_phase_primitive':I,
                     'recurrence_primitive_gap':I-direct_i,
                     'recurrence_endpoint_moment_gap':endpoint_gap,
                     'quantization_only_RP_bound':cumulative_bound})
    return ctx, rows, indices, dt_ld


def diagnostics(a, setting, ctx, rows, indices, dt):
    ia,ib=int(indices[0]),int(indices[-1])
    hist=a['history_ledger_integrand']
    S=simpson_to(hist,ib,dt)-simpson_to(hist,ia,dt) # subtraction remains binary80.
    reset=simpson_to(hist[ia:ib+1],ib-ia,dt)
    doubled=None
    if setting == 'fine':
        doubled=simpson_to(hist[::2],ib//2,2*dt)-simpson_to(hist[::2],ia//2,2*dt)
    rr,pp=QUANTITIES.index('rho'),QUANTITIES.index('p')
    for index,row in enumerate(rows):
        stored_R=[represented(ctx,v) for v in a['history_values'][indices,index,rr]]
        stored_P=[represented(ctx,v) for v in a['history_values'][indices,index,pp]]
        stored_F=[represented(ctx,v) for v in hist[indices,index]]
        deltaR=represented(ctx,a['history_values'][ib,index,rr])-represented(ctx,a['history_values'][ia,index,rr])
        row.update({'stored_R':stored_R,'stored_P':stored_P,'stored_F':stored_F,
                    'S_ab':represented(ctx,S[index]),'S_direct_reset':represented(ctx,reset[index]),
                    'DeltaR':deltaR})
        row['S_global_minus_reset']=row['S_ab']-row['S_direct_reset']
        row['D_S']=deltaR-row['S_ab']
        row['D_cont']=deltaR-row['I_ab']
        row['E_Q']=row['I_ab']-row['S_ab']
        row['decomposition_error']=row['D_S']-row['D_cont']-row['E_Q']
        row['flow_operator_closure']=row['D_cont']-row['E_flow']
        row['R_profile_max_error']=max(abs(x-y) for x,y in zip(row['R'],stored_R))
        row['P_profile_max_error']=max(abs(x-y) for x,y in zip(row['P'],stored_P))
        row['F_profile_max_error']=max(abs(x-y) for x,y in zip(row['F'],stored_F))
        if doubled is not None:
            row['S_ab_double_global']=represented(ctx,doubled[index])
            row['S_native_minus_double']=row['S_ab']-row['S_ab_double_global']
    return rows


def serialize(ctx, rows):
    """All nonmetadata scientific MP scalar/profile entries use exact-decimal checks."""
    maximum=Fraction(0)
    def conv(x):
        nonlocal maximum
        if isinstance(x,list):
            return [conv(v) for v in x]
        require(ctx.isfinite(x),'Nonfinite reported scientific value')
        s=ctx.nstr(x,n=ctx.dps+12,strip_zeros=False)
        error=abs(fraction_of_mpf(x)-Fraction(Decimal(s)))
        maximum=max(maximum,error)
        return s
    out=[]
    for row in rows:
        out.append({key:(value if key=='K' else conv(value)) for key,value in row.items()})
    require(maximum <= Fraction(GATES['serialization']),'Serialization precision gate failed')
    return out,str(Decimal(maximum.numerator)/Decimal(maximum.denominator))


def gap_and_gates(level80, level100):
    require(len(level80)==len(level100)==3 and tuple(r['K'] for r in level80)==CUTOFFS
            and tuple(r['K'] for r in level100)==CUTOFFS,'Exactly three ordered K rows required at each precision')
    gap=Fraction(0)
    failures=[]; fatal=[]
    def compare(a,b):
        nonlocal gap
        if isinstance(a,list):
            require(len(a)==len(b),'Precision profile lengths differ')
            for x,y in zip(a,b): compare(x,y)
        else:
            gap=max(gap,abs(Fraction(Decimal(a))-Fraction(Decimal(b))))
    for r80,r100 in zip(level80,level100):
        require(r80.keys()==r100.keys() and r80['K']==r100['K'],'Precision output schema differs')
        for key in r80:
            if key!='K': compare(r80[key],r100[key])
        for row,dps in ((r80,80),(r100,100)):
            for key in ('decomposition_error','recurrence_primitive_gap'):
                if abs(Fraction(Decimal(row[key])))>Fraction(GATES['decomposition']):
                    fatal.append({'K':row['K'],'decimal_dps':dps,'gate':key,
                                  'value':row[key],'threshold':GATES['decomposition']})
        for key,gate in (('R_profile_max_error','profile'),('P_profile_max_error','profile'),
                         ('D_cont','D_cont'),('E_flow','E_flow')):
            if abs(Fraction(Decimal(r100[key])))>Fraction(GATES[gate]):
                item={'K':r100['K'],'gate':key,'value':r100[key],'threshold':GATES[gate]}
                failures.append(item)
    if gap>Fraction(GATES['precision_gap']): fatal.append({'gate':'precision_gap','exact_gap':str(gap)})
    return gap,failures,fatal


def run_arrays(inputs, started, progress=None):
    records=[]; failures=[]; fatal_failures=[]; max_gap=Fraction(0)
    for source in SOURCES:
        for setting in SETTINGS:
            a=inputs(source,setting)
            zero_signs={}
            for key,value in a.items():
                if value.dtype.kind == 'f':
                    zero_signs[key]=int(np.count_nonzero((value == 0) & np.signbit(value)))
                elif value.dtype.kind == 'c':
                    zero_signs[key]={'real':int(np.count_nonzero((value.real == 0) & np.signbit(value.real))),
                                     'imag':int(np.count_nonzero((value.imag == 0) & np.signbit(value.imag)))}
            serialized=[]; ser_errors=[]
            for dps,p in LEVELS:
                ctx,rows,indices,dt=moment_kernel(a,setting,dps,p,started)
                rows=diagnostics(a,setting,ctx,rows,indices,dt)
                encoded,error=serialize(ctx,rows)
                serialized.append(encoded); ser_errors.append(error)
                if progress is not None:
                    progress({'source':source,'setting':setting,'decimal_dps':dps,'event':'precision_complete',
                              'elapsed_seconds':time.monotonic()-started})
            gap,bad,fatal=gap_and_gates(*serialized)
            max_gap=max(max_gap,gap)
            failures.extend({'source':source,'setting':setting,**f} for f in bad)
            fatal_failures.extend({'source':source,'setting':setting,**f} for f in fatal)
            records.append({'source':source,'setting':setting,'levels':{'80':serialized[0],'100':serialized[1]},
                            'input_negative_zero_counts':zero_signs,
                            'serialization_max_error':{'80':ser_errors[0],'100':ser_errors[1]},
                            'precision_gap_exact_rational':str(gap)})
    elapsed=time.monotonic()-started
    rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if elapsed>900 or rss>262144: fatal_failures.append({'gate':'shared_resource_budget','seconds':elapsed,'rss_kib':rss})
    attribution=[]
    if not failures and not fatal_failures:
        for record in records:
            for row in record['levels']['100']:
                ds=abs(Fraction(Decimal(row['D_S'])))
                if ds>Fraction(GATES['attribution']) and abs(Fraction(Decimal(row['D_cont'])))<=ds/10 \
                   and abs(Fraction(Decimal(row['E_Q'])))>=9*ds/10:
                    attribution.append({'source':record['source'],'setting':record['setting'],'K':row['K']})
    return {'schema_version':1,'route':'primary_integer_three_moment',
            'status':'COMPLETED_SAVED_DATA_DIAGNOSTIC' if not fatal_failures else 'FATAL_DIAGNOSTIC_FAILURE',
            'classification':('CONSISTENCY_FAILURE' if failures or fatal_failures else
                              ('LEDGER_ERROR_DEMONSTRATED' if attribution else 'NO_GATE_SCALE_ATTRIBUTION')),
            'old_metric_status':'FAIL',
            'prior_scientific_calibration':'FAIL_UNCHANGED','fixed_cases':12,
            'attribution':'POST_SUPPORT_LEDGER_ERROR_DEMONSTRATED' if attribution else 'NOT_DEMONSTRATED',
            'attribution_cases':attribution,'failures':failures,'fatal_failures':fatal_failures,'records':records,
            'precision_gap_max_exact_rational':str(max_gap),
            'resources':{'seconds':elapsed,'peak_rss_kib':rss,'scope':'all twelve cases and both decimal precision passes',
                         'limit_seconds':900,'limit_peak_rss_kib':262144},
            'arithmetic':{'levels_decimal_dps':[80,100],'binary_fractional_scales':[266,333],
                          'guard_digits':0,'rounding':'nearest ties-even','phase_refresh':'none',
                          'phase_work_decimal_dps':[80,100],
                          'quantization_bound_scope':'integer projections/recurrence only, conditional on work-context transcendental seeds; not total arithmetic or physical error',
                          'gap_universe':'every non-K scalar/profile value under records.levels'},
            'gates':GATES}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--inputs',type=Path)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--checkpoint-root',type=Path)
    parser.add_argument('--registration-sha256')
    parser.add_argument('--freeze-commit')
    parser.add_argument('--output-dir',type=Path)
    parser.add_argument('--preflight-only',action='store_true')
    args=parser.parse_args()
    started=time.monotonic()
    auth=None
    if args.checkpoint_root is not None:
        require(args.inputs is None and args.output is None,'Registered CLI has fixed checkpoint input/output paths')
        require(args.output_dir is not None and not args.output_dir.exists(),'Fresh output directory required')
        require(sys.version_info[:3]==(3,12,14) and np.__version__=='2.2.6' and mpmath.__version__=='1.3.0',
                'Pinned Python3.12.14/NumPy2.2.6/mpmath1.3.0 runtime required')
        verify_frozen=authenticated_helper(args.checkpoint_root,args.registration_sha256,args.freeze_commit)
        auth=verify_frozen(args.checkpoint_root,args.registration_sha256,args.freeze_commit)
        args.inputs=args.checkpoint_root/'inputs'
        args.output=args.output_dir/'diagnostic.json'
        args.output_dir.mkdir(parents=True)
    else:
        require(args.preflight_only and args.inputs is not None and args.output is not None,
                'Physical saved-data execution requires frozen checkpoint authentication')
    require(not args.output.exists(),'Output must be new; old evidence cannot be overwritten')
    manifest,caps=input_metadata(args.inputs)
    if args.preflight_only:
        result={'status':'PASS_METADATA_ONLY','arrays_loaded':False,'input_manifest_sha256':sha(args.inputs/'INPUT_MANIFEST.json'),
                'capsule_count':len(caps),'main_execution_budget_scope':'12 cases including both80/100 decimal passes'}
    else:
        def load(source,setting):
            with np.load(caps[(source,setting)],allow_pickle=False) as z:
                return {key:z[key] for key in MEMBERS}
        def progress(record):
            with (args.output_dir/'progress.jsonl').open('a') as f:
                f.write(json.dumps(record,sort_keys=True)+'\n')
        result=run_arrays(load,started,progress)
        result['input_manifest_sha256']=sha(args.inputs/'INPUT_MANIFEST.json')
        result['input_lineage_public_commit']=manifest['published_lineage_commit']
        result['registration_sha256']=args.registration_sha256
        result['freeze_commit']=args.freeze_commit
        verify_frozen(args.checkpoint_root,args.registration_sha256,args.freeze_commit)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':result['status'],'output':str(args.output),'sha256':sha(args.output)}))
    return 0 if not result.get('fatal_failures') else 2


if __name__=='__main__':
    sys.exit(main())
