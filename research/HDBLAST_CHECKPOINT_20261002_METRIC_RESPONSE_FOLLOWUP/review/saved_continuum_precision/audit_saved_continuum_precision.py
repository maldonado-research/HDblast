#!/usr/bin/env python3
"""Optional post-run saved-data audit, separate from the frozen experiment.

Loads saved decimals and floats only. Never imports the physical producer,
source evaluator, mode code, or quadrature routine. No physical evaluations.
"""
from __future__ import annotations
import argparse,copy,hashlib,json,math,platform,re,sys,traceback
from pathlib import Path
import mpmath as mp
from mpmath.libmp import dps_to_prec

CONTRACT={'method':'mpmath_tanh_sinh_endpoint_subtraction_squared','decimal_precisions':[50,70],'maxdegree':10,'panels':['0','0.5','1'],'reported_precision':70,'error_policy':'max_precision_gap_quad_estimates_and_float_conversion','source_arithmetic':'analytic_integer_polynomials_mpmath_only'}
DECIMAL=re.compile(r'^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?$',re.I)

def require(ok,message):
    if not ok:raise RuntimeError(message)
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(path,data):Path(path).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')

def decimal_box(text,dps):
    require(isinstance(text,str) and DECIMAL.fullmatch(text) is not None,'Malformed saved decimal')
    value=mp.mpf(text);require(mp.isfinite(value),'Nonfinite saved decimal')
    if value==0:return value,mp.mpf(0)
    # mp.nstr(x,dps) rounds to at most dps significant decimal digits.
    exponent=int(mp.floor(mp.log10(abs(value))))
    radius=mp.mpf('0.5000000000000001')*mp.power(10,exponent-dps+1)
    return value,radius

def absolute_difference_box(a,b):
    center=abs(a[0]-b[0]);radius=a[1]+b[1]
    return max(mp.mpf(0),center-radius),center+radius

def interval_contains_saved(bounds,saved,label):
    lo,hi=bounds;value,radius=saved
    require(lo<=value+radius and value-radius<=hi,label+' outside saved-decimal rounding range')

def exact_float_from_box(saved,value,label):
    center,radius=saved
    require(type(value) is float and math.isfinite(value),label+' is not a finite float')
    require(float(center-radius)==value==float(center+radius),label+' does not match unique rounding of saved decimal')


def audit_rows(data,experiment):
    require(data['continuum_algorithm']==CONTRACT==experiment['primary']['continuum'],'Algorithm configuration mismatch')
    expected=[(source,eta) for source in experiment['source']['ids'] for eta in experiment['observation_eta']]
    require([(r['source'],r['eta']) for r in data['rows']]==expected,'Source/observation coverage mismatch')
    require(len(expected)==12,'Expected twelve original primary rows')
    counts={key:0 for key in ('precision_runs','qjet_formula_components','quadrature_normalizations','precision_gap_ranges','float_conversion_ranges','scalar_memory_float_bindings','raw_memory_float_bindings','endpoint_jet_comparisons','endpoint_precision_pair_comparisons','saved_forcing_reconstructions','allowance_float_bindings','upward_allowance_checks','propagated_error_fields','metric_variance_sums','metric_stress_sums','current_equalities','trace_arithmetic')}
    records=[]
    with mp.workdps(110):
      for row in data['rows']:
        eta=mp.mpf(str(row['eta']));L=-1/eta;lf=-1/row['eta'];pref=1/(8*mp.pi**2);den=8*math.pi**2
        block=row['continuum'];diag=block['continuum_precision_diagnostics']
        for key in ('method','decimal_precisions','maxdegree','panels','reported_precision','error_policy'):
            require(diag[key]==CONTRACT[key],'Row algorithm declaration differs: '+key)
        require('complete qjet' in diag['saved_memory_error_scope'] and 'true raw mp.quad estimates' in diag['saved_memory_error_scope'],'Effective allowance scope is missing')
        require(block['log_history_error_scope']==diag['saved_memory_error_scope'],'Saved memory allowance scopes differ')
        require('not interval error certificates' in diag['uncertainty'],'Empirical uncertainty qualification missing')
        require([r['dps'] for r in diag['precision_runs']]==[50,70],'Precision run order/coverage mismatch')
        all_runs=[];formula_ratios=[]
        for run in diag['precision_runs']:
            dps=run['dps'];counts['precision_runs']+=1
            require(all(len(run[name])==n for name,n in [('qjet_decimal',3),('log_history_decimal',3),('quadrature_estimates_decimal',3),('raw_log_history_quad_errors_decimal',3),('endpoint_forcing_jet_decimal',4)]),'Saved precision vector shape mismatch')
            q=[decimal_box(x,dps) for x in run['qjet_decimal']]
            F=[decimal_box(x,dps) for x in run['log_history_decimal']]
            g=[decimal_box(x,dps) for x in run['endpoint_forcing_jet_decimal']]
            e=[decimal_box(x,dps) for x in run['quadrature_estimates_decimal']]
            raw=[decimal_box(x,dps) for x in run['raw_log_history_quad_errors_decimal']]
            # Reconstruct only algebra on saved values, without any source call.
            inside=[F[0][0],F[1][0]+L*g[0][0],F[2][0]+2*L*g[1][0]+L*L*g[0][0]]
            input_radii=[F[0][1],F[1][1]+abs(L)*g[0][1],F[2][1]+2*abs(L)*g[1][1]+L*L*g[0][1]]
            absolute_terms=[abs(F[0][0]),abs(F[1][0])+abs(L*g[0][0]),abs(F[2][0])+2*abs(L*g[1][0])+abs(L*L*g[0][0])]
            for n in range(3):
                # This is a saved-formula consistency tolerance for decimal and
                # context rounding, not a physical quadrature-error enclosure.
                arithmetic=64*mp.power(2,-dps_to_prec(dps))*abs(pref)*absolute_terms[n]
                tolerance=q[n][1]+abs(pref)*input_radii[n]+arithmetic
                difference=abs(q[n][0]+pref*inside[n])
                require(difference<=tolerance,'Saved qjet/contact formula mismatch')
                formula_ratios.append(float(difference/tolerance) if tolerance else 0.0)
                counts['qjet_formula_components']+=1
                require(raw[n][0]>=0 and e[n][0]>=0,'Negative saved quadrature estimate')
                tolerance_e=e[n][1]+abs(pref)*raw[n][1]+16*mp.power(2,-dps_to_prec(dps))*abs(pref*raw[n][0])
                require(abs(e[n][0]-pref*raw[n][0])<=tolerance_e,'Raw quadrature estimate normalization mismatch')
                counts['quadrature_normalizations']+=1
            all_runs.append({'q':q,'F':F,'g':g,'e':e})
        low,high=all_runs
        gaps=[decimal_box(x,70) for x in diag['qjet_precision_gaps_decimal']]
        conversions=[decimal_box(x,70) for x in diag['qjet_float_conversion_errors_decimal']]
        allowances=[decimal_box(x,70) for x in diag['effective_raw_memory_allowances_decimal']]
        require(len(gaps)==len(conversions)==len(allowances)==3,'Saved diagnostic vector shape mismatch')
        scalar_names=('scalar_memory_q','scalar_memory_q_prime','scalar_memory_q_second')
        memory=block['log_history_integrals_over_epsilon']
        require([m['forcing_derivative_order'] for m in memory]==[1,2,3],'Saved raw memory derivative order mismatch')
        qerrors=[];row_gap_bounds=[]
        for n in range(3):
            require(gaps[n][0]>=0 and conversions[n][0]>=0 and allowances[n][0]>=0,'Negative saved precision diagnostic')
            gap_box=absolute_difference_box(high['q'][n],low['q'][n])
            interval_contains_saved(gap_box,gaps[n],'Declared 50/70 precision gap')
            row_gap_bounds.append([mp.nstr(x,30) for x in gap_box]);counts['precision_gap_ranges']+=1
            actual=block['components'][scalar_names[n]]
            exact_float_from_box(high['q'][n],actual,'Scalar memory qjet float')
            counts['scalar_memory_float_bindings']+=1
            interval_contains_saved(absolute_difference_box((mp.mpf(actual),mp.mpf(0)),high['q'][n]),conversions[n],'Float conversion diagnostic')
            counts['float_conversion_ranges']+=1
            exact_float_from_box(high['F'][n],memory[n]['value'],'Raw memory float')
            counts['raw_memory_float_bindings']+=1
            exact_float_from_box(allowances[n],memory[n]['quadrature_error_estimate'],'Effective raw allowance float')
            counts['allowance_float_bindings']+=1
            normalized=memory[n]['quadrature_error_estimate']/den;qerrors.append(normalized)
            # Every possible saved decimal rounding of these diagnostics is
            # bounded by the serialized allowance or overlaps its boundary.
            upper=max(gaps[n][0]+gaps[n][1],conversions[n][0]+conversions[n][1],low['e'][n][0]+low['e'][n][1],high['e'][n][0]+high['e'][n][1])
            require(mp.mpf(normalized)>=upper,'Serialized allowance rounded downward relative to saved diagnostics')
            counts['upward_allowance_checks']+=1
            require(block['quadrature_error_estimates'][('q','q_prime','q_second')[n]]==normalized,'Saved error does not equal legacy exact division')
        endpoint_differences=[];endpoint_precision_differences=[];endpoint_precision_units=[]
        hjet=[row['source_jet_over_epsilon']['h'+str(n)] for n in range(6)]
        require(all(type(x) is float and math.isfinite(x) for x in hjet),'Malformed saved h jets')
        for n in range(4):
            hi=high['g'][n][0];lo=low['g'][n][0]
            double=row['canonical_forcing_jet_over_epsilon']['g'+str(n)]
            require(type(double) is float and math.isfinite(double),'Nonfinite saved endpoint jet')
            endpoint_differences.append(float(abs(hi-mp.mpf(double))));counts['endpoint_jet_comparisons']+=1
            delta=abs(hi-lo);endpoint_precision_differences.append(float(delta))
            unit=mp.power(2,-dps_to_prec(50))*max(1,abs(hi),abs(lo))
            endpoint_precision_units.append(float(delta/unit))
            counts['endpoint_precision_pair_comparisons']+=1
            # Pure recurrence on already saved h jets, without source sampling.
            reconstructed=4*sum(math.comb(n,j)*math.factorial(n-j+1)*lf**(n-j+2)*hjet[j] for j in range(n+1))-2*sum(math.comb(n,j)*math.factorial(n-j)*lf**(n-j+1)*hjet[j+1] for j in range(n+1))-hjet[n+2]
            require(double==reconstructed,'Saved double forcing recurrence mismatch')
            counts['saved_forcing_reconstructions']+=1
        e0,e1,e2=qerrors
        error_map={'q':e0,'q_prime':e1,'q_second':e2,'rho':(3*lf*lf*e0+lf*e1)/2,'p':(e2+3*lf*e1+3*lf*lf*e0)/6,'Q0':0.0,'rho0':0.0,'p0':0.0,'current':e0}
        require(block['quadrature_error_estimates']==error_map,'Saved complete propagated error map mismatch')
        counts['propagated_error_fields']+=len(error_map)
        require(block['density_derivative']['quadrature_error_estimate']==3*lf**3*e0+lf*lf*e1+lf*e2/2,'Saved density derivative error propagation mismatch')
        counts['propagated_error_fields']+=1
        for n,(name,local) in enumerate([('q','q_local'),('q_prime','q_local_prime'),('q_second','q_local_second')]):
            require(block['values'][name]==block['components'][scalar_names[n]]+block['local_coefficients'][local],'Saved metric variance sum mismatch')
            counts['metric_variance_sums']+=1
        for name,scalar,local in [('rho','scalar_mass_density','full_metric_local_density'),('p','scalar_mass_pressure','full_metric_local_pressure')]:
            require(block['values'][name]==block['components'][scalar]+block['components'][local],'Saved metric stress sum mismatch')
            counts['metric_stress_sums']+=1
        require(block['values']['current']==block['values']['q'],'Fixed-phi current does not equal saved q');counts['current_equalities']+=1
        require(block['trace_from_direct_stresses']==-block['values']['rho']+3*block['values']['p'],'Saved direct trace arithmetic mismatch');counts['trace_arithmetic']+=1
        if row['eta']<=-5 or row['eta']>=-3:
            require(all(value[0]==0 for run in all_runs for value in run['g']),'Source-free saved endpoint jets are not zero')
        records.append({'source':row['source'],'eta':row['eta'],'max_saved_formula_rounding_ratio':max(formula_ratios),'precision_gap_rounding_ranges':row_gap_bounds,'max_declared_precision_gap':max(float(x[0]) for x in gaps),'max_qjet_float_conversion_error':max(float(x[0]) for x in conversions),'max_serialized_normalized_allowance':max(qerrors),'max_saved_endpoint_double_difference':max(endpoint_differences),'max_saved_50_70_endpoint_decimal_difference':max(endpoint_precision_differences),'max_saved_endpoint_precision_difference_in_50_dps_units':max(endpoint_precision_units),'endpoint_double_comparison_scope':'Descriptive saved-data comparison only; no source is evaluated and no new source acceptance gate is imposed.'})
    return counts,records


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--checkpoint',type=Path,required=True);parser.add_argument('--results',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    checkpoint=args.checkpoint.resolve();results=args.results.resolve();output=args.output.resolve()
    require(not output.exists(),'Refusing to overwrite any earlier audit')
    require(not output.is_relative_to(checkpoint),'Audit output must be outside frozen checkpoint')
    output.mkdir(parents=True,exist_ok=False)
    source=Path(__file__).resolve();write(output/'STARTED.json',{'source_sha256':sha(source),'command':sys.argv,'scope':'Optional post-run saved-data audit; zero physical evaluations'})
    try:
        registration_path=checkpoint/'FULL_REGISTRATION.json';experiment_path=checkpoint/'EXPERIMENT.json'
        registration=json.loads(registration_path.read_text());experiment=json.loads(experiment_path.read_text());data=json.loads(results.read_text())
        before={name:sha(checkpoint/name) for name in registration['files']}
        require(before==registration['files'],'Frozen registration files changed before audit')
        require(registration['frozen_configuration']==experiment,'Frozen experiment configuration differs')
        provenance=data['provenance'];require(data['status']=='completed','Primary did not complete')
        require(provenance['registration_sha256']==sha(registration_path),'Primary registration SHA mismatch')
        require(provenance['experiment_sha256']==sha(experiment_path),'Primary experiment SHA mismatch')
        for key,path in [('producer_sha256','code/metric_primary.py'),('continuum_algorithm_sha256','code/continuum_metric_mp.py'),('source_jet_sha256','code/source_jet.py'),('contact_coefficients_sha256','code/metric_contact_coefficients.py')]:
            require(provenance[key]==registration['files'][path],'Primary source provenance mismatch: '+key)
        require(re.fullmatch('[0-9a-f]{40}',provenance['public_freeze_commit']) is not None,'Public freeze declaration is not a full commit')
        require(data['continuum_algorithm']==provenance['continuum_algorithm'],'Result/start algorithm declarations differ')
        started=json.loads((results.parent/'started.json').read_text());execution=json.loads((results.parent/'EXECUTION.json').read_text())
        require(started==provenance,'Saved started provenance differs from result')
        require(execution['status']=='completed' and execution['results_sha256']==sha(results),'Execution result SHA mismatch')
        counts,records=audit_rows(data,experiment)
        # Mutations affect copies of saved metadata only; no physical evaluations.
        controls=[]
        def reject(name,mutator):
            changed=copy.deepcopy(data);mutator(changed)
            try:audit_rows(changed,experiment)
            except RuntimeError as error:controls.append({'name':name,'rejected':True,'reason':str(error)})
            else:raise RuntimeError('Saved-metadata mutation was not rejected: '+name)
        reject('wrong precision declaration',lambda d:d['rows'][1]['continuum']['continuum_precision_diagnostics']['precision_runs'][0].update(dps=49))
        reject('corrupt raw memory decimal',lambda d:d['rows'][1]['continuum']['continuum_precision_diagnostics']['precision_runs'][1]['log_history_decimal'].__setitem__(0,'999'))
        reject('corrupt reported precision gap',lambda d:d['rows'][1]['continuum']['continuum_precision_diagnostics']['qjet_precision_gaps_decimal'].__setitem__(0,'1'))
        reject('corrupt float conversion diagnostic',lambda d:d['rows'][1]['continuum']['continuum_precision_diagnostics']['qjet_float_conversion_errors_decimal'].__setitem__(0,'1'))
        reject('downward raw allowance serialization',lambda d:d['rows'][1]['continuum']['log_history_integrals_over_epsilon'][0].update(quadrature_error_estimate=0.0))
        require({name:sha(checkpoint/name) for name in registration['files']}==before,'Frozen source changed during audit')
        receipt={'status':'PASS','scope':'Optional post-run saved-data consistency audit; not one of the frozen numerical experiment commands','physical_evaluations':0,'new_source_samples':0,'new_mode_evolutions':0,'new_quadratures':0,'source_sha256':sha(source),'results_sha256':sha(results),'registration_sha256':sha(registration_path),'experiment_sha256':sha(experiment_path),'public_freeze_commit_declaration':provenance['public_freeze_commit'],'frozen_source_files_verified':len(before),'frozen_source_bytes_unchanged':True,'python':platform.python_version(),'mpmath':mp.__version__,'optimized':bool(sys.flags.optimize),'counts':counts,'rows':records,'saved_metadata_mutations':controls,'independently_reconstructible':['Algebraic qjet/contact consistency from saved memories and endpoint jets within declared decimal/context rounding','Normalized raw quadrature estimate arithmetic','Reported precision gaps consistent with saved decimal rounding ranges','Final binary64 qjet and raw-memory conversion','Saved float-conversion diagnostics consistent with 70-digit qjet rounding ranges','Upward effective allowance and exact legacy error propagation','Saved metric component sums, current equality and trace arithmetic','Original input/result SHA bindings and source-byte preservation'],'source_bound_declarations_only':['Actual tanh-sinh nodes,degree reached and integration convergence are not independently replayed','Underlying B/uB endpoint values are not newly evaluated','Exact70-digit precision gap cannot generally be recovered from only50-digit low-precision qjet strings; audit verifies the allowed rounding range','Raw quadrature estimates and precision agreement are empirical diagnostics, not certified physical error bounds','Public remote commit publication timing is declared and source-bound, not re-queried by this offline audit']}
        write(output/'SAVED_CONTINUUM_PRECISION_AUDIT.json',receipt)
        print(json.dumps({'status':'PASS','rows':len(records),'precision_runs':counts['precision_runs'],'qjet_formula_components':counts['qjet_formula_components'],'saved_metadata_mutations':len(controls),'physical_evaluations':0}))
    except Exception as error:
        write(output/'FAILURE.json',{'status':'FAIL','source_sha256':sha(source),'error':str(error),'traceback':traceback.format_exc(),'physical_evaluations':0})
        raise
if __name__=='__main__':main()
