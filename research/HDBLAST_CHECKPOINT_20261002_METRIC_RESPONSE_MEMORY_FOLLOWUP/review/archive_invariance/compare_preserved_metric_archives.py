#!/usr/bin/env python3
"""Post-run archive comparison, loading one matching member pair at a time.

No producer, source, mode equation, quadrature or scientific validator is
executed. Exact ndarray value equality is separate from storage-byte equality.
"""
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import platform
import resource
import numpy as np

def require(ok,msg):
    if not ok:raise RuntimeError(msg)
def read(path):return json.loads(path.read_text())
def sha(path):
    digest=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):digest.update(block)
    return digest.hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prior',type=Path,required=True)
    parser.add_argument('--current',type=Path,required=True)
    parser.add_argument('--completion-receipt',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    require(not args.output.exists(),'Refusing overwrite')
    control=read(args.completion_receipt)
    require(control['status']=='PASS_EXPECTED_SCIENTIFIC_FAILURE_CONTROL' and control['control_status']=='EXPECTED_SCIENTIFIC_FAILURE_CONFIRMED' and control['underlying_scientific_status']=='FAIL' and control['registered_metric_calibration_passed'] is False,'Third control must finish before archived-data comparison')
    require(control['source_verification_before']==control['source_verification_after']=='PASS','Completed source authentication required')
    prior=read(args.prior/'results.json');current=read(args.current/'results.json')
    require(prior['status']==current['status']=='failed','Both historical producer statuses must remain failed')
    require(len(prior['runs'])==2 and len(current['runs'])==4,'Expected preserved partial and complete-run coverage')
    require([(r['source'],r['setting']) for r in current['runs']]==[('positive_B','coarse'),('positive_B','fine'),('signed_uB','coarse'),('signed_uB','fine')],'Third producer did not complete four ordered runs')
    new_lookup={(r['source'],r['setting']):r for r in current['runs']}
    records=[];total_arrays=total_elements=numeric_arrays=text_arrays=0
    unequal_arrays=storage_unequal_arrays=nonfinite_arrays=0
    max_member_pair_bytes=0
    for oldrun in prior['runs']:
        identity=(oldrun['source'],oldrun['setting']);newrun=new_lookup[identity]
        before=args.prior/oldrun['archive']['path'];after=args.current/newrun['archive']['path']
        old_hash,new_hash=sha(before),sha(after)
        require(old_hash==oldrun['archive']['sha256'] and new_hash==newrun['archive']['sha256'],'Run archive digest differs from saved producer report')
        fields=[]
        with np.load(before,allow_pickle=False) as old,np.load(after,allow_pickle=False) as new:
            missing=sorted(set(old.files)-set(new.files));extra=sorted(set(new.files)-set(old.files))
            member_names_equal=not missing and not extra
            for name in sorted(set(old.files)&set(new.files)):
                # Only this pair of arrays is materialized at a time.
                a=old[name];b=new[name]
                same_dtype=a.dtype==b.dtype;same_shape=a.shape==b.shape
                same_values=same_dtype and same_shape and bool(np.array_equal(a,b))
                numeric=a.dtype.kind in 'biufc';text=a.dtype.kind in 'SU'
                finite_a=bool(np.all(np.isfinite(a))) if a.dtype.kind in 'fc' else True
                finite_b=bool(np.all(np.isfinite(b))) if b.dtype.kind in 'fc' else True
                # This diagnostic can differ for longdouble padding even when
                # numerical values, shapes and dtypes are exactly equal.
                same_storage=bool(a.tobytes(order='A')==b.tobytes(order='A')) if same_dtype and same_shape else False
                info={'name':name,'prior_dtype':a.dtype.str,'current_dtype':b.dtype.str,'prior_shape':list(a.shape),'current_shape':list(b.shape),'dtype_equal':same_dtype,'shape_equal':same_shape,'values_exactly_equal':same_values,'array_storage_bytes_equal':same_storage,'numeric':numeric,'text':text,'prior_finite':finite_a,'current_finite':finite_b,'prior_elements':int(a.size),'current_elements':int(b.size)}
                if not same_values and same_dtype and same_shape:
                    unequal=a!=b;indices=np.argwhere(unequal)
                    info['unequal_element_count']=int(np.count_nonzero(unequal))
                    info['first_16_unequal_indices']=indices[:16].tolist()
                    if numeric:
                        with np.errstate(all='ignore'):
                            difference=np.abs(a-b)
                            info['maximum_absolute_difference_repr']=str(np.max(difference))
                        del difference
                    del unequal,indices
                fields.append(info)
                total_arrays+=1;total_elements+=int(a.size)
                numeric_arrays+=int(numeric);text_arrays+=int(text)
                unequal_arrays+=int(not same_values);storage_unequal_arrays+=int(not same_storage)
                nonfinite_arrays+=int(not finite_a or not finite_b)
                max_member_pair_bytes=max(max_member_pair_bytes,a.nbytes+b.nbytes)
                del a,b
            records.append({'source':identity[0],'setting':identity[1],'prior_file':str(before),'current_file':str(after),'prior_sha256':old_hash,'current_sha256':new_hash,'raw_zip_sha256_equal':old_hash==new_hash,'prior_bytes':before.stat().st_size,'current_bytes':after.stat().st_size,'prior_member_count':len(old.files),'current_member_count':len(new.files),'member_names_equal':member_names_equal,'missing_members':missing,'extra_members':extra,'fields':fields})
    membership_ok=all(r['member_names_equal'] and r['prior_member_count']==r['current_member_count']==631 for r in records)
    status='PASS_EXACT_SAVED_ARRAY_INVARIANCE' if membership_ok and unequal_arrays==0 else 'FAIL_SAVED_ARRAY_INVARIANCE'
    report={'status':status,'scope':'Read-only comparison of historical positive_B coarse/fine raw arrays after third-run completion. No new physical evaluation.','new_physical_evaluations':0,'mode_or_source_functions_executed':False,'prior_scientific_status':'FAIL','current_scientific_status':'FAIL','metric_calibration_passed':False,'archive_pair_count':len(records),'compared_member_pairs':total_arrays,'numeric_member_pairs':numeric_arrays,'text_member_pairs':text_arrays,'compared_elements':total_elements,'unequal_array_pairs':unequal_arrays,'array_storage_byte_inequality_pairs':storage_unequal_arrays,'nonfinite_array_pairs':nonfinite_arrays,'maximum_loaded_member_pair_bytes':max_member_pair_bytes,'comparison_process_peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'python':platform.python_version(),'numpy':np.__version__,'recorded_utc':datetime.now(timezone.utc).isoformat(),'completion_receipt_sha256':sha(args.completion_receipt),'prior_producer_report_sha256':sha(args.prior/'results.json'),'current_producer_report_sha256':sha(args.current/'results.json'),'comparison_source_sha256':sha(Path(__file__)),'storage_note':'NPZ compression/container metadata and unused longdouble/complex-longdouble padding bytes are not numerical equality gates. This report records raw hashes/storage differences separately and requires exact ndarray dtype/shape/value equality.','archives':records}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:report[k] for k in ('status','compared_member_pairs','numeric_member_pairs','text_member_pairs','compared_elements','unequal_array_pairs','array_storage_byte_inequality_pairs','comparison_process_peak_rss_kib')}))
    require(status=='PASS_EXACT_SAVED_ARRAY_INVARIANCE','All mismatches preserved; saved archive invariance failed')

if __name__=='__main__':main()
