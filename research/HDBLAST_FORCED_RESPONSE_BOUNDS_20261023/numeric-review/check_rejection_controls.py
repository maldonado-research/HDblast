#!/usr/bin/env python3
"""Targeted rejection checks for scientifically consequential false certificates."""
from pathlib import Path
from fractions import Fraction as Q
import argparse, copy, hashlib, importlib.util, json, tempfile

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('independent_response_audit',HERE/'check_cone_response.py')
review=importlib.util.module_from_spec(spec);spec.loader.exec_module(review)

def changed(x,case):
    x=copy.deepcopy(x)
    if case=='double_K_normalization':
        x['evaluated_Duhamel_mean_interval']=[str(2*Q(a)) for a in x['evaluated_Duhamel_mean_interval']]
    elif case=='omit_full_feedback_error':
        x['actual_W_yy_probability_mean_interval']=x['evaluated_Duhamel_mean_interval'][:]
    elif case=='wrong_coupled_coefficient_sum':
        x['full_coupling']['sigma_upper']='703/1000000'
    elif case=='omit_positive_area_cell':
        x['actual_source_cell_enclosures'].pop()
    elif case=='omit_scalar_coupled_sector':
        del x['actual_source_cell_enclosures'][0]['sources']['scalar_forcing']
    elif case=='wrong_clipped_triangle_area':
        row=x['actual_source_cell_enclosures'][0];row['exact_area']=str(2*Q(row['exact_area']))
    elif case=='wrong_rectangle_derivative_weight':
        row=x['continuous_rectangle_corollary'];row['mean_variation_upper']=str(3*Q(row['mean_variation_upper']))
    elif case=='wrong_inherited_source_authentication':
        x['source_authentication']['manifest_sha256']='0'*64
    elif case=='wrong_evaluator_pin':
        x['evaluator_sha256']='0'*64
    else:raise ValueError('UNKNOWN_CONTROL')
    return x

def run(report,source,producer):
    original=review.load(report)
    baseline=review.run(report,source,producer,False)
    cases=['double_K_normalization','omit_full_feedback_error','wrong_coupled_coefficient_sum',
           'omit_positive_area_cell','omit_scalar_coupled_sector','wrong_clipped_triangle_area',
           'wrong_rectangle_derivative_weight','wrong_inherited_source_authentication','wrong_evaluator_pin']
    results=[]
    with tempfile.TemporaryDirectory(prefix='hdblast-cone-rejections-',dir='/tmp') as folder:
        root=Path(folder)
        for case in cases:
            path=root/(case+'.json');path.write_text(json.dumps(changed(original,case),indent=2)+'\n')
            try:review.run(path,source,producer,False)
            except ValueError as error:results.append({'case':case,'rejected':True,'error':str(error)})
            else:raise ValueError('UNREJECTED_FALSE_CERTIFICATE_'+case)
    return {'status':'PASS_TARGETED_FALSE_CERTIFICATE_REJECTIONS','baseline_passed':True,
            'baseline_report_sha256':baseline['producer_report_sha256'],
            'reviewer_source_sha256':review.digest(HERE/'check_cone_response.py'),
            'control_source_sha256':review.digest(Path(__file__)),
            'rejections':results,'count':len(results),
            'scope':'Nine targeted metadata, geometry, source-sector, normalization and error-ledger corruptions; no claim that arbitrary fabricated source rows are independently disproved by this arithmetic audit.'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--report',type=Path,required=True);p.add_argument('--source-root',type=Path,required=True)
    p.add_argument('--producer-script',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    review.need(not a.output.exists(),'REFUSE_OUTPUT_OVERWRITE')
    r=run(a.report,a.source_root,a.producer_script);a.output.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))

if __name__=='__main__':main()
