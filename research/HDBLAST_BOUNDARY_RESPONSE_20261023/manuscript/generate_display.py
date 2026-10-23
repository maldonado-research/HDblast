#!/usr/bin/env python3
"""Authenticate exact response/source reports and generate manuscript displays.

No numerical source evaluation. All decisions and decimal rounding are exact.
Output paths and input locations do not enter the deterministic result bytes.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys

MEAN_SHA='efb890b316af82a5d2c5235770446baf160c9c5b30b495a9c54151498a1efaf5'
SOURCE_SHA='6caf16c304fb6db940899d58d09b648066eb2356385189ccb79eb2b0a6879675'

def need(ok,message):
    if not ok:raise ValueError(message)

def load(path,pin):
    data=path.read_bytes()
    need(hashlib.sha256(data).hexdigest()==pin,'exact input report authentication')
    return json.loads(data)

def format_ticks(ticks,places=6):
    sign='-' if ticks<0 else ''
    whole,frac=divmod(abs(ticks),10**places)
    return sign+str(whole)+'.'+str(frac).zfill(places)

def decimal_bounds(value,places=6):
    z=value*10**places
    return [format_ticks(z.__floor__(),places),format_ticks(z.__ceil__(),places)]

def grouped(value):
    whole,dot,frac=value.partition('.')
    return format(int(whole),',')+(dot+frac if dot else '')

def calculate(mean,source):
    need(mean['status']=='PASS_CONDITIONAL_MEAN_BLOCK_POSITIVE_POINT_CERTIFICATE','mean certificate status')
    need(source['status']=='CONDITIONAL_ACTUAL_BOUNDARY_AUXILIARY_MEAN_ENCLOSURE','source certificate status')
    need(source['partition']['n']==16,'selected source mesh')
    need(source['domain']['a']=='3/2' and source['domain']['b']=='1/4','selected apex')
    slo,shi=map(Q,source['signed_integral_ledger']['S_equals_KF_plus_LG_mean'])
    eb=Q(mean['bulk_feedback_value_bound']);er=Q(mean['Robin_feedback_value_bound'])
    e=Q(mean['total_feedback_value_bound']);lo,hi=map(Q,mean['actual_Wyy_probability_mean_interval'])
    need(e==eb+er and e>0,'complete feedback allowance')
    need(mean['evaluated_source_integral_interval']==[str(slo),str(shi)],'source join identity')
    need((lo,hi)==(slo-e,shi+e) and lo>0,'exact positive interval')
    simple_lo,simple_hi=Q(55843882),Q(737593144521)
    need(simple_lo<lo<=hi<simple_hi,'strict safe integer display')
    q=Q(mean['unweighted_max_row_contraction'])
    need(0<=q<Q(853,1000000),'visible contraction bound')
    values={'SourceLower':slo,'SourceUpper':shi,'BulkFeedback':eb,
            'RobinFeedback':er,'Feedback':e,'MeanLower':lo,'MeanUpper':hi}
    brackets={key:decimal_bounds(value) for key,value in values.items()}
    display={key:grouped(bounds[0 if key.endswith('Lower') else 1]) for key,bounds in brackets.items()}
    # Feedback fields are displayed as upper bounds; all lower endpoints use floor.
    display.update({'SimpleLower':grouped(str(simple_lo)),'SimpleUpper':grouped(str(simple_hi)),
                    'ContractionUpper':'0.000853'})
    tex='% Generated from authenticated exact fractions by generate_display.py.\n'
    for name,value in display.items():
        tex+='\\newcommand{\\'+name+'}{'+value.replace(',','{,}')+'}\n'
    plain=('A certified auxiliary response after shell contact\n'
        '23 October 2026; point (x,r)=(3/2,1/4).\n'
        'S mean interval: ['+display['SourceLower']+', '+display['SourceUpper']+'].\n'
        'Bulk feedback upper bound: '+display['BulkFeedback']+'.\n'
        'Robin feedback upper bound: '+display['RobinFeedback']+'.\n'
        'Complete mean feedback upper bound: '+display['Feedback']+'.\n'
        'Actual auxiliary W_yy probability mean: ['+display['MeanLower']+', '+display['MeanUpper']+'].\n'
        'Strict safe display: '+display['SimpleLower']+' < mean W_yy < '+display['SimpleUpper']+'.\n'
        'Decimal interval endpoints are rounded outward with exact rational decisions.\n'
        'Conditional on the retained actual first-solution/operator/source premises; one point only.\n')
    result={'status':'PASS_AUTHENTICATED_EXACT_MANUSCRIPT_DISPLAYS',
        'input_sha256':{'mean_certificate':MEAN_SHA,'source_integral':SOURCE_SHA},
        'generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'exact_values':{key:str(value) for key,value in values.items()},
        'outward_decimal_brackets':brackets,'visible_values':display,
        'strict_positive_margin_exact':str(lo),'safe_integer_inequality_verified':True,
        'single_point_only':True,'source_evaluations':0,
        'generated_tex_sha256':hashlib.sha256(tex.encode()).hexdigest(),
        'generated_text_sha256':hashlib.sha256(plain.encode()).hexdigest()}
    return result,tex,plain

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for arg in ['mean-report','source-report','json-output','tex-output','text-output']:
        p.add_argument('--'+arg,type=Path,required=True)
    a=p.parse_args()
    need(sys.flags.isolated and sys.flags.no_site and sys.flags.dont_write_bytecode,'run python -I -S -B')
    outs=[a.json_output,a.tex_output,a.text_output]
    need(len({str(x.resolve()) for x in outs})==3 and all(not x.exists() for x in outs),'new distinct output files required')
    result,tex,plain=calculate(load(a.mean_report,MEAN_SHA),load(a.source_report,SOURCE_SHA))
    for path,data in zip(outs,[json.dumps(result,indent=2)+'\n',tex,plain]):path.write_text(data)
    print(result['status'])
