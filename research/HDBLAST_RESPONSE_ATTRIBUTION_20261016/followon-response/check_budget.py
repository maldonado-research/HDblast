#!/usr/bin/env python3
"""Small exact conditional follow-on arithmetic; no response computation."""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import hashlib
import json

ROOT=Path(__file__).resolve().parent

def need(test,message):
    if not test:raise ValueError(message)

def calculate():
    pins=json.loads((ROOT/'INPUT_PINS.json').read_text())
    for name,pin in pins.items():
        raw=(ROOT/name).read_bytes()
        need(len(raw)==pin['bytes'] and hashlib.sha256(raw).hexdigest()==pin['sha256'],'INPUT_AUTHENTICATION:'+name)
    refined=json.loads((ROOT/'inputs/REFINED_RESPONSE_NORMAL.json').read_text())
    actual=json.loads((ROOT/'inputs/ACTUAL_SOURCES_NORMAL.json').read_text())
    signs=json.loads((ROOT/'inputs/EVOLVED_SCALAR_SIGNS_NORMAL.json').read_text())
    need(signs['input_sha256']['refined']==pins['inputs/REFINED_RESPONSE_NORMAL.json']['sha256'],'MATCHED_SIGN_AND_REFINED_PREMISES')
    q=Q(refined['retained_join']['q']);B=Q(refined['response']['X_simple_upper'])
    FH=Q(actual['global_cover']['summary']['positive_bundle_bulk_forcing_H_bound']['F_H_upper'])
    GH=Q(refined['retained_join']['G_H'])
    need(q==Q(2109,10**6) and B==Q(103,100)*10**33,'REGISTERED_CONTRACTION_AND_RESPONSE')
    need(FH==17387513156573 and GH==3851400,'FULL_BULK_AND_SHELL_CAPS')
    S=3*FH+GH;tail=q*q*B;total_approx=Q(10**29)
    mean_error=3*(tail+total_approx)
    margins=[Q(c['signed_magnitude_exceeds_common']) for c in signs['persistence']['cells']]
    need(len(margins)==2 and margins[0]==margins[1]==Q(185,100)*10**30,'ACCEPTED_EVOLVED_MEAN_MARGIN')
    margin=margins[0]
    need(mean_error<margin,'NONVACUOUS_TOTAL_MEAN_TARGET')
    eb=Q(4*10**28);tau=Q(4,10**5);btilde=Q(104,100)*10**33;ealg=Q(10**28)
    need((1+q)*B+eb<btilde,'APPROXIMATE_BASE_NORM_TARGET')
    ledger=(1+q)*eb+tau*btilde+ealg
    need(ledger<total_approx,'ILLUSTRATIVE_LEDGER_WITHIN_TARGET')
    need(3*10**28+2*S<eb,'CONSERVATIVE_SOURCE_ALLOWANCE_FITS_BASE_TARGET')
    WB=S/(1-q);Wmean=3*WB;Wratio=Wmean/margin;Wtail=q*q*WB
    need(WB<Q(5228,1000)*10**13,'FORCED_INCREMENT_ROUNDED_BOUND')
    need(Wratio<Q(8477,1000)*Q(1,10**17),'FORCED_INCREMENT_TO_SIGN_MARGIN')
    need(Wtail<Q(2326,1000)*10**8,'FORCED_TWO_TERM_TAIL')
    need(Wtail<10**10 and 2*S>10**10,'FORCED_TARGET_REQUIRES_SOURCE_REFINEMENT')
    result={'status':'PASS_PROSPECTIVE_BUDGET_ARITHMETIC_ONLY',
            'registered':{'q':str(q),'response_X_upper':str(B),'F_H':str(FH),'full_G_H':str(GH)},
            'full_source_ledger':{'S_X_upper':str(S),'two_cap_source_approximation_allowance':str(2*S),
                                  'one_correction_source_error_amplification':str((1+q)*S)},
            'total_response_target':{'exact_two_term_X_tail':str(tail),'hypothetical_numerical_X_error':str(total_approx),
                                     'combined_X_error':str(tail+total_approx),'mean_value_error':str(mean_error),
                                     'existing_mean_margin':str(margin),'remaining_margin':str(margin-mean_error),
                                     'illustrative_approximation_ledger':str(ledger),
                                     'epsilon_b_target':str(eb),'operator_error_target':str(tau),
                                     'btilde_X_target':str(btilde),'final_numerical_X_target':str(ealg)},
            'forced_increment_comparison':{'W_X_upper':str(WB),'W_mean_upper':str(Wmean),
                                           'ratio_to_existing_mean_margin':str(Wratio),
                                           'two_term_X_tail':str(Wtail),
                                           'illustrative_absolute_X_target':str(10**10),
                                           'coarse_source_error_exceeds_smaller_target':True,
                                           'auxiliary_reduced_operator_split_not_separate_constraint_solutions':True},
            'scope':{'computed_candidate':False,'continuous_residual_certified':False,
                     'source_derivative_estimate_used':False,'parent_scientific_proofs_replayed':False,
                     'no_waveform_or_nonlinear_claim':True}}
    return result

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path);args=ap.parse_args()
    if args.output and args.output.exists():raise ValueError('OUTPUT_EXISTS')
    raw=json.dumps(calculate(),indent=2,sort_keys=True)+'\n'
    if args.output:args.output.write_text(raw)
    else:print(raw,end='')

if __name__=='__main__':main()
