#!/usr/bin/env python3
"""Exact read-only verification of the manuscript's reported numerical claims."""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,json
ROOT=Path(__file__).resolve().parents[1]

def need(ok,label):
    if not ok:raise ValueError(label)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(name):return json.loads((ROOT/name).read_text())

def check(pdf):
    table=read('manuscript/TABLE_VERIFICATION_RECORD.json')
    source=read('numeric/ACTUAL_SOURCES_NORMAL.json')
    norm=read('numeric-review/SOURCE_REPORT_REVIEW_NORMAL.json')
    response=read('response-signs/EVOLVED_SCALAR_SIGNS_NORMAL.json')
    need(table['source_report_sha256']==sha(ROOT/'numeric/ACTUAL_SOURCES_NORMAL.json'),'TABLE_SOURCE_HASH')
    need(table['table_sha256']==sha(ROOT/'manuscript/NUMERIC_TABLE.tex'),'TABLE_TEX_HASH')
    actual={'exterior_QE_yy_mean':source['bulk_witness']['sources']['QE']['yy'][0],
        'reflected_Qphi_mean':source['reflected_witness']['sources']['Qphi'][0],
        'shell_Qp_mean':source['shell_witness']['sources']['one_copy_outward_hessian']['scalar'][0]}
    for entry in table['entries']:
        need(entry['exact']==actual[entry['name']],'TABLE_EXACT_ENTRY_'+entry['name'])
        lo,hi=map(F,entry['exact']);dlo,dhi=map(F,entry['outward_display'])
        need(dlo<=lo<=hi<=dhi,'TABLE_OUTWARD_DISPLAY_'+entry['name'])
    need(source['global_cover']['cells']==81 and source['evolved_radial_profile']['cells']==64
         and source['shell_whole_time_cover']['cells']==6,'COMPLETE_CELL_COUNTS')
    need(all(F(row['sources'][name]['yy'][0][1])<0
             for row in source['evolved_radial_profile']['cell_enclosures']
             for name in ('QE','QEred')),'ALL_64_PROFILE_NATURAL_REDUCED_MEANS_NEGATIVE')
    fh=F(table['bulk_F_H_upper']);oldfh=F(table['parent_bulk_F_H_upper'])
    need(fh==F(norm['independent_positive_forcing_H_upper'])==17387513156573,'FULL_FORCING_H_BOUND')
    need(oldfh==171610068644000 and oldfh/fh>F(98,10),'UPPER_BOUND_IMPROVEMENT_RATIO')
    need(response['mean_projection']['functional_norm_exact']=='1','MEAN_FUNCTIONAL_NORM_ONE')
    bound=F(response['response_join']['scalar_mean_x_derivative_strict_upper'])
    need(bound==F('1.03e33'),'RETAINED_FIRST_DERIVATIVE_BOUND')
    drift=F(1,200)*bound; margin=7*10**30-drift
    need(drift==F('5.15e30') and margin==F('1.85e30'),'DRIFT_AND_PERSISTENCE_MARGIN')
    targetcells=[(F(73,64),F(585,512)),(F(951,512),F(119,64))]
    for cell,target in zip(response['persistence']['cells'],targetcells):
        need(tuple(map(F,cell['scaled_r_closed_interval']))==target,'AUTHENTICATED_RADIAL_CELL')
        lo,hi=map(F,cell['initial_scalar_mean_interval'])
        need(lo>7*10**30 if cell['sign']=='positive' else hi< -7*10**30,'INITIAL_STRICT_MEAN_MARGIN')
        need(F(1,200)+target[1]<=8 and target[0]>=0,'FULL_RECTANGLE_CAUSAL_INCLUSION')
    need((F(1,10**6)*F(1,1000),F(1,10**6)*F(1,200))==(F(1,10**9),F(5,10**9)),'PHYSICAL_TIME_CONVERSION')
    need(response['zero_consequence']['open_scaled_radius_bracket']==['585/512','951/512']
         and not response['zero_consequence']['continuous_zero_branch_claimed'],'RADIAL_ZERO_SCOPE')
    paths=['manuscript/GEOMETRIC_BULK_SOURCES_NOTE.tex','manuscript/NUMERIC_TABLE.tex',
       'manuscript/TABLE_VERIFICATION_RECORD.json','numeric/ACTUAL_SOURCES_NORMAL.json',
       'numeric-review/SOURCE_REPORT_REVIEW_NORMAL.json','response-signs/EVOLVED_SCALAR_SIGNS_NORMAL.json',
       'projection/GEOMETRIC_COEFFICIENT_MAPS.json','review/FINAL_REVIEW_RECEIPT.json',
       'manuscript-review/check_manuscript_arithmetic.py']
    return {'status':'PASS_INDEPENDENT_MANUSCRIPT_ARITHMETIC_AND_BINDINGS',
       'files':{name:{'sha256':sha(ROOT/name),'bytes':(ROOT/name).stat().st_size} for name in paths},
       'pdf':{'name':pdf.name,'sha256':sha(pdf),'bytes':pdf.stat().st_size},
       'checks':['three exact source entries and outward table endpoints','81 causal/64 profile/6 shell complete cells',
           'all 64 profile natural and reduced yy means strictly negative','accepted full bulk forcing H bound and ratio above9.8',
           'mean projection norm1 and derivative bound1.03e33','strict persistence margin1.85e30 on both complete rectangles',
           'exact physical time conversion and per-time radial IVT scope'],
       'scope':'Arithmetic and byte binding only; formula/proof/scope reading and PDF visual inspection are recorded separately by the manuscript reviewer.'}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--pdf',type=Path,default=ROOT/'manuscript/GEOMETRIC_BULK_SOURCES_NOTE.pdf');args=ap.parse_args()
    need(not args.output.exists(),'OUTPUT_ALREADY_EXISTS')
    out=check(args.pdf);args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(out['status'])

if __name__=='__main__':main()
