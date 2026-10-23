#!/usr/bin/env python3
"""Independent exact geometry, Gram, and response-transfer audit.

The source reports are accepted actual-source enclosure premises. Optional
oracle spot replay verifies one clipped and one full cell through the frozen
all-source adapter, separate from the producer's forcing-only call path.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse, hashlib, importlib.util, json, math, sys

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('independent_cone_geometry',HERE/'cone_geometry.py')
geometry=importlib.util.module_from_spec(spec);spec.loader.exec_module(geometry)
SOURCE_MANIFEST='3ff1f95c151d2da6a1ab5e9b54718c2bb325c10d9c62a7ed9f97ef2bb36f0bb2'
SOURCE_REPORT='45a5550d2e4647bdfb84f9c4a36bf9861ef500a1de89a76f978c5afe9c51a60f'
MAP_SHA='34b10b62d0b36a4506a35407a591a9fdff319971994f25e0b5f031f48eb4f615'
EXPECTED={'metric_forcing/yy':(3,1),'metric_forcing/y0':(3,2),'metric_forcing/00':(3,1),
          'metric_forcing/ya':(2,2),'metric_forcing/0a':(2,2),'metric_forcing/ab':(7,1),'scalar_forcing':(3,1)}

def need(v,message):
    if not v:raise ValueError(message)

def no_dupes(pairs):
    out={}
    for k,v in pairs:
        need(k not in out,'DUPLICATE_JSON_KEY');out[k]=v
    return out

def load(path):return json.loads(path.read_text(),object_pairs_hook=no_dupes)
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def run(report_path,source,producer_script,oracle_spots=False):
    need(digest(source/'MANIFEST.json')==SOURCE_MANIFEST,'SOURCE_MANIFEST')
    manifest=load(source/'MANIFEST.json');listed={r['path']:r for r in manifest['files']}
    for name,row in listed.items():
        p=source/name
        need(p.stat().st_size==row['bytes'] and digest(p)==row['sha256'],'SOURCE_PAYLOAD_'+name)
    need(digest(source/'numeric/ACTUAL_SOURCES_NORMAL.json')==SOURCE_REPORT,'ACCEPTED_SOURCE_REPORT')
    need(digest(source/'projection/GEOMETRIC_COEFFICIENT_MAPS.json')==MAP_SHA,'ACCEPTED_SOURCE_MAP')
    accepted=load(source/'numeric/ACTUAL_SOURCES_NORMAL.json')
    grams=accepted['global_cover']['summary']['positive_bundle_bulk_forcing_H_bound']['block_Gram_certificates']
    r=load(report_path)
    need(r['source_authentication']['manifest_sha256']==SOURCE_MANIFEST and r['source_authentication']['files_authenticated']==169,'REPORT_SOURCE_AUTHENTICATION')
    need(r['evaluator_sha256']==digest(producer_script),'PRODUCER_CODE_PIN')
    part=r['partition'];n=part['n'];T=Q(1,4);R=Q(5,4);h=T/n
    need(isinstance(n,int) and n>=1,'MESH_COUNT')
    need(part['h']==str(h),'MESH_WIDTH')
    cells=r['actual_source_cell_enclosures'];seen=set();area=Q(0);moments=[Q(0)]*3
    lo=hi=Q(0);caps={};two_spots=[]
    for cell in cells:
        i,k=cell['i'],cell['k'];need((i,k) not in seen,'DUPLICATE_CELL');seen.add((i,k))
        need(0<=i<n and i<=k<2*n-i,'CELL_INDEX')
        xx=list(map(Q,cell['x']));rr=list(map(Q,cell['r']))
        need(xx==[i*h,(i+1)*h] and rr==[R-T+k*h,R-T+(k+1)*h],'CELL_ENDPOINTS')
        need(0<=xx[0]<=xx[1]<=3 and 0<=rr[0]<=rr[1] and xx[1]+rr[1]<=8,'SOURCE_RECTANGLE_DOMAIN')
        need(rr[0]>=xx[1],'EXTERIOR_ORACLE_BRANCH')
        mm=geometry.rectangle_cone_moments(*xx,*rr,T,R)
        need(mm[0]>0 and Q(cell['exact_area'])==mm[0],'INDEPENDENT_CLIPPED_AREA')
        moments=[a+b for a,b in zip(moments,mm)];area+=mm[0]
        vals=cell['sources'];need(set(vals)==set(EXPECTED),'COMPLETE_23_SOURCE_BLOCKS')
        for field,(count,weight) in EXPECTED.items():
            need(len(vals[field])==count,'SOURCE_BLOCK_LENGTH')
            for j,pair in enumerate(vals[field]):
                a,b=map(Q,pair);need(a<=b,'SOURCE_INTERVAL_ORDER')
                key=field+'/'+str(j);caps[key]=max(caps.get(key,Q(0)),abs(a),abs(b))
        a,b=map(Q,vals['metric_forcing/yy'][0]);lo+=mm[0]*a/2;hi+=mm[0]*b/2
        if not two_spots:two_spots.append(cell)
        if len(two_spots)==1 and mm[0]==h*h:two_spots.append(cell)
    expected_indices={(i,k) for i in range(n) for k in range(i,2*n-i)}
    need(seen==expected_indices and len(cells)==n*(n+1),'GAP_FREE_CONE_PARTITION')
    need(tuple(moments)==(T*T,T**3/3,R*T*T),'EXACT_CONE_MOMENTS')
    need(Q(part['area'])==area==Q(1,16),'CONE_AREA')
    need({k:Q(v) for k,v in r['forcing_coefficient_caps'].items()}==caps,'FORCING_CAPS')
    need(list(map(Q,r['evaluated_Duhamel_mean_interval']))==[lo,hi],'NORMALIZED_DUHAMEL_INTERVAL')
    square=Q(0)
    for prefix,(count,weight) in EXPECTED.items():
        field='scalar' if prefix=='scalar_forcing' else prefix.split('/')[1]
        gg=[[Q(z) for z in row] for row in grams[field]['gram']]
        need(len(gg)==count and all(len(z)==count for z in gg),'GRAM_DIMENSIONS')
        if field=='yy':need(gg[0]==[Q(1),Q(0),Q(0)],'PROBABILITY_MEAN_CONSTANT_BASIS')
        need(grams[field]['multiplicity']==weight,'REGISTERED_FROBENIUS_WEIGHT')
        new=r['forcing_H_bound']['blocks'][field]
        need(new['gram']==grams[field]['gram'] and new['multiplicity']==weight,'GRAM_BASIS_IDENTITY')
        cc=[caps[prefix+'/'+str(j)] for j in range(count)]
        ss=weight*sum((abs(gg[i][j])*cc[i]*cc[j] for i in range(count) for j in range(count)),Q(0))
        need(Q(new['squared_upper'])==ss,'BLOCK_H_BOUND');square+=ss
    FH=Q(r['forcing_H_bound']['F_H_upper'])
    need(FH.denominator==1 and (FH-1)**2<=square<=FH**2,'EXACT_SQRT_ROUNDUP')
    need(Q(r['forcing_H_bound']['squared_upper'])==square,'TOTAL_H_BOUND')
    global_F=Q(accepted['global_cover']['summary']['positive_bundle_bulk_forcing_H_bound']['F_H_upper'])
    need(global_F==17387513156573 and Q(r['accepted_global_F_H_upper'])==global_F,'GLOBAL_SOURCE_H_BOUND')
    F_used=min(FH,global_F);need(Q(r['F_H_upper_used'])==F_used,'LOCAL_GLOBAL_SOURCE_CAP_MINIMUM')
    sigma=Q(203,1000000);q=T*sigma;Wcap=T*F_used/(1-q)
    # Independent simplified expression from t0=1/4 and sigma=203/1e6.
    correction=203*F_used/Q(127993504)
    need(correction==T*T*sigma*Wcap/2,'POINT_TRANSFER_FACTOR')
    coupling=r['full_coupling']
    for k,v in [('A_upper',Q(1,10000)),('B_upper',Q(1,10000)),('D0_upper',Q(1,1000000)),
                ('sigma_upper',sigma),('K_X_operator_upper',T),('q_C_upper',q),('W_X_C_upper',Wcap)]:
        need(Q(coupling[k])==v,'COUPLING_CONSTANT_'+k)
    need(Q(r['coupling_point_value_error_upper'])==correction,'POINT_ERROR')
    final=[lo-correction,hi+correction]
    need(list(map(Q,r['actual_W_yy_probability_mean_interval']))==final,'FULL_RESPONSE_INTERVAL')
    need(r['strict_positive']==(final[0]>0) and r['strict_negative']==(final[1]<0),'SIGN_FLAGS')
    need(final[0]>0,'STRICT_POSITIVE_FULL_RESPONSE')
    rect=r['continuous_rectangle_corollary'];epsilon=Q(1,1000000);Tplus=T+2*epsilon
    need(list(map(Q,rect['x']))==[T-epsilon,T+epsilon] and list(map(Q,rect['r']))==[R-epsilon,R+epsilon],'RECTANGLE_ENDPOINTS')
    need(Q(rect['epsilon'])==epsilon and Q(rect['enlarged_cone_apex_x'])==Tplus and Q(rect['enlarged_cone_apex_r'])==R,'ENLARGED_CONE_APEX')
    need(Q(rect['enlarged_cone_min_r_minus_x'])==R-Tplus>0 and Tplus+R<8 and Tplus<3,'ENLARGED_EXTERIOR_DOMAIN')
    for x in [T-epsilon,T+epsilon]:
        for r0 in [R-epsilon,R+epsilon]:need(abs(r0-R)<=Tplus-x,'WHOLE_RECTANGLE_IN_CONVEX_ENLARGED_CONE')
    qplus=Tplus*sigma;Bplus=Tplus*global_F/(1-qplus);variation=2*epsilon*Bplus
    need(Q(rect['global_F_H_upper_used'])==global_F and Q(rect['q_enlarged_cone_upper'])==qplus,'ENLARGED_CONE_SOURCE_AND_COUPLING')
    need(Q(rect['W_gradient_H_upper'])==Bplus and Q(rect['point_to_rectangle_l1_distance_upper'])==2*epsilon,'DERIVATIVE_WEIGHT_ONE_AND_L1_DISTANCE')
    need(Q(rect['mean_variation_upper'])==variation,'RECTANGLE_MEAN_VARIATION')
    rect_final=[final[0]-variation,final[1]+variation]
    need(list(map(Q,rect['actual_W_yy_probability_mean_interval']))==rect_final,'WHOLE_RECTANGLE_RESPONSE_INTERVAL')
    need(rect['strict_positive'] and rect_final[0]>0,'STRICT_POSITIVE_COMPLETE_RECTANGLE')
    oracle=[]
    if oracle_spots:
        sys.path.insert(0,str(source/'numeric'));sys.path.insert(0,str(source/'projection'))
        import evaluate_actual_sources as frozen
        for cell in two_spots:
            v=frozen.cell(frozen.I(*map(Q,cell['x'])),frozen.I(*map(Q,cell['r'])),source/'projection/GEOMETRIC_COEFFICIENT_MAPS.json',audit_actual=True)
            vals=v['metric_forcing']
            got={'metric_forcing/'+field:[list(map(str,(z.lo,z.hi))) for z in vec] for field,vec in vals.items()}
            got['scalar_forcing']=[list(map(str,(z.lo,z.hi))) for z in v['scalar_forcing']]
            need(got==cell['sources'],'INDEPENDENT_ACTUAL_ORACLE_SPOT')
            oracle.append({'i':cell['i'],'k':cell['k'],'all23_actual_forcing_endpoints_match':True})
    return {'status':'PASS_INDEPENDENT_CAUSAL_RESTRICTED_RESPONSE_AUDIT','producer_report_sha256':digest(report_path),
            'producer_evaluator_sha256':digest(producer_script),'reviewer_source_sha256':digest(Path(__file__)),
            'geometry_source_sha256':digest(HERE/'cone_geometry.py'),
            'source_manifest_sha256':SOURCE_MANIFEST,'source_report_sha256':SOURCE_REPORT,'source_map_sha256':MAP_SHA,
            'cells_checked':len(cells),'independent_cone_moments':list(map(str,moments)),
            'Duhamel_mean_interval':list(map(str,(lo,hi))),'F_H_upper':str(FH),'q_C':str(q),
            'full_coupling_mean_error':str(correction),'full_W_yy_mean_interval':list(map(str,final)),
            'rectangle_mean_variation_upper':str(variation),'whole_rectangle_full_W_yy_mean_interval':list(map(str,rect_final)),
            'strict_positive':True,'actual_oracle_spots':oracle,
            'review_scope':'Independent rational vertical-slice geometry and moments; source rows/23 amplitudes/Frobenius Gram/coupling arithmetic; accepted source proof is a premise, not replayed in full.',
            'limits':'One coordinate metric sphere mean of the full auxiliary forced response before shell contact; no global waveform, physical force, separately constrained solution, or nonlinear existence claim.'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--report',required=True,type=Path);p.add_argument('--source-root',required=True,type=Path)
    p.add_argument('--producer-script',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path);p.add_argument('--oracle-spots',action='store_true');a=p.parse_args()
    need(not a.output.exists(),'REFUSE_OUTPUT_OVERWRITE')
    r=run(a.report,a.source_root,a.producer_script,a.oracle_spots);a.output.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))

if __name__=='__main__':main()
