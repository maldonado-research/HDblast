#!/usr/bin/env python3
"""Portable exact report arithmetic/coverage checker; no upstream imports."""
import argparse,copy,hashlib,json,math
from fractions import Fraction as F
from pathlib import Path
OLD_COEFFICIENT_PIN='25ba98e5835960e98248f85cf4493ecbd263180722d8d871cf599922b36209e6'
RECENTERED_COEFFICIENT_PIN='739240283c2bee53afa559c4ce71c1e3356cf2d0ad4e4b0f18b9621c077e20d3'
BRANCH_PIN='024acc3720696322e9fa2f864de4c75b17c88d81ec1a6b80c020e928cefd9cd0'
def need(ok,msg):
 if not ok:raise ValueError(msg)
def sqrtup(q):
 g=1<<112;return F(math.isqrt((q*g*g).__ceil__())+1,g)
def norm(row):return sqrtup(sum(w*x*x for w,x in zip((1,2,1,3,1),row)))
def check(d):
 need(d['means']==['K0','L0','N0','J0','w0'] and d['gram']==[1,2,1,3,1],'registered means/Gram')
 need(F(d['ell'])==F(1,10**6),'physical derivative scale')
 rows=d['rows'];start,stop=d['partition_indices'];need(128<=start<stop<=224 and len(rows)==stop-start==d['whole_cells'],'finite partition')
 complete=start==128 and stop==224;need(d['full_footprint_covered']==complete,'coverage label')
 caps=[[[F(0)]*5 for _ in range(3)] for _ in range(3)]
 for index,row in zip(range(start,stop),rows):
  need(row['index']==index,'partition order');need(list(map(F,row['r']))==[F(index,128),F(index+1,128)],'whole closed cell')
  cc=row['channel_component_error_caps'];need(len(cc)==3 and all(len(c)==3 and all(len(v)==5 for v in c) for c in cc),'channel dimensions')
  for c in range(3):
   for k in range(3):
    vals=list(map(F,cc[c][k]));need(all(v>=0 for v in vals),'negative error cap');caps[c][k]=[max(a,b) for a,b in zip(caps[c][k],vals)]
 if 'channel_component_error_caps' in d:need(caps==[[list(map(F,row)) for row in c] for c in d['channel_component_error_caps']],'component maxima')
 hs=[[norm(row) for row in chan] for chan in caps]
 need(hs==[list(map(F,row)) for row in d['channel_H_caps']],'Gram norm caps')
 total=[sum(hs[c][k] for c in range(3)) for k in range(3)];need(total==list(map(F,d['total_H_caps'])),'three-channel total')
 free=max((total[0]+F(3,2)*total[2])/3,total[1]+total[2]);need(free==F(d['free_wave_X_cap']),'free-wave exact formula')
 numeric=free<=10**28
 if 'numeric_cap_below_allocation' in d:need(d['numeric_cap_below_allocation']==numeric,'numerical allocation comparison')
 branch=d['case']=='branch';gate=complete and branch and numeric
 need(d['passes_1e28']==gate,'whole certified initial gate label')
 if branch:need(d['branch_input_certificate']['sha256']==BRANCH_PIN,'independent beta input pin')
 need(not d['actual_error_lower_bound_proved'],'upper bound scope')
 need(d['reauthenticated_after'],'complete after authentication')
 expected={'P':(107,'8bf173293b5cecf7f78e138ac5004abd9c7b1821ef29addb7ecfbc2a6409b4f6'),'N':(111,'928f3e4e58d222a3c7e1274babf9e94930ec501b74cf6ba04209f97c6bbded44'),'G':(168,'3ff1f95c151d2da6a1ab5e9b54718c2bb325c10d9c62a7ed9f97ef2bb36f0bb2'),'C':(48,'72924328a4e79cc77f0da90999b8dd56956bfdebd5e614e15c3c0498c3a8ae76')}
 for k,(n,pin) in expected.items():need(d['authentication_summary'][k]['files_verified']==n and d['authentication_summary'][k]['manifest_sha256']==pin,'upstream authentication summary '+k)
 need(d['current_authentication']['files_including_manifest']==49 and d['current_authentication']['manifest_sha256']=='cffb213a69f7f9a27932e990ffdecb7569d1007e5e5baeda708d8dc39d0cce72','49-file dependency pin')
 coeff=d['representation_background_coefficients_ascending'];need(set(coeff)=={'rho','eta','p','a','alpha'} and all(len(v)==5 for v in coeff.values()),'single degree4 nominal representation')
 for v in coeff.values():list(map(F,v))
 cp=hashlib.sha256(json.dumps(coeff,sort_keys=True,separators=(',',':')).encode()).hexdigest();need(cp==d['representation_background_coefficients_sha256'],'coefficient content hash')
 need(cp==(RECENTERED_COEFFICIENT_PIN if d['candidate_changed'] else OLD_COEFFICIENT_PIN),'independently fixed candidate coefficient pin')
 design=[[F(5*10**27),F(2*10**27),F(5*10**26)],[F(10**28),F(4*10**27),F(10**27)],[F(5*10**27),F(2*10**27),F(5*10**26)]]
 failures=[{'channel':c,'quantity':k,'actual_upper':str(hs[c][k]),'illustrative_design_upper':str(design[c][k])} for c in range(3) for k in range(3) if hs[c][k]>design[c][k]]
 total_design=[F(2*10**28),F(8*10**27),F(2*10**27)]
 total_failures=[{'quantity':k,'actual_upper':str(total[k]),'illustrative_design_upper':str(total_design[k])} for k in range(3) if total[k]>total_design[k]]
 return {'status':'PASS_PORTABLE_DIFFERENCE_LEDGER','whole_cells':len(rows),'candidate_changed':d['candidate_changed'],'free_wave_X_cap':str(free),'initial_gate_passes':gate,'illustrative_design_component_failures':failures,'illustrative_design_total_failures':total_failures,'scope':'Checks reported exact arithmetic, partition, pins and pass/fail honesty; does not reexecute native interval formulas or reauthenticate unavailable dependency trees.'}
def controls(d):
 mutations=[]
 def rejected(name,edit):
  q=copy.deepcopy(d);edit(q)
  try:check(q)
  except (ValueError,KeyError,TypeError,IndexError):mutations.append(name);return
  raise ValueError('undetected mutation '+name)
 rejected('wrong Gram weight',lambda q:q['gram'].__setitem__(3,2))
 rejected('missing whole cell',lambda q:q['rows'].pop())
 rejected('shrunken free cap',lambda q:q.__setitem__('free_wave_X_cap','1'))
 rejected('lost physical ell',lambda q:q.__setitem__('ell','1'))
 rejected('zeroed packet channel norm',lambda q:q['channel_H_caps'][1].__setitem__(1,'0'))
 rejected('changed nominal coefficient',lambda q:q['representation_background_coefficients_ascending']['p'].__setitem__(0,'-1/10000'))
 if d['case']=='branch':rejected('changed beta provenance',lambda q:q['branch_input_certificate'].__setitem__('sha256','0'*64))
 rejected('false gate decision',lambda q:q.__setitem__('passes_1e28',not q['passes_1e28']))
 return mutations
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--report',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();d=json.loads(args.report.read_text());r=check(d);r['mutation_controls_rejected']=controls(d);args.output.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r))
if __name__=='__main__':main()
