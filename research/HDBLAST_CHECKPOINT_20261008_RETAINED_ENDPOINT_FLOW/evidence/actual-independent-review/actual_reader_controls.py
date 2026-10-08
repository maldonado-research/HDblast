"""Export-only negative controls and exact representation boundary controls."""
from pathlib import Path
from fractions import Fraction as Q
import gzip
import hashlib
import json
import sys
BASE=Path(__file__).resolve().parent
raw=(BASE/'actual_allnode_readback.py').read_bytes();namespace={'__name__':'control_reader','__file__':str(BASE/'actual_allnode_readback.py')}
exec(compile(raw,namespace['__file__'],'exec'),namespace)
check=namespace['represented_binary80'];exact=namespace['exact'];EPS=namespace['EPS']

def main(directory,label):
 directory=Path(directory);out=BASE/label;out.mkdir();cases=[]
 positives=[('zero',Q(0),False),('negative_zero',Q(0),True),('minimum_subnormal',Q(1,1<<16445),False),
 ('maximum_subnormal',Q((1<<63)-1,1<<16445),False),('minimum_normal',Q(1,1<<16382),False),
 ('maximum_finite',Q(((1<<64)-1)*(1<<16320)),False),('negative_finite',Q(-3,8),False)]
 for name,value,flag in positives:
  check(value,flag);cases.append({'case':name,'passed':True,'expected':'accept exact representable scalar'})
 negatives=[('below_subnormal',Q(1,1<<16446),False),('above_maximum',Q(1<<16384),False),
 ('excess_significand',Q((1<<65)+1,1<<64),False),('nondyadic',Q(1,3),False),
 ('nonzero_negative_zero',Q(1),True),('double_normalized',Q(1)/EPS,False)]
 for name,value,flag in negatives:
  try:check(value,flag)
  except RuntimeError as exc:cases.append({'case':name,'passed':True,'error':str(exc)})
  else:raise RuntimeError('negative representation control accepted:'+name)
 for text in ['2/2','-0/1','1.0','+1/1','01/1']:
  try:exact(text)
  except RuntimeError:cases.append({'case':'noncanonical_'+text,'passed':True})
  else:raise RuntimeError('noncanonical ratio accepted')
 data=(directory/'DATA.json').read_bytes();summary=(directory/'SCIENCE_SUMMARY.json').read_bytes()
 with gzip.open(directory/'NODE_CERTIFICATES.jsonl.gz','rb') as stream:first=json.loads(stream.readline())
 plans={'double_epsilon':'not dyadic','negative_zero_nonzero':'negative zero flag on nonzero',
        'noncanonical_ratio':'noncanonical exact ratio','overwide_target':'target width gate',
        'canonical_drift':'canonical drift differs','delta_tamper':'saved-minus-target serialization differs',
        'wrong_prefix':'cutoff/prefix mask differs','manufactured_scope':'actual exact represented-anchor endpoint scope'}
 for name,expected in plans.items():
  fixture=out/name;fixture.mkdir();row=json.loads(json.dumps(first));document=json.loads(data)
  if name=='double_epsilon':row['incoming_U'][0]=namespace['canonical'](Q(row['incoming_U'][0])/EPS)
  if name=='negative_zero_nonzero':row['negative_zero_flags'][2]=True
  if name=='noncanonical_ratio':
   q=Q(row['k']);row['k']=str(q.numerator*2)+'/'+str(q.denominator*2)
  if name=='overwide_target':
   pair=row['endpoints']['-4']['target_dyadic96']['U']['real'];pair[1]=pair[0]+(1<<96)//10**10
  if name=='canonical_drift':row['endpoints']['-4']['canonical_drift']=namespace['canonical'](Q(row['endpoints']['-4']['canonical_drift'])+1)
  if name=='delta_tamper':row['endpoints']['-4']['delta_U'][0][0]=namespace['canonical'](Q(row['endpoints']['-4']['delta_U'][0][0])+1)
  if name=='wrong_prefix':row['k']='100/1'
  if name=='manufactured_scope':document['manufactured']=True
  (fixture/'DATA.json').write_text(json.dumps(document));(fixture/'SCIENCE_SUMMARY.json').write_bytes(summary)
  with gzip.open(fixture/'NODE_CERTIFICATES.jsonl.gz','wb') as stream:stream.write((json.dumps(row)+'\n').encode())
  try:namespace['run'](fixture,label+'/'+name+'/FORBIDDEN_SUCCESS.json')
  except RuntimeError as exc:
   message=str(exc)
   if expected not in message:raise RuntimeError('wrong rejection branch '+name+': '+message)
   cases.append({'case':name,'passed':True,'error':message})
  else:raise RuntimeError('export mutation accepted:'+name)
 result={'status':'PASS_INDEPENDENT_ACTUAL_READER_CONTROLS','reader_sha256':hashlib.sha256(raw).hexdigest(),
 'checks':len(cases),'cases':cases,'retained_arrays_opened':0,'physical_source_calls':0,'production_modules_imported':0,
 'input_source':'Existing exported first node only; no retained bytes or second evolution.'}
 (out/'RESULT.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');print(json.dumps(result,sort_keys=True))
if __name__=='__main__':main(*sys.argv[1:])
