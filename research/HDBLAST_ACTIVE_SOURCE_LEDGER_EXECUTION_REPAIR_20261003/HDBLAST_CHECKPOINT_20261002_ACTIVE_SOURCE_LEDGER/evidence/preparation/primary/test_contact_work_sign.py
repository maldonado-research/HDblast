"""Fabricated contact-work truth and deliberately wrong-sign rejection."""
from pathlib import Path
from decimal import Decimal
from fractions import Fraction
import json,sys
root=Path(__file__).parent
payload=json.loads((root/'SYNTHETIC_SCHEMA_FINAL.json').read_text())
checks=[]
for control in payload['case']['controls']:
 for dps,rows in control['levels'].items():
  for row in rows:
   value=lambda key:Fraction(Decimal(row[key]))
   reconstructed=value('E_reconstruction')
   exact_formula=(value('canonical_constant_drift')+Fraction(Decimal(row['contact_R'][-1]))
                    -Fraction(Decimal(row['contact_R'][0]))+value('analytic_source_moment')*value('source_work_integral')
                    -value('full_contact_ledger_integral'))
   if abs(exact_formula-reconstructed)>Fraction('1e-12'):raise ValueError('Fabricated direct contact primitive bookkeeping differs')
   wrong_sign=exact_formula+2*value('full_contact_ledger_integral')
   if abs(exact_formula)>Fraction('2e-7'):raise ValueError('Fabricated contact primitive fails reconstruction gate')
   if abs(wrong_sign)<=Fraction('2e-7'):raise ValueError('Deliberately changed contact-work sign escaped reconstruction gate')
   checks.append({'order':control['quadrature_order'],'dps':dps,'K':row['K'],
                  'correct_reconstruction_exact_rational':str(exact_formula),'wrong_sign_exact_rational':str(wrong_sign)})
result={'status':'PASS','check_count':len(checks),'python_optimization':sys.flags.optimize,
        'producer_sha256':payload['producer_sha256'],'study_source_or_arrays_used':False,'checks':checks}
out=root/('SYNTHETIC_CONTACT_WORK_SIGN_'+('OPTIMIZED' if sys.flags.optimize else 'NORMAL')+'_'+payload['producer_sha256'][:12]+'.json')
if out.exists():raise ValueError('Output must be fresh')
out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='checks'}))
