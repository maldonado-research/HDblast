#!/usr/bin/env python3
"""Compare complete reports, excluding exactly two location-only fields."""
import argparse,copy,hashlib,json
from pathlib import Path
def normalize(d):
 d=copy.deepcopy(d)
 # These locations differ between native externally pinned and portable
 # embedded-hash authentication. Their manifest/content pins stay compared.
 d['current_authentication'].pop('external_pin')
 if 'branch_input_certificate' in d:d['branch_input_certificate'].pop('path')
 return d
def main():
 ap=argparse.ArgumentParser();ap.add_argument('reports',type=Path,nargs='+');ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
 if len(a.reports)<2:raise ValueError('two reports required')
 data=[normalize(json.loads(p.read_text())) for p in a.reports]
 if any(d!=data[0] for d in data[1:]):raise ValueError('report mismatch outside the two declared location fields')
 r={'status':'PASS_COMPLETE_REPORT_COMPARISON','reports':[{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in a.reports],'ignored_only':['current_authentication.external_pin','branch_input_certificate.path'],'normalized_sha256':hashlib.sha256(json.dumps(data[0],sort_keys=True,separators=(',',':')).encode()).hexdigest()};a.output.write_text(json.dumps(r,indent=2)+'\n');print(r['status'])
if __name__=='__main__':main()
