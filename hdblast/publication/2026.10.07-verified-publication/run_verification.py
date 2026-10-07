"""Run the unchanged, hash-pinned Notes verifier with a portable packet input path."""
import argparse, hashlib, importlib.util, json
from pathlib import Path
ROOT=Path(__file__).resolve().parent
VERIFIER_PIN='869e30cb4ad79138d8b528a076cfeac75dca90b28ecb2a104c6ad3607b380668'
def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--packet-directory',type=Path,default=ROOT.parent/'2026.10.03-verified-source-operator')
 p.add_argument('--saved-response',type=Path,required=True)
 p.add_argument('--public-response',type=Path,required=True)
 p.add_argument('--output',type=Path,required=True)
 a=p.parse_args();source=ROOT/'verify_reviewed_html_serialization.py'
 if hashlib.sha256(source.read_bytes()).hexdigest()!=VERIFIER_PIN:raise RuntimeError('Reviewed verifier pin differs.')
 spec=importlib.util.spec_from_file_location('byte_preserved_notes_verifier',source)
 v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
 v.PACKET=a.packet_directory.resolve()
 report=v.verify(json.loads(a.saved_response.read_text()),published=True,public_saved=json.loads(a.public_response.read_text()))
 a.output.open('x').write(json.dumps(report,indent=2,sort_keys=True)+'\n')
 print(json.dumps({'status':report['status'],'original_literal_guard':report['literal_frozen_guard']['status'],'publication_verified':report['publication_verified']}))
 return 0 if report['status']=='PASS_REVIEWED_HTML_SERIALIZATION_PUBLISHED_RECORD' else 1
if __name__=='__main__':raise SystemExit(main())
