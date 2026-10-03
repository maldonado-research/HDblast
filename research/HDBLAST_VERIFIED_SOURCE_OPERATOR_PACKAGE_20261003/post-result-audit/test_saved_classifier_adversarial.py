#!/usr/bin/env python3
"""Post-result numerical/classifier mutations of copies; no source calls."""
from copy import deepcopy
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import review_saved_operator_results as own


def main():
    base=Path('/workspace/hdblast-research-work/verified-integration-20261003')
    p=own.load((base/'actual/primary/OUTPUT.json').read_bytes())['payload']
    i=own.load((base/'actual/independent/OUTPUT.json').read_bytes())['payload']
    outer=own.load((base/'actual/primary/EXECUTION.json').read_bytes())
    command=outer['command'];root=Path(command[command.index('--root')+1])
    source=root/'protocol/operator_probe_contract.py'
    registration=own.load((root/'FULL_REGISTRATION.json').read_bytes())
    expected_pin=registration['files']['protocol/operator_probe_contract.py']['sha256']
    own.require(own.sha(source.read_bytes())==expected_pin,'Frozen classifier source pin differs')
    spec=importlib.util.spec_from_file_location('frozen_output_validator_only',source)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    checks=[];negative=[]
    def rejected(name,mutate):
        value=deepcopy(i);mutate(value)
        try:module.validate_payload(value)
        except (ValueError,TypeError,KeyError):negative.append(name)
        else:raise RuntimeError('Invalid execution/schema mutation survived: '+name)
    own.require(module.validate_pair(p,i)['status']=='PASS_BOTH_REGISTERED_PROBE_OUTPUTS','Actual pair classifier differs')
    checks.append('actual saved pair classified PASS')
    rejected('truncated output universe',lambda x:x['panel_rows'].pop())
    rejected('bool archive count',lambda x:x.__setitem__('archive_arrays_decoded',False))
    rejected('bool schema version',lambda x:x.__setitem__('schema_version',True))
    rejected('wrong rational momentum',lambda x:x['whole_rows'][1].__setitem__('momentum','1/1099511627775'))
    rejected('incorrect complete radius',lambda x:x['whole_rows'][1]['total_absolute_radii'].__setitem__('Mu','0/1'))
    rejected('noncanonical rational endpoint',lambda x:x['whole_rows'][1]['moments']['Mu']['real'].__setitem__('lo','0'))
    rejected('wrong row order',lambda x:x['panel_rows'].reverse())
    def nonreal(x):
        row=x['whole_rows'][1]
        row['moments']['M0']['imag']={'lo':'1/1','hi':'1/1'}
    rejected('nonreal M0',nonreal)
    # Valid wide enclosures preserve schema; they are a scientific unresolved result.
    wide=deepcopy(i);row=wide['whole_rows'][1];box=row['moments']['Mu']
    lo=Q(box['real']['lo'])-Q(1,10**20);hi=Q(box['real']['hi'])+Q(1,10**20)
    box['real']={'lo':own.canonical(lo),'hi':own.canonical(hi)}
    row['total_absolute_radii']['Mu']=own.canonical(module.absolute_radius(box))
    own.require(module.validate_pair(p,wide)['status']=='UNRESOLVED_PROBE_RADIUS','Wide valid result incorrectly passes/classifies')
    checks.append('valid wide enclosure classified UNRESOLVED_PROBE_RADIUS (completed reviewer maps to UNRESOLVED_CERTIFICATE)')
    shifted=deepcopy(i);row=shifted['whole_rows'][1];box=row['moments']['Mu']
    box['real']={name:own.canonical(Q(value)+Q(1,10**8)) for name,value in box['real'].items()}
    own.require(module.validate_pair(p,shifted)['status']=='CERTIFICATE_CONSISTENCY_FAILURE',
                'Valid disjoint enclosure incorrectly passes/classifies')
    checks.append('valid disjoint enclosure classified CERTIFICATE_CONSISTENCY_FAILURE')
    # An exact output of the own audit is recomputed without the frozen helper.
    pr,_=own.inspect_payload(p,'primary');ir,_=own.inspect_payload(i,'independent')
    for field in pr:
        for key in pr[field]:
            for axis in ('real','imag'):
                a,b=pr[field][key][axis],ir[field][key][axis]
                own.require(max(a[0],b[0])<=min(a[1],b[1]),'Independent component intersection failed')
    checks.append('all saved rectangle intersections independently checked')
    for q in (Q(-3,7),Q(0),Q(1,3),Q(1,10**45),Q(-1,10**45)):
        lo=Q(own.outward_decimal(q,lower=True));hi=Q(own.outward_decimal(q,lower=False))
        own.require(lo<=q<=hi,'Inward decimal number formatting')
    checks.append('signed and nearzero outward decimal formatting exact')
    report={'status':'PASS_POST_RESULT_ADVERSARIAL_CLASSIFIER_REVIEW',
            'checks':checks,'invalid_schema_mutations_rejected':negative,
            'frozen_classifier_sha256':expected_pin,'source_callbacks':0,'retained_arrays_decoded':0,
            'modified_original_result_files':False,'fixture_scope':'mutated in-memory copies of already saved numerical output only',
            'own_reviewer_source_sha256':own.sha(Path(own.__file__).read_bytes()),
            'test_source_sha256':own.sha(Path(__file__).read_bytes())}
    output=Path(sys.argv[1]);own.require(not output.exists(),'Fresh receipt required')
    output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':report['status'],'checks':len(checks),'rejected_mutations':len(negative)}))


if __name__=='__main__':main()
