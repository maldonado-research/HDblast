#!/usr/bin/env python3
"""Synthetic offline saved-record guard checks, plus one authorized saved draft."""
import argparse
import copy
import hashlib
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

sys.dont_write_bytecode = True


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + '\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--verifier', type=Path, required=True)
    parser.add_argument('--results', type=Path, required=True)
    parser.add_argument('--fixtures-root', type=Path, required=True)
    parser.add_argument('--saved-draft', type=Path, required=True)
    args = parser.parse_args()
    args.fixtures_root.mkdir(parents=True, exist_ok=True)
    candidate = args.verifier.parent
    sys.path.insert(0, str(candidate))
    spec = importlib.util.spec_from_file_location('verifier_under_review', args.verifier)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    manifest_path = candidate/'FILE_MANIFEST.json'
    manifest = json.loads(manifest_path.read_text())
    request = json.loads((candidate/manifest['metadata_request']['path']).read_text())
    selected_id = 90000001
    rows = manifest['inherited_files'] + manifest['new_files']
    complete = {'id': selected_id, 'conceptrecid':'17088132','conceptdoi':'10.5281/zenodo.17088132','submitted':False,'state':'unsubmitted','metadata':copy.deepcopy(request['metadata']), 'files':[{'filename':f['filename'],'filesize':f['bytes'],'checksum':f['md5']} for f in rows]}
    complete['metadata']['prereserve_doi'] = {'doi':'10.5281/zenodo.90000001','recid':selected_id}
    results = []

    def check(name, data, status='FAIL_SAVED_RECORD_VALIDATION', published=False, expected_id=selected_id, expected_manifest=None):
        report, exception = None, None
        try:
            report = module.verify(data, expected_manifest or manifest, request, expected_id, published)
        except Exception as exc:
            exception = f'{type(exc).__name__}: {exc}'
        passed = exception is None and report['status']==status
        if passed and status.startswith('PASS_'):
            passed = report['observed_files']==19 and not report['failures'] and not report['metadata_mismatch_fields'] and report['publication_performed_by_this_tool'] is False and 'no network, upload or publish action' in report['scope']
        if passed and status=='FAIL_SAVED_RECORD_VALIDATION':
            passed = bool(report['failures'])
        results.append({'case':name,'passed':passed,'expected_status':status,'report':report,'exception':exception})
        return report

    def mutated(mutator):
        data = copy.deepcopy(complete)
        mutator(data)
        return data

    check('raw_complete_draft_with_reserved_doi', complete, 'PASS_COMPLETE_SAVED_DRAFT')
    check('receipt_wrapper_complete_draft', {'status':200,'data':complete}, 'PASS_COMPLETE_SAVED_DRAFT')
    check('integer_conceptrecid', mutated(lambda d:d.update(conceptrecid=17088132)), 'PASS_COMPLETE_SAVED_DRAFT')
    alternate = mutated(lambda d:[f.update(size=f.pop('filesize'),checksum='MD5:'+f['checksum'].upper()) for f in d['files']])
    check('alternate_size_key_and_prefixed_uppercase_md5', alternate, 'PASS_COMPLETE_SAVED_DRAFT')
    reordered = mutated(lambda d:[d['metadata'][key].reverse() for key in ['keywords','references','related_identifiers']])
    check('setlike_metadata_list_reordering', reordered, 'PASS_COMPLETE_SAVED_DRAFT')
    published = mutated(lambda d:d.update(submitted=True, state='done', doi='10.5281/zenodo.'+str(selected_id)))
    check('raw_submitted_done_published', published, 'PASS_VERIFIED_PUBLISHED_RECORD', published=True)
    check('receipt_wrapper_submitted_done_published', {'data':published}, 'PASS_VERIFIED_PUBLISHED_RECORD', published=True)
    check('reserved_doi_draft_not_published', complete, published=True)
    check('published_requires_explicit_mode', published)
    missing_published_doi = copy.deepcopy(published)
    missing_published_doi.pop('doi')
    missing_published_doi['metadata'].pop('prereserve_doi')
    check('published_missing_top_level_doi', missing_published_doi, published=True)
    mismatched_published_doi = copy.deepcopy(published)
    mismatched_published_doi['doi'] = '10.5281/zenodo.'+str(selected_id+1)
    check('published_mismatched_top_level_doi', mismatched_published_doi, published=True)
    reserved_only_published = copy.deepcopy(published)
    reserved_only_published.pop('doi')
    check('published_reserved_metadata_doi_only', reserved_only_published, published=True)
    for name, changes in [
        ('wrong_family_conceptdoi',{'conceptdoi':'10.5281/zenodo.22922927'}),
        ('wrong_family_conceptrecid',{'conceptrecid':'22922927'}),
        ('missing_conceptdoi',{'conceptdoi':None}),
        ('wrong_record_id',{'id':selected_id+1}),
        ('string_record_id',{'id':str(selected_id)}),
        ('missing_record_id',{'id':None}),
        ('wrong_draft_state',{'state':'done'}),
        ('draft_submitted_true',{'submitted':True}),
        ('draft_submitted_null',{'submitted':None}),
        ('draft_submitted_zero',{'submitted':0}),
        ('empty_metadata',{'metadata':{}}),
        ('null_metadata',{'metadata':None}),
        ('list_metadata',{'metadata':[]}),
        ('missing_files',{'files':[]}),
        ('null_files',{'files':None}),
        ('object_files',{'files':{}}),
    ]:
        check(name, mutated(lambda d, changes=changes:d.update(changes)))
    for name, changes in [
        ('published_submitted_false_done',{'submitted':False,'state':'done'}),
        ('published_submitted_true_unsubmitted',{'submitted':True,'state':'unsubmitted'}),
        ('published_submitted_one_done',{'submitted':1,'state':'done'}),
        ('published_submitted_true_processing',{'submitted':True,'state':'processing'}),
    ]:
        published_changes = copy.deepcopy(published)
        published_changes.update(changes)
        check(name, published_changes, published=True)
    check('missing_one_remote_file', mutated(lambda d:d['files'].pop()))
    check('extra_remote_file', mutated(lambda d:d['files'].append({'filename':'unlisted.bin','filesize':1,'checksum':'0'*32})))
    check('duplicate_remote_filename', mutated(lambda d:d['files'].append(copy.deepcopy(d['files'][0]))))
    check('malformed_remote_file', mutated(lambda d:d['files'].append('invalid')))
    check('nonstring_remote_filename', mutated(lambda d:d['files'][0].update(filename=123)))
    check('boolean_byte_count', mutated(lambda d:d['files'][0].update(filesize=True)))
    check('sha256_is_not_md5', mutated(lambda d:d['files'][0].update(checksum='sha256:'+'0'*64)))
    for i, row in enumerate(rows):
        classification = 'inherited' if i<10 else 'new'
        check(f'{classification}_{i:02}_changed_bytes', mutated(lambda d,i=i:d['files'][i].update(filesize=d['files'][i]['filesize']+1)))
        check(f'{classification}_{i:02}_changed_md5', mutated(lambda d,i=i:d['files'][i].update(checksum='0'*32)))
    for key in request['metadata']:
        check('metadata_field_missing_'+key, mutated(lambda d,key=key:d['metadata'].pop(key)))
    check('nonobject_saved_response', [])
    check('wrapper_nonobject_deposition', {'data':[]})
    duplicate_manifest = copy.deepcopy(manifest)
    duplicate_manifest['new_files'][0]['filename']=duplicate_manifest['inherited_files'][0]['filename']
    check('nonunique_expected_inventory', complete, expected_manifest=duplicate_manifest)

    # Exercise command-line loading, pinned local metadata/manifest, output,
    # and strict --published state checks without copying a scientific file.
    def cli_case(name, data, published_mode, expected_status, expected_exit, wrong_pin=False, existing=False):
        saved_path = args.fixtures_root/(name+'-saved.json')
        output_path = args.fixtures_root/(name+'-report.json')
        write_json(saved_path,data)
        if existing:
            output_path.write_text('preserve existing review evidence\n')
        command=[sys.executable,'-B']
        if sys.flags.optimize:
            command.append('-O')
        command += [str(args.verifier),'--saved-response',str(saved_path),'--expected-record-id',str(selected_id),'--expected-manifest-sha256','0'*64 if wrong_pin else sha(manifest_path),'--output',str(output_path)]
        if published_mode:
            command.append('--published')
        process=subprocess.run(command,capture_output=True,text=True)
        report=json.loads(output_path.read_text()) if output_path.exists() and not existing else None
        passed=process.returncode==expected_exit and (report['status']==expected_status if expected_status else report is None)
        if existing:
            passed=passed and output_path.read_text()=='preserve existing review evidence\n'
        if wrong_pin:
            passed=passed and not output_path.exists()
        results.append({'case':name,'passed':passed,'expected_status':expected_status,'exit_code':process.returncode,'report':report,'stdout':process.stdout.strip(),'existing_output_preserved':not existing or output_path.read_text()=='preserve existing review evidence\n'})
    cli_case('cli_raw_complete_draft',complete,False,'PASS_COMPLETE_SAVED_DRAFT',0)
    cli_case('cli_wrapper_complete_draft',{'data':complete},False,'PASS_COMPLETE_SAVED_DRAFT',0)
    cli_case('cli_reserved_doi_not_published',complete,True,'FAIL_SAVED_RECORD_VALIDATION',1)
    cli_case('cli_submitted_done_published',published,True,'PASS_VERIFIED_PUBLISHED_RECORD',0)
    cli_case('cli_wrong_manifest_pin',complete,False,None,1,wrong_pin=True)
    cli_case('cli_existing_output',complete,False,None,1,existing=True)
    cli_case('cli_null_metadata',mutated(lambda d:d.update(metadata=None)),False,'FAIL_SAVED_RECORD_VALIDATION',1)

    saved_draft = json.loads(args.saved_draft.read_text())
    draft_data = saved_draft.get('data',saved_draft)
    actual_report = check('actual_saved_main_draft_is_incomplete',saved_draft,expected_id=draft_data['id'])
    if actual_report and (actual_report['observed_files']!=10 or len(actual_report['metadata_mismatch_fields'])!=7):
        results[-1]['passed']=False
    report={'status':'PASS_FOCUSED_SAVED_RECORD_VERIFIER_REVIEW' if all(r['passed'] for r in results) else 'FAIL_FOCUSED_SAVED_RECORD_VERIFIER_REVIEW','verifier_sha256':sha(args.verifier),'preparation_helper_sha256':sha(candidate/'prepare_candidate.py'),'manifest_sha256':sha(manifest_path),'metadata_sha256':sha(candidate/manifest['metadata_request']['path']),'optimization':sys.flags.optimize,'cases':len(results),'passed':sum(r['passed'] for r in results),'failures':[r['case'] for r in results if not r['passed']],'results':results,'actual_saved_draft_sha256':sha(args.saved_draft),'scope':'Saved-response comparisons and synthetic CLI checks only; no network, science, upload or publication.'}
    write_json(args.results,report)
    print(json.dumps({key:report[key] for key in ['status','verifier_sha256','optimization','cases','passed','failures']}))
    return 0 if not report['failures'] else 1


if __name__=='__main__':
    raise SystemExit(main())
