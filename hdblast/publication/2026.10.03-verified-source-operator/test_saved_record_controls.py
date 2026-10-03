"""Synthetic publication-guard boundary controls; no network or scientific run."""
import argparse
import copy
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location('candidate_saved_response_guard', ROOT / 'verify_saved_record.py')
guard = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(guard)

def require(value, message):
    if not value:
        raise RuntimeError(message)

def capture(files, metadata, published=False, public=False):
    endpoint = ('https://zenodo.org/api/records/' if public else 'https://zenodo.org/api/deposit/depositions/') + '23114217'
    data = {'id':23114217, 'conceptrecid':'17088132', 'conceptdoi':'10.5281/zenodo.17088132', 'owner':1386319,
            'submitted':published, 'state':'done' if published else 'unsubmitted', 'metadata':copy.deepcopy(metadata), 'files':[]}
    for x in files:
        data['files'].append({'key' if public else 'filename':x['filename'], 'size' if public else 'filesize':x['bytes'], 'checksum':'md5:'+x['md5']})
    if published:
        data['doi']='10.5281/zenodo.23114217'
    if public:
        data['links']={'self':endpoint,'parent':'https://zenodo.org/api/records/17088132',
                       'parent_doi':'https://doi.org/10.5281/zenodo.17088132',
                       'latest':'https://zenodo.org/api/records/23114217/versions/latest'}
        data['metadata']['license']={'id':data['metadata']['license']}
        data['metadata']['resource_type']={'type':data['metadata'].pop('upload_type')}
    return {'method':'GET','http_status':200,'requested_url':endpoint,'final_url':endpoint,'redirects':[],
            'authenticated_request':not public,'data':data,'synthetic_control':True}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    manifest=guard.read_json(ROOT/'FILE_MANIFEST.json')
    request=guard.read_json(ROOT/'METADATA.json')
    guard.validate_candidate(manifest,ROOT)
    inherited=manifest['inherited_files']
    all_files=inherited+manifest['new_files']
    meta=capture(inherited,request['metadata'])
    complete=capture(all_files,request['metadata'])
    done=capture(all_files,request['metadata'],True)
    public=capture(all_files,request['metadata'],True,True)
    records=[]
    for name,saved,options,status in (
        ('valid_metadata_and_inherited',meta,{'metadata_only':True},'PASS_EXACT_SAVED_DRAFT_METADATA'),
        ('valid_complete_23_file_draft',complete,{},'PASS_COMPLETE_SAVED_DRAFT'),
        ('valid_done_and_direct_public_witness',done,{'published':True,'public_saved':public},'PASS_VERIFIED_PUBLISHED_RECORD')):
        report=guard.verify(saved,manifest,request,23114217,**options)
        require(report['status']==status,name+' was rejected')
        records.append({'name':name,'expected':status,'observed':report['status']})
    mutations=[
        ('no_authenticated_witness',lambda x:x.update(authenticated_request=False)),
        ('wrong_owner',lambda x:x['data'].update(owner=1386320)),
        ('bool_owner',lambda x:x['data'].update(owner=True)),
        ('wrong_record',lambda x:x['data'].update(id=22347452)),
        ('bool_record',lambda x:x['data'].update(id=True)),
        ('wrong_concept',lambda x:x['data'].update(conceptrecid='22922927',conceptdoi='10.5281/zenodo.22922927')),
        ('empty_metadata',lambda x:x['data'].update(metadata={})),
        ('wrong_title',lambda x:x['data']['metadata'].update(title='Unreviewed title')),
        ('wrong_semantic_version',lambda x:x['data']['metadata'].update(version='2026.10.02-ledger-checkpoint')),
        ('wrong_publication_date',lambda x:x['data']['metadata'].update(publication_date='2026-10-02')),
        ('missing_scientific_limit',lambda x:x['data']['metadata'].update(description='Unbounded claim')),
        ('old_inherited_only_draft',lambda x:x['data'].update(files=copy.deepcopy(meta['data']['files']))),
        ('old_22_file_inventory',lambda x:x['data']['files'].pop()),
        ('extra_file',lambda x:x['data']['files'].append({'filename':'extra.zip','filesize':1,'checksum':'md5:'+'0'*32})),
        ('duplicate_filename',lambda x:x['data']['files'].append(copy.deepcopy(x['data']['files'][0]))),
        ('wrong_size',lambda x:x['data']['files'][0].update(filesize=1)),
        ('bool_size',lambda x:x['data']['files'][0].update(filesize=True)),
        ('wrong_checksum',lambda x:x['data']['files'][0].update(checksum='md5:'+'0'*32)),
        ('missing_checksum',lambda x:x['data']['files'][0].pop('checksum')),
        ('inherited_file_replaced',lambda x:x['data']['files'][0].update(filename='replaced.zip')),
        ('malformed_file',lambda x:x['data']['files'].append([])),
        ('reserved_doi_and_wrong_state',lambda x:x['data'].update(state='done',doi='10.5281/zenodo.23114217')),
        ('wrong_http_status',lambda x:x.update(http_status=201)),
        ('redirected_witness',lambda x:x.update(redirects=[{'status':302,'url':x['final_url']}])),
        ('wrong_endpoint',lambda x:x.update(requested_url='https://zenodo.org/api/deposit/depositions/22347452')),
        ('tombstone',lambda x:x['data'].update(tombstone={'reason':'removed'})),
        ('error_even_with_http200',lambda x:x.update(error='Incomplete response')),
    ]
    for name,mutation in mutations:
        saved=copy.deepcopy(complete)
        mutation(saved)
        report=guard.verify(saved,manifest,request,23114217)
        require(report['status']=='FAIL_SAVED_RECORD_VALIDATION',name+' was accepted')
        records.append({'name':name,'expected':'FAIL_SAVED_RECORD_VALIDATION','observed':report['status']})
    for name,mutation in (
        ('public_wrong_family_link',lambda x:x['data']['links'].update(parent='https://zenodo.org/api/records/22922927')),
        ('public_removed',lambda x:x.update(http_status=410)),
        ('public_missing_file',lambda x:x['data']['files'].pop()),
        ('public_wrong_metadata',lambda x:x['data']['metadata'].update(description='Other result')),
        ('public_unsubmitted',lambda x:x['data'].update(submitted=False,state='unsubmitted'))):
        witness=copy.deepcopy(public)
        mutation(witness)
        report=guard.verify(done,manifest,request,23114217,published=True,public_saved=witness)
        require(report['status']=='FAIL_SAVED_RECORD_VALIDATION',name+' was accepted')
        records.append({'name':name,'expected':'FAIL_SAVED_RECORD_VALIDATION','observed':report['status']})
    report=guard.verify(done,manifest,request,23114217,published=True,metadata_only=True,public_saved=public)
    require(report['status']=='FAIL_SAVED_RECORD_VALIDATION','Metadata-only check accepted as publication')
    records.append({'name':'metadata_only_cannot_verify_publication','expected':'FAIL_SAVED_RECORD_VALIDATION','observed':report['status']})
    report={'status':'PASS_SYNTHETIC_SAVED_RESPONSE_GUARD_CONTROLS','controls':len(records),'records':records,
            'guard_sha256':guard.sha256(ROOT/'verify_saved_record.py'),
            'manifest_sha256':guard.sha256(ROOT/'FILE_MANIFEST.json'),
            'metadata_sha256':guard.sha256(ROOT/'METADATA.json'),
            'optimized':not __debug__,
            'scope':'Synthetic boundary witnesses only; no network, remote metadata/file change, publication or source evaluation.',
            'scientific_result_recomputed':False,'publication_verified_by_these_controls':False}
    if args.output:
        args.output.open('x').write(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':report['status'],'controls':len(records),'remote_mutations':0}))

if __name__=='__main__':
    main()
