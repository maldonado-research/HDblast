"""Strict Notes serialization and unchanged-guard controls; synthetic only."""
import copy
import importlib.util
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent
SPEC=importlib.util.spec_from_file_location('reviewed_html',ROOT/'verify_reviewed_html_serialization.py')
review=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(review)

def require(value,message):
    if not value:raise RuntimeError(message)

def main():
    frozen,addendum,proof,manifest,guard=review.load_verified()
    controls_spec=importlib.util.spec_from_file_location('controls',review.PACKET/'test_saved_record_controls.py')
    controls=importlib.util.module_from_spec(controls_spec);controls_spec.loader.exec_module(controls)
    all_files=manifest['inherited_files']+manifest['new_files']
    saved=controls.capture(all_files,addendum['metadata'])
    done=controls.capture(all_files,addendum['metadata'],published=True)
    public=controls.capture(all_files,addendum['metadata'],published=True,public=True)
    records=[]
    result=review.verify(saved)
    require(result['status']=='PASS_REVIEWED_HTML_SERIALIZATION_COMPLETE_SAVED_DRAFT','Exact reviewed complete witness rejected')
    require(result['literal_frozen_guard']['status']=='FAIL_SAVED_RECORD_VALIDATION' and result['literal_frozen_guard']['metadata_mismatch_fields']==['notes'],'Literal discrepancy not disclosed')
    require(result['literal_frozen_guard']['failures']==['Saved metadata differs in: notes'],'Literal failure has unexpected scope')
    require(result['publication_ready'] is True and result['publication_verified'] is False,'Draft/publication flags wrong')
    records.append('exact_serialized_draft_passes_all_unchanged_guard_fields_and_discloses_literal_notes_FAIL')
    result=review.verify(done,published=True,public_saved=public)
    require(result['status']=='PASS_REVIEWED_HTML_SERIALIZATION_PUBLISHED_RECORD' and result['publication_verified'] is True,'Two published witnesses rejected')
    require(result['literal_frozen_guard']['status']=='FAIL_SAVED_RECORD_VALIDATION','Published literal discrepancy concealed')
    records.append('exact_serialized_authenticated_and_direct_public_witnesses_pass')
    notes=addendum['metadata']['notes']
    mutations=[
        ('ASCII_apostrophe',notes.replace('&rsquo;',"'")),
        ('different_entity_same_character',notes.replace('&rsquo;','&#8217;')),
        ('literal_curly_same_character',notes.replace('&rsquo;','’')),
        ('missing_separator',notes.replace('</p>\n<p>','</p><p>',1)),
        ('extra_separator',notes.replace('</p>\n<p>','</p>\n\n<p>',1)),
        ('changed_words',notes.replace('Semantic version','Forged version',1)),
        ('changed_tag',notes.replace('<p>','<div>',1)),
        ('changed_attribute',notes.replace('<p>','<p class="other">',1)),
        ('nested_tag',notes.replace('ZIP&rsquo;s','<span>ZIP&rsquo;s</span>',1)),
        ('outer_content',notes+'unreviewed'),
        ('HTML_comment',notes+'<!--new-->'),
        ('leading_whitespace',' '+notes),
    ]
    for name,variant in mutations:
        require(variant!=notes,name+' mutation did not change string')
        try:
            review.certify_pair(frozen['metadata']['notes'],variant)
            raise RuntimeError(name+' certified')
        except review.Invalid:pass
        changed=copy.deepcopy(saved);changed['data']['metadata']['notes']=variant
        require(review.verify(changed)['status']=='FAIL_REVIEWED_HTML_SERIALIZATION_RECORD',name+' record accepted')
        records.append(name+'_rejected_by_pair_and_exact_record_guard')
    for name,mutation in (
        ('other_metadata',lambda x:x['data']['metadata'].update(title='Unreviewed')),
        ('wrong_owner',lambda x:x['data'].update(owner=1)),
        ('wrong_family',lambda x:x['data'].update(conceptrecid='22922927')),
        ('wrong_record',lambda x:x['data'].update(id=22347452)),
        ('published_draft_mismatch',lambda x:x['data'].update(submitted=True,state='done')),
        ('missing_file',lambda x:x['data']['files'].pop()),
        ('duplicate_file',lambda x:x['data']['files'].append(copy.deepcopy(x['data']['files'][0]))),
        ('wrong_file_bytes',lambda x:x['data']['files'][0].update(filesize=1)),
        ('wrong_MD5',lambda x:x['data']['files'][0].update(checksum='md5:'+'0'*32)),
        ('no_authenticated_GET',lambda x:x.update(authenticated_request=False)),
        ('redirected_GET',lambda x:x.update(redirects=[{'status':302}])),
        ('tombstone',lambda x:x['data'].update(tombstone={'reason':'removed'})),
    ):
        changed=copy.deepcopy(saved);mutation(changed)
        require(review.verify(changed)['status']=='FAIL_REVIEWED_HTML_SERIALIZATION_RECORD',name+' accepted')
        records.append(name+'_rejected_by_unchanged_guard')
    badpublic=copy.deepcopy(public);badpublic['data']['metadata']['notes']=notes.replace('&rsquo;',"'")
    require(review.verify(done,published=True,public_saved=badpublic)['status']=='FAIL_REVIEWED_HTML_SERIALIZATION_RECORD','Public ASCII quote accepted')
    records.append('public_ASCII_apostrophe_rejected')
    result={'status':'PASS_STRICT_REVIEWED_HTML_SERIALIZATION_CONTROLS','controls':len(records),'records':records,
            'verifier_sha256':review.sha((ROOT/'verify_reviewed_html_serialization.py').read_bytes()),
            'reviewed_request_sha256':review.ADDENDUM_PIN,'frozen_request_sha256':review.METADATA_PIN,
            'remote_requests':0,'remote_writes':0,'scientific_result_recomputed':False,
            'scope':'Synthetic witnesses and exact pinned HTML pair only; no live metadata/publication claim.'}
    (ROOT/'OFFLINE_HTML_SERIALIZATION_CONTROLS.json').open('x').write(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:result[k] for k in ('status','controls','remote_requests','remote_writes')}))

if __name__=='__main__':main()
