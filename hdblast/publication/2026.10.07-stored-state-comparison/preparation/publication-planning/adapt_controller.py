from pathlib import Path
import hashlib
OLD=Path('/workspace/hdblast-research-work/continuation-next-20261007/release-planning/new_edition_controller.py')
OUT=Path(__file__).resolve().parent
s=OLD.read_text();assert hashlib.sha256(s.encode()).hexdigest()=='7dade987e77ff166e5849db06cd47d019f6ab16d67bd18a22c2120564c480110'
s=s.replace('This is not the old 23114217 publisher. All writable inputs must be frozen by','Prospective 27 inherited + 2 added edition. No mutation has been run.\nThe completed 23225288 publisher and its state remain untouched.\nAll writable inputs must be frozen by')
s=s.replace('PRIOR_ID = "23114217"','PRIOR_ID = "23225288"\nINHERITED_COUNT = 27\nADDITION_COUNT = 2\nFINAL_COUNT = INHERITED_COUNT + ADDITION_COUNT\nPRIOR_SEMANTIC_VERSION = "2026.10.07-bd-prehistory-target"\nSCIENTIFIC_REPLAY_STATUS = "PASS_COMPLETE_STORED_STATE_COMPARISON_REPLAY"')
s=s.replace("PRIOR_RECEIPT = Path('/workspace/hdblast-research-work/zenodo-file-bytes-20261007-02/ALL23_CONTENT_CHECKSUMS.json')", "PRIOR_RECEIPT = HERE/'read-only-observation/PRIOR27_CONTENT_RECEIPT.json'")
s=s.replace("'f3cb8d3c5280fe90307c98988bb2bf97de252bc2bebfbc20aef296612deb96ff'", "'2664ec6a867d806e6146847fcf677d2e5c20061a400262acb8563d3ca552529c'")
s=s.replace("'6c54a28fe50097a53df7dcf5803c443284bf9a7edbf58313e641f23ac7a94169'", "'94ec06030315e99c6465c784c3dcae741ed028dd89e276335b656482e11c63b4'")
s=s.replace("'817c30d77b464876b7e302aa96cfc6a8eefc413534b7662aee77d0a92280ad2c'", "'261d81302396c1683f5c0e1185b071189cf83bbb5fde5fb4a098ca64e8fe47ca'")
s=s.replace("HERE/'CURRENT_PUBLISHED_INVENTORY.json'", "HERE/'read-only-observation/CURRENT_PUBLISHED_INVENTORY.json'")
s=s.replace("HERE/'PUBLISHED_MODERN_VENDOR.json'", "HERE/'read-only-observation/CURRENT_PRIOR.raw.json'")
s=s.replace("require(len(self.old)==23 and sum(p['bytes'] for p in self.old)==440222994\n                and self.baseline['concept_id']==17088132 and self.baseline['current_record']==23114217, 'BASELINE_IDENTITY_DIFFERS')", "require(len(self.old)==INHERITED_COUNT and type(self.baseline.get('expected_bytes')) is int\n                and self.baseline['expected_bytes']>0\n                and sum(p['bytes'] for p in self.old)==self.baseline['expected_bytes']\n                and self.baseline['concept_id']==int(PARENT_ID) and self.baseline['current_record']==int(PRIOR_ID), 'BASELINE_IDENTITY_DIFFERS')")
s=s.replace("r['status']=='PASS_COMPLETE_REMOTE_BYTES_MD5_SHA256'", "r['status']=='PASS_PRIOR_COMPLETE_BYTES_WITH_FRESH_IMMUTABLE_IDENTITY'")
old="self.state=json.loads(self.state_path.read_text()) if self.state_path.exists() else {'pending':None,'draft_id':None,'published':False}"
new="""binding={'prior_record_id':PRIOR_ID,'parent_id':PARENT_ID,'controller_sha256':sha(Path(__file__).read_bytes()),
                 'baseline_sha256':sha(canonical_bytes(self.baseline)),
                 'published_baseline_sha256':sha(canonical_bytes(self.published_baseline)),
                 'prior_receipt_sha256':sha(canonical_bytes(self.prior))}
        if self.state_path.exists():
            self.state=json.loads(self.state_path.read_text())
            require(self.state.get('continuation_binding')==binding,'JOURNAL_BASELINE_OR_CONTROLLER_BINDING_DIFFERS')
        else:
            require(not (self.dir/'JOURNAL.jsonl').exists(),'ORPHAN_JOURNAL_REQUIRES_MANUAL_REVIEW')
            self.state={'pending':None,'draft_id':None,'published':False,'continuation_binding':binding}"""
assert old in s;s=s.replace(old,new)
s=s.replace("require(len(old)==23 and len(self.add)==4, 'INVENTORY_COUNTS_DIFFER')", "require(len(old)==INHERITED_COUNT and len(self.add)==ADDITION_COUNT, 'INVENTORY_COUNTS_DIFFER')\n        require(i.get('scientific_replay_status')==SCIENTIFIC_REPLAY_STATUS,'COMPLETED_SCIENTIFIC_REPLAY_REQUIRED')")
s=s.replace("len(names)==len(set(names))==27", "len(names)==len(set(names))==FINAL_COUNT")
s=s.replace("i.get('expected_total_files_if_four_additions_finalized')", "i.get('expected_total_files_final')")
s=s.replace("i['expected_total_files_if_four_additions_finalized']==27", "i['expected_total_files_final']==FINAL_COUNT")
s=s.replace("(None,'','2026.10.03-verified-source-operator')", "(None,'',PRIOR_SEMANTIC_VERSION)")
s=s.replace("(None,'','2026.10.03-verified-source-operator',self.metadata['metadata']['version'])", "(None,'',PRIOR_SEMANTIC_VERSION,self.metadata['metadata']['version'])")
needle="require(not any(k in self.metadata for k in ('id','parent','pids','files','versions','links')), 'SERVER_ASSIGNED_METADATA_FIELDS')"
replacement=needle+"""
        # Preserve all inherited metadata except explicitly reviewed update fields.
        prior=self.published_baseline
        for key in ('custom_fields','access'):
            expected=prior.get(key,{})
            if key=='access':expected={k:v for k,v in expected.items() if k!='status'}
            require(value_equal(self.metadata[key],expected,(key,)),'INHERITED_'+key.upper()+'_DIFFERS')
        old_meta=prior['metadata'];new_meta=self.metadata['metadata']
        require(set(new_meta)==set(old_meta),'INHERITED_METADATA_FIELD_SET_DIFFERS')
        editable={'title','description','additional_descriptions','publication_date','version','related_identifiers'}
        for key in set(old_meta)-editable:
            require(value_equal(new_meta[key],old_meta[key],('metadata',key)),'INHERITED_METADATA_'+key.upper()+'_DIFFERS')
        before=old_meta.get('additional_descriptions',[]);after=new_meta.get('additional_descriptions',[])
        require(isinstance(before,list) and isinstance(after,list) and len(after)==len(before)>0,'NOTES_SHAPE_DIFFERS')
        require(value_equal(after[1:],before[1:],('metadata','additional_descriptions')),'HISTORICAL_NOTES_DIFFERS')
        require(set(after[0])==set(before[0]) and all(value_equal(after[0][k],before[0][k],('metadata','additional_descriptions',0,k)) for k in set(before[0])-{'description'}),'CURRENT_NOTES_TYPE_OR_LANGUAGE_DIFFERS')
        previous=old_meta.get('related_identifiers',[]);related=new_meta.get('related_identifiers',[])
        require(isinstance(previous,list) and isinstance(related,list) and len(related)>=len(previous),'INHERITED_REFERENCES_REMOVED')
        require(value_equal(related[:len(previous)],previous,('metadata','related_identifiers')),'INHERITED_REFERENCES_DIFFERS')
        allowed_append={'identifier':'10.5281/zenodo.'+PRIOR_ID,'relation_type':{'id':'references'},'resource_type':{'id':'software'},'scheme':'doi'}
        require(len(related)==len(previous) or (len(related)==len(previous)+1 and value_equal(related[-1],allowed_append,('metadata','related_identifiers',len(previous)))),'UNREVIEWED_REFERENCE_ADDITION')
"""
assert needle in s;s=s.replace(needle,replacement)
# Public and draft family ownership are separate from the exhausted owner-list guard.
needle="metadata_equal(o,{k:self.published_baseline[k] for k in ('metadata','custom_fields','access')})"
s=s.replace(needle, "require(str(o.get('parent',{}).get('access',{}).get('owned_by',{}).get('user'))==str(EXPECTED_OWNER),'CURRENT_PRIOR_OWNER_DIFFERS')\n        "+needle)
needle="require(o['links']['self']==ORIGIN+'/api/records/'+str(id_)+'/draft','DRAFT_SELF_LINK_DIFFERS')"
s=s.replace(needle,"require(str(o.get('parent',{}).get('access',{}).get('owned_by',{}).get('user'))==str(EXPECTED_OWNER),'DRAFT_OWNER_DIFFERS')\n        "+needle)
start=s.index("        content=[]\n",s.index('    def verify('));end=s.index("        self.capture('VERIFIED_PUBLIC'",start)
s=s[:start]+"""        content=[]
        for p in self.old+self.add:
            f=files[p['filename']]
            prefix=ORIGIN+'/api/records/'+id_+('/files/' if published else '/draft/files/')+parse.quote(p['filename'],safe='')+'/content'
            require(f['links']['content']==prefix,'CONTENT_LINK_DIFFERS_FROM_VERIFIED_RECORD_AND_FILENAME')
            self.stream_check(f['links']['content'],p,public=published)
            content.append({'filename':p['filename'],'sha256':p['sha256'],'basis':'FRESH_COMPLETE_CONTENT_STREAM',
                            'remote_file_id':f.get('file_id',f.get('id')),'remote_version_id':f.get('version_id')})
        require(len(content)==FINAL_COUNT and all(x['basis']=='FRESH_COMPLETE_CONTENT_STREAM' for x in content),'ALL_EDITION_FILES_REQUIRE_FRESH_STREAMS')
        receipt={'status':'PASS_COMPLETE_'+str(FINAL_COUNT)+'_FILE_EDITION_READBACK','utc':utc(),'draft_id':id_,'parent_id':PARENT_ID,'inventory_sha256':self.inventory_sha,'metadata_sha256':self.metadata_sha,'published':published,'file_count':FINAL_COUNT,'total_bytes':self.inventory['total_bytes_final'],'content':content,'etag':next((v for k,v in h.items() if k.lower()=='etag'),None)}
"""+s[end:]
s=s.replace("        require(not self.state.get('published'),'EDITION_ALREADY_PUBLISHED_NO_REPEAT_POST')\n        require(not self.state.get('published'),'EDITION_ALREADY_PUBLISHED_NO_REPEAT_POST')", "        require(not self.state.get('published'),'EDITION_ALREADY_PUBLISHED_NO_REPEAT_POST')")
revision=OUT/'controller-revisions'/hashlib.sha256(s.encode()).hexdigest();revision.mkdir(parents=True,exist_ok=True)
(revision/'new_edition_controller.py').write_text(s)
print('candidate_sha256',hashlib.sha256(s.encode()).hexdigest())
