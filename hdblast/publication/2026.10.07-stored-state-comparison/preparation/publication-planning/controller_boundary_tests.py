#!/usr/bin/env python3
"""Offline controls for critical metadata, inventory, and credential boundaries."""
import copy
import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parent
SPEC=importlib.util.spec_from_file_location('controller_boundary_subject',ROOT/'new_edition_controller.py')
C=importlib.util.module_from_spec(SPEC);sys.modules[SPEC.name]=C;SPEC.loader.exec_module(C)

class BoundaryControls(unittest.TestCase):
    def expected(self):
        return {'metadata':{'title':'Frozen title','description':'<p>First &amp; second.</p>\n<p>Third.</p>',
                'resource_type':{'id':'publication-preprint'},'subjects':[{'subject':'physics'}]},
                'custom_fields':{'example':112,'code:programmingLanguage':[{'id':'python'}]},
                'access':{'record':'public','files':'public','embargo':{'active':False,'reason':None}}}
    def rejected_metadata(self,change):
        e=self.expected();a=copy.deepcopy(e);change(a)
        with self.assertRaises(C.Stop):C.metadata_equal(a,e)
    def test_identical_complete_metadata_passes(self):
        self.assertTrue(C.metadata_equal(self.expected(),self.expected()))
    def test_html_entities_with_preserved_interior_whitespace_pass(self):
        e=self.expected();a=copy.deepcopy(e);a['metadata']['description']='<p>First &#38; second.</p>\n<p>Third.</p>'
        self.assertTrue(C.metadata_equal(a,e))
    def test_html_inside_paragraph_whitespace_change_rejected(self):
        self.rejected_metadata(lambda a:a['metadata'].__setitem__('description','<p>First  &amp; second.</p>\n<p>Third.</p>'))
    def test_html_changed_paragraph_rejected(self):
        self.rejected_metadata(lambda a:a['metadata'].__setitem__('description','<p>First &amp; second.</p>\n<p>Different.</p>'))
    def test_custom_float_substitution_rejected(self):
        self.rejected_metadata(lambda a:a['custom_fields'].__setitem__('example',112.0))
    def test_custom_boolean_substitution_rejected(self):
        e=self.expected();e['custom_fields']['example']=1;a=copy.deepcopy(e);a['custom_fields']['example']=True
        with self.assertRaises(C.Stop):C.metadata_equal(a,e)
    def test_access_boolean_integer_substitution_rejected(self):
        self.rejected_metadata(lambda a:a['access']['embargo'].__setitem__('active',0))
    def test_extra_metadata_field_rejected(self):
        self.rejected_metadata(lambda a:a['metadata'].__setitem__('new_unreviewed_field','value'))
    def test_missing_metadata_field_rejected(self):
        self.rejected_metadata(lambda a:a['metadata'].pop('title'))
    def test_missing_custom_field_rejected(self):
        self.rejected_metadata(lambda a:a['custom_fields'].pop('example'))
    def test_known_vocabulary_decoration_passes(self):
        e=self.expected();a=copy.deepcopy(e);a['metadata']['resource_type']['title']={'en':'Preprint'};a['custom_fields']['code:programmingLanguage'][0]['title']={'en':'Python'}
        self.assertTrue(C.metadata_equal(a,e))
    def test_vocabulary_id_substitution_rejected(self):
        self.rejected_metadata(lambda a:a['metadata']['resource_type'].__setitem__('id','dataset'))
    def test_unrecognized_id_decoration_not_dropped(self):
        e=self.expected();e['metadata']['subjects']=[{'id':'physics'}];a=copy.deepcopy(e);a['metadata']['subjects'][0]['title']='unreviewed'
        with self.assertRaises(C.Stop):C.metadata_equal(a,e)
    def test_unknown_vocabulary_property_rejected(self):
        self.rejected_metadata(lambda a:a['metadata']['resource_type'].__setitem__('evil_extra',True))
    def test_known_server_access_status_passes(self):
        e=self.expected();a=copy.deepcopy(e);a['access']['status']='open';self.assertTrue(C.metadata_equal(a,e))
    def test_unknown_access_field_rejected(self):
        self.rejected_metadata(lambda a:a['access'].__setitem__('unreviewed_policy','hidden'))
    def test_filename_size_md5_complete_pass(self):
        p={'filename':'a.zip','bytes':1,'md5':'a'*32};self.assertTrue(C.assert_pins({'a.zip':{'size':1,'checksum':'md5:'+'a'*32,'status':'completed'}},[p]))
    def test_float_and_boolean_file_size_rejected(self):
        p={'filename':'a.zip','bytes':1,'md5':'a'*32}
        for size in (1.0,True):
            with self.subTest(size=size),self.assertRaises(C.Stop):C.assert_pins({'a.zip':{'size':size,'checksum':'a'*32}},[p])
    def test_extra_file_rejected(self):
        p={'filename':'a.zip','bytes':1,'md5':'a'*32}
        with self.assertRaises(C.Stop):C.assert_pins({'a.zip':{'size':1,'checksum':'a'*32},'b.zip':{}},[p])
    def test_wrong_md5_rejected(self):
        p={'filename':'a.zip','bytes':1,'md5':'a'*32}
        with self.assertRaises(C.Stop):C.assert_pins({'a.zip':{'size':1,'checksum':'b'*32}},[p])
    def test_incomplete_file_rejected(self):
        p={'filename':'a.zip','bytes':1,'md5':'a'*32}
        with self.assertRaises(C.Stop):C.assert_pins({'a.zip':{'size':1,'checksum':'a'*32,'status':'pending'}},[p])
    def test_unsafe_destination_and_credential_query_rejected(self):
        for url in ('http://zenodo.org/api/records/1','https://evil.example/api/records/1','https://token@zenodo.org/api/records/1','https://zenodo.org/api/records/1?access_token=fake','https://zenodo.org/api/records/1#fragment'):
            with self.subTest(url=url),self.assertRaises(C.Stop):C.safe_url(url,query=True)
    def test_bound_owner_query_allowed(self):
        u='https://zenodo.org/api/deposit/depositions?page=1&size=100&sort=mostrecent';self.assertEqual(C.safe_url(u,query=True),u)
    def test_fake_token_and_query_redacted(self):
        fake='OFFLINE_FAKE_TOKEN_ABC';o=C.redact({'token':fake,'url':'https://zenodo.org/api/records/1?access_token='+fake,'message':'contains '+fake},fake)
        self.assertNotIn(fake,json.dumps(o));self.assertNotIn('access_token=',json.dumps(o))
    def test_embedded_error_url_query_and_malformed_url_are_redacted(self):
        value=C.redact('Provider says https://zenodo.org/api/records/1?secret=OTHER_PRIVATE_VALUE and https://[bad?secret=OTHER_PRIVATE_VALUE')
        self.assertNotIn('OTHER_PRIVATE_VALUE',value)
        self.assertIn('[REDACTED_MALFORMED_URL]',value)
    def test_error_json_structure_and_nonfinite_limits(self):
        for body in (b'['*70+b'0'+b']'*70,b'{"x":NaN}',b'{"x":1e309}'):
            with self.subTest(body=body[:20]),self.assertRaises(Exception):C.bounded_error_json(body)
    def test_canonical_quoted_etag_converts_to_matching_decimal_revision(self):
        self.assertEqual(C.revision_if_match('"7"',{'revision_id':7}),'7')
    def test_weak_negative_nondecimal_unquoted_or_noncanonical_etags_rejected(self):
        for value in ('W/"7"','"-7"','"seven"','7','"07"','"7 "','"７"',None):
            with self.subTest(value=value),self.assertRaisesRegex(C.Stop,'ETAG_NOT_CANONICAL_QUOTED_DECIMAL'):
                C.revision_if_match(value,{'revision_id':7})
    def test_etag_requires_integer_matching_record_revision(self):
        for revision in (True,7.0,8,None,-7):
            with self.subTest(revision=revision),self.assertRaisesRegex(C.Stop,'ETAG_REVISION_ID_DIFFERS'):
                C.revision_if_match('"7"',{'revision_id':revision})
    def test_transport_http400_reads_bounded_body_and_preserves_status(self):
        body=b'{"message":"Provider validation error"}'
        error=C.error.HTTPError('https://zenodo.org/api/records/1/draft',400,'error',{},io.BytesIO(body))
        class FakeOpener:
            def open(self,*args,**kwargs):raise error
        transport=C.Transport('OFFLINE_FAKE_TOKEN_ONLY');transport.opener=FakeOpener()
        result=transport.request('PUT','https://zenodo.org/api/records/1/draft',data=b'{}')
        self.assertEqual(result.status,400);self.assertEqual(result.body,body);self.assertFalse(result.body_truncated)
    def test_transport_error_oversize_is_bounded(self):
        error=C.error.HTTPError('https://zenodo.org/api/records/1/draft',400,'error',{},io.BytesIO(b'x'*(C.MAX_JSON+100)))
        class FakeOpener:
            def open(self,*args,**kwargs):raise error
        transport=C.Transport('OFFLINE_FAKE_TOKEN_ONLY');transport.opener=FakeOpener()
        result=transport.request('PUT','https://zenodo.org/api/records/1/draft',data=b'{}')
        self.assertEqual(result.status,400);self.assertEqual(len(result.body),C.MAX_JSON+1);self.assertTrue(result.body_truncated)
    def test_cli_mutation_flags_reject_before_transport(self):
        with tempfile.TemporaryDirectory() as evidence:
            for stage in ('create','prepare','publish'):
                with self.subTest(stage=stage),self.assertRaisesRegex(C.Stop,'MUTATION_STAGE_REQUIRES_EXPLICIT_MUTATE_FLAG'):
                    C.main(['--stage',stage,'--evidence-dir',evidence])
    def test_cli_publish_requires_second_explicit_flag(self):
        with tempfile.TemporaryDirectory() as evidence,self.assertRaisesRegex(C.Stop,'PUBLISH_STAGE_REQUIRES_EXPLICIT_PUBLISH_FLAG'):
            C.main(['--stage','publish','--mutate','--evidence-dir',evidence])
    def test_cli_publish_flag_cannot_hide_in_other_stage(self):
        with tempfile.TemporaryDirectory() as evidence,self.assertRaisesRegex(C.Stop,'PUBLISH_FLAG_ONLY_VALID_FOR_PUBLISH_STAGE'):
            C.main(['--stage','reconcile','--publish','--evidence-dir',evidence])
    def test_cli_406_retry_requires_explicit_create_mutation(self):
        with tempfile.TemporaryDirectory() as evidence,self.assertRaisesRegex(C.Stop,'CREATE_406_RETRY_FLAG_ONLY_FOR_EXPLICIT_CREATE_MUTATION'):
            C.main(['--stage','reconcile','--retry-create-406','--evidence-dir',evidence])
    def test_frozen_input_digest_change_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'input.json';p.write_bytes(b'{"reviewed":true}\n');pin=C.sha(p.read_bytes());p.write_bytes(b'{"reviewed":false}\n')
            with self.assertRaisesRegex(C.Stop,'FROZEN_INPUT_SHA256_DIFFERS'):C.frozen_json(p,pin)
    def test_evidence_process_lock_excludes_second_controller(self):
        with tempfile.TemporaryDirectory() as evidence:
            with C.EvidenceLock(evidence):
                with self.assertRaisesRegex(C.Stop,'EVIDENCE_DIRECTORY_ALREADY_LOCKED'):
                    with C.EvidenceLock(evidence):pass
            with C.EvidenceLock(evidence):pass

if __name__=='__main__':unittest.main(verbosity=2)
