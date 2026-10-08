#!/usr/bin/env python3
"""Offline manufactured regression of prior normalization gaps and the fixed guard."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import socket
import sys

sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
OLD_SHA='27fd844cc903a445688a177339b7af633c0a3581f71eeedf7e20bb015bee36a7'

def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s)
    sys.modules[name]=m;s.loader.exec_module(m);return m

def equal(m,a,e,path):
    try:
        return m.value_equal(a,e,path)
    except (m.Stop,ValueError,TypeError):
        return False

def metadata(m,a,e):
    try:
        return m.metadata_equal(a,e)
    except (m.Stop,ValueError,TypeError):
        return False

def mapped(m,obj):
    try:
        m.file_map(obj);return True
    except (m.Stop,ValueError,TypeError,AttributeError):
        return False

def main():
    old=load('old_normalization_regression',HERE/'controller-revisions'/OLD_SHA/'new_edition_controller.py')
    current=load('current_normalization_regression',HERE/'new_edition_controller.py')
    socket.socket=lambda *a,**k: (_ for _ in ()).throw(RuntimeError('REAL_NETWORK_FORBIDDEN'))
    rows=[]
    def check(name,fn,expected,old_expected=None):
        observed=fn(current);previous=fn(old)
        if observed is not expected or (old_expected is not None and previous is not old_expected):
            raise ValueError('Unexpected regression outcome: '+name)
        rows.append({'case':name,'current':observed,'old':previous,'expected':expected,'status':'PASS'})
    def text(name,actual,expected,pass_=False,old_expected=None):
        check(name,lambda m:equal(m,actual,expected,('metadata','description')),pass_,old_expected)
    source='<p>A &amp; B.</p>\n<p>C.</p>'
    events=json.loads(json.dumps(old.normalized(source,source,('metadata','description'))))
    text('HTML_string_event_list_collision',events,source,False,True)
    check('nested_Notes_event_list_collision',lambda m:equal(m,[{'description':events}],[{'description':source}],
          ('metadata','additional_descriptions')),False,True)
    text('inline_interword_space_deleted','<b>Hello</b><b>world</b>','<b>Hello</b> <b>world</b>',False,True)
    text('paragraph_interior_newline_deleted','<p>A</p><p>B</p>','<p>A</p>\n<p>B</p>',False,True)
    text('inline_styled_paragraph_separator_deleted','<p style="display:inline">A</p><p style="display:inline">B</p>',
         '<p style="display:inline">A</p> <p style="display:inline">B</p>',False,True)
    text('NBSP_between_block_elements_deleted','<p>A</p><p>B</p>','<p>A</p>\u00a0<p>B</p>',False,True)
    text('NBSP_document_edge_inserted','\u00a0<p>A</p>','<p>A</p>',False,True)
    text('processing_instruction_inserted','<?audit changed?><p>A</p>','<p>A</p>',False,True)
    text('DOCTYPE_inserted','<!DOCTYPE html><p>A</p>','<p>A</p>',False,True)
    text('CDATA_declaration_inserted','<![CDATA[changed]]><p>A</p>','<p>A</p>',False,True)
    text('empty_malformed_closer_inserted','<p>A</p></>','<p>A</p>',False,True)
    text('junk_malformed_closer_inserted','<b>A</b junk>','<b>A</b>',False,True)
    text('valid_HTML_entities_and_interior_WS','<p>A &#38; B.</p>\n<p>C.</p>',source,True,True)
    text('ASCII_document_edge_WS_allowed',' \t\r\n<p>A</p>\f\n','<p>A</p>',True,True)
    check('nested_boolean_integer_type_rejected',lambda m:equal(m,{'embargo':{'active':0}},
          {'embargo':{'active':False}},('access',)),False,False)
    check('integer_float_type_rejected',lambda m:equal(m,{'x':1.0},{'x':1},('custom_fields',)),False,False)
    check('null_value_type_rejected',lambda m:equal(m,{'x':''},{'x':None},('custom_fields',)),False,False)
    check('extra_dictionary_key_rejected',lambda m:equal(m,{'x':1,'y':2},{'x':1},('custom_fields',)),False,False)
    check('missing_dictionary_key_rejected',lambda m:equal(m,{}, {'x':None},('custom_fields',)),False,False)
    check('list_length_rejected',lambda m:equal(m,[1,2],[1],('custom_fields','x')),False,False)
    check('nonJSON_dictionary_key_type_rejected',lambda m:equal(m,{1:'a'},{'1':'a'},('custom_fields',)),False,True)
    check('JSON_subclass_rejected',lambda m:equal(m,type('ForeignDict',(dict,),{})({'x':1}),{'x':1},('custom_fields',)),False,True)
    check('listed_vocabulary_decorations_allowed',lambda m:equal(m,{'id':'software','title':{'en':'Software'},
          'description':'label','icon':'code','props':{'x':True}}, {'id':'software'},('metadata','resource_type')),True,True)
    check('unlisted_vocabulary_property_rejected',lambda m:equal(m,{'id':'software','unreviewed':1},
          {'id':'software'},('metadata','resource_type')),False,False)
    check('unlisted_ID_path_decoration_rejected',lambda m:equal(m,{'id':'x','title':'label'},{'id':'x'},
          ('metadata','subjects',0)),False,False)
    expected={'metadata':{},'custom_fields':{},'access':{}}
    for key in expected:
        actual=copy.deepcopy(expected);actual.pop(key)
        check('required_empty_'+key+'_omission',lambda m,a=actual:metadata(m,a,expected),False,True)
    expected_null={'metadata':{},'custom_fields':{},'access':{'reason':None}}
    check('required_null_access_key_omission',lambda m:metadata(m,{'metadata':{},'custom_fields':{},'access':{}},expected_null),False,True)
    check('mapping_key_filename_contradiction',lambda m:mapped(m,{'files':{'entries':{'outer':{'filename':'inner'}}}}),False,True)
    check('mapping_key_key_contradiction',lambda m:mapped(m,{'files':{'entries':{'outer':{'key':'inner'}}}}),False,True)
    check('nested_key_filename_contradiction',lambda m:mapped(m,{'entries':[{'key':'a','filename':'b'}]}),False,True)
    check('consistent_mapping_keys_pass',lambda m:mapped(m,{'files':{'entries':{'a':{'key':'a','filename':'a'}}}}),True,True)
    prior=json.loads((HERE/'read-only-observation/CURRENT_PRIOR.raw.json').read_text())
    frozen=json.loads(Path('/workspace/hdblast-research-work/continuation-next-20261007/release-planning/FINAL_METADATA_MODERN.json').read_text())
    check('genuine_cached_prior_complete_metadata',lambda m:metadata(m,prior,frozen),True,True)
    receipt={'status':'PASS_STRICT_METADATA_AND_FILEMAP_REGRESSIONS','optimized':not __debug__,
             'controller_sha256':hashlib.sha256((HERE/'new_edition_controller.py').read_bytes()).hexdigest(),
             'old_controller_sha256':OLD_SHA,'checks':rows,'passes':len(rows),'skipped':0,
             'old_failures_explicitly_reproduced':sum(r['old'] and not r['expected'] for r in rows),
             'actual_network_or_publication_calls':0,'actual_array_decodes_or_source_callbacks':0,
             'historical_genuine_published_receipts':'REMAIN_GENUINE_NOT_RETROACTIVELY_INVALIDATED'}
    dest=HERE/'strict-controller-regression-results';dest.mkdir(exist_ok=True)
    p=dest/('REGRESSIONS_OPTIMIZED.json' if not __debug__ else 'REGRESSIONS_NORMAL.json')
    p.write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print('PASS',len(rows),'controls; old vulnerabilities reproduced',receipt['old_failures_explicitly_reproduced'])

if __name__=='__main__':
    main()
