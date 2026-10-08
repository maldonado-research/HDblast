#!/usr/bin/env python3
"""Manufactured saved-witness tests only; no network or numeric source work."""
from pathlib import Path
import argparse
import copy
import hashlib
import importlib.util
import json
import socket
import sys
import tempfile

sys.dont_write_bytecode = True
BASE = Path('/workspace/hdblast-research-work/continuation-next-20261007')
HERE = Path(__file__).absolute().parent
spec = importlib.util.spec_from_file_location('offline_witness_verifier', HERE/'verify_publication.py')
v = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v)

# A network attempt from any tested path fails even before a socket is opened.
def no_network(*args, **kwargs):
    raise RuntimeError('Network forbidden in manufactured witness controls')
socket.socket = no_network
socket.create_connection = no_network

INVENTORY_RAW = (BASE/'release-planning/FINAL_UPLOAD_INVENTORY.json').read_bytes()
METADATA_RAW = (BASE/'release-planning/FINAL_METADATA_MODERN.json').read_bytes()
PRIOR_RAW = Path('/workspace/hdblast-research-work/zenodo-file-bytes-20261007-02/ALL23_CONTENT_CHECKSUMS.json').read_bytes()
INV_SHA = v.sha(INVENTORY_RAW)
META_SHA = v.sha(METADATA_RAW)
v.require(INV_SHA == 'a78c80160a48e0171e8bee61e0455aa9b283b8c20e4c17e6e461c5c5fdb45153', "Manufactured positive check failed")
v.require(META_SHA == '72415e84db788b8dec0c421012c1d264fc669e6f8e5d08d93e80f0075e3485ba', "Manufactured positive check failed")
v.require(v.sha(PRIOR_RAW) == v.PRIOR_STREAM_PIN, "Manufactured positive check failed")
ID = '90000001'  # Explicit fabricated identity; never a publication assertion.
PREFIX = 'https://zenodo.org/api/records/' + ID

def fixture():
    inventory = v.load(INVENTORY_RAW)
    metadata = v.load(METADATA_RAW)
    prior = v.load(PRIOR_RAW)
    pri = {row['filename']: row for row in prior['rows']}
    entries = []
    all_pins = inventory['inherited_files'] + inventory['new_files']
    for i, pin in enumerate(all_pins):
        name = pin['filename']
        entries.append({'key': name, 'size': pin['bytes'], 'checksum': 'md5:'+pin['md5'],
                        'file_id': pri[name]['remote_file_id'] if name in pri else 'fabricated-file-'+str(i),
                        'version_id': pri[name]['remote_version_id'] if name in pri else 'fabricated-version-'+str(i),
                        'status': 'completed',
                        'links': {'content': PREFIX+'/files/'+v.quote(name, safe='')+'/content'}})
    record = dict(copy.deepcopy(metadata), id=ID, parent={'id': v.PARENT},
                  is_published=True, is_draft=False,
                  pids={'doi': {'identifier': '10.5281/zenodo.'+ID}},
                  files={'entries': {f['key']: copy.deepcopy(f) for f in entries}},
                  links={'self': PREFIX, 'files': PREFIX+'/files'})
    files = {'entries': copy.deepcopy(entries)}
    verified = {'status':'PASS_COMPLETE_27_FILE_EDITION_READBACK', 'published':True,
                'draft_id':ID, 'parent_id':v.PARENT, 'inventory_sha256':INV_SHA,
                'metadata_sha256':META_SHA, 'file_count':27,
                'total_bytes':inventory['total_bytes_final'], 'content': []}
    journal = [{'kind':'GET','url':PREFIX+'/files','status':200,'authenticated_request':False}]
    for pin in all_pins:
        name = pin['filename']
        basis = ('SAME_IMMUTABLE_FILE_AND_CONTENT_VERSION_AS_PRIOR_FULL_STREAM' if name in pri
                 else 'FRESH_COMPLETE_CONTENT_STREAM')
        verified['content'].append({'filename':name,'sha256':pin['sha256'],'basis':basis})
        if name not in pri:
            journal.append({'kind':'FULL_CONTENT_STREAM_VERIFIED','filename':name,
                            'observed':{k:pin[k] for k in ('bytes','md5','sha256')}})
    old = json.loads((BASE/'release-planning/PUBLISHED_MODERN_VENDOR.json').read_text())
    companion = json.loads((BASE/'publication-prepare/zenodo_companion.json').read_text())
    fabricated_earlier = {'id':22347452, 'conceptrecid':v.PARENT, 'submitted':True,'state':'done',
                         'metadata':{'doi':'10.5281/zenodo.22347452','title':'Fabricated earlier-record fixture'},
                         'files':[{'filename':'fabricated-previous.bin','filesize':17,'checksum':'md5:'+'1'*32}]}
    preserved = {}
    for rid, source in [('23114217',old),('22347452',fabricated_earlier),('23111008',companion)]:
        saved_files = ({'entries': list(source['files']['entries'].values())}
                       if type(source['files']) is dict else {'files': source['files']})
        preserved[rid] = {'record':copy.deepcopy(source),'baseline':copy.deepcopy(source),
                          'files':copy.deepcopy(saved_files)}
    return [inventory,metadata,record,files,copy.deepcopy(record),verified,prior,journal,preserved,
            ID,INV_SHA,META_SHA]

def checked(args):
    result=v.verify_values(*args)
    v.require(result['status']=='PASS_OFFLINE_27_FILE_PUBLICATION_WITNESSES', "Manufactured positive check failed")
    v.require(result['file_count']==27 and result['total_bytes']==448387915, "Manufactured positive check failed")
    v.require(result['historical_original_literal_Notes_guard']=='FAIL_PRESERVED', "Manufactured positive check failed")
    v.require(result['network_calls']==result['numeric_imports']==result['source_callbacks']==0, "Manufactured positive check failed")
    return result

def set_all(args, field, value):
    for index in (2,4): args[index][field]=copy.deepcopy(value)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',required=True);ns=parser.parse_args()
    positives=[];rejections=[]
    initial=fixture();r=checked(initial);positives.append('23 immutable prior streams + four fabricated fresh streams')
    v.require(r['fresh_complete_public_streams']==4 and r['prior_same_immutable_stream_reuses']==23, "Manufactured positive check failed")
    allfresh=fixture()
    for row in allfresh[5]['content']:
        row['basis']='FRESH_COMPLETE_CONTENT_STREAM'
    allfresh[7]=allfresh[7][:1]+[{'kind':'FULL_CONTENT_STREAM_VERIFIED','filename':pin['filename'],
        'observed':{k:pin[k] for k in ('bytes','md5','sha256')}}
        for pin in allfresh[0]['inherited_files']+allfresh[0]['new_files']]
    v.require(checked(allfresh)['fresh_complete_public_streams']==27, "Manufactured positive check failed")
    positives.append('all27 fabricated fresh stream witnesses')
    decorated=fixture()
    for index in (2,4):
        decorated[index]['metadata']['resource_type']['title']={'en':'Software'}
        decorated[index]['metadata']['creators'][0]['role']['props']={'server decoration':True}
        decorated[index]['metadata']['description']=decorated[index]['metadata']['description'].replace('&amp;','&#38;')
        decorated[index]['access']['status']='open'
    checked(decorated);positives.append('whitelisted vocabulary decoration/HTML character-reference normalization')

    def reject(name, mutate):
        args=fixture();mutate(args)
        try: checked(args)
        except (ValueError,TypeError,KeyError):rejections.append(name)
        else:raise AssertionError('Accepted mutant: '+name)
    reject('wrong new ID',lambda a:a.__setitem__(9,'90000002'))
    reject('reuse previous identity',lambda a:a.__setitem__(9,'23114217'))
    reject('wrong family',lambda a:a[2]['parent'].__setitem__('id','22922927'))
    reject('public not published',lambda a:a[2].__setitem__('is_published',False))
    reject('public marked draft',lambda a:a[2].__setitem__('is_draft',True))
    reject('draft marker omitted',lambda a:a[2].pop('is_draft'))
    reject('wrong DOI',lambda a:a[2]['pids']['doi'].__setitem__('identifier','10.5281/zenodo.90000002'))
    reject('latest stale',lambda a:a[4].__setitem__('id','23114217'))
    reject('public files link wrong',lambda a:a[2]['links'].__setitem__('files',PREFIX+'/draft/files'))
    reject('missing published file',lambda a:a[3]['entries'].pop())
    reject('duplicate published file',lambda a:a[3]['entries'].append(copy.deepcopy(a[3]['entries'][0])))
    reject('wrong published size',lambda a:a[3]['entries'][0].__setitem__('size',1))
    reject('size changes type',lambda a:a[3]['entries'][0].__setitem__('size',float(a[3]['entries'][0]['size'])))
    reject('wrong published MD5',lambda a:a[3]['entries'][0].__setitem__('checksum','md5:'+'0'*32))
    reject('unfinished published file',lambda a:a[3]['entries'][0].__setitem__('status','pending'))
    reject('embedded file identity differs',lambda a:next(iter(a[2]['files']['entries'].values())).__setitem__('version_id','other'))
    reject('missing SHA basis',lambda a:a[5]['content'].pop())
    reject('duplicate SHA basis',lambda a:a[5]['content'].append(copy.deepcopy(a[5]['content'][0])))
    reject('wrong SHA basis hash',lambda a:a[5]['content'][0].__setitem__('sha256','0'*64))
    reject('unknown SHA basis',lambda a:a[5]['content'][0].__setitem__('basis','SELF_ATTESTED'))
    reject('new file uses inherited basis',lambda a:a[5]['content'][-1].__setitem__('basis','SAME_IMMUTABLE_FILE_AND_CONTENT_VERSION_AS_PRIOR_FULL_STREAM'))
    reject('inherited immutable identity differs',lambda a:a[3]['entries'][0].__setitem__('file_id','other'))
    reject('inherited prior stream hash differs',lambda a:a[6]['rows'][0]['observed'].__setitem__('sha256','0'*64))
    reject('missing prior stream row',lambda a:a[6]['rows'].pop())
    reject('missing terminal public GET',lambda a:a[7].pop(0))
    reject('terminal GET authenticated',lambda a:a[7][0].__setitem__('authenticated_request',True))
    reject('missing fresh stream',lambda a:a[7].pop())
    reject('fresh stream truncated',lambda a:a[7][-1]['observed'].__setitem__('bytes',1))
    reject('fresh stream MD5 altered',lambda a:a[7][-1]['observed'].__setitem__('md5','0'*32))
    reject('fresh stream SHA altered',lambda a:a[7][-1]['observed'].__setitem__('sha256','0'*64))
    reject('fresh stream before terminal GET',lambda a:a[7].append(copy.deepcopy(a[7][0])))
    reject('duplicate fresh stream',lambda a:a[7].append(copy.deepcopy(a[7][-1])))
    reject('fresh content URL wrong',lambda a:a[3]['entries'][-1]['links'].__setitem__('content',PREFIX+'/draft/files/other/content'))
    reject('verification is draft',lambda a:a[5].__setitem__('published',False))
    reject('verification wrong inventory pin',lambda a:a[5].__setitem__('inventory_sha256','0'*64))
    reject('verification wrong metadata pin',lambda a:a[5].__setitem__('metadata_sha256','0'*64))
    reject('verification wrong total',lambda a:a[5].__setitem__('total_bytes',448387916))
    reject('metadata field omitted',lambda a:a[2]['metadata'].pop('languages'))
    reject('metadata unexpected field',lambda a:a[2]['metadata'].__setitem__('unregistered',True))
    reject('ORCID altered',lambda a:a[2]['metadata']['creators'][0]['person_or_org']['identifiers'][0].__setitem__('identifier','0000-0000-0000-0000'))
    reject('references omitted',lambda a:a[2]['metadata'].__setitem__('references',[]))
    reject('description meaning altered',lambda a:a[2]['metadata'].__setitem__('description',a[2]['metadata']['description']+'NOVELTY PROVED'))
    reject('metadata type altered',lambda a:a[2]['metadata'].__setitem__('version',True))
    reject('vocabulary ID altered',lambda a:a[2]['metadata']['resource_type'].__setitem__('id','dataset'))
    reject('access altered',lambda a:a[2]['access'].__setitem__('files','restricted'))
    reject('custom fields omitted',lambda a:a[2].__setitem__('custom_fields',{}))
    reject('missing preserved record',lambda a:a[8].pop('22347452'))
    reject('preserved record metadata altered',lambda a:a[8]['22347452']['record']['metadata'].__setitem__('title','changed'))
    reject('preserved record file omitted',lambda a:a[8]['22347452']['files']['files'].pop())
    reject('preserved main embedded files omitted',lambda a:a[8]['23114217']['record'].__setitem__('files',[]))
    reject('preserved earlier embedded files omitted',lambda a:a[8]['22347452']['record'].__setitem__('files',[]))
    reject('preserved companion embedded files omitted',lambda a:a[8]['23111008']['record'].__setitem__('files',[]))
    reject('preserved earlier embedded MD5 changed',lambda a:a[8]['22347452']['record']['files'][0].__setitem__('checksum','md5:'+'0'*32))
    reject('preserved main embedded identity changed',lambda a:next(iter(a[8]['23114217']['record']['files']['entries'].values())).__setitem__('version_id','other'))
    reject('preserved record size altered',lambda a:a[8]['22347452']['files']['files'][0].__setitem__('filesize',18))
    reject('companion family changed',lambda a:a[8]['23111008']['record'].__setitem__('conceptrecid','17088132'))
    reject('inventory unsealed',lambda a:a[0].__setitem__('sealed_upload_manifest',False))
    reject('inventory wrong total',lambda a:a[0].__setitem__('total_bytes_final',1))
    reject('inventory pin type altered',lambda a:a[0]['new_files'][0].__setitem__('bytes',True))

    with tempfile.TemporaryDirectory(prefix='offline-publication-controls-') as tmp:
        root=Path(tmp)
        args=fixture();names=['inventory','metadata','new_public_record','new_public_files','latest_public_record',
                             'verified_public','prior_full_stream']
        vals=dict(zip(names,args[:7])); roles={}; files={}
        for role,value in vals.items():
            name=role+'.json';raw=PRIOR_RAW if role=='prior_full_stream' else canonical(value)
            (root/name).write_bytes(raw);roles[role]=name;files[name]={'bytes':len(raw),'sha256':v.sha(raw)}
        raw=b''.join(canonical(event)+b'\n' for event in args[7]);name='journal.jsonl'
        (root/name).write_bytes(raw);roles['controller_journal']=name;files[name]={'bytes':len(raw),'sha256':v.sha(raw)}
        groups={}
        for rid,group in args[8].items():
            groups[rid]={}
            for role,value in group.items():
                name=rid+'-'+role+'.json';raw=canonical(value);(root/name).write_bytes(raw)
                groups[rid][role]=name;files[name]={'bytes':len(raw),'sha256':v.sha(raw)}
        manifest={'schema_version':1,'files':files,'witnesses':roles,'preserved_records':groups,
                  'historical_original_literal_Notes_guard':'FAIL_PRESERVED'}
        mraw=canonical(manifest);mpath=root/'WITNESS_MANIFEST.json';mpath.write_bytes(mraw)
        # Directory fixture retains actual externally pinned inventory/metadata bytes.
        for role,raw in [('inventory',INVENTORY_RAW),('metadata',METADATA_RAW)]:
            (root/roles[role]).write_bytes(raw);files[roles[role]]={'bytes':len(raw),'sha256':v.sha(raw)}
        mraw=canonical(manifest);mpath.write_bytes(mraw);mpin=v.sha(mraw)
        def directory():return v.verify_directory(root,mpath,mpin,ID,INV_SHA,META_SHA)
        v.require(directory()['status']=='PASS_OFFLINE_27_FILE_PUBLICATION_WITNESSES', "Manufactured positive check failed")
        positives.append('complete externally pinned manufactured witness directory')
        def directory_reject(name,operation):
            try: operation()
            except (ValueError,TypeError,KeyError):rejections.append(name)
            else:raise AssertionError('Accepted directory mutant: '+name)
        directory_reject('external witness manifest pin wrong',lambda:v.verify_directory(root,mpath,'0'*64,ID,INV_SHA,META_SHA))
        directory_reject('external inventory pin wrong',lambda:v.verify_directory(root,mpath,mpin,ID,'0'*64,META_SHA))
        directory_reject('external metadata pin wrong',lambda:v.verify_directory(root,mpath,mpin,ID,INV_SHA,'0'*64))
        target=root/roles['verified_public'];original=target.read_bytes();target.write_bytes(original+b' ')
        directory_reject('saved witness bytes changed',directory);target.write_bytes(original)
        real=root/'real.json';real.write_bytes(original);target.unlink();target.symlink_to(real)
        directory_reject('symlink witness',directory);target.unlink();target.write_bytes(original)
        edited=copy.deepcopy(manifest);edited['historical_original_literal_Notes_guard']='PASS'
        eraw=canonical(edited);mpath.write_bytes(eraw)
        directory_reject('original Notes failure promoted',lambda:v.verify_directory(root,mpath,v.sha(eraw),ID,INV_SHA,META_SHA))
        mpath.write_bytes(mraw)
        directory_reject('duplicate JSON keys',lambda:v.load('{"x":1,"x":2}'))
        directory_reject('nonfinite JSON value',lambda:v.load('{"x":NaN}'))
        directory_reject('unsafe relative path',lambda:v.safe_file(root,'../escape.json'))
    receipt={'status':'PASS_MANUFACTURED_OFFLINE_PUBLICATION_CONTROLS','positive_checks':len(positives),
             'rejected_mutants':len(rejections),'positives':positives,'rejections':rejections,
             'fixture_scope':'FABRICATED_NEW_RECORD_WITNESSES_ONLY_WITH_PINNED_OLD_INVENTORY_AND_STREAM_RECEIPT',
             'fabricated_record_id':ID,'verifier_sha256':v.sha((HERE/'verify_publication.py').read_bytes()),
             'test_sha256':v.sha(Path(__file__).read_bytes()),'network_calls':0,'numeric_imports':0,
             'source_callbacks':0,'array_decodes':0,'remote_writes':0}
    Path(ns.output).write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps({k:receipt[k] for k in ('status','positive_checks','rejected_mutants','network_calls')}))

def canonical(value):return json.dumps(value,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()
if __name__=='__main__':main()
