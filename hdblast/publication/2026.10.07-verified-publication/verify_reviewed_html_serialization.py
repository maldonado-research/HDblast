"""Certify ONE exact Notes HTML variant with the unchanged saved-record guard.

The frozen request and literal guard remain unchanged. Four declared outside-
paragraph LF nodes and one named character entity are the entire allowed pair.
No general HTML normalization, remote request, metadata write or publish exists.
"""
import argparse
import hashlib
from html.parser import HTMLParser
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PACKET = Path('/workspace/HDblast/hdblast/publication/2026.10.03-verified-source-operator')
MANIFEST_PIN = '00f0ed73c27e846603680fdc52667392935e59bbe18ac942f30e8e66248937bb'
METADATA_PIN = 'd96514fa7d8723706f2bd7fd3dcf8a5928ce58505c09bdd911b98f68c3c772af'
GUARD_PIN = 'a05af063d490d208f23292d11da68bb023fd69c639cff7b0f31a91666b8d47db'
ADDENDUM_PIN = '66fea4d6a5666b441963e63372977221917dd3239335b660887fea8869821cb3'
FROZEN_NOTES_PIN = '74cff9d8938341e3d6643102c01a41bd73c6f364310b45c60a392a7ed7f95688'
SERIALIZED_NOTES_PIN = '77668fcd9ef479c303f67e30dfc7f7c770e9c85d968217a2e9ff981fcc8e4c98'

class Invalid(RuntimeError):
    pass

def require(value, message):
    if not value:
        raise Invalid(message)

def sha(raw):
    return hashlib.sha256(raw).hexdigest()

def read_pinned(path, expected):
    raw = path.read_bytes()
    require(sha(raw) == expected, 'Pinned verification input differs.')
    return json.loads(raw)

class Paragraphs(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.inside = False
        self.paragraphs = []
        self.elements = []
        self.outside = []
    def handle_starttag(self, tag, attrs):
        require(tag == 'p' and attrs == [] and self.inside is False, 'Unexpected tag, attribute or nesting.')
        self.inside = True
        self.paragraphs.append('')
        self.elements.append(('start', tag, attrs))
    def handle_endtag(self, tag):
        require(tag == 'p' and self.inside is True, 'Unexpected closing tag.')
        self.inside = False
        self.elements.append(('end', tag))
    def handle_data(self, value):
        if self.inside:
            self.paragraphs[-1] += value
        else:
            self.outside.append((len(self.paragraphs), value))
    def handle_startendtag(self, tag, attrs):
        raise Invalid('Self-closing tag rejected.')
    def handle_comment(self, text):
        raise Invalid('Comment rejected.')
    def handle_decl(self, text):
        raise Invalid('Declaration rejected.')
    def handle_pi(self, text):
        raise Invalid('Processing instruction rejected.')
    def unknown_decl(self, text):
        raise Invalid('Unknown declaration rejected.')

def signature(notes):
    parser = Paragraphs()
    parser.feed(notes)
    parser.close()
    require(parser.inside is False and len(parser.paragraphs) == 5 and len(parser.elements) == 10,
            'Expected five complete paragraph elements.')
    return parser

def certify_pair(frozen_notes, serialized_notes):
    require(sha(frozen_notes.encode()) == FROZEN_NOTES_PIN, 'Frozen Notes content pin differs.')
    require(sha(serialized_notes.encode()) == SERIALIZED_NOTES_PIN, 'Only the one exact reviewed Notes serialization is allowed.')
    require(frozen_notes.count('ZIP’s') == 1 and frozen_notes.count('</p><p>') == 4,
            'Registered entity/separator locations differ.')
    derived = frozen_notes.replace('ZIP’s', 'ZIP&rsquo;s', 1).replace('</p><p>', '</p>\n<p>')
    require(serialized_notes == derived and serialized_notes.count('&rsquo;') == 1,
            'Changes exceed the registered one-entity/four-separator transformation.')
    frozen = signature(frozen_notes)
    serialized = signature(serialized_notes)
    require(frozen.elements == serialized.elements and frozen.paragraphs == serialized.paragraphs,
            'Paragraph elements/attributes or decoded Unicode text differ.')
    require(frozen.outside == [] and serialized.outside == [(i, '\n') for i in range(1, 5)],
            'Outside text differs from the four declared inter-paragraph LF nodes.')
    return {'status':'PASS_EXACT_REVIEWED_HTML_SERIALIZATION_PAIR',
            'frozen_notes_sha256':FROZEN_NOTES_PIN, 'serialized_notes_sha256':SERIALIZED_NOTES_PIN,
            'paragraph_elements':5, 'paragraph_attributes':0, 'identical_decoded_Unicode_paragraph_text':True,
            'identical_element_sequence':True, 'one_registered_entity':'&rsquo; -> U+2019',
            'only_outside_text_difference':'Four LF text nodes between block paragraph elements at boundaries 1,2,3,4.',
            'scope':'Exact pinned pair; full raw DOM has the four disclosed whitespace nodes. Paragraph structure and Unicode text are identical.'}

def load_verified():
    frozen = read_pinned(PACKET/'METADATA.json', METADATA_PIN)
    addendum = read_pinned(ROOT/'REVIEWED_HTML_SERIALIZATION_METADATA_REQUEST.json', ADDENDUM_PIN)
    require(isinstance(addendum,dict) and set(addendum)=={'metadata'}, 'Reviewed metadata wrapper differs.')
    check = json.loads(json.dumps(addendum))
    check['metadata']['notes'] = frozen['metadata']['notes']
    require(check == frozen, 'Reviewed metadata differs beyond exact Notes serialization.')
    proof = certify_pair(frozen['metadata']['notes'], addendum['metadata']['notes'])
    manifest = read_pinned(PACKET/'FILE_MANIFEST.json', MANIFEST_PIN)
    require(sha((PACKET/'verify_saved_record.py').read_bytes()) == GUARD_PIN, 'Frozen guard changed.')
    spec = importlib.util.spec_from_file_location('unchanged_saved_record_guard',PACKET/'verify_saved_record.py')
    guard = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(guard)
    guard.validate_candidate(manifest,PACKET)
    return frozen,addendum,proof,manifest,guard

def verify(saved,published=False,public_saved=None):
    frozen,addendum,proof,manifest,guard=load_verified()
    literal=guard.verify(saved,manifest,frozen,23114217,published=published,public_saved=public_saved)
    reviewed=guard.verify(saved,manifest,addendum,23114217,published=published,public_saved=public_saved)
    wanted='PASS_VERIFIED_PUBLISHED_RECORD' if published else 'PASS_COMPLETE_SAVED_DRAFT'
    passed=reviewed['status']==wanted
    status=('PASS_REVIEWED_HTML_SERIALIZATION_PUBLISHED_RECORD' if published else
            'PASS_REVIEWED_HTML_SERIALIZATION_COMPLETE_SAVED_DRAFT') if passed else 'FAIL_REVIEWED_HTML_SERIALIZATION_RECORD'
    return {'status':status,'record_id':23114217,'concept_doi':'10.5281/zenodo.17088132',
            'frozen_metadata_sha256':METADATA_PIN,'reviewed_metadata_sha256':ADDENDUM_PIN,
            'frozen_guard_sha256':GUARD_PIN,'manifest_sha256':MANIFEST_PIN,'semantic_pair_proof':proof,
            'literal_frozen_guard':literal,'unchanged_guard_against_reviewed_expected_metadata':reviewed,
            'literal_frozen_guard_unchanged':True,'literal_failure_disclosed':literal['status'].startswith('FAIL_'),
            'publication_ready':passed and not published,'publication_verified':passed and published,
            'remote_requests':0,'remote_mutations':0,'scientific_result_recomputed':False,
            'scope':'Separately pinned reviewed HTML serialization, exact other metadata/inventory/owner/family/state; original literal result retained.'}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--saved-response',type=Path,required=True)
    parser.add_argument('--published',action='store_true')
    parser.add_argument('--public-response',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    _,_,_,_,guard=load_verified()
    saved=guard.read_json(args.saved_response)
    public=guard.read_json(args.public_response) if args.public_response else None
    report=verify(saved,published=args.published,public_saved=public)
    args.output.open('x').write(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':report['status'],'literal_frozen_guard_status':report['literal_frozen_guard']['status'],
                      'publication_ready':report['publication_ready'],'publication_verified':report['publication_verified']}))
    return 0 if report['status'].startswith('PASS_') else 1

if __name__=='__main__':
    raise SystemExit(main())
