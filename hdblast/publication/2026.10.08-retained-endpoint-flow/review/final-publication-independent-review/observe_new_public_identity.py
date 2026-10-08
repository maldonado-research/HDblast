"""Independent anonymous public JSON GETs; use pinned verifier for metadata semantics."""
from pathlib import Path
import hashlib
import json
import runpy

HERE = Path(__file__).resolve().parent
GET = runpy.run_path(str(HERE / 'observe_public_preservation.py'))
require = GET['require']
out = HERE / 'new-public-observation-001'
out.mkdir()
record, record_get = GET['get'](out, '23244754', 'record')
listing, listing_get = GET['get'](out, '23244754', 'files')
source_path = HERE / 'verifier-closure/verify_publication.py'
raw = source_path.read_bytes()
require(hashlib.sha256(raw).hexdigest() == 'a6d7e0db5f8a68cb893c37e24fde2aa2e37d12dad96ec13ad3c9983570b9e9bf', 'verifier source pin')
module = {'__name__': 'captured_read_only_publication_verifier', '__file__': str(source_path)}
exec(compile(raw, str(source_path), 'exec'), module)
inventory_raw = (HERE / 'SEALED_UPLOAD_INVENTORY.json').read_bytes()
metadata_raw = (HERE / 'SEALED_METADATA_MODERN.json').read_bytes()
require(hashlib.sha256(inventory_raw).hexdigest() == '2c6d0a3d4c4e227c3a8aa5f6e48dd3ff8494ac4e2f1b9045c603cf35617ced3c', 'inventory pin')
require(hashlib.sha256(metadata_raw).hexdigest() == '43feec5ca6bd9ab8c5c6a2cd925f6f42ff13418bd90b334de7ab604367dc385c', 'metadata pin')
inventory, metadata = json.loads(inventory_raw), json.loads(metadata_raw)
module['new_identity'](record, '23244754')
module['metadata_equal'](record, metadata)
rows = inventory['inherited_files'] + inventory['new_files']
pins = {row['filename']: row for row in rows}
require(len(pins) == 31, '31 inventory names')
listed = module['file_map'](listing)
embedded = module['file_map'](record)
module['assert_pins'](listed, pins)
module['assert_pins'](embedded, pins)
identities = {}
for name, item in listed.items():
    file_id = module['immutable_file_id'](item)
    version_id = module['immutable_version_id'](item)
    require(module['immutable_file_id'](embedded[name]) == file_id, 'embedded identity differs')
    if 'version_id' in embedded[name]:
        require(module['immutable_version_id'](embedded[name]) == version_id, 'embedded version differs')
    identities[name] = {'bytes': item['size'], 'checksum': item['checksum'], 'file_id': file_id, 'version_id': version_id}
receipt = {
    'status': 'PASS_INDEPENDENT_ANONYMOUS_NEW_PUBLIC_RECORD_AND_FILE_IDENTITIES',
    'record_id': '23244754', 'concept_id': '17088132', 'owner': '1386319',
    'version_doi': '10.5281/zenodo.23244754', 'concept_doi': '10.5281/zenodo.17088132',
    'published_true_draft_false_latest_true': True,
    'full_editable_metadata_matches_external_pin_with_reviewed_normalization': True,
    'inventory_sha256': hashlib.sha256(inventory_raw).hexdigest(),
    'metadata_sha256': hashlib.sha256(metadata_raw).hexdigest(),
    'verifier_sha256': hashlib.sha256(raw).hexdigest(),
    'file_count': len(identities), 'total_file_bytes': sum(row['bytes'] for row in identities.values()),
    'file_identities': identities, 'record_GET': record_get, 'files_GET': listing_get,
    'GET_requests': 2, 'authenticated_requests': 0, 'remote_writes': 0, 'attachment_streams': 0,
    'new_attachment_content_SHA256_verification': False,
    'content_attestation_scope': 'Identity/listing observation only; complete content hashes remain bound to the separate root-owned source-authenticated publisher journal.',
    'internal_independent_review_not_external_peer_review': True,
}
GET['dump'](out / 'NEW_PUBLIC_IDENTITY_REVIEW.json', receipt)
print(json.dumps({key: receipt[key] for key in ('status', 'record_id', 'file_count', 'total_file_bytes', 'remote_writes', 'attachment_streams')}))
