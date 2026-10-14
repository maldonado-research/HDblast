"""Write MANIFEST.json: path, bytes and SHA-256 for every file in the capsule
(excluding the manifest itself, caches and transient logs of killed runs)."""
import hashlib, json, os
HERE = os.path.dirname(os.path.abspath(__file__))
entries = []
for root, dirs, files in os.walk(HERE):
    dirs[:] = sorted(d for d in dirs if d not in ('__pycache__', 'site_patch'))
    for f in sorted(files):
        if f in ('MANIFEST.json', '.DS_Store') or f.endswith(('.pyc', '.pkl')):
            continue
        p = os.path.join(root, f)
        b = open(p, 'rb').read()
        entries.append({'path': os.path.relpath(p, HERE), 'bytes': len(b),
                        'sha256': hashlib.sha256(b).hexdigest()})
json.dump({'capsule': 'HDBLAST_SECOND_ORDER_LOCALIZED_20261014', 'files': entries},
          open(os.path.join(HERE, 'MANIFEST.json'), 'w'), indent=2)
print(len(entries), 'files')
