"""Install a newly frozen capsule and dedicated replay CI; preserve old CI bytes."""
from hashlib import sha256
import json
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parent.parent
CAP = 'HDBLAST_UNIFORM_FIRST_VARIATION_WAVEFORM_20261015'
OLD = 'HDBLAST_SHELL_RESPONSE_ERROR_CERTIFICATE_20261015'
WT = ROOT / 'publish-HDblast'


def require(condition, reason):
    if not condition:
        raise RuntimeError(reason)


receipt = json.loads((ROOT / 'frozen-source/SOURCE_FREEZE_RECEIPT.json').read_text())
source = ROOT / 'frozen-source' / CAP
target = WT / 'research' / CAP
require(not target.exists(), 'Fresh capsule path required')
shutil.copytree(source, target)
shutil.copyfile(ROOT / 'frozen-source' / (CAP + '.zip'), WT / 'research' / (CAP + '.zip'))
old = (WT / '.github/workflows/hdblast_shell_response_certificate.yml').read_text()
new = old.replace(OLD, CAP).replace('Shell-response quantitative certificate replay',
                                    'Uniform first-variation waveform enclosure replay')
new = new.replace('hdblast_shell_response_certificate.yml', 'hdblast_uniform_first_variation_waveform.yml')
new = new.replace('shell-response-error-certificate', 'uniform-first-variation-waveform')
pins = {
    'MANIFEST_SHA256': receipt['manifest']['sha256'],
    'ARCHIVE_SHA256': receipt['archive']['sha256'],
    'RUNNER_SHA256': sha256((source / 'reproduce.py').read_bytes()).hexdigest(),
    'VERIFIER_SHA256': sha256((source / 'verify_manifest.py').read_bytes()).hexdigest(),
}
import re
for key, value in pins.items():
    new, count = re.subn(r"(" + key + r": ')[0-9a-f]{64}(')",
                         lambda match: match[1] + value + match[2], new)
    require(count == 1, 'Unique workflow pin required: ' + key)
new = new.replace('Finite identity checks do not certify an analytic response theorem.',
                  'Finite checks do not formally certify the conditional analytic PDE proof.')
workflow = WT / '.github/workflows/hdblast_uniform_first_variation_waveform.yml'
require(not workflow.exists(), 'New workflow required')
workflow.write_text(new)
names = subprocess.check_output(['git', '-C', str(WT), 'ls-tree', '-r', '--name-only',
                                 'HEAD', '.github/workflows'], text=True).splitlines()
for name in names:
    require((WT / name).read_bytes() == subprocess.check_output(
        ['git', '-C', str(WT), 'show', 'HEAD:' + name]), 'Existing workflow altered: ' + name)
for path in source.rglob('*'):
    if path.is_file():
        require(path.read_bytes() == (target / path.relative_to(source)).read_bytes(),
                'Installed source bytes differ')
result = {'status': 'PASS', 'capsule': CAP, 'source_files': receipt['source_files'],
          'archive_members': receipt['archive_members'], 'old_workflows_preserved': len(names),
          'workflow_path': workflow.relative_to(WT).as_posix(),
          'workflow_sha256': sha256(workflow.read_bytes()).hexdigest(), 'pins': pins}
(ROOT / 'publication/CI_INSTALL_PREFLIGHT.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
