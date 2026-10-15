"""Add the waveform checkpoint and prove recovery of all previous summary bytes.

No Git mutations are performed. Only nine owned summaries and two new PNG
assets are written; source and historical checkpoint files are left untouched.
"""
import hashlib
import html
import json
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parent.parent
MARK = 'uniform-first-variation-waveform-20261015'
CAPSULE = 'HDBLAST_UNIFORM_FIRST_VARIATION_WAVEFORM_20261015'
BASE = 'https://github.com/maldonado-research/HDblast/blob/main/research/' + CAPSULE + '/'
TREE = 'https://github.com/maldonado-research/HDblast/tree/main/research/' + CAPSULE
PROOF = 'UNIFORM_FIRST_VARIATION_WAVEFORM.md'
ASSET = 'shell_first_variation_waveform_20261015.png'
OLD_VERSION = '2026.10.15-shell-response-error-certificate (GitHub checkpoint)'
NEW_VERSION = '2026.10.15-uniform-first-variation-waveform (GitHub checkpoint)'
OLD_TREE = 'https://github.com/maldonado-research/HDblast/tree/main/research/HDBLAST_PHYSICAL_SHELL_RESPONSE_20261015'
EN = {
    'title': 'Conditional first-variation shell waveform',
    'result': 'For the same specified narrow aligned quadrupole (delta = 1/1000, ell = 10^-10 in model coordinates), a continuous 4098-cell enclosure resolves the physical scalar’s linear first variation. The normalized trace H differs from its explicit comparator by less than 1.0000004 × 10^-7 on 0 ≤ x ≤ 3; the negative physical gain reverses its sign. The leading proper-time response coefficient now satisfies 6.06 × 10^15 < C < 7.44 × 10^15 in R_epsilon = epsilon^2 C + o(epsilon^2).',
    'limits': 'The retained background, exact-data, gauge, regularity and nonlinear support hypotheses remain explicit. The post-pulse estimate |H| < 4.752 × 10^-14 on 2 ≤ x ≤ 3 is a small-tail bound. No finite-amplitude nonlinear waveform, observed time scale or higher-dimensional Big Bang origin is established. Next: a wider state with a validated background and an effective nonlinear remainder estimate. Internal AI-assisted review and finite checks are distinct from external peer review. Earlier entries retain their original scope. No new Zenodo version or APS submission was made.',
    'description': 'Conditional continuous-cell enclosure of a specified first-variation shell waveform and a sharper leading response coefficient; no finite-amplitude nonlinear waveform or observed time scale.',
    'proof': 'Waveform theorem, normalization and assumptions',
    'checks': 'Reproducible calculations and independent reviews',
    'data': 'Continuous-cell enclosure data',
    'caption': 'Conditional physical first-variation envelope in model coordinates. The dashed curve is an illustrative comparator at midpoint gain; the band encloses the first variation. The inset bounds a small tail, whose shape remains unresolved.',
    'alt': 'Conditional shell scalar first-variation envelope with opposite-signed lobes and a bounded post-pulse tail',
    'short': 'For the same narrow state, a continuous 4098-cell first-variation envelope and a negative physical gain give 6.06 × 10^15 < C < 7.44 × 10^15. This conditional result supplies no finite-amplitude nonlinear waveform or observed time scale. The next target is an effective nonlinear remainder.',
}
ES = {
    'title': 'Forma de onda condicional de la primera variación de la capa',
    'result': 'Para el mismo cuadrupolo alineado estrecho especificado (delta = 1/1000, ell = 10^-10 en las coordenadas del modelo), una envolvente continua de 4098 celdas resuelve la primera variación lineal del escalar físico. La traza normalizada H difiere de su comparador explícito en menos de 1.0000004 × 10^-7 para 0 ≤ x ≤ 3; la ganancia física negativa invierte su signo. El coeficiente principal de respuesta en tiempo propio ahora satisface 6.06 × 10^15 < C < 7.44 × 10^15 en R_epsilon = epsilon^2 C + o(epsilon^2).',
    'limits': 'Se conservan las hipótesis explícitas sobre el fondo, los datos exactos, el gauge, la regularidad y el soporte no lineal. La estimación posterior al pulso |H| < 4.752 × 10^-14 para 2 ≤ x ≤ 3 acota una cola pequeña. No se establece una forma de onda no lineal de amplitud finita, una escala temporal observada ni un origen del Big Bang en dimensiones superiores. Próximo objetivo: un estado más ancho con un fondo validado y una estimación efectiva del resto no lineal. Las revisiones internas asistidas por IA y las verificaciones finitas son distintas de la revisión externa. Las entradas anteriores conservan su alcance original. No se creó una nueva versión de Zenodo ni se realizó un envío a APS.',
    'description': 'Envolvente condicional continua de la primera variación de una forma de onda de la capa y una cota más precisa del coeficiente principal; sin forma de onda no lineal de amplitud finita ni escala temporal observada.',
    'proof': 'Teorema de la forma de onda, normalización e hipótesis',
    'checks': 'Cálculos reproducibles y revisiones independientes',
    'data': 'Datos de la envolvente continua por celdas',
    'caption': 'Envolvente condicional de la primera variación física en las coordenadas del modelo. La curva discontinua es un comparador ilustrativo con ganancia intermedia; la banda encierra la primera variación. El recuadro acota una cola pequeña cuya forma sigue sin resolverse.',
    'alt': 'Envolvente condicional de la primera variación del escalar de la capa con lóbulos de signos opuestos y una cola posterior al pulso acotada',
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def baseline(folder, name):
    return subprocess.check_output(['git', '-C', str(ROOT / folder), 'show', 'HEAD:' + name]).decode('utf-8')


def markdown(words):
    return ('<!-- ' + MARK + ' -->\n**' + words['title'] + ':** ' + words['result']
            + '\n\n' + words['limits'] + '\n\n[' + words['proof'] + ']('
            + BASE + PROOF + ') · [' + words['checks'] + '](' + BASE
            + 'README.md) · [' + words['data'] + '](' + BASE
            + 'computed/WAVEFORM_ENCLOSURE.json).\n\n')


def block(words, image_path):
    e = html.escape
    return ('<div data-hdb-first-variation-waveform="' + MARK + '"><p><strong>'
            + e(words['title']) + ':</strong> ' + e(words['result']) + '</p><p>'
            + e(words['limits']) + '</p><figure><a href="' + BASE + PROOF
            + '"><img src="' + image_path + '" alt="' + e(words['alt'])
            + '" loading="lazy" width="2000" height="1480" style="display:block;max-width:100%;height:auto"></a>'
            + '<figcaption>' + e(words['caption']) + '</figcaption></figure><p><a href="'
            + BASE + PROOF + '">' + e(words['proof']) + '</a> · <a href="'
            + BASE + 'README.md">' + e(words['checks']) + '</a> · <a href="'
            + BASE + 'computed/WAVEFORM_ENCLOSURE.json">' + e(words['data'])
            + '</a>.</p></div>\n\n')


edits, records = [], []
for folder, paths in [
        ('publish-HDblast', ['README.md', 'hdblast/README.md', 'hdblast/CURRENT_STATUS.md', 'hdblast/PROGRAM_MAP.md']),
        ('publish-shared-site', ['README.md'])]:
    for name in paths:
        old = baseline(folder, name)
        require(MARK not in old, 'Checkpoint marker already exists in baseline: ' + name)
        anchor = '<!-- shell-response-error-certificate-20261015 -->'
        require(old.count(anchor) == 1, 'Expected one previous certificate anchor: ' + name)
        addition = markdown(EN)
        new = old.replace(anchor, addition + anchor, 1)
        require(new.replace(addition, '', 1) == old, 'Historical Markdown changed: ' + name)
        edits.append((folder, name, old, new))
        records.append({'repository': folder, 'path': name, 'all_prior_bytes_recovered': True})

for folder, name, words, image_path in [
        ('publish-HDblast', 'site/index.html', EN, 'assets/' + ASSET),
        ('publish-HDblast', 'site/es/index.html', ES, '../assets/' + ASSET),
        ('publish-shared-site', 'projects/hdblast/index.html', EN, 'assets/' + ASSET)]:
    old = baseline(folder, name)
    require(MARK not in old, 'Checkpoint marker already exists in baseline: ' + name)
    anchor = '<div data-hdb-response-certificate="shell-response-error-certificate-20261015">'
    require(old.count(anchor) == 1, 'Expected one previous HTML certificate anchor: ' + name)
    addition = block(words, image_path)
    new = old.replace(anchor, addition + anchor, 1)
    match = re.search(r'<meta name="description" content="([^"]+)"', old)
    require(match is not None, 'Missing HTML description: ' + name)
    old_description = match.group(1)
    require(old.count(old_description) == 4, 'Expected four synchronized description fields: ' + name)
    new = new.replace(old_description, words['description'])
    version_count = old.count(OLD_VERSION)
    require(version_count in (1, 2), 'Unexpected current-version count: ' + name)
    new = new.replace(OLD_VERSION, NEW_VERSION)
    old_suffix = 'conditional leading-coefficient enclosure in model units; numerical waveform and effective nonlinear amplitude bounds remain open.'
    new_suffix = 'conditional continuous first-variation waveform in model coordinates; effective nonlinear amplitude and remainder bounds remain open.'
    if version_count == 2:
        require(new.count(old_suffix) == 1, 'Missing current version-line scope: ' + name)
        new = new.replace(old_suffix, new_suffix, 1)
    require(new.count('"isBasedOn": "' + OLD_TREE + '"') == 1, 'Unexpected current JSON-LD source: ' + name)
    new = new.replace('"isBasedOn": "' + OLD_TREE + '"', '"isBasedOn": "' + TREE + '"', 1)
    subject = json.dumps({'@type': 'CreativeWork', 'url': BASE + PROOF,
                          'name': words['proof']}, ensure_ascii=False) + ', '
    require(new.count('"subjectOf": [') == 1, 'Unexpected JSON-LD subject count: ' + name)
    new = new.replace('"subjectOf": [', '"subjectOf": [' + subject, 1)
    recovered = new.replace(addition, '', 1).replace(words['description'], old_description).replace(NEW_VERSION, OLD_VERSION)
    recovered = recovered.replace('"isBasedOn": "' + TREE + '"', '"isBasedOn": "' + OLD_TREE + '"', 1)
    recovered = recovered.replace('"subjectOf": [' + subject, '"subjectOf": [', 1)
    if version_count == 2:
        recovered = recovered.replace(new_suffix, old_suffix, 1)
    require(recovered == old, 'Historical HTML bytes changed outside registered current fields: ' + name)
    for payload in re.findall(r'<script type="application/ld\+json">(.*?)</script>', new, re.S):
        json.loads(payload)
    resolved_image = (ROOT / folder / name).parent / image_path
    expected_image = (ROOT / folder / ('site/assets/' if folder == 'publish-HDblast' else 'projects/hdblast/assets/') / ASSET)
    require(resolved_image.resolve() == expected_image.resolve(), 'Local image path does not resolve: ' + name)
    edits.append((folder, name, old, new))
    records.append({'repository': folder, 'path': name, 'all_prior_body_bytes_recovered': True,
                    'description_fields_updated': 4, 'version_fields_updated': version_count,
                    'current_jsonld_source_updated': 1, 'jsonld_historical_subjects_preserved': True,
                    'local_figure_path': image_path})

folder, name = 'publish-shared-site', 'index.html'
old = baseline(folder, name)
require(MARK not in old, 'Checkpoint marker already exists in shared-home baseline')
anchor = '<p data-hdb-response-certificate="shell-response-error-certificate-20261015">'
require(old.count(anchor) == 1, 'Expected one previous shared-home certificate anchor')
addition = ('<p data-hdb-first-variation-waveform="' + MARK + '"><strong>'
            + html.escape(EN['title']) + ':</strong> ' + html.escape(EN['short'])
            + ' <a href="' + BASE + PROOF + '">' + html.escape(EN['proof']) + '</a>.</p>')
new = old.replace(anchor, addition + anchor, 1)
require(new.replace(addition, '', 1) == old, 'Other projects or previous shared-home content changed')
edits.append((folder, name, old, new))
records.append({'repository': folder, 'path': name, 'all_prior_bytes_recovered': True,
                'other_projects_and_shared_home_jsonld_preserved': True})

# Validate all intended writes before changing any existing file.
require(len(edits) == 9, 'Owned summary inventory differs from nine')
for folder, name, old, new in edits:
    require((ROOT / folder / name).read_text() in (old, new), 'Unexpected preexisting user edit: ' + name)
source_figure = ROOT / 'figures/SHELL_FIRST_VARIATION_WAVEFORM.png'
figure_raw = source_figure.read_bytes()
require(figure_raw.startswith(b'\x89PNG\r\n\x1a\n'), 'Figure is not PNG')
asset_paths = [('publish-HDblast', 'site/assets/' + ASSET),
               ('publish-shared-site', 'projects/hdblast/assets/' + ASSET)]
for folder, name in asset_paths:
    target = ROOT / folder / name
    tracked = subprocess.run(['git', '-C', str(ROOT / folder), 'cat-file', '-e', 'HEAD:' + name], capture_output=True)
    require(tracked.returncode != 0, 'Asset would overwrite an existing tracked file: ' + name)
    require(not target.exists() or target.read_bytes() == figure_raw, 'Asset would overwrite an unrelated local file: ' + name)

for folder, name, old, new in edits:
    path = ROOT / folder / name
    path.write_text(new)
    record = next(x for x in records if x['repository'] == folder and x['path'] == name)
    record.update({'baseline_sha256': digest(old.encode()), 'updated_sha256': digest(new.encode()),
                   'baseline_bytes': len(old.encode()), 'updated_bytes': len(new.encode())})
for folder, name in asset_paths:
    target = ROOT / folder / name
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(figure_raw)

preserved = []
for folder in ['publish-HDblast', 'publish-shared-site']:
    wt = ROOT / folder
    owned = {r['path'] for r in records if r['repository'] == folder}
    names = subprocess.check_output(['git', '-C', str(wt), 'ls-tree', '-r', '--name-only', 'HEAD'], text=True).splitlines()
    changed = set(subprocess.check_output(['git', '-C', str(wt), 'diff', '--name-only', 'HEAD'], text=True).splitlines())
    require(changed == owned, 'An existing file outside the owned summaries changed: ' + folder + ' ' + repr(changed ^ owned))
    preserved.append({'repository': folder,
                      'base_sha': subprocess.check_output(['git', '-C', str(wt), 'rev-parse', 'HEAD'], text=True).strip(),
                      'other_existing_files_preserved': len(names) - len(owned)})
    subprocess.run(['git', '-C', str(wt), 'diff', '--check', '--'] + sorted(owned), check=True)

receipt = {'status': 'PASS', 'checkpoint': MARK, 'capsule': CAPSULE,
           'changed_owned_summaries': records, 'preservation': preserved,
           'new_assets': [{'repository': folder, 'path': name, 'bytes': len(figure_raw),
                           'sha256': digest(figure_raw), 'identical_to_source_figure': True}
                          for folder, name in asset_paths],
           'wording': {'en': EN, 'es': ES},
           'scope': 'All historical body bytes and other-project content preserved. Only new checkpoint entries, registered current description/version/source metadata and two new identical figure assets were added. No Git mutations or science edits.'}
(ROOT / 'publication/SUMMARY_INTEGRATION_PREFLIGHT.json').write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + '\n')
(ROOT / 'publication/APPROVED_SUMMARY_WORDING.json').write_text(json.dumps({'en': EN, 'es': ES}, ensure_ascii=False, indent=2) + '\n')
figure_receipt = {'status': 'PASS', 'checkpoint': MARK,
                  'source': {'path': 'figures/SHELL_FIRST_VARIATION_WAVEFORM.png',
                             'bytes': len(figure_raw), 'sha256': digest(figure_raw)},
                  'site_assets': [
                      {'repository': 'maldonado-research/HDblast',
                       'path': 'site/assets/' + ASSET,
                       'url': 'https://maldonado-research.github.io/HDblast/assets/' + ASSET,
                       'bytes': len(figure_raw), 'sha256': digest(figure_raw)},
                      {'repository': 'maldonado-research/maldonado-research.github.io',
                       'path': 'projects/hdblast/assets/' + ASSET,
                       'url': 'https://maldonado-research.github.io/projects/hdblast/assets/' + ASSET,
                       'bytes': len(figure_raw), 'sha256': digest(figure_raw)}],
                  'all_asset_bytes_equal_to_source': True,
                  'html_usage': [{'repository': r['repository'], 'path': r['path'],
                                  'image_path': r['local_figure_path']}
                                 for r in records if 'local_figure_path' in r]}
(ROOT / 'publication/FIGURE_INTEGRATION_RECEIPT.json').write_text(json.dumps(figure_receipt, indent=2) + '\n')
print(json.dumps({'status': 'PASS', 'summaries': len(records), 'identical_png_assets': 2,
                  'other_existing_files_preserved': sum(x['other_existing_files_preserved'] for x in preserved)}))
