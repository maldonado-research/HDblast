"""Add the dated quadratic material tension result to nine existing summaries.

Only the nine OWNED summaries and this helper's three JSON receipts are written.
Figure installation, workflow changes and Git mutations are outside this tool.
"""
import argparse
import hashlib
import html
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import subprocess


ROOT = Path(__file__).resolve().parent
MARK = 'quadratic-shell-forcing-20261015'
PRIOR_MARK = 'inertial-shell-waveform-20261015'
BASE = 'https://github.com/maldonado-research/HDblast/blob/main/research/HDBLAST_QUADRATIC_SHELL_FORCING_20261015/'
NOTE = 'QUADRATIC_MATERIAL_TENSION_FORCING.md'
EVIDENCE = 'https://github.com/maldonado-research/HDblast/blob/main/hdblast/continuation/2026.10.15-quadratic-shell-forcing/PUBLICATION_RECEIPT.json'
OLD_VERSION = '2026.10.15-inertial-shell-waveform (GitHub checkpoint)'
NEW_VERSION = '2026.10.15-quadratic-shell-forcing (GitHub checkpoint)'
OLD_SUFFIX = 'conditional inertial first-variation waveform and signed post-pulse tail in model coordinates; effective nonlinear amplitude and remainder bounds remain open.'
NEW_SUFFIX = 'enclosed quadratic material tension contribution and negative mean tail in model units; the complete coordinate source and full second-order response remain unevaluated.'
ASSET = 'quadratic_material_forcing_20261015.png'
EXPECTED_BASES = {
    'hdblast-worktree': '6d4070dd7fc85dfc61c3e89484b083077bb3b472',
    'shared-worktree': '751a79d93c1ac9969ab49b727f509dce2852c74f',
}
OWNED = {
    'hdblast-worktree': ['README.md', 'hdblast/README.md', 'hdblast/CURRENT_STATUS.md',
                         'hdblast/PROGRAM_MAP.md', 'site/index.html', 'site/es/index.html'],
    'shared-worktree': ['README.md', 'index.html', 'projects/hdblast/index.html'],
}
EN = {
    'title': '15 October 2026: quadratic material tension forcing',
    'result': 'For the previously certified wider quadrupole (delta = 1/1000, ell = 10^-6 in model coordinates), the explicit material scalar-junction tension contribution is Q_phi = -4 psi_1^2 in the outward-normal SUM convention. It has scalar degrees 0, 2 and 4. Exact continuous enclosures give -2.3 × 10^-19 < ell^2 mean(Q_phi) < -6.7 × 10^-20 throughout 2 ≤ x ≤ 3 and 5.778 × 10^21 < ||Q_phi||_L2(dmu_0)^2 < 8.644 × 10^21 over 0 ≤ x ≤ 3, using the background proper-time measure and model units. The angular mean stays negative throughout the tail; the field vanishes on Y = 0.',
    'limits': 'This isolates the explicit material tension contribution; the complete quadratic source for bulk-coordinate unknowns and full U2 are not numerically evaluated. The earlier conditional construction of the complete second coefficient is preserved under its stated data and evolution assumptions. No effective finite amplitude, nonlinear remainder, observational scale or cosmological origin is established. Independent internal exact checks are distinct from external peer review and formal PDE verification. Earlier dated entries retain their scope. No new Zenodo version or APS submission was made.',
    'description': 'Enclosed quadratic material tension contribution at the shell, with scalar degrees 0, 2 and 4 and a negative continuous mean tail; the complete coordinate source and full second-order response remain unevaluated.',
    'proof': 'Material tension theorem, normalization and assumptions',
    'checks': 'Reproducible calculations and independent internal checks',
    'evidence': 'Publication checks remain pending until documented in the separate final additive evidence receipt.',
    'evidence_label': 'Final additive publication evidence',
    'alt': 'Enclosures of the quadratic material tension forcing, including its degree-zero mean and continuous post-pulse tail.',
    'caption': 'Explicit material tension forcing for the specified wider state in model units. The continuous angular mean is strictly negative on 2 ≤ x ≤ 3; the angular field is nonpositive and vanishes where Y = 0. The complete coordinate source and full second-order response are not evaluated here.',
}
ES = {
    'title': '15 de octubre de 2026: fuente cuadrática de la tensión material',
    'result': 'Para el cuadrupolo de mayor anchura ya certificado (delta = 1/1000, ell = 10^-6 en coordenadas del modelo), la contribución explícita de la tensión a la condición de unión escalar material es Q_phi = -4 psi_1^2 en la convención SUM de normales exteriores. Contiene los grados escalares 0, 2 y 4. Las cotas continuas exactas dan -2.3 × 10^-19 < ell^2 media(Q_phi) < -6.7 × 10^-20 en todo 2 ≤ x ≤ 3 y 5.778 × 10^21 < ||Q_phi||_L2(dmu_0)^2 < 8.644 × 10^21 sobre 0 ≤ x ≤ 3, con la medida de tiempo propio del fondo y unidades del modelo. La media angular permanece negativa en toda la cola; el campo se anula donde Y = 0.',
    'limits': 'Se aísla la contribución explícita de la tensión material; no se evalúan numéricamente la fuente cuadrática completa para las incógnitas en coordenadas del volumen ni U2 completo. Se conserva la construcción condicional anterior del segundo coeficiente completo bajo sus hipótesis sobre los datos y la evolución. No se establece una amplitud finita efectiva, un resto no lineal, una escala observacional ni un origen cosmológico. Las verificaciones internas exactas e independientes son distintas de la revisión externa y de la verificación formal de las EDP. Las entradas anteriores conservan su alcance fechado. No se creó una nueva versión de Zenodo ni se realizó un envío a APS.',
    'description': 'Contribución cuadrática acotada de la tensión material en la capa, con grados escalares 0, 2 y 4 y una cola continua de media negativa; la fuente completa en coordenadas y la respuesta completa de segundo orden quedan sin evaluar.',
    'proof': 'Teorema de la tensión material, normalización e hipótesis',
    'checks': 'Cálculos reproducibles y verificaciones internas independientes',
    'evidence': 'Las verificaciones de publicación quedan pendientes hasta que se documenten en el registro final aditivo independiente.',
    'evidence_label': 'Registro final aditivo de publicación',
    'alt': 'Cotas de la fuente cuadrática de la tensión material, incluidas su media de grado cero y su cola continua posterior al pulso.',
    'caption': 'Fuente explícita de la tensión material del estado de mayor anchura especificado, en unidades del modelo. La media angular continua es estrictamente negativa en 2 ≤ x ≤ 3; el campo angular no es positivo y se anula donde Y = 0. Aquí no se evalúan la fuente completa en coordenadas ni la respuesta completa de segundo orden.',
}


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def git(wt, *args):
    """Read-only Git access; this helper performs no Git mutations."""
    return subprocess.check_output(['git', '-C', str(wt), *args])


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def md(words):
    return ('<!-- '+MARK+' -->\n**'+words['title']+':** '+words['result']+'\n\n'+words['limits']
            +'\n\n['+words['proof']+']('+BASE+NOTE+') · ['+words['checks']+']('
            +BASE+'README.md). '+words['evidence']+' ['+words['evidence_label']+']('+EVIDENCE+').\n\n')


def block(words, image):
    e = html.escape
    return ('<div data-hdb-quadratic-forcing="'+MARK+'"><p><strong>'+e(words['title'])+':</strong> '
            +e(words['result'])+'</p><p>'+e(words['limits'])+'</p><p><a href="'+BASE
            +NOTE+'">'+e(words['proof'])+'</a> · <a href="'+BASE+'README.md">'
            +e(words['checks'])+'</a>.</p><p>'+e(words['evidence'])+' <a href="'+EVIDENCE+'">'
            +e(words['evidence_label'])+'</a>.</p><figure><a href="'+BASE+'computed/QUADRATIC_MATERIAL_FORCING.png">'
            +'<img src="'+image+'" alt="'+e(words['alt'], quote=True)
            +'" loading="lazy" style="max-width:100%;height:auto"></a><figcaption>'+e(words['caption'])
            +'</figcaption></figure></div>\n\n')


def homepage_addition():
    return ('<p data-hdb-quadratic-forcing="'+MARK+'"><strong>'+html.escape(EN['title'])
            +':</strong> The explicit material tension contribution Q_phi = -4 psi_1^2 has scalar degrees 0, 2 and 4, a strictly negative continuous mean tail, and an enclosed squared L2 norm in model units. The complete coordinate source and full U2 remain unevaluated; the earlier conditional second-coefficient construction is preserved. <a href="'
            +BASE+NOTE+'">Result and assumptions</a>. Publication evidence remains pending until the <a href="'
            +EVIDENCE+'">final additive receipt</a>.</p>')


def detailed(old, words, image, shared=False):
    anchor = '<div data-hdb-inertial-waveform="'+PRIOR_MARK+'">'
    require(old.count(anchor) == 1 and MARK not in old, 'Unique prior inertial insertion anchor required')
    addition = block(words, image)
    new = old.replace(anchor, addition+anchor, 1)
    match = re.search(r'<meta name="description" content="([^"]+)"', old)
    require(match is not None, 'Prior description is missing')
    desc = match.group(1)
    require(old.count(desc) == 4, 'Exactly four prior description fields required')
    new = new.replace(desc, words['description'])
    version_count = old.count(OLD_VERSION)
    require(version_count == (2 if shared else 1), 'Expected current version fields differ')
    new = new.replace(OLD_VERSION, NEW_VERSION)
    if shared:
        require(new.count(OLD_SUFFIX) == 1, 'Shared current-version suffix differs')
        new = new.replace(OLD_SUFFIX, NEW_SUFFIX, 1)
        dates = re.findall(r'"@type": "WebPage".*?"dateModified": "([^"]+)"', old)
        require(dates == ['2026-10-15'], 'Shared current WebPage date differs')
    return new, {'description': desc, 'addition': addition, 'version_count': version_count}


def transform(folder, name, old):
    require(MARK not in old, 'Checkpoint already present in baseline: '+name)
    if name.endswith('.md'):
        addition = md(EN)
        anchor = '<!-- '+PRIOR_MARK+' -->'
        require(old.count(anchor) == 1, 'Unique prior Markdown insertion anchor required: '+name)
        return old.replace(anchor, addition+anchor, 1), {'addition': addition}
    if folder == 'shared-worktree' and name == 'index.html':
        anchor = '<p data-hdb-inertial-waveform="'+PRIOR_MARK+'">'
        require(old.count(anchor) == 1, 'Unique shared HDBLAST card anchor required')
        addition = homepage_addition()
        return old.replace(anchor, addition+anchor, 1), {'addition': addition}
    words = ES if name == 'site/es/index.html' else EN
    image = ('../assets/' if name == 'site/es/index.html' else 'assets/')+ASSET
    new, detail = detailed(old, words, image, folder == 'shared-worktree')
    detail.update({'words': words, 'image': image})
    return new, detail


def reverse(new, detail, shared=False):
    """Undo only the documented addition and current metadata replacements."""
    require(new.count(detail['addition']) == 1, 'Unique reviewed addition required for reversal')
    reverted = new.replace(detail['addition'], '', 1)
    if 'description' in detail:
        require(reverted.count(detail['words']['description']) == 4,
                'Exactly four replacement descriptions required for reversal')
        reverted = reverted.replace(detail['words']['description'], detail['description'])
        require(reverted.count(NEW_VERSION) == detail['version_count'], 'Version reversal differs')
        reverted = reverted.replace(NEW_VERSION, OLD_VERSION)
        if shared:
            require(reverted.count(NEW_SUFFIX) == 1, 'Shared suffix reversal differs')
            reverted = reverted.replace(NEW_SUFFIX, OLD_SUFFIX, 1)
    return reverted


class Parser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.meta, self.links, self.jsonld = {}, [], []
        self.in_jsonld = False

    def handle_starttag(self, tag, attrs):
        values = dict(attrs)
        if tag == 'meta':
            key = values.get('name') or values.get('property')
            self.meta.setdefault(key, []).append(values.get('content'))
        if tag == 'a' and 'href' in values:
            self.links.append(values['href'])
        if tag == 'script' and values.get('type') == 'application/ld+json':
            self.in_jsonld = True

    def handle_endtag(self, tag):
        if tag == 'script':
            self.in_jsonld = False

    def handle_data(self, data):
        if self.in_jsonld:
            self.jsonld.append(json.loads(data))


def review(bases, require_assets=False):
    records, preservation = [], []
    for folder, names in OWNED.items():
        wt, base = ROOT/folder, bases[folder]
        old_names = set(git(wt, 'ls-tree', '-r', '--name-only', base).decode().splitlines())
        changed = set(git(wt, 'diff', '--name-only', base, '--').decode().splitlines()) & old_names
        require(changed == set(names), 'Unexpected modification to existing tracked files: '+folder)
        for name in names:
            old_bytes = git(wt, 'show', base+':'+name)
            new_bytes = (wt/name).read_bytes()
            old, new = old_bytes.decode(), new_bytes.decode()
            expected, detail = transform(folder, name, old)
            require(new == expected and new.count(MARK) == 1,
                    'Summary differs from the unique reviewed insertion: '+name)
            reverted = reverse(new, detail, folder == 'shared-worktree').encode()
            require(reverted == old_bytes, 'Reversal failed to recover all baseline bytes: '+name)
            record = {'repository': folder, 'path': name, 'baseline_sha256': sha256(old_bytes),
                      'updated_sha256': sha256(new_bytes), 'reversed_sha256': sha256(reverted),
                      'exact_byte_reversal_passed': True}
            if name.endswith('.md'):
                record['all_prior_bytes_preserved'] = True
                records.append(record)
                continue
            parser, old_parser = Parser(), Parser()
            parser.feed(new)
            old_parser.feed(old)
            require(re.findall(r'<nav\b[^>]*>.*?</nav>', old, re.S)
                    == re.findall(r'<nav\b[^>]*>.*?</nav>', new, re.S), 'Navigation changed: '+name)
            require(set(old_parser.links) <= set(parser.links), 'Historical link removed: '+name)
            record['historical_links_navigation_retained'] = True
            if folder == 'shared-worktree' and name == 'index.html':
                record.update({'all_prior_bytes_preserved': True, 'only_hdblast_card_modified': True})
                records.append(record)
                continue
            values = []
            for key in ('description', 'og:description', 'twitter:description'):
                require(len(parser.meta.get(key, [])) == 1, 'Unique description meta field required: '+key)
                values += parser.meta[key]
            json_desc = [node['description'] for data in parser.jsonld for node in data['@graph']
                         if node.get('@type') in ('Dataset', 'CreativeWork') and 'description' in node]
            require(len(json_desc) == 1 and values+json_desc == [detail['words']['description']]*4,
                    'Exactly four matching description fields required: '+name)
            image_path = ((wt/name).parent/detail['image']).resolve()
            if require_assets:
                require(image_path.is_file(), 'Figure asset not yet installed: '+str(image_path))
                require(image_path.read_bytes() == (ROOT.parent/'computed/QUADRATIC_MATERIAL_FORCING.png').read_bytes(),
                        'Installed site figure differs from the authoritative image: '+str(image_path))
            record.update({'all_historical_body_bytes_preserved': True,
                           'description_fields_consistent': 4,
                           'current_version_fields_updated': detail['version_count'],
                           'current_modified_date_preserved': True,
                           'image_repository_path': image_path.relative_to(wt).as_posix(),
                           'image_present': image_path.is_file()})
            records.append(record)
        preservation.append({'repository': folder, 'base': base,
                             'other_existing_files_preserved': len(old_names)-len(names)})
    assets = [{'repository': 'maldonado-research/HDblast', 'path': 'site/assets/'+ASSET,
               'url': 'https://maldonado-research.github.io/HDblast/assets/'+ASSET},
              {'repository': 'maldonado-research/maldonado-research.github.io',
               'path': 'projects/hdblast/assets/'+ASSET,
               'url': 'https://maldonado-research.github.io/projects/hdblast/assets/'+ASSET}]
    note_path = ROOT.parent/NOTE
    return {'status': 'PASS' if require_assets else 'PASS_SUMMARY_CONTENT_FIGURE_ASSETS_PENDING',
            'checkpoint': MARK, 'base_commits': bases, 'changed_owned_summaries': records,
            'preservation': preservation, 'wording': {'en': EN, 'es': ES},
            'authoritative_note': {'path': str(note_path), 'sha256': sha256(note_path.read_bytes())},
            'planned_site_assets': assets, 'figure_assets_installed_by_this_tool': False,
            'publication_evidence_status': 'PENDING_FINAL_ADDITIVE_RECEIPT',
            'site_build': 'Existing GitHub Pages workflows deploy the static sites; no source generator or deployment workflow changed.',
            'scope': 'Nine additive summary edits plus current metadata replacements. All old body bytes, links, navigation, other projects and existing workflows preserved. Capsule, asset copying and Git operations remain with root.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--review-only', action='store_true')
    parser.add_argument('--require-assets', action='store_true')
    args = parser.parse_args()
    baseline_path = ROOT/'SITE_BASELINES.json'
    if baseline_path.exists():
        bases = json.loads(baseline_path.read_text())
        require(bases == EXPECTED_BASES, 'Saved site baselines differ from the requested checkpoints')
    else:
        require(not args.review_only, 'Saved site baselines required for review')
        bases = {folder: git(ROOT/folder, 'rev-parse', 'HEAD').decode().strip() for folder in OWNED}
        require(bases == EXPECTED_BASES, 'Requested isolated base commits differ')
    if not args.review_only:
        planned = []
        for folder, names in OWNED.items():
            wt, base = ROOT/folder, bases[folder]
            for name in names:
                old_bytes = git(wt, 'show', base+':'+name)
                new, detail = transform(folder, name, old_bytes.decode())
                new_bytes = new.encode()
                require(reverse(new, detail, folder == 'shared-worktree').encode() == old_bytes,
                        'Pre-write exact byte reversal failed: '+name)
                path = wt/name
                require(path.read_bytes() in (old_bytes, new_bytes),
                        'Existing independent summary edit must be preserved: '+name)
                planned.append((path, new_bytes))
        baseline_path.write_text(json.dumps(bases, indent=2)+'\n')
        for path, new_bytes in planned:
            path.write_bytes(new_bytes)
    receipt = review(bases, args.require_assets)
    for name in ('SITE_INTEGRATION_REVIEW.json', 'SUMMARY_INTEGRATION_PREFLIGHT.json'):
        (ROOT/name).write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({'status': receipt['status'], 'summaries': len(receipt['changed_owned_summaries']),
                      'other_existing_files_preserved': sum(r['other_existing_files_preserved'] for r in receipt['preservation']),
                      'exact_byte_reversals': len(receipt['changed_owned_summaries']),
                      'figure_assets': 'checked against authoritative PNG' if args.require_assets else 'pending root installation'}, indent=2))


if __name__ == '__main__':
    main()
