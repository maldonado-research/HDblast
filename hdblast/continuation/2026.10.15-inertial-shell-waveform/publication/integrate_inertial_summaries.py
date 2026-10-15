"""Insert the dated inertial first-variation result into nine existing summaries."""
import argparse
import hashlib
import html
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parent
MARK = 'inertial-shell-waveform-20261015'
BASE = 'https://github.com/maldonado-research/HDblast/blob/main/research/HDBLAST_INERTIAL_SHELL_WAVEFORM_20261015/'
EVIDENCE = 'https://github.com/maldonado-research/HDblast/blob/main/hdblast/continuation/2026.10.15-inertial-shell-waveform/PUBLICATION_RECEIPT.json'
OLD_VERSION = '2026.10.15-uniform-first-variation-waveform (GitHub checkpoint)'
NEW_VERSION = '2026.10.15-inertial-shell-waveform (GitHub checkpoint)'
OLD_SUFFIX = 'conditional continuous first-variation waveform in model coordinates; effective nonlinear amplitude and remainder bounds remain open.'
NEW_SUFFIX = 'conditional inertial first-variation waveform and signed post-pulse tail in model coordinates; effective nonlinear amplitude and remainder bounds remain open.'
OWNED = {
    'hdblast-worktree': ['README.md', 'hdblast/README.md', 'hdblast/CURRENT_STATUS.md', 'hdblast/PROGRAM_MAP.md', 'site/index.html', 'site/es/index.html'],
    'shared-worktree': ['README.md', 'index.html', 'projects/hdblast/index.html'],
}
EN = {
    'title': 'Conditional inertial shell waveform',
    'result': 'For a newly specified wider aligned quadrupole (delta = 1/1000, ell = 10^-6 in model coordinates, T = 3ell and Z = 8ell), exact interval bounds resolve the shell response at first variation. The normalized trace H(x) = ell D(ell x) satisfies 3.2 × 10^-8 < H < 5.3 × 10^-8 throughout 2 ≤ x ≤ 3. The negative physical gain G makes the corresponding scaled physical first variation ell S negative there. At least one zero lies in each of (1.49998, 1.49999) and (1.96, 1.97); uniqueness is not asserted. The leading proper-time response coefficient satisfies 6.06 × 10^11 < C < 7.44 × 10^11 in R_epsilon = epsilon^2 C + o(epsilon^2), in the stated model normalization.',
    'limits': 'The actual-background, exact-data, gauge, regularity and whole nonlinear correction-support hypotheses remain explicit. This first-variation result does not close the full second-order nonlinear waveform or provide an effective epsilon or nonlinear remainder bound. No observational scale or higher-dimensional Big Bang origin is established. Internal AI-assisted review and finite checks are distinct from external peer review and formal PDE verification. Earlier entries retain their dated scope. No new Zenodo version or APS submission was made.',
    'description': 'Conditional inertial first-variation shell waveform for a specified wider quadrupole, with a signed post-pulse tail and two zero brackets; no finite-amplitude nonlinear waveform or cosmological-origin claim.',
    'proof': 'Waveform theorem, normalization and assumptions',
    'checks': 'Reproducible calculations and independent reviews',
    'evidence': 'Publication checks remain pending until documented in the separate final additive evidence receipt.',
    'evidence_label': 'Final additive publication evidence',
    'alt': 'Normalized first-variation inertial shell waveform, with its post-pulse interval enclosure.',
    'caption': 'First variation for the specified wider state in model coordinates. The post-pulse trace H is positive; the negative physical gain reverses the sign of the scaled physical scalar response ell S. This is not a finite-amplitude nonlinear waveform.',
}
ES = {
    'title': 'Forma de onda inercial condicional de la capa',
    'result': 'Para un nuevo cuadrupolo alineado de mayor anchura (delta = 1/1000, ell = 10^-6 en las coordenadas del modelo, T = 3ell y Z = 8ell), cotas de intervalos exactas resuelven la respuesta de la capa en primera variación. La traza normalizada H(x) = ell D(ell x) satisface 3.2 × 10^-8 < H < 5.3 × 10^-8 en todo 2 ≤ x ≤ 3. La ganancia física negativa G hace negativa allí la primera variación física escalada ell S. Existe al menos un cero en cada intervalo (1.49998, 1.49999) y (1.96, 1.97); no se afirma unicidad. El coeficiente principal de respuesta en tiempo propio satisface 6.06 × 10^11 < C < 7.44 × 10^11 en R_epsilon = epsilon^2 C + o(epsilon^2), con la normalización indicada del modelo.',
    'limits': 'Se conservan las hipótesis explícitas sobre el fondo real, los datos exactos, el gauge, la regularidad y el soporte de toda la corrección no lineal. Este resultado de primera variación no cierra la forma de onda no lineal completa de segundo orden ni proporciona una cota efectiva de epsilon o del resto no lineal. No se establece una escala observacional ni un origen del Big Bang en dimensiones superiores. La revisión interna asistida por IA y las verificaciones finitas son distintas de la revisión externa y de la verificación formal de las EDP. Las entradas anteriores conservan su alcance fechado. No se creó una nueva versión de Zenodo ni se realizó un envío a APS.',
    'description': 'Forma de onda inercial condicional de la primera variación de la capa para un cuadrupolo de mayor anchura, con cola de signo definido y dos intervalos de ceros; sin forma de onda no lineal de amplitud finita ni afirmación sobre el origen cosmológico.',
    'proof': 'Teorema de la forma de onda, normalización e hipótesis',
    'checks': 'Cálculos reproducibles y revisiones independientes',
    'evidence': 'Las verificaciones de publicación quedan pendientes hasta que se documenten en el registro final aditivo independiente.',
    'evidence_label': 'Registro final aditivo de publicación',
    'alt': 'Forma de onda inercial normalizada de primera variación de la capa, con intervalo certificado posterior al pulso.',
    'caption': 'Primera variación del estado de mayor anchura especificado, en coordenadas del modelo. La traza H después del pulso es positiva; la ganancia física negativa invierte el signo de la respuesta escalar física escalada ell S. No es una forma de onda no lineal de amplitud finita.',
}


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def git(wt, *args):
    return subprocess.check_output(['git', '-C', str(wt), *args])


def md(words):
    return ('<!-- '+MARK+' -->\n**'+words['title']+':** '+words['result']+'\n\n'+words['limits']
            +'\n\n['+words['proof']+']('+BASE+'INERTIAL_SHELL_WAVEFORM.md) · ['+words['checks']+']('
            +BASE+'README.md). '+words['evidence']+' ['+words['evidence_label']+']('+EVIDENCE+').\n\n')


def block(words, image):
    e = html.escape
    return ('<div data-hdb-inertial-waveform="'+MARK+'"><p><strong>'+e(words['title'])+':</strong> '
            +e(words['result'])+'</p><p>'+e(words['limits'])+'</p><p><a href="'+BASE
            +'INERTIAL_SHELL_WAVEFORM.md">'+e(words['proof'])+'</a> · <a href="'+BASE+'README.md">'
            +e(words['checks'])+'</a>.</p><p>'+e(words['evidence'])+' <a href="'+EVIDENCE+'">'
            +e(words['evidence_label'])+'</a>.</p><figure><a href="'+BASE+'computed/INERTIAL_SHELL_WAVEFORM.png">'
            +'<img src="'+image+'" alt="'+e(words['alt'], quote=True)
            +'" loading="lazy" style="max-width:100%;height:auto"></a><figcaption>'+e(words['caption'])
            +'</figcaption></figure></div>\n\n')


def detailed(old, words, image, shared=False):
    anchor = '<div data-hdb-first-variation-waveform="uniform-first-variation-waveform-20261015">'
    require(old.count(anchor) == 1 and MARK not in old, 'Unique prior waveform insertion anchor required')
    addition = block(words, image)
    new = old.replace(anchor, addition+anchor, 1)
    desc = re.search(r'<meta name="description" content="([^"]+)"', old).group(1)
    require(old.count(desc) == 4, 'Four prior description fields required')
    new = new.replace(desc, words['description'])
    version_count = old.count(OLD_VERSION)
    require(version_count == (2 if shared else 1), 'Expected current version fields differ')
    new = new.replace(OLD_VERSION, NEW_VERSION)
    if shared:
        require(new.count(OLD_SUFFIX) == 1, 'Shared current-version suffix differs')
        new = new.replace(OLD_SUFFIX, NEW_SUFFIX, 1)
        pattern = r'("@type": "WebPage".*?"dateModified": ")2026-10-14(")'
        new, changed_dates = re.subn(pattern, r'\g<1>2026-10-15\g<2>', new, count=1)
        require(changed_dates == 1, 'Shared current WebPage date differs')
    return new, {'description': desc, 'addition': addition, 'version_count': version_count}


class Parser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.meta, self.links, self.jsonld = {}, [], []
        self.in_jsonld = False

    def handle_starttag(self, tag, attrs):
        values = dict(attrs)
        if tag == 'meta':
            self.meta[values.get('name') or values.get('property')] = values.get('content')
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
            old = git(wt, 'show', base+':'+name).decode()
            new = (wt/name).read_text()
            require(new.count(MARK) == 1, 'Exactly one new checkpoint entry required: '+name)
            if name.endswith('.md'):
                addition = md(EN)
                require(new.count(addition) == 1 and new.replace(addition, '', 1) == old,
                        'Historical Markdown bytes changed: '+name)
                records.append({'repository': folder, 'path': name, 'all_prior_bytes_preserved': True})
                continue
            shared_home = folder == 'shared-worktree' and name == 'index.html'
            if shared_home:
                anchor = '<p data-hdb-first-variation-waveform="uniform-first-variation-waveform-20261015">'
                addition = ('<p data-hdb-inertial-waveform="'+MARK+'"><strong>'+html.escape(EN['title'])
                            +':</strong> A specified wider quadrupole has a signed post-pulse first-variation tail, two zero brackets and a leading coefficient interval in model units. No zero uniqueness, finite-amplitude nonlinear waveform or effective epsilon is asserted. <a href="'
                            +BASE+'INERTIAL_SHELL_WAVEFORM.md">Result and assumptions</a>. Publication evidence remains pending until the <a href="'
                            +EVIDENCE+'">final additive receipt</a>.</p>')
                require(new.replace(addition, '', 1) == old, 'Shared homepage historical bytes changed')
                records.append({'repository': folder, 'path': name, 'all_prior_bytes_preserved': True,
                                'only_hdblast_card_modified': True})
                continue
            words = ES if name == 'site/es/index.html' else EN
            image = ('../assets/inertial_shell_waveform_20261015.png' if name == 'site/es/index.html'
                     else 'assets/inertial_shell_waveform_20261015.png')
            expected, details = detailed(old, words, image, folder == 'shared-worktree')
            require(new == expected, 'Detailed page differs from reviewed insertion/metadata replacements: '+name)
            parser = Parser()
            parser.feed(new)
            values = [parser.meta[k] for k in ('description', 'og:description', 'twitter:description')]
            json_desc = [node['description'] for data in parser.jsonld for node in data['@graph']
                         if node.get('@type') in ('Dataset', 'CreativeWork')]
            require(len(json_desc) == 1 and len(set(values+json_desc)) == 1, 'Four descriptions disagree: '+name)
            require(re.findall(r'<nav\b[^>]*>.*?</nav>', old, re.S) == re.findall(r'<nav\b[^>]*>.*?</nav>', new, re.S),
                    'Navigation changed: '+name)
            old_parser = Parser()
            old_parser.feed(old)
            require(set(old_parser.links) <= set(parser.links), 'Historical link removed: '+name)
            image_path = ((wt/name).parent/image).resolve()
            if require_assets:
                require(image_path.is_file(), 'Pending figure asset not yet installed: '+str(image_path))
            records.append({'repository': folder, 'path': name, 'all_historical_body_bytes_preserved': True,
                            'description_fields_consistent': 4, 'current_version_fields_updated': details['version_count'],
                            'current_modified_date_updated': folder == 'shared-worktree',
                            'historical_links_navigation_retained': True,
                            'image_repository_path': image_path.relative_to(wt).as_posix(),
                            'image_present': image_path.is_file()})
        preservation.append({'repository': folder, 'base': base, 'other_existing_files_preserved': len(old_names)-len(names)})
    assets = [{'repository': 'maldonado-research/HDblast', 'path': 'site/assets/inertial_shell_waveform_20261015.png',
               'url': 'https://maldonado-research.github.io/HDblast/assets/inertial_shell_waveform_20261015.png'},
              {'repository': 'maldonado-research/maldonado-research.github.io', 'path': 'projects/hdblast/assets/inertial_shell_waveform_20261015.png',
               'url': 'https://maldonado-research.github.io/projects/hdblast/assets/inertial_shell_waveform_20261015.png'}]
    receipt = {'status': 'PASS_SUMMARY_CONTENT_FIGURE_ASSETS_PENDING' if not require_assets else 'PASS',
               'checkpoint': MARK, 'base_commits': bases, 'changed_owned_summaries': records, 'preservation': preservation,
               'wording': {'en': EN, 'es': ES}, 'planned_site_assets': assets,
               'figure_assets_installed_by_this_tool': False, 'publication_evidence_status': 'PENDING_FINAL_ADDITIVE_RECEIPT',
               'site_build': 'Existing GitHub Pages workflow deploys static site directly; no source generator or deployment workflow changed.',
               'scope': 'Nine additive summary edits plus current metadata replacements. Old research, links, other projects and workflows preserved. Capsule/workflow/asset copying and Git operations remain with root.'}
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--review-only', action='store_true')
    parser.add_argument('--require-assets', action='store_true')
    args = parser.parse_args()
    baseline_path = ROOT/'SITE_BASELINES.json'
    if baseline_path.exists():
        bases = json.loads(baseline_path.read_text())
    else:
        require(not args.review_only, 'Saved site baselines required for review')
        bases = {folder: git(ROOT/folder, 'rev-parse', 'HEAD').decode().strip() for folder in OWNED}
        require(bases['hdblast-worktree'] == 'c490f17139f1266426995595f5ba0eaefb14c636'
                and bases['shared-worktree'].startswith('7246'), 'Requested isolated base commits differ')
        baseline_path.write_text(json.dumps(bases, indent=2)+'\n')
    if not args.review_only:
        for folder, names in OWNED.items():
            wt, base = ROOT/folder, bases[folder]
            for name in names:
                old = git(wt, 'show', base+':'+name).decode()
                require(MARK not in old, 'Checkpoint already present in baseline')
                if name.endswith('.md'):
                    addition = md(EN)
                    anchor = '<!-- uniform-first-variation-waveform-20261015 -->'
                    if old.count(anchor) == 1:
                        new = old.replace(anchor, addition+anchor, 1)
                    else:
                        require(name == 'hdblast/PROGRAM_MAP.md' and old.count(anchor) == 0, 'Markdown insertion anchor differs')
                        first, rest = old.split('\n\n', 1)
                        new = first+'\n\n'+addition+rest
                elif folder == 'shared-worktree' and name == 'index.html':
                    anchor = '<p data-hdb-first-variation-waveform="uniform-first-variation-waveform-20261015">'
                    require(old.count(anchor) == 1, 'Shared HDBLAST card anchor differs')
                    addition = ('<p data-hdb-inertial-waveform="'+MARK+'"><strong>'+html.escape(EN['title'])
                                +':</strong> A specified wider quadrupole has a signed post-pulse first-variation tail, two zero brackets and a leading coefficient interval in model units. No zero uniqueness, finite-amplitude nonlinear waveform or effective epsilon is asserted. <a href="'
                                +BASE+'INERTIAL_SHELL_WAVEFORM.md">Result and assumptions</a>. Publication evidence remains pending until the <a href="'
                                +EVIDENCE+'">final additive receipt</a>.</p>')
                    new = old.replace(anchor, addition+anchor, 1)
                else:
                    words = ES if name == 'site/es/index.html' else EN
                    image = '../assets/inertial_shell_waveform_20261015.png' if name == 'site/es/index.html' else 'assets/inertial_shell_waveform_20261015.png'
                    new, _ = detailed(old, words, image, folder == 'shared-worktree')
                path = wt/name
                require(path.read_text() in (old, new), 'Existing independent summary edit must be preserved: '+name)
                path.write_text(new)
    receipt = review(bases, args.require_assets)
    (ROOT/'SITE_INTEGRATION_REVIEW.json').write_text(json.dumps(receipt, indent=2)+'\n')
    # The installer consumes this established filename as well.
    (ROOT/'SUMMARY_INTEGRATION_PREFLIGHT.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({'status': receipt['status'], 'summaries': len(receipt['changed_owned_summaries']),
                      'other_existing_files_preserved': sum(r['other_existing_files_preserved'] for r in receipt['preservation']),
                      'figure_assets': 'checked' if args.require_assets else 'pending root installation'}, indent=2))


if __name__ == '__main__':
    main()
