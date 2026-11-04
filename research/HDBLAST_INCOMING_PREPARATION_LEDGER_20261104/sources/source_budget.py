#!/usr/bin/env python3
"""Exact postprocessing of retained continuous source boxes; no project imports."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
PINS_SHA = 'fd623a2189857ecd23eb203ec36e3b679d5ba7822d98a1cf42913b1b2ff9e02b'
DESIGN_SHA = 'cac609399486af8c97e01e4f74b1eb84882241a25ea692756b49a48787862e7f'
CLOSURE_SHA = 'dd2384c76c3cfb6e115373b39f1319e624f226f5817a8300f87a793b9438d87b'
PACKET_SHA = '591a70b5a483c9e250ed41635fd67c8f932f587ed52121c7475b833d72253a8a'
GEO = 'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/'
COMMON = 'research/HDBLAST_COMMON_GAUGE_FORCING_20261016/'
PREP = 'research/HDBLAST_PREPARED_SHELL_VOLUME_20261101/'
T = Q(5, 4)
ELL = Q(1, 10**6)
GRAM = (1, 2, 1, 3, 1)


def need(ok, reason):
    if not ok:
        raise ValueError(reason)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def read_regular(path):
    need(not any(p.is_symlink() for p in (path, *path.parents)), 'SYMLINK_PATH')
    need(path.is_file(), 'REGULAR_INPUT_REQUIRED')
    return path.read_bytes()


def authenticate_portable(repo):
    # No project module is imported; every scientific input is authenticated
    # before any source arithmetic, including the selected packet joins.
    design_raw = read_regular(HERE/'FROZEN_DESIGN.json')
    need(sha(design_raw) == DESIGN_SHA, 'DESIGN_AUTHENTICATION')
    design = json.loads(design_raw)
    need(Q(design['ell']) == ELL and Q(design['T']) == T
         and design['mean_gram'] == list(GRAM), 'DESIGN_SCALE_AND_NORM')
    closure_raw = read_regular(HERE/'SOURCE_CLOSURE_INPUTS.json')
    need(sha(closure_raw) == CLOSURE_SHA, 'CLOSURE_PIN_AUTHENTICATION')
    closure = json.loads(closure_raw)['files']
    for rel, pin in closure.items():
        need(not Path(rel).is_absolute() and '..' not in Path(rel).parts,
             'REPOSITORY_RELATIVE_SOURCE_PATH')
        raw = read_regular(repo/rel)
        need(sha(raw) == pin['sha256'] and len(raw) == pin['bytes'],
             'SOURCE_AUTHENTICATION:' + rel)
    packet_raw = read_regular(HERE/'PACKET_SOURCE_JOIN.json')
    need(sha(packet_raw) == PACKET_SHA, 'PACKET_JOIN_AUTHENTICATION')
    packet = json.loads(packet_raw)['source_files_checked']
    need(len(packet) == 7, 'PACKET_JOIN_COUNT')
    for item in packet:
        need(closure[item['path']] == {'bytes':item['bytes'],'sha256':item['sha256']},
             'PACKET_SOURCE_PIN')
        declared = json.loads(read_regular(repo/item['declared_pin_record']))
        if 'input_sha256' in declared:
            accepted = declared['input_sha256'].get(item['path']) == item['sha256']
        elif 'source_files' in declared or 'sources' in declared:
            candidates = declared.get('source_files',declared.get('sources',[]))
            accepted = any(x.get('path') == item['path'] and x.get('sha256') == item['sha256']
                           and x.get('bytes') == item['bytes'] for x in candidates)
        else:
            direct = declared.get(item['path'])
            if direct is None:
                parent = Path(item['declared_pin_record']).parent
                relative = (str(Path(item['path']).relative_to(parent))
                            if Path(item['path']).is_relative_to(parent) else '')
                direct = declared.get(relative)
            accepted = direct == {'bytes':item['bytes'],'sha256':item['sha256']}
        need(accepted, 'PACKET_DECLARED_PIN_JOIN')
    return design, closure


def authenticate(repo):
    raw = (HERE / 'SOURCE_PINS.json').read_bytes()
    need(sha(raw) == PINS_SHA, 'LOCAL_PIN_AUTHENTICATION')
    pins = json.loads(raw)['files']
    for rel, pin in pins.items():
        content = (repo / rel).read_bytes()
        need(len(content) == pin['bytes'] and sha(content) == pin['sha256'],
             'SOURCE_AUTHENTICATION:' + rel)
    # The prepared gauge proof's published pins actually bind the old operator
    # and the current exact-family identification, rather than an inferred name.
    gauge = json.loads((repo / PREP / 'gauge/SOURCE_PINS.json').read_bytes())
    for item in gauge['sources']:
        if item['key'] in ('S1', 'S2', 'S3'):
            raw = (repo / item['path']).read_bytes()
            need(sha(raw) == item['sha256'] and len(raw) == item['bytes'],
                 'PREPARED_PUBLISHED_JOIN:' + item['key'])
    old = json.loads((repo / GEO / 'numeric/PREDECESSOR_PINS.json').read_bytes())
    for rel in ('inputs/ACTUAL_FIRST_FIELDS_PREMISE.md',
                'inputs/INHERITED_APPROXIMATION_ERROR_PROOF.md'):
        raw = (repo / GEO / 'numeric' / rel).read_bytes()
        need(sha(raw) == old[rel]['sha256'], 'INHERITED_FIRST_FIELD_JOIN:' + rel)
    current_first = (repo / COMMON / 'first-fields/FIRST_HARMONIC_FIELDS_UNPUBLISHED.md').read_bytes()
    need(sha(current_first) == old['inputs/ACTUAL_FIRST_FIELDS_PREMISE.md']['sha256'],
         'EXACT_SAME_FIRST_FIELD_PROOF_BYTES')
    return pins


def interval(row):
    lo, hi = map(Q, row)
    need(lo <= hi, 'ORDERED_SOURCE_INTERVAL')
    return lo, hi


def mag(row):
    return max(map(abs, interval(row)))


def linear_interval(terms):
    lo = hi = Q(0)
    for coefficient, row in terms:
        a, b = interval(row)
        lo += coefficient * (a if coefficient >= 0 else b)
        hi += coefficient * (b if coefficient >= 0 else a)
    return lo, hi


def clip(poly, level, keep_below):
    """Exact half-plane clipping against x+r=level."""
    out = []
    if not poly:
        return out
    def inside(p):
        return (p[0] + p[1] <= level if keep_below
                else p[0] + p[1] >= level)
    for a, b in zip(poly, poly[1:] + poly[:1]):
        ia, ib = inside(a), inside(b)
        if ia:
            out.append(a)
        if ia != ib:
            fraction = (level - sum(a)) / (sum(b) - sum(a))
            out.append(tuple(a[k] + fraction * (b[k] - a[k]) for k in (0, 1)))
    return out


def area(poly):
    if len(poly) < 3:
        return Q(0)
    return abs(sum(a[0] * b[1] - a[1] * b[0]
                   for a, b in zip(poly, poly[1:] + poly[:1]))) / 2


def cell_area(x, r, active):
    a, b = map(Q, x)
    c, d = map(Q, r)
    poly = clip([(a, c), (b, c), (b, d), (a, d)], T, True)
    return area(clip(poly, Q(1), False) if active else poly)


def mean_square(sources):
    components = [sources['metric_forcing'][key][0]
                  for key in ('yy', 'y0', '00', 'ab')]
    components.append(sources['scalar_forcing'][0])
    return sum(Q(w) * mag(v)**2 for w, v in zip(GRAM, components))


def full_square(sources, certificates):
    total = Q(0)
    for field, certificate in certificates.items():
        values = (sources['scalar_forcing'] if field == 'scalar'
                  else sources['metric_forcing'][field])
        caps = list(map(mag, values))
        gram = [[Q(v) for v in row] for row in certificate['gram']]
        need(len(gram) == len(caps) and all(len(row) == len(caps) for row in gram),
             'ANGULAR_GRAM_SHAPE')
        total += Q(certificate['multiplicity'])*sum(
            abs(gram[i][j])*caps[i]*caps[j]
            for i in range(len(caps)) for j in range(len(caps)))
    return total


def run(repo):
    pins = authenticate(repo)
    raw = (repo / GEO / 'numeric/ACTUAL_SOURCES_NORMAL.json').read_bytes()
    need(raw == (repo / GEO / 'numeric/ACTUAL_SOURCES_OPTIMIZED.json').read_bytes(),
         'RETAINED_NORMAL_OPTIMIZED_IDENTITY')
    actual = json.loads(raw)
    need(actual['status'] == 'PASS_CONTINUOUS_ACTUAL_EVOLVED_GEOMETRIC_SOURCE_CELL_ENCLOSURES',
         'ACTUAL_CONTINUOUS_SOURCE_STATUS')
    need(actual['domain']['causal_triangle'] == '0<=x<=3,r>=0,x+r<=8'
         and Q(actual['domain']['ell']) == ELL, 'SOURCE_DOMAIN_AND_SCALE')
    need(actual['normalization'] == 'Second parameter derivatives along the affine actual-first-pair path; metric_forcing=-2 ell^2 rho^2 QEred; scalar_forcing=ell^2 rho^2 Qphi. QC is the actual-metric-lowered quadratic gauge Hessian.',
         'ACTUAL_SOURCE_HESSIAN_NORMALIZATION')
    cover = actual['global_cover']
    d = cover['division_per_unit']
    rows = cover['cell_enclosures']
    expected = [(i, k) for i in range(3*d) for k in range(8*d) if i+k < 8*d]
    need([tuple(row['index']) for row in rows] == expected and len(rows) == cover['cells'],
         'COMPLETE_OLD_CAUSAL_COVER')
    for row in rows:
        i, k = row['index']
        need(list(map(Q, row['x'])) == [Q(i, d), Q(i+1, d)]
             and list(map(Q, row['r'])) == [Q(k, d), Q(k+1, d)]
             and row['intersection'] == 'x+r<=8', 'EXACT_CELL_GEOMETRY')
    bound = cover['summary']['positive_bundle_bulk_forcing_H_bound']
    FH = Q(bound['F_H_upper'])
    need(FH == 17387513156573 and FH**2 >= Q(bound['squared_upper']),
         'REGISTERED_FULL_H_FORCING_BOUND')
    join = json.loads((repo / COMMON / 'evolution/SECOND_RESPONSE_FINAL_NORMAL.json').read_bytes())['source_join']
    GH = Q(join['G_H_upper'])
    need(GH == 3851400 and GH == ELL * 131 * (8*Q(15*10**8) + 44*10**8 + 13*10**9),
         'BOUNDARY_DERIVATIVE_SOURCE_JOIN')

    bulk_rows = []
    SF_full = SF_active = full_H_cell_active = total_area = active_area = Q(0)
    for row in rows:
        a, aa = cell_area(row['x'], row['r'], False), cell_area(row['x'], row['r'], True)
        total_area += a
        active_area += aa
        if a:
            sq = mean_square(row['sources'])
            fullsq = full_square(row['sources'], bound['block_Gram_certificates'])
            # Either the directly enclosed invariant slots or the independently
            # certified full angular H cap bounds the projected source.
            chosen = min(sq, fullsq, FH**2)
            SF_full += a*chosen
            SF_active += aa*chosen
            full_H_cell_active += aa*min(fullsq, FH**2)
            bulk_rows.append({'old_cell_index': row['index'], 'cone_area': str(a),
                              'active_area': str(aa), 'direct_mean_square_cap': str(sq),
                              'full_angular_H_square_cap': str(fullsq),
                              'chosen_square_cap': str(chosen)})
    need(total_area == Q(25, 32) and active_area == Q(9, 32), 'EXACT_CONE_AREAS')

    shell_rows = []
    SG_full = SG_active = shell_length = active_length = Q(0)
    for shell in actual['shell_whole_time_cover']['cell_enclosures']:
        xa, xb = map(Q, shell['x'])
        width = max(Q(0), min(xb, T) - max(xa, Q(0)))
        awidth = max(Q(0), min(xb, T) - max(xa, Q(1)))
        shell_length += width
        active_length += awidth
        if not width:
            continue
        matches = [row for row in rows if row['x'] == shell['x'] and row['r'][0] == '0']
        need(len(matches) == 1, 'BOUNDARY_QC_CELL_ALIGNMENT')
        physical = shell['sources']['one_copy_outward_hessian']
        qr00, qrJ, qp = physical['00'][0], physical['ab'][0], physical['scalar'][0]
        qcy = matches[0]['sources']['QC']['y'][0]
        gn = linear_interval([(Q(2), qr00), (Q(-6), qrJ), (Q(-2), qcy)])
        # Five-mean G=( -ell*rho*GN,0,2ell*rho*QR00,
        #              2ell*rho*QRJ,ell*rho*QP ); rho_b <=131.
        sq = (ELL*131)**2*(mag(gn)**2 + 4*mag(qr00)**2 + 12*mag(qrJ)**2 + mag(qp)**2)
        chosen = min(sq, GH**2)
        SG_full += width*chosen
        SG_active += awidth*chosen
        shell_rows.append({'old_time_cell': shell['x'], 'length': str(width),
                           'active_length': str(awidth), 'GN_mean_interval': list(map(str, gn)),
                           'direct_mean_square_cap': str(sq), 'chosen_square_cap': str(chosen)})
    need(shell_length == T and active_length == Q(1,4), 'EXACT_BOUNDARY_LENGTHS')
    Cf, Cg = Q(1,2), Q(5002001,4000000)
    budgets = {}
    for name, f, g in (
            ('full_cone_global_caps', Q(25,32)*FH**2, T*GH**2),
            ('wedge_zero_global_caps', Q(9,32)*FH**2, Q(1,4)*GH**2),
            ('full_cone_retained_cell_caps', SF_full, SG_full),
            ('wedge_zero_retained_cell_caps', SF_active, SG_active)):
        budget = Cf*f + Cg*g
        budgets[name] = {'bulk_L2_squared_upper': str(f), 'boundary_L2_squared_upper': str(g),
                         'energy_source_ledger_upper': str(budget),
                         'conditional_zero_initial_energy_upper': str(44*budget),
                         'conditional_zero_initial_mean_trace_squared_upper': str(176*budget)}
    need(SF_active <= Q(9,32)*FH**2 and SG_active <= Q(1,4)*GH**2,
         'CELL_BUDGET_NO_WORSE_THAN_GLOBAL_BUDGET')
    need(authenticate(repo) == pins, 'SOURCE_AUTHENTICATION_AFTER')
    return {
        'schema': 'PRIVATE_RETAINED_ACTUAL_SOURCE_L2_POSTPROCESS_V1',
        'status': 'PASS_EXACT_ZERO_ORACLE_SOURCE_LEDGER',
        'source_pins_sha256': PINS_SHA,
        'input_files_authenticated': len(pins),
        'retained_actual_source_report_sha256': sha(raw),
        'registered_first_state_binding': 'Same full first Cauchy jet and common-harmonic solution via the pinned Oct31 exact-family and Nov1 prepared gauge derivations; the Oct31 native velocity change is a second-response approximation candidate only.',
        'analytic_identity_premises_reproved': False,
        'mean_order': ['K0','L0','N0','J0','w0'], 'mean_gram': list(GRAM),
        'forcing_convention': 'Second parameter derivatives; S0=P0(-2ell^2*rho^2*QEred,ell^2*rho^2*Qphi), G0=(-ell*rho*GN,0,2ell*rho*QR00,2ell*rho*QRJ,ell*rho*QP).',
        'registered_caps': {'F_H': str(FH), 'G_H': str(GH)},
        'exact_geometry': {'T': str(T), 'cone_area': str(total_area), 'active_area': str(active_area),
                           'boundary_length': str(shell_length), 'active_boundary_length': str(active_length)},
        'energy_coefficients_conditional_on_accompanying_forced_proof':
            {'C': '301/100', 'Cf': str(Cf), 'Cg': str(Cg), 'exp_CT_upper': '44'},
        'budgets': budgets,
        'additional_full_angular_H_cellwise_bulk_L2_squared_upper': str(full_H_cell_active),
        'retained_bulk_cells_intersecting_cone': bulk_rows,
        'retained_shell_cells_intersecting_boundary': shell_rows,
        'scientific_activity': {'project_imports':0,'actual_first_jet_calls':0,'source_oracle_calls':0,
                                'background_queries':0,'PDE_solves':0,'bump_queries':0},
        'scope': [
            'Upper bounds for continuous actual source integrals on the existing future cone, using retained complete-cell enclosures; no new source evaluation or retarded quadrature.',
            'Wedge-zero versions additionally retain the analytic exact support/finite-propagation premise x+r<=1.',
            'Energy and trace consequences require the accompanying forced-energy proof and an independently specified compatible initial energy.',
            'Zero-initial numbers describe an auxiliary reduced response only; they do not declare zero physical second data admissible.',
            'No source certification on x<0, changed first states, or a proposed past-retarded preparation is supplied.',
            'No absolute proper-observable sign follows without the retained first-only clock, measure, and Gauss terms.'
        ]}


def main():
    class Parser(argparse.ArgumentParser):
        def error(self, message):
            raise ValueError('CLI_ARGUMENTS')
    parser = Parser(add_help=False)
    parser.add_argument('--repo-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    need(args.repo_root.exists() and args.repo_root.is_dir(), 'REPOSITORY_REQUIRED')
    need(not any(p.is_symlink() for p in (args.repo_root, *args.repo_root.parents)), 'SYMLINK_PATH')
    repo = args.repo_root.resolve()
    need(args.output.is_absolute(), 'ABSOLUTE_OUTPUT_REQUIRED')
    need(not any(p.is_symlink() for p in (args.output, *args.output.parents)), 'SYMLINK_PATH')
    out = args.output.resolve()
    need(not out.exists(), 'OUTPUT_EXISTS')
    need(not out.is_relative_to(repo) and not out.is_relative_to(HERE), 'OUTPUT_INSIDE_SOURCE')
    need(out.parent.exists(), 'OUTPUT_PARENT_MISSING')
    design, closure = authenticate_portable(repo)
    result = run(repo)
    exclusions = design['scientific_projection_excluded_top_keys']
    projected = {k:v for k,v in result.items() if k not in exclusions}
    projection_sha = sha(json.dumps(projected,sort_keys=True,separators=(',',':')).encode())
    need(projection_sha == design['scientific_projection_sha256'], 'FROZEN_SCIENTIFIC_PROJECTION')
    result['schema'] = 'HDBLAST_INCOMING_PREPARATION_SOURCE_LEDGER_V1'
    result['scope'] = [
        'Exact upper bounds for continuous actual sources on the retained future cone, from existing complete-cell enclosures; no new source evaluation or retarded quadrature.',
        'Wedge-zero versions retain the analytic exact support and finite-propagation premise x+r<=1.',
        'Forced-energy and trace consequences require the accompanying forced proof and a separately specified compatible initial energy.',
        'The zero-initial trace columns concern an auxiliary reduced response, which is not alone a constraint-compatible physical second state.',
        'No source certification on x<0, changed first states, or past-retarded selection is supplied.',
        'The overall trace estimate remains looser than the retained Duhamel bound; no numerical prediction improvement or absolute proper-observable sign is inferred.'
    ]
    result['public_provenance'] = {
        'original_private_canonical_sha256':design['original_private_canonical_sha256'],
        'original_scientific_source_sha256':'db49d85349f35c008e6b5950acde70912fb6c92961089192106f6db48efac354',
        'scientific_projection_excluded_top_keys':exclusions,
        'scientific_projection_sha256':projection_sha,
        'source_closure_inputs_sha256':CLOSURE_SHA,
        'distinct_repository_sources_authenticated':len(closure),
        'selected_packet_join_sha256':PACKET_SHA,
        'selected_packet_pin_comparisons':7,
        'frozen_design_sha256':DESIGN_SHA,
        'local_source_sha256':sha(read_regular(HERE/'source_budget.py')),
        'provenance_paths':'Repository-relative source names and component-local input roles only.'
    }
    need(authenticate_portable(repo) == (design,closure), 'PORTABLE_AUTHENTICATION_AFTER')
    encoded = (json.dumps(result,sort_keys=True,indent=2)+'\n').encode()
    with out.open('xb') as stream:
        stream.write(encoded)
    print(json.dumps({'status':result['status'],'output_sha256':sha(encoded)},sort_keys=True))


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, ArithmeticError, KeyError, TypeError, IndexError, AttributeError, UnicodeError) as exc:
        reason = str(exc) if isinstance(exc,ValueError) else type(exc).__name__
        print(json.dumps({'status':'REJECTED','witness':reason},sort_keys=True))
        raise SystemExit(2)
