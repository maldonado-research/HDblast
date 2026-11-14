"""Authenticate fixed source snapshots and replay finite temporal/heat ledgers."""
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import PurePosixPath
import re

REF = '2592e90650cb440cab2ce74f726e3dcf9b9d1247'
REPOSITORY = 'https://github.com/maldonado-research/HDblast'
SOURCE_PATHS = {'ACTUAL_POTENTIAL_JET_BOUNDS.json': 'research/HDBLAST_UNIFORM_CAUCHY_DATA_20261113/inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json', 'BOUNDARY_FLOW_DATA.json': 'research/HDBLAST_BOUNDARY_FLOW_PATCH_20261106/DATA.json', 'BOUNDARY_FLOW_MANIFEST.json': 'research/HDBLAST_BOUNDARY_FLOW_PATCH_20261106/MANIFEST.json', 'BOUNDARY_FLOW_PROOF.md': 'research/HDBLAST_BOUNDARY_FLOW_PATCH_20261106/BOUNDARY_FLOW_PROOF.md', 'CANONICAL_CONSTRAINTS.md': 'research/HDBLAST_UNIFORM_CAUCHY_DATA_20261113/inputs/CANONICAL_CONSTRAINTS.md', 'COEFFICIENTS.md': 'research/HDBLAST_EVOLUTION_BOUNDS_20261106/COEFFICIENTS.md', 'CURRENT_CAUCHY_DATA.json': 'research/HDBLAST_UNIFORM_CAUCHY_DATA_20261113/DATA.json', 'CURRENT_CAUCHY_MANIFEST.json': 'research/HDBLAST_UNIFORM_CAUCHY_DATA_20261113/MANIFEST.json', 'CURRENT_CAUCHY_README.md': 'research/HDBLAST_UNIFORM_CAUCHY_DATA_20261113/README.md', 'CURRENT_CAUCHY_SOURCE_PINS.json': 'research/HDBLAST_UNIFORM_CAUCHY_DATA_20261113/SOURCE_PINS.json', 'CURRENT_CONTINUATION.md': 'research/HDBLAST_UNIFORM_CAUCHY_DATA_20261113/CONTINUATION.md', 'CYLINDER_EMBEDDING.md': 'research/HDBLAST_ORDERED_EVOLUTION_20261106/QUADRATIC_FORCING.md', 'EVOLUTION_DATA.json': 'research/HDBLAST_EVOLUTION_BOUNDS_20261106/DATA.json', 'EVOLUTION_MANIFEST.json': 'research/HDBLAST_EVOLUTION_BOUNDS_20261106/MANIFEST.json', 'METRIC_ENERGY_CERTIFICATE.md': 'research/HDBLAST_BOUNDARY_FLOW_PATCH_20261106/METRIC_ENERGY_CERTIFICATE.md', 'ORDERED_BACKGROUND.md': 'research/HDBLAST_ORDERED_EVOLUTION_20261106/ORDERED_BACKGROUND.md', 'ORDERED_DATA.json': 'research/HDBLAST_ORDERED_EVOLUTION_20261106/DATA.json', 'ORDERED_MANIFEST.json': 'research/HDBLAST_ORDERED_EVOLUTION_20261106/MANIFEST.json', 'PROPER_GEOMETRY_MANIFEST.json': 'research/HDBLAST_PROPER_SECOND_GEOMETRY_20261113/MANIFEST.json', 'QUASILINEAR_VARIABLE_MIXED_EVOLUTION.md': 'research/HDBLAST_QUASILINEAR_SHELL_EVOLUTION_20261015/astra/QUASILINEAR_VARIABLE_MIXED_EVOLUTION.md', 'REVERSE_TIME.md': 'research/HDBLAST_PROPER_SECOND_GEOMETRY_20261113/REVERSE_TIME.md', 'SPHERE_NORM_BRIDGE.md': 'research/HDBLAST_ALL_MODE_INVERSE_20261113/WEIGHTED_NORM_BRIDGE.md', 'UNIFORM_CAUCHY_DATA.md': 'research/HDBLAST_UNIFORM_CAUCHY_DATA_20261113/UNIFORM_CAUCHY_DATA.md'}
EXPECTED_SOURCE_PINS_SHA256 = 'f43aa808cc57e1627188dc5187dc54539635fced9c64171c5e19658103bb2dfe'
MODULES = {'arithmetic': 'initialization_arithmetic', 'temporal': 'temporal_recursion', 'heat': 'heat_jet_oracle'}

def need(condition, message):
    if not condition:
        raise ValueError(message)

def load_json(raw):
    def unique(items):
        result = {}
        for key, value in items:
            need(key not in result, "Duplicate JSON key: " + key)
            result[key] = value
        return result

    def nonfinite(value):
        raise ValueError("Nonfinite JSON constant: " + value)

    return json.loads(raw.decode("utf-8"), object_pairs_hook=unique,
                      parse_float=str, parse_constant=nonfinite)

def digest(raw):
    return hashlib.sha256(raw).hexdigest()

def public_sources(root, check):
    raw_pins = (root / "SOURCE_PINS.json").read_bytes()
    check("accepted_source_index_sha256", digest(raw_pins) == EXPECTED_SOURCE_PINS_SHA256)
    pins = load_json(raw_pins)
    need(isinstance(pins, dict) and set(pins) == {"schema", "repository", "ref", "files"},
         "Unexpected public pin keys")
    need(pins["schema"] == "hdblast-temporal-seed-pins-v1",
         "Public pin schema mismatch")
    need(pins["repository"] == REPOSITORY and pins["ref"] == REF,
         "Public source repository/ref mismatch")
    files = pins["files"]
    need(isinstance(files, list) and len(files) == len(SOURCE_PATHS),
         "Original public input count differs")
    names = set()
    by_local = {}
    for row in files:
        need(isinstance(row, dict) and set(row) ==
             {"local", "path", "ref", "url", "bytes", "sha256", "git_blob_sha1"},
             "Unexpected source pin keys")
        need(isinstance(row["local"], str) and isinstance(row["path"], str),
             "Source paths must be strings")
        local, source = PurePosixPath(row["local"]), PurePosixPath(row["path"])
        need(not local.is_absolute() and len(local.parts) == 2 and local.parts[0] == "inputs"
             and local.name not in names and "\\" not in row["local"]
             and local.as_posix() == row["local"] and ".." not in local.parts,
             "Unsafe or repeated local source pin")
        need(local.name in SOURCE_PATHS and row["path"] == SOURCE_PATHS[local.name],
             "Original source roster differs")
        need(not source.is_absolute() and source.parts[0] == "research"
             and "\\" not in row["path"] and ".." not in source.parts
             and source.as_posix() == row["path"], "Unsafe original source path")
        need(row["ref"] == REF and row["url"] ==
             REPOSITORY + "/blob/" + REF + "/" + source.as_posix(),
             "Immutable public URL mismatch")
        raw = (root / local).read_bytes()
        check("source_size:" + local.name, type(row["bytes"]) is int
              and row["bytes"] >= 0 and len(raw) == row["bytes"])
        check("source_sha256:" + local.name, isinstance(row["sha256"], str)
              and re.fullmatch(r"[0-9a-f]{64}", row["sha256"]) is not None
              and digest(raw) == row["sha256"])
        blob = b"blob " + str(len(raw)).encode("ascii") + b"\0" + raw
        check("source_git_blob:" + local.name, isinstance(row["git_blob_sha1"], str)
              and re.fullmatch(r"[0-9a-f]{40}", row["git_blob_sha1"]) is not None
              and hashlib.sha1(blob).hexdigest() == row["git_blob_sha1"])
        names.add(local.name)
        by_local[local.name] = row
    need({item.name for item in (root / "inputs").iterdir()} == names,
         "Public input roster mismatch")
    need([row["local"] for row in files] == sorted(row["local"] for row in files),
         "Public input order differs")

    provenance = load_json((root / "INPUT_PROVENANCE.json").read_bytes())
    need(isinstance(provenance, dict) and set(provenance) ==
         {"schema", "repository", "ref", "snapshot_method", "snapshot_producers_executed",
          "source_pins_sha256", "snapshot_count", "files"}, "Input provenance keys differ")
    check("provenance_fixed_base", provenance["schema"] ==
          "HDBLAST_TEMPORAL_JETS_SEED_INPUT_PROVENANCE_V1"
          and provenance["repository"] == REPOSITORY and provenance["ref"] == REF)
    check("provenance_snapshot_only", provenance["snapshot_method"] ==
          "git cat-file of immutable base blobs" and provenance["snapshot_producers_executed"] is False)
    check("provenance_source_index", provenance["source_pins_sha256"] == digest(raw_pins)
          and type(provenance["snapshot_count"]) is int
          and provenance["snapshot_count"] == len(files))
    check("provenance_same_blob_roster", provenance["files"] == [
        {"local": row["local"], "git_blob_sha1": row["git_blob_sha1"], "sha256": row["sha256"]}
        for row in files])

    ancestors = {'CURRENT_CAUCHY_MANIFEST.json': {'CURRENT_CAUCHY_README.md': 'README.md', 'CURRENT_CAUCHY_DATA.json': 'DATA.json', 'CURRENT_CAUCHY_SOURCE_PINS.json': 'SOURCE_PINS.json', 'UNIFORM_CAUCHY_DATA.md': 'UNIFORM_CAUCHY_DATA.md', 'CURRENT_CONTINUATION.md': 'CONTINUATION.md', 'ACTUAL_POTENTIAL_JET_BOUNDS.json': 'inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json', 'CANONICAL_CONSTRAINTS.md': 'inputs/CANONICAL_CONSTRAINTS.md'}, 'EVOLUTION_MANIFEST.json': {'EVOLUTION_DATA.json': 'DATA.json', 'COEFFICIENTS.md': 'COEFFICIENTS.md'}, 'ORDERED_MANIFEST.json': {'ORDERED_DATA.json': 'DATA.json', 'ORDERED_BACKGROUND.md': 'ORDERED_BACKGROUND.md', 'CYLINDER_EMBEDDING.md': 'QUADRATIC_FORCING.md'}, 'BOUNDARY_FLOW_MANIFEST.json': {'BOUNDARY_FLOW_DATA.json': 'DATA.json', 'BOUNDARY_FLOW_PROOF.md': 'BOUNDARY_FLOW_PROOF.md', 'METRIC_ENERGY_CERTIFICATE.md': 'METRIC_ENERGY_CERTIFICATE.md'}, 'PROPER_GEOMETRY_MANIFEST.json': {'REVERSE_TIME.md': 'REVERSE_TIME.md'}}
    for manifest, members in ancestors.items():
        previous=load_json((root/"inputs"/manifest).read_bytes())
        need(isinstance(previous,dict) and isinstance(previous.get("files"),dict), "Ancestor manifest schema")
        for local, member in members.items():
            row=by_local[local]
            check("ancestor_member:"+local, previous["files"].get(member)=={"bytes":row["bytes"],"sha256":row["sha256"]})
    return len(names)


def source_models(root, check):
    cauchy=load_json((root/'inputs/CURRENT_CAUCHY_DATA.json').read_bytes())
    check('current_cauchy_schema',cauchy['schema']=='HDBLAST_UNIFORM_CAUCHY_DATA_V1')
    ca=cauchy['arithmetic']; algebra=cauchy['algebra']; chart=cauchy['chart']
    check('retained_original_initial_interval',ca['initial_amplitude_exponent']==-6402 and ca['initial_K_exponent']==6400)
    check('retained_uniform_initial_pair',ca['completed_first_derivative_exponent']==4200 and ca['completed_Hessian_exponent']==10400 and ca['direct_initial_size_exponent']==-2200)
    check('retained_initial_topologies',ca['flow_total_order']==18 and ca['retained_copy_count']==2)
    check('retained_half_action_density',F(algebra['K_conversion_factor'])==2 and F(algebra['normal_scalar_conversion_factor'])==1)
    check('retained_nonzero_state_chain',F(chart['velocity_second_P_first_Q_first_factor'])==2 and F(chart['velocity_second_P_second_Q_factor'])==1 and F(chart['velocity_second_third_differential_factor'])==1)
    check('retained_reference_source_variation',F(algebra['reference_source_variation_factor'])==1)
    evolution=load_json((root/'inputs/EVOLUTION_DATA.json').read_bytes())
    check('evolution_schema',evolution['schema']=='HDBLAST_EVOLUTION_BOUNDS_DATA_V1')
    ec=evolution['coefficients']; ee=evolution['energy']
    check('actual_ordinary_coefficient_inputs',ec['coefficient_jet_order']==16 and ec['final_coefficient_jet_decade']==150000 and ec['inverse_jet_decade']==17000 and ec['flow_derivative_decade']==126)
    check('actual_all_state_energy_member',F(ee['ell'])==F(1,10**6) and F(ee['candidate_time'])==F(13,8*10**6) and F(ee['coercivity'])==F(1,18000))
    ordered=load_json((root/'inputs/ORDERED_DATA.json').read_bytes())
    check('ordered_schema',ordered['schema']=='HDBLAST_ORDERED_EVOLUTION_DATA_V1')
    ob=ordered['bounds']; of=ordered['forcing']
    check('retained_sharp_principal_branch',ob['inverse_aggregate_decade']==ec['inverse_jet_decade']==17000 and ob['inverse_ordinary_derivative_decade']==18000)
    check('retained_connection_and_substitution',ob['connection_derivative_decade']==60 and ob['substituted_ordinary_decade']==150200 and ob['state_dimension']==16)
    check('retained_full_tensor_embedding',of['H3_embedding_decade']==7 and F(of['collar_length'])==F(6,10**6))
    flow=load_json((root/'inputs/BOUNDARY_FLOW_DATA.json').read_bytes())
    check('flow_schema',flow['schema']=='HDBLAST_BOUNDARY_FLOW_PATCH_DATA_V1')
    bf=flow['flow']; bm=flow['metric']
    check('same_sixteen_state_flow',bf['state_coordinates']==16 and bf['jet_total_order']==18 and F(bf['jet_ceiling'])==10**126)
    check('same_flow_real_patch',F(bf['ell'])==F(1,10**6) and F(bf['inner_initial_radius'])==F(1,1000) and bf['real_time_multiple']==6)
    check('same_all_state_principal_margin',F(bm['safe_coercivity'])==F(ee['coercivity']) and F(bm['safe_inverse_upper'])==F(1,16000) and F(bm['candidate_time'])==F(ee['candidate_time']))
    actual=load_json((root/'inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json').read_bytes())
    potential=[F(v) for v in actual['actual_model']['U_coefficients']]
    check('same_actual_polynomial_model',potential==[F(v) for v in ca['potential_coefficients']]==[F(-2,27),F(0),F(14,9),F(50,27),F(-1,6),F(-4,9),F(-2,27)])
    check('same_authenticated_background', [F(v) for v in actual['background_intervals']['rho']['exact']]==[129,131] and [F(v) for v in actual['background_intervals']['p']['exact']]==[F(-7,50000),F(-3,25000)])
    anchors={
      'UNIFORM_CAUCHY_DATA.md':['10^4200','10^10400','10^-6402','full nonzero-velocity seed'],
      'COEFFICIENTS.md':['q^{ij} ∂_i∂_j X = Q(x,X,∂X).','Ψ_VV(P_i,P_j)','∂_iΩ_j+Ω_iΩ_j−γ^k_ijΩ_k'],
      'QUASILINEAR_VARIABLE_MIXED_EVOLUTION.md':['a u_tt − 2 b^i D_i u_t − c^{ij} D_iD_j u','No normal commutation of the boundary conditions is performed.','The exact nonlinear initial recursion defines J_j'],
      'CYLINDER_EMBEDDING.md':['E3=10^7'],
      'BOUNDARY_FLOW_PROOF.md':['F_z(0)=0','m+k<=18','It is not the physical bulk radial ODE'],
      'METRIC_ENERGY_CERTIFICATE.md':['a=−g^{tt}','1/18000','full spatial inverse block'],
      'REVERSE_TIME.md':['V^-(s,z,Omega)=J V(-s,z,Omega).','(V1^-)_s(0)=-J (V1)_t(0)','S J4 S=J4'],
      'CANONICAL_CONSTRAINTS.md':['H   = R/2 + (2/det γ)','M_i = −2 D_j π^j_i']}
    for name,phrases in anchors.items():
        text=' '.join((root/'inputs'/name).read_text().split())
        check('retained_text_anchor:'+name,all(' '.join(phrase.split()) in text for phrase in phrases))
    return {'initial_first_derivative_exponent':ca['completed_first_derivative_exponent'],
      'old_initial_amplitude_exponent':ca['initial_amplitude_exponent'],'old_initial_size_exponent':ca['direct_initial_size_exponent'],
      'raw_lower_exponent':ec['final_coefficient_jet_decade'],'q_aggregate_exponent':ec['inverse_jet_decade'],
      'q_ordinary_exponent':ob['inverse_ordinary_derivative_decade'],'ordinary_coefficient_order':ec['coefficient_jet_order'],
      'connection_exponent':ob['connection_derivative_decade'],'first_jet_substitution_exponent':ob['substituted_ordinary_decade'],
      'state_dimension':ob['state_dimension'],'embedding_exponent':of['H3_embedding_decade'],'ell':str(F(ee['ell'])),
      'cylinder_length':str(F(of['collar_length'])),'flow_jet_order':bf['jet_total_order'],'flow_jet_exponent':126,
      'all_state_a_lower':str(F(ee['coercivity'])),'candidate_time':str(F(ee['candidate_time'])),
      'snapshot_producers_executed':False}


def load_module(root,name):
    spec=importlib.util.spec_from_file_location('temporal_seed_'+name,root/'checks'/(name+'.py'))
    need(spec is not None and spec.loader is not None,'Evaluator loader unavailable')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def evaluate(root):
    data=load_json((root/'DATA.json').read_bytes())
    need(isinstance(data,dict) and set(data)=={'schema','scope',*MODULES} and data['schema']=='HDBLAST_TEMPORAL_JETS_SEED_DATA_V1' and isinstance(data['scope'],str) and data['scope'],'Data schema mismatch')
    suites={name:[] for name in ('sources','source_models',*MODULES,'interfaces')}
    def checker(suite):
        def check(name,condition):
            need(isinstance(name,str) and name,'Nonempty check name required')
            need(bool(condition),'Arithmetic failed: '+suite+':'+name)
            need(name not in suites[suite],'Duplicate check: '+suite+':'+name)
            suites[suite].append(name)
        return check
    count=public_sources(root,checker('sources'))
    exact={'source_models':source_models(root,checker('source_models'))}
    for suite,name in MODULES.items():
        need(isinstance(data[suite],dict),'Claims must be objects: '+suite)
        module=load_module(root,name)
        exact[suite]=module.evaluate(data[suite],checker(suite))
        need(isinstance(exact[suite],dict) and exact[suite],'Exact suite report required')
    joined=checker('interfaces'); a=data['arithmetic']; source=exact['source_models']
    for field in ('initial_first_derivative_exponent','old_initial_amplitude_exponent','old_initial_size_exponent',
       'raw_lower_exponent','q_aggregate_exponent','q_ordinary_exponent','ordinary_coefficient_order',
       'connection_exponent','first_jet_substitution_exponent','state_dimension','embedding_exponent'):
        joined('same_authenticated_'+field,a[field]==source[field])
    joined('same_scale_and_cylinder',F(a['ell'])==F(source['ell']) and a['cylinder_length_in_ell']*F(a['ell'])==F(source['cylinder_length']))
    joined('actual_all_state_lapse_margin',F(a['a_reciprocal_zeroth_ceiling'])>=1/F(source['all_state_a_lower']))
    joined('flow_derivative_reserve',source['flow_jet_order']==18 and source['flow_jet_exponent']==126 and a['ordinary_coefficient_order']<=source['flow_jet_order'])
    joined('unchanged_time_patch',F(source['candidate_time'])==F(13,8)*F(a['ell']))
    joined('temporal_and_heat_orders',data['temporal']['max_recursion_order']==a['max_order']-2 and data['heat']['max_jet_order']==data['heat']['max_time_derivative']==a['max_order'] and a['time_derivative_count']==a['max_order']+1)
    joined('actual_negative_source',F(data['temporal']['source_sign_from_Fbar'])==-1)
    joined('horizontal_state_connection_retained',F(data['temporal']['horizontal_state_transport_factor'])==-1 and a['all_state_horizontal_conversion'] is True)
    joined('genuine_forward_heat_lift',F(data['heat']['heat_exponential_sign'])==-1 and data['heat']['heat_time_power']==2 and data['heat']['counterterms']=='full')
    joined('distinct_lift_and_data_norm_inputs',a['seed_from_data_exponent']==a['C_J_exponent']+a['lift_exponent'])
    joined('finite_Gaussian_degree_and_spectral_scope',exact['heat']['finite_exact']['maximum_polynomial_degree']==28 and exact['heat']['finite_exact']['arbitrary_jet_count']==15 and exact['heat']['finite_exact']['spectral_weight_identities']==960)
    exact['interfaces']={'C_J_exponent':a['C_J_exponent'],'C_lift_exponent':a['lift_exponent'],
        'C_seed_from_Cauchy_pair_exponent':a['seed_from_data_exponent'],'a_seed_exponent':a['a_seed_exponent'],
        'seed_size_exponent':a['seed_size_exponent'],'seed_point_exponent':a['seed_point_exponent']}
    names=[suite+':'+name for suite,entries in suites.items() for name in entries]
    return {'schema':'hdblast-temporal-jets-seed-v1','status':'PASS','named_check_count':len(names),
       'named_checks':names,'suite_counts':{suite:len(entries) for suite,entries in suites.items()},
       'source_count':count,'exact':exact,'scope':data['scope']}
