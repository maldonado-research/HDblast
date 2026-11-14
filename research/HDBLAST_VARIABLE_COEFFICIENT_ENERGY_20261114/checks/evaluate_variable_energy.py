"""Authenticate fixed source snapshots and replay finite variable-energy ledgers."""
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import PurePosixPath
import re

REF = 'd27cfd999aa2e57995c08b1228a73d3da7c973d1'
REPOSITORY = 'https://github.com/maldonado-research/HDblast'
SOURCE_PATHS = {'BOUNDARY_FLOW_DATA.json': 'research/HDBLAST_BOUNDARY_FLOW_PATCH_20261106/DATA.json', 'BOUNDARY_FLOW_MANIFEST.json': 'research/HDBLAST_BOUNDARY_FLOW_PATCH_20261106/MANIFEST.json', 'BOUNDARY_FLOW_PROOF.md': 'research/HDBLAST_BOUNDARY_FLOW_PATCH_20261106/BOUNDARY_FLOW_PROOF.md', 'CAUCHY_DATA.json': 'research/HDBLAST_UNIFORM_CAUCHY_DATA_20261113/DATA.json', 'CAUCHY_MANIFEST.json': 'research/HDBLAST_UNIFORM_CAUCHY_DATA_20261113/MANIFEST.json', 'COEFFICIENTS.md': 'research/HDBLAST_EVOLUTION_BOUNDS_20261106/COEFFICIENTS.md', 'CURRENT_CONTINUATION.md': 'research/HDBLAST_TEMPORAL_JETS_SEED_20261114/CONTINUATION.md', 'CURRENT_SEED_DATA.json': 'research/HDBLAST_TEMPORAL_JETS_SEED_20261114/DATA.json', 'CURRENT_SEED_MANIFEST.json': 'research/HDBLAST_TEMPORAL_JETS_SEED_20261114/MANIFEST.json', 'CURRENT_SEED_README.md': 'research/HDBLAST_TEMPORAL_JETS_SEED_20261114/README.md', 'CURRENT_SEED_SOURCE_PINS.json': 'research/HDBLAST_TEMPORAL_JETS_SEED_20261114/SOURCE_PINS.json', 'CYLINDER_EMBEDDING.md': 'research/HDBLAST_ORDERED_EVOLUTION_20261106/QUADRATIC_FORCING.md', 'EVOLUTION_DATA.json': 'research/HDBLAST_EVOLUTION_BOUNDS_20261106/DATA.json', 'EVOLUTION_MANIFEST.json': 'research/HDBLAST_EVOLUTION_BOUNDS_20261106/MANIFEST.json', 'METRIC_ENERGY_CERTIFICATE.md': 'research/HDBLAST_BOUNDARY_FLOW_PATCH_20261106/METRIC_ENERGY_CERTIFICATE.md', 'MIXED_INITIAL_JET_PRODUCTS.md': 'research/HDBLAST_TEMPORAL_JETS_SEED_20261114/MIXED_INITIAL_JET_PRODUCTS.md', 'NATURAL_SYMMETRY.md': 'research/HDBLAST_HIGH_ORDER_FIRST_WAVE_20261112/SYMMETRY.md', 'NONLINEAR_INITIAL_JETS.md': 'research/HDBLAST_TEMPORAL_JETS_SEED_20261114/NONLINEAR_INITIAL_JETS.md', 'ORDERED_BACKGROUND.md': 'research/HDBLAST_ORDERED_EVOLUTION_20261106/ORDERED_BACKGROUND.md', 'ORDERED_DATA.json': 'research/HDBLAST_ORDERED_EVOLUTION_20261106/DATA.json', 'ORDERED_MANIFEST.json': 'research/HDBLAST_ORDERED_EVOLUTION_20261106/MANIFEST.json', 'PROPER_GEOMETRY_MANIFEST.json': 'research/HDBLAST_PROPER_SECOND_GEOMETRY_20261113/MANIFEST.json', 'QUANTITATIVE_COMPATIBLE_SEED.md': 'research/HDBLAST_TEMPORAL_JETS_SEED_20261114/QUANTITATIVE_COMPATIBLE_SEED.md', 'QUASILINEAR_VARIABLE_MIXED_EVOLUTION.md': 'research/HDBLAST_QUASILINEAR_SHELL_EVOLUTION_20261015/astra/QUASILINEAR_VARIABLE_MIXED_EVOLUTION.md', 'REVERSE_TIME.md': 'research/HDBLAST_PROPER_SECOND_GEOMETRY_20261113/REVERSE_TIME.md', 'UNIFORM_CAUCHY_DATA.md': 'research/HDBLAST_UNIFORM_CAUCHY_DATA_20261113/UNIFORM_CAUCHY_DATA.md'}
EXPECTED_SOURCE_PINS_SHA256 = '6073858dfb015b5bd751e9dc94031aaa86bac457299d274f930fbe63fc1a8a62'
MODULES = {'arithmetic': 'energy_arithmetic', 'oracle': 'evolution_oracles'}

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
    need(pins["schema"] == "hdblast-variable-energy-pins-v1",
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
          "HDBLAST_VARIABLE_COEFFICIENT_ENERGY_INPUT_PROVENANCE_V1"
          and provenance["repository"] == REPOSITORY and provenance["ref"] == REF)
    check("provenance_snapshot_only", provenance["snapshot_method"] ==
          "git cat-file of immutable base blobs" and provenance["snapshot_producers_executed"] is False)
    check("provenance_source_index", provenance["source_pins_sha256"] == digest(raw_pins)
          and type(provenance["snapshot_count"]) is int
          and provenance["snapshot_count"] == len(files))
    check("provenance_same_blob_roster", provenance["files"] == [
        {"local": row["local"], "git_blob_sha1": row["git_blob_sha1"], "sha256": row["sha256"]}
        for row in files])

    ancestors = {'CURRENT_SEED_MANIFEST.json': {'CURRENT_SEED_README.md': 'README.md', 'CURRENT_SEED_DATA.json': 'DATA.json', 'CURRENT_SEED_SOURCE_PINS.json': 'SOURCE_PINS.json', 'NONLINEAR_INITIAL_JETS.md': 'NONLINEAR_INITIAL_JETS.md', 'MIXED_INITIAL_JET_PRODUCTS.md': 'MIXED_INITIAL_JET_PRODUCTS.md', 'QUANTITATIVE_COMPATIBLE_SEED.md': 'QUANTITATIVE_COMPATIBLE_SEED.md', 'CURRENT_CONTINUATION.md': 'CONTINUATION.md'}, 'CAUCHY_MANIFEST.json': {'CAUCHY_DATA.json': 'DATA.json', 'UNIFORM_CAUCHY_DATA.md': 'UNIFORM_CAUCHY_DATA.md'}, 'EVOLUTION_MANIFEST.json': {'EVOLUTION_DATA.json': 'DATA.json', 'COEFFICIENTS.md': 'COEFFICIENTS.md'}, 'ORDERED_MANIFEST.json': {'ORDERED_DATA.json': 'DATA.json', 'ORDERED_BACKGROUND.md': 'ORDERED_BACKGROUND.md', 'CYLINDER_EMBEDDING.md': 'QUADRATIC_FORCING.md'}, 'BOUNDARY_FLOW_MANIFEST.json': {'BOUNDARY_FLOW_DATA.json': 'DATA.json', 'BOUNDARY_FLOW_PROOF.md': 'BOUNDARY_FLOW_PROOF.md', 'METRIC_ENERGY_CERTIFICATE.md': 'METRIC_ENERGY_CERTIFICATE.md'}, 'PROPER_GEOMETRY_MANIFEST.json': {'REVERSE_TIME.md': 'REVERSE_TIME.md'}}
    for manifest, members in ancestors.items():
        previous=load_json((root/"inputs"/manifest).read_bytes())
        need(isinstance(previous,dict) and isinstance(previous.get("files"),dict), "Ancestor manifest schema")
        for local, member in members.items():
            row=by_local[local]
            check("ancestor_member:"+local, previous["files"].get(member)=={"bytes":row["bytes"],"sha256":row["sha256"]})
    return len(names)


def source_models(root, check):
    seed=load_json((root/'inputs/CURRENT_SEED_DATA.json').read_bytes())
    check('seed_schema',seed['schema']=='HDBLAST_TEMPORAL_JETS_SEED_DATA_V1')
    sa=seed['arithmetic']; st=seed['temporal']; sh=seed['heat']
    check('accepted_seed_input_scope',sa['review_status']=='accepted_conditional_analytic_inputs' and sa['max_order']==14 and sa['ordinary_coefficient_order']==16)
    check('accepted_all_state_covariant_coefficients',sa['all_state_coefficient_exponent']==500000 and sa['all_state_horizontal_conversion'] is True and sa['coefficient_state_patch_exponent']==-3)
    check('accepted_temporal_source_and_horizontal_signs',F(st['source_sign_from_Fbar'])==-1 and F(st['mixed_principal_factor'])==2 and F(st['horizontal_state_transport_factor'])==-1)
    check('accepted_nonzero_seed_and_linear_lift',sa['retained_nonzero_first_velocity'] is True and sa['bounded_linear_lift'] is True and sa['initial_curve_zero_at_parameter_zero'] is True)
    check('accepted_jet_and_seed_operator_exponents',sa['C_J_exponent']==6503000 and sa['lift_exponent']==255 and sa['seed_from_data_exponent']==6503255)
    check('accepted_unchanged_seed_interval',sa['a_seed_exponent']==-6507485 and sa['initial_first_derivative_exponent']==4200)
    check('accepted_seed_size_hierarchy',sa['seed_data_exponent']==-6503285 and sa['seed_jet_exponent']==-285 and sa['seed_size_exponent']==-30 and sa['seed_point_exponent']==-23)
    check('accepted_seed_full_spectrum_and_split',sa['full_angular_spectrum'] is True and sa['D_components']==4 and sa['N_components']==12 and sa['state_dimension']==16 and sa['retained_copy_count']==2)
    check('accepted_fixed_collars_and_full_slot_symmetry',sa['initial_pair_fixed_zero_collars'] is True and sa['Gamma_natural_full_slot_action'] is True and sa['Gamma_preserves_boundary_split'] is True and sa['input_jet_table_Gamma_invariant'] is True)
    check('accepted_heat_counterterms',sh['heat_exponential_sign']=='-1' and sh['heat_time_power']==2 and sh['counterterms']=='full' and sa['cutoff_after_heat'] is True)
    cauchy=load_json((root/'inputs/CAUCHY_DATA.json').read_bytes())
    check('cauchy_schema',cauchy['schema']=='HDBLAST_UNIFORM_CAUCHY_DATA_V1')
    ca=cauchy['arithmetic']; algebra=cauchy['algebra']; chart=cauchy['chart']
    check('same_original_cauchy_interval',ca['initial_amplitude_exponent']==-6402 and ca['initial_K_exponent']==6400 and ca['completed_first_derivative_exponent']==sa['initial_first_derivative_exponent']==4200)
    check('retained_current_canonical_conversion',F(algebra['K_conversion_factor'])==2 and F(algebra['normal_scalar_conversion_factor'])==1)
    check('retained_nonzero_state_velocity_chain',F(chart['velocity_second_P_first_Q_first_factor'])==2 and F(chart['velocity_second_P_second_Q_factor'])==1 and F(chart['velocity_second_third_differential_factor'])==1)
    evolution=load_json((root/'inputs/EVOLUTION_DATA.json').read_bytes())
    check('evolution_schema',evolution['schema']=='HDBLAST_EVOLUTION_BOUNDS_DATA_V1')
    ec=evolution['coefficients']; ee=evolution['energy']
    check('actual_ordinary_coefficient_inputs',ec['coefficient_jet_order']==sa['ordinary_coefficient_order']==16 and ec['final_coefficient_jet_decade']==sa['raw_lower_exponent']==150000 and ec['inverse_jet_decade']==sa['q_aggregate_exponent']==17000 and ec['flow_derivative_decade']==126)
    check('actual_energy_scale_and_time',F(ee['ell'])==F(1,10**6) and F(ee['candidate_time'])==F(13,8*10**6) and F(ee['coercivity'])==F(1,18000))
    ordered=load_json((root/'inputs/ORDERED_DATA.json').read_bytes())
    check('ordered_schema',ordered['schema']=='HDBLAST_ORDERED_EVOLUTION_DATA_V1')
    ob=ordered['bounds']; of=ordered['forcing']
    check('retained_ordered_coefficients',ob['inverse_aggregate_decade']==ec['inverse_jet_decade']==sa['q_aggregate_exponent']==17000 and ob['inverse_ordinary_derivative_decade']==sa['q_ordinary_exponent']==18000 and ob['connection_derivative_decade']==sa['connection_exponent']==60 and ob['substituted_ordinary_decade']==sa['first_jet_substitution_exponent']==150200)
    check('same_full_tensor_embedding',of['H3_embedding_decade']==sa['embedding_exponent']==7 and ob['state_dimension']==sa['state_dimension']==16 and F(of['collar_length'])==F(6,10**6))
    flow=load_json((root/'inputs/BOUNDARY_FLOW_DATA.json').read_bytes())
    check('flow_schema',flow['schema']=='HDBLAST_BOUNDARY_FLOW_PATCH_DATA_V1')
    bf=flow['flow']; bm=flow['metric']
    check('same_all_state_boundary_flow',bf['state_coordinates']==16 and bf['jet_total_order']==18 and F(bf['jet_ceiling'])==10**126 and F(bf['inner_initial_radius'])==F(1,1000))
    check('same_retained_time_patch',F(bf['ell'])==F(ee['ell']) and bf['real_time_multiple']==6 and F(bm['candidate_time'])==F(ee['candidate_time']))
    check('same_original_metric_margin',F(bm['safe_coercivity'])==F(ee['coercivity']) and F(bm['safe_inverse_upper'])==F(1,16000))
    anchors={
      'NONLINEAR_INITIAL_JETS.md':['A=2b/a, B=c/a, C=-cbar/a, E=-bbar/a, G=-dbar/a.','a(x,U) U_ss=2b^i(x,U)D_iU_s+c^{ij}(x,U)H_ijU-Fbar(x,U,DU),','U(s)=J_E V(-s)','J_0=J_E V(0), J_1=-J_E V_t(0).'],
      'MIXED_INITIAL_JET_PRODUCTS.md':['10^500000','- L_U[Omega_i U].'],
      'QUANTITATIVE_COMPATIBLE_SEED.md':['10^6503000','10^6503255','10^-6507485'],
      'QUASILINEAR_VARIABLE_MIXED_EVOLUTION.md':['No normal commutation of the boundary conditions is performed.','P_a u = a u_tt − 2 b^i D_i u_t − c^{ij} D_iD_j u,','b^z=0, c^{zA}=0 at z=0, A angular.'],
      'CYLINDER_EMBEDDING.md':['E3=10^7'],
      'BOUNDARY_FLOW_PROOF.md':['F_z(0)=0','m+k<=18'],
      'REVERSE_TIME.md':['V^-(s,z,Omega)=J V(-s,z,Omega).','(V1^-)_s(0)=-J (V1)_t(0)','S J4 S=J4'],
      'COEFFICIENTS.md':['q^{ij} ∂_i∂_j X = Q(x,X,∂X).','∂_iΩ_j+Ω_iΩ_j−γ^k_ijΩ_k']}
    for name,phrases in anchors.items():
        text=' '.join((root/'inputs'/name).read_text().split())
        check('retained_text_anchor:'+name,all(' '.join(phrase.split()) in text for phrase in phrases))
    return {'Cstar_exponent':sa['all_state_coefficient_exponent'],'coefficient_total_order':sa['ordinary_coefficient_order'],
       'C_J_exponent':sa['C_J_exponent'],'C_lift_exponent':sa['lift_exponent'],'C_seed_exponent':sa['seed_from_data_exponent'],
       'data_d_slope_exponent':sa['initial_first_derivative_exponent'],'seed_amplitude_exponent':sa['a_seed_exponent'],
       'initial_jet_guard_exponent':sa['jet_data_guard_exponent'],'seed_data_exponent':sa['seed_data_exponent'],
       'seed_I14_exponent':sa['seed_jet_exponent'],'seed_X14_exponent':sa['seed_size_exponent'],
       'retained_copy_count':sa['retained_copy_count'],'component_count':sa['state_dimension'],
       'mixed_order':sa['max_order'],'max_original_tensor_rank':sa['original_rank_ceiling'],
       'source_C_zero':sa['exact_zero_state_source'],
       'embedding_exponent':sa['embedding_exponent'],'ell':str(F(ee['ell'])),'T_factor':str(F(ee['candidate_time'])/F(ee['ell'])),
       'seed_jet_exponent':sa['seed_jet_exponent'],'seed_size_exponent':sa['seed_size_exponent'],
       'seed_point_exponent':sa['seed_point_exponent'],'seed_max_order':sa['max_order'],'copies':sa['retained_copy_count'],
       'D_components':sa['D_components'],'N_components':sa['N_components'],'state_dimension':sa['state_dimension'],
       'snapshot_producers_executed':False}


def load_module(root,name):
    spec=importlib.util.spec_from_file_location('variable_energy_'+name,root/'checks'/(name+'.py'))
    need(spec is not None and spec.loader is not None,'Evaluator loader unavailable')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def evaluate(root):
    data=load_json((root/'DATA.json').read_bytes())
    need(isinstance(data,dict) and set(data)=={'schema','scope',*MODULES} and data['schema']=='HDBLAST_VARIABLE_COEFFICIENT_ENERGY_DATA_V1' and isinstance(data['scope'],str) and data['scope'],'Data schema mismatch')
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
    joined=checker('interfaces');a=data['arithmetic'];o=data['oracle'];source=exact['source_models']
    for field in ('Cstar_exponent','coefficient_total_order','C_J_exponent','C_lift_exponent','C_seed_exponent','data_d_slope_exponent','seed_amplitude_exponent','embedding_exponent',
       'initial_jet_guard_exponent','seed_data_exponent','seed_I14_exponent','seed_X14_exponent',
       'retained_copy_count','component_count','mixed_order','max_original_tensor_rank','source_C_zero'):
        joined('same_authenticated_'+field,a[field]==source[field])
    joined('same_scale_and_time',F(a['ell'])==F(source['ell']) and F(a['T_factor'])==F(source['T_factor'])==F(13,8))
    joined('same_complete_jet_order',o['maximum_field_order']==source['seed_max_order']==14)
    joined('full_source_orders',o['energy_source_order']==13 and o['normal_source_order']==12)
    joined('normalized_energy_without_lapse_mass',F(o['energy_spurious_a_time_factor'])==0)
    joined('retained_interior_mixed_principal',F(o['normal_Az_factor'])==-1 and F(o['normal_Bzy_factor'])==-2 and o['normal_divide_before_differentiate_mode']=='unisolated')
    joined('genuine_predecessor_difference',o['difference_second_derivative_target']=='predecessor' and F(o['difference_A_factor'])==F(o['difference_B_factor'])==1)
    joined('preserved_ordered_word_action',o['commutator_word_mode']=='preserve' and F(o['commutator_sign'])==1)
    exact['interfaces']={'source_jet_order':source['seed_max_order'],'source_jet_size_exponent':source['seed_jet_exponent'],
       'source_seed_size_exponent':source['seed_size_exponent'],'source_seed_point_exponent':source['seed_point_exponent']}
    names=[suite+':'+name for suite,entries in suites.items() for name in entries]
    return {'schema':'hdblast-variable-energy-v1','status':'PASS','named_check_count':len(names),
       'named_checks':names,'suite_counts':{suite:len(entries) for suite,entries in suites.items()},
       'source_count':count,'exact':exact,'scope':data['scope']}
