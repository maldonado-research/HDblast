"""Authenticate fixed source snapshots and replay finite nonlinear-limit ledgers."""
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import PurePosixPath
import re

REF = '9b74d57b6caef415775abc069838a4bcea303831'
REPOSITORY = 'https://github.com/maldonado-research/HDblast'
SOURCE_PATHS = {'BOUNDARY_FLOW_DATA.json': 'research/HDBLAST_BOUNDARY_FLOW_PATCH_20261106/DATA.json', 'BOUNDARY_FLOW_MANIFEST.json': 'research/HDBLAST_BOUNDARY_FLOW_PATCH_20261106/MANIFEST.json', 'BOUNDARY_FLOW_PROOF.md': 'research/HDBLAST_BOUNDARY_FLOW_PATCH_20261106/BOUNDARY_FLOW_PROOF.md', 'CAUCHY_DATA.json': 'research/HDBLAST_UNIFORM_CAUCHY_DATA_20261113/DATA.json', 'CAUCHY_MANIFEST.json': 'research/HDBLAST_UNIFORM_CAUCHY_DATA_20261113/MANIFEST.json', 'CONDITIONAL_PICARD_ADMISSION.md': 'research/HDBLAST_VARIABLE_COEFFICIENT_ENERGY_20261114/CONDITIONAL_PICARD_ADMISSION.md', 'CURRENT_CONTINUATION.md': 'research/HDBLAST_VARIABLE_COEFFICIENT_ENERGY_20261114/CONTINUATION.md', 'CURRENT_EVOLUTION_DATA.json': 'research/HDBLAST_VARIABLE_COEFFICIENT_ENERGY_20261114/DATA.json', 'CURRENT_EVOLUTION_MANIFEST.json': 'research/HDBLAST_VARIABLE_COEFFICIENT_ENERGY_20261114/MANIFEST.json', 'CURRENT_EVOLUTION_README.md': 'research/HDBLAST_VARIABLE_COEFFICIENT_ENERGY_20261114/README.md', 'CURRENT_EVOLUTION_SOURCE_PINS.json': 'research/HDBLAST_VARIABLE_COEFFICIENT_ENERGY_20261114/SOURCE_PINS.json', 'CYLINDER_EMBEDDING.md': 'research/HDBLAST_ORDERED_EVOLUTION_20261106/QUADRATIC_FORCING.md', 'FROZEN_ENERGY_COMMUTATORS.md': 'research/HDBLAST_VARIABLE_COEFFICIENT_ENERGY_20261114/FROZEN_ENERGY_COMMUTATORS.md', 'INTERIOR_NORMAL_RECOVERY.md': 'research/HDBLAST_VARIABLE_COEFFICIENT_ENERGY_20261114/INTERIOR_NORMAL_RECOVERY.md', 'METRIC_ENERGY_CERTIFICATE.md': 'research/HDBLAST_BOUNDARY_FLOW_PATCH_20261106/METRIC_ENERGY_CERTIFICATE.md', 'MIXED_INITIAL_JET_PRODUCTS.md': 'research/HDBLAST_TEMPORAL_JETS_SEED_20261114/MIXED_INITIAL_JET_PRODUCTS.md', 'NATURAL_SYMMETRY.md': 'research/HDBLAST_HIGH_ORDER_FIRST_WAVE_20261112/SYMMETRY.md', 'NONLINEAR_INITIAL_JETS.md': 'research/HDBLAST_TEMPORAL_JETS_SEED_20261114/NONLINEAR_INITIAL_JETS.md', 'NORMALIZED_METRIC_MARGINS.md': 'research/HDBLAST_VARIABLE_COEFFICIENT_ENERGY_20261114/NORMALIZED_METRIC_MARGINS.md', 'ORDERED_DATA.json': 'research/HDBLAST_ORDERED_EVOLUTION_20261106/DATA.json', 'ORDERED_MANIFEST.json': 'research/HDBLAST_ORDERED_EVOLUTION_20261106/MANIFEST.json', 'PROPER_GEOMETRY_MANIFEST.json': 'research/HDBLAST_PROPER_SECOND_GEOMETRY_20261113/MANIFEST.json', 'QUANTITATIVE_COMPATIBLE_SEED.md': 'research/HDBLAST_TEMPORAL_JETS_SEED_20261114/QUANTITATIVE_COMPATIBLE_SEED.md', 'QUASILINEAR_VARIABLE_MIXED_EVOLUTION.md': 'research/HDBLAST_QUASILINEAR_SHELL_EVOLUTION_20261015/astra/QUASILINEAR_VARIABLE_MIXED_EVOLUTION.md', 'REVERSE_TIME.md': 'research/HDBLAST_PROPER_SECOND_GEOMETRY_20261113/REVERSE_TIME.md', 'ROTATION_MIXED_PRODUCTS.md': 'research/HDBLAST_VARIABLE_COEFFICIENT_ENERGY_20261114/ROTATION_MIXED_PRODUCTS.md', 'TEMPORAL_SEED_DATA.json': 'research/HDBLAST_TEMPORAL_JETS_SEED_20261114/DATA.json', 'TEMPORAL_SEED_MANIFEST.json': 'research/HDBLAST_TEMPORAL_JETS_SEED_20261114/MANIFEST.json', 'UNIFORM_CAUCHY_DATA.md': 'research/HDBLAST_UNIFORM_CAUCHY_DATA_20261113/UNIFORM_CAUCHY_DATA.md'}
EXPECTED_SOURCE_PINS_SHA256 = 'c0909cd0ff0414560795c78aa08514d424149d3d986537165d01a6efaf40fde6'
MODULES = {'arithmetic': 'convergence_arithmetic', 'oracle': 'convergence_oracles', 'extension': 'extension_oracles'}

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
    need(pins["schema"] == "hdblast-nonlinear-limit-pins-v1",
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
          "HDBLAST_NONLINEAR_PICARD_LIMIT_INPUT_PROVENANCE_V1"
          and provenance["repository"] == REPOSITORY and provenance["ref"] == REF)
    check("provenance_snapshot_only", provenance["snapshot_method"] ==
          "git cat-file of immutable base blobs" and provenance["snapshot_producers_executed"] is False)
    check("provenance_source_index", provenance["source_pins_sha256"] == digest(raw_pins)
          and type(provenance["snapshot_count"]) is int
          and provenance["snapshot_count"] == len(files))
    check("provenance_same_blob_roster", provenance["files"] == [
        {"local": row["local"], "git_blob_sha1": row["git_blob_sha1"], "sha256": row["sha256"]}
        for row in files])

    ancestors = {'CURRENT_EVOLUTION_MANIFEST.json': {'BOUNDARY_FLOW_DATA.json': 'inputs/BOUNDARY_FLOW_DATA.json', 'BOUNDARY_FLOW_MANIFEST.json': 'inputs/BOUNDARY_FLOW_MANIFEST.json', 'BOUNDARY_FLOW_PROOF.md': 'inputs/BOUNDARY_FLOW_PROOF.md', 'CAUCHY_DATA.json': 'inputs/CAUCHY_DATA.json', 'CAUCHY_MANIFEST.json': 'inputs/CAUCHY_MANIFEST.json', 'CONDITIONAL_PICARD_ADMISSION.md': 'CONDITIONAL_PICARD_ADMISSION.md', 'CURRENT_CONTINUATION.md': 'CONTINUATION.md', 'CURRENT_EVOLUTION_DATA.json': 'DATA.json', 'CURRENT_EVOLUTION_README.md': 'README.md', 'CURRENT_EVOLUTION_SOURCE_PINS.json': 'SOURCE_PINS.json', 'CYLINDER_EMBEDDING.md': 'inputs/CYLINDER_EMBEDDING.md', 'FROZEN_ENERGY_COMMUTATORS.md': 'FROZEN_ENERGY_COMMUTATORS.md', 'INTERIOR_NORMAL_RECOVERY.md': 'INTERIOR_NORMAL_RECOVERY.md', 'METRIC_ENERGY_CERTIFICATE.md': 'inputs/METRIC_ENERGY_CERTIFICATE.md', 'MIXED_INITIAL_JET_PRODUCTS.md': 'inputs/MIXED_INITIAL_JET_PRODUCTS.md', 'NATURAL_SYMMETRY.md': 'inputs/NATURAL_SYMMETRY.md', 'NONLINEAR_INITIAL_JETS.md': 'inputs/NONLINEAR_INITIAL_JETS.md', 'NORMALIZED_METRIC_MARGINS.md': 'NORMALIZED_METRIC_MARGINS.md', 'ORDERED_DATA.json': 'inputs/ORDERED_DATA.json', 'ORDERED_MANIFEST.json': 'inputs/ORDERED_MANIFEST.json', 'PROPER_GEOMETRY_MANIFEST.json': 'inputs/PROPER_GEOMETRY_MANIFEST.json', 'QUANTITATIVE_COMPATIBLE_SEED.md': 'inputs/QUANTITATIVE_COMPATIBLE_SEED.md', 'QUASILINEAR_VARIABLE_MIXED_EVOLUTION.md': 'inputs/QUASILINEAR_VARIABLE_MIXED_EVOLUTION.md', 'REVERSE_TIME.md': 'inputs/REVERSE_TIME.md', 'ROTATION_MIXED_PRODUCTS.md': 'ROTATION_MIXED_PRODUCTS.md', 'TEMPORAL_SEED_DATA.json': 'inputs/CURRENT_SEED_DATA.json', 'TEMPORAL_SEED_MANIFEST.json': 'inputs/CURRENT_SEED_MANIFEST.json', 'UNIFORM_CAUCHY_DATA.md': 'inputs/UNIFORM_CAUCHY_DATA.md'}, 'TEMPORAL_SEED_MANIFEST.json': {'TEMPORAL_SEED_DATA.json': 'DATA.json', 'NONLINEAR_INITIAL_JETS.md': 'NONLINEAR_INITIAL_JETS.md', 'MIXED_INITIAL_JET_PRODUCTS.md': 'MIXED_INITIAL_JET_PRODUCTS.md', 'QUANTITATIVE_COMPATIBLE_SEED.md': 'QUANTITATIVE_COMPATIBLE_SEED.md'}, 'CAUCHY_MANIFEST.json': {'CAUCHY_DATA.json': 'DATA.json', 'UNIFORM_CAUCHY_DATA.md': 'UNIFORM_CAUCHY_DATA.md'}, 'ORDERED_MANIFEST.json': {'ORDERED_DATA.json': 'DATA.json', 'CYLINDER_EMBEDDING.md': 'QUADRATIC_FORCING.md'}, 'BOUNDARY_FLOW_MANIFEST.json': {'BOUNDARY_FLOW_DATA.json': 'DATA.json', 'BOUNDARY_FLOW_PROOF.md': 'BOUNDARY_FLOW_PROOF.md', 'METRIC_ENERGY_CERTIFICATE.md': 'METRIC_ENERGY_CERTIFICATE.md'}, 'PROPER_GEOMETRY_MANIFEST.json': {'REVERSE_TIME.md': 'REVERSE_TIME.md'}}
    for manifest, members in ancestors.items():
        previous=load_json((root/"inputs"/manifest).read_bytes())
        need(isinstance(previous,dict) and isinstance(previous.get("files"),dict), "Ancestor manifest schema")
        for local, member in members.items():
            row=by_local[local]
            check("ancestor_member:"+manifest+":"+local, previous["files"].get(member)=={"bytes":row["bytes"],"sha256":row["sha256"]})
    return len(names)


def source_models(root,check):
    current=load_json((root/'inputs/CURRENT_EVOLUTION_DATA.json').read_bytes())
    check('current_evolution_schema',current['schema']=='HDBLAST_VARIABLE_COEFFICIENT_ENERGY_DATA_V1')
    ca=current['arithmetic'];co=current['oracle']
    check('current_accepted_analytic_interfaces',ca['review_status']=='accepted_conditional_analytic_inputs' and ca['full_tensor_derivative_slots'] is True and ca['unknown_ball_required'] is False and ca['inverse_derivatives_required'] is False)
    check('current_all_state_radius_and_coefficient_order',ca['Cstar_exponent']==500000 and ca['coefficient_total_order']==16 and ca['predecessor_radius_exponent']==-12 and F(ca['point_state_radius'])==F(1,1000))
    check('current_complete_source_and_predecessor_conventions',ca['source_C_zero'] is True and co['difference_second_derivative_target']=='predecessor' and F(co['difference_A_factor'])==F(co['difference_B_factor'])==1 and co['energy_source_order']==13 and co['normal_source_order']==12)
    check('current_six_evaluated_outputs',[ca[k] for k in ('Y_initial_exponent','energy_unknown_exponent','energy_predecessor_exponent','normal_Y_exponent','normal_unknown_exponent','normal_predecessor_exponent')]==[30,500101,500100,110,500208,500208])
    check('current_closed_smooth_iteration_admission',ca['picard_C0_exponent']==7003210 and ca['picard_H_exponent']==500212 and ca['iterated_amplitude_power_exponent']==500209 and ca['majorant_size_exponent']==-15 and ca['smooth_approximation_factor']==2)
    check('current_same_elapsed_slab',F(ca['ell'])==F(1,10**6) and F(ca['T_factor'])==F(13,8))
    check('current_metric_margins',F(ca['A_ceiling'])==F(1,19) and F(ca['B_lower'])==F(931,1280) and F(ca['B_upper'])==F(20,19) and F(ca['Bzz_lower'])==F(19,20))
    check('current_uncommuted_energy_constants',ca['energy_production_coefficient']==21 and ca['energy_aggregate_production_coefficient']==16 and ca['energy_lower_denominator']==4 and F(co['energy_spurious_a_time_factor'])==0)
    check('current_inner_boundary_split_and_scalar_principal',co['shell_D_time_trace'] is True and co['shell_N_normal_trace'] is True and F(co['shell_A_normal'])==F(co['shell_B_mixed'])==F(co['shell_principal_fiber_coupling'])==0)
    seed=load_json((root/'inputs/TEMPORAL_SEED_DATA.json').read_bytes())
    check('seed_schema',seed['schema']=='HDBLAST_TEMPORAL_JETS_SEED_DATA_V1')
    sa=seed['arithmetic'];st=seed['temporal'];sh=seed['heat']
    check('same_full_array_coefficient_certificate',sa['all_state_coefficient_exponent']==ca['Cstar_exponent'] and sa['ordinary_coefficient_order']==ca['coefficient_total_order'] and sa['all_state_horizontal_conversion'] is True)
    check('same_completed_nonlinear_jets_and_lift',sa['C_J_exponent']==ca['C_J_exponent']==6503000 and sa['lift_exponent']==ca['C_lift_exponent']==255 and sa['seed_from_data_exponent']==ca['C_seed_exponent']==6503255)
    check('same_seed_interval_and_guard',sa['a_seed_exponent']==ca['seed_amplitude_exponent']==-6507485 and sa['jet_data_guard_exponent']==ca['initial_jet_guard_exponent']==-6503020)
    check('same_seed_size_hierarchy',sa['seed_data_exponent']==ca['seed_data_exponent']==-6503285 and sa['seed_jet_exponent']==ca['seed_I14_exponent']==-285 and sa['seed_size_exponent']==ca['seed_X14_exponent']==-30)
    check('actual_untruncated_jet_order',sa['max_order']==st['max_recursion_order']+2==sh['max_jet_order']==sh['max_time_derivative']==14 and sa['full_angular_spectrum'] is True)
    check('actual_source_normalization_and_horizontal_derivative',F(st['source_sign_from_Fbar'])==-1 and F(st['mixed_principal_factor'])==2 and F(st['horizontal_state_transport_factor'])==-1)
    check('same_global_copy_and_component_norms',sa['I14_is_global_union_norm'] is True and F(sa['seed_union_norm_factor'])==1 and sa['retained_copy_count']==ca['retained_copy_count']==2 and sa['state_dimension']==ca['component_count']==16 and sa['base_dimension']==5)
    check('same_shell_ranks',sa['D_components']==4 and sa['N_components']==12 and sa['D_components']+sa['N_components']==sa['state_dimension'])
    check('same_fixed_collars_and_symmetry_action',sa['initial_pair_fixed_zero_collars'] is True and sa['input_jet_table_Gamma_invariant'] is True and sa['Gamma_natural_full_slot_action'] is True and sa['Gamma_preserves_boundary_split'] is True)
    check('same_zero_state_source',sa['exact_zero_state_source'] is ca['source_C_zero'] is True)
    cauchy=load_json((root/'inputs/CAUCHY_DATA.json').read_bytes())
    check('cauchy_schema',cauchy['schema']=='HDBLAST_UNIFORM_CAUCHY_DATA_V1')
    ia=cauchy['arithmetic'];algebra=cauchy['algebra'];chart=cauchy['chart']
    check('same_original_parameter_and_completed_map',ia['initial_amplitude_exponent']==-6402 and ia['completed_first_derivative_exponent']==sa['initial_first_derivative_exponent']==ca['data_d_slope_exponent']==4200 and ia['retained_copy_count']==sa['retained_copy_count'])
    check('current_half_action_density_convention',F(algebra['K_conversion_factor'])==2 and F(algebra['normal_scalar_conversion_factor'])==1 and F(algebra['reference_source_variation_factor'])==1)
    check('retained_nonzero_state_chart_chain',F(chart['velocity_second_P_first_Q_first_factor'])==2 and F(chart['velocity_second_P_second_Q_factor'])==1 and F(chart['velocity_second_third_differential_factor'])==1)
    ordered=load_json((root/'inputs/ORDERED_DATA.json').read_bytes())
    check('ordered_schema',ordered['schema']=='HDBLAST_ORDERED_EVOLUTION_DATA_V1')
    ob=ordered['bounds'];of=ordered['forcing']
    check('same_ordered_coefficient_and_frame_branches',ob['inverse_aggregate_decade']==sa['q_aggregate_exponent']==17000 and ob['inverse_ordinary_derivative_decade']==sa['q_ordinary_exponent']==18000 and ob['connection_derivative_decade']==sa['connection_exponent']==60 and ob['substituted_ordinary_decade']==sa['first_jet_substitution_exponent']==150200)
    check('same_full_slot_embedding_and_cylinder',of['H3_embedding_decade']==sa['embedding_exponent']==ca['embedding_exponent']==7 and F(of['collar_length'])==6*F(ca['ell']) and ob['state_dimension']==16)
    flow=load_json((root/'inputs/BOUNDARY_FLOW_DATA.json').read_bytes())
    check('flow_schema',flow['schema']=='HDBLAST_BOUNDARY_FLOW_PATCH_DATA_V1')
    bf=flow['flow'];bm=flow['metric']
    check('same_actual_flow_patch',bf['state_coordinates']==16 and bf['jet_total_order']==18 and F(bf['jet_ceiling'])==10**126 and F(bf['inner_initial_radius'])==F(ca['point_state_radius']) and bf['real_time_multiple']==6)
    check('same_actual_time_and_metric_reference',F(bf['ell'])==F(ca['ell']) and F(bm['candidate_time'])==F(ca['T_factor'])*F(ca['ell']) and F(bm['safe_coercivity'])==F(1,18000))
    anchors={
      'FROZEN_ENERGY_COMMUTATORS.md':['P_v u := u_ss-A^i(x,v)D_i u_s-B^ij(x,v)H_ij u=F(x,v,Dv),','a=-q^ss, A^i=2q^si/a, B^ij=q^ij/a,','C=-cbar/a, E=-bbar/a, G=-dbar/a.','U(s)=J_E V(-s)','A^z=B^zA=0'],
      'INTERIOR_NORMAL_RECOVERY.md':['10^110','10^500208'],
      'CONDITIONAL_PICARD_ADMISSION.md':['10^(-10^500209)','R/4','10^6503255'],
      'NONLINEAR_INITIAL_JETS.md':['A=2b/a, B=c/a, C=-cbar/a, E=-bbar/a, G=-dbar/a.','J_0=J_E V(0), J_1=-J_E V_t(0).'],
      'NATURAL_SYMMETRY.md':['V_D=0, D_z V_N=0 at z=0; V=0 at z=Z.','natural action on every old and new derivative slot','dy²+rho(y)²[-dt²+cosh²(t)kappa]'],
      'BOUNDARY_FLOW_PROOF.md':['B(0)=0 also implies F_z(0)=0 exactly','physical perturbations vanish','g^{za}=0'],
      'QUASILINEAR_VARIABLE_MIXED_EVOLUTION.md':['No normal commutation of the boundary conditions is performed.'],
      'REVERSE_TIME.md':['V^-(s,z,Omega)=J V(-s,z,Omega).','(V1^-)_s(0)=-J (V1)_t(0)','S J4 S=J4']}
    for name,phrases in anchors.items():
        text=' '.join((root/'inputs'/name).read_text().split())
        check('retained_text_anchor:'+name,all(' '.join(phrase.split()) in text for phrase in phrases))
    names=('Cstar_exponent','coefficient_total_order','embedding_exponent','point_state_radius','Y_initial_exponent','energy_unknown_exponent','energy_predecessor_exponent','normal_Y_exponent','normal_unknown_exponent','normal_predecessor_exponent','C_J_exponent','C_lift_exponent','C_seed_exponent','initial_jet_guard_exponent','seed_amplitude_exponent','D0_safe_exponent','picard_C0_exponent','picard_H_exponent','smooth_approximation_factor','iterated_amplitude_power_exponent','growth_power_exponent','majorant_combined_linear_exponent','nested_gap_coefficient','majorant_size_exponent','T_factor','energy_lower_denominator','A_ceiling','B_upper','source_C_zero','full_tensor_derivative_slots')
    result={name:ca[name] for name in names}
    result.update({'R_exponent':ca['predecessor_radius_exponent'],'copy_count':ca['retained_copy_count'],
       'data_slope_exponent':ca['data_d_slope_exponent'],'component_count':ca['component_count'],
       'spatial_dimension':sa['base_dimension']-1,'first_jet_directions':sa['base_dimension'],
       'energy_growth_multiplier':ca['energy_aggregate_production_coefficient'],'energy_production_multiplier':ca['energy_production_coefficient'],
       'ell_exponent':-6,'max_order':sa['max_order'],'D_components':sa['D_components'],'N_components':sa['N_components'],
       'global_copy_norms':sa['I14_is_global_union_norm'] and F(sa['seed_union_norm_factor'])==1,
       'snapshot_producers_executed':False})
    return result


def load_module(root,name):
    spec=importlib.util.spec_from_file_location('nonlinear_limit_'+name,root/'checks'/(name+'.py'))
    need(spec is not None and spec.loader is not None,'Evaluator loader unavailable')
    loaded=importlib.util.module_from_spec(spec);spec.loader.exec_module(loaded)
    return loaded


def evaluate(root):
    data=load_json((root/'DATA.json').read_bytes())
    need(isinstance(data,dict) and set(data)=={'schema','scope',*MODULES} and data['schema']=='HDBLAST_NONLINEAR_PICARD_LIMIT_DATA_V1' and isinstance(data['scope'],str) and data['scope'],'Data schema mismatch')
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
        exact[suite]=load_module(root,name).evaluate(data[suite],checker(suite))
        need(isinstance(exact[suite],dict) and exact[suite],'Exact suite report required')
    a=data['arithmetic'];o=data['oracle'];e=data['extension'];source=exact['source_models'];joined=checker('interfaces')
    for name in source:
        if name in ('max_order','D_components','N_components','snapshot_producers_executed'):continue
        need(name in a,'Required authenticated claim missing: '+name)
        joined('same_authenticated_'+name,F(a[name])==F(source[name]) if isinstance(a[name],str) and name!='review_status' else a[name]==source[name])
    joined('full_source_state_derivative_charged',a['FV_group_count']==3 and all(F(o[k])==1 for k in ('F_V_C_factor','F_V_E_factor','F_V_G_factor')))
    joined('full_source_firstjet_derivative_charged',a['FP_cross_count']==2==F(o['F_P_G_factor']) and F(o['F_P_E_factor'])==F(o['elapsed_firstjet_variation_factor'])==1)
    joined('both_principal_differences_positive',a['principal_difference_count']==2 and a['principal_difference_sign']==F(o['difference_A_factor'])==F(o['difference_B_factor'])==1)
    joined('Picard_predecessor_and_uniqueness_second_solution',o['Picard_second_derivative_target']=='preceding_solution' and o['uniqueness_second_derivative_target']=='second_solution')
    joined('seed_starts_with_norm_bound_not_PDE',a['recurrence_start']==1 and a['first_difference_has_PDE'] is o['seed_PDE_assumed'] is False)
    joined('stability_keeps_initial_energy',F(o['nonzero_initial_energy_factor'])==1)
    joined('one_growth_conjugation',F(o['outer_growth_weight_sign'])==-1 and F(o['kernel_growth_difference_factor'])==F(o['inner_predecessor_growth_sign'])==1)
    joined('factorial_control_without_small_time',o['simplex_factorial_mode']=='factorial' and a['requires_small_betaT'] is False and a['extra_amplitude_restriction'] is False)
    joined('actual_finite_witness_scope',o['maximum_iteration']==14 and exact['oracle']['finite_exact']['tail_cases']==55)
    joined('strong_lower_rows_and_weak_top_row',o['lower_mixed_total_order']==o['strong_row_maximum']+2==13 and o['weak_top_row']==12 and o['weak_top_spatial_order']==-1 and o['coefficient_iterate_shift']==-1)
    joined('actual_extension_spatial_and_fiber_arrays',exact['extension']['finite_exact']['spatial_directions']==a['spatial_dimension']==4 and exact['extension']['finite_exact']['fiber_components']==a['component_count']==16)
    joined('actual_extension_D_N_ranks',exact['extension']['finite_exact']['D_components']==source['D_components']==4 and exact['extension']['finite_exact']['N_components']==source['N_components']==12)
    joined('correct_firstorder_wave_gradient_and_symmetrizer',F(e['firstorder_gradient_factor'])==F(e['firstorder_A_factor'])==F(e['characteristic_Bzz_factor'])==1 and e['symmetrizer_gradient_mode']=='B')
    joined('actual_base_reflection_interior_terms',F(e['normal_A_reflection_sign'])==F(e['mixed_B_reflection_sign'])==-1 and F(e['interior_normal_A_factor'])==F(e['interior_mixed_B_factor'])==1)
    joined('mixed_inner_and_odd_outer_field_parity',F(e['inner_D_parity'])==F(e['outer_field_parity'])==-1 and F(e['inner_N_parity'])==1)
    joined('ordinary_fourZ_closure',e['closure_multiple']==4 and F(e['Neumann_twoZ_sign'])==-1)
    joined('odd_coefficient_trace_zero_both_ends',F(e['inner_odd_coefficient_trace'])==F(e['outer_odd_coefficient_trace'])==0)
    joined('ordinary_traces_below_top_N_row',e['D_ordinary_trace_maximum']+1==source['max_order']==14 and e['N_ordinary_trace_maximum']+1==e['D_ordinary_trace_maximum']==13 and e['top_N_ordinary_trace_required'] is False)
    joined('causal_PDE_argument_is_written_premise',a['finite_propagation_is_text_premise'] is True)
    exact['interfaces']={'source_mixed_order':source['max_order'],'D_components':source['D_components'],'N_components':source['N_components'],'ordinary_extension_period':e['closure_multiple']}
    named=[suite+':'+name for suite,entries in suites.items() for name in entries]
    return {'schema':'hdblast-nonlinear-limit-v1','status':'PASS','named_check_count':len(named),'named_checks':named,
       'suite_counts':{suite:len(entries) for suite,entries in suites.items()},'source_count':count,'exact':exact,'scope':data['scope']}
