"""Bind completed source-free actual audits into a portable compact manifest."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib
import json
HERE=Path(__file__).resolve().parent
BASE=HERE.parent
REG='33eded898f47a96142f85002b896edd4c294180653f50f2b38e83f22d5a75796'
GO='edded4f67840bc54e95673c5f43b71f6804e373d13211c90ba48ca99f72a8389'

def require(ok,message):
 if not ok:raise ValueError(message)
def pin(path):
 raw=path.read_bytes();return {'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()}
def load(path):return json.loads(path.read_bytes())
def ratio(value):return str(value.numerator)+'/'+str(value.denominator)

def main():
 target=HERE/'FINAL_ACTUAL_REVIEW_RECEIPT.json';manifest=HERE/'FINAL_ACTUAL_REVIEW_MANIFEST.json'
 require(not target.exists() and not manifest.exists(),'actual review already sealed')
 node=load(HERE/'ACTUAL_NORMAL_ALLNODE_REVIEW.json')
 pair=load(HERE/'receipt-audit/ACTUAL_PAIR_RECEIPT_AUDIT.json')
 seal=load(HERE/'receipt-audit/RECEIPT_AUDIT_SEAL.json')
 zero=load(HERE/'SIGNED_INTERVAL_ZERO_EXCLUSIONS_48.json')
 require(node['status']=='PASS_INDEPENDENT_ACTUAL_EXPORTED_ALLNODE_ARITHMETIC','all-node status')
 require(pair['status']=='PASS_INDEPENDENT_ACTUAL_RECEIPTS_AND_EXPORTED_SOURCE_ACCOUNTING' and pair['paired_scientific_identity_checked'] is True,'pair status')
 require(pair['controls']['registration']['sha256']==REG and pair['controls']['GO']['sha256']==GO,'actual control pins')
 stable=pair['runs']['normal']['stable_scientific_files']
 require(stable==pair['runs']['optimized']['stable_scientific_files'] and len(stable)==9,'nine identical scientific files')
 for name,digest in node['input_output_files'].items():require(stable[name]['sha256']==digest,'all-node/custody binding')
 for name,item in seal['artifacts'].items():require(pin(HERE/'receipt-audit'/name)==item,'subreview changed:'+name)
 for mode in ('normal','optimized'):
  controls=load(HERE/('reader-controls-'+mode)/'RESULT.json')
  require(controls['status']=='PASS_INDEPENDENT_ACTUAL_READER_CONTROLS' and controls['checks']==26 and
          controls['reader_sha256']==pin(HERE/'actual_allnode_readback.py')['sha256'],'reader controls')
  run=pair['runs'][mode]
  for name,item in run['file_pins'].items():require(pin(Path(run['directory'])/name)==item,'actual output changed:'+mode+'/'+name)
 require(node['nodes']==49152 and node['comparisons']==98304 and node['prefixes']==24 and
         node['normalized_stored_scalar_components_checked']==589824 and node['exported_scalar_binary80_representability_checks']==688128,'exact readback coverage')
 require(zero['source_exact_review_sha256']==pin(HERE/'ACTUAL_NORMAL_ALLNODE_REVIEW.json')['sha256'] and zero['intervals']==48 and zero['zero_exclusions']==48 and zero['includes_zero']==0,'signed zero exclusions')
 for field in ('physical_source_calls','retained_arrays_opened','production_modules_imported'):require(node[field]==0,'no extra review operations')
 for mode,run in pair['runs'].items():
  require(run['source_accounting']['source_attempts_verified']==128 and run['input_authentication']['selected_array_attempts_verified']==36 and
          run['input_authentication']['authenticated_member_receipts']==76 and run['input_authentication']['planned_real_slots_derived_from_shapes']==688152,'actual operation universe')
 # Copies of the already issued public authorization are provenance artifacts,
 # never a new authorization issued by this audit.
 for source,name,digest in [(BASE/'root/public-readback005/PUBLIC_GO.json','PUBLIC_GO.json',GO),
                           (BASE/'manufactured-registration-005.json','FULL_REGISTRATION.json',REG)]:
  require(pin(source)['sha256']==digest,'external control changed')
  destination=HERE/name;require(not destination.exists(),'provenance copy already exists');destination.write_bytes(source.read_bytes())
 global_bounds={}
 for group,names in [('maxima',('U_L1_upper','W_L1_upper','first_order_wronskian_drift_abs')),
                    ('sums',('weighted_U_L1_upper','weighted_W_L1_upper','R_triangle_upper','P_triangle_upper'))]:
  for name in names:
   maximum=max(Q(row[group][name]) for row in node['computed_finite_aggregates'])
   global_bounds[name]={'exact_outward_upper':ratio(maximum),'attaining_reported_cases':[row['case'] for row in node['computed_finite_aggregates'] if Q(row[group][name])==maximum]}
 selected={}
 for path in HERE.iterdir():
  if path.is_file() and path.name not in (target.name,manifest.name):selected[path.name]=pin(path)
 selected['receipt-audit/RECEIPT_AUDIT_SEAL.json']=pin(HERE/'receipt-audit/RECEIPT_AUDIT_SEAL.json')
 for name,item in seal['artifacts'].items():selected['receipt-audit/'+name]=item
 for mode in ('normal','optimized'):
  name='reader-controls-'+mode+'/RESULT.json';selected[name]=pin(HERE/name)
 proof_names=['provenance_decoder/exact_binary80.py','provenance_decoder/STATIC_PROVENANCE_SUMMARY.md',
 'provenance_decoder/LAYOUT_PROVENANCE_BASIS.json','theory/LATER_TRAJECTORY_THEOREM.md',
 'theory/SOURCE_AND_ENDPOINT_ENGINE_PROOF.md','theory/DYADIC_SOURCE_ERROR_UPDATE.md']
 report={'schema_version':1,'status':'PASS_INDEPENDENT_ACTUAL_ENDPOINT_EXPORTED_CERTIFICATE_REVIEW',
 'scope':'FINITE_SAVED_ENDPOINT_TOTAL_LATER_ERROR_FROM_EXACT_REPRESENTED_INCOMING_ANCHOR',
 'public_freeze_commit':pair['controls']['freeze_commit'],'registration_sha256':REG,'public_GO_sha256':GO,
 'registered_file_pins':pair['controls']['source_files'],
 'independent_allnode_arithmetic_pass':True,'source_budget_journal_input_metadata_custody_pass':True,
 'nine_scientific_files_identical':True,'scientific_file_pins':stable,
 'coverage_per_run':{'capsules':4,'nodes':49152,'endpoint_comparisons':98304,'finite_aggregates':24,
 'source_callbacks':128,'selected_array_decodes':36,'authenticated_members':76,'decoded_real_slots':688152},
 'independent_export_coverage':{'normalized_saved_scalar_components':589824,'binary80_representability_checks':688128,
 'negative_zero_flags_checked':688128,'negative_zero_flags_true':0,'observation_time_values_exported':False,
 'normal_stream_checked_directly':True,'optimized_stream_covered_by_exact_byte_identity':True},
 'target_export_gate':{'exact_L1_radius_gate':'1/1000000000000000000','exact_maximum_complete_radii':node['maximum_radius'],
 'all_targets_within_gate':all(Q(value)<=Q(1,10**18) for value in node['maximum_radius'].values()),
 'gate_applies_to_target_uncertainty_not_measured_saved_discrepancy':True},
 'exact_global_finite_error_upper_bounds':global_bounds,
 'signed_interval_zero_exclusions':zero,
 'exact_aggregate_authority':{'path':'ACTUAL_NORMAL_ALLNODE_REVIEW.json',**pin(HERE/'ACTUAL_NORMAL_ALLNODE_REVIEW.json')},
 'bounded_actual_runs':{mode:{'wall_seconds':run['wall_seconds'],'cpu_seconds_wait4':run['cpu_seconds_wait4'],
 'peak_rss_kib_wait4':run['peak_rss_kib_wait4'],'worker_output':run['worker_output'],'all_output_file_pins':run['file_pins']} for mode,run in pair['runs'].items()},
 'controls_each_python_mode':{'reader_representation_and_export_mutations':26,'source_accounting_and_authorization_mutations':8},
 'audit_operations':{'production_modules_imported':0,'retained_arrays_opened':0,'retained_values_decoded':0,
 'physical_source_callbacks':0,'producer_executions':0,'network_requests':0,'new_GO_issued':False},
 'preserved_history':{'sealed_preflight_artifacts_rehashed_unchanged':184,'source003_source004_rejections_preserved':True,
 'prior_metric_calibration':'FAIL_UNCHANGED','original_unsaved_continuous_solver_trajectory':'UNRECOVERABLE_FROM_RETAINED_MODE_SNAPSHOTS',
 'continuous_momentum_pressure_contact_time_UV_certificate':'UNRESOLVED'},
 'limits':['Original represented values are not independently reauthenticated against archived arrays; duplicate retained decoding was prohibited.',
 'Physical source coefficients, their construction-error premise and target evolution were not reconstructed; containment relies on the frozen proof and implementation.',
 'Exact inverse-epsilon/binary80 representability is a consistency check, not independent archived-value identity.',
 'The 24 observation-time values are checked by the authenticated worker but absent from the node export.',
 'The bound covers total subsequent error at the two saved endpoints; incoming preparation error and unsaved per-step error attribution are separate.',
 'Normal/optimized identity is a reproducibility control, not an independent analytic proof.'],
 'proof_and_decoder_provenance_pins':{name:pair['controls']['source_files'][name] for name in proof_names},
 'compact_review_evidence_pins':selected,'authorizes_new_physical_or_retained_operations':False}
 require(report['target_export_gate']['all_targets_within_gate'],'target gate')
 target.write_text(json.dumps(report,sort_keys=True,indent=2)+'\n');selected[target.name]=pin(target)
 record={'schema_version':1,'status':'SEALED_PORTABLE_COMPACT_ACTUAL_REVIEW','files':selected,
 'file_count':len(selected),'total_bytes':sum(item['bytes'] for item in selected.values()),
 'scope':'Relative paths only; original saved scientific outputs remain external and are bound by exact hashes.',
 'excluded':'Copied mutation fixtures, large node stream/source exports, development receipt-audit versions; preserved locally or in complete saved-run assets.'}
 manifest.write_text(json.dumps(record,sort_keys=True,indent=2)+'\n')
 for name,item in selected.items():require(pin(HERE/name)==item,'sealed compact artifact changed:'+name)
 print(json.dumps({'receipt':pin(target),'manifest':pin(manifest),'files':len(selected),'bytes':record['total_bytes'],'signed_zero_exclusions':48},sort_keys=True))
if __name__=='__main__':main()
