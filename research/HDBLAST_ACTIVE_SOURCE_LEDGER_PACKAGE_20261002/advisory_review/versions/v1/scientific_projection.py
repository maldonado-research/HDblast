"""Stdlib-only explicit science projection; imports no producer or array reader."""
from __future__ import annotations
import copy,hashlib,json,re
from pathlib import Path


def need(ok,message):
 if not ok:raise ValueError(message)


def policy():
 p=Path(__file__).with_name('SCIENTIFIC_PROJECTION_ALLOWLIST.json')
 return json.loads(p.read_text())


def project(result,kind):
 rule=policy()
 selected=rule[kind+'_scientific_allowlist'];ignored=rule[kind+'_excluded_top_keys']
 need(type(result) is dict and set(result)==set(selected)|set(ignored),kind+' top-level field universe changed; update review policy explicitly')
 return {name:result[name] for name in selected}


def scientific_digest(result,kind):
 payload=json.dumps(project(result,kind),sort_keys=True,separators=(',',':'),ensure_ascii=False,allow_nan=False).encode('utf-8')
 return hashlib.sha256(payload).hexdigest()


def provenance(result,kind):
 if kind=='primary':
  need(result['source_unchanged_during_execution'] is True,'Executed primary source was changed')
  return {k:result[k] for k in ('producer_sha256','registration_sha256','freeze_commit','input_manifest_sha256','source_unchanged_during_execution')},None
 if kind=='independent':
  p=copy.deepcopy(result['provenance'])
  need(p['post_run_frozen_verification']=='PASS','Independent frozen source verification absent')
  manifest=p['full_frozen_verification'].pop('package_manifest_sha256')
  return {'producer_sha256':result['producer_sha256'],'input_manifest_sha256':result['input_manifest_sha256'],'provenance':p},manifest
 need(kind=='validation','Unknown result kind')
 return {k:result[k] for k in ('registration_sha256','public_freeze_commit')},None


def compare(original,fresh,kind):
 left,left_manifest=provenance(original,kind);right,right_manifest=provenance(fresh,kind)
 need(left==right,kind+' frozen provenance or producer/input pins differ')
 transition=None
 if kind=='independent' and left_manifest!=right_manifest:
  need(left_manifest is None and isinstance(right_manifest,str) and re.fullmatch('[0-9a-f]{64}',right_manifest) is not None,
       'Only the documented original-preZIP None to packaged-manifest hash transition is permitted')
  transition={'original':None,'fresh':right_manifest,'scope':'Package-manifest protection added; all other frozen provenance must be identical.'}
 lproject,rproject=project(original,kind),project(fresh,kind)
 return {'kind':kind,'science_equal':lproject==rproject,'original_scientific_sha256':scientific_digest(original,kind),
         'fresh_scientific_sha256':scientific_digest(fresh,kind),'stable_frame_equal':stable_frame(original,kind)==stable_frame(fresh,kind),
         'original_stable_full_frame_sha256':stable_diagnostic_digest(original,kind) if kind!='validation' else None,
         'fresh_stable_full_frame_sha256':stable_diagnostic_digest(fresh,kind) if kind!='validation' else None,
         'allowed_package_manifest_transition':transition,
         'ignored_top_level_fields':policy()[kind+'_excluded_top_keys']}


def stable_frame(result,kind):
    """Preserve the full declared frame except exact nondeterministic paths."""
    project(result,kind)  # Exact top-level membership is mandatory.
    value=copy.deepcopy(result)
    if kind in ('primary','independent'):
        time_key='seconds' if kind=='primary' else 'elapsed_seconds'
        expected={time_key,'peak_rss_kib','scope','limit_seconds','limit_peak_rss_kib'}
        need(set(value['resources'])==expected,'Resource field universe changed')
        need(value['resources']['limit_seconds']==900 and value['resources']['limit_peak_rss_kib']==262144,'Resource limits changed')
        need(0<value['resources'][time_key]<=900 and 0<value['resources']['peak_rss_kib']<=262144,'Actual resources outside fixed bounds')
        value['resources'].pop(time_key);value['resources'].pop('peak_rss_kib')
        if kind=='independent':value['provenance']['full_frozen_verification'].pop('package_manifest_sha256')
    else:
        for execution in value['executions']:
            execution.pop('elapsed_seconds');execution.pop('peak_rss_kib')
        value.pop('result_sha256');value.pop('python_optimization')
    return value


def stable_diagnostic_digest(result,route):
    need(route in ('primary','independent'),'Stable diagnostic route must be primary or independent')
    frame={'scientific_hash_policy':'HDBLAST_ACTIVE_SOURCE_SCIENTIFIC_PROJECTION_V1','route':route,'diagnostic':stable_frame(result,route)}
    return hashlib.sha256(json.dumps(frame,sort_keys=True,separators=(',',':'),ensure_ascii=False,allow_nan=False).encode('utf-8')).hexdigest()
