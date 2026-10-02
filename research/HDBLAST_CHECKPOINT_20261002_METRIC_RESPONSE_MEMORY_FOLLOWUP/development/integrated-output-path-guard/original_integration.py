from pathlib import Path
import json,hashlib,subprocess,datetime,ast
root=Path('/workspace/HDblast/research/HDBLAST_CHECKPOINT_20261002_METRIC_RESPONSE')
root.joinpath('REPRODUCE.md').write_text('''# Reproduce the registered metric calibration

Use Python 3.12 and install `requirements-replay.txt` in an isolated environment. The standalone checkpoint contains all preserved reference inputs; it needs no Git credentials or parent checkout.

```sh
python -m pip install -r requirements-replay.txt
python code/replay_metric.py --output /tmp/hdblast-metric-fresh-replay
```

Choose an unused external output directory. `--plan-only` verifies all source hashes and lists the 32 commands without running calculations. A final package verifies its complete membership as well as the original prospective registration and post-freeze receipt. Physical producers run normally once each; both validators also run under optimized Python. Every command has separate logs and exit evidence, and any failure stops the replay with its outputs preserved.

The output retains independent complex mode archives, direct stresses, complete conservation histories, finite-cutoff comparisons, derivative certificates, saved-data summaries and nine figure artifacts. The conservative UV envelopes cover omitted bands only. Finite quadrature and refinement estimates do not certify total numerical error.

`FULL_REGISTRATION.json` authenticates the prospective bytes. `FREEZE_RECEIPT.json` is added only after the public commit has been checked remotely. Original execution evidence and failures are added under `outputs/` and `evidence/original/`; the final ZIP replay receipt is stored beside this checkpoint in the repository to avoid a package self-reference.

This is a homogeneous conformal metric test at the special plane-wave reference, with fixed physical mass, reference subtraction, scalar and incoming BD state. General lapse, actual shifted-root, bulk/shell and state response remain necessary before coupled initial-data or stability work.
''')
# Capture integrated guards without invoking physical callables.
python='/workspace/hdblast-cloud-setup/venv-frw/bin/python'
env=dict(__import__('os').environ,PYTHONDONTWRITEBYTECODE='1',MPLBACKEND='Agg',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
out=root/'evidence/preparation/integrated';out.mkdir(parents=True,exist_ok=True)
records=[]
for optimized in (False,True):
 suffix='_optimized' if optimized else ''
 commands=[('validator_guards','test_metric_validator_guards.py',[]),('raw_baseline_algebra','verify_raw_metric_baselines.py',[]),('inherited_source_pins','verify_inherited_sources.py',['--output',str(out/('INHERITED_SOURCES'+suffix+'.json'))])]
 for label,script,args in commands:
  cmd=[python,*(['-O'] if optimized else []),str(root/'code'/script),*args]
  start=datetime.datetime.now(datetime.timezone.utc).isoformat()
  r=subprocess.run(cmd,env=env,capture_output=True)
  log=out/(label+suffix+'.stdout.log');log.write_bytes(r.stdout);(out/(label+suffix+'.stderr.log')).write_bytes(r.stderr)
  record={'check':label+suffix,'command':cmd,'started_utc':start,'exit_code':r.returncode,'physical_evaluations':0,'stdout_sha256':hashlib.sha256(r.stdout).hexdigest()}
  records.append(record)
  if r.returncode:raise RuntimeError(record)
(out/'EXECUTION.json').write_text(json.dumps(records,indent=2,sort_keys=True)+'\n')
print(json.dumps({'status':'PASS','integrated_pure_checks':len(records),'physical_evaluations':0}))
