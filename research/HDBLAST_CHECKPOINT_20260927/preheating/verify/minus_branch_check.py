#!/usr/bin/env python3
"""Audit script 6: where does the archived minus-branch trajectory reverse (h crosses 0)?  Writes MINUS_BRANCH_CHECK.json."""
import io, json, zipfile
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
Z = ('/home/user/unified-theory-maldonado/new-files/latest-work/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922/'
     'HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922/source_audit/inputs/CHAT14_COMPLETED_REFERENCE.zip')
pre = 'HDBLAST_CHAT14_REGISTERED_DETUNING_ROLLOFF_20260922/runs/v2_t1e3_minus_c4_finecoarse_timeseries.npz'
rec = np.load(io.BytesIO(zipfile.ZipFile(Z).read(pre)))['rec']
t, phi, h, s = rec[:, 0], rec[:, 1], rec[:, 2], rec[:, 9]
i0 = int(np.argmax(h < 0)) if np.any(h < 0) else None
i3 = int(np.argmax(h < 0.3))
ic = int(np.argmax(phi < -0.5))
out = {'s_h_below_0.3': float(s[i3]), 't_h_below_0.3': float(t[i3]),
       's_h_first_negative': float(s[i0]) if i0 is not None else None, 't_h_first_negative': float(t[i0]) if i0 is not None else None,
       's_phi_crosses_-0.5_nodes': [float(s[ic-1]), float(s[ic])], 'data_end_s': float(s[-1]), 'data_end_h': float(h[-1])}
(HERE/'MINUS_BRANCH_CHECK.json').write_text(json.dumps(out, indent=1) + '\n'); print(out)
