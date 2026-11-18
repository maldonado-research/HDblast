"""Recompute all 16 wider-band density contrasts and frequency envelopes with standard Python."""
import csv
import io
import json
import math
from pathlib import Path
import zipfile

HERE=Path(__file__).resolve().parent
def close(a,b):
    if not math.isclose(a,b,rel_tol=2e-9,abs_tol=5e-10):
        raise ValueError((a,b))


def read(z,name):
    with z.open(name+'.csv') as f:
        return [{k:float(v) for k,v in row.items()} for row in csv.DictReader(io.TextIOWrapper(f,encoding='ascii'))]


def main():
    report=json.loads((HERE/'BAND_RESULTS.json').read_bytes())
    checked=0
    with zipfile.ZipFile(HERE/'BAND_SAMPLES.zip') as data:
        for group,cases in report.items():
            edge=int(group.split('/K')[1])
            pre=read(data,group+'/pre')
            for arm,r in cases.items():
                samples=pre+read(data,group+'/'+arm)
                keep=[v for v in samples if v['s']<2.1]
                b=samples[len(keep)]
                a=keep[-1]
                w=(2.1-a['s'])/(b['s']-a['s'])
                end={k:v+w*(b[k]-v) for k,v in a.items()}
                keep.append(end)
                maximum=max(math.hypot(edge*math.exp(-v['log_a']),100*(v['phi']-.5)) for v in keep)
                close(maximum,r['frequency_envelope_bound'])
                if maximum>6.49e-6**(-1/3):raise ValueError('Mode screen')
                for key in ('rad_hat','massive_hat','rad_fraction'):close(end[key],r[key])
                close(max(v['rad_fraction'] for v in keep),r['max_rad_fraction_through_cut'])
                close(100*(end['rad_hat']/cases['quiet']['rad_hat']-1),r['contrast_percent'])
                checked+=1
    if checked!=16:raise ValueError('Incomplete control matrix')
    print('Recomputed 16 wider-band native-shell endpoints and all-segment frequency screens.')


if __name__=='__main__':
    main()
