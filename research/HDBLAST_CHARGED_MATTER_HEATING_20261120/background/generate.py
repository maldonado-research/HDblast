"""Densely sample the published neutral wall. Only the output time schedule changes."""
import argparse
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import subprocess
import tempfile

from scipy.integrate import quad
from scipy.special import expit

HERE = Path(__file__).resolve().parent
COLUMNS = "t a H psi pi rad sheet e V source_proper_time source_emitted source_returned source_radius K Bmag X wall_mass".split()
ORIGINAL = """    vector<double> times = {0., .001, .002, .003, .004, .005, .006, .008, .01, .02, .05};
    for (int j = 1; j <= static_cast<int>(10 * end); j++)
        times.push_back(j / 10.);"""
DENSE = """    vector<double> times;
    for (int j = 0; j <= 4000; ++j)
        times.push_back(j / 50000.);
    for (int j = 161; j <= static_cast<int>(2000 * end); ++j)
        times.push_back(j / 2000.);"""


def source_constants():
    stop = 3 + math.log(39) / 4 + .0625

    def smooth(x):
        if x <= 0:
            return 0.
        if x >= 1:
            return 1.
        return expit(1/(1-x) - 1/x)

    def power(t):
        w = expit(4*(t-3))
        return w * expit(-4*(t-3))**2 * (smooth(t/.35) * smooth((stop-t)/.125))**2

    integral, _ = quad(power, 0., stop, points=[.35, 3., stop-.125], epsabs=2.e-14, epsrel=4.e-13)
    return stop, math.sqrt(2*integral)


def generate(destination):
    source = (HERE/'forward.cpp').read_bytes()
    expected = (HERE/'SOURCE_SHA256').read_text().split()[0]
    if hashlib.sha256(source).hexdigest() != expected:
        raise ValueError('Published forward solver has changed')
    code = source.decode()
    if code.count(ORIGINAL) != 1:
        raise ValueError('Expected the published reporting schedule exactly once')
    stop, norm = source_constants()
    destination.mkdir(exist_ok=True, parents=True)
    with tempfile.TemporaryDirectory() as tmp:
        scratch = Path(tmp)
        cpp = scratch/'dense.cpp'
        cpp.write_text(code.replace(ORIGINAL, DENSE))
        exe = scratch/'wave'
        subprocess.run(['g++','-O3','-std=c++17',str(cpp),'-o',str(exe)],check=True)
        args = [64,.002,.123,16.,norm,stop,.08,6,.08,13200.,1]
        output = subprocess.check_output([str(exe), *map(str,args)],text=True)
    series, diagnostics = [], {}
    for line in output.splitlines():
        name, *values = line.split()
        if name == '@':
            if len(values) != len(COLUMNS):
                raise ValueError('Incomplete sample')
            series.append(values)
        else:
            diagnostics[name] = float(values[0])
    times = [float(row[0]) for row in series]
    if any(b <= a for a,b in zip(times,times[1:])) or times[0] != 0 or times[-1] < 12:
        raise ValueError('Invalid times')
    if not math.isclose(diagnostics['longest_em90_expanding_efolds'],.12140056,abs_tol=2e-7):
        raise ValueError('Parent solution failed reproduction')
    if not math.isclose(diagnostics['peak_exit_sheet_fraction'],.91973925,abs_tol=2e-7):
        raise ValueError('Unexpected parent energy share')
    with (destination/'trace.csv.gz').open('wb') as raw:
        with gzip.GzipFile(filename='',fileobj=raw,mode='wb',compresslevel=6,mtime=0) as stream:
            with io.TextIOWrapper(stream) as f:
                writer = csv.writer(f,lineterminator='\n')
                writer.writerow(COLUMNS)
                writer.writerows(series)
    selected = ['bounce_tau','bounce_radius','exit_complete_tau','maximum_tau','maximum_radius',
                'first_em90_radius','peak_exit_sheet_fraction','longest_em90_expanding_efolds','friedmann_error']
    summary = dict(pinned_source_sha256=expected, samples=len(times), stop=stop, norm=norm,
                   final_time=times[-1], landmarks={name:diagnostics[name] for name in selected})
    (destination/'TRACE.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=HERE)
    generate(parser.parse_args().output_dir)

