"""Check closure, phase-independent bounds and agreement with the recorded reference."""
import argparse
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def close(a, b, label, rtol=3e-6, atol=1e-5):
    if not math.isclose(a, b, rel_tol=rtol, abs_tol=atol):
        raise ValueError(f'{label}: {a} versus {b}')


def verify(report, reference):
    if set(report) != {'coarse_32','fine_32','fine_48'}:
        raise ValueError('Unexpected grid settings')
    for setting, grid in report.items():
        if set(grid['arms']) != {'quiet','near_positive','near_negative','far_positive'}:
            raise ValueError('Missing history')
        for arm, case in grid['arms'].items():
            other = reference[setting]['arms'][arm]
            cut = case['cut']
            for label in ('early','original','quarter_turn'):
                ledger = case[label]
                if max(abs(ledger[k]) for k in ('exc_defect','raw_defect')) > 1e-4:
                    raise ValueError('Scalar/pressure work residual')
                if abs(ledger['potential_cnumber_defect'])>2e-8 or ledger['wronskian_defect']>1e-8:
                    raise ValueError('Potential or state norm')
                close(ledger['mass_work']+ledger['expansion_work'],ledger['energy_change'],label,atol=1e-4)
                close(ledger['ground_mass_work']+ledger['ground_expansion_work'],ledger['ground_energy_change'],label,atol=1e-4)
                for key in ('mass_work','expansion_work','energy_change','end_exc_energy'):
                    close(ledger[key],other[label][key],key)
            if max(abs(cut[k]) for k in ('max_mode_n_gap','energy_gap_at_rotation','purity_gap')) > 1e-9:
                raise ValueError('Rotation changed spectrum or purity')
            if not cut['force_quarter_turn'] < cut['force_ground'] < cut['force_original']:
                raise ValueError('Rotation lost the ground-relative force reversal')
            close((cut['force_original']+cut['force_quarter_turn'])/2,cut['force_diagonal'],'Force midpoint')
            for key in ('force_original','force_quarter_turn','force_ground','sharp_force_half_range'):
                close(cut[key],other['cut'][key],key)
            env = case['future_envelopes']
            if env['ward_matrix_defect']>2e-6 or env['transfer_symplectic_defect']>2e-8:
                raise ValueError('Quadratic energy ledger or symplectic normalization')
            for key in ('mass_work','expansion_work','energy_change'):
                group = env[key]
                if not group['minimum'] < group['quarter_turn'] < group['actual'] < group['maximum']:
                    raise ValueError('Work supports')
                close(group['center']-group['half_range'],group['minimum'],key)
                close(group['center']+group['half_range'],group['maximum'],key)
                close((group['actual']+group['quarter_turn'])/2,group['center'],key)
                if group['max_direct_replay_error']>2e-4:
                    raise ValueError('Matrix and complex mode transports differ')
                for v in ('minimum','maximum','actual','quarter_turn'):
                    close(group[v],other['future_envelopes'][key][v],key)
            for k in ('original_final_inst_number','rotated_final_inst_number'):
                close(case[k], other[k], k)
        for arm, contrast in grid['pulse_minus_quiet'].items():
            for key, value in contrast.items():
                close(value,reference[setting]['pulse_minus_quiet'][arm][key],key,rtol=1e-3,
                      atol=1e-8 if key.startswith('potential_') else 2e-4)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report', type=Path, nargs='?', default=HERE/'PHASE_WORK.json')
    args = parser.parse_args()
    verify(json.loads(args.report.read_bytes()),json.loads((HERE/'PHASE_WORK.json').read_bytes()))
    print('Verified work balances, equal-number rotations and sharp Gaussian phase ranges on four shell histories.')
