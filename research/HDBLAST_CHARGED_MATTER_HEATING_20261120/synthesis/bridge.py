"""Convert finite harmonic occupancies into dimensionless wall energy and carrier counts."""
import argparse
import csv
import gzip
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load(path):
    return json.loads(path.read_bytes())


def interpolate_em(trace, target):
    previous = None
    with gzip.open(trace, 'rt') as data:
        for row in csv.DictReader(data):
            point = tuple(float(row[name]) for name in ('t','sheet','e'))
            if previous and previous[0] <= target <= point[0]:
                w = (target-previous[0])/(point[0]-previous[0])
                return [a+w*(b-a) for a,b in zip(previous[1:], point[1:])]
            previous = point
    raise ValueError('Event outside the simulated interval')


def calculate(inputs):
    pairs = load(inputs/'quantum/PAIRS.json')
    event = pairs['events']['zero_crossing']
    radius = event['radius']
    sheet, total = interpolate_em(inputs/'background/trace.csv.gz', event['t'])
    model = pairs['protocols']['zero_crossing']['q1']
    quantum_energy = model['conformal_excitation_energy_finite_band']
    factor = 1/(2*math.pi**2*radius**4)
    energy_ratio = quantum_energy*factor/sheet
    assert math.isclose(sheet/total,event['em_fraction'],rel_tol=1e-7)
    # This is one explicitly named kinematic comparison, not a fitted particle.
    # The data table includes values per S0, so it can be used for any choice of wall scale.
    primary = load(inputs/'coupled/CHARGED.json')['parameters']['primary']
    q, charge_to_mass, rest, initial_radius = 1., primary['Q'], primary['D'], 2.
    mu = q/(math.sqrt(3)*charge_to_mass)
    carrier_number = 2*math.pi**2*initial_radius**3*rest/mu
    # Lowest angular harmonic for the chosen conformal scalar on the initial
    # stationary sphere with zero gauge field, including its curvature gap.
    ground_ratio=math.sqrt(1.+1./(initial_radius*mu)**2)
    allowed_excess=.1
    min_mu=1./(initial_radius*math.sqrt((1+allowed_excess)**2-1))
    min_q=math.sqrt(3)*charge_to_mass*min_mu
    min_S0=min_q**2/(4*math.pi*.001)
    out_bounds = {}
    for name in ('taper_short','taper_long'):
        result = pairs['protocols'][name]['q1']
        mean_pairs = result['total_quanta_finite_band']/2
        prob = result['probability_of_any_pair_finite_band']
        assert 0 < prob <= mean_pairs < .002
        out_bounds[name] = dict(mean_pairs_single_species=mean_pairs,
                                probability_at_least_one_pair_single_species=prob,
                                mean_charged_quanta_three_independent_species=6*mean_pairs)
    return dict(
        matching_example=dict(total_rest_density_in_tension_units_at_a2=rest, effective_charge_to_mass=charge_to_mass,
            quantum_charge=q, mass_in_inverse_curvature_units=mu, total_seed_carriers_per_unit_S0=carrier_number,
            early_K=.001, late_K=13200.04,
            local_early_alpha_times_S0=q**2/(4*math.pi*.001),
            local_late_alpha_times_S0=q**2/(4*math.pi*13200.04),
            initial_free_scalar_ground_energy_over_rest=ground_ratio,
            ground_energy_max_fractional_excess=allowed_excess,
            minimum_quantum_mass_for_near_rest_ground_state=min_mu,
            minimum_charge_at_fixed_classical_Q=min_q,
            minimum_S0_if_initial_local_alpha_below_one_and_near_rest_ground_state=min_S0),
        finite_band=dict(n_max=pairs['parameters']['nmax_exact'],q=q,radius_at_equal_energy_crossing=radius,
            wall_time_at_crossing=event['t'], finite_band_comoving_excitation_energy=quantum_energy,
            single_species_density_in_tension_units_times_S0=quantum_energy*factor,
            single_species_energy_to_sheet_ratio_times_S0=energy_ratio, sheet_density=sheet),
        external_field_readouts=out_bounds)


if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--inputs',type=Path,default=ROOT)
    cli.add_argument('--output',type=Path,required=True)
    args=cli.parse_args()
    args.output.write_text(json.dumps(calculate(args.inputs),indent=2,allow_nan=False)+'\n')
