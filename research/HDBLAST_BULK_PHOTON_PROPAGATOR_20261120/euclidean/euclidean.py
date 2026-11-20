"""Static transverse Maxwell current response on a tuned HDblast half-leaf."""
import argparse
import json
import math
from pathlib import Path

from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq


THROAT_SLOPE = 5/9
PHI = (0., .5, -.5)
MOMENTA = (.001, .003, .01, .03, .1, .3, 1., 3., 10., 32., 128., 512.)


def primitive(y):
    logcosh = abs(y) + math.log1p(math.exp(-2*abs(y))) - math.log(2.)
    return y + 2*logcosh/3 + math.tanh(y)**2/6


class Leaf:
    def __init__(self, phi):
        self.b = -math.atanh(phi)
        self.warp_offset = primitive(self.b)
        self.norm = (1+phi)**2 * (2-phi)/3
        self.shell_weight = (1-phi**2)**2

    def warp(self, y):
        return (self.warp_offset - primitive(y))/3

    def weight(self, y):
        t = math.exp(-2*abs(y))
        return (4*t/(1+t)**2)**2

    def slopes(self, y):
        t = math.tanh(y)
        return -(1+t-t**3/3)/3, -4*t

    def optical_end(self, q, depth):
        # A fixed dimensionless optical depth in the asymptotic region suffices.
        # Dirichlet/Neumann comparison below tests the whole remaining tail.
        def distance(y):
            value = quad(lambda t: math.exp(-self.warp(t)), self.b, y)[0]
            return q*value-depth
        upper = self.b + 1.
        while distance(upper) < 0:
            upper = self.b + 2*(upper-self.b)
        return brentq(distance, self.b, upper, xtol=1e-14)

    def response(self, q, depth=28, rtol=3e-12):
        end = self.optical_end(q, depth)

        def rhs(y, state):
            warp_slope, log_weight_slope = self.slopes(y)
            k = q*math.exp(-self.warp(y))
            r = state[0]
            return [k*(r*r-1) - (warp_slope+log_weight_slope)*r]

        solution = solve_ivp(rhs, (end,self.b), [1.], method='DOP853',
                             rtol=rtol, atol=rtol*1e-4*q/(1+q))
        if not solution.success:
            raise ArithmeticError(solution.message)
        return self.shell_weight*q*solution.y[0,-1]

    def endpoint_bracket(self, q, depth=19):
        end = self.optical_end(q, depth)
        # Propagate the original second-order equation for the two independent
        # final conditions. The Dirichlet derivative normalization is arbitrary.
        def rhs(y, state):
            warp_slope, log_weight_slope = self.slopes(y)
            k2 = q*q*math.exp(-2*self.warp(y))
            drag = 2*warp_slope + log_weight_slope
            u, du, v, dv = state
            return [du, k2*u-drag*du, dv, k2*v-drag*dv]
        k_end = q*math.exp(-self.warp(end))
        sol = solve_ivp(rhs, (end,self.b), [1.,0.,0.,-k_end], method='DOP853',
                        rtol=6e-13, atol=2e-16*min(1.,q*q))
        if not sol.success:
            raise ArithmeticError(sol.message)
        u,du,v,dv = sol.y[:,-1]
        return -self.shell_weight*du/u, -self.shell_weight*dv/v

    def low_moment(self):
        # (int_y^inf G)/G = (2+tanh y)/(3(1+tanh y)^2).
        def integrand(y):
            t = math.tanh(y)
            return ((2+t)/(3*(1+t)**2))**2*self.weight(y)*math.exp(-2*self.warp(y))
        end = 32.
        value, error = quad(integrand,self.b,end,epsabs=2e-14,epsrel=2e-13)
        decay = 4-2*THROAT_SLOPE
        tail = (64/9)*math.exp(-2*THROAT_SLOPE*self.b-decay*end)/decay
        return value, error, tail

    def fraction(self, q):
        return 1-self.response(q)/(self.norm*q*q)

    def results(self):
        moment,quad_error,tail_bound = self.low_moment()
        moment_ratio = moment/self.norm
        independent_norm = quad(self.weight,self.b,32,epsabs=2e-14,epsrel=2e-13)[0]
        if abs(independent_norm-self.norm)>3e-14:
            raise ArithmeticError('Independent norm disagrees with the analytic result')
        records = {}
        fractions = []
        for q in MOMENTA:
            kernel = self.response(q)
            lo,hi = self.endpoint_bracket(q)
            discrepancy = max(abs(kernel/lo-1),abs(kernel/hi-1))
            if discrepancy>8e-10:
                raise ArithmeticError(f'Direct radial solutions disagree at {q}: {discrepancy}')
            fraction = 1-kernel/(self.norm*q*q)
            upper = q*q*moment_ratio/(1+q*q*moment_ratio)
            if not 0 < fraction <= upper+3e-11:
                raise ArithmeticError(f'Spectral fraction bound failed at {q}')
            records[str(q)] = dict(dtn=kernel, brane_propagator=1/kernel,
                                   extra_response_fraction=fraction, spectral_upper_bound=upper,
                                   finite_neumann_dtn=lo,finite_dirichlet_dtn=hi,
                                   relative_bracket_and_integration_difference=discrepancy)
            fractions.append(fraction)
        if not all(right>left for left,right in zip(fractions,fractions[1:])):
            raise ArithmeticError('Spectral fraction must increase with spacelike momentum')
        low_check = fractions[1]/(MOMENTA[1]**2*moment_ratio)
        A1,G1 = self.slopes(self.b)
        high_q = MOMENTA[-1]
        high_offset = (A1+G1)/2
        high_error = records[str(high_q)]['dtn']/self.shell_weight-high_q-high_offset
        if abs(low_check-1)>5e-5 or abs(high_error)>0.008:
            raise ArithmeticError(f'Asymptotic mismatch: {low_check}, {high_error}')
        thresholds = {}
        for target in (.01,.10,.50):
            q = brentq(lambda t:self.fraction(t)-target,.001,128.,xtol=3e-13)
            thresholds[str(target)] = q
        return dict(y_b=self.b,one_copy_norm=self.norm,shell_weight=self.shell_weight,
                    independent_norm_integral=independent_norm,first_dtn_correction=moment,
                    first_correction_quadrature_error_estimate=quad_error,
                    first_correction_omitted_tail_bound=tail_bound,correction_over_norm=moment_ratio,
                    finite_propagator_after_zero_mode_at_zero=moment/self.norm**2,
                    momentum_at_response_fraction=thresholds,
                    sufficient_distance_one_percent_static_potential=math.sqrt(4*math.exp(-2)*moment_ratio/.01),
                    low_q_leading_ratio_at_0p003=low_check,
                    high_q_dtn_over_weight_minus_q_limit=high_offset,
                    high_q_absolute_offset_error_at_512=high_error,momenta=records)


class CriticalLeaf(Leaf):
    """A Maxwell weight proportional to warp squared; normalize it to one at the wall."""
    def __init__(self, phi):
        super().__init__(phi)
        self.norm = self.tail_ratio(self.b)
        self.shell_weight = 1.

    def weight(self, y):
        return math.exp(2*self.warp(y))

    def slopes(self, y):
        warp_slope, _ = super().slopes(y)
        return warp_slope, 2*warp_slope

    def tail_ratio(self, y):
        """Exact T/w from a scaled quadrature without a small tail normalization."""
        reference = primitive(y)
        return quad(lambda d:math.exp(-2*(primitive(y+d)-reference)/3),
                    0.,80.,epsabs=2e-13,epsrel=2e-13)[0]

    def subtracted(self, q, depth=27., rtol=8e-12):
        # For t=q*exp(-A) integrate U=T/w and B=(T-M)/q².
        # Starting deep in the throat at a mixed boundary gives M/w=1/t.
        end = self.optical_end(q, depth)
        t_end = q*math.exp(-self.warp(end))
        u_end = self.tail_ratio(end)
        b_end = (u_end-1/t_end)/t_end**2
        def rhs(y,state):
            u,b = state
            t2 = (q*math.exp(-self.warp(y)))**2
            warp_slope,_ = self.slopes(y)
            return [-1-2*warp_slope*u,-(u-t2*b)**2]
        sol = solve_ivp(rhs,(end,self.b),[u_end,b_end],method='DOP853',
                        rtol=rtol,atol=rtol*.001)
        if not sol.success:
            raise ArithmeticError(sol.message)
        u,b = sol.y[:,-1]
        if abs(u-self.norm)>2e-10:
            raise ArithmeticError('Tail integral differs from the wall norm')
        return b

    def results(self):
        log_response = 1/(4*THROAT_SLOPE**3*self.norm)
        small = {}
        previous = None
        for q in (1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8):
            ratio = self.subtracted(q)/self.norm
            refined = self.subtracted(q,depth=33,rtol=1e-12)/self.norm
            if abs(ratio-refined)>2e-9:
                raise ArithmeticError('Subtracted small-Q calculation differs')
            row = dict(extra_response_fraction=q*q*ratio,
                       extra_response_fraction_over_q_squared=ratio,
                       subleading_remainder=ratio-log_response*math.log(1/q),
                       refinement_absolute_in_scaled_fraction=abs(ratio-refined))
            if previous is not None:
                row['slope_per_log_inverse_q'] = (ratio-previous)/math.log(10)
            previous = ratio
            small[str(q)] = row
        if abs(small[str(1e-8)]['slope_per_log_inverse_q']/log_response-1)>3e-7:
            raise ArithmeticError('Critical logarithm coefficient differs')
        measurements = {}
        for q in (.003,.01,.03,.1,.3,1.,3.,10.,32.,128.):
            d = self.response(q)
            f = 1-d/(q*q*self.norm)
            lo,hi = self.endpoint_bracket(q)
            err = max(abs(d/lo-1),abs(d/hi-1))
            if f<=0 or f>=1 or err>2e-9:
                raise ArithmeticError('Critical response or endpoint check failed')
            if q<=.3:
                other = q*q*self.subtracted(q)/self.norm
                if abs(f-other)>9e-12:
                    raise ArithmeticError('Direct and subtracted fractions disagree')
            measurements[str(q)]=dict(dtn=d,extra_response_fraction=f,
                                      relative_bracket_and_integration_difference=err)
        thresholds = {str(target):brentq(lambda t:self.fraction(t)-target,.001,128.,
                                          xtol=4e-13) for target in (.01,.1,.5)}
        global_wall_weight = math.exp(-2*primitive(self.b)/3)
        phi = -math.tanh(self.b)
        scalar_formula = (1+phi)**(5/9)*(1-phi)**(-1/9)*math.exp(-phi**2/9)
        if abs(global_wall_weight-scalar_formula)>3e-14:
            raise ArithmeticError('Global weight and its scalar expression disagree')
        global_norm = global_wall_weight*self.norm
        central_norm = quad(lambda y:math.exp(-2*primitive(y)/3),0.,80.,
                            epsabs=2e-13,epsrel=2e-13)[0]
        return dict(y_b=self.b,wall_weight=1.,one_copy_norm=self.norm,
                    global_wall_weight=global_wall_weight,global_one_copy_norm=global_norm,
                    log_I_over_center_for_fixed_global_weight=.5*math.log(global_norm/central_norm),
                    dtn_logarithmic_coefficient=1/(4*THROAT_SLOPE**3),
                    response_fraction_logarithmic_coefficient=log_response,
                    response_fraction_polynomial_coefficient_diverges=True,
                    momentum_at_response_fraction=thresholds,small_momenta=small,momenta=measurements)



def calculate():
    return dict(model='One flat tuned half-leaf with G=(1-phi^2)^2 and no brane kinetic term',
                units='Wall proper momentum in inverse L0; lengths in L0; g5^2 omitted; A(y_b)=0',
                boundary='Finite action in the throat; transverse current at the wall',
                leaves={str(phi):Leaf(phi).results() for phi in PHI},
                selected_weight_leaves={str(phi):CriticalLeaf(phi).results() for phi in PHI})


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('EUCLIDEAN.json'))
    arguments=parser.parse_args()
    arguments.output.write_text(json.dumps(calculate(),indent=2,allow_nan=False)+'\n')
