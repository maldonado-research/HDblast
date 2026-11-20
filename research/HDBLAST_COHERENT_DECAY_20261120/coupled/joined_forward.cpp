// Neutral coherent parent excitation on a finite Gaussian-normal annulus. Each element
// has four Lobatto nodes and an upwind numerical flux at each interface.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
using std::array;
using std::vector;

constexpr int P = 4;
double late_gauge = -1.;
double fall_phi=.36, fall_width=.02;
double prepare_phi=.18, prepare_width=.015, prepare_gauge=0., prepare_height=0.;
double radius_limit=120.;
double daughter_mass=1., amplitude=.01, support_begin=.03, support_width=.15;
double width_scale=0.,charge_mass=0.,drag_rate=0.;
double threshold_phi=.53,field_width=.003,mass_low=.95,mass_high=1.05;
constexpr double width = .001, cap = 40., calibration = .12, height = .111;
constexpr double lob[P] = {-1., -.4472135954999579393, .4472135954999579393, 1.};
constexpr double weight[P] = {1. / 6, 5. / 6, 5. / 6, 1. / 6};
enum {
    A_W,
    H,
    MASS,
    C,
    MOM,
    PHI,
    PI,
    RAD,
    SU,
    EMITTED,
    RETURNED,
    SHELL_IN,
    SHELL_OUT,
    SHELL_WORK,
    RADIATION_COMOV,
    VOLUME,
    INN_MASS,
    LATE_RETURN,
    WAVE_TO_SHEET,
    SCALAR_TO_SHEET,
    LATE_SCALAR_TO_SHEET,
    CHI, CHI_VEL, CHI_DRIVE_WORK, CHI_MASS_WORK, CHI_PRESS_WORK,
    DAUGHTER_REST, DIRECTED, LAB_EXCESS, ELECTRIC_WORK, ELECTRIC_WORK_A4, GAS_PRESS_WORK,
    DECAY_ENERGY, BIRTH_KINETIC, PAIR_COUNT, REST_WORK_MIX, REST_WORK_DRAG,
    PSI_EM_WORK, PSI_PRESS_WORK, FIELD_EXPOSURE, FORMATION_PHASE, OPTICAL_DEPTH,
    BIRTH_CHI2_OVER_K, BIRTH_AMPLITUDE2_OVER_K,
    COUNT
};
enum { I, O, GAUGE };
struct Geo {
    double r, rx, u, n, nx;
};
struct Fields {
    double i, o, a;
};
struct Values {
    double sheet, e, p, b, x, pressure, flux, k, dk, dv, gamma, acc, extra_dk;
};
struct Diag {
    double r, u, n, mass, i, o, minmass = 1e100, minn = 1e100, minr = 1e100, minf = 1e100,
                                maxrx = -1e100, ax_error = 0., ax_scale = 0., az_integral = 0.,
                                a_integral_error = 0., a_wall_error = 0., fried = 0., ew = 0., pw = 0.,
                                sh = 0.;
};

double logistic(double z) { return z > 0. ? 1 / (1 + std::exp(-z)) : std::exp(z) / (1 + std::exp(z)); }
array<double, 2> step(double x) {
    if (x <= 0)
        return {0, 0};
    if (x >= 1)
        return {1, 0};
    double z = 1 / (1 - x) - 1 / x, f = logistic(z);
    return {f, logistic(-std::abs(z)) * logistic(std::abs(z)) * (1 / (x * x) + 1 / ((1 - x) * (1 - x)))};
}
double incoming(double tau, double dose, double norm, double stop) {
    double s = tau / width;
    if (s <= 0 || s >= stop)
        return 0.;
    double z = 4 * (s - 3.);
    return std::sqrt(dose / width * logistic(z)) * logistic(-z) * step(s / .35)[0] *
           step((stop - s) / .125)[0] / norm;
}
struct Oscillator {
    double minimum, minimum_slope, mass, mass_slope, delta, energy, pressure, drive_slope, mass_gradient, mass_ratio, width, power;
};
Oscillator oscillator(const vector<double>& y) {
    const double psi=y[PHI], half=.5*support_width;
    const auto up=step((psi-support_begin)/half), down=step((psi-support_begin-half)/half);
    const auto opening=step((psi-threshold_phi)/field_width+.5);
    Oscillator o{};
    o.minimum=amplitude*(up[0]-down[0]);
    o.minimum_slope=amplitude*(up[1]-down[1])/half;
    o.mass=2.*daughter_mass*(mass_low+(mass_high-mass_low)*opening[0]);
    o.mass_slope=2.*daughter_mass*(mass_high-mass_low)*opening[1]/field_width;
    o.mass_ratio=o.mass/(2.*daughter_mass);
    const double beta2=std::max(0.,1.-1./(o.mass_ratio*o.mass_ratio));
    o.width=width_scale*o.mass_ratio*beta2*std::sqrt(beta2);
    o.power=o.width*y[CHI_VEL]*y[CHI_VEL];
    o.delta=y[CHI]-o.minimum;
    const double kinetic=.5*y[CHI_VEL]*y[CHI_VEL], potential=.5*std::pow(o.mass*o.delta,2);
    o.energy=kinetic+potential;
    o.pressure=kinetic-potential;
    o.drive_slope=-o.mass*o.mass*o.delta*o.minimum_slope;
    o.mass_gradient=o.mass*o.mass_slope*o.delta*o.delta;
    return o;
}

struct Gas {
    double comoving_rest, dilution;
    double cold, excess, lab_pressure, gas_momentum, energy, pressure, velocity, gamma, theta;
    double current, random_lab, intrinsic_random, rest, u, reconstruction_error;
};

Gas charged_gas(const vector<double> &y) {
    Gas g{};
    g.comoving_rest=y[DAUGHTER_REST];
    g.dilution=1./std::pow(y[A_W],3);
    g.gamma=1.;
    if(g.comoving_rest<=0.) return g;
    g.rest=g.comoving_rest*g.dilution;
    const double moment=y[DIRECTED], extra=y[LAB_EXCESS];
    const double cold=std::hypot(g.comoving_rest,moment), E=cold+extra;
    g.cold=cold; g.excess=extra; g.gas_momentum=moment;
    // Invert (h-Theta)(h-4Theta)=1 for nonnegative rest-frame pressure P.
    // Normalize by E. Compute E^2-Z^2 from U to retain the cold limit.
    const double z=g.comoving_rest/E, zz=z*z, surplus=extra/E*((2.*cold+extra)/E);
    const double A=zz+surplus, rhs=A*surplus, quad=3.*A+2.+zz, lin=A+2.*zz;
    double p=0.;
    if(rhs>0.) {
        p=std::min({1./3.,rhs/lin,std::sqrt(rhs/quad)});
        for(int j=0;j<10;++j) {
            double dp=(((3.*p+quad)*p+lin)*p-rhs)/((9.*p+2.*quad)*p+lin);
            p-=dp;
            if(std::abs(dp)<4e-15*p) break;
        }
    }
    g.lab_pressure=E*p; // Comoving rest-frame pressure, before the isotropic stream sum.
    g.velocity=(moment/E)/(1.+p);
    g.gamma=(1.+p)/std::sqrt(A+p*(2.+p));
    g.u=g.gamma*g.velocity;
    g.theta=p*g.gamma/z;
    double root=std::hypot(1.,1.5*g.theta);
    double thermal=1.5*g.theta+(2.25*g.theta*g.theta)/(1.+root);
    g.energy=g.dilution*E;
    g.pressure=g.dilution*(g.lab_pressure+moment*g.velocity/3.);
    g.current=g.rest*charge_mass*g.velocity;
    g.intrinsic_random=g.rest*thermal;
    g.random_lab=g.rest*((thermal+g.theta)*g.gamma-g.theta/g.gamma);
    g.reconstruction_error=g.energy-(g.rest*g.gamma+g.random_lab);
    return g;
}

Values wall(const vector<double> &y, double i) {
    const double a = y[A_W], psi = y[PHI], pi = y[PI], strength = std::sqrt(8 / calibration) / width;
    const double l = std::log(cap), s = std::abs(strength * psi), q = (s - l + .2) / .4;
    double f = s, der = 1.;
    if (q >= 1) {
        f = l;
        der = 0;
    } else if (q > 0) {
        f = l - .2 + .4 * (q - 2.5 * std::pow(q, 4) + 3 * std::pow(q, 5) - std::pow(q, 6));
        der = std::pow(1 - q, 3) * (1 + 3 * q + 6 * q * q);
    }
    double k = width * std::exp(std::copysign(f, psi)), dk = k * strength * der;
    auto rise = step((psi - .004) / .008), fall = step((psi - fall_phi) / fall_width);
    // The preparation interval ends before the final release starts.
    const auto prep=step((psi-prepare_phi)/prepare_width);
    const double swing=late_gauge-prepare_gauge;
    const double k_pre=k+prepare_gauge*prep[0], dk_pre=dk+prepare_gauge*prep[1]/prepare_width;
    double release_dk=swing*fall[1]/fall_width;
    // Magnetic energy is linear in K at frozen flux and electric energy is
    // linear in 1/K at frozen displacement. Use the latter for a downward step.
    if(prepare_gauge>late_gauge && fall[0]>0.) {
        const double k_after=k+late_gauge, s=fall[0];
        const double reciprocal_swing=1./k_after-1./k_pre;
        k=1./((1.-s)/k_pre+s/k_after);
        release_dk=-k*k*fall[1]*reciprocal_swing/fall_width;
        dk=k*k*((1.-s)*dk_pre/(k_pre*k_pre)+s*dk/(k_after*k_after))+release_dk;
    } else {
        k=k_pre+swing*fall[0];
        dk=dk_pre+release_dk;
    }
    double v=height*rise[0]-prepare_height*prep[0]-(height-prepare_height)*fall[0];
    double dv=height*rise[1]/.008-prepare_height*prep[1]/prepare_width-(height-prepare_height)*fall[1]/fall_width;
    const double gamma=0.;
    double b = 2 * y[C] / a, x = y[MOM] / k, o = x - i, sheet = .5 * k * (x * x + b * b);
    double e = .5 * pi * pi + v + y[RAD] + sheet, p = .5 * pi * pi - v + (y[RAD] + sheet) / 3;
    const auto chi=oscillator(y);
    const auto gas=charged_gas(y);
    e+=chi.energy+gas.energy;
    p+=chi.pressure+gas.pressure;
    double press = i * i + o * o - .5 * b * b, flux = i * i - o * o;
    return {sheet,
            e,
            p,
            b,
            x,
            press,
            flux,
            k,
            dk,
            dv,
            gamma,
            (1 + e) * (1 - 2 * e - 3 * p) - 1 - y[MASS] / std::pow(a, 4) - 2 * press,
            release_dk};
}
double localmass(Geo g) { return g.r * g.r * (1 + g.r * g.r + g.u * g.u - g.rx * g.rx); }
Geo radial(Geo g, Fields w) {
    double e = w.i + w.o, z = w.i - w.o, b = 2 * w.a / (g.r * g.r),
           mr = g.r * (1 + g.r * g.r + g.u * g.u - g.rx * g.rx);
    // mr is m/r. Keep T_nn independent of flux for mixed characteristic data.
    return {g.rx, g.r + mr / (g.r * g.r) - g.r * (e * e + z * z + b * b), -2 * g.r * e * z, g.nx,
            g.n * (1 - 3 * mr / (g.r * g.r * g.r) + 5 * b * b - z * z + e * e)};
}
Geo shifted(Geo g, Geo d, double dx) {
    return {g.r + dx * d.r, g.rx + dx * d.rx, g.u + dx * d.u, g.n + dx * d.n, g.nx + dx * d.nx};
}
Geo integrate(Geo g, Fields a, Fields m, Fields b, double dx) {
    Geo k1 = radial(g, a), k2 = radial(shifted(g, k1, dx / 2), m),
        k3 = radial(shifted(g, k2, dx / 2), m), k4 = radial(shifted(g, k3, dx), b);
    return shifted(g,
                   {k1.r + 2 * k2.r + 2 * k3.r + k4.r, k1.rx + 2 * k2.rx + 2 * k3.rx + k4.rx,
                    k1.u + 2 * k2.u + 2 * k3.u + k4.u, k1.n + 2 * k2.n + 2 * k3.n + k4.n,
                    k1.nx + 2 * k2.nx + 2 * k3.nx + k4.nx},
                   dx / 6);
}

struct Solver {
    int cells;
    double dx, dose, norm, stop, D[P][P]{}, mid[P - 1][P]{};
    Solver(int cells_, double d, double dose_, double norm_, double stop_)
        : cells(cells_), dx(d / cells), dose(dose_), norm(norm_), stop(stop_) {
        for (int k = 0; k < P; k++)
            for (int j = 0; j < P; j++) {
                if (k != j) {
                    double c = 1 / (lob[j] - lob[k]);
                    for (int l = 0; l < P; l++)
                        if (l != j && l != k)
                            c *= (lob[k] - lob[l]) / (lob[j] - lob[l]);
                    D[k][j] = c;
                } else
                    for (int l = 0; l < P; l++)
                        if (l != j)
                            D[k][j] += 1 / (lob[j] - lob[l]);
            }
        for (int k = 0; k < P - 1; k++)
            for (int j = 0; j < P; j++) {
                double x = .5 * (lob[k] + lob[k + 1]), c = 1;
                for (int l = 0; l < P; l++)
                    if (l != j)
                        c *= (x - lob[l]) / (lob[j] - lob[l]);
                mid[k][j] = c;
            }
    }
    int index(int n, int channel) { return COUNT + 3 * n + channel; }
    void rhs(const vector<double> &y, vector<double> &dy, Diag *diag = nullptr) {
        double a = y[A_W], h = y[H], pi = y[PI], iwall = y[index(0, I)];
        Values w = wall(y, iwall);
        Geo g = {a, -a * (1 + w.e), a * h, 1, -(1 - 2 * w.e - 3 * w.p)};
        double emitted = incoming(y[SU], dose, norm, stop), outgoing = 0.;
        for (int cell = 0; cell < cells; cell++) {
            Geo geom[P];
            Fields field[P];
            const int base = P * cell;
            for (int k = 0; k < P; k++) {
                int j = index(base + k, 0);
                field[k] = {y[j + I], y[j + O], y[j + GAUGE]};
            }
            geom[0] = g;
            for (int k = 0; k < P - 1; k++) {
                Fields m{0., 0., 0.};
                for (int j = 0; j < P; j++) {
                    m.i += mid[k][j] * field[j].i;
                    m.o += mid[k][j] * field[j].o;
                    m.a += mid[k][j] * field[j].a;
                }
                g = integrate(g, field[k], m, field[k + 1], dx * (lob[k + 1] - lob[k]) / 2);
                geom[k + 1] = g;
            }
            for (int k = 0; k < P; k++) {
                double di = 0., doo = 0., da = 0.;
                for (int j = 0; j < P; j++) {
                    di += D[k][j] * field[j].i;
                    doo += D[k][j] * field[j].o;
                    if (diag)
                        da += D[k][j] * field[j].a;
                }
                // Upwind SAT in the strong nodal equations for I_t=N I_x and O_t=-N O_x.
                int node = base + k;
                if (k == P - 1) {
                    double right = cell + 1 == cells ? emitted : y[index(node + 1, I)];
                    di += (right - field[k].i) / weight[k];
                }
                if (k == 0) {
                    double left = cell == 0 ? w.x - iwall : y[index(node - 1, O)];
                    doo -= (left - field[k].o) / weight[k];
                }
                const Geo &z = geom[k];
                double e = field[k].i + field[k].o, bz = field[k].i - field[k].o;
                double se = (z.nx + 2 * z.n * z.rx / z.r) * bz - 2 * z.n * z.u * e / z.r -
                            4 * z.n * field[k].a / std::pow(z.r, 3);
                double sz = (z.nx + z.n * z.rx / z.r) * e - z.n * z.u * bz / z.r;
                int j = index(node, 0);
                dy[j + I] = 2 * z.n * di / dx + .5 * (se + sz);
                dy[j + O] = -2 * z.n * doo / dx + .5 * (se - sz);
                dy[j + GAUGE] = z.n * z.r * e;
                if (diag) {
                    diag->minmass = std::min(diag->minmass, localmass(z));
                    diag->minr = std::min(diag->minr, z.r);
                    diag->minf = std::min(diag->minf, z.rx * z.rx - z.u * z.u);
                    diag->minn = std::min(diag->minn, z.n);
                    diag->maxrx = std::max(diag->maxrx, z.rx);
                    diag->ax_error = std::max(diag->ax_error, std::abs(2 * da / dx - z.r * bz));
                    diag->ax_scale = std::max(diag->ax_scale, std::abs(z.r * bz));
                    diag->az_integral += dx * .5 * weight[k] * z.r * bz;
                }
            }
            outgoing = field[P - 1].o;
        }
        double ow = w.x - iwall, mu = y[MASS], r = g.r,
               bmag = 2 * y[index(P * cells - 1, GAUGE)] / (r * r);
        double fin = emitted * emitted - outgoing * outgoing,
               pr = emitted * emitted + outgoing * outgoing - .5 * bmag * bmag;
        dy[A_W] = a * h;
        dy[H] = w.acc - h * h;
        dy[MASS] = -4 * std::pow(a, 4) * (w.pressure * h + w.flux * (1 + w.e));
        dy[C] = w.x - h * y[C];
        dy[MOM] = -2 * h * y[MOM] - 4 * w.k * y[C] / (a * a) - 2 * (ow - iwall);
        const double a3=a*a*a;
        dy[PHI] = pi;
        dy[PI] = .5 * w.dk * (w.x * w.x - w.b * w.b) - w.dv - (3 * h + w.gamma) * pi;
        const auto chi=oscillator(y);
        dy[CHI]=y[CHI_VEL];
        dy[CHI_VEL]=-(3.*h+chi.width)*y[CHI_VEL]-chi.mass*chi.mass*chi.delta;
        dy[PI]-=chi.drive_slope+chi.mass_gradient;
        dy[CHI_DRIVE_WORK]=a3*pi*chi.drive_slope;
        dy[CHI_MASS_WORK]=a3*pi*chi.mass_gradient;
        dy[CHI_PRESS_WORK]=3.*h*a3*chi.pressure;
        const auto gas=charged_gas(y);
        const double moment=y[DIRECTED];
        const double vc=gas.cold>0. ? moment/gas.cold : 0.;
        const double inv_gamma_c=gas.cold>0. ? gas.comoving_rest/gas.cold : 1.;
        const double dv=gas.cold>0. ? -vc*(gas.excess+gas.lab_pressure)/(gas.cold+gas.excess+gas.lab_pressure) : 0.;
        dy[MOM]-=gas.current;
        dy[DECAY_ENERGY]=a3*chi.power;
        dy[DAUGHTER_REST]=dy[DECAY_ENERGY]/chi.mass_ratio;
        dy[BIRTH_KINETIC]=dy[DECAY_ENERGY]*(1.-1./chi.mass_ratio);
        dy[PAIR_COUNT]=dy[DECAY_ENERGY]/chi.mass;
        dy[FIELD_EXPOSURE]=gas.comoving_rest*std::abs(w.x);
        dy[FORMATION_PHASE]=std::max(0.,chi.mass-2.*daughter_mass);
        dy[OPTICAL_DEPTH]=chi.width;
        dy[BIRTH_CHI2_OVER_K]=dy[DAUGHTER_REST]*std::pow(y[CHI],2)/w.k;
        dy[BIRTH_AMPLITUDE2_OVER_K]=dy[BIRTH_CHI2_OVER_K]+dy[DAUGHTER_REST]*std::pow(y[CHI_VEL]/chi.mass,2)/w.k;
        dy[DIRECTED]=charge_mass*gas.comoving_rest*w.x-(h+drag_rate)*moment;
        // d(E-sqrt(D²+S²)) remains well conditioned for newly born cold gas.
        dy[LAB_EXCESS]=charge_mass*gas.comoving_rest*w.x*dv+h*(-moment*dv-3.*gas.lab_pressure)
                     +drag_rate*moment*vc+(chi.mass_ratio-inv_gamma_c)*dy[DAUGHTER_REST];
        dy[ELECTRIC_WORK]=a3*gas.current*w.x;
        dy[ELECTRIC_WORK_A4]=a*dy[ELECTRIC_WORK];
        dy[GAS_PRESS_WORK]=3.*a3*h*gas.pressure;
        dy[REST_WORK_DRAG]=drag_rate*gas.u*moment;
        dy[REST_WORK_MIX]=(gas.gamma-1.)*dy[DECAY_ENERGY];
        // A separate driver-scalar ledger records the opposite energy transfer.
        const double driver_energy=w.e-y[RAD]-w.sheet-chi.energy-gas.energy;
        dy[PSI_EM_WORK]=.5*a3*pi*w.dk*(w.x*w.x-w.b*w.b);
        dy[PSI_PRESS_WORK]=3.*h*a3*(pi*pi-driver_energy);
        dy[RAD] = w.gamma * pi * pi - 4 * h * y[RAD];
        dy[SU] = g.n;
        dy[EMITTED] = 2 * g.n * r * r * r * emitted * emitted;
        dy[RETURNED] = 2 * g.n * r * r * r * outgoing * outgoing;
        dy[SHELL_IN] = 2 * a * a * a * iwall * iwall;
        dy[SHELL_OUT] = 2 * a * a * a * ow * ow;
        dy[SHELL_WORK] = 3 * a * a * a * h * w.p;
        dy[RADIATION_COMOV] = std::pow(a, 4) * w.gamma * pi * pi;
        dy[VOLUME] = g.n * g.u;
        dy[INN_MASS] = -4 * r * r * r * g.n * (pr * g.u - fin * g.rx);
        dy[LATE_RETURN] = y[SU] > .1 ? dy[RETURNED] : 0.;
        double drive = .5 * std::pow(a, 4) * pi * (w.b * w.b - w.x * w.x);
        dy[WAVE_TO_SHEET] = 2 * std::pow(a, 4) * w.flux;
        dy[SCALAR_TO_SHEET] = drive * w.dk;
        dy[LATE_SCALAR_TO_SHEET] = drive * w.extra_dk;
        if (diag) {
            diag->r = r;
            diag->u = g.u;
            diag->n = g.n;
            diag->mass = localmass(g);
            diag->i = emitted;
            diag->o = outgoing;
            diag->ew = w.e;
            diag->pw = w.p;
            diag->sh = w.sheet;
            diag->a_wall_error = std::abs(y[index(0, GAUGE)] - a * y[C]);
            diag->a_integral_error =
                std::abs(y[index(P * cells - 1, GAUGE)] - a * y[C] - diag->az_integral);
            diag->fried = h * h + 1 / (a * a) - w.e * (2 + w.e) - mu / std::pow(a, 4);
        }
    }
};

int main(int argc, char **argv) {
    if (argc != 30) {
        std::cerr << "Run with joined.py or pass the documented parameters.\n";
        return 1;
    }
    int cells = std::stoi(argv[1]);
    double depth = std::stod(argv[2]), dose = std::stod(argv[3]), end = std::stod(argv[4]);
    double norm = std::stod(argv[5]), stop = std::stod(argv[6]), cfl = std::stod(argv[7]);
    int late_cells = std::stoi(argv[8]);
    double transfer_time = std::stod(argv[9]);
    late_gauge = std::stod(argv[10]);
    threshold_phi=std::stod(argv[12]);
    field_width=std::stod(argv[13]);
    mass_low=std::stod(argv[14]);
    mass_high=std::stod(argv[15]);
    fall_phi=std::stod(argv[16]); fall_width=std::stod(argv[17]);
    prepare_phi=std::stod(argv[18]); prepare_width=std::stod(argv[19]);
    prepare_gauge=std::stod(argv[20]); prepare_height=std::stod(argv[21]); radius_limit=std::stod(argv[22]);
    daughter_mass=std::stod(argv[23]); amplitude=std::stod(argv[24]); support_begin=std::stod(argv[25]); support_width=std::stod(argv[26]);
    width_scale=std::stod(argv[27]);charge_mass=std::stod(argv[28]);drag_rate=std::stod(argv[29]);
    if(width_scale<0. || charge_mass<0. || drag_rate<0.) throw std::runtime_error("Nonnegative width, charge and drag required");
    if(daughter_mass<=0. || support_width<=0. || support_begin<0. || support_begin+support_width>=threshold_phi-field_width/2.)
        throw std::runtime_error("Drive support must end before the mass ramp");
    if (fall_phi<=.012 || fall_width<=0. || prepare_width<=0. || prepare_gauge<0. || prepare_height<0. || prepare_height>=height ||
        ((prepare_gauge>0. || prepare_height>0.) && (prepare_phi<.012 || prepare_phi+prepare_width>fall_phi)) || radius_limit<=2.)
        throw std::runtime_error("Invalid staged scalar or gauge interaction");
    if (!(threshold_phi>=field_width/2. && field_width>0. &&
          mass_low>0. && mass_low<1. && mass_high>1. && std::abs(.5*(mass_low+mass_high)-1.)<1.e-14))
        throw std::runtime_error("Invalid decay parameters");
    if (late_gauge<0.)
        throw std::runtime_error("Invalid parent, daughter or gauge parameter");
    if (cells < 1 || late_cells < 1 || late_cells > cells || depth <= 0 || cfl <= 0)
        throw std::runtime_error("Invalid mesh");
    Solver s(cells, depth, dose, norm, stop);
    int n = COUNT + 3 * cells * P;
    vector<double> y(n, 0.), dy(n), stage(n), acc(n);
    y[A_W] = 2;
    // Choose mu using the exact initial Friedmann constraint at a=2,H=0.
    y[MASS] = 4.;
    y[INN_MASS] = y[MASS];
    Diag initial;
    s.rhs(y, dy, &initial);
    y[VOLUME] = initial.r;
    y[INN_MASS] = initial.mass;
    double t = 0., dt0 = cfl * s.dx;
    int sampling = std::stoi(argv[11]);
    bool want_samples = sampling!=0;
    vector<double> times = {0., .001, .002, .003, .004, .005, .006, .008, .01, .02, .05};
    for (int j = 1; j <= static_cast<int>(10 * end); j++)
        times.push_back(j / 10.);
    if(want_samples) {
        for(int i=0;i<=std::min(2.,end)*25;i++) times.push_back(i*.04);
        // The primary trace has dense coverage of decay and first turnaround.
        int samples_per_time=(sampling==2?1000:50);
        for(int i=5*samples_per_time;i<=9*samples_per_time;i++) times.push_back(double(i)/samples_per_time);
        std::sort(times.begin(),times.end());
        times.erase(std::unique(times.begin(),times.end(),[](double a,double b){return std::abs(a-b)<1e-9;}),times.end());
    }
    size_t next_sample = 0;
    auto sample = [&](double time, const vector<double> &state) {
        Values mat = wall(state, state[s.index(0, I)]);
        const auto chi=oscillator(state);
        const auto gas=charged_gas(state);
        double potential =
            height*step((state[PHI]-.004)/.008)[0]-prepare_height*step((state[PHI]-prepare_phi)/prepare_width)[0]-(height-prepare_height)*step((state[PHI]-fall_phi)/fall_width)[0];
        std::cout << "@ " << std::setprecision(16) << time << ' ' << state[A_W] << ' ' << state[H] << ' '
                  << state[PHI] << ' ' << state[PI] << ' ' << state[RAD] << ' ' << mat.sheet << ' '
                  << mat.e << ' ' << potential << ' ' << state[SU] << ' ' << state[EMITTED] << ' '
                  << state[RETURNED] << ' ' << state[VOLUME] << ' ' << mat.k << ' ' << mat.b << ' '
                  << mat.x << ' ' << state[MASS] << ' ' << mat.p << ' ' << state[WAVE_TO_SHEET] << ' ' << state[SCALAR_TO_SHEET]
                  << ' ' << state[CHI] << ' ' << state[CHI_VEL] << ' ' << chi.minimum << ' ' << chi.mass
                  << ' ' << chi.energy << ' ' << chi.pressure << ' ' << state[CHI_DRIVE_WORK] << ' ' << state[CHI_MASS_WORK]
                  << ' ' << state[CHI_PRESS_WORK] << ' ' << chi.drive_slope << ' ' << chi.mass_gradient << ' ' << state[PSI_EM_WORK] << ' ' << state[PSI_PRESS_WORK] << ' ' << chi.width << ' ' << chi.power
                  << ' ' << state[DAUGHTER_REST] << ' ' << state[DIRECTED] << ' ' << state[LAB_EXCESS]
                  << ' ' << gas.energy << ' ' << gas.pressure << ' ' << gas.intrinsic_random << ' ' << gas.theta
                  << ' ' << gas.velocity << ' ' << gas.current << ' ' << state[DECAY_ENERGY] << ' ' << state[PAIR_COUNT]
                  << ' ' << state[BIRTH_KINETIC] << ' ' << state[ELECTRIC_WORK] << ' ' << state[ELECTRIC_WORK_A4]
                  << ' ' << state[GAS_PRESS_WORK] << ' ' << state[REST_WORK_DRAG] << ' ' << state[REST_WORK_MIX] << ' ' << state[FIELD_EXPOSURE] << ' ' << state[FORMATION_PHASE]
                  << ' ' << state[OPTICAL_DEPTH] << ' ' << state[BIRTH_CHI2_OVER_K] << ' ' << state[BIRTH_AMPLITUDE2_OVER_K] << '\n';
    };
    std::map<std::string, double> out{{"min_lapse", 1.},
                                      {"min_areal_radius", 2.},
                                      {"min_radial_gradient_norm", 100.},
                                      {"min_bulk_mass", 4.},
                                      {"max_inward_radius_derivative", -1e9}};
    auto peak = [&](const std::string &key, double x) { out[key] = std::max(out[key], std::abs(x)); };
    auto event = [&](const std::string &name, double t0, double a0, double h0, double t1, double a1,
                     double h1) {
        if (out.count(name + "_tau"))
            return;
        double frac = -h0 / (h1 - h0);
        out[name + "_tau"] = t0 + (t1 - t0) * frac;
        out[name + "_radius"] = a0 + (a1 - a0) * frac;
    };
    int steps = 0;
    double hot_start=-1.;
    bool transferred = false;
    double last_t = 0., last_r = initial.r, last_Nu = initial.n * initial.u;
    while (t + 1e-13 < end) {
        if (!transferred && t >= transfer_time - 1e-14 && cells != late_cells) {
            Diag before;
            s.rhs(y, dy, &before);
            int oldcells = cells;
            vector<double> old = y;
            cells = late_cells;
            s = Solver(cells, depth, dose, norm, stop);
            n = COUNT + 3 * cells * P;
            y.resize(n);
            dy.resize(n);
            stage.resize(n);
            acc.resize(n);
            for (int cell = 0; cell < cells; cell++)
                for (int k = 0; k < P; k++) {
                    double oldcoord = (cell + (1 + lob[k]) * .5) * oldcells / cells;
                    int from = std::min(static_cast<int>(oldcoord), oldcells - 1);
                    double xi = 2 * (oldcoord - from) - 1;
                    for (int channel = 0; channel < 3; channel++) {
                        double interpolated = 0.;
                        for (int j = 0; j < P; j++) {
                            double l = 1.;
                            for (int m = 0; m < P; m++)
                                if (j != m)
                                    l *= (xi - lob[m]) / (lob[j] - lob[m]);
                            interpolated += l * old[COUNT + 3 * (P * from + j) + channel];
                        }
                        y[s.index(cell * P + k, channel)] = interpolated;
                    }
                }
            Diag after;
            s.rhs(y, dy, &after);
            out["remesh_tube_radius_jump"] = std::abs(before.r - after.r);
            out["remesh_tube_mass_jump"] = std::abs(before.mass - after.mass);
            dt0 = cfl * s.dx;
            transferred = true;
        }
        Diag z;
        bool check = steps % 32 == 0 || t < .025;
        s.rhs(y, dy, check ? &z : nullptr);
        if (check) {
            out["min_lapse"] = std::min(out["min_lapse"], z.minn);
            out["min_areal_radius"] = std::min(out["min_areal_radius"], z.minr);
            out["min_radial_gradient_norm"] = std::min(out["min_radial_gradient_norm"], z.minf);
            out["min_bulk_mass"] = std::min(out["min_bulk_mass"], z.minmass);
            out["max_inward_radius_derivative"] = std::max(out["max_inward_radius_derivative"], z.maxrx);
            peak("friedmann_error", z.fried);
            peak("tube_radius_error", z.r - y[VOLUME]);
            peak("tube_mass_error", z.mass - y[INN_MASS]);
            peak("gauge_gradient_error", z.ax_error);
            peak("gauge_gradient_scale", z.ax_scale);
            peak("gauge_radial_integral_error", z.a_integral_error);
            if (t > last_t)
                peak("tube_radius_velocity_residual",
                     (z.r - last_r) / (t - last_t) - .5 * (last_Nu + z.n * z.u));
            last_t = t;
            last_r = z.r;
            last_Nu = z.n * z.u;
            peak("boundary_A_error", z.a_wall_error);
            const auto chi=oscillator(y);
            const auto gas=charged_gas(y);
            const double driver_energy=z.ew-y[RAD]-z.sh-chi.energy-gas.energy;
            peak("driver_balance_error",std::pow(y[A_W],3)*driver_energy+y[PSI_PRESS_WORK]-y[PSI_EM_WORK]+y[CHI_DRIVE_WORK]+y[CHI_MASS_WORK]);
            peak("oscillator_balance_error",std::pow(y[A_W],3)*chi.energy+y[CHI_PRESS_WORK]+y[DECAY_ENERGY]-y[CHI_DRIVE_WORK]-y[CHI_MASS_WORK]);
            peak("surface_balance_error",
                 std::pow(y[A_W], 3) * z.ew + y[SHELL_WORK] - y[SHELL_IN] + y[SHELL_OUT]);
            peak("sheet_balance_error",
                 std::pow(y[A_W], 4) * z.sh - y[WAVE_TO_SHEET] - y[SCALAR_TO_SHEET]+y[ELECTRIC_WORK_A4]);
            peak("gas_balance_error",std::pow(y[A_W],3)*gas.energy+y[GAS_PRESS_WORK]-y[ELECTRIC_WORK]-y[DECAY_ENERGY]);
            peak("birth_rest_error",y[DECAY_ENERGY]-y[BIRTH_KINETIC]-y[DAUGHTER_REST]);
            peak("pair_number_error",2.*daughter_mass*y[PAIR_COUNT]-y[DAUGHTER_REST]);
            peak("gas_primitive_error",gas.reconstruction_error);
            if(y[DAUGHTER_REST]<0. || y[LAB_EXCESS]<-1e-14 || !std::isfinite(gas.theta) || std::abs(gas.velocity)>1.)
                throw std::runtime_error("Nonphysical daughter state");
            peak("fluid_balance_error", std::pow(y[A_W], 4) * y[RAD] - y[RADIATION_COMOV]);
            if (!std::isfinite(z.mass) || z.minn <= 0)
                throw std::runtime_error("Gaussian chart failed");
        }
        const std::array<double,3> old_work{y[CHI_DRIVE_WORK],y[CHI_MASS_WORK],y[CHI_PRESS_WORK]};
        const auto previous_gas=charged_gas(y);
        const auto previous_wall=wall(y,y[s.index(0,I)]);
        const std::array<double,7> old_sources{y[DAUGHTER_REST],y[DECAY_ENERGY],y[BIRTH_KINETIC],y[FIELD_EXPOSURE],
                                                y[ELECTRIC_WORK],y[GAS_PRESS_WORK],y[FORMATION_PHASE]};
        const double old_margin=previous_gas.intrinsic_random-.9*previous_wall.e;
        const double old_chi=y[CHI],old_chi_vel=y[CHI_VEL],oldpi=y[PI];
        double olda = y[A_W], oldh = y[H], oldphi = y[PHI];
        double dt = std::min(dt0, end - t);
        if (!transferred && cells != late_cells && t < transfer_time - 1e-14)
            dt = std::min(dt, transfer_time - t);
        if(charge_mass>0. && previous_gas.comoving_rest>0.) {
            const double enthalpy=2.5*previous_gas.theta+std::hypot(1.,1.5*previous_gas.theta);
            const double frequency=charge_mass*std::sqrt(previous_gas.rest/(previous_wall.k*enthalpy*previous_gas.gamma));
            if (frequency>0.) dt=std::min(dt,.5/frequency);
            const double size=std::max(1e-26,previous_gas.comoving_rest);
            if(dy[DIRECTED]!=0.)dt=std::min(dt,.05*(size+std::abs(y[DIRECTED]))/std::abs(dy[DIRECTED]));
            if(dy[LAB_EXCESS]!=0.)dt=std::min(dt,.05*(size+std::abs(y[LAB_EXCESS]))/std::abs(dy[LAB_EXCESS]));
        }
        bool due = want_samples && next_sample < times.size() && times[next_sample] <= t + dt;
        vector<double> old_sample;
        if (due)
            old_sample = y;
        for (int k = 0; k < n; k++) {
            acc[k] = dy[k];
            stage[k] = y[k] + .5 * dt * dy[k];
        }
        s.rhs(stage, dy);
        for (int k = 0; k < n; k++) {
            acc[k] += 2 * dy[k];
            stage[k] = y[k] + .5 * dt * dy[k];
        }
        s.rhs(stage, dy);
        for (int k = 0; k < n; k++) {
            acc[k] += 2 * dy[k];
            stage[k] = y[k] + dt * dy[k];
        }
        s.rhs(stage, dy);
        for (int k = 0; k < n; k++)
            y[k] += dt / 6 * (acc[k] + dy[k]);
        steps++;
        if (due) {
            double f = (times[next_sample] - t) / dt;
            for (int k = 0; k < n; k++)
                old_sample[k] += (y[k] - old_sample[k]) * f;
            sample(times[next_sample++], old_sample);
        }
        Values wnew = wall(y, y[s.index(0, I)]);
        const double a = y[A_W], h = y[H];
        const auto capture=[&](const std::string& name, double field) {
            if (!(oldphi<field && y[PHI]>=field) || out.count(name+"_action")) return;
            const double z=(field-oldphi)/(y[PHI]-oldphi);
            vector<double> ev=y;
            ev[PHI]=field;
            ev[CHI]=y[CHI]+(1.-z)*(old_chi-y[CHI]);
            ev[CHI_VEL]=y[CHI_VEL]+(1.-z)*(old_chi_vel-y[CHI_VEL]);
            auto c=oscillator(ev);
            const double radius=olda+z*(a-olda), hubble=oldh+z*(h-oldh);
            const double action=std::pow(radius,3)*c.energy/c.mass;
            out[name+"_time"]=t+z*dt; out[name+"_radius"]=radius;out[name+"_H"]=hubble;
            out[name+"_mass"]=c.mass;out[name+"_density"]=c.energy;out[name+"_pressure"]=c.pressure;
            out[name+"_action"]=action;out[name+"_target_ratio"]=action/(.0008/(2.*daughter_mass*mass_low));
            out[name+"_mass_WKB"]=std::abs(c.mass_slope*(oldpi+z*(y[PI]-oldpi)))/(c.mass*c.mass);
            out[name+"_H_over_M"]=std::abs(hubble)/c.mass;
            out[name+"_drive_work"]=old_work[0]+z*(y[CHI_DRIVE_WORK]-old_work[0]);
            out[name+"_mass_work"]=old_work[1]+z*(y[CHI_MASS_WORK]-old_work[1]);
            out[name+"_pressure_work"]=old_work[2]+z*(y[CHI_PRESS_WORK]-old_work[2]);
        };
        capture("drive_begin",support_begin);
        capture("drive_end",support_begin+support_width);
        capture("before_mass",threshold_phi-.5*field_width);
        capture("opening",threshold_phi);
        capture("mass_complete",threshold_phi+.5*field_width);
        if ((prepare_gauge>0. || prepare_height>0.) && oldphi < prepare_phi && y[PHI]>=prepare_phi && !out.count("prepare_radius")) {
            out["prepare_radius"]=a;out["prepare_time"]=t+dt;out["prepare_B"]=wnew.b;out["prepare_X"]=wnew.x;out["prepare_pi"]=y[PI];
        }
        if (oldphi < fall_phi && y[PHI] >= fall_phi && !out.count("first_exit_Bmag")) {
            out["first_exit_Bmag"] = wnew.b;
            out["first_exit_X"] = wnew.x;
            out["first_exit_K"]=wnew.k;out["first_exit_aH"]=a*h;out["first_exit_pi"]=y[PI];
            out["first_exit_magnetic_budget_K"] = (height-prepare_height) / (.5 * wnew.b * wnew.b);
            out["first_exit_electric_budget_K"] = std::pow(y[MOM],2)/(2.*(height-prepare_height));
        }
        if (oldphi < fall_phi+fall_width && y[PHI] >= fall_phi+fall_width && !out.count("exit_complete_tau")) {
            out["exit_complete_aH"]=a*h;out["exit_complete_field_density"]=wnew.sheet;
            out["exit_complete_electric_fraction"]=wnew.x*wnew.x/(wnew.x*wnew.x+wnew.b*wnew.b+1.e-300);
            event("exit_complete", t, olda, oldphi-(fall_phi+fall_width), t + dt, a, y[PHI]-(fall_phi+fall_width));
        }
        const auto gas=charged_gas(y);
        const double margin=gas.intrinsic_random-.9*wnew.e;
        if(oldphi>=fall_phi+fall_width && oldh>0. && out.count("bounce_tau")) {
            const auto value=[](double x,double y,double z){return x+(y-x)*z;};
            auto store_budget=[&](const std::string& point,double fraction){
                const std::array<int,7> index{DAUGHTER_REST,DECAY_ENERGY,BIRTH_KINETIC,FIELD_EXPOSURE,ELECTRIC_WORK,GAS_PRESS_WORK,FORMATION_PHASE};
                const std::array<const char*,7> key{"rest","birth_energy","birth_kinetic","exposure","electric_work","gas_pressure_work","formation_phase"};
                for(size_t j=0;j<key.size();++j)out[point+"_"+key[j]]=value(old_sources[j],y[index[j]],fraction);
            };
            if(h>0. && y[PHI]>=fall_phi+fall_width && gas.intrinsic_random/wnew.e>out["peak_hot_fraction"]) {
                out["peak_hot_fraction"]=gas.intrinsic_random/wnew.e;out["peak_hot_time"]=t+dt;
                out["peak_hot_radius"]=a;out["peak_hot_theta"]=gas.theta;out["peak_hot_pressure_ratio"]=gas.pressure/gas.energy;
                out["peak_electric_work"]=y[ELECTRIC_WORK];out["peak_birth_energy"]=y[DECAY_ENERGY];
                out["peak_gas_rest"]=y[DAUGHTER_REST];out["peak_pair_count"]=y[PAIR_COUNT];
            }
            double limit=1.;
            if(h<=0.) limit=std::min(limit,oldh/(oldh-h));
            if(y[PHI]<fall_phi+fall_width) limit=std::min(limit,(oldphi-fall_phi-fall_width)/(oldphi-y[PHI]));
            const double last_margin=value(old_margin,margin,limit);
            if(hot_start<0. && (old_margin>=0. || last_margin>=0.)) {
                const double f=old_margin>=0. ? 0. : old_margin/(old_margin-margin);
                hot_start=value(olda,a,f);
                out["hot_intervals"]+=1.;
                if(!out.count("first_hot_time")) {
                    store_budget("first_hot",f);
                    out["first_hot_time"]=t+dt*f;out["first_hot_radius"]=hot_start;out["first_hot_aH"]=hot_start*value(oldh,h,f);
                }
            }
            if(hot_start>0.) {
                if(last_margin<0.) limit=std::min(limit,old_margin/(old_margin-margin));
                const double end_radius=value(olda,a,limit);
                out["hot_efolds"]=std::max(out["hot_efolds"],std::log(end_radius/hot_start));
                if(olda<2.*hot_start && end_radius>=2.*hot_start && !out.count("doubling_time")){
                    const double f=(2.*hot_start-olda)/(a-olda);
                    store_budget("doubling",f);
                    out["doubling_time"]=t+f*dt;out["doubling_radius"]=2.*hot_start;
                    out["doubling_hot_fraction"]=value(previous_gas.intrinsic_random,gas.intrinsic_random,f)/value(previous_wall.e,wnew.e,f);
                    out["doubling_aH"]=2.*hot_start*value(oldh,h,f);
                    out["doubling_theta"]=value(previous_gas.theta,gas.theta,f);
                    out["doubling_w"]=value(previous_gas.pressure,gas.pressure,f)/value(previous_gas.energy,gas.energy,f);
                }
                if(margin<0. || h<=0. || y[PHI]<fall_phi+fall_width) hot_start=-1.;
            }
        }
        if (oldh < 0 && h >= 0)
            event("bounce", t, olda, oldh, t + dt, a, h);
        if (oldh > 0 && h <= 0 && t > .001)
            event("maximum", t, olda, oldh, t + dt, a, h);
        for (auto q : {std::pair<const char *, double>{"plateau", .012}, {"exit", fall_phi}})
            if (oldphi < q.second && y[PHI] >= q.second)
                event(q.first, t, olda, oldphi - q.second, t + dt, a, y[PHI] - q.second);
        t += dt;
        if ((out.count("maximum_tau") && h < -.005) || a < 1.6 || a > radius_limit)
            break;
    }
    if (want_samples)
        sample(t, y);
    Diag z;
    s.rhs(y, dy, &z);
    out["stop_turnaround"]=out.count("maximum_tau") && y[H]<-.005;
    out["stop_radius_limit"]=y[A_W]>radius_limit;
    out["stop_contraction_limit"]=y[A_W]<1.6;
    const auto osc=oscillator(y);
    const auto gas=charged_gas(y);
    out["end_field_exposure"]=y[FIELD_EXPOSURE];out["end_formation_phase"]=y[FORMATION_PHASE];
    out["end_optical_depth"]=y[OPTICAL_DEPTH];out["end_birth_chi2_over_K"]=y[BIRTH_CHI2_OVER_K];
    out["end_birth_amplitude2_over_K"]=y[BIRTH_AMPLITUDE2_OVER_K];
    out["end_gas_density"]=gas.energy;out["end_gas_random"]=gas.intrinsic_random;out["end_rest"]=y[DAUGHTER_REST];
    out["end_decay_energy"]=y[DECAY_ENERGY];out["end_birth_kinetic"]=y[BIRTH_KINETIC];out["end_pair_count"]=y[PAIR_COUNT];
    out["end_electric_work"]=y[ELECTRIC_WORK];out["end_gas_pressure_work"]=y[GAS_PRESS_WORK];out["end_electric_work_a4"]=y[ELECTRIC_WORK_A4];
    out["end_parent_osc_density"]=osc.energy;
    out["end_parent_osc_pressure"]=osc.pressure;
    out["end_osc_action"]=std::pow(y[A_W],3)*osc.energy/osc.mass;
    out["end_chi_drive_work"]=y[CHI_DRIVE_WORK];out["end_chi_mass_work"]=y[CHI_MASS_WORK];out["end_chi_pressure_work"]=y[CHI_PRESS_WORK];
    out["end_driver_em_work"]=y[PSI_EM_WORK];
    out["end_driver_pressure_work"]=y[PSI_PRESS_WORK];
    out["end_tau"] = t;
    out["end_phi"]=y[PHI]; out["end_pi"]=y[PI];
    out["initial_bulk_mass"]=4.;
    out["tube_proper_time"] = y[SU];
    out["end_radius"] = y[A_W];
    out["end_H"] = y[H];
    out["end_mu"] = y[MASS];
    out["end_sheet_energy_density"] = z.sh;
    out["end_scalar_velocity"] = y[PI];
    out["end_scalar"] = y[PHI];
    out["end_density"] = z.ew;
    out["end_radiation"] = y[RAD];
    out["tube_incident_proper_energy"] = y[EMITTED];
    out["tube_returned_proper_energy"] = y[RETURNED];
    out["tube_returned_after_proper_time_0p1"] = y[LATE_RETURN];
    out["wall_incident_proper_energy"] = y[SHELL_IN];
    out["wall_outgoing_proper_energy"] = y[SHELL_OUT];
    out["wall_pressure_work"] = y[SHELL_WORK];
    out["fluid_comoving_radiation_produced"] = y[RADIATION_COMOV];
    out["comoving_wave_to_sheet"] = y[WAVE_TO_SHEET];
    out["comoving_scalar_to_sheet"] = y[SCALAR_TO_SHEET];
    out["comoving_exit_coupling_transfer"] = y[LATE_SCALAR_TO_SHEET];
    out["end_comoving_sheet_energy"] = std::pow(y[A_W], 4) * z.sh;
    out["end_tube_radius"] = z.r;
    out["end_tube_mass"] = z.mass;
    out["steps"] = steps;
    for (auto &[key, val] : out)
        std::cout << key << ' ' << std::setprecision(16) << val << '\n';
}
