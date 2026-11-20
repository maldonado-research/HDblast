// Parity-even two-stage release on a finite Gaussian-normal annulus. Each element
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
// Rest energy density summed over the six streams at a=2, q/(sqrt(3)*m), drag rate.
double rest_at_start = 0., charge_mass = 0., collision_rate = 0.;
double threshold_phi=.77, field_width=.03, mass_low=.95, mass_high=1.05, width_scale=100.;
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
    FRACTION, GAS_MOM, LAB_EXCESS, ELECTRIC_WORK, ELECTRIC_WORK_A4, GAS_PRESSURE_WORK,
    DECAY_ENERGY, DECAY_KINETIC, MASS_WORK, DRAG_HEAT, MIXING_HEAT, OPTICAL_DEPTH,
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
// Each opposite-sign pair has opposite velocity on the same spatial axis;
// three identical copies give exactly isotropic total stress. The spatial
// velocity component in units of c is u/sqrt(1+u^2).
struct Gas {
    double f, dz, parent_mass, mass_slope, rate, parent_density, stable_mass_density;
    double cold, excess, lab_pressure, gas_momentum, energy, pressure, velocity, gamma, theta;
    double current, random_lab, intrinsic_random, rest, u, reconstruction_error;
};

Gas charged_gas(const vector<double> &y) {
    Gas g{};
    const auto mass_step=step((y[PHI]-threshold_phi)/field_width+.5);
    g.parent_mass=mass_low+(mass_high-mass_low)*mass_step[0]; // M_chi / (2m)
    g.mass_slope=(mass_high-mass_low)*mass_step[1]/field_width;
    double beta_squared=std::max(0.,1.-1./(g.parent_mass*g.parent_mass));
    g.rate=width_scale*g.parent_mass*beta_squared*std::sqrt(beta_squared);
    g.f=y[FRACTION];
    g.dz=g.rate*(1.-g.f);
    g.stable_mass_density=rest_at_start/mass_low*std::pow(2./y[A_W],3);
    g.parent_density=g.stable_mass_density*g.parent_mass*(1.-g.f);
    g.gamma=1.;
    if(g.f<=0.) return g;
    g.rest=g.stable_mass_density*g.f;
    const double moment=y[GAS_MOM], extra=y[LAB_EXCESS];
    const double cold=std::hypot(g.f,moment), E=cold+extra;
    g.cold=cold; g.excess=extra; g.gas_momentum=moment;
    // Invert (h-Theta)(h-4Theta)=1 for nonnegative rest-frame pressure P.
    // Normalize by E. Compute E^2-Z^2 from U to retain the cold limit.
    const double z=g.f/E, zz=z*z, surplus=extra/E*((2.*cold+extra)/E);
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
    g.lab_pressure=E*p; // Rest-frame pressure divided by B, not isotropic wall pressure.
    g.velocity=(moment/E)/(1.+p);
    g.gamma=(1.+p)/std::sqrt(A+p*(2.+p));
    g.u=g.gamma*g.velocity;
    g.theta=p*g.gamma/z;
    double root=std::hypot(1.,1.5*g.theta);
    double thermal=1.5*g.theta+(2.25*g.theta*g.theta)/(1.+root);
    g.energy=g.stable_mass_density*E;
    g.pressure=g.stable_mass_density*(g.lab_pressure+moment*g.velocity/3.);
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
    auto gas=charged_gas(y);
    e+=gas.parent_density+gas.energy;
    p+=gas.pressure;
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
        const auto gas=charged_gas(y);
        dy[MOM] -= gas.current;
        // Stream groups receive equal and opposite daughter charges with zero mean
        // injected momentum. Rest-frame daughter energy per birth is M/(2m).
        const double moment=y[GAS_MOM], a3=a*a*a;
        const double vc=gas.cold>0. ? moment/gas.cold : 0.;
        const double inv_gamma_c=gas.cold>0. ? gas.f/gas.cold : 1.;
        const double dv=gas.cold>0. ? -vc*(gas.excess+gas.lab_pressure)/
                            (gas.cold+gas.excess+gas.lab_pressure) : 0.;
        dy[FRACTION]=gas.dz;
        dy[GAS_MOM]=charge_mass*gas.f*w.x-(h+collision_rate)*moment;
        dy[LAB_EXCESS]=charge_mass*gas.f*w.x*dv+h*(-moment*dv-3.*gas.lab_pressure)
                       +collision_rate*moment*vc+(gas.parent_mass-inv_gamma_c)*gas.dz;
        dy[ELECTRIC_WORK]=a3*gas.current*w.x;
        dy[ELECTRIC_WORK_A4]=a*dy[ELECTRIC_WORK];
        dy[GAS_PRESSURE_WORK]=3.*a3*h*gas.pressure;
        dy[DECAY_ENERGY]=a3*gas.stable_mass_density*gas.parent_mass*gas.dz;
        dy[DECAY_KINETIC]=a3*gas.stable_mass_density*(gas.parent_mass-1.)*gas.dz;
        dy[MASS_WORK]=a3*gas.stable_mass_density*(1.-gas.f)*gas.mass_slope*pi;
        dy[DRAG_HEAT]=a3*gas.stable_mass_density*collision_rate*gas.u*moment;
        dy[MIXING_HEAT]=a3*gas.stable_mass_density*(gas.gamma-1.)*gas.parent_mass*gas.dz;
        dy[OPTICAL_DEPTH]=gas.rate;
        dy[PHI] = pi;
        dy[PI] = .5 * w.dk * (w.x * w.x - w.b * w.b) - w.dv - (3 * h + w.gamma) * pi
                 -gas.stable_mass_density*(1.-gas.f)*gas.mass_slope;
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
    if (argc != 27) {
        std::cerr << "Run with growth.py or pass the 26 parameters in its stated order.\n";
        return 1;
    }
    int cells = std::stoi(argv[1]);
    double depth = std::stod(argv[2]), dose = std::stod(argv[3]), end = std::stod(argv[4]);
    double norm = std::stod(argv[5]), stop = std::stod(argv[6]), cfl = std::stod(argv[7]);
    int late_cells = std::stoi(argv[8]);
    double transfer_time = std::stod(argv[9]);
    late_gauge = std::stod(argv[10]);
    rest_at_start=std::stod(argv[12]);
    charge_mass=std::stod(argv[13]);
    collision_rate=std::stod(argv[14]);
    threshold_phi=std::stod(argv[15]);
    field_width=std::stod(argv[16]);
    width_scale=std::stod(argv[17]);
    mass_low=std::stod(argv[18]);
    mass_high=std::stod(argv[19]);
    fall_phi=std::stod(argv[20]); fall_width=std::stod(argv[21]);
    prepare_phi=std::stod(argv[22]); prepare_width=std::stod(argv[23]);
    prepare_gauge=std::stod(argv[24]); prepare_height=std::stod(argv[25]); radius_limit=std::stod(argv[26]);
    if (fall_phi<=.012 || fall_width<=0. || prepare_width<=0. || prepare_gauge<0. || prepare_height<0. || prepare_height>=height ||
        ((prepare_gauge>0. || prepare_height>0.) && (prepare_phi<.012 || prepare_phi+prepare_width>fall_phi)) || radius_limit<=2.)
        throw std::runtime_error("Invalid staged scalar or gauge interaction");
    if (!(threshold_phi>=field_width/2. && field_width>0. && width_scale>=0. &&
          mass_low>0. && mass_low<=1. && mass_high>=mass_low))
        throw std::runtime_error("Invalid decay parameters");
    if (rest_at_start<0. || charge_mass<0. || collision_rate<0. || late_gauge<0.)
        throw std::runtime_error("Invalid parent, daughter or gauge parameter");
    if (cells < 1 || late_cells < 1 || late_cells > cells || depth <= 0 || cfl <= 0)
        throw std::runtime_error("Invalid mesh");
    Solver s(cells, depth, dose, norm, stop);
    int n = COUNT + 3 * cells * P;
    vector<double> y(n, 0.), dy(n), stage(n), acc(n);
    y[A_W] = 2;
    // Choose mu using the exact initial Friedmann constraint at a=2,H=0.
    y[MASS] = 4.-16.*rest_at_start*(2.+rest_at_start);
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
        // The primary trace has dense coverage of decay and first turnaround.
        int samples_per_time=(sampling==2?1000:50);
        for(int i=4*samples_per_time;i<=13*samples_per_time;i++) times.push_back(double(i)/samples_per_time);
        std::sort(times.begin(),times.end());
        times.erase(std::unique(times.begin(),times.end(),[](double a,double b){return std::abs(a-b)<1e-9;}),times.end());
    }
    size_t next_sample = 0;
    auto sample = [&](double time, const vector<double> &state) {
        Values mat = wall(state, state[s.index(0, I)]);
        auto gas=charged_gas(state);
        double potential =
            height*step((state[PHI]-.004)/.008)[0]-prepare_height*step((state[PHI]-prepare_phi)/prepare_width)[0]-(height-prepare_height)*step((state[PHI]-fall_phi)/fall_width)[0];
        std::cout << "@ " << std::setprecision(16) << time << ' ' << state[A_W] << ' ' << state[H] << ' '
                  << state[PHI] << ' ' << state[PI] << ' ' << state[RAD] << ' ' << mat.sheet << ' '
                  << mat.e << ' ' << potential << ' ' << state[SU] << ' ' << state[EMITTED] << ' '
                  << state[RETURNED] << ' ' << state[VOLUME] << ' ' << mat.k << ' ' << mat.b << ' '
                  << mat.x << ' ' << state[MASS] << ' ' << gas.rest << ' ' << gas.u << ' ' << gas.theta
                  << ' ' << gas.energy << ' ' << gas.pressure << ' ' << gas.random_lab
                  << ' ' << gas.intrinsic_random << ' ' << state[ELECTRIC_WORK] << ' ' << state[DRAG_HEAT] << ' ' << gas.f
                  << ' ' << gas.parent_mass << ' ' << gas.parent_density << ' ' << gas.rate
                  << ' ' << state[DECAY_ENERGY] << ' ' << state[DECAY_KINETIC] << ' ' << state[MASS_WORK]
                  << ' ' << state[MIXING_HEAT] << ' ' << state[OPTICAL_DEPTH] << ' ' << gas.current
                  << ' ' << mat.p << ' ' << state[WAVE_TO_SHEET] << ' ' << state[SCALAR_TO_SHEET]
                  << ' ' << state[ELECTRIC_WORK_A4] << ' ' << state[GAS_PRESSURE_WORK]
                  << ' ' << state[GAS_MOM] << ' ' << state[LAB_EXCESS] << '\n';
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
    double em_start = -1.;
    int em_intervals = 0;
    std::array<double,2> intrinsic_start{-1.,-1.};
    const std::array<double,2> thresholds{.5,.9};
    const std::array<std::string,2> labels{"half","ninety"};
    int steps = 0;
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
            peak("surface_balance_error",
                 std::pow(y[A_W], 3) * z.ew + y[SHELL_WORK] - y[SHELL_IN] + y[SHELL_OUT] - 8.*rest_at_start);
            peak("sheet_balance_error",
                 std::pow(y[A_W], 4) * z.sh - y[WAVE_TO_SHEET] - y[SCALAR_TO_SHEET] + y[ELECTRIC_WORK_A4]);
            peak("fluid_balance_error", std::pow(y[A_W], 4) * y[RAD] - y[RADIATION_COMOV]);
            const auto gas=charged_gas(y);
            const double a3=std::pow(y[A_W],3);
            peak("gas_balance_error",a3*gas.energy+y[GAS_PRESSURE_WORK]-y[ELECTRIC_WORK]-y[DECAY_ENERGY]);
            peak("parent_balance_error",a3*gas.parent_density-8.*rest_at_start-y[MASS_WORK]+y[DECAY_ENERGY]);
            peak("particle_number_error",gas.f+std::expm1(-y[OPTICAL_DEPTH]));
            peak("daughter_injected_rest_error",y[DECAY_ENERGY]-y[DECAY_KINETIC]-a3*gas.rest);
            peak("gas_primitive_error",gas.reconstruction_error);
            if(gas.f< -1.e-11 || gas.f>1.+1.e-10 || y[LAB_EXCESS]< -1.e-12)
                throw std::runtime_error("Unphysical conservative gas state");
            if(gas.theta < -1e-8 || std::abs(gas.velocity)>1. || !std::isfinite(gas.theta))
                throw std::runtime_error("Nonphysical gas state");
            if (y[H] > 0 && z.ew > 1e-12)
                peak("peak_expanding_radiation_share", y[RAD] / z.ew);
            if (!std::isfinite(z.mass) || z.minn <= 0)
                throw std::runtime_error("Gaussian chart failed");
        }
        double olda = y[A_W], oldh = y[H], oldphi = y[PHI], oldrad = y[RAD];
        Values wold = wall(y, y[s.index(0, I)]);
        double old_e = wold.e;
        const auto prev_gas=charged_gas(y);
        const double old_intrinsic=prev_gas.intrinsic_random;
        const double old_charge_work=y[ELECTRIC_WORK], old_drag=y[DRAG_HEAT];
        double dt = std::min(dt0, end - t);
        if (!transferred && cells != late_cells && t < transfer_time - 1e-14)
            dt = std::min(dt, transfer_time - t);
        // Resolve rapid initial acceleration after a very large charge is activated.
        if (charge_mass>0. && prev_gas.f>0.) {
            // The cold carrier plasma can be stiff even at a zero of its current.
            const double enthalpy=2.5*prev_gas.theta+std::hypot(1.,1.5*prev_gas.theta);
            const double frequency=charge_mass*std::sqrt(prev_gas.rest/(wold.k*enthalpy*prev_gas.gamma));
            if (frequency>0.) dt=std::min(dt,.5/frequency);
            const double size=std::max(1.e-16,prev_gas.f);
            if (dy[GAS_MOM]!=0.) dt=std::min(dt,.05*(size+std::abs(y[GAS_MOM]))/std::abs(dy[GAS_MOM]));
            if (dy[LAB_EXCESS]!=0.) dt=std::min(dt,.05*(size+std::abs(y[LAB_EXCESS]))/std::abs(dy[LAB_EXCESS]));
        }
        if (prev_gas.rate>0.) dt=std::min(dt,.20/prev_gas.rate);
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
        const double a = y[A_W], h = y[H], e = wnew.e;
        const auto now_gas=charged_gas(y);
        if (prev_gas.parent_mass<=1. && now_gas.parent_mass>1. && !out.count("threshold_tau")) {
            const double z=(1.-prev_gas.parent_mass)/(now_gas.parent_mass-prev_gas.parent_mass);
            double aevent=olda+(a-olda)*z;
            event("threshold",t,olda,prev_gas.parent_mass-1.,t+dt,a,now_gas.parent_mass-1.);
            out["threshold_K"]=wold.k+(wnew.k-wold.k)*z;
            out["threshold_electric_share"]=(wold.x*wold.x)/(wold.x*wold.x+wold.b*wold.b+1.e-300);
            out["threshold_em_a4"]=std::pow(aevent,4)*(wold.sheet+z*(wnew.sheet-wold.sheet));
            out["threshold_bulk_sheet_flux"]=y[WAVE_TO_SHEET];
            out["threshold_scalar_sheet_flux"]=y[SCALAR_TO_SHEET];
        }
        if(prev_gas.f<.5 && now_gas.f>=.5) event("half_decayed",t,olda,prev_gas.f-.5,t+dt,a,now_gas.f-.5);
        if(prev_gas.f<.9 && now_gas.f>=.9) event("ninety_decayed",t,olda,prev_gas.f-.9,t+dt,a,now_gas.f-.9);
        if(oldphi>=fall_phi+fall_width && oldh>0. && wold.b*wnew.b<0. && !out.count("postexit_magnetic_zero_tau"))
            event("postexit_magnetic_zero", t, olda, wold.b, t+dt, a, wnew.b);
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
        if (h > 0 && out.count("bounce_tau")) {
            out["maximum_psi_after_bounce"] = std::max(out["maximum_psi_after_bounce"], y[PHI]);
            const auto gas=charged_gas(y);
            if (gas.theta>out["maximum_temperature_in_expansion"]) out["maximum_temperature_in_expansion"]=gas.theta;
            if (gas.random_lab/e>out["peak_random_fraction_in_expansion"]) {
                out["peak_random_fraction_in_expansion"]=gas.random_lab/e;
                out["random_peak_radius"]=a;
                out["random_peak_time"]=t+dt;
                out["random_peak_theta"]=gas.theta;
                out["random_peak_velocity"]=gas.velocity;
                out["random_peak_intrinsic_fraction"]=gas.intrinsic_random/e;
                out["random_peak_comoving_work"]=y[ELECTRIC_WORK];
            }
            if (y[PHI] >= fall_phi+fall_width && e > 1e-12) {
                if(gas.intrinsic_random/e>out["maximum_intrinsic_fraction"]) {
                    out["peak_wall_density"]=e;
                    out["peak_daughter_rest_density"]=gas.rest;
                    out["peak_parent_density"]=gas.parent_density;
                    out["maximum_intrinsic_fraction"]=gas.intrinsic_random/e;
                    out["maximum_intrinsic_tau"]=t+dt;
                    out["maximum_intrinsic_radius"]=a;
                    out["peak_daughter_fraction"]=gas.f;
                    out["peak_daughter_theta"]=gas.theta;
                    out["peak_daughter_velocity"]=gas.velocity;
                    out["peak_decay_energy_a3"]=y[DECAY_ENERGY];
                    out["peak_injected_kinetic_a3"]=y[DECAY_KINETIC];
                    out["peak_mass_work_a3"]=y[MASS_WORK];
                    out["peak_daughter_work_a3"]=y[ELECTRIC_WORK];
                    out["peak_gas_pressure_work"]=y[GAS_PRESSURE_WORK];
                    out["peak_gas_total_a3"]=std::pow(a,3)*gas.energy;
                    out["peak_postexit_em_a4"]=std::pow(a,4)*wnew.sheet;
                    out["peak_work_a4"]=y[ELECTRIC_WORK_A4];
                    out["peak_drag_restframe_source"]=y[DRAG_HEAT];
                    out["peak_mixing_restframe_source"]=y[MIXING_HEAT];
                }
                const double random_fraction=gas.random_lab/e;
                if(random_fraction>out["postexit_peak_random_fraction"]){
                    out["postexit_peak_random_fraction"]=random_fraction;
                    out["postexit_peak_radius"]=a;
                    out["postexit_peak_time"]=t+dt;
                    out["postexit_peak_theta"]=gas.theta;
                    out["postexit_peak_velocity"]=gas.velocity;
                    out["postexit_peak_intrinsic_fraction"]=gas.intrinsic_random/e;
                    out["postexit_peak_sheet_fraction"]=wnew.sheet/e;
                    out["postexit_peak_gas_pressure_to_density"]=gas.pressure/gas.energy;
                    out["postexit_peak_work_a4"]=y[ELECTRIC_WORK_A4];
                    out["postexit_peak_bulk_sheet_flux"]=y[WAVE_TO_SHEET];
                    out["postexit_peak_scalar_sheet_flux"]=y[SCALAR_TO_SHEET];
                    out["postexit_peak_em_a4"]=std::pow(a,4)*wnew.sheet;
                }
                for (size_t j=0;j<thresholds.size();++j) {
                    const double f_old=old_intrinsic-thresholds[j]*old_e;
                    const double f_new=gas.intrinsic_random-thresholds[j]*e;
                    const std::string &label=labels[j];
                    if(f_new>=0. && f_old<0.) {
                        const double frac=-f_old/(f_new-f_old);
                        intrinsic_start[j]=olda+(a-olda)*frac;
                        out[label+"_intrinsic_intervals"] += 1;
                        if(!out.count("first_"+label+"_intrinsic_radius")) {
                            out["first_"+label+"_intrinsic_radius"]=intrinsic_start[j];
                            out["first_"+label+"_intrinsic_time"]=t+dt*frac;
                            out["first_"+label+"_intrinsic_aH"]=intrinsic_start[j]*(oldh+frac*(h-oldh));
                        }
                    }
                    if(intrinsic_start[j]>0.) {
                        double limit=a;
                        if(f_old>=0. && f_new<0.) limit=olda-(a-olda)*f_old/(f_new-f_old);
                        out[label+"_intrinsic_expanding_efolds"]=std::max(out[label+"_intrinsic_expanding_efolds"],std::log(limit/intrinsic_start[j]));
                        out["last_"+label+"_intrinsic_radius"]=limit;
                        if(f_new<0.) intrinsic_start[j]=-1.;
                    }
                }
                if (out.count("first_ninety_intrinsic_radius") && !out.count("first_hot_doubling_tau") &&
                    intrinsic_start[1]>0.) {
                    const double target=2.*intrinsic_start[1];
                    if(olda<target && a>=target) {
                        const double f=(target-olda)/(a-olda);
                        const auto interp=[&](double v0,double v1){return v0+f*(v1-v0);};
                        out["first_hot_doubling_tau"]=t+f*dt;
                        out["first_hot_doubling_radius"]=target;
                        out["doubling_aH"]=target*interp(oldh,h);
                        out["doubling_intrinsic_fraction"]=interp(old_intrinsic,gas.intrinsic_random)/interp(old_e,e);
                        out["doubling_pressure_over_gas_density"]=interp(prev_gas.pressure,gas.pressure)/interp(prev_gas.energy,gas.energy);
                        out["doubling_theta"]=interp(prev_gas.theta,gas.theta);
                    }
                }
                if(random_fraction > .5) {
                    if(!out.count("postexit_first_half_random_radius")) {
                        out["postexit_first_half_random_radius"]=a;
                        out["postexit_first_half_random_tau"]=t+dt;
                    }
                    out["postexit_last_half_random_radius"]=a;
                }
                if(random_fraction > .9) {
                    if(!out.count("postexit_first_ninety_random_radius")) out["postexit_first_ninety_random_radius"]=a;
                    out["postexit_last_ninety_random_radius"]=a;
                }
            }
            if (y[PHI] >= fall_phi && e > 1e-12) {
                double sh = wnew.sheet / e;
                if (sh > out["peak_exit_sheet_fraction"]) {
                    out["peak_exit_sheet_fraction"] = sh;
                    out["at_peak_sheet_K"] = wnew.k;
                    out["at_peak_sheet_electric_density"] = .5 * wnew.k * wnew.x * wnew.x;
                    out["at_peak_sheet_magnetic_density"] = .5 * wnew.k * wnew.b * wnew.b;
                    out["at_peak_sheet_radius"] = a;
                    out["at_peak_sheet_total_density"] = e;
                    out["at_peak_comoving_wave_to_sheet"] = y[WAVE_TO_SHEET];
                    out["at_peak_comoving_scalar_to_sheet"] = y[SCALAR_TO_SHEET];
                    out["at_peak_comoving_late_scalar_to_sheet"] = y[LATE_SCALAR_TO_SHEET];
                    out["at_peak_comoving_sheet_energy"] = std::pow(a, 4) * wnew.sheet;
                }
            }
        }
        double f0 = wold.sheet - .9 * old_e, f1 = wnew.sheet - .9 * e;
        if (em_start < 0 && f0 < 0 && f1 >= 0 && h > 0 && y[PHI] >= fall_phi && out.count("bounce_tau")) {
            em_start = olda + (a - olda) * (-f0 / (f1 - f0));
            em_intervals++;
            if (em_intervals == 1)
                out["first_em90_radius"] = em_start;
        }
        if (em_start > 0) {
            double abound = a;
            if (f1 < 0 && f0 >= 0)
                abound = olda + (a - olda) * (-f0 / (f1 - f0));
            out["longest_em90_expanding_efolds"] =
                std::max(out["longest_em90_expanding_efolds"], std::log(abound / em_start));
            if (f1 < 0 || h <= 0)
                em_start = -1;
        }
        if (oldphi < fall_phi+fall_width && y[PHI] >= fall_phi+fall_width && !out.count("exit_complete_tau")) {
            const double frac=(fall_phi+fall_width-oldphi)/(y[PHI]-oldphi);
            out["exit_complete_aH"]=a*h;out["exit_complete_field_density"]=wnew.sheet;
            out["exit_complete_electric_fraction"]=wnew.x*wnew.x/(wnew.x*wnew.x+wnew.b*wnew.b+1.e-300);
            out["exit_charge_work"]=old_charge_work+frac*(y[ELECTRIC_WORK]-old_charge_work);
            out["exit_heat_source_integral"]=old_drag+frac*(y[DRAG_HEAT]-old_drag);
            event("exit_complete", t, olda, oldphi-(fall_phi+fall_width), t + dt, a, y[PHI]-(fall_phi+fall_width));
        }
        if (oldh < 0 && h >= 0)
            event("bounce", t, olda, oldh, t + dt, a, h);
        if (oldh > 0 && h <= 0 && t > .001)
            event("maximum", t, olda, oldh, t + dt, a, h);
        for (auto q : {std::pair<const char *, double>{"plateau", .012}, {"exit", fall_phi}})
            if (oldphi < q.second && y[PHI] >= q.second)
                event(q.first, t, olda, oldphi - q.second, t + dt, a, y[PHI] - q.second);
        if (h > 0 && oldrad < .9 * old_e && y[RAD] >= .9 * e && y[RAD] > 1e-10)
            event("rad90", t, olda, oldrad - .9 * old_e, t + dt, a, y[RAD] - .9 * e);
        if (h > 0 && y[RAD] >= .9 * e && y[RAD] > 1e-10)
            out["last_rad90_expanding_radius"] = a;
        t += dt;
        if ((out.count("maximum_tau") && h < -.005) || a < 1.6 || a > radius_limit)
            break;
    }
    if (want_samples)
        sample(t, y);
    Diag z;
    s.rhs(y, dy, &z);
    const auto gas=charged_gas(y);
    out["end_tau"] = t;
    out["end_phi"]=y[PHI]; out["end_pi"]=y[PI];
    out["initial_bulk_mass"]=4.-16.*rest_at_start*(2.+rest_at_start);
    out["end_gas_density"]=gas.energy;
    out["end_temperature"]=gas.theta;
    out["end_daughter_fraction"]=gas.f;
    out["end_optical_depth"]=y[OPTICAL_DEPTH];
    out["end_parent_density"]=gas.parent_density;
    out["end_parent_mass_ratio"]=gas.parent_mass;
    out["end_decay_energy_a3"]=y[DECAY_ENERGY];
    out["end_decay_kinetic_a3"]=y[DECAY_KINETIC];
    out["end_scalar_mass_work_a3"]=y[MASS_WORK];
    out["end_mixing_restframe_source"]=y[MIXING_HEAT];
    out["end_thermal_density"]=gas.random_lab;
    out["end_gas_pressure"]=gas.pressure;
    out["end_intrinsic_random_density"]=gas.intrinsic_random;
    out["end_charge_work"]=y[ELECTRIC_WORK];
    out["end_charge_work_a4"]=y[ELECTRIC_WORK_A4];
    out["end_heat_source_integral"]=y[DRAG_HEAT];
    out["end_gas_pressure_work"]=y[GAS_PRESSURE_WORK];
    out["tube_proper_time"] = y[SU];
    out["end_radius"] = y[A_W];
    out["end_H"] = y[H];
    out["end_mu"] = y[MASS];
    out["end_sheet_energy_density"] = z.sh;
    out["end_scalar_velocity"] = y[PI];
    out["end_scalar"] = y[PHI];
    out["end_density"] = z.ew;
    out["end_radiation"] = y[RAD];
    out["postexit_em90_intervals"] = em_intervals;
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
    if (out.count("rad90_radius"))
        out["radiation_expanding_efolds"] =
            std::log(out.at("last_rad90_expanding_radius") / out.at("rad90_radius"));
    out["steps"] = steps;
    for (auto &[key, val] : out)
        std::cout << key << ' ' << std::setprecision(16) << val << '\n';
}
