// Forward initial-boundary evolution on a Gaussian-normal annulus. Each element
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
    auto rise = step((psi - .004) / .008), fall = step((psi - .36) / .02);
    if (late_gauge >= 0) {
        k += late_gauge * fall[0];
        dk += late_gauge * fall[1] / .02;
    }
    double v = height * (rise[0] - fall[0]), dv = height * (rise[1] / .008 - fall[1] / .02),
           gamma = late_gauge >= 0 ? 0. : 4 * fall[0];
    double b = 2 * y[C] / a, x = y[MOM] / k, o = x - i, sheet = .5 * k * (x * x + b * b);
    double e = .5 * pi * pi + v + y[RAD] + sheet, p = .5 * pi * pi - v + (y[RAD] + sheet) / 3;
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
            late_gauge >= 0 ? late_gauge * fall[1] / .02 : 0.};
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
        dy[PHI] = pi;
        dy[PI] = .5 * w.dk * (w.x * w.x - w.b * w.b) - w.dv - (3 * h + w.gamma) * pi;
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
    if (argc != 12)
        return 1;
    int cells = std::stoi(argv[1]);
    double depth = std::stod(argv[2]), dose = std::stod(argv[3]), end = std::stod(argv[4]);
    double norm = std::stod(argv[5]), stop = std::stod(argv[6]), cfl = std::stod(argv[7]);
    int late_cells = std::stoi(argv[8]);
    double transfer_time = std::stod(argv[9]);
    late_gauge = std::stod(argv[10]);
    if (cells < 1 || late_cells < 1 || late_cells > cells || depth <= 0 || cfl <= 0)
        throw std::runtime_error("Invalid mesh");
    Solver s(cells, depth, dose, norm, stop);
    int n = COUNT + 3 * cells * P;
    vector<double> y(n, 0.), dy(n), stage(n), acc(n);
    y[A_W] = 2;
    y[MASS] = 4;
    y[INN_MASS] = 4;
    Diag initial;
    s.rhs(y, dy, &initial);
    y[VOLUME] = initial.r;
    double t = 0., dt0 = cfl * s.dx;
    bool want_samples = std::stoi(argv[11]) != 0;
    vector<double> times = {0., .001, .002, .003, .004, .005, .006, .008, .01, .02, .05};
    for (int j = 1; j <= static_cast<int>(10 * end); j++)
        times.push_back(j / 10.);
    size_t next_sample = 0;
    auto sample = [&](double time, const vector<double> &state) {
        Values mat = wall(state, state[s.index(0, I)]);
        double potential =
            height * (step((state[PHI] - .004) / .008)[0] - step((state[PHI] - .36) / .02)[0]);
        std::cout << "@ " << std::setprecision(16) << time << ' ' << state[A_W] << ' ' << state[H] << ' '
                  << state[PHI] << ' ' << state[PI] << ' ' << state[RAD] << ' ' << mat.sheet << ' '
                  << mat.e << ' ' << potential << ' ' << state[SU] << ' ' << state[EMITTED] << ' '
                  << state[RETURNED] << ' ' << state[VOLUME] << ' ' << mat.k << ' ' << mat.b << ' '
                  << mat.x << ' ' << state[MASS] << '\n';
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
                 std::pow(y[A_W], 3) * z.ew + y[SHELL_WORK] - y[SHELL_IN] + y[SHELL_OUT]);
            peak("sheet_balance_error",
                 std::pow(y[A_W], 4) * z.sh - y[WAVE_TO_SHEET] - y[SCALAR_TO_SHEET]);
            peak("fluid_balance_error", std::pow(y[A_W], 4) * y[RAD] - y[RADIATION_COMOV]);
            if (y[H] > 0 && z.ew > 1e-12)
                peak("peak_expanding_radiation_share", y[RAD] / z.ew);
            if (!std::isfinite(z.mass) || z.minn <= 0)
                throw std::runtime_error("Gaussian chart failed");
        }
        double olda = y[A_W], oldh = y[H], oldphi = y[PHI], oldrad = y[RAD];
        Values wold = wall(y, y[s.index(0, I)]);
        double old_e = wold.e;
        double dt = std::min(dt0, end - t);
        if (!transferred && cells != late_cells && t < transfer_time - 1e-14)
            dt = std::min(dt, transfer_time - t);
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
        if (oldphi < .36 && y[PHI] >= .36 && !out.count("first_exit_Bmag")) {
            out["first_exit_Bmag"] = wnew.b;
            out["first_exit_X"] = wnew.x;
            out["first_exit_target_K"] = height / (.5 * wnew.b * wnew.b);
        }
        if (h > 0 && out.count("bounce_tau")) {
            out["maximum_psi_after_bounce"] = std::max(out["maximum_psi_after_bounce"], y[PHI]);
            if (y[PHI] >= .36 && e > 1e-12) {
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
        if (em_start < 0 && f0 < 0 && f1 >= 0 && h > 0 && y[PHI] >= .36 && out.count("bounce_tau")) {
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
        if (oldphi < .38 && y[PHI] >= .38)
            event("exit_complete", t, olda, oldphi - .38, t + dt, a, y[PHI] - .38);
        if (oldh < 0 && h >= 0)
            event("bounce", t, olda, oldh, t + dt, a, h);
        if (oldh > 0 && h <= 0 && t > .001)
            event("maximum", t, olda, oldh, t + dt, a, h);
        for (auto q : {std::pair<const char *, double>{"plateau", .012}, {"exit", .36}})
            if (oldphi < q.second && y[PHI] >= q.second)
                event(q.first, t, olda, oldphi - q.second, t + dt, a, y[PHI] - q.second);
        if (h > 0 && oldrad < .9 * old_e && y[RAD] >= .9 * e && y[RAD] > 1e-10)
            event("rad90", t, olda, oldrad - .9 * old_e, t + dt, a, y[RAD] - .9 * e);
        if (h > 0 && y[RAD] >= .9 * e && y[RAD] > 1e-10)
            out["last_rad90_expanding_radius"] = a;
        t += dt;
        if ((out.count("maximum_tau") && h < -.005) || a < 1.6 || a > 16.)
            break;
    }
    if (want_samples)
        sample(t, y);
    Diag z;
    s.rhs(y, dy, &z);
    out["end_tau"] = t;
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
