// One row per source cell gives its contribution to the two physical fields at birth.
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

struct Map {
    double ee, eb, be, bb;
};

Map propagate(double dt, double a, double K, double k) {
    // Matrix exponential of the two source-free canonical wall equations,
    // evaluated with midpoint coefficients. The critical-damping limit is regular.
    const double d = dt / K;
    const double v = dt * k / a;
    const double s = d*d - v*v;
    double c, f;
    if (std::abs(s) < 1.e-8) {
        c = 1. + s*(.5+s/24.);
        f = 1. + s*(1./6.+s/120.);
    } else if (s > 0.) {
        const double r = std::sqrt(s);
        c = std::cosh(r);
        f = std::sinh(r)/r;
    } else {
        const double r = std::sqrt(-s);
        c = std::cos(r);
        f = std::sin(r)/r;
    }
    const double decay = std::exp(-d);
    return {decay*(c-d*f), -decay*v*K*f, decay*v*f/K, decay*(c+d*f)};
}

Map compose(Map l, Map r) {
    return {l.ee*r.ee+l.eb*r.be, l.ee*r.eb+l.eb*r.bb,
            l.be*r.ee+l.bb*r.be, l.be*r.eb+l.bb*r.bb};
}

int main(int argc, char** argv) {
    if (argc != 4) return 2;
    const double k = std::stod(argv[1]);
    std::ifstream input(argv[2], std::ios::binary);
    std::ofstream output(argv[3], std::ios::binary);
    if (!input || !output || !(k > 0.)) throw std::runtime_error("invalid job");
    using Row = std::array<double, 3>;  // time, radius, gauge coefficient
    std::vector<Row> table;
    Row row;
    while (input.read(reinterpret_cast<char*>(row.data()), sizeof row)) table.push_back(row);
    if (!input.eof() || input.gcount() || table.size()<2) throw std::runtime_error("bad trace");
    const double a_end=table.back()[1], K_end=table.back()[2];
    Map left{1./(a_end*a_end*K_end), 0., 0., 1./(a_end*a_end)};
    std::vector<std::array<double,2>> columns(table.size()-1);
    for (size_t j=columns.size(); j-- > 0;) {
        const auto& lo=table[j];
        const auto& hi=table[j+1];
        const double dt=hi[0]-lo[0], a=.5*(hi[1]+lo[1]), K=.5*(lo[2]+hi[2]);
        if (dt <= 0. || a <= 0. || K <= 0.) throw std::runtime_error("bad coefficients");
        const Map half=propagate(.5*dt,a,K,k);
        const Map mid=compose(left,half);
        columns[j]={4.*dt*a*a*mid.ee, 4.*dt*a*a*mid.be};
        left=compose(mid,half);
    }
    output.write(reinterpret_cast<const char*>(columns.data()), columns.size()*sizeof(columns[0]));
    return !output;
}
