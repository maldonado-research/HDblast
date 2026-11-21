// Read archived accepted substeps written by replay.py. The emitted birth law is
// isotropic in u=cos(theta), and within a step the signed field is linear in time.
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

#include "sweep.hpp"

using Real = long double;

struct Angular {
    Real unfilled = 0;    // Mean of the positive directional shortfall divided by demand.
    Real deficient = 0;   // Solid-angle fraction having any shortfall.
    Angular operator+(const Angular& b) const { return {unfilled+b.unfilled, deficient+b.deficient}; }
    Angular operator*(Real x) const { return {unfilled*x, deficient*x}; }
};

// c and d are absolute speeds in units of the per-area request. On the
// interval [c-d,c+d] the missing weight is a triangle with unit half-width.
// Compute the overlap separately on either side of its tip.
Angular at_field(Real c, Real d) {
    if (d == 0) return {std::max(0.L, 1-c), c < 1 ? 1.L : 0.L};
    if (d >= c+1) return {0.5L/d, 1/d};
    if (c-d >= 1) return {};
    if (c >= d && d <= 1-c) return {1-c, 1};
    const Real left = c-d, right = c+d;
    const Real pos_low = std::max(0.L, left), pos_high = std::min(1.L, right);
    const Real neg_width = std::clamp(-left, 0.L, 1.L);
    const Real pos_width = std::max(0.L, pos_high-pos_low);
    const Real area = pos_width*((1-pos_low)+(1-pos_high))/2 + neg_width*(1-neg_width/2);
    return {area/(2*d), (pos_width+neg_width)/(2*d)};
}

constexpr std::array<Real,4> gx = {
    .18343464249564980494L,.52553240991632898582L,.79666647741362673959L,.96028985649753623168L};
constexpr std::array<Real,4> gw = {
    .36268378337836198297L,.31370664587788728734L,.22238103445337447054L,.10122853629037625915L};

// Integrate over nonnegative field magnitude. Split at every triangle endpoint
// and tip event. Outside the last event the answer is proportional to 1/d and
// has an exact log integral, including a very broad outer field range.
Angular positive_integral(Real lo, Real hi, Real c) {
    if (lo == hi) return {};
    if (hi <= c-1) return {};
    auto tail = [](Real a, Real b) -> Angular {
        const Real logarithm = std::log1p((b-a)/a);
        return {logarithm/2, logarithm};
    };
    if (lo >= c+1) return tail(lo,hi);
    Angular integral;
    const Real cplus = c+1;
    if (hi > cplus) {
        integral = tail(cplus,hi);
        hi = cplus;
    }
    std::array<Real,5> knots = {lo,hi,std::abs(c-1),c,cplus};
    std::sort(knots.begin(),knots.end());
    for (unsigned j=1; j<knots.size(); ++j) {
        Real a = std::max(lo,knots[j-1]), b = std::min(hi,knots[j]);
        if (b <= a || b <= c-1) continue;
        // With c=1 the origin is a removable singularity. Otherwise a nonzero
        // knot starts any segment with a logarithmic term. Limit its ratio so
        // low-order Gauss integration retains accuracy down to the inner layer.
        while (a < b) {
            const Real end = a>0 ? std::min(b,2*a) : b;
            const Real half = (end-a)/2, mid = a+half;
            Angular pair;
            for (unsigned i=0;i<gx.size();++i) {
                pair = pair + (at_field(c,mid-half*gx[i]) + at_field(c,mid+half*gx[i]))*gw[i];
            }
            integral = integral + pair*half;
            a = end;
        }
    }
    return integral;
}

Angular field_average(Real c, Real d0, Real d1) {
    if (d0 == d1) return at_field(c,std::abs(d0));
    if (d0>d1) std::swap(d0,d1);
    if (d0>=0) return positive_integral(d0,d1,c)*(1/(d1-d0));
    if (d1<=0) return positive_integral(-d1,-d0,c)*(1/(d1-d0));
    return (positive_integral(0,-d0,c)+positive_integral(0,d1,c))*(1/(d1-d0));
}

void test_close(Real a, Real b) {
    if (std::abs(a-b)>2e-13L) throw std::runtime_error("Analytic angular control failed");
}

void self_test() {
    test_close(at_field(0,2).unfilled,.25L);
    test_close(at_field(1,1).unfilled,.25L);
    test_close(at_field(0,0).unfilled,1);
    test_close(at_field(2,1).unfilled,0);
    test_close(field_average(0,-2,2).unfilled,(.75L+.5L*std::log(2.L))/2);
    test_close(field_average(0,-2,2).deficient,(1+std::log(2.L))/2);
    test_close(field_average(1,-2,2).unfilled,(.75L-.5L*std::log(2.L))/2);
    test_close(field_average(1,-2,2).deficient,.5L);
    test_close(field_average(2,1,2).unfilled,(std::log(2.L)-.5L)/4);
    test_close(field_average(1.001L,-2.e8L,4.e8L).unfilled,
               (field_average(1.001L,-2.e8L,0).unfilled+2*field_average(1.001L,0,4.e8L).unfilled)/3);
    test_close(at_field(.8L,std::sqrt(1.8L*1.8L-2)).unfilled,
               (1.8L-std::sqrt(1.8L*1.8L-2))/2);
}

// Select the field magnitude which maximizes the instantaneous angular mean.
// For c>1/2 the maximum is interior to c<d<c+1 at sqrt((1+c)^2-2).
// Otherwise the maximum occurs on a plateau containing d=0.
Real peak_field(Real c, Real d0, Real d1) {
    const Real optimum = c > .5L ? std::sqrt((c+1)*(c+1)-2) : 0;
    Real chosen = d0;
    for (Real candidate : {d1, std::clamp(optimum,std::min(d0,d1),std::max(d0,d1)),
                                std::clamp(-optimum,std::min(d0,d1),std::max(d0,d1))}) {
        if (at_field(c,std::abs(candidate)).unfilled > at_field(c,std::abs(chosen)).unfilled) chosen=candidate;
    }
    return chosen;
}

int main(int argc, char** argv) {
    if (argc==2 && std::string(argv[1])=="--test") {self_test(); return 0;}
    if (argc!=4) throw std::runtime_error("Expected substep file, reference end rest, action scale");
    const Real final_rest = std::stold(argv[2]), scale = std::stold(argv[3]);
    if (final_rest<=0 || scale<=0) throw std::runtime_error("Positive normalization required");
    std::ifstream in(argv[1],std::ios::binary);
    if (!in) throw std::runtime_error("Could not open substeps");
    Angular totals;
    Real peak=-1, peak_time=0, peak_fraction=0, peak_c=0, peak_d=0, delivered=0, tangent_weight=0;
    Real max_pointwise=0, aggregate_unfilled=0, aggregate_peak=0;
    unsigned long count=0;
    std::array<double,12> record;
    while (in.read(reinterpret_cast<char*>(record.data()),sizeof(record))) {
        const auto [t,dt,k0,k1,field0,field1,drest,f0,f1,dk,k,target] = record;
        if (!(target>0 && dk>=0 && dt>0 && drest>0)) throw std::runtime_error("Invalid substep");
        const Real request = scale*target, c=dk/request, d0=field0/request, d1=field1/request;
        delivered+=drest;
        totals = totals + field_average(c,d0,d1)*drest;
        const auto aggregate = sweep_integral(field0,field1,dk,static_cast<double>(request));
        aggregate_unfilled += drest*aggregate.unmet/request;
        const double min_field = field0*field1 <= 0 ? 0 : std::min(std::abs(field0),std::abs(field1));
        aggregate_peak = std::max(aggregate_peak,request/directional_mean(min_field,dk));
        const Real d = peak_field(c,d0,d1);
        const Real here = at_field(c,std::abs(d)).unfilled;
        if (here>peak) {
            const Real s = d0==d1 ? 0 : (d-d0)/(d1-d0);
            peak=here; peak_time=t+s*dt; peak_fraction=f0+s*(f1-f0);
            peak_c=c; peak_d=d;
        }
        const Real lo=std::min(d0,d1), hi=std::max(d0,d1);
        max_pointwise=std::max(max_pointwise,1-std::max(0.L,c-std::max(std::abs(lo),std::abs(hi))));
        const Real no_tangent = hi==lo ? (std::abs(lo)<c ? 1 : 0) :
            std::max(0.L,std::min(hi,c)-std::max(lo,-c))/(hi-lo);
        tangent_weight+=drest*(1-no_tangent);
        ++count;
    }
    if (!in.eof() || in.gcount() || !count) throw std::runtime_error("Invalid substep file");
    std::cout << std::setprecision(17) << "{\n"
        << "\"steps\":" << count << ",\n"
        << "\"aggregate_peak_demand_over_capacity\":" << aggregate_peak << ",\n"
        << "\"aggregate_unfilled_fraction_of_total\":" << aggregate_unfilled/final_rest << ",\n"
        << "\"fraction_tested\":" << delivered/final_rest << ",\n"
        << "\"unfilled_fraction_of_total\":" << totals.unfilled/final_rest << ",\n"
        << "\"unfilled_fraction_of_tested\":" << totals.unfilled/delivered << ",\n"
        << "\"in_deficient_directions_fraction_of_total\":" << totals.deficient/final_rest << ",\n"
        << "\"tangent_present_fraction_of_total\":" << tangent_weight/final_rest << ",\n"
        << "\"peak_instantaneous_unfilled_fraction\":" << peak << ",\n"
        << "\"peak_tau\":" << peak_time << ",\n"
        << "\"peak_cumulative_production\":" << peak_fraction << ",\n"
        << "\"peak_radial_speed_over_request\":" << peak_c << ",\n"
        << "\"peak_signed_field_over_request\":" << peak_d << ",\n"
        << "\"largest_single_direction_unfilled_fraction\":" << max_pointwise << "}\n";
}
