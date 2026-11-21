// Independent angular audit. Average the clipped deficit over proper time
// at fixed angle before doing the angular quadrature. Speeds are divided by
// the requested source per shell area; the remaining deficit is a fraction.
namespace reverse_order {
using Real = long double;

Real time_mean(Real x, Real c, Real d0, Real d1) {
    const Real a = c-d0*x, b = c-d1*x;
    const Real lo = std::min(a,b), hi = std::max(a,b);
    if (lo >= 1 || hi <= -1) return 0;
    if (lo == hi) return 1-std::abs(lo);
    Real area = 0;
    const Real left = std::max(lo,Real(-1)), right = std::min(hi,Real(1));
    const Real negative_end = std::min(Real(0),right);
    if (negative_end > left)
        area += (negative_end-left)*(1+(left+negative_end)/2);
    const Real positive_start = std::max(Real(0),left);
    if (right > positive_start)
        area += (right-positive_start)*(1-(positive_start+right)/2);
    return area/(hi-lo);
}

Real panel(Real lo, Real hi, Real c, Real d0, Real d1) {
    // Positive nodes/weights of order-eight Gauss-Legendre quadrature.
    constexpr Real x[] = {0.18343464249564980494L,0.52553240991632898582L,
                          0.79666647741362673959L,0.96028985649753623168L};
    constexpr Real w[] = {0.36268378337836198297L,0.31370664587788728734L,
                          0.22238103445337447054L,0.10122853629037625915L};
    const Real mid = (lo+hi)/2, half = (hi-lo)/2;
    Real sum=0;
    for (int i=0;i<4;++i) {
        const Real shift=half*x[i];
        sum+=w[i]*(time_mean(mid-shift,c,d0,d1)+time_mean(mid+shift,c,d0,d1));
    }
    return half*sum;
}

Real angular_mean(Real c, Real d0, Real d1) {
    Real total=0;
    for (int sign : {-1,1}) {
        // Work on the positive magnitude of the direction cosine. Endpoint
        // preimages of the tent kinks (-1,0,1) split every change of formula.
        const Real e0=sign*d0,e1=sign*d1;
        std::array<Real,8> cuts{0,1};
        unsigned count=2;
        for (const Real end : {e0,e1})
            if (end!=0)
                for (int kink : {-1,0,1}) {
                    const Real x=(c-kink)/end;
                    if (x>0 && x<1) cuts[count++]=x;
                }
        std::sort(cuts.begin(),cuts.begin()+count);
        for (unsigned i=1;i<count;++i) {
            Real lo=cuts[i-1];
            const Real hi=cuts[i];
            if (lo==hi) continue;
            // In the first panel the integrand is linear. Elsewhere it has
            // form A/x+B+Cx; geometric panels resolve the reciprocal term.
            while (lo>0 && 2*lo<hi) {
                total+=panel(lo,2*lo,c,e0,e1);
                lo*=2;
            }
            total+=panel(lo,hi,c,e0,e1);
        }
    }
    return total/2;
}
}

void check_sweep_integral() {
    using reverse_order::angular_mean;
    // d0=d1=2, c=0: angular integral of the central tent is 1/4.
    // d varies uniformly on [-2,2], c=0: its product with a uniform angular
    // cosine gives 3/8+log(2)/4 for the expected clipped deficit.
    if (std::abs(angular_mean(0,2,2)-.25L)>1e-12L ||
        std::abs(angular_mean(0,-2,2)-(.375L+std::log(2.L)/4))>1e-12L ||
        std::abs(angular_mean(.3L,0,0)-.7L)>1e-12L ||
        angular_mean(10,-3,3)!=0)
        throw std::runtime_error("Reverse-order angular controls failed");
}

struct SweepCollector {
    long double unmet=0, delivered=0;
    void collect(double, double, double d0, double d1, double c, double k,
                 double source, double drest, double, double) {
        if (!(source>0 && k>0 && drest>=0))
            throw std::runtime_error("Invalid source or birth shell");
        const long double scale=12.L*k*k/source;
        const auto deficit=reverse_order::angular_mean(c*scale,d0*scale,d1*scale);
        if (!(deficit>=-1.e-12L && deficit<=1+1.e-12L))
            throw std::runtime_error("Angular deficit outside its physical bounds");
        unmet+=drest*deficit;
        delivered+=drest;
    }
    void save(std::map<std::string,double>& out,const std::string& prefix,double final_rest) {
        out[prefix+"unmet_fraction"]=double(unmet/final_rest);
        out[prefix+"fraction_tested"]=double(delivered/final_rest);
    }
};
