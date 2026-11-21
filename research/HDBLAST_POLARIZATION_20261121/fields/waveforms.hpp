// Tube proper-time profiles. Compact profiles are normalized first in d tau.
// The reported modal responses are subsequently normalized to equal emitted
// energy, 2 integral r_tube^3 i^2 d tau, on each numerical background.
double shaped(double tau, int profile, double early) const {
    if (profile == 0) return early;
    if (profile == 1) return incoming(tau-4.989, dose, norm, stop);
    constexpr double centers[] = {6.2100, 6.2446, 6.2546, 6.2446, 6.2486};
    constexpr double durations[] = {.008, .006, .006, .026, .022};
    const double duration=durations[profile-2];
    const double u=(tau-centers[profile-2])/duration+.5;
    if (u<=0. || u>=1.) return 0.;
    const double pi=3.14159265358979323846;
    const double sine=std::sin(pi*u);
    if (profile == 6) return std::sqrt(32./(5.*duration))*sine*sine*std::sin(2.*pi*u);
    return std::sqrt(8./(3.*duration))*sine*sine;
}
