"""Certified fixed-frequency-ratio distance coefficients for two source clocks."""
import csv
import hashlib
import io
import json
from fractions import Fraction as F
from pathlib import Path

PARENT = Path(__file__).resolve().parent.parent / "HDBLAST_JOINT_CURVATURE_20261118"
PINS = {
    "inputs/FUTURE_CELLS.csv": "2f76fd1ae2195b445582a29a352bd76658dea17c2233cbb9acc395b255172bf7",
    "inputs/SHELL_BOX.json": "8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8",
    "METRIC.json": "d7ce7f5c602bd5f7df159ccdd7e2e099c03c4b8740bfb4397c178ab3f2addf42",
}
ELL = F(1,10**6)
X0, X1 = F(49,32), F(99,64)
SPAN = X1-X0


def inputs(parent_dir):
    files = {}
    for name, pin in PINS.items():
        content = (parent_dir/name).read_bytes()
        if hashlib.sha256(content).hexdigest() != pin:
            raise ValueError(f"input digest mismatch: {name}")
        files[name] = content
    return files


def multiply(a, b):
    products = [x*y for x in a for y in b]
    return min(products), max(products)


def shifted_corners(a, error):
    return a[0]-error, a[1]+error


def combined_integrals(csv_data):
    # For unit acceleration at s, z0(X1)=t and integral z0=t^2/2, t=X1-s.
    # The combined kernel for 2*integral z0 - SPAN*z0(X1) is t*(t-SPAN).
    # Its change of sign at s=X0 falls on a certified cell boundary.
    edge, count = F(1), 0
    intervals = [[F(0),F(0)] for _ in range(3)]
    for row in csv.DictReader(io.StringIO(csv_data.decode("ascii"))):
        a = max(F(1),F(int(row["left_grid"]),4096))
        b = min(X1,F(int(row["right_grid"]),4096))
        if b<=a:
            continue
        if a!=edge or a<X0<b:
            raise ValueError("coverage or kernel sign")
        acc = (F(int(row["H_lo"]),2**112)-F(3,10**6),
               F(int(row["H_hi"]),2**112)+F(3,10**6))
        if acc[1]<acc[0]:
            raise ValueError("acceleration interval")
        ta,tb = X1-a,X1-b
        end_weight = (ta**2-tb**2)/2
        integral_weight = (ta**3-tb**3)/6
        joint_weight = 2*integral_weight-SPAN*end_weight
        for total, weight in zip(intervals,(end_weight,integral_weight,joint_weight)):
            low,high = sorted((weight*acc[0],weight*acc[1]))
            total[0]+=low
            total[1]+=high
        edge=b
        count+=1
    if edge!=X1:
        raise ValueError("source time is outside the input")
    return count, intervals


def compute(parent_dir):
    data=inputs(parent_dir)
    metric=json.loads(data["METRIC.json"])
    box=json.loads(data["inputs/SHELL_BOX.json"])["boxes"]
    rho_low,rho_high=map(F,box["rho_b"]["exact"])
    b_low,b_high=map(F,box["B_b"]["exact"])
    physical_scale=(ELL*b_low/rho_high,ELL*b_high/rho_low)
    count,(endpoint,integral,combined)=combined_integrals(data["inputs/FUTURE_CELLS.csv"])
    e=F(metric["absolute_mean_error"])  # Uniform supremum error from the Volterra estimate.
    bound=F(metric["z_bound"])  # Bounds both |z| and |M|=sech(q)^(3/2)|z|.
    M=shifted_corners(endpoint,e)
    integral_M=shifted_corners(integral,(X1-1)*e)
    c=1+F(1,10**11)
    q0,q1=ELL*X0,ELL*X1
    if 1/(1-q1**2)>=c:
        raise ValueError("cosh enclosure")
    # Normalization: U=(ell B_b/rho_b)*M. For supported sources the
    # bracket is -(sin(theta)/ell)*M1; geodesic sources add 2*int M /cosh q0.
    # d <= (sinh(q1)-sinh(q0))/ell <= C*d; 1 <= cosh(qi) <= C.
    sin_over_ell=(SPAN/c**2, SPAN*c)
    supported=multiply((-sin_over_ell[1],-sin_over_ell[0]),M)
    jacobian_error=(2*(X1-1)*(1-1/c)+SPAN*(c**2-1))*bound
    joint_error=(2*(X1-1)+SPAN)*e+jacobian_error
    geodesic=shifted_corners(combined,joint_error)
    # cos(theta) >= 1 - sin(theta)^2; sinh(q1)/q1 in [1,C].
    common=((1-(ELL*sin_over_ell[1])**2)/(c*X1),1/X1)
    if common[0]<=0 or supported[1]>=0 or geodesic[0]<=0:
        raise ValueError("clock or distance signs")
    # 2*sinh((q1+q0)/2)*sinh((q1-q0)/2)=cosh(q1)-cosh(q0).
    redshift_lead=(q1**2-q0**2)/2
    target_z=(-redshift_lead*c**2,-redshift_lead/c)
    doppler_ratio=(2*integral_M[0]/(c*sin_over_ell[1]*M[1]),
                   2*integral_M[1]/(sin_over_ell[0]*M[0]))
    if not doppler_ratio[0] > 17:
        raise ValueError("velocity term must exceed the supported term")
    modes={}
    for name,bracket in (("geodesic",geodesic),("supported",supported)):
        normalized=multiply(common,bracket)
        physical=multiply(multiply(normalized,physical_scale),(F(18,25),F(18,25)))
        modes[name]={"axis_per_ell_B_over_rho":list(map(str,normalized)),
                     "fractional_distance_cap_half_difference":list(map(str,physical)),
                     "display_only":list(map(float,physical))}
    return {"source_sha256":PINS, "receiver_x":str(X0), "source_background_x":str(X1),
            "source_redshift_background":list(map(str,target_z)), "acceleration_cells":count,
            "M_at_source":list(map(str,M)), "integral_M_dx_to_source":list(map(str,integral_M)),
            "ell_B_over_rho":list(map(str,physical_scale)),"geodesic_bracket_error":str(joint_error),
            "half_difference_convention":"(positive cap minus negative cap)/(2*background distance), differentiated at zero",
            "doppler_to_abs_supported_term":list(map(str,doppler_ratio)), "comparisons":modes}


if __name__=="__main__":
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("joint_capsule",nargs="?",type=Path,default=PARENT)
    print(json.dumps(compute(parser.parse_args().joint_capsule),indent=2))
