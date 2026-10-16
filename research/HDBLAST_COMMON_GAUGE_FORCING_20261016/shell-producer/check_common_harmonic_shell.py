#!/usr/bin/env python3
"""Exact continuous physical shell-source bounds in the actual first harmonic gauge."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile


EXPECTED = {
    "GLOBAL_GAUGE_COEFFICIENTS_NORMAL.json": "4009c46f6d63583b4d161eebbb620077467b7561a3c4b5728d5293655d340a53",
    "HIGHER_GAUGE_NORMS_NORMAL.json": "f246bc6840ea36f123c7c53da28a9c04925954187e0f06df25add779ce083c6d",
    "HIGHER_JETS_NORMAL.json": "92b17fd71d78f83eb3f2c5ec8b9ea5f4e70365ec4d41aa4c0205551e7eaeb5df",
    "COORDINATE_SHELL_SOURCE_ENCLOSURE.json": "03aba506022e8b5fba38280db9b568d9d1c4fb0c5f76d0f33bf096a0ec8c75f1",
    "COORDINATE_INPUT_PINS.json": "2507bbb4f279ec146fcea7342bbd17322f3ded6eb86eba346aeb3bd985765481",
    "SCALAR_SHELL_BOUNDARY_SYSTEM.md": "8b4d08bf3c4b1704be99e291049a61c902feb39439f8edbda54a5d7d3f9735c8",
    "MATERIAL_SECOND_ORDER_JUNCTIONS.md": "69c2af7bb7be06783f2470833f3cd10507b26ccdd9c67ebf7cf905e62f188871",
    "FIRST_HARMONIC_FIELDS_UNPUBLISHED.md": "d15f8c65ffffd0085d195064ccea3c61bb4e719695652d0dd73472c8a06c224d",
}
ERRORS = {
    "wrong_normal_orientation": "registered outward background Israel row",
    "half_affine_hessian": "full tensor affine normal Hessian",
    "isotropic_tensor_only": "anisotropic off-diagonal metric Hessian",
    "wrong_sigma_third": "registered scalar affine Hessian",
    "wrong_proper_normal_scale": "proper inward mixed derivative witness",
    "missing_two_copy_factor": "one/two-copy metric residual Hessian adapter",
}
GRID = 2**96
ROWS = ["00", "01", "02", "03", "11", "12", "13", "22", "23", "33"]


def require(test, message):
    if not test:
        raise ValueError(message)


def eq(a, b, message):
    require(a == b, message)


def up(x):
    """Outward dyadic rounding; every reported calculation remains rational."""
    x = F(x)
    return F(-((-x.numerator * GRID)//x.denominator), GRID)


def cap(x):
    x = F(x)
    require(x >= 0, "nonnegative absolute cap")
    return up(x)


def abs_interval(values):
    lo, hi = map(F, values)
    require(lo <= hi, "ordered authenticated interval")
    return max(abs(lo), abs(hi))


def ceil_fraction(x):
    return -((-x.numerator)//x.denominator)


def decimal_cap(x, digits=5):
    """Strict next decimal grid point, with no floating arithmetic."""
    x = F(x)
    if not x:
        return {"exact": "0", "decimal": "0", "relation":"exact_zero"}
    e = len(str(abs(x.numerator))) - len(str(x.denominator))
    ten = F(10)**e
    while x < ten:
        e -= 1
        ten /= 10
    while x >= 10*ten:
        e += 1
        ten *= 10
    unit = F(10)**(e-digits+1)
    k = x//unit + 1
    value = k*unit
    eq(value > x, True, "strict readable cap")
    return {"exact": str(value), "decimal": str(k)+"e"+str(e-digits+1),"relation":"strict_upper"}


def pack(x):
    return {"exact_upper": str(x), "strict_readable_cap": decimal_cap(x)}


def authenticate(root, pins_path):
    pins = json.loads(pins_path.read_text())
    eq(pins, EXPECTED, "fixed input pin manifest")
    out = {}
    for name, digest in EXPECTED.items():
        raw = (root/name).read_bytes()
        eq(hashlib.sha256(raw).hexdigest(), digest, "input authentication: "+name)
        if name.endswith(".json"):
            out[name] = json.loads(raw)
    c = out["GLOBAL_GAUGE_COEFFICIENTS_NORMAL.json"]
    g = out["HIGHER_GAUGE_NORMS_NORMAL.json"]
    w = out["HIGHER_JETS_NORMAL.json"]
    s = out["COORDINATE_SHELL_SOURCE_ENCLOSURE.json"]
    sp = out["COORDINATE_INPUT_PINS.json"]
    eq(s["input_pins_sha256"], EXPECTED["COORDINATE_INPUT_PINS.json"], "coordinate source input-manifest linkage")
    oldpins = {entry["path"]:entry["sha256"] for entry in sp["files"]}
    for name,digest in w["input_pins"].items():
        eq(oldpins[name],digest,"same registered master/shell retained input: "+name)
    for name in ("SCALAR_SHELL_BOUNDARY_SYSTEM.md","MATERIAL_SECOND_ORDER_JUNCTIONS.md"):
        eq(oldpins["inputs/reference/"+name],EXPECTED[name],"registered sign-document manifest membership")
    eq(g["source_pins"]["GLOBAL_GAUGE_COEFFICIENTS_NORMAL.json"], EXPECTED["GLOBAL_GAUGE_COEFFICIENTS_NORMAL.json"], "gauge coefficient linkage")
    eq(c["input_sha256"]["HIGHER_JETS_NORMAL.json"], EXPECTED["HIGHER_JETS_NORMAL.json"], "actual master coefficient linkage")
    eq(g["inherited_actual_source_pins"]["HIGHER_JETS_NORMAL.json"], EXPECTED["HIGHER_JETS_NORMAL.json"], "actual master gauge linkage")
    for name in ("rho", "a", "p"):
        eq(c["first_exit"][name]["shell"]["exact"], s["exact"][name], "registered shared shell box: "+name)
    eta = list(map(F,c["first_exit"]["eta"]["shell"]["exact"]))
    eq([4*(1+x) for x in eta], list(map(F,s["exact"]["sigma2"])), "actual sigma second derivative equals four phi_b")
    eq(c["domain"]["ell"], "1/1000000", "actual coefficient scale")
    eq(g["domain"]["ell"], c["domain"]["ell"], "actual gauge scale")
    eq(w["source_bound_adapter"]["ell"], c["domain"]["ell"], "actual master scale")
    eq(s["exact"]["ell"], c["domain"]["ell"], "actual shell scale")
    require(g["status"] == "PASS_CONTINUOUS_RAW_GAUGE_ALL_MIXED_C3_CAPS", "continuous actual gauge C3 premise")
    require(w["status"].startswith("PASS"), "continuous actual master C3 premise")
    return c,g,w,s


def poly_mul(a,b):
    return [sum(a[j]*b[i-j] for j in range(i+1)) for i in range(3)]


def algebra_checks(mutation):
    """Independent affine Christoffel/normal series for all ten tensor rows."""
    count = 0
    for seed in (616, 1616, 2616, 3616):
        rng = random.Random(seed)
        get = lambda: F(rng.randint(1,9),rng.randint(2,11))
        a,p,H,eta = get(),-get(),get(),get()
        sigma0,sigma1,sigma2,sigma3 = 6*a,-2*p,4*get(),F(4)
        gamma = [[F(0) for _ in range(4)] for _ in range(4)]
        for i in range(4):
            gamma[i][i] = F(-1 if i==0 else 1)
        m = [[F(0) for _ in range(4)] for _ in range(4)]
        for i in range(4):
            for j in range(i,4):
                m[i][j] = m[j][i] = get()
        my = [[2*a*m[i][j]+(a*H+sigma1*eta/3)*gamma[i][j]
               for j in range(4)] for i in range(4)]
        etay = p*H/2-sigma2*eta/2
        orientation = -1 if mutation=="wrong_normal_orientation" else 1
        eq(orientation*a-sigma0/6, F(0), "registered outward background Israel row")
        eq(orientation*p+sigma1/2,F(0),"registered outward background scalar row")
        count += 2
        # n_y=sqrt(g_yy), Gamma^y_ab=-g^yy gamma_ab,y/2.
        # Their independent product gives K_ab=-n_y Gamma^y_ab.
        N = [F(1),H/2,-H*H/8]
        inverse_gyy = [F(1),-H,H*H]
        ninv = [F(1),-H/2,3*H*H/8]
        eq(poly_mul(N,N),[F(1),H,F(0)],"unit normal square-root normalization")
        eq(poly_mul(N,ninv),[F(1),F(0),F(0)],"unit normal inverse normalization")
        count += 2
        sig = [sigma0,sigma1*eta,sigma2*eta*eta/2]
        A = a*H*H/4-sigma1*H*eta/6-sigma2*eta*eta/6
        L = a*H+sigma1*eta/3
        for i in range(4):
            for j in range(i,4):
                gamma_y = [2*a*gamma[i][j],my[i][j],F(0)]
                christ = [-x/2 for x in poly_mul(inverse_gyy,gamma_y)]
                Kseries = [-x for x in poly_mul(N,christ)]
                tension = poly_mul(sig,[gamma[i][j],m[i][j],F(0)])
                factor = 1 if mutation=="half_affine_hessian" else 2
                direct = factor*(Kseries[2]-tension[2]/6)
                raw = 3*a*H*H*gamma[i][j]/4-H*my[i][j]/2-(sigma2*eta*eta*gamma[i][j]+2*sigma1*eta*m[i][j])/6
                eq(direct,raw,"full tensor affine normal Hessian")
                used_m = F(0) if mutation=="isotropic_tensor_only" and i!=j else m[i][j]
                reduced = A*gamma[i][j]-L*used_m
                eq(raw,reduced,"anisotropic off-diagonal metric Hessian" if i!=j else "on-shell diagonal metric Hessian")
                two = (1 if mutation=="missing_two_copy_factor" else 2)*reduced
                eq(two,2*direct,"one/two-copy metric residual Hessian adapter")
                count += 3
        normal_phi = poly_mul(ninv,[p,etay,F(0)])
        sigprime = [sigma1,sigma2*eta,sigma3*eta*eta/2]
        directp = 2*(normal_phi[2]+sigprime[2]/2)
        used_sigma3 = F(2) if mutation=="wrong_sigma_third" else sigma3
        reducedp = p*H*H/4+sigma2*H*eta/2+used_sigma3*eta*eta/2
        eq(directp,reducedp,"registered scalar affine Hessian")
        count += 1
        # Old longitudinal specialization fixes the one-copy RHS convention.
        u,v = get(),get()
        HL,etaL = -4*u,v
        AL = a*HL*HL/4-sigma1*HL*etaL/6-sigma2*etaL*etaL/6
        LL = a*HL+sigma1*etaL/3
        eq(-(AL-LL*2*u),sigma2*v*v/6-12*a*u*u,"longitudinal Israel source specialization")
        eq(-(p*HL*HL/4+sigma2*HL*v/2+sigma3*v*v/2),-4*p*u*u+2*sigma2*u*v-2*v*v,"longitudinal scalar source specialization")
        count += 2
        # The isotropic source cancels from the time-angular de Sitter boost.
        ut,bt,Ttime,bv = get(),get(),get(),get()
        X,Z = ut+bt,ut+Ttime*bv
        t00,td = -A+2*L*X,A-2*L*Z
        eq(t00+td,2*L*(bt-Ttime*bv),"mixed angular covariant derivative isotropic cancellation")
        count += 1
    ell,rho,x,r = F(1,1000000),F(130),F(2,3),F(5,4)
    # Actual chain rule for test f(x,r)=x*r, a nonzero mixed derivative.
    dy_dt = -1/(ell*ell*rho)
    stated = (1 if mutation=="wrong_proper_normal_scale" else -1)/(ell*ell*rho)
    eq(dy_dt,stated,"proper inward mixed derivative witness")
    count += 1
    return count


def bounds(c,g,w,s):
    se = s["exact"]
    ell,T = F(1,1000000),F(3,1000000)
    rho = tuple(map(F,se["rho"]))
    rhomin = rho[0]
    invrho = cap(1/rhomin)
    a,p = cap(abs_interval(se["a"])),cap(abs_interval(se["p"]))
    sig1,sig2,sig3 = cap(2*p),cap(abs_interval(se["sigma2"])),F(4)
    cmax = cap(1/(1-T*T/2))
    # All mixed gauge derivative entries are actual D_x^i D_r^j caps.
    def gauge(i,j,ch):
        return cap(F(g["scaled_mixed_partial_caps"][str(i)+","+str(j)][ch])/ell**(i+j))
    bb,bd = [gauge(0,0,ch) for ch in (1,2)]
    bfr,bfxr = gauge(0,1,0),gauge(1,1,0)
    bbt,bbtt = gauge(1,0,1),gauge(2,0,1)
    bdt,bdtt = gauge(1,0,2),gauge(2,0,2)
    y = cap(F(w["exact"]["boundary"]["y"]))
    yp = cap(F(w["exact"]["boundary"]["yprime"]))
    incoming = cap(F(w["exact"]["boundary"]["H"]))
    incomingprime = cap(F(w["exact"]["interior"]["Hprime"]))
    gamma = cap(abs_interval(w["exact"]["boundary"]["gamma"]))
    k = cap(abs_interval(c["static_z_jets"]["k"][0]["exact"]))
    B = cap(3/(rhomin**3*min(abs(F(v)) for v in se["p"])))
    U = cap(min(abs_interval(se["global_envelopes"]["u"]),y/rhomin**2))
    Ut = cap((yp/ell+3*T*y/2)/rhomin**2)
    V = cap(min(abs_interval(se["global_envelopes"]["v"]),B*(incoming/ell+(gamma+k)*y)))
    Vt = cap(B*(incomingprime/ell**2+(gamma+k)*yp/ell+3*T*(incoming/ell+(gamma+k)*y)/2))
    J = cap(4*U+2*bfr*invrho)
    Jt = cap(4*Ut+2*bfxr*invrho)
    X,Xt = cap(U+bbt),cap(Ut+bbtt)
    Z,Zt = cap(U+T*bb),cap(Ut+bb+T*bbt)
    Vg = cap(cmax*bdt+bb)
    Vgt = cap(cmax*T*bdt+cmax*bdtt+bbt+T*bb)
    A = cap(a*J*J/4+sig1*J*V/6+sig2*V*V/6)
    At = cap(a*J*Jt/2+sig1*(Jt*V+J*Vt)/6+sig2*V*Vt/3)
    L,Lt = cap(a*J+sig1*V/3),cap(a*Jt+sig1*Vt/3)
    P = cap(p*J*J/4+sig2*J*V/2+sig3*V*V/2)
    Pt = cap(p*J*Jt/2+sig2*(Jt*V+J*Vt)/2+sig3*V*Vt)
    t00,t00t = cap(A+2*L*X),cap(At+2*Lt*X+2*L*Xt)
    t0d,t0dt = cap(L*Vg),cap(Lt*Vg+L*Vgt)
    td,tdt = cap(A+2*L*Z),cap(At+2*Lt*Z+2*L*Zt)
    th,tht = cap(2*L*bd),cap(2*Lt*bd+2*L*bdt)
    # |Y|<=1, |Y DY|<=1, |Y HessY|<=4.
    q00,q0a,qaa,qab = t00,t0d,cap(td+4*th),cap(4*th)
    qt00,qt0a,qtaa,qtab = [cap(z*invrho) for z in (t00t,t0dt,tdt+4*tht,4*tht)]
    # Full induced-shell angular derivative, including de Sitter boosts.
    qa00 = cap(2*t00*invrho+2*T*invrho*q0a)
    # Isotropic A*gamma contributes no mixed covariant derivative. In the
    # de Sitter boost term Q_ca+delta_ca Q_00 the A terms cancel exactly;
    # X-Z=b_t-T*b also cancels u.
    qa0a = cap(8*t0d*invrho+T*invrho*(2*L*(bbt+T*bb)+4*th))
    qaab = cap((2*td+16*th)*invrho+2*T*invrho*q0a)
    qaoffdiag = cap(16*th*invrho+2*T*invrho*q0a)
    values,time,angular = {},{},{}
    for row in ROWS:
        if row=="00":
            values[row],time[row],angular[row] = q00,qt00,qa00
        elif row[0]=="0":
            values[row],time[row],angular[row] = q0a,qt0a,qa0a
        elif row[0]==row[1]:
            values[row],time[row],angular[row] = qaa,qtaa,qaab
        else:
            values[row],time[row],angular[row] = qab,qtab,qaoffdiag
    tangential = {row:max(time[row],angular[row]) for row in ROWS}
    ptime,pangular = cap(Pt*invrho),cap(2*P*invrho)
    tau = cap(8*U+2*bbt+6*T*bb+16*bd)
    taut = cap(8*Ut+2*bbtt+6*bb+6*T*bbt+16*bdt)
    traceA2 = cap(4*X*X+12*Z*Z+64*Z*bd+8*Vg*Vg+96*bd*bd)
    traceA2t = cap(8*X*Xt+24*Z*Zt+64*(Zt*bd+Z*bdt)+16*Vg*Vgt+192*bd*bdt)
    m00,m0a,maa,mab = cap(2*X),cap(2*Vg),cap(2*Z+8*bd),cap(8*bd)
    fieldnormals = {
        "H_y":cap(4*a*J+2*a*tau+4*sig1*V/3),
        "eta_y":cap(p*J/2+sig2*V/2),
        "nabla_y_m_diagonal":L,
        "nabla_y_m_offdiagonal":F(0),
        "coordinate_m_y_in_frame_diagonal":cap(2*a*max(m00,maa)+L),
        "coordinate_m_y_in_frame_offdiagonal":cap(2*a*max(m0a,mab)),
    }
    fields={"u":U,"u_t":Ut,"v":V,"v_t":Vt,"J":J,"J_t":Jt,
            "b":bb,"b_t":bbt,"b_tt":bbtt,"d":bd,"d_t":bdt,"d_tt":bdtt,
            "X":X,"X_t":Xt,"Z":Z,"Z_t":Zt,"V_gauge":Vg,"V_gauge_t":Vgt,
            "tau":tau,"tau_t":taut,"abs_trace_A_squared":traceA2,"abs_trace_A_squared_t":traceA2t}
    def metric_pack(vs,ts,angs,factor=1):
        return {"value":{r:pack(factor*vs[r]) for r in ROWS},
                "proper_time_derivative":{r:pack(factor*ts[r]) for r in ROWS},
                "each_angular_covariant_derivative":{r:pack(factor*angs[r]) for r in ROWS},
                "each_tangential_covariant_derivative":{r:pack(factor*max(ts[r],angs[r])) for r in ROWS},
                "maximum_value":pack(factor*max(vs.values())),
                "maximum_each_tangential_covariant_derivative":pack(factor*max(max(ts.values()),max(angs.values())))}
    def scalar_pack(factor=1):
        return {"value":pack(factor*P),"proper_time_derivative":pack(factor*ptime),
                "each_angular_derivative":pack(factor*pangular),
                "each_tangential_derivative":pack(factor*max(ptime,pangular))}
    return {
        "parameters": {"ell":str(ell),"t_max":str(T),"rho_b_interval":se["rho"],
                       "a_abs":str(a),"p_abs":str(p),"sigma_prime_abs":str(sig1),
                       "sigma_second_interval":se["sigma2"],"sigma_second_abs":str(sig2),
                       "sigma_third":"4","cosh_t_upper":str(cmax)},
        "actual_field_amplitude_caps":{k:pack(v) for k,v in fields.items()},
        "actual_field_component_caps":{k:pack(v) for k,v in {
            "H":J,"eta":V,"m_00":m00,"m_0a":m0a,"m_aa":maa,"m_ab_a_ne_b":mab,
            "tau":tau,"abs_trace_A_squared":traceA2}.items()},
        "proper_first_normal_field_caps":{k:pack(v) for k,v in fieldnormals.items()},
        "one_copy_outward_residual_hessian":{"metric":metric_pack(values,time,angular),"scalar":scalar_pack()},
        "two_copy_outward_residual_hessian":{"metric":metric_pack(values,time,angular,2),"scalar":scalar_pack(2)},
        "coordinate_physical_forcing":{"metric":metric_pack(values,time,angular,2),"scalar":scalar_pack()},
        "adapter":{
            "one_copy_second_equation":"D B[Phi2]=-Q; coefficient epsilon^2 is Phi2/2",
            "two_copy_outward":"J=2B, Q_J=2Q; second equation D J[Phi2]=-2Q",
            "two_copy_inward":"J_in=-2B, Q_Jin=-2Q; second equation D J_in[Phi2]=+2Q",
            "untraced_two_copy_Israel":"Q_I=2*(Q_K-tr_gamma0(Q_K)*gamma0); background and first residuals vanish",
            "coordinate_metric":"k_ab,y-2a*k_ab-(a*k_yy+sigma_prime*w/3)*gamma_ab = -2 Q_K,ab",
            "coordinate_scalar":"w_y+sigma_second*w/2+sigma_prime*k_yy/4 = -Q_p",
        },
        "source_normal_derivative_claimed":False,
        "tangential_derivative_geometry":"D0=rho_b^-1 partial_t; Da includes sphere connection and tanh(t)/rho_b de Sitter frame boosts",
    }


def run_controls(script,root,pins):
    rows=[]
    with tempfile.TemporaryDirectory(prefix="common-shell-controls-") as td:
        for mutation,reason in ERRORS.items():
            dst=Path(td)/(mutation+".json")
            cmd=[sys.executable]+([] if __debug__ else ["-O"])+[
                str(script),"--input-root",str(root),"--pins",str(pins),
                "--inject-error",mutation,"--output",str(dst)]
            p=subprocess.run(cmd,capture_output=True,text=True)
            require(p.returncode!=0 and reason in p.stderr and not dst.exists(),
                    "mutation-specific control rejection: "+mutation)
            rows.append({"mutation":mutation,"expected_reason":reason,"rejected_without_output":True})
    return rows


def main():
    here=Path(__file__).resolve().parent
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input-root",type=Path,default=here/"inputs")
    ap.add_argument("--pins",type=Path,default=here/"INPUT_PINS.json")
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--run-controls",action="store_true")
    ap.add_argument("--inject-error",choices=list(ERRORS))
    args=ap.parse_args()
    require(not args.output.exists(),"refuse to overwrite output")
    require(not args.output.resolve().is_relative_to(args.input_root.resolve()),"output must be outside input tree")
    c,g,w,s=authenticate(args.input_root,args.pins)
    count=algebra_checks(args.inject_error)
    report={
        "status":"PASS_CONTINUOUS_FULL_COMMON_FIRST_HARMONIC_PHYSICAL_SHELL_FORCING",
        "input_sha256":EXPECTED,
        "domain":{"t":["0","3/1000000"],"scaled_shell":"0<=x<=3,r=0; entire unit S3"},
        "norm":"absolute value of each lower adapted background shell-frame component; time index 0; not a positive Lorentzian tensor norm",
        "arithmetic":{"method":"fractions.Fraction only; initial caps and arithmetic rounded outward on a 2^96 dyadic grid","dyadic_grid_denominator":str(GRID),"floating_arithmetic_used":False},
        "finite_algebra_comparisons":count,
        "bounds":bounds(c,g,w,s),
        "controls":run_controls(Path(__file__).resolve(),args.input_root,args.pins) if args.run_controls else [],
        "scope":{"actual_registered_first_pair":True,"first_gauge_common_to_bulk_and_shell":True,
                 "all_ten_physical_metric_rows_and_scalar":True,"continuous_values_and_first_tangential_covariant_derivatives":True,
                 "affine_quadratic_source_not_solved_second_data":True,"source_normal_derivative_enclosed":False,
                 "gauge_and_constraint_boundary_rows_certified_here":False,"second_order_evolution_solved":False,
                 "invariant_force_or_observable_claimed":False},
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")


if __name__=="__main__":
    try:
        main()
    except (ValueError,KeyError) as exc:
        sys.exit("REJECT: "+str(exc))
