"""Explanatory geometry and exact flat-test-scalar flux, not a model simulation."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

out = Path(__file__).resolve().parent
out.mkdir(exist_ok=True)
plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False})
fig, (ax,bx) = plt.subplots(1,2,figsize=(11.2,4.3),layout="constrained")

T = np.linspace(-1.2,1.2,600)
shellR = np.sqrt(1+T*T)
ax.fill_betweenx(T,np.abs(T),shellR,color="#edf3f8",label="Exterior radial patch")
ax.plot(np.abs(T),T,color="#2878a5",lw=2)
ax.plot(shellR,T,color="#53616c",lw=1.5,label=r"Illustrative shell: $r=1$")
for rc in (0.12,0.30,0.55):
    ax.plot(np.sqrt(rc*rc+T*T),T,color="#9babb8",lw=.8,alpha=.8)
rayR=np.linspace(0,1.0,60)
ax.plot(rayR,np.tanh(.6)*rayR,"--",color="#7047a3",lw=1.8)
ax.annotate("Fixed finite $t$\nreaches the vertex",xy=(.05,.027),xytext=(.47,-.43),
            color="#7047a3",arrowprops={"arrowstyle":"->","color":"#7047a3"})
u=np.linspace(.03,.60,100); v=1.6
ax.plot((v+u)/2,(v-u)/2,color="#d7652c",lw=2)
ax.scatter([.8],[.8],color="#d7652c",s=30,zorder=5)
ax.annotate("Fixed $v=R+T>0$\nreaches a null-cone point",xy=(.8,.8),xytext=(1.0,.16),
            color="#b95022",arrowprops={"arrowstyle":"->","color":"#d7652c"})
ax.scatter([0],[0],color="#222222",s=18,zorder=6)
ax.text(.15,1.06,"Future cone: $u=R-T=0$",color="#2878a5")
ax.set(xlim=(-.03,1.9),ylim=(-1.2,1.2),xlabel="Ambient spatial radius $R$",
       ylabel="Ambient time $T$",title="A. Two different endpoint limits")
ax.set_aspect("equal",adjustable="box")
ax.legend(loc="lower right",fontsize=8,frameon=False)

eps=np.geomspace(1e-4,.3,180)
fixed=eps**4*(np.cosh(1.)**4-1.)/4
null=((2+eps**2/2)**4-(1+eps**2)**4)/64
bx.loglog(eps,fixed,color="#7047a3",lw=2,label=r"Fixed $0\leq t\leq1$")
bx.loglog(eps,null,color="#d7652c",lw=2,label=r"Fixed $1\leq v\leq2$")
bx.axhline(15/64,color="#53616c",lw=1,ls="--",label=r"Null limit $15/64$")
bx.set(xlabel=r"Regulator radius $y=\epsilon$",ylabel="Integrated symplectic flux per solid angle",
       title="B. Exact flat-cone test-scalar example")
bx.grid(which="major",alpha=.18)
bx.legend(loc="lower right",frameon=False,fontsize=9)
fig.savefig(out/"cone-limits.pdf")
fig.savefig(out/"cone-limits.png",dpi=170)
plt.close(fig)
with (out/"cone-limits.csv").open("w",newline="") as f:
    w=csv.writer(f);w.writerow(["epsilon","fixed_t_0_to_1_flux","fixed_v_1_to_2_flux","null_limit"])
    w.writerows(zip(eps,fixed,null,np.full_like(eps,15/64)))
(out/"README.md").write_text("""# Cone limits: explanatory figure

Panel A uses the exact local conformal geometry `g=Omega(r²)² eta`,
`r²=R²-T²`, and an illustrative shell normalized to `r=1`. It is a causal
coordinate diagram, not a numerical solution of the registered background.
Gray curves have fixed conformal radius. A bounded de Sitter time approaches
the vertex; fixed `v=R+T>0` approaches a finite null-cone cut.

Panel B plots exact fluxes for the flat massless test-scalar pair
`phi1=1`, `phi2=T`, per unit solid angle, oriented toward increasing y.
They can be realized locally by compact smooth finite-energy initial data.
This is a comparator, not a coupled HDBLAST solution or positive energy flux.

`I_fixed = epsilon^4 (cosh(1)^4-1)/4 -> 0`.
`I_null = [(2+epsilon²/2)^4-(1+epsilon²)^4]/64 -> 15/64`.

The CSV contains floating-point evaluations of these exact formulas for
plotting. The independently derived symbolic limits are in the domain audit;
plot points are not additional scientific tests.
""")
print(str(out/"cone-limits.pdf"))
