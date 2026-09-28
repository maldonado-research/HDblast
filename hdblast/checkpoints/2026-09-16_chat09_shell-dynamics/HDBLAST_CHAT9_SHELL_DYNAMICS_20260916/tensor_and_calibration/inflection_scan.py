import numpy as np, math
import shell_inflation as si
Y,PHI,F,FY,FYY,PY,SK,TH,TG = si.Y,si.PHI,si.F,si.FY,si.FYY,si.PY,si.SK,si.TH,si.TG
# inflection conditions at y_i: (num/f^2)_y = 0, (num/f^2)_yy = 0 ; num = 1 + c p + d p^2/2
for phi_i in (-0.9,-0.8,-0.6,-0.4,-0.2,-0.1,0.0,0.1,0.2,0.4,0.6,0.8):
    i = int(np.argmin(np.abs(PHI-phi_i))); p=PHI[i]; py=PY[i]; pyy = 2*p*(1-p*p)   # phi_yy = d/dy(-(1-p^2)) = 2 p p_y = -2p(1-p^2) ... check sign
    pyy = 2*p*py
    f,fy,fyy = F[i],FY[i],FYY[i]
    # V = num f^-2 ; V_y = (num_y - 2 num fy/f) f^-2 ; set a = fy/f
    a = fy/f; ay = fyy/f - a*a
    # num_y - 2 a num = 0 ; (num_y - 2 a num)_y = 0 (given first) -> num_yy - 2 ay num - 2 a num_y = 0
    # num = 1 + c p + d p^2/2 ; num_y = (c + d p) py ; num_yy = d py^2 + (c + d p) pyy
    A = np.array([[py - 2*a*p, p*py - a*p*p],
                  [pyy - 2*ay*p - 2*a*py, py*py + p*pyy - ay*p*p - 2*a*p*py]])
    b = np.array([2*a, 2*ay])
    c,d = np.linalg.solve(A,b)
    lnV,G,GT = si.tables(c,d); GTT=np.gradient(GT,TG); j=int(np.argmin(np.abs(TG-TH[i])))
    num = 1+c*PHI+d*PHI**2/2
    print("phi_i=%+.2f Theta_i=%+.4f: c=%.6f d=%.6f | G=%.1e GT=%.1e gamma=%+.3f | num min=%.3f | max eps_V throat side of i: %.3f, other side: %.3f | sign G beyond (throat side) %+d, (plus side) %+d"%(
        p,TH[i],c,d,G[j],GT[j],GTT[j]+3*G[j]*GT[j],num.min(),(G[j+50:-200]**2/2).max(),(G[200:j-50]**2/2).max(), np.sign(G[min(j+2000,len(G)-300)]), np.sign(G[max(j-2000,0)])))
