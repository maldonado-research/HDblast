import sys, numpy as np
sys.path.insert(0,'..'); sys.path.insert(0,'.')
import a1_analyze as A
s,ts=A.load(sys.argv[1]); d=A.derived(s,ts)
ie=A.reliable_end(d); print('reliable end idx',ie,'tau',d['H0tau'][ie], 'n',len(d['T']))
idx=np.unique(np.linspace(0,len(d['T'])-1,int(sys.argv[2]) if len(sys.argv)>2 else 60).astype(int))
print('%7s %7s %7s %8s %9s %9s %9s %9s %8s %8s %8s %8s %9s %8s'%('T','tau','lna','H','Wa4','Ra4','r','Om_r','Om_vac','lapse','phi_b','v','Ydot/4HR','Hnear'))
for i in idx:
    Y=s['params']['Y']; v=d['v_over_H0'][i]; H=d['H_over_H0'][i]; R=d['R'][i]
    src=Y*v*v/s['rho_b']**0 ; 
    print('%7.3f %7.2f %7.3f %8.4f %9.3e %9.3e %9.4f %8.3f %8.3f %8.3f %8.5f %9.2e %9.2e %8.1e'%(d['T'][i],d['H0tau'][i],d['ln_a'][i],H,d['Wa4'][i],d['Ra4'][i],d['r'][i],d['Omega_r'][i],d['Omega_vac'][i],d['lapse'][i],d['phi_b'][i],v,(Y*v*v/s['rho_b'])/(4*H*R) if R>0 and H!=0 else np.nan, d['Hmax_near'][i]))
