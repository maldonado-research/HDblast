"""Independent fast-time source solution for the heavy, narrow family."""
import math


def bump(u, width):
    def step(x):
        if x <= 0 or x >= 1:
            return float(x >= 1)
        z=1/(1-x)-1/x
        ex=math.exp(-abs(z))
        return 1/(1+ex) if z>=0 else ex/(1+ex)
    up = step(2*u/width)
    down = step(2*u/width-1)
    return up-down


def solve(incoming, steps):
    # u is the independent variable. Eliminate its positive velocity using energy conservation.
    width, amplitude, mass = .24, .008, 1.9
    h = width/steps
    c, momentum = 0., 0.
    def rhs(u, c, momentum):
        delta = c-amplitude*bump(u,width)
        velocity=math.sqrt(incoming**2-momentum**2-(mass*delta)**2)
        return momentum/velocity, -mass**2*delta/velocity
    for i in range(steps):
        u = h*i
        k1,l1=rhs(u,c,momentum)
        k2,l2=rhs(u+h/2,c+h*k1/2,momentum+h*l1/2)
        k3,l3=rhs(u+h/2,c+h*k2/2,momentum+h*l2/2)
        k4,l4=rhs(u+h,c+h*k3,momentum+h*l3)
        c+=h*(k1+2*k2+2*k3+k4)/6
        momentum+=h*(l1+2*l2+2*l3+l4)/6
    return .5*(momentum**2+(mass*c)**2)


def audit(joined, runs):
    primary=joined["runs"]["primary"]
    samples=(dict(zip(joined["trace_columns"],row)) for row in primary["trace"])
    conserved=[s["pi"]*s["a"]**3 for s in samples
               if .012<s["psi"]<.03 and s["pi"]>0. and s["chi"]==0.]
    # Flat pre-drive scalar plateau obeys d(a^3 psi_dot)/dt=0. The primary
    # mesh sets the incoming datum; other spatial meshes have small differences.
    incoming=(sum(conserved)/len(conserved))/primary["drive_begin_radius"]**3
    e=solve(incoming,2048)
    return dict(incoming_velocity=incoming,pre_drive_momentum_relative_range=max(conserved)/min(conserved)-1,
                outgoing_parent_energy=e,half_resolution_relative=solve(incoming,1024)/e-1,
                match_energy_relative={n:r['handover_action']*r['drive_begin_mass']/
                                       (r['drive_begin_radius']**3*e)-1
                                       for n,r in runs.items() if n.startswith(('mass_4e10','mass_4e12'))})
