"""A finite comoving band of minimally coupled Gaussian scalar modes on the B1 wall."""
import math
import numpy as np

G = 100.0
PHI_STAR = 0.5
B = 6.49e-6


class QuantumShell:
    kind = "quantum"
    pure_radiation = False
    Y = 0.0

    def __init__(self, rb, loga_ref, nodes, feedback, decay):
        x, w = np.polynomial.legendre.leggauss(nodes)
        self.k = 4*(1+x)
        self.weight = 4*w*self.k**2/(2*np.pi**2)
        self.rb = rb
        self.ref = loga_ref
        self.mult = B/rb
        self.feedback = feedback
        self.decay = decay
        self.initial_corner = None

    def init_state(self):
        # Prior to the declared preparation slice the quantum sector is absent.
        return np.zeros(1)

    def background(self, Ab, phi):
        a = math.exp(Ab-self.ref)
        p2 = (self.k/a)**2
        delta = phi-PHI_STAR
        m2 = G**2*delta**2
        e = np.sqrt(p2+m2)
        return a, p2, delta, m2, e

    def prepare(self, Ab, phi):
        a, _, _, _, e = self.background(Ab, phi)
        return np.concatenate(([0.,0.,0.], .5/(a**3*e), a**3*e*.5, np.zeros_like(e)))

    def view(self, state):
        return state[3:].reshape(3,-1)

    def quantities(self, state, Ab, phi):
        a, p2, d, m2, e = self.background(Ab, phi)
        x, pp, xp = self.view(state)
        x0 = .5/(a**3*e)
        pp0 = .5*a**3*e
        dx = x-x0
        dp = (pp-pp0)/a**6
        rho = np.dot(self.weight, .5*(dp+(p2+m2)*dx))
        pr = np.dot(self.weight, .5*(dp-(m2+p2/3)*dx))
        j = G**2*d*np.dot(self.weight,dx)
        return (rho,pr,j), (a,p2,d,m2,e,x,pp,xp,x0,pp0,dx,dp)

    def totals(self, m, Ab, phib, **kwargs):
        if len(m)==1 or not self.feedback:
            return 0.,0.,0.
        (rho,pr,j), _ = self.quantities(m,Ab,phib)
        return self.mult*(m[0]+rho), self.mult*(m[0]/3+pr), self.mult*j

    def rhs(self, m, Ab, phib, AT, phiT, eB, T=None):
        if len(m)==1:
            return np.zeros(1)
        (rho,pr,j), z = self.quantities(m,Ab,phib)
        a,p2,d,m2,e,x,pp,xp,x0,pp0,dx,dp = z
        lapse=eB/self.rb
        gamma = self.decay**2*m2/(8*np.pi*e)
        power = np.dot(self.weight,.5*gamma*(dp+e**2*dx))
        rad = lapse*power - 4*AT*m[0]
        W = a**3*j*phiT
        pressure = 3*a**3*AT*(pr+m[0]/3)
        return np.concatenate(([rad,W,pressure], lapse*(2*xp/a**3-gamma*dx),
                               lapse*(-2*a**3*e**2*xp-gamma*(pp-pp0)),
                               lapse*(pp/a**3-a**3*e**2*x-gamma*xp)))

    def totals_T(self, m, dm, Ab, phi, AT, phiT, eB, T=None):
        if len(m)==1 or not self.feedback:
            return 0.,0.,0.
        (rho,pr,j),z=self.quantities(m,Ab,phi)
        a,p2,d,m2,e,x,pp,xp,x0,pp0,dx,dp = z
        xd,ppd,_ = self.view(dm)
        m2d=2*G**2*d*phiT
        p2d=-2*AT*p2
        gd=3*AT+.5*(m2d+p2d)/e**2
        dxd=xd+x0*gd
        dpd=(ppd-pp0*gd)/a**6-6*AT*dp
        rd=.5*np.dot(self.weight,dpd+(m2d+p2d)*dx+e**2*dxd)
        pd=.5*np.dot(self.weight,dpd-(m2d+p2d/3)*dx-(m2+p2/3)*dxd)
        jd=G**2*(phiT*np.dot(self.weight,dx)+d*np.dot(self.weight,dxd))
        return self.mult*(dm[0]+rd), self.mult*(dm[0]/3+pd), self.mult*jd

    def extras(self, m, Ab, phib, AT, phiT, eB, T=None):
        if len(m)==1:
            return {}
        (rho,pr,j),z=self.quantities(m,Ab,phib)
        a,p2,d,m2,e,x,pp,xp,x0,pp0,dx,dp = z
        energy_modes=.5*(dp+e**2*dx)
        n=energy_modes*a**3/e
        gamma=self.decay**2*m2/(8*np.pi*e)
        return dict(rho_hat=rho,p_hat=pr,j_hat=j,n_hat=float(np.dot(self.weight,n))/a**3,comoving_number=float(np.dot(self.weight,n)),m_hat=math.sqrt(m2),
                    Q_hat=float(np.dot(self.weight,gamma*energy_modes)),rad_hat=m[0],
                    scalar_work=m[1],pressure_work=m[2],comoving_energy=a**3*(m[0]+rho),
                    max_physical_frequency=math.hypot(8/a,math.sqrt(m2)),
                    min_uncertainty=float(np.min(x*pp-xp*xp)),min_occupation=float(n.min()))
