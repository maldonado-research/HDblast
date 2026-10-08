"""Exact signed finite-weight endpoint error operators; no I/O or source."""
from fractions import Fraction as Q

PI = Q(14488038916154245685,4611686018427387904)
EPSILON = Q(3777893186295716171,37778931862957161709568)
SCALE = 1<<96

def require(ok,text):
    if not ok:raise ValueError(text)

def add(a,b):return (a[0]+b[0],a[1]+b[1])
def scale(a,q):return (q*a[0],q*a[1]) if q>=0 else (q*a[1],q*a[0])
def absolute(a):return max(abs(a[0]),abs(a[1]))
def intersect(a,b):
    out=(max(a[0],b[0]),min(a[1],b[1]))
    require(out[0]<=out[1],'canonical correlation/target intervals disjoint')
    return out

def difference(saved,rectangle):
    require(type(saved) is Q,'exact normalized saved component required')
    require(type(rectangle) is list and len(rectangle)==2 and all(type(v) is int for v in rectangle)
            and rectangle[0]<=rectangle[1],'ordered dyadic96 interval required')
    return (saved-Q(rectangle[1],SCALE),saved-Q(rectangle[0],SCALE))

def node_terms(k,weight,endpoint,saved_u,saved_w,target,incoming_u,incoming_w):
    require(type(k) is Q and 0<k<256 and type(weight) is Q and weight>0,'positive exact k/weight required')
    require(type(endpoint) is Q and endpoint in (Q(-4),Q(-7,2)),'registered endpoint required')
    require(all(type(v) is Q for v in (*saved_u,*saved_w,*incoming_u,*incoming_w)),
            'exact normalized saved/incoming state required')
    du=[difference(saved_u[j],target['U'][axis]) for j,axis in enumerate(('real','imag'))]
    dw=[difference(saved_w[j],target['W'][axis]) for j,axis in enumerate(('real','imag'))]
    L=-1/endpoint; q=weight*k/(4*PI*PI);mu=2*q*k
    R=add(add(scale(du[0],q*(2*k*k+3*L*L)),scale(dw[0],-q*L)),scale(dw[1],-q*k))
    P=add(add(scale(du[0],q*(2*k*k/3-L*L)),scale(dw[0],-q*L)),scale(dw[1],-q*k))
    drift=(saved_u[0]-saved_w[1]/(2*k))-(incoming_u[0]-incoming_w[1]/(2*k))
    intersect(du[0],add((drift,drift),scale(dw[1],Q(1)/(2*k))))
    q0=weight/(4*PI*PI)
    RC=add(add((q0*k*(2*k*k+3*L*L)*drift,)*2,scale(dw[1],q0*3*L*L/2)),scale(dw[0],-q0*k*L))
    PC=add(add((q0*k*(2*k*k/3-L*L)*drift,)*2,scale(dw[1],-q0*(2*k*k/3+L*L/2))),scale(dw[0],-q0*k*L))
    R=intersect(R,RC);P=intersect(P,PC)
    return {'delta_U':du,'delta_W':dw,'R_error':R,'P_error':P,
            'weighted_U_L1_upper':mu*sum(map(absolute,du)),
            'weighted_W_L1_upper':mu*sum(map(absolute,dw)),
            'R_triangle_upper':absolute(R),'P_triangle_upper':absolute(P),
            'U_L1_upper':sum(map(absolute,du)),'W_L1_upper':sum(map(absolute,dw)),
            'canonical_drift':drift,'canonical_drift_weighted_signed':mu*drift,
            'canonical_drift_weighted_abs':mu*abs(drift),
            'first_order_wronskian_drift_abs':2*EPSILON*abs(drift)}

class Accumulator:
    def __init__(self):
        self.count=0;self.last=Q(0);self.sums={k:Q(0) for k in ('weight_sum','weighted_U_L1_upper',
            'weighted_W_L1_upper','R_triangle_upper','P_triangle_upper',
            'canonical_drift_weighted_signed','canonical_drift_weighted_abs')}
        self.intervals={'R_error':(Q(0),Q(0)),'P_error':(Q(0),Q(0))}
        self.maxima={'U_L1_upper':Q(0),'W_L1_upper':Q(0),'first_order_wronskian_drift_abs':Q(0)}
    def add(self,k,weight,terms):
        require(k>self.last,'monotonically ordered independent capsule required')
        self.last=k;self.count+=1;self.sums['weight_sum']+=weight
        for field in self.sums:
            if field!='weight_sum':self.sums[field]+=terms[field]
        for field in self.intervals:self.intervals[field]=add(self.intervals[field],terms[field])
        for field in self.maxima:self.maxima[field]=max(self.maxima[field],terms[field])
    def snapshot(self,capsule,K,endpoint):
        return {'case':f'{capsule}/{K}/{endpoint}','capsule':capsule,'K':str(K),'endpoint':str(endpoint),
                'node_count':self.count,'last_k':str(self.last),'sums':dict(self.sums),
                'signed_intervals':dict(self.intervals),'maxima':dict(self.maxima),
                'incoming_state_error':'ZERO_BY_EXACT_REPRESENTED_ANCHOR_DEFINITION',
                'contact_difference':'EXACT_ZERO_UNDER_IDENTICAL_COMPLETE_CONTACTS',
                'full_momentum_and_continuous_time_integrals':'UNRESOLVED'}
