"""Universal sphere-quotient checks for all quadratic angular source products."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
sys.path.insert(0,str(Path(__file__).resolve().parent))
import angular_exact as e
import source_adapters as a
from fractions import Fraction as Q

def validate():
    DY=e.grad(e.Y);HY=e.hess(e.Y)
    basis=[e.P,e.tmul(e.Z,e.P),e.tmul(e.F,e.P),*e.HESS,*e.make_tt()]
    square={(i,j):e.add(*(e.mul(HY[i,k],HY[k,j]) for k in range(4))) for i in range(4) for j in range(4)}
    rows=(('DY_DY',e.outer(DY,DY),a.DYDY),('Y_HESS',e.tmul(e.Y,HY),a.YHESS),('HESS_SQUARE',square,a.HESS_SQUARE))
    for label,target,coeff in rows:e.expect(target,e.tadd(*(e.tscale(c,b) for c,b in zip(coeff,basis))),'SOURCE_'+label)
    scalars=(('Y2',e.mul(e.Y,e.Y),a.Y2),('GRAD2',e.dot(DY,DY),a.GRAD2),('HESS2',e.dot(HY,HY),a.HESS2))
    for label,target,coeff in scalars:e.expect(e.scalar(target),e.tadd(*(e.tscale(c,b) for c,b in zip(coeff,e.SCALARS))),'SOURCE_'+label)
    hdy={(i,):e.add(*(e.mul(HY[i,j],DY[j,]) for j in range(4))) for i in range(4)}
    e.expect(e.tmul(e.Y,DY),e.tadd(e.tscale(Q(1,4),e.GRADS[0]),e.tscale(Q(1,2),e.GRADS[1])),'SOURCE_Y_DY')
    e.expect(hdy,e.tadd(e.GRADS[0],e.tscale(-2,e.GRADS[1])),'SOURCE_HESS_DY')
    D3=e.covd(HY)
    rhs={(c,i,j):e.add(e.scale(-1,e.mul(e.P[c,i],DY[j,])),e.scale(-1,e.mul(e.P[c,j],DY[i,])),e.scale(-2,e.mul(e.P[i,j],DY[c,]))) for c in range(4) for i in range(4) for j in range(4)}
    e.expect(D3,rhs,'SOURCE_THIRD_ANGULAR_JET')
    D4=e.covd(D3)
    rhs={(d,c,i,j):e.add(e.scale(-1,e.mul(e.P[c,i],HY[d,j])),e.scale(-1,e.mul(e.P[c,j],HY[d,i])),e.scale(-2,e.mul(e.P[i,j],HY[d,c]))) for d in range(4) for c in range(4) for i in range(4) for j in range(4)}
    e.expect(D4,rhs,'SOURCE_FOURTH_ANGULAR_JET')
    print('PASS: 10 universal angular identities, including third/fourth derivative tensors.')
    return 10
if __name__=='__main__':validate()
