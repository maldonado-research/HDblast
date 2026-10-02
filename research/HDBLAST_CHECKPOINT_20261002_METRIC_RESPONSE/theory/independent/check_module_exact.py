import sympy as s
import metric_wkb as mw
mw.HALF=s.Rational(1,2)
mw.ONE=s.Integer(1)
L,w=s.symbols('L w',positive=True)
h=s.symbols('h0:5',real=True)
a=tuple(s.factorial(n)*L**(n+1) for n in range(5))
da=tuple(sum(s.binomial(n,j)*a[n-j]*h[j] for j in range(n+1)) for n in range(5))
x=mw.generic_subtractions(w*w-2*L*L,a,da)
for name in ['R','P','S','W2','W4']:
 print(name, type(x[name]['value']),s.factor(x[name]['value']))
print('done')
