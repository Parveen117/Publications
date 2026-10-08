"""ECP1 exact controls. Requires Python 3.12 and sympy."""
import sys,json,hashlib
from pathlib import Path
from fractions import Fraction as Q
from math import factorial
import sympy as S
assert sys.version_info[:2]==(3,12)
# Free associative words; t,s degrees remain central and explicit.
def mul(p,q,cap=4):
 out={}
 for (i,j,w),a in p.items():
  for (k,l,v),b in q.items():
   if i+j+k+l<=cap:
    key=(i+k,j+l,w+v);out[key]=out.get(key,Q(0))+a*b
 return {k:v for k,v in out.items() if v}
def exp(letter,axis,sign):
 return {(n if axis==0 else 0,n if axis==1 else 0,letter*n):Q(sign**n,factorial(n)) for n in range(5)}
loop={(0,0,''):Q(1)}
for letter,axis,sign in [('A',0,1),('B',1,1),('A',0,-1),('B',1,-1)]:loop=mul(loop,exp(letter,axis,sign))
lead={k:v for k,v in loop.items() if k[0]+k[1]<=2}
assert lead=={(0,0,''):Q(1),(1,1,'AB'):Q(1),(1,1,'BA'):Q(-1)}
assert all(i and j for i,j,w in loop if i+j>0)
A=S.Matrix([[0,1],[0,0]]);B=S.Matrix([[0,0],[1,0]]);I=S.eye(2)
t,s,x,y,a=S.symbols('t s x y a',real=True)
L=(I+t*A)*(I+s*B)*(I-t*A)*(I-s*B)
assert S.simplify(L-S.Matrix([[1+s*t+s*s*t*t,-s*t*t],[s*s*t,1-s*t]]))==S.zeros(2)
assert A*B-B*A==S.diag(1,-1) # rational coefficients can be curved
for n in range(1,9):
 h=S.Rational(1,n); assert L.subs({t:h,s:h})!=I
# Scalar response witness and a blind reading.
DX=lambda f:S.diff(f,x)
DY=lambda f:S.exp(a*x)*S.diff(f,y)
for f in [y,x,x*y,x*x+y*y]:
 assert S.simplify(DX(DY(f))-DY(DX(f))-a*DY(f))==0
assert S.simplify(DX(DY(x))-DY(DX(x)))==0
assert S.simplify(DX(DY(y))-DY(DX(y)))==a*S.exp(a*x)
# Pure gauge: noncommutation canceled by derivative terms.
g=(I+x*A)*(I+y*B);Xs=g.inv()*g.diff(x);Ys=g.inv()*g.diff(y);comm=Xs*Ys-Ys*Xs
assert comm!=S.zeros(2)
assert S.simplify(Ys.diff(x)-Xs.diff(y)+comm)==S.zeros(2)
q=S.symbols('q');F=q*(Ys.diff(x)-Xs.diff(y))+q*q*comm
assert S.simplify(F-q*(q-1)*comm)==S.zeros(2)
# The draft's actual affine words; no manufactured connection.
y0=S.Function('y0')(x);f=S.Function('f')(x)
T=lambda z,c:(S.diff(z,x)+c*y0,0)
Dh=lambda z,c:(S.diff(z,x),0)
Sh=lambda z,c:(z+c*y0,c)
abc=T(*T(f,1))[0]+y0;adc=T(*Dh(*Sh(f,1)))[0]+y0
assert S.simplify(abc-adc)==S.diff(y0,x)-S.diff(y0,x,2)
C1,C2=S.symbols('C1 C2');sol=C1*S.exp(x)+C2
assert S.diff(sol,x)-S.diff(sol,x,2)==0
# Negative controls must fail their proposed identities.
false=[lead=={(0,0,''):Q(1)},comm==S.zeros(2),S.simplify(F-q*q*comm)==S.zeros(2),A*B-B*A==S.zeros(2),S.simplify(DX(DY(y))-DY(DX(y)))==0,S.simplify(abc-adc)==0]
assert not any(false)
root=Path(__file__).parent
result={'status':'PASS_SCOPED_ECP1','python':sys.version.split()[0],'sympy':S.__version__,'groups':7,'rational_loop_controls':8,'rejected_false_identities':len(false),'sha256':{p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in ['THEOREM.md','verify.py']},'limitations':'Written general proofs plus exact finite/formal controls; no proof-assistant, independent review or physical validation.'}
(root/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
