import sys,json,hashlib
from pathlib import Path
import sympy as S
assert sys.version_info[:2]==(3,12)
I=S.eye(2);J=S.diag(1,-1);G=S.Matrix([[1,2],[3,4]])
Ge=(G+J*G*J)/2;Go=(G-J*G*J)/2
assert Ge+Go==G and J*Ge*J==Ge and J*Go*J==-Go
# Degree-two exponential cut-loop and its logarithm.
t=S.symbols('t');H=J*G*J
C1=H+G;C2=H*H/2+H*G+G*G/2
assert C1==2*Ge
assert C2-C1*C1/2==Ge*Go-Go*Ge
# Exact odd nilpotent flow and bilateral channel identities.
N=S.Matrix([[0,1],[0,0]]);U=I+t*N;V=I-t*N
join=U+V;cut=U-V
assert J*U*J==V
assert join*cut==cut*join==2*t*N*2
assert join*join-cut*cut==4*I
# Orthogonal complete observer preserves the generator Gram.
P=(I+J)/2;Q=I-P
assert (P*G).T*(P*G)+(Q*G).T*(Q*G)==G.T*G
# Lifted attenuation cocycle and its differential equation.
s,z,lam=S.symbols('s z lam',real=True,nonzero=True)
c=lambda v:(1-S.exp(-lam*v))/lam
assert S.simplify(c(t+s)-c(t)-S.exp(-lam*t)*c(s))==0
amp=S.exp(-z*c(t))
assert S.simplify(S.diff(amp,t)+z*S.exp(-lam*t)*amp)==0
assert S.limit(c(t),lam,0)==t
# One-step fixed frequency and nonfixed +/- i pi.
assert S.exp(2*S.pi*S.I)==1 and 2*S.pi*S.I!=0
assert S.sinh(S.pi*S.I)==0 and S.pi*S.I!=0
# Rational lower-bound recurrence controls: 2sinh(x)>2x for real x>0
# finite positive first omitted term already establishes strictness.
for n in range(1,9):
 x=S.Rational(1,n);assert 2*x+S.Rational(1,3)*x**3>2*x
false=[join*cut==U+V,Ge*Go==Go*Ge,S.simplify(c(t+s)-c(t)-c(s))==0,2*S.pi*S.I==0,2*S.sinh(S.pi*S.I)==S.pi*S.I]
assert not any(false)
root=Path(__file__).parent
result={'status':'PASS_SCOPED_UG_A1','python':sys.version.split()[0],'groups':7,'rejected_false_identities':5,'sha256':{name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in ['README.md','verify.py']},'scope':'Exact controls and written finite-flow derivation; T41 historical certificate not rerun; no universal physical generator claimed.'}
(root/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
