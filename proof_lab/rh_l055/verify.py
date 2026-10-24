"""Directed Decimal enclosures of a fixed polynomial RH source trial.

No quadrature, floating-point eigenvalue, or hidden-mode truncation is used
in the certificate.  See the accompanying proof note for analytic errors.
"""
from decimal import Decimal, Context, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction as F
from math import comb, factorial
from pathlib import Path
import json, sys, hashlib

PREC = int(sys.argv[1]) if len(sys.argv)>1 else 180
if PREC<100: raise ValueError('use at least 100 digits for this fixed polynomial packet')
DN=Context(prec=PREC,rounding=ROUND_FLOOR)
UP=Context(prec=PREC,rounding=ROUND_CEILING)
NEAR=Context(prec=PREC)

class I:
    __slots__=('lo','hi')
    def __init__(self,x=0,hi=None):
        if isinstance(x,I): self.lo,self.hi=x.lo,x.hi; return
        if isinstance(x,F):
            self.lo=DN.divide(Decimal(x.numerator),Decimal(x.denominator))
            self.hi=UP.divide(Decimal(x.numerator),Decimal(x.denominator)); return
        if isinstance(x,float): raise TypeError('float forbidden in certificate')
        self.lo=Decimal(x);self.hi=Decimal(hi) if hi is not None else self.lo
        if self.lo>self.hi: raise ValueError('unordered interval')
    def __add__(self,b):
        b=I(b);return I(DN.add(self.lo,b.lo),UP.add(self.hi,b.hi))
    __radd__=__add__
    def __neg__(self):return I(self.hi.copy_negate(),self.lo.copy_negate())
    def __sub__(self,b):return self+-I(b)
    def __rsub__(self,b):return I(b)+-self
    def __mul__(self,b):
        b=I(b);pairs=[(a,c) for a in (self.lo,self.hi) for c in (b.lo,b.hi)]
        return I(min(DN.multiply(a,c) for a,c in pairs),max(UP.multiply(a,c) for a,c in pairs))
    __rmul__=__mul__
    def __truediv__(self,b):
        b=I(b)
        if b.lo<=0<=b.hi:raise ZeroDivisionError('interval contains zero')
        return self*I(DN.divide(Decimal(1),b.hi),UP.divide(Decimal(1),b.lo))
    def __rtruediv__(self,b):return I(b)/self
    def __pow__(self,n):
        assert isinstance(n,int) and n>=0
        ans=I(1);base=self
        while n:
            if n&1:ans=ans*base
            base=base*base;n//=2
        return ans
    def ln(self):
        if self.lo<=0:raise ValueError('log domain')
        if self.lo==self.hi==1:return I(0)
        return I(NEAR.ln(self.lo).next_minus(DN),NEAR.ln(self.hi).next_plus(UP))
    def sqrt(self):
        if self.lo<0:raise ValueError('sqrt domain')
        return I(max(Decimal(0),NEAR.sqrt(self.lo).next_minus(DN)),NEAR.sqrt(self.hi).next_plus(UP))
    def json(self):return [str(self.lo),str(self.hi)]
    def __repr__(self):return f'I({self.lo},{self.hi})'

def add(p,q):
    out=list(p)+[0]*max(0,len(q)-len(p))
    for j,v in enumerate(q):out[j]=out[j]+v
    return out
def scale(p,a):return [v*a for v in p]
def mul(p,q):
    out=[0]*(len(p)+len(q)-1)
    for j,v in enumerate(p):
        for k,w in enumerate(q):out[j+k]=out[j+k]+v*w
    return out
def shift(p,a):
    out=[I(0) for _ in p];powers=[I(1)]
    for _ in range(len(p)-1):powers.append(powers[-1]*a)
    for j,v in enumerate(p):
        for k in range(j+1):out[k]+=I(v)*comb(j,k)*powers[j-k]
    return out
def legendre(n):return [F((-1)**(n+k)*comb(n,k)*comb(n+k,k)) for k in range(n+1)]
def inner_poly(p,q):return 2*L*sum((I(v)/ (k+1) for k,v in enumerate(mul(p,q))),I(0))

def atan_small(q,n=110):
    x=F(1,q);s=sum(((-1)**j*x**(2*j+1)/F(2*j+1) for j in range(n)),F(0))
    tail=x**(2*n+1)/F(2*n+1)
    return I(s)+I(-tail.numerator,tail.numerator)/tail.denominator

def bernoulli(n):
    A=[F(0)]*(n+1);out=[]
    for m in range(n+1):
        A[m]=F(1,m+1)
        for j in range(m,0,-1):A[j-1]=j*(A[j-1]-A[j])
        out.append(A[0])
    return out

L=F(11,20); PI=16*atan_small(5)-4*atan_small(239)
B=bernoulli(66); nn=1024
H=sum((F(1,k) for k in range(1,nn+1)),F(0))
GAMMA=I(H)-I(nn).ln()-F(1,2*nn)
GAMMA+=sum((I(B[2*k]/F(2*k*nn**(2*k))) for k in range(1,33)),I(0))
gerr=abs(B[66]/F(66*nn**66));GAMMA+=I(-gerr.numerator,gerr.numerator)/gerr.denominator
LN2=I(2).ln(); LN3=I(3).ln();LN4=2*LN2
c=-(GAMMA+(2*PI*L).ln())/2
u,v=LN2,LN3;a,b=u/I(2).sqrt(),v/I(3).sqrt()
assert u.lo>I(L).hi and v.hi<I(2*L).lo and I(2*L).hi<(2*u).lo

trial_path=Path(__file__).with_name('trial.json')
trial=json.loads(trial_path.read_text());M,N=trial['M'],trial['N']
assert N==5 and M==12
den=trial['denominator']
coeff=[[F(x,den) for x in row] for row in trial['legendre_numerators']]
assert all(coeff[j][k]==(1 if j==k else 0) for j in range(N) for k in range(N))
Us=[legendre(2*j) for j in range(M)]
Z=[];JZ=[]
for k in range(N):
    z=[F(0)];jz=[F(0)]
    for j in range(M):
        z=add(z,scale(Us[j],coeff[j][k]))
        hj=sum((F(1,l) for l in range(1,2*j+1)),F(0))/2
        jz=add(jz,scale(Us[j],coeff[j][k]*hj))
    Z.append(z);JZ.append(jz)

# Exact rational Taylor coefficients of j(z)=e^(z/2)/(2 sinh z)-1/(2z).
K=48; ff=[]
for n in range(K+2):
    f=F(1,2**n*factorial(n))
    for q in range(1,n//2+1):f-=F(2,factorial(2*q+1))*ff[n-2*q]
    ff.append(f/2)
jcoef=ff[1:]
assert jcoef[:5]==[F(1,4),-F(1,48),-F(1,32),F(7,11520),F(5,1536)]

def kernel_poly(p):
    out=[F(0)]*(len(p)+K+1)
    for m,jm in enumerate(jcoef):
        fac=-L*jm*(2*L)**m
        for j,pj in enumerate(p):
            if not pj:continue
            f=fac*pj
            for r in range(m+1):out[r]+=f*F((-1)**r*comb(m,r),m-r+j+1)
            if m%2:out[m+j+1]+=f*F(2*factorial(m)*factorial(j),factorial(m+j+1))
    return out

cosh_poly=[F(0)]
for r in range(9):
    n=2*r
    term=[F(comb(n,k))*L**n*(-F(1,2))**(n-k)/factorial(n) for k in range(n+1)]
    cosh_poly=add(cosh_poly,term)

P=[]
for k,z in enumerate(Z):
    p=add(JZ[k],kernel_poly(z))
    p=add([I(x) for x in p],scale(z,c))
    p=add(p,scale(cosh_poly,inner_poly(cosh_poly,z)))
    P.append(p)

def hm(n,p=1):return sum((F(1,k**p) for k in range(1,n+1)),F(0))
maxdegree=2*max(map(len,P))+2
Wmom=[];W2mom=[]
for r in range(maxdegree+1):
    n=r+1;h=hm(n);h2=hm(n,2)
    sy=-I(F(1,n*n));st=-I(h/n)
    yy=I(F(2,n**3));tt=I((h*h+h2)/n)
    yt=I(h/F(n*n))+I(h2/n)-PI*PI/F(6*n)
    Wmom.append(-(LN4/n+sy+st)/4)
    W2mom.append((LN4*LN4/n+2*LN4*(sy+st)+yy+tt+2*yt)/16)

def weighted(p,kind):
    moments=Wmom if kind==1 else W2mom
    return 2*L*sum((I(v)*moments[j] for j,v in enumerate(p)),I(0))

def primitives(x,d):
    powers=[I(1)]
    for _ in range(d+1):powers.append(powers[-1]*x)
    ordinary=[powers[j+1]/(j+1) for j in range(d+1)]
    iszero=x.lo==x.hi==0;isone=x.lo==x.hi==1
    lx=I(0) if iszero else x.ln()
    l1=I(0) if isone else (1-x).ln()
    accum=I(0);w=[]
    for j in range(d+1):
        n=j+1;accum+=powers[n]/n
        logy=I(0) if iszero else powers[n]*(lx/n-F(1,n*n))
        log1=-I(hm(n)/n) if isone else ((powers[n]-1)*l1-accum)/n
        w.append(-(LN4*ordinary[j]+logy+log1)/4)
    return ordinary,w

edges=[I(0),1-v/(2*L),1-u/(2*L),u/(2*L),v/(2*L),I(1)]
assert all(x.hi<y.lo for x,y in zip(edges,edges[1:]))
prims=[primitives(x,maxdegree) for x in edges]

def interval_int(p,seg,weight=0):
    left=prims[seg][weight];right=prims[seg+1][weight]
    return 2*L*sum((I(z)*(right[k]-left[k]) for k,z in enumerate(p)),I(0))

pieces=[]
for seg in range(len(edges)-1):
    row=[]
    for k,z in enumerate(Z):
        p=list(P[k])
        if seg<=1:p=add(p,scale(shift(z,u/(2*L)),-a/2))
        if seg==0:p=add(p,scale(shift(z,v/(2*L)),-b/2))
        if seg>=3:p=add(p,scale(shift(z,-u/(2*L)),-a/2))
        if seg==4:p=add(p,scale(shift(z,-v/(2*L)),-b/2))
        row.append(p)
    pieces.append(row)

def pairing(poly,k):
    return weighted(mul(poly,Z[k]),1)+sum((interval_int(mul(poly,pieces[s][k]),s) for s in range(5)),I(0))

A=[[pairing(Z[i],j) for j in range(N)] for i in range(N)]
project=[[pairing(Us[i],j) for j in range(N)] for i in range(N)]
Gram=[[I(0) for _ in range(N)] for _ in range(N)]
for i in range(N):
    for j in range(i,N):
        value=weighted(mul(Z[i],Z[j]),2)
        for s in range(5):
            value+=interval_int(add(mul(Z[i],pieces[s][j]),mul(Z[j],pieces[s][i])),s,1)
            value+=interval_int(mul(pieces[s][i],pieces[s][j]),s)
        for k in range(N):value-=project[k][i]*project[k][j]/F(2*L,4*k+1)
        Gram[i][j]=Gram[j][i]=value

tz=sum((inner_poly(z,z) for z in Z),I(0))
tr=sum((Gram[k][k] for k in range(N)),I(0))
assert tz.lo>0 and tr.lo>=0
ej=I(F(100,9)*F(11,20)**(K+1))
ec=I(2*F(11,40)**18/F(factorial(18)))
epsilon=L*ej+8*L*ec
kappa=F(6300,1871)
delta=epsilon*tz+kappa*(2*epsilon*(tz*tr).sqrt()+epsilon*epsilon*tz)
# Freeze one rational outward scalar allowance; interval inputs cannot shrink it.
delta=I(delta.hi)
lower=[[((A[i][j]+A[j][i])/2-kappa*Gram[i][j]-(delta if i==j else 0)) for j in range(N)] for i in range(N)]

def ldlt(mat,shift=F(0)):
    n=len(mat);ll=[[I(0) for _ in range(n)] for _ in range(n)];dd=[]
    for i in range(n):
        d=mat[i][i]-shift-sum((ll[i][k]*ll[i][k]*dd[k] for k in range(i)),I(0))
        if d.lo<=0:return False,dd+[d]
        dd.append(d);ll[i][i]=I(1)
        for j in range(i+1,n):ll[j][i]=(mat[j][i]-sum((ll[j][k]*ll[i][k]*dd[k] for k in range(i)),I(0)))/d
    return True,dd

passed,pivots=ldlt(lower)
safe_shift=F(1,10**10)
shift_pass,shift_pivots=ldlt(lower,safe_shift)

# A second sign check uses only exact Fraction arithmetic on a rational
# matrix dominated by every matrix in the directed entry enclosure.
grid=10**18
rat=[[F(0) for _ in range(N)] for _ in range(N)]
for i in range(N):
    for j in range(i,N):
        midpoint=(F(lower[i][j].lo)+F(lower[i][j].hi))/2
        z=(midpoint*grid).__floor__()
        rat[i][j]=rat[j][i]=F(z,grid)
row_errors=[]
for i in range(N):
    row_errors.append(sum((max(abs(F(lower[i][j].lo)-rat[i][j]),abs(F(lower[i][j].hi)-rat[i][j])) for j in range(N)),F(0)))
round_allow=F(1,10**17)
assert max(row_errors)<round_allow

def exact_ldlt(mat):
    n=len(mat);ll=[[F(0) for _ in range(n)] for _ in range(n)];dd=[]
    for i in range(n):
        d=mat[i][i]-sum((ll[i][k]**2*dd[k] for k in range(i)),F(0))
        if d<=0:return False,dd+[d]
        dd.append(d);ll[i][i]=1
        for j in range(i+1,n):ll[j][i]=(mat[j][i]-sum((ll[j][k]*ll[i][k]*dd[k] for k in range(i)),F(0)))/d
    return True,dd

rat_test=[[rat[i][j]-(round_allow+safe_shift if i==j else 0) for j in range(N)] for i in range(N)]
exact_pass,exact_pivots=exact_ldlt(rat_test)
assert exact_pass and passed and shift_pass

# All-function transfer: E=Z-C^(-1)R has norm below 8/5.
assert tz.hi<I(F(17,10)).lo and tr.hi<I(F(1,250)).lo
assert (tr.sqrt()+epsilon*tz.sqrt()).hi<I(F(8,125)).lo
assert I(kappa).hi<I(F(7,2)).lo
assert F(4,3)+F(7,2)*F(8,125)<F(8,5)
assert F(128,25)*F(1,10**11)<safe_shift
assert F(2,10**11)<F(1871,6300)

# Adversarial controls for the sign/uncertainty decisions.
controls={}
controls['positive_diagonal_indefinite_matrix_rejected']=not exact_ldlt([[F(1),F(2)],[F(2),F(1)]])[0]
controls['visible_only_false_positive_rejected']=not exact_ldlt([[F(1)-F(2)**2]])[0]
controls['zero_pivot_rejected']=not exact_ldlt([[F(0)]])[0]
bad=[row[:] for row in rat_test];bad[0][0]=-1
controls['negative_direction_rejected']=not exact_ldlt(bad)[0]
controls['inflated_uncertainty_rejected']=not ldlt(lower,F(1,100000))[0]
try:I(1)/I(-1,1)
except ZeroDivisionError:controls['singular_interval_inverse_rejected']=True
else:controls['singular_interval_inverse_rejected']=False
try:I(0.1)
except TypeError:controls['binary_float_input_rejected']=True
else:controls['binary_float_input_rejected']=False
assert all(controls.values())

result={'status':'PASS_FIXED_WINDOW_LOWER_MATRIX' if passed else 'INCONCLUSIVE',
    'python':sys.version.split()[0],'decimal_precision':PREC,'L':'11/20','retained_modes':N,'trial_modes':M,
    'kernel_degree':K,'cosh_degree':16,'trial_sha256':hashlib.sha256(trial_path.read_bytes()).hexdigest(),
    'operator_error':epsilon.json(),'lower_scalar_allowance':delta.json(),'trial_norm_trace':tz.json(),
    'residual_gram_trace':tr.json(),'lower_matrix':[[x.json() for x in row] for row in lower],
    'ldlt_pivots':[x.json() for x in pivots],
    'shift':'1/10000000000','shift_pass':shift_pass,'shifted_pivots':[x.json() for x in shift_pivots],
    'rational_grid':grid,'rational_matrix_numerators':[[int(x*grid) for x in row] for row in rat],
    'rational_rounding_allowance':str(round_allow),'rational_row_error_upper':[str(x) for x in row_errors],
    'exact_rational_shift_pass':exact_pass,'exact_shifted_pivots':[str(x) for x in exact_pivots],
    'full_family_window':'0<L<=11/20','full_family_lower_bound':'1/100000000000',
    'negative_controls':controls,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'claim':'Fixed arithmetic source at L=11/20 only; not all-width positivity or RH.'}
out=Path(__file__).with_name(f'result-{PREC}.json');out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['status','decimal_precision','operator_error','lower_scalar_allowance','shift_pass']}),flush=True)
print(json.dumps({'exact_rational_shift_pass':exact_pass,'full_family_lower_bound':'1/100000000000','negative_controls':controls}),flush=True)
