"""Directed Decimal enclosures of a fixed polynomial RH source trial at a
declared support width L with every active prime power retained.

General-width version of the L=11/20 packet (proof_lab/rh_l055): the
support half-width L, the set of active prime powers, the number N of
retained Legendre modes, the kernel and cosh degrees and the
complement comparison levels are read from config.json. The trial is
read from trial.json. No quadrature, floating-point eigenvalue or
hidden-mode truncation is used in the certificate. See the accompanying
proof note for the analytic error bounds and the complement comparison.

Usage: python3 verify.py [precision] --packet DIR
           certificate run for the packet directory DIR, which holds
           config.json and trial.json; result-<precision>.json is written there
       python3 verify.py [precision] --packet DIR --rayleigh
           two-sided enclosure of the Rayleigh quotients of the trial columns
       python3 verify.py [precision] --packet DIR --matrices OUT.json
           operator matrices for an identity trial (exploration aid, not a
           certificate)
Without --packet the directory of this file is used.
"""
from decimal import Decimal, Context, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction as F
from math import comb, factorial
from pathlib import Path
import json, sys, hashlib

argv=list(sys.argv[1:])
PACKET=Path(__file__).parent
if '--packet' in argv:
    i=argv.index('--packet');PACKET=Path(argv[i+1]).resolve();del argv[i:i+2]
if '--matrices' in argv:
    i=argv.index('--matrices');MATRIX_OUT=argv[i+1];del argv[i:i+2]
args=[a for a in argv if not a.startswith('--')]
flags=[a for a in argv if a.startswith('--')]
PREC=int(args[0]) if args else 180
MATRIX_MODE='MATRIX_OUT' in globals()
RAYLEIGH_MODE='--rayleigh' in flags
if PREC<100 and not MATRIX_MODE: raise ValueError('use at least 100 digits for a certificate run')
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
    def exp(self):
        return I(NEAR.exp(self.lo).next_minus(DN),NEAR.exp(self.hi).next_plus(UP))
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

def hm(n,p=1):return sum((F(1,k**p) for k in range(1,n+1)),F(0))

# ---------------------------------------------------------------- config
here=PACKET
config=json.loads((here/'config.json').read_text())
L=F(config['L']); assert 0<L<1
N=int(config['retained_modes']); K=int(config['kernel_degree']); COSH_TERMS=int(config['cosh_terms'])
complement=config['complement']
HLEV=int(complement['explicit_levels']); assert HLEV>=0
MARGIN=F(config['margin']); assert MARGIN>0

def inner_poly(p,q):return 2*L*sum((I(v)/ (k+1) for k,v in enumerate(mul(p,q))),I(0))

PI=16*atan_small(5)-4*atan_small(239)
B=bernoulli(66); nn=1024
H=sum((F(1,k) for k in range(1,nn+1)),F(0))
GAMMA=I(H)-I(nn).ln()-F(1,2*nn)
GAMMA+=sum((I(B[2*k]/F(2*k*nn**(2*k))) for k in range(1,33)),I(0))
gerr=abs(B[66]/F(66*nn**66));GAMMA+=I(-gerr.numerator,gerr.numerator)/gerr.denominator
LN4=2*I(2).ln()
c=-(GAMMA+(2*PI*L).ln())/2

# Active prime powers: every n=p^k with log n < 2L, weight Lambda(n)/sqrt(n)=log(p)/sqrt(n).
def prime_powers(limit):
    out=[]
    for n in range(2,limit+1):
        p=min(d for d in range(2,n+1) if n%d==0)
        m=n
        while m%p==0:m//=p
        if m==1:out.append((n,p))
    return out
TWO_L=2*I(L)
active=[];excluded_min=None
for n,p in prime_powers(64):
    ln_n=I(n).ln()
    if ln_n.hi<TWO_L.lo: active.append((n,p,ln_n))
    elif ln_n.lo>TWO_L.hi:
        excluded_min=n if excluded_min is None else excluded_min
    else: raise ValueError(f'log {n} not separated from 2L')
assert excluded_min is not None and I(64).ln().lo>TWO_L.hi
assert [n for n,_,_ in active]==list(config['active_prime_powers'])
PRIMES=[(n,ln_n,I(p).ln()/I(n).sqrt()) for n,p,ln_n in active]

# ---------------------------------------------------------------- trial
trial_path=here/'trial.json'
if MATRIX_MODE:
    M=int(config['matrix_modes']);NCOL=M
    coeff=[[F(1 if j==k else 0) for k in range(NCOL)] for j in range(M)]
else:
    trial=json.loads(trial_path.read_text());M,NCOL=trial['M'],trial['N']
    den=trial['denominator']
    coeff=[[F(x,den) for x in row] for row in trial['legendre_numerators']]
    assert len(coeff)==M and all(len(r)==NCOL for r in coeff)
    if RAYLEIGH_MODE:
        N=NCOL;HLEV=0
    else:
        assert NCOL==N
        assert all(coeff[j][k]==(1 if j==k else 0) for j in range(N) for k in range(N))
NROWS=max(M,N+HLEV)
Us=[legendre(2*j) for j in range(NROWS)]
Z=[];JZ=[]
for k in range(NCOL):
    z=[F(0)];jz=[F(0)]
    for j in range(M):
        if coeff[j][k]==0:continue
        z=add(z,scale(Us[j],coeff[j][k]))
        hj=sum((F(1,l) for l in range(1,2*j+1)),F(0))/2
        jz=add(jz,scale(Us[j],coeff[j][k]*hj))
    Z.append(z);JZ.append(jz)

# Exact rational Taylor coefficients of j(z)=e^(z/2)/(2 sinh z)-1/(2z).
ff=[]
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
for r in range(COSH_TERMS):
    n=2*r
    term=[F(comb(n,k))*L**n*(-F(1,2))**(n-k)/factorial(n) for k in range(n+1)]
    cosh_poly=add(cosh_poly,term)

P=[]
for k,z in enumerate(Z):
    p=add(JZ[k],kernel_poly(z))
    p=add([I(x) for x in p],scale(z,c))
    p=add(p,scale(cosh_poly,inner_poly(cosh_poly,z)))
    P.append(p)

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

# Segment structure in y=(x+L)/(2L): each active shift s'=log n/(2L) creates the
# edges s' (where z(y-s') switches on) and 1-s' (where z(y+s') switches off).
edge_list=[(I(0),'zero',None),(I(1),'one',None)]
for idx,(n,ln_n,a_n) in enumerate(PRIMES):
    s=ln_n/(2*L)
    edge_list.append((s,'on',idx));edge_list.append((1-s,'off',idx))
edge_list.sort(key=lambda e:e[0].lo)
assert all(x[0].hi<y[0].lo for x,y in zip(edge_list,edge_list[1:]))
edges=[e[0] for e in edge_list]
on_index={e[2]:i for i,e in enumerate(edge_list) if e[1]=='on'}
off_index={e[2]:i for i,e in enumerate(edge_list) if e[1]=='off'}
prims=[primitives(x,maxdegree) for x in edges]
NSEG=len(edges)-1

def interval_int(p,seg,weight=0):
    left=prims[seg][weight];right=prims[seg+1][weight]
    return 2*L*sum((I(z)*(right[k]-left[k]) for k,z in enumerate(p)),I(0))

pieces=[]
for seg in range(NSEG):
    row=[]
    for k,z in enumerate(Z):
        p=list(P[k])
        for idx,(n,ln_n,a_n) in enumerate(PRIMES):
            s=ln_n/(2*L)
            if on_index[idx]<=seg:p=add(p,scale(shift(z,-s),-a_n/2))
            if off_index[idx]>=seg+1:p=add(p,scale(shift(z,s),-a_n/2))
        row.append(p)
    pieces.append(row)

def pairing(poly,k):
    return weighted(mul(poly,Z[k]),1)+sum((interval_int(mul(poly,pieces[s][k]),s) for s in range(NSEG)),I(0))

A=[[pairing(Z[i],j) for j in range(NCOL)] for i in range(NCOL)]
project=[[pairing(Us[i],j) for j in range(NCOL)] for i in range(NROWS)]
full=[[I(0) for _ in range(NCOL)] for _ in range(NCOL)]
for i in range(NCOL):
    for j in range(i,NCOL):
        value=weighted(mul(Z[i],Z[j]),2)
        for s in range(NSEG):
            value+=interval_int(add(mul(Z[i],pieces[s][j]),mul(Z[j],pieces[s][i])),s,1)
            value+=interval_int(mul(pieces[s][i],pieces[s][j]),s)
        full[i][j]=full[j][i]=value
norms=[F(2*L,4*k+1) for k in range(NROWS)]

if MATRIX_MODE:
    out={'L':str(L),'modes':M,'active_prime_powers':[n for n,_,_ in PRIMES],
         'A':[[x.json() for x in row] for row in A],
         'project':[[x.json() for x in row] for row in project],
         'full_gram':[[x.json() for x in row] for row in full],
         'basis_norms':[str(x) for x in norms],'note':'exploration output, not a certificate'}
    Path(MATRIX_OUT).write_text(json.dumps(out)+'\n')
    print('matrices written',MATRIX_OUT);sys.exit(0)

if RAYLEIGH_MODE:
    # Rigorous two-sided enclosure of the Rayleigh quotient q[z]/||z||^2 of each
    # fixed column: q[z] = <z,T z> and |<z,(T-bar T)z>| <= epsilon ||z||^2.
    ej=I(5*L**(K+1)/(1-L))
    ec=I(2*(L/2)**(2*COSH_TERMS)/F(factorial(2*COSH_TERMS)))
    epsilon=L*ej+8*L*ec
    report=[]
    for k,z in enumerate(Z):
        n2=inner_poly(z,z);q=A[k][k]
        upper=(q+epsilon*n2).hi;lower=(q-epsilon*n2).lo
        report.append({'column':k,'norm_squared':n2.json(),'operator_value':q.json(),
            'q_upper':str(upper),'q_lower':str(lower),
            'rayleigh_upper':str(UP.divide(upper,n2.lo)),'rayleigh_lower':str(DN.divide(lower,n2.hi))})
    result={'status':'RAYLEIGH_ENCLOSURE','python':sys.version.split()[0],'decimal_precision':PREC,'L':str(L),
        'active_prime_powers':[n for n,_,_ in PRIMES],'trial_modes':M,'columns':NCOL,'kernel_degree':K,'cosh_terms':COSH_TERMS,
        'operator_error':epsilon.json(),'trial_sha256':hashlib.sha256(trial_path.read_bytes()).hexdigest(),
        'config_sha256':hashlib.sha256((here/'config.json').read_bytes()).hexdigest(),
        'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'columns_report':report,
        'claim':f'Two-sided enclosure of q_L[z]/||z||^2 for the fixed columns at L={L}; not a positivity or RH statement.'}
    out=here/f'rayleigh-{PREC}.json';out.write_text(json.dumps(result,indent=2)+'\n')
    for r in report:print(json.dumps({k:r[k] for k in ['column','rayleigh_lower','rayleigh_upper']}))
    sys.exit(0)

Gram=[[full[i][j]-sum((project[k][i]*project[k][j]/norms[k] for k in range(N)),I(0)) for j in range(N)] for i in range(N)]
tz=sum((inner_poly(z,z) for z in Z),I(0))
tr=sum((Gram[k][k] for k in range(N)),I(0))
assert tz.lo>0 and tr.lo>=0

# Complement comparison: C_N >= D on Q_N with D phi_(2j) = d_j phi_(2j),
# d_j = D0 + (H_(2j)-H_(2 J0))/2 for j >= J0 (strictly increasing in j).
# Either D0 is a declared proved constant (rule "declared"), or it is computed
# by the two-prime boundary-split rule (rule "two_prime_boundary_split"):
#   q[h] >= G_L N_h + J_L[h] + W_L[h] + l(h)^2 - a C_h(log2) - b C_h(log3)
#   for mean-zero h (T50 (203)-(204) with nu=L), |C_h(u)| <= N_h/2 on
#   disjoint caps (T56 (222)), prime 3 split between the norm and the
#   boundary energy W_L with total boundary spend at most 1, and
#   J_L[h] >= H_(2N) N_h/2 on Q_N (ASM (14)).
def atan_series(x):
    # 0<x<1 alternating series, remainder bounded by the first omitted term
    x=I(x);s=I(0);term=x;n=0
    while True:
        s=s+(term/(2*n+1) if n%2==0 else -(term/(2*n+1)))
        term=term*x*x;n+=1
        if term.hi<Decimal(10)**(-PREC-5):break
    return s+I(-term.hi,term.hi)
complement_rule=complement.get('rule','declared')
if complement_rule=='declared':
    D0=F(complement['base_value']);J0=int(complement['base_level'])
elif complement_rule=='boundary_split':
    # q[h] >= G_L N_h + J_L[h] + W_L[h] + l(h)^2 - sum_n a_n C_h(log n) for
    # mean-zero h (T50 (203)-(204), nu=L).  For a shift u with L<u<2L the two
    # overlap caps are disjoint and |C_h(u)| <= N_h/2 <= W_L-cap integral/(2 w_min)
    # (T56 (222)); otherwise |C_h(u)| <= N_h.  The boundary energy W_L (coefficient
    # one) is spent greedily on the caps with the largest w_min.  Then
    # J_L[h] >= H_(2N) N_h/2 on Q_N (ASM (14)).
    x=(-I(L)).exp()
    F2L=((1+x)/(1-x)).ln()/2+atan_series(x)
    OMEGA0=-GAMMA-PI/2-3*I(2).ln()-PI.ln()
    nu=I(L)
    G_L=OMEGA0/2+F2L+I(2).ln()/2+I(L)/2-nu*L/48-I(L)**3/12-I(L)**3/(36*nu)
    G_L=G_L+((2*I(L)-nu)**3-(I(L)-nu)**3)/(576*nu)
    cap_info=[]
    for n,ln_n,a_n in PRIMES:
        if ln_n.lo>I(L).hi:
            assert ln_n.hi<TWO_L.lo
            wmin=(I(L)*L/(ln_n*(2*I(L)-ln_n))).ln()/4
            assert wmin.lo>0
            cap_info.append((n,a_n,wmin,True))
        else:
            cap_info.append((n,a_n,None,False))
    order=sorted([ci for ci in cap_info if ci[3]],key=lambda ci:-ci[2].lo)
    budget=F(1);charges={};spend={}
    for n,a_n,wmin,_ in order:
        rho=a_n/(2*wmin)          # boundary coefficient if fully charged to W_L
        if F(rho.hi)<=budget: f=F(1);used=F(rho.hi)
        elif budget<=0: f=F(0);used=F(0)
        else: f=budget/F(rho.hi);used=budget   # f*rho <= f*rho.hi = budget exactly
        spend[n]=used
        budget=budget-used
        charges[n]=((1-f)*a_n/2).hi
    for n,a_n,wmin,disjoint in cap_info:
        if not disjoint: charges[n]=a_n.hi
    assert sum(spend.values(),F(0))<=1 and budget>=0
    J0=N
    D0=F(G_L.lo)+hm(2*N)/2-sum((F(v) for v in charges.values()),F(0))
    complement_report={'G_L':G_L.json(),'F_2L':F2L.json(),
        'caps':[{'n':n,'disjoint_caps':disjoint,'w_min':(wmin.json() if disjoint else None),
                 'norm_charge':str(charges[n]),'boundary_spend':str(spend.get(n,0))} for n,a_n,wmin,disjoint in cap_info],
        'd_N':str(D0)}
else: raise ValueError('unknown complement rule')
def dlevel(j):
    assert j>=J0
    return D0+(hm(2*j)-hm(2*J0))/2
assert dlevel(N)>0
kappa_max=1/dlevel(N)
levels=[(k,dlevel(k)) for k in range(N,N+HLEV)]
tail_weight=1/dlevel(N+HLEV)
Ehier=[[I(0) for _ in range(N)] for _ in range(N)]
for i in range(N):
    for j in range(N):
        explicit=sum((project[k][i]*project[k][j]/norms[k] for k,_ in levels),I(0))
        value=sum((project[k][i]*project[k][j]/(norms[k]*d) for k,d in levels),I(0))
        value+=tail_weight*(Gram[i][j]-explicit)
        Ehier[i][j]=value

# Analytic operator approximation: j by its degree-K Taylor polynomial on
# 0<=s<=2L (|j|<5 on |z|=2, Cauchy, geometric tail ratio L), cosh(x/2) by
# COSH_TERMS even terms (positive tail, successive ratio below 1/2).
ej=I(5*L**(K+1)/(1-L))
ec=I(2*(L/2)**(2*COSH_TERMS)/F(factorial(2*COSH_TERMS)))
epsilon=L*ej+8*L*ec
delta=epsilon*tz+kappa_max*(2*epsilon*(tz*tr).sqrt()+epsilon*epsilon*tz)
delta=I(delta.hi)  # one frozen outward rational allowance
lower=[[((A[i][j]+A[j][i])/2-Ehier[i][j]-(delta if i==j else 0)) for j in range(N)] for i in range(N)]

def ldlt(mat,shift=F(0)):
    n=len(mat);ll=[[I(0) for _ in range(n)] for _ in range(n)];dd=[]
    for i in range(n):
        d=mat[i][i]-shift-sum((ll[i][k]*ll[i][k]*dd[k] for k in range(i)),I(0))
        if d.lo<=0:return False,dd+[d]
        dd.append(d);ll[i][i]=I(1)
        for j in range(i+1,n):ll[j][i]=(mat[j][i]-sum((ll[j][k]*ll[i][k]*dd[k] for k in range(i)),I(0)))/d
    return True,dd

passed,pivots=ldlt(lower)
safe_shift=MARGIN
shift_pass,shift_pivots=ldlt(lower,safe_shift)

# Second sign check: exact Fraction arithmetic on a rational matrix dominating
# every matrix in the directed entry enclosure.
grid=10**int(config['rational_grid_digits'])
rat=[[F(0) for _ in range(N)] for _ in range(N)]
for i in range(N):
    for j in range(i,N):
        midpoint=(F(lower[i][j].lo)+F(lower[i][j].hi))/2
        z=(midpoint*grid).__floor__()
        rat[i][j]=rat[j][i]=F(z,grid)
row_errors=[]
for i in range(N):
    row_errors.append(sum((max(abs(F(lower[i][j].lo)-rat[i][j]),abs(F(lower[i][j].hi)-rat[i][j])) for j in range(N)),F(0)))
round_allow=F(config['rational_rounding_allowance'])
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

# All-function transfer. With E=Z-C^(-1)R and C^(-1) <= D^(-1) in form order,
# ||C^(-1)R|| <= ||C^(-1/2)|| ||C^(-1/2)R|| <= sqrt(kappa_max) sqrt(||R^*C^(-1)R||)
# and ||R^*C^(-1)R|| <= tr E_hier(R) <= tr E_hier(bar R) + kappa_max(2 eps sqrt(t_z t_r)+eps^2 t_z).
# Then q[g] >= MARGIN|alpha|^2 + d_N ||z||^2 and ||g||^2 <= 2||E||^2|alpha|^2+2||z||^2.
trE=sum((Ehier[k][k] for k in range(N)),I(0))
inverse_residual_bound=trE+kappa_max*(2*epsilon*(tz*tr).sqrt()+epsilon*epsilon*tz)
Ebound=(tz.sqrt()+I(kappa_max).sqrt()*inverse_residual_bound.sqrt()).hi
Ebound=F(Ebound)
family_bound=min(MARGIN/(2*Ebound*Ebound),dlevel(N)/2)
assert family_bound>0
declared_family=F(config['declared_family_bound'])
assert declared_family<=family_bound

controls={}
controls['positive_diagonal_indefinite_matrix_rejected']=not exact_ldlt([[F(1),F(2)],[F(2),F(1)]])[0]
controls['visible_only_false_positive_rejected']=not exact_ldlt([[F(1)-F(2)**2]])[0]
controls['zero_pivot_rejected']=not exact_ldlt([[F(0)]])[0]
bad=[row[:] for row in rat_test];bad[0][0]=-1
controls['negative_direction_rejected']=not exact_ldlt(bad)[0]
controls['inflated_uncertainty_rejected']=not ldlt(lower,F(1,100))[0]
try:I(1)/I(-1,1)
except ZeroDivisionError:controls['singular_interval_inverse_rejected']=True
else:controls['singular_interval_inverse_rejected']=False
try:I(0.1)
except TypeError:controls['binary_float_input_rejected']=True
else:controls['binary_float_input_rejected']=False
# Dropping the complete residual debt must not be what makes the sign pass:
controls['residual_debt_is_nonzero']=tr.lo>0
# Lowering every complement level to its tail value only weakens the bound:
weak=[[((A[i][j]+A[j][i])/2-kappa_max*Gram[i][j]-(delta if i==j else 0)) for j in range(N)] for i in range(N)]
controls['uniform_allowance_run_recorded']=ldlt(weak,safe_shift)[0]
assert all(v for k,v in controls.items() if k!='uniform_allowance_run_recorded')

result={'status':'PASS_FIXED_WINDOW_LOWER_MATRIX' if passed else 'INCONCLUSIVE',
    'python':sys.version.split()[0],'decimal_precision':PREC,'L':str(L),'active_prime_powers':[n for n,_,_ in PRIMES],
    'retained_modes':N,'trial_modes':M,'kernel_degree':K,'cosh_terms':COSH_TERMS,
    'complement_rule':complement_rule,'complement_base_level':J0,'complement_base_value':str(D0),'complement_report':complement_report if complement_rule!='declared' else None,'complement_levels':[[k,str(d)] for k,d in levels],'complement_tail_level':N+HLEV,'complement_tail_value':str(dlevel(N+HLEV)),
    'kappa_max':str(kappa_max),'trial_sha256':hashlib.sha256(trial_path.read_bytes()).hexdigest(),
    'config_sha256':hashlib.sha256((here/'config.json').read_bytes()).hexdigest(),
    'operator_error':epsilon.json(),'lower_scalar_allowance':delta.json(),'trial_norm_trace':tz.json(),
    'residual_gram_trace':tr.json(),'lower_matrix':[[x.json() for x in row] for row in lower],
    'ldlt_pivots':[x.json() for x in pivots],
    'shift':str(safe_shift),'shift_pass':shift_pass,'shifted_pivots':[x.json() for x in shift_pivots],
    'rational_grid':grid,'rational_matrix_numerators':[[int(x*grid) for x in row] for row in rat],
    'rational_rounding_allowance':str(round_allow),'rational_row_error_upper':[str(x) for x in row_errors],
    'exact_rational_shift_pass':exact_pass,'exact_shifted_pivots':[str(x) for x in exact_pivots],
    'lift_norm_bound':str(Ebound),'inverse_residual_trace_bound':inverse_residual_bound.json(),'full_family_window':f'0<L<={L}','full_family_lower_bound':str(declared_family),
    'full_family_lower_bound_computed':str(family_bound),
    'negative_controls':controls,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'claim':f'Fixed arithmetic source at L={L} only; not all-width positivity or RH.'}
out=here/f'result-{PREC}.json';out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['status','decimal_precision','operator_error','lower_scalar_allowance','shift_pass','full_family_lower_bound']}),flush=True)
print(json.dumps({'exact_rational_shift_pass':exact_pass,'negative_controls':controls}),flush=True)
