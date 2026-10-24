"""Ball-arithmetic engine for the fixed polynomial RH source trials.

Second implementation of verify.py's certificate with the same inputs
(config.json, trial.json), the same analytic pipeline and the same
decisions, but with every enclosure carried by FLINT/Arb balls
(python-flint) instead of directed Decimal endpoints, and with
polynomial products done in Arb's compiled polynomial arithmetic.
Every ball is rigorous: a midpoint with a radius that is rounded
outward by the library, and the lower/upper bounds extracted for the
exact rational sign check are exact dyadic numbers below/above the
whole ball. Constants (pi, Euler's constant, logarithms, exp, sqrt,
arctan) are Arb's certified enclosures. The exact rational LDL^T, the
rational rounding allowance and every control are unchanged.

Purpose: speed. The Decimal verifier needs hours for the widths beyond
3/4; this engine needs minutes. Both engines are kept: verify.py is the
standard-library reference, verify_arb.py the fast implementation; a
certificate run with either is a run of the same declared certificate.

Usage: python3 verify_arb.py [precision_digits] --packet DIR [--rayleigh]
       python3 verify_arb.py [precision_digits] --packet DIR --matrices OUT.json
Writes result-arb-<digits>.json (or rayleigh-arb-<digits>.json).
Requires python-flint (pip install python-flint).
"""
from decimal import Decimal, Context, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction as F
from math import comb, factorial, ceil, log2
from pathlib import Path
import json, sys, hashlib
from flint import arb, arb_poly, fmpq, ctx

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
BITS=int(ceil(PREC*log2(10)))+64
ctx.prec=BITS
DN=Context(prec=PREC,rounding=ROUND_FLOOR)
UP=Context(prec=PREC,rounding=ROUND_CEILING)

def _exact(a):
    # exact dyadic value of a zero-radius arb as a Fraction
    assert a.is_exact()
    m,e=a.man_exp();m=int(m);e=int(e)
    return F(m)*(F(2)**e) if e>=0 else F(m,2**(-e))
def _dec(f,ctxt):
    return ctxt.divide(Decimal(f.numerator),Decimal(f.denominator))

class I:
    """Rigorous ball with the interface of verify.py's interval class."""
    __slots__=('b',)
    def __init__(self,x=0,hi=None):
        if isinstance(x,I): self.b=x.b; return
        if isinstance(x,arb): self.b=x; return
        if isinstance(x,float): raise TypeError('float forbidden in certificate')
        if isinstance(x,F): self.b=arb(fmpq(x.numerator,x.denominator))
        elif isinstance(x,int): self.b=arb(x)
        elif isinstance(x,(str,Decimal)): self.b=arb(fmpq(F(str(x)).numerator,F(str(x)).denominator))
        else: raise TypeError(type(x))
        if hi is not None:
            h=I(hi).b
            self.b=self.b.union(h)
    @property
    def lo(self): return _exact(self.b.lower())
    @property
    def hi(self): return _exact(self.b.upper())
    def __add__(self,o): return I(self.b+I(o).b)
    __radd__=__add__
    def __neg__(self): return I(-self.b)
    def __sub__(self,o): return I(self.b-I(o).b)
    def __rsub__(self,o): return I(I(o).b-self.b)
    def __mul__(self,o): return I(self.b*I(o).b)
    __rmul__=__mul__
    def __truediv__(self,o):
        d=I(o).b
        if d.contains(0): raise ZeroDivisionError('interval contains zero')
        return I(self.b/d)
    def __rtruediv__(self,o): return I(o)/self
    def __pow__(self,n):
        # binary powering by multiplication: ball products stay finite when
        # the ball contains zero (arb's real power would not)
        assert isinstance(n,int) and n>=0
        ans=arb(1);base=self.b
        while n:
            if n&1:ans=ans*base
            base=base*base;n//=2
        return I(ans)
    def ln(self):
        if not (self.b>0): raise ValueError('log domain')
        return I(self.b.log())
    def exp(self): return I(self.b.exp())
    def sqrt(self):
        if self.b<0: raise ValueError('sqrt domain')
        return I(self.b.sqrt())
    def json(self): return [str(_dec(self.lo,DN)),str(_dec(self.hi,UP))]
    def __repr__(self): return 'I('+self.b.str(20,radius=True)+')'

def A_(x): return I(x).b           # arb from Fraction/int/I
def poly(coeffs): return arb_poly([A_(c) for c in coeffs])
def pcoeffs(p): return p.coeffs()
def pshift(p,a):
    # p(y+a) by Horner in the polynomial ring
    cs=p.coeffs();lin=arb_poly([A_(a),1]);out=arb_poly([])
    for c in reversed(cs): out=out*lin+arb_poly([c])
    return out
def pscale(p,a): return p*arb_poly([A_(a)])
def legendre(n): return [F((-1)**(n+k)*comb(n,k)*comb(n+k,k)) for k in range(n+1)]
def hm(n,p=1): return sum((F(1,k**p) for k in range(1,n+1)),F(0))
def dot(cs,vec):
    s=arb(0)
    for k,c in enumerate(cs): s+=c*vec[k]
    return s

# ---------------------------------------------------------------- config
here=PACKET
config=json.loads((here/'config.json').read_text())
L=F(config['L']); assert 0<L<1
N=int(config['retained_modes']); K=int(config['kernel_degree']); COSH_TERMS=int(config['cosh_terms'])
complement=config['complement']
HLEV=int(complement['explicit_levels']); assert HLEV>=0
MARGIN=F(config['margin']); assert MARGIN>0
Lb=A_(L); TWO_Lb=A_(2*L)

def inner_poly(p,q):
    cs=(p*q).coeffs()
    return I(TWO_Lb*sum((c/(k+1) for k,c in enumerate(cs)),arb(0)))

PI=I(arb.pi()); GAMMA=I(arb.const_euler())
LN4=2*I(2).ln()
c=-(GAMMA+(2*PI*L).ln())/2

def prime_powers(limit):
    out=[]
    for n in range(2,limit+1):
        p=min(d for d in range(2,n+1) if n%d==0)
        m=n
        while m%p==0:m//=p
        if m==1:out.append((n,p))
    return out
TWO_L=I(2*L)
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
Us_q=[legendre(2*j) for j in range(NROWS)]
Us=[poly(u) for u in Us_q]
Z=[];JZ=[]
for k in range(NCOL):
    zq=[F(0)];jzq=[F(0)]
    for j in range(M):
        if coeff[j][k]==0:continue
        hj=hm(2*j)/2
        u=Us_q[j]
        zq=zq+[F(0)]*max(0,len(u)-len(zq)); jzq=jzq+[F(0)]*max(0,len(u)-len(jzq))
        for t,v in enumerate(u): zq[t]+=coeff[j][k]*v; jzq[t]+=coeff[j][k]*hj*v
    Z.append(poly(zq));JZ.append(poly(jzq))
Zq=[]   # exact rational columns, used for the kernel column
for k in range(NCOL):
    zq=[F(0)]
    for j in range(M):
        if coeff[j][k]==0:continue
        u=Us_q[j]; zq=zq+[F(0)]*max(0,len(u)-len(zq))
        for t,v in enumerate(u): zq[t]+=coeff[j][k]*v
    Zq.append(zq)

# Exact rational Taylor coefficients of j(z)=e^(z/2)/(2 sinh z)-1/(2z).
ff=[]
for n in range(K+2):
    f=F(1,2**n*factorial(n))
    for q in range(1,n//2+1):f-=F(2,factorial(2*q+1))*ff[n-2*q]
    ff.append(f/2)
jcoef=ff[1:]
assert jcoef[:5]==[F(1,4),-F(1,48),-F(1,32),F(7,11520),F(5,1536)]
jcoef_b=[A_(x) for x in jcoef]
TWOL_pow=[arb(1)]
for m in range(1,K+1): TWOL_pow.append(TWOL_pow[-1]*TWO_Lb)

def kernel_poly(pq):
    # pq: exact rational coefficient list; returns arb list of -(1/2)K_j applied to it
    out=[arb(0) for _ in range(len(pq)+K+1)]
    pb=[A_(x) for x in pq]
    for m in range(K+1):
        fac=-Lb*jcoef_b[m]*TWOL_pow[m]
        binom=[comb(m,r) for r in range(m+1)]
        for j,pj in enumerate(pq):
            if not pj:continue
            f=fac*pb[j]
            for r in range(m+1):
                out[r]+=f*((-binom[r] if r%2 else binom[r])/arb(m-r+j+1))
            if m%2:out[m+j+1]+=f*arb(fmpq(2*factorial(m)*factorial(j),factorial(m+j+1)))
    return out

cosh_q=[F(0)]
for r in range(COSH_TERMS):
    n=2*r
    term=[F(comb(n,k))*L**n*(-F(1,2))**(n-k)/factorial(n) for k in range(n+1)]
    cosh_q=cosh_q+[F(0)]*max(0,len(term)-len(cosh_q))
    for t,v in enumerate(term): cosh_q[t]+=v
cosh_poly=poly(cosh_q)

P=[]
for k in range(NCOL):
    kp=arb_poly(kernel_poly(Zq[k]))
    p=JZ[k]+kp+pscale(Z[k],c)+pscale(cosh_poly,inner_poly(cosh_poly,Z[k]))
    P.append(p)

maxdegree=2*max(p.degree()+1 for p in P)+2
Wmom=[];W2mom=[]
PIb=PI.b;LN4b=LN4.b
for r in range(maxdegree+1):
    n=r+1;h=hm(n);h2=hm(n,2)
    sy=-A_(F(1,n*n));st=-A_(h/n)
    yy=A_(F(2,n**3));tt=A_((h*h+h2)/n)
    yt=A_(h/F(n*n))+A_(h2/n)-PIb*PIb/(6*n)
    Wmom.append(-(LN4b/n+sy+st)/4)
    W2mom.append((LN4b*LN4b/n+2*LN4b*(sy+st)+yy+tt+2*yt)/16)

def weighted(p,kind):
    moments=Wmom if kind==1 else W2mom
    return I(TWO_Lb*dot(p.coeffs(),moments))

def primitives(x,d):
    xb=x.b
    powers=[arb(1)]
    for _ in range(d+1):powers.append(powers[-1]*xb)
    ordinary=[powers[j+1]/(j+1) for j in range(d+1)]
    iszero=x.lo==x.hi==0;isone=x.lo==x.hi==1
    lx=arb(0) if iszero else xb.log()
    l1=arb(0) if isone else (1-xb).log()
    accum=arb(0);w=[]
    for j in range(d+1):
        n=j+1;accum+=powers[n]/n
        logy=arb(0) if iszero else powers[n]*(lx/n-arb(fmpq(1,n*n)))
        log1=-A_(hm(n)/n) if isone else ((powers[n]-1)*l1-accum)/n
        w.append(-(LN4b*ordinary[j]+logy+log1)/4)
    return ordinary,w

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
# differences of primitives per segment and weight, as plain lists
diffs=[[[prims[s+1][w][k]-prims[s][w][k] for k in range(maxdegree+1)] for w in (0,1)] for s in range(NSEG)]

def interval_int(p,seg,weight=0):
    return I(TWO_Lb*dot(p.coeffs(),diffs[seg][weight]))

pieces=[]
for seg in range(NSEG):
    row=[]
    for k in range(NCOL):
        p=P[k]
        for idx,(n,ln_n,a_n) in enumerate(PRIMES):
            s=ln_n/(2*L)
            if on_index[idx]<=seg:p=p+pscale(pshift(Z[k],-s),-a_n/2)
            if off_index[idx]>=seg+1:p=p+pscale(pshift(Z[k],s),-a_n/2)
        row.append(p)
    pieces.append(row)

def pairing(polyobj,k):
    return weighted(polyobj*Z[k],1)+sum((interval_int(polyobj*pieces[s][k],s) for s in range(NSEG)),I(0))

A=[[pairing(Z[i],j) for j in range(NCOL)] for i in range(NCOL)]
project=[[pairing(Us[i],j) for j in range(NCOL)] for i in range(NROWS)]
full=[[I(0) for _ in range(NCOL)] for _ in range(NCOL)]
for i in range(NCOL):
    for j in range(i,NCOL):
        value=weighted(Z[i]*Z[j],2)
        for s in range(NSEG):
            value+=interval_int(Z[i]*pieces[s][j]+Z[j]*pieces[s][i],s,1)
            value+=interval_int(pieces[s][i]*pieces[s][j],s)
        full[i][j]=full[j][i]=value
norms=[F(2*L,4*k+1) for k in range(NROWS)]

if MATRIX_MODE:
    out={'L':str(L),'modes':M,'active_prime_powers':[n for n,_,_ in PRIMES],
         'A':[[x.json() for x in row] for row in A],
         'project':[[x.json() for x in row] for row in project],
         'full_gram':[[x.json() for x in row] for row in full],
         'basis_norms':[str(x) for x in norms],'engine':'arb','note':'exploration output, not a certificate'}
    Path(MATRIX_OUT).write_text(json.dumps(out)+'\n')
    print('matrices written',MATRIX_OUT);sys.exit(0)

ej=I(5*L**(K+1)/(1-L))
ec=I(2*(L/2)**(2*COSH_TERMS)/F(factorial(2*COSH_TERMS)))
epsilon=L*ej+8*L*ec

if RAYLEIGH_MODE:
    report=[]
    for k in range(NCOL):
        n2=inner_poly(Z[k],Z[k]);q=A[k][k]
        upper=(q+epsilon*n2).hi;lower=(q-epsilon*n2).lo
        report.append({'column':k,'norm_squared':n2.json(),'operator_value':q.json(),
            'q_upper':str(_dec(upper,UP)),'q_lower':str(_dec(lower,DN)),
            'rayleigh_upper':str(_dec(upper/n2.lo,UP)),'rayleigh_lower':str(_dec(lower/n2.hi,DN))})
    result={'status':'RAYLEIGH_ENCLOSURE','engine':'arb','arb_precision_bits':BITS,'python':sys.version.split()[0],'decimal_precision':PREC,'L':str(L),
        'active_prime_powers':[n for n,_,_ in PRIMES],'trial_modes':M,'columns':NCOL,'kernel_degree':K,'cosh_terms':COSH_TERMS,
        'operator_error':epsilon.json(),'trial_sha256':hashlib.sha256(trial_path.read_bytes()).hexdigest(),
        'config_sha256':hashlib.sha256((here/'config.json').read_bytes()).hexdigest(),
        'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'columns_report':report,
        'claim':f'Two-sided enclosure of q_L[z]/||z||^2 for the fixed columns at L={L}; not a positivity or RH statement.'}
    out=here/f'rayleigh-arb-{PREC}.json';out.write_text(json.dumps(result,indent=2)+'\n')
    for r in report:print(json.dumps({k:r[k] for k in ['column','rayleigh_lower','rayleigh_upper']}))
    sys.exit(0)

Gram=[[full[i][j]-sum((project[k][i]*project[k][j]/norms[k] for k in range(N)),I(0)) for j in range(N)] for i in range(N)]
tz=sum((inner_poly(z,z) for z in Z),I(0))
tr=sum((Gram[k][k] for k in range(N)),I(0))
assert tz.lo>0 and tr.lo>=0

complement_rule=complement.get('rule','declared')
if complement_rule=='declared':
    D0=F(complement['base_value']);J0=int(complement['base_level'])
elif complement_rule=='boundary_split':
    x=(-I(L)).exp()
    F2L=((1+x)/(1-x)).ln()/2+I(x.b.atan())
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
        rho=a_n/(2*wmin)
        if F(rho.hi)<=budget: f=F(1);used=F(rho.hi)
        elif budget<=0: f=F(0);used=F(0)
        else: f=budget/F(rho.hi);used=budget
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

delta=epsilon*tz+kappa_max*(2*epsilon*(tz*tr).sqrt()+epsilon*epsilon*tz)
delta=I(delta.hi)
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

grid=10**int(config['rational_grid_digits'])
rat=[[F(0) for _ in range(N)] for _ in range(N)]
lo_hi=[[(lower[i][j].lo,lower[i][j].hi) for j in range(N)] for i in range(N)]
for i in range(N):
    for j in range(i,N):
        midpoint=(lo_hi[i][j][0]+lo_hi[i][j][1])/2
        z=(midpoint*grid).__floor__()
        rat[i][j]=rat[j][i]=F(z,grid)
row_errors=[]
for i in range(N):
    row_errors.append(sum((max(abs(lo_hi[i][j][0]-rat[i][j]),abs(lo_hi[i][j][1]-rat[i][j])) for j in range(N)),F(0)))
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

trE=sum((Ehier[k][k] for k in range(N)),I(0))
inverse_residual_bound=trE+kappa_max*(2*epsilon*(tz*tr).sqrt()+epsilon*epsilon*tz)
Ebound=(tz.sqrt()+I(kappa_max).sqrt()*inverse_residual_bound.sqrt()).hi
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
controls['residual_debt_is_nonzero']=tr.lo>0
weak=[[((A[i][j]+A[j][i])/2-kappa_max*Gram[i][j]-(delta if i==j else 0)) for j in range(N)] for i in range(N)]
controls['uniform_allowance_run_recorded']=ldlt(weak,safe_shift)[0]
assert all(v for k,v in controls.items() if k!='uniform_allowance_run_recorded')

result={'status':'PASS_FIXED_WINDOW_LOWER_MATRIX' if passed else 'INCONCLUSIVE','engine':'arb','arb_precision_bits':BITS,
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
out=here/f'result-arb-{PREC}.json';out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['status','engine','decimal_precision','operator_error','lower_scalar_allowance','shift_pass','full_family_lower_bound']}),flush=True)
print(json.dumps({'exact_rational_shift_pass':exact_pass,'negative_controls':controls}),flush=True)
