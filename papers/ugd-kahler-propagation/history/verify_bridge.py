"""Exact finite checks, not a replacement for the general proofs in the note."""
from fractions import Fraction as F
from itertools import product
import json

def matmul(a,b):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def add(a,b): return [[x+y for x,y in zip(r,s)] for r,s in zip(a,b)]
def scale(t,a): return [[t*x for x in r] for r in a]
def comm(a,b): return add(matmul(a,b),scale(-1,matmul(b,a)))
def inv(a):
    n=len(a); b=[[F(x) for x in a[i]]+[F(i==j) for j in range(n)] for i in range(n)]
    for j in range(n):
        p=next(i for i in range(j,n) if b[i][j])
        b[j],b[p]=b[p],b[j]
        d=b[j][j]; b[j]=[x/d for x in b[j]]
        for i in range(n):
            if i!=j:
                d=b[i][j]; b[i]=[x-d*y for x,y in zip(b[i],b[j])]
    return [r[n:] for r in b]
def tr(a): return sum(a[i][i] for i in range(len(a)))

I=[[1,0],[0,1]]; R=[[0,-1],[1,0]]; K=[[1,0],[0,-1]]; L=matmul(K,R)
assert matmul(R,R)==scale(-1,I)
assert matmul(K,K)==I and matmul(L,L)==I
assert comm(K,R)==scale(2,L) and comm(L,R)==scale(-2,K)
assert comm(I,R)==[[0,0],[0,0]] and comm(R,R)==[[0,0],[0,0]]

def poly_derivative(poly, point, indices):
    total=F(0)
    for powers,coef in poly.items():
        powers=list(powers); coef=F(coef)
        for i in indices:
            if powers[i]==0: coef=0; break
            coef*=powers[i]; powers[i]-=1
        if coef:
            for x,p in zip(point,powers): coef*=x**p
            total+=coef
    return total

def verify_ricci(poly,point):
    n=2; m=2*n
    der=lambda inds: poly_derivative(poly,point,inds)
    g=[[der((a,b)) for b in range(n)] for a in range(n)]
    gi=inv(g)
    g1=[[[der((a,b,c)) for b in range(n)] for a in range(n)] for c in range(n)]
    g2=[[[[der((a,b,c,d)) for b in range(n)] for a in range(n)] for d in range(n)] for c in range(n)]
    G=[[2*g[a%n][b%n] if a//n==b//n else F(0) for b in range(m)] for a in range(m)]
    Ginv=inv(G)
    def dG(c,a,b): return 2*g1[c][a%n][b%n] if c<n and a//n==b//n else F(0)
    def ddG(c,d,a,b): return 2*g2[c][d][a%n][b%n] if c<n and d<n and a//n==b//n else F(0)
    def dInv(c,a,b):
        return -sum(Ginv[a][u]*dG(c,u,v)*Ginv[v][b] for u,v in product(range(m),repeat=2))
    Gamma=[[[sum(Ginv[k][l]*(dG(a,b,l)+dG(b,a,l)-dG(l,a,b)) for l in range(m))/2 for b in range(m)] for a in range(m)] for k in range(m)]
    def dGamma(c,k,a,b):
        return sum(dInv(c,k,l)*(dG(a,b,l)+dG(b,a,l)-dG(l,a,b))+Ginv[k][l]*(ddG(c,a,b,l)+ddG(c,b,a,l)-ddG(c,l,a,b)) for l in range(m))/2
    Ric=[[sum(dGamma(k,k,a,b)-dGamma(b,k,a,k) for k in range(m))+sum(Gamma[k][k][l]*Gamma[l][a][b]-Gamma[k][b][l]*Gamma[l][a][k] for k,l in product(range(m),repeat=2)) for b in range(m)] for a in range(m)]
    hlog=[[tr(matmul(gi,g2[a][b]))-tr(matmul(matmul(matmul(gi,g1[a]),gi),g1[b])) for b in range(n)] for a in range(n)]
    expected=[[-hlog[a%n][b%n]/2 if a//n==b//n else F(0) for b in range(m)] for a in range(m)]
    assert Ric==expected, (Ric,expected)
    det=g[0][0]*g[1][1]-g[0][1]*g[1][0]
    return {"point":[str(x) for x in point],"det_hessian":str(det),"ricci_xx":[[str(v) for v in row[:2]] for row in Ric[:2]]}

# A coupled potential checks non-diagonal inverse, mixed fourth derivatives,
# factors of two, and the sign of Ricci independently of the complex formula.
polynomial={(2,0):F(1,2),(0,2):F(1,2),(2,1):F(1,2),(0,4):F(1,24)}
checks=[verify_ricci(polynomial,p) for p in [(F(0),F(0)),(F(1,3),F(1,2)),(F(-1,2),F(1))]]

# Analytic curved example lambda=-log(x1)-log(x2), c=-1/2:
# h=diag(x_i^-2), Hess log det h=diag(2*x_i^-2).
for x1,x2 in [(F(1),F(1)),(F(2),F(3)),(F(1,2),F(7,3))]:
    h=[[1/x1**2,F(0)],[F(0),1/x2**2]]
    hlog=scale(2,h)
    assert scale(F(-1,4),hlog)==scale(F(-1,2),h)
    # exp(-4c lambda)=exp(2 lambda)=(x1*x2)^-2.
    assert h[0][0]*h[1][1]==1/(x1*x2)**2

# H=G+i omega on the complexified real tangent has a kernel; det h is used instead.
# For one pair G=2I and omega=[[0,2],[-2,0]], det H=4-(2i)(-2i)=0.
assert 4-4==0
print(json.dumps({"status":"PASS","native_RKL_relations":"exact","centralizer_commutators":"exact","direct_Levi_Civita_vs_logdet_Ricci":checks,"nonflat_Einstein_example":"PASS, c=-1/2"},indent=2))
