"""Exact finite controls for the written correction proofs. No physical verdict."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib,json,sys
HERE=Path(__file__).resolve().parent

def spectral_upper(mu, tail):
    if tail < 0: raise ValueError('tail must be a certified nonnegative upper bound')
    return (max(mu,F(0)) if mu is not None else F(0))+tail

def pluecker(c):
    if len(c)!=6: raise ValueError('six independent skew coordinates required')
    return c[0]*c[5]-c[1]*c[4]+c[2]*c[3]

def wedge(u,v): return tuple(u[i]*v[j]-v[i]*u[j] for i in range(4) for j in range(i+1,4))

def pluecker_upper(coordinate_norm_upper, eta):
    if coordinate_norm_upper<0 or eta<0: raise ValueError('nonnegative bounds required')
    return coordinate_norm_upper*eta+F(3,2)*eta**2

def build():
    checks={}
    checks['old_spectral_bound_false_accepts'] = -2+2<=1 and 2>1
    checks['corrected_spectral_bound_rejects']=spectral_upper(F(-2),F(2))==2>1
    checks['zero_rank_projection']=spectral_upper(None,F(2))==2
    checks['full_projection_negative_is_safe']=spectral_upper(F(-2),F(0))==0>=-2
    checks['threshold_equality_allowed']=spectral_upper(F(1),F(0))==1
    checks['diagonal_exhaustion']=all(max(a,d)<=spectral_upper(F(a),F(abs(d)))
        for a,d in product(range(-5,6),repeat=2))
    # K=[[-2,1],[1,2]], P=diag(1,0): tail Frobenius square 6 <= (5/2)^2.
    # The upper-minus-K has nonnegative principal minors, checked exactly.
    u=spectral_upper(F(-2),F(5,2))
    checks['coupled_complement_bound']=F(6)<=F(25,4) and u+2>=0 and u-2>=0 and (u+2)*(u-2)-1>=0
    checks['wedge_identity_finite_grid']=all(pluecker(wedge(x[:4],x[4:]))==0
        for x in product((F(-1),F(0),F(1)),repeat=8))
    c=(F(1),F(0),F(0),F(0),F(0),F(1))
    checks['independent_bracket_can_violate']=pluecker(c)==1
    checks['nonzero_error_negative_control']=sum(x*x for x in c)<=F(3,2)**2 and pluecker(c)>pluecker_upper(F(3,2),F(1,10))
    # Decomposable B with B12=1, perturb only C34 by 1/10; known error=1/10.
    c2=(F(1),F(0),F(0),F(0),F(0),F(1,10))
    checks['perturbed_independent_positive_control']=sum(x*x for x in c2)<=F(11,10)**2 and abs(pluecker(c2))<=pluecker_upper(F(11,10),F(1,10))
    return {'schema':'native-prerequisite-corrections-v1','status':'PASS_FINITE_CONTROLS' if all(checks.values()) else 'FAIL',
      'scope':'Exact declared finite controls only. General proofs are in README.md and corrected sources. No Lean, peer-review, empirical, RH, YM or quantum-gravity promotion.',
      'checks':checks}

def canonical(payload): return (json.dumps(payload,indent=2,sort_keys=True)+'\n').encode()
if __name__=='__main__':
    body=canonical(build());digest=hashlib.sha256(body).hexdigest()
    if '--check' in sys.argv:
        assert body==(HERE/'CERTIFICATE.json').read_bytes()
        assert digest==(HERE/'EXPECTED.sha256').read_text().strip()
    else:
        (HERE/'CERTIFICATE.json').write_bytes(body);(HERE/'EXPECTED.sha256').write_text(digest+'\n')
    print(build()['status'],digest)
    raise SystemExit(0 if all(build()['checks'].values()) else 1)
