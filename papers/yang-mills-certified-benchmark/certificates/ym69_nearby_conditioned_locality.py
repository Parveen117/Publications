"""YM69: uniform finite near pins and averaged infinite distant exteriors.

Python 3.12 only. Exact posterior identities and directed estimates
accompany the written proofs; no interacting vacuum is discretized.
--check/default are read-only; --write updates this stage's evidence.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import itertools
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
import ym68_averaged_exterior_locality as previous

pointwise = previous.previous
native, p, y, spatial = previous.native, previous.p, previous.y, previous.spatial
Iv = previous.Iv
RESULT = HERE/'YM69_RESULT.json'
PIN = HERE/'EXPECTED_YM69.sha256'
SOURCES = HERE/'YM69_SOURCE_PINS.json'


@lru_cache(maxsize=None)
def logarithms():
    two = spatial.parameters(0)['log_two']
    # 40=32*(5/4) keeps the logarithm series argument near one.
    five_fourths = spatial.model.blocks.old.log_iv(Iv(F(5,4)),80)
    forty = p._r(Iv(5)*two+five_fourths)
    return two,forty


def joint_pointwise_budget(eta, distance, joint_links, block_factor):
    """sup_x E_Z |p_(X|Z)(x)-q_X(x)|, without an L1 cap of two.

    Combine the exponential with the ball-volume factor BEFORE rounding;
    otherwise a tiny covariance rounding floor would be magnified.
    """
    cell = pointwise.parameters(eta)
    for n in (distance,joint_links,block_factor):
        native.width(n)
    if not cell['eta']:
        return dict(bound=F(0),radius=None,method='exact_zero',exponent_lower=None)
    rho = pointwise.previous.parameters(eta)['rho'].lo
    log_two,log_forty = logarithms()
    s = joint_links
    K = pointwise.power_upper(cell['kappa'].hi,s)
    cap = pointwise.up(K-1/K)
    curvature = s*cell['density_second_ratio'].hi*K
    prefactor = 1+F(2*block_factor*s,3)+1/(3*rho*distance)
    candidates = [(cap,-1,None)]
    for order in range(81):
        exponent = rho*distance-s*log_forty.hi-3*s*order*log_two.hi
        if exponent < 0:
            continue  # Omitting a candidate can only weaken the upper bound.
        radius = F(1,2**order)
        first = prefactor*pointwise.decay_upper(exponent)
        value = pointwise.up(first+curvature*radius*radius)
        candidates.append((value,order,exponent))
    value,order,exponent = min(candidates,key=lambda item:item[0])
    return dict(bound=value,radius=None if order<0 else F(1,2**order),
                method='density_range' if order<0 else 'combined_exponential',
                exponent_lower=exponent)


def conditioned_budget(eta, distance, target_links, near_links, block_factor):
    cell = pointwise.parameters(eta)
    native.width(target_links)
    native.width(near_links,0)
    for n in (distance,block_factor):
        native.width(n)
    m,n = target_links,near_links
    joint = joint_pointwise_budget(eta,distance,m+n,block_factor)
    K_target = pointwise.power_upper(cell['kappa'].hi,m)
    K_near = pointwise.power_upper(cell['kappa'].hi,n)
    normalization = 2*K_near if n else F(1)
    cap = min(F(2),pointwise.up(K_target-1/K_target))
    error = min(cap,pointwise.up(normalization*joint['bound']))
    # With no near pins, retain the earlier (sometimes sharper) arithmetic.
    if not n:
        old = previous.exterior_budget(eta,distance,m,block_factor)
        error = min(error,old['expected_L1'])
    score,second = cell['score'].hi,cell['log_second'].hi
    log_cap = 8*m*spatial.pi_interval().hi*cell['wave_log_gradient'].hi
    coefficient = 6*m*second+2*m*score*score
    return dict(distance=distance,target_links=m,near_links=n,
                block_factor=block_factor,joint_links=m+n,
                fine_joint_links=block_factor*(m+n),joint_pointwise=joint,
                near_density_floor=1/K_near,posterior_normalization_upper=normalization,
                uniform_near_expected_L1=error,uniform_near_expected_TV=pointwise.up(error/2),
                conditional_mutual_information=pointwise.up(min(log_cap,log_cap*error)),
                conditional_relative_Fisher=pointwise.up(min(4*m*score*score,coefficient*error)),
                Fisher_coefficient=pointwise.up(coefficient),
                analytic_exponent=F(2,3*(m+n)+2))


def square_root_growth_threshold(eta):
    cell = pointwise.parameters(eta)
    if not cell['eta']:
        raise ValueError('The free case is exactly independent, with no finite threshold')
    rho = pointwise.previous.parameters(eta)['rho']
    log_K = p._r(Iv(4)*spatial.pi_interval()*cell['wave_log_gradient'])
    if log_K.lo <= 0:
        raise ArithmeticError('Positive log K is below interval precision')
    return p._r(rho/(Iv(3)*log_K))


def discrete_posterior_controls():
    """Independent finite probability law; not a native-vacuum approximation."""
    weights = (1,1,1,1,1,3,3,9)
    law = dict(zip(itertools.product(range(2),repeat=3),map(lambda v:F(v,sum(weights)),weights)))
    pz = {z:sum(v for (y,a,k),v in law.items() if k==z) for z in range(2)}
    rows = []
    for a in range(2):
        qa = sum(v for (y,n,z),v in law.items() if n==a)
        qy = {y:sum(law[y,a,z] for z in range(2))/qa for y in range(2)}
        actual = canceled = residual = wrong = F(0)
        posterior = []
        for z in range(2):
            az = sum(law[y,a,z] for y in range(2))
            nu = az/qa
            posterior.append(nu)
            difference = sum(abs(law[y,a,z]/az-qy[y]) for y in range(2))
            actual += nu*difference
            wrong += pz[z]*difference
            A = {y:law[y,a,z]/pz[z]-qa*qy[y] for y in range(2)}
            H = sum(A.values())
            canceled += pz[z]/qa*sum(abs(A[y]-qy[y]*H) for y in range(2))
            residual += pz[z]/qa*sum(abs(A[y]) for y in range(2))
        assert sum(posterior)==1 and canceled==actual and actual<=2*residual
        rows.append(dict(near_value=a,near_probability=qa,posterior_exterior=posterior,
                         conditional_L1=actual,canceled_identity=canceled,
                         joint_residual_over_near_density=residual,
                         wrong_unconditional_exterior_average=wrong))
    assert rows[0]['conditional_L1']>rows[0]['joint_residual_over_near_density']
    assert rows[0]['conditional_L1']!=rows[0]['wrong_unconditional_exterior_average']
    return dict(scope='Finite positive comparison law, not a YM discretization',cases=rows)


def native_posterior_controls():
    # w(x,a,z)=1+alpha*a*z+epsilon*x*a*z on three native links.
    alpha,epsilon = F(3,5),F(1,10)
    root = F(4,5)
    G = 2/(1+root)  # Integral 1/(1+alpha*z), using the native semicircle law.
    lower = sum(alpha**(2*k)*previous.coordinate_even_moment(k) for k in range(80))
    upper = lower+alpha**160/(1-alpha*alpha)
    assert lower <= G <= upper
    xz = {(1,0,0,0,1,0,0,0):F(1)}
    zpoly = {(0,0,0,0,1,0,0,0):F(1)}
    xpoly = {(1,0,0,0,0,0,0,0):F(1)}
    w = p.add({(0,)*8:F(1)},p.scale(zpoly,alpha),p.scale(xz,epsilon))
    assert native.phi(w,2)==1
    assert native.phi(p.mul(w,xpoly),2)==0
    cases = []
    for a,z in itertools.product((F(-1),F(0),F(3,5),F(1)),repeat=2):
        posterior = 1+alpha*a*z
        amplitude = epsilon*a*z/posterior
        assert posterior>0 and abs(amplitude)<1
        cases.append(dict(near_coordinate=a,remote_coordinate=z,
                          posterior_density=posterior,target_amplitude=amplitude,
                          score_at_target_i=amplitude/2))
    wrong_mean = epsilon/(4*alpha)*(1-G)
    assert wrong_mean==F(-1,216)
    return dict(scope='Smooth three-link native comparison density, not the YM vacuum',
                alpha=alpha,epsilon=epsilon,density_floor=1-alpha-epsilon,
                exact_resolvent_integral=G,resolvent_series_enclosure=Iv(lower,upper),
                correct_conditional_target_mean=F(0),wrong_unweighted_target_mean=wrong_mean,
                cases=cases)


def growing_pin_controls(eta):
    critical = square_root_growth_threshold(eta)
    assert critical.lo>3
    rows = [conditioned_budget(eta,n*n,4,n,2) for n in (128,256,512)]
    assert all(a['uniform_near_expected_L1']>b['uniform_near_expected_L1']
               for a,b in zip(rows,rows[1:]))
    return dict(critical_c_squared=critical,chosen_c=F(1),cases=rows,
                scope='n=sqrt(R) with fixed four-link target; not a filled spatial collar')


def source_checks():
    manifest = json.loads(SOURCES.read_text())
    checked = {}
    for path,expected in manifest['upstream_sha256'].items():
        actual = native.digest(ROOT/path)
        if actual!=expected:
            raise ValueError('upstream source changed: '+path)
        checked[path]=actual
    for path in manifest['local_inputs']:
        checked[path]=native.digest(ROOT/path)
    checked[str(SOURCES.relative_to(ROOT))]=native.digest(SOURCES)
    return checked


def run():
    native.previous.runtime_check()
    eta = F(1,16000)
    return spatial.encode(dict(
        certificate_type='YM69_EXACT_NEAR_PINS_AND_AVERAGED_DISTANT_EXTERIORS',
        verdict='PASS',runtime_policy='Python 3.12 only',input_sha256=source_checks(),
        fixed_near_budgets=[conditioned_budget(eta,R,4,4,2) for R in (256,512,1024,2048)],
        finite_posterior=discrete_posterior_controls(),native_posterior=native_posterior_controls(),
        growing_near_pins=growing_pin_controls(eta),
        model=dict(actual_YM62_vacuum_retained=True,unbounded_link_spaces_retained=True,
                   uniform_over_all_fixed_finite_near_pin_values=True,
                   distant_exterior_cardinality_independent=True,
                   countably_infinite_distant_exterior=True,
                   posterior_null_set_common_to_all_near_values=True,
                   conditioned_score_and_diagonal_drift_control=True,
                   sparse_growing_near_sets_controlled=True,
                   arbitrary_growing_near_collar_uniformity=False,
                   worst_case_distant_exterior_configuration_locality=False,
                   full_Schur_memory_kernel_locality=False,time_zero_Markov_identified=False,
                   continuum_R4_QFT=False,physical_mass_prediction=False),
        claim_status='EXACT_FINITE_NEAR_PINS__AVERAGED_FAR_EXTERIOR__FULL_COLLAR_AND_SCHUR_OPEN',
        proof_scope='Written Bayes and native regularity proofs with exact posterior and interval checks'))


def check(cert,result=RESULT,pin=PIN):
    sha=p.canonical_sha(cert)
    if json.loads(result.read_text())!=cert or pin.read_text().strip()!=sha:
        raise ValueError('YM69 fresh certificate/pin mismatch')
    return sha


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--write',action='store_true')
    group.add_argument('--check',action='store_true')
    args=parser.parse_args()
    cert=run()
    sha=p.canonical_sha(cert)
    if args.write:
        RESULT.write_text(json.dumps(cert,sort_keys=True,indent=2)+'\n')
        PIN.write_text(sha+'\n')
    else:
        check(cert)
    print('YM69 PASS',sha)
    print(json.dumps(dict(runtime=sys.version.split()[0],uniform_exact_near_pins=True,
                          averaged_infinite_far_exterior=True,full_collar_locality=False,
                          full_memory_locality=False,continuum_R4_QFT=False)))


if __name__=='__main__':
    main()
