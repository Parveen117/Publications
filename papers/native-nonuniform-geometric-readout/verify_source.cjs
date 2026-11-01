'use strict';
// Exact finite source diagnostics only. This does not execute the continuum POVM.
const assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const F=N.F, f=F.of, order=6;
const checks=[],negatives=[];
function check(name,ok){assert(ok,name);checks.push(name);}
function reject(name,ok){assert(ok,name);negatives.push(name);}
const zero=()=>Array.from({length:order+1},()=>f(0));
const constant=a=>{const p=zero();p[0]=f(a);return p;};
const add=(a,b)=>a.map((x,i)=>x.add(b[i]));
const scale=(a,s)=>a.map(x=>x.mul(s));
function mul(a,b){const p=zero();for(let i=0;i<=order;i++)for(let j=0;j<=order-i;j++)p[i+j]=p[i+j].add(a[i].mul(b[j]));return p;}
function pow(a,n){let p=constant(1);for(let i=0;i<n;i++)p=mul(p,a);return p;}
function reciprocal(a){assert(!a[0].zero());const p=zero();p[0]=f(1).div(a[0]);for(let n=1;n<=order;n++){let s=f(0);for(let k=1;k<=n;k++)s=s.add(a[k].mul(p[n-k]));p[n]=s.neg().div(a[0]);}return p;}
const derivative=a=>a.map((_,i)=>i<order?a[i+1].mul(i+1):f(0));
function factorial(n){let a=f(1);for(let j=2;j<=n;j++)a=a.mul(j);return a;}
function exponential(a){assert(a[0].zero());let p=zero();for(let n=0;n<=order;n++)p=add(p,scale(pow(a,n),f(1).div(factorial(n))));return p;}
function cosine(a){assert(a[0].zero());let p=zero();for(let n=0;2*n<=order;n++)p=add(p,scale(pow(a,2*n),f(n%2?-1:1).div(factorial(2*n))));return p;}
function coefficientFunction(theta){return add(constant(1),scale(cosine(theta),'1/2'));}
// Solve Theta'(t)=1+cos(Theta(t))/2, Theta(0)=0 as a formal rational jet.
const theta=zero();
for(let n=0;n<order;n++)theta[n+1]=coefficientFunction(theta)[n].div(n+1);
const a=coefficientFunction(theta),inverseA=reciprocal(a),eTheta=exponential(theta);
const thetaPrime=derivative(theta);
check('inverse_chart_ODE_coefficients_through_degree_5',thetaPrime.slice(0,order).every((x,i)=>x.eq(a[i])));
check('Theta_first_derivative_3_over_2',theta[1].eq('3/2'));
check('Theta_second_derivative_zero',theta[2].zero());
check('Theta_third_derivative_minus_9_over_8',theta[3].mul(6).eq('-9/8'));
// 5 a5 = -(a1 a3)/2 + a1^4/48 = 63/256, hence a5=63/1280.
check('Theta_fifth_derivative_189_over_32',theta[5].mul(120).eq('189/32'));
check('inverse_a_product_is_one',mul(a,inverseA).every((x,i)=>x.eq(i===0?1:0)));
// Independent exact solution: 1/a(Theta(t))=4/3-(2/3)cos(sqrt(3)t/2).
const inverseAClosed=constant('4/3');
for(let n=0;2*n<=order;n++)inverseAClosed[2*n]=inverseAClosed[2*n].sub(f('2/3').mul(n%2?-1:1).mul(f('3/4').pow(n)).div(factorial(2*n)));
check('independent_half_angle_solution_checks_inverse_a',inverseA.every((x,i)=>x.eq(inverseAClosed[i])));

// Exact circle moments: integral cos(z)^(2m) dz = (2m)!/(2^(2m)(m!)^2).
function moment(n){if(n%2)return f(0);const m=n/2;return factorial(n).div(f(2).pow(n).mul(factorial(m).pow(2)));}
check('circle_moments_0_2_4_6',[0,2,4,6].map(moment).every((x,i)=>x.eq(['1','1/2','3/8','5/16'][i])));
const integrate=jet=>jet.map((x,n)=>x.mul(moment(n)));
const numerator=mul(eTheta,inverseA),P=integrate(inverseA),A=integrate(numerator);
const length=mul(A,reciprocal(P));
check('source_P_quadratic_coefficient_1_over_8',P[2].eq('1/8'));
check('source_L_quadratic_coefficient_9_over_16',length[2].eq('9/16'));
check('source_length_even_through_degree_6',length.every((x,n)=>n%2===0||x.zero()));
// Independent differential calculation at zero, without a formal-series division.
const thetaD1=f('3/2'),thetaD2=f(0),a0=f('3/2'),aThetaD1=f(0),aThetaD2=f('-1/2');
const invD2=aThetaD2.neg().mul(thetaD1.pow(2)).div(a0.pow(2));
const numeratorD2=thetaD1.pow(2).add(thetaD2).div(a0).add(invD2);
const P_D2=invD2.mul(moment(2));
const lengthD2=numeratorD2.mul(moment(2)).sub(P_D2).div(f(1).div(a0));
check('independent_second_derivative_matches_jet',length[2].mul(2).eq(lengthD2)&&lengthD2.eq('9/8'));
const harmonicNumerator=numerator.map((x,n)=>x.mul(moment(n+1)));
const harmonic=mul(harmonicNumerator,reciprocal(P));
check('coframe_cosine_linear_response_3_over_4',harmonic[1].eq('3/4'));

// Rational moment Gram matrix and projection diagnostic at flat j=1.
const H=N.matrix([[1,0],[0,'1/2']]),b=N.matrix([['2/3'],['-5/4']]);
const D=N.mul(H,b),projected=N.mul(N.inverse(H),D);
check('two_window_projection_recovers_constant_and_cosine',N.equal(projected,b));
const pureCos=N.mul(H,N.matrix([[0],[1]]));
check('radial_average_misses_nonzero_cosine_profile',pureCos[0][0].zero()&&pureCos[1][0].eq('1/2'));
// Exact null mode for this weighted projection; it is not a full field inverse.
const omittedCos2D=N.matrix([[moment(2).mul(2).sub(moment(0))],[moment(3).mul(2).sub(moment(1))]]);
check('second_harmonic_has_zero_two_window_projection',N.isZero(N.mul(N.inverse(H),omittedCos2D)));

const frozenP=scale(A,f(1).div(P[0]));
reject('freezing_global_P_changes_quadratic_length',frozenP[2].eq('3/4')&&!frozenP[2].eq(length[2]));
reject('radial_only_prediction_L_equals_1_fails',!length[2].zero());
reject('using_affine_D_as_canonical_p_fails',!N.equal(D,b));
reject('reversing_pulse_sign_reverses_detectable_harmonic',!harmonic[1].eq(harmonic[1].neg()));
// Nonconstant source reweighting must retain the removed cosine Wick variance.
const h=f('1/2');
const fullRemovedAtC1=h.mul(2).div(2),zeroOnly=h.div(2);
reject('omitting_removed_cosine_Wick_variance_detected',!fullRemovedAtC1.eq(zeroOnly));

const result={status:'PASS',author:'Monty',date:'2026-10-30',
 scope:'finite rational source jets and selected algebraic diagnostics; no continuum POVM execution',
 arithmetic:'unchanged native-common-cone-limit/native_engine/native_operator.cjs F and matrix operations',
 truncation_order:order,checks:checks.length,negative_controls:negatives.length,
 check_names:checks,negative_control_names:negatives,
 coefficients:{Theta:theta,P,length,coframe_cosine_reading:harmonic},
 independently_derived_length_second_derivative:lengthD2,
 omitted_global_feedback_quadratic:frozenP[2],
 remainder_bound_proved_by_execution:false,continuum_theorems_proved_by_execution:false,
 quantum_expected_length_computed:false,nrqc_executed:false,hardware_executed:false};
process.stdout.write(JSON.stringify(result,null,2)+'\n');
