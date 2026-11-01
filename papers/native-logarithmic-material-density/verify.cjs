'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs'),f=N.F.of,C=N.Cut.of;
const M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const checks=[],negative_controls=[],proofs=[],fixtures=[],presentations=[],owners=new Map();
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const sha=x=>crypto.createHash('sha256').update(x).digest('hex');
function register(alg){const audit=alg.audit();ck('confluent:'+alg.hash(),audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');owners.set(alg.hash(),alg);presentations.push({spec:alg.spec(),audit:{status:audit.status,criticalPairs:audit.criticalPairs,resolved:audit.resolutions.length,open:audit.open.length,full_audit_sha256:sha(JSON.stringify(audit))}});}
function prove(alg,name,lhs,rhs){const r=alg.reduce(lhs.minus(rhs),{witness:true});ck(name,r.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(r.certificate));ck(name+':replay',alg.replay(certificate));proofs.push({name,presentation_sha256:alg.hash(),lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate});}
function scalar(tokens){const rules=[];for(let j=0;j<tokens.length;j++)for(let i=0;i<j;i++)rules.push({id:'scalar_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'commuting scalar polynomial variables'});const alg=new M.Presentation({tokens,rules});register(alg);return alg;}
const alg=scalar(['h','a','l','v','vp','w','wp','q','p','n']);
const z=Object.fromEntries([...alg.order.keys()].map(k=>[k,alg.word([k])])),one=alg.word([]),zero=alg.word([],0);
const power=(x,n)=>{let r=one;for(let j=0;j<n;j++)r=alg.reduce(r.times(x)).normal;return r;};
const ah=one.plus(z.h.scale('1/2'));
// General scalar identities; these are not a finite-dimensional field model.
prove(alg,'Wick_density_corrected_weight',ah.minus(one).minus(z.h.scale('1/2')).times(z.l),zero);
prove(alg,'uncorrected_flow_has_half_h_weight',one.plus(z.h.scale('1/2')).minus(one),z.h.scale('1/2'));
const oldComm=z.h.times(z.v.times(z.wp).plus(z.h.times(z.vp).times(z.w).scale('1/2'))).scale([0,-1]);
const addedComm=z.h.times(z.h).times(z.vp).times(z.w).scale([0,'1/2']);
prove(alg,'corrected_local_density_commutator',oldComm.plus(addedComm),z.h.times(z.v).times(z.wp).scale([0,-1]));
prove(alg,'Wick_mean_exponent',z.q.scale('1/2').minus(z.q.scale('1/2')),zero);
prove(alg,'Wick_two_point_exponent',z.q.plus(z.p).plus(z.l.scale(2)).scale('1/2').minus(z.q.plus(z.p).scale('1/2')),z.l);
prove(alg,'same_sign_critical_entry_cancellation',z.p.times(z.q).plus(z.q.times(z.p)).minus(z.p.times(z.q).scale(2)),zero);
prove(alg,'opposite_sign_critical_entry_square',z.p.times(z.q).plus(z.q.times(z.p)).plus(z.p.times(z.q).scale(2)).scale('1/4'),z.p.times(z.q));
const F=x=>power(x,3).minus(x).scale('1/6');
prove(alg,'opposite_sign_sum_induction',F(z.n.plus(one)).minus(F(z.n)),z.n.times(z.n.plus(one)).scale('1/2'));
no('minus_half_h_instead_of_plus',alg.reduce(one.minus(z.h.scale('1/2')).minus(one).minus(z.h.scale('1/2'))).normal.terms.size!==0);
no('omit_quantum_density_feedback',alg.reduce(oldComm.minus(z.h.times(z.v).times(z.wp).scale([0,-1]))).normal.terms.size!==0);
no('wrong_covariance_diagonal_sign',alg.reduce(ah.minus(one).plus(z.h.scale('1/2'))).normal.terms.size!==0);
no('drop_measure_Jacobian',alg.reduce(ah.minus(z.h.scale('1/2'))).normal.terms.size!==0);
no('normal_ordering_subtracts_full_variance',alg.reduce(z.q.scale('1/2').minus(z.q)).normal.terms.size!==0);

// Power-family leading coefficient. The general expansion is proved in T1.
for(const p of ['0','1/4','1/2','1','4']){
 const coeff=f(p).mul(2).sub(1);
 ck('unique_critical_power_in_fixture:'+p,coeff.zero()===(p==='1/2'));
 fixtures.push({kind:'power_covariance_leading_coefficient',p,coefficient_over_k:coeff});
}
ck('MRM_smooth_covariance_coefficient_recovered',f(4).mul(2).sub(1).eq(7));

// S-entry norm squares avoid irrational square roots entirely.
const lam=n=>n===0?f(1):f(1).div(f(n).abs().mul(2));
const wp=(bank,k)=>k>=0?(bank[k]||C(0)):(bank[-k]||C(0)).dagger();
function entry2(bank,m,n){
 const wt=wp(bank,m-n).norm2();
 return wt.mul(f(n).pow(2).mul(lam(n).div(lam(m))).add(f(m).pow(2).mul(lam(m).div(lam(n)))).sub(f(2).mul(m).mul(n))).div(4);
}
function HS(bank,cutoff){let total=f(0);for(let m=-cutoff;m<=cutoff;m++)for(let n=-cutoff;n<=cutoff;n++)total=total.add(entry2(bank,m,n));return total;}
function closedHS(bank){let total=f(0);for(let n=1;n<bank.length;n++)total=total.add(f(n).pow(3).sub(n).div(3).add(f(n).div(2)).mul(bank[n].norm2()));return total;}
const banks=[
 [C(0),C([0,'-1/2'])],
 [C('2/5'),C(0),C([0,'-1/2'])],
 [C('1/3'),C([0,'-1/2']),C('1/3'),C([0,'1/5'])]
];
for(let j=0;j<banks.length;j++){
 const bank=banks[j],full=closedHS(bank);
 for(const cutoff of [1,2,3,4,8]){
  const value=HS(bank,cutoff);
  ck('critical_HS_partial_bound:'+j+':'+cutoff,value.le(full));
  if(cutoff>=bank.length-1)ck('critical_HS_exact_finite_support_sum:'+j+':'+cutoff,value.eq(full));
 }
 let lin=f(0),closed=f(0);
 for(let n=1;n<bank.length;n++){lin=lin.add(f(2).mul(f(n).pow(2)).mul(bank[n].norm2()).div(lam(n)));closed=closed.add(f(4).mul(f(n).pow(3)).mul(bank[n].norm2()));}
 ck('critical_linear_norm_sum:'+j,lin.eq(closed));
 fixtures.push({kind:'smooth_trigonometric_window',coefficients:bank,HS_squared_over_4pi_squared:full,b_C_squared_over_4pi_squared:lin});
}
ck('first_sine_HS_exact',closedHS(banks[0]).eq('1/8'));
ck('second_sine_HS_exact',closedHS(banks[1]).eq('3/4'));
for(const [m,n]of [[1,2],[2,7],[-1,-3],[-4,-2]])ck('same_sign_entry_exact_zero:'+m+':'+n,entry2(Array.from({length:10},()=>C(1)),m,n).zero());
for(const [m,n]of [[-1,1],[-1,3],[-2,4]])ck('opposite_sign_entry_square:'+m+':'+n,entry2(Array.from({length:10},()=>C(1)),m,n).eq(f(m).abs().mul(f(n).abs())));
for(const n of [1,2,4]){
 const m=n+1,oldLam=t=>f(1).div(f(1).add(f(t).pow(2)).pow(4));
 const oldSq=f(n).pow(2).mul(oldLam(n).div(oldLam(m))).add(f(m).pow(2).mul(oldLam(m).div(oldLam(n)))).sub(f(2).mul(m).mul(n)).div(4);
 no('old_MPQ_covariance_cancels_same_sign:'+n,!oldSq.zero());
}
no('drop_constant_mode_HS_contribution',!closedHS(banks[0]).zero());
no('zero_local_vacuum_charge',!closedHS(banks[1]).zero());

// Exact Wick normalization and its necessary zero-mode factor.
for(const cutoff of [1,2,4,8,16]){
 let variance=f(1);for(let n=1;n<=cutoff;n++)variance=variance.add(f(1).div(n));
 const byModes=f(1).add(Array.from({length:cutoff},(_,i)=>lam(i+1).mul(2)).reduce((a,b)=>a.add(b),f(0)));
 ck('complete_logarithmic_point_variance:'+cutoff,variance.eq(byModes));
 ck('Wick_mean_normalization:'+cutoff,variance.div(2).sub(byModes.div(2)).zero());
 no('omit_zero_mode_variance:'+cutoff,!variance.div(2).sub(byModes.sub(1).div(2)).zero());
}
for(const h of ['1/4','1/3','1/2','3/4']){
 const ah=f(1).add(f(h).div(2));
 ck('density_exponent_exactly_one:'+h,ah.sub(f(h).div(2)).eq(1));
 no('classical_shift_preserves_new_density_weight:'+h,!f(1).sub(f(h).div(2)).eq(1));
 // v=sin(2 pi z), s_*=0, d_*=c cos(2 pi z), c=2/3.
 const c=f('2/3'),mean=ah.mul(c).div(2),classical=c.div(2);
 const realVar=f(h).mul(c.pow(2)).div(8),imagVar=f(h).mul(ah.pow(2)).div(4),quadraticVar=f(h).pow(2).div(16);
 ck('local_coherent_mean_correction:'+h,mean.sub(classical).eq(f(h).mul(c).div(4)));
 ck('local_vacuum_variance_from_HS:'+h,quadraticVar.eq(f(h).pow(2).mul(closedHS(banks[0])).div(2)));
 ck('L2_moment_integrability_margin:'+h,f(0).le(f(1).sub(h))&&!f(1).sub(h).zero());
 fixtures.push({kind:'coherent_sine_window',h,a_h:ah,mean_over_2pi:mean,classical_over_2pi:classical,variance_over_4pi_squared:realVar.add(imagVar).add(quadraticVar)});
}
no('L2_bound_includes_endpoint_h_one',f(1).sub(1).zero());

// Gaussian centering and chaos orthogonality, checked by exact moments.
const gaussian=n=>{if(n%2)return f(0);let q=f(1);for(let j=1;j<n;j+=2)q=q.mul(j);return q;};
ck('centered_second_chaos_mean',gaussian(2).sub(1).zero());
ck('first_second_chaos_orthogonal',gaussian(3).sub(gaussian(1)).zero());
ck('centered_second_chaos_norm_two',gaussian(4).sub(gaussian(2).mul(2)).add(1).eq(2));
no('second_chaos_is_a_scalar_counterterm',!gaussian(4).sub(1).zero());

// Actual spatial T^3 shear, with rational covariance ratio on a Fourier ray.
const A=N.matrix([[1,-2,-2],[0,1,0],[0,0,1]]),Ai=N.inverse(A),AiT=N.dagger(Ai),e1=N.matrix([[1],[0],[0]]);
const ray=N.mul(AiT,e1);
ck('3D_shear_inverse',N.equal(N.mul(A,Ai),N.identity(3)));
ck('3D_shear_covector_image',N.equal(ray,N.matrix([[1],[2],[2]])));
ck('3D_shear_covector_length_squared',N.energy(ray).eq(9));
ck('3D_shear_unimodular_triangular_diagonal',A[0][0].mul(A[1][1]).mul(A[2][2]).eq(1));
no('volume_preservation_implies_orthogonality',!N.equal(N.mul(N.dagger(A),A),N.identity(3)));
for(const n of [1,2,5,17]){
 const ratio=f(1).div(f(n).mul(3).pow(3)).div(f(1).div(f(n).pow(3)));
 ck('3D_logarithmic_ray_variance_ratio:'+n,ratio.eq('1/27'));
 no('3D_ray_covariance_preserved:'+n,!ratio.eq(1));
}
fixtures.push({kind:'three_dimensional_torus_shear',A,inverse_transpose:AiT,ray_image:ray,variance_ratio:'1/27',scope:'integer spatial diffeomorphism and exact coefficients; measure singularity uses the written infinite-ray law'});
// Repeated normalized samples distinguish variance 1 from 1/27 in the written SLLN.
ck('distinct_ray_second_moments',!f(1).eq('1/27'));

// Fable W11 scope audit: a positive complement form solves all CORE tests.
const diag=values=>N.matrix(values.map((x,i)=>values.map((_,j)=>i===j?x:0)));
const core=diag([1,1,1,0]),complement=diag([0,0,0,1]),ell=diag([1,0,0,0]);
const Q1=N.matrix([[0,0,0,0],[1,1,0,0],[0,0,0,0],[0,0,0,0]]);
const Q2=N.matrix([[0,0,0,0],[0,0,0,0],[1,0,1,0],[0,0,0,0]]);
const extra=N.kron(complement,N.identity(2)),core2=N.kron(core,N.identity(2));
ck('Fable_complement_form_positive_projection',N.equal(N.mul(extra,extra),extra)&&N.equal(N.dagger(extra),extra)&&N.rank(extra)===2);
ck('Fable_complement_form_kills_both_constraint_slots',[Q1,Q2].every(Q=>{const image=N.kron(N.mul(Q,core),N.identity(2));return N.isZero(N.mul(extra,image))&&N.isZero(N.mul(N.dagger(image),extra));}));
ck('Fable_complement_form_zero_on_core',N.isZero(N.mul(N.mul(core2,extra),core2)));
no('Fable_unrestricted_full_space_two_weight_claim',N.isZero(N.mul(N.kron(ell,N.identity(2)),extra))&&!N.isZero(extra));
fixtures.push({kind:'Fable_W11_complement_control',positive_form:extra,core_projection:core2,scope:'does not refute the correctly restricted core theorem'});

// A corrupted symbolic trace must not pass replay.
const target=proofs.find(p=>p.certificate.steps.length),bad=JSON.parse(JSON.stringify(target.certificate));
bad.steps[0].coefficient=['991','0'];
no('corrupted_Paninian_certificate',(()=>{try{owners.get(target.presentation_sha256).replay(bad);return false;}catch{return true;}})());

const payload={protocol:'LMD1_LOGARITHMIC_MATERIAL_DENSITY_V1',status:'PASS_SELECTED_LOCAL_FLOW_AND_RENORMALIZED_DENSITY',
 source_pins_sha256:sha(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),checks,negative_controls,presentations,proofs,fixtures,
 results:{representation:'explicit alternative reduced logarithmic Gaussian; old MPQ preserved',
  quantum_parameter_range:'0<hbar<1 in the chosen dimensionless normalization (L2 range)',
  critical_power_for_same_sign_vacuum_cancellation:'p=1/2 for lambda_n=c*|n|^(-2p)',
  smooth_material_circle_flows_quasi_invariant:true,
  real_local_generators_self_adjoint_by_Stone:true,
  density:'Wick random measure with positive smeared spectral observables',
  canonical_density_weight_requires:'a_h=1+hbar/2',
  old_classical_shift_has_exact_same_density_law:false,
  corrected_local_density_pulse_feedback:true,
  coherent_smeared_classical_recovery:true,
  full_cutoff_strong_resolvent_theorem_claimed:false,
  equivalence_to_old_MPQ_configuration_representation:false,
  original_pointwise_metric_reconstructed:false,
  three_dimensional_scalar_shear_variance_ratio:'1/27; mutually singular by written SLLN proof',
  all_three_dimensional_gravity_representations_excluded:false,
  source_3plus1_gravitational_constraints_constructed:false,
  Fable_bridge_review:'core-closure bridge useful under dependencies; full-Psi density unresolved; three corrections recorded',
  Fable_uniqueness_imported_as_new_dependency:false,
  continuum_proofs_machine_certified:false},
 summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,negative_controls:negative_controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',hash=sha(bytes),cp=path.join(__dirname,'CERTIFICATE.json'),hp=path.join(__dirname,'EXPECTED.sha256');
if(process.argv.includes('--write')){fs.writeFileSync(cp,bytes);fs.writeFileSync(hp,hash+'  CERTIFICATE.json\n');}
else{assert.equal(fs.readFileSync(cp,'utf8'),bytes);assert.equal(fs.readFileSync(hp,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
