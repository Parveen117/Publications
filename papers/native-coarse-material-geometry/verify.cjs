'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs'),f=N.F.of,C=N.Cut.of;
const M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const checks=[],negative_controls=[],proofs=[],fixtures=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const sha=x=>crypto.createHash('sha256').update(x).digest('hex');
const tokens=['t','s','A','B','q','x','n'];
const rules=[];for(let j=0;j<tokens.length;j++)for(let i=0;i<j;i++)rules.push({id:'commute_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'commuting Gaussian generating parameters and positive scalar filter bounds'});
const alg=new M.Presentation({tokens,rules}),audit=alg.audit();
ck('scalar_presentation_confluent',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');
const presentations=[{spec:alg.spec(),audit:{status:audit.status,criticalPairs:audit.criticalPairs,resolved:audit.resolutions.length,open:audit.open.length,full_audit_sha256:sha(JSON.stringify(audit))}}];
const V=Object.fromEntries(tokens.map(t=>[t,alg.word([t])])),one=alg.word([]),zero=alg.word([],0);
const pow=(p,n)=>{let r=one;for(let i=0;i<n;i++)r=r.times(p);return r;};
function prove(name,lhs,rhs){const r=alg.reduce(lhs.minus(rhs),{witness:true});ck(name,r.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(r.certificate));ck(name+':replay',alg.replay(certificate));proofs.push({name,presentation_sha256:alg.hash(),lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate});}
prove('complex_Gaussian_generating_exponent',pow(V.t.plus(V.s),2).minus(pow(V.t.minus(V.s),2)).scale('1/4'),V.t.times(V.s));
let geom=zero;for(let k=0;k<7;k++)geom=geom.plus(pow(V.B,6-k).times(pow(V.A,k)));
prove('positive_power_difference_factorization',pow(V.B,7).minus(pow(V.A,7)),V.B.minus(V.A).times(geom));
prove('filter_scale_comparison',one.plus(V.q.times(V.x)).minus(V.q.times(one.plus(V.x))),one.minus(V.q));
prove('tail_comparison_numerator',pow(V.n,2).minus(V.n.times(V.n.minus(one))),V.n);
// If x=n^2-1>=0, the following has nonnegative coefficients.
const binom=(n,k)=>{let r=f(1);for(let i=1;i<=k;i++)r=r.mul(n-i+1).div(i);return r;};
const tailCoeffs=Array.from({length:8},(_,k)=>binom(7,k).mul(f(128).sub(f(2).pow(7-k))));
const tailExpr=tailCoeffs.reduce((r,c,k)=>r.plus(pow(V.x,k).scale(c.toString())),zero);
prove('all_n_tail_numerator_bound',pow(one.plus(V.x),7).scale(128).minus(pow(one.scale(2).plus(V.x),7)),tailExpr);
for(let k=0;k<tailCoeffs.length;k++)ck('nonnegative_tail_bound_coefficient:'+k,f(0).le(tailCoeffs[k]));
fixtures.push({kind:'general_tail_bound_polynomial',substitution:'x=n^2-1>=0',coefficients:tailCoeffs.map(x=>x.toString()),conclusion:'(1+n^2)^7 <= 128 n^14 for n>=1'});
no('wrong_complex_Gaussian_normalization',alg.reduce(pow(V.t.plus(V.s),2).minus(pow(V.t.minus(V.s),2)).scale('1/2').minus(V.t.times(V.s))).normal.terms.size!==0);
no('wrong_filter_scale_sign',alg.reduce(one.plus(V.q.times(V.x)).minus(V.q.times(one.plus(V.x))).minus(V.q.minus(one))).normal.terms.size!==0);

// Exact circular Gaussian moments from independent real variance 1/2.
const factorial=n=>{let r=f(1);for(let k=2;k<=n;k++)r=r.mul(k);return r;};
const moment=(n,var_)=>{if(n%2)return f(0);let r=f(1);for(let k=1;k<n;k+=2)r=r.mul(k).mul(var_);return r;};
const cpow=(z,n)=>{let r=N.ONE;for(let k=0;k<n;k++)r=r.mul(z);return r;};
function circular(a,b,var_){let r=N.ZERO;for(let k=0;k<=a;k++)for(let l=0;l<=b;l++){
 const c=binom(a,k).mul(binom(b,l)).mul(moment(a+b-k-l,var_)).mul(moment(k+l,var_));
 r=r.add(cpow(N.IOTA,k).mul(cpow(C([0,-1]),l)).mul(c));
}return r;}
for(let a=0;a<=6;a++)for(let b=0;b<=6;b++)ck('Hermite_Bargmann_Gram:'+a+':'+b,circular(a,b,'1/2').eq(a===b?factorial(a):0));
for(const a of [1,2,3,4])no('complex_variance_one_per_real_axis:'+a,!circular(a,a,1).eq(factorial(a)));
ck('two_mode_marginal_orthogonality',circular(2,2,'1/2').mul(circular(3,3,'1/2')).eq(factorial(2).mul(factorial(3))));
ck('unused_mode_offdiagonal_vanishes',circular(2,3,'1/2').zero());

// The selected rational resolvent response uses native F throughout.
const filter=(ell,n)=>f(1).div(f(1).add(f(ell).pow(2).mul(f(n).pow(2))).pow(4));
const base=n=>f(1).add(f(n).pow(2));
for(const ell of ['1','1/2','1/4']){
 let vH3=f(1),uH3=f(1);
 for(let n=1;n<=24;n++){
  const b=base(n),t=filter(ell,n),vterm=b.pow(7).mul(t.pow(2)),uterm=t.pow(2).div(b);
  vH3=vH3.add(vterm.mul(2));uH3=uH3.add(uterm.mul(2));
  ck('filter_tail_bound:'+ell+':'+n,vterm.le(f(128).div(f(ell).pow(16).mul(f(n).pow(2)))));
  ck('classical_scale_bound:'+ell+':'+n,vterm.le(f(1).div(f(ell).pow(14).mul(f(1).add(f(ell).pow(2).mul(f(n).pow(2)))))));
  if(n<=4){
   ck('strict_contraction:'+ell+':'+n,f(0).le(t)&&t.le(1)&&!t.eq(1));
   const lambda=f(1).div(b.pow(4)),qu=f(2).mul(lambda).mul(t.pow(2)),qv=f(1).div(lambda.mul(2)).mul(t.pow(2));
   ck('filtered_joint_noise_product:'+ell+':'+n,qu.mul(qv).eq(t.pow(4)));
  }
 }
 fixtures.push({kind:'finite_H3_covariance_budget_h1',ell,cutoff:24,configuration:uH3.mul(2).toString(),conjugate:vH3.div(2).toString(),scope:'finite sums only; infinite tails bounded analytically in theorem'});
 ck('zero_mode_retained:'+ell,filter(ell,0).eq(1));
}
for(let n=0;n<=6;n++){
 const lambda=f(1).div(base(n).pow(4)),t=filter(1,n);
 ck('ell_one_filtered_momentum_covariance:'+n,t.pow(2).div(lambda.mul(2)).eq(lambda.div(2)));
 ck('ell_one_filtered_configuration_covariance:'+n,lambda.mul(t.pow(2)).mul(2).eq(lambda.pow(3).mul(2)));
}
no('smoothing_is_canonical_coordinate_change',!filter(1,1).pow(2).eq(1));
// For u=cos(2pi z), constant coefficients of R(u^2) and (Ru)^2 differ.
no('filter_commutes_with_nonlinear_inverse',!f('1/2').eq(filter(1,1).pow(2).div(2)));
no('omit_inverse_covariance_in_momentum_noise',!base(2).pow(4).eq(f(1).div(base(2).pow(4))));

// Positive weighted circumference and the retained global normalization.
for(const sign of [-1,1]){
 const mass=[f(2).mul(sign),f(3).mul(sign),f(5).mul(sign)],eTheta=[f('1/2'),f(3),f('7/4')];
 const P=mass.reduce((s,x)=>s.add(x),f(0)),weights=mass.map(x=>x.div(P));
 const L=weights.reduce((s,w,i)=>s.add(w.mul(eTheta[i])),f(0));
 ck('length_weights_normalized:'+sign,weights.reduce((s,x)=>s.add(x),f(0)).eq(1));
 ck('length_weights_positive:'+sign,weights.every(x=>f(0).le(x)));
 ck('length_between_theta_extrema:'+sign,f('1/2').le(L)&&L.le(3));
 no('omit_global_P_from_length:'+sign,!mass.reduce((s,x,i)=>s.add(x.mul(eTheta[i])),f(0)).eq(L));
}
for(const delta of ['1/4','1/2','1'])for(const j of ['1','2','5']){
 const y=f(3).div(j),bound=f(3).div(delta);
 ck('accepted_inverse_density_bound:'+delta+':'+j,y.le(bound));
}
no('threshold_bound_without_acceptance',!f(3).div('1/8').le(f(3).div('1/4')));
// Native rate witness: ell=1/k, h=ell^16 gives h ell^-15=1/k.
for(const k of [2,4,8,16]){
 const ell=f(1).div(k),h=ell.pow(16);
 ck('joint_classical_scale:'+k,h.div(ell.pow(15)).eq(ell));
 no('fixed_hbar_removes_scale:'+k,!f(1).div(ell.pow(15)).le(1));
}

// Finite dilation regression: two success outcomes (lengths 2 and 5)
// and a separately recorded rejection. This is not a cutoff continuum proof.
const W=N.matrix([['3/5',0],[0,1],['4/5',0]]),Wd=N.dagger(W);
const diag=vs=>N.matrix(vs.map((v,i)=>vs.map((_,j)=>i===j?v:0))),compress=vs=>N.mul(N.mul(Wd,diag(vs)),W);
const success=compress([1,1,0]),reject=compress([0,0,1]),mean=compress([2,5,0]),second=compress([4,25,0]);
ck('detector_analysis_isometry',N.equal(N.mul(Wd,W),N.identity(2)));
ck('success_plus_rejection_normalized',N.equal(N.add(success,reject),N.identity(2)));
ck('length_mean_selfadjoint',N.equal(mean,N.dagger(mean)));
ck('length_mean_positive_diagonal',mean[0][0].eq('18/25')&&mean[1][1].eq(5));
const p0=f('9/25'),p1=f(1),l0=f('18/25'),l1=f(5);
ck('accepted_conditional_mean_first_state',l0.div(p0).eq(2));
ck('accepted_conditional_mean_second_state',l1.div(p1).eq(5));
const mixedConditional=l0.add(l1).div(p0.add(p1));
ck('mixed_state_conditional_reading',mixedConditional.eq('143/34'));
no('conditional_mean_is_linear_observable',!mixedConditional.eq(l0.div(p0).add(l1.div(p1)).div(2)));
no('discard_rejection_probability',!N.equal(success,N.identity(2)));
no('mean_operator_spectrum_is_POVM',!N.equal(second,N.mul(mean,mean)));
fixtures.push({kind:'finite_flagged_detector',analysis:W,success,reject,length_mean:mean,length_second_moment:second,mixed_conditional_mean:mixedConditional.toString()});
const bad=JSON.parse(JSON.stringify(proofs.find(p=>p.certificate.steps.length).certificate));bad.steps[0].coefficient=['997','0'];
no('corrupted_native_replay',(()=>{try{alg.replay(bad);return false;}catch{return true;}})());

const payload={protocol:'CMG1_COARSE_MATERIAL_GEOMETRY_V1',status:'PASS_CONTINUUM_COARSE_GEOMETRY_WITH_DECLARED_THRESHOLD',source_pins_sha256:sha(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),checks,negative_controls,presentations,proofs,fixtures,results:{physical_space:'unchanged selected RGR/MPQ quotient',countable_POVM_and_finite_marginals:true,filter:'(I+ell^2 Lambda)^(-4)',threshold:'min exp(u_ell)>=delta with rejection probability retained',fixed_hbar_coarse_C2_profiles:true,all_normal_states_have_profile_outcomes_or_rejection:true,positive_length_mean_operator_dense_domain:true,vacuum_accepted_length_moments_finite:true,all_normal_states_have_finite_length_moments:false,ungated_vacuum_mean_length:'infinite even after smoothing',fixed_hbar_zero_response_scale_point_limit:false,classical_mean_length_recovery:'hbar ell^(-15) -> 0, fixed threshold below classical minimum |j|',response_scale_and_threshold_derived_from_source:false,sharp_continuum_r_operator:false,spacetime_curvature_operator:false,residual_gauge_completeness_assumed:false,physical_three_plus_one_gravity_closed:false,continuum_proofs_machine_certified:false},summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,negative_controls:negative_controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',hash=sha(bytes),cp=path.join(__dirname,'CERTIFICATE.json'),hp=path.join(__dirname,'EXPECTED.sha256');
if(process.argv.includes('--write')){fs.writeFileSync(cp,bytes);fs.writeFileSync(hp,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cp,'utf8'),bytes);assert.equal(fs.readFileSync(hp,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
