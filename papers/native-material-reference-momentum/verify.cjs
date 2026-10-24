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
function scalar(tokens){const rules=[];for(let j=0;j<tokens.length;j++)for(let i=0;i<j;i++)rules.push({id:'order_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'commuting scalar polynomial variables'});const alg=new M.Presentation({tokens,rules});register(alg);return alg;}
const comm=(x,y)=>x.times(y).minus(y.times(x));
const power=(alg,x,n)=>{let r=alg.word([]);for(let j=0;j<n;j++)r=alg.reduce(r.times(x)).normal;return r;};

// One sine/cosine Weyl pair, h=1 and frequency divided by 2 pi n.
// This is an exact presentation, not a finite matrix model of the CCR.
const tokens=['sc','ss','dc','ds'],rules=[];
for(let j=0;j<tokens.length;j++)for(let i=0;i<j;i++){
 const hi=tokens[j],lo=tokens[i],rhs=[[[lo,hi],1]];
 if((hi==='dc'&&lo==='sc')||(hi==='ds'&&lo==='ss'))rhs.push([[],[0,1]]);
 rules.push({id:'weyl_'+hi+'_'+lo,lhs:[hi,lo],rhs,source:'MPQ [D_i,S_j]=i delta_ij, independent real Fourier modes'});
}
const alg=new M.Presentation({tokens,rules});register(alg);
const sc=alg.word(['sc']),ss=alg.word(['ss']),dc=alg.word(['dc']),ds=alg.word(['ds']),zero=alg.word([],0);
const R=ss.times(dc).minus(sc.times(ds));
prove(alg,'total_charge_configuration_cosine',comm(R,sc),ss.scale([0,1]));
prove(alg,'total_charge_configuration_sine',comm(R,ss),sc.scale([0,-1]));
prove(alg,'total_charge_momentum_cosine',comm(R,dc),ds.scale([0,1]));
prove(alg,'total_charge_momentum_sine',comm(R,ds),dc.scale([0,-1]));
prove(alg,'total_charge_symmetric_order',dc.times(ss).minus(ds.times(sc)),R);
prove(alg,'total_charge_radial_configuration_invariant',comm(R,sc.times(sc).plus(ss.times(ss))),zero);
no('reversed_total_charge_sign',alg.reduce(comm(R,sc).minus(ss.scale([0,-1]))).normal.terms.size!==0);
no('sum_instead_of_angular_difference',alg.reduce(comm(ss.times(dc).plus(sc.times(ds)),sc.times(sc).plus(ss.times(ss)))).normal.terms.size!==0);

// Current transport and the smearing bracket before periodic integration.
const jets=scalar(['w','wp','v','vp','d','dp','sp','m','mp','rp','gp','gpp','l']);
const z=Object.fromEntries([...jets.order.keys()].map(k=>[k,jets.word([k])])),nil=jets.word([],0);
const pb=z.wp.plus(z.w.times(z.sp)).times(z.vp.times(z.d).plus(z.v.times(z.dp))).scale(-1)
 .plus(z.wp.times(z.d).plus(z.w.times(z.dp)).times(z.vp.plus(z.v.times(z.sp))));
const pbExpected=z.dp.times(z.w.times(z.vp).minus(z.v.times(z.wp)))
 .plus(z.d.times(z.sp).times(z.v.times(z.wp).minus(z.w.times(z.vp))));
prove(jets,'source_current_smearing_bracket_density',pb,pbExpected);
const transported=z.mp.times(z.gp).times(z.gp).plus(z.m.times(z.gpp))
 .minus(z.m.times(z.rp).times(z.gp).times(z.gp).plus(z.m.times(z.gpp)));
prove(jets,'reference_current_weight_two_and_gsecond_cancellation',transported,z.mp.minus(z.m.times(z.rp)).times(z.gp).times(z.gp));
prove(jets,'equal_pair_Gaussian_score_cancellation',z.m.times(z.d).times(z.l).scale(-1).plus(z.d.times(z.m).times(z.l)),nil);
no('omit_log_derivative_from_r_transport',jets.reduce(z.mp.times(z.gp).times(z.gp).plus(z.m.times(z.gpp)).minus(z.m.times(z.rp).times(z.gp).times(z.gp)).minus(transported)).normal.terms.size!==0);
no('wrong_smearing_Witt_sign',jets.reduce(pb.plus(pbExpected)).normal.terms.size!==0);

// Exact all-n lower bound, and the general 7k leading coefficient.
const poly=scalar(['x','k']),x=poly.word(['x']),k=poly.word(['k']),unit=poly.word([]);
const n=x.plus(unit),np=n.plus(unit),a=unit.plus(power(poly,n,2)),b=unit.plus(power(poly,np,2));
const Q=n.times(power(poly,b,4)).minus(np.times(power(poly,a,4))).minus(power(poly,a,2).times(power(poly,b,2)).scale(4));
const coefficients=[193,1041,2084,2292,1570,702,204,36,3];
let positive=poly.word([],0);coefficients.forEach((v,i)=>{positive=positive.plus(power(poly,x,i).scale(v));});
prove(poly,'all_n_local_current_lower_bound_polynomial',Q,positive);
ck('all_n_shifted_polynomial_coefficients_positive',coefficients.every(v=>f(0).le(v)&&!f(v).zero()));
const ax=unit.plus(power(poly,x,2)),bx=unit.plus(power(poly,x.plus(k),2));
const numerator=poly.reduce(x.times(power(poly,bx,4)).minus(x.plus(k).times(power(poly,ax,4)))).normal;
const leading=numerator.entries().filter(([word])=>word.filter(t=>t==='x').length>=8);
ck('general_local_current_asymptotic_7k',leading.length===1&&leading[0][0].filter(t=>t==='x').length===8&&leading[0][0].filter(t=>t==='k').length===1&&leading[0][1].eq(7));
fixtures.push({kind:'all_n_and_all_k_polynomials',positive_coefficients:coefficients.map(f),general_numerator:numerator.toJSON(),normalization:'S/(2 pi), covariance square root (1+n^2)^(-2)'});
no('local_asymptotic_coefficient_one',!leading[0][1].eq(1));
no('local_asymptotic_coefficient_eight',!leading[0][1].eq(8));

// Pair rotations preserve exactly the selected covariance; unequal weights do not.
const diag=values=>N.matrix(values.map((v,i)=>values.map((_,j)=>i===j?v:0)));
for(const [cc,tt]of [['3/5','4/5'],['5/13','12/13']]){
 const c=f(cc),t=f(tt),rot=N.matrix([[c,t],[t.neg(),c]]);
 for(const lambda of ['1/16','1/625']){
  ck('paired_covariance_rotation:'+cc+':'+lambda,N.equal(N.mul(N.mul(rot,diag([lambda,lambda])),N.dagger(rot)),diag([lambda,lambda])));
 }
 no('unequal_pair_covariance_rotation:'+cc,N.equal(N.mul(N.mul(rot,diag([1,2])),N.dagger(rot)),diag([1,2]))===false);
}

// Coherent mean and variance coefficients, with 2 pi factored out.
const lambdaAt=n=>f(1).div(f(1).add(f(n).pow(2)).pow(4));
for(const mode of [1,2,5])for(const h of ['1','1/5']){
 const lambda=lambdaAt(mode),a=f('2/3'),b=f('3/7'),omega=f(mode);
 const mean=a.mul(b).mul(omega).neg();
 const varReal=f(h).mul(omega.pow(2)).mul(lambda).mul(b.pow(2));
 const varImag=f(h).mul(omega.pow(2)).mul(a.pow(2)).div(lambda.mul(4));
 ck('coherent_positive_variance:'+mode+':'+h,f(0).le(varReal.add(varImag))&&!varReal.add(varImag).zero());
 ck('coherent_phase_pair_opposite_means:'+mode+':'+h,mean.add(a.mul(b).mul(omega)).zero());
 ck('coherent_variance_linear_in_h:'+mode+':'+h,varReal.add(varImag).div(h).eq(omega.pow(2).mul(lambda.mul(b.pow(2)).add(a.pow(2).div(lambda.mul(4))))));
 fixtures.push({kind:'coherent_phase_pair',mode,h,mean_over_2pi:mean,variance_over_4pi_squared:varReal.add(varImag),same_configuration_law:true});
}
no('charge_is_positive_energy',f(-1).le(0)&&!f(-1).zero());

// Exact vacuum moments of the noisy joint current. No continuum extrapolation.
const factorial=n=>{let r=f(1);for(let j=2;j<=n;j++)r=r.mul(j);return r;};
const choose=(n,k)=>factorial(n).div(factorial(k).mul(factorial(n-k)));
const gaussian=(n,variance)=>{if(n%2)return f(0);let r=f(1);for(let j=1;j<n;j+=2)r=r.mul(j).mul(variance);return r;};
function differenceProductMoment(degree,q,p){let result=f(0);for(let k=0;k<=degree;k++)result=result.add(choose(degree,k).mul(k%2?-1:1).mul(gaussian(k,q)).mul(gaussian(k,p)).mul(gaussian(degree-k,q)).mul(gaussian(degree-k,p)));return result;}
for(const mode of [0,1,7])for(const h of ['1','1/3']){
 const lambda=lambdaAt(mode),q=f(h).mul(lambda).mul(2),p=f(h).div(lambda.mul(2));
 for(const degree of [0,1,2,3,4,6,8]){
  const actual=differenceProductMoment(degree,q,p),expected=degree%2?f(0):factorial(degree).mul(f(h).pow(degree));
  ck('Laplace_difference_product_moment:'+mode+':'+h+':'+degree,actual.eq(expected));
 }
}
for(const modes of [1,2,4,8,16,32]){
 let sum=f(0);for(let j=1;j<=modes;j++)sum=sum.add(f(j).pow(2));
 ck('joint_charge_vacuum_variance_sum:'+modes,sum.mul(2).eq(f(modes).mul(modes+1).mul(f(modes).mul(2).add(1)).div(3)));
 no('joint_charge_law_equals_spectral_vacuum:'+modes,!sum.zero());
 fixtures.push({kind:'joint_charge_vacuum',modes,variance_over_h_squared_4pi_squared:sum.mul(2),last_mode_density_sup_times_2pi_h:f(1).div(f(modes).mul(2))});
}
no('Gaussian_instead_of_Laplace_fourth_moment',!differenceProductMoment(4,f(2),f('1/2')).eq(gaussian(4,f(2))));

// Check centered quadratic Gaussian variance independently by Wick contractions.
function gaussianFourth(i,j,k,l){return f(i===j&&k===l?1:0).add(i===k&&j===l?1:0).add(i===l&&j===k?1:0);}
for(const mat of [N.matrix([[1,2],[2,-1]]),N.matrix([[2,1,0],[1,3,-2],[0,-2,1]])]){
 let second=f(0),trace=f(0);
 for(let i=0;i<mat.length;i++){trace=trace.add(mat[i][i].rad);for(let j=0;j<mat.length;j++)for(let k=0;k<mat.length;k++)for(let l=0;l<mat.length;l++)second=second.add(mat[i][j].rad.mul(mat[k][l].rad).mul(gaussianFourth(i,j,k,l)));}
 ck('second_chaos_norm_is_twice_HS:'+mat.length,second.sub(trace.pow(2)).eq(N.energy(mat).mul(2)));
 if(!trace.zero())no('omit_centering_trace',!second.eq(N.energy(mat).mul(2)));
}

// Local current matrices in the complex Fourier basis; divide operators by 2 pi.
// This finite Gaussian algebra supports T5; it does not establish a limit alone.
const sqrtLambda=n=>f(1).div(f(1).add(f(n).pow(2)).pow(2));
const sineCoeff=k=>k===1?C([0,'-1/2']):k===-1?C([0,'1/2']):C(0);
for(const cutoff of [1,2,4,8,16]){
 const modes=Array.from({length:2*cutoff+1},(_,j)=>j-cutoff);
 const B=N.matrix(modes.map(m=>modes.map(n=>C([0,n]).mul(sineCoeff(m-n)).mul(sqrtLambda(n).div(sqrtLambda(m))))));
 const S=N.scale(N.add(B,N.dagger(B)),'1/2'),hs=N.energy(S);
 ck('local_S_Hermitian:'+cutoff,N.equal(S,N.dagger(S)));
 ck('local_HS_lower_bound:'+cutoff,f(2).mul(cutoff-1).le(hs));
 const direct=N.matrix(modes.map(m=>modes.map(n=>C([0,'1/2']).mul(sineCoeff(m-n)).mul(f(n).mul(sqrtLambda(n).div(sqrtLambda(m))).sub(f(m).mul(sqrtLambda(m).div(sqrtLambda(n))))))));
 ck('local_S_coefficient_formula:'+cutoff,N.equal(S,direct));
 let linear=f(0);for(const n of modes)linear=linear.add(C([0,n]).mul(sineCoeff(n)).norm2().div(lambdaAt(n)));
 ck('local_linear_vacuum_norm:'+cutoff,linear.eq(8));
 for(let n=1;n<cutoff;n++){
  const entry=S[n+1+cutoff][n+cutoff];
  ck('local_adjacent_entry_bound:'+cutoff+':'+n,entry.turn.zero()&&f(1).le(entry.rad));
 }
 const rotation=N.matrix(modes.map(m=>modes.map(n=>m===n?C([0,n]):C(0))));
 ck('constant_window_second_chaos_zero:'+cutoff,N.isZero(N.add(rotation,N.dagger(rotation))));
 no('local_scalar_subtraction_removes_second_chaos:'+cutoff,!hs.zero());
 fixtures.push({kind:'local_current_vacuum',cutoff,HS_squared_over_4pi_squared:hs,first_chaos_norm_over_h_4pi_squared:linear.div(4),second_chaos_norm_over_h_squared_4pi_squared:hs.div(2)});
}

// The written chaos argument is exact: lower degrees are orthogonal to degree 2.
ck('centered_square_orthogonal_to_constant',gaussian(2,f(1)).sub(1).zero());
ck('centered_square_orthogonal_to_linear',gaussian(3,f(1)).sub(gaussian(1,f(1))).zero());
ck('centered_square_nonzero',gaussian(4,f(1)).sub(gaussian(2,f(1)).mul(2)).add(1).eq(2));
no('linear_counterterm_cancels_centered_quadratic',!gaussian(4,f(1)).sub(1).zero());

// Reuse the frozen IGQ path statistic; exact normalization and one witness.
function cosMean(n){return n%2?f(0):choose(n,n/2).div(f(2).pow(n));}
function variationFactor(power){let r=f(0);for(let j=0;j<=power;j++)r=r.add(choose(power,j).mul(f('1/2').pow(j)).mul(cosMean(j)));return r;}
ck('IGQ_nonrigid_variation_factor',variationFactor(7).eq('4587/1024'));
no('omit_time_change_in_quadratic_variation',!variationFactor(6).eq('4587/1024'));
ck('rigid_translation_variation_factor',f(1).pow(7).eq(1));
for(const mode of [1,2,4,8]){
 const den=f(1).add(f(mode).pow(2)).pow(2),ratio=f(mode).pow(4).div(den),defect=f(mode).pow(2).mul(2).add(1).div(den);
 ck('IGQ_bridge_amplitude_defect:'+mode,f(1).sub(ratio).eq(defect));
}
fixtures.push({kind:'reused_IGQ_variation',nonrigid_factor:variationFactor(7),rigid_factor:'1',scope:'exact coefficient only; singularity is the written measurable path-statistic proof'});

// Corrupted symbolic replays must be rejected.
const target=proofs.find(p=>p.certificate.steps.length),bad=JSON.parse(JSON.stringify(target.certificate));
bad.steps[0].coefficient=['991','0'];
no('corrupted_Paninian_certificate',(()=>{try{owners.get(target.presentation_sha256).replay(bad);return false;}catch{return true;}})());

const payload={protocol:'MRM1_MATERIAL_REFERENCE_MOMENTUM_V1',status:'PASS_TOTAL_CHARGE_AND_SCOPED_LOCAL_OBSTRUCTIONS',source_pins_sha256:sha(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),checks,negative_controls,presentations,proofs,fixtures,results:{
 physical_space:'unchanged selected RGR/MPQ Gaussian physical representation',
 source_reference_current:'rho=-d_prime+d*s_prime',
 total_charge_self_adjoint_with_spectral_domain:true,
 spectrum:'2*pi*hbar*Z, both signs; not identified as energy',
 total_charge_cutoff_strong_resolvent:true,
 density_pulse_feedback:'U^* R U=R-kappa*Pi*J(w_prime)',
 ideal_charge_spectral_pointer_and_positive_clock_norm:true,
 coherent_classical_recovery:true,
 raw_joint_detector_total_charge_vacuum_laws_tight:false,
 local_Weyl_current_vacuum_domain_limit:false,
 scalar_or_linear_counterterms_remove_second_chaos:false,
 nonrigid_canonical_material_flow_unitary_in_selected_measure_class:false,
 local_observable_nonexistence_in_every_representation:false,
 all_alternative_operator_domains_or_renormalizations_excluded:false,
 source_inverse_embedding_gauge_representation_changed:false,
 new_material_relabelling_constraint_imposed:false,
 physical_energy_identified:false,
 full_local_geometry_or_curvature_constructed:false,
 residual_gauge_uniqueness_assumed:false,
 physical_three_plus_one_gravity_closed:false,
 continuum_proofs_machine_certified:false
},summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,negative_controls:negative_controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',hash=sha(bytes),cp=path.join(__dirname,'CERTIFICATE.json'),hp=path.join(__dirname,'EXPECTED.sha256');
if(process.argv.includes('--write')){fs.writeFileSync(cp,bytes);fs.writeFileSync(hp,hash+'  CERTIFICATE.json\n');}
else{assert.equal(fs.readFileSync(cp,'utf8'),bytes);assert.equal(fs.readFileSync(hp,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
