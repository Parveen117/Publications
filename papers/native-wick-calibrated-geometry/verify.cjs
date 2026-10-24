'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const W=require('./wick_geometry_adapter.cjs');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const P=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const {F,Cut}=N;
const argv=process.argv.slice(2);assert.ok(argv.length===1&&['--write','--check'].includes(argv[0]));
const checks=[],negative_controls=[],fixtures=[],presentations=[],proofs=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const rejects=f=>{try{f();return false;}catch{return true;}};
const sum=xs=>xs.reduce((s,x)=>s.add(x),new F(0));
const comm=(a,b)=>N.sub(N.mul(a,b),N.mul(b,a));
const det2=M=>M[0][0].mul(M[1][1]).sub(M[0][1].mul(M[1][0]));
function positive2(M){return N.equal(M,N.dagger(M))&&M[0][0].turn.zero()&&M[0][0].rad.n>0n&&det2(M).turn.zero()&&det2(M).rad.n>0n;}
function register(S,name){const audit=S.audit();ck('words:confluent:'+name,audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name,spec:S.spec(),audit});}
function prove(S,name,a,b){const r=S.reduce(a.minus(b),{witness:true});ck('words:'+name,r.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(r.certificate));ck('words:replay:'+name,S.replay(certificate));proofs.push({name,presentation_hash:S.hash(),lhs:a.toJSON(),rhs:b.toJSON(),certificate});}
function commutative(tokens){const rules=[];for(let j=1;j<tokens.length;j++)for(let i=0;i<j;i++)rules.push({id:'commute_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'commuting real scalar polynomial variables'});return new P.Presentation({tokens,rules});}
function power(S,x,n){let y=S.word([]);for(let k=0;k<n;k++)y=y.times(x);return y;}

// Finite Gaussian coefficient arithmetic. The continuum Wick limits and
// conditional-expectation theorem are analytic claims in THEOREM.md.
for(const h of ['1/8','1/4'])for(const cutoff of [0,1,4]){
  const z=W.calibration(h,cutoff,new F(cutoff-1,3),8),tag=h+':'+cutoff;
  ck('covariance:cutoff_includes_zero_mode:'+tag,z.c.eq(W.lambda(0).add(sum(Array.from({length:cutoff},(_,k)=>W.lambda(k+1).mul(2))))));
  for(let n=0;n<=8;n++){
    ck('wick:conditional_equals_sharp:degree_'+n+':'+tag,z.conditional[n].eq(z.sharp[n]));
    ck('wick:measured_vacuum_normalization:degree_'+n+':'+tag,z.vacuum[n].eq(n===0?1:0));
  }
  const halfWrong=W.wickCoefficient(2,0,z.nu.mul(2),z.nu.div(2));
  no('wick:half_measured_counterterm:'+tag,halfWrong.eq(z.nu.div(2))&&!halfWrong.zero());
  no('wick:half_conditional_counterterm:'+tag,!W.wickCoefficient(2,z.mean,z.nu,z.nu.div(2)).eq(z.sharp[2]));
  const omittedZero=W.wickCoefficient(2,0,z.nu.mul(2),z.h.mul(z.c.sub(1)));
  no('wick:omitted_zero_mode_counterterm:'+tag,omittedZero.eq(z.h)&&!omittedZero.zero());
  fixtures.push({kind:'finite_gaussian_wick_coefficients',...z,wrongHalfVacuumDegreeTwo:halfWrong,omittedZeroDegreeTwo:omittedZero});
}
ck('covariance:lambda_zero_is_one',W.lambda(0).eq(1));
ck('covariance:nonzero_LMD_modes',[1,2,5].every(n=>W.lambda(n).eq(new F(1,2*n))));
{
  const h=new F(1,8),c=W.cutoffVariance(3),k=new F(1,2),nu=h.mul(c),sharpCross=h.mul(k),measuredCross=sharpCross.mul(2),coefficients=[];
  for(let a=0;a<=4;a++)for(let b=0;b<=4;b++){
    const measured=W.jointWickCoefficient(a,b,nu.mul(2),measuredCross,nu),sharp=W.jointWickCoefficient(a,b,nu,sharpCross,nu.div(2));
    ck('wick:measured_cross_exponent:'+a+','+b,measured.eq(a===b?measuredCross.pow(a).div(W.factorial(a)):0));
    ck('wick:sharp_cross_exponent:'+a+','+b,sharp.eq(a===b?sharpCross.pow(a).div(W.factorial(a)):0));
    coefficients.push({a,b,measured,sharp});
  }
  no('wick:measured_and_sharp_cross_exponents_equal',!W.jointWickCoefficient(1,1,nu.mul(2),measuredCross,nu).eq(W.jointWickCoefficient(1,1,nu,sharpCross,nu.div(2))));
  fixtures.push({kind:'two_point_cross_covariance_polynomial',h,c,k,coefficients});
}
{
  const probabilities=[new F(1,3),new F(2,3)],values=[new F(2),new F(5)],mean=sum(values.map((x,i)=>x.mul(probabilities[i]))),phaseVariance=sum(values.map((x,i)=>x.sub(mean).pow(2).mul(probabilities[i])));
  ck('current:declared_phase_mass_variance_fixture',mean.eq(4)&&phaseVariance.eq(2));
  for(const h of ['1/8','1/4'])for(const n of [0,1,2,5]){
    const z=W.modeCovariance(h,n,phaseVariance),tag=h+':'+n;
    ck('measurement:configuration_covariance_doubled:'+tag,z.totalConfiguration.eq(z.h.mul(z.lambda).mul(2)));
    ck('measurement:added_current_noise:'+tag,z.addedCurrent.eq(z.h.div(z.lambda.mul(4))));
    ck('measurement:vacuum_current_doubled:'+tag,z.totalVacuumCurrent.eq(z.h.div(z.lambda.mul(2))));
    ck('measurement:phase_current_full_variance:'+tag,z.measuredPhaseCurrent.eq(z.h.div(z.lambda.mul(2)).add(phaseVariance)));
    ck('measurement:sharp_phase_current_variance:'+tag,z.sharpPhaseCurrent.eq(z.h.div(z.lambda.mul(4)).add(phaseVariance)));
    ck('measurement:added_noise_product:'+tag,z.addedConfiguration.mul(z.addedCurrent).eq(z.h.pow(2).div(4)));
    fixtures.push({kind:'mode_covariance',...z});
  }
}

// Actual compression of coherent Gaussian integrals on span{1,z}; these
// finite matrices are compressed moments, never a representation of the CCR.
for(const [h,r] of [['1/9','1/3'],['1/4','1/2']]){
  const z=W.zeroModeCompression(h,r),half=new F(1,2),zero=new F(0);
  const gx=n=>W.gaussianMoment(n,half),tilt=n=>W.shiftedGaussianMoment(n,z.r,half);
  const gram=N.matrix([[1,new Cut(gx(1),gx(1))],[new Cut(gx(1),gx(1).neg()),gx(2).add(gx(2))]]);
  ck('coherent:two_Hermite_analysis_Gram:'+h,N.equal(gram,z.gram));
  const q=N.matrix([[1,new Cut(tilt(1),gx(1))],[new Cut(tilt(1),gx(1).neg()),tilt(2).add(gx(2))]]);
  ck('coherent:actual_Wick_density_integrals:'+h,N.equal(q,z.wickDensity));
  // Multiplication by exp(S-h/2), S~N(0,h), on 1,S/sqrt(h).
  const sharp=N.matrix([[1,z.h.div(z.r)],[z.h.div(z.r),W.shiftedGaussianMoment(2,z.h,z.h).div(z.h)]]);
  ck('coherent:conditional_sharp_density_compression:'+h,N.equal(sharp,q));
  ck('coherent:positive_Wick_density_compression:'+h,positive2(q)&&det2(q).eq(1));
  const first=N.matrix([[zero,new Cut(0,z.vScale.mul(gx(2)))],[new Cut(0,z.vScale.mul(gx(2)).neg()),0]]);
  const second=N.matrix([[z.vScale.pow(2).mul(gx(2)),0],[0,z.vScale.pow(2).mul(gx(2).pow(2).add(gx(4)))]]);
  ck('coherent:actual_current_first_integral:'+h,N.equal(first,z.currentFirstMoment));
  ck('coherent:actual_current_second_integral:'+h,N.equal(second,z.currentSecondMoment));
  ck('coherent:sharp_current_second_oscillator_compression:'+h,N.equal(z.sharpCurrentSecondMoment,N.scale(N.matrix([[1,0],[0,3]]),z.h.div(4))));
  ck('coherent:current_noise_is_h_over_four:'+h,N.equal(N.sub(second,z.sharpCurrentSecondMoment),N.scale(N.identity(2),z.h.div(4))));
  no('coherent:compressed_density_and_current_commute:'+h,N.isZero(comm(q,first))===false);
  no('coherent:second_moment_equals_square_of_first:'+h,!N.equal(second,N.mul(first,first)));
  const configuration=N.matrix([[0,z.r],[z.r,0]]);
  no('coherent:finite_compression_satisfies_sharp_CCR:'+h,!N.equal(comm(first,configuration),N.scale(N.identity(2),N.IOTA.mul(z.h))));
  fixtures.push({kind:'actual_two_Hermite_coherent_compression',...z});
}

// A finite normalized measurement records simultaneous classical labels.
// Nonlinear geometry is evaluated only in the exact theta=0 chart. Rejected
// outcomes retain their probability and raw density/current labels.
{
  const A=N.matrix([['3/5',0],[0,'3/5'],['12/25','16/25'],['-16/25','12/25']]);
  ck('finite_geometry:analysis_isometry',N.equal(N.mul(N.dagger(A),A),N.identity(2)));
  const inputs=[{q:[2,4],current:[0,0]},{q:[3,1],current:[0,0]},{q:['1/8',2],current:[1,2]},{q:[2,'1/10'],current:[-2,1]}];
  const outcomes=inputs.map(x=>W.zeroAngleProfile(x.q,x.current));
  ck('finite_geometry:threshold_acceptance_pattern',outcomes.map(x=>x.accepted).join(',')==='true,true,false,false');
  ck('finite_geometry:actual_nonlinear_metric_labels',outcomes[0].boundedMetric[0].eq(new F(4,13))&&outcomes[1].boundedMetric[0].eq(new F(9,13)));
  ck('finite_geometry:accepted_positive_length',outcomes.slice(0,2).every(x=>x.length.eq(1)));
  const success=W.compressLabels(A,outcomes.map(x=>x.accepted?1:0)),rejected=W.compressLabels(A,outcomes.map(x=>x.accepted?0:1));
  ck('finite_geometry:all_probability_preserved',N.equal(N.add(success,rejected),N.identity(2)));
  ck('finite_geometry:explicit_success_probability',N.equal(success,N.scale(N.identity(2),'9/25')));
  ck('finite_geometry:rejection_mass_preserved',N.equal(rejected,N.scale(N.identity(2),'16/25')));
  ck('finite_geometry:rejected_raw_currents_retained',outcomes[2].current[0].eq(1)&&outcomes[3].current[0].eq(-2));
  const densityLabels=outcomes.map(x=>sum(x.q).div(x.q.length)),currentLabels=outcomes.map(x=>x.current[0]);
  const density=W.compressLabels(A,densityLabels),current=W.compressLabels(A,currentLabels);
  ck('finite_geometry:positive_density_first_moment',positive2(density));
  ck('finite_geometry:current_first_moment_Hermitian',N.equal(current,N.dagger(current)));
  no('finite_geometry:simultaneous_labels_force_commuting_first_moments',!N.isZero(comm(density,current)));
  const product=W.compressLabels(A,densityLabels.map((x,i)=>x.mul(currentLabels[i])));
  no('finite_geometry:compression_is_multiplicative',!N.equal(product,N.mul(density,current)));
  const boundedLabels=outcomes.map(x=>x.accepted?x.boundedMetric[0]:new F(0)),bounded=W.compressLabels(A,boundedLabels),boundedSecond=W.compressLabels(A,boundedLabels.map(x=>x.pow(2)));
  ck('finite_geometry:bounded_metric_effect_positive',positive2(bounded));
  ck('finite_geometry:bounded_metric_below_success_effect',positive2(N.sub(success,bounded)));
  const length=W.compressLabels(A,outcomes.map(x=>x.accepted?x.length:0));
  ck('finite_geometry:accepted_length_effect',N.equal(length,success));
  ck('finite_geometry:normalized_accepted_length_one',N.equal(N.scale(length,'25/9'),N.identity(2)));
  no('finite_geometry:nonlinear_metric_sharp_after_compression',!N.equal(boundedSecond,N.mul(bounded,bounded)));
  no('finite_geometry:discard_rejected_mass',!N.equal(success,N.identity(2)));
  fixtures.push({kind:'finite_label_geometry_dilation_not_sharp_joint_spectrum',analysis:A,outcomes,density,current,success,rejected,boundedMetric:bounded,acceptedLength:length,contract:'finite normalized measurement regression; literal coherent Gaussian compression is certified separately'});
}
{
  const q=[2,3,5],a=['1/2',1,'3/2'],expTheta=['1/2',3,'7/4'],ps=[1,-1].map(sign=>W.profileWeights(q,a,sign));
  for(const p of ps){const length=W.weightedLength(p,expTheta);ck('length:normalized_positive_weights:'+p.sign,sum(p.normalized).eq(1)&&p.normalized.every(x=>x.n>0n));ck('length:global_P_sign:'+p.sign,p.P.n*BigInt(p.sign)>0n);ck('length:positive_convex_bounds:'+p.sign,new F(1,2).le(length)&&length.le(3));ck('length:global_mass_range:'+p.sign,p.T.mul(new F(2,3)).le(p.P.abs())&&p.P.abs().le(p.T.mul(2)));fixtures.push({kind:'positive_weight_bound_comparison_not_full_inverse',profile:p,declaredExpTheta:expTheta,length,contract:'independent positive inputs check convex weights and global normalization only; only zeroAngleProfile evaluates the actual source inverse'});}
  ck('length:orientation_separate_geometry_equal',ps[0].pointFactors.every((x,i)=>x.eq(ps[1].pointFactors[i]))&&W.weightedLength(ps[0],expTheta).eq(W.weightedLength(ps[1],expTheta)));
  no('length:omit_global_normalization',!sum(ps[0].masses.map((x,i)=>x.abs().mul(ps[0].quadrature[i]))).eq(1));
  for(const a0 of ['1/2',1,'3/2'])for(const ratio of ['2/3',1,2]){const z=W.thresholdFactor(2,5,a0,F.of(ratio).mul(5));ck('length:threshold_bound:'+a0+':'+ratio,z.lower.le(z.actual)&&z.actual.le(z.upper));}
  ck('length:sharp_lower_threshold_endpoint',W.thresholdFactor(2,5,'3/2',10).actual.eq(new F(2,15)));
  ck('length:sharp_upper_threshold_endpoint',W.thresholdFactor(2,5,'1/2','10/3').actual.eq(new F(6,5)));
}

for(const ell of ['1/2','1/4','1/8','1/16']){
  const l=F.of(ell),h=l.pow(9),z=W.rateBudget({h,ell,densityTraceConstant:2,currentTraceConstant:3,tolerance:'1/2'});
  ck('rates:h_ell9_density:'+ell,z.density.eq(l.pow(2).mul(2)));
  ck('rates:h_ell9_current:'+ell,z.current.eq(l.mul(3)));
  ck('rates:joint_Markov_budget:'+ell,z.markov.eq(z.joint.mul(4)));
  const sr=W.sourceRateBudget({h,ell,densityKernelConstant:2,currentPhaseConstant:'3/5'});
  ck('rates:source_trace_constants:'+ell,sr.densityTraceConstant.eq(4170)&&sr.currentTraceConstant.eq(new F(10785,2)));
  ck('rates:source_current_vacuum_survives_zero_phase:'+ell,W.sourceRateBudget({h,ell,densityKernelConstant:2}).currentTraceConstant.eq(new F(8283,2)));
  fixtures.push({kind:'analytic_trace_bound_rational_budget',budget:z,sourceBudget:sr});
}
{
  const fine=W.rateBudget({h:'1/8',ell:'1/4',densityTraceConstant:1,currentTraceConstant:1}),coarse=W.rateBudget({h:'1/8',ell:'1/2',densityTraceConstant:1,currentTraceConstant:1});
  no('rates:fixed_h_decreases_current_error',coarse.current.le(fine.current)&&!coarse.current.eq(fine.current));
  no('rates:h_ell8_gives_vanishing_current',W.rateBudget({h:new F(1,4).pow(8),ell:'1/4',densityTraceConstant:1,currentTraceConstant:1}).current.eq(1));
}
ck('heat:factorial_majorant_constant_90',W.factorial(6).mul(8).div(new F(2).pow(6)).eq(90));
ck('heat:configuration_trace_constant_2085',new F(2048).add(36).add(1).eq(2085));
ck('heat:current_trace_constant_8283',new F(8192).add(90).add(1).eq(8283));
for(const [ell,L] of [['1/2',2],['1/4',4],['2/5',3],['1/3',3]]){
  const z=W.heatTraceSplit(ell,L),l=F.of(ell),tag=ell+':'+L;
  ck('heat:configuration_split_bound:'+tag,z.A.le(z.boundA));
  ck('heat:current_split_bound:'+tag,z.B.le(z.boundB));
  ck('heat:configuration_low_and_high:'+tag,z.lowA.le(new F(2048).div(l.pow(7)))&&z.highA.le(new F(36).div(l.pow(7))));
  ck('heat:current_low_and_high:'+tag,z.lowB.le(new F(8192).div(l.pow(8)))&&z.highB.le(new F(90).div(l.pow(8))));
  fixtures.push({kind:'exact_heat_trace_split_arithmetic',...z});
}
for(const ell of ['1/2','1/4'])for(const n of [1,2,5,10]){
  const z=W.heatMajorant(n,ell),den=F.of(ell).pow(12).mul(new F(n).pow(6));
  ck('heat:polynomial_weight_majorant:'+ell+':'+n,z.weight.le(z.weightBound));
  ck('heat:weighted_sixth_majorant:'+ell+':'+n,z.weightedMajorant.le(new F(90).div(den)));
}

// Symbolic conditional Wick coefficient identity with arbitrary commuting
// mean s and variance nu, proved from Gaussian pairings through degree eight.
{
  const S=commutative(['s','nu']);register(S,'conditional_Wick_polynomials');const s=S.word(['s']),nu=S.word(['nu']);
  const moment=n=>{let r=S.zero();for(let k=0;k<=n;k+=2){let pair=new F(1);for(let j=1;j<k;j+=2)pair=pair.mul(j);r=r.plus(power(S,s,n-k).times(power(S,nu,k/2)).scale(W.choose(n,k).mul(pair)));}return r;};
  const E=P.formalExp(s,8),G=P.formalExp(nu.scale('-1/2'),4),G2=Array.from({length:9},(_,n)=>n%2?S.zero():G[n/2]);
  const target=P.seriesMultiply(E,G2,8);
  for(let n=0;n<=8;n++){let actual=S.zero();for(let k=0;2*k<=n;k++)actual=actual.plus(power(S,nu.scale(-1),k).times(moment(n-2*k)).scale(new F(1).div(W.factorial(k).mul(W.factorial(n-2*k)))));prove(S,'conditional_Wick_degree_'+n,actual,target[n]);}
}
{
  const tokens=['h','s','a','d','b'],rules=[];
  for(let j=1;j<tokens.length;j++)for(let i=0;i<j;i++){const rhs=[[[tokens[i],tokens[j]],1]];if((tokens[j]==='d'&&tokens[i]==='s')||(tokens[j]==='b'&&tokens[i]==='a'))rhs.push([['h'],N.IOTA]);rules.push({id:'apparatus_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs,source:'two independent canonical pairs with matching commutators'});}
  const S=new P.Presentation({tokens,rules});register(S,'symbolic_commuting_apparatus_only');const [h,s,a,d,b]=tokens.map(t=>S.word([t]));
  prove(S,'original_CCR_nonzero',P.bracket(d,s),h.scale(N.IOTA));
  prove(S,'apparatus_CCR_matches',P.bracket(b,a),h.scale(N.IOTA));
  prove(S,'joint_apparatus_outputs_commute',P.bracket(d.minus(b),s.plus(a)),S.zero());
  no('words:wrong_apparatus_current_sign',S.reduce(P.bracket(d.plus(b),s.plus(a))).normal.terms.size>0);
  no('words:original_sharp_pair_commutes',S.reduce(P.bracket(d,s)).normal.terms.size>0);
}
{
  const S=commutative(['a','m','mp','y','yp','jp']);register(S,'source_inverse_cleared_relations');const [a,m,mp,y,yp,jp]=['a','m','mp','y','yp','jp'].map(t=>S.word([t]));
  const thetaPrime=a.times(yp),j=a.times(m),d=j.times(y),dPrime=jp.times(y).plus(j.times(yp)),mRPrime=m.times(thetaPrime).plus(mp),rho=mp.minus(mRPrime);
  prove(S,'source_density_equals_minus_j_yprime',rho,j.times(yp).scale(-1));
  prove(S,'source_current_inverse_without_division',rho.times(j),j.times(dPrime).scale(-1).plus(d.times(jp)));
  no('words:source_inverse_opposite_sign',S.reduce(rho.minus(j.times(yp))).normal.terms.size>0);
}
{
  const S=commutative(['x']);register(S,'positive_heat_weight_polynomial');const x=S.word(['x']),one=S.word([]);
  prove(S,'heat_weight_gap_positive_coefficients',power(S,one.plus(x),3).scale(8).minus(power(S,one.scale(2).plus(x),3)),x.scale(12).plus(power(S,x,2).scale(18)).plus(power(S,x,3).scale(7)));
}

no('domain:h_zero',rejects(()=>W.calibration(0,0)));
no('domain:h_half_measured_L2_boundary',rejects(()=>W.calibration('1/2',0)));
no('domain:negative_shifted_degree',rejects(()=>W.shiftedGaussianMoment(-1,0,1)));
no('domain:negative_joint_Wick_degree',rejects(()=>W.jointWickCoefficient(-1,2,1,0,1)));
no('domain:negative_variance',rejects(()=>W.gaussianMoment(2,-1)));
no('domain:nonpositive_joint_covariance',rejects(()=>W.jointGaussianMoment(2,2,1,2)));
no('domain:incorrect_exact_root',rejects(()=>W.zeroModeCompression('1/4','1/3')));
no('domain:unnormalized_analysis',rejects(()=>W.compressLabels([[2,0],[0,1]],[1,2])));
no('domain:negative_density',rejects(()=>W.profileWeights([1,-1],[1,1])));
no('domain:mixed_orientation_not_a_global_sign',rejects(()=>W.profileWeights([1,2],[1,1],[1,-1])));
no('domain:arbitrary_exact_Y_inverse_not_supplied',rejects(()=>W.zeroAngleProfile([2,3],[1,0])));
no('domain:outside_threshold_bounds',rejects(()=>W.thresholdFactor(1,1,2,1)));
no('domain:bad_heat_split',rejects(()=>W.heatTraceSplit('1/4',2)));
no('domain:negative_phase_variance',rejects(()=>W.modeCovariance('1/4',0,-1)));
no('domain:negative_source_current_constant',rejects(()=>W.sourceRateBudget({h:'1/4',ell:'1/2',densityKernelConstant:1,currentPhaseConstant:-1})));

const payload={protocol:'WCG1_WICK_CALIBRATED_GEOMETRY_V1',status:'PASS_WCG1_WICK_CALIBRATED_GEOMETRY',source_parent:'82b362b94f460cc4679bdc51ebc9cee5af9a4beb',checks,negative_controls,presentations,proofs,fixtures,
  scope:{Gaussian_coefficients:'finite exact Gaussian moments and formal coefficients through degree eight; no sampling or numerical exponentials',coherent_compression:'actual first-two-Hermite compression of positive Wick-density and current output moments',finite_geometry:'normalized finite dilation with nonlinear theta=0 geometry labels and explicit rejection probability',source_inverse:'polynomial certificates conditional on analytic derivative identities; no arbitrary Y inverse or cosine evaluator',currents:'symbolic canonical apparatus algebra and actual compressed moments; no finite CCR representation',rates:'rational consequences of separately proved analytic heat and GMC trace estimates',configuration_trace_constant:2085,current_trace_constant:8283,measured_Wick_L2_range:'0<h<1/2',finite_fixtures_prove_continuum_Wick_limit:false,finite_fixtures_prove_current_domains:false,finite_fixtures_prove_heat_trace_bounds:false,finite_joint_labels_are_sharp_joint_spectrum:false,spacetime_or_quantum_gravity_claim:false},
  summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,negative_controls:negative_controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',hash=crypto.createHash('sha256').update(bytes).digest('hex'),cp=path.join(__dirname,'CERTIFICATE.json'),hp=path.join(__dirname,'EXPECTED.sha256');
if(argv[0]==='--write'){fs.writeFileSync(cp,bytes);fs.writeFileSync(hp,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cp,'utf8'),bytes);assert.equal(fs.readFileSync(hp,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
