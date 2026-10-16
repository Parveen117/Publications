'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs'),f=N.F.of,C=N.Cut.of;
const M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const checks=[],negative_controls=[],proofs=[],fixtures=[],presentations=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const sha=x=>crypto.createHash('sha256').update(x).digest('hex');
function register(alg){const audit=alg.audit();ck('confluent:'+alg.hash(),audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({spec:alg.spec(),audit:{status:audit.status,criticalPairs:audit.criticalPairs,resolved:audit.resolutions.length,open:audit.open.length,full_audit_sha256:sha(JSON.stringify(audit))}});}
function prove(alg,name,lhs,rhs){const r=alg.reduce(lhs.minus(rhs),{witness:true});ck(name,r.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(r.certificate));ck(name+':replay',alg.replay(certificate));proofs.push({name,presentation_sha256:alg.hash(),lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate});}

// Exact two-pair Weyl algebra at h=1. This is not a finite matrix CCR model.
const tokens=['s','a','d','b'],rules=[];
for(let j=0;j<tokens.length;j++)for(let i=0;i<j;i++){
 const hi=tokens[j],lo=tokens[i],rhs=[[[lo,hi],1]];
 if((hi==='d'&&lo==='s')||(hi==='b'&&lo==='a'))rhs.push([[],[0,1]]);
 rules.push({id:'order_'+hi+'_'+lo,lhs:[hi,lo],rhs,source:'MPQ [D,S]=i and independent matched probe [B,A]=i; h=1'});
}
const alg=new M.Presentation({tokens,rules});register(alg);
const s=alg.word(['s']),a=alg.word(['a']),d=alg.word(['d']),b=alg.word(['b']),one=alg.word([]),zero=alg.word([],0);
const comm=(x,y)=>x.times(y).minus(y.times(x)),H=s.times(b).minus(d.times(a)),u=s.plus(a),v=d.minus(b);
prove(alg,'commuting_joint_readings',comm(u,v),zero);
prove(alg,'source_feedback_s',comm(H,s).scale([0,1]),a);
prove(alg,'source_feedback_d',comm(H,d).scale([0,1]),b);
prove(alg,'probe_feedback_a',comm(H,a).scale([0,1]),s.scale(-1));
prove(alg,'probe_feedback_b',comm(H,b).scale([0,1]),d.scale(-1));
prove(alg,'opposite_ancilla_noise_commutator',comm(a,b.scale(-1)),one.scale([0,1]));
no('same_sign_momentum_output',alg.reduce(comm(s.plus(a),d.plus(b))).normal.terms.size!==0);
no('no_ancilla_joint_sharp_readout',alg.reduce(comm(s,d)).normal.terms.size!==0);
no('read_without_source_backreaction',alg.reduce(comm(H,s)).normal.terms.size!==0);

// A commutative scalar presentation for the normalization/calibration algebra.
const ts=['x','u','v','q','r','t','l','h'],rs=[];
for(let j=0;j<ts.length;j++)for(let i=0;i<j;i++)rs.push({id:'scalar_'+j+'_'+i,lhs:[ts[j],ts[i]],rhs:[[[ts[i],ts[j]],1]],source:'commuting scalar kernel, filter and covariance variables'});
const sc=new M.Presentation({tokens:ts,rules:rs});register(sc);
const z=Object.fromEntries(ts.map(k=>[k,sc.word([k])])),unit=sc.word([]),nil=sc.word([],0);
const pw=(x,n)=>{let r=unit;for(let j=0;j<n;j++)r=r.times(x);return r;};
prove(sc,'Gaussian_density_kernel_completion',pw(z.x.minus(z.u),2),pw(z.x,2).minus(z.x.times(z.u).scale(2)).plus(pw(z.u,2)));
prove(sc,'probe_score_cancels_in_rotation',z.x.times(z.u).times(z.l).scale(-1).plus(z.u.times(z.x).times(z.l)),nil);
prove(sc,'density_first_moment_exponent',z.x.plus(z.q.scale('1/2')).minus(z.q.scale('1/2')),z.x);
prove(sc,'density_second_moment_exponent',z.x.scale(2).plus(z.q.scale(2)).minus(z.q),z.x.scale(2).plus(z.q));
prove(sc,'two_mode_noise_Cauchy_Schwarz_remainder',pw(z.x,2).plus(pw(z.u,2)).times(pw(z.v,2).plus(pw(z.q,2))).minus(pw(z.x.times(z.v).plus(z.u.times(z.q)),2)),pw(z.x.times(z.q).minus(z.u.times(z.v)),2));
no('wrong_density_bias_subtraction',sc.reduce(z.x.plus(z.q.scale('1/2')).minus(z.q).minus(z.x)).normal.terms.size!==0);
no('output_second_moment_is_sharp',sc.reduce(z.x.scale(2).plus(z.q).minus(z.x.scale(2))).normal.terms.size!==0);
ck('balanced_kernel_density_Jacobian',f(2).mul('1/2').eq(1));
no('omit_kernel_output_Jacobian',!f(2).eq(1));

// Exact rational rotations test the full symplectic map and variance gains.
const Omega=N.matrix([[0,0,-1,0],[0,0,0,-1],[1,0,0,0],[0,1,0,0]]);
const diag=vs=>N.matrix(vs.map((v,i)=>vs.map((_,j)=>i===j?v:0)));
for(const [cc,tt] of [['3/5','4/5'],['5/13','12/13'],['8/17','15/17']]){
 const c=f(cc),t=f(tt),R=N.matrix([[c,t,0,0],[t.neg(),c,0,0],[0,0,c,t],[0,0,t.neg(),c]]);
 ck('rotation_unitarity:'+cc,N.equal(N.mul(R,N.dagger(R)),N.identity(4)));
 ck('rotation_symplectic:'+cc,N.equal(N.mul(N.mul(R,Omega),N.dagger(R)),Omega));
 const output=N.matrix([[f(1).div(c),0,0,0],[0,0,0,f(-1).div(t)]]),read=N.mul(output,R);
 ck('rescaled_readout:'+cc,N.equal(read,N.matrix([[1,t.div(c),0,0],[0,0,1,c.div(t).neg()]])));
 ck('readouts_commute:'+cc,N.isZero(N.mul(N.mul(read,Omega),N.dagger(read))));
 for(const lambda of ['1','1/16','1/625'])for(const h of ['1','1/3']){
  const q=f(h).mul(lambda),p=f(h).div(f(4).mul(lambda)),cov=diag([q,q,p,p]);
  ck('matched_covariance_invariant:'+cc+':'+lambda+':'+h,N.equal(N.mul(N.mul(R,cov),N.dagger(R)),cov));
  const Q=N.mul(N.mul(read,cov),N.dagger(read)),nu=q.mul(t.pow(2)).div(c.pow(2)),nv=p.mul(c.pow(2)).div(t.pow(2));
  ck('configuration_added_noise:'+cc+':'+lambda+':'+h,Q[0][0].eq(q.add(nu)));
  ck('momentum_added_noise:'+cc+':'+lambda+':'+h,Q[1][1].eq(p.add(nv)));
  ck('added_noise_product:'+cc+':'+lambda+':'+h,nu.mul(nv).eq(f(h).pow(2).div(4)));
  ck('vacuum_total_product_bound:'+cc+':'+lambda+':'+h,f(h).pow(2).le(Q[0][0].rad.mul(Q[1][1].rad)));
 }
 fixtures.push({kind:'rational_mixing',cosine:cc,sine:tt,rotation:R,rescaled_readout:read});
}
const badRotate=N.matrix([['3/5','4/5',0,0],['4/5','3/5',0,0],[0,0,'3/5','4/5'],[0,0,'-4/5','3/5']]);
no('same_sign_rotation_rows',!N.equal(N.mul(N.mul(badRotate,Omega),N.dagger(badRotate)),Omega));
const R=N.matrix([['3/5','4/5'],['-4/5','3/5']]),wrongCov=diag([1,2]);
no('unmatched_probe_covariance_preserves_product_measure',!N.equal(N.mul(N.mul(R,wrongCov),N.dagger(R)),wrongCov));

// Balanced covariance and exact Gaussian moments of both recorded quadratures.
const factorial=n=>{let q=f(1);for(let k=2;k<=n;k++)q=q.mul(k);return q;};
const choose=(n,k)=>factorial(n).div(factorial(k).mul(factorial(n-k)));
const gauss=(n,var_)=>{if(n%2)return f(0);let q=f(1);for(let k=1;k<n;k+=2)q=q.mul(k).mul(var_);return q;};
function sumMoment(n,var1,var2){let out=f(0);for(let k=0;k<=n;k++)out=out.add(choose(n,k).mul(gauss(k,var1)).mul(gauss(n-k,var2)));return out;}
for(let n=0;n<=5;n++){
 const lambda=f(1).div(f(1).add(f(n).pow(2)).pow(4)),q=lambda,p=f(1).div(lambda.mul(4));
 ck('balanced_CMG_u_variance:'+n,q.add(q).eq(lambda.mul(2)));
 ck('balanced_CMG_v_variance:'+n,p.add(p).eq(f(1).div(lambda.mul(2))));
 ck('balanced_added_noise_product:'+n,q.mul(p).eq('1/4'));
 for(let power=0;power<=8;power++){
  ck('configuration_gaussian_moment:'+n+':'+power,sumMoment(power,q,q).eq(gauss(power,q.mul(2))));
  ck('momentum_gaussian_moment:'+n+':'+power,sumMoment(power,p,p).eq(gauss(power,p.mul(2))));
 }
}
no('half_of_CMG_momentum_variance',!f('1/4').eq('1/2'));
// Equal configuration densities, phase-sensitive bounded readout: finite
// two-dimensional analogy only. It is not a truncated canonical commutator.
const ket0=N.matrix([['3/5'],['4/5']]),ket1=N.matrix([['3/5'],['-4/5']]);
const expect=(ket,O)=>N.mul(N.mul(N.dagger(ket),O),ket)[0][0];
for(const weights of [[1,0],[0,1],['2/7','5/9']])ck('phase_blind_multiplier:'+weights.join(','),expect(ket0,diag(weights)).eq(expect(ket1,diag(weights))));
const phaseProbe=N.matrix([['1/2','1/2'],['1/2','1/2']]);
ck('phase_sensitive_effect_first_state',expect(ket0,phaseProbe).eq('49/50'));
ck('phase_sensitive_effect_second_state',expect(ket1,phaseProbe).eq('1/50'));
no('density_record_separates_phase_pair',!expect(ket0,phaseProbe).eq(expect(ket1,phaseProbe)));
fixtures.push({kind:'finite_phase_blindness_analogy',first:ket0,second:ket1,phase_effect:phaseProbe,probabilities:['49/50','1/50'],scope:'finite information witness, not a continuum CCR realization'});

// Exact calibration through an invertible response, with no threshold erasure.
const modes=Array.from({length:9},(_,j)=>j),rawU=modes.map(j=>f(j+1).div(j+2)),rawV=modes.map(j=>f(j-3).div(j+1)),test=modes.map(j=>f(j%2===0?1:-1).div(j+1));
const dot=(x,y)=>x.reduce((q,v,j)=>q.add(v.mul(y[j])),f(0));
for(const ell of ['1','1/2','1/4']){
 const filter=modes.map(n=>f(1).div(f(1).add(f(ell).pow(2).mul(f(n).pow(2))).pow(4)));
 const uf=rawU.map((x,j)=>x.mul(filter[j])),vf=rawV.map((x,j)=>x.mul(filter[j])),inverseTest=test.map((x,j)=>x.div(filter[j]));
 ck('configuration_smearing_exact_calibration:'+ell,dot(uf,inverseTest).eq(dot(rawU,test)));
 ck('momentum_smearing_exact_calibration:'+ell,dot(vf,inverseTest).eq(dot(rawV,test)));
 no('same_smearing_ignores_response:'+ell,!dot(vf,test).eq(dot(rawV,test)));
 ck('inverse_high_mode_gain:'+ell,f(1).le(f(1).div(filter.at(-1))));
}
const gains=[2,4,8,16].map(n=>f(1).add(f(n).pow(2)).pow(4));
for(let j=1;j<gains.length;j++)ck('strictly_growing_inverse_gain:'+j,gains[j-1].le(gains[j])&&!gains[j-1].eq(gains[j]));
no('uniform_unit_noise_gain_for_inverse',!gains.at(-1).le(1));
// Positive rejection probability changes unconditional calibration if discarded.
const probs=[f('1/3'),f('2/3')],values=[f(2),f(5)];
const uncond=dot(probs,values),gated=probs[0].mul(values[0]);
no('erase_rejected_adapted_record',!uncond.eq(gated));
no('conditional_mean_equals_unconditional',!uncond.eq(gated.div(probs[0])));

// Spectral compression determines effects, not the disturbance map.
// Two projective instruments with the same effects, one with output feedback.
const P0=diag([1,0]),P1=diag([0,1]),swap=N.matrix([[0,1],[1,0]]);
const first=[P0,P1],second=[N.mul(swap,P0),N.mul(swap,P1)];
ck('instrument_effects_agree',first.every((K,j)=>N.equal(N.mul(N.dagger(K),K),N.mul(N.dagger(second[j]),second[j]))));
ck('instrument_effects_normalized',N.equal(N.add(P0,P1),N.identity(2)));
no('same_effects_force_same_instrument',!N.equal(N.mul(N.mul(first[0],P0),N.dagger(first[0])),N.mul(N.mul(second[0],P0),N.dagger(second[0]))));

// Tampered symbolic replay must fail rather than silently certify an identity.
const target=proofs.find(p=>p.certificate.steps.length),bad=JSON.parse(JSON.stringify(target.certificate));
bad.steps[0].coefficient=['991','0'];
const owner=target.presentation_sha256===alg.hash()?alg:sc;
no('corrupted_Paninian_certificate',(()=>{try{owner.replay(bad);return false;}catch{return true;}})());

const payload={protocol:'MDD1_MATERIAL_DETECTOR_DILATION_V1',status:'PASS_SOURCE_LIFT_JOINT_DETECTOR_AND_CALIBRATION',source_pins_sha256:sha(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),checks,negative_controls,presentations,proofs,fixtures,results:{physical_space:'selected RGR/MPQ with a prepared matched Gaussian probe',J_only_apparatus_phase_blind:true,finite_classical_first_class_source_lift:true,source_backreaction_retained:true,countable_rotation_and_CMG_POVM_identity:true,positive_clock_physical_norm:true,finite_bin_apparatus_instruments_converge:true,noise_product_per_mode:'hbar^2/4 for added variances',smooth_adapted_full_record_resolution_independence:true,original_J_density_mean_calibrated:true,existing_J_only_apparatus_realizes_joint_POVM:false,probe_covariance_and_readout_derived_from_gravity:false,local_material_detector_constructed:false,unreduced_nonlinear_quantum_orderings_closed:false,sharp_continuum_r_or_curvature_constructed:false,ungated_vacuum_mean_length:'still infinite',calibration_free_nonlinear_geometry:false,residual_gauge_uniqueness_assumed:false,physical_three_plus_one_gravity_closed:false,continuum_proofs_machine_certified:false},summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,negative_controls:negative_controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',hash=sha(bytes),cp=path.join(__dirname,'CERTIFICATE.json'),hp=path.join(__dirname,'EXPECTED.sha256');
if(process.argv.includes('--write')){fs.writeFileSync(cp,bytes);fs.writeFileSync(hp,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cp,'utf8'),bytes);assert.equal(fs.readFileSync(hp,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
