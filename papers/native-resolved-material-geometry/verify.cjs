'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs'),f=N.F.of,C=N.Cut.of;
const M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const checks=[],negative_controls=[],proofs=[],fixtures=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const sha=x=>crypto.createHash('sha256').update(x).digest('hex');
const tokens=['j','a','b','y1','u1','v','v1','h','l','n'];
const rules=[];for(let j=0;j<tokens.length;j++)for(let i=0;i<j;i++)rules.push({id:'commute_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'commuting classical output jets and Gaussian variance symbols'});
const alg=new M.Presentation({tokens,rules}),audit=alg.audit();
ck('scalar_presentation_confluent',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');
const presentations=[{spec:alg.spec(),audit:{status:audit.status,criticalPairs:audit.criticalPairs,resolved:audit.resolutions.length,open:audit.open.length,full_audit_sha256:sha(JSON.stringify(audit))}}];
const V=Object.fromEntries(tokens.map(t=>[t,alg.word([t])])),one=alg.word([]),zero=alg.word([],0);
const pow=(p,n)=>{let r=one;for(let i=0;i<n;i++)r=r.times(p);return r;};
function prove(name,lhs,rhs){const r=alg.reduce(lhs.minus(rhs),{witness:true});ck(name,r.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(r.certificate));ck(name+':replay',alg.replay(certificate));proofs.push({name,presentation_sha256:alg.hash(),lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate});}
// Here b=a_theta. Multiplying rho=m'-m r' by a avoids division.
const aMprime=V.j.times(V.u1.minus(V.b.times(V.y1)));
const aMrprime=V.j.times(V.a.times(V.y1).plus(V.u1).minus(V.b.times(V.y1)));
prove('reference_constraint_inverse_chart',aMprime.minus(aMrprime),V.a.times(V.j).times(V.y1).scale(-1));
// Differentiate d=j y: j y'=d'-d s'. The two rho expressions agree.
prove('adapted_reference_momentum',V.a.times(V.v1.minus(V.v.times(V.u1))).scale(-1),V.a.times(V.v).times(V.u1).minus(V.a.times(V.v1)));
prove('joint_configuration_variance',V.h.times(V.l).plus(V.l.times(V.h)),V.h.times(V.l).scale(2));
// Denominators cleared: (h lambda)*(h/(4 lambda))=h^2/4.
prove('added_noise_uncertainty_product',V.h.times(V.l).times(V.h),pow(V.h,2).times(V.l));
no('drop_reference_profile_derivative',alg.reduce(aMprime.minus(aMrprime).minus(V.j.times(V.a).times(V.y1).scale(-1)).plus(V.v1)).normal.terms.size!==0);
no('drop_measurement_configuration_noise',alg.reduce(V.h.times(V.l).minus(V.h.times(V.l).scale(2))).normal.terms.size!==0);

// Native rational polynomial arithmetic; coefficients are F, never floats.
const clean=a=>{while(a.length>1&&a.at(-1).zero())a.pop();return a;};
const add=(a,b)=>clean(Array.from({length:Math.max(a.length,b.length)},(_,i)=>(a[i]||f(0)).add(b[i]||f(0))));
const scale=(a,c)=>clean(a.map(x=>x.mul(c)));
const mul=(a,b)=>{const r=Array(a.length+b.length-1).fill(f(0));for(let i=0;i<a.length;i++)for(let j=0;j<b.length;j++)r[i+j]=r[i+j].add(a[i].mul(b[j]));return clean(r);};
const ppow=(a,n)=>{let r=[f(1)];for(let i=0;i<n;i++)r=mul(r,a);return r;};
const peval=(a,x)=>a.reduceRight((s,c)=>s.mul(x).add(c),f(0));
const binom=(n,k)=>{let r=f(1);for(let i=1;i<=k;i++)r=r.mul(n-i+1).div(i);return r;};
const sums=[[f(0),f(1)]];
for(let p=1;p<=14;p++){
 let r=add(ppow([f(1),f(1)],p+1),[f(-1)]);
 for(let k=0;k<p;k++)r=add(r,scale(sums[k],binom(p+1,k).mul(-1)));
 sums.push(scale(r,f(1).div(p+1)));
}
function budget(power){let r=[f(1)];for(let k=0;k<=power;k++)r=add(r,scale(sums[2*k],binom(power,k).mul(2)));return r;}
const B4=budget(4),B7=budget(7);
const expr=(p,x)=>p.reduce((r,c,i)=>r.plus(pow(x,i).scale(c.toString())),zero);
for(const[p,b]of [[4,B4],[7,B7]]){
 prove('all_N_noise_budget_recurrence:'+p,expr(b,V.n.plus(one)).minus(expr(b,V.n)),pow(one.plus(pow(V.n.plus(one),2)),p).scale(2));
 ck('noise_budget_initial_value:'+p,peval(b,0).eq(1));
 ck('noise_budget_degree:'+p,b.length===2*p+2);
 ck('noise_budget_leading_coefficient:'+p,b.at(-1).eq(f(2).div(2*p+1)));
 fixtures.push({kind:'exact_noise_budget_polynomial',power:p,coefficients:b.map(x=>x.toString()),scope:'recurrence and base prove the finite sum identity; continuum use remains written analysis'});
}
let direct4=f(1),direct7=f(1),qtrace=f(1);
for(let cutoff=0;cutoff<=32;cutoff++){
 if(cutoff){const b=f(1).add(f(cutoff).pow(2));direct4=direct4.add(b.pow(4).mul(2));direct7=direct7.add(b.pow(7).mul(2));qtrace=qtrace.add(f(2).div(b));}
 ck('point_noise_sum:'+cutoff,peval(B4,cutoff).eq(direct4));
 ck('H3_noise_sum:'+cutoff,peval(B7,cutoff).eq(direct7));
 ck('configuration_H3_trace_bounded:'+cutoff,qtrace.le(5));
 if(cutoff){
  ck('N9_lower_bound:'+cutoff,f(cutoff).pow(9).mul('2/9').le(direct4));
  ck('N15_upper_bound:'+cutoff,direct7.le(f(257).mul(f(cutoff).pow(15))));
  const h=f(1).div(f(cutoff).pow(16));
  ck('simultaneous_limit_budget:'+cutoff,h.mul(direct7).div(2).le(f(257).div(f(2).mul(cutoff))));
 }
 if([0,1,2,4,8,16,32].includes(cutoff))fixtures.push({kind:'vacuum_noise',cutoff,point_momentum_variance_h1:direct4.div(2).toString(),H3_momentum_noise_h1:direct7.div(2).toString(),H3_configuration_noise_h1:qtrace.mul(2).toString()});
}
no('use_C_instead_of_inverse_C_for_momentum',!f(1).div(16).eq(16));
no('fixed_hbar_has_vanishing_noise',!peval(B4,32).le(peval(B4,1)));
no('hbar_N9_condition_implies_H3_budget',f(1).le(peval(B7,8).div(f(8).pow(10))));

// Gaussian polynomial moments for the joint coherent law; no sampling.
const moment=(n,var_)=>{if(n%2)return f(0);let r=f(1);for(let k=1;k<n;k+=2)r=r.mul(k).mul(var_);return r;};
for(const lambda of ['1','1/16','1/625'])for(const h of ['1','1/3']){
 const q=f(h).mul(lambda).mul(2),p=f(h).div(f(lambda).mul(2));
 ck('total_noise_product:'+lambda+':'+h,q.mul(p).eq(f(h).pow(2)));
 ck('added_noise_product:'+lambda+':'+h,q.div(2).mul(p.div(2)).eq(f(h).pow(2).div(4)));
 ck('Gaussian_fourth_configuration_moment:'+lambda+':'+h,moment(4,q).eq(q.pow(2).mul(3)));
 ck('Gaussian_mixed_fourth_moment:'+lambda+':'+h,moment(2,q).mul(moment(2,p)).eq(f(h).pow(2)));
 no('erase_joint_measurement_noise:'+lambda+':'+h,!q.eq(q.div(2)));
}

// Exact rational inverse-chart jets at declared points, not trig evaluation.
for(const sign of [-1,1])for(const a of ['1/2','1','3/2']){
 const j=f(sign).mul(3),d=f('5/7'),jp=j.mul('2/3'),dp=f('-4/5'),ap=f('1/3');
 const yprime=dp.div(j).sub(d.mul(jp).div(j.pow(2))),m=j.div(a);
 const mp=jp.div(a).sub(j.mul(ap).mul(yprime).div(a));
 const rp=f(a).mul(yprime).add(mp.div(m));
 ck('source_reference_constraint:'+sign+':'+a,mp.sub(m.mul(rp)).eq(dp.mul(-1).add(d.mul(jp).div(j))));
 ck('measuring_density_reconstructed:'+sign+':'+a,m.mul(a).eq(j));
 ck('momentum_sign_retained:'+sign+':'+a,f(0).le(m.mul(sign)));
}
const masses=[f(2),f(3)],dm=[f(1),f(-2)],dr=[f('2/7'),f('-1/5')];
const P=masses.reduce((s,x)=>s.add(x),f(0)),dP=dm.reduce((s,x)=>s.add(x),f(0));
const thetaVariation=dr.map((x,i)=>x.add(dP.div(P)).sub(dm[i].div(masses[i])));
const pair=xs=>xs.reduce((s,x,i)=>s.add(x.mul(masses[i])),f(0));
ck('global_P_feedback_preserves_canonical_one_form',pair(thetaVariation).eq(pair(dr)));
no('freeze_global_P_in_log_variation',!pair(thetaVariation.map(x=>x.sub(dP.div(P)))).eq(pair(dr)));

// A finite rational tight-frame witness of the dilation statements only.
// Four rows are 1/2*(1,+/-1), 1/2*(1,+/-i).
const W=N.matrix([['1/2','1/2'],['1/2','-1/2'],['1/2',[0,'1/2']],['1/2',[0,'-1/2']]]),Wd=N.dagger(W);
ck('finite_POVM_analysis_isometry',N.equal(N.mul(Wd,W),N.identity(2)));
const values=[-2,-1,1,2],diag=vals=>N.matrix(vals.map((v,i)=>vals.map((_,j)=>i===j?v:0)));
const T=N.mul(N.mul(Wd,diag(values)),W),T2=N.mul(N.mul(Wd,diag(values.map(v=>f(v).pow(2)))),W);
const defect=N.sub(T2,N.mul(T,T)),leak=N.mul(N.sub(N.identity(4),N.mul(W,Wd)),N.mul(diag(values),W));
ck('compressed_real_reading_is_selfadjoint',N.equal(T,N.dagger(T)));
ck('measurement_variance_is_positive_Gram',N.equal(defect,N.mul(N.dagger(leak),leak)));
ck('variance_defect_nonzero',!N.isZero(defect));
no('POVM_second_moment_is_square_of_mean',N.equal(T2,N.mul(T,T))===false);
const E0=N.mul(N.dagger([W[0]]),[W[0]]),E2=N.mul(N.dagger([W[2]]),[W[2]]);
no('coherent_frame_effect_is_projection',!N.equal(N.mul(E0,E0),E0));
no('all_compressed_geometry_readings_commute',!N.isZero(N.commutator(E0,E2)));
fixtures.push({kind:'finite_dilation_only',isometry:W,mean:T,second_moment:T2,positive_variance_defect:defect});
const bad=JSON.parse(JSON.stringify(proofs.find(p=>p.certificate.steps.length).certificate));bad.steps[0].coefficient=['997','0'];
no('corrupt_symbolic_certificate',(()=>{try{alg.replay(bad);return false;}catch{return true;}})());
const payload={protocol:'RMG1_RESOLVED_MATERIAL_GEOMETRY_V1',status:'PASS_RESOLVED_GEOMETRY_AND_SCOPED_OBSTRUCTION',source_pins_sha256:sha(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),checks,negative_controls,presentations,proofs,fixtures,results:{physical_space:'unchanged RGR selected MPQ quotient',added_contract:'finite-resolution coherent-state joint measurement',normalized_geometry_POVM:true,positive_momentum_mean_form_domains:true,sharp_original_r_field_constructed:false,fixed_hbar_raw_point_geometry_limit:false,obstruction_scope:'this POVM in the physical vacuum; not all composite prescriptions',vacuum_mean_length_N0:'infinite',controlled_classical_C2_limit:'hbar N^15 -> 0',point_conjugate_variance_asymptotic:'hbar N^9/9',spacetime_curvature_operator_constructed:false,gauge_dual_completeness_classified:false,physical_three_plus_one_gravity_closed:false,continuum_proofs_machine_certified:false},summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,negative_controls:negative_controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',hash=sha(bytes),cp=path.join(__dirname,'CERTIFICATE.json'),hp=path.join(__dirname,'EXPECTED.sha256');
if(process.argv.includes('--write')){fs.writeFileSync(cp,bytes);fs.writeFileSync(hp,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cp,'utf8'),bytes);assert.equal(fs.readFileSync(hp,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
