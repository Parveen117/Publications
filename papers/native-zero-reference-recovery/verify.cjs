'use strict';
// Native rational arithmetic and proof replay; no numerical fitting or quadrature.
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const f=N.F.of,mm=(...a)=>a.reduce(N.mul),sc=N.scale,eq=N.equal,zero=N.isZero,cm=N.commutator;
const checks=[],negative_controls=[],proofs=[],presentations=[],systems=new Map(),results={};
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const real=z=>{assert.ok(z.turn.zero());return z.rad;},tr=a=>a.reduce((s,row,i)=>s.add(row[i]),N.ZERO);
const I=N.identity(2),X=N.matrix([[0,1],[1,0]]),Z=N.matrix([[1,0],[0,-1]]),P=N.matrix([[1,0],[0,0]]),Q=N.sub(I,P),V=N.matrix([[1,0],[0,2]]);
const frames=[I,N.matrix([['3/5','-4/5'],['4/5','3/5']]),N.matrix([[['3/5','0'],['0','4/5']],[['0','4/5'],['3/5','0']]])];
const hams=[X,Z,N.matrix([[2,['3/5','4/5']],[['3/5','-4/5'],-1]])];
for(let i=0;i<frames.length;i++)for(let j=0;j<hams.length;j++){
 const u=frames[i],p=mm(u,P,N.dagger(u)),q=N.sub(I,p),h=hams[j],b=mm(q,h,p),a=mm(p,h,p),defect=N.sub(mm(p,h,h,p),mm(a,a));
 ck('closure:frame_unitarity:'+i+':'+j,eq(mm(u,N.dagger(u)),I));
 ck('closure:dagger_square:'+i+':'+j,eq(defect,mm(N.dagger(b),b)));
 ck('closure:commutator_criterion:'+i+':'+j,zero(b)===zero(cm(h,p)));
 const shifted=N.add(h,sc(I,'17/3'));
 ck('zero:energy_shift:'+i+':'+j,eq(mm(q,shifted,p),b));
 const defect2=N.sub(mm(p,shifted,shifted,p),mm(p,shifted,p,shifted,p));
 ck('zero:second_moment_shift:'+i+':'+j,eq(defect2,defect));
}
ck('zero:scalar_reading_shift',eq(mm(N.sub(V,I),P),mm(N.sub(N.sub(V,sc(I,'7/3')),sc(I,'-4/3')),P)));
ck('stationary:density',zero(cm(Z,P)));
ck('stationary:quantum_variance',real(tr(mm(P,X,X))).sub(real(tr(mm(P,X))).pow(2)).eq(1));
no('stationary:implies_commutative_observables',!zero(cm(X,Z)));
ck('sharp:initial_lapse',eq(mm(V,P),P));
no('sharp:implies_generator_closure',!zero(mm(Q,X,P)));
ck('closed:scalar_sector',zero(cm(Z,P))&&eq(mm(V,P),P));
no('closed:alone_implies_scalar_metric',zero(cm(Z,I))&&!eq(V,sc(I,'3/2')));
const u=N.add(sc(I,'3/5'),sc(X,N.IOTA.mul('-4/5'))),compressed=mm(P,u,P),survival=real(tr(mm(P,N.dagger(u),P,u,P))),loss=real(tr(mm(P,N.dagger(u),Q,u,P)));
ck('return:exact_unitarity',eq(mm(N.dagger(u),u),I));
ck('return:retained_amplitude',eq(compressed,sc(P,'3/5')));
ck('return:survival',survival.eq('9/25'));ck('return:loss',loss.eq('16/25'));
ck('return:probability_conservation',survival.add(loss).eq(1));
no('return:initial_sharpness_removes_memory',!eq(compressed,P));
results.finite={phase_cosine:'3/5',phase_sine:'4/5',survival,loss,initial_lapse_variance:'0',stationary_X_variance:'1',scope:'Exact unitary phase fixture; the finite-time inequalities have written proofs.'};

// Finite Laurent coefficient vectors represent exact circle polynomials.
// Shifts extend the support instead of imposing artificial finite boundaries.
function choose(n,k){if(k<0||k>n)return f(0);let r=f(1);for(let j=1;j<=k;j++)r=r.mul(n-k+j).div(j);return r;}
const add=(a,b)=>{const c=new Map(a);for(const[n,z]of b)c.set(n,(c.get(n)||N.ZERO).add(z));return c;};
const scale=(a,z)=>new Map([...a].map(([n,x])=>[n,x.mul(z)]));
const shift=(a,k)=>new Map([...a].map(([n,x])=>[n+k,x]));
const dot=(a,b)=>[...a].reduce((s,[n,z])=>s.add(z.dagger().mul(b.get(n)||N.ZERO)),N.ZERO);
const rOp=a=>add(scale(a,'1/2'),scale(add(shift(a,1),shift(a,-1)),'-1/4'));
const sinOp=a=>scale(add(shift(a,1),scale(shift(a,-1),-1)),N.IOTA.mul('-1/2'));
const number=a=>new Map([...a].map(([n,z])=>[n,z.mul(n)]));
const Hg=(a,eta)=>add(add(scale(number(number(a)),eta.pow(2)),scale(a,8)),scale(add(shift(a,2),shift(a,-2)),-1));
const packets=[];
for(const m of [1,2,3,4,8,12]){
 const a=new Map(Array.from({length:2*m+1},(_,i)=>[i-m,N.Cut.of(choose(2*m,i))])),norm=real(dot(a,a)),eta=f(1).div(m);
 const ex=op=>real(dot(a,op(a))).div(norm),q=ex(rOp),r2=ex(x=>rOp(rOp(x))),l2=ex(x=>number(number(x))),sin2=ex(x=>sinOp(sinOp(x))),epsilon=ex(x=>Hg(x,eta)).sub(6),variance=r2.sub(q.pow(2));
 const expectedE=eta.pow(2).mul(f(m).pow(2)).div(4*m-1).add(f(2*(4*m+1)).div((m+1)*(2*m+1)));
 ck('packet:normalization:'+m,norm.eq(choose(4*m,2*m)));
 ck('packet:first_moment:'+m,q.eq(f(1).div(4*m+2)));
 ck('packet:second_moment:'+m,r2.eq(f(3).div((4*m+2)*(4*m+4))));
 ck('packet:number_square:'+m,l2.eq(f(m).pow(2).div(4*m-1)));
 ck('packet:sine_square:'+m,sin2.eq(q.sub(r2).mul(4)));
 ck('packet:energy:'+m,epsilon.eq(expectedE));
 ck('packet:variance_positive:'+m,f(0).le(variance)&&!variance.zero());
 ck('packet:energy_budget:'+m,epsilon.le(f(13).div(3*m)));
 ck('packet:initial_rate_zero:'+m,real(dot(a,add(Hg(rOp(a),eta),scale(rOp(Hg(a,eta)),-1))).mul(N.IOTA)).zero());
 const translated=new Map([...a].map(([n,z])=>[n,z.mul(n%2===0?1:-1)]));
 const qOther=real(dot(translated,rOp(translated))).div(norm),eOther=real(dot(translated,Hg(translated,eta))).div(norm).sub(6);
 ck('well:complement_moment:'+m,qOther.eq(f(1).sub(q)));
 ck('well:same_energy:'+m,eOther.eq(epsilon));
 no('well:small_energy_selects_reference:'+m,!qOther.eq(q));
 packets.push({m,eta,normalization_squared:norm,q,r2,number_square:l2,epsilon,variance,other_well_q:qOther});
}
results.packets=packets;
const p2=packets.find(x=>x.m===2),af=f('3/10');
ck('readout:variable_profile_mean_at_x1',f(1).add(af.mul(p2.q)).eq('103/100'));
ck('readout:variable_profile_variance_at_x1',af.pow(2).mul(p2.variance).eq('27/20000'));
ck('readout:mean_curvature_at_origin',f('-12/5').mul(p2.q).eq('-6/25'));
no('readout:finite_width_exactly_flat',!p2.q.zero());
const kinetic1=f(1).div(3),kinetic12=f(144).div(47);
no('width:narrowing_at_fixed_eta_has_zero_cost',!kinetic12.le(kinetic1));

const budgets=[];
for(const n of [1,10,100,1000]){
 const m=n*n,eta=f(1).div(m),q=f(1).div(4*m+2),eps=eta.pow(2).mul(f(m).pow(2)).div(4*m-1).add(f(2*(4*m+1)).div(f(m+1).mul(2*m+1)));
 const qRoot=f(1).div(2*n),eRoot=f(21).div(10*n),a=f('3/5'),K=f('3/2'),error=a.mul(K).mul(qRoot.add(eRoot.div(2)));
 ck('budget:q_square:'+n,q.le(qRoot.pow(2)));ck('budget:energy_square:'+n,eps.le(eRoot.pow(2)));
 ck('budget:energy_majorant:'+n,eps.le(f(13).div(3*m)));
 ck('budget:state_error:'+n,error.eq(f(279).div(200*n)));
 ck('budget:observable_error:'+n,error.mul(2).eq(f(279).div(100*n)));
 budgets.push({n,m,eta,q,epsilon:eps,q_root_upper:qRoot,energy_root_upper:eRoot,state_error:error,observable_error:error.mul(2)});
}
results.budgets={a:'3/5',g:'0',s_max:'1',matter_derivative_norm_max:'1',K_max:'3/2',cases:budgets,scope:'Analytic family; large-m moments and bounds evaluated exactly without constructing large matrices.'};

function register(name,spec){const s=new M.Presentation(spec),audit=s.audit();ck(name+':confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name,spec:s.spec(),audit});systems.set(s.hash(),s);return{s,w:Object.fromEntries(spec.tokens.map(t=>{const n=typeof t==='string'?t:t.name;return[n,s.word([n])];}))};}
function prove(name,lhs,rhs){const s=lhs.system,r=s.reduce(lhs.minus(rhs),{witness:true});ck(name,r.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(r.certificate));ck(name+':replay',s.replay(certificate));proofs.push({name,lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate});}
{
 const{s,w:{p,h,v,z}}=register('observation_projection',{tokens:['p','h','v','z'],rules:[
 {id:'cut',lhs:['p','p'],rhs:[[['p'],1]],source:'Retained cut idempotence'},
 ...['p','h','v'].map(t=>({id:'central_'+t,lhs:['z',t],rhs:[[[t,'z'],1]],source:'Scalar reference shift'}))]});
 const one=s.one(),q=one.minus(p),b=q.times(h).times(p),a=p.times(h).times(p);
 prove('closure:second_moment_defect',p.times(h).times(h).times(p).minus(a.times(a)),p.times(h).times(q).times(h).times(p));
 prove('closure:complement_orthogonality',q.times(p),s.zero());
 prove('zero:energy_origin_invariance',q.times(h.plus(z)).times(p),b);
 prove('zero:geometry_reference_invariance',v.minus(z).minus(one.minus(z)).times(p),v.minus(one).times(p));
 prove('closure:commutator_blocks',h.times(p).minus(p.times(h)),b.minus(p.times(h).times(q)));
}
{
 const source=JSON.parse(fs.readFileSync(path.join(__dirname,'../native-action-geometry-dynamics/CERTIFICATE.json'),'utf8'));
 const spec=source.presentations.find(x=>x.name==='circle_momentum').spec;
 const{s,w:{u,v,h,p}}=register('inherited_circle_momentum',spec),one=s.one(),r=one.scale('1/2').minus(u.plus(v).scale('1/4')),sin=u.minus(v).scale(N.IOTA.mul('-1/2')),rp=sin.scale('1/2'),hg=p.times(p).plus(one.scale(8)).minus(u.times(u)).minus(v.times(v));
 prove('circle:positive_localization_identity',rp.times(rp),r.times(one.minus(r)));
 prove('circle:potential_two_wells',hg.minus(p.times(p)).minus(one.scale(6)),r.times(one.minus(r)).scale(16));
 prove('circle:momentum_reading_commutator',p.times(r).minus(r.times(p)),h.times(u.minus(v)).scale('-1/4'));
 prove('circle:localization_rate',hg.times(r).minus(r.times(hg)),h.times(p.times(rp).plus(rp.times(p))).scale(N.IOTA.neg()));
}
{
 const{s,w:{m}}=register('packet_budget_polynomials',{tokens:['m'],rules:[]}),one=s.one();
 prove('budget:potential_majorant_numerator',m.plus(one).times(m.scale(2).plus(one)).scale(4).minus(m.times(m.scale(4).plus(one)).scale(2)),m.scale(10).plus(one.scale(4)));
 prove('budget:kinetic_majorant_numerator',m.scale(4).minus(one).minus(m.scale(3)),m.minus(one));
}
const chosen=proofs.find(x=>x.certificate.steps.length),bad=JSON.parse(JSON.stringify(chosen.certificate));bad.steps[0].coefficient=['991','0'];
no('integrity:corrupted_word_certificate',(()=>{try{systems.get(chosen.certificate.presentation).replay(bad);return false;}catch{return true;}})());
const hash=x=>crypto.createHash('sha256').update(x).digest('hex');
const payload={protocol:'ZCR1_FIXED_REFERENCE_RECOVERY_V1',status:'PASS_NATIVE_CLOSURE_PACKET_MOMENTS_AND_RECOVERY_BUDGET',scope:'Written finite-time flat-metric recovery for the full variable-profile AGD operator; native algebra and exact rational controls. Not a full classical-gravity or measurement-outcome theorem.',source_pins_sha256:hash(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),presentations,proofs,checks,negative_controls,results,summary:{exact_checks:checks.length,negative_controls:negative_controls.length,native_presentations:presentations.length,symbolic_certificates:proofs.length,serialized_replays:proofs.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',sha=hash(bytes),file=path.join(__dirname,'CERTIFICATE.json');
if(process.argv.includes('--write')){fs.writeFileSync(file,bytes);fs.writeFileSync(path.join(__dirname,'EXPECTED.sha256'),sha+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(file,'utf8'),bytes);assert.equal(fs.readFileSync(path.join(__dirname,'EXPECTED.sha256'),'utf8'),sha+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:sha,state_error_budget:'279/200000',observable_error_budget:'279/100000'}));
