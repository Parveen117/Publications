'use strict';
// GCP-1: original native rational-cut arithmetic and serialized word rewriting.
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const f=N.F.of,cut=N.Cut.of,I=N.identity(2),Z=N.matrix([[1,0],[0,-1]]),X=N.matrix([[0,1],[1,0]]),R=N.matrix([[0,-1],[1,0]]),Y=N.scale(R,N.IOTA),I4=N.identity(4);
const mm=(...a)=>a.reduce(N.mul),tr=a=>a.reduce((s,r,i)=>s.add(r[i]),N.ZERO),real=z=>{assert.ok(z.turn.zero());return z.rad;};
const sc=N.scale,eq=N.equal,zero=N.isZero,cm=N.commutator,kron=N.kron;
const P0=sc(N.add(I,Z),'1/2'),P1=N.sub(I,P0),plus=sc(N.add(I,X),'1/2');
const checks=[],controls=[],proofs=[],presentations=[],systems=new Map(),results={};
const ck=(name,v)=>{assert.ok(v,name);checks.push({name,passed:true});};
const no=(name,v)=>{assert.ok(v,name);controls.push({name,rejected:true});};
const conj=(u,r)=>mm(u,r,N.dagger(u));
function partial(r,keep){return N.matrix(Array.from({length:2},(_,i)=>Array.from({length:2},(_,j)=>[0,1].reduce((s,k)=>s.add(keep==='matter'?r[2*i+k][2*j+k]:r[2*k+i][2*k+j]),N.ZERO))));}
const purity=r=>real(tr(mm(r,r))),expect=(r,a)=>real(tr(mm(r,a))),diag=(a,b)=>N.matrix([[a,0],[0,b]]);
ck('native:cut_square',eq(mm(Z,Z),I));ck('native:turn_square',eq(mm(R,R),sc(I,-1)));ck('native:anticommutation',zero(N.add(mm(Z,R),mm(R,Z))));

// Exact rational-function differentiation, using only native F arithmetic.
const trim=p=>{p=p.map(f);while(p.length>1&&p.at(-1).zero())p.pop();return p;};
const padd=(p,q)=>trim(Array.from({length:Math.max(p.length,q.length)},(_,i)=>(p[i]||f(0)).add(q[i]||0)));
const pscale=(p,a)=>trim(p.map(v=>f(v).mul(a)));
function pmul(p,q){let out=Array(p.length+q.length-1).fill(f(0));p.forEach((a,i)=>q.forEach((b,j)=>out[i+j]=out[i+j].add(f(a).mul(b))));return trim(out);}
const pdiff=p=>p.length===1?[f(0)]:trim(p.slice(1).map((a,i)=>a.mul(i+1)));
const peval=(p,x)=>p.reduceRight((a,c)=>a.mul(x).add(c),f(0));
const rf=(n,d=[1])=>({n:trim(n),d:trim(d)});
const rdiff=a=>rf(padd(pmul(pdiff(a.n),a.d),pscale(pmul(a.n,pdiff(a.d)),-1)),pmul(a.d,a.d));
const req=(a,b)=>padd(pmul(a.n,b.d),pscale(pmul(b.n,a.d),-1)).every(z=>z.zero());
const reval=(a,x)=>peval(a.n,x).div(peval(a.d,x));
const aa=f('3/5'),v=rf([1,0,f(1).add(aa)],[1,0,1]),vp=rdiff(v),vpp=rdiff(vp);
ck('profile:first_derivative_identity',req(vp,rf([0,aa.mul(2)],[1,0,2,0,1])));
ck('profile:second_derivative_identity',req(vpp,rf([aa.mul(2),0,aa.mul(-6)],[1,0,3,0,3,0,1])));
const profile=[];
for(const xs of ['-3','-1','0','1/2','1','2','5']){const x=f(xs),speed=reval(v,x),first=reval(vp,x),second=reval(vpp,x),curvature=second.mul(-2).div(speed);
 ck('profile:positive_bounded_fixture:'+xs,f(1).le(speed)&&speed.le('8/5'));
 const torsionCoefficient=first.neg().add(first);ck('profile:torsion_free_fixture:'+xs,torsionCoefficient.zero());
 profile.push({x,speed,first_derivative:first,second_derivative:second,scalar_curvature:curvature});
}
ck('profile:curvature_at_origin',profile.find(p=>p.x.eq(0)).scalar_curvature.eq('-12/5'));
ck('profile:curvature_at_one',profile.find(p=>p.x.eq(1)).scalar_curvature.eq('6/13'));
no('profile:same_value_and_slope_imply_same_curvature',reval(v,0).eq(1)&&reval(vp,0).eq(0)&&!reval(vpp,0).zero());
no('transport:omit_half_derivative',!reval(vp,1).zero());
results.curved_profile={a:aa,formula:'v1(x)=1+a*x^2/(1+x^2)',points:profile,uniform_bound:'1 <= v_j(x) <= 8/5',curvature_convention:'R=-2 v_double_prime/v'};
const gamma0=R,gamma1=mm(Z,R);
ck('coframe:clifford_time',eq(mm(gamma0,gamma0),sc(I,-1)));
ck('coframe:clifford_space',eq(mm(gamma1,gamma1),I));
ck('coframe:clifford_cross',zero(N.add(mm(gamma0,gamma1),mm(gamma1,gamma0))));
for(const p of profile){
 const multiplier=sc(R,p.speed),spin=sc(Z,p.first_derivative.div(-2));
 ck('coframe:time_coefficient:'+p.x,eq(mm(multiplier,sc(R,f(1).div(p.speed))),sc(I,-1)));
 ck('coframe:spatial_coefficient:'+p.x,eq(mm(multiplier,gamma1),sc(Z,p.speed)));
 ck('coframe:connection_match:'+p.x,eq(mm(multiplier,sc(R,f(1).div(p.speed)),spin),sc(Z,p.first_derivative.div(2))));
 ck('coframe:curvature_match:'+p.x,eq(sc(Z,p.second_derivative.div(2)),sc(Z,p.speed.mul(p.scalar_curvature).div(-4))));
}
ck('coframe:nonzero_native_curvature_at_origin',eq(sc(Z,reval(vpp,0).div(2)),sc(Z,'3/5')));

// Matrix-level principal leakage, measured on the retained geometry cut.
const V=diag(1,2),A=kron(Z,V),P=kron(I,plus),Q=N.sub(I4,P),vbar=f('3/2'),variance=f('1/4');
ck('compression:retained_cut',eq(mm(P,P),P));
ck('compression:exact_principal_variance',eq(N.sub(mm(P,A,A,P),mm(P,A,P,A,P)),sc(P,variance)));
ck('compression:lost_channel_square',eq(mm(P,A,Q,A,P),sc(P,variance)));
no('compression:mean_only_closure',!zero(mm(Q,A,P)));
ck('compression:scalar_lapse_zero_leakage',zero(mm(Q,kron(Z,sc(I,'3/2')),P)));
ck('compression:curved_point_variance',f('13/10').sub(1).pow(2).div(4).eq('9/400'));
for(const [ts,xs]of [['1','1'],['2','1'],['3/2','1'],['0','0'],['7/3','2/5']]){
 const t=f(ts),x=f(xs),symbol=N.sub(sc(I4,t),sc(A,x)),det=symbol.reduce((acc,row,i)=>acc.mul(row[i]),N.ONE);
 const fact=t.pow(2).sub(x.pow(2)).mul(t.pow(2).sub(x.pow(2).mul(4)));
 ck('symbol:two_branch_determinant:'+ts+':'+xs,det.eq(fact));
}
no('symbol:mean_cone_matches_branches',!f(1).sub(vbar.pow(2)).zero());
results.compression={mean_speed:vbar,second_moment:f('5/2'),lost_principal_coefficient:variance,curved_point_lost_coefficient:'9/400'};

// Exact finite-mode feedback and energy exchange.
const omega=f('2/3'),delta=f('-1/2'),G=k=>N.add(sc(kron(Z,V),N.IOTA.mul(k)),sc(kron(I,R),omega));
const initial=kron(P0,P0),observable=kron(I,P1),response=[];
for(const k of [0,1,2]){
 const generator=G(k);ck('mode:anti_dagger:'+k,eq(N.dagger(generator),sc(generator,-1)));
 let derivative=observable,values=[];
 for(let n=1;n<=4;n++){derivative=cm(derivative,generator);values.push(expect(initial,derivative));}
 const h2=omega.pow(2).add(f(k).pow(2).mul(delta.pow(2)));
 ck('feedback:first_population_derivative:'+k,values[0].eq(0));
 ck('feedback:second_population_derivative:'+k,values[1].eq(omega.pow(2).mul(2)));
 ck('feedback:third_population_derivative:'+k,values[2].eq(0));
 ck('feedback:fourth_population_derivative:'+k,values[3].eq(omega.pow(2).mul(h2).mul(-8)));
 response.push({k,derivatives:values,h_squared:h2});
}
no('feedback:second_derivative_is_complete',response[0].derivatives[1].eq(response[1].derivatives[1])&&!response[0].derivatives[3].eq(response[1].derivatives[3]));
const Hkin=sc(kron(Z,V),-1),Hgeo=sc(kron(I,R),N.IOTA.mul(omega)),total=N.add(Hkin,Hgeo),exchange=cm(Hkin,G(1));
ck('exchange:kinetic_identity',eq(exchange,sc(kron(Z,X),f(2).mul(omega).mul(delta))));
ck('exchange:joint_conservation',zero(cm(total,G(1))));
ck('exchange:equal_opposite',eq(cm(Hgeo,G(1)),sc(exchange,-1)));
const exchangeInput=kron(P0,plus);ck('exchange:kinetic_reading',expect(exchangeInput,exchange).eq('-2/3'));ck('exchange:geometry_reading',expect(exchangeInput,cm(Hgeo,G(1))).eq('2/3'));
results.feedback={population_derivatives:response,kinetic_exchange:'-2/3',geometry_exchange:'2/3'};

// The exact exponential at h*s=acos(3/5), after a removable local direction phase.
const c=f('3/5'),b=f('16/25'),d=f('12/25');
const U=N.add(sc(I4,c),N.add(sc(kron(I,R),b),sc(kron(Z,Z),N.IOTA.mul(d).neg())));
ck('demo:unitarity',eq(mm(N.dagger(U),U),I4));
const input=kron(plus,P0),out=conj(U,input),rm=partial(out,'matter'),rg=partial(out,'geometry'),entropy=f(1).sub(purity(rg));
ck('demo:input_pure',eq(mm(input,input),input));ck('demo:output_pure',eq(mm(out,out),out));
ck('demo:trace_one',tr(out).eq(1));
ck('demo:geometry_population',expect(rg,P1).eq('256/625'));
ck('demo:exact_entanglement',entropy.eq('73728/390625'));
ck('demo:matching_reduced_purities',purity(rm).eq(purity(rg)));
const amp=N.matrix([[cut([c,d.neg()])],[b],[cut([c,d])],[b]]);
ck('demo:amplitude_formula',eq(sc(mm(amp,N.dagger(amp)),'1/2'),out));
const mixed=conj(U,kron(sc(I,'1/2'),P0));
const gp=N.matrix([[cut([c,d.neg()])],[b]]),gm=N.matrix([[cut([c,d])],[b]]);
const separable=sc(N.add(kron(P0,mm(gp,N.dagger(gp))),kron(P1,mm(gm,N.dagger(gm)))),'1/2');
ck('control:mixed_input_separable_decomposition',eq(mixed,separable));
ck('control:same_geometry_marginal',eq(partial(mixed,'geometry'),rg));
no('control:geometry_entropy_alone_certifies_entanglement',eq(partial(mixed,'geometry'),rg)&&!eq(mixed,out));
const sameMetric=kron(I,N.add(sc(I,'3/5'),sc(R,'4/5'))),sameOut=conj(sameMetric,input);
ck('control:scalar_metric_factorizes',purity(partial(sameOut,'geometry')).eq(1));
const noMix=N.add(sc(I4,'3/5'),sc(kron(Z,Z),N.IOTA.mul('4/5'))),noMixOut=conj(noMix,input);
ck('control:no_mixing_eigengeometry_product',purity(partial(noMixOut,'geometry')).eq(1));
results.demo={unitary:U,initial_state:input,final_state:out,matter_reduction:rm,geometry_reduction:rg,geometry_population:'256/625',linear_entropy:entropy,mixed_control:mixed,scope:'Pure direction-geometry preparation at one Fourier mode. This is an exact model output, not an experiment.'};

// Homogeneous exact Schur memory control on the same finite carrier.
const T=kron(I,N.add(sc(I,'3/5'),sc(R,'4/5'))),Grot=mm(N.dagger(T),G(1),T),z=f(3);
const keep=[0,2],lost=[1,3],block=(a,rows,cols)=>N.matrix(rows.map(i=>cols.map(j=>a[i][j])));
const aaBlock=block(Grot,keep,keep),bbBlock=block(Grot,keep,lost),ccBlock=block(Grot,lost,keep),ddBlock=block(Grot,lost,lost);
const ret=mm(bbBlock,N.inverse(N.sub(sc(I,z),ddBlock)),ccBlock),eff=N.inverse(N.sub(N.sub(sc(I,z),aaBlock),ret));
ck('memory:schur_return_identity',eq(block(N.inverse(N.sub(sc(I4,z),Grot)),keep,keep),eff));
no('memory:discard_return_is_exact',!eq(eff,N.inverse(N.sub(sc(I,z),aaBlock))));

// Exact frame-jet check. T_x=T Jx and T_s=T Js enforce the unitary jet relations.
const T0=kron(N.add(sc(I,'3/5'),sc(R,'4/5')),N.add(sc(I,'5/13'),sc(R,'12/13'))),Jx=kron(R,Z),Js=kron(Z,R),Tx=mm(T0,Jx),Ts=mm(T0,Js);
const Aprime=kron(Z,diag(0,'3/10')),K=N.add(sc(kron(R,V),'2/5'),sc(kron(I,R),omega)),Tdag=N.dagger(T0),Txd=N.dagger(Tx);
const Anew=mm(Tdag,A,T0),AnewPrime=N.add(mm(Txd,A,T0),N.add(mm(Tdag,Aprime,T0),mm(Tdag,A,Tx)));
const Knew=N.sub(N.add(mm(Tdag,K,T0),sc(N.sub(mm(Tdag,A,Tx),mm(Txd,A,T0)),'1/2')),mm(Tdag,Ts));
const directLower=N.sub(N.add(mm(Tdag,A,Tx),mm(Tdag,N.add(sc(Aprime,'1/2'),K),T0)),mm(Tdag,Ts));
ck('frame:complete_lower_term',eq(directLower,N.add(sc(AnewPrime,'1/2'),Knew)));
ck('frame:anti_dagger_connection',eq(N.dagger(Knew),sc(Knew,-1)));
ck('frame:hermitian_principal',eq(N.dagger(Anew),Anew));
no('frame:drop_connection_terms',!eq(Knew,mm(Tdag,K,T0)));

// Native word proofs: each exact identity has a serialized replayable witness.
function system(name,tokens,rules=[],commute=false){if(commute)for(let j=1;j<tokens.length;j++)for(let i=0;i<j;i++)rules.push({id:'central_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'Declared central coefficients'});const s=new M.Presentation({tokens,rules}),audit=s.audit();ck(name+':confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name,spec:s.spec(),audit});systems.set(s.hash(),s);return{s,w:Object.fromEntries(tokens.map(t=>[t,s.word([t])]))};}
function prove(name,lhs,rhs){const s=lhs.system,r=s.reduce(lhs.minus(rhs),{witness:true});ck(name,r.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(r.certificate));ck(name+':replay',s.replay(certificate));proofs.push({name,lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate});}
{
 const{s,w:{p,v,w}}=system('geometry_variance',['p','v','w'],[],true),one=s.one(),mean=p.times(v).plus(one.minus(p).times(w)),second=p.times(v).times(v).plus(one.minus(p).times(w).times(w));
 prove('variance:two_branch_identity',second.minus(mean.times(mean)),p.times(one.minus(p)).times(v.minus(w)).times(v.minus(w)));
 const v2=v.times(v),w2=w.times(w);
 prove('cone:repeated_root_obstruction',v2.plus(w2).times(v2.plus(w2)).minus(v2.times(w2).scale(4)),v2.minus(w2).times(v2.minus(w2)));
}
{
 const{s,w:{p,a}}=system('retained_principal',['p','a'],[{id:'projector',lhs:['p','p'],rhs:[[['p'],1]],source:'Retained cut'}]);
 prove('memory:principal_lost_channel',p.times(a).times(a).times(p).minus(p.times(a).times(p).times(a).times(p)),p.times(a).times(s.one().minus(p)).times(a).times(p));
}
{
 const{s,w:{r,z,b,d}}=system('traceless_mode',['r','z','b','d'],[
  {id:'r_square',lhs:['r','r'],rhs:[[[],-1]],source:'Native turn'},
  {id:'z_square',lhs:['z','z'],rhs:[[[],1]],source:'Native cut'},
  {id:'anti',lhs:['z','r'],rhs:[[['r','z'],-1]],source:'Cut/turn anticommutation'},
  ...['b','d'].flatMap(t=>['r','z'].map(u=>({id:t+'_'+u,lhs:[t,u],rhs:[[[u,t],1]],source:'Central coefficient'}))),
  {id:'db',lhs:['d','b'],rhs:[[['b','d'],1]],source:'Central coefficients'}
 ]);
 const h=r.times(b).plus(z.times(d).scale(N.IOTA));
 prove('mode:traceless_generator_square',h.times(h),b.times(b).plus(d.times(d)).scale(-1));
}
{
 const{s,w:{c,b,d}}=system('pure_joint_entropy',['c','b','d'],[{id:'normalization',lhs:['d','d'],rhs:[[[],1],[['c','c'],-1],[['b','b'],-1]],source:'Normalized two branch vectors'}],true);
 const realOverlap=c.times(c).minus(d.times(d)).plus(b.times(b)),imagOverlap=c.times(d).scale(2);
 prove('entropy:overlap_identity',s.one().minus(realOverlap.times(realOverlap)).minus(imagOverlap.times(imagOverlap)),b.times(b).times(d).times(d).scale(4));
 prove('entropy:linear_entropy',s.one().minus(realOverlap.times(realOverlap)).minus(imagOverlap.times(imagOverlap)).scale('1/2'),b.times(b).times(d).times(d).scale(2));
}
const chosen=proofs.find(p=>p.certificate.steps.length),bad=JSON.parse(JSON.stringify(chosen.certificate));bad.steps[0].coefficient=['997','0'];
no('integrity:corrupted_word_proof',(()=>{try{systems.get(chosen.certificate.presentation).replay(bad);return false;}catch{return true;}})());
const hash=x=>crypto.createHash('sha256').update(x).digest('hex');
const payload={protocol:'GCP1_GEOMETRY_CONTROLLED_PROPAGATION_V1',status:'PASS_NATIVE_ALGEBRA_AND_EXACT_QUANTUM_DEMO',scope:'Written analytic existence/locality and branch metric proofs; native serialized polynomial identities and finite exact fixtures. Selected operator-valued lapse and geometry mixing. Not a physical quantum-gravity completion.',source_pins_sha256:hash(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),presentations,proofs,checks,negative_controls:controls,results,summary:{exact_checks:checks.length,negative_controls:controls.length,native_presentations:presentations.length,symbolic_certificates:proofs.length,serialized_replays:proofs.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',sha=hash(bytes),file=path.join(__dirname,'CERTIFICATE.json');
if(process.argv.includes('--write')){fs.writeFileSync(file,bytes);fs.writeFileSync(path.join(__dirname,'EXPECTED.sha256'),sha+'  CERTIFICATE.json\n');}
else{assert.equal(fs.readFileSync(file,'utf8'),bytes);assert.equal(fs.readFileSync(path.join(__dirname,'EXPECTED.sha256'),'utf8'),sha+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:sha,linear_entropy:entropy,geometry_population:'256/625'}));
