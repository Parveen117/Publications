'use strict';
// Native arithmetic only. No alternative scalar, polynomial or matrix engine.
const fs=require('node:fs'), path=require('node:path'), crypto=require('node:crypto');
const assert=require('node:assert/strict'), root=__dirname, args=process.argv.slice(2);
const arg=k=>{const i=args.indexOf(k);assert.ok(i>=0&&args[i+1],k);return args[i+1];};
const nativeRoot=path.resolve(arg('--native-root')), upstream=JSON.parse(fs.readFileSync(arg('--upstream-report'),'utf8'));
const write=args.includes('--write');assert.ok(write!==args.includes('--check'),'Choose --write or --check');
const hash=x=>crypto.createHash('sha256').update(x).digest('hex');
const pins=JSON.parse(fs.readFileSync(path.join(root,'SOURCE_PINS.json'),'utf8'));
for(const [file,expected] of Object.entries(pins.native_engine.files))assert.equal(hash(fs.readFileSync(path.join(nativeRoot,file))),expected,file);
const N=require(path.join(nativeRoot,'native_operator.cjs')), M=require(path.join(nativeRoot,'paninian_operator.cjs'));
const {F,Cut,IOTA}=N, f=x=>F.of(x), I=N.identity(2), O=N.zeros(2);
const R=N.matrix([[0,-1],[1,0]]), J=N.scale(N.matrix([[0,1],[1,0]]),IOTA), L=N.mul(R,J);
const basis=[R,J,L], mm=(...xs)=>xs.reduce(N.mul), trace=A=>A[0][0].add(A[1][1]);
const real=x=>{assert.ok(x.turn.zero(),'Real native scalar required');return x.rad;};
const sc=A=>real(trace(A)).div(2), pair=(A,B)=>sc(N.mul(A,B)).neg();
const combo=coeffs=>coeffs.reduce((a,c,i)=>N.add(a,N.scale(basis[i],c)),O);
const checks=[],controls=[],certificates=[],presentations=[],systems=new Map();
const check=(name,ok,detail={})=>{assert.ok(ok,name);checks.push({name,passed:true,...detail});};
const control=(name,ok,detail={})=>{assert.ok(ok,name);controls.push({name,false_inference_rejected:true,...detail});};
const rejects=fn=>{try{fn();return false;}catch{return true;}};
function make(name,tokens,{commutative=false,rules=[],star=null}={}){
  const commute=[];if(commutative)for(let j=1;j<tokens.length;j++)for(let i=0;i<j;i++)commute.push({id:tokens[j]+'_past_'+tokens[i],lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'Real central native coefficients'});
  const S=new M.Presentation({tokens,rules:commute.concat(rules)}),w=Object.fromEntries(tokens.map(t=>[t,S.word([t])])),audit=S.audit();
  check(name+':confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');
  const images=star?star(S,w):commutative?w:null,dagger=images?M.daggerGate(S,images):null;
  if(dagger)check(name+':dagger',dagger.status==='DAGGER_DESCENDS');
  presentations.push({name,spec:S.spec(),audit,dagger,dagger_images:images?Object.fromEntries(Object.entries(images).map(([k,e])=>[k,e.toJSON()])):null});
  systems.set(S.hash(),S);return {S,w};
}
function prove(name,lhs,rhs){
  const S=lhs.system,r=M.proveEquality(lhs,rhs);check(name,r.status==='EQUAL_IN_DECLARED_QUOTIENT');
  const c=JSON.parse(JSON.stringify(r.certificate));assert.ok(S.fromJSON(c.input).equals(lhs.minus(rhs)));assert.ok(S.fromJSON(c.output).equals(S.zero(lhs.from,lhs.to)));
  check(name+':serialized_replay',S.replay(c));certificates.push({name,lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate:c});
}
const times=(...xs)=>xs.reduce((a,b)=>a.times(b)),sq=x=>x.times(x);
function sym(S){
  const one=S.one(),zero=S.zero(),constant=A=>A.map(row=>row.map(c=>one.scale(c)));
  const add=(A,B)=>A.map((row,i)=>row.map((x,j)=>x.plus(B[i][j]))),sub=(A,B)=>A.map((row,i)=>row.map((x,j)=>x.minus(B[i][j])));
  const scale=(A,x)=>A.map(row=>row.map(a=>a.times(x))),mul=(A,B)=>A.map(row=>B[0].map((_,j)=>row.reduce((z,x,k)=>z.plus(x.times(B[k][j])),zero)));
  const tr=A=>A[0][0].plus(A[1][1]),eq=(name,A,B)=>{for(let i=0;i<2;i++)for(let j=0;j<2;j++)prove(name+':'+i+j,A[i][j],B[i][j]);};
  const sb=basis.map(constant),vector=c=>c.reduce((A,x,i)=>add(A,scale(sb[i],x)),constant(O));
  return {constant,add,sub,scale,mul,tr,eq,vector,id:constant(I),zero:constant(O)};
}

// Inherited cut identities expose the actual force in the second cut jet.
{
  const tokens=['p','v','h'],rules=[
    {id:'cut',lhs:['p','p'],rhs:[[['p'],1]],source:'P^2=P'},
    {id:'horizontal',lhs:['v','p'],rhs:[[['v'],1],[['p','v'],-1]],source:'PV+VP=V'}];
  const {S,w:{p,v,h}}=make('R1_native_force_extraction',tokens,{rules,star:(_,w)=>({p:w.p,v:w.v.scale(-1),h:w.h})});
  const br=M.bracket,Fc=br(p,h),X=br(v,p),Z=br(Fc,p).plus(br(v,X));
  prove('R1:force_from_second_cut_derivative',br(Z,p),Fc);
  prove('R1:rest_cut_acceleration_is_tangent',p.times(br(Fc,p)).plus(br(Fc,p).times(p)),br(Fc,p));
  prove('R1:parameter_acceleration_recovers_velocity',br(X,p),v);
  prove('R1:cut_diagonal_gradient_is_unobservable',br(p,times(p,h,p).plus(times(S.one().minus(p),h,S.one().minus(p)))),S.zero());
}
{
  const {S,w:{a,b,c,u,v,k,r}}=make('R2_native_Gram_reconstruction',['a','b','c','u','v','k','r'],{commutative:true});
  const delta=a.times(c).minus(sq(b)),read1=k.times(a.times(u).plus(b.times(v))),read2=k.times(b.times(u).plus(c.times(v)));
  const nA=c.times(read1).minus(b.times(read2)),nB=a.times(read2).minus(b.times(read1));
  prove('R2:information_coefficient_numerator',nA,k.times(delta).times(u));
  prove('R2:curvature_coefficient_numerator',nB,k.times(delta).times(v));
  prove('R2:unknown_parameter_speed_cancels',nA.times(v).minus(nB.times(u)),S.zero());
  const aa=sq(k).times(a.times(sq(u)).plus(b.times(u).times(v).scale(2)).plus(c.times(sq(v))));
  const ff=a.plus(b.times(r).scale(2)).plus(c.times(sq(r)));
  const af=k.times(a.times(u).plus(b.times(u.times(r).plus(v))).plus(c.times(v).times(r)));
  prove('R2:independent_prediction_residual',aa.times(ff).minus(sq(af)),sq(k).times(sq(u.times(r).minus(v))).times(delta));
}

// Native matrix calculus for the control family (s R,J,L), with a sphere cut.
{
  const {S,w:{s,x,y,z}}=make('R3_deformed_control_cut',['s','x','y','z'],{commutative:true,rules:[
    {id:'unit_sphere',lhs:['z','z'],rhs:[[[],1],[['x','x'],-1],[['y','y'],-1]],source:'Rank-one cut constraint'}]});
  const raw=sym(S),norm=A=>A.map(row=>row.map(x=>S.reduce(x).normal));
  const m={...raw,add:(A,B)=>norm(raw.add(A,B)),sub:(A,B)=>norm(raw.sub(A,B)),mul:(A,B)=>norm(raw.mul(A,B)),tr:A=>S.reduce(raw.tr(A)).normal};
  const sm=(A,c)=>A.map(row=>row.map(x=>x.scale(c))),br=(A,B)=>m.sub(m.mul(A,B),m.mul(B,A));
  const sigma=[m.constant(N.matrix([[0,1],[1,0]])),m.constant(N.scale(R,IOTA)),m.constant(N.matrix([[1,0],[0,-1]]))];
  const vec=cs=>cs.reduce((A,c,i)=>m.add(A,m.scale(sigma[i],c)),m.zero),P=sm(m.add(m.id,vec([x,y,z])),'1/2');
  const source=B=>{
    const X=B.map(b=>br(b,P)),w=[1,1,2],pairs=[[0,1,1],[0,2,2],[1,2,2]];
    let HI=m.zero,HC=m.zero,I0=S.zero(),C0=S.zero();
    X.forEach((q,a)=>{HI=m.sub(HI,sm(br(B[a],q),w[a]));I0=I0.plus(m.tr(m.mul(q,q)).scale(f(w[a]).div(2)));});
    pairs.forEach(([a,b,h])=>{const W=br(X[a],X[b]),C=m.mul(m.mul(P,W),P);
      HC=m.add(HC,sm(m.add(sm(m.add(m.mul(W,C),m.mul(C,W)),-1),m.add(br(B[a],br(X[b],C)),br(B[b],br(C,X[a])))),h));
      C0=C0.minus(m.tr(m.mul(C,C)).scale(f(h).div(2)));
    });return {I0,C0,FI:br(P,HI),FC:br(P,HC)};
  };
  const B=basis.map(m.constant),plain=source(B),deformed=source([m.scale(B[0],s),B[1],B[2]]),u=sq(s);
  prove('R3:symmetric_action_blindness',plain.C0.minus(plain.I0.scale(2)),S.one().scale(-2));
  m.eq('R3:symmetric_force_blindness',plain.FC,sm(plain.FI,2));
  prove('R3:deformed_information_functional',deformed.I0,u.plus(S.one().scale(3)).minus(sq(x)).minus(u.times(sq(y))).minus(sq(z).scale(2)));
  prove('R3:deformed_curvature_functional',deformed.C0,u.times(sq(x)).scale(4).plus(sq(y).scale(4)).plus(u.times(sq(z)).scale(2)));
  const HI=sm(vec([x,u.times(y),z.scale(2)]),-2),HC=sm(vec([u.times(x).scale(4),y.scale(4),u.times(z).scale(2)]),2);
  m.eq('R3:information_force_form',deformed.FI,br(P,HI));
  m.eq('R3:curvature_force_form',deformed.FC,br(P,HC));
  const pairing=(A,B)=>m.tr(m.mul(A,B)).scale(-1),g11=pairing(deformed.FI,deformed.FI),g12=pairing(deformed.FI,deformed.FC),g22=pairing(deformed.FC,deformed.FC);
  const determinant=g11.times(g22).minus(sq(g12)),expected=sq(u.minus(S.one())).times(sq(u.minus(S.one().scale(2)))).times(sq(x)).times(sq(y)).times(sq(z)).scale(256);
  prove('R3:factorized_force_Gram_determinant',determinant,expected);
  prove('R3:second_blind_control_setting',deformed.C0.minus(deformed.I0.scale(4)).plus(S.one().scale(8)),u.minus(S.one().scale(2)).times(sq(z)).scale(-2));
}

const br=N.commutator,fullPair=(A,B)=>real(trace(mm(A,B))).neg(),sum=xs=>xs.reduce((s,x)=>s.add(x),f(0));
const sx=N.matrix([[0,1],[1,0]]),sy=N.scale(R,IOTA),sz=N.matrix([[1,0],[0,-1]]);
const cut=n=>{n=n.map(f);assert.ok(sum(n.map(x=>x.pow(2))).eq(1));return N.scale(n.reduce((A,x,i)=>N.add(A,N.scale([sx,sy,sz][i],x)),I),'1/2');};
const trans=(G,A)=>mm(N.dagger(G),A,G),anti=A=>N.equal(N.dagger(A),N.scale(A,-1)),herm=A=>N.equal(N.dagger(A),A);
function forces(P,B){
  assert.ok(herm(P)&&N.equal(mm(P,P),P));assert.ok(B.every(anti));
  const X=B.map(b=>br(b,P)),w=[1,1,2],pairs=[[0,1,1],[0,2,2],[1,2,2]];
  let HI=O,HC=O,I0=f(0),C0=f(0);
  X.forEach((x,a)=>{HI=N.sub(HI,N.scale(br(B[a],x),w[a]));I0=I0.add(N.energy(x).mul(w[a]).div(2));});
  pairs.forEach(([a,b,h])=>{const W=br(X[a],X[b]),C=mm(P,W,P),H=N.add(N.scale(N.add(mm(W,C),mm(C,W)),-1),N.add(br(B[a],br(X[b],C)),br(B[b],br(C,X[a]))));HC=N.add(HC,N.scale(H,h));C0=C0.add(N.energy(C).mul(h).div(2));});
  const FI=br(P,HI),FC=br(P,HC),g11=fullPair(FI,FI),g12=fullPair(FI,FC),g22=fullPair(FC,FC),delta=g11.mul(g22).sub(g12.pow(2));
  return {P,FI,FC,I0,C0,G:[[g11,g12],[g12,g22]],delta};
}
function recover(P,B,first,second){
  assert.ok(N.isZero(first),'Initial covariant cut derivative must vanish');
  assert.ok(herm(second)&&N.equal(N.add(mm(P,second),mm(second,P)),second),'Initial-rest second derivative must be cut tangent');
  const r=forces(P,B);assert.ok(!r.delta.le(0),'Independent force directions required');
  const A=br(second,P),b=[fullPair(r.FI,A),fullPair(r.FC,A)],[[a,d],[,c]]=r.G;
  const alpha=c.mul(b[0]).sub(d.mul(b[1])).div(r.delta),beta=a.mul(b[1]).sub(d.mul(b[0])).div(r.delta);
  assert.ok(f(0).le(alpha)&&f(0).le(beta),'Nonnegative action coefficients required');
  assert.ok(!alpha.add(beta).zero(),'A zero response identifies no nonzero ratio');
  const residual=N.sub(A,N.add(N.scale(r.FI,alpha),N.scale(r.FC,beta)));assert.ok(N.isZero(residual),'Observed force is outside the declared span');
  return {...r,A,readings:b,alpha,beta,branch:alpha.zero()?'pure_curvature':'finite_ratio',ratio:alpha.zero()?null:beta.div(alpha),residual_norm_squared:N.energy(residual)};
}
function response(P,B,alpha,beta,q='1',secondParameter='0'){
  const r=forces(P,B),Fsource=N.add(N.scale(r.FI,alpha),N.scale(r.FC,beta)),first=O;
  const second=N.add(N.scale(br(Fsource,P),f(q).pow(2)),N.scale(first,secondParameter));
  return {source:r,Fsource,first,second,parameter_speed:f(q),parameter_acceleration:f(secondParameter),reading:recover(P,B,first,second)};
}
function diagnostic(A,Fp){
  const aa=N.energy(A),ff=N.energy(Fp);assert.ok(!aa.zero()&&!ff.zero(),'Nonzero observation and prediction required');
  const af=fullPair(A,Fp),gap=aa.mul(ff).sub(af.pow(2));assert.ok(f(0).le(gap));
  return {gap,normalized:gap.div(aa.mul(ff)),orientation:af,accepted:gap.zero()&&!af.le(0)};
}
const nstar=['2/3','2/3','1/3'],Pstar=cut(nstar),Bstar=[N.scale(R,2),J,L],fixture=response(Pstar,Bstar,2,1);
check('R4:exact_G_and_determinant',fixture.reading.G[0][0].eq('1312/81')&&fixture.reading.G[0][1].eq('5120/81')&&fixture.reading.G[1][1].eq('20992/81')&&fixture.reading.delta.eq('16384/81'));
check('R4:exact_two_reading_fixture',fixture.reading.readings[0].eq('7744/81')&&fixture.reading.readings[1].eq('31232/81')&&fixture.reading.ratio.eq('1/2'));
const observed=[],frames=[],blindness=[],predictions=[];
for(const [alpha,beta] of [['2','1'],['3','1/2'],['4','0'],['0','2']]){
  const old=forces(Pstar,basis),Fsource=N.add(N.scale(old.FI,alpha),N.scale(old.FC,beta));
  check('R4:old_protocol_same_force:'+alpha+','+beta,N.equal(Fsource,N.scale(old.FI,4))&&old.delta.zero());
  blindness.push({alpha:f(alpha),beta:f(beta),Fsource,potential:old.I0.mul(alpha).add(old.C0.mul(beta))});
  for(const [q,q2] of [['1/3','5'],['1','0'],['-2','-3'],['5/2','7/4']]){
    const r=response(Pstar,Bstar,alpha,beta,q,q2),actual=r.reading;
    check('R4:initial_rest_parameter_cancellation:'+alpha+','+beta+','+q,actual.alpha.eq(f(alpha).mul(f(q).pow(2)))&&actual.beta.eq(f(beta).mul(f(q).pow(2)))&&(f(alpha).zero()?actual.ratio===null:actual.ratio.eq(f(beta).div(alpha))));
    observed.push({alpha:f(alpha),beta:f(beta),...r});
  }
}
for(const t of ['1/7','1/3','-2/5']){
  const G=N.mul(N.cayley(R,t),N.cayley(J,'1/4')),P=trans(G,Pstar),B=Bstar.map(b=>trans(G,b)),Z=trans(G,fixture.second),r=recover(P,B,O,Z);
  check('R5:unit_frame_invariance:'+t,r.ratio.eq('1/2')&&r.delta.eq(fixture.reading.delta)&&r.readings.every((x,i)=>x.eq(fixture.reading.readings[i])));
  const connection=N.add(R,N.scale(J,'1/3')),connectionDerivative=L,ordinaryFirst=N.scale(br(connection,P),-1);
  const ordinarySecond=N.add(N.sub(Z,br(connectionDerivative,P)),br(connection,br(connection,P)));
  const covFirst=N.add(ordinaryFirst,br(connection,P));
  const covSecond=N.add(N.add(ordinarySecond,br(connectionDerivative,P)),N.add(N.scale(br(connection,ordinaryFirst),2),br(connection,br(connection,P))));
  const corrected=recover(P,B,covFirst,covSecond);
  check('R5:variable_frame_covariant_rest_reconstruction:'+t,corrected.ratio.eq(r.ratio)&&N.equal(covSecond,Z));
  control('R5:ordinary_derivative_does_not_certify_covariant_rest:'+t,rejects(()=>recover(P,B,ordinaryFirst,ordinarySecond)));
  frames.push({parameter:t,reading:r,corrected_ratio:corrected.ratio});
}
const targetRatio=fixture.reading.ratio;
for(const [n,s,q] of [[['1/3','2/3','2/3'],'2','3/2'],[['2/7','3/7','6/7'],'3','1/5'],[['-2/3','1/3','2/3'],'3/2','-3']]){
  const P=cut(n),B=[N.scale(R,s),J,L],r=forces(P,B),gain=f(q).pow(2),A=N.scale(N.add(N.scale(r.FI,2),r.FC),gain);
  const predicted=N.add(r.FI,N.scale(r.FC,targetRatio)),good=diagnostic(A,predicted),wrongRatio=f('1/6'),bad=diagnostic(A,N.add(r.FI,N.scale(r.FC,wrongRatio)));
  check('R6:independent_preparation_prediction:'+n.join(','),good.accepted&&good.normalized.zero());
  check('R6:exact_wrong_ratio_residual:'+n.join(','),bad.gap.eq(gain.pow(2).mul(f(2).mul(wrongRatio).sub(1).pow(2)).mul(r.delta))&&!bad.gap.zero());
  control('R6:wrong_ratio_rejected_without_refitting:'+n.join(','),!bad.accepted);
  predictions.push({coordinates:n,control_amplitude:f(s),parameter_speed:f(q),inferred_ratio_used:targetRatio,good,bad});
}
control('R7:blind_protocol_cannot_reconstruct_coefficients',rejects(()=>response(Pstar,basis,2,1)));
control('R7:zero_coordinate_is_blind_even_with_amplitude_change',rejects(()=>response(cut(['4/5','0','3/5']),Bstar,2,1)));
const rotated=N.cayley(J,'1/3'),blindRotated=forces(trans(rotated,Pstar),basis.map(b=>trans(rotated,b)));
control('R7:frame_change_does_not_repair_force_blindness',blindRotated.delta.zero()&&N.equal(blindRotated.FC,N.scale(blindRotated.FI,2)));
control('R7:negative_action_coefficient_rejected',rejects(()=>response(Pstar,Bstar,-1,2)));
control('R7:zero_response_has_no_ratio',rejects(()=>response(Pstar,Bstar,0,0)));
control('R7:singular_parameter_rejected',rejects(()=>response(Pstar,Bstar,2,1,'0','3')));
const Fp=N.add(fixture.reading.FI,N.scale(fixture.reading.FC,targetRatio)),reversed=diagnostic(N.scale(Fp,-1),Fp);
control('R7:squared_alignment_without_orientation_is_insufficient',reversed.normalized.zero()&&!reversed.accepted);
control('R7:zero_prediction_denominator_rejected',rejects(()=>diagnostic(fixture.reading.A,O)));
// Away from rest a nonlinear reparameterization adds f'' V to the force.
const V=fixture.reading.FC,X=br(V,Pstar),Fsource=fixture.Fsource,q=f(2),q2=f(1);
const originalSecond=N.add(br(Fsource,Pstar),br(V,br(V,Pstar)));
const changedSecond=N.add(N.scale(originalSecond,q.pow(2)),N.scale(X,q2)),changedForce=br(changedSecond,Pstar);
check('R7:nonrest_parameter_term',N.equal(changedForce,N.add(N.scale(Fsource,q.pow(2)),N.scale(V,q2))));
control('R7:initial_rest_cannot_be_dropped',rejects(()=>recover(Pstar,Bstar,N.scale(X,q),changedSecond))&&!N.equal(changedForce,N.scale(Fsource,q.pow(2))));
const wrongReadings=[fullPair(fixture.reading.FI,changedForce),fullPair(fixture.reading.FC,changedForce)],[[g11,g12],[,g22]]=fixture.reading.G;
const wrongAlpha=g22.mul(wrongReadings[0]).sub(g12.mul(wrongReadings[1])).div(fixture.reading.delta);
const wrongBeta=g11.mul(wrongReadings[1]).sub(g12.mul(wrongReadings[0])).div(fixture.reading.delta),wrongRat=wrongBeta.div(wrongAlpha);
control('R7:nonrest_speed_change_can_change_inferred_ratio',wrongAlpha.eq(8)&&wrongBeta.eq(5)&&wrongRat.eq('5/8')&&!wrongRat.eq(targetRatio));

const output={schema:'native-coupling-reconstruction.v1',status:'PASS_NATIVE_COUPLING_RECONSTRUCTION',source_commit:pins.source_commit,
  scope:'Exact symmetric-action blindness; native Gram recovery of information/induced-curvature coupling ratio from initial-rest covariant cut acceleration; regular parameter and frame invariance; factorized control-design criterion; independent-preparation diagnostic. Synthetic fixtures are not physical measurements or universal constants.',
  source_pins_sha256:hash(fs.readFileSync(path.join(root,'SOURCE_PINS.json'))),native_engine:pins.native_engine,upstream,presentations,certificates,checks,negative_controls:controls,
  examples:{fixture,symmetric_blindness:blindness,parameter_changes:observed,frame_changes:frames,independent_predictions:predictions,nonrest_counterexample:{V,X,changedForce,wrongReadings,wrongAlpha,wrongBeta,correct_ratio:targetRatio,apparent_ratio:wrongRat}},
  summary:{symbolic_certificates:certificates.length,serialized_replays:certificates.length,native_presentations:presentations.length,exact_checks:checks.length,negative_controls:controls.length,parameter_cases:observed.length,frame_cases:frames.length,independent_preparations:predictions.length}};
const serialized=JSON.stringify(output,null,2)+'\n',digest=hash(serialized),file=path.join(root,'CERTIFICATE.json'),checksum=path.join(root,'EXPECTED.sha256');
if(write){fs.writeFileSync(file,serialized);fs.writeFileSync(checksum,digest+'  CERTIFICATE.json\n');}
else{assert.equal(fs.readFileSync(file,'utf8'),serialized,'Frozen certificate differs from native replay');assert.equal(fs.readFileSync(checksum,'utf8'),digest+'  CERTIFICATE.json\n','Certificate checksum mismatch');}
console.log(JSON.stringify({status:output.status,mode:write?'explicit generation':'read-only replay',...output.summary,certificate_sha256:digest},null,2));
