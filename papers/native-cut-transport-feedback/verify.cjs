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

// Native word certificates. Trace conclusions use the displayed commutators.
{
  const tokens=['p','v','h'],rules=[
    {id:'cut',lhs:['p','p'],rhs:[[['p'],1]],source:'P squared equals P'},
    {id:'horizontal',lhs:['v','p'],rhs:[[['v'],1],[['p','v'],-1]],source:'PV+VP=V'}];
  const {S,w:{p,v,h}}=make('E1_cut_and_transport',tokens,{rules,star:(_,w)=>({p:w.p,v:w.v.scale(-1),h:w.h})});
  const br=M.bracket,x=br(v,p),force=br(p,h),one=S.one();
  prove('E1:cut_constraint_derivative',x.times(p).plus(p.times(x)).minus(x),S.zero());
  prove('E1:horizontal_constraint_derivative',x.times(v).plus(v.times(x)).plus(p.times(force)).plus(force.times(p)).minus(force),S.zero());
  prove('E1:recover_generator',br(x,p),v);
  prove('E1:kinetic_operator_identity',sq(x),sq(v).scale(-1));
  prove('E1:kinetic_trace_equals_cut_information',sq(v).scale(-1).plus(times(p,v,v,p).scale(2)),br(p.times(v),v.times(p)));
  const accel=br(force,p).plus(br(v,x));
  prove('E1:stationary_cut_equation',br(p,accel.plus(h)),S.zero());
  prove('E1:normal_acceleration_commutes_with_cut',br(p,br(v,x)),S.zero());
  prove('E1:force_is_cut_off_diagonal',p.times(force).plus(force.times(p)),force);
}
{
  const {S,w:{a,b,c,p,v,h}}=make('E2_cyclic_variation',['a','b','c','p','v','h']);
  const br=M.bracket;
  prove('E2:trace_integration_by_parts_residual',a.times(br(b,c)).minus(br(a,b).times(c)),br(b,a.times(c)));
  prove('E2:energy_derivative_is_sum_of_commutators',v.times(br(p,h)).scale(-1).plus(h.times(br(v,p))),br(h,v.times(p)).plus(br(v,h.times(p))));
  prove('E2:force_covariance_jacobi',br(a,br(p,h)),br(br(a,p),h).plus(br(p,br(a,h))));
}

// All algebra below uses native Expr coefficients and native matrix constants.
// Central sphere relation is certified by checked critical pairs and dagger.
{
  const {S,w:{x,y,z}}=make('E3_rank_one_two_mode',['x','y','z'],{commutative:true,rules:[
    {id:'unit_sphere',lhs:['z','z'],rhs:[[[],1],[['x','x'],-1],[['y','y'],-1]],source:'Native rank-one cut constraint x squared+y squared+z squared=1'}]});
  const raw=sym(S),normalize=A=>A.map(row=>row.map(x=>S.reduce(x).normal));
  const m={...raw,add:(A,B)=>normalize(raw.add(A,B)),sub:(A,B)=>normalize(raw.sub(A,B)),mul:(A,B)=>normalize(raw.mul(A,B)),tr:A=>S.reduce(raw.tr(A)).normal};
  const br=(A,B)=>m.sub(m.mul(A,B),m.mul(B,A)),sm=(A,c)=>A.map(r=>r.map(a=>a.scale(c)));
  const sx=N.matrix([[0,1],[1,0]]),sy=N.scale(R,IOTA),sz=N.matrix([[1,0],[0,-1]]);
  const P=sm(m.add(m.id,m.add(m.scale(m.constant(sx),x),m.add(m.scale(m.constant(sy),y),m.scale(m.constant(sz),z)))),'1/2');
  const B=basis.map(m.constant),X=B.map(b=>br(b,P)),weights=[1,1,2],pairs=[[0,1,1],[0,2,2],[1,2,2]];
  let HI=m.zero,HC=m.zero,phiI=S.zero(),phiC=S.zero();
  X.forEach((q,a)=>{HI=m.add(HI,sm(br(B[a],q),-2*weights[a]));phiI=phiI.plus(m.tr(m.mul(q,q)).scale(weights[a]));});
  for(const [a,b,w] of pairs){
    const W=br(X[a],X[b]),C=m.mul(m.mul(P,W),P);
    const H=m.add(sm(m.add(m.mul(W,C),m.mul(C,W)),-1),m.add(br(B[a],br(X[b],C)),br(B[b],br(C,X[a]))));
    HC=m.add(HC,sm(H,w));phiC=phiC.minus(m.tr(m.mul(C,C)).scale(f(w).div(2)));
  }
  m.eq('E3:rank_one_cut',m.mul(P,P),P);
  prove('E3:information_functional',phiI,S.one().scale(6).minus(sq(z).scale(2)));
  prove('E3:curvature_functional',phiC,S.one().scale(4).minus(sq(z).scale(2)));
  const axis=br(P,br(B[2],br(B[2],P)));
  m.eq('E3:information_force',br(P,HI),sm(axis,-2));
  m.eq('E3:curvature_force',br(P,HC),sm(axis,-2));
  const W=sm(B[2],2),V=br(br(W,P),P),dP=br(W,P),dV=br(W,V);
  m.eq('E3:closed_cut_orbit',br(V,P),dP);
  m.eq('E3:closed_generator_orbit',br(P,m.add(HI,HC)),dV);
  prove('E3:transport_information',m.tr(m.mul(V,V)).scale('-1/2'),S.one().minus(sq(z)).scale(4));
  prove('E3:conserved_orbit_total',m.tr(m.mul(V,V)).scale('-1/2').plus(phiI).plus(phiC),S.one().scale(14).minus(sq(z).scale(8)));
}

// The returning pair has a nontrivial fifth-root full-carrier transport.
{
  const tokens=['p','q'],rules=[
    {id:'cut',lhs:['p','p'],rhs:[[['p'],1]],source:'Retained cut'},
    {id:'phase_is_central',lhs:['q','p'],rhs:[[['p','q'],1]],source:'Native scalar phase'},
    {id:'fifth_root',lhs:['q','q','q','q'],rhs:[[[],-1],[['q'],-1],[['q','q'],-1],[['q','q','q'],-1]],source:'Primitive fifth-root polynomial 1+q+q^2+q^3+q^4=0'}];
  const {S,w:{p,q}}=make('E4_return_transport',['p','q'],{rules,star:(S,w)=>({p:w.p,q:times(w.q,w.q,w.q,w.q)})});
  const qp=times(q,q,q,q),U=qp.times(p).plus(q.times(S.one().minus(p))),Ud=q.times(p).plus(qp.times(S.one().minus(p)));
  prove('E4:phase_order_five',times(q,q,q,q,q),S.one());
  prove('E4:transport_is_unit',Ud.times(U),S.one());
  prove('E4:visible_cut_return',times(U,p,Ud),p);
  prove('E4:five_cycle_transport_return',times(U,U,U,U,U),S.one());
  let power=S.one();for(let n=1;n<5;n++){power=S.reduce(power.times(U)).normal;const residual=S.reduce(power.minus(S.one())).normal;control('E4:pair_return_does_not_close_transport_cycle_'+n,residual.terms.size>0,{nonzero_native_normal_form:residual.toJSON()});}
}

const sum=xs=>xs.reduce((s,x)=>s.add(x),f(0)),br=N.commutator;
const daggerEq=A=>N.equal(N.dagger(A),A),antiEq=A=>N.equal(N.dagger(A),N.scale(A,-1));
const trReal=A=>real(trace(A)),Pstar=N.matrix([['4/5','2/5'],['2/5','1/5']]);
const transform=(G,A)=>mm(N.dagger(G),A,G);
const phiWeights={alpha:f(2),beta:f(1),w:[f(1),f(1),f(2)],pairs:[[0,1,f(1)],[0,2,f(2)],[1,2,f(2)]]};
function functional(P,B,weights=phiWeights){
  const {alpha,beta,w,pairs}=weights,X=B.map(b=>br(b,P));
  let HI=O,HC=O,phiI=f(0),phiC=f(0);const curvature=[];
  X.forEach((x,a)=>{HI=N.sub(HI,N.scale(br(B[a],x),alpha.mul(w[a])));phiI=phiI.add(N.energy(x).mul(alpha).mul(w[a]).div(2));});
  pairs.forEach(([a,b,h])=>{
    const W=br(X[a],X[b]),C=mm(P,W,P),hC=N.add(N.scale(N.add(mm(W,C),mm(C,W)),-1),N.add(br(B[a],br(X[b],C)),br(B[b],br(C,X[a]))));
    HC=N.add(HC,N.scale(hC,beta.mul(h)));phiC=phiC.add(N.energy(C).mul(beta).mul(h).div(2));curvature.push({a,b,h,W,C});
  });
  assert.ok(daggerEq(HI)&&daggerEq(HC));
  const H=N.add(HI,HC),Q=X.map(a=>X.map(b=>trace(mm(P,a,b,P))));
  return {P,X,HI,HC,H,phiI,phiC,phi:phiI.add(phiC),curvature,Q};
}
function directVariation(P,B,dP,weights=phiWeights){
  // Differentiate the defining functional by native Leibniz products,
  // independently of the closed gradient formula.
  const r=functional(P,B,weights),dX=B.map(b=>br(b,dP)),{alpha,beta,w}=weights;
  let dI=f(0),dC=f(0);
  r.X.forEach((x,a)=>{dI=dI.add(trReal(mm(x,dX[a])).mul(alpha).mul(w[a]));});
  r.curvature.forEach(({a,b,h,W,C})=>{
    const dW=N.add(br(dX[a],r.X[b]),br(r.X[a],dX[b]));
    const dZ=N.add(N.add(mm(dP,W,P),mm(P,W,dP)),mm(P,dW,P));
    dC=dC.sub(trReal(mm(C,dZ)).mul(beta).mul(h));
  });
  return {dI,dC,total:dI.add(dC),gradient:trReal(mm(r.H,dP))};
}
function coupled(P,V,B,weights=phiWeights){
  assert.ok(daggerEq(P)&&N.equal(mm(P,P),P)&&antiEq(V)&&N.equal(N.add(mm(P,V),mm(V,P)),V));
  const r=functional(P,B,weights),dP=br(V,P),dV=br(P,r.H),dv=directVariation(P,B,dP,weights);
  const dT=trReal(mm(V,dV)).neg(),kinetic=N.energy(V).div(2);
  assert.ok(dv.total.eq(dv.gradient));assert.ok(dT.add(dv.total).zero());
  const dX=B.map(b=>br(b,dP));
  const dQ=r.X.map((xa,a)=>r.X.map((xb,b)=>trace(mm(P,N.add(mm(dX[a],xb),mm(xa,dX[b])),P))));
  const directQ=r.X.map((xa,a)=>r.X.map((xb,b)=>trace(N.add(N.add(mm(dP,xa,xb,P),mm(P,N.add(mm(dX[a],xb),mm(xa,dX[b])),P)),mm(P,xa,xb,dP)))));
  assert.ok(dQ.every((row,a)=>row.every((x,b)=>x.eq(directQ[a][b]))));
  return {...r,V,dP,dV,dT,dPhi:dv.total,kinetic,energy:kinetic.add(r.phi),tensor_derivative:dQ};
}
const phases=['0','1/7','1/3','1/2','1','-2/3'],orbit=[],gradientCases=[],frameCases=[];
for(const t of phases){
  const G=N.cayley(L,t),P=mm(G,Pstar,N.dagger(G)),W=N.scale(L,2),V=br(br(W,P),P),r=coupled(P,V,basis),tag=t;
  check('E5:exact_orbit:'+tag,N.equal(r.dP,br(W,P))&&N.equal(r.dV,br(W,V))&&!N.isZero(r.dV));
  check('E5:constant_selected_totals:'+tag,r.kinetic.eq('64/25')&&r.phiI.eq('132/25')&&r.phiC.eq('82/25')&&r.energy.eq('278/25'));
  check('E5:quarter_turn_between_cut_and_generator:'+tag,V[0][1].eq(P[0][1].mul(new Cut(0,'12/5'))));
  check('E5:both_feedback_terms_nonzero:'+tag,!N.isZero(br(P,r.HI))&&N.equal(br(P,r.HI),br(P,r.HC)));
  const hermDirections=[N.matrix([[1,0],[0,-1]]),N.matrix([[0,1],[1,0]]),N.scale(R,IOTA),N.add(P,N.matrix([[0,'1/3'],['1/3',0]]))];
  hermDirections.forEach((dP,j)=>{const v=directVariation(P,basis,dP);check('E5:ambient_gradient_variation:'+tag+':'+j,v.total.eq(v.gradient));gradientCases.push({phase:t,direction:j,...v});});
  const frame=N.mul(N.cayley(R,'1/3'),N.cayley(J,t)),moved=coupled(transform(frame,P),transform(frame,V),basis.map(b=>transform(frame,b)));
  check('E5:frame_covariance:'+tag,moved.energy.eq(r.energy)&&N.equal(moved.dP,transform(frame,r.dP))&&N.equal(moved.dV,transform(frame,r.dV))&&moved.Q.every((row,a)=>row.every((x,b)=>x.eq(r.Q[a][b]))));
  const connection=N.add(R,N.scale(J,'2/7')),ordinaryP=N.sub(moved.dP,br(connection,moved.P)),ordinaryV=N.sub(moved.dV,br(connection,moved.V));
  check('E5:variable_frame_correction:'+tag,N.equal(N.add(ordinaryP,br(connection,moved.P)),moved.dP)&&N.equal(N.add(ordinaryV,br(connection,moved.V)),moved.dV));
  frameCases.push({phase:t,energy:moved.energy,connection});orbit.push({phase:t,...r});
}
check('E5:information_components_really_evolve',orbit.some(r=>r.tensor_derivative.some(row=>row.some(x=>!x.rad.zero()))));
check('E5:curvature_components_really_evolve',orbit.some(r=>r.tensor_derivative.some(row=>row.some(x=>!x.turn.zero()))));
const force=functional(Pstar,basis),exchangeV=br(Pstar,force.H),exchange=coupled(Pstar,exchangeV,basis);
check('E6:nonzero_information_curvature_energy_exchange',!exchange.dPhi.zero()&&exchange.dPhi.eq(exchange.dT.neg()));
const mixedControls=[N.add(R,N.scale(J,'1/3')),N.add(J,N.scale(L,'2/5')),N.add(L,N.scale(R,'-1/4'))];
const mixed=coupled(Pstar,br(br(N.add(R,J),Pstar),Pstar),mixedControls,{alpha:f('3/2'),beta:f('2/3'),w:[f(2),f(1),f('1/2')],pairs:[[0,1,f('1/3')],[1,2,f('5/4')]]});
check('E6:nonaxial_model_conserves_total',mixed.dT.add(mixed.dPhi).zero());
control('E6:freezing_generator_violates_feedback',!N.isZero(orbit[0].dV));
control('E6:force_sign_is_load_bearing',!exchange.dT.neg().add(exchange.dPhi).zero());
control('E6:kinetic_term_alone_is_not_conserved',!exchange.dT.zero());
control('E6:curvature_force_cannot_be_omitted',!N.equal(br(Pstar,force.HI),br(Pstar,force.H)));
control('E6:frame_change_requires_control_change',!functional(transform(N.cayley(R,'1/3'),Pstar),basis).phi.eq(force.phi));
control('E6:cut_diagonal_generator_freedom_is_fixed',rejects(()=>coupled(Pstar,N.add(orbit[0].V,N.scale(N.sub(N.scale(Pstar,2),I),IOTA)),basis)));
const staticCut=N.matrix([[1,0],[0,0]]),staticState=coupled(staticCut,O,basis);
control('E6:coupling_does_not_force_every_orbit_to_be_helical',N.isZero(staticState.dP)&&N.isZero(staticState.dV));
check('E6:reversal_symmetry',N.equal(coupled(Pstar,N.scale(orbit[0].V,-1),basis).dP,N.scale(orbit[0].dP,-1))&&N.equal(coupled(Pstar,N.scale(orbit[0].V,-1),basis).dV,orbit[0].dV));

const output={schema:'native-cut-transport-feedback.v1',status:'PASS_NATIVE_CUT_TRANSPORT_FEEDBACK',source_commit:pins.source_commit,
  scope:'Selected finite positive cut action derives coupled cut/generator evolution, exact information-curvature force, conserved total, co-rotating example and five-cycle transport memory. Native algebra is executed; analytic variational/ordered-exponential arguments have written proofs. No universal action or gravitational law is inferred.',
  source_pins_sha256:hash(fs.readFileSync(path.join(root,'SOURCE_PINS.json'))),native_engine:pins.native_engine,upstream,presentations,certificates,checks,negative_controls:controls,
  examples:{orbit,gradient_variations:gradientCases,frames:frameCases,nonzero_exchange:exchange,nonaxial_case:mixed,static_case:staticState},
  summary:{symbolic_certificates:certificates.length,serialized_replays:certificates.length,native_presentations:presentations.length,exact_checks:checks.length,negative_controls:controls.length,orbit_cases:orbit.length,gradient_variations:gradientCases.length,frame_cases:frameCases.length}};
const serialized=JSON.stringify(output,null,2)+'\n',digest=hash(serialized),file=path.join(root,'CERTIFICATE.json'),checksum=path.join(root,'EXPECTED.sha256');
if(write){fs.writeFileSync(file,serialized);fs.writeFileSync(checksum,digest+'  CERTIFICATE.json\n');}
else{assert.equal(fs.readFileSync(file,'utf8'),serialized,'Frozen certificate differs from native replay');assert.equal(fs.readFileSync(checksum,'utf8'),digest+'  CERTIFICATE.json\n','Certificate checksum mismatch');}
console.log(JSON.stringify({status:output.status,mode:write?'explicit generation':'read-only replay',...output.summary,certificate_sha256:digest},null,2));
