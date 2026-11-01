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

// F1: change of observation base, followed by reference-correct composition.
{
  const names=['g','h','t','p','v1','v2','v3','u1','u2','u3'];
  const tokens=names.flatMap(x=>[x,x+'d']),rules=names.flatMap(x=>[
    {id:x+'_right',lhs:[x,x+'d'],rhs:[[[],1]],source:'Native two-sided dagger inverse'},
    {id:x+'_left',lhs:[x+'d',x],rhs:[[[],1]],source:'Native two-sided dagger inverse'}]);
  const {S,w}=make('F1_based_unit_words',tokens,{rules,star:(_,w)=>Object.fromEntries(names.flatMap(x=>[[x,w[x+'d']],[x+'d',w[x]]]))});
  const T=times(w.hd,w.t,w.g),Td=times(w.gd,w.td,w.h),P=times(w.hd,w.p,w.h);
  prove('F1:bridge_frame_covariance',times(Td,P,T),times(w.gd,w.td,w.p,w.t,w.g));
  prove('F1:bridge_is_reversible',times(w.t,w.td,w.p,w.t,w.td),w.p);
  const q1=times(w.v1d,w.u1),q2=times(w.v1d,w.v2d,w.u2,w.v1);
  const q3=times(w.v1d,w.v2d,w.v3d,w.u3,w.v2,w.v1);
  prove('F1:two_reference_prefixes',times(q2,q1),times(w.v1d,w.v2d,w.u2,w.u1));
  prove('F1:three_reference_prefixes',times(q3,q2,q1),times(w.v1d,w.v2d,w.v3d,w.u3,w.u2,w.u1));
}

// F2: three-loop scalar reconstruction, with an oriented third-order reading.
{
  const {S,w}=make('F2_three_native_vectors',['a','b','c','x','y','z','u','v','w','r','s','t'],{commutative:true});
  const {a,b,c,x,y,z,u,v,w:wv,r,s,t}=w,m=sym(S),vs=[[x,y,z],[u,v,wv],[r,s,t]];
  const dot=(A,B)=>A.reduce((sum,q,i)=>sum.plus(q.times(B[i])),S.zero());
  const det=times(x,v,t).plus(times(y,wv,r)).plus(times(z,u,s)).minus(times(z,v,r)).minus(times(y,u,t)).minus(times(x,wv,s));
  const V=vs.map(m.vector),Q=V.map((v,i)=>m.add(m.scale(m.id,[a,b,c][i]),v));
  const g=vs.map(v=>vs.map(u=>dot(v,u)));
  const gramDet=times(g[0][0],g[1][1],g[2][2]).plus(times(g[0][1],g[1][2],g[2][0]).scale(2))
    .minus(g[0][0].times(sq(g[1][2]))).minus(g[1][1].times(sq(g[0][2]))).minus(g[2][2].times(sq(g[0][1])));
  prove('F2:oriented_native_triple',m.tr(m.mul(V[2],m.mul(V[1],V[0]))).scale('1/2'),det);
  prove('F2:triple_orientation_reversal',m.tr(m.mul(V[0],m.mul(V[1],V[2]))).scale('1/2'),det.scale(-1));
  prove('F2:three_loop_scalar_decoder',m.tr(m.mul(Q[2],m.mul(Q[1],Q[0]))).scale('1/2'),times(a,b,c).minus(c.times(g[0][1])).minus(b.times(g[0][2])).minus(a.times(g[1][2])).plus(det));
  prove('F2:orientation_square_equals_Gram',sq(det),gramDet);
}
// F3/F4: native telescoping, positive square certificate, and budget merging.
{
  const {S,w:{p,q,r}}=make('F3_ordered_telescope',['p','q','r']);
  prove('F3:three_loop_telescope',S.one().minus(times(r,q,p)),S.one().minus(r).plus(r.times(S.one().minus(q))).plus(times(r,q,S.one().minus(p))));
}
{
  const {w:{x,y,z}}=make('F4_positive_records',['x','y','z'],{commutative:true});
  prove('F4:three_term_sum_of_squares',sq(x).plus(sq(y)).plus(sq(z)).scale(3).minus(sq(x.plus(y).plus(z))),sq(x.minus(y)).plus(sq(y.minus(z))).plus(sq(z.minus(x))));
}
// F5: the fused small-cell curvature is a sum; retained records retain squares.
{
  const names=['A','B','C'],tokens=['e',...names],rules=[{id:'e2',lhs:['e','e'],rhs:[],source:'First native cell-area jet e^2=0'}];
  for(const x of names)rules.push({id:x+'_past_e',lhs:[x,'e'],rhs:[[['e',x],1]],source:'Central area increment'});
  const {S,w:{e,A,B,C}}=make('F5_fused_curvature_jet',tokens,{rules});
  prove('F5:three_curvature_sum',times(S.one().minus(e.times(C)),S.one().minus(e.times(B)),S.one().minus(e.times(A))),S.one().minus(e.times(A.plus(B).plus(C))));
}

// CID-D1: cut tangency, induced curvature and the joint tensor evolution.
{
  const tangents=['x','y','z'],tokens=['p',...tangents,'h','j','f'];
  const rules=[{id:'cut_idempotent',lhs:['p','p'],rhs:[[['p'],1]],source:'Native observation cut P^2=P'},
    ...tangents.map(x=>({id:x+'_cut_tangent',lhs:[x,'p'],rhs:[[[x],1],[['p',x],-1]],source:'D(P^2)=DP: XP+PX=X'}))];
  const {S,w:{p,x,y,z,h,j,f:Fcurv}}=make('D1_native_moving_cut',tokens,{rules,star:(_,w)=>Object.fromEntries(tokens.map(t=>[t,t==='f'?w[t].scale(-1):w[t]]))});
  const one=S.one(),K=p.scale(2).minus(one),br=M.bracket,Ax=br(p,x),Ay=br(p,y);
  prove('D1:cut_involution',sq(K),one);
  prove('D1:tangent_has_no_retained_diagonal',times(p,x,p),S.zero());
  prove('D1:recover_minimal_cut_generator',br(br(x,p),p),x);
  prove('D2:minimal_connection_parallelizes_cut',x.plus(br(Ax,p)),S.zero());
  const full=Fcurv.plus(br(x,y).scale(2)).plus(br(p,br(Fcurv,p))).plus(br(Ax,Ay));
  const expected=times(p,Fcurv,p).plus(times(one.minus(p),Fcurv,one.minus(p))).plus(br(x,y));
  prove('D2:full_cut_connection_curvature',full,expected);
  prove('D2:cut_connection_curvature_preserves_cut',br(expected,p),S.zero());
  prove('D2:induced_retained_curvature',times(p,expected,p),times(p,Fcurv,p).plus(times(p,br(x,y),p)));
  prove('D2:pure_cut_score_equation',p.times(x.scale(2)).plus(x.scale(2).times(p)),x.scale(2));
  const derivative=times(z,x,y,p).plus(times(p,h,y,p)).plus(times(p,x,j,p)).plus(times(p,x,y,z));
  prove('D3:covariant_tensor_evolution_corner',times(p,derivative,p),p.times(h.times(y).plus(x.times(j))).times(p));
}
{
  const {w:{a,b,p,u,v,c,d,w}}=make('D3_covariant_derivative_commutator',['a','b','p','u','v','c','d','w']);
  const br=M.bracket,first=w.plus(br(c,p)).plus(br(b,u)).plus(br(a,v.plus(br(b,p))));
  const second=w.plus(br(d,p)).plus(br(a,v)).plus(br(b,u.plus(br(a,p))));
  prove('D3:curvature_correction_to_cut_acceleration',first.minus(second),br(c.minus(d).plus(br(a,b)),p));
}
{
  const {w:{x,y,u,v}}=make('D4_information_curvature_area',['x','y','u','v'],{commutative:true});
  prove('D4:two_mode_information_curvature_duality',sq(x).plus(sq(y)).times(sq(u).plus(sq(v))).minus(sq(x.times(u).plus(y.times(v)))),sq(x.times(v).minus(y.times(u))));
}

function unit(Q){
  assert.ok(N.equal(mm(N.dagger(Q),Q),I),'Unit native transport required');
  const a=sc(Q),v=basis.map(B=>pair(B,Q));
  assert.ok(N.equal(Q,N.add(N.scale(I,a),combo(v))),'Real native unit sector required');
  assert.ok(v.reduce((s,x)=>s.add(x.pow(2)),a.pow(2)).eq(1));return {a,v};
}
const product=Qs=>Qs.reduce((p,q)=>N.mul(q,p),I),conj=(T,P)=>mm(N.dagger(T),P,T);
function defect(Q){unit(Q);const c=f(1).sub(sc(Q)).div(2);assert.ok(f(0).le(c)&&c.le(1));assert.ok(c.eq(N.energy(N.sub(I,Q)).div(8)));return c;}
const sum=xs=>xs.reduce((s,x)=>s.add(x),f(0));
const dot=(a,b)=>sum(a.map((x,i)=>x.mul(b[i])));
const det3=g=>g[0][0].mul(g[1][1]).mul(g[2][2]).add(g[0][1].mul(g[1][2]).mul(g[2][0])).add(g[0][2].mul(g[1][0]).mul(g[2][1]))
  .sub(g[0][2].mul(g[1][1]).mul(g[2][0])).sub(g[0][1].mul(g[1][0]).mul(g[2][2])).sub(g[0][0].mul(g[1][2]).mul(g[2][1]));
function triple(Qs){
  assert.equal(Qs.length,3);const parts=Qs.map(unit),a=parts.map(x=>x.a),V=Qs.map((Q,i)=>N.sub(Q,N.scale(I,a[i])));
  const g=parts.map(x=>parts.map(y=>dot(x.v,y.v))),tau=sc(product(V));
  const cross=[[0,1],[0,2],[1,2]].map(([i,j])=>sc(mm(N.dagger(Qs[i]),Qs[j])));
  const scalar=a[0].mul(a[1]).mul(a[2]).sub(a[2].mul(g[0][1])).sub(a[1].mul(g[0][2])).sub(a[0].mul(g[1][2])).add(tau);
  assert.ok(scalar.eq(sc(product(Qs))),'Oriented decoder disagrees with direct native product');
  assert.ok(tau.pow(2).eq(det3(g)),'Gram/orientation identity failed');
  const decoded=a[0].mul(a[1]).mul(a[2]).sub(a[2].mul(cross[0].sub(a[0].mul(a[1])))).sub(a[1].mul(cross[1].sub(a[0].mul(a[2])))).sub(a[0].mul(cross[2].sub(a[1].mul(a[2])))).add(tau);
  assert.ok(decoded.eq(scalar));return {scalar_readings:a,cross_readings:cross,tau,gram:g,fused_scalar:scalar,fused_defect:f(1).sub(scalar).div(2),retained_defect:sum(Qs.map(defect))};
}
function fusion(U,V,T){
  assert.ok(U.length>0&&U.length===V.length&&V.length===T.length,'Matched nonempty source/reference/bridge records required');
  U.concat(V,T).forEach(unit);
  const u=U.map((x,i)=>conj(T[i],x)),v=V.map((x,i)=>conj(T[i],x));
  let prefix=I;const local=[],dressed=[];
  for(let i=0;i<u.length;i++){const p=mm(N.dagger(v[i]),u[i]);local.push(p);dressed.push(conj(prefix,p));prefix=mm(v[i],prefix);}
  const direct=mm(N.dagger(product(v)),product(u)),fused=product(dressed);
  assert.ok(N.equal(direct,fused),'Reference-correct ordered fusion failed');
  assert.ok(sum(local.map(defect)).eq(sum(dressed.map(defect))));
  return {local,dressed,direct,fused,retained_defect:sum(local.map(defect)),fused_defect:defect(fused),naive:product(local)};
}
function positiveBound(Qs){
  assert.ok(Qs.length>0);Qs.forEach(unit);let suffix=I;const d=[];
  for(let i=Qs.length-1;i>=0;i--){d.push(mm(suffix,N.sub(I,Qs[i])));suffix=mm(suffix,Qs[i]);}
  const total=d.reduce(N.add,O);assert.ok(N.equal(total,N.sub(I,product(Qs))));
  const residual=sum(d.flatMap((x,i)=>d.slice(i+1).map(y=>N.energy(N.sub(x,y)))));
  const retained=sum(Qs.map(defect)),fused=defect(product(Qs));
  assert.ok(retained.mul(Qs.length).sub(fused).mul(8).eq(residual));
  return {count:Qs.length,retained,fused,positive_residual:residual};
}
function record(Qs,budgets){
  assert.ok(Qs.length>0&&Qs.length===budgets.length,'Matched nonempty records required');
  budgets=budgets.map(f);assert.ok(budgets.every(x=>!x.le(0)),'Positive native squared areas required');
  const defect_sum=sum(Qs.map(defect)),budget_sum=sum(budgets),density=defect_sum.div(budget_sum);
  const child_density=Qs.map((q,i)=>defect(q).div(budgets[i]));
  const weights=budgets.map(b=>b.div(budget_sum));assert.ok(density.eq(sum(weights.map((w,i)=>w.mul(child_density[i])))));
  return {defect_sum,budget_sum,density,child_density,weights};
}
function blockRecord(Qs){const b=N.zeros(2*Qs.length);Qs.forEach((Q,k)=>Q.forEach((r,i)=>r.forEach((x,j)=>b[2*k+i][2*k+j]=x)));return b;}
function cell(P,t){
  const x=N.cayley(R,f(t).div(-2)),y=N.cayley(J,f(t).div(-2)),top=mm(x,N.dagger(P),N.dagger(y));
  const loop=mm(N.dagger(y),N.dagger(top),I,x);assert.ok(N.equal(loop,P));
  const readout=Q=>N.scale(mm(N.sub(I,Q),N.inverse(N.add(I,Q))),2);
  const dx=readout(x),dy=readout(y),area=pair(dx,dx).mul(pair(dy,dy)).sub(pair(dx,dy).pow(2));
  assert.ok(area.eq(f(t).pow(4)));return {links:{x,y,r:I,t:top},area_squared:area,shape:f(1),defect:defect(P)};
}

const Gs=[N.cayley(combo([1,2,-1]),'1/5'),mm(N.cayley(R,'1/3'),N.cayley(J,'2/5'))];
const tripleCases=[];let frameCases=0;
for(const x of ['1/3','1/2','-2/3'])for(const y of ['1/4','2/5'])for(const z of ['1/5','-3/7']){
  const Qs=[N.cayley(R,x),N.cayley(J,y),N.cayley(L,z)],r=triple(Qs),tag=[x,y,z].join(',');
  check('F2:exact_three_loop_reconstruction:'+tag,r.fused_defect.eq(defect(product(Qs))));
  Gs.forEach((G,i)=>{const s=triple(Qs.map(q=>conj(G,q)));check('F2:common_frame_invariance:'+tag+':'+i,s.tau.eq(r.tau)&&s.fused_defect.eq(r.fused_defect)&&s.cross_readings.every((x,j)=>x.eq(r.cross_readings[j])));frameCases++;});
  tripleCases.push({parameters:[x,y,z],readout:r,bound:positiveBound(Qs)});
}
const planar=triple([N.cayley(R,'1/3'),N.cayley(J,'1/4'),N.cayley(N.add(R,J),'1/5')]);
check('F2:planar_triple_needs_no_orientation_sign',planar.tau.zero());
const positive=triple([R,J,L]),negative=triple([R,J,N.scale(L,-1)]);
control('F2:six_pairwise_scalars_do_not_determine_fusion',positive.scalar_readings.every((x,i)=>x.eq(negative.scalar_readings[i]))&&positive.cross_readings.every((x,i)=>x.eq(negative.cross_readings[i]))&&positive.tau.eq(1)&&negative.tau.eq(-1)&&positive.fused_defect.zero()&&negative.fused_defect.eq(1));
control('F2:order_is_observation_data',defect(product([R,J,L])).zero()&&defect(product([L,J,R])).eq(1));
control('F1:connector_choice_is_not_a_frame_symmetry',defect(product([R,J,conj(R,L)])).eq(1)&&defect(product([R,J,L])).zero());

const source=[N.cayley(R,'1/3'),N.cayley(J,'2/5'),N.cayley(L,'-1/4')];
const reference=[N.cayley(J,'1/5'),N.cayley(L,'-1/6'),N.cayley(R,'2/7')];
const bridges=[N.cayley(L,'1/4'),N.cayley(R,'1/7'),N.cayley(J,'-1/8')];
const corrected=fusion(source,reference,bridges);
check('F1:nonflat_reference_fusion',N.equal(corrected.direct,corrected.fused));
control('F1:omitting_reference_prefixes_changes_the_answer',!N.equal(corrected.naive,corrected.direct)&&!defect(corrected.naive).eq(corrected.fused_defect));
Gs.forEach((G,i)=>{
  const sites=[N.cayley(R,'1/9'),N.cayley(J,'-2/9'),N.cayley(L,'1/6')];
  const moved=fusion(source.map((q,j)=>conj(sites[j],q)),reference.map((q,j)=>conj(sites[j],q)),bridges.map((q,j)=>mm(N.dagger(sites[j]),q,G)));
  check('F1:independent_site_and_base_frames:'+i,N.equal(moved.fused,conj(G,corrected.fused))&&moved.retained_defect.eq(corrected.retained_defect));frameCases++;
});

const finiteBounds=[];
for(let m=1;m<=7;m++){
  const Qs=Array.from({length:m},(_,i)=>N.cayley(basis[i%3],f(i+1).div(m+3)));
  const r=positiveBound(Qs);check('F3:finite_telescoping_square_bound:'+m,f(0).le(r.positive_residual)&&r.fused.le(r.retained.mul(m)));finiteBounds.push(r);
  const b=blockRecord(Qs);check('F4:labelled_block_retains_all_defects:'+m,N.energy(N.sub(N.identity(2*m),b)).div(8).eq(r.retained));
}
const merging=record(source,['1/4','1/9','1/16']);
check('F4:budget_merge_associative',sum(source.slice(0,2).map(defect)).add(defect(source[2])).eq(merging.defect_sum)&&f('1/4').add('1/9').add('1/16').eq(merging.budget_sum));
check('F4:weighted_mean_stays_in_child_range',merging.child_density.some(x=>x.le(merging.density))&&merging.child_density.some(x=>merging.density.le(x)));
check('F4:zero_record_means_all_unit_returns',record([I,I,I],[1,1,1]).defect_sum.zero()&&!record([I,R,I],[1,1,1]).defect_sum.zero());

const cancellation=[],alignment=[],population=[];
for(let n=1;n<=6;n++){
  const t=f(1).div(f(2).pow(n)),t2=t.pow(2),t4=t.pow(4),q=N.cayley(R,t.neg()),r=N.cayley(J,t.neg());
  const Qs=[q,r,N.dagger(mm(r,q))],data=record(Qs,[t4,t4,t4]),bound=positiveBound(Qs),cells=Qs.map(p=>cell(p,t));
  const retained_density=data.defect_sum.div(t4),formula=f(2).mul(f(2).add(t2)).div(t2.mul(f(1).add(t2).pow(2)));
  check('F6:cancellation_family_exact:'+n,N.equal(product(Qs),I)&&data.defect_sum.eq(f(2).mul(t2).mul(f(2).add(t2)).div(f(1).add(t2).pow(2)))&&retained_density.eq(formula));
  check('F6:cancellation_family_nonzero_leaves:'+n,cells.every(c=>!c.defect.le(0)&&c.area_squared.eq(t4)&&c.shape.eq(1)));
  check('F6:retained_concentration_lower_bound:'+n,f('64/25').div(t2).le(retained_density));
  check('F6:budget_normalized_density_retains_concentration:'+n,data.density.eq(retained_density.div(3)));
  cancellation.push({n,t,retained_density,normalized_density:data.density,record:data,bound,cells});
  const aligned=[q,q,q],al=positiveBound(aligned),efficiency=al.fused.div(al.retained.mul(3));
  const expected=t2.mul(f(3).sub(t2).pow(2)).div(f(1).add(t2).pow(3));
  check('F3:alignment_sharpness_formula:'+n,al.fused.eq(expected)&&efficiency.eq(f(3).sub(t2).pow(2).div(f(9).mul(f(1).add(t2).pow(2)))));
  alignment.push({n,t,...al,efficiency});
  const count=2**n,regular=N.cayley(L,t2.neg()),many=record(Array(count).fill(regular),Array(count).fill(t4));
  check('F5:population_count_is_separate_budget:'+n,many.density.eq(f(1).div(f(1).add(t4)))&&many.defect_sum.div(t4).eq(f(count).div(f(1).add(t4))));
  population.push({n,t,count,normalized_density:many.density,total_density:many.defect_sum.div(t4),budget_ratio:many.budget_sum.div(t4)});
}
control('F6:fused_flatness_does_not_exclude_retained_concentration',cancellation.every(r=>r.bound.fused.zero())&&cancellation[0].retained_density.le(cancellation.at(-1).retained_density)&&!cancellation[0].retained_density.eq(cancellation.at(-1).retained_density));
control('F5:growing_record_count_is_not_a_local_singularity',population.every(r=>r.normalized_density.le(1))&&population.at(-1).total_density.le(population[0].total_density)===false);
control('domain:negative_budget_rejected',rejects(()=>record([R],[-1])));
control('domain:zero_budget_rejected',rejects(()=>record([R],[0])));
control('domain:empty_records_rejected',rejects(()=>record([],[])));
control('domain:unmatched_records_rejected',rejects(()=>record([R,J],[1])));
control('domain:nonunit_bridge_rejected',rejects(()=>fusion(source,reference,[N.scale(I,2),I,I])));
control('domain:non_native_central_phase_rejected',rejects(()=>defect(N.scale(I,IOTA))));
control('native:floating_point_scalar_rejected',rejects(()=>f(0.25)));
const witness=certificates.find(c=>c.certificate.steps.length>0),S=systems.get(witness.certificate.presentation),bad=JSON.parse(JSON.stringify(witness.certificate));
bad.steps[0].coefficient=['997','0'];control('certificate:mutated_rewrite_step_rejected',rejects(()=>S.replay(bad)));
const wrong=JSON.parse(JSON.stringify(witness.certificate));wrong.presentation='0'.repeat(64);control('certificate:wrong_presentation_rejected',rejects(()=>S.replay(wrong)));

const P0=N.matrix([[1,0],[0,0]]),K0=N.sub(N.scale(P0,2),I),SX=N.matrix([[0,1],[1,0]]),SY=N.scale(R,IOTA);
function cutTensor(P,X,Y){
  assert.ok(N.equal(mm(P,P),P)&&N.equal(N.dagger(P),P)&&trace(P).eq(1),'Trace-one Hermitian cut required');
  for(const Z of [X,Y])assert.ok(N.equal(N.dagger(Z),Z)&&N.equal(N.add(mm(P,Z),mm(Z,P)),Z),'Hermitian cut tangent required');
  const tangents=[X,Y],Q=tangents.map(x=>tangents.map(y=>trace(mm(P,x,y,P)))),g=Q.map(row=>row.map(x=>x.rad)),omega=Q[0][1].turn;
  const gram=g[0][0].mul(g[1][1]).sub(g[0][1].pow(2)),Ac=tangents.map(x=>N.commutator(P,x)),curvature=N.commutator(X,Y);
  assert.ok(gram.eq(omega.pow(2)),'Two-mode information/curvature area mismatch');
  assert.ok(trace(mm(P,curvature,P)).eq(new Cut(0,omega.mul(2))));
  assert.ok(pair(Ac[0],Ac[0]).eq(g[0][0])&&pair(Ac[0],Ac[1]).eq(g[0][1])&&pair(Ac[1],Ac[1]).eq(g[1][1]),'Minimal-connection Gram must equal cut information metric');
  for(let i=0;i<2;i++)for(let j=0;j<2;j++)assert.ok(trace(mm(P,N.scale(tangents[i],2),N.scale(tangents[j],2))).eq(Q[i][j].mul(4)),'QTH-type score normalization failed');
  return {Q,information_metric:g,omega,gram_determinant:gram,minimal_connection:Ac,curvature,normalized_curvature:gram.zero()?null:pair(curvature,curvature).div(gram.mul(4))};
}
function cayleyDerivative(G,s){s=f(s);const den=f(1).add(s.pow(2)).pow(2);return N.add(N.scale(I,s.mul(-4).div(den)),N.scale(G,f(2).mul(f(1).sub(s.pow(2))).div(den)));}
function movingCut(lambda,x,y){
  [lambda,x,y]=[lambda,x,y].map(f);const a=f(1).add(lambda),s=a.mul(x).div(2),t=y.div(2),U1=N.cayley(R,s),U2=N.cayley(J,t),U=mm(U1,U2);
  const dx=mm(N.scale(cayleyDerivative(R,s),a.div(2)),U2),dy=mm(U1,N.scale(cayleyDerivative(J,t),'1/2'));
  const tangent=dU=>N.add(mm(dU,P0,N.dagger(U)),mm(U,P0,N.dagger(dU)));
  return {P:mm(U,P0,N.dagger(U)),X:tangent(dx),Y:tangent(dy)};
}
const cutCases=[],cutEvolution=[];
for(const [lambda,x,y] of [['0','1/3','1/4'],['1/2','-1/3','1/5'],['1','1/4','-1/6'],['-1/2','2/5','1/3']]){
  const cut=movingCut(lambda,x,y),r=cutTensor(cut.P,cut.X,cut.Y),tag=[lambda,x,y].join(',');
  check('D4:finite_moving_cut_duality:'+tag,r.normalized_curvature!==null&&r.normalized_curvature.eq(1));
  Gs.forEach((G,i)=>{const moved=cutTensor(conj(G,cut.P),conj(G,cut.X),conj(G,cut.Y));check('D2:moving_cut_frame_invariance:'+tag+':'+i,moved.Q.every((row,k)=>row.every((z,j)=>z.eq(r.Q[k][j]))));frameCases++;});
  cutCases.push({parameters:{lambda,x,y},...r});
}
for(const lambda of ['-1/2','0','1/2','1','2']){
  const a=f(1).add(lambda),X=N.scale(SX,a),Y=SY,r=cutTensor(P0,X,Y);
  const derivative=(A,B,dA,dB)=>trace(mm(P0,N.add(mm(dA,B),mm(A,dB)),P0));
  const qxx=derivative(X,X,SX,SX),qxy=derivative(X,Y,SX,O),qyy=derivative(Y,Y,O,O);
  check('D3:joint_information_curvature_evolution:'+lambda,qxx.eq(a.mul(2))&&qxy.eq(IOTA)&&qyy.zero()&&r.Q[0][0].eq(a.pow(2))&&r.Q[0][1].eq(new Cut(0,a)));
  const derivativeCurvature=N.sub(N.commutator(X,Y),N.commutator(Y,X)),frozen=N.commutator(r.minimal_connection[0],r.minimal_connection[1]);
  check('D2:full_minimal_curvature_from_derivatives:'+lambda,N.equal(N.add(derivativeCurvature,frozen),r.curvature)&&N.equal(r.curvature,N.scale(L,a.mul(-2))));
  cutEvolution.push({lambda:f(lambda),a,...r,tensor_derivative:[[qxx,qxy],[qxy.dagger(),qyy]]});
}
control('D2:frozen_connection_bracket_misses_cut_derivatives',N.equal(N.commutator(cutEvolution[1].minimal_connection[0],cutEvolution[1].minimal_connection[1]),N.scale(cutEvolution[1].curvature,-1))&&!N.isZero(cutEvolution[1].curvature));
const cutOnly=N.cayley(K0,'1/3');
control('D1:exponentiating_the_cut_does_not_move_that_cut',N.equal(mm(cutOnly,P0,N.inverse(cutOnly)),P0)&&!N.isZero(N.commutator(R,P0)));
control('D1:cut_diagonal_generator_is_not_recoverable_from_cut_motion',N.equal(N.commutator(N.add(R,N.scale(K0,IOTA)),P0),N.commutator(R,P0))&&!N.equal(N.add(R,N.scale(K0,IOTA)),R));
control('D2:non_tangent_direction_rejected',rejects(()=>cutTensor(P0,K0,SY)));

const output={schema:'native-loop-concentration.v2',status:'PASS_NATIVE_LOOP_CONCENTRATION',source_commit:pins.source_commit,
  scope:'Native common-base/reference-correct fusion, orientation reconstruction, retained concentration and a covariant cut-information-curvature tensor/evolution equation. Native algebra and finite rational cases are executed; general induction, differentiation and refinement arguments have written proofs. No physical focusing law, selected dynamics or gravitational singularity is derived.',
  source_pins_sha256:hash(fs.readFileSync(path.join(root,'SOURCE_PINS.json'))),native_engine:pins.native_engine,upstream,
  presentations,certificates,checks,negative_controls:controls,
  examples:{three_loop_cases:tripleCases,planar,orientation_blindness:{positive,negative},reference_correct_fusion:corrected,finite_bounds:finiteBounds,merging,cancellation_family:cancellation,aligned_family:alignment,population_control:population,cut_information_curvature:cutCases,cut_tensor_evolution:cutEvolution},
  summary:{symbolic_certificates:certificates.length,serialized_replays:certificates.length,native_presentations:presentations.length,exact_checks:checks.length,negative_controls:controls.length,three_loop_cases:tripleCases.length,frame_cases:frameCases,finite_fusion_sizes:finiteBounds.length,cancellation_refinements:cancellation.length,alignment_refinements:alignment.length,population_refinements:population.length,moving_cut_cases:cutCases.length,cut_evolution_cases:cutEvolution.length}};
const serialized=JSON.stringify(output,null,2)+'\n',digest=hash(serialized),file=path.join(root,'CERTIFICATE.json'),checksum=path.join(root,'EXPECTED.sha256');
if(write){fs.writeFileSync(file,serialized);fs.writeFileSync(checksum,digest+'  CERTIFICATE.json\n');}
else{assert.equal(fs.readFileSync(file,'utf8'),serialized,'Frozen certificate differs from native replay');assert.equal(fs.readFileSync(checksum,'utf8'),digest+'  CERTIFICATE.json\n','Certificate checksum mismatch');}
console.log(JSON.stringify({status:output.status,mode:write?'explicit generation':'read-only replay',...output.summary,certificate_sha256:digest},null,2));
