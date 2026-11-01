'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs'),M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const f=N.F.of,mm=(...a)=>a.reduce(N.mul),sc=N.scale,eq=N.equal,zero=N.isZero,cm=N.commutator;
const checks=[],negative_controls=[],proofs=[],presentations=[],systems=new Map(),results={};
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const tr=a=>a.reduce((s,row,i)=>s.add(row[i]),N.ZERO),real=z=>{assert.ok(z.turn.zero());return z.rad;};
const hs=(a,b)=>real(tr(mm(N.dagger(a),b))),proj=v=>mm(v,N.dagger(v)),inner=(a,b)=>mm(N.dagger(a),b)[0][0];
const id=N.identity(2),X=N.matrix([[0,1],[1,0]]),Z=N.matrix([[1,0],[0,-1]]),R=N.matrix([[0,-1],[1,0]]),I=N.IOTA;
const vs=[N.matrix([[1],[0]]),N.matrix([['3/5'],['4/5']]),N.matrix([['3/5'],[N.Cut.of([0,'4/5'])]])],ps=vs.map(proj),mirror=proj(N.matrix([['3/5'],[N.Cut.of([0,'-4/5'])]])),PM=[ps[0],ps[1],mirror];
const overlap=(p,q)=>real(tr(mm(p,q))),dist2=(p,q)=>hs(N.sub(p,q),N.sub(p,q)).div(2),cycle=arr=>tr(mm(...arr));
for(let i=0;i<3;i++){
 ck('cuts:rank_one:'+i,eq(mm(ps[i],ps[i]),ps[i])&&eq(N.dagger(ps[i]),ps[i])&&tr(ps[i]).eq(1));
 for(let j=0;j<3;j++){
  ck('geometry:distance_overlap:'+i+':'+j,dist2(ps[i],ps[j]).eq(f(1).sub(overlap(ps[i],ps[j]))));
  ck('geometry:mirror_pair_agreement:'+i+':'+j,overlap(ps[i],ps[j]).eq(overlap(PM[i],PM[j])));
 }
}
ck('geometry:overlap01',overlap(ps[0],ps[1]).eq('9/25'));
ck('geometry:overlap02',overlap(ps[0],ps[2]).eq('9/25'));
ck('geometry:overlap12',overlap(ps[1],ps[2]).eq('337/625'));
ck('geometry:positive_cycle',cycle(ps).eq(N.Cut.of(['81/625','144/625'])));
ck('geometry:mirror_cycle',cycle(PM).eq(N.Cut.of(['81/625','-144/625'])));
no('geometry:distances_classify_unitary_cut_orbit',!cycle(ps).eq(cycle(PM)));
const W=N.matrix([['3/5','-4/5'],['4/5','3/5']]),phases=[N.Cut.of(['5/13','12/13']),I,N.Cut.of(['3/5','4/5'])],target=vs.map((v,i)=>sc(mm(W,v),phases[i]));
const gauged=target.map((v,j)=>j===0?v:sc(v,inner(target[0],v).dagger().div('3/5')));
const gram=arr=>N.matrix(arr.map(a=>arr.map(b=>inner(a,b))));
ck('reconstruct:tree_gauge_gram',eq(gram(gauged),gram(vs)));
ck('reconstruct:non_tree_edge_from_cycle',cycle(ps).div(f('3/5').pow(2)).eq(inner(vs[1],vs[2])));
ck('reconstruct:span_rank',N.rank(gram(vs))===2);
const columns=(a,b)=>N.matrix(a.map((row,i)=>[row[0],b[i][0]])),U=mm(columns(gauged[0],gauged[1]),N.inverse(columns(vs[0],vs[1])));
ck('reconstruct:unitary_lift',eq(mm(N.dagger(U),U),id));
for(let j=0;j<3;j++)ck('reconstruct:all_cuts:'+j,eq(mm(U,ps[j],N.dagger(U)),proj(target[j])));
no('reconstruct:mirror_non_tree_edge_agrees',!cycle(PM).div(f('3/5').pow(2)).eq(inner(vs[1],vs[2])));
results.reconstruction={overlaps:['9/25','9/25','337/625'],oriented_cycle:cycle(ps),mirror_cycle:cycle(PM),recovered_unitary:U,span_dimension:2};
const frameGenerators=[sc(X,I),R,sc(Z,I)],localFrames=frameGenerators.map(k=>N.cayley(k,'1/5'));
for(let j=0;j<3;j++){
 const gj=localFrames[j],local=mm(gj,ps[j],N.dagger(gj)),link=mm(W,N.dagger(gj)),based=mm(link,local,N.dagger(link));
 ck('frame:independent_links:'+j,eq(based,mm(W,ps[j],N.dagger(W))));
}
no('frame:ignore_comparison_link',!eq(mm(localFrames[0],ps[0],N.dagger(localFrames[0])),mm(W,ps[0],N.dagger(W))));
const pathCut=mm(Z,ps[1],Z);
ck('path:alternate_distance',dist2(ps[1],pathCut).eq('576/625'));
no('path:arbitrary_loop_is_path_independent',!zero(cm(Z,ps[1]))&&!eq(pathCut,ps[1]));
ck('path:commuting_loop',eq(mm(Z,ps[0],Z),ps[0]));
results.paths={same_path_distance_squared:'0',changed_path_distance_squared:'576/625'};

// Full two-mode anti-Hermitian algebra modulo its scalar iI gauge.
const basis=frameGenerators;
function fit(cuts,tangents){
 const images=basis.map(k=>cuts.map(p=>cm(k,p))),normal=N.matrix(basis.map((_,a)=>basis.map((_,b)=>images[a].reduce((s,x,j)=>s.add(hs(x,images[b][j])),f(0))))),rhs=N.matrix(basis.map((_,a)=>[images[a].reduce((s,x,j)=>s.add(hs(x,tangents[j])),f(0))]));
 const coefficients=mm(N.inverse(normal),rhs),K=basis.reduce((s,k,a)=>N.add(s,sc(k,coefficients[a][0])),N.zeros(2)),residuals=cuts.map((p,j)=>N.sub(tangents[j],cm(K,p))),epsilon2=residuals.reduce((s,x)=>s.add(hs(x,x)),f(0));
 return {K,coefficients,normal,rhs,residuals,epsilon2,images};
}
const known=N.add(N.add(sc(basis[0],'2/3'),sc(basis[1],'3/5')),sc(basis[2],'-4/7')),tangents=ps.map(p=>cm(known,p)),rigid=fit(ps,tangents);
ck('fit:common_generator_recovered',eq(rigid.K,known));
ck('fit:common_residual_zero',rigid.epsilon2.zero());
ck('fit:normal_equation',eq(mm(rigid.normal,rigid.coefficients),rigid.rhs));
ck('fit:three_visible_generator_directions',N.rank(rigid.normal)===3);
ck('fit:scalar_gauge_invisible',ps.every(p=>zero(cm(sc(id,I),p))));
ck('fit:source_reading_changes',!real(tr(mm(ps[1],tangents[0]))).zero());
for(let i=0;i<3;i++)for(let j=i+1;j<3;j++)ck('fit:rigid_overlap_derivative:'+i+':'+j,real(tr(N.add(mm(tangents[i],ps[j]),mm(ps[i],tangents[j])))).zero());
const derivative=(cuts,tans)=>cuts.reduce((s,_,a)=>s.add(tr(mm(...cuts.map((p,j)=>a===j?tans[j]:p)))),N.ZERO);
ck('fit:rigid_cycle_derivative',derivative(ps,tangents).zero());
const rvs=[vs[0],vs[1],N.matrix([['5/13'],['12/13']])],rps=rvs.map(proj),x2=N.matrix([[0,N.Cut.of([0,'-60/169'])],[N.Cut.of([0,'60/169']),0]]),tx=[N.zeros(2),N.zeros(2),x2],strained=fit(rps,tx);
for(let j=0;j<3;j++)ck('strain:valid_cut_tangent:'+j,eq(N.add(mm(tx[j],rps[j]),mm(rps[j],tx[j])),tx[j])&&eq(N.dagger(tx[j]),tx[j]));
for(let i=0;i<3;i++)for(let j=i+1;j<3;j++){
 const rate=real(tr(N.add(mm(tx[i],rps[j]),mm(rps[i],tx[j])))),rr=real(tr(N.add(mm(strained.residuals[i],rps[j]),mm(rps[i],strained.residuals[j]))));
 ck('strain:pair_blind:'+i+':'+j,rate.zero());
 ck('strain:residual_drives_overlap:'+i+':'+j,rate.eq(rr));
 ck('strain:budget:'+i+':'+j,rate.pow(2).le(strained.epsilon2.mul(2)));
}
ck('strain:cycle_changes',derivative(rps,tx).eq(N.Cut.of([0,'720/4225'])));
ck('strain:residual_drives_cycle',derivative(rps,tx).eq(derivative(rps,strained.residuals)));
ck('strain:normal_equation',eq(mm(strained.normal,strained.coefficients),strained.rhs));
for(let a=0;a<3;a++)ck('strain:orthogonal_residual:'+a,strained.images[a].reduce((s,x,j)=>s.add(hs(x,strained.residuals[j])),f(0)).zero());
no('strain:pair_derivatives_zero_implies_common_transport',!strained.epsilon2.zero());
const firstTwo=fit(rps.slice(0,2),tx.slice(0,2));
ck('strain:fixed_two_force_scalar_gauge',N.rank(firstTwo.normal)===3&&zero(firstTwo.K));
no('strain:individual_transport_means_joint_transport',tx.every((x,j)=>eq(cm(cm(x,rps[j]),rps[j]),x))&&!strained.epsilon2.zero());
const movingFrameRate=sc(basis[0],'1/2'),gps=rps.map(p=>mm(W,p,N.dagger(W))),gtx=tx.map((x,j)=>N.add(mm(W,x,N.dagger(W)),cm(movingFrameRate,gps[j]))),gf=fit(gps,gtx);
ck('strain:time_dependent_frame_generator',eq(gf.K,N.add(mm(W,strained.K,N.dagger(W)),movingFrameRate)));
ck('strain:frame_invariant_residual',gf.epsilon2.eq(strained.epsilon2));
for(let j=0;j<3;j++)ck('strain:frame_covariant_components:'+j,eq(gf.residuals[j],mm(W,strained.residuals[j],N.dagger(W))));
const cldDefects=tx.map((x,j)=>N.add(x,sc(cm(Z,rps[j]),I))),relativeFit=fit(rps,cldDefects);
ck('CLD:known_source_shift',eq(relativeFit.K,N.add(strained.K,sc(Z,I))));
ck('CLD:source_shift_same_residual',relativeFit.epsilon2.eq(strained.epsilon2));
const changing=ps.map((p,j)=>j===1?cm(R,p):N.zeros(2)),changedFit=fit(ps,changing);
no('strain:single_cut_motion_is_rigid_family',!changedFit.epsilon2.zero());
for(let i=0;i<3;i++)for(let j=i+1;j<3;j++){
 const rate=real(tr(N.add(mm(changing[i],ps[j]),mm(ps[i],changing[j])))),rr=real(tr(N.add(mm(changedFit.residuals[i],ps[j]),mm(ps[i],changedFit.residuals[j]))));
 ck('changing:residual_drives_distance:'+i+':'+j,rate.eq(rr));
 ck('changing:nontrivial_budget:'+i+':'+j,rate.pow(2).le(changedFit.epsilon2.mul(2)));
 if(i===0&&j===1)ck('changing:nonzero_distance_rate',!rate.zero());
}
ck('changing:cycle_derivative',derivative(ps,changing).eq(derivative(ps,changedFit.residuals)));
results.common_transport={rigid_generator:rigid.K,rigid_residual_squared:rigid.epsilon2,distance_blind_generator_fit:strained.K,distance_blind_residual_squared:strained.epsilon2,distance_blind_cycle_derivative:derivative(rps,tx),normal_matrix:strained.normal};

const e=a=>N.add(mm(ps[0],a,ps[0]),mm(N.sub(id,ps[0]),a,N.sub(id,ps[0]))),delta=N.sub(ps[1],ps[2]),diag=e(delta),off=N.sub(delta,diag);
ck('measurement:pinching_idempotent',eq(e(e(ps[1])),e(ps[1])));
ck('measurement:orthogonal_split',hs(diag,off).zero());
ck('measurement:pythagoras',hs(delta,delta).eq(hs(diag,diag).add(hs(off,off))));
ck('measurement:initial_separation',dist2(ps[1],ps[2]).eq('288/625'));
ck('measurement:erased_separation',dist2(e(ps[1]),e(ps[2])).zero());
ck('measurement:identical_outputs',eq(e(ps[1]),N.matrix([['9/25',0],[0,'16/25']]))&&eq(e(ps[1]),e(ps[2])));
no('measurement:outputs_are_still_pure_cuts',!eq(mm(e(ps[1]),e(ps[1])),e(ps[1])));
no('measurement:cut_dependence_implies_geometry_invariance',!dist2(ps[1],ps[2]).eq(dist2(e(ps[1]),e(ps[2]))));
results.measurement={before_distance_squared:'288/625',after_distance_squared:'0',common_output:e(ps[1])};

function register(name,spec){const s=new M.Presentation(spec),audit=s.audit();ck(name+':confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name,spec:s.spec(),audit});systems.set(s.hash(),s);return{s,w:Object.fromEntries(spec.tokens.map(t=>{const n=typeof t==='string'?t:t.name;return[n,s.word([n])];}))};}
function prove(name,lhs,rhs){const s=lhs.system,r=s.reduce(lhs.minus(rhs),{witness:true});ck(name,r.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(r.certificate));ck(name+':replay',s.replay(certificate));proofs.push({name,lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate});}
{
 const{s,w:{p,q,k,a}}=register('cut_relations',{tokens:['p','q','k','a'],rules:[{id:'pp',lhs:['p','p'],rhs:[[['p'],1]],source:'First cut idempotence'},{id:'qq',lhs:['q','q'],rhs:[[['q'],1]],source:'Second cut idempotence'}]}),one=s.one(),c=one.minus(p),comm=(x,y)=>x.times(y).minus(y.times(x)),pinch=x=>p.times(x).times(p).plus(c.times(x).times(c));
 prove('words:distance_square',p.minus(q).times(p.minus(q)),p.plus(q).minus(p.times(q)).minus(q.times(p)));
 prove('words:double_commutator',comm(comm(k,p),p),k.times(p).plus(p.times(k)).minus(p.times(k).times(p).scale(2)));
 prove('words:pinching_idempotence',pinch(pinch(a)),pinch(a));
 prove('words:pinching_offdiagonal_zero',pinch(a.minus(pinch(a))),s.zero());
 prove('words:complete_block_record',pinch(a).plus(p.times(a).times(c)).plus(c.times(a).times(p)),a);
 prove('words:common_generator_pair',comm(k,p).times(q).plus(p.times(comm(k,q))),comm(k,p.times(q)));
 prove('words:common_generator_cycle',comm(k,p).times(q).times(a).plus(p.times(comm(k,q)).times(a)).plus(p.times(q).times(comm(k,a))),comm(k,p.times(q).times(a)));
}
{
 const{s,w:{p,q,u,v}}=register('unitary_comparison',{tokens:['p','q','u','v'],rules:[{id:'uv',lhs:['u','v'],rhs:[[[],1]],source:'Unitary and inverse'},{id:'vu',lhs:['v','u'],rhs:[[[],1]],source:'Inverse and unitary'}]});
 prove('words:transport_product',u.times(p).times(v).times(u).times(q).times(v),u.times(p).times(q).times(v));
 prove('words:transport_commutator',u.times(p).times(v).times(u).times(q).times(v).minus(u.times(q).times(v).times(u).times(p).times(v)),u.times(p.times(q).minus(q.times(p))).times(v));
}
const chosen=proofs.find(x=>x.certificate.steps.length),bad=JSON.parse(JSON.stringify(chosen.certificate));bad.steps[0].coefficient=['991','0'];
no('integrity:corrupted_word_proof',(()=>{try{systems.get(chosen.certificate.presentation).replay(bad);return false;}catch{return true;}})());
const hash=x=>crypto.createHash('sha256').update(x).digest('hex'),payload={protocol:'RCG1_RELATIONAL_CUT_SPACE_AND_COMMON_TRANSPORT_V1',status:'PASS_NATIVE_RELATIONAL_CUT_GEOMETRY_AND_TRANSPORT_RESIDUAL',scope:'Finite exact native matrices and word witnesses with separate written reconstruction and transport proofs. Candidate operational geometry, not a physical spacetime or gravity derivation.',source_pins_sha256:hash(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),presentations,proofs,checks,negative_controls,results,summary:{exact_checks:checks.length,negative_controls:negative_controls.length,native_presentations:presentations.length,symbolic_certificates:proofs.length,serialized_replays:proofs.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',sha=hash(bytes),file=path.join(__dirname,'CERTIFICATE.json');
if(process.argv.includes('--write')){fs.writeFileSync(file,bytes);fs.writeFileSync(path.join(__dirname,'EXPECTED.sha256'),sha+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(file,'utf8'),bytes);assert.equal(fs.readFileSync(path.join(__dirname,'EXPECTED.sha256'),'utf8'),sha+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:sha,distance_blind_residual_squared:strained.epsilon2}));
