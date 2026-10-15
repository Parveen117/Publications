'use strict';
const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const m=N.matrix,mm=(...xs)=>xs.reduce(N.mul),dg=N.dagger,eq=N.equal,sc=N.scale,f=N.F.of;
const checks=[],negative_controls=[],proofs=[],presentations=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name:'binary:'+name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name:'binary:'+name,rejected:true});};
const rejects=fn=>{try{fn();return false;}catch{return true;}};
const diag=values=>{const D=N.zeros(values.length);values.forEach((v,j)=>D[j][j]=N.Cut.of(v));return D;};
const form=(H,x)=>mm(dg(x),H,x)[0][0];
// Constructive exact completion of the square, without spectral fitting or roots.
function factor(H){
 H=m(H);assert.ok(eq(H,dg(H)),'Hermitian form required');
 let residual=H;const rows=[],weights=[],pivots=[];
 while(!N.isZero(residual)){
  let p=-1;
  for(let j=0;j<H.length;j++){
   const d=residual[j][j];assert.ok(d.turn.zero()&&f(0).le(d.rad),'Negative/nonreal pivot');
   if(!d.zero()&&p===-1)p=j;
  }
  assert.ok(p>=0,'Nonzero residual with all zero diagonals is not positive');
  const a=residual[p][p],row=[residual[p].map(v=>v.div(a))];
  residual=N.sub(residual,sc(mm(dg(row),row),a));
  rows.push(row[0]);weights.push(a);pivots.push(p);
  assert.ok(rows.length<=H.length,'Rank termination');
 }
 return {L:rows.length?m(rows):null,D:rows.length?diag(weights):null,weights,pivots,rank:rows.length};
}
function tree(values){
 if(values.length===0)return {sum:f(0),leaves:0,merges:0,depth:0};
 if(values.length===1)return {sum:values[0],leaves:1,merges:0,depth:0};
 const mid=(values.length+1)>>1,left=tree(values.slice(0,mid)),right=tree(values.slice(mid));
 return {sum:left.sum.add(right.sum),leaves:left.leaves+right.leaves,merges:left.merges+right.merges+1,depth:1+Math.max(left.depth,right.depth)};
}
const H=m([[2,1],[1,5]]),I=N.identity(2),P=m([[1,0],[0,0]]),Q=N.sub(I,P);
const source=JSON.parse(fs.readFileSync(path.join(__dirname,'../yang-mills-certified-benchmark/certificates/YM54_RESULT.json'),'utf8')).thermo_fixture;
ck('published_thermo_source',eq(H,m(source.H)));
const oldB=m([[2,1],[0,3]]);
ck('HCM_factor_reused',eq(mm(dg(oldB),oldB),sc(H,2)));
const v=m([['3/5'],[new N.Cut(0,'4/5')]]),complexH=N.add(sc(P,3),mm(v,dg(v)));
const rows3=m([[1,'1/2',0,0],[0,1,'1/3',0],[0,0,1,0]]),rank3=mm(dg(rows3),diag([2,3,5]),rows3);
const fixtures=[['zero',N.zeros(2)],['rank_one',P],['thermo',H],['complex',complexH],['rank_three',rank3],['rank_five',N.identity(5)]];
const factors=[];
for(const [name,metric]of fixtures){
 const result=factor(metric),{L,D,rank}=result;
 ck(name+':rank_minimal',rank===N.rank(metric));
 const recovered=rank?mm(dg(L),D,L):N.zeros(metric.length);
 ck(name+':exact_reconstruction',eq(recovered,metric));
 ck(name+':positive_weights',result.weights.every(a=>a.turn.zero()&&f(0).le(a.rad)&&!a.zero()));
 if(rank){
  ck(name+':independent_channels',N.rank(L)===rank);
  const S=mm(dg(L),N.inverse(mm(L,dg(L))));
  ck(name+':quotient_section',eq(mm(L,S),N.identity(rank)));
  ck(name+':kernel_is_metric_kernel',eq(mm(metric,N.sub(N.identity(metric.length),mm(S,L))),N.zeros(metric.length)));
 }
 const x=m(Array.from({length:metric.length},(_,j)=>[new N.Cut(j+1,j%2)]));
 const y=rank?mm(L,x):null;
 const values=rank?result.weights.map((a,j)=>a.rad.mul(y[j][0].norm2())):[];
 const t=tree(values),expected=form(metric,x);
 ck(name+':binary_total',expected.turn.zero()&&t.sum.eq(expected.rad));
 ck(name+':binary_minimum_merges',t.merges===Math.max(0,rank-1));
 let capacity=1,depth=0;while(capacity<rank){capacity*=2;depth++;}
 ck(name+':binary_minimum_depth',t.depth===depth);
 factors.push({name,...result,tree:t});
}
const hf=factor(H);
ck('thermo:actual_two_channels',eq(hf.L,m([[1,'1/2'],[0,1]]))&&eq(hf.D,diag([2,'9/2'])));
no('indefinite_form_is_positive_sum',rejects(()=>factor(m([[1,0],[0,-1]]))));
no('zero_diagonal_hides_nonpositive_cross_terms',rejects(()=>factor(m([[0,1],[1,0]]))));
no('nonhermitian_form_admitted',rejects(()=>factor(m([[1,2],[0,1]]))));
// Exact invariant-preserving transport; the operation is declared, not physical time.
const J=m([[0,-1],[1,0]]),G=mm(N.inverse(H),J),T=N.cayley(G,'1/3');
ck('transport:H_skew_generator',N.isZero(N.add(mm(dg(G),H),mm(H,G))));
ck('transport:H_invariant',eq(mm(dg(T),H,T),H));
const S=mm(dg(hf.L),N.inverse(mm(hf.L,dg(hf.L)))),U=mm(hf.L,T,S);
ck('transport:exact_descent',eq(mm(hf.L,T),mm(U,hf.L)));
ck('transport:weighted_isometry',eq(mm(dg(U),hf.D,U),hf.D));
ck('transport:inverse_recovery',eq(mm(N.inverse(U),U),I));
ck('transport:composition',eq(mm(hf.L,T,T,S),mm(U,U)));
no('transport:Euclidean_metric_silently_substituted',!eq(mm(dg(T),T),I));
const H3=m([[2,1,0],[1,5,0],[0,0,0]]),h3=factor(H3);
const T3=m([[T[0][0],T[0][1],0],[T[1][0],T[1][1],0],[2,3,0]]);
ck('quotient:source_isometry_with_lost_kernel',eq(mm(dg(T3),H3,T3),H3));
const U3=N.descendedAction(T3,h3.L);
ck('quotient:descended_weighted_isometry',U3!==null&&eq(mm(dg(U3),h3.D,U3),h3.D));
no('quotient:invertible_full_source_from_invertible_record',N.rank(T3)===2&&N.rank(U3)===2);
const e3=m([[0],[0],[1]]);
no('quotient:kernel_is_recovered',N.isZero(mm(h3.L,e3))&&!N.isZero(e3));
// The old operation-closure engine supplies the smallest future-complete record.
const u=m([['3/5'],['4/5']]),B=mm(u,dg(u)),L0=m([[1,0]]),e2=m([[0],[1]]);
const closure=N.rowClosure(L0,[B]);
ck('future:rank_one_grows_to_two',JSON.stringify(closure.ranks)==='[1,2]');
no('future:initial_single_channel_is_autonomous',N.descendedAction(B,L0)===null);
const descended=N.descendedAction(B,closure.observer);
ck('future:minimal_repaired_intertwiner',descended!==null&&eq(mm(closure.observer,B),mm(descended,closure.observer)));
ck('future:initial_hidden_reading_zero',form(P,e2).zero());
ck('future:after_cut_hidden_reading',form(P,mm(B,e2)).eq('144/625'));
// Natural fixed branch metric versus the fully transported observable.
const W=Q.concat(P),branchMetric=N.kron(I,H),encodedH=mm(W,H,dg(W));
const fixedPullback=mm(dg(W),branchMetric,W),cross=N.add(mm(P,H,Q),mm(Q,H,P));
ck('metric:fixed_branch_pullback',eq(fixedPullback,N.add(mm(P,H,P),mm(Q,H,Q))));
ck('metric:double_commutator',eq(N.sub(H,fixedPullback),N.commutator(P,N.commutator(P,H))));
ck('metric:missing_cross_record',eq(N.sub(H,fixedPullback),cross));
ck('metric:complete_observable_recovery',eq(mm(dg(W),encodedH,W),H));
no('metric:norm_isometry_preserves_every_fixed_branch_metric',!eq(fixedPullback,H));
ck('metric:commuting_metric_is_conserved',eq(mm(dg(W),N.kron(I,diag([2,5])),W),diag([2,5])));
const plus=m([['3/5'],['4/5']]),minus=m([['3/5'],['-4/5']]);
for(const [label,x,value,loss]of [['plus',plus,'122/25','24/25'],['minus',minus,'74/25','-24/25']]){
 ck('metric:'+label+':source',form(H,x).eq(value));
 ck('metric:'+label+':classical_branches',form(fixedPullback,x).eq('98/25'));
 ck('metric:'+label+':cross_information',form(cross,x).eq(loss));
 ck('metric:'+label+':coherent_value_retained',form(encodedH,mm(W,x)).eq(value));
}
no('amount:total_norm_encodes_next_cut',form(I,plus).eq(form(I,minus))&&!form(B,plus).eq(form(B,minus)));
ck('amount:distinguishing_cut_values',form(B,plus).eq(1)&&form(B,minus).eq('49/625'));
// Root descriptors stay exact; no approximate root scalar is introduced.
const square2=m([[0,2],[1,0]]),cube2=m([[0,0,2],[1,0,0],[0,1,0]]);
ck('roots:square_companion',eq(N.power(square2,2),sc(I,2)));
ck('roots:cube_companion',eq(N.power(cube2,3),sc(N.identity(3),2)));
ck('roots:cubic_minimal_polynomial_degree',N.rank([N.flatten(N.identity(3)),N.flatten(cube2),N.flatten(N.power(cube2,2))])===3);
for(const candidate of [-2,-1,1,2])ck('roots:cubic_rational_candidate_rejected:'+candidate,!f(candidate).pow(3).sub(2).zero());
ck('roots:sqrt2_isolating_interval',f(1).pow(2).le(2)&&!f(1).pow(2).eq(2)&&f(2).le(f('3/2').pow(2))&&!f('3/2').pow(2).eq(2));
ck('roots:cuberoot2_isolating_interval',f(1).pow(3).le(2)&&!f(1).pow(3).eq(2)&&f(2).le(f(2).pow(3))&&!f(2).pow(3).eq(2));
no('roots:rational_approximation_is_exact_sqrt2',!f('7/5').pow(2).eq(2));
// One commutative square-completion identity and two noncommutative metric identities.
function prove(system,name,lhs,rhs){
 const r=system.reduce(lhs.minus(rhs),{witness:true});ck(name,r.normal.terms.size===0);
 const certificate=JSON.parse(JSON.stringify(r.certificate));ck(name+':replay',system.replay(certificate));
 proofs.push({name:'binary:'+name,lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate});
}
function register(name,s){const audit=s.audit();ck(name+':confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name,spec:s.spec(),audit});return s;}
const scalar=register('real_source_coordinates',new M.Presentation({tokens:['a','b'],rules:[{id:'ba',lhs:['b','a'],rhs:[[['a','b'],1]],source:'Commuting real scalar coordinates of the inherited source'}]}));
const a=scalar.word(['a']),b=scalar.word(['b']),ell=a.plus(b.scale('1/2'));
prove(scalar,'weighted_Pythagorean_source',a.times(a).scale(2).plus(a.times(b).scale(2)).plus(b.times(b).scale(5)),ell.times(ell).scale(2).plus(b.times(b).scale('9/2')));
const ops=register('cut_and_fixed_metric',new M.Presentation({tokens:['p','h'],rules:[{id:'pp',lhs:['p','p'],rhs:[[['p'],1]],source:'Orthogonal cut idempotence; h is unrestricted'}]}));
const p=ops.word(['p']),h=ops.word(['h']),q=ops.one().minus(p),com=p.times(h).minus(h.times(p));
const missing=h.minus(p.times(h).times(p).plus(q.times(h).times(q)));
prove(ops,'fixed_metric_double_commutator',missing,p.times(com).minus(com.times(p)));
prove(ops,'fixed_metric_cross_record',missing,p.times(h).times(q).plus(q.times(h).times(p)));
no('words:cut_commutes_with_every_metric',ops.reduce(com).normal.terms.size!==0);
module.exports={checks,negative_controls,proofs,presentations,results:{factors,source_metric:H,transport:T,descended_transport:U,future_closure:closure,fixed_branch_pullback:fixedPullback,full_coherent_observable:encodedH,root_companions:{square2,cube2},scope:{class:'finite positive quadratic forms over Q+iota Q',minimum:'linear scalar channels with retained sign/phase',arbitrary_roots_by_square_towers:false,root_descriptors_are_not_finite_positional_expansions:true,full_UGD_homomorphism:false,physical_gravity_completed:false}}};
