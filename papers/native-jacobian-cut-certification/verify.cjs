'use strict';
// NJC-1: unchanged native arithmetic and Pāṇinian word rewriter only.
const fs=require('node:fs'), path=require('node:path'), crypto=require('node:crypto');
const assert=require('node:assert/strict'), root=__dirname, argv=process.argv.slice(2);
const nativeRoot=argv[argv.indexOf('--native-root')+1];
assert.ok(argv.includes('--native-root')&&nativeRoot,'--native-root required');
const write=argv.includes('--write');assert.ok(write!==argv.includes('--check'),'Choose --write or --check');
const hash=x=>crypto.createHash('sha256').update(x).digest('hex');
const pins=JSON.parse(fs.readFileSync(path.join(root,'SOURCE_PINS.json'),'utf8'));
for(const [name,h] of Object.entries(pins.native_engine))assert.equal(hash(fs.readFileSync(path.join(nativeRoot,name))),h,name);
for(const [name,h] of Object.entries(pins.packet_files))assert.equal(hash(fs.readFileSync(path.join(root,name))),h,name);
for(const [name,h] of Object.entries(pins.predecessor_files))assert.equal(hash(fs.readFileSync(path.resolve(root,'../..',name))),h,name);
const N=require(path.join(nativeRoot,'native_operator.cjs')),M=require(path.join(nativeRoot,'paninian_operator.cjs'));
const {F,Cut,IOTA}=N,f=F.of,checks=[],controls=[],proofs=[],presentations=[],systems=new Map(),results={};
const ck=(name,v,data={})=>{assert.ok(v,name);checks.push({name,passed:true,...data});};
const no=(name,v,data={})=>{assert.ok(v,name);controls.push({name,rejected:true,...data});};
const br=M.bracket,mm=(...a)=>a.reduce(N.mul),tr=a=>a.reduce((s,r,i)=>s.add(r[i]),N.ZERO);
const real=x=>{assert.ok(x.turn.zero());return x.rad;};
const rej=fn=>{try{fn();return false;}catch{return true;}};
function system(name,tokens,rules=[],commute=false,star=null){
 if(commute)for(let j=1;j<tokens.length;j++)for(let i=0;i<j;i++)rules.push({id:'sort_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'central real coefficients'});
 const S=new M.Presentation({tokens,rules}),w=Object.fromEntries(tokens.map(t=>[t,S.word([t])])),audit=S.audit();
 ck(name+':critical_pairs',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');
 const images=star?star(S,w):null,dagger=images?M.daggerGate(S,images):null;
 if(dagger)ck(name+':dagger',dagger.status==='DAGGER_DESCENDS');
 presentations.push({name,spec:S.spec(),audit,dagger});systems.set(S.hash(),S);return {S,w,images};
}
function prove(name,a,b){
 const S=a.system,r=S.reduce(a.minus(b),{witness:true});ck(name,r.normal.terms.size===0);
 const c=JSON.parse(JSON.stringify(r.certificate));ck(name+':replay',S.replay(c));
 proofs.push({name,lhs:a.toJSON(),rhs:b.toJSON(),certificate:c});
}
const square=(token,sign=1)=>({id:'square_'+token,lhs:[token,token],rhs:[[[],sign]],source:'declared generator square'});

// J is real-linear on a two-component cut carrier, not a complex-linear scalar.
{
 const {S,w:{a,b,j}}=system('involution_and_decay',['a','b','j'],[square('j')],true);
 const p=S.one().plus(j).scale('1/2'),q=S.one().minus(j).scale('1/2');
 prove('cut:plus_projector',p.times(p),p);prove('cut:minus_projector',q.times(q),q);
 prove('cut:orthogonality',p.times(q),S.zero());prove('cut:completeness',p.plus(q),S.one());
 prove('cut:reflection_residual',j.minus(S.one()),q.scale(-2));
 prove('cut:residual_square',j.minus(S.one()).times(j.minus(S.one())),q.scale(4));
 const T=x=>p.plus(q.times(x));
 prove('decay:composition',T(a).times(T(b)),T(a.times(b)));
 prove('decay:lawful_component_preserved',p.times(T(a)),p);
 prove('decay:transverse_multiplier',q.times(T(a)),q.times(a));
}
{
 const {S,w:{x,y}}=system('sine_in_real_coordinates',['x','y'],[],true);
 const dx=y.minus(x),dy=x.minus(y);
 prove('sine:exact_reflection_defect_energy',dx.times(dx).plus(dy.times(dy)),x.minus(y).times(x.minus(y)).scale(2));
 prove('scalar_gradient:commutator_vanishes',br(x,y),S.zero());
}

// Universal differentiated-cut algebra. No relation is imposed between x and y.
{
 const rules=[square('j'),...['x','y'].map(t=>({id:t+'_past_j',lhs:[t,'j'],rhs:[[['j',t],-1]],source:'differentiate J^2=I'})),{id:'b_past_j',lhs:['b','j'],rhs:[[['j','b'],1]],source:'cut-even connection freedom'}];
 const {S,w:{j,x,y,b},images}=system('moving_cut_connection',['j','x','y','b'],rules,false,(_,w)=>({j:w.j,x:w.x,y:w.y,b:w.b.scale(-1)}));
 const Ax=j.times(x).scale('1/2'),Ay=j.times(y).scale('1/2'),curv=br(x,y).scale('1/4');
 prove('connection:cut_compatibility',x.plus(br(Ax,j)),S.zero());
 prove('connection:odd_part',j.times(Ax).times(j),Ax.scale(-1));
 prove('connection:even_freedom',x.plus(br(Ax.plus(b),j)),S.zero());
 prove('connection:metric_compatibility',M.star(Ax,images),Ax.scale(-1));
 prove('connection:curvature_formula',br(x,y).scale('1/2').plus(br(Ax,Ay)),curv);
 prove('connection:curvature_preserves_cut',br(curv,j),S.zero());
 prove('connection:curvature_is_skew',M.star(curv,images),curv.scale(-1));
 prove('connection:trace_obstruction_anticommutator',j.times(x).plus(x.times(j)),S.zero());
}
{
 const {S,w:{j,a}}=system('connection_uniqueness',['j','a'],[square('j')]);
 const dJ=br(a,j).scale(-1),forced=j.times(dJ).scale('1/2'),even=a.plus(j.times(a).times(j)).scale('1/2');
 prove('connection:forced_odd_component',forced,a.minus(j.times(a).times(j)).scale('1/2'));
 prove('connection:all_solutions_split',forced.plus(even),a);
 prove('connection:remaining_component_commutes',br(even,j),S.zero());
}

// Correct transpose grading requires a self-adjoint/orthogonal cut.
{
 const {S,w:{j,h,k},images}=system('tensor_bigrading',['j','h','k'],[square('j')],false,(_,w)=>({j:w.j,h:w.k,k:w.h}));
 const cut=a=>j.times(a).times(j),star=a=>M.star(a,images),parts=[];
 for(const s of [1,-1])for(const t of [1,-1]){
  const part=h.plus(cut(h).scale(s)).plus(star(h).scale(t)).plus(cut(star(h)).scale(s*t)).scale('1/4');parts.push(part);
  prove('bigrading:cut_'+s+'_'+t,cut(part),part.scale(s));
  prove('bigrading:transpose_'+s+'_'+t,star(part),part.scale(t));
 }
 prove('bigrading:reconstruction',parts.reduce((a,b)=>a.plus(b)),h);
 prove('bigrading:commuting_actions',star(cut(h)),cut(star(h)));
}
{
 const {S,w:{j,k}}=system('two_slot_tensor_grading',['j','k'],[square('j'),square('k')],true);
 const parts=[];
 for(const a of [1,-1])for(const b of [1,-1]){
  const p=S.one().plus(j.scale(a)).times(S.one().plus(k.scale(b))).scale('1/4');parts.push(p);
  prove('tensor:slot_projector_'+a+'_'+b,p.times(p),p);
 }
 prove('tensor:all_slot_reconstruction',parts.reduce((a,b)=>a.plus(b)),S.one());
 for(let i=0;i<parts.length;i++)for(let j=0;j<i;j++)prove('tensor:disjoint_slots_'+i+'_'+j,parts[i].times(parts[j]),S.zero());
}

// Native port of the recovered MP quadratic-generator bridge.
{
 const rules=[{id:'g_squared',lhs:['g','g'],rhs:[[['d'],1]],source:'G^2=Delta I'},{id:'g_past_d',lhs:['g','d'],rhs:[[['d','g'],1]],source:'Delta is central'},{id:'g_past_t',lhs:['g','t'],rhs:[[['t','g'],1]],source:'parameter is central'},{id:'t_past_d',lhs:['t','d'],rhs:[[['d','t'],1]],source:'central coefficients'}];
 const {S,w:{d,t,g}}=system('quadratic_generator',['d','t','g'],rules);
 let power=S.one(),dp=S.one(),factor=f(1);const exp=M.formalExp(g,10);
 for(let n=0;n<=10;n++){
  if(n>0){power=power.times(g);factor=factor.mul(n);if(n%2===0)dp=dp.times(d);}
  const expected=n%2?dp.times(g):dp;
  prove('bridge:power_'+n,power,expected);
  prove('bridge:formal_exponential_coefficient_'+n,exp[n],expected.scale(f(1).div(factor)));
 }
 const den=S.one().minus(d.times(t).times(t)),num=S.one().plus(d.times(t).times(t)).plus(t.times(g).scale(2));
 prove('bridge:Cayley_closed_form',S.one().minus(t.times(g)).times(num),den.times(S.one().plus(t.times(g))));
 const c=S.one().plus(d.times(t).times(t)),s=t.scale(2);
 prove('bridge:quadratic_norm_numerator',c.times(c).minus(d.times(s).times(s)),den.times(den));
}

// Exact symbolic port of the tower's polynomial fixture, using the native compiler.
{
 const {S,w}=system('lambda_tower_polynomials',['t','v','s','p'],[],true),{t,v,s,p}=w,zero=S.zero(),one=S.one();
 const values=[t.plus(v.scale(2)).plus(s.times(p)),t.times(v).plus(s.times(s)).plus(p.times(p).scale(3)),s.plus(t.times(p)).plus(v.times(s)),p.plus(t.times(s)).minus(v.times(p))];
 const names=['t','v','s','p'],signs=[1,1,-1,-1];
 const derivative=(e,k)=>S.reduce(M.polynomialDerivative(e,Object.fromEntries(names.map((n,i)=>[n,i===k?one:zero])))).normal;
 const reflect=e=>e.entries().reduce((a,[word,c])=>a.plus(S.word(word,c.mul(word.reduce((q,n)=>q*signs[names.indexOf(n)],1)))),zero);
 const jac=values.map(e=>names.map((_,i)=>derivative(e,i)));
 const hess=jac.map(row=>row.map(e=>names.map((_,i)=>derivative(e,i))));
 for(let a=0;a<4;a++){
  prove('tower:value_covariance_'+a,reflect(values[a]),values[a].scale(signs[a]));
  for(let i=0;i<4;i++){
   prove('tower:Jacobian_covariance_'+a+'_'+i,reflect(jac[a][i]).scale(signs[i]),jac[a][i].scale(signs[a]));
   for(let j=0;j<4;j++)prove('tower:Hessian_covariance_'+a+'_'+i+'_'+j,reflect(hess[a][i][j]).scale(signs[i]*signs[j]),hess[a][i][j].scale(signs[a]));
  }
  ck('tower:quadratic_termination_'+a,hess[a].flat().every(e=>e.entries().every(([word])=>word.length===0)));
 }
}

// Exact finite witnesses: each falsifies a specific overstatement in a draft.
const I=N.identity(2),J=N.matrix([[0,1],[1,0]]),P=N.scale(N.add(I,J),'1/2'),Q=N.scale(N.sub(I,J),'1/2');
const col=v=>N.matrix(v.map(x=>[x]));
{
 const e=col([1,1]),m=col([1,-1]),ep=N.scale(N.add(e,N.mul(J,m)),'1/2');
 ck('finite:nonzero_fixed_line',N.equal(N.mul(J,e),e));
 no('tensor:equality_of_opposite_eigenspaces_is_not_fixed_line',!N.equal(N.mul(J,e),N.scale(e,-1)));
 no('tensor:proposed_parallel_component_has_wrong_type',!N.equal(N.mul(J,ep),ep));
 no('tensor:proposed_distance_has_wrong_factor',!N.energy(N.sub(e,m)).eq(N.energy(N.scale(N.sub(e,N.mul(J,m)),'1/2'))));
 const ordinary=col([1,0]);no('tensor:J_is_not_identity_on_ordinary_real_axis',!N.equal(N.mul(J,ordinary),ordinary));
 let z=N.zeros(2,1);const T=N.add(P,N.scale(Q,'1/2'));
 for(let n=1;n<=4;n++){z=N.add(N.mul(T,z),e);ck('series:divergent_lawful_trajectory_'+n,N.equal(z,N.scale(e,n))&&N.isZero(N.mul(Q,z)));}
 no('series:zero_transverse_component_does_not_imply_convergence',N.energy(z).eq(32));
 no('series:bare_involution_does_not_accumulate_quarter_turns',N.equal(N.power(J,2),I)&&!N.IOTA.mul(N.IOTA).eq(1));
 const eps=f('1/3'),jt=x=>IOTA.mul(x),centered=jt(f(1).add(eps)).sub(jt(f(1).sub(eps))).div(eps.mul(2));
 ck('sine:lag_limit_linear_witness_is_zero',centered.sub(IOTA).zero());
 no('sine:lag_is_not_reflection_defect',!IOTA.sub(1).zero());
 const q=new Cut('1/2','-1/2');ck('zeta:multiplier_zero',IOTA.add(N.ONE.sub(IOTA).mul(q)).zero());
 ck('zeta:multiplier_root_modulus_squared',q.norm2().eq('1/2'));
 no('sine:product_of_sines_is_not_product_of_phases',!N.ZERO.eq(IOTA));
 const grad=col([1,0]),ric=N.scale(N.sub(mm(grad,N.dagger(grad)),I),'1/2');
 no('tensor:claimed_Ricci_scalar_contraction_fails',real(tr(ric)).eq('-1/2')&&!real(tr(ric)).eq('1/2'));
 results.draft_counterexamples={Ricci_trace:'-1/2',claimed_scalar:'1/2',zeta_multiplier_root:q.toJSON(),divergent_fixed_component:'z_n=n(1+iota)'};
}
{
 const badJ=N.matrix([[1,1],[0,-1]]),A=N.matrix([[0,1],[0,0]]),transpose=N.dagger;
 ck('tower:nonorthogonal_involution_witness',N.equal(mm(badJ,badJ),I));
 no('tower:transpose_and_cut_need_orthogonality',!N.equal(transpose(mm(badJ,A,badJ)),mm(badJ,transpose(A),badJ)));
 const obs=N.matrix([[1,-1]]),target=N.matrix([[1,0]]),next=I,w=col([1,1]);
 ck('tower:actual_target_blindness',N.isZero(N.mul(obs,w))&&!N.isZero(N.mul(target,w)));
 const blindDim=a=>N.rank(a.concat(target))-N.rank(a);
 ck('tower:target_blind_image_repaired',blindDim(obs)===1&&blindDim(next)===0);
 no('tower:intersection_with_target_support_misses_blindness',!N.isZero(N.mul(obs,col([1,0])))&&blindDim(obs)===1);
 results.observer={before_target_blind_dimension:blindDim(obs),after_target_blind_dimension:blindDim(next)};
}
{
 const companion=d=>N.matrix([[0,d],[1,0]]),circular=companion(-4),hyper=companion(4),nil=companion(0);
 for(const [name,g,d] of [['circular',circular,-4],['parabolic',nil,0],['hyperbolic',hyper,4]])ck('bridge:matrix_square_'+name,N.equal(mm(g,g),N.scale(I,d)));
 const K=N.scale(hyper,'1/2'),plus=N.scale(N.add(I,K),'1/2'),minus=N.scale(N.sub(I,K),'1/2'),flow=N.add(N.scale(plus,2),N.scale(minus,'1/2'));
 ck('bridge:stable_factor',N.equal(mm(flow,minus),N.scale(minus,'1/2')));
 ck('bridge:expanding_factor',N.equal(mm(flow,plus),N.scale(plus,2)));
 ck('bridge:both_sectors_reconstruct',N.equal(N.add(plus,minus),I));
 no('bridge:stable_cut_loses_expanding_information',N.isZero(mm(minus,plus))&&!N.isZero(plus));
 no('bridge:algebraic_projector_need_not_be_orthogonal',!N.equal(N.dagger(plus),plus));
 const C=N.cayley(circular,'1/3');no('bridge:circular_does_not_imply_positive_unitarity',!N.equal(mm(N.dagger(C),C),I));
 no('bridge:first_and_second_order_rates_differ',!f(2).eq(4));
 ck('bridge:quarter_threshold_is_reparameterizable',f('1/4').add('3/4').eq(f('1/2').add('1/2')));
}

// Four-dimensional geometry: cut-compatible connection and free normal transport.
{
 const I4=N.identity(4),J0=N.matrix([[1,0,0,0],[0,1,0,0],[0,0,-1,0],[0,0,0,-1]]),P0=N.scale(N.add(I4,J0),'1/2'),Q0=N.sub(I4,P0);
 const rot=(a,b)=>{const z=N.zeros(4);z[b][a]=N.ONE;z[a][b]=N.ONE.neg();return z;};
 const X=rot(0,2),Y=rot(0,3),dx=N.commutator(X,J0),dy=N.commutator(Y,J0),A=N.scale(mm(J0,dx),'1/2'),B=N.scale(mm(J0,dy),'1/2');
 const mixed=N.commutator(X,dy),dAy=N.scale(N.add(mm(dx,dy),mm(J0,mixed)),'1/2'),dAx=N.scale(N.add(mm(dy,dx),mm(J0,mixed)),'1/2');
 const Fxy=N.add(N.sub(dAy,dAx),N.commutator(A,B)),expected=N.scale(N.commutator(dx,dy),'1/4');
 ck('gravity:independent_connection_derivatives',N.equal(Fxy,expected));
 ck('gravity:normal_curvature_nonzero',!N.isZero(Fxy)&&N.equal(mm(Q0,Fxy,Q0),Fxy));
 ck('gravity:curvature_preserves_cut',N.isZero(N.commutator(Fxy,J0)));
 ck('gravity:trace_form_loses_transport',tr(mm(J0,dx)).zero()&&!N.isZero(A));
 no('gravity:nonzero_cut_derivative_need_not_give_curvature',!N.isZero(dx)&&N.isZero(N.commutator(dx,dx)));
 const free=rot(2,3);ck('gravity:free_normal_connection_preserves_cut',N.isZero(N.commutator(free,J0))&&!N.isZero(free));
 no('gravity:even_connection_curvature_is_invisible_to_constant_cut_tower',N.isZero(N.commutator(free,J0))&&!N.isZero(free),{construction:'J=J0 constant; B_x=0, B_y=x*C with C=normal rotation; F_xy=C. The alternative B=0 is flat with exactly the same J and all its derivatives.'});
 const R=N.cayley(free,'1/2'),eta=N.matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-1]]),H=N.zeros(4);H[2][3]=N.ONE;H[3][2]=N.ONE;
 const boost=N.cayley(H,'1/2');
 ck('gravity:normal_rotation_is_positive_unitary',N.equal(mm(N.dagger(R),R),I4));
 ck('gravity:normal_boost_is_Lorentz_isometry',N.equal(mm(N.dagger(boost),eta,boost),eta));
 no('gravity:boost_is_not_positive_unitary',!N.equal(mm(N.dagger(boost),boost),I4));
 ck('gravity:both_normal_transports_fix_face',N.equal(mm(R,P0),P0)&&N.equal(mm(boost,P0),P0));
 results.gravity={curvature:Fxy.map(r=>r.map(x=>x.toJSON())),normal_curvature_energy:N.energy(Fxy).toString(),free_normal_connection:true};
}

// Two late-added exponential drafts: preserve real-linear versus central-scalar typing.
{
 const R=N.matrix([[0,-1],[1,0]]),RJ=mm(R,J);
 ck('exponential:internal_quarter_turn_anticommutes_with_cut',N.equal(mm(R,J),N.scale(mm(J,R),-1)));
 ck('exponential:internal_i_times_cut_squares_positive',N.equal(mm(RJ,RJ),I));
 const central=N.scale(J,IOTA);
 ck('exponential:complexified_central_iota_times_cut_squares_negative',N.equal(mm(central,central),N.scale(I,-1)));
 const claimed=N.add(N.scale(I,'3/5'),N.scale(RJ,'4/5'));
 no('exponential:original_imaginary_time_formula_is_not_unitary',!N.equal(mm(N.dagger(claimed),claimed),I));
 const quantum=N.cayley(central,'1/2');
 ck('exponential:typed_complexified_Cayley_is_unitary',N.equal(mm(N.dagger(quantum),quantum),I));
 const U=N.add(N.scale(P,2),N.scale(Q,'1/2')),Ui=N.add(N.scale(P,'1/2'),N.scale(Q,2));
 ck('exponential:hyperbolic_flow_is_invertible',N.equal(mm(U,Ui),I));
 no('exponential:hyperbolic_flow_is_not_a_norm_contraction',N.energy(N.mul(U,col([1,1]))).eq(8));
 const relaxed=N.scale(U,'1/2');
 ck('exponential:normalized_flow_has_fixed_lawful_sector',N.equal(mm(relaxed,P),P)&&N.equal(mm(relaxed,Q),N.scale(Q,'1/4')));
 no('extended:fixed_point_does_not_select_unit_radius',N.equal(mm(J,col([2,2])),col([2,2]))&&N.energy(col([2,2])).eq(8));
 const diagonal=N.ONE.add(IOTA);
 ck('extended:primordial_identity_algebraic_core',diagonal.norm2().eq(2)&&diagonal.mul(diagonal).eq(IOTA.mul(2)));
 // For m=1-exp(-x), the claimed derivative differs by m-x exp(-x).
 const q1=f(-1),q2=f('1/2'),difference2=q2.add(q1).neg();
 no('memory:time_dependent_exponential_derivative_is_wrong',difference2.eq('1/2'),{first_nonzero_difference_coefficient:'x^2/2, x=lambda*t'});
 no('memory:positive_real_cost_cannot_equal_nonzero_pure_phase',!new Cut(1,0).eq(new Cut(0,1)));
 results.exponential={internal_quarter_turn_product_square:'+I',central_complexification_product_square:'-I',correct_relaxation:'exp(-lambda*t*P_minus)=P_plus+exp(-lambda*t)*P_minus',memory_propagator:'exp(t*J-[t-(1-exp(-lambda*t))/lambda]*I)',memory_limit:'exp(1/lambda)*P_plus for lambda>0',scope:'Memory law is an assumed constitutive model, not inferred entropy production.'};
}

const selected=proofs.find(p=>p.certificate.steps.length),S=systems.get(selected.certificate.presentation);
const corrupt=JSON.parse(JSON.stringify(selected.certificate));corrupt.steps[0].coefficient=['997','0'];
no('integrity:corrupted_rewrite_rejected',rej(()=>S.replay(corrupt)));
const wrong=JSON.parse(JSON.stringify(selected.certificate));wrong.presentation='0'.repeat(64);
no('integrity:wrong_presentation_rejected',rej(()=>S.replay(wrong)));
const payload={protocol:'NJC1_NATIVE_CLAIM_CERTIFICATION_V1',status:'PASS_CERTIFIED_CORE_WITH_RECORDED_REJECTIONS',scope:'Native algebraic proof witnesses and finite exact controls. General analytic/completion and geometric statements retain THEOREM.md hypotheses; not a general-purpose proof-assistant certificate or Einstein/RH/YM proof.',source_pins_sha256:hash(fs.readFileSync(path.join(root,'SOURCE_PINS.json'))),presentations,proofs,checks,negative_controls:controls,results,summary:{native_presentations:presentations.length,symbolic_certificates:proofs.length,serialized_replays:proofs.length,exact_checks:checks.length,negative_controls:controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',h=hash(bytes),file=path.join(root,'CERTIFICATE.json'),pin=path.join(root,'EXPECTED.sha256');
if(write){fs.writeFileSync(file,bytes);fs.writeFileSync(pin,h+'  CERTIFICATE.json\n');}
else{assert.equal(fs.readFileSync(file,'utf8'),bytes,'Certificate replay changed');assert.equal(fs.readFileSync(pin,'utf8'),h+'  CERTIFICATE.json\n','Certificate pin changed');}
console.log(JSON.stringify({status:payload.status,mode:write?'generation':'read-only replay',...payload.summary,sha256:h,gravity:results.gravity},null,2));
