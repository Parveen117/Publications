'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const T=require('./transport_adapter.cjs');
const L=require('../native-tvsp-lambda-matrices/tvsp_adapter.cjs');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const P=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const argv=process.argv.slice(2);assert.ok(argv.length===1&&['--write','--check'].includes(argv[0]));
const checks=[],negative_controls=[],fixtures=[],presentations=[],proofs=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const rejects=f=>{try{f();return false;}catch{return true;}};
const eq=N.equal,mul=N.mul,scale=N.scale,mat=N.matrix,cut=N.Cut.of,zero=N.isZero;
const p=T.project,c0=x=>T.coefficient(x,0),m0=L.matrixOuter;
const tr=M=>M.reduce((s,row,i)=>s.add(row[i]),N.ZERO);
const H0=mat([[2,1],[1,5]]),P0=mat([[1,'1/2'],[0,0]]),Q=mat([[1,0],['-1/5',0]]),I=N.identity(2),Z=N.zeros(2);
ck('fixture:rank_one_projector',eq(mul(P0,P0),P0)&&tr(P0).eq(1));
ck('fixture:H0_orthogonal_projector',eq(mul(N.dagger(P0),H0),mul(H0,P0)));
no('fixture:Euclidean_orthogonal_projector',!eq(P0,N.dagger(P0)));
for(const mu of ['-5','0','1/2','1','5']){
  const f=L.buildFixture(mu,4),curvature=T.sourceGeometry(f).Fsv;
  // alpha_i=R_i/sqrt(Delta), R_i=(b_i-(b/a)a_i)/2. Clearing
  // Delta^(3/2) allows a fully rational check of the exact NT-4 phase curl.
  const ba=L.project(L.div(f.b,f.a),3);
  const Rs=L.scale(L.sub(L.derivative(f.b,'s'),L.mul(ba,L.derivative(f.a,'s'))),'1/2');
  const Rv=L.scale(L.sub(L.derivative(f.b,'v'),L.mul(ba,L.derivative(f.a,'v'))),'1/2');
  const phaseNumerator=L.sub(L.mul(L.project(f.Delta,2),L.sub(L.derivative(Rv,'s'),L.derivative(Rs,'v'))),
    L.scale(L.sub(L.mul(L.project(Rv,2),L.project(L.derivative(f.Delta,'s'),2)),L.mul(L.project(Rs,2),L.project(L.derivative(f.Delta,'v'),2))),'1/2'));
  ck('orthonormal_phase:cleared_source_curl:'+mu,L.equal(phaseNumerator,L.scalar(phaseNumerator.context,cut(mu).mul(mu).mul('-5/4'))));
  for(const [name,pathSpec] of [['s',{a:1,b:0}],['v',{a:0,b:1}],['diagonal',{a:1,b:1}]]){
    const {H,A,U,r}=T.sourcePath(f,pathSpec),tag=mu+':'+name,C=U.context,depth=C.order,Ui=T.inverse(U);
    ck('path:metric_to_connection_depth:'+tag,A.context.order===depth-1);
    ck('path:ordered_ODE:'+tag,T.equal(T.derivative(U),T.scale(T.mul(A,p(U,depth-1)),-1)));
    ck('path:inverse:'+tag,T.equal(T.mul(U,Ui),T.identity(C,2))&&T.equal(T.mul(Ui,U),T.identity(C,2)));
    ck('path:affine_binomial_solution:'+tag,T.equal(U,T.affineBinomialTransport(H)));
    ck('path:metric_compatibility:'+tag,T.equal(T.derivative(H),T.add(T.mul(T.dagger(A),p(H,depth-1)),T.mul(p(H,depth-1),A))));
    ck('path:transported_metric:'+tag,T.equal(T.mul(T.dagger(U),T.mul(H,U)),T.constant(C,H0)));
    const observer=T.observer(U,mat([[2,3]])),w=T.mul(U,T.constant(C,mat([[1],[-2]])));
    ck('observer:co_moving_readout:'+tag,T.equal(T.mul(observer,w),T.constant(C,mat([[-4]]))));
    ck('observer:differential_cancellation:'+tag,T.equal(T.derivative(observer),T.mul(p(observer,depth-1),A)));
    const Pt=T.mul(U,T.mul(T.constant(C,P0),Ui));
    ck('cut:idempotence:'+tag,T.equal(T.mul(Pt,Pt),Pt));
    ck('cut:metric_adjoint:'+tag,T.equal(T.mul(T.dagger(Pt),H),T.mul(H,Pt)));
    ck('cut:rank_one_trace:'+tag,T.coefficients(Pt).every((M,k)=>tr(M).eq(k===0?1:0)));
    if(name==='s')ck('response:s_axis_constant_five:'+mu,T.equal(r,T.constant(C,mat([[5]]))));
    if(name==='v'){
      ck('response:v_axis_changes:'+mu,T.coefficient(r,1)[0][0].eq(5));
      const k=cut(mu).mul('5/9');let alpha=new N.F(1);const closed=[I];
      for(let n=1;n<=depth;n++){alpha=alpha.mul(new N.F(1-2*n,2*n));closed.push(scale(Q,cut(alpha).mul(new N.Cut(k.rad.pow(n)))));}
      ck('path:v_axis_projector_closed_form:'+mu,T.equal(U,T.series(C,closed)));
      if(mu!=='0')no('path:Euclidean_unitarity:'+mu,!T.equal(T.mul(T.dagger(U),U),T.identity(C,2)));
    }
    fixtures.push({kind:'straight_source_path',mu,path:name,U:T.serialize(U),response:T.serialize(r)});
  }
  const curved=T.sourcePath(f,{a:1,b:1,p:1,q:2}),ordinary=T.ordinaryExponential(T.scale(T.integrateZero(curved.A),-1));
  const defect=T.sub(curved.U,ordinary),F0=m0(curvature);
  ck('ordering:source_curve_equal_through_degree_two:'+mu,[0,1,2].every(k=>zero(T.coefficient(defect,k))));
  ck('ordering:source_curve_degree_three_F_over_six:'+mu,eq(T.coefficient(defect,3),scale(F0,'1/6')));
  if(mu!=='0')no('ordering:source_ordinary_exponential_replaces_ordering:'+mu,!zero(T.coefficient(defect,3)));
  const loop=T.smallLoop(f),loopU=loop.holonomy,LC=loopU.context;
  ck('loop:first_order_cancels:'+mu,zero(T.coefficient(loopU,1)));
  ck('loop:source_derivatives_plus_commutator:'+mu,eq(T.coefficient(loopU,2),scale(F0,-1)));
  const co=T.observer(loopU,mat([[1,0]])),section=T.mul(loopU,T.constant(LC,mat([[1],[0]])));
  ck('loop:co_moving_observer_cancels_holonomy:'+mu,T.equal(T.mul(co,section),T.constant(LC,mat([[1]]))));
  const fixed=T.mul(T.constant(LC,mat([[1,0]])),section);
  const transportedCut=T.mul(loopU,T.mul(T.constant(LC,P0),T.inverse(loopU)));
  ck('loop:cut_commutator:'+mu,eq(T.coefficient(transportedCut,2),scale(N.commutator(F0,P0),-1)));
  if(mu!=='0'){
    no('loop:fixed_observer_is_always_blind:'+mu,!zero(T.coefficient(fixed,2)));
    no('loop:co_moving_observer_reads_curvature:'+mu,zero(T.coefficient(T.mul(co,section),2)));
    no('loop:reverse_oriented_sign:'+mu,!eq(T.coefficient(loopU,2),F0));
  }
  const commuting=T.mul(loopU,T.mul(T.identity(LC,2),T.inverse(loopU)));
  ck('loop:commuting_cut_can_be_blind:'+mu,T.equal(commuting,T.identity(LC,2)));
  fixtures.push({kind:'curved_source_path',mu,U:T.serialize(curved.U),ordering_defect:T.serialize(defect),loop_coefficient:T.coefficient(loopU,2)});
}

// An exact endpoint lying in the stated positive local interval; no numeric
// Taylor approximation or finite-jet-to-endpoint inference is used here.
{
  const t=new N.F(369,10000),H=N.add(H0,scale(mat([[5,0],[0,0]]),t)),U=N.sub(I,scale(Q,'1/21')),Ui=N.inverse(U);
  ck('endpoint:exact_square_ratio',new N.F(1).add(new N.F(25,9).mul(t)).eq(new N.F(21,20).pow(2)));
  ck('endpoint:inside_positive_source_interval',t.le(new N.F(1,20))&&!t.eq(new N.F(1,20)));
  ck('endpoint:exact_metric_transport',eq(mul(N.dagger(U),mul(H,U)),H0));
  const Pt=mul(U,mul(P0,Ui));
  ck('endpoint:cut_idempotent',eq(mul(Pt,Pt),Pt));
  ck('endpoint:cut_H_adjoint',eq(mul(N.dagger(Pt),H),mul(H,Pt)));
  ck('endpoint:cut_trace_one',tr(Pt).eq(1));
  no('endpoint:U_is_Euclidean_unitary',!eq(mul(N.dagger(U),U),I));
  ck('endpoint:actual_response_changes',new N.F(5).mul(new N.F(1).add(t)).eq(new N.F(10369,2000)));
  fixtures.push({kind:'exact_rational_endpoint',mu:'5',t:t.toString(),H,U,cut:Pt,response:'10369/2000'});
}

// Ordered versus ordinary exponential for a genuinely noncommuting A=X+tY.
{
  const C=T.createContext(2,'independent ordering witness'),X=mat([[0,1],[0,0]]),Y=mat([[0,0],[1,0]]);
  const A=T.series(C,[X,Y,Z]),U=T.transport(A),ordinary=T.ordinaryExponential(T.scale(T.integrateZero(A),-1));
  ck('ordering:general_degree_three_commutator',eq(T.coefficient(T.sub(U,ordinary),3),scale(N.commutator(Y,X),'1/12')));
  no('ordering:wrong_commutator_sign',!eq(T.coefficient(T.sub(U,ordinary),3),scale(N.commutator(X,Y),'1/12')));
}

// Finite-depth blindness is shown with an actual integrable source potential,
// not an arbitrary symmetric matrix field incorrectly called a Hessian.
for(const depth of [0,1,2,3]){
  const sourceDepth=depth+3,ctx=L.createContext(sourceDepth,'blindness source '+depth),s=L.coordinate(ctx,'s'),v=L.coordinate(ctx,'v');
  const base=L.add(L.sub(L.scale(s,10),L.scale(v,10)),L.scale(L.add(L.power(s,2),L.power(v,2)),'1/2'));
  const u=L.add(base,L.scale(L.power(s,depth+3),new N.F(1,(depth+2)*(depth+3))));
  function actualResponses(potential){
    const Tfull=L.derivative(potential,'s'),Pfull=L.scale(L.derivative(potential,'v'),-1);
    const a=L.derivative(Tfull,'s'),b=L.derivative(Tfull,'v'),c=L.scale(L.derivative(Pfull,'v'),-1),C=a.context;
    const S=L.add(L.scalar(C,1),L.coordinate(C,'s')),V=L.add(L.scalar(C,1),L.coordinate(C,'v'));
    const T0=L.project(Tfull,C.order),P0=L.project(Pfull,C.order),Delta=L.sub(L.mul(a,c),L.mul(b,b));
    const responses={lambdaP:L.scale(L.div(L.mul(S,Delta),c),-1),lambdaV:L.scale(L.mul(S,a),-1),lambdaS:L.scale(L.mul(V,c),-1),
      zT:L.scale(L.div(L.mul(V,Delta),a),-1),zS:L.div(P0,c),lambdaT:L.div(L.mul(P0,a),Delta)};
    return {a,b,c,H:L.matrix([[a,b],[b,c]]),T:T0,P:P0,responses};
  }
  const perturbed=actualResponses(u),unperturbed=actualResponses(base),H=perturbed.H;
  ck('blindness:potential_Hessian:'+depth,L.equal(perturbed.a,L.add(L.scalar(perturbed.a.context,1),L.power(L.coordinate(perturbed.a.context,'s'),depth+1))));
  ck('blindness:baseline_T_and_P:'+depth,L.outer(perturbed.T).eq(10)&&L.outer(perturbed.P).eq(10));
  ck('blindness:same_retained_metric:'+depth,L.matrixEqual(L.matrixProject(H,depth),L.matrixProject(unperturbed.H,depth)));
  for(const name of Object.keys(perturbed.responses))ck('blindness:same_actual_response_Lambda_'+name+':'+depth,
    L.equal(L.project(perturbed.responses[name],depth),L.project(unperturbed.responses[name],depth)));
  const path=T.restrictMatrix(H,{a:1,b:0}),U=T.transport(T.connectionFromMetric(path));
  ck('blindness:same_transport_through_known_depth:'+depth,T.equal(T.project(U,depth),T.identity(T.lowerContext(U.context,depth),2)));
  no('blindness:next_transport_coefficient_is_determined:'+depth,T.coefficient(U,depth+1)[0][0].eq('-1/2'));
  const t=new N.F(1,2),metric11=new N.F(1).add(t.pow(depth+1));
  ck('blindness:distinct_analytic_endpoint_metric:'+depth,!metric11.eq(1));
  no('blindness:finite_jet_fixes_exact_endpoint:'+depth,!new N.F(1).div(metric11).eq(1));
}

// Actual covariant source layers are transported with their tensor slots.
// This checks the variation bridge, not merely inverse-frame cancellation.
{
  const f=L.buildFixture(1,4),pathSpec={a:1,b:1,p:1,q:2},z=T.sourcePath(f,pathSpec),C=z.U.context;
  const velocity=T.series(C,Array.from({length:5},(_,k)=>k===0?mat([[1],[0]]):k===1?mat([[0],[2]]):N.zeros(2,1)));
  const xi=T.mul(T.inverse(z.U),velocity),tower=L.responseTower(f,2);
  for(const name of ['lambdaS','lambdaV']){
    const gradient=tower[1].components[name],hessian=tower[2].components[name];
    const g=T.restrictMatrix(L.matrix([[gradient['0'],gradient['1']]]),pathSpec);
    const K=T.restrictMatrix(L.matrix([[hessian['00'],hessian['01']],[hessian['10'],hessian['11']]]),pathSpec);
    const beta=T.mul(g,T.project(z.U,3)),U2=T.project(z.U,2),Kbar=T.mul(T.dagger(U2),T.mul(K,U2));
    const scalarPath=T.restrictScalar(f.responses[name],pathSpec);
    ck('tensor_bridge:scalar_variation_equals_beta_xi:'+name,T.equal(T.derivative(scalarPath),T.mul(beta,T.project(xi,3))));
    ck('tensor_bridge:gradient_variation_equals_pulled_Hessian:'+name,
      T.equal(T.derivative(beta),T.mul(T.dagger(T.project(xi,2)),Kbar)));
    ck('tensor_bridge:explicit_depth_consumption:'+name,beta.context.order===3&&T.derivative(beta).context.order===2&&Kbar.context.order===2);
    const rawHessian=L.matrix([[L.derivative(gradient['0'],'s'),L.derivative(gradient['1'],'s')],
      [L.derivative(gradient['0'],'v'),L.derivative(gradient['1'],'v')]]);
    const rawK=T.restrictMatrix(rawHessian,pathSpec),rawKbar=T.mul(T.dagger(U2),T.mul(rawK,U2));
    no('tensor_bridge:omit_connection_corrections:'+name,!T.equal(T.derivative(beta),T.mul(T.dagger(T.project(xi,2)),rawKbar)));
  }
}

// Typed and representation controls: discarded layers are never silently zero.
{
  const f=L.buildFixture(1,3),z=T.sourcePath(f),U=z.U;
  no('depth:unknown_coefficient',rejects(()=>T.coefficient(U,4)));
  no('depth:invent_higher_layers',rejects(()=>T.project(U,4)));
  no('depth:mixed_depth_product',rejects(()=>T.mul(U,T.project(U,2))));
  no('depth:zero_depth_derivative',rejects(()=>T.derivative(T.project(U,0))));
  no('depth:insufficient_coefficient_array',rejects(()=>T.series(U.context,[I])));
  no('depth:source_origin_must_be_retained',rejects(()=>T.restrictMatrix(f.H,{a:1,b:1,p:0,q:1})));
  no('depth:independent_context_merge',rejects(()=>T.add(U,T.constant(T.createContext(3,U.context.label),I))));
  no('exponential:nonzero_constant_exponent',rejects(()=>T.ordinaryExponential(T.identity(U.context,2))));
  ck('immutability:records_and_coefficients',Object.isFrozen(U)&&Object.isFrozen(T.coefficients(U))&&Object.isFrozen(T.coefficient(U,0)[0]));
  const ctx=L.createContext(2),s=L.coordinate(ctx,'s'),R=L.regularMatrix(L.matrix([[s]]));
  ck('formal:coordinate_is_selfadjoint',L.equal(L.dagger(s),s));
  ck('formal:regular_coordinate_nonzero_nilpotent',!zero(R)&&zero(N.power(R,3)));
  no('formal:regular_coordinate_is_positive_Hilbert_adjoint',!eq(R,N.dagger(R)));
}

// Universal algebraic cancellations are separately serialized and replayed
// using the unchanged Paninian engine. These do not certify analytic domains.
function register(S,name){const audit=S.audit();ck('words:confluent:'+name,audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name,spec:S.spec(),audit});}
function prove(S,name,a,b){const r=S.reduce(a.minus(b),{witness:true});ck('words:'+name,r.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(r.certificate));ck('words:replay:'+name,S.replay(certificate));proofs.push({name,presentation_hash:S.hash(),lhs:a.toJSON(),rhs:b.toJSON(),certificate});}
{
  const S=new P.Presentation({tokens:['h','j','k'],rules:[{id:'left_inverse',lhs:['j','h'],rhs:[[[],1]],source:'declared metric inverse'},{id:'right_inverse',lhs:['h','j'],rhs:[[[],1]],source:'declared metric inverse'}]});register(S,'selected_metric_compatibility');
  const h=S.word(['h']),j=S.word(['j']),k=S.word(['k']),A=j.times(k).scale('1/2'),Ad=k.times(j).scale('1/2');
  prove(S,'Astar_H_plus_H_A_equals_Hprime',Ad.times(h).plus(h.times(A)),k);
}
{
  const S=new P.Presentation({tokens:['o','a','u','w'],rules:[]});register(S,'observer_differential_cancellation');const o=S.word(['o']),a=S.word(['a']),u=S.word(['u']),w=S.word(['w']);
  prove(S,'observer_times_transport_derivative_zero',o.times(a).times(u).minus(o.times(a.times(u))),S.zero());
  prove(S,'observer_reading_derivative_zero',o.times(a).times(w).minus(o.times(a.times(w))),S.zero());
}
{
  const S=new P.Presentation({tokens:['x','y','xs','xv','ys','yv'],rules:[]});register(S,'oriented_loop');
  const [x,y,xs,xv,ys,yv]=['x','y','xs','xv','ys','yv'].map(t=>S.word([t])),one=S.one(),z=S.zero();
  const product=(A,B)=>Array.from({length:3},(_,n)=>{let out=z;for(let j=0;j<=n;j++)out=out.plus(A[j].times(B[n-j]));return out;});
  const E1=[one,x.scale(-1),x.times(x).minus(xs).scale('1/2')];
  const E2=[one,y.scale(-1),y.times(y).minus(yv).scale('1/2').minus(ys)];
  const E3=[one,x,x.times(x).plus(xs).scale('1/2').plus(xv)];
  const E4=[one,y,y.times(y).plus(yv).scale('1/2')];
  const loop=product(E4,product(E3,product(E2,E1))),F=ys.minus(xv).plus(P.bracket(x,y));
  prove(S,'loop_first_order_zero',loop[1],z);prove(S,'oriented_loop_second_order_minus_F',loop[2],F.scale(-1));
  no('words:constant_commutator_alone_is_curvature',S.reduce(loop[2].plus(P.bracket(x,y))).normal.terms.size>0);
}
{
  const S=new P.Presentation({tokens:['x','y'],rules:[]});register(S,'ordering_defect');const x=S.word(['x']),y=S.word(['y']);
  const ordered=x.times(x).times(x).scale('-1/6').plus(x.times(y).scale('1/6')).plus(y.times(x).scale('1/3'));
  const ordinary=x.times(x).times(x).scale('-1/6').plus(x.times(y).plus(y.times(x)).scale('1/4'));
  prove(S,'third_order_defect_YX_minus_XY_over_twelve',ordered.minus(ordinary),P.bracket(y,x).scale('1/12'));
}
{
  const tokens=['a','b','m'],rules=[];
  for(let j=1;j<tokens.length;j++)for(let i=0;i<j;i++)rules.push({id:'commute_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'central thermodynamic source coefficients'});
  const S=new P.Presentation({tokens,rules});register(S,'orthonormal_phase_curl');
  const [a,b,m]=tokens.map(t=>S.word([t])),delta=a.scale(5).minus(b.times(b)),m2=m.times(m);
  const cleared=m2.times(delta).scale('-1/2').minus(m2.times(b).times(b).scale('1/2')).plus(m2.times(a).scale('5/4'));
  prove(S,'a_times_phase_curl_numerator_equals_minus_five_mu_squared_a_over_four',cleared,m2.times(a).scale('-5/4'));
}
const payload={protocol:'TET1_TVSP_EXPONENTIAL_TRANSPORT_V1',status:'PASS_TET1_TVSP_EXPONENTIAL_TRANSPORT',source_parent:'babb4fb1922bd91c2eadc82abd4b9538a19745d0',checks,negative_controls,presentations,proofs,fixtures,
  scope:{source:'TLM six-response Hessian source with selected strain connection and declared local paths',transport:'finite exact ordered coefficients and exact rational endpoint fixture',metric_transport:true,co_moving_observer_cancellation:true,source_curved_path_ordering_defect:true,oriented_source_loop_sign:true,finite_depth_determines_analytic_endpoint:false,Euclidean_unitarity:false,formal_regular_jets_positive_star_representation:false,physical_space_quantized:false,quantum_gravity_closed:false,universal_claims_from_finite_tests:false},
  summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,negative_controls:negative_controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',hash=crypto.createHash('sha256').update(bytes).digest('hex'),cp=path.join(__dirname,'CERTIFICATE.json'),hp=path.join(__dirname,'EXPECTED.sha256');
if(argv[0]==='--write'){fs.writeFileSync(cp,bytes);fs.writeFileSync(hp,hash+'  CERTIFICATE.json\n');}
else{assert.equal(fs.readFileSync(cp,'utf8'),bytes);assert.equal(fs.readFileSync(hp,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
