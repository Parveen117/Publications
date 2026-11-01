'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const R=require('./response_adapter.cjs');
const L=require('../native-tvsp-lambda-matrices/tvsp_adapter.cjs');
const T=require('../native-tvsp-exponential-transport/transport_adapter.cjs');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const P=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const argv=process.argv.slice(2);assert.ok(argv.length===1&&['--write','--check'].includes(argv[0]));
const checks=[],negative_controls=[],fixtures=[],presentations=[],proofs=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const rejects=f=>{try{f();return false;}catch{return true;}};
const atoms=[[1,1],[2,1],[1,2]],mx=L.matrixMul,md=L.matrixDagger,mp=L.matrixProject,meq=L.matrixEqual,m0=L.matrixOuter;
const tr=M=>M.reduce((s,row,i)=>s.add(row[i]),N.ZERO);
const primaryAmplitudes=['1/3','2/3','2/3'];
for(const [label,amplitudes] of [['primary',primaryAmplitudes],['permuted',['2/3','1/3','2/3']]]){
  const z=R.finiteSource(amplitudes,atoms,3,'QRT '+label),g=R.responseGeometry(z),depth=z.context.order,low=depth-1;
  const mom=R.moments(z.weights,atoms),E=z.E,Ed=md(E),psi=z.psi,H=z.H,Q=g.Q;
  ck('source:state_normalization:'+label,meq(mx(md(psi),psi),L.identity(z.context,1)));
  ck('source:raw_state_norm_is_partition:'+label,L.equal(mx(md(z.unnormalized),z.unnormalized).rows[0][0],z.partition));
  no('source:delete_partition_normalization:'+label,!meq(mx(md(z.unnormalized),z.unnormalized),L.identity(z.context,1)));
  ck('source:centered_frame_orthogonal_to_state:'+label,meq(mx(md(psi),E),L.matrix([[L.scalar(z.context),L.scalar(z.context)]])));
  ck('source:response_Gram_is_covariance:'+label,meq(mx(Ed,E),H));
  ck('source:exact_point_covariance:'+label,N.equal(m0(H),mom.H));
  ck('source:point_means:'+label,z.means.every((x,i)=>L.outer(x).eq(mom.means[i])));
  no('source:covariance_is_raw_second_moment:'+label,!N.equal(mom.H,mom.raw));
  ck('source:strict_covariance_determinant:'+label,!mom.Delta.rad.le(0));
  const third=[L.derivative(z.a,'s'),L.derivative(z.a,'v'),L.derivative(z.c,'s'),L.derivative(z.c,'v')];
  ck('source:third_derivatives_are_negative_cumulants:'+label,third.every((x,i)=>L.outer(x).eq(mom.third[i].neg())));
  for(const coordinate of ['s','v']){
    const dE=L.matrixDerivative(E,coordinate),dH=L.matrixDerivative(H,coordinate),A=coordinate==='s'?g.As:g.Av;
    ck('source:derived_connection_half_factor:'+label+':'+coordinate,meq(mx(mp(Ed,low),dE),L.matrixScale(dH,'1/2')));
    ck('source:connection_from_response_frame:'+label+':'+coordinate,meq(mx(mp(g.Hinv,low),mx(mp(Ed,low),dE)),A));
    no('source:drop_connection_half:'+label+':'+coordinate,!meq(mx(mp(Ed,low),dE),dH));
    no('source:reverse_connection_sign:'+label+':'+coordinate,!meq(mx(mp(Ed,low),dE),L.matrixScale(dH,'-1/2')));
    const dpsi=L.matrixDerivative(psi,coordinate),index=coordinate==='s'?0:1;
    ck('source:normalized_state_derivative:'+label+':'+coordinate,meq(dpsi,L.matrixScale(L.matrix(mp(E,low).rows.map(row=>[row[index]])),'-1/2')));
    const berry=mx(mp(md(psi),low),dpsi);
    ck('source:base_line_Berry_connection_zero:'+label+':'+coordinate,meq(berry,L.matrix([[L.scalar(berry.context)]])));
  }
  ck('projection:idempotence:'+label,meq(mx(Q,Q),Q));
  ck('projection:self_adjoint:'+label,meq(md(Q),Q));
  ck('projection:fixes_response_frame:'+label,meq(mx(Q,E),E));
  ck('projection:orthogonal_to_base_state:'+label,meq(mx(Q,psi),L.matrix(psi.rows.map(()=>[L.scalar(z.context)]))));
  ck('projection:categorical_complement_only_in_three_atom_fixture:'+label,meq(Q,L.matrixSub(L.identity(z.context,3),mx(psi,md(psi)))));
  ck('projection:rank_two_trace:'+label,tr(m0(Q)).eq(2));
  const HF=mx(mp(H,depth-2),g.Fsv),curvature=L.div(HF.rows[0][1],L.project(z.Delta,depth-2));
  ck('curvature:NT_commutator_equals_cumulant_formula:'+label,L.outer(curvature).eq(mom.curvature));
  ck('curvature:whole_retained_categorical_family_quarter:'+label,L.equal(curvature,L.scalar(curvature.context,'1/4')));
  no('curvature:zero_base_Berry_forces_zero_response_curvature:'+label,!N.isZero(m0(g.Fsv)));
  const pathData=R.responsePath(z,{a:1,b:1,p:1,q:2}),{K,Kgen,U,A}=pathData,C=K.context,pE=T.project(pathData.E,low),pQ=T.project(pathData.Q,low);
  ck('Kato:projection_derivative_tangent:'+label,T.equal(T.add(T.mul(pathData.Qprime,pQ),T.mul(pQ,pathData.Qprime)),pathData.Qprime));
  ck('Kato:generator_anti_adjoint:'+label,T.equal(T.dagger(Kgen),T.scale(Kgen,-1)));
  ck('Kato:generator_moves_projection:'+label,T.equal(T.sub(T.mul(Kgen,pQ),T.mul(pQ,Kgen)),pathData.Qprime));
  ck('Kato:response_frame_horizontal_identity:'+label,T.equal(T.mul(Kgen,pE),T.sub(T.derivative(pathData.E),T.mul(pE,A))));
  ck('Kato:ambient_ODE:'+label,T.equal(T.derivative(K),T.mul(Kgen,T.project(K,low))));
  ck('Kato:ambient_unitarity:'+label,T.equal(T.mul(T.dagger(K),K),T.identity(C,3)));
  ck('Kato:ambient_inverse_unitarity:'+label,T.equal(T.mul(K,T.dagger(K)),T.identity(C,3)));
  ck('Kato:exact_response_intertwiner:'+label,T.equal(T.mul(K,T.constant(C,m0(E))),T.mul(pathData.E,U)));
  ck('Kato:projector_intertwiner:'+label,T.equal(T.mul(K,T.mul(T.constant(C,m0(Q)),T.dagger(K))),pathData.Q));
  ck('Kato:response_metric_transport:'+label,T.equal(T.mul(T.dagger(U),T.mul(pathData.H,U)),T.constant(C,m0(H))));
  no('Kato:coordinate_transport_is_Euclidean_unitary:'+label,!T.equal(T.mul(T.dagger(U),U),T.identity(C,2)));
  no('Kato:response_space_is_extra_qubit:'+label,K.height===atoms.length&&U.height===2&&E.height===atoms.length);
  const loop=T.smallLoop(z);
  ck('curvature:actual_source_loop_sign:'+label,N.equal(T.coefficient(loop.holonomy,2),N.scale(m0(g.Fsv),-1)));
  no('curvature:actual_source_loop_flat:'+label,!N.isZero(T.coefficient(loop.holonomy,2)));
  fixtures.push({kind:'finite_exact_spectral_comparison_only',label,amplitudes,atoms,weights:z.weights,
    point_moments:mom,point_response_frame:m0(E),point_response_projection:m0(Q),K:T.serialize(K),U:T.serialize(U),loop_coefficient:T.coefficient(loop.holonomy,2)});
}

// Probability models need no square-root amplitudes for exact cumulants.
for(const probabilities of [['1/3','1/3','1/3'],['1/9','4/9','4/9'],['4/9','1/9','4/9'],['1/6','1/3','1/2'],['1/10','3/10','3/5']]){
  const m=R.moments(probabilities,atoms),tag=probabilities.join(',');
  ck('categorical:curvature_quarter:'+tag,m.curvature.eq('1/4'));
  ck('categorical:determinant_identity:'+tag,m.D.eq(m.Delta.mul(m.Delta).neg()));
  if(probabilities[0]==='1/3'){
    ck('categorical:equal_weight_H',N.equal(m.H,N.scale(N.matrix([[2,-1],[-1,2]]),'1/9')));
    ck('categorical:equal_weight_u3',m.u3.every((x,i)=>x.eq(new N.F([-2,1,1,-2][i],27))));
    ck('categorical:equal_weight_Delta',m.Delta.eq('1/27'));
    ck('categorical:equal_weight_D',m.D.eq('-1/729'));
  }
  fixtures.push({kind:'exact_probability_cumulant_model',probabilities,moments:m});
}

// A four-atom comparison prevents the curvature test from checking only the
// special categorical constant-curvature family.
{
  const support=[[1,1],[2,1],[1,2],[3,4]],z=R.finiteSource(['1/2','1/2','1/2','1/2'],support,2,'generic four-atom source');
  const g=R.responseGeometry(z),m=R.moments(z.weights,support),HF=mx(mp(z.H,0),g.Fsv),kappa=L.outer(HF.rows[0][1]).div(L.outer(z.Delta));
  ck('generic_support:response_Gram',meq(mx(md(z.E),z.E),z.H));
  ck('generic_support:covariance_cumulant_curvature',kappa.eq(m.curvature));
  no('generic_support:all_positive_sources_have_quarter_curvature',!m.curvature.eq('1/4'));
  fixtures.push({kind:'generic_four_atom_curvature_comparison',support,moments:m,connection_curvature:kappa});
}

// Fixed orthonormal response-plane detectors read signed orientation. This
// is the already-derived response plane, with no added tensor factor.
{
  const Pplus=N.scale(N.matrix([[1,1],[1,1]]),'1/2'),P2=N.matrix([[0,0],[0,1]]),initial=N.matrix([[1],[0]]),c=N.Cut.of('3/5');
  ck('observer:forty_five_degree_projector',N.equal(N.mul(Pplus,Pplus),Pplus)&&N.equal(N.dagger(Pplus),Pplus));
  const readings=[];
  for(const sign of [1,-1]){
    const s=N.Cut.of(new N.F(4*sign,5)),W=N.matrix([[c,s.neg()],[s,c]]),state=N.mul(W,initial);
    const pplus=N.mul(N.dagger(state),N.mul(Pplus,state))[0][0],p2=N.mul(N.dagger(state),N.mul(P2,state))[0][0];
    ck('observer:declared_response_rotation_unitary:'+sign,N.equal(N.mul(N.dagger(W),W),N.identity(2)));
    ck('observer:pplus_equals_half_plus_cs:'+sign,pplus.eq(c.mul(s).add('1/2')));
    ck('observer:exact_pplus:'+sign,pplus.eq(sign===1?'49/50':'1/50'));
    ck('observer:exact_p2:'+sign,p2.eq('16/25'));
    readings.push({sign,c,s,pplus,p2});
  }
  no('observer:p2_recovers_orientation_sign',readings[0].p2.eq(readings[1].p2));
  ck('observer:fixed_diagonal_detector_distinguishes_sign',!readings[0].pplus.eq(readings[1].pplus));
  fixtures.push({kind:'orthonormal_response_plane_fixed_detector',projector:Pplus,readings});
}

// These rational inequalities certify the propagation of analytically proved
// support/preparation error bounds. No continuum distribution is discretized.
{
  const b=R.continuumBounds(),F=N.F.of;
  ck('unbinned:half_width',b.delta.eq('1/10000'));
  ck('unbinned:covariance_error',b.covarianceError.eq(b.delta.mul(4)));
  ck('unbinned:third_central_error',b.thirdError.eq(b.delta.mul(6)));
  ck('unbinned:covariance_entry_margin',F('2/9').add(b.covarianceError).le(F('1/4')));
  ck('unbinned:third_central_entry_margin',F('2/27').add(b.thirdError).le(F('1/12')));
  ck('unbinned:positive_covariance_determinant_margin',!b.deltaLower.le(0));
  ck('unbinned:negative_cumulant_determinant_margin',!b.determinantLower.le(0));
  ck('unbinned:curvature_lower_strictly_above_eighth',F('1/8').le(b.curvatureLower)&&!F('1/8').eq(b.curvatureLower));
  ck('unbinned:curvature_upper_strictly_below_three_eighths',b.curvatureUpper.le(F('3/8'))&&!b.curvatureUpper.eq(F('3/8')));
  no('unbinned:unbounded_window_size_is_certified',rejects(()=>R.continuumBounds('1/10')));
  fixtures.push({kind:'exact_rational_bound_for_analytic_unbinned_preparation',assumptions:['three positive-mass spectral cells','coordinate half-width delta','equal normalized cell probabilities','analytically established moment-error bounds'],bounds:b});
}

// Domain failures and distinctions are controls, not omitted hypotheses.
{
  const redundant=R.finiteSource(primaryAmplitudes,[[1,2],[2,4],[3,6]],2,'redundant observable control');
  ck('domain:redundant_observable_covariance_singular',L.outer(redundant.Delta).zero());
  no('domain:singular_response_inverse',rejects(()=>R.responseGeometry(redundant)));
  no('domain:unnormalized_initial_amplitudes',rejects(()=>R.finiteSource([1,1,1],atoms,2)));
  no('domain:nonpositive_observable',rejects(()=>R.moments(['1/3','1/3','1/3'],[[0,1],[2,1],[1,2]])));
  no('domain:invent_source_jet_depth',rejects(()=>L.project(redundant.a,3)));
  no('domain:zero_depth_derivative',rejects(()=>L.derivative(L.project(redundant.a,0),'s')));
}

function register(S,name){const audit=S.audit();ck('words:confluent:'+name,audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name,spec:S.spec(),audit});}
function prove(S,name,a,b){const r=S.reduce(a.minus(b),{witness:true});ck('words:'+name,r.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(r.certificate));ck('words:replay:'+name,S.replay(certificate));proofs.push({name,presentation_hash:S.hash(),lhs:a.toJSON(),rhs:b.toJSON(),certificate});}
{
  const S=new P.Presentation({objects:['response','ambient'],tokens:[
    {name:'e',from:'response',to:'ambient'},{name:'f',from:'ambient',to:'response'},
    {name:'h',from:'response',to:'response'},{name:'j',from:'response',to:'response'}],rules:[
    {id:'response_Gram',lhs:['f','e'],rhs:[[['h'],1]],source:'f=E adjoint and H=E adjoint E'},
    {id:'metric_right_inverse',lhs:['h','j'],rhs:[[[],1]],source:'j=H inverse'},
    {id:'metric_left_inverse',lhs:['j','h'],rhs:[[[],1]],source:'j=H inverse'}]});
  register(S,'typed_response_projection');const e=S.word(['e']),f=S.word(['f']),j=S.word(['j']),q=e.times(j).times(f);
  prove(S,'response_projection_idempotent',q.times(q),q);
  prove(S,'response_projection_fixes_frame',q.times(e),e);
}
{
  const S=new P.Presentation({tokens:['q','d'],rules:[
    {id:'projection',lhs:['q','q'],rhs:[[['q'],1]],source:'Q squared equals Q'},
    {id:'projection_tangent',lhs:['d','q'],rhs:[[['d'],1],[['q','d'],-1]],source:'derivative of Q squared equals Q'}]});
  register(S,'Kato_projection_tangent');const q=S.word(['q']),d=S.word(['d']),k=P.bracket(d,q);
  prove(S,'projection_sandwich_derivative_zero',q.times(d).times(q),S.zero());
  prove(S,'Kato_double_commutator_is_projection_derivative',P.bracket(k,q),d);
  prove(S,'Kato_generator_anti_adjoint',P.star(k,{q,d}),k.scale(-1));
}
{
  const S=new P.Presentation({tokens:['q','r'],rules:[{id:'central_probabilities',lhs:['r','q'],rhs:[[['q','r'],1]],source:'commuting categorical probabilities'}]});
  register(S,'categorical_curvature');const q=S.word(['q']),r=S.word(['r']),one=S.one();
  const a=q.times(one.minus(q)),b=q.times(r).scale(-1),c=r.times(one.minus(r));
  const h111=a.times(one.minus(q.scale(2))).scale(-1),h112=q.times(r).times(one.minus(q.scale(2))),h122=q.times(r).times(one.minus(r.scale(2))),h222=c.times(one.minus(r.scale(2))).scale(-1);
  const delta=a.times(c).minus(b.times(b));
  const determinant=a.times(h112.times(h222).minus(h122.times(h122))).minus(b.times(h111.times(h222).minus(h122.times(h112)))).plus(c.times(h111.times(h122).minus(h112.times(h112))));
  prove(S,'categorical_Delta_equals_pqr',delta,q.times(r).times(one.minus(q).minus(r)));
  prove(S,'categorical_D_equals_minus_Delta_squared',determinant,delta.times(delta).scale(-1));
}
const payload={protocol:'QRT1_QUANTUM_RESPONSE_TRANSPORT_V1',status:'PASS_QRT1_QUANTUM_RESPONSE_TRANSPORT',source_parent:'9a0621637b2cc12fb57eb6fd36e56daef8cf5f34',checks,negative_controls,presentations,proofs,fixtures,
  scope:{exact_examples:'finite positive spectral comparison models with rational amplitudes or probabilities',continuum:'only rational propagation of separately proved unbinned support/preparation bounds',normalized_actual_response_Gram:true,derived_connection_factor_one_half:true,Kato_unitary_response_intertwiner:true,base_line_Berry_distinct_from_response_curvature:true,extra_qubit_added:false,coordinate_U_Euclidean_unitary:false,finite_models_prove_GMC_support:false,spacetime_curvature_established:false,quantum_gravity_closed:false,universal_claims_from_finite_checks:false},
  summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,negative_controls:negative_controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',hash=crypto.createHash('sha256').update(bytes).digest('hex'),cp=path.join(__dirname,'CERTIFICATE.json'),hp=path.join(__dirname,'EXPECTED.sha256');
if(argv[0]==='--write'){fs.writeFileSync(cp,bytes);fs.writeFileSync(hp,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cp,'utf8'),bytes);assert.equal(fs.readFileSync(hp,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
