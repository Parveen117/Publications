'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const PCR=require('./phase_current_adapter.cjs');
const QRT=require('../native-quantum-response-transport/response_adapter.cjs');
const L=require('../native-tvsp-lambda-matrices/tvsp_adapter.cjs');
const T=require('../native-tvsp-exponential-transport/transport_adapter.cjs');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const P=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const argv=process.argv.slice(2);assert.ok(argv.length===1&&['--write','--check'].includes(argv[0]));
const checks=[],negative_controls=[],fixtures=[],presentations=[],proofs=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const rejects=f=>{try{f();return false;}catch{return true;}};
const tr=M=>M.reduce((sum,row,i)=>sum.add(row[i]),N.ZERO),eq=N.equal,mm=N.mul,md=N.dagger,c0=x=>T.coefficient(x,0);
function permutations(xs,prefix=[]){if(!xs.length)return [prefix];return xs.flatMap((x,i)=>permutations(xs.filter((_,j)=>j!==i),prefix.concat(x)));}
function det4(M){let result=N.ZERO;for(const p of permutations([0,1,2,3])){let v=N.ONE,inversions=0;for(let i=0;i<4;i++){v=v.mul(M[i][p[i]]);for(let j=i+1;j<4;j++)if(p[i]>p[j])inversions++;}result=result.add(inversions%2?v.neg():v);}return result;}
function checkGeometry(H,means,label){
  for(const h of ['1/3','1/2','3/4']){
    const z=PCR.quantumGeometry(H,h),tag=label+':'+h,hf=N.F.of(h),I=N.identity(2);
    ck('QGT:Hermitian:'+tag,eq(z.QGT,md(z.QGT)));
    ck('QGT:source_factorization:'+tag,eq(z.QGT,mm(md(z.factor),mm(H,z.factor))));
    ck('QGT:complex_rank_two:'+tag,N.rank(z.QGT)===2);
    no('QGT:four_independent_complex_directions:'+tag,N.rank(z.QGT)!==4);
    const nullspace=N.scale(I,N.IOTA.mul(new N.F(-2).div(hf))).concat(I);
    ck('QGT:explicit_complex_nullspace:'+tag,N.isZero(mm(z.QGT,nullspace))&&N.rank(nullspace)===2);
    ck('metric:real_rank_four:'+tag,N.rank(z.g)===4);
    ck('metric:positive_first_covariance_block:'+tag,!z.g[0][0].rad.le(0)&&!z.g[0][0].mul(z.g[1][1]).sub(z.g[0][1].mul(z.g[1][0])).rad.le(0));
    ck('metric:positive_second_covariance_block:'+tag,!z.g[2][2].rad.le(0)&&!z.g[2][2].mul(z.g[3][3]).sub(z.g[2][3].mul(z.g[3][2])).rad.le(0));
    ck('complex_structure:square_minus_identity:'+tag,eq(mm(z.J,z.J),N.scale(N.identity(4),-1)));
    ck('complex_structure:metric_compatibility:'+tag,eq(mm(md(z.J),mm(z.g,z.J)),z.g));
    ck('symplectic:omega_equals_minus_half_Berry_F:'+tag,eq(z.omega,N.scale(z.BerryF,'-1/2')));
    ck('symplectic:Berry_form_nondegenerate:'+tag,N.rank(z.BerryF)===4);
    ck('symplectic:antisymmetry:'+tag,eq(md(z.BerryF),N.scale(z.BerryF,-1)));
    const minusTwiceImaginary=z.QGT.map(row=>row.map(x=>new N.Cut(x.turn.mul(-2))));
    ck('symplectic:Berry_F_is_minus_twice_Im_QGT:'+tag,eq(z.BerryF,minusTwiceImaginary));
    no('symplectic:opposite_Berry_sign:'+tag,!eq(z.omega,N.scale(z.BerryF,'1/2')));
    const readout=PCR.readoutGeometry(means,H,['3/5','-2/7'],h);
    ck('readout:canonical_pullback_Berry_F:'+tag,eq(readout.pullback,z.BerryF));
    ck('readout:mp_Jacobian_determinant:'+tag,det4(readout.mpJacobian).eq(readout.mpDeterminant));
    ck('readout:qp_Jacobian_determinant:'+tag,det4(readout.qpJacobian).eq(readout.qpDeterminant));
    ck('readout:positive_local_invertibility:'+tag,!readout.mpDeterminant.rad.le(0));
    ck('readout:phase_recovery_on_fixed_family:'+tag,readout.betaRecovered.every((x,i)=>x.eq(readout.beta[i])));
    fixtures.push({kind:'exact_point_phase_density_geometry',label,h,H,means,geometry:z,readout});
  }
}
const probabilityModel=QRT.moments(['1/3','1/3','1/3'],[[1,1],[2,1],[1,2]]);
checkGeometry(probabilityModel.H,probabilityModel.means,'equal_categorical_probabilities');
const cases=[
  {label:'three_atom_phase_family',amplitudes:['1/3','2/3','2/3'],atoms:[[1,1],[2,1],[1,2]],order:3,h:'1/2'},
  {label:'four_atom_phase_family',amplitudes:['1/2','1/2','1/2','1/2'],atoms:[[1,1],[2,1],[1,2],[3,4]],order:2,h:'1/3'}
];
for(const f of cases){
  const source=QRT.finiteSource(f.amplitudes,f.atoms,f.order,f.label),H0=L.matrixOuter(source.H),m0=source.means.map(L.outer),E0=L.matrixOuter(source.E);
  checkGeometry(H0,m0,f.label);
  const point=PCR.quantumGeometry(H0,f.h),tangents=E0.map(row=>row.map(x=>x.mul('-1/2')).concat(row.map(x=>x.mul(N.IOTA).mul(new N.F(-1).div(f.h)))));
  ck('source:actual_horizontal_tangent_Gram:'+f.label,eq(mm(md(tangents),tangents),point.QGT));
  for(let i=0;i<2;i++)for(let j=0;j<2;j++)ck('source:Berry_curvature_from_mean_derivatives:'+f.label+':'+i+j,
    L.equal(L.derivative(source.means[i],j===0?'s':'v'),L.scale(L.project(source.H.rows[i][j],f.order-1),-1)));
  const z=PCR.phaseCurrentPath(source,{h:f.h,c:'3/5',d:'-2/7'}),C=z.lift.context,n=f.atoms.length,depth=f.order-1,h=N.F.of(f.h),psi=T.project(z.psi,depth);
  ck('phase:unitarity:'+f.label,T.equal(T.mul(T.dagger(z.V),z.V),T.identity(C,n)));
  ck('phase:state_normalization:'+f.label,T.equal(T.mul(T.dagger(z.psi),z.psi),T.identity(C,1)));
  ck('phase:Berry_connection_actual:'+f.label,T.equal(z.BerryA,z.expectedBerry));
  ck('phase:Berry_connection_is_real:'+f.label,T.equal(z.BerryA,T.dagger(z.BerryA)));
  no('phase:Berry_connection_stays_zero:'+f.label,N.isZero(c0(z.BerryA))===false);
  ck('phase:horizontal_norm_equals_real_metric:'+f.label,T.equal(T.mul(T.dagger(z.horizontal),z.horizontal),z.horizontalNormExpected));
  ck('phase:horizontal_derivative_orthogonal:'+f.label,T.equal(T.mul(T.dagger(psi),z.horizontal),T.constant(z.BerryA.context,[[0]])));
  for(let i=0;i<2;i++){
    const Xi=N.zeros(n);for(let k=0;k<n;k++)Xi[k][k]=source.atoms[k][i];
    const densityMean=T.mul(T.dagger(z.psi),T.mul(T.constant(C,Xi),z.psi));
    ck('phase:density_readout_unchanged:'+f.label+':'+i,T.equal(densityMean,z.meanPaths[i]));
    ck('phase:specified_current_readout_beta_times_m:'+f.label+':'+i,T.equal(z.pReadouts[i],T.mul(z.betaPaths[i],densityMean)));
    ck('phase:recover_beta_from_mean_current_pair:'+f.label+':'+i,T.equal(z.betaRecovered[i],z.betaPaths[i]));
    // Tr([P,X])=0 for every finite pair, while Tr(i h X) is nonzero here.
    no('current:finite_matrix_affine_representation:'+f.label+':'+i,!tr(N.scale(Xi,N.IOTA.mul(h))).zero());
  }
  const Ephase=T.mul(z.V,z.base.E),Qphase=T.mul(z.V,T.mul(z.base.Q,T.dagger(z.V)));
  ck('phase:response_Gram_unchanged:'+f.label,T.equal(T.mul(T.dagger(Ephase),Ephase),z.base.H));
  ck('lift:unitarity:'+f.label,T.equal(T.mul(T.dagger(z.lift),z.lift),T.identity(C,n)));
  ck('lift:normalized_state_transport:'+f.label,T.equal(T.mul(z.lift,T.constant(C,c0(z.base.psi))),z.psi));
  ck('lift:response_frame_transport:'+f.label,T.equal(T.mul(z.lift,T.constant(C,c0(z.base.E))),T.mul(Ephase,z.base.U)));
  ck('lift:response_projector_transport:'+f.label,T.equal(T.mul(z.lift,T.mul(T.constant(C,c0(z.base.Q)),T.dagger(z.lift))),Qphase));
  ck('lift:Hamiltonian_Hermitian:'+f.label,T.equal(z.liftGenerator,T.dagger(z.liftGenerator)));
  ck('lift:Schrodinger_ODE:'+f.label,T.equal(T.derivative(z.lift),T.scale(T.mul(z.liftGenerator,T.project(z.lift,depth)),N.IOTA.mul(new N.F(-1).div(h)))));
  ck('lift:state_Schrodinger_ODE:'+f.label,T.equal(T.derivative(z.psi),T.scale(T.mul(z.liftGenerator,psi),N.IOTA.mul(new N.F(-1).div(h)))));
  const energy=T.mul(T.dagger(psi),T.mul(z.liftGenerator,psi));
  ck('lift:energy_mean_is_h_Berry_connection:'+f.label,T.equal(energy,T.scale(z.BerryA,h)));
  no('lift:wrong_Hamiltonian_sign:'+f.label,!T.equal(T.derivative(z.lift),T.scale(T.mul(z.liftGenerator,T.project(z.lift,depth)),N.IOTA.mul(new N.F(1).div(h)))));
  ck('horizontal_phase:state_norm_preserved:'+f.label,T.equal(T.mul(T.dagger(z.horizontalState),z.horizontalState),T.identity(C,1)));
  ck('horizontal_phase:Berry_connection_removed:'+f.label,T.equal(T.mul(T.dagger(T.project(z.horizontalState,depth)),T.derivative(z.horizontalState)),T.constant(z.BerryA.context,[[0]])));
  const wrongPhase=T.ordinaryExponential(T.scale(T.integrateZero(z.BerryA),N.IOTA.neg())),wrongState=T.mul(PCR.scalarIdentity(wrongPhase,n),z.psi);
  no('horizontal_phase:wrong_scalar_phase_sign:'+f.label,!T.equal(T.mul(T.dagger(T.project(wrongState,depth)),T.derivative(wrongState)),T.constant(z.BerryA.context,[[0]])));
  fixtures.push({kind:'finite_density_and_state_phase_path_not_current_representation',...f,phase:T.serialize(z.V),BerryA:T.serialize(z.BerryA),horizontalNorm:T.serialize(z.horizontalNormExpected),lift:T.serialize(z.lift),Hamiltonian:T.serialize(z.liftGenerator),currentReadoutContract:z.currentContract});
}

no('domain:h_zero',rejects(()=>PCR.quantumGeometry(N.identity(2),0)));
no('domain:h_one',rejects(()=>PCR.quantumGeometry(N.identity(2),1)));
no('domain:singular_covariance',rejects(()=>PCR.quantumGeometry([[1,1],[1,1]],'1/2')));
no('domain:nonpositive_mean_logarithm',rejects(()=>PCR.readoutGeometry([0,1],N.identity(2),[0,0])));
no('domain:complex_phase_parameter',rejects(()=>PCR.readoutGeometry([1,1],N.identity(2),[N.IOTA,0])));

function register(S,name){const audit=S.audit();ck('words:confluent:'+name,audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name,spec:S.spec(),audit});}
function prove(S,name,a,b){const r=S.reduce(a.minus(b),{witness:true});ck('words:'+name,r.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(r.certificate));ck('words:replay:'+name,S.replay(certificate));proofs.push({name,presentation_hash:S.hash(),lhs:a.toJSON(),rhs:b.toJSON(),certificate});}
{
  // Infinite affine-current algebra contract. No finite current matrices are
  // supplied. a,b denote beta/h, and the central h remains symbolic.
  const tokens=['a','b','h','x','y','p','q'],rules=[];
  for(let j=1;j<tokens.length;j++)for(let i=0;i<Math.min(j,3);i++)rules.push({id:'central_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'central real phase and normalization parameters'});
  rules.push({id:'commuting_densities',lhs:['y','x'],rhs:[[['x','y'],1]],source:'commuting positive density observables'});
  rules.push({id:'current_p_x',lhs:['p','x'],rhs:[[['x','p'],1],[['h','x'],N.IOTA]],source:'[P1,X1]=i h X1 on the analytic common core'});
  rules.push({id:'current_p_y',lhs:['p','y'],rhs:[[['y','p'],1]],source:'disjoint current-density windows'});
  rules.push({id:'current_q_x',lhs:['q','x'],rhs:[[['x','q'],1]],source:'disjoint current-density windows'});
  rules.push({id:'current_q_y',lhs:['q','y'],rhs:[[['y','q'],1],[['h','y'],N.IOTA]],source:'[P2,X2]=i h X2 on the analytic common core'});
  const S=new P.Presentation({tokens,rules});register(S,'symbolic_affine_currents_only');
  const [a,b,h,x,y,p,q]=tokens.map(t=>S.word([t])),phase=a.times(x).plus(b.times(y)).scale(N.IOTA);
  prove(S,'affine_current_diagonal_p',P.bracket(p,x),h.times(x).scale(N.IOTA));
  prove(S,'affine_current_diagonal_q',P.bracket(q,y),h.times(y).scale(N.IOTA));
  prove(S,'affine_current_cross_p_y',P.bracket(p,y),S.zero());
  prove(S,'affine_current_cross_q_x',P.bracket(q,x),S.zero());
  const shifts=[h.times(a).times(x),h.times(b).times(y)];
  for(const [i,current] of [p,q].entries()){
    prove(S,'phase_conjugation_first_commutator_'+i,P.bracket(phase,current),shifts[i]);
    prove(S,'phase_conjugation_nested_commutator_zero_'+i,P.bracket(phase,P.bracket(phase,current)),S.zero());
    const degree=3,E=P.formalExp(phase,degree),Ei=P.formalExp(phase.scale(-1),degree),middle=[current,S.zero(),S.zero(),S.zero()];
    const conjugated=P.seriesMultiply(E,P.seriesMultiply(middle,Ei,degree),degree);
    prove(S,'phase_current_shift_coefficient_'+i,conjugated[1],shifts[i]);
    prove(S,'phase_current_shift_no_quadratic_'+i,conjugated[2],S.zero());
    prove(S,'phase_current_shift_no_cubic_'+i,conjugated[3],S.zero());
  }
  no('words:affine_commutator_is_constant_CCR',S.reduce(P.bracket(p,x).minus(h.scale(N.IOTA))).normal.terms.size>0);
  no('words:phase_conjugation_shift_has_opposite_sign',S.reduce(P.bracket(phase,p).plus(shifts[0])).normal.terms.size>0);
}
const payload={protocol:'PCR1_PHASE_CURRENT_RESPONSE_V1',status:'PASS_PCR1_PHASE_CURRENT_RESPONSE',source_parent:'7664681c901f4d063b09e0d021d1d74b0a615ca1',checks,negative_controls,presentations,proofs,fixtures,
  scope:{finite_matrices:'commuting density spectra, normalized states, response frames and phase lifts only',currents:'symbolic affine algebra plus analytically specified mean readout; no finite current or CCR representation',QGT_complex_rank:2,real_parameter_metric_rank:4,Berry_sign:'F_theta_beta=-H/h; omega=-F/2',phase_mean_recovery:'fixed source family with positive means and invertible covariance',phase_lift_generator:'dot(beta) dot X + i h V B V adjoint',finite_fixtures_prove_current_domains:false,finite_fixtures_prove_continuum_support:false,four_parameters_are_spacetime:false,extra_qubit_added:false,quantum_gravity_closed:false},
  summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,negative_controls:negative_controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',hash=crypto.createHash('sha256').update(bytes).digest('hex'),cp=path.join(__dirname,'CERTIFICATE.json'),hp=path.join(__dirname,'EXPECTED.sha256');
if(argv[0]==='--write'){fs.writeFileSync(cp,bytes);fs.writeFileSync(hp,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cp,'utf8'),bytes);assert.equal(fs.readFileSync(hp,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
