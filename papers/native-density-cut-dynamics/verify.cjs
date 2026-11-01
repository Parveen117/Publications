'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const DCD=require('./density_cut_adapter.cjs');
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
const comm=(a,b)=>T.sub(T.mul(a,b),T.mul(b,a)),c0=x=>T.coefficient(x,0),tr=M=>M.reduce((s,row,i)=>s.add(row[i]),N.ZERO);
const cases=[
  {name:'three_atom_saturated',amplitudes:['1/3','2/3','2/3'],atoms:[[1,1],[2,1],[1,2]],order:3},
  {name:'four_atom_leakage',amplitudes:['1/2','1/2','1/2','1/2'],atoms:[[1,1],[2,1],[1,2],[3,4]],order:3},
  {name:'five_atom_leakage',amplitudes:['1/7','2/7','2/7','2/7','6/7'],atoms:[[1,1],[2,1],[1,2],[3,4],[4,2]],order:2}
];
for(const f of cases){
  const source=QRT.finiteSource(f.amplitudes,f.atoms,f.order,f.name),z=DCD.densityCutPath(source,{a:1,b:1,p:1,q:2}),C=z.K.context,n=source.atoms.length,d=f.order-1;
  const p=x=>T.project(x,d),Pi=p(z.Pi),Q=p(z.Q),R=p(z.remainder),D=p(z.D),psi=p(z.psi),E=p(z.E),I=T.identity(Pi.context,n),Z=T.constant(Pi.context,N.zeros(n));
  ck('source:exact_density_generator:'+f.name,T.equal(z.B,z.Bprojection));
  ck('source:cut_orthogonality:'+f.name,T.equal(T.mul(Pi,Q),Z)&&T.equal(T.mul(Q,Pi),Z));
  ck('source:three_sector_resolution:'+f.name,T.equal(T.add(T.add(Pi,Q),R),I));
  ck('source:remainder_projection:'+f.name,T.equal(T.mul(R,R),R));
  ck('source:density_maps_state_into_base_plus_response:'+f.name,T.equal(T.mul(R,T.mul(D,psi)),T.constant(Pi.context,N.zeros(n,1))));
  ck('source:base_remainder_density_block_zero:'+f.name,T.equal(T.mul(Pi,T.mul(D,R)),Z));
  ck('generator:anti_adjoint:'+f.name,T.equal(T.dagger(z.B),T.scale(z.B,-1)));
  ck('generator:normalized_state_derivative:'+f.name,T.equal(T.mul(z.B,psi),T.derivative(z.psi)));
  ck('generator:response_derivative:'+f.name,T.equal(T.mul(z.B,E),T.sub(T.derivative(z.E),T.mul(E,z.A))));
  ck('generator:base_state_transport:'+f.name,T.equal(T.mul(z.K,T.constant(C,c0(z.psi))),z.psi));
  ck('generator:equals_Kato_transport:'+f.name,T.equal(z.K,z.Kprojection));
  ck('generator:full_unitarity:'+f.name,T.equal(T.mul(T.dagger(z.K),z.K),T.identity(C,n)));
  ck('generator:response_frame_transport:'+f.name,T.equal(T.mul(z.K,T.constant(C,c0(z.E))),T.mul(z.E,z.U)));
  ck('generator:projection_transport:'+f.name,T.equal(T.mul(z.K,T.mul(T.constant(C,c0(z.Q)),T.dagger(z.K))),z.Q));
  const W=p(z.W);
  ck('cut:Hermitian_weighted_operator:'+f.name,T.equal(W,T.dagger(W)));
  ck('cut:three_exact_weights:'+f.name,T.equal(T.mul(W,T.mul(T.sub(W,T.scale(I,'1/2')),T.sub(W,I))),Z));
  ck('cut:identity_shift_does_not_change_commutator:'+f.name,T.equal(comm(T.add(W,T.scale(I,'3/7')),D),z.B));
  const alpha=N.Cut.of('7/5'),beta=N.Cut.of('2/7'),weighted=DCD.weightedGenerator(z,alpha,beta);
  ck('weights:base_response_block:'+f.name,T.equal(T.mul(Pi,T.mul(weighted,Q)),T.scale(T.mul(Pi,T.mul(D,Q)),alpha.sub(beta))));
  ck('weights:response_remainder_block:'+f.name,T.equal(T.mul(Q,T.mul(weighted,R)),T.scale(T.mul(Q,T.mul(D,R)),beta)));
  ck('weights:remainder_response_block:'+f.name,T.equal(T.mul(R,T.mul(weighted,Q)),T.scale(T.mul(R,T.mul(D,Q)),beta.neg())));
  no('weights:wrong_base_speed:'+f.name,!T.equal(T.mul(DCD.weightedGenerator(z,1,0),psi),T.derivative(z.psi)));
  no('weights:reversed_commutator_sign:'+f.name,!T.equal(T.mul(T.scale(z.B,-1),psi),T.derivative(z.psi)));
  const wrong=DCD.weightedGenerator(z,'1/2',0),wrongK=T.transport(T.scale(wrong,-1));
  ck('weights:state_only_test_cannot_fix_response_weight:'+f.name,T.equal(T.mul(wrong,psi),T.derivative(z.psi)));
  ck('weights:wrong_cut_still_transports_base_state:'+f.name,T.equal(T.mul(wrongK,T.constant(C,c0(z.psi))),z.psi));
  const leakage=T.mul(R,T.mul(D,E)),leakageEnergy=N.energy(c0(leakage));
  if(n===3){
    ck('saturation:remainder_zero',T.equal(R,Z));
    ck('saturation:wrong_weights_hidden',T.equal(wrong,z.B));
    ck('saturation:wrong_weights_transport_hidden',T.equal(wrongK,z.K));
    ck('saturation:normal_leakage_zero',leakageEnergy.zero());
  }else{
    ck('leakage:positive_exact_normal_energy:'+f.name,!leakageEnergy.le(0));
    ck('leakage:remainder_dimension:'+f.name,tr(c0(R)).eq(n-3));
    no('leakage:Pi_over_two_is_full_generator:'+f.name,!T.equal(wrong,z.B));
    no('leakage:state_transport_implies_response_transport:'+f.name,!T.equal(T.mul(wrong,E),T.sub(T.derivative(z.E),T.mul(E,z.A))));
    no('leakage:wrong_cut_preserves_response_intertwiner:'+f.name,!T.equal(T.mul(wrongK,T.constant(C,c0(z.E))),T.mul(z.E,z.U)));
  }
  const W0=c0(z.W),D0=c0(z.D),B0=c0(z.B),pulse=DCD.pulseJet(W0,D0,3),opposite=DCD.pulseJet(W0,D0,3,-1);
  ck('pulse:constant_identity:'+f.name,N.equal(T.coefficient(pulse.product,0),N.identity(n)));
  ck('pulse:first_order_cancels:'+f.name,N.isZero(T.coefficient(pulse.product,1)));
  ck('pulse:ordered_second_coefficient_plus_B:'+f.name,N.equal(T.coefficient(pulse.product,2),B0));
  const cubic=N.scale(N.commutator(N.sub(D0,W0),B0),N.IOTA.mul('1/2'));
  ck('pulse:cubic_commutator_formula:'+f.name,N.equal(T.coefficient(pulse.product,3),cubic));
  ck('pulse:finite_unitarity_coefficients:'+f.name,T.equal(T.mul(T.dagger(pulse.product),pulse.product),T.identity(pulse.context,n)));
  no('pulse:reordered_product_keeps_sign:'+f.name,!N.equal(T.coefficient(opposite.product,2),B0));
  ck('pulse:opposite_order_minus_B:'+f.name,N.equal(T.coefficient(opposite.product,2),N.scale(B0,-1)));
  const ideal=T.ordinaryExponential(T.series(pulse.context,[N.zeros(n),N.zeros(n),B0,N.zeros(n)]));
  no('pulse:four_pulses_are_exact_B_exponential:'+f.name,!N.isZero(T.coefficient(T.sub(pulse.product,ideal),3)));
  ck('pulse:finite_fixture_density_norm_bound:'+f.name,source.atoms.every(atom=>atom[0].rad.le(4)));
  ck('pulse:finite_fixture_generator_bound:'+f.name,N.energy(B0).le(64));
  fixtures.push({kind:'finite_exact_density_cut_comparison',...f,leakageEnergy,point_W:W0,point_D:D0,point_B:B0,
    density_generator:T.serialize(z.B),K:T.serialize(z.K),pulse:T.serialize(pulse.product),cubicDefect:cubic,
    source_contract:'finite spectral model; no continuum domain or nonatomic proof supplied by finite replay'});
}

// Conditional rational resource arithmetic for the reviewer-proved norm bound.
for(const [steps,delta] of [[100,'1/10'],[10000,'1/100'],[1000000,'1/1000']]){
  const bound=DCD.pulseResourceBound({delta,w:1,d:4,steps}),F=N.F.of;
  ck('resource:fixed_effective_duration:'+steps,bound.effectiveDuration.eq(1));
  ck('resource:commutator_bound_two_d:'+steps,bound.b.eq(8));
  ck('resource:exact_cubic_plus_quartic:'+steps,bound.local.eq(F(delta).pow(3).mul(40).add(F(delta).pow(4).mul(32))));
  ck('resource:telescoping_arithmetic:'+steps,bound.total.eq(bound.local.mul(steps)));
  if(steps===1000000){ck('resource:one_million_pulse_blocks_bound',bound.total.eq('1251/31250'));ck('resource:bound_below_one_twentieth',bound.total.le('1/20'));}
  fixtures.push({kind:'conditional_fixed_generator_resource_bound',...bound,analytic_hypotheses:'norm theorem and bounded commutators proved separately; arithmetic does not establish a continuum product limit'});
}
no('resource:negative_step_size',rejects(()=>DCD.pulseResourceBound({delta:'-1/10',d:4})));
no('pulse:nonHermitian_coherent_control',rejects(()=>DCD.pulseJet(N.matrix([[0,1],[0,0]]),N.identity(2))));
no('pulse:unknown_orientation',rejects(()=>DCD.pulseJet(N.identity(2),N.identity(2),3,0)));

function register(S,name){const audit=S.audit();ck('words:confluent:'+name,audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name,spec:S.spec(),audit});}
function prove(S,name,a,b){const r=S.reduce(a.minus(b),{witness:true});ck('words:'+name,r.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(r.certificate));ck('words:replay:'+name,S.replay(certificate));proofs.push({name,presentation_hash:S.hash(),lhs:a.toJSON(),rhs:b.toJSON(),certificate});}
{
  // Three decomposition sectors, not three physical coordinates or a tensor
  // product. Free entries preserve arbitrary internal dimensions of each block.
  const S=new P.Presentation({tokens:['a','b','c','d','e','f','g'],rules:[]});register(S,'source_projection_blocks');
  const [a,b,c,d,e,f,g]=['a','b','c','d','e','f','g'].map(t=>S.word([t])),z=S.zero(),one=S.one();
  const D=[[a,b,z],[c,d,e],[z,f,g]],Q=[[z,z,z],[z,one,z],[z,z,z]],W=[[one,z,z],[z,one.scale('1/2'),z],[z,z,z]];
  const Qprime=[[z,b.scale('1/2'),z],[c.scale('1/2'),z,e.scale('-1/2')],[z,f.scale('-1/2'),z]];
  const mul=(A,B)=>A.map(row=>B[0].map((_,j)=>row.reduce((out,x,k)=>out.plus(x.times(B[k][j])),z)));
  const sub=(A,B)=>A.map((row,i)=>row.map((x,j)=>x.minus(B[i][j]))),br=(A,B)=>sub(mul(A,B),mul(B,A));
  const source=br(Qprime,Q),density=br(W,D);
  for(let i=0;i<3;i++)for(let j=0;j<3;j++)prove(S,'density_cut_block_'+i+j,source[i][j],density[i][j]);
  no('words:Pi_over_two_has_response_remainder_block',source[1][2].terms.size>0);
}
{
  const S=new P.Presentation({tokens:['w','d'],rules:[]});register(S,'ordered_density_cut_pulses');
  const w=S.word(['w']),d=S.word(['d']),C=P.bracket(w,d),i=N.IOTA,degree=3;
  const E=x=>P.formalExp(x,degree),mul=(a,b)=>P.seriesMultiply(a,b,degree);
  const pulse=mul(E(w.scale(i.neg())),mul(E(d.scale(i)),mul(E(w.scale(i)),E(d.scale(i.neg())))));
  prove(S,'pulse_first_coefficient_zero',pulse[1],S.zero());
  prove(S,'pulse_second_coefficient_plus_density_commutator',pulse[2],C);
  prove(S,'pulse_third_coefficient',pulse[3],P.bracket(d.minus(w),C).scale(i.mul('1/2')));
  no('words:pulse_cubic_defect_vanishes',pulse[3].terms.size>0);
}
{
  // Algebraic normalizer obstruction on Sym^3 of an oriented plane. It does
  // not establish spectral support or continuum measurement normalizers.
  const S=new P.Presentation({tokens:['c','s'],rules:[
    {id:'commuting_rotation_coefficients',lhs:['s','c'],rhs:[[['c','s'],1]],source:'central real rotation coefficients'},
    {id:'unit_circle',lhs:['s','s'],rhs:[[[],1],[['c','c'],-1]],source:'c squared plus s squared equals one'}]});
  register(S,'cubic_rotation_normalizer');const c=S.word(['c']),s=S.word(['s']),one=S.one(),z=S.zero();
  const c2=c.times(c),s2=s.times(s),c3=c2.times(c),s3=s2.times(s),c2s=c2.times(s),cs2=c.times(s2);
  const R3=[[c3,c2s,cs2,s3],[c2s.scale(-3),c3.minus(cs2.scale(2)),c2s.scale(2).minus(s3),cs2.scale(3)],
    [cs2.scale(3),s3.minus(c2s.scale(2)),c3.minus(cs2.scale(2)),c2s.scale(3)],[s3.scale(-1),cs2,c2s.scale(-1),c3]];
  const M=R3.map((row,i)=>row.map((x,j)=>(i===j?one:z).minus(x)));
  let determinant=z;
  function permutations(xs,prefix=[]){if(!xs.length)return [prefix];return xs.flatMap((x,i)=>permutations(xs.filter((_,j)=>j!==i),prefix.concat(x)));}
  for(const p of permutations([0,1,2,3])){let inversions=0,term=one;for(let i=0;i<4;i++){term=term.times(M[i][p[i]]);for(let j=i+1;j<4;j++)if(p[i]>p[j])inversions++;}determinant=determinant.plus(term.scale(inversions%2?-1:1));}
  const expected=one.minus(c).times(one.minus(c)).times(c.scale(2).plus(one)).times(c.scale(2).plus(one)).scale(4);
  prove(S,'Sym3_rotation_fixed_tensor_determinant',determinant,expected);
  const cr=new N.F(3,5),sr=new N.F(4,5),cr2=cr.pow(2),sr2=sr.pow(2),cr3=cr.pow(3),sr3=sr.pow(3),cr2sr=cr2.mul(sr),crsr2=cr.mul(sr2);
  const nativeR=N.matrix([[cr3,cr2sr,crsr2,sr3],[cr2sr.mul(-3),cr3.sub(crsr2.mul(2)),cr2sr.mul(2).sub(sr3),crsr2.mul(3)],
    [crsr2.mul(3),sr3.sub(cr2sr.mul(2)),cr3.sub(crsr2.mul(2)),cr2sr.mul(3)],[sr3.neg(),crsr2,cr2sr.neg(),cr3]]);
  ck('normalizer:nonresonant_rational_rotation_has_no_fixed_cubic',N.rank(N.sub(N.identity(4),nativeR))===4);
  ck('normalizer:determinant_rational_marker',new N.F(4).mul(new N.F(1).sub(cr).pow(2)).mul(cr.mul(2).add(1).pow(2)).eq('1936/625'));
}
const payload={protocol:'DCD1_DENSITY_CUT_DYNAMICS_V1',status:'PASS_DCD1_DENSITY_CUT_DYNAMICS',source_parent:'59c68caac773a6c086c16d3c5c565cd7275cb266',checks,negative_controls,presentations,proofs,fixtures,
  scope:{finite_models:'three, four and five positive spectral atoms with exact normalized amplitudes',generator:'B=[Qprime,Q]=[Pi+Q/2,D]',state_and_response_transport:true,higher_response_leakage_tested:true,three_atom_weight_blindness:true,pulses:'exact coefficient identities through cubic order',pulse_norm_bounds:'conditional rational resource arithmetic for separate analytic theorem',normalizer_obstruction:'algebraic Sym3 rotation identity only',finite_checks_prove_continuum_domain_theorem:false,finite_checks_prove_nonatomic_no_go:false,finite_pulses_equal_exact_transport:false,extra_qubit_added:false,spacetime_curvature_established:false,quantum_gravity_closed:false},
  summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,negative_controls:negative_controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',hash=crypto.createHash('sha256').update(bytes).digest('hex'),cp=path.join(__dirname,'CERTIFICATE.json'),hp=path.join(__dirname,'EXPECTED.sha256');
if(argv[0]==='--write'){fs.writeFileSync(cp,bytes);fs.writeFileSync(hp,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cp,'utf8'),bytes);assert.equal(fs.readFileSync(hp,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
