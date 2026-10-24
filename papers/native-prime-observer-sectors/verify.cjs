'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const O=require('./observer_adapter.cjs'),N=require('../native-common-cone-limit/native_engine/native_operator.cjs'),P=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const TET=require('../native-tvsp-exponential-transport/transport_adapter.cjs');
const argv=process.argv.slice(2);assert.ok(argv.length===1&&['--write','--check'].includes(argv[0]));
const checks=[],negative_controls=[],fixtures=[],presentations=[],proofs=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});},no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const rejects=f=>{try{f();return false;}catch{return true;}},same=(a,b)=>a.length===b.length&&a.every((x,i)=>x===b[i]),key=x=>x.join(','),setEqual=(a,b)=>a.length===b.length&&a.every(x=>new Set(b.map(key)).has(key(x)));
const z=O.sourceMatrices(),I=z.I,mm=N.mul,eq=N.equal,scale=N.scale;
function rawAction(G,x){return G.map(row=>row[0].rad.n*x[0]+row[1].rad.n*x[1]);}
function register(S,name){const audit=S.audit();ck('words:confluent:'+name,audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name,spec:S.spec(),audit});}
function prove(S,name,a,b){const r=S.reduce(a.minus(b),{witness:true});ck('words:'+name,r.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(r.certificate));ck('words:replay:'+name,S.replay(certificate));proofs.push({name,presentation_hash:S.hash(),lhs:a.toJSON(),rhs:b.toJSON(),certificate});}

ck('source:rotation_and_reflection_relations',eq(mm(z.R,z.R),scale(I,-1))&&eq(mm(z.K,z.K),I)&&eq(mm(z.R,z.K),scale(mm(z.K,z.R),-1)));
ck('source:nilpotent_seam',N.isZero(mm(z.N,z.N))&&eq(N.add(I,z.N),z.T));
ck('source:native_matrix_units',eq(mm(z.N,z.R),scale(N.add(I,z.K),'1/2'))&&eq(mm(z.R,z.N),scale(N.sub(I,z.K),'1/2'))&&eq(mm(z.R,mm(z.N,z.R)),N.add(z.R,z.N)));
ck('source:adapted_integer_shear',eq(z.U,[[1,1],[0,1]])&&eq(O.adapted(z.R),z.R)&&eq(z.Kadapted,[[1,0],[0,-1]]));
const canonical={e11:[[1,0],[0,0]],e12:[[0,1],[0,0]],e21:[[0,0],[1,0]],e22:[[0,0],[0,1]]};
ck('source:all_four_matrix_units_conjugate',Object.entries(z.matrixUnits).every(([name,M])=>eq(O.adapted(M),canonical[name])));
ck('source:native_basis_index_two',O.det(z.S).eq(2));
ck('lattice:native_S_basis_preserved',O.integral(O.adapted(z.T))&&O.integral(O.adapted(z.R)));
no('lattice:ambient_Z2_is_native_lattice',!O.integral(z.T));
{
  const labels=O.labels(2),ambient=labels.map(x=>O.action(z.S,x,2)),s1=mm(z.S,[[1],[0]]),s2=mm(z.S,[[0],[1]]);
  ck('lattice:intrinsic_mod_two_has_four_labels',labels.length===4);
  no('lattice:ambient_mod_two_faithful',new Set(ambient.map(key)).size===2&&same(ambient[1],ambient[2]));
  ck('lattice:colliding_native_columns_distinct',!eq(s1,s2));
  const nativeEven=mm(z.T,[[2],[0]]);
  no('lattice:ambient_reduction_is_T_equivariant',nativeEven[0][0].eq(1)&&nativeEven[1][0].eq(-1));
  fixtures.push({kind:'intrinsic_versus_ambient_parity',nativeBasis:z.S,intrinsicLabels:labels,ambientLabels:ambient,evenInput:[[2],[0]],nativeTOutput:nativeEven});
}
{
  const bases=O.hnfBases(24),results=bases.map(B=>O.latticeInvariance(B)),invariant=results.filter(x=>x.invariant);
  ck('lattice:all_HNF_through_index_24',results.every(x=>x.invariant===x.scalarLattice));
  ck('lattice:only_scalar_HNF_indices',invariant.map(x=>x.index).sort((a,b)=>Number(a)-Number(b)).join(',')==='1,4,9,16');
  ck('lattice:unimodular_basis_changes_do_not_change_class',O.latticeInvariance([[2,4],[0,2]]).invariant&&O.latticeInvariance([[2,4],[0,2]]).scalarLattice);
  const checker=O.latticeInvariance([[2,1],[0,1]],[z.R,z.Kadapted]);
  ck('lattice:checkerboard_preserves_quartic_generators',checker.invariant);
  no('lattice:quartic_generators_alone_force_scalar',!checker.scalarLattice&&!O.latticeInvariance([[2,1],[0,1]]).invariant);
  fixtures.push({kind:'bounded_column_HNF_classification',maxIndex:24,basesChecked:bases.length,invariant,checkerboard:checker,scope:'finite enumeration supplements the written matrix-unit classification theorem'});
}
{
  const samples=[[0n,0n],[6n,9n],[-12n,18n],[5n,0n]];
  ck('content:positive_gcd_invariant',samples.every(x=>[z.R,z.U].every(G=>O.content(rawAction(G,x))===O.content(x))));
  ck('content:zero_mod_prime_power_iff_divisibility',samples.every(x=>[2n,4n,8n,3n,9n].every(q=>same(O.action(I,x,q),[0n,0n])===(O.content(x)%q===0n))));
  for(const n of [2,3,4,5,6,7,8]){
    const labels=O.labels(n),results=labels.map(x=>{const d=O.gcd(O.content(x),n),closure=O.invariantClosure(n,[x]),expected=labels.filter(y=>y.every(v=>v%d===0n));return {seed:x,contentModulus:d,closureSize:closure.length,valid:setEqual(closure,expected)};});
    ck('module:all_seed_additive_closures:'+n,results.every(x=>x.valid));
    ck('orbit:gcd_classification:'+n,labels.every(x=>setEqual(O.orbit(n,x),labels.filter(y=>O.gcd(O.content(y),n)===O.gcd(O.content(x),n)))));
    if([2,3,5].includes(n))ck('module:prime_is_elementary:'+n,results.every(x=>same(x.seed,[0n,0n])?x.closureSize===1:x.closureSize===n*n));
    if([4,6].includes(n))no('module:composite_is_elementary:'+n,results.some(x=>x.closureSize>1&&x.closureSize<n*n));
    fixtures.push({kind:'finite_module_closure_and_orbit',n,seedsChecked:labels.length,closureSizes:[...new Set(results.map(x=>x.closureSize))],scope:'submodule simplicity of the natural modular module; not complex Hilbert irreducibility'});
  }
}
{
  const points=O.labels(12);
  ck('CRT:all_12_to_4_and_3_roundtrips',points.every(x=>{const a=O.crt(x,4,3);return same(O.crtInverse(a.left,a.right,4,3),x);}));
  ck('CRT:equivariant_for_both_generators',points.every(x=>[z.R,z.U].every(G=>{const a=O.crt(x,4,3),b=O.crt(O.action(G,x,12),4,3);return same(b.left,O.action(G,a.left,4))&&same(b.right,O.action(G,a.right,3));})));
  no('CRT:noncoprime_product_identification',rejects(()=>O.crt([1,1],4,6)));
  fixtures.push({kind:'CRT',modulus:12,factors:[4,3],labelsChecked:points.length});
}
for(const [p,r] of [[2,3],[3,2]]){
  const n=p**r,points=O.labels(n),tag=p+'^'+r;
  ck('tower:digit_roundtrip:'+tag,points.every(x=>same(O.undigits(O.digits(x,p,r),p),x)));
  ck('tower:carry_recurrence_matches_full_action:'+tag,points.every(x=>[z.R,z.U].every(G=>same(O.carryAction(G,O.digits(x,p,r),p).value,O.action(G,x,n)))));
  ck('tower:each_successive_layer_has_p_squared_classes:'+tag,Array.from({length:r},(_,k)=>k).every(k=>points.filter(x=>x.every(v=>v%BigInt(p**k)===0n)).length===p**(2*(r-k))));
  no('tower:additively_split_into_exponent_p_layers:'+tag,BigInt(p)%BigInt(n)!==0n);
  fixtures.push({kind:'prime_power_carry_tower',p,r,modulus:n,labelsChecked:points.length,nonSplitWitness:{element:[1,0],order:n,layerExponent:p}});
}
{
  const input=[[1,1],[0,0]],c=O.carryAction(z.U,input,2),factorwise=input.map(x=>O.action(z.U,x,2));
  ck('carry:explicit_low_digit_and_carry',same(c.digits[0],[0n,1n])&&same(c.carries[1],[1n,0n])&&same(c.digits[1],[1n,0n]));
  no('carry:independent_factorwise_digit_action',!same(O.undigits(factorwise,2),c.value));
  fixtures.push({kind:'binary_shear_carry',input,output:c,factorwise});
}

for(const n of [2,3,4]){
  const PR=O.permutation(z.R,n),PU=O.permutation(z.U,n),P0=O.uniformProjector(n),size=n*n;
  ck('coherent:permutation_unitarity:'+n,eq(mm(N.dagger(PR),PR),N.identity(size))&&eq(mm(N.dagger(PU),PU),N.identity(size)));
  ck('coherent:composition_law:'+n,eq(O.permutation(mm(z.R,z.U),n),mm(PR,PU)));
  ck('coherent:uniform_projector_invariant:'+n,eq(mm(P0,P0),P0)&&eq(N.dagger(P0),P0)&&eq(mm(PR,P0),P0)&&eq(mm(PU,P0),P0));
  no('coherent:complex_lift_is_irreducible:'+n,N.rank(P0)===1&&N.rank(N.sub(N.identity(size),P0))===size-1);
  fixtures.push({kind:'exact_finite_permutation_lift',modulus:n,dimension:size,PR,PU,uniformProjector:P0});
}
{
  const W=O.coarseRecord(2,2,1),edgeKeys=W.edges.map(e=>key(e.coarse)+'|'+key(e.environment));
  ck('coarse_record:sparse_isometry',new Set(edgeKeys).size===W.edges.length&&W.edges.length===16);
  ck('coarse_record:equivariance',W.edges.every(e=>[z.R,z.U].every(G=>{const fine=O.action(G,e.input,4),coarse=O.action(G,e.coarse,2);return same(O.action(I,fine,2),coarse)&&same(fine,O.action(G,e.environment,4));})));
  const reduced=N.zeros(4);for(const e of W.edges){const i=Number(e.coarse[0])*2+Number(e.coarse[1]);reduced[i][i]=reduced[i][i].add('1/16');}
  ck('coarse_record:uniform_fine_state_reduced_probabilities',eq(reduced,scale(N.identity(4),'1/4')));
  no('coarse_record:discard_environment_preserves_pure_coherence',N.rank(reduced)===4);
  fixtures.push({kind:'sparse_coherent_coarse_record',...W,reducedUniformState:reduced});
  // Classicalize the coarse label while retaining the entire source. Equality
  // is checked block by block, so this is distinct from tracing the source.
  const coarseProjectors=O.labels(2).map(a=>N.matrix(W.edges.map((e,i)=>W.edges.map((_,j)=>i===j&&same(e.coarse,a)?1:0))));
  ck('classical_record:coarse_projector_resolution',eq(coarseProjectors.reduce((sum,Q)=>N.add(sum,Q),N.zeros(16)),N.identity(16))&&coarseProjectors.every((Q,i)=>coarseProjectors.every((R,j)=>eq(mm(Q,R),i===j?Q:N.zeros(16)))));
  const plus=N.zeros(16,1),minus=N.zeros(16,1);plus[0][0]=new N.Cut('3/5');plus[1][0]=new N.Cut('4/5');minus[0][0]=new N.Cut('3/5');minus[1][0]=new N.Cut('-4/5');
  const rhoPlus=mm(plus,N.dagger(plus)),rhoMinus=mm(minus,N.dagger(minus)),blocksPlus=coarseProjectors.map(Q=>mm(Q,mm(rhoPlus,Q))),blocksMinus=coarseProjectors.map(Q=>mm(Q,mm(rhoMinus,Q)));
  ck('classical_record:normalized_distinct_inputs',eq(mm(N.dagger(plus),plus),[[1]])&&eq(mm(N.dagger(minus),minus),[[1]])&&!eq(rhoPlus,rhoMinus));
  ck('classical_record:retained_source_outputs_equal',blocksPlus.every((Q,i)=>eq(Q,blocksMinus[i])));
  ck('classical_record:output_probability_retained',blocksPlus.reduce((sum,Q)=>sum.add(Q.reduce((tr,row,i)=>tr.add(row[i]),N.ZERO)),N.ZERO).eq(1));
  no('classical_record:retaining_source_restores_erased_phase',!eq(rhoPlus,rhoMinus)&&blocksPlus.every((Q,i)=>eq(Q,blocksMinus[i])));
  fixtures.push({kind:'classical_coarse_record_with_source_retained',occupiedLabels:[[0,0],[0,1]],amplitudesPlus:['3/5','4/5'],amplitudesMinus:['3/5','-4/5'],rhoPlus,rhoMinus,classicalSourceBlocks:blocksPlus,contract:'the output is the classical direct sum of P_a rho P_a; phase loss persists with the full source retained'});
  const J=O.layerUniform(2,1),uSmall=N.matrix(Array.from({length:4},()=>['1/2'])),uFine=N.matrix(Array.from({length:16},()=>['1/4']));
  ck('uniform_layer:isometry',eq(mm(N.dagger(J.J),J.J),N.identity(4)));
  ck('uniform_layer:both_generators_intertwine',[z.R,z.U].every(G=>eq(mm(O.permutation(G,4),J.J),mm(J.J,O.permutation(G,2)))));
  ck('uniform_layer:uniform_vacua_match',eq(mm(J.J,uSmall),uFine));
  ck('uniform_layer:chosen_Hamiltonians_intertwine',eq(mm(O.hamiltonian(4,'7/3').H,J.J),mm(J.J,O.hamiltonian(2,'7/3').H)));
  fixtures.push({kind:'uniform_prime_layer_isometry',...J});
  const range=mm(J.J,N.dagger(J.J)),detail=N.sub(N.identity(16),range),oldExcitation=N.sub(range,O.uniformProjector(4));
  const fineH=N.add(scale(oldExcitation,'1/2'),scale(detail,'1/4')),coarseH=O.hamiltonian(2,'1/2').H;
  ck('energy_tower:orthogonal_detail_projections',eq(mm(detail,detail),detail)&&eq(mm(oldExcitation,oldExcitation),oldExcitation)&&N.isZero(mm(detail,oldExcitation)));
  ck('energy_tower:nonzero_detail_dimensions',N.rank(detail)===12&&N.rank(oldExcitation)===3);
  ck('energy_tower:same_coarse_intertwining',eq(mm(fineH,J.J),mm(J.J,coarseH)));
  const witness=N.zeros(16,1);witness[0][0]=N.ONE;witness[2][0]=new N.Cut(-1);
  ck('energy_tower:quarter_energy_detail_witness',N.isZero(mm(N.dagger(J.J),witness))&&eq(mm(fineH,witness),scale(witness,'1/4')));
  ck('energy_tower:quarter_lower_bound_on_excited_projection',eq(N.sub(fineH,scale(N.sub(N.identity(16),O.uniformProjector(4)),'1/4')),scale(oldExcitation,'1/4')));
  no('energy_tower:coarse_intertwining_forces_same_gap',!eq(mm(fineH,fineH),scale(fineH,'1/2')));
  fixtures.push({kind:'independent_positive_detail_energy',coarseEnergy:'1/2',newDetailEnergy:'1/4',coarseH,fineH,detailWitness:witness});
}
for(const n of [2,4])for(const kappa of ['0','1/3','7/2']){
  const h=O.hamiltonian(n,kappa),size=n*n,ones=N.matrix(Array.from({length:size},()=>[1])),tag=n+':'+kappa;
  ck('energy:Hermitian_and_projection_polynomial:'+tag,eq(h.H,N.dagger(h.H))&&eq(mm(h.H,h.H),scale(h.H,kappa)));
  ck('energy:uniform_vacuum_zero:'+tag,N.isZero(mm(h.H,ones)));
  ck('energy:rank_from_declared_scale:'+tag,N.rank(h.H)===(kappa==='0'?0:size-1));
  ck('energy:symmetry_invariant:'+tag,[z.R,z.U].every(G=>eq(mm(O.permutation(G,n),h.H),mm(h.H,O.permutation(G,n)))));
  if(kappa==='0')no('energy:positive_excitation_bound_for_zero_model:'+n,N.isZero(h.H));
  fixtures.push({kind:'chosen_energy_scale_not_derived_mass_gap',modulus:n,...h});
}
no('energy:observer_modulus_fixes_energy_scale',!eq(O.hamiltonian(2,'1/3').H,O.hamiltonian(2,'7/2').H));

// TET lives in the adapted diagonal-reflection frame. Only S-conjugation
// identifies it with the NM1 swap-reflection frame; U_TET S is not that seam.
{
  const b=O.transportFrames(1,3),C=b.path.U.context,S=TET.constant(C,z.S),Si=TET.constant(C,z.Sinverse),Ui=TET.inverse(b.path.U),frameInverse=TET.inverse(b.nativeFrame);
  ck('transport:faithful_native_frame_seam',TET.equal(b.nativeFrame,TET.mul(S,b.path.U))&&TET.equal(b.nativeTransport,TET.mul(S,TET.mul(b.path.U,Si))));
  for(const [name,G] of [['R',z.R],['T',z.T]]){
    const moving=TET.mul(b.nativeTransport,TET.mul(TET.constant(C,G),TET.inverse(b.nativeTransport))),inMoving=TET.mul(frameInverse,TET.mul(moving,b.nativeFrame));
    ck('transport:moving_lattice_generator_coordinates:'+name,TET.equal(inMoving,TET.constant(C,O.adapted(G))));
  }
  const loop=b.loop,adaptedCoeff=TET.coefficient(loop.holonomy,2),nativeCoeff=mm(z.S,mm(adaptedCoeff,z.Sinverse)),F0=scale([[5,25],[-10,-5]],'1/324');
  ck('transport:actual_adapted_loop_minus_curvature',eq(loop.curvature,F0)&&eq(adaptedCoeff,scale(F0,-1)));
  ck('transport:native_loop_conjugacy',eq(nativeCoeff,scale(mm(z.S,mm(F0,z.Sinverse)),-1)));
  no('transport:opposite_loop_orientation_sign',!eq(adaptedCoeff,F0));
  no('transport:uncorrected_frame_order',!TET.equal(b.nativeFrame,TET.mul(b.path.U,S)));
  const e=O.rationalEndpoint();
  ck('transport:exact_rational_endpoint_closed_form',new N.F(1).add(new N.F(25,9).mul(e.t)).eq(new N.F(21,20).pow(2))&&e.t.le(new N.F(1,20)));
  ck('transport:endpoint_metric_compatibility',eq(mm(N.dagger(e.U),mm(e.H,e.U)),e.H0));
  ck('transport:endpoint_native_conjugate',eq(e.nativeU,scale([[102,-3],[-2,103]],'1/105')));
  ck('transport:endpoint_determinant',O.det(e.U).eq('20/21')&&O.det(e.nativeU).eq('20/21'));
  no('transport:fixed_adapted_lattice_preserved',!O.integral(e.U));
  no('transport:fixed_native_lattice_preserved',!O.integral(O.adapted(e.nativeU)));
  no('transport:any_fixed_full_rank_lattice_inclusion',O.det(e.U).rad.d!==1n);
  ck('transport:distinct_moving_lattice_is_valid',eq(mm(N.inverse(mm(z.S,e.U)),mm(e.nativeU,z.S)),I));
  fixtures.push({kind:'faithful_TET_frame_bridge',adaptedPath:TET.serialize(b.path.U),nativeFrame:TET.serialize(b.nativeFrame),adaptedCurvature:F0,nativeLoopCoefficient:nativeCoeff,endpoint:e,scope:'fixed-lattice obstruction; moving lattices remain legitimate'});
}

{
  const tokens=['e11','e12','e21','e22'],rules=[];for(const a of tokens)for(const b of tokens){const i=a[1],j=a[2],k=b[1],l=b[2];rules.push({id:a+'_'+b,lhs:[a,b],rhs:j===k?[[['e'+i+l],1]]:[],source:'canonical integer matrix-unit multiplication'});}
  const S=new P.Presentation({tokens,rules});register(S,'canonical_matrix_units');const [e11,e12,e21,e22]=tokens.map(t=>S.word([t])),one=e11.plus(e22),R=e21.minus(e12);
  prove(S,'matrix_unit_identity',one.times(one),one);
  prove(S,'rotation_square',R.times(R),one.scale(-1));
  prove(S,'nilpotent_matrix_unit',e12.times(e12),S.zero());
  prove(S,'NR_plus_RN_identity',e12.times(R).plus(R.times(e12)),one);
  prove(S,'RNR_source_seam',R.times(e12).times(R),R.plus(e12));
  prove(S,'shear_integral_inverse',one.plus(e12).times(one.minus(e12)),one);
}
{
  const S=new P.Presentation({tokens:['a','p'],rules:[{id:'central_scale',lhs:['p','a'],rhs:[[['a','p'],1]],source:'declared nonnegative energy scalar is central'},{id:'projection',lhs:['p','p'],rhs:[[['p'],1]],source:'uniform orthogonal projector'}]});register(S,'chosen_projection_Hamiltonian');const a=S.word(['a']),p=S.word(['p']),one=S.word([]),H=a.times(one.minus(p));
  prove(S,'chosen_Hamiltonian_polynomial',H.times(H),a.times(H));
  prove(S,'vacuum_annihilation',H.times(p),S.zero());
}
no('guard:scientific_or_floating_modulus',rejects(()=>O.action(z.U,[1,1],'2e1'))&&rejects(()=>O.action(z.U,[1,1],2.5)));
no('guard:noninteger_adapted_action',rejects(()=>O.action(z.T,[1,1],3)));
no('guard:singular_lattice',rejects(()=>O.latticeInvariance([[1,1],[1,1]])));
no('guard:rational_lattice_basis',rejects(()=>O.latticeInvariance([['1/2',0],[0,1]])));
no('guard:nonunimodular_lattice_generator',rejects(()=>O.latticeInvariance(I,[[[2,0],[0,2]]])));
no('guard:composite_digit_base',rejects(()=>O.digits([1,1],4,2)));
no('guard:noncanonical_digits',rejects(()=>O.undigits([[2,0]],2)));
no('guard:noninvertible_modular_permutation',rejects(()=>O.permutation([[2,0],[0,1]],4)));
no('guard:negative_energy_scale',rejects(()=>O.hamiltonian(2,-1)));
no('guard:unbounded_dense_model',rejects(()=>O.permutation(z.U,5)));
no('guard:wrong_coarse_depth',rejects(()=>O.coarseRecord(2,2,3)));

const payload={protocol:'POS1_PRIME_OBSERVER_SECTORS_V1',status:'PASS_POS1_PRIME_OBSERVER_SECTORS',source_parent:'7f1a8b9704d3ddc57f7ec079e55bf0efb34f2a0a',source_pins_sha256:crypto.createHash('sha256').update(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))).digest('hex'),checks,negative_controls,presentations,proofs,fixtures,
  scope:{native_frame:'NM1 Lambda=S Z^2; labels always use intrinsic adapted coordinates',lattice:'exact integrality tests and bounded HNF checks supplement the written universal matrix-unit proof',prime:'elementary natural modular module; complex permutation lift is reducible',digits:'canonical integer digits with carries, not an additive direct product of layers',coherent:'finite permutation representations, coarse-record isometry, and uniform layer embeddings',TET:'adapted diagonal frame with native transport S U_TET S^-1',energy:'external nonnegative scale kappa; observer arithmetic does not derive a mass gap',zero_model:'no positive excitation lower bound relative to the declared uniform vacuum',finite_checks_prove_universal_classification:false,prime_modulus_selected_by_dynamics:false,finite_quotient_quantizes_physical_space:false},
  summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,negative_controls:negative_controls.length}};
const bytes=JSON.stringify(payload,(_,value)=>typeof value==='bigint'?value.toString():value,2)+'\n',hash=crypto.createHash('sha256').update(bytes).digest('hex'),cp=path.join(__dirname,'CERTIFICATE.json'),hp=path.join(__dirname,'EXPECTED.sha256');
if(argv[0]==='--write'){fs.writeFileSync(cp,bytes);fs.writeFileSync(hp,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cp,'utf8'),bytes);assert.equal(fs.readFileSync(hp,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
