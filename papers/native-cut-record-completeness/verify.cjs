'use strict';
// CRC-1: use the existing native field, row closure and proof compiler unchanged.
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const root=__dirname,args=process.argv.slice(2),nativeRoot=args[args.indexOf('--native-root')+1],write=args.includes('--write');
assert.ok(args.includes('--native-root')&&nativeRoot);assert.ok(write!==args.includes('--check'));
const hash=x=>crypto.createHash('sha256').update(x).digest('hex');
const pins=JSON.parse(fs.readFileSync(path.join(root,'SOURCE_PINS.json'),'utf8'));
for(const[n,h]of Object.entries(pins.native_engine))assert.equal(hash(fs.readFileSync(path.join(nativeRoot,n))),h,n);
for(const[n,h]of Object.entries(pins.packet_files))assert.equal(hash(fs.readFileSync(path.join(root,n))),h,n);
for(const[n,h]of Object.entries(pins.predecessor_files))assert.equal(hash(fs.readFileSync(path.resolve(root,'../..',n))),h,n);
const N=require(path.join(nativeRoot,'native_operator.cjs')),M=require(path.join(nativeRoot,'paninian_operator.cjs'));
const f=N.F.of,C=N.Cut.of,mm=(...a)=>a.reduce(N.mul),real=z=>{assert.ok(z.turn.zero());return z.rad;};
const col=a=>N.matrix(a.map(z=>[z])),cols=a=>N.matrix(a[0].map((_,i)=>a.map(c=>c[i][0]))),tr=a=>a.reduce((s,row,i)=>s.add(row[i]),N.ZERO);
const checks=[],controls=[],proofs=[],presentations=[],systems=new Map(),results={};
const ck=(name,value,detail={})=>{assert.ok(value,name);checks.push({name,passed:true,...detail});};
const no=(name,value,detail={})=>{assert.ok(value,name);controls.push({name,rejected:true,...detail});};
const positive=x=>!x.le(0),I=N.identity(2),Z=N.matrix([[1,0],[0,-1]]),X=N.matrix([[0,1],[1,0]]),R=N.matrix([[0,-1],[1,0]]),Y=N.scale(R,N.IOTA),B=[I,Z,X,Y];
const cut=(z,x,y=0)=>N.scale(N.add(I,N.add(N.scale(Z,z),N.add(N.scale(X,x),N.scale(Y,y)))),'1/2');
const Q=cut(1,0),Qminus=cut(-1,0),Qx=cut(0,1),Qy=cut(0,0,1);
const rotation=(c,s)=>N.add(N.scale(I,c),N.scale(R,s)),U=rotation('3/5','4/5'),V=N.add(N.scale(I,'3/5'),N.scale(X,N.IOTA.mul('4/5')));
const evolve=(u,rho)=>mm(u,rho,N.dagger(u)),read=(q,u,rho)=>real(tr(mm(q,evolve(u,rho)))),encode=rho=>col(B.map(b=>real(tr(mm(b,rho)))));
const decode=m=>N.scale(B.reduce((sum,b,i)=>N.add(sum,N.scale(b,m[i][0])),N.zeros(2)),'1/2');
const transfer=u=>N.matrix(B.map(a=>B.map(b=>tr(mm(a,u,b,N.dagger(u))).mul('1/2'))));
const Tu=transfer(U),Tv=transfer(V),seed=N.matrix([[1,0,0,0],['1/2','1/2',0,0]]),realRecord=N.matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0]]);
for(let i=0;i<4;i++)for(let j=0;j<4;j++)ck('carrier:trace_pair_'+i+'_'+j,tr(mm(B[i],B[j])).eq(i===j?2:0));
for(const[name,u]of[['U',U],['V',V]])ck('transport:'+name+':unitary',N.equal(mm(N.dagger(u),u),I));
const closureReal=N.rowClosure(seed,[Tu]),closureFull=N.rowClosure(seed,[Tu,Tv]);
ck('closure:real_three_coordinates',closureReal.rank===3&&N.equal(closureReal.observer,realRecord));
ck('closure:extra_phase_four_coordinates',closureFull.rank===4&&N.equal(closureFull.observer,N.identity(4)));
no('closure:initial_two_readings_are_not_closed',N.descendedAction(Tu,seed)===null);
no('closure:real_record_not_closed_under_extra_transport',N.descendedAction(Tv,realRecord)===null);
ck('closure:real_update_descends',N.equal(mm(N.descendedAction(Tu,realRecord),realRecord),mm(realRecord,Tu)));
for(const[name,t]of[['U',Tu],['V',Tv]])ck('closure:full_update_descends_'+name,N.equal(N.descendedAction(t,N.identity(4)),t));

// Retain actual operation words for each independent readout. The arithmetic
// closure is still N.rowClosure; this adapter supplies its finite witnesses.
function wordWitnesses(initial,generators){
 const kept=[];
 for(let i=0;i<initial.length;i++)if(!kept.length||N.rank(kept.map(a=>a.row).concat([initial[i]]))>kept.length)kept.push({row:initial[i],seed:i,word:[]});
 for(let index=0;index<kept.length;index++)for(let g=0;g<generators.length;g++){
  const row=mm([kept[index].row],generators[g])[0];
  if(N.rank(kept.map(a=>a.row).concat([row]))>kept.length)kept.push({row,seed:kept[index].seed,word:[g,...kept[index].word]});
 }
 return kept;
}
const wordsReal=wordWitnesses(seed,[Tu]),wordsFull=wordWitnesses(seed,[Tu,Tv]);
for(const[name,records,generators,closed]of[['real',wordsReal,[Tu],closureReal],['full',wordsFull,[Tu,Tv],closureFull]]){
 ck('closure:'+name+':witness_basis_rank',records.length===closed.rank&&N.rank(records.map(a=>a.row).concat(closed.observer))===closed.rank);
 for(let i=0;i<records.length;i++){
  let action=N.identity(4);for(const j of records[i].word)action=mm(generators[j],action);
  ck('closure:'+name+':witness_word_'+i,N.equal([records[i].row],mm([seed[records[i].seed]],action)));
 }
}
const pPlus=cut('3/5','4/5'),pMinus=cut('3/5','-4/5'),initialPair=[pPlus,pMinus].map(p=>mm(seed,encode(p))),afterPair=[pPlus,pMinus].map(p=>read(Q,U,p));
ck('repair:missing_direction_same_initial_record',N.equal(initialPair[0],initialPair[1]));
ck('repair:missing_direction_distinguished_by_one_rotation',afterPair[0].eq('4/125')&&afterPair[1].eq('4/5'));
no('repair:amount_and_initial_cut_reading_are_sufficient',!afterPair[0].eq(afterPair[1]));
const yPlus=cut(0,0,1),yMinus=cut(0,0,-1),phasePair=[yPlus,yMinus].map(p=>read(Q,V,p));
ck('repair:phase_same_real_record',N.equal(mm(realRecord,encode(yPlus)),mm(realRecord,encode(yMinus))));
ck('repair:phase_new_transport_distinguishes',phasePair[0].eq('1/50')&&phasePair[1].eq('49/50'));
no('repair:real_record_is_universal',!phasePair[0].eq(phasePair[1]));
for(const p of[pPlus,pMinus,yPlus,yMinus])ck('carrier:projector_fixture_'+real(tr(mm(Z,p)))+'_'+real(tr(mm(X,p)))+'_'+real(tr(mm(Y,p))),N.equal(mm(p,p),p)&&real(tr(p)).eq(1));

const sources=[pPlus,pMinus,yPlus,yMinus,N.scale(I,'3/2'),N.scale(pPlus,7)],samples=[];
for(let i=0;i<sources.length;i++){
 const rho=sources[i],m=encode(rho);
 ck('readout:'+i+':exact_reconstruction',N.equal(decode(m),rho));
 for(const[name,u,t]of[['U',U,Tu],['V',V,Tv],['VU',mm(V,U),mm(Tv,Tu)]]){
  ck('readout:'+i+':transfer_'+name,N.equal(encode(evolve(u,rho)),mm(t,m)));
  const out=read(Q,u,rho),pred=real(mm([seed[1]],t,m)[0][0]);
  ck('readout:'+i+':prediction_'+name,out.eq(pred));
  const W=mm(V,U),rw=evolve(W,rho),qw=evolve(W,Q),uw=evolve(W,u);
  ck('readout:'+i+':joint_frame_invariance_'+name,read(qw,uw,rw).eq(out));
 }
 const n=real(m[0][0]),z=real(m[1][0]),x=real(m[2][0]),y=real(m[3][0]);
 const readings=[read(Q,I,rho),read(Qminus,I,rho),read(Qx,I,rho),read(Qy,I,rho)];
 ck('readout:'+i+':cut_reconstruction',readings[0].add(readings[1]).eq(n)&&readings[0].sub(readings[1]).eq(z)&&readings[2].mul(2).sub(n).eq(x)&&readings[3].mul(2).sub(n).eq(y));
 samples.push({source:rho,record:m,cut_readings:readings});
}
no('scale:normalization_discards_source_amount',N.equal(N.scale(sources[5],'1/7'),pPlus)&&!read(Q,I,sources[5]).eq(read(Q,I,pPlus)));
results.observation={basis:B,coordinate_order:['N=Tr(rho)','z=Tr(Z*rho)','x=Tr(X*rho)','y=Tr(Y*rho)'],seed,transports:{U,V},transfer_matrices:{Tu,Tv},real_closure:closureReal,full_closure:closureFull,real_word_witnesses:wordsReal,full_word_witnesses:wordsFull,one_step_direction_witness:{initial:initialPair,after:afterPair},one_step_phase_witness:{after:phasePair},samples,scope:'Minimal dimension among linear real readout records containing the seed and closed under the declared fixed transports. Not a minimum of arbitrary nonlinear encodings or a number of one-shot experimental outcomes.'};

// Reflection products, signed rotation, and the distinction between a state
// readout target and a path/count target.
const J1=Z,J2=N.add(N.scale(Z,'3/5'),N.scale(X,'4/5')),U12=mm(J2,J1),U21=mm(J1,J2),px=Qx;
ck('rotation:two_reflections',N.equal(U12,U)&&N.equal(U21,N.dagger(U)));
ck('rotation:reflection_involutions',N.equal(mm(J1,J1),I)&&N.equal(mm(J2,J2),I));
const orderReadings=[read(Q,U12,px),read(Q,U21,px)];
ck('rotation:order_witness',orderReadings[0].eq('1/50')&&orderReadings[1].eq('49/50'));
no('rotation:count_and_unsigned_trace_determine_next_reading',tr(U12).eq(tr(U21))&&!orderReadings[0].eq(orderReadings[1]));
const historyRecords=[{label:'empty',word:[],count:0},{label:'ZZ',word:[Z,Z],count:2},{label:'XX',word:[X,X],count:2}].map(h=>({...h,endpoint:h.word.reduce((u,j)=>mm(j,u),I)}));
ck('history:distinct_counts_same_endpoint',N.equal(historyRecords[0].endpoint,historyRecords[1].endpoint));
ck('history:same_count_and_endpoint_distinct_intermediates',N.equal(historyRecords[1].endpoint,historyRecords[2].endpoint)&&!N.equal(Z,X));
for(let i=0;i<sources.length;i++)ck('history:identical_endpoint_is_sufficient_for_this_readout_'+i,read(Q,mm(V,historyRecords[0].endpoint),sources[i]).eq(read(Q,mm(V,historyRecords[1].endpoint),sources[i])));
no('history:transport_endpoint_recovers_operation_count',historyRecords[0].count!==historyRecords[1].count&&N.equal(historyRecords[0].endpoint,historyRecords[1].endpoint));
no('history:count_and_endpoint_recover_full_word',historyRecords[1].count===historyRecords[2].count&&N.equal(historyRecords[1].endpoint,historyRecords[2].endpoint)&&!N.equal(Z,X));
ck('history:four_quarter_turns_return',N.equal(N.power(R,4),I));
ck('history:cut_already_returns_after_half_turn',N.equal(evolve(N.power(R,2),pPlus),pPlus));
results.history={reflection_pair:{J1,J2,forward:U12,reverse:U21,readings:orderReadings},records:historyRecords,lift_contract:'For declared continuous real rotations, retain a real angle lift when the target depends on accumulated angle; U alone determines it only modulo 2*pi. A cut orientation is modulo pi. Discrete count alone specifies neither lift nor path.',relevance:'The same endpoint source state gives the same readings after identical memoryless continuations. Count/history is necessary only for declared count/path-dependent targets or dynamics.'};

// Reuse the certified ten-edge tensor. Change only an explicitly declared
// quadratic source-law coefficient, preserving source and first derivative.
const egr2=JSON.parse(fs.readFileSync(path.resolve(root,'../native-nonlinear-gravity-continuation/CERTIFICATE.json'),'utf8'));
const egr1=JSON.parse(fs.readFileSync(path.resolve(root,'../native-exponential-gravity-response/CERTIFICATE.json'),'utf8')).results.geometry,g=egr2.results.geometry;
ck('gravity:predecessor_pass',egr2.status==='PASS_EXACT_NONLINEAR_JET_AND_LOCAL_BRANCH_HYPOTHESES');
const A=N.matrix(g.mode_readout),S=N.matrix(g.section),G=N.matrix(g.gauge),K=N.matrix(g.stiffness),T=g.third_derivative_tensor.map(pl=>pl.map(row=>row.map(C))),H=mm(N.dagger(A),K,A);
const sCols=[0,1].map(j=>S.map(row=>[row[j]])),gCols=Array.from({length:8},(_,j)=>G.map(row=>[row[j]])),B10=cols([...sCols,...gCols]),Gamma=N.inverse(B10).slice(2);
const source=col(['3/5','4/5']),d0=mm(S,N.inverse(K),source),base=egr2.results.response.plus,M0=N.matrix(base.null_hessian_first_coefficient),eta0=N.matrix(base.selected_flat_displacement);
ck('gravity:displacement_coordinate',N.equal(mm(Gamma,G),N.identity(8))&&N.isZero(mm(Gamma,S)));
function contract(a,b){return col(T.map(pl=>pl.reduce((sum,row,j)=>sum.add(row.reduce((s,z,k)=>s.add(z.mul(a[j][0]).mul(b[k][0])),N.ZERO)),N.ZERO)));}
function volumeRate(d){const idx=Object.fromEntries(g.interior_edges.map((e,i)=>[e,i])),dv=(a,b)=>{const key=[a,b].sort().join(':');return idx[key]===undefined?f(0):real(d[idx[key]][0]);};return egr1.cells.map((cell,i)=>{const gd=N.matrix(cell.slice(1).map(a=>cell.slice(1).map(b=>dv(cell[0],a).add(dv(cell[0],b)).sub(a===b?0:dv(a,b)).div(2))));return real(tr(mm(N.inverse(N.matrix(egr1.cell_records[i].gram)),gd))).div(2);});}
const b0=mm(N.dagger(G),contract(d0,d0)),couplings=[];
for(const kappa of[0,2]){
 const k=f(kappa),eta=N.scale(eta0,f(1).div(f(1).sub(k))),d=N.add(d0,mm(G,eta)),F2=N.scale(mm(N.dagger(Gamma),M0,Gamma),k),residual=N.sub(N.scale(contract(d,d),'1/2'),mm(F2,d));
 const z2=N.scale(mm(N.inverse(K),N.dagger(S),residual),-1),e=mm(S,z2),lift=N.sub(mm(N.dagger(G),cols(gCols.map(v=>contract(d,v)))),mm(N.dagger(G),F2,G));
 const expectedLift=N.scale(M0,f(1).sub(k)),block=N.matrix([[lift[0][0],lift[0][4]],[lift[4][0],lift[4][4]]]),det=real(block[0][0].mul(block[1][1]).sub(block[0][1].mul(block[1][0])));
 ck('gravity:kappa_'+kappa+':same_first_order_equations',N.equal(mm(H,d),mm(N.dagger(A),source)));
 ck('gravity:kappa_'+kappa+':same_first_order_curvature_readout',N.equal(mm(A,d),mm(N.inverse(K),source)));
 ck('gravity:kappa_'+kappa+':eight_compatibility_equations',N.isZero(mm(N.dagger(G),residual)));
 ck('gravity:kappa_'+kappa+':all_ten_second_order_equations',N.isZero(N.add(mm(H,e),residual)));
 ck('gravity:kappa_'+kappa+':exact_lift',N.equal(lift,expectedLift));
 ck('gravity:kappa_'+kappa+':implicit_jacobian_nondegenerate',N.rank(lift)===8&&det.eq(48));
 ck('gravity:kappa_'+kappa+':block_sign',positive(real(block[0][0]).mul(kappa===0?-1:1))&&positive(det));
 const rates=volumeRate(d);
 if(kappa===0){ck('gravity:original_response_reproduced',N.equal(d,N.matrix(base.edge_derivative))&&N.equal(z2,N.matrix(base.readout_second_coefficient)));ck('gravity:original_volume_rates_reproduced',rates.every((v,i)=>v.eq(base.volume_log_derivative_at_zero[i])));}
 couplings.push({kappa,source_gradient:mm(N.dagger(A),source),source_hessian:F2,flat_displacement:eta,edge_derivative:d,readout_first:mm(A,d),readout_second:z2,one_second_order_edge_jet:e,null_hessian_first_coefficient:lift,volume_log_derivative:rates,inertia_small_positive_epsilon:kappa===0?{positive:2,negative:8,zero:0}:{positive:10,negative:0,zero:0}});
}
ck('gravity:all_source_first_derivatives_equal',N.equal(couplings[0].source_gradient,couplings[1].source_gradient));
no('gravity:source_and_first_curvature_response_fix_second_response',!N.equal(couplings[0].readout_second,couplings[1].readout_second));
no('gravity:source_and_first_curvature_response_fix_stability',couplings[0].inertia_small_positive_epsilon.negative!==couplings[1].inertia_small_positive_epsilon.negative);
ck('gravity:kappa_two_mixed_volume_signs',couplings[1].volume_log_derivative.every((v,i)=>positive(i===4?v.neg():v)));
const criticalLift=N.scale(M0,0);
no('gravity:kappa_one_has_an_analytic_linear_size_branch',N.isZero(criticalLift)&&!N.isZero(b0),{reason:'Second-order null equation reduces to b0/2=0, contradicted by the exact nonzero b0. This excludes the differentiable O(epsilon) continuation at q0, not all possible nonanalytic or remote solutions.'});
results.gravity={source_cut:pPlus,source_readout:source,source_family:'F_kappa(h)=r^T*A*h+(kappa/2)*(Gamma*h)^T*M(r)*(Gamma*h); h=q-q0',Gamma,old_null_lift:M0,null_obstruction:b0,cases:couplings,new_null_lift:'(1-kappa)*M(r)',additional_record:'The specified nonlinear source map, or its derivatives to the required response order. F second derivative is independent model data, not an extra component secretly implied by rho=N*P.',scope:'Same nine-cell geometry, N, P and source gradient; distinct declared source laws. Written EGR-2 implicit-function proof applies for kappa != 1 and r1*r2 != 0. This is an identifiability witness, not selection of a physical coupling.'};

function system(name,tokens,rules=[],commute=false){if(commute)for(let j=1;j<tokens.length;j++)for(let i=0;i<j;i++)rules.push({id:'central_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'Declared commuting real record coefficients'});const S=new M.Presentation({tokens,rules}),w=Object.fromEntries(tokens.map(t=>[t,S.word([t])])),audit=S.audit();ck(name+':confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name,spec:S.spec(),audit});systems.set(S.hash(),S);return{S,w};}
function prove(name,lhs,rhs){const S=lhs.system,out=S.reduce(lhs.minus(rhs),{witness:true});ck(name,out.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(out.certificate));ck(name+':replay',S.replay(certificate));proofs.push({name,lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate});}
{
 const{S,w:{p}}=system('cut_to_reflection',['p'],[{id:'cut',lhs:['p','p'],rhs:[[['p'],1]],source:'P squared equals P'}]);
 const j=p.scale(2).minus(S.one()),complement=S.one().minus(p);
 prove('cut:reflection_involution',j.times(j),S.one());
 prove('cut:complementary_projector',complement.times(complement),complement);
 prove('cut:complementary_reflection',complement.scale(2).minus(S.one()),j.scale(-1));
}
{
 const{S,w:{z,x}}=system('native_reflection_pair',['z','x'],[{id:'z_square',lhs:['z','z'],rhs:[[[],1]],source:'Z squared equals I'},{id:'x_square',lhs:['x','x'],rhs:[[[],1]],source:'X squared equals I'},{id:'anticommute',lhs:['x','z'],rhs:[[['z','x'],-1]],source:'XZ=-ZX'}]);
 const r=x.times(z),ja=z.scale('3/5').plus(x.scale('4/5')),jb=z.scale('5/13').plus(x.scale('12/13'));
 prove('rotation:quarter_turn_square',r.times(r),S.one().scale(-1));
 prove('rotation:exact_relative_reflection_product',jb.times(ja),S.one().scale('63/65').plus(r.scale('16/65')));
 prove('rotation:reversed_product',ja.times(jb),S.one().scale('63/65').minus(r.scale('16/65')));
 prove('rotation:cancelled_cut_pair',ja.times(ja),S.one());
}
{
 const{S,w:{n,z,x,y,c,s,a,b}}=system('carrier_readout',['n','z','x','y','c','s','a','b'],[],true),one=S.one();
 const c2=c.times(c).minus(s.times(s)),s2=c.times(s).scale(2),zp=c2.times(z).minus(s2.times(x)),xp=s2.times(z).plus(c2.times(x));
 prove('record:rotation_preserves_planar_square_before_unit_constraint',zp.times(zp).plus(xp.times(xp)),c.times(c).plus(s.times(s)).times(c.times(c).plus(s.times(s))).times(z.times(z).plus(x.times(x))));
 const qp=n.plus(z).scale('1/2'),qm=n.minus(z).scale('1/2'),qx=n.plus(x).scale('1/2');
 prove('record:recover_amount',qp.plus(qm),n);
 prove('record:recover_z',qp.minus(qm),z);
 prove('record:recover_x',qx.scale(2).minus(qp).minus(qm),x);
 prove('record:common_scale_invisible_to_ratio_cross_product',s.times(a).times(b).minus(a.times(s.times(b))),S.zero());
 // a,b here are receiver direction coefficients, c,s are transport coefficients.
 const observed=n.plus(a.times(zp)).plus(b.times(xp)).scale('1/2');
 prove('record:relative_reading_formula',observed,n.plus(a.times(c2).plus(b.times(s2)).times(z)).plus(b.times(c2).minus(a.times(s2)).times(x)).scale('1/2'));
}
{
 const{S,w:{k,m,o,p,r,s}}=system('coupling_second_derivative',['k','m','o','p','r','s'],[],true),one=S.one();
 const a=r.scale(-10).minus(s.scale('2/5')),b=s.scale(2),d=s.scale(-10);
 prove('coupling:general_null_block_determinant',a.times(d).minus(b.times(b)),r.times(s).scale(100));
 prove('coupling:lift_changes_by_source_hessian',m.minus(k.times(m)),one.minus(k).times(m));
 const q=o.times(o).times(a).plus(o.times(p).times(b).scale(2)).plus(p.times(p).times(d));
 prove('coupling:quadratic_source_completion',q.minus(k.times(q)),one.minus(k).times(q));
}
const chosen=proofs.find(p=>p.certificate.steps.length),sys=systems.get(chosen.certificate.presentation),bad=JSON.parse(JSON.stringify(chosen.certificate));bad.steps[0].coefficient=['997','0'];
no('integrity:corrupted_proof_rejected',(()=>{try{sys.replay(bad);return false;}catch{return true;}})());
const payload={protocol:'CRC1_NATIVE_CUT_RECORD_COMPLETENESS_V1',status:'PASS_DECLARED_OPERATION_CLOSURE_AND_SOURCE_LAW_WITNESSES',scope:'Exact finite operator closure, cut/readout reconstruction, reflection/history witnesses and second-order gravity source-law comparison. Written proofs cover all words and local analytic continuation under stated contracts. No universal finite information carrier, physical coupling selection or full quantum-gravity theory.',source_pins_sha256:hash(fs.readFileSync(path.join(root,'SOURCE_PINS.json'))),presentations,proofs,checks,negative_controls:controls,results,summary:{native_presentations:presentations.length,symbolic_certificates:proofs.length,serialized_replays:proofs.length,exact_checks:checks.length,negative_controls:controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',h=hash(bytes),cert=path.join(root,'CERTIFICATE.json'),expected=path.join(root,'EXPECTED.sha256');
if(write){fs.writeFileSync(cert,bytes);fs.writeFileSync(expected,h+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cert,'utf8'),bytes);assert.equal(fs.readFileSync(expected,'utf8'),h+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,mode:write?'generation':'read-only replay',...payload.summary,sha256:h,closure_ranks:{real:closureReal.ranks,extra_transport:closureFull.ranks},gravity:couplings.map(c=>({kappa:c.kappa,inertia:c.inertia_small_positive_epsilon,volume_rates:c.volume_log_derivative}))},null,2));
