'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const checks=[],negative_controls=[],proofs=[],presentations=[],systems=new Map();
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const rejects=fn=>{try{fn();return false;}catch{return true;}};
const f=N.F.of,mm=(...xs)=>xs.reduce(N.mul),eq=N.equal,sc=N.scale,dg=N.dagger;
const I=N.identity(2),G=N.matrix([[0,1],[0,0]]),C=N.matrix([[1,-1]]),x=N.matrix([[1],[1]]);
const Hp=N.matrix([[0,0],[0,1]]),Hm=sc(Hp,-1),e1=N.matrix([[1],[0]]),e2=N.matrix([[0],[1]]);
ck('exponential:nilpotent_generator',N.isZero(mm(G,G)));
for(const t of ['-2','-1/2','0','1/3','3']){
 const U=N.add(I,sc(G,t)),Ui=N.sub(I,sc(G,t)),Ct=mm(C,Ui),xt=mm(U,x),Ht=mm(dg(Ui),Hm,Ui);
 ck('exponential:inverse:'+t,eq(mm(U,Ui),I));
 ck('cut:transported_reading:'+t,eq(mm(Ct,xt),mm(C,x)));
 ck('cut:kernel_vector_moves:'+t,N.isZero(mm(Ct,xt)));
 ck('cut:covariant_derivative:'+t,eq(sc(mm(C,G),-1),sc(mm(Ct,G),-1)));
 ck('zero:numerical_vector_fixed:'+t,N.isZero(mm(U,N.zeros(2,1))));
 const r=N.matrix([[2],[-1]]),b=N.matrix([[t],[1]]),rt=N.add(mm(U,r),b),y=N.add(mm(U,N.add(r,x)),b);
 ck('affine:comparison_zero:'+t,N.isZero(mm(C,Ui,N.sub(y,rt))));
 const yp=mm(U,e1),yn=mm(U,e2),Hpt=mm(dg(Ui),Hp,Ui);
 ck('sign:same_positive_kernel:'+t,N.isZero(mm(Hpt,yp)));
 ck('sign:same_negative_kernel:'+t,N.isZero(mm(Ht,yp)));
 ck('sign:negative_witness_retained:'+t,mm(dg(yn),Ht,yn)[0][0].eq(-1));
 ck('sign:positive_witness_retained:'+t,mm(dg(yn),Hpt,yn)[0][0].eq(1));
 const Hdot=sc(N.add(mm(dg(G),Hm,Ui),mm(dg(Ui),Hm,G)),-1);
 ck('frame:covariant_form_derivative:'+t,N.isZero(N.add(Hdot,N.add(mm(dg(G),Ht),mm(Ht,G)))));
 if(t!=='0')no('cut:fixed_observer_keeps_balance:'+t,!N.isZero(mm(C,xt)));
}
const U=N.matrix([[1,2],[0,1]]),Ui=N.inverse(U),D=N.matrix([[1,0],[0,-1]]),Dt=mm(U,D,Ui),metric=mm(dg(Ui),Ui),form=mm(dg(Ui),D,Ui);
ck('metric:operator_form_relation',eq(mm(metric,Dt),form));
ck('metric:self_adjoint_in_retained_metric',eq(mm(dg(Dt),metric),mm(metric,Dt)));
no('metric:nonunitary_similarity_is_Euclidean_self_adjoint',!eq(Dt,dg(Dt)));
no('zero:kernel_determines_sign',!eq(Hp,Hm));
const P=N.matrix([[1,0],[0,0]]);
ck('loss:projection_hides_negative_direction',N.isZero(mm(P,Hm,P)));
no('loss:projection_is_faithful_transport',rejects(()=>N.inverse(P)));
const reciprocal=[];
for(const t of ['9/16','16/9']){
 const z=f(t),delta=f(1).sub(z).div(f(1).add(z)),K=z.div(f(1).add(z).pow(2));
 ck('RH:diagonal_curvature_identity:'+t,delta.pow(2).add(K.mul(4)).eq(1));
 reciprocal.push({gain:t,diagonal:delta.toString(),curvature:K.toString()});
}
ck('RH:reciprocal_curvature_agreement',reciprocal[0].curvature===reciprocal[1].curvature);
no('RH:curvature_determines_signed_diagonal',reciprocal[0].diagonal!==reciprocal[1].diagonal);
const A=N.matrix([[2,1],[1,3]]),b=N.matrix([[1],[2]]),Ai=N.inverse(A);
const schur=matrix=>{const old=matrix.slice(0,2).map(r=>r.slice(0,2)),cross=matrix.slice(0,2).map(r=>[r[2]]);return matrix[2][2].sub(mm(dg(cross),N.inverse(old),cross)[0][0]);};
const innovation_results=[];
for(const d of ['2','1','7/5'])for(const a of ['2','-3','1/2']){
 const matrix=N.matrix([[2,1,1],[1,3,2],[1,2,d]]),W=N.matrix([[1,2,3],[0,2,-1],[0,0,a]]),newM=mm(dg(W),matrix,W),s=schur(matrix),sp=schur(newM);
 ck('innovation:old_block_invertible:'+d+':'+a,N.rank(newM.slice(0,2).map(r=>r.slice(0,2)))===2);
 ck('innovation:margin_scale:'+d+':'+a,sp.eq(s.mul(f(a).pow(2))));
 ck('innovation:zero_preserved:'+d+':'+a,s.zero()===sp.zero());
 ck('innovation:sign_preserved:'+d+':'+a,f(0).le(s.rad)===f(0).le(sp.rad));
 innovation_results.push({d,a,margin:s.toJSON(),transported_margin:sp.toJSON()});
}
// Positive old diagonal and positive new diagonal alone do not suffice.
no('innovation:positive_diagonal_implies_positive_extension',!f(0).le(schur(N.matrix([[2,1,1],[1,3,2],[1,2,1]])).rad));
const nullComplement=N.matrix([[2,1,1],[1,3,2],[1,2,'7/5']]),v=N.matrix([['-1/5'],['-3/5'],[1]]);
ck('innovation:threshold_null_witness',N.isZero(mm(nullComplement,v)));
ck('innovation:completed_square_center',eq(mm(A,Ai,b),b));
// MZC-2: circular operators retain oriented sectors, not only 0/1 membership.
const cayley=(H,a=1)=>mm(N.add(sc(I,a),sc(H,N.IOTA)),N.inverse(N.sub(sc(I,a),sc(H,N.IOTA))));
const uncayley=(V,a=1)=>sc(mm(N.sub(V,I),N.inverse(N.add(V,I))),N.IOTA.neg().mul(a));
const imag=V=>sc(N.sub(V,dg(V)),N.IOTA.neg().div(2));
const circular_cases=[{name:'positive',H:N.matrix([[2,0],[0,1]])},{name:'mixed',H:N.matrix([[1,2],[2,1]])},{name:'zero_and_positive',H:Hp}];
for(const {name,H}of circular_cases)for(const a of ['1','3/2']){
 const V=cayley(H,a),R=N.inverse(N.sub(sc(I,a),sc(H,N.IOTA)));
 ck('circle:unitary:'+name+':'+a,eq(mm(dg(V),V),I));
 ck('circle:inverse:'+name+':'+a,eq(uncayley(V,a),H));
 ck('circle:oriented_form:'+name+':'+a,eq(imag(V),sc(mm(dg(R),H,R),f(a).mul(2))));
}
ck('circle:marked_zero',N.isZero(mm(N.sub(cayley(Hp),I),e1)));
const Vm=cayley(D),Pminus=N.matrix([[0,0],[0,1]]),rotation=N.matrix([['3/5','-4/5'],['4/5','3/5']]);
const Pmoved=mm(rotation,Pminus,dg(rotation)),Vmoved=mm(rotation,Vm,dg(rotation));
ck('sector:unitary_frame',eq(mm(dg(rotation),rotation),I));
ck('sector:projector_transported',eq(mm(Pmoved,Pmoved),Pmoved));
ck('sector:reducing_sector',eq(mm(Pmoved,Vmoved),mm(Vmoved,Pmoved)));
ck('sector:lower_phase_retained',eq(mm(Vmoved,Pmoved),sc(Pmoved,N.IOTA.neg())));
ck('sector:lower_rank_retained',N.rank(Pmoved)===1);
const cross=N.matrix([[1,2],[2,1]]),negative=N.matrix([[1],[-1]]);
ck('sector:first_diagonal_positive',mm(dg(e1),cross,e1)[0][0].eq(1));
ck('sector:second_diagonal_positive',mm(dg(e2),cross,e2)[0][0].eq(1));
no('sector:positive_compressions_imply_positive_full_form',mm(dg(negative),cross,negative)[0][0].eq(-2));
no('sector:cut_sectors_automatically_reduce_source',!N.isZero(N.commutator(Pminus,cross)));
const Upositive=cayley(sc(I,2)),Uproduct=mm(Upositive,Upositive),Usum=cayley(sc(I,4));
ck('circle:positive_initial_phase',eq(imag(Upositive),sc(I,'4/5')));
ck('circle:product_remains_unitary',eq(mm(dg(Uproduct),Uproduct),I));
ck('circle:product_negative_phase',eq(imag(Uproduct),sc(I,'-24/25')));
ck('circle:product_source_reading',eq(uncayley(Uproduct),sc(I,'-4/3')));
ck('circle:source_sum_positive_phase',eq(imag(Usum),sc(I,'8/17')));
no('circle:positive_arc_closed_under_composition',!eq(Uproduct,Usum));
no('circle:chart_boundary_is_finite_zero',rejects(()=>uncayley(sc(I,-1))));
// MZC-3: combine the declared sector defect and cut, then reuse CRC row closure.
const I3=N.identity(3),Ps=N.matrix([[1,0,0],[0,1,0],[0,0,0]]),Qs=N.sub(I3,Ps);
const classifier=sc(N.sub(Ps,Qs),N.IOTA),c3=N.matrix([[0,0,1]]);
const g12=N.matrix([[0,-1,0],[1,0,0],[0,0,0]]),g23=N.matrix([[0,0,0],[0,0,-1],[0,1,0]]);
const col1=N.matrix([[1],[0],[0]]),col2=N.matrix([[0],[1],[0]]);
ck('common:classifier_unitary',eq(mm(dg(classifier),classifier),I3));
ck('common:allowed_upper_arc',eq(mm(classifier,Ps),sc(Ps,N.IOTA)));
ck('common:excluded_lower_arc',eq(mm(classifier,Qs),sc(Qs,N.IOTA.neg())));
const common_cases=[
 {name:'internal_circle',C:c3,G:g12,initial:2,persistent:2,x:col1,basis:N.matrix([[1,0],[0,1],[0,0]])},
 {name:'one_leaking_direction',C:c3,G:g23,initial:2,persistent:1,x:col2,basis:col1},
 {name:'coupled_planes',C:c3,G:N.add(g12,g23),initial:2,persistent:0,x:col1},
 {name:'cut_changes_inside_allowed_arc',C:N.matrix([[1,0,0]]),G:g12,initial:1,persistent:0,x:col2}
];
const common_equilibrium_results=[];
for(const fixture of common_cases){
 const {name,C,G,x}=fixture,defect=C.concat(Qs),closure=N.rowClosure(defect,[G]),O=closure.observer;
 const B=N.descendedAction(G,O),W=N.cayley(G,'1/3'),wx=mm(W,x),Ct=mm(C,dg(W)),Pt=mm(W,Ps,dg(W));
 ck('common:skew_generator:'+name,N.isZero(N.add(G,dg(G))));
 ck('common:initial_dimension:'+name,3-N.rank(defect)===fixture.initial);
 ck('common:persistent_dimension:'+name,3-closure.rank===fixture.persistent);
 ck('common:descended_record_exists:'+name,B!==null);
 ck('common:descended_identity:'+name,eq(mm(O,G),mm(B,O)));
 ck('common:initial_cut_and_sector:'+name,N.isZero(mm(defect,x)));
 ck('common:rational_step_unitary:'+name,eq(mm(dg(W),W),I3));
 ck('common:moving_cut_zero:'+name,N.isZero(mm(Ct,wx)));
 ck('common:moving_sector_admissible:'+name,eq(mm(Pt,wx),wx));
 const powers=defect.concat(mm(defect,G),mm(defect,G,G));
 ck('common:Krylov_span_matches_native_closure:'+name,N.rank(powers)===closure.rank&&N.rank(powers.concat(O))===closure.rank);
 if(fixture.basis){
  const E=fixture.basis,Pe=mm(E,dg(E));
  ck('common:basis_dimension:'+name,N.rank(E)===fixture.persistent);
  ck('common:complete_kernel_basis:'+name,N.isZero(mm(O,E)));
  ck('common:generator_preserves_carrier:'+name,eq(mm(Pe,G,E),mm(G,E)));
  ck('common:finite_step_preserves_carrier:'+name,eq(mm(Pe,W,E),mm(W,E)));
 }
 if(name==='internal_circle'){
  ck('common:fixed_frame_zero_with_internal_motion',N.isZero(mm(defect,wx))&&!eq(wx,x));
  ck('common:entire_initial_carrier_persists',N.descendedAction(G,N.rref(defect).basis)!==null);
 }else{
  no('common:initial_equilibrium_implies_persistence:'+name,!N.isZero(mm(defect,wx)));
  no('common:entire_initial_carrier_is_invariant:'+name,N.descendedAction(G,N.rref(defect).basis)===null);
 }
 if(name==='coupled_planes'){
  ck('common:first_jet_zero',N.isZero(mm(defect,G,x)));
  no('common:first_jet_suffices',!N.isZero(mm(defect,G,G,x)));
 }
 if(name==='cut_changes_inside_allowed_arc'){
  ck('common:allowed_sector_preserved_but_cut_changes',N.isZero(N.commutator(G,Ps))&&eq(mm(Ps,wx),wx));
  no('common:sector_preservation_implies_cut_equilibrium',!N.isZero(mm(C,wx)));
 }
 common_equilibrium_results.push({name,initial_dimension:fixture.initial,persistent_dimension:fixture.persistent,closure_ranks:closure.ranks,defect_record:defect,closed_record:O,descended_generator:B,kernel_basis:fixture.basis||null,finite_step:W,fixed_defect_after_step:mm(defect,wx)});
}
function register(name,spec){const s=new M.Presentation(spec),audit=s.audit();ck(name+':confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name,spec:s.spec(),audit});systems.set(s.hash(),s);return s;}
function prove(name,lhs,rhs){const s=lhs.system,r=s.reduce(lhs.minus(rhs),{witness:true});ck(name,r.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(r.certificate));ck(name+':replay',s.replay(certificate));proofs.push({name,lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate});}
const inv=register('invertible_transport',{tokens:['c','u','v','x'],rules:[{id:'uv',lhs:['u','v'],rhs:[[[],1]],source:'U inverse V'},{id:'vu',lhs:['v','u'],rhs:[[[],1]],source:'V inverse U'}]});
prove('words:transported_cut_reading',inv.word(['c','v','u','x']),inv.word(['c','x']));
prove('words:inverse_kernel_pullback',inv.word(['u','v','x']),inv.word(['x']));
const tokens=['a','b','d','p','q','s'],rules=[];
for(let j=1;j<tokens.length;j++)for(let i=0;i<j;i++)rules.push({id:'commute_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'Real central scalar specialization of Schur block identity'});
const ring=register('real_Schur_coefficients',{tokens,rules});
const [a,bb,d,p,q,s]=tokens.map(t=>ring.word([t])),sq=z=>z.times(z);
const Ap=a.times(sq(p)),bp=p.times(a.times(q).plus(bb.times(s))),dp=a.times(sq(q)).plus(bb.times(q).times(s).scale(2)).plus(d.times(sq(s)));
prove('words:Schur_denominator_cleared',Ap.times(dp).minus(sq(bp)),sq(p).times(sq(s)).times(a.times(d).minus(sq(bb))));
const chosen=proofs[0],bad=JSON.parse(JSON.stringify(chosen.certificate));bad.steps[0].coefficient=['991','0'];
no('integrity:tampered_proof_rejected',rejects(()=>systems.get(chosen.certificate.presentation).replay(bad)));
const hash=x=>crypto.createHash('sha256').update(x).digest('hex');
const payload={protocol:'MZC1_MOVING_ZERO_COVARIANCE_V1',status:'PASS_NATIVE_MOVING_ZERO_COVARIANCE',source_pins_sha256:hash(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),presentations,proofs,checks,negative_controls,results:{reciprocal,innovation_results,common_equilibrium_results,scope:{examples:'exact synthetic controls; no new zeta data',RH_proved:false,QG_complete:false,analytic_RH_bridge_certified:false,general_proofs:['THEOREM.md','CIRCULAR_SECTORS.md','SECTOR_EQUILIBRIUM.md'],source_RH_certificate_rerun:false,independent_RH_generator_identified:false}},summary:{exact_checks:checks.length,negative_controls:negative_controls.length,native_presentations:presentations.length,symbolic_certificates:proofs.length,serialized_replays:proofs.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',sha=hash(bytes),file=path.join(__dirname,'CERTIFICATE.json');
if(process.argv.includes('--write')){fs.writeFileSync(file,bytes);fs.writeFileSync(path.join(__dirname,'EXPECTED.sha256'),sha+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(file,'utf8'),bytes);assert.equal(fs.readFileSync(path.join(__dirname,'EXPECTED.sha256'),'utf8'),sha+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:sha}));
