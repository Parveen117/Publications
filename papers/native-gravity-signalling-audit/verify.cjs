'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const m=N.matrix,mm=(...xs)=>xs.reduce(N.mul),dg=N.dagger,eq=N.equal,sc=N.scale,kron=N.kron,f=N.F.of;
const checks=[],negative_controls=[],proofs=[],presentations=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const rejects=fn=>{try{fn();return false;}catch{return true;}};
const I=N.identity(2),I4=N.identity(4),P0=m([[1,0],[0,0]]),P1=N.sub(I,P0),Z=N.sub(P0,P1),X=m([[0,1],[1,0]]);
const tr=A=>A.reduce((s,row,j)=>s.add(row[j]),N.ZERO),conj=(U,rho)=>mm(U,rho,dg(U));
const sum=(xs,dim)=>xs.reduce(N.add,N.zeros(dim));
function partialB(rho,da,db=2){
 const out=N.zeros(db);assert.equal(rho.length,da*db);
 for(let i=0;i<db;i++)for(let j=0;j<db;j++)for(let a=0;a<da;a++)out[i][j]=out[i][j].add(rho[a*db+i][a*db+j]);
 return out;
}
const cut=(z,x)=>sc(N.add(I,N.add(sc(Z,z),sc(X,x))),'1/2');
const rhoAB=sc(m([[4,0,0,2],[0,0,0,0],[0,0,0,0],[2,0,0,1]]),'1/5');
ck('steering:source_pure_positive_cut',eq(rhoAB,dg(rhoAB))&&eq(mm(rhoAB,rhoAB),rhoAB)&&tr(rhoAB).eq(1));
const common=partialB(rhoAB,2),Qplus=sc(N.add(I,X),'1/2'),Qminus=N.sub(I,Qplus);
const zBranches=[P0,P1].map(K=>partialB(conj(kron(K,I),rhoAB),2));
const xBranches=[Qplus,Qminus].map(K=>partialB(conj(kron(K,I),rhoAB),2));
const normalize=branches=>branches.map(rho=>{const weight=tr(rho);assert.ok(!weight.zero());return {weight,state:sc(rho,weight.inv())};});
const ez=normalize(zBranches),ex=normalize(xBranches);
ck('steering:Z_weights',ez[0].weight.eq('4/5')&&ez[1].weight.eq('1/5'));
ck('steering:X_weights',ex.every(e=>e.weight.eq('1/2')));
ck('steering:Z_conditional_cuts',eq(ez[0].state,P0)&&eq(ez[1].state,P1));
ck('steering:X_conditional_cuts',eq(ex[0].state,cut('3/5','4/5'))&&eq(ex[1].state,cut('3/5','-4/5')));
ck('steering:unconditional_marginals_equal',eq(sum(zBranches,2),common)&&eq(sum(xBranches,2),common));
const qgc=JSON.parse(fs.readFileSync(path.join(__dirname,'../native-quantum-gravity-consistency/CERTIFICATE.json'),'utf8')).results.ensemble;
ck('source:QGC_common_density',eq(common,m(qgc.common_density)));
for(const [label,entries,source]of [['Z',ez,qgc.ensembles.ea],['X',ex,qgc.ensembles.eb]])entries.forEach((e,j)=>ck('source:QGC_steering_ensemble:'+label+':'+j,e.weight.eq(source[j].weight)&&eq(e.state,m(source[j].state))));
no('postselection:conditional_difference_is_unconditional_signal',!eq(ez[0].state,ex[0].state)&&eq(sum(zBranches,2),sum(xBranches,2)));
no('postselection:unit_probability_outcomes',!ez[0].weight.eq(1)&&!ex[0].weight.eq(1));
// Local null controls on a full operator basis, not only selected positive states.
const filter=[m([[1,0],[0,'3/5']]),m([[0,0],[0,'4/5']])];
const damping=[m([[1,0],[0,'3/5']]),m([[0,'4/5'],[0,0]])];
const u=m([['3/5'],['4/5']]),P=mm(u,dg(u)),Q=N.sub(I,P);
const history=[mm(Q,P1),mm(P,P1),mm(Q,P0),mm(P,P0)];
const instruments=[['Z_measurement',[P0,P1]],['X_measurement',[Qplus,Qminus]],['QGC_filter',filter],['damping',damping],['CHD_history',history]];
for(const [name,Ks]of instruments){
 ck('local:'+name+':complete',eq(sum(Ks.map(K=>mm(dg(K),K)),2),I));
 for(let i=0;i<4;i++)for(let j=0;j<4;j++){
  const e=N.zeros(4);e[i][j]=N.ONE;
  const out=sum(Ks.map(K=>conj(kron(K,I),e)),4);
  ck('local:'+name+':basis_'+i+'_'+j,eq(partialB(out,2),partialB(e,2)));
 }
}
const W=history.flatMap(K=>K),WA=kron(W,I),coherent=conj(WA,rhoAB);
ck('local:coherent_history_isometry',eq(mm(dg(W),W),I));
ck('local:coherent_receiver_unchanged',eq(partialB(coherent,W.length),common));
ck('local:coherent_receiver_effect_intertwiner',eq(mm(dg(WA),kron(N.identity(W.length),P0),WA),kron(I,P0)));
no('local:coherent_retention_alone_sends_a_bit',eq(partialB(coherent,W.length),common));
// Actual inherited nonlinear candidate; do not replace it by a new fitted law.
const g=f('37/470');
const F=P=>{
 const z=tr(mm(Z,P)),x=tr(mm(X,P)),H=sc(N.add(sc(Z,z.div(10)),sc(X,x.div(47))),-1);
 return N.commutator(N.commutator(P,H),P);
};
const avgF=entries=>sum(entries.map(e=>sc(F(e.state),e.weight)),2);
const fz=avgF(ez),fx=avgF(ex),difference=N.sub(fx,fz);
ck('candidate:stationary_Z_ensemble',N.isZero(fz));
ck('candidate:source_velocity_reproduced',eq(fx,m(qgc.average_velocity.eb)));
ck('candidate:Z_contrast_rate',tr(mm(Z,difference)).eq(qgc.Z_readout_velocity_difference)&&tr(mm(Z,difference)).eq('1776/29375'));
ck('candidate:probability_contrast_rate',tr(mm(P0,difference)).eq('888/29375'));
no('candidate:affine_density_only_extension',!eq(fz,fx)&&eq(sum(zBranches,2),sum(xBranches,2)));
const trajectories=[];
for(const [zz,xx,uu]of [['3/5','4/5','1'],['4/5','3/5','81/256'],['12/13','5/13','25/256']]){
 const z=f(zz),x=f(xx),u=f(uu),den=f(9).add(f(16).mul(u));
 ck('trajectory:'+zz+':unit_circle',z.pow(2).add(x.pow(2)).eq(1));
 ck('trajectory:'+zz+':closed_form_z',z.pow(2).mul(den).eq(9));
 ck('trajectory:'+zz+':closed_form_x',x.pow(2).mul(den).eq(f(16).mul(u)));
 const zdot=f(2).mul(g).mul(z).mul(x.pow(2)),xdot=f(-2).mul(g).mul(z.pow(2)).mul(x);
 for(const sign of [1,-1]){
  const state=cut(z,x.mul(sign)),value=F(state);
  ck('trajectory:'+zz+':actual_source_derivative:'+sign,eq(value,sc(N.add(sc(Z,zdot),sc(X,xdot.mul(sign))),'1/2')));
  ck('trajectory:'+zz+':trace_zero:'+sign,tr(value).zero());
 }
 const zFromClosed=f(32).mul(g).mul(u).mul(z).div(den);
 ck('trajectory:'+zz+':derivative_closed_form',zdot.eq(zFromClosed));
 trajectories.push({z,x,u,zdot,xdot});
}
ck('trajectory:exact_count_prefactor',f(4).mul(g).mul('235/74').eq(1));
const endX=sc(N.add(cut('4/5','3/5'),cut('4/5','-3/5')),'1/2'),endZ=common;
ck('candidate:finite_X_density',eq(endX,m([['9/10',0],[0,'1/10']])));
ck('candidate:finite_Z_density',eq(endZ,m([['4/5',0],[0,'1/5']])));
const finiteContrast=tr(mm(P0,N.sub(endX,endZ)));
ck('candidate:finite_probability_contrast',finiteContrast.eq('1/10'));
ck('candidate:linearized_response_contrast',tr(mm(Z,N.sub(endX,endZ))).div(10).eq('1/50'));
// An affine repair applied to either initial ensemble gives exactly the same result.
const channel=rho=>sum(filter.map(K=>conj(K,rho)),2);
const avgChannel=e=>sum(e.map(a=>sc(channel(a.state),a.weight)),2);
ck('candidate:affine_repair_has_zero_contrast',eq(avgChannel(ez),avgChannel(ex))&&eq(avgChannel(ez),channel(common)));
// Positive control: the actual MGH matter/geometry interaction, with no speed assigned.
const mgh=JSON.parse(fs.readFileSync(path.join(__dirname,'../native-matter-geometry-holonomy/CERTIFICATE.json'),'utf8')).results.quantum_model;
const U=m(mgh.joint_operation),V=m([[new N.Cut('9/25','12/25'),new N.Cut('16/25','-12/25')],[new N.Cut('16/25','-12/25'),new N.Cut('9/25','12/25')]]);
ck('joint:published_operation_reproduced',eq(U,N.add(kron(P0,I),kron(P1,V))));
ck('joint:unitary',eq(mm(dg(U),U),I4));
const flip=kron(X,I),effect=kron(I,P1),O=mm(dg(U),effect,U),signal=N.sub(mm(flip,O,flip),O);
const ref=kron(P0,P0),refFlipped=conj(flip,ref);
ck('joint:no_signal_before_interaction',tr(mm(effect,ref)).eq(0)&&tr(mm(effect,refFlipped)).eq(0));
ck('joint:reference_output_zero',tr(mm(effect,conj(U,ref))).eq(0));
ck('joint:reference_output_one',tr(mm(effect,conj(U,refFlipped))).eq('16/25'));
ck('joint:Heisenberg_signal_formula',tr(mm(ref,signal)).eq('16/25'));
ck('joint:exact_signal_norm_square',eq(mm(signal,signal),sc(I4,'16/25')));
const optimal=m([['9/10',new N.Cut(0,'3/10')],[new N.Cut(0,'-3/10'),'1/10']]);
ck('joint:optimal_geometry_positive_cut',eq(optimal,dg(optimal))&&eq(mm(optimal,optimal),optimal)&&tr(optimal).eq(1));
ck('joint:attained_maximum_contrast',tr(mm(kron(P0,optimal),signal)).eq('4/5'));
no('joint:interaction_is_sender_only',!N.isZero(signal));
no('joint:reference_contrast_is_optimal',!f('16/25').eq('4/5'));
// Native word proofs with explicit commutation contracts.
function register(name,s){const audit=s.audit();ck(name+':confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name,spec:s.spec(),audit});return s;}
const systems=new Map();
function prove(s,name,lhs,rhs){const r=s.reduce(lhs.minus(rhs),{witness:true});ck(name,r.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(r.certificate));ck(name+':replay',s.replay(certificate));proofs.push({name,lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate});systems.set(s.hash(),s);}
const local=register('separate_registers',new M.Presentation({tokens:['p','b'],rules:[{id:'pp',lhs:['p','p'],rhs:[[['p'],1]],source:'Sender orthogonal cut'},{id:'bp',lhs:['b','p'],rhs:[[['p','b'],1]],source:'Receiver effect acts on a separate tensor register'}]}));
const p=local.word(['p']),b=local.word(['b']),q=local.one().minus(p);
prove(local,'words:local_cut_receiver_null',p.times(b).times(p).plus(q.times(b).times(q)),b);
const tokens=['g','z','x'],rules=[];
for(let j=1;j<tokens.length;j++)for(let i=0;i<j;i++)rules.push({id:'commute_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'Commuting real branch coordinates and fixed source coefficient'});
rules.push({id:'circle',lhs:['x','x'],rhs:[[[],1],[['z','z'],-1]],source:'Pure real cut unit circle'});
const scalar=register('source_pure_circle',new M.Presentation({tokens,rules})),gg=scalar.word(['g']),zz=scalar.word(['z']),xx=scalar.word(['x']);
const zdot=gg.times(zz).times(xx).times(xx).scale(2),xdot=gg.times(zz).times(zz).times(xx).scale(-2);
prove(scalar,'words:pure_circle_preserved',zz.times(zdot).scale(2).plus(xx.times(xdot).scale(2)),scalar.zero());
prove(scalar,'words:logistic_square_equation',zz.times(zdot).scale(2),gg.times(zz).times(zz).times(scalar.one().minus(zz.times(zz))).scale(4));
const bad=JSON.parse(JSON.stringify(proofs[0].certificate));bad.steps[0].coefficient=['991','0'];
no('integrity:corrupted_symbolic_witness',rejects(()=>systems.get(bad.presentation).replay(bad)));
const hash=x=>crypto.createHash('sha256').update(x).digest('hex');
const payload={protocol:'GSA1_CONDITIONAL_GRAVITY_SIGNALLING_AUDIT_V1',status:'PASS_CONDITIONAL_SIGNALLING_AUDIT_NOT_PHYSICAL_FTL',source_pins_sha256:hash(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),presentations,proofs,checks,negative_controls,results:{source:rhoAB,receiver_common_density:common,steering:{Z:ez,X:ex},candidate:{initial_probability_contrast_rate:'888/29375',trajectories,finite_count:'(235/74)*log(256/81)',finite_X_density:endX,finite_Z_density:endZ,finite_probability_contrast:finiteContrast,linearized_response_contrast:'1/50',required_extra_premise:'componentwise nonlinear physical evolution of remotely prepared branches without delivered labels'},joint_control:{operation:U,signal_operator:signal,reference_contrast:'16/25',maximal_effect_contrast:'4/5',optimal_geometry:optimal},scope:{local_affine_remote_signal:false,nonlinear_steering_extension_has_mathematical_signal:true,nonlinear_law_admitted_as_affine_quantum_dynamics:false,physical_distance_supplied:false,physical_clock_calibrated:false,physical_FTL_demonstrated:false,consciousness_mechanism_derived:false,quantum_gravity_complete:false}},summary:{exact_checks:checks.length,negative_controls:negative_controls.length,native_presentations:presentations.length,symbolic_certificates:proofs.length,serialized_replays:proofs.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',sha=hash(bytes),file=path.join(__dirname,'CERTIFICATE.json');
if(process.argv.includes('--write')){fs.writeFileSync(file,bytes);fs.writeFileSync(path.join(__dirname,'EXPECTED.sha256'),sha+'  CERTIFICATE.json\n');}
else{assert.equal(fs.readFileSync(file,'utf8'),bytes);assert.equal(fs.readFileSync(path.join(__dirname,'EXPECTED.sha256'),'utf8'),sha+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:sha}));
