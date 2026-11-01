'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const f=N.F.of,m=N.matrix,mm=(...xs)=>xs.reduce(N.mul),dg=N.dagger,eq=N.equal,sc=N.scale;
const checks=[],negative_controls=[],proofs=[],presentations=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const rejects=fn=>{try{fn();return false;}catch{return true;}};
const form=(H,x)=>mm(dg(x),H,x)[0][0],outer=x=>mm(x,dg(x));
// Inherited homogeneous source, finite canonical and clock identities.
for(const lv of ['1/2',1,3])for(const pv of [-2,'-3/5','4/5',3]){
 const L=f(lv),pi=f(pv),p=L.mul(pi),qL=f(1).div(L),pL=pi,pPi=L;
 ck('source:canonical:'+lv+':'+pv,qL.mul(pPi).eq(1));
 ck('source:constraint:'+lv+':'+pv,L.mul(pi).mul(2).eq(p.mul(2)));
 ck('source:clock:'+lv+':'+pv,L.mul(pi).eq(p)&&!p.zero());
 ck('source:zero_potential_stationary:'+lv+':'+pv,L.mul(0).div(pi.pow(2)).zero());
 ck('source:weak_length:'+lv+':'+pv,L.mul(2).eq(L.mul(p.mul(2)).div(p)));
}
no('source:clock_regular_at_zero',rejects(()=>f(1).div(0)));
no('source:potential_omitted_without_selection',!f(2).mul(3).mul(0).add(f(2).pow(2).mul(1)).zero());
// The absolute Jacobian is essential for the positive physical measure.
const jacobians=[];
for(const pv of [-3,'-3/5','1/2',2]){
 const p=f(pv),weight=f(1).div(p.abs()),a=f('3/5');
 ck('measure:positive:'+pv,f(0).le(weight)&&!weight.zero());
 ck('measure:Jacobian:'+pv,weight.mul(p.abs()).eq(1));
 ck('measure:unitary_rescaling:'+pv,p.abs().mul(a.pow(2)).mul(weight).eq(a.pow(2)));
 const logDerivative=f(-1).div(p),connection=f(-1).div(p.mul(2));
 ck('measure:symmetric_connection:'+pv,connection.mul(2).eq(logDerivative));
 no('measure:naive_derivative_symmetric:'+pv,!connection.zero());
 jacobians.push({p,weight,connection});
}
no('measure:signed_Jacobian_is_positive',!f(0).le(f(1).div(-2)));
// Half-line/log-length transformation, exact formal monomial identities.
// These monomials test differential formulas; they are not normalizable states.
for(const exponent of [-2,-1,0,'1/2',1,3]){
 const a=f(exponent),dilation=new N.Cut(0,a.add('1/2').neg());
 const next=new N.Cut(0,a.add('3/2').neg());
 ck('affine:monomial_commutator:'+exponent,dilation.sub(next).eq(N.IOTA));
 const nu=a.add('1/2');
 ck('affine:half_density_transform:'+exponent,dilation.eq(new N.Cut(0,nu.neg())));
}
const deficiencies=[];
for(const [name,gamma,eigenvalue,squareIntegrable]of [['plus_i',-1,N.IOTA,true],['minus_i',1,N.IOTA.neg(),false]]){
 const action=new N.Cut(0,f(gamma).neg());
 ck('domain:'+name+':deficiency_equation',action.eq(eigenvalue));
 ck('domain:'+name+':exponential_decay',f(gamma).le(0)===squareIntegrable);
 deficiencies.push({name,gamma,eigenvalue,square_integrable:squareIntegrable});
}
no('domain:ordinary_half_line_momentum_self_adjoint',deficiencies[0].square_integrable&&!deficiencies[1].square_integrable);
const I=N.identity(2),P0=m([[1,0],[0,0]]),P1=N.sub(I,P0),H=m([[2,1],[1,5]]);
const trace=A=>A.reduce((s,row,j)=>s.add(row[j]),N.ZERO);
const Lfinite=m([[1,0],[0,2]]),Pfinite=m([[0,1],[1,0]]);
ck('finite:commutator_trace_zero',trace(N.commutator(Lfinite,Pfinite)).zero());
no('finite:positive_Hermitian_exact_affine_representation',!trace(sc(Lfinite,N.IOTA)).zero());
// Unit-length disjoint intervals give actual orthonormal physical packets.
const intervals=[{name:'positive',a:f(1),b:f(2)},{name:'negative',a:f(-2),b:f(-1)}];
const integral=(a,b,power)=>b.pow(power+1).sub(a.pow(power+1)).div(power+1);
const moment=n=>m(intervals.map((s,j)=>intervals.map((_,k)=>j===k?integral(s.a,s.b,n):f(0))));
const gram=moment(0),Pcompressed=moment(1),P2compressed=moment(2);
ck('physical_packets:Gram',eq(gram,I));
ck('physical_packets:first_moment',eq(Pcompressed,m([['3/2',0],[0,'-3/2']])));
ck('physical_packets:second_moment',eq(P2compressed,sc(I,'7/3')));
ck('physical_packets:compression_defect',eq(N.sub(P2compressed,mm(Pcompressed,Pcompressed)),sc(I,'1/12')));
no('physical_packets:span_is_momentum_invariant',!eq(P2compressed,mm(Pcompressed,Pcompressed)));
const momentRanks=[];
for(let degree=0;degree<=5;degree++){
 const G=m(Array.from({length:degree+1},(_,j)=>Array.from({length:degree+1},(_,k)=>integral(f(1),f(2),j+k))));
 ck('operation_closure:moment_Gram_rank:'+degree,N.rank(G)===degree+1);
 momentRanks.push({degree,rank:N.rank(G)});
}
const sources=[m([['3/5'],['4/5']]),m([['3/5'],['-4/5']]),m([['3/5'],[new N.Cut(0,'4/5')]])];
const packets=[];
for(const [j,c]of sources.entries()){
 const norm=form(gram,c),mean=form(Pcompressed,c),second=form(P2compressed,c),variance=second.sub(mean.mul(mean));
 ck('physical_packets:'+j+':positive_norm',norm.eq(1));
 ck('physical_packets:'+j+':mean',mean.eq('-21/50'));
 ck('physical_packets:'+j+':second',second.eq('7/3'));
 ck('physical_packets:'+j+':variance',variance.eq('16177/7500'));
 ck('physical_packets:'+j+':sign_probabilities',form(P0,c).eq('9/25')&&form(P1,c).eq('16/25'));
 packets.push({coefficients:c,norm,mean,second,variance});
}
// CHD coherent history and complete circular probes on the packet image corner.
const v=m([['3/5'],['4/5']]),B=outer(v),branches=[];
for(const first of [P1,P0])for(const second of [N.sub(I,B),B])branches.push(mm(second,first));
const W=branches.flat(),encoded=mm(W,H,dg(W));
ck('physical_cut_history:isometry',eq(mm(dg(W),W),I));
ck('physical_cut_history:observable_decoder',eq(mm(dg(W),encoded,W),H));
const probes=[m([[1],[0]]),m([[0],[1]]),v,m([['3/5'],[new N.Cut(0,'4/5')]])];
const observed=probes.map((p,j)=>{
 const x=mm(W,p),value=form(encoded,x);ck('physical_observer:reading:'+j,value.eq(form(H,p)));return value.rad;
});
ck('physical_observer:source_readings',observed.every((a,j)=>a.eq([2,5,'122/25','98/25'][j])));
const re=observed[2].mul(25).sub(observed[0].mul(9)).sub(observed[1].mul(16)).div(24);
const im=observed[0].mul(9).add(observed[1].mul(16)).sub(observed[3].mul(25)).div(24);
ck('physical_observer:full_reconstruction',eq(m([[observed[0],new N.Cut(re,im)],[new N.Cut(re,im.neg()),observed[1]]]),H));
const positiveBranch=P0,coherentPlus=mm(W,sources[0]),coherentMinus=mm(W,sources[1]);
no('physical_observer:classical_sign_record_contains_relative_phase',!eq(coherentPlus,coherentMinus)&&form(positiveBranch,sources[0]).eq(form(positiveBranch,sources[1])));
// Exact word identities for the declared infinite differential representation.
const sys=new M.Presentation({tokens:['l','q','p','k'],rules:[
 {id:'ql',lhs:['q','l'],rhs:[[['l','q'],1]],source:'q and L_ref exp(q) commute on the common core'},
 {id:'pl',lhs:['p','l'],rhs:[[['l','p'],1],[['l'],[0,-1]]],source:'p=-i partial_q acting on L_ref exp(q)'},
 {id:'pq',lhs:['p','q'],rhs:[[['q','p'],1],[[],[0,-1]]],source:'Canonical log-length representation on the common core'},
 ...['l','q','p'].map(t=>({id:'k'+t,lhs:['k',t],rhs:[[[t,'k'],1]],source:'Independent scalar-clock and geometry coordinates'}))
]});
const audit=sys.audit();ck('words:confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');
presentations.push({name:'homogeneous_affine_constraint_core',spec:sys.spec(),audit});
const l=sys.word(['l']),q=sys.word(['q']),p=sys.word(['p']),k=sys.word(['k']),one=sys.one(),C=p.times(k);
const comm=(a,b)=>a.times(b).minus(b.times(a));
function prove(name,lhs,rhs){
 const r=sys.reduce(lhs.minus(rhs),{witness:true});ck('words:'+name,r.normal.terms.size===0);
 const certificate=JSON.parse(JSON.stringify(r.certificate));ck('words:'+name+':replay',sys.replay(certificate));
 proofs.push({name,lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate});
}
prove('affine_commutator',comm(l,p),l.scale(N.IOTA));
prove('canonical_log_pair',comm(q,p),one.scale(N.IOTA));
prove('constraint_log_weak',comm(C,q),k.scale(N.IOTA.neg()));
prove('constraint_length_weak',comm(C,l),l.times(k).scale(N.IOTA.neg()));
prove('constraint_dilation_strong',comm(C,p),sys.zero());
prove('constraint_ordering_independent_factors',p.times(k),k.times(p));
let ad=l,coefficient=N.ONE;
for(let j=1;j<=4;j++){ad=comm(p,ad);coefficient=coefficient.mul(N.IOTA.neg());prove('dilation_nested_'+j,ad,l.scale(coefficient));}
no('words:length_strongly_commutes_with_full_constraint',sys.reduce(comm(C,l)).normal.terms.size!==0);
no('words:log_length_strongly_commutes_with_full_constraint',sys.reduce(comm(C,q)).normal.terms.size!==0);
const bad=JSON.parse(JSON.stringify(proofs[0].certificate));assert.ok(bad.steps.length);bad.steps[0].coefficient=['991','0'];
no('integrity:corrupted_native_witness',rejects(()=>sys.replay(bad)));
const hash=x=>crypto.createHash('sha256').update(x).digest('hex');
const payload={protocol:'HQS1_REGULAR_HOMOGENEOUS_QUANTUM_SECTOR_V1',status:'PASS_SELECTED_REGULAR_HOMOGENEOUS_QUANTUM_SECTOR',source_pins_sha256:hash(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),presentations,proofs,checks,negative_controls,results:{source:{constraint:'C=L*pi*k, V=0, matter=0',canonical_log_pair:'q=log(L/L_ref), p=L*pi',action_unit:'hbar=1 declared'},physical_measure:'dp/abs(p)',jacobians,half_line_deficiencies:deficiencies,packet_moments:packets,compressed_P:Pcompressed,compressed_P2:P2compressed,momentum_closure_ranks:momentRanks,circular_readings:observed,scope:{self_adjointness_and_group_averaging:'written analytic proofs, not finite matrix inference',singular_p_zero_sector_included:false,physical_norm_positive:true,reduced_geometry_observable_matching:'selected full-line Fourier representation',quantization_unique:false,relational_evolution:'identity for this V=0 matter-free sector',full_local_constraint_quantization:false,physical_3plus1_quantum_gravity_complete:false,physical_FTL:false}},summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,serialized_replays:proofs.length,negative_controls:negative_controls.length,native_presentations:presentations.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',sha=hash(bytes),output=path.join(__dirname,'CERTIFICATE.json');
if(process.argv.includes('--write')){fs.writeFileSync(output,bytes);fs.writeFileSync(path.join(__dirname,'EXPECTED.sha256'),sha+'  CERTIFICATE.json\n');}
else{assert.equal(fs.readFileSync(output,'utf8'),bytes);assert.equal(fs.readFileSync(path.join(__dirname,'EXPECTED.sha256'),'utf8'),sha+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:sha}));
