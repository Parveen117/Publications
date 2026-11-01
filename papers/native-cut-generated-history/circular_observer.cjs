'use strict';
// CHD-3: exact native scalar/matrix arithmetic only; no approximate fitting.
const assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const m=N.matrix,mm=(...xs)=>xs.reduce(N.mul),dg=N.dagger,eq=N.equal,sc=N.scale,f=N.F.of;
const checks=[],negative_controls=[],proofs=[],presentations=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name:'circular:'+name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name:'circular:'+name,rejected:true});};
const rejects=fn=>{try{fn();return false;}catch{return true;}};
const form=(H,x)=>mm(dg(x),H,x)[0][0],outer=x=>mm(x,dg(x));
function probes(n){
 const ps=[];
 for(let i=0;i<n;i++){const p=N.zeros(n,1);p[i][0]=N.ONE;ps.push({kind:'axis',i,p});}
 for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)for(const kind of ['real','imag']){
  const p=N.zeros(n,1);p[i][0]=N.Cut.of('3/5');p[j][0]=kind==='real'?N.Cut.of('4/5'):new N.Cut(0,'4/5');
  ps.push({kind,i,j,p});
 }
 return ps;
}
function read(H,ps=probes(H.length)){
 assert.ok(eq(H,dg(H)),'Hermitian observation contract');
 return ps.map(({p})=>{const z=form(H,p);assert.ok(z.turn.zero());return z.rad;});
}
function reconstruct(n,values){
 assert.equal(values.length,n*n,'Complete fixed probe record required');
 const v=values.map(f),H=N.zeros(n);for(let i=0;i<n;i++)H[i][i]=N.Cut.of(v[i]);
 let k=n;
 for(let i=0;i<n;i++)for(let j=i+1;j<n;j++){
  const baseline=v[i].mul(9).add(v[j].mul(16));
  H[i][j]=new N.Cut(v[k++].mul(25).sub(baseline).div(24),baseline.sub(v[k++].mul(25)).div(24));
  H[j][i]=H[i][j].dagger();
 }
 return H;
}
function basis(n){
 const out=[];
 for(let i=0;i<n;i++){const E=N.zeros(n);E[i][i]=N.ONE;out.push(E);}
 for(let i=0;i<n;i++)for(let j=i+1;j<n;j++){
  const R=N.zeros(n),J=N.zeros(n);R[i][j]=R[j][i]=N.ONE;J[i][j]=N.IOTA;J[j][i]=N.IOTA.neg();out.push(R,J);
 }
 return out;
}
const families=[];
for(const n of [1,2,3,5]){
 const ps=probes(n),bs=basis(n),observation=bs.map(H=>read(H,ps));
 ck('n'+n+':probe_count',ps.length===n*n);
 ck('n'+n+':full_rank_observation_map',N.rank(observation)===n*n);
 for(const [k,{p}]of ps.entries()){
  ck('n'+n+':unit_position:'+k,form(N.identity(n),p).eq(1));
  const P=outer(p);ck('n'+n+':binary_cut:'+k,eq(dg(P),P)&&eq(mm(P,P),P));
 }
 for(const [k,H]of bs.entries())ck('n'+n+':basis_roundtrip:'+k,eq(reconstruct(n,read(H,ps)),H));
 const duals=[];
 for(let k=0;k<n*n;k++){
  const e=Array(n*n).fill(f(0));e[k]=f(1);const D=reconstruct(n,e),values=read(D,ps);
  ck('n'+n+':dual_probe:'+k,values.every((v,j)=>v.eq(j===k?1:0)));
  no('n'+n+':omit_probe_'+k+'_still_complete',!N.isZero(D)&&values.every((v,j)=>j===k||v.zero()));
  duals.push(D);
 }
 families.push({n,probe_count:ps.length,observation_rank:N.rank(observation),probes:ps,dual_witnesses:duals});
}
const ps=probes(2),I=N.identity(2),P=m([[1,0],[0,0]]),Q=N.sub(I,P),H=m([[2,1],[1,5]]);
const complexH=m([[2,new N.Cut(1,2)],[new N.Cut(1,-2),5]]);
const thermo=read(H),complex=read(complexH);
ck('thermo:readings',thermo.every((x,j)=>x.eq([2,5,'122/25','98/25'][j])));
ck('thermo:reconstruct',eq(reconstruct(2,thermo),H));
ck('complex:readings',complex.every((x,j)=>x.eq([2,5,'122/25',2][j])));
ck('complex:reconstruct',eq(reconstruct(2,complex),complexH));
no('input:nonhermitian_contract',rejects(()=>read(m([[1,1],[0,1]]))));
no('input:incomplete_record',rejects(()=>reconstruct(2,thermo.slice(0,3))));
const xp=m([['3/5'],[new N.Cut(0,'4/5')]]),xm=m([['3/5'],[new N.Cut(0,'-4/5')]]);
const rp=outer(xp),rm=outer(xm),vp=read(rp),vm=read(rm);
no('phase:real_positions_are_complete',vp.slice(0,3).every((x,j)=>x.eq(vm[j]))&&!eq(rp,rm));
ck('phase:imaginary_position_separates',!vp[3].eq(vm[3]));
ck('phase:Gram_reconstruction',eq(reconstruct(2,vp),rp)&&eq(reconstruct(2,vm),rm));
no('phase:quadratic_readings_recover_global_phase',!eq(xp,sc(xp,N.IOTA))&&eq(rp,outer(sc(xp,N.IOTA))));
// Actual CHD coherent history, with two noncommuting source cuts.
const u=m([['3/5'],['4/5']]),B=outer(u),W=Q.concat(P),branches=[];
for(const A0 of [Q,P])for(const A1 of [N.sub(I,B),B])branches.push(mm(A1,A0));
const W2=branches.flat();
for(const [name,C]of [['single',W],['double',W2]]){
 ck('history:'+name+':isometry',eq(mm(dg(C),C),I));
 for(const [k,O]of [...basis(2),H,complexH].entries()){
  const J=mm(C,O,dg(C)),values=ps.map(({p})=>form(J,mm(C,p)).rad);
  ck('history:'+name+':complete_observer_recovery:'+k,eq(reconstruct(2,values),O));
  ck('history:'+name+':compression:'+k,eq(mm(dg(C),J,C),O));
 }
 const complement=N.sub(N.identity(C.length),mm(C,dg(C)));
 no('history:'+name+':source_probes_determine_unused_complement',!N.isZero(complement)&&N.isZero(mm(dg(C),complement,C)));
}
const fixed=N.kron(I,H),defect=N.sub(mm(dg(W),fixed,W),H),defectRead=read(defect);
ck('branch:exact_missing_information',eq(defect,m([[0,-1],[-1,0]])));
ck('branch:observer_defect_recovery',eq(reconstruct(2,defectRead),defect));
ck('branch:real_circle_failure',defectRead[2].eq('-24/25'));
no('branch:axis_only_certification',defectRead[0].zero()&&defectRead[1].zero()&&!N.isZero(defect));
// Reuse CHD-2's exact H-isometric Cayley operation, no new evolution law.
const J=m([[0,-1],[1,0]]),T=N.cayley(mm(N.inverse(H),J),'1/3');
const transports=[['identity',I,H,H],['source_Cayley',T,H,H],['history',W,I,N.identity(4)],['bad_filter',m([[1,0],[0,'3/5']]),I,I]];
const transportResults=[];
for(const [name,A,gin,gout]of transports){
 const D=N.sub(mm(dg(A),gout,A),gin),values=read(D),pass=values.every(x=>x.zero());
 ck('transport:'+name+':finite_iff',pass===N.isZero(D));
 ck('transport:'+name+':full_defect_recovery',eq(reconstruct(2,values),D));
 if(pass){const decoder=mm(N.inverse(gin),dg(A),gout);ck('transport:'+name+':left_decoder',eq(mm(decoder,A),I));}
 else no('transport:'+name+':accepted_as_invariant',!pass);
 transportResults.push({name,pass,defect:D,readings:values});
}
ck('transport:composition',eq(mm(dg(mm(T,T)),H,T,T),H));
// One invariant can survive a nonlinear erasure; the linear hypothesis is essential.
const e1=m([[1],[0]]),e2=m([[0],[1]]);
no('transport:norm_alone_implies_complete_information',form(I,e1).eq(form(I,e2))&&!eq(e1,e2));
const S=m([[2,1],[0,3]]),Si=N.inverse(S),Hp=mm(dg(Si),complexH,Si),Gp=mm(dg(Si),Si);
for(const [k,{p}]of ps.entries()){
 ck('frame:reading:'+k,form(Hp,mm(S,p)).eq(form(complexH,p)));
 ck('frame:circle:'+k,form(Gp,mm(S,p)).eq(1));
}
no('frame:metric_left_untransformed',!form(I,mm(S,e1)).eq(1));
no('frame:calibration_ignored',!form(complexH,mm(S,e1)).eq(form(complexH,e1)));
const U=m([['3/5','-4/5'],['4/5','3/5']]),V=m([[1,0],[0,N.IOTA]]),K=mm(dg(V),dg(U),V,U);
ck('loop:isometric',eq(mm(dg(K),K),I));
const loopDefect=N.sub(mm(dg(K),P,K),P);
no('loop:norm_preservation_implies_identity',!eq(K,I)&&!N.isZero(loopDefect));
ck('loop:complete_observer_detects',read(loopDefect).some(x=>!x.zero())&&eq(reconstruct(2,read(loopDefect)),loopDefect));
// Deterministic reconstruction error, all vertices of the four-reading error box.
const epsilon=f('1/100'),bound=epsilon.mul(f('31/6'));
for(let mask=0;mask<16;mask++){
 const errors=Array.from({length:4},(_,j)=>epsilon.mul(mask&(1<<j)?1:-1));
 const D=reconstruct(2,errors),rows=D.map(row=>row.reduce((s,z)=>s.add(z.gauge()),f(0)));
 ck('error:box_vertex:'+mask,rows.every(v=>v.le(bound)));
}
// Native word proofs of the scalar polarization pair, including its phase sign.
const tokens=['a','b','d','e','r','s'],rules=[];
for(let i=1;i<tokens.length;i++)for(let j=0;j<i;j++)rules.push({id:'commute_'+tokens[i]+tokens[j],lhs:[tokens[i],tokens[j]],rhs:[[[tokens[j],tokens[i]],1]],source:'Declared commuting scalar coordinates of a Hermitian two-channel form'});
const sys=new M.Presentation({tokens,rules}),audit=sys.audit();
ck('words:confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');
presentations.push({name:'circular_observer_scalar_coordinates',spec:sys.spec(),audit});
const [a,b,d,e,r,s]=tokens.map(t=>sys.word([t])),z=r.plus(s.scale(N.IOTA)),zb=r.minus(s.scale(N.IOTA));
const baseline=a.times(a).times(d).plus(b.times(b).times(e));
const real=a.times(d).times(a).plus(a.times(z).times(b)).plus(b.times(zb).times(a)).plus(b.times(e).times(b));
const imag=a.times(d).times(a).plus(a.times(z).times(b).scale(N.IOTA)).plus(b.times(zb).times(a).scale(N.IOTA.neg())).plus(b.times(e).times(b));
function prove(name,lhs,rhs){
 const out=sys.reduce(lhs.minus(rhs),{witness:true});ck('words:'+name,out.normal.terms.size===0);
 const certificate=JSON.parse(JSON.stringify(out.certificate));ck('words:'+name+':replay',sys.replay(certificate));
 proofs.push({name:'circular:'+name,lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate});
}
prove('real_polarization',real,baseline.plus(a.times(b).times(r).scale(2)));
prove('imaginary_polarization',imag,baseline.minus(a.times(b).times(s).scale(2)));
const qr=d.scale('9/25').plus(e.scale('16/25')).plus(r.scale('24/25'));
const qi=d.scale('9/25').plus(e.scale('16/25')).minus(s.scale('24/25'));
prove('real_exact_inverse',qr.scale(25).minus(d.scale(9)).minus(e.scale(16)).scale('1/24'),r);
prove('imaginary_exact_inverse',d.scale(9).plus(e.scale(16)).minus(qi.scale(25)).scale('1/24'),s);
const corrupt=JSON.parse(JSON.stringify(proofs[0].certificate));assert.ok(corrupt.steps.length);corrupt.steps[0].coefficient=['919','0'];
no('integrity:corrupted_polarization_witness',rejects(()=>sys.replay(corrupt)));
no('phase:wrong_imaginary_sign',sys.reduce(imag.minus(baseline.plus(a.times(b).times(s).scale(2)))).normal.terms.size!==0);
module.exports={checks,negative_controls,proofs,presentations,results:{families,thermo_readings:thermo,complex_readings:complex,branch_defect:defect,branch_defect_readings:defectRead,transport:transportResults,loop:K,loop_observable_defect:loopDefect,error_constant_n2:'31/6',scope:{target:'Complete finite quadratic-observer verification',closed_within_declared_contract:true,general_proofs:'written; native certificates verify specified finite bases and identities',fixed_probe_minimum:'n^2 for unrestricted Hermitian observable; n^2-1 independent data if trace known',global_phase_recovered:false,single_copy_tomography:false,physical_dynamics_selected:false,quantum_gravity_complete:false}}};
