'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs'),M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const f=N.F.of,checks=[],negative_controls=[],proofs=[],presentations=[],systems=new Map(),results={};
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const rejects=fn=>{try{fn();return false;}catch{return true;}};
const mm=(...xs)=>xs.reduce(N.mul),sc=N.scale,eq=N.equal,tr=a=>a.reduce((s,row,i)=>s.add(row[i]),N.ZERO);
const I=N.identity(2),K=N.matrix([[1,0],[0,-1]]),R=N.matrix([[0,-1],[1,0]]),J=N.mul(R,K),P=sc(N.add(I,J),'1/2'),Q=N.sub(I,P);
const Hplus=N.matrix([['5/2',1],[1,'3/2']]),Hminus=N.matrix([['3/2',1],[1,'5/2']]);
const determinant=h=>h[0][0].mul(h[1][1]).sub(h[0][1].mul(h[1][0]));
const positive2=h=>h[0][0].turn.zero()&&f(0).le(h[0][0].rad)&&!h[0][0].rad.zero()&&f(0).le(determinant(h).rad)&&!determinant(h).rad.zero();
ck('native:J_squared',eq(mm(J,J),I));
ck('native:diagonal_projector',eq(mm(P,P),P));
ck('native:axis_response_reconstruction',eq(Hplus,N.add(N.add(mm(P,Hplus,P),mm(P,Hplus,Q)),N.add(mm(Q,Hplus,P),mm(Q,Hplus,Q)))));
for(const [j,h]of [Hplus,Hminus].entries()){
 ck('axis:positive:'+j,positive2(h));
 ck('axis:determinant:'+j,determinant(h).eq('11/4'));
 ck('axis:ratio:'+j,determinant(h).div(h[0][0].mul(h[1][1])).eq('11/15'));
 ck('axis:diagonal_compression:'+j,eq(mm(P,h,P),sc(P,3)));
}
ck('transport:axis_swap',eq(mm(J,Hplus,J),Hminus));
for(const [j,U]of [I,J].entries()){
 ck('transport:unitary:'+j,eq(mm(N.dagger(U),U),I));
 ck('transport:intertwines_diagonal:'+j,eq(mm(P,U),mm(U,P)));
 ck('transport:same_retained_action:'+j,eq(mm(U,P),P));
}
no('axis:diagonal_determines_full_response',!eq(Hplus,Hminus));
no('transport:retained_action_selects_hidden_transport',!eq(mm(I,Q),mm(J,Q)));
const cut=N.matrix([[1,0]]),descended=N.matrix([[2]]),lifted=mm(N.dagger(cut),descended,cut),blind=N.matrix([[0],[1]]);
ck('descent:annihilates_kernel',N.isZero(mm(lifted,blind)));
ck('descent:bilinear_reconstruction',eq(lifted,N.matrix([[2,0],[0,0]])));
no('descent:symmetry_suffices',!N.isZero(mm(Hplus,blind)));
// Matrix polynomials. Only loop indices/degrees use JS integers; coefficients use N.
const trim=p=>{while(p.length>1&&N.isZero(p.at(-1)))p.pop();return p;};
const ev=(p,t)=>p.reduceRight((s,m)=>N.add(sc(s,t),m),N.zeros(p[0].length));
const pd=p=>p.length===1?[N.zeros(p[0].length)]:p.slice(1).map((m,k)=>sc(m,k+1));
const peq=(p,q)=>{const n=p[0].length;return Array.from({length:Math.max(p.length,q.length)},(_,k)=>eq(p[k]||N.zeros(n),q[k]||N.zeros(n))).every(Boolean);};
function dress(C,w){const n=w.length,out=Array.from({length:C.length+2*Math.max(...w)},()=>N.zeros(n));C.forEach((m,k)=>m.forEach((row,i)=>row.forEach((z,j)=>out[k+w[i]+w[j]][i][j]=z)));return trim(out);}
function extend(H,w){const n=w.length;assert.equal(H[0].length,n);const out=Array.from({length:H.length},()=>N.zeros(n));for(let k=0;k<H.length;k++)for(let i=0;i<n;i++)for(let j=0;j<n;j++){const degree=k-w[i]-w[j],z=H[k][i][j];if(degree<0)assert.ok(z.zero(),'Non-removable negative power');else out[degree][i][j]=z;}return trim(out);}
const diag=xs=>N.matrix(xs.map((x,i)=>xs.map((_,j)=>i===j?x:0)));
const coeff=(H,k,i,j)=>H[k]?.[i]?.[j]||N.ZERO;
const cases=[{name:'one_collapsed_axis',w:[0,1],C:[[[2,1],[1,2]]]},{name:'variable_response',w:[0,1],C:[[[2,1],[1,2]],[[1,2],[2,0]],[[0,1],[1,1]]]},{name:'three_axis_grades',w:[0,1,2],C:[[[3,1,0],[1,3,1],[0,1,3]],[[1,0,2],[0,1,0],[2,0,1]]]}];
results.graded_cases=[];
for(const c of cases){
 const C=c.C.map(N.matrix),H=dress(C,c.w),recovered=extend(H,c.w),n=c.w.length;
 ck('graded:exact_polynomial_extension:'+c.name,peq(C,recovered));
 for(let i=0;i<n;i++)for(let j=0;j<n;j++)for(let k=0;k<C.length;k++)ck('graded:jet_rule:'+c.name+':'+i+':'+j+':'+k,coeff(recovered,k,i,j).eq(coeff(H,k+c.w[i]+c.w[j],i,j)));
 for(const t0 of ['-2','-1/3','1/2','3']){
  const t=f(t0),A=diag(c.w.map(w=>t.pow(w))),Ai=N.inverse(A),Ht=ev(H,t),Ct=ev(C,t),Hp=ev(pd(H),t),Cp=ev(pd(C),t),Gamma=diag(c.w.map(w=>f(w).div(t)));
  ck('graded:direct_congruence:'+c.name+':'+t0,eq(mm(N.dagger(Ai),Ht,Ai),Ct));
  ck('graded:symmetric:'+c.name+':'+t0,eq(Ct,N.dagger(Ct)));
  ck('frame:derivative_connection:'+c.name+':'+t0,eq(N.sub(N.sub(mm(N.dagger(Ai),Hp,Ai),mm(N.dagger(Gamma),Ct)),mm(Ct,Gamma)),Cp));
  const x=N.matrix(c.w.map((_,i)=>[i+1])),y=mm(A,x);
  ck('potential:preserved:'+c.name+':'+t0,eq(mm(N.dagger(x),Ht,x),mm(N.dagger(y),Ct,y)));
 }
 const bad=H.map(N.matrix),j=c.w.findIndex(w=>w>0),k=2*c.w[j]-1;
 bad[k][j][j]=bad[k][j][j].add(1);
 no('graded:lower_order_contamination:'+c.name,rejects(()=>extend(bad,c.w)));
 results.graded_cases.push({name:c.name,weights:c.w,raw_coefficients:H,recovered_coefficients:recovered,required_boundary_depth:2*Math.max(...c.w)});
}
const Cplus=N.matrix([[2,1],[1,2]]),Cminus=N.matrix([[2,-1],[-1,2]]),Cother=N.matrix([[2,1],[1,3]]),Bplus=dress([Cplus],[0,1]),Bminus=dress([Cminus],[0,1]),Bother=dress([Cother],[0,1]);
ck('singular:same_collapsed_point',eq(Bplus[0],Bminus[0]));
ck('singular:positive_reconstructed',positive2(extend(Bplus,[0,1])[0]));
ck('singular:ratio_sign_blind',determinant(Cplus).div(Cplus[0][0].mul(Cplus[1][1])).eq(determinant(Cminus).div(Cminus[0][0].mul(Cminus[1][1]))));
no('singular:value_recovers_cross_sign',!eq(extend(Bplus,[0,1])[0],extend(Bminus,[0,1])[0]));
ck('singular:first_jet_same_second_axis_ambiguity',eq(Bplus[0],Bother[0])&&eq(Bplus[1],Bother[1]));
no('singular:first_jet_suffices',!eq(extend(Bplus,[0,1])[0],extend(Bother,[0,1])[0]));
no('singular:arbitrary_symmetric_response_extends',rejects(()=>extend([I],[0,1])));
no('singular:weights_are_optional',!eq(extend(Bplus,[0,0])[0],extend(Bplus,[0,1])[0]));
const degeneratingC=[N.matrix([[1,0],[0,0]]),N.zeros(2),N.matrix([[0,0],[0,1]])];
ck('positivity:away_from_zero',positive2(ev(degeneratingC,f(1))));
no('positivity:positive_boundary_automatic',determinant(extend(dress(degeneratingC,[0,1]),[0,1])[0]).zero());
const t=f(2),Ai=diag([1,'1/2']);
no('frame:omit_connection',!N.isZero(mm(Ai,ev(pd(Bplus),t),Ai)));
// Partial Legendre change preserves the signed TVSP/EM potential differential.
for(const [j,raw]of [[2,1,2],[3,-1,2]].entries()){
 const [a,b,c]=raw.map(f),S=N.matrix([[a.sub(b.pow(2).div(c)),b.div(c)],[b.div(c),f(-1).div(c)]]);
 ck('adapter:reciprocity:'+j,eq(S,N.dagger(S)));
 for(const [D0,H0]of [['1','3'],['-2','1/2']]){
  const D=f(D0),Hm=f(H0),B=Hm.sub(b.mul(D)).div(c),E=a.mul(D).add(b.mul(B)),grad=mm(S,N.matrix([[D],[Hm]]));
  ck('adapter:gradient_E:'+j+':'+D0,grad[0][0].eq(E));
  ck('adapter:gradient_minus_B:'+j+':'+D0,grad[1][0].eq(B.neg()));
  if(!B.zero())no('adapter:unsigned_conjugate_relabel:'+j+':'+D0,!grad[1][0].eq(B));
  const w=a.mul(D.pow(2)).add(b.mul(D).mul(B).mul(2)).add(c.mul(B.pow(2))).div(2),psi=w.sub(B.mul(Hm)),q=N.matrix([[D],[Hm]]);
  ck('adapter:Legendre_potential:'+j+':'+D0,mm(N.dagger(q),S,q)[0][0].eq(psi.mul(2)));
 }
 no('adapter:energy_positivity_survives_partial_Legendre:'+j,!f(0).le(S[1][1].rad));
}
// Finite oriented triangle complex: a chain-map control, not a physical EM grid.
const d0=N.matrix([[-1,1,0],[0,-1,1],[1,0,-1]]),d1=N.matrix([[1,1,1]]),T0=N.matrix([[0,1,0],[1,0,0],[0,0,1]]),T1=N.matrix([[-1,0,0],[0,0,-1],[0,-1,0]]),T2=N.matrix([[-1]]);
ck('chain:boundary_squared_zero',N.isZero(mm(d1,d0)));
ck('chain:first_intertwiner',eq(mm(d0,T0),mm(T1,d0)));
ck('chain:second_intertwiner',eq(mm(d1,T1),mm(T2,d1)));
const chi=N.matrix([[2,1,0],[1,3,1],[0,1,4]]),chiPrime=mm(T1,chi,N.inverse(T1)),F=mm(d0,N.matrix([[1],[2],[4]])),G=mm(chi,F),current=mm(d1,G),Fp=mm(T1,F),Gp=mm(T1,G);
ck('chain:constitutive_intertwiner',eq(mm(chiPrime,T1),mm(T1,chi)));
ck('chain:closed_field_transported',N.isZero(mm(d1,Fp)));
ck('chain:source_transported',eq(mm(d1,Gp),mm(T2,current)));
ck('chain:constitutive_field_transported',eq(mm(chiPrime,Fp),Gp));
no('chain:rename_without_chi_transport',!eq(mm(chi,Fp),Gp));
const O=diag([1,1,0]),compressed=mm(O,chi,O);
ck('chain:compressed_response_symmetric',eq(compressed,N.dagger(compressed)));
no('chain:symmetry_preserves_flux_equation',!N.isZero(mm(d1,O,F)));
// Native symbolic identities, with exact serialized witnesses.
function register(name,tokens){const rules=[];for(let j=1;j<tokens.length;j++)for(let i=0;i<j;i++)rules.push({id:'commute_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'Central radial coefficients of the native response element'});const s=new M.Presentation({tokens,rules}),audit=s.audit();ck(name+':confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name,spec:s.spec(),audit});systems.set(s.hash(),s);return{s,w:Object.fromEntries(tokens.map(t=>[t,s.word([t])]))};}
function prove(name,lhs,rhs){const s=lhs.system,r=s.reduce(lhs.minus(rhs),{witness:true});ck(name,r.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(r.certificate));ck(name+':replay',s.replay(certificate));proofs.push({name,lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate});}
{
 const{s,w:{a,b,c,t,x,y,d}}=register('graded_response_coefficients',['a','b','c','t','x','y','d']);
 prove('words:graded_determinant',a.times(c.times(t).times(t)).minus(b.times(t).times(b).times(t)),t.times(t).times(a.times(c).minus(b.times(b))));
 prove('words:quadratic_potential',a.times(x).times(x).plus(b.times(t).times(x).times(y).scale(2)).plus(c.times(t).times(t).times(y).times(y)),x.times(a.times(x).plus(b.times(t.times(y)))).plus(t.times(y).times(b.times(x).plus(c.times(t.times(y))))));
 prove('words:frame_second_axis_cancellation',c.times(t).scale(2).plus(t.times(d).times(t)).minus(t.times(c).scale(2)),t.times(t).times(d));
 prove('words:Legendre_denominator_cleared',a.times(c).minus(b.times(b)).times(x).plus(b.times(y)),c.times(a.times(x)).plus(b.times(y.minus(b.times(x)))));
 prove('words:corner_differential',a.times(x).plus(b.times(y)).minus(x.times(a)).minus(d.times(c)),b.times(y).minus(c.times(d)));
}
const chosen=proofs.find(q=>q.certificate.steps.length),bad=JSON.parse(JSON.stringify(chosen.certificate));bad.steps[0].coefficient=['991','0'];
no('integrity:corrupt_native_proof',rejects(()=>systems.get(chosen.certificate.presentation).replay(bad)));
results.scope={chart:'declared analytic scaling family',quantum_gravity_complete:false,physical_spacetime_singularity_resolved:false,weights_derived_from_cut_dynamics:false,full_Maxwell_PDE_transport_constructed:false,analytic_proof:'THEOREM.md section 3'};
const hash=x=>crypto.createHash('sha256').update(x).digest('hex'),payload={protocol:'SRT1_GRADED_SINGULAR_RESPONSE_TRANSPORT_V1',status:'PASS_NATIVE_GRADED_SINGULAR_RESPONSE_TRANSPORT',source_pins_sha256:hash(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),presentations,proofs,checks,negative_controls,results,summary:{exact_checks:checks.length,negative_controls:negative_controls.length,native_presentations:presentations.length,symbolic_certificates:proofs.length,serialized_replays:proofs.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',sha=hash(bytes),file=path.join(__dirname,'CERTIFICATE.json');
if(process.argv.includes('--write')){fs.writeFileSync(file,bytes);fs.writeFileSync(path.join(__dirname,'EXPECTED.sha256'),sha+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(file,'utf8'),bytes);assert.equal(fs.readFileSync(path.join(__dirname,'EXPECTED.sha256'),'utf8'),sha+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:sha}));
