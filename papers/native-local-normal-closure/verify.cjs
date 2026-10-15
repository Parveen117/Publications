'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs'),M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const f=N.F.of,I=N.identity(2),Z=N.matrix([[1,0],[0,-1]]),R=N.matrix([[0,-1],[1,0]]);
const mm=(...a)=>a.reduce(N.mul),sc=N.scale,eq=N.equal,zero=N.isZero,cm=N.commutator;
const tr=a=>a.reduce((s,row,i)=>s.add(row[i]),N.ZERO),real=z=>{assert.ok(z.turn.zero());return z.rad;};
const checks=[],negative_controls=[],proofs=[],presentations=[],systems=new Map(),results={};
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
// Coefficient helper reused verbatim from AGD-1 apart from A_N's coframe factor.
// Polynomial differential operators in x; all coefficients are native 2x2 matrices.
const zmat=()=>N.zeros(2),ptrim=p=>{while(p.length>1&&zero(p.at(-1)))p.pop();return p;};
const pa=(p,r)=>ptrim(Array.from({length:Math.max(p.length,r.length)},(_,i)=>N.add(p[i]||zmat(),r[i]||zmat())));
const ps=(p,c)=>ptrim(p.map(x=>sc(x,c)));
function pm(p,r){const out=Array.from({length:p.length+r.length-1},zmat);p.forEach((a,i)=>r.forEach((b,j)=>out[i+j]=N.add(out[i+j],mm(a,b))));return ptrim(out);}
const pd=p=>p.length===1?[zmat()]:p.slice(1).map((x,i)=>sc(x,i+1));
const peq=(p,r)=>zero(pa(p,ps(r,-1)).reduce((s,x)=>N.add(s,N.matrix(x.map(row=>row.map(v=>v.norm2())))),zmat()));
const scalar=p=>p.map(x=>sc(I,x)),left=(m,p)=>p.map(x=>mm(m,x));
function compose(F,G){const[a,b]=F,[c,d]=G;return[pm(a,c),pa(pa(pm(a,pd(c)),pm(a,d)),pm(b,c)),pa(pm(a,pd(d)),pm(b,d))];}
const commOp=(F,G)=>{const x=compose(F,G),y=compose(G,F);return x.map((p,i)=>pa(p,ps(y[i],-1)));};
const Dop=f=>[f,ps(pd(f),'1/2')],Aop=(n,e,g)=>[left(Z,n),pa(ps(left(Z,pd(n)),'1/2'),ps(left(R,pm(n,e)),g))];
const lift=op=>[[zmat()],...op],opEq=(a,b)=>a.every((p,i)=>peq(p,b[i]));
const br=(a,b)=>pa(pm(a,pd(b)),ps(pm(b,pd(a)),-1)),pzero=[zmat()],opAdd=(a,b)=>a.map((p,i)=>pa(p,b[i])),opSub=(a,b)=>a.map((p,i)=>pa(p,ps(b[i],-1))),mult=p=>[pzero,p],g='2/3';
const sum=(...xs)=>xs.reduce(pa,pzero),sub=(a,b)=>pa(a,ps(b,-1)),prod=(...xs)=>xs.reduce(pm),integ=p=>p.reduce((v,a,i)=>v.add(real(a[0][0]).div(i+1)),f(0));
const pb=(a,b)=>sum(pm(a[0],b[1]),ps(pm(a[1],b[0]),-1),pm(a[2],b[3]),ps(pm(a[3],b[2]),-1));
function hgGrad(n,q,p,phi,k,U,Up,a='1',b='1',c='-1'){
 return [sum(ps(pd(pm(n,pd(phi))),f(b).neg()),ps(pm(n,U),2)),ps(pm(n,k),a),sum(ps(pd(pm(n,pd(q))),f(b).neg()),ps(pd(pd(n)),c),pm(n,Up)),ps(pm(n,p),a)];
}
const dgGrad=(ff,q,p,phi,k)=>[ps(pd(pm(ff,p)),-1),pa(pm(ff,pd(q)),pd(ff)),ps(pd(pm(ff,k)),-1),pm(ff,pd(phi))];
const hd=(q,p,phi,k,U,a='1',b='1',c='-1')=>sum(ps(pm(p,k),a),ps(pm(pd(q),pd(phi)),b),ps(pd(pd(phi)),c),U);
const dd=(q,p,phi,k)=>sum(pm(p,pd(q)),pm(k,pd(phi)),ps(pd(p),-1));
const nn=scalar([0,0,1,-2,1]),mm0=scalar([0,0,0,1,-3,3,-1]),ff=scalar([0,0,2,-3,0,1]),omega=br(nn,mm0);
const configs=[{q:[1,2,1],p:[2,-1,1],phi:[1,1,2],k:[1,0,1],potential:'zero'},{q:[0],p:[1,2],phi:[1,-1,1],k:[2,1],potential:'quadratic'},{q:[0],p:[2,0,1],phi:[1,2],k:[1,-1],potential:'linear'}];
for(let z=0;z<configs.length;z++){
 const cfg=configs[z],q=scalar(cfg.q),p=scalar(cfg.p),phi=scalar(cfg.phi),k=scalar(cfg.k),U=cfg.potential==='zero'?pzero:cfg.potential==='linear'?phi:pm(phi,phi),Up=cfg.potential==='zero'?pzero:cfg.potential==='linear'?scalar([1]):ps(phi,2),Dden=dd(q,p,phi,k),Hden=hd(q,p,phi,k,U),G=sm=>hgGrad(sm,q,p,phi,k,U,Up),DG=sm=>dgGrad(sm,q,p,phi,k),normal=pb(G(nn),G(mm0));
 ck('geometry:normal_total_derivative:'+z,peq(sub(normal,pm(omega,Dden)),pd(pm(omega,p))));
 ck('geometry:normal_integrated:'+z,integ(normal).eq(integ(pm(omega,Dden))));
 ck('geometry:spatial_integrated:'+z,integ(pb(DG(ff),DG(mm0))).eq(integ(pm(br(ff,mm0),Dden))));
 ck('geometry:mixed_integrated:'+z,integ(pb(DG(ff),G(nn))).eq(integ(pm(br(ff,nn),Hden))));
 const dq=pa(pm(ff,pd(q)),pd(ff)),dp=pd(pm(ff,p)),dphi=pm(ff,pd(phi)),dk=pd(pm(ff,k)),dU=pa(pm(ff,pd(U)),ps(pm(pd(ff),U),2));
 const deltaH=sum(pm(dp,k),pm(p,dk),pm(pd(dq),pd(phi)),pm(pd(q),pd(dphi)),ps(pd(pd(dphi)),-1),dU);
 ck('geometry:density_weight_two:'+z,peq(deltaH,pa(pm(ff,pd(Hden)),ps(pm(pd(ff),Hden),2))));
 for(const [a,b,c]of [['1','1','-1'],['2','1/2','-1/2'],['1','1','0'],['2','1','-1']]){
  const label=[z,a,b,c].join(':'),x=hgGrad(nn,q,p,phi,k,U,Up,a,b,c),y=hgGrad(mm0,q,p,phi,k,U,Up,a,b,c),actual=integ(pb(x,y)),expected=integ(pa(ps(pm(omega,Dden),f(a).mul(b)),ps(pm(omega,pd(p)),f(a).mul(f(b).add(c)))));
  ck('ansatz:general_normal_formula:'+label,actual.eq(expected));
  if(f(a).mul(b).eq(1)&&f(b).add(c).zero())ck('ansatz:normalized_closure:'+label,actual.eq(integ(pm(omega,Dden))));
  else if(z===0)no('ansatz:wrong_coefficients_close:'+label,!actual.eq(integ(pm(omega,Dden))));
 }
 const badDensity=sum(pm(dp,k),pm(p,dk),pm(pd(dq),pd(phi)),pm(pd(q),pd(dphi)),dU),badH=sum(pm(p,k),pm(pd(q),pd(phi)),U),densityDefect=sub(badDensity,pa(pm(ff,pd(badH)),ps(pm(pd(ff),badH),2)));
 ck('ansatz:omitted_improvement_defect:'+z,peq(densityDefect,pm(pd(pd(ff)),pd(phi))));
 no('ansatz:omit_improvement_covariance:'+z,!peq(densityDefect,pzero));
}
// Nonconstant coframe and inherited matter blocks: geometry cancels the mixed source.
for(const [z,ee,kk]of [[0,[1,0,1],[1,1]],[1,[2,1,1],[2,0,1]],[2,[1,1,2],[1,-1]]]){
 const e=scalar(ee),k=scalar(kk),A=Aop(nn,e,g),B=Aop(mm0,e,g),massVariation=mult(ps(left(R,pm(nn,pd(pm(ff,e)))),g));
 ck('matter:normal_bracket:'+z,opEq(commOp(A,B),lift(Dop(omega))));
 ck('matter:total_mixed_bracket:'+z,opEq(opSub(commOp(Dop(ff),A),lift(massVariation)),lift(Aop(br(ff,nn),e,g))));
 no('matter:freeze_geometry_mixed:'+z,!opEq(commOp(Dop(ff),A),lift(Aop(br(ff,nn),e,g))));
 const j=scalar([0,'4/3']),crossNM=ps(prod(nn,mm0,e,k,j),-1),crossMN=prod(mm0,nn,e,k,j);
 ck('source:antisymmetric_normal_cross_cancels:'+z,peq(pa(crossNM,crossMN),pzero));
 no('source:individual_cross_is_zero:'+z,!peq(crossNM,pzero));
 const de=scalar([1,1]),difference=sub(prod(nn,pa(e,de),j),prod(nn,e,j));
 ck('source:coframe_variation_retained:'+z,peq(difference,prod(nn,j,de)));
}
// Normal Hamilton-vector-field commutator in the positive-coframe chart.
// qprime and qdoubleprime are supplied exact local jets; matter does not enter kdot.
for(const [z,q1,q2]of [[0,scalar([1,1]),scalar([1])],[1,scalar([0,2]),scalar([2])]]){
 const N1=pd(nn),M1=pd(mm0),N2=pd(N1),M2=pd(M1),Uphi=scalar([2,1]),kdotN=sum(pm(N1,q1),pm(nn,q2),N2,ps(pm(nn,Uphi),-1)),kdotM=sum(pm(M1,q1),pm(mm0,q2),M2,ps(pm(mm0,Uphi),-1));
 const relative=sub(pm(mm0,kdotN),pm(nn,kdotM));
 ck('normal:field_direction_restored:'+z,peq(relative,ps(pa(pm(omega,q1),pd(omega)),-1)));
}
// Two exactly constrained initial jets with different selected potentials.
const e=f(1),pi=f(0),phi=f(0),k=f(0),N0=f(1);
const alternatives=[];
for(const [label,V,Vprime]of [['zero','0','0'],['linear','0','1']]){
 const h=e.mul(pi).mul(k).add(e.pow(2).mul(V)),edot=N0.mul(e).mul(k),kdot=N0.mul(e.pow(2)).mul(Vprime).neg(),eddot=edot.mul(k).add(e.mul(kdot));
 ck('selection:initial_constraint:'+label,h.zero());
 ck('selection:initial_velocity:'+label,edot.zero());
 ck('selection:acceleration:'+label,eddot.eq(label==='zero'?0:-1));
 alternatives.push({potential:label,H_initial:h,D_initial:'0',e_dot:edot,p_phi_dot:kdot,e_ddot:eddot});
}
no('selection:constraint_algebra_selects_potential',!alternatives[0].e_ddot.eq(alternatives[1].e_ddot));
const negativeDensity=hd(scalar([0]),scalar([1]),scalar([0]),scalar([-1]),pzero);
no('positivity:kinematic_hamiltonian_positive',real(negativeDensity[0][0][0]).eq(-1));
results.selection={models:alternatives,kinematic_negative_density:'-1'};
results.constraints={geometry_pair_count:2,spatial_directions:1,brackets:['DD=D','DH=H','HH=D'],matter:'classical commuting complex canonical fields',quantum_closure_claimed:false};
function register(name,tokens){const rules=[];for(let j=1;j<tokens.length;j++)for(let i=0;i<j;i++)rules.push({id:'central_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'Commuting canonical scalar jets'});const s=new M.Presentation({tokens,rules}),audit=s.audit();ck(name+':confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name,spec:s.spec(),audit});systems.set(s.hash(),s);return{s,w:Object.fromEntries(tokens.map(t=>[t,s.word([t])]))};}
function prove(name,lhs,rhs){const s=lhs.system,r=s.reduce(lhs.minus(rhs),{witness:true});ck(name,r.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(r.certificate));ck(name+':replay',s.replay(certificate));proofs.push({name,lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate});}
{
 const{s,w:{n,m,dn,dm,t,u,k}}=register('normal_bracket_scalar_sector',['n','m','dn','dm','t','u','k']),gn=dn.times(t).scale(-1).plus(n.times(u)),gm=dm.times(t).scale(-1).plus(m.times(u)),omega=n.times(dm).minus(m.times(dn));
 prove('words:normal_sector_derivative',gn.times(m).times(k).minus(n.times(k).times(gm)),omega.times(t).times(k));
 prove('words:potential_antisymmetry',n.times(u).times(m).times(k).minus(n.times(k).times(m).times(u)),s.zero());
 prove('words:matter_geometry_cross_cancel',n.times(m).times(t).times(u).times(k).scale(-1).plus(m.times(n).times(t).times(u).times(k)),s.zero());
}
{
 const{s,w:{f,df,ddf,q1,q2,t1,t2,t3}}=register('improved_density_jets',['f','df','ddf','q1','q2','t1','t2','t3']);
 const dq1=f.times(q2).plus(df.times(q1)).plus(ddf),dt1=f.times(t2).plus(df.times(t1)),dt2=f.times(t3).plus(df.times(t2).scale(2)).plus(ddf.times(t1));
 prove('words:improvement_weight_two',dq1.times(t1).plus(q1.times(dt1)).minus(dt2),f.times(q2.times(t1).plus(q1.times(t2)).minus(t3)).plus(df.times(q1.times(t1).minus(t2)).scale(2)));
}
const chosen=proofs.find(x=>x.certificate.steps.length),bad=JSON.parse(JSON.stringify(chosen.certificate));bad.steps[0].coefficient=['991','0'];
no('integrity:corrupted_word_proof',(()=>{try{systems.get(chosen.certificate.presentation).replay(bad);return false;}catch{return true;}})());
const hash=x=>crypto.createHash('sha256').update(x).digest('hex'),payload={protocol:'LNC1_CLASSICAL_LOCAL_NORMAL_CONSTRAINT_CLOSURE_V1',status:'PASS_NATIVE_CLASSICAL_LOCAL_CONSTRAINT_ALGEBRA',scope:'Selected smooth classical 1+1 canonical coframe-scalar-matter model. Written all-smearing constraint proof, native algebra and exact polynomial certificates. No quantum anomaly/positivity or physical 3+1 completion claim.',source_pins_sha256:hash(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),presentations,proofs,checks,negative_controls,results,summary:{exact_checks:checks.length,negative_controls:negative_controls.length,native_presentations:presentations.length,symbolic_certificates:proofs.length,serialized_replays:proofs.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',sha=hash(bytes),file=path.join(__dirname,'CERTIFICATE.json');
if(process.argv.includes('--write')){fs.writeFileSync(file,bytes);fs.writeFileSync(path.join(__dirname,'EXPECTED.sha256'),sha+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(file,'utf8'),bytes);assert.equal(fs.readFileSync(path.join(__dirname,'EXPECTED.sha256'),'utf8'),sha+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:sha}));
