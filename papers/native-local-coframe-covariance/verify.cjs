'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs'),M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const f=N.F.of,I=N.identity(2),Z=N.matrix([[1,0],[0,-1]]),S=N.matrix([[0,1],[1,0]]),R=N.matrix([[0,-1],[1,0]]);
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
const data=[[[0,1],[1],[0,1],[1,0,1],[1,1],[1]],[[1,0,1],[1,1],[0,1,1],[2,1,1],[0,1],[0,1]],[[1,2,1],[1,-1],[1,0,1],[3,0,2],[1,1,1],[1,0,1]]];
for(let z=0;z<data.length;z++){
 const [ff,nn,mm0,ee,bb,tt]=data[z],f0=scalar(ff),n=scalar(nn),m=scalar(mm0),e=scalar(ee),beta=scalar(bb),ft=scalar(tt),one=scalar([1]),de=pd(pm(f0,e)),h=m,b=br(n,m),A=Aop(n,e,g),D=Dop(f0);
 ck('algebra:DD:'+z,opEq(commOp(D,Dop(h)),lift(Dop(br(f0,h)))));
 ck('algebra:AA_variable_coframe:'+z,opEq(commOp(A,Aop(m,e,g)),lift(Dop(b))));
 const extra=mult(ps(left(R,pm(n,de)),g)),raw=commOp(D,A),expected=lift(opAdd(Aop(br(f0,n),e,g),extra));
 ck('algebra:mixed_residual:'+z,opEq(raw,expected));
 ck('algebra:mixed_repaired:'+z,opEq(opSub(raw,lift(extra)),lift(Aop(br(f0,n),e,g))));
 no('algebra:freeze_coframe:'+z,!opEq(raw,lift(Aop(br(f0,n),e,g))));
 const wrongDe=pm(f0,pd(e)),wrong=opSub(raw,lift(mult(ps(left(R,pm(n,wrongDe)),g)))),weightResidual=mult(ps(left(R,pm(pm(n,e),pd(f0))),g));
 ck('weight:wrong_weight_residual:'+z,opEq(wrong,lift(opAdd(Aop(br(f0,n),e,g),weightResidual))));
 no('weight:scalar_coframe_suffices:'+z,!opEq(wrong,lift(Aop(br(f0,n),e,g))));
 ck('field:V_bracket_sign:'+z,peq(pa(pd(pm(h,de)),ps(pd(pm(f0,pd(pm(h,e)))),-1)),ps(pd(pm(br(f0,h),e)),-1)));
 const deltaN=br(f0,n),deltaBeta=pa(ft,br(f0,beta)),deltaB=opAdd(opAdd(Aop(deltaN,e,g),extra),Dop(deltaBeta)),B=opAdd(A,Dop(beta));
 ck('action:time_dependent_covariance:'+z,opEq(lift(deltaB),opAdd(lift(Dop(ft)),commOp(D,B))));
 ck('lapse:alpha_scalar:'+z,peq(pa(pm(deltaN,e),pm(n,de)),pm(f0,pd(pm(n,e)))));
 ck('metric:density_weight_two:'+z,peq(ps(pm(e,de),2),pa(pm(f0,pd(pm(e,e))),ps(pm(pd(f0),pm(e,e)),2))));
 const missing=pd(pm(b,e));
 no('normal:full_field_closure:'+z,!peq(missing,pzero));
 ck('massless:old_mixed_special_case:'+z,opEq(commOp(D,Aop(n,e,0)),lift(Aop(br(f0,n),e,0))));
 const old=commOp(D,Aop(n,one,g)),oldextra=mult(ps(left(R,pm(n,pd(f0))),g));
 ck('inherit:AGD_residual:'+z,opEq(old,lift(opAdd(Aop(br(f0,n),one,g),oldextra))));
}
const n0=scalar([1]),m0=scalar([0,1]),e0=scalar([1,0,1]),b0=br(n0,m0),missing0=pd(pm(b0,e0));
ck('normal:explicit_missing_field_direction',peq(missing0,scalar([0,2])));
results.constraint={mixed_residual_cancelled:true,closed_brackets:['spatial-spatial','spatial-normal'],normal_normal_missing_field_direction:'V_(N Mprime - M Nprime)',example:{N:'1',M:'x',e:'1+x^2',missing_on_e:'2x'}};

const circle=[];
for(const [cs,ss]of [['1','0'],['0','1'],['3/5','4/5'],['5/13','12/13']])for(const speed of ['0','2','3/5']){
 const P=sc(N.add(I,N.add(sc(Z,cs),sc(S,ss))),'1/2'),pt=sc(N.add(sc(Z,f(ss).neg()),sc(S,cs)),'1/2'),px=sc(pt,speed),W=real(tr(mm(px,px))).div(2),tag=[cs,ss,speed].join(':');
 ck('cut:circle_projection:'+tag,eq(mm(P,P),P));
 ck('cut:metric_coefficient:'+tag,W.eq(f(speed).pow(2).div(4)));
 ck('cut:density_squared:'+tag,f(speed).div(2).pow(2).eq(W));
 if(speed==='0')no('cut:constant_cut_nondegenerate:'+tag,W.zero());
 circle.push({cosine:cs,sine:ss,theta_prime:speed,W});
}
results.cut_metric={scope:'Parallel reference chart; positive orientation patch for the square root',circle};
// Exact local bilinear J for psi=(1,i*x), with a compact cutoff understood.
const a=scalar([1]),bb=ps(scalar([0,1]),N.IOTA),dag=p=>p.map(N.dagger),j=ps(pa(ps(pm(dag(a),bb),-1),pm(dag(bb),a)),N.IOTA.mul(g)),n=scalar([1,1]),force=ps(pd(pm(n,j)),'1/2');
ck('source:bilinear_real',j.every(x=>eq(x,N.dagger(x))));
ck('source:J_polynomial',peq(j,scalar([0,'4/3'])));
ck('source:cut_force',peq(force,scalar(['2/3','4/3'])));
const integral=p=>p.reduce((s,a,i)=>s.add(real(a[0][0]).div(i+1)),f(0));
for(const [z,variation]of [[0,[0,1,-1]],[1,[0,0,1,-1]],[2,[0,1,0,-1]]]){
 const t=scalar(variation),lhs=integral(pm(pm(n,j),ps(pd(t),'1/2'))),rhs=integral(ps(pm(force,t),-1));
 ck('source:variation_by_parts:'+z,lhs.eq(rhs));
 no('source:drop_gradient_force:'+z,!rhs.zero());
}
const f0=scalar([1,0,1]),deltaJ=pa(pm(f0,pd(j)),pm(pd(f0),j)),deltaN=br(f0,n),deltaNJ=pa(pm(deltaN,j),pm(n,deltaJ));
ck('source:covariant_scalar_pairing',peq(deltaNJ,pm(f0,pd(pm(n,j)))));
results.source={eta:'1',g,J:'4*x/3',N:'1+x',matter_cut_momentum_source:'2/3+4*x/3',variation_contract:'fixed N; theta_prime>0; compactly supported variations'};

function register(name,spec){const s=new M.Presentation(spec),audit=s.audit();ck(name+':confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name,spec:s.spec(),audit});systems.set(s.hash(),s);return{s,w:Object.fromEntries(spec.tokens.map(t=>[t,s.word([t])]))};}
function prove(name,lhs,rhs){const s=lhs.system,r=s.reduce(lhs.minus(rhs),{witness:true});ck(name,r.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(r.certificate));ck(name+':replay',s.replay(certificate));proofs.push({name,lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate});}
{
 const{s,w:{c,r}}=register('EMK_direction_block',{tokens:['c','r'],rules:[{id:'cc',lhs:['c','c'],rhs:[[[],1]],source:'EMK involution'},{id:'rr',lhs:['r','r'],rhs:[[[],-1]],source:'EMK circular generator'},{id:'rc',lhs:['r','c'],rhs:[[['c','r'],-1]],source:'EMK mixed parity'}]});
 prove('words:mixed_mass_parity',c.times(r).plus(r.times(c)),s.zero());
 prove('words:transport_principal_square',c.times(c),s.one());
 prove('words:mixed_involution',c.times(r).times(c).times(r),s.one());
}
{
 const tokens=['f','df','n','dn','m','dm','e','de','j','dj'],rules=[];
 for(let j=1;j<tokens.length;j++)for(let i=0;i<j;i++)rules.push({id:'central_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'Scalar local jets commute'});
 const{s,w:{f,df,n,dn,m,dm,e,de,j,dj}}=register('coframe_scalar_jets',{tokens,rules}),b=f.times(dn).minus(n.times(df)),deltaE=f.times(de).plus(df.times(e));
 prove('words:mixed_residual_cancel',f.times(dn.times(e).plus(n.times(de))).minus(b.times(e)).minus(n.times(deltaE)),s.zero());
 prove('words:normal_cross_mass_cancel',n.times(dm.times(e).plus(m.times(de))).plus(dn.times(m).times(e)).minus(m.times(dn.times(e).plus(n.times(de)))).minus(dm.times(n).times(e)),s.zero());
 prove('words:alpha_scalar',b.times(e).plus(n.times(deltaE)),f.times(dn.times(e).plus(n.times(de))));
 prove('words:metric_weight_two',e.times(deltaE).scale(2),f.times(e).times(de).scale(2).plus(df.times(e).times(e).scale(2)));
 prove('words:source_pairing_scalar',b.times(j).plus(n.times(f.times(dj).plus(df.times(j)))),f.times(dn.times(j).plus(n.times(dj))));
}
const chosen=proofs.find(x=>x.certificate.steps.length),bad=JSON.parse(JSON.stringify(chosen.certificate));bad.steps[0].coefficient=['991','0'];
no('integrity:corrupted_word_proof',(()=>{try{systems.get(chosen.certificate.presentation).replay(bad);return false;}catch{return true;}})());
const hash=x=>crypto.createHash('sha256').update(x).digest('hex'),payload={protocol:'LCF1_LOCAL_COFRAME_COVARIANCE_AND_CUT_SOURCE_V1',status:'PASS_NATIVE_MIXED_LOCAL_REPAIR_WITH_NORMAL_OBLIGATION',scope:'Native polynomial and word certification plus written local covariance/source proofs; no full gravitational constraint algebra, unbounded-domain evolution or local quantum field completion.',source_pins_sha256:hash(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),presentations,proofs,checks,negative_controls,results,summary:{exact_checks:checks.length,negative_controls:negative_controls.length,native_presentations:presentations.length,symbolic_certificates:proofs.length,serialized_replays:proofs.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',sha=hash(bytes),file=path.join(__dirname,'CERTIFICATE.json');
if(process.argv.includes('--write')){fs.writeFileSync(file,bytes);fs.writeFileSync(path.join(__dirname,'EXPECTED.sha256'),sha+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(file,'utf8'),bytes);assert.equal(fs.readFileSync(path.join(__dirname,'EXPECTED.sha256'),'utf8'),sha+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:sha}));
