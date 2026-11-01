'use strict';
// AGD-1: all scientific arithmetic uses the unchanged native cut engine.
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs'),M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const f=N.F.of,I=N.identity(2),Z=N.matrix([[1,0],[0,-1]]),X=N.matrix([[0,1],[1,0]]),R=N.matrix([[0,-1],[1,0]]),J=N.scale(X,N.IOTA),L=N.mul(R,J);
const mm=(...a)=>a.reduce(N.mul),sc=N.scale,eq=N.equal,zero=N.isZero,cm=N.commutator;
const tr=a=>a.reduce((s,r,i)=>s.add(r[i]),N.ZERO),real=z=>{assert.ok(z.turn.zero());return z.rad;},expect=(rho,A)=>real(tr(mm(rho,A)));
const checks=[],controls=[],proofs=[],presentations=[],systems=new Map(),results={};
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);controls.push({name,rejected:true});};
const B=[R,J,L],weights=[1,1,2],pairWeights=[[0,1,1],[0,2,2],[1,2,2]],circle=[];
for(const[cs,ss]of [['1','0'],['-1','0'],['0','1'],['3/5','4/5'],['5/13','12/13'],['-3/5','4/5']]){
 const c=f(cs),s=f(ss),P=sc(N.add(I,N.add(sc(Z,c),sc(X,s))),'1/2'),dP=sc(N.add(sc(Z,s.neg()),sc(X,c)),'1/2');
 ck('action:cut:'+cs+':'+ss,eq(mm(P,P),P)&&eq(N.dagger(P),P)&&tr(P).eq(1));
 ck('action:inertia:'+cs+':'+ss,N.energy(dP).eq('1/2'));
 const tangents=B.map(b=>cm(b,P));
 const phiI=tangents.reduce((v,t,i)=>v.add(N.energy(t).mul(weights[i])),f(0));
 const phiC=pairWeights.reduce((v,[i,j,h])=>v.add(N.energy(mm(P,cm(tangents[i],tangents[j]),P)).mul(h).div(2)),f(0));
 ck('action:information:'+cs+':'+ss,phiI.eq(f(6).sub(c.pow(2).mul(2))));
 ck('action:curvature:'+cs+':'+ss,phiC.eq(f(4).sub(c.pow(2).mul(2))));
 ck('action:potential:'+cs+':'+ss,phiI.add(phiC).eq(f(10).sub(c.pow(2).mul(4))));
 const Q1=sc(N.sub(I,Z),'1/2'),overlap=expect(P,Q1);
 ck('lapse:cut_overlap:'+cs+':'+ss,overlap.eq(f(1).sub(c).div(2)));
 circle.push({cosine:c,sine:s,kinetic_metric:'1/2',information_potential:phiI,curvature_potential:phiC,overlap});
}
results.action={inertia:'1/2',potential:'8-2*cos(2*theta)',circle_fixtures:circle};

function rotor(radius){const count=2*radius+1,id=N.identity(count),u=N.zeros(count),number=N.zeros(count),parity=N.zeros(count);
 for(let j=0;j<count;j++){number[j][j]=N.Cut.of(j-radius);parity[j][j]=N.Cut.of((j-radius)%2===0?1:-1);if(j+1<count)u[j+1][j]=N.ONE;}
 const ud=N.dagger(u),cos=sc(N.add(u,ud),'1/2'),sin=sc(N.sub(u,ud),N.IOTA.mul('-1/2'));
 const hg=N.sub(N.add(mm(number,number),sc(id,8)),N.add(mm(u,u),mm(ud,ud)));
 return{radius,count,id,u,ud,number,parity,cos,sin,hg};
}
const small=rotor(3),top=N.zeros(7),bottom=N.zeros(7);top[6][6]=N.ONE;bottom[0][0]=N.ONE;
ck('finite:open_shift_commutator',eq(cm(small.number,small.u),small.u));
ck('finite:upper_boundary_defect',eq(mm(small.ud,small.u),N.sub(small.id,top)));
ck('finite:lower_boundary_defect',eq(mm(small.u,small.ud),N.sub(small.id,bottom)));
no('finite:open_shift_is_unitary',!eq(mm(small.ud,small.u),small.id));
const cyc=N.add(small.u,N.matrix(Array.from({length:7},(_,i)=>Array.from({length:7},(_,j)=>i===0&&j===6?1:0)))),edge=N.zeros(7);edge[0][6]=N.Cut.of(-7);
ck('finite:cyclic_shift_unitarity',eq(mm(N.dagger(cyc),cyc),small.id));
ck('finite:cyclic_commutator_defect',eq(N.sub(cm(small.number,cyc),cyc),edge));
no('finite:cyclic_shift_exact_canonical_relation',!eq(cm(small.number,cyc),cyc));

const a=f('3/5'),q=rotor(5),v=N.sub(sc(q.id,f(1).add(a.div(2))),sc(q.cos,a.div(2))),k=f(2),hm=sc(v,k.neg()),H=N.add(q.hg,hm);
const ket=N.zeros(q.count,1);ket[q.radius][0]=N.ONE;ket[q.radius+1][0]=N.IOTA;
const rho=sc(mm(ket,N.dagger(ket)),'1/2');
ck('source:normalized_pure_preparation',tr(rho).eq(1)&&eq(mm(rho,rho),rho));
ck('source:joint_Hermitian',eq(N.dagger(H),H));
const sin2=sc(N.sub(mm(q.u,q.u),mm(q.ud,q.ud)),N.IOTA.mul('-1/2'));
const force=sc(cm(H,q.number),N.IOTA),expectedForce=N.add(sc(sin2,-4),sc(q.sin,k.mul(a).div(2)));
ck('source:force_operator_identity',eq(force,expectedForce));
ck('source:first_sine',expect(rho,q.sin).eq('-1/2'));ck('source:second_sine',expect(rho,sin2).eq(0));
ck('source:momentum_rate',expect(rho,force).eq('-3/10'));
const geoRate=sc(cm(H,q.hg),N.IOTA),matterRate=sc(cm(H,hm),N.IOTA);
ck('source:geometry_energy_rate',expect(rho,geoRate).eq('-3/10'));
ck('source:matter_energy_rate',expect(rho,matterRate).eq('3/10'));
ck('source:total_energy_conservation',zero(N.add(geoRate,matterRate)));
const vMean=expect(rho,v),vVar=expect(rho,mm(v,v)).sub(vMean.pow(2)),pMean=expect(rho,q.number),pVar=expect(rho,mm(q.number,q.number)).sub(pMean.pow(2)),vtheta=sc(q.sin,a.div(2)),vthetaMean=expect(rho,vtheta),uncertaintyLower=vthetaMean.pow(2).div(4);
ck('variance:mean_lapse',vMean.eq('13/10'));ck('variance:lapse',vVar.eq('9/200'));ck('variance:momentum',pVar.eq('1/4'));
ck('variance:source_derivative',vthetaMean.eq('-3/20'));ck('variance:uncertainty_product',pVar.mul(vVar).eq('9/800')&&uncertaintyLower.eq('9/1600')&&uncertaintyLower.le(pVar.mul(vVar)));
results.source={eta:'1',k:'2',direction:'1',a,initial_density:rho,geometry_H:q.hg,matter_H:hm,force_operator:force,momentum_rate:'-3/10',geometry_energy_rate:'-3/10',matter_energy_rate:'3/10',mean_lapse:vMean,lapse_variance:vVar,momentum_variance:pVar,uncertainty_lower:uncertaintyLower};

// Uncoupled parity is an exact symmetry; it is not a preserved coupled constraint.
ck('constraint:geometry_parity_symmetry',zero(cm(q.parity,q.hg)));
no('constraint:coupled_parity_preserved',!zero(cm(q.parity,H)));
const even0=N.zeros(q.count,1);even0[q.radius][0]=N.ONE;
const odd=sc(N.sub(q.id,q.parity),'1/2'),leak=mm(odd,H,even0);
ck('constraint:exact_parity_leak',N.energy(leak).eq('9/50'));
results.constraints={uncoupled_parity_preserved:true,coupled_parity_preserved:false,even_initial_leak_squared:'9/50',clock_constraint:'global p_s+H only; not local gravitational constraints'};

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
const Dop=f=>[f,ps(pd(f),'1/2')],Aop=(n,g)=>[left(Z,n),pa(ps(left(Z,pd(n)),'1/2'),ps(left(R,n),g))];
const lift=op=>[[zmat()],...op],opEq=(a,b)=>a.every((p,i)=>peq(p,b[i]));
for(const [ff,nn,hh]of [[[0,1],[1],[1,0,1]],[[1,0,1],[0,1],[1,1]],[[1,2,1],[1,-1],[0,0,1]]]){
 const f0=scalar(ff),n0=scalar(nn),h0=scalar(hh),br=(f,g)=>pa(pm(f,pd(g)),ps(pm(g,pd(f)),-1));
 ck('local_constraint:DD:'+ff.join(','),opEq(commOp(Dop(f0),Dop(h0)),lift(Dop(br(f0,h0)))));
 ck('local_constraint:AA:'+ff.join(','),opEq(commOp(Aop(f0,'2/3'),Aop(n0,'2/3')),lift(Dop(br(f0,n0)))));
 const expected=lift(Aop(br(f0,n0),'2/3'));expected[2]=pa(expected[2],ps(left(R,pm(n0,pd(f0))),'2/3'));
 ck('local_constraint:DA_with_residual:'+ff.join(','),opEq(commOp(Dop(f0),Aop(n0,'2/3')),expected));
 no('local_constraint:massive_naive_closure:'+ff.join(','),!opEq(commOp(Dop(f0),Aop(n0,'2/3')),lift(Aop(br(f0,n0),'2/3'))));
 ck('local_constraint:massless_special_case:'+ff.join(','),opEq(commOp(Dop(f0),Aop(n0,0)),lift(Aop(br(f0,n0),0))));
}

// Band-path finite-section witnesses: compare powers on a safely larger carrier.
const big=rotor(9),Vbig=N.sub(sc(big.id,'13/10'),sc(big.cos,'3/10')),Hbig=N.sub(big.hg,Vbig),cutoff=5;
const mask=N.matrix(Array.from({length:big.count},(_,i)=>Array.from({length:big.count},(_,j)=>i===j&&Math.abs(i-big.radius)<=cutoff?1:0))),HL=mm(mask,Hbig,mask),seed=N.zeros(big.count,1);seed[big.radius+1][0]=N.ONE;
let full=seed,finite=seed;
for(let order=0;order<=3;order++){
 if(order<3)ck('truncation:matching_power_before_exit:'+order,eq(full,finite));else no('truncation:claim_no_boundary_error',!eq(full,finite));
 full=mm(Hbig,full);finite=mm(HL,finite);
}
const t=f('23/100'),bound=t.pow(6).mul(6).div(720);
ck('truncation:rational_error_ceiling',bound.eq('148035889/120000000000000')&&bound.le('1/800000'));
results.truncation={eta:'1',g:'0',k:'1',a:'3/5',s:'1/10',initial_radius:1,retained_radius:11,minimum_exit_insertions:6,bounded_offdiagonal_norm:'23/10',t,certified_norm_error:bound,scope:'Homogeneous fixed spatial mode; written Dyson-band proof supplies the infinite comparison.'};

// Native serialized algebraic proofs, separate from the analytic domain proofs.
function system(name,tokens,rules=[],commute=false){if(commute)for(let j=1;j<tokens.length;j++)for(let i=0;i<j;i++)rules.push({id:'central_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'Declared central coefficients'});const s=new M.Presentation({tokens,rules}),audit=s.audit();ck(name+':confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name,spec:s.spec(),audit});systems.set(s.hash(),s);return{s,w:Object.fromEntries(tokens.map(t=>[t,s.word([t])]))};}
function prove(name,lhs,rhs){const s=lhs.system,r=s.reduce(lhs.minus(rhs),{witness:true});ck(name,r.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(r.certificate));ck(name+':replay',s.replay(certificate));proofs.push({name,lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate});}
// The circle algebra uses u inverse v; p is its unbounded angular momentum word.
{
 const{s,w:{u,v,h,p}}=system('circle_momentum',['u','v','h','p'],[
  {id:'uv',lhs:['u','v'],rhs:[[[],1]],source:'Unitary circle shift'},
  {id:'vu',lhs:['v','u'],rhs:[[[],1]],source:'Unitary circle shift'},
  {id:'hu',lhs:['h','u'],rhs:[[['u','h'],1]],source:'Central action unit'},
  {id:'hv',lhs:['h','v'],rhs:[[['v','h'],1]],source:'Central action unit'},
  {id:'ph',lhs:['p','h'],rhs:[[['h','p'],1]],source:'Central action unit'},
  {id:'pu',lhs:['p','u'],rhs:[[['u','p'],1],[['u','h'],1]],source:'p u = u(p+eta)'},
  {id:'pv',lhs:['p','v'],rhs:[[['v','p'],1],[['v','h'],-1]],source:'p u_dagger = u_dagger(p-eta)'}
 ]);
 const one=s.one(),cos=u.plus(v).scale('1/2'),potential=one.scale(10).minus(cos.times(cos).scale(4)),Hg=p.times(p).plus(potential);
 prove('action:circle_potential',potential,one.scale(8).minus(u.times(u)).minus(v.times(v)));
 prove('geometry:momentum_source',Hg.times(p).minus(p.times(Hg)),h.times(u.times(u).minus(v.times(v))).scale(2));
 prove('geometry:kinetic_shift_commutator',p.times(p).times(u).minus(u.times(p).times(p)),u.times(h).times(p.scale(2).plus(h)));
 prove('geometry:inverse_shift_commutator',p.times(p).times(v).minus(v.times(p).times(p)),v.times(h).times(p.scale(-2).plus(h)));
 const lapse=one.scale('13/10').minus(u.plus(v).scale('3/20'));
 prove('source:lapse_momentum_commutator',lapse.times(p).minus(p.times(lapse)),h.times(u.minus(v)).scale('3/20'));
}
{
 const{s,w:{c,sn}}=system('cut_circle_metric',['c','sn'],[{id:'circle',lhs:['sn','sn'],rhs:[[[],1],[['c','c'],-1]],source:'Real cut circle'}],true);
 prove('action:kinetic_metric',c.times(c).plus(sn.times(sn)).scale('1/2'),s.one().scale('1/2'));
 prove('action:potential_split',s.one().scale(6).minus(c.times(c).scale(2)).plus(s.one().scale(4).minus(c.times(c).scale(2))),s.one().scale(10).minus(c.times(c).scale(4)));
}
{
 const{s,w:{p}}=system('legendre',['p']);prove('action:legendre_inertia',p.times(p.scale(2)).minus(p.scale(2).times(p.scale(2)).scale('1/4')),p.times(p));
}
const chosen=proofs.find(p=>p.certificate.steps.length),bad=JSON.parse(JSON.stringify(chosen.certificate));bad.steps[0].coefficient=['991','0'];
no('integrity:corrupted_word_proof',(()=>{try{systems.get(chosen.certificate.presentation).replay(bad);return false;}catch{return true;}})());
const hash=x=>crypto.createHash('sha256').update(x).digest('hex'),payload={protocol:'AGD1_NATIVE_ACTION_GEOMETRY_DYNAMICS_V1',status:'PASS_ACTION_DERIVATION_NATIVE_PROOFS_AND_FINITE_CONTROLS',scope:'Inherited cut-action reduction under a declared periodic canonical quantization contract. Written full-domain self-adjointness, locality and finite-section error proofs. Exact native symbolic/finite checks and explicit constraint rejections; not a completed physical quantum-gravity theory.',source_pins_sha256:hash(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),presentations,proofs,checks,negative_controls:controls,results,summary:{exact_checks:checks.length,negative_controls:controls.length,native_presentations:presentations.length,symbolic_certificates:proofs.length,serialized_replays:proofs.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',sha=hash(bytes),file=path.join(__dirname,'CERTIFICATE.json');
if(process.argv.includes('--write')){fs.writeFileSync(file,bytes);fs.writeFileSync(path.join(__dirname,'EXPECTED.sha256'),sha+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(file,'utf8'),bytes);assert.equal(fs.readFileSync(path.join(__dirname,'EXPECTED.sha256'),'utf8'),sha+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:sha,source_rate:'-3/10',parity_leak:'9/50',finite_section_error:bound}));
