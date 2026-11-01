'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs'),M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const f=N.F.of,checks=[],negative_controls=[],proofs=[],presentations=[],systems=new Map(),results={};
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const rejected=fn=>{try{fn();return false;}catch{return true;}};
// Rational first derivatives in canonical order (L,pi,phi,k).
const pb=(a,b)=>a[0].mul(b[1]).sub(a[1].mul(b[0])).add(a[2].mul(b[3])).sub(a[3].mul(b[2]));
const fixtures=[['6','3','1','0'],['2','-2','3/5','2/7'],['5/3','4/5','-2/3','1/2']];
results.constrained_jets=[];
for(const [j,raw]of fixtures.entries()){
 const [L,p,V,Vp]=raw.map(f),k=L.mul(V).div(p).neg(),C=L.mul(p).mul(k).add(L.pow(2).mul(V));
 const dC=[p.mul(k).add(L.mul(V).mul(2)),L.mul(k),L.pow(2).mul(Vp),L.mul(p)];
 const dK=[f(1).div(p),L.div(p.pow(2)).neg(),f(0),f(0)],dMass=[f(0),p.mul(2),V.mul(2),f(0)];
 ck('constraint:on_surface:'+j,C.zero());
 ck('clock:regular:'+j,!dC[3].zero());
 ck('observable:K_weak:'+j,pb(dK,dC).zero());
 ck('observable:M_weak:'+j,pb(dMass,dC).zero());
 ck('observable:canonical_K_M:'+j,pb(dK,dMass).eq(2));
 const Ltau=L.mul(V).div(p.pow(2)).neg(),ptau=V.div(p).neg(),ktau=L.mul(Vp).div(p).neg();
 ck('invariant:ratio_derivative:'+j,Ltau.mul(p).sub(L.mul(ptau)).zero());
 ck('invariant:mass_derivative:'+j,p.mul(ptau).mul(2).add(V.mul(2)).zero());
 const kFromSolved=Ltau.mul(V).div(p).neg().sub(L.mul(Vp).div(p)).add(L.mul(V).mul(ptau).div(p.pow(2)));
 ck('constraint:remaining_momentum_equation:'+j,kFromSolved.eq(ktau));
 for(const lapse of ['1','2/3','-5/2']){
  const n=f(lapse),clock=n.mul(L).mul(p),Ldot=n.mul(L).mul(k),pdot=n.mul(dC[0]).neg(),kdot=n.mul(L.pow(2)).mul(Vp).neg();
  ck('clock:lapse_cancels_L:'+j+':'+lapse,Ldot.div(clock).eq(Ltau));
  ck('clock:lapse_cancels_pi:'+j+':'+lapse,pdot.div(clock).eq(ptau));
  ck('clock:lapse_cancels_k:'+j+':'+lapse,kdot.div(clock).eq(ktau));
 }
 // Nonconstraint jet verifies proportionality, not merely zero on chosen states.
 const offk=k.add(1),offC=L.mul(p).mul(offk).add(L.pow(2).mul(V)),offdC=[p.mul(offk).add(L.mul(V).mul(2)),L.mul(offk),L.pow(2).mul(Vp),L.mul(p)];
 ck('observable:off_surface_K_factor:'+j,pb(dK,offdC).eq(offC.mul(2).div(p.pow(2))));
 ck('observable:off_surface_M_factor:'+j,pb(dMass,offdC).eq(offC.mul(p).mul(-2).div(L)));
 no('observable:K_strongly_commutes:'+j,!pb(dK,offdC).zero());
 results.constrained_jets.push({L,p,V,Vprime:Vp,k,C,L_tau:Ltau,pi_tau:ptau,k_tau:ktau});
}
// Polynomial coefficients are exact native rationals, in ascending clock powers.
const deriv=a=>a.length===1?[f(0)]:a.slice(1).map((c,j)=>c.mul(j+1));
const evalp=(a,t)=>a.reduceRight((s,c)=>s.mul(t).add(c),f(0));
const models=[{name:'A',V:[1],W:[0,1],p2:[9,-2]},{name:'B',V:[0,'1/2'],W:[0,0,'1/4'],p2:[9,0,'-1/2']}].map(o=>({...o,V:o.V.map(f),W:o.W.map(f),p2:o.p2.map(f)}));
results.event_models=[];
for(const m of models){
 const recV=deriv(m.p2).map(c=>c.mul('-1/2')),dW=deriv(m.W);
 ck('inverse:polynomial_potential:'+m.name,JSON.stringify(recV)===JSON.stringify(m.V));
 ck('inverse:primitive:'+m.name,JSON.stringify(dW)===JSON.stringify(m.V));
 for(const t of [0,2,4])ck('observable:M_on_curve:'+m.name+':'+t,evalp(m.p2,f(t)).add(evalp(m.W,f(t)).mul(2)).eq(9));
 ck('events:initial_pi_squared:'+m.name,evalp(m.p2,f(0)).eq(9));
 ck('events:final_pi_squared:'+m.name,evalp(m.p2,f(4)).eq(1));
 ck('inverse:endpoint_integral:'+m.name,evalp(m.W,f(4)).sub(evalp(m.W,f(0))).eq(f(9).sub(1).div(2)));
 results.event_models.push({name:m.name,V:m.V,W:m.W,pi_squared:m.p2,initial:{phi:'0',L:'6',pi:'3'},final:{phi:'4',L:'2',pi:'1'},interior:{phi:'2',pi_squared:evalp(m.p2,f(2)),L_squared:evalp(m.p2,f(2)).mul(4)}});
}
no('inverse:endpoints_select_interior',!evalp(models[0].p2,f(2)).eq(evalp(models[1].p2,f(2))));
no('inverse:endpoints_select_initial_k',!evalp(models[0].V,f(0)).eq(evalp(models[1].V,f(0))));
// Finite rational event on Model A: tau=5/2, pi=2, L=4; no square root approximated.
const tau=f('5/2'),L=f(4),p=f(2),L0=f(6),p0=f(3),a=L.div(L0),b=p.div(p0);
ck('event:exact_square_branch',p.pow(2).eq(evalp(models[0].p2,tau)));
ck('event:normalized_diagonal',a.eq(b));
ck('event:normalized_invariant',a.div(b).eq(1));
// Complete-observable gradients at the initial point L=6, pi=3, phi=0, k=-2,
// evaluated at clock target tau=5/2 (Pi_tau=2, L_tau=4).
const dPiTarget=['0','3/2','1/2','0'].map(f),dLTarget=['2/3','5/3','1','0'].map(f),dCInitial=['6','-12','0','18'].map(f);
ck('event:complete_observable_canonical',pb(dLTarget,dPiTarget).eq(1));
ck('event:complete_Pi_weak',pb(dPiTarget,dCInitial).zero());
ck('event:complete_L_weak',pb(dLTarget,dCInitial).zero());
no('calibration:independent_axis_scales_preserve_slope_one',!a.mul(3).div(b.mul(2)).eq(1));
no('rate:diagonal_implies_constant_rate',!f(-1).div(9).eq(f(-1).div(4)));
ck('baseline:zero_potential_L_static',L.mul(0).div(p.pow(2)).zero());
ck('baseline:zero_potential_pi_static',f(0).div(p).zero());
ck('baseline:zero_potential_clock_runs',L.mul(p).eq(8));
no('clock:zero_momentum_division',rejected(()=>f(1).div(0)));
ck('clock:zero_momentum_can_satisfy_constraint',f(1).mul(0).mul(0).add(f(1).pow(2).mul(0)).zero());
no('clock:constraint_implies_clock_regular',f(1).mul(0).zero());
no('topology:fixed_length_allows_nonzero_V',!L.mul(1).div(p.pow(2)).neg().zero());
// Closed noncommutative algebra and state-dependent uncertainty, all exact.
const I=N.identity(2),X=N.matrix([[0,1],[1,0]]),Y=N.matrix([[0,[0,-1]],[[0,1],0]]),Z=N.matrix([[1,0],[0,-1]]);
const basis=[I,X,Y,Z],flatten=N.flatten;
ck('closure:matrix_basis_dimension',N.rank(basis.map(flatten))===4);
ck('closure:X_Z_generate_full_algebra',N.operatorAlgebra([X,Z]).dimension===4);
for(let i=0;i<4;i++)for(let j=0;j<4;j++)ck('closure:product:'+i+':'+j,N.rank([...basis.map(flatten),flatten(N.mul(basis[i],basis[j]))])===4);
const rhoY=N.scale(N.add(I,Y),'1/2'),rhoZ=N.scale(N.add(I,Z),'1/2');
const tr=m=>m.reduce((s,row,j)=>s.add(row[j]),N.ZERO),expect=(rho,A)=>tr(N.mul(rho,A));
const variance=(rho,A)=>expect(rho,N.mul(A,A)).sub(expect(rho,A).mul(expect(rho,A)));
ck('closure:rho_Y_projector',N.equal(N.mul(rhoY,rhoY),rhoY));
ck('uncertainty:var_X',variance(rhoY,X).eq(1));
ck('uncertainty:var_Z',variance(rhoY,Z).eq(1));
ck('uncertainty:nonzero_commutator_bound',expect(rhoY,N.commutator(X,Z)).norm2().div(4).eq(1));
no('closure:noncommutativity_forbids_closure',!N.isZero(N.commutator(X,Z)));
ck('uncertainty:Z_eigenstate_zero_variance',variance(rhoZ,Z).zero());
no('closure:closure_forces_nonzero_variance',variance(rhoZ,Z).zero());
function register(name,tokens){const rules=[];for(let j=1;j<tokens.length;j++)for(let i=0;i<j;i++)rules.push({id:'commute_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'Commuting scalar symbols in the classical homogeneous canonical sector'});const s=new M.Presentation({tokens,rules}),audit=s.audit();ck(name+':confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name,spec:s.spec(),audit});systems.set(s.hash(),s);return{s,w:Object.fromEntries(tokens.map(t=>[t,s.word([t])]))};}
function prove(name,lhs,rhs){const s=lhs.system,r=s.reduce(lhs.minus(rhs),{witness:true});ck(name,r.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(r.certificate));ck(name+':replay',s.replay(certificate));proofs.push({name,lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate});}
{
 const{s,w:{l,p,k,v,dv}}=register('homogeneous_clock',['l','p','k','v','dv']);
 const C=l.times(p).times(k).plus(l.times(l).times(v)),CL=p.times(k).plus(l.times(v).scale(2)),Cp=l.times(k),Ck=l.times(p);
 prove('words:K_weak_factor',p.times(Cp).plus(l.times(CL)),C.scale(2));
 prove('words:M_weak_factor',l.times(p.times(CL).scale(-2).plus(v.times(Ck).scale(2))),p.times(C).scale(-2));
 prove('words:ratio_conservation',p.times(l).times(v).scale(-1).plus(l.times(v).times(p)),s.zero());
 prove('words:remaining_momentum_equation',l.times(v).times(v).minus(l.times(dv).times(p).times(p)).minus(v.times(l).times(v)),l.times(dv).times(p).times(p).scale(-1));
}
const chosen=proofs.find(q=>q.certificate.steps.length),bad=JSON.parse(JSON.stringify(chosen.certificate));bad.steps[0].coefficient=['991','0'];
no('integrity:corrupted_native_proof',rejected(()=>systems.get(chosen.certificate.presentation).replay(bad)));
results.scope={classical:true,homogeneous:true,matter:'zero',clock:'phi',domain:'L>0; pi nonzero; fixed sign branch; radicand>0',global_single_angle_lift_assumed:false,full_local_quantization:false,continuum_proofs:'THEOREM.md'};
const hash=x=>crypto.createHash('sha256').update(x).digest('hex'),payload={protocol:'REC1_RELATIONAL_EVENT_CLOCK_V1',status:'PASS_NATIVE_RELATIONAL_EVENT_CLOCK',scope:'Selected homogeneous classical LNC reduction; finite exact checks and symbolic denominator-cleared identities. No full local quantum gravity claim.',source_pins_sha256:hash(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),presentations,proofs,checks,negative_controls,results,summary:{exact_checks:checks.length,negative_controls:negative_controls.length,native_presentations:presentations.length,symbolic_certificates:proofs.length,serialized_replays:proofs.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',sha=hash(bytes),file=path.join(__dirname,'CERTIFICATE.json');
if(process.argv.includes('--write')){fs.writeFileSync(file,bytes);fs.writeFileSync(path.join(__dirname,'EXPECTED.sha256'),sha+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(file,'utf8'),bytes);assert.equal(fs.readFileSync(path.join(__dirname,'EXPECTED.sha256'),'utf8'),sha+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:sha}));
