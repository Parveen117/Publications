'use strict';
const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict');
const N=require('./native_engine/native_operator.cjs'),f=N.F.of;
const I=N.identity(2),C=N.matrix([[1,0],[0,-1]]),R=N.matrix([[0,-1],[1,0]]),X=N.matrix([[0,1],[1,0]]);
const mm=(...x)=>x.reduce(N.mul),checks=[],controls=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);controls.push({name,rejected:true});};
const eq=N.equal,zero=N.isZero,scale=N.scale;
ck('cut_square',eq(mm(C,C),I));
ck('turn_square',eq(mm(R,R),scale(I,-1)));
ck('anticommutation',zero(N.add(mm(C,R),mm(R,C))));
ck('cut_dagger',eq(N.dagger(C),C));
ck('turn_anti_dagger',eq(N.dagger(R),scale(R,-1)));
const fixtures=[['0','0'],['0','3/2'],['2/3','0'],['3/4','5/7'],['-4/5','7/3'],['2','-3']];
for(const [qs,gs]of fixtures){
 const q=f(qs),g=f(gs),A=scale(C,N.IOTA.mul(q)),B=scale(R,g),G=N.add(A,B),sum=q.pow(2).add(g.pow(2));
 ck(`generator_square:${qs}:${gs}`,eq(mm(G,G),scale(I,sum.neg())));
 ck(`generator_anti_dagger:${qs}:${gs}`,eq(N.dagger(G),scale(G,-1)));
 const comm=N.commutator(A,B);
 ck(`commutator_norm_square:${qs}:${gs}`,eq(mm(N.dagger(comm),comm),scale(I,q.pow(2).mul(g.pow(2)).mul(4))));
 const correction=q.pow(2).mul(g.pow(2)).div(-3);
 const lhs=correction.div(-2).add(sum.pow(2).div(24));
 const rhs=q.pow(4).add(g.pow(4)).div(24).add(q.pow(2).mul(g.pow(2)).div(4));
 ck(`dispersion_fourth_order:${qs}:${gs}`,lhs.eq(rhs));
}
const ch=f('5/3'),sh=f('4/3'),J=C,L=N.matrix([[ch,sh],[sh,ch]]);
ck('rational_boost_preserves_form',eq(mm(N.dagger(L),J,L),J));
ck('boosted_rest_example',f('10/3').pow(2).sub(f('8/3').pow(2)).eq(4));
const boostGenerator=X;
ck('boost_generator_compatibility',zero(N.add(mm(N.dagger(boostGenerator),J),mm(J,boostGenerator))));
no('dilation_does_not_preserve_form',!eq(mm(scale(I,2),J,scale(I,2)),J));
no('incompatible_generator',!zero(N.add(mm(I,J),mm(J,I))));
// Native real Cl(3,1) coefficient module used by MG/SM.
const K=C,Lcut=mm(K,R);
const gamma=[N.kron(R,I),N.kron(K,K),N.kron(K,Lcut),N.kron(Lcut,I)],I4=N.identity(4),eta=[-1,1,1,1];
for(let a=0;a<4;a++)for(let b=a;b<4;b++)
 ck(`clifford:${a}:${b}`,eq(N.add(mm(gamma[a],gamma[b]),mm(gamma[b],gamma[a])),scale(I4,a===b?2*eta[a]:0)));
for(const vec of [[1,1,0,0],[2,0,0,0],['3/2','1/2','2/3','-4/5'],[1,0,1,0]]){
 let S=N.zeros(4),metric=f(0);
 vec.forEach((v,a)=>{S=N.add(S,scale(gamma[a],v));metric=metric.add(f(v).pow(2).mul(eta[a]));});
 ck('symbol_square:'+vec.join(','),eq(mm(S,S),scale(I4,metric)));
 ck('symbol_singular_iff_null:'+vec.join(','),(N.rank(S)<4)===metric.eq(0));
}
// Unequal address shifts preserve exponential/unitary structure but change the cone.
const q=f(1),g=f('1/2'),v=f(2),Gv=N.add(scale(C,N.IOTA.mul(v)),scale(R,g));
ck('unequal_speed_actual_square',eq(mm(Gv,Gv),scale(I,v.pow(2).add(g.pow(2)).neg())));
no('unequal_speed_common_shell',!eq(mm(Gv,Gv),scale(I,q.pow(2).add(g.pow(2)).neg())));
// MGH local coupling: direction tensor matter tensor geometry; no spatial derivatives.
const P=N.scale(N.sub(I,C),'1/2'),Q=N.scale(N.sub(I,X),'1/2');
const V=N.add(scale(N.kron(R,N.identity(4)),'3/5'),scale(N.kron(I,N.kron(P,Q)),N.IOTA.mul('7/4')));
ck('local_quantum_geometry_anti_dagger',eq(N.dagger(V),scale(V,-1)));
ck('local_quantum_geometry_preserves_source_moment',zero(N.commutator(V,N.kron(I,N.kron(P,I)))));
// An exact local-coin lattice support fixture; this supplements the written induction proof.
let sites=new Map([[0,N.matrix([[1],[0]])]]);
const coin=N.add(scale(I,'3/5'),scale(R,'4/5'));
for(let step=1;step<=8;step++){
 const next=new Map();
 for(const [j,psi]of sites){const turned=mm(coin,psi);for(const [dest,value]of [[j-1,N.matrix([[turned[0][0]],[0]])],[j+1,N.matrix([[0],[turned[1][0]]])]])next.set(dest,N.add(next.get(dest)||N.zeros(2,1),value));}
 sites=next;
 ck('finite_support_step:'+step,[...sites.keys()].every(j=>Math.abs(j)<=step));
 ck('lattice_norm_step:'+step,[...sites.values()].reduce((s,z)=>s.add(N.energy(z)),f(0)).eq(1));
}
const numerical=[0.2,0.1,0.05,0.025].map(e=>{
 const q=1.2,g=0.7,omega=Math.acos(Math.cos(e*g)*Math.cos(e*q))/e,base=q*q+g*g;
 const corrected=base-e*e*g*g*q*q/3;
 return{epsilon:e,omega_squared:omega*omega,continuum_squared:base,error:omega*omega-base,corrected_residual:omega*omega-corrected};
});
const result={protocol:'CCL1_SUPPLEMENTARY_CHECKS_V1',scope:'Native exact finite fixtures plus separately labeled floating-point illustrations. Written continuum proofs are not proof-assistant certified.',exact_checks:checks.length,negative_controls:controls.length,checks,negative_control_details:controls,numerical_illustration:numerical,passed:true};
const target=path.join(__dirname,'CHECKS.json');
if(process.argv.includes('--write'))fs.writeFileSync(target,JSON.stringify(result,null,2)+'\n');
else{const saved=JSON.parse(fs.readFileSync(target,'utf8'));assert.deepEqual(result,saved);}
console.log(JSON.stringify({passed:true,exact_checks:checks.length,negative_controls:controls.length,numerical_illustration:numerical}));
