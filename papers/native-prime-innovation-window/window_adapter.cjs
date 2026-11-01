'use strict';
/** PIW-1 retained-source scalar window; no positivity correction is added. */
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const POC=require('../native-prime-observation-closure/observation_adapter.cjs');
const {F}=N,{Interval:I}=POC;
const dot=(a,b)=>a.reduce((s,x,i)=>s.add(I.of(x).mul(b[i])),new I(0));
function forward(f,b){if(!f.positive||b.length!==f.d.length)throw new RangeError('Certified factor and matching vector required');const w=[];for(let i=0;i<b.length;i++){let x=I.of(b[i]);for(let j=0;j<i;j++)x=x.sub(f.L[i][j].mul(w[j]));w.push(x);}return w;}
function solve(f,b){const w=forward(f,b).map((x,i)=>x.div(f.d[i])),v=Array(b.length);for(let i=b.length-1;i>=0;i--){let x=w[i];for(let j=i+1;j<b.length;j++)x=x.sub(f.L[j][i].mul(v[j]));v[i]=x;}return v;}
function inverseForm(f,b){return forward(f,b).reduce((s,x,i)=>s.add(x.square().div(f.d[i])),new I(0));}
function form(f,v){return f.d.reduce((s,d,i)=>{let w=new I(0);for(let j=i;j<v.length;j++)w=w.add(f.L[j][i].mul(v[j]));return s.add(d.mul(w.square()));},new I(0));}
function validate(bank,x){
 x=F.of(x);if(x.le(0))throw new RangeError('Positive new node required');const n=bank.nodes.length,ts=bank.nodes.map(F.of);
 if(!n||bank.values.length!==n||[bank.L,bank.J].some(M=>M.length!==n||M.some(r=>r.length!==n)))throw new TypeError('Pencil dimensions must agree');
 if(ts.some(t=>t.le(0))||new Set(ts.map(String)).size!==n||ts.some(t=>t.eq(x)))throw new RangeError('Distinct positive old nodes and a genuinely new node required');
 return {x,ts};
}
function coefficients(bank,x){
 const valid=validate(bank,x);x=valid.x;const ts=valid.ts,l=POC.intervalLDL(bank.L),j=POC.intervalLDL(bank.J);
 if(!l.positive||!j.positive)throw new RangeError('Both old source forms must be certified positive');
 const pencil=bank.J.map((row,i)=>row.map((z,k)=>I.of(z).add(I.of(bank.L[i][k]).mul(x)))),p=POC.intervalLDL(pencil);
 if(!p.positive)throw new RangeError('Response pencil not certified positive');
 const a=ts.map(t=>new F(1).div(x.sub(t))),v=solve(p,bank.values),r0=inverseForm(p,bank.values),d0=form(l,v).neg();
 const A=inverseForm(l,a),C=inverseForm(j,a),b=dot(a,v),kappa=A.add(C.mul(x)).mul(x);
 const determinant=f=>f.d.reduce((s,z)=>s.mul(z),new I(1)),nodeProduct=ts.reduce((s,t)=>s.mul(x.sub(t)),new F(1)),h=determinant(l).mul(nodeProduct).div(determinant(p));
 return {x,ts,a,v,r0,d0,A,C,b,kappa,valueUpper:r0.add(new I(1).div(kappa)),h,g:h.square(),factors:{L:l,J:j,P:p}};
}
function acceptance(bank,x,r,d){
 const c=coefficients(bank,x);x=c.x;r=I.of(r);d=I.of(d);const delta=r.sub(c.r0),deltaPrime=d.sub(c.d0),delta2=delta.square();
 const alpha=deltaPrime.neg().add(c.b.mul(delta).mul(2)).sub(c.A.mul(delta2));
 const beta=delta.add(deltaPrime.mul(x)).sub(c.b.mul(delta).mul(x.mul(2))).sub(c.C.mul(delta2).mul(x.pow(2)));
 const ell=c.ts.map((t,i)=>I.of(bank.values[i]).sub(r).div(x.sub(t))),j=c.ts.map((t,i)=>r.mul(x).sub(I.of(bank.values[i]).mul(t)).div(x.sub(t)));
 const rawAlpha=d.neg().sub(inverseForm(c.factors.L,ell)),rawBeta=r.add(d.mul(x)).sub(inverseForm(c.factors.J,j));
 const valueFraction=c.kappa.mul(delta),budget=delta.mul(new I(1).sub(valueFraction));
 const theta=budget.positive()?alpha.mul(x).div(budget):null;
 const normalizedResidual=delta.div(c.g),normalizedDerivative=deltaPrime.sub(c.b.mul(delta).mul(2)).div(c.g);
 const normalizedLower=normalizedResidual.div(x).neg().add(c.C.mul(c.g).mul(normalizedResidual.square()).mul(x));
 const normalizedUpper=c.A.mul(c.g).mul(normalizedResidual.square()).neg();
 return {x,r,d,coefficients:{r0:c.r0,d0:c.d0,A:c.A,C:c.C,b:c.b,kappa:c.kappa,valueUpper:c.valueUpper,h:c.h,g:c.g},delta,deltaPrime,alpha,beta,rawAlpha,rawBeta,valueFraction,budget,theta,normalizedResidual,normalizedDerivative,normalizedLower,normalizedUpper,accepted:rawAlpha.positive()&&rawBeta.positive()&&alpha.positive()&&beta.positive()};
}
function append(bank,jet){
 const {x,ts}=validate(bank,jet.t),r=I.of(jet.q),d=I.of(jet.derivative),ell=ts.map((t,i)=>I.of(bank.values[i]).sub(r).div(x.sub(t))),j=ts.map((t,i)=>r.mul(x).sub(I.of(bank.values[i]).mul(t)).div(x.sub(t)));
 return {nodes:[...bank.nodes,x.toString()],values:[...bank.values,r],L:[...bank.L.map((row,i)=>[...row,ell[i]]),[...ell,d.neg()]],J:[...bank.J.map((row,i)=>[...row,j[i]]),[...j,r.add(d.mul(x))]]};
}
function exactPencil(nodes,values,derivatives){
 const ts=nodes.map(F.of),u=values.map(F.of),d=derivatives.map(F.of);if(ts.length!==u.length||u.length!==d.length||!ts.length||new Set(ts.map(String)).size!==ts.length||ts.some(t=>t.le(0)))throw new TypeError('Compatible distinct positive response jets required');
 return {nodes:ts,values:u,L:ts.map((t,i)=>ts.map((x,j)=>i===j?d[i].neg():u[i].sub(u[j]).div(x.sub(t)))),J:ts.map((t,i)=>ts.map((x,j)=>i===j?u[i].add(t.mul(d[i])):x.mul(u[j]).sub(t.mul(u[i])).div(x.sub(t))))};
}
function determinant(M){
 const A=N.matrix(M).map(r=>r.map(x=>{if(!x.turn.zero())throw new TypeError('Real determinant required');return x.rad;})),n=A.length;if(A.some(r=>r.length!==n))throw new TypeError('Square determinant required');let out=new F(1);
 for(let k=0;k<n;k++){const j=A.findIndex((r,i)=>i>=k&&!r[k].zero());if(j<0)return new F(0);if(j!==k){[A[k],A[j]]=[A[j],A[k]];out=out.neg();}const p=A[k][k];out=out.mul(p);for(let i=k+1;i<n;i++){const q=A[i][k].div(p);for(let m=k;m<n;m++)A[i][m]=A[i][m].sub(q.mul(A[k][m]));}}
 return out;
}
function exactCoefficients(bank,x){
 const valid=validate(bank,x);x=valid.x;const ts=valid.ts,u=N.matrix(bank.values.map(z=>[z])),a=N.matrix(ts.map(t=>[new F(1).div(x.sub(t))])),p=N.add(bank.J,N.scale(bank.L,x)),inv=N.inverse(p),v=N.mul(inv,u),quad=(M,z)=>N.mul(N.dagger(z),N.mul(N.inverse(M),z))[0][0].rad;
 const r0=N.mul(N.dagger(u),v)[0][0].rad,d0=N.mul(N.dagger(v),N.mul(bank.L,v))[0][0].rad.neg(),A=quad(bank.L,a),C=quad(bank.J,a),b=N.mul(N.dagger(a),v)[0][0].rad,kappa=x.mul(A.add(x.mul(C))),h=determinant(bank.L).mul(ts.reduce((s,t)=>s.mul(x.sub(t)),new F(1))).div(determinant(p));
 return {x,a,v,p,r0,d0,A,C,b,kappa,h,g:h.pow(2),valueUpper:r0.add(new F(1).div(kappa)),traceInverseL:N.mul(inv,bank.L).reduce((s,r,i)=>s.add(r[i].rad),new F(0))};
}
module.exports=Object.freeze({forward,solve,inverseForm,form,coefficients,acceptance,append,exactPencil,determinant,exactCoefficients});
