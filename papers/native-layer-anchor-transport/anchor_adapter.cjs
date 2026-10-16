'use strict';
const {F}=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const O=require('../native-prime-observation-closure/observation_adapter.cjs');
const L=require('../native-infinite-layer-cut/layer_adapter.cjs');
const {Interval:I}=O;
const power=(x,n)=>{let p=new I(1);for(let i=0;i<n;i++)p=p.mul(x);return p;};
function hierarchy(raw){
 const c=raw.map(I.of);if(c.length<2)throw new RangeError('Two layers required');
 return c.slice(0,-1).map((z,n)=>z.mul(c[1].add(n)).sub(c[n+1].mul(n+1)));
}
function tail(K,j,z){
 z=F.of(z);if(!Number.isSafeInteger(K)||!Number.isSafeInteger(j)||j<0||j>K||z.le(0)&&!z.zero())throw new RangeError('Complete nonnegative tail indices and radius required');
 const start=K+1-j,ratio=z.mul(K+2).div(K+2-j);
 if(!ratio.le(1)||ratio.eq(1))throw new RangeError('Tail ratio must be below one');
 return new F(L.choose(K+1,j)).mul(z.pow(start)).div(new F(1).sub(ratio));
}
function errorBounds(q0,K,d,R='15/16',includeModel=true){
 q0=I.of(q0);d=F.of(d);R=F.of(R);
 if(!q0.positive()||d.le(0)||!d.le(R)||d.eq(R)||!R.eq('15/16'))throw new RangeError('Positive normalization and 0<d<R=15/16 required');
 const M=new F('221/228').div(q0.lo),errors=[];
 for(let j=0;j<=K;j++){
  const ratio=d.div(R).mul(K+2).div(K+2-j);
  if(!ratio.le(1)||ratio.eq(1))break;
  let e=M.div(R.pow(j)).mul(tail(K,j,d.div(R)));
  if(includeModel)e=e.add(tail(K,j,d));
  errors.push(new F(1).add(d).pow(j).mul(e));
 }
 return {K,d,R,M,errors,includeModel};
}
function corridor(raw,q0,K,m,d){
 if(!Number.isSafeInteger(K)||K<1||K%2!==1||raw.length<K+1||!Number.isSafeInteger(m)||m<1||2*m>K+1)throw new RangeError('Complete odd-degree source bank and finite target forms required');
 const model=L.forms(raw,(K+1)/2);
 if(!model.factors.A.positive||!model.factors.B.positive||!model.factors.H.positive)throw new RangeError('Complete positive finite surrogate required');
 const f=L.forms(raw,m),e=errorBounds(q0,K,d),one=new F(1);d=e.d;
 if(e.errors.length<2*m)throw new RangeError('Insufficient convergent derivative tails');
 const max=xs=>xs.reduce((a,b)=>a.le(b)?b:a,new F(0));
 const epsA=max(e.errors.slice(1,2*m)),epsB=max(Array.from({length:2*m-1},(_,j)=>e.errors[j].add(e.errors[j+1])));
 const tauA=epsA.mul(m).mul(one.add(d).pow(2*m)).div(one.sub(d).pow(2*m-1));
 const tauB=epsB.mul(m).mul(one.add(d).pow(2*m)).div(one.sub(d).pow(2*m-2));
 const A=O.intervalLDL(f.A.map((r,i)=>r.map((z,j)=>z.sub(i===j?tauA:0)))),B=O.intervalLDL(f.B.map((r,i)=>r.map((z,j)=>z.sub(i===j?tauB:0))));
 const normalizationLower=new F(1).div(one.add(d)).sub(e.errors[0]);
 return {K,m,d,bounds:e,epsA,epsB,tauA,tauB,normalizationLower,model_order:(K+1)/2,factors:{A,B},certified:A.positive&&B.positive&&!normalizationLower.le(0)};
}
function forecast(raw,q0,rho,degree){
 const c=raw.map(I.of),K=c.length-1;rho=F.of(rho);const delta=rho.sub(1),d=delta.abs();
 if(rho.le(0)||!Number.isSafeInteger(degree)||degree<0||degree>K)throw new RangeError('Positive anchor ratio and available degree required');
 if(d.zero())return {rho,raw:c.slice(0,degree+1),normalized:c.slice(0,degree+1),q:I.of(q0)};
 const e=errorBounds(q0,K,d,'15/16',false);if(e.errors.length<=degree)throw new RangeError('Requested forecast tail is not convergent');
 const values=Array.from({length:degree+1},(_,j)=>{
  let z=new I(0);for(let k=0;k+j<=K;k++)z=z.add(c[j+k].mul(L.choose(j+k,j)).mul(delta.neg().pow(k)));
  return z.mul(rho.pow(j)).add(new I(e.errors[j].neg(),e.errors[j]));
 });
 if(!values[0].positive())throw new RangeError('Forecast normalization not certified positive');
 return {rho,raw:values,normalized:values.map((z,j)=>j?z.div(values[0]):new I(1)),q:values[0].mul(q0),tails:e};
}
function exactMove(atoms,weights,rho){
 atoms=atoms.map(F.of);weights=weights.map(F.of);rho=F.of(rho);
 if(rho.le(0)||atoms.length!==weights.length||!atoms.length||atoms.some(x=>x.le(0)||new F(1).le(x))||weights.some(x=>x.le(0))||!weights.reduce((a,b)=>a.add(b),new F(0)).eq(1))throw new RangeError('Strict finite probability source and positive anchor ratio required');
 const denominators=atoms.map(y=>new F(1).add(rho.sub(1).mul(y))),normalization=weights.reduce((z,w,j)=>z.add(w.div(denominators[j])),new F(0));
 return {rho,atoms:atoms.map((y,j)=>rho.mul(y).div(denominators[j])),weights:weights.map((w,j)=>w.div(denominators[j]).div(normalization)),normalization};
}
function polynomialTransport(m,rho){
 rho=F.of(rho);if(!Number.isSafeInteger(m)||m<1||rho.le(0))throw new RangeError('Positive degree and anchor ratio required');
 return Array.from({length:m},(_,i)=>Array.from({length:m},(_,j)=>i<j?new F(0):new F(L.choose(m-1-j,i-j)).mul(rho.pow(j)).mul(rho.sub(1).pow(i-j))));
}
module.exports=Object.freeze({hierarchy,tail,errorBounds,corridor,forecast,exactMove,polynomialTransport,power});
