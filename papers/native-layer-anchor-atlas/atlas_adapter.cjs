'use strict';
const {F}=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const O=require('../native-prime-observation-closure/observation_adapter.cjs');
const L=require('../native-infinite-layer-cut/layer_adapter.cjs');
const {Interval:I}=O,R=new F('15/16');
const positive=x=>!F.of(x).le(0),max=xs=>xs.reduce((a,b)=>a.le(b)?b:a,new F(0));
function exactTail(K,j,z){
 z=F.of(z);
 if(!Number.isSafeInteger(K)||K<0||!Number.isSafeInteger(j)||j<0||j>K||z.le(0)&&!z.zero()||!z.le(1)||z.eq(1))throw new RangeError('Retained derivative and 0<=z<1 required');
 let total=new F(0);
 for(let ell=0;ell<=j;ell++)total=total.add(new F(L.choose(K+1,j-ell)).mul(z.pow(K+1-j+ell)).div(new F(1).sub(z).pow(ell+1)));
 return total;
}
function sourceContract(jet){
 const s=jet.s,w=F.of(jet.w),a=F.of(jet.a),q0=I.of(jet.q0),c=jet.moments.map(I.of),K=c.length-1;
 if(!Number.isSafeInteger(s)||s<5||s>501||s%2!==1||!w.eq(new F(s).sub('1/2'))||!a.eq(w.pow(2))||!q0.positive()||K<1||K>31||jet.degree!==K||!c[0].lo.eq(1)||!c[0].hi.eq(1))throw new RangeError('Complete retained actual-source jet contract required');
 const majorant=new F(41).div(w.mul(2)).add('1/4'),M=majorant.div(q0.lo);
 return {s,w,a,q0,c,K,R,majorant,M};
}
function errors(jet,d,includeModel=false){
 const base=sourceContract(jet);d=F.of(d);
 if(d.le(0)&&!d.zero()||!d.le(R)||d.eq(R))throw new RangeError('0<=d<R required');
 const values=Array.from({length:base.K+1},(_,j)=>{
  let e=base.M.div(R.pow(j)).mul(exactTail(base.K,j,d.div(R)));
  if(includeModel)e=e.add(exactTail(base.K,j,d));
  return new F(1).add(d).pow(j).mul(e);
 });
 return {K:base.K,d,R,majorant:base.majorant,M:base.M,includeModel,values};
}
function matrix(K,rho){
 rho=F.of(rho);if(!Number.isSafeInteger(K)||K<0||K>31||!positive(rho))throw new RangeError('Finite retained degree and positive ratio required');
 return Array.from({length:K+1},(_,j)=>Array.from({length:K+1},(_,n)=>n<j?new F(0):new F(L.choose(n,j)).mul(rho.pow(j)).mul(new F(1).sub(rho).pow(n-j))));
}
function apply(raw,rho){
 const c=raw.map(I.of),P=matrix(c.length-1,rho);
 return P.map(row=>row.reduce((v,x,j)=>v.add(c[j].mul(x)),new I(0)));
}
function forecast(jet,ratios,degree){
 const base=sourceContract(jet);
 if(!Array.isArray(ratios)||!Number.isSafeInteger(degree)||degree<0||degree>base.K||ratios.some(x=>!positive(x)))throw new RangeError('Positive path and retained degree required');
 let rho=new F(1);const path=ratios.map(x=>{rho=rho.mul(x);if(!rho.sub(1).abs().le(R)||rho.sub(1).abs().eq(R))throw new RangeError('Every prefix must remain inside the original source disk');return rho;});
 const d=rho.sub(1).abs(),budget=errors(jet,d),polynomial=apply(base.c,rho).slice(0,degree+1),raw=polynomial.map((v,j)=>v.add(new I(budget.values[j].neg(),budget.values[j])));
 if(!raw[0].positive())throw new RangeError('Actual target normalization is not certified positive');
 return {source_s:base.s,source_a:base.a,source_degree:base.K,target_a:base.a.mul(rho),path,rho,budget,polynomial,raw,normalized:raw.map((v,j)=>j?v.div(raw[0]):new I(1)),q:base.q0.mul(raw[0])};
}
function corridor(jet,m,d){
 const base=sourceContract(jet),{K,c}=base;
 if(K%2!==1||!Number.isSafeInteger(m)||m<1||2*m>K+1)throw new RangeError('Complete odd-degree source and available target order required');
 const model=L.forms(c,(K+1)/2);
 if(!model.factors.A.positive||!model.factors.B.positive||!model.factors.H.positive)throw new RangeError('Certified complete positive model required');
 const bounds=errors(jet,d,true),e=bounds.values,one=new F(1);d=bounds.d;
 if(!positive(d))throw new RangeError('Positive corridor radius required');
 const f=L.forms(c,m),epsA=max(e.slice(1,2*m)),epsB=max(Array.from({length:2*m-1},(_,j)=>e[j].add(e[j+1])));
 const tauA=epsA.mul(m).mul(one.add(d).pow(2*m)).div(one.sub(d).pow(2*m-1)),tauB=epsB.mul(m).mul(one.add(d).pow(2*m)).div(one.sub(d).pow(2*m-2));
 const A=O.intervalLDL(f.A.map((row,i)=>row.map((v,j)=>v.sub(i===j?tauA:0)))),B=O.intervalLDL(f.B.map((row,i)=>row.map((v,j)=>v.sub(i===j?tauB:0))));
 const normalizationLower=one.div(one.add(d)).sub(e[0]);
 return {s:base.s,a:base.a,K,m,d,range:[base.a.mul(one.sub(d)),base.a.mul(one.add(d))],bounds,epsA,epsB,tauA,tauB,normalizationLower,model_order:(K+1)/2,factors:{A,B},certified:A.positive&&B.positive&&positive(normalizationLower)};
}
function glue(charts,m){
 if(!Array.isArray(charts)||!charts.length||!Number.isSafeInteger(m)||m<1||charts.some(x=>!x.certified||x.m<m))throw new RangeError('Certified charts of the requested order required');
 const sorted=charts.slice().sort((a,b)=>F.of(a.range[0]).eq(b.range[0])?0:F.of(a.range[0]).le(b.range[0])?-1:1),overlaps=[];
 const lower=F.of(sorted[0].range[0]);let upper=F.of(sorted[0].range[1]);
 for(let i=1;i<sorted.length;i++){
  const lo=F.of(sorted[i].range[0]),hi=F.of(sorted[i].range[1]);
  if(!lo.le(upper))throw new RangeError('Uncertified gap between source charts');
  overlaps.push([lo,hi.le(upper)?hi:upper]);if(upper.le(hi))upper=hi;
 }
 return {m,anchors:sorted.map(x=>x.s),range:[lower,upper],overlaps,connected:true};
}
module.exports=Object.freeze({exactTail,sourceContract,errors,matrix,apply,forecast,corridor,glue});
