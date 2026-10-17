'use strict';
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const O=require('../native-prime-observation-closure/observation_adapter.cjs');
const S=require('../native-prime-innovation-window/window_adapter.cjs');
const L=require('../native-infinite-layer-cut/layer_adapter.cjs');
const T=require('../native-layer-observer-transport/transport_adapter.cjs');
const {F}=N,{Interval:I}=O;
function integer(n,min=0){if(!Number.isSafeInteger(n)||n<min)throw new RangeError('Nonnegative safe index required');return n;}
function trim(p){p=p.map(F.of);if(!p.length)p=[new F(0)];while(p.length>1&&p[p.length-1].zero())p.pop();return p;}
function add(a,b){return trim(Array.from({length:Math.max(a.length,b.length)},(_,i)=>F.of(a[i]||0).add(b[i]||0)));}
function scale(a,v){return trim(a.map(x=>F.of(x).mul(v)));}
function multiply(a,b){const c=Array.from({length:a.length+b.length-1},()=>new F(0));a.forEach((x,i)=>b.forEach((y,j)=>{c[i+j]=c[i+j].add(F.of(x).mul(y));}));return trim(c);}
function action(pair,kind,steps=1){
 integer(steps);if(!Array.isArray(pair)||pair.length!==2||pair.some(p=>!Array.isArray(p)))throw new TypeError('Two real coefficient polynomials required');
 let [x,z]=pair.map(trim);if(!['wave','inverse','generator'].includes(kind))throw new RangeError('Known source action required');
 for(let k=0;k<steps;k++){
  if(kind==='wave'){const u=multiply([0,2],add(x,z));[x,z]=[add(u,add(scale(x,-1),scale(z,-2))),add(u,scale(z,-1))];}
  else if(kind==='inverse'){const u=multiply([0,2],add(x,scale(z,-1)));[x,z]=[add(u,add(scale(x,-1),scale(z,2))),add(scale(u,-1),scale(z,-1))];}
  else [x,z]=[multiply([-2,2],z),multiply([0,2],x)];
 }
 return [x,z];
}
function pairColumn(pair,n){integer(n,1);const p=pair.map(trim);if(p.length!==2||p.some(x=>x.length>n))throw new RangeError('Polynomial does not fit requested carrier');return N.matrix(p.flatMap(a=>Array.from({length:n},(_,j)=>[a[j]||new F(0)])));}
function correlation(u,v){return add(multiply([0,1],multiply(u[0],v[0])),multiply([1,-1],multiply(u[1],v[1])));}
function functional(raw,poly){poly=trim(poly);if(raw.length<poly.length)throw new RangeError('Unobserved source coefficient required');return poly.reduce((s,v,j)=>s.add(I.of(raw[j]).mul(v)),new I(0));}
function factorial(n){integer(n);let f=new F(1);for(let j=2;j<=n;j++)f=f.mul(j);return f;}
function tailBound(r,t){integer(r,1);const tau=F.of(t).abs();if(new F(r+1).le(tau))throw new RangeError('Tail ratio must be below one');return tau.pow(r).div(factorial(r)).div(new F(1).sub(tau.div(r+1)));}
function budgets(n,m,t){integer(n,1);integer(m,n);const d=m-n;return {n,m,margin:d,t:F.of(t),relativeStateBound:tailBound(d+1,t).mul(2),relativeCorrelationBound:tailBound(2*d+1,t).mul(2)};}
function generator(cut){return N.sub(N.mul(cut.cut,cut.Z),N.mul(cut.Z,cut.cut));}
function finiteReturns(raw,m,maxStep,kind='wave'){
 if(!['wave','generator'].includes(kind))throw new RangeError('Known finite return action required');
 integer(maxStep);const f=T.strict(raw,m),columns=f.A.map((_,j)=>S.solve(f.factors.H,f.A.map(row=>row[j]))),Y=columns.map((_,i)=>columns.map(col=>col[i])),seed=Array.from({length:m},(_,j)=>new I(j===0?1:0)),dot=(a,b)=>a.reduce((s,v,j)=>s.add(v.mul(b[j])),new I(0)),apply=v=>Y.map(row=>dot(row,v)),reading=v=>dot(f.A[0],v),values=[reading(seed)];
 if(kind==='generator'){
  let x=seed,z=seed.map(()=>new I(0));
  for(let k=1;k<=maxStep;k++){const yz=apply(z),yx=apply(x);[x,z]=[yz.map((v,j)=>v.sub(z[j]).mul(2)),yx.map(v=>v.mul(2))];values.push(reading(x));}
  return {m,kind,values};
 }
 let previous=seed,current=apply(seed).map((v,j)=>v.mul(2).sub(seed[j]));if(maxStep>=1)values.push(reading(current));
 for(let k=2;k<=maxStep;k++){const next=apply(current).map((v,j)=>v.mul(4).sub(current[j].mul(2)).sub(previous[j]));values.push(reading(next));previous=current;current=next;}
 return {m,kind,values};
}
function returnPolynomial(k){integer(k);return correlation([[1],[0]],action([[1],[0]],'wave',k));}
function exponentialReading(raw,m,t){
 T.strict(raw,m);t=F.of(t);let partial=new I(0);const terms=[];
 for(let j=0;j<m;j++){const moment=L.difference(raw,j+1,j),coefficient=t.pow(2).mul(-4).pow(j).div(factorial(2*j)),term=moment.mul(coefficient);partial=partial.add(term);terms.push({j,moment,coefficient,term});}
 const normalizedRadius=t.abs().pow(2*m).div(factorial(2*m)),radius=I.of(raw[1]).hi.mul(normalizedRadius),enclosure=partial.add(new I(radius.neg(),radius));
 return {m,t,terms,partial,radius,normalizedRadius,enclosure,normalized:enclosure.div(raw[1])};
}
module.exports=Object.freeze({trim,add,scale,multiply,action,pairColumn,correlation,functional,factorial,tailBound,budgets,generator,finiteReturns,returnPolynomial,exponentialReading});
