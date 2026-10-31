'use strict';
// Scientific arithmetic stays in the frozen native rational interval engine.
const assert=require('node:assert/strict');
const R=require('../native-lmd-certified-pulse-solver/interval.cjs');
const {f,iv,point,add,sub,mul,scale,div,exp,normalCDF,contains}=R;

function bessel0(x,N=32){
 x=f(x);assert(x.n>=0n&&x.le(1));assert(Number.isInteger(N)&&N>=8&&N<=64);
 const a=x.pow(2).div(4);let term=point(1),sum=point(1);
 for(let n=1;n<=N;n++){term=scale(term,a.div(n*n));sum=add(sum,term);}
 const tail=term.hi.mul(a).div((N+1)**2).div(f(1).sub(a.div((N+2)**2)));
 return add(sum,iv(0,tail));
}

function gammaTail(x,N){
 if(x===Infinity)return Array.from({length:N+1},()=>point(0));
 x=f(x);assert(x.n>=0n&&x.le(16));
 let term=point(1),sum=point(1);const ex=exp(point(x.neg())),values=[ex];
 for(let n=1;n<=N;n++){term=scale(term,x.div(n));sum=add(sum,term);values.push(mul(ex,sum));}
 return values;
}

function unitInterval(x){
 const lo=x.lo.n<0n?f(0):x.lo,hi=f(1).le(x.hi)?f(1):x.hi;
 assert(lo.le(hi));return iv(lo,hi);
}

function summarize(rows){
 let totalP=point(0),totalA=point(0),second=point(0);
 for(const row of rows){
  assert(row.p.lo.n>0n&&row.first.lo.n>0n);
  row.mass=div(row.first,row.p);
  totalP=add(totalP,row.p);totalA=add(totalA,row.first);
  second=add(second,div(mul(row.first,row.first),row.p));
 }
 assert(contains(totalP,1)&&contains(totalA,1),'partition lost probability or mass');
 return{rows,totalP,totalA,second};
}

function gaussianPartition(g,bins,radius=f(2)){
 g=f(g);radius=f(radius);assert(g.n>0n&&g.le('3/4'));
 assert(Number.isInteger(bins)&&bins>=2&&bins<=128&&radius.eq(2));
 const endpoints=[null];
 for(let j=0;j<=bins;j++)endpoints.push(radius.neg().add(radius.mul(2*j).div(bins)));
 endpoints.push(Infinity);
 const rows=[];
 for(let j=0;j<endpoints.length-1;j++){
  const left=endpoints[j],right=endpoints[j+1];
  const shift=x=>x===null||x===Infinity?x:x.sub(g);
  rows.push({left,right,p:sub(normalCDF(g,right),normalCDF(g,left)),
    first:sub(normalCDF(g,shift(right)),normalCDF(g,shift(left)))});
 }
 return summarize(rows);
}

function radialPartition(g,bins,end=f(16),N=32){
 g=f(g);end=f(end);assert(g.n>0n&&g.le('3/4'));
 assert(Number.isInteger(bins)&&bins>=2&&bins<=32&&end.eq(16));
 assert(Number.isInteger(N)&&N>=8&&N<=64);
 const lam=g.div(2),prefactor=exp(point(lam.neg()));
 let weight=prefactor;const weights=[weight];
 for(let n=1;n<=N;n++){weight=scale(weight,lam.div(n));weights.push(weight);}
 const poissonRemainder=weight.hi.mul(lam).div(N+1).div(f(1).sub(lam.div(N+2)));
 const edges=Array.from({length:bins+1},(_,j)=>end.mul(j).div(bins));edges.push(Infinity);
 const tails=edges.map(x=>gammaTail(x,N)),rows=[];
 for(let j=0;j<edges.length-1;j++){
  const p=unitInterval(sub(tails[j][0],tails[j+1][0]));let first=point(0);
  for(let n=0;n<=N;n++){
   const probability=unitInterval(sub(tails[j][n],tails[j+1][n]));
   first=add(first,mul(weights[n],probability));
  }
  first=add(first,iv(0,poissonRemainder));
  rows.push({left:edges[j],right:edges[j+1],p,first});
 }
 return{...summarize(rows),poissonRemainder};
}

function productPartition(g,bins,sourceI){
 g=f(g);const q=gaussianPartition(g,bins.q_bins),radial=radialPartition(g,bins.radial_bins);
 const exactQ=exp(point(g)),exactR=bessel0(g),idealSecond=mul(exactQ,exactR);
 const second=mul(q.second,radial.second),residual=sub(sourceI,second);
 assert(residual.lo.n>0n,'nonpositive source variance budget');
 const floors={zero:sub(sourceI,exactQ),first:sub(sourceI,idealSecond)};
 const decomposition={higher_modes:floors.first,
  q_bins:mul(sub(exactQ,q.second),exactR),radial_bins:mul(q.second,sub(exactR,radial.second))};
 assert(Object.values(decomposition).every(x=>x.lo.n>0n));
 const rows=[];
 for(let i=0;i<q.rows.length;i++)for(let j=0;j<radial.rows.length;j++){
  rows.push({q_index:i,radial_index:j,p:mul(q.rows[i].p,radial.rows[j].p),mass:mul(q.rows[i].mass,radial.rows[j].mass)});
 }
 return{q,radial,rows,second,residual,floors,decomposition,idealSecond};
}

module.exports={...R,bessel0,gammaTail,gaussianPartition,radialPartition,productPartition};
