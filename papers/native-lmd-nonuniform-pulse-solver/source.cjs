'use strict';
const assert=require('node:assert/strict');
const R=require('../native-lmd-exact-mass-marginal/source.cjs');
const {f,iv,point,add,sub,mul,scale,div,exp,sqrt,widen,absmax,contains,bessel0,gamma,gammaCDF,gaussianPartition,radialPartition}=R;
const radialCache=new Map(),qCache=new Map(),rowCache=new Map();
const clampUnit=x=>iv(x.lo.n<0n?0:x.lo,f(1).le(x.hi)?1:x.hi);

function bessel1(x,N=32){
 x=f(x);assert(x.n>=0n&&x.le(1));assert(N===32);
 const A=x.pow(2).div(4);let term=point(x.div(2)),sum=term;
 for(let n=1;n<=N;n++){term=scale(term,A.div(n*(n+1)));sum=add(sum,term);}
 const tail=term.hi.mul(A).div((N+1)*(N+2)).div(f(1).sub(A.div((N+2)*(N+3))));
 return add(sum,iv(0,tail));
}

function halfGammaTails(x,N=32){
 assert(N===32);
 if(x===Infinity)return Array.from({length:N+1},()=>point(0));
 x=f(x);assert(x.n>=0n&&x.le(16));
 if(x.zero())return Array.from({length:N+1},()=>point(1));
 let q=sub(point(1),gammaCDF('1/2',x));
 let inc=div(mul(exp(point(x.neg())),sqrt(point(x))),gamma('3/2'));
 const values=[];
 for(let n=0;n<=N;n++){
  q=add(q,inc);values.push(clampUnit(q));
  inc=scale(inc,x.div(f(n).add('3/2')));
 }
 return values;
}

function radialWeighted(h,bins){
 h=f(h);assert(h.n>0n&&h.le('3/8'));assert([16,32].includes(bins));
 const key=h+':'+bins;if(radialCache.has(key))return radialCache.get(key);
 const base=radialPartition(h,bins),lambda=h.div(2),N=32;
 let weight=mul(mul(exp(point(lambda.neg())),sqrt(point(lambda))),gamma('3/2'));
 const weights=[weight],ratio=n=>lambda.mul(f(n).add('3/2')).div((n+1)*(n+2));
 for(let n=0;n<N;n++){weight=scale(weight,ratio(n));weights.push(weight);}
 const remainder=weight.hi.mul(ratio(N)).div(f(1).sub(ratio(N+1)));
 const edges=base.rows.map(r=>r.left).concat([Infinity]),tails=edges.map(x=>halfGammaTails(x));
 let second1=point(0),first1Total=point(0);
 const rows=base.rows.map((row,j)=>{
  let first1=point(0);
  for(let n=0;n<=N;n++)first1=add(first1,mul(weights[n],clampUnit(sub(tails[j][n],tails[j+1][n]))));
  first1=add(first1,iv(0,remainder));
  first1=iv(first1.lo.n<0n?0:first1.lo,first1.hi.le(row.first.hi)?first1.hi:row.first.hi);
  assert(first1.lo.n>0n&&first1.hi.le(row.first.lo),'positive weighted cell dominance');
  second1=add(second1,div(mul(first1,first1),row.p));first1Total=add(first1Total,first1);
  return{...row,first1,mass1:div(first1,row.p)};
 });
 const result={...base,rows,second1,first1Total,mixtureRemainder:remainder,
  coefficientTotal:add(weights.reduce(add,point(0)),iv(0,remainder))};
 assert(base.second.hi.le(bessel0(h).lo)&&second1.hi.le(bessel1(h).lo));
 radialCache.set(key,result);return result;
}

function productPartition(h,e,bins,sourceI){
 h=f(h);e=f(e);assert(e.abs().le('3/4'));assert([64,128].includes(bins.q_bins));
 const qkey=h+':'+bins.q_bins;if(!qCache.has(qkey))qCache.set(qkey,gaussianPartition(h,bins.q_bins));
 const q=qCache.get(qkey),radial=radialWeighted(h,bins.radial_bins),e2=e.pow(2).div(2);
 const idealRadial=add(bessel0(h),scale(bessel1(h),e2)),capturedRadial=add(radial.second,scale(radial.second1,e2));
 const fullSecond=scale(sourceI,f(1).add(e2.mul(h).div(f(2).sub(h))));
 const fullQ=exp(point(h)),idealSecond=mul(fullQ,idealRadial),second=mul(q.second,capturedRadial);
 const residual=sub(fullSecond,second),floor=sub(fullSecond,idealSecond);
 const decomposition={higher_modes:floor,q_bins:mul(sub(fullQ,q.second),idealRadial),
  radial_bins:mul(q.second,sub(idealRadial,capturedRadial))};
 assert(residual.lo.n>0n&&Object.values(decomposition).every(x=>x.lo.n>0n));
 const rowkey=qkey+':'+bins.radial_bins;
 if(!rowCache.has(rowkey)){
  const rows=[];
  for(let i=0;i<q.rows.length;i++)for(let j=0;j<radial.rows.length;j++)rows.push({q_index:i,radial_index:j,
   p:mul(q.rows[i].p,radial.rows[j].p),mass0:mul(q.rows[i].mass,radial.rows[j].mass),mass1:mul(q.rows[i].mass,radial.rows[j].mass1)});
  rowCache.set(rowkey,rows);
 }
 return{q,radial,rows:rowCache.get(rowkey),fullSecond,idealSecond,second,residual,floor,decomposition};
}

function angular0(x){
 const m=absmax(x);
 if(m.zero())return{value:point(1),tail:f(0),terms:0,fallback:false};
 if(!m.le(32))return{value:iv(-1,1),tail:f(1),terms:0,fallback:true};
 const N=m.le(1)?16:m.le(4)?24:m.le(16)?40:64;
 const a=scale(mul(x,x),'-1/4'),A=m.pow(2).div(4);let term=point(1),sum=term;
 for(let n=1;n<=N;n++){term=scale(mul(term,a),f(1).div(n*n));sum=add(sum,term);}
 const tail=absmax(term).mul(A).div((N+1)**2).div(f(1).sub(A.div((N+2)**2)));
 const enclosed=widen(sum,tail),lo=f(-1).le(enclosed.lo)?enclosed.lo:f(-1),hi=enclosed.hi.le(1)?enclosed.hi:f(1);
 assert(lo.le(hi));return{value:iv(lo,hi),tail,terms:N,fallback:false};
}
module.exports={...R,bessel1,halfGammaTails,radialWeighted,productPartition,angular0};
