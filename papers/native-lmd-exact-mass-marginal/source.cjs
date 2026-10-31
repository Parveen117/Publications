'use strict';
const assert=require('node:assert/strict');
const R=require('../native-lmd-first-mode-pulse-solver/source.cjs');
const {f,iv,point,add,sub,mul,scale,div,exp,sqrt,contains,gaussianPartition}=R;
const gammaCache=new Map(),lowerCache=new Map(),qCache=new Map(),eCache=new Map();

function power(a,n){
 assert(Number.isInteger(n)&&n>=0&&n<=32);let result=point(1),base=a;
 while(n){if(n%2)result=mul(result,base);n=Math.floor(n/2);if(n)base=mul(base,base);}
 return result;
}

function dyadicPower(x,s){
 x=f(x);s=f(s);assert(x.n>=0n&&s.n>0n&&s.le(2));
 assert(s.d<=16n&&(s.d&(s.d-1n))===0n,'supported dyadic exponent required');
 if(x.zero())return point(0);
 let root=point(x),den=s.d;
 while(den>1n){root=sqrt(root);den/=2n;}
 return power(root,Number(s.n));
}

function gammaLower(s,x){
 s=f(s);x=f(x);assert(s.n>0n&&s.le(2)&&x.n>=0n&&x.le(64));
 if(x.zero())return point(0);
 const key=s+':'+x;if(lowerCache.has(key))return lowerCache.get(key);
 let term=point(1),sum=point(1);const N=384;
 for(let n=1;n<=N;n++){term=scale(term,x.div(s.add(n)));sum=add(sum,term);}
 const tail=term.hi.mul(x).div(s.add(N+1)).div(f(1).sub(x.div(s.add(N+2))));
 const prefactor=scale(mul(dyadicPower(x,s),exp(point(x.neg()))),f(1).div(s));
 const value=mul(prefactor,add(sum,iv(0,tail)));lowerCache.set(key,value);return value;
}

function gamma(s){
 s=f(s);assert(s.n>0n&&s.le(2));const key=s.toString();if(gammaCache.has(key))return gammaCache.get(key);
 const x=f(64),lower=gammaLower(s,x),powerAtEndpoint=scale(dyadicPower(x,s),f(1).div(x));
 const exponent=s.sub(1),denominator=exponent.n>0n?f(1).sub(exponent.div(x)):f(1);
 const tail=scale(mul(exp(point(x.neg())),powerAtEndpoint),f(1).div(denominator));
 const value=add(lower,iv(0,tail.hi));gammaCache.set(key,value);return value;
}

function unitInterval(value){
 const lo=value.lo.n<0n?f(0):value.lo,hi=f(1).le(value.hi)?f(1):value.hi;
 assert(lo.le(hi));return iv(lo,hi);
}

function gammaCDF(s,x){
 if(x===Infinity)return point(1);
 if(f(x).zero())return point(0);
 return unitInterval(div(gammaLower(s,x),gamma(s)));
}

function exponentialPartition(h,subdivisions){
 h=f(h);assert(h.n>0n&&h.le('3/8'));
 assert([2,4].includes(subdivisions));
 const key=h+':'+subdivisions;if(eCache.has(key))return eCache.get(key);
 const shape=f(1).sub(h.div(2)),normalizer=gamma(shape),edges=[f(0),f(1).div(f(2).pow(16))];
 for(let j=-16;j<4;j++){
  const left=j<0?f(1).div(f(2).pow(-j)):f(2).pow(j);
  for(let m=1;m<=subdivisions;m++)edges.push(left.mul(f(1).add(f(m).div(subdivisions))));
 }
 edges.push(Infinity);
 const cdf=edges.map(x=>gammaCDF(shape,x)),survival=edges.map(x=>x===Infinity?point(0):exp(point(x.neg())));
 let totalP=point(0),totalA=point(0),second=point(0);const rows=[];
 for(let j=0;j<edges.length-1;j++){
  const p=unitInterval(sub(survival[j],survival[j+1])),first=unitInterval(sub(cdf[j+1],cdf[j]));
  assert(p.lo.n>0n&&first.lo.n>0n,'conditional cell lost strict positivity');
  rows.push({left:edges[j],right:edges[j+1],p,first,mass:div(first,p)});
  totalP=add(totalP,p);totalA=add(totalA,first);second=add(second,div(mul(first,first),p));
 }
 assert(contains(totalP,1)&&contains(totalA,1));
 const fullSecond=div(gamma(f(1).sub(h)),mul(normalizer,normalizer));
 assert(second.hi.le(fullSecond.lo));
 const result={rows,totalP,totalA,second,fullSecond,normalizer,shape};eCache.set(key,result);return result;
}

function productPartition(h,settings){
 h=f(h);assert([32,64].includes(settings.q_bins));
 const qkey=h+':'+settings.q_bins;
 if(!qCache.has(qkey))qCache.set(qkey,gaussianPartition(h,settings.q_bins));
 const q=qCache.get(qkey),exponential=exponentialPartition(h,settings.subdivisions_per_octave);
 const fullQ=exp(point(h)),fullSecond=mul(fullQ,exponential.fullSecond),second=mul(q.second,exponential.second);
 const residual=sub(fullSecond,second);
 const decomposition={q_bins:mul(sub(fullQ,q.second),exponential.fullSecond),
  exponential_bins:mul(q.second,sub(exponential.fullSecond,exponential.second))};
 assert(residual.lo.n>0n&&Object.values(decomposition).every(x=>x.lo.n>0n));
 const rows=[];
 for(let i=0;i<q.rows.length;i++)for(let j=0;j<exponential.rows.length;j++){
  rows.push({q_index:i,exponential_index:j,p:mul(q.rows[i].p,exponential.rows[j].p),mass:mul(q.rows[i].mass,exponential.rows[j].mass)});
 }
 return{q,exponential,rows,fullSecond,second,residual,decomposition};
}

module.exports={...R,power,dyadicPower,gammaLower,gamma,gammaCDF,exponentialPartition,productPartition};
