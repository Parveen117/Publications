'use strict';
// Certified elementary functions over the frozen native rational field.
const {F}=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const assert=require('node:assert/strict'),f=F.of,BITS=192,S=1n<<BigInt(BITS);
const floor=x=>{x=f(x);let q=x.n/x.d;if(x.n<0n&&x.n%x.d)q--;return q;};
const ceil=x=>-floor(f(x).neg());
const down=x=>new F(floor(f(x).mul(S)),S),up=x=>new F(ceil(f(x).mul(S)),S);
const iv=(a,b=a)=>({lo:f(a),hi:f(b)}),point=x=>iv(x);
const min=xs=>xs.reduce((a,b)=>a.le(b)?a:b),max=xs=>xs.reduce((a,b)=>a.le(b)?b:a);
const add=(a,b)=>iv(down(a.lo.add(b.lo)),up(a.hi.add(b.hi)));
const neg=a=>iv(a.hi.neg(),a.lo.neg()),sub=(a,b)=>add(a,neg(b));
const mul=(a,b)=>{const v=[a.lo.mul(b.lo),a.lo.mul(b.hi),a.hi.mul(b.lo),a.hi.mul(b.hi)];return iv(down(min(v)),up(max(v)));};
const scale=(a,x)=>mul(a,point(x));
const inv=a=>{assert(a.lo.n>0n||a.hi.n<0n,'division interval contains zero');return iv(down(f(1).div(a.hi)),up(f(1).div(a.lo)));};
const div=(a,b)=>mul(a,inv(b));
const mid=a=>a.lo.add(a.hi).div(2),width=a=>a.hi.sub(a.lo),absmax=a=>max([a.lo.abs(),a.hi.abs()]);
const widen=(a,e)=>{e=f(e);assert(e.n>=0n);return iv(down(a.lo.sub(e)),up(a.hi.add(e)));};
const overlap=(a,b)=>a.lo.le(b.hi)&&b.lo.le(a.hi),contains=(a,x)=>a.lo.le(x)&&f(x).le(a.hi);
function isqrt(n){assert(n>=0n);if(n<2n)return n;let x=n,y=(x+1n)/2n;while(y<x){x=y;y=(x+n/x)/2n;}return x;}
function sqrt(a){assert(a.lo.n>=0n);return iv(new F(isqrt(a.lo.n*S*S/a.lo.d),S),new F(isqrt(a.hi.n*S*S/a.hi.d)+1n,S));}
function expPoint(x){
 x=f(x);if(x.n<0n)return inv(expPoint(x.neg()));assert(x.le(128));
 let r=x,squarings=0;while(!r.le(1)){r=r.div(2);squarings++;}
 let term=point(1),sum=point(1);
 for(let n=1;n<=96;n++){term=scale(term,r.div(n));sum=add(sum,term);}
 const tail=term.hi.mul(r).div(97).div(f(1).sub(r.div(98)));
 sum=iv(sum.lo,up(sum.hi.add(tail)));for(let n=0;n<squarings;n++)sum=mul(sum,sum);return sum;
}
const exp=a=>iv(expPoint(a.lo).lo,expPoint(a.hi).hi);
const sinh=a=>scale(sub(exp(a),exp(neg(a))),'1/2');
function atanReciprocal(q){
 const x=f(1).div(q);let term=point(x),sum=point(x);
 for(let n=1;n<=160;n++){term=scale(term,x.pow(2).neg().mul(2*n-1).div(2*n+1));sum=add(sum,term);}
 return widen(sum,absmax(term).mul(x.pow(2)).mul(321).div(323));
}
const PI=sub(scale(atanReciprocal(5),16),scale(atanReciprocal(239),4));
const cdfCache=new Map();
function normalCDF(gamma,x){
 gamma=f(gamma);if(x===null)return point(0);if(x===Infinity)return point(1);
 x=f(x);const key=gamma+':'+x;if(cdfCache.has(key))return cdfCache.get(key);
 if(x.n<0n){const r=sub(point(1),normalCDF(gamma,x.neg()));cdfCache.set(key,r);return r;}
 const A=x.pow(2).div(gamma.mul(2));assert(A.le(64));
 let term=point(x),sum=point(x);
 for(let n=1;n<=256;n++){term=scale(term,A.neg().mul(2*n-1).div(f(n).mul(2*n+1)));sum=add(sum,term);}
 const next=absmax(term).mul(A).mul(513).div(f(257).mul(515));
 const r=add(point('1/2'),div(widen(sum,next),sqrt(scale(PI,gamma.mul(2)))));
 assert(!r.lo.le('-1/100000000000000000000')&&r.hi.le('100000000000000000001/100000000000000000000'));
 cdfCache.set(key,r);return r;
}
function trig(a){
 const tau=scale(PI,2),period=floor(mid(div(a,tau)).add('1/2'));
 const x=sub(a,scale(tau,period)),x2=mul(x,x);assert(absmax(x).le(4));
 let sin=x,cos=point(1),st=x,ct=point(1);
 for(let n=1;n<=80;n++){
  st=scale(mul(st,x2),f(-1).div((2*n)*(2*n+1)));sin=add(sin,st);
  ct=scale(mul(ct,x2),f(-1).div((2*n-1)*(2*n)));cos=add(cos,ct);
 }
 // Bound the absolute tail, without relying on alternating interval signs.
 const m=absmax(x),m2=m.pow(2);
 const se=absmax(st).mul(m2).div(162*163).div(f(1).sub(m2.div(164*165)));
 const ce=absmax(ct).mul(m2).div(161*162).div(f(1).sub(m2.div(163*164)));
 return{sin:widen(sin,se),cos:widen(cos,ce)};
}
function decimal(s){if(!s.includes('.'))return f(s);const sign=s.startsWith('-')?-1n:1n,t=s.replace(/^[-+]/,''),[a,b]=t.split('.');return new F(sign*BigInt(a+b),10n**BigInt(b.length));}
function fixed(x,places=9,direction=0){const s=10n**BigInt(places),v=f(x).mul(s),n=direction<0?floor(v):direction>0?ceil(v):v.n/v.d,abs=n<0n?-n:n;return(n<0n?'-':'')+String(abs/s)+'.'+String(abs%s).padStart(places,'0');}
const pack=a=>({lower:a.lo.toString(),upper:a.hi.toString(),display_lower:fixed(a.lo,9,-1),display_upper:fixed(a.hi,9,1)});
module.exports={F,f,BITS,iv,point,add,neg,sub,mul,scale,inv,div,mid,width,absmax,widen,overlap,contains,sqrt,exp,sinh,PI,normalCDF,trig,decimal,fixed,pack,floor};
