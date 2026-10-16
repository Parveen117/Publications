'use strict';
/** ILC-1: complete single-point xi layers; exact outward native arithmetic. */
const {F}=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const {Interval:I}=require('../native-prime-observation-closure/observation_adapter.cjs');
const S=require('../native-prime-innovation-window/full_source.cjs');
const factorial=n=>{let r=1n;for(let k=2;k<=n;k++)r*=BigInt(k);return r;};
const rising=(r,n)=>{let z=1n;for(let k=0;k<n;k++)z*=BigInt(r+k);return z;};
function integer(n,lo,hi){if(!Number.isSafeInteger(n)||n<lo||n>hi)throw new RangeError('Integer outside declared layer range');return n;}
function power(x,n){x=I.of(x);let r=new I(1);for(let k=0;k<n;k++)r=r.mul(x);return r;}
function product(a,b,degree){
 const out=Array.from({length:degree+1},()=>new I(0));
 for(let i=0;i<a.length&&i<=degree;i++)for(let j=0;j<b.length&&i+j<=degree;j++)out[i+j]=out[i+j].add(I.of(a[i]).mul(b[j]));
 return out;
}
function compose(a,b,degree){
 if(!I.of(b[0]).lo.zero()||!I.of(b[0]).hi.zero())throw new RangeError('Composition requires zero inner constant');
 let out=Array.from({length:degree+1},()=>new I(0)),p=[new I(1)];
 for(let k=0;k<a.length&&k<=degree;k++){for(let j=0;j<p.length;j++)out[j]=out[j].add(p[j].mul(a[k]));p=product(p,b,degree);}
 return out;
}
function binomialSeries(exponent,degree){
 exponent=F.of(exponent);let c=new F(1);const out=[c];
 for(let n=1;n<=degree;n++){c=c.mul(exponent.sub(n-1)).div(n);out.push(c);}return out;
}
function hurwitz(r,v,{shift=512,terms=32}={}){
 integer(r,2,32);integer(shift,1,1024);integer(terms,1,49);v=F.of(v);if(v.le(0))throw new RangeError('Positive Hurwitz parameter required');
 const y=v.add(shift),B=S.bernoulli(2*terms);let finite=new I(0);
 for(let k=0;k<shift;k++)finite=finite.add(new F(1).div(v.add(k).pow(r)));
 let tail=new I(new F(1).div(y.pow(r-1).mul(r-1)).add(new F(1).div(y.pow(r).mul(2))));
 let last;
 for(let k=1;k<=terms;k++){last=B[2*k].mul(rising(r,2*k-1)).div(factorial(2*k)).div(y.pow(r+2*k-1));tail=tail.add(last);}
 // Euler--Maclaurin includes the last printed term. Its periodic-Bernoulli
 // integral is bounded by the magnitude of that term, not silently dropped.
 const remainder=last.abs();return {value:finite.add(tail).add(new I(remainder.neg(),remainder)),r,v,shift,terms,remainder};
}
function eulerian(degree){
 integer(degree,0,31);const rows=[[1n]];
 for(let k=0;k<degree;k++){
  const a=rows[k],b=Array(a.length+1).fill(0n);
  for(let j=0;j<a.length;j++){b[j]+=BigInt(j+1)*a[j];b[j+1]+=BigInt(k-j)*a[j];}
  while(b.length>1&&b.at(-1)===0n)b.pop();rows.push(b);
 }return rows;
}
function logarithmicIntegral(degree,s,cutoff){
 const l=S.logRational(cutoff),a=new F(s-1);let out=new I(0);
 for(let j=0;j<=degree;j++)out=out.add(power(l,degree-j).mul(new F(factorial(degree),factorial(degree-j))).div(a.pow(j+1)));
 return out.mul(new F(1,BigInt(cutoff)**BigInt(s-1)));
}
function eulerLayers(s,degree,cutoff=4096){
 integer(s,3,501);integer(degree,0,31);integer(cutoff,2,65536);
 if(S.logRational(cutoff).lo.mul(s).le(degree+1))throw new RangeError('Omitted log-power integrands must be decreasing');
 const rows=eulerian(degree),primes=S.primesThrough(cutoff),sums=Array.from({length:degree+1},()=>new I(0));
 for(const p of primes){
  const P=BigInt(p)**BigInt(s),den=P-1n,log=S.logRational(p);let lp=new I(1);
  for(let k=0;k<=degree;k++){
   lp=lp.mul(log);let numerator=0n;
   for(let j=0;j<rows[k].length;j++)numerator+=rows[k][j]*P**BigInt(k-j);
   const weight=new F(numerator,den**BigInt(k+1)*factorial(k));
   sums[k]=sums[k].add(lp.mul(weight));
  }
 }
 const c=new F(1).sub(new F(1,BigInt(cutoff+1)**BigInt(s))),tails=[],coefficients=[];
 for(let k=0;k<=degree;k++){
  const tail=logarithmicIntegral(k+1,s,cutoff).div(c.pow(k+1)).hi;
  tails.push(tail);const positive=sums[k].add(new I(0,tail));coefficients.push(k%2?positive:positive.neg());
 }
 return {coefficients,retained_positive_sums:sums,tails,cutoff,retained_prime_count:primes.length,last_retained_prime:primes.at(-1),eulerian_rows:rows};
}
function sourceLayers(s,degree=15,{cutoff=4096,shift=512,terms=32}={}){
 integer(s,3,501);integer(degree,1,31);if(s%2!==1)throw new RangeError('Odd original source argument required');
 const w=new F(2*s-1,2),a=w.pow(2),v=new F(s,2),cc=S.constants(),psi=S.digamma(v),euler=eulerLayers(s,degree,cutoff),gamma=[],poles=[],Q=[],hurwitz_records=[];
 for(let k=0;k<=degree;k++){
  const sign=k%2?-1:1;const pole=new F(sign).mul(new F(1).div(new F(s).pow(k+1)).add(new F(1).div(new F(s-1).pow(k+1))));
  let g;if(k===0)g=psi.value.sub(cc.logpi).div(2);else{const h=hurwitz(k+1,v,{shift,terms});hurwitz_records.push(h);g=h.value.mul(k%2?1:-1).div(2n**BigInt(k+1));}
  poles.push(pole);gamma.push(g);Q.push(g.add(pole).add(euler.coefficients[k]));
 }
 const root=binomialSeries('1/2',degree),inverseRoot=binomialSeries('-1/2',degree),deltaS=root.map((c,k)=>k?c.mul(w):new F(0));
 const qLayers=product(compose(Q,deltaS,degree),inverseRoot,degree).map(c=>c.div(w.mul(2))),q0=qLayers[0];
 if(!q0.positive())throw new RangeError('Positive normalization not certified');
 const moments=qLayers.map((c,k)=>k?c.mul(k%2?-1:1).div(q0):new I(1));
 return {s,w,a,degree,q0,q_layers:qLayers,moments,Q_layers:Q,components:{poles,gamma,psi,hurwitz:hurwitz_records,euler},protocol:{lambda:'(t-a)/a',moments:'(-1)^n [lambda^n] q(a(1+lambda))/q(a)',rounding_decimal_places:110,source:'complete xi; all source tails retained; no assumed zero measure'}};
}
module.exports=Object.freeze({factorial,rising,power,product,compose,binomialSeries,hurwitz,eulerian,logarithmicIntegral,eulerLayers,sourceLayers});
