'use strict';
/** PIW-1 complete positive-real xi jets. Native exact outward arithmetic. */
const {F}=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const {Interval:I}=require('../native-prime-observation-closure/observation_adapter.cjs');
const LOG_TERMS=128,ATAN_TERMS=96,PSI_SHIFT=256,PSI_TERMS=24,PRIME_CUTOFF=4096;
const cache=new Map();
function boundedInt(n,lo,hi){if(!Number.isSafeInteger(n)||n<lo||n>hi)throw new RangeError('Exact integer outside declared source range');return n;}
function atanhLog(x,terms=LOG_TERMS){
 x=F.of(x);boundedInt(terms,1,512);if(x.le(0)||!x.le(2)||new F(1).le(x)===false)throw new RangeError('Reduced logarithm requires 1 <= x <= 2');
 const z=new I(x.sub(1).div(x.add(1))),zz=z.square();let power=z,sum=new I(0);
 for(let k=0;k<terms;k++){sum=sum.add(power.div(2*k+1));power=power.mul(zz);}
 const tail=power.mul(2).div(2*terms+1).div(new I(1).sub(zz));
 return sum.mul(2).add(new I(0,tail.hi));
}
function logRational(x){
 x=F.of(x);if(x.le(0))throw new RangeError('Logarithm requires a positive rational');
 const key='log:'+x.toString();if(cache.has(key))return cache.get(key);
 let y=x,k=0;while(!y.le(2)||y.eq(2)){y=y.div(2);if(++k>2048)throw new RangeError('Logarithm scale bound');}
 while(!new F(1).le(y)){y=y.mul(2);if(--k < -2048)throw new RangeError('Logarithm scale bound');}
 if(!cache.has('reduced_log2'))cache.set('reduced_log2',atanhLog(2));
 const result=atanhLog(y).add(cache.get('reduced_log2').mul(k));cache.set(key,result);return result;
}
function atanReciprocal(n,terms=ATAN_TERMS){
 boundedInt(n,2,10000);boundedInt(terms,1,512);const x=new I(new F(1,n)),xx=x.square();let power=x,sum=new I(0);
 for(let k=0;k<terms;k++){const term=power.div(2*k+1);sum=k%2?sum.sub(term):sum.add(term);power=power.mul(xx);}
 const tail=power.div(2*terms+1).hi;return sum.add(new I(tail.neg(),tail));
}
function constants(){
 if(cache.has('constants'))return cache.get('constants');
 const pi=atanReciprocal(5).mul(16).sub(atanReciprocal(239).mul(4));
 const result={pi,logpi:new I(logRational(pi.lo).lo,logRational(pi.hi).hi),log2:logRational(2)};
 cache.set('constants',result);return result;
}
function choose(n,k){let b=1n;for(let j=1;j<=k;j++)b=b*BigInt(n-j+1)/BigInt(j);return b;}
function bernoulli(n){
 boundedInt(n,0,100);const key='bernoulli:'+n;if(cache.has(key))return cache.get(key);
 const B=[new F(1)];for(let m=1;m<=n;m++){let s=new F(0);for(let k=0;k<m;k++)s=s.add(B[k].mul(choose(m+1,k)));B.push(s.neg().div(m+1));}
 cache.set(key,B);return B;
}
function digamma(v,shift=PSI_SHIFT,terms=PSI_TERMS){
 v=F.of(v);if(v.le(0))throw new RangeError('Positive digamma argument required');boundedInt(shift,0,1024);boundedInt(terms,1,49);
 const y=v.add(shift),B=bernoulli(2*terms+2);let approx=logRational(y).sub(new F(1).div(y.mul(2)));
 for(let k=1;k<=terms;k++)approx=approx.sub(B[2*k].div(y.pow(2*k).mul(2*k)));
 const firstOmitted=B[2*terms+2].neg().div(y.pow(2*terms+2).mul(2*terms+2));
 const remainder=firstOmitted.le(0)?new I(firstOmitted,0):new I(0,firstOmitted);
 let recurrence=new F(0);for(let k=0;k<shift;k++)recurrence=recurrence.add(new F(1).div(v.add(k)));
 return {value:approx.add(remainder).sub(recurrence),shift,terms,shifted_argument:y,first_omitted:firstOmitted,remainder};
}
function trigammaOdd(s){
 boundedInt(s,3,501);if(s%2!==1)throw new RangeError('Odd original s required');
 let sum=new F(0);for(let k=1;k<=(s-1)/2;k++)sum=sum.add(new F(1,BigInt(2*k-1)**2n));
 return constants().pi.square().div(2).sub(sum.mul(4));
}
function primesThrough(n){boundedInt(n,2,65536);const marked=Array(n+1).fill(false),primes=[];for(let p=2;p<=n;p++)if(!marked[p]){primes.push(p);for(let m=p*p;m<=n;m+=p)marked[m]=true;}return primes;}
function eulerJet(s,cutoff=PRIME_CUTOFF){
 boundedInt(s,3,501);boundedInt(cutoff,2,65536);const primes=primesThrough(cutoff);let sum0=new I(0),sum1=new I(0);
 // Form the exact rational weight before interval multiplication. Splitting
 // ps / (ps-1)^2 into a huge product and a sub-grid reciprocal is valid but
 // destroys useful precision on the fixed outward grid.
 for(const p of primes){const lp=logRational(p),ps=new F(BigInt(p)**BigInt(s)),den=ps.sub(1);sum0=sum0.add(lp.mul(new F(1).div(den)));sum1=sum1.add(lp.square().mul(ps.div(den.pow(2))));}
 // Bound all omitted primes and their complete powers by decreasing
 // integer sums, then the exact logarithmic integrals from cutoff to infinity.
 const l=logRational(cutoff),a=new F(s-1),power=new F(1,BigInt(cutoff)**BigInt(s-1)),factor=new F(1).sub(new F(1,BigInt(cutoff+1)**BigInt(s)));
 const tail0=l.div(a).add(new F(1).div(a.pow(2))).mul(power).div(factor).hi;
 const tail1=l.square().div(a).add(l.mul(2).div(a.pow(2))).add(new F(2).div(a.pow(3))).mul(power).div(factor.pow(2)).hi;
 return {log_derivative:sum0.neg().add(new I(tail0.neg(),0)),log_derivative_prime:sum1.add(new I(0,tail1)),retained_positive_sum:sum0,retained_derivative_sum:sum1,tail0,tail1,cutoff,retained_prime_count:primes.length,last_retained_prime:primes.at(-1)};
}
function sourceJet(s,{cutoff=PRIME_CUTOFF,shift=PSI_SHIFT,terms=PSI_TERMS}={}){
 boundedInt(s,3,501);if(s%2!==1)throw new RangeError('This source adapter requires odd original s');
 const w=new F(2*s-1,2),t=w.pow(2),psi=digamma(new F(s,2),shift,terms),psi1=trigammaOdd(s),euler=eulerJet(s,cutoff),cc=constants();
 const pole=new F(1,s).add(new F(1,s-1)),polePrime=new F(-1,s*s).sub(new F(1,(s-1)*(s-1)));
 const gamma=psi.value.sub(cc.logpi).div(2),gammaPrime=psi1.div(4),Q=gamma.add(pole).add(euler.log_derivative),Qprime=gammaPrime.add(polePrime).add(euler.log_derivative_prime);
 const q=Q.div(w.mul(2)),derivative=Qprime.mul(w).sub(Q).div(w.pow(3).mul(4));
 return {s,w,t,q,derivative,Q,Qprime,components:{pole,polePrime,gamma,gammaPrime,psi,psi1,euler},protocol:{log_terms:LOG_TERMS,atan_terms:ATAN_TERMS,rounding_decimal_places:110,source:'completed xi; all omitted Euler terms enclosed; no RH or zero list'}};
}
module.exports=Object.freeze({atanhLog,logRational,atanReciprocal,constants,bernoulli,digamma,trigammaOdd,primesThrough,eulerJet,sourceJet});
