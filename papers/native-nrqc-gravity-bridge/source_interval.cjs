'use strict';
// Outward rational enclosures using the frozen native arithmetic engine.
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const {F}=N, f=F.of, assert=require('node:assert/strict');
const SCALE=1n<<80n, K=1024, TERMS=80;
function down(x){x=f(x);assert(x.n>=0n);return new F(x.n*SCALE/x.d,SCALE);}
function up(x){x=f(x);assert(x.n>=0n);return new F((x.n*SCALE+x.d-1n)/x.d,SCALE);}
const interval=(lo,hi=lo)=>({lo:f(lo),hi:f(hi)});
const add=(a,b)=>interval(down(a.lo.add(b.lo)),up(a.hi.add(b.hi)));
const mul=(a,b)=>interval(down(a.lo.mul(b.lo)),up(a.hi.mul(b.hi)));
function zeta(n){
 let lower=0n,upper=0n;
 for(let k=1;k<=K;k++){
  const d=BigInt(k)**BigInt(n);
  lower+=SCALE/d; upper+=(SCALE+d-1n)/d;
 }
 // Integral test for the tail after k=K, in both directions.
 return interval(down(new F(lower,SCALE).add(f(1).div(f(n-1).mul(f(K+1).pow(n-1))))),
                 up(new F(upper,SCALE).add(f(1).div(f(n-1).mul(f(K).pow(n-1))))));
}
function expInterval(x){
 assert(x.lo.n>=0n && x.hi.le(4));
 let term=interval(1),sum=interval(1);
 for(let n=1;n<=120;n++){term=mul(term,mul(x,interval(f(1).div(n))));sum=add(sum,term);}
 const next=up(term.hi.mul(x.hi).div(121));
 const tail=up(next.div(f(1).sub(x.hi.div(122))));
 return interval(sum.lo,up(sum.hi.add(tail)));
}
function isqrt(n){assert(n>=0n);if(n<2n)return n;let x=n,y=(x+1n)/2n;while(y<x){x=y;y=(x+n/x)/2n;}return x;}
function sqrtFloor(x,bits=22){x=f(x);assert(x.n>0n);const s=1n<<BigInt(bits);return new F(isqrt(x.n*s*s/x.d),s);}
function fixed(x,places=12){x=f(x);const s=10n**BigInt(places),neg=x.n<0n;let n=(neg?-x.n:x.n)*s/x.d;return (neg?'-':'')+String(n/s)+'.'+String(n%s).padStart(places,'0');}
function bounds(x){return {lower:x.lo.toString(),upper:x.hi.toString(),display_lower:fixed(x.lo),display_upper:fixed(x.hi.add('1/1000000000000'))};}
function build(){
 const begin=process.hrtime.bigint(), zs=new Map();
 for(let n=2;n<=TERMS;n++)zs.set(n,zeta(n));
 const cases=['1/4','1/2','3/4'].map(hs=>{
  const h=f(hs);let logC=interval(0);
  for(let n=2;n<=TERMS;n++){
   const weight=h.pow(n).mul(f(1).sub(f(1).div(f(2).pow(n-1)))).div(n);
   logC=add(logC,mul(zs.get(n),interval(weight)));
  }
  // zeta(n)<2 and 0<1-2^(1-n)<1.
  const tail=f(2).mul(h.pow(TERMS+1)).div(f(TERMS+1).mul(f(1).sub(h)));
  logC.hi=up(logC.hi.add(tail));
  const B=mul(expInterval(add(interval(h),logC)),interval(h.div(f(2).mul(f(2).sub(h)))));
  const A=h.mul(f(1).add(h.div(2)).pow(2)).div(4).add(h.pow(2).div(16)),c=h.div(4);
  const alpha=sqrtFloor(A), beta=down(c.div(alpha));
  // Quantize beta to the same declared 22-bit grid.
  const den=1n<<22n,betaHat=new F(beta.n*den/beta.d,den);
  const gamma=sqrtFloor(B.lo.add(B.hi).div(2).sub(betaHat.pow(2)));
  const AHat=alpha.pow(2),BHat=betaHat.pow(2).add(gamma.pow(2)),cHat=alpha.mul(betaHat);
  const errB1=BHat.sub(B.lo).abs(),errB2=BHat.sub(B.hi).abs();
  const errors={A:AHat.sub(A).abs(),B:errB1.le(errB2)?errB2:errB1,c:cHat.sub(c).abs()};
  for(const e of Object.values(errors))assert(e.le('1/100000'));
  assert(!A.mul(B.lo).sub(c.pow(2)).le(0));
  return {h:hs,A,c,B:bounds(B),alpha,beta:betaHat,gamma,
    encoded:{A:AHat,B:BHat,c:cHat},absolute_error_bounds:errors,
    determinant_lower:A.mul(B.lo).sub(c.pow(2)),
    leading_only_error_A:A.sub(h.div(4)),
    display:{A:fixed(A),B_mid:fixed(B.lo.add(B.hi).div(2)),c:fixed(c),
     error_A:fixed(errors.A),error_B:fixed(errors.B),error_c:fixed(errors.c)}};
 });
 return {method:'native rational outward log-Gamma series, zeta integral tails, exponential Taylor tail',
  h_range:'only the three declared points, dimensionless',coefficient_bits:22,
  fixed_interval_bits:80,zeta_terms:K,log_gamma_terms:TERMS,exp_terms:120,
  source_preparation_seconds:Number(process.hrtime.bigint()-begin)/1e9,cases};
}
if(require.main===module)process.stdout.write(JSON.stringify(build(),null,2)+'\n');
module.exports={build,fixed};
