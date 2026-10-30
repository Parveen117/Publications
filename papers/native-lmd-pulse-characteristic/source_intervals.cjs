'use strict';
// Adapt the frozen NRQC pilot's interval method; arithmetic remains in its
// existing native rational field. No floating special functions are called.
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const {F}=N,f=F.of,assert=require('node:assert/strict');
const SCALE=1n<<80n,K=2048,TERMS=80,BITS=28;
const iv=(lo,hi=lo)=>({lo:f(lo),hi:f(hi)});
function down(x){x=f(x);assert(x.n>=0n);return new F(x.n*SCALE/x.d,SCALE);}
function up(x){x=f(x);assert(x.n>=0n);return new F((x.n*SCALE+x.d-1n)/x.d,SCALE);}
const add=(a,b)=>iv(down(a.lo.add(b.lo)),up(a.hi.add(b.hi)));
const mul=(a,b)=>iv(down(a.lo.mul(b.lo)),up(a.hi.mul(b.hi)));
function zeta(n){
 let lo=0n,hi=0n;
 for(let k=1;k<=K;k++){const d=BigInt(k)**BigInt(n);lo+=SCALE/d;hi+=(SCALE+d-1n)/d;}
 return iv(down(new F(lo,SCALE).add(f(1).div(f(n-1).mul(f(K+1).pow(n-1))))),
           up(new F(hi,SCALE).add(f(1).div(f(n-1).mul(f(K).pow(n-1))))));
}
function expPositive(x){
 assert(x.lo.n>=0n&&x.hi.le(4));
 let term=iv(1),sum=iv(1);
 for(let n=1;n<=120;n++){term=mul(term,mul(x,iv(f(1).div(n))));sum=add(sum,term);}
 const next=up(term.hi.mul(x.hi).div(121));
 return iv(sum.lo,up(sum.hi.add(next.div(f(1).sub(x.hi.div(122))))));
}
function isqrt(n){assert(n>=0n);if(n<2n)return n;let x=n,y=(x+1n)/2n;while(y<x){x=y;y=(x+n/x)/2n;}return x;}
function sqrtFloor(x){x=f(x);assert(x.n>0n);const s=1n<<BigInt(BITS);return new F(isqrt(x.n*s*s/x.d),s);}
function gridFloor(x){x=f(x);assert(x.n>=0n);const s=1n<<BigInt(BITS);return new F(x.n*s/x.d,s);}
function fixed(x,places=9){x=f(x);const s=10n**BigInt(places),neg=x.n<0n,n=(neg?-x.n:x.n)*s/x.d;return(neg?'-':'')+String(n/s)+'.'+String(n%s).padStart(places,'0');}
function bounds(a){return{lower:a.lo.toString(),upper:a.hi.toString()};}
function maxError(x,a){const lo=x.sub(a.lo).abs(),hi=x.sub(a.hi).abs();return lo.le(hi)?hi:lo;}
function cutoffSecond(h,n){
 // Constant Laurent coefficient of exp(h sum_{j=1}^n cos(j theta)/j).
 // All Laurent coefficients of the exponent are nonnegative. Rounding each
 // Taylor-convolution contribution down/up therefore encloses the result.
 const eh=expPositive(iv(h));if(n===0)return eh;
 let low=new Map([[0,SCALE]]),high=new Map([[0,SCALE]]),sum=iv(1),x=f(0),majorant=f(1);
 for(let j=1;j<=n;j++)x=x.add(h.div(j));
 const count=48;
 for(let m=1;m<=count;m++){
  const nextLow=new Map(),nextHigh=new Map();
  for(const [offset,v] of low)for(let j=1;j<=n;j++)for(const sign of [-1,1]){
   const key=offset+sign*j,den=h.d*BigInt(2*j*m),num=v*h.n;
   nextLow.set(key,(nextLow.get(key)||0n)+num/den);
  }
  for(const [offset,v] of high)for(let j=1;j<=n;j++)for(const sign of [-1,1]){
   const key=offset+sign*j,den=h.d*BigInt(2*j*m),num=v*h.n;
   nextHigh.set(key,(nextHigh.get(key)||0n)+(num+den-1n)/den);
  }
  low=nextLow;high=nextHigh;
  sum=add(sum,iv(new F(low.get(0)||0n,SCALE),new F(high.get(0)||0n,SCALE)));
  majorant=majorant.mul(x).div(m);
 }
 assert(x.le(count+1));
 const tail=majorant.mul(x).div(count+1).div(f(1).sub(x.div(count+2)));
 sum.hi=up(sum.hi.add(tail));return mul(eh,sum);
}
function build(){
 const start=process.hrtime.bigint(),zs=new Map();for(let n=2;n<=TERMS;n++)zs.set(n,zeta(n));
 const cache=new Map();
 function massSecond(g){
  g=f(g);assert(g.n>0n&&!f(1).le(g));let log=iv(g);
  for(let n=2;n<=TERMS;n++)log=add(log,mul(zs.get(n),iv(g.pow(n).mul(f(1).sub(f(1).div(f(2).pow(n-1)))).div(n))));
  log.hi=up(log.hi.add(f(2).mul(g.pow(TERMS+1)).div(f(TERMS+1).mul(f(1).sub(g)))));
  const result=expPositive(log);cache.set(g.toString(),bounds(result));return result;
 }
 const cases=['1/8','1/4','3/8'].map(hs=>{
  const h=f(hs),I=massSecond(h),I2=massSecond(h.mul(2)),B=iv(I.lo.sub(1),I.hi.sub(1)),A=h.div(4),c=h.div(2);
  const alpha=sqrtFloor(A),beta=gridFloor(c.div(alpha)),gamma=sqrtFloor(B.lo.add(B.hi).div(2).sub(beta.pow(2)));
  const ah=alpha.pow(2),bh=beta.pow(2).add(gamma.pow(2)),ch=alpha.mul(beta);
  const errors={A:ah.sub(A).abs(),B:maxError(bh,B),c:ch.sub(c).abs()};
  for(const x of Object.values(errors))assert(x.le('1/1000000'));
  assert(!A.mul(B.lo).sub(c.pow(2)).le(0));
  const gap=iv(I2.lo.sub(I.hi),I2.hi.sub(I.lo));assert(!gap.lo.le(0));
  return{h:hs,A,c,I:bounds(I),I2:bounds(I2),B:bounds(B),alpha,beta,gamma,
   encoded:{A:ah,B:bh,c:ch},absolute_error_bounds:errors,
   determinant_lower:A.mul(B.lo).sub(c.pow(2)),variance_gap_at_t1:bounds(gap),
   display:{actual_variance_t1:fixed(h.div(2).add(I.lo.add(I.hi).div(2)).sub(1)),
    naive_variance_t1:fixed(h.div(2).add(I2.lo.add(I2.hi).div(2)).sub(1)),
    gap_t1:fixed(gap.lo.add(gap.hi).div(2))}};
 });
 // Certified first-time response of the bounded sine reading at k=2:
 // exp(-h)*(2/h*sinh(h)-2), with sinh computed as a positive series.
 const bounded_response=cases.map(row=>{
  const h=f(row.h),eh=expPositive(iv(h));let term=h,sum=h;
  for(let n=1;n<=40;n++){term=term.mul(h.pow(2)).div((2*n)*(2*n+1));sum=sum.add(term);}
  const next=term.mul(h.pow(2)).div(82*83),tail=next.div(f(1).sub(h.pow(2).div(84*85)));
  const diff=iv(f(2).mul(sum).div(h).sub(2),f(2).mul(sum.add(tail)).div(h).sub(2));
  const gap=iv(diff.lo.div(eh.hi),diff.hi.div(eh.lo));assert(!gap.lo.le(0));
  return{h:row.h,k:'2',sine_response_slope_gap:bounds(gap),display_gap:fixed(gap.lo.add(gap.hi).div(2))};
 });
 const cutoff_moments=cases.map(row=>({h:row.h,rows:[0,1,2,4,8,16].map(n=>{
  const value=cutoffSecond(f(row.h),n),gap=iv(f(row.I.lower).sub(value.hi),f(row.I.upper).sub(value.lo));
  assert(!gap.lo.le(0));return{N:n,I_N:bounds(value),mass_L2_error_squared:bounds(gap),
   display_I_N:fixed(value.lo.add(value.hi).div(2)),display_L2_error_squared:fixed(gap.lo.add(gap.hi).div(2))};
 })}));
 return{scope:'conditional source moments plus finite encoding coefficients; not a continuum simulator',
  method:'native outward log-Gamma/zeta/exponential intervals',interval_bits:80,coefficient_bits:BITS,
  zeta_terms:K,log_gamma_terms:TERMS,exp_terms:120,times:['-1','-1/2','0','1/2','1'],
  cases,bounded_response,cutoff_moments,cutoff_taylor_terms:48,integrals:Object.fromEntries(cache),
  source_preparation_seconds:Number(process.hrtime.bigint()-start)/1e9};
}
if(require.main===module)process.stdout.write(JSON.stringify(build(),null,2)+'\n');
module.exports={build,fixed};
