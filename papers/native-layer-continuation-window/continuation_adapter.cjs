'use strict';
// LCW-1: sharp next-layer windows in the inherited rational/interval engines.
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const O=require('../native-prime-observation-closure/observation_adapter.cjs');
const W=require('../native-prime-innovation-window/window_adapter.cjs');
const L=require('../native-infinite-layer-cut/layer_adapter.cjs');
const {F}=N,{Interval:I}=O;
function evenWindow(raw,n){
 const f=L.forms(raw,n),c=raw.map(I.of);
 if(!f.factors.A.positive||!f.factors.B.positive||!f.factors.H.positive)throw new RangeError('Certified strict old layer forms required');
 const h=c.slice(n,2*n),C=Array.from({length:n-1},(_,i)=>Array.from({length:n-1},(_,j)=>c[i+j+1].sub(c[i+j+2]))),d=Array.from({length:n-1},(_,i)=>c[n+i].sub(c[n+i+1]));
 const cf=n===1?null:O.intervalLDL(C);
 if(cf&&!cf.positive)throw new RangeError('Mixed old form not certified positive');
 const lower=W.inverseForm(f.factors.H,h),upper=c[2*n-1].sub(n===1?0:W.inverseForm(cf,d)),width=upper.sub(lower);
 if(!width.positive())throw new RangeError('Even window not certified open');
 return {order:n,lower,upper,width,h,C,d,mixed_pivots:cf?cf.d:[]};
}
function oddWindow(raw,n,x){
 const f=L.forms(raw,n),c=raw.map(I.of);x=I.of(x);
 if(!f.factors.A.positive||!f.factors.B.positive)throw new RangeError('Certified strict old layer forms required');
 const v=[...c.slice(n+1,2*n),x],w=Array.from({length:n},(_,i)=>c[n+i].sub(i===n-1?x:c[n+i+1]));
 const lower=W.inverseForm(f.factors.A,v),upper=x.sub(W.inverseForm(f.factors.B,w));
 return {order:n,x,lower,upper,width:upper.sub(lower),v,w};
}
function acceptance(raw,n,x,y){
 const even=evenWindow(raw,n),odd=oddWindow(raw,n,x);x=I.of(x);y=I.of(y);
 const evenLower=x.sub(even.lower),evenUpper=even.upper.sub(x),oddLower=y.sub(odd.lower),oddUpper=odd.upper.sub(y);
 const evenFraction=evenLower.div(even.width),factorizedWidth=evenLower.mul(evenUpper).div(even.width);
 const oddFraction=odd.width.positive()?oddLower.div(odd.width):null;
 return {order:n,x,y,even,odd,evenLower,evenUpper,oddLower,oddUpper,evenFraction,oddFraction,factorizedWidth,accepted:[evenLower,evenUpper,oddLower,oddUpper].every(z=>z.positive())};
}
const quad=(M,b)=>{const v=N.matrix(b.map(z=>[z]));return N.mul(N.dagger(v),N.mul(N.inverse(M),v))[0][0].rad;};
function exactWindow(raw,n,x){
 const c=raw.map(F.of),f=L.exactForms(c,n);x=F.of(x);
 const h=c.slice(n,2*n),C=Array.from({length:n-1},(_,i)=>Array.from({length:n-1},(_,j)=>c[i+j+1].sub(c[i+j+2]))),d=Array.from({length:n-1},(_,i)=>c[n+i].sub(c[n+i+1]));
 const lower=quad(f.H,h),upper=c[2*n-1].sub(n===1?0:quad(C,d)),v=[...c.slice(n+1,2*n),x],w=Array.from({length:n},(_,i)=>c[n+i].sub(i===n-1?x:c[n+i+1]));
 const oddLower=quad(f.A,v),oddUpper=x.sub(quad(f.B,w));
 return {lower,upper,width:upper.sub(lower),x,oddLower,oddUpper,oddWidth:oddUpper.sub(oddLower)};
}
module.exports=Object.freeze({evenWindow,oddWindow,acceptance,exactWindow});
