'use strict';
const {F}=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const {Interval:I}=require('../native-prime-observation-closure/observation_adapter.cjs');
const L=require('../native-infinite-layer-cut/layer_adapter.cjs');
const {factorial}=require('../native-layer-wave-propagation/propagation_adapter.cjs');
const R=new F(15,16);
function index(n,min=0){if(!Number.isSafeInteger(n)||n<min)throw new RangeError('Safe nonnegative index required');return n;}
function majorant(M,r=R){M=F.of(M);r=F.of(r);if(M.le(0)||r.le(0))throw new RangeError('Positive coefficient majorant required');return {M,R:r,C:M.div(r),A:r.add(1).mul(4).div(r.pow(2))};}
function sourceMajorant(source){
 const w=F.of(source.s).sub('1/2'),a=F.of(source.a),q0=I.of(source.q0);
 if(w.le(4)||!a.eq(w.pow(2))||!q0.positive())throw new RangeError('LAA source anchor and positive q0 required');
 const Q=new F(41).div(w.mul(2)).add('1/4');
 return {...majorant(Q.div(q0.lo)),w,a,q0,Q};
}
function coefficient(j,t){index(j);return F.of(t).pow(2).mul(-4).pow(j).div(factorial(2*j));}
function fullBound(j,b){index(j);return b.C.mul(b.A.div(4).pow(j));}
function tail(J,t,b){
 index(J,1);const x=b.A.mul(F.of(t).pow(2)),ratio=x.div((2*J+1)*(2*J+2));
 if(new F(1).le(ratio))throw new RangeError('Increase cutoff until the geometric ratio is below one');
 return {J,x,ratio,radius:b.C.mul(x.pow(J)).div(factorial(2*J)).div(new F(1).sub(ratio))};
}
function missing(K,j,b){
 index(K,1);index(j);let total=new F(0);
 for(let ell=Math.max(0,K-j);ell<=j;ell++)total=total.add(new F(L.choose(j,ell)).div(b.R.pow(ell)));
 return b.M.mul(total).div(b.R.pow(j+1));
}
function known(raw,t){
 const K=raw.length-1;index(K,1);let value=new I(0);const terms=[];
 for(let j=0;j<K;j++){
  let moment=new I(0);
  for(let ell=0;ell<=Math.min(j,K-j-1);ell++)moment=moment.add(I.of(raw[j+1+ell]).mul(new F(L.choose(j,ell)).mul(ell%2?-1:1)));
  const weight=coefficient(j,t),term=moment.mul(weight);value=value.add(term);terms.push({j,moment,weight,term});
 }
 return {K,t:F.of(t),value,terms};
}
function kernel(n,t){
 index(n,1);let value=new F(0);
 for(let j=Math.ceil((n-1)/2);j<=n-1;j++){const ell=n-j-1;value=value.add(coefficient(j,t).mul(new F(L.choose(j,ell))).mul(ell%2?-1:1));}
 return value;
}
function checkRetained(raw,b){
 if(!Array.isArray(raw)||raw.length<2)throw new RangeError('At least c0 and c1 required');
 if(!I.of(raw[0]).contains(1))throw new RangeError('Normalized c0 required');
 // This is a consistency guard, not a proof of the unknown coefficient bound.
 for(let n=0;n<raw.length;n++)if(!I.of(raw[n]).absUpper().le(b.M.div(b.R.pow(n))))throw new RangeError('Retained interval exceeds supplied coefficient budget');
}
function enclosure(raw,t,b,J=20){
 checkRetained(raw,b);const K=raw.length-1;index(J,K);
 const p=known(raw,t),rest=tail(J,t,b),missingTerms=[];let budget=new F(0);
 for(let j=0;j<J;j++){const bound=missing(K,j,b),term=coefficient(j,t).abs().mul(bound);budget=budget.add(term);if(!term.zero())missingTerms.push({j,bound,term});}
 const radius=budget.add(rest.radius),value=p.value.add(new I(radius.neg(),radius));
 return {K,J,t:F.of(t),partial:p.value,missingTerms,finiteMissingBudget:budget,tail:rest,radius,enclosure:value};
}
function coarse(raw,t,b){
 checkRetained(raw,b);const J=Math.floor(raw.length/2);let partial=new I(0);
 for(let j=0;j<J;j++)partial=partial.add(L.difference(raw,j+1,j).mul(coefficient(j,t)));
 const rest=tail(J,t,b);return {J,t:F.of(t),partial,tail:rest,radius:rest.radius,enclosure:partial.add(new I(rest.radius.neg(),rest.radius))};
}
module.exports=Object.freeze({R,majorant,sourceMajorant,coefficient,fullBound,tail,missing,known,kernel,checkRetained,enclosure,coarse});
