'use strict';
const {F}=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const O=require('../native-prime-observation-closure/observation_adapter.cjs');
const W=require('../native-prime-innovation-window/window_adapter.cjs');
const L=require('../native-infinite-layer-cut/layer_adapter.cjs');
const {Interval:I}=O,one=new F(1);
const dot=(a,b)=>a.reduce((s,v,j)=>s.add(I.of(v).mul(b[j])),new I(0));
const det=f=>{if(!f.positive)throw new RangeError('Certified positive factor required');return f.d.reduce((s,v)=>s.mul(v),new I(1));};
function record(raw,m){
 if(!Number.isSafeInteger(m)||m<0||raw.length<=m)throw new RangeError('Available nonnegative record length required');
 const c=raw.map(I.of);if(!c[0].lo.eq(1)||!c[0].hi.eq(1))throw new RangeError('Unit source normalization required');
 const words=Array.from({length:m+1},(_,j)=>L.difference(c,j,m-j)),weights=words.map((v,j)=>v.mul(L.choose(m,j)));
 if(!words.every(v=>v.positive()))throw new RangeError('Strict source record probabilities required');
 return {m,words,weights,total:weights.reduce((s,v)=>s.add(v),new I(0))};
}
function coarsen(weights){
 if(!Array.isArray(weights)||weights.length<2)throw new RangeError('At least a one-bit record required');
 const n=weights.length-1,w=weights.map(I.of);
 return Array.from({length:n},(_,j)=>w[j].mul(new F(n-j,n)).add(w[j+1].mul(new F(j+1,n))));
}
function adjacent(raw,m){
 if(!Number.isSafeInteger(m)||m<1)throw new RangeError('Positive shorter record length required');
 const lower=record(raw,m),upper=record(raw,m+1),coarse=coarsen(upper.weights),joint=[];let meanSquare=new I(0);
 for(let j=0;j<=m+1;j++){
  for(const k of [j-1,j]){
   if(k<0||k>m)continue;
   const conditional=k===j?new F(m+1-j,m+1):new F(j,m+1),mass=upper.weights[j].mul(conditional),difference=new F(k,m).sub(new F(j,m+1));
   meanSquare=meanSquare.add(mass.mul(difference.pow(2)));joint.push({old_count:k,new_count:j,conditional,mass});
  }
 }
 const expected=I.of(raw[1]).sub(raw[2]).div(m*(m+1));
 return {m,lower,upper,coarse,joint,meanSquare,expected};
}
function nested(raw,n){
 if(!Number.isSafeInteger(n)||n<1||raw.length<2*n+2)throw new RangeError('Two complete nested source banks required');
 const small=L.forms(raw,n),large=L.forms(raw,n+1);
 if(![small,large].every(f=>['H','A','B'].every(k=>f.factors[k].positive)))throw new RangeError('Strict nested source forms required');
 return {c:raw.map(I.of),small,large,n};
}
function boundary(raw,n){
 const {c,small,large}=nested(raw,n),h=c.slice(n,2*n),alpha=W.solve(small.factors.H,h),monic=[...alpha.map(v=>v.neg()),new I(1)],norm=large.factors.H.d[n],schurNorm=c[2*n].sub(W.inverseForm(small.factors.H,h)),previousNorm=small.factors.H.d[n-1],beta=norm.div(previousNorm);
 const orthogonality=Array.from({length:n},(_,i)=>dot(c.slice(i,i+n+1),monic));
 return {n,monic,norm,schurNorm,previousNorm,beta,orthogonality,detH:det(small.factors.H),detHNext:det(large.factors.H)};
}
function innovation(raw,n,lambda){
 lambda=F.of(lambda);if(lambda.le(-1))throw new RangeError('Real displacement lambda>-1 required');
 const {c,small,large}=nested(raw,n),pencil=f=>O.intervalLDL(f.H.map((row,i)=>row.map((v,j)=>v.add(f.A[i][j].mul(lambda))))),p=pencil(small),pNext=pencil(large);
 if(!p.positive||!pNext.positive)throw new RangeError('Positive response pencils required');
 const h=c.slice(n,2*n),k=c.slice(n+1,2*n+1),v=h.map((x,j)=>x.add(k[j].mul(lambda))),residual=c[n].sub(dot(v,W.solve(p,small.u))),schur=c[2*n].add(c[2*n+1].mul(lambda)).sub(W.inverseForm(p,v)),d=det(p),dNext=det(pNext),hNext=det(large.factors.H);
 const norm=large.factors.H.d[n],oldRatio=det(small.factors.H).div(d),newRatio=hNext.div(dNext);
 const formula=norm.mul(oldRatio).mul(newRatio).mul(lambda.pow(2*n)),residualFormula=norm.mul(oldRatio).mul(lambda.neg().pow(n)),before=W.inverseForm(p,small.u),after=W.inverseForm(pNext,large.u),difference=after.sub(before);
 return {n,lambda,before,after,difference,residual,residualFormula,schur,schurRatio:dNext.div(d),formula,norm,oldRatio,newRatio,detPencil:d,detPencilNext:dNext,detHNext:hNext};
}
function falling(n,r){let value=1n;for(let j=0;j<r;j++)value*=BigInt(n-j);return value;}
function reading(weights,r){
 const m=weights.length-1;if(!Number.isSafeInteger(r)||r<0||r>m)throw new RangeError('Available factorial reading required');
 return weights.reduce((s,v,j)=>s.add(I.of(v).mul(new F(falling(j,r),falling(m,r)))),new I(0));
}
module.exports=Object.freeze({record,coarsen,adjacent,boundary,innovation,reading,falling,det});
