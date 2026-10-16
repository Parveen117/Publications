'use strict';
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const O=require('../native-prime-observation-closure/observation_adapter.cjs');
const W=require('../native-prime-innovation-window/window_adapter.cjs');
const {F}=N,{Interval:I}=O;
function forms(moments,order){
 if(!Number.isSafeInteger(order)||order<1||moments.length<2*order)throw new RangeError('Two full moment layers per requested order required');
 const c=moments.map(I.of);if(!c[0].lo.eq(1)||!c[0].hi.eq(1))throw new RangeError('Exact unit source normalization required');
 const H=Array.from({length:order},(_,i)=>Array.from({length:order},(_,j)=>c[i+j])),A=H.map((r,i)=>r.map((_,j)=>c[i+j+1])),B=H.map((r,i)=>r.map((x,j)=>x.sub(A[i][j])));
 return {order,H,A,B,u:c.slice(0,order),factors:{H:O.intervalLDL(H),A:O.intervalLDL(A),B:O.intervalLDL(B)}};
}
function response(moments,order,lambda){
 lambda=F.of(lambda);if(lambda.le(-1))throw new RangeError('Positive original t requires lambda > -1');
 const f=forms(moments,order);if(!f.factors.A.positive||!f.factors.B.positive||!f.factors.H.positive)throw new RangeError('Both layer forms must be certified positive');
 const pencil=f.H.map((r,i)=>r.map((x,j)=>x.add(f.A[i][j].mul(lambda)))),factor=O.intervalLDL(pencil);
 if(!factor.positive)throw new RangeError('Layer response pencil not certified positive');
 return {value:W.inverseForm(factor,f.u),lambda,order,pencil_pivots:factor.d};
}
function exactForms(moments,order){
 if(!Number.isSafeInteger(order)||order<1||moments.length<2*order)throw new RangeError('Insufficient exact moments');
 const c=moments.map(F.of);if(!c[0].eq(1))throw new RangeError('Exact unit source normalization required');
 const H=N.matrix(Array.from({length:order},(_,i)=>Array.from({length:order},(_,j)=>c[i+j]))),A=N.matrix(Array.from({length:order},(_,i)=>Array.from({length:order},(_,j)=>c[i+j+1]))),B=N.sub(H,A),Y=N.mul(N.inverse(H),A);
 return {H,A,B,Y,u:N.matrix(c.slice(0,order).map(x=>[x]))};
}
function block(a,b,c,d){return [...a.map((r,i)=>[...r,...b[i]]),...c.map((r,i)=>[...r,...d[i]])];}
function cutWave(Y,S,H=N.identity(Y.length)){
 Y=N.matrix(Y);S=N.matrix(S);H=N.matrix(H);const id=N.identity(Y.length),zero=N.zeros(Y.length);
 const positive=M=>O.intervalLDL(M.map(row=>row.map(x=>{if(!x.turn.zero())throw new TypeError('Real acquired pairing required');return x.rad;}))).positive;
 if(!N.equal(H,N.dagger(H))||!N.equal(N.mul(N.dagger(Y),H),N.mul(H,Y))||!N.equal(N.mul(N.dagger(S),H),N.mul(H,S)))throw new RangeError('Acquired self-adjoint source operators required');
 if(![H,N.mul(H,Y),N.mul(H,N.sub(id,Y)),N.mul(H,S)].every(positive))throw new RangeError('Strict finite contraction and positive square root required');
 if(!N.equal(N.mul(Y,S),N.mul(S,Y))||!N.equal(N.mul(S,S),N.mul(Y,N.sub(id,Y))))throw new RangeError('Commuting source square root required');
 const cut=block(Y,S,S,N.sub(id,Y)),reflection=N.sub(N.scale(cut,2),N.identity(2*Y.length)),fixedReflection=block(id,zero,zero,N.scale(id,-1)),wave=N.mul(reflection,fixedReflection);
 return {cut,reflection,fixedReflection,wave};
}
function choose(n,k){let v=1n;for(let j=1;j<=k;j++)v=v*BigInt(n-j+1)/BigInt(j);return v;}
function difference(moments,r,k){
 if(!Number.isSafeInteger(r)||!Number.isSafeInteger(k)||r<0||k<0||r+k>=moments.length)throw new RangeError('Complete nonnegative layer indices required');
 let d=new I(0);for(let j=0;j<=k;j++)d=d.add(I.of(moments[r+j]).mul(choose(k,j)*(j%2?-1n:1n)));return d;
}
function exactDifference(moments,r,k){
 if(!Number.isSafeInteger(r)||!Number.isSafeInteger(k)||r<0||k<0||r+k>=moments.length)throw new RangeError('Complete nonnegative layer indices required');
 let d=new F(0);for(let j=0;j<=k;j++)d=d.add(F.of(moments[r+j]).mul(choose(k,j)*(j%2?-1n:1n)));return d;
}
module.exports=Object.freeze({forms,response,exactForms,block,cutWave,choose,difference,exactDifference});
