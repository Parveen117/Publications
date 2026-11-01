'use strict';
// LOT-1. All scientific arithmetic stays in the frozen native engines.
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const O=require('../native-prime-observation-closure/observation_adapter.cjs');
const W=require('../native-prime-innovation-window/window_adapter.cjs');
const L=require('../native-infinite-layer-cut/layer_adapter.cjs');
const C=require('../native-layer-continuation-window/continuation_adapter.cjs');
const D=require('../native-layer-depth-compatibility/depth_adapter.cjs');
const {F}=N,{Interval:I}=O;
function strict(raw,n){
 const f=L.forms(raw,n);
 if(!Object.values(f.factors).every(t=>t.positive))throw new RangeError('Strict acquired source pairings required');
 return f;
}
function inclusion(n,m){
 if(!Number.isSafeInteger(n)||!Number.isSafeInteger(m)||n<1||m<n)throw new RangeError('Positive nested carrier dimensions required');
 return N.matrix(Array.from({length:m},(_,i)=>Array.from({length:n},(_,j)=>i===j?1:0)));
}
function weightedCut(raw,n){
 strict(raw,n);const f=L.exactForms(raw,n),id=N.identity(n),zero=N.zeros(n),other=N.sub(id,f.Y);
 const G=L.block(f.A,zero,zero,f.B),cut=L.block(f.Y,other,f.Y,other),Z=L.block(id,zero,zero,N.scale(id,-1)),wave=N.mul(N.sub(N.scale(cut,2),N.identity(2*n)),Z);
 return {...f,G,cut,Z,wave};
}
function adjoint(T,domainMetric,rangeMetric){return N.mul(N.inverse(domainMetric),N.mul(N.dagger(T),rangeMetric));}
function exactTransport(raw,n,m=n+1){
 if(m<=n)throw new RangeError('Strict depth extension required');
 const old=weightedCut(raw,n),next=weightedCut(raw,m),E=inclusion(n,m),zero=N.zeros(m,n),V=L.block(E,zero,zero,E),star=adjoint(V,old.G,next.G),d=N.matrix([...Array.from({length:n},(_,j)=>[j===n-1?1:0]),...Array.from({length:n},(_,j)=>[j===n-1?-1:0])]);
 const hcol=N.matrix(raw.slice(n,2*n).map(x=>[F.of(x)])),alpha=N.mul(N.inverse(old.H),hcol),monic=N.matrix([...alpha.map(row=>[row[0].neg()]),[1],...Array.from({length:m-n-1},()=>[0])]),w=N.matrix([...monic,...monic]),norm=N.mul(N.dagger(monic),N.mul(next.H,monic))[0][0].rad;
 const K=N.scale(N.mul(N.inverse(old.G),N.mul(d,N.dagger(d))),norm);
 return {n,m,old,next,E,V,star,monic,w,norm,K,p: norm.mul(N.mul(N.dagger(d),N.mul(N.inverse(old.G),d))[0][0].rad),boundary:N.scale(N.mul(w,N.mul(N.dagger(w),next.G)),new F(1).div(norm)),compressed:N.mul(star,N.mul(next.cut,V)),compressedWave:N.mul(star,N.mul(next.wave,V))};
}
function sourceTransport(raw,n){
 const old=strict(raw,n);strict(raw,n+1);
 const boundary=D.boundary(raw,n),window=C.evenWindow(raw,n),ell=Array.from({length:n},(_,j)=>new I(j===n-1?1:0)),dualA=W.solve(old.factors.A,ell),dualB=W.solve(old.factors.B,ell),inverseWidth=dualA[n-1].add(dualB[n-1]),p=boundary.norm.mul(inverseWidth),windowPosition=boundary.norm.div(window.width),oneMinus=new I(1).sub(p);
 if(!p.positive()||!oneMinus.positive())throw new RangeError('Strict canonical boundary position not certified');
 const cutDefect=p.mul(oneMinus),waveLoss=cutDefect.mul(4),waveDiscrepancy=p.mul(4);
 return {n,boundaryNorm:boundary.norm,windowWidth:window.width,inverseWidth,position:p,windowPosition,cutDefect,waveLoss,waveDiscrepancy,dualA,dualB};
}
module.exports=Object.freeze({strict,inclusion,weightedCut,adjoint,exactTransport,sourceTransport});
