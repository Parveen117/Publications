'use strict';
/** POC-1: exact native equal-reading centres and outward source-bank replay. */
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const O=require('../native-prime-observer-sectors/observer_adapter.cjs');
const {F,Cut}=N;
const int=x=>O.integer(x);
function real(x){x=Cut.of(x);if(!x.turn.zero())throw new TypeError('Real source data required');return x.rad;}
function realMatrix(A){A=N.matrix(A);A.forEach(r=>r.forEach(real));return A;}
function column(x){if(!Array.isArray(x)||x.length!==2)throw new TypeError('Two coordinates required');return realMatrix(x.map(v=>[v]));}
function quadratic(H,x){return N.mul(N.dagger(x),N.mul(H,x))[0][0];}
function centre(points,H=[[2,1],[1,5]]){
 if(!Array.isArray(points)||points.length!==3)throw new TypeError('Three source points required');
 H=realMatrix(H);if(H.length!==2||H[0].length!==2||!N.equal(H,N.dagger(H)))throw new TypeError('Symmetric 2x2 metric required');
 if(!new F(0).le(real(H[0][0]))||real(H[0][0]).zero()||!new F(0).le(real(O.det(H)))||O.det(H).zero())throw new RangeError('Positive metric required');
 const xs=points.map(column),q0=quadratic(H,xs[0]);
 const A=xs.slice(1).map(x=>N.scale(N.mul(N.dagger(N.sub(x,xs[0])),H),2)[0]);
 const b=xs.slice(1).map(x=>[quadratic(H,x).sub(q0)]);
 if(O.det(A).zero())throw new RangeError('Noncollinear points required');
 const x=N.mul(N.inverse(A),b),radii=xs.map(a=>quadratic(H,N.sub(x,a)));
 return {H,points:xs,A,b,x,radii,determinant:O.det(A),contract:'quadratic readings in the retained tangent metric; no prime modulus selected'};
}
function integerSystem(A,b){A=realMatrix(A);b=realMatrix(b);if(A.length!==2||A[0].length!==2||b.length!==2||b[0].length!==1||!O.integral(A)||!O.integral(b))throw new TypeError('Integral 2x2 system with one column required');return {A,b};}
function residues(A,b,n){({A,b}=integerSystem(A,b));n=int(n);if(n<2n||n>256n)throw new RangeError('Finite residue modulus 2..256 required');const out=[];for(let a=0n;a<n;a++)for(let c=0n;c<n;c++)if(A.every((row,i)=>(row[0].rad.n*a+row[1].rad.n*c-b[i][0].rad.n)%n===0n))out.push([a.toString(),c.toString()]);return out;}
function extended(a,b){let u=1n,v=0n,x=0n,y=1n;while(b){const q=a/b;[a,b]=[b,a-q*b];[u,x]=[x,u-q*x];[v,y]=[y,v-q*y];}return [a,u,v];}
const mod=(a,n)=>(a%n+n)%n;
function rationalMod(x,n){x=real(x);n=int(n);if(n<2n)throw new RangeError('Modulus at least two required');let [g,u]=extended(x.d,n);if(g!==1n)throw new RangeError('Denominator not invertible');return mod(x.n*u,n).toString();}
function uniqueResidue(A,b,n){({A,b}=integerSystem(A,b));n=int(n);if(n<2n)throw new RangeError('Modulus at least two required');const d=real(O.det(A));if(O.gcd(d.n,n)!==1n)throw new RangeError('Unique residue requires unit determinant');return N.mul(N.inverse(A),b).map(r=>rationalMod(r[0],n));}
function decimal(x){if(x instanceof F)return x;if(typeof x==='number'||typeof x==='bigint')return F.of(x);if(typeof x!=='string')throw new TypeError('Exact decimal or rational string required');if(x.includes('/'))return new F(x);const m=x.match(/^([+-]?)(\d+)(?:\.(\d*))?(?:[eE]([+-]?\d+))?$/);if(!m)throw new TypeError('Invalid decimal');const digits=m[2]+(m[3]||''),e=Number(m[4]||0)-(m[3]||'').length;if(!Number.isSafeInteger(e)||Math.abs(e)>2000)throw new RangeError('Decimal exponent bound');const sign=m[1]==='-'?-1n:1n;return e>=0?new F(sign*BigInt(digits)*10n**BigInt(e)):new F(sign*BigInt(digits),10n**BigInt(-e));}
const SCALE=10n**110n;
function floorDiv(a,b){let q=a/b;if(a<0n&&a%b)q--;return q;}
function low(x){return new F(floorDiv(x.n*SCALE,x.d),SCALE);}
function high(x){return new F(-floorDiv(-x.n*SCALE,x.d),SCALE);}
class Interval{
 constructor(lo,hi=lo){lo=decimal(lo);hi=decimal(hi);if(!lo.le(hi))throw new RangeError('Reversed interval');this.lo=low(lo);this.hi=high(hi);Object.freeze(this);}
 static of(x){return x instanceof Interval?x:Array.isArray(x)?new Interval(x[0],x[1]):new Interval(x);}
 add(x){x=Interval.of(x);return new Interval(this.lo.add(x.lo),this.hi.add(x.hi));}
 neg(){return new Interval(this.hi.neg(),this.lo.neg());}
 sub(x){return this.add(Interval.of(x).neg());}
 mul(x){x=Interval.of(x);const a=[this.lo.mul(x.lo),this.lo.mul(x.hi),this.hi.mul(x.lo),this.hi.mul(x.hi)];return new Interval(a.reduce((p,q)=>p.le(q)?p:q),a.reduce((p,q)=>p.le(q)?q:p));}
 square(){if(this.lo.le(0)&&new F(0).le(this.hi)){const a=this.lo.pow(2),b=this.hi.pow(2);return new Interval(0,a.le(b)?b:a);}return this.mul(this);}
 div(x){x=Interval.of(x);if(x.lo.le(0)&&new F(0).le(x.hi))throw new RangeError('Interval divisor contains zero');return this.mul(new Interval(new F(1).div(x.hi),new F(1).div(x.lo)));}
 contains(x){x=decimal(x);return this.lo.le(x)&&x.le(this.hi);}
 positive(){return !this.lo.le(0);}
 absUpper(){const a=this.lo.abs(),b=this.hi.abs();return a.le(b)?b:a;}
 toJSON(){return [this.lo.toString(),this.hi.toString()];}
}
function intervalLDL(raw){if(!Array.isArray(raw)||!raw.length||raw.some(r=>r.length!==raw.length))throw new TypeError('Square interval matrix required');const A=raw.map(r=>r.map(Interval.of)),n=A.length;
 for(let i=0;i<n;i++)for(let j=0;j<i;j++)if(!A[i][j].lo.eq(A[j][i].lo)||!A[i][j].hi.eq(A[j][i].hi))throw new TypeError('Symmetric interval enclosure required');
 const L=Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>new Interval(i===j?1:0))),d=[];
 for(let j=0;j<n;j++){let p=A[j][j];for(let k=0;k<j;k++)p=p.sub(L[j][k].square().mul(d[k]));if(!p.positive())return {positive:false,failedAt:j+1,pivot:p};d.push(p);for(let i=j+1;i<n;i++){let v=A[i][j];for(let k=0;k<j;k++)v=v.sub(L[i][k].mul(L[j][k]).mul(d[k]));L[i][j]=v.div(p);}}
 return {positive:true,L,d,minimumLower:d.map(x=>x.lo).reduce((a,b)=>a.le(b)?a:b)};
}
function response(rawL,rawJ,rawU,z){
 if(!Array.isArray(rawU)||!rawU.length||[rawL,rawJ].some(A=>!Array.isArray(A)||A.length!==rawU.length||A.some(r=>!Array.isArray(r)||r.length!==rawU.length)))throw new TypeError('Response dimensions must agree');
 z=decimal(z);if(z.le(0))throw new RangeError('Positive real response point required');
 const A=rawJ.map((r,i)=>r.map((x,j)=>Interval.of(x).add(Interval.of(rawL[i][j]).mul(z)))),factor=intervalLDL(A);
 if(!factor.positive)throw new RangeError('Response pencil not certified positive');
 const v=[];for(let i=0;i<rawU.length;i++){let q=Interval.of(rawU[i]);for(let j=0;j<i;j++)q=q.sub(factor.L[i][j].mul(v[j]));v.push(q);}
 let q=new Interval(0);for(let i=0;i<v.length;i++)q=q.add(v[i].square().div(factor.d[i]));return {value:q,factor};
}
function primeBank(source){if(source.dimension!==24||source.status!=='PASS_FINITE_SOURCE_AND_SCOPED_CONVERGENCE_CONTRACT')throw new TypeError('Pinned SE49 archive required');const primes=[3,5,7,11,13,17,19,23],indices=primes.map(p=>p-3);for(let j=0;j<primes.length;j++)if(!new F(source.nodes[indices[j]]).eq(new F(2*primes[j]-1,2).pow(2)))throw new Error('Prime node/source mismatch');return {primes,indices,nodes:indices.map(i=>source.nodes[i]),values:indices.map(i=>source.values[i]),L:indices.map(i=>indices.map(j=>source.L[i][j])),J:indices.map(i=>indices.map(j=>source.J[i][j]))};}
function schwarzProduct(ws,w,a=0){w=Cut.of(w);a=decimal(a);if(!a.le(w.rad)||a.eq(w.rad))throw new RangeError('Point must be inside half-plane');let product=new F(1);for(const raw of ws){const b=decimal(raw);if(b.le(a))throw new RangeError('Probe outside half-plane');const numerator=w.sub(b).norm2(),denominator=w.add(b.sub(a.mul(2))).norm2();product=product.mul(numerator.div(denominator));}return product;}
module.exports=Object.freeze({centre,residues,uniqueResidue,rationalMod,decimal,Interval,intervalLDL,response,primeBank,schwarzProduct});
