'use strict';
/** Exact finite operator calculus over the rational cut field Q + iota Q.
 * The scalar multiplication is defined by iota^2=-1; no floating point is used.
 * Matrix coordinates realize finite path composition, not a physical spacetime.
 */
function gcd(a,b) { a=a<0n?-a:a; b=b<0n?-b:b; while(b){[a,b]=[b,a%b];} return a; }
function integer(x) {
  if(typeof x==='bigint') return x;
  if(typeof x==='number' && Number.isSafeInteger(x)) return BigInt(x);
  if(typeof x==='string' && /^[-+]?\d+$/.test(x)) return BigInt(x);
  throw new TypeError('Use exact integers or rational strings, not floating-point inputs');
}
class F {
  constructor(n=0,d=1) {
    if(n instanceof F && d===1) {this.n=n.n;this.d=n.d;Object.freeze(this);return;}
    if(typeof n==='string' && n.includes('/') && d===1){const a=n.split('/');if(a.length!==2)throw new TypeError('Invalid fraction');[n,d]=a;}
    n=integer(n);d=integer(d);if(!d)throw new RangeError('Zero denominator');if(d<0n){n=-n;d=-d;}
    const k=gcd(n,d);this.n=n/k;this.d=d/k;Object.freeze(this);
  }
  static of(x){return x instanceof F?x:new F(x);}
  add(x){x=F.of(x);return new F(this.n*x.d+x.n*this.d,this.d*x.d);}
  neg(){return new F(-this.n,this.d);}
  sub(x){return this.add(F.of(x).neg());}
  mul(x){x=F.of(x);return new F(this.n*x.n,this.d*x.d);}
  div(x){x=F.of(x);return new F(this.n*x.d,this.d*x.n);}
  abs(){return new F(this.n<0n?-this.n:this.n,this.d);}
  eq(x){x=F.of(x);return this.n===x.n&&this.d===x.d;}
  le(x){x=F.of(x);return this.n*x.d<=x.n*this.d;}
  zero(){return this.n===0n;}
  pow(k){if(!Number.isSafeInteger(k)||k<0)throw new RangeError('Nonnegative integer power required');return new F(this.n**BigInt(k),this.d**BigInt(k));}
  toString(){return this.d===1n?`${this.n}`:`${this.n}/${this.d}`;}
  toJSON(){return this.toString();}
}
class Cut {
  constructor(rad=0,turn=0){this.rad=F.of(rad);this.turn=F.of(turn);Object.freeze(this);}
  static of(x){return x instanceof Cut?x:Array.isArray(x)&&x.length===2?new Cut(x[0],x[1]):new Cut(x);}
  add(x){x=Cut.of(x);return new Cut(this.rad.add(x.rad),this.turn.add(x.turn));}
  neg(){return new Cut(this.rad.neg(),this.turn.neg());}
  sub(x){return this.add(Cut.of(x).neg());}
  mul(x){x=Cut.of(x);return new Cut(this.rad.mul(x.rad).sub(this.turn.mul(x.turn)),this.rad.mul(x.turn).add(this.turn.mul(x.rad)));}
  dagger(){return new Cut(this.rad,this.turn.neg());}
  norm2(){return this.rad.mul(this.rad).add(this.turn.mul(this.turn));}
  gauge(){return this.rad.abs().add(this.turn.abs());}
  inv(){const n=this.norm2();if(n.zero())throw new RangeError('Zero cut scalar has no inverse');return new Cut(this.rad.div(n),this.turn.neg().div(n));}
  div(x){return this.mul(Cut.of(x).inv());}
  zero(){return this.rad.zero()&&this.turn.zero();}
  eq(x){x=Cut.of(x);return this.rad.eq(x.rad)&&this.turn.eq(x.turn);}
  toJSON(){return [this.rad.toString(),this.turn.toString()];}
}
const ZERO=new Cut(), ONE=new Cut(1), IOTA=new Cut(0,1);
function dimension(n){if(!Number.isSafeInteger(n)||n<1)throw new RangeError('Positive integer dimension required');}
function matrix(a){if(!Array.isArray(a)||!a.length||!Array.isArray(a[0])||!a[0].length)throw new TypeError('Nonempty rectangular matrix required');const n=a[0].length;if(a.some(r=>!Array.isArray(r)||r.length!==n))throw new TypeError('Ragged matrix');return a.map(r=>r.map(Cut.of));}
function zeros(m,n=m){dimension(m);dimension(n);return Array.from({length:m},()=>Array(n).fill(ZERO));}
function identity(n){const a=zeros(n);for(let i=0;i<n;i++)a[i][i]=ONE;return a;}
function shape(a){return [a.length,a[0].length];}
function sameShape(a,b){return a.length===b.length&&a[0].length===b[0].length;}
function add(a,b){a=matrix(a);b=matrix(b);if(!sameShape(a,b))throw new RangeError('Addition shape mismatch');return a.map((r,i)=>r.map((x,j)=>x.add(b[i][j])));}
function scale(a,s){s=Cut.of(s);return matrix(a).map(r=>r.map(x=>x.mul(s)));}
function sub(a,b){return add(a,scale(b,-1));}
function dagger(a){a=matrix(a);return a[0].map((_,j)=>a.map(r=>r[j].dagger()));}
function mul(a,b){a=matrix(a);b=matrix(b);if(a[0].length!==b.length)throw new RangeError('Composition shape mismatch');return a.map(r=>b[0].map((_,j)=>r.reduce((z,x,k)=>z.add(x.mul(b[k][j])),ZERO)));}
function equal(a,b){a=matrix(a);b=matrix(b);return sameShape(a,b)&&a.every((r,i)=>r.every((x,j)=>x.eq(b[i][j])));}
function mass(a){return matrix(a).flat().reduce((s,x)=>s.add(x.gauge()),new F());}
function energy(a){return matrix(a).flat().reduce((s,x)=>s.add(x.norm2()),new F());}
function isZero(a){return matrix(a).flat().every(x=>x.zero());}
function power(a,n){a=matrix(a);dimension(a.length);if(a.length!==a[0].length||!Number.isSafeInteger(n)||n<0)throw new RangeError('Square matrix and nonnegative power required');let r=identity(a.length);for(let k=0;k<n;k++)r=mul(r,a);return r;}
function kron(a,b){a=matrix(a);b=matrix(b);return a.flatMap(ar=>b.map(br=>ar.flatMap(x=>br.map(y=>x.mul(y)))));}
function rref(a){a=matrix(a);let p=0;const pivots=[];for(let j=0;j<a[0].length&&p<a.length;j++){const k=a.findIndex((r,i)=>i>=p&&!r[j].zero());if(k<0)continue;[a[p],a[k]]=[a[k],a[p]];const d=a[p][j];a[p]=a[p].map(x=>x.div(d));for(let i=0;i<a.length;i++)if(i!==p){const q=a[i][j];if(!q.zero())a[i]=a[i].map((x,t)=>x.sub(q.mul(a[p][t])));}pivots.push(j);p++;}return {matrix:a,pivots,rank:p,basis:a.slice(0,p)};}
function rank(a){return rref(a).rank;}
function inverse(a){a=matrix(a);const n=a.length;if(n!==a[0].length)throw new RangeError('Square inverse required');const id=identity(n);const r=rref(a.map((row,i)=>row.concat(id[i])));if(!equal(r.matrix.map(row=>row.slice(0,n)),id))throw new RangeError('Singular matrix');return r.matrix.map(row=>row.slice(n));}
function cayley(g,t){g=matrix(g);const id=identity(g.length);return mul(inverse(sub(id,scale(g,t))),add(id,scale(g,t)));}
function commutator(a,b){return sub(mul(a,b),mul(b,a));}
function flatten(a){return matrix(a).flat();}
function unflatten(v,n){if(v.length!==n*n)throw new RangeError('Wrong vectorized dimension');return Array.from({length:n},(_,i)=>v.slice(i*n,(i+1)*n));}
function validateGenerators(g,n){return g.map(a=>{a=matrix(a);if(a.length!==n||a[0].length!==n)throw new RangeError('Generator carrier mismatch');return a;});}
function rowClosure(seed,generators,{includeDagger=false}={}){
  seed=matrix(seed);const n=seed[0].length;let g=validateGenerators(generators,n);if(includeDagger)g=g.concat(g.map(dagger));
  let b=rref(seed).basis;const ranks=[b.length];if(!b.length)return {observer:[Array(n).fill(ZERO)],rank:0,ranks,extra:0};
  const initial=b.length;
  while(true){const c=rref(b.concat(...g.map(a=>mul(b,a)))).basis;if(c.length===b.length)break;b=c;ranks.push(b.length);if(b.length>n)throw new Error('Rank saturation failed');}
  return {observer:b,rank:b.length,ranks,extra:b.length-initial};
}
function operatorAlgebra(generators,{includeDagger=true}={}){
  if(!generators.length)throw new RangeError('At least one generator is required to specify dimension');const n=matrix(generators[0]).length;let g=validateGenerators(generators,n);if(includeDagger)g=g.concat(g.map(dagger));
  let b=[flatten(identity(n))];const ranks=[1];while(true){const added=b.flatMap(v=>g.map(a=>flatten(mul(unflatten(v,n),a))));const c=rref(b.concat(added)).basis;if(c.length===b.length)break;b=c;ranks.push(b.length);if(b.length>n*n)throw new Error('Algebra saturation failed');}
  return {basis:b.map(v=>unflatten(v,n)),dimension:b.length,ranks};
}
function descendedAction(g,c){g=matrix(g);c=matrix(c);const r=rref(c);if(r.rank!==c.length)throw new RangeError('Readout must have independent rows');const n=c[0].length,m=c.length;const pivotBlock=c.map(row=>r.pivots.map(j=>row[j]));const inv=inverse(pivotBlock);const section=zeros(n,m);r.pivots.forEach((j,i)=>section[j]=inv[i]);const t=mul(mul(c,g),section);return equal(mul(t,c),mul(c,g))?t:null;}
function regularAction(g,side){g=matrix(g);const n=g.length;if(n!==g[0].length||!['left','right'].includes(side))throw new TypeError('Square generator and left/right action required');const cols=[];for(let a=0;a<n;a++)for(let b=0;b<n;b++){const e=zeros(n);e[a][b]=ONE;cols.push(flatten(side==='left'?mul(g,e):mul(e,g)));}return Array.from({length:n*n},(_,i)=>cols.map(c=>c[i]));}
function idealReadoutClosureFullMatrix(seed,generators){
  const n=matrix(generators[0]).length;if(operatorAlgebra(generators).dimension!==n*n)throw new RangeError('This executable ideal adapter requires a generated full matrix algebra');
  seed=matrix(seed);if(seed[0].length!==n*n)throw new RangeError('Algebra readout acts on n^2 coordinates');let b=rref(seed).basis;const start=b.length,ranks=[start];if(!start)return {rank:0,extra:0,ranks};
  const g=generators.map(matrix).flatMap(a=>[a,dagger(a)]), actions=g.flatMap(a=>[regularAction(a,'left'),regularAction(a,'right')]);
  while(true){const star=b.map(row=>Array.from({length:n*n},(_,k)=>row[(k%n)*n+Math.floor(k/n)].dagger()));const c=rref(b.concat(star,...actions.map(a=>mul(b,a)))).basis;if(c.length===b.length)break;b=c;ranks.push(b.length);}
  return {observer:b,rank:b.length,extra:b.length-start,ranks};
}
function resolvent(a,z){a=matrix(a);return inverse(sub(scale(identity(a.length),z),a));}
function neumann(a,z,n){a=matrix(a);if(!Number.isSafeInteger(n)||n<0)throw new RangeError('Nonnegative truncation order required');z=Cut.of(z);const zi=z.inv(),q=mass(a).mul(zi.gauge());if(!q.le(1)||q.eq(1))throw new RangeError('Mass-Neumann bound requires q<1');const k=scale(a,zi);let sum=identity(a.length),term=identity(a.length);for(let j=1;j<=n;j++){term=mul(term,k);sum=add(sum,term);}return {approximation:scale(sum,zi),q,tail:zi.gauge().mul(q.pow(n+1)).div(new F(1).sub(q))};}
module.exports={F,Cut,ZERO,ONE,IOTA,matrix,zeros,identity,add,sub,scale,mul,dagger,equal,mass,energy,isZero,power,kron,rref,rank,inverse,cayley,commutator,flatten,unflatten,rowClosure,operatorAlgebra,descendedAction,idealReadoutClosureFullMatrix,resolvent,neumann};
