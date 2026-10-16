'use strict';
/** TET-1 publication adapter: exact finite path/transport coefficients over the
 * frozen native Cut matrix arithmetic and source TLM quotient calculus.
 * Ordinary jets cannot be lifted. ODE integration supplies one new coefficient
 * from a stated initial value; differentiation consumes one retained order.
 */
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const L=require('../native-tvsp-lambda-matrices/tvsp_adapter.cjs');
const {F,Cut,ZERO,ONE}=N;
const contexts=new WeakMap(),records=new WeakMap(),paths=new WeakMap();
const nonnegative=(x,name)=>{if(!Number.isSafeInteger(x)||x<0)throw new RangeError(name+' must be a nonnegative safe integer');};
function contextData(c){const d=contexts.get(c);if(!d)throw new TypeError('Native path context required');return d;}
function data(x){const d=records.get(x);if(!d)throw new TypeError('Native path coefficient record required');return d;}
function contextAt(chart,order){nonnegative(order,'Order');if(chart.cache.has(order))return chart.cache.get(order);const c=Object.freeze({kind:'TET_PATH_CONTEXT',order,label:chart.label,path:chart.path});contexts.set(c,{chart});chart.cache.set(order,c);return c;}
function createContext(order,label='declared path coefficients'){return contextAt({label,path:null,cache:new Map()},order);}
function lowerContext(c,order){const d=contextData(c);nonnegative(order,'Projected order');if(order>c.order)throw new RangeError('Cannot invent higher source layers');return contextAt(d.chart,order);}
function same(a,b){data(a);data(b);if(a.context!==b.context)throw new TypeError('Path context/depth mismatch; project explicitly');}
function freezeMatrix(a){return Object.freeze(N.matrix(a).map(r=>Object.freeze(r)));}
function series(c,coefficients){contextData(c);if(!Array.isArray(coefficients)||coefficients.length!==c.order+1)throw new TypeError('Supply exactly order+1 declared coefficients');const z=coefficients.map(freezeMatrix),height=z[0].length,width=z[0][0].length;if(z.some(a=>a.length!==height||a[0].length!==width))throw new TypeError('Coefficient matrix shape mismatch');const x=Object.freeze({kind:'TET_MATRIX_PATH_JET',context:c,height,width});records.set(x,Object.freeze(z));return x;}
function coefficients(x){return data(x);}
function coefficient(x,k){const z=data(x);nonnegative(k,'Coefficient order');if(k>x.context.order)throw new RangeError('Coefficient lies beyond retained order');return z[k];}
function constant(c,a){a=N.matrix(a);return series(c,Array.from({length:c.order+1},(_,k)=>k?N.zeros(a.length,a[0].length):a));}
function identity(c,n){return constant(c,N.identity(n));}
function project(x,order){const c=lowerContext(x.context,order);return c===x.context?x:series(c,data(x).slice(0,order+1));}
function add(a,b){same(a,b);return series(a.context,data(a).map((x,k)=>N.add(x,data(b)[k])));}
function scale(a,c){return series(a.context,data(a).map(x=>N.scale(x,Cut.of(c))));}
function sub(a,b){return add(a,scale(b,-1));}
function mul(a,b){same(a,b);if(a.width!==b.height)throw new RangeError('Path matrix composition shape mismatch');return series(a.context,Array.from({length:a.context.order+1},(_,k)=>{let x=N.zeros(a.height,b.width);for(let j=0;j<=k;j++)x=N.add(x,N.mul(data(a)[j],data(b)[k-j]));return x;}));}
function dagger(x){return series(x.context,data(x).map(N.dagger));}
function equal(a,b){same(a,b);return a.height===b.height&&a.width===b.width&&data(a).every((x,k)=>N.equal(x,data(b)[k]));}
function derivative(x){data(x);if(x.context.order===0)throw new RangeError('Depth-zero data does not determine a derivative');const c=lowerContext(x.context,x.context.order-1);return series(c,data(x).slice(1).map((a,k)=>N.scale(a,k+1)));}
function inverse(x){data(x);if(x.height!==x.width)throw new RangeError('Square inverse required');const z=data(x),out=[N.inverse(z[0])];for(let n=1;n<=x.context.order;n++){let sum=N.zeros(x.height);for(let k=1;k<=n;k++)sum=N.add(sum,N.mul(z[k],out[n-k]));out.push(N.scale(N.mul(out[0],sum),-1));}return series(x.context,out);}
function power(x,k){data(x);nonnegative(k,'Power');if(x.height!==x.width)throw new RangeError('Square power required');let y=identity(x.context,x.height);for(let n=0;n<k;n++)y=mul(y,x);return y;}
function integrateZero(x){data(x);const c=contextAt(contextData(x.context).chart,x.context.order+1);return series(c,[N.zeros(x.height,x.width),...data(x).map((a,k)=>N.scale(a,new F(1,k+1)))]);}
function transport(A,initial=N.identity(A.height)){
  data(A);if(A.height!==A.width)throw new RangeError('Connection must be square');initial=N.matrix(initial);if(initial.length!==A.height)throw new RangeError('Initial section carrier mismatch');
  const out=[initial],z=data(A);for(let n=0;n<=A.context.order;n++){let sum=N.zeros(A.height,initial[0].length);for(let j=0;j<=n;j++)sum=N.add(sum,N.mul(z[j],out[n-j]));out.push(N.scale(sum,new F(-1,n+1)));}
  return series(contextAt(contextData(A.context).chart,A.context.order+1),out);
}
function ordinaryExponential(B){data(B);if(B.height!==B.width||!N.isZero(data(B)[0]))throw new RangeError('Finite ordinary exponential requires square zero-constant exponent');let out=identity(B.context,B.height),term=out;for(let k=1;k<=B.context.order;k++){term=scale(mul(term,B),new F(1,k));out=add(out,term);}return out;}
function observer(U,C0){data(U);const C=constant(U.context,C0);return mul(C,inverse(U));}
function pathSpec(raw={}){const a=Cut.of(raw.a===undefined?0:raw.a),b=Cut.of(raw.b===undefined?1:raw.b),p=raw.p===undefined?1:raw.p,q=raw.q===undefined?1:raw.q;nonnegative(p,'s path degree');nonnegative(q,'v path degree');if(!p||!q||!a.turn.zero()||!b.turn.zero())throw new TypeError('Real zero-origin paths s=a t^p,v=b t^q require positive integer p,q');return Object.freeze({a,b,p,q});}
function pathContext(source,raw){const path=pathSpec(raw),key=JSON.stringify([path.a,path.b,path.p,path.q]);let chart;
  for(let d=source.order;d>=0;d--){const lower=L.lowerContext(source,d),found=paths.get(lower)?.get(key);if(found){chart=contextData(found).chart;break;}}
  if(!chart)chart={label:source.label+' restricted path',path,cache:new Map()};
  for(let d=0;d<=source.order;d++){const lower=L.lowerContext(source,d);if(!paths.has(lower))paths.set(lower,new Map());paths.get(lower).set(key,contextAt(chart,d));}
  return contextAt(chart,source.order);
}
function cutPower(x,n){let z=ONE;for(let k=0;k<n;k++)z=z.mul(x);return z;}
function restrictScalar(x,raw={}){const c=pathContext(x.context,raw),p=c.path,z=Array.from({length:c.order+1},()=>[[ZERO]]);for(const [a,b,value] of L.entries(x)){const k=a*p.p+b*p.q;if(k<=c.order)z[k][0][0]=z[k][0][0].add(value.mul(cutPower(p.a,a)).mul(cutPower(p.b,b)));}return series(c,z);}
function restrictMatrix(H,raw={}){const c=pathContext(H.context,raw),parts=H.rows.map(r=>r.map(x=>restrictScalar(x,raw)));return series(c,Array.from({length:c.order+1},(_,k)=>parts.map(r=>r.map(x=>coefficient(x,k)[0][0]))));}
function connectionFromMetric(H){data(H);if(H.context.order<1)throw new RangeError('Metric depth at least one required');if(!equal(H,dagger(H)))throw new TypeError('Hermitian metric path required');const d=H.context.order-1;return scale(mul(project(inverse(H),d),derivative(H)),new F(1,2));}
function sourcePath(fixture,raw={}){const H=restrictMatrix(fixture.H,raw),A=connectionFromMetric(H),U=transport(A);return Object.freeze({H,A,U,r:restrictScalar(fixture.r,raw),path:H.context.path});}
function affineBinomialTransport(H){data(H);if(H.height!==H.width||data(H).slice(2).some(a=>!N.isZero(a)))throw new TypeError('An affine retained metric path is required');const B=N.mul(N.inverse(coefficient(H,0)),H.context.order?coefficient(H,1):N.zeros(H.height));let binomial=new F(1),term=N.identity(H.height);const out=[term];for(let k=1;k<=H.context.order;k++){binomial=binomial.mul(new F(1-2*k,2*k));term=N.mul(term,B);out.push(N.scale(term,binomial));}return series(H.context,out);}
function sourceGeometry(fixture){const z=fixture,d=z.context.order;if(d<2)throw new RangeError('Geometry needs source depth at least two');const Hinv=L.matrix([[L.div(z.c,z.Delta),L.scale(L.div(z.b,z.Delta),-1)],[L.scale(L.div(z.b,z.Delta),-1),L.div(z.a,z.Delta)]]);const As=L.matrixScale(L.matrixMul(L.matrixProject(Hinv,d-1),L.matrixDerivative(z.H,'s')),'1/2'),Av=L.matrixScale(L.matrixMul(L.matrixProject(Hinv,d-1),L.matrixDerivative(z.H,'v')),'1/2');const DsAv=L.matrixDerivative(Av,'s'),DvAs=L.matrixDerivative(As,'v'),comm=L.matrixProject(L.matrixSub(L.matrixMul(As,Av),L.matrixMul(Av,As)),d-2);const Fsv=L.matrixAdd(L.matrixSub(DsAv,DvAs),comm);return Object.freeze({As,Av,Fsv,DsAv,DvAs});}
function smallLoop(fixture){const g=sourceGeometry(fixture),S=L.matrixOuter(g.As),V=L.matrixOuter(g.Av),Ss=L.matrixOuter(L.matrixDerivative(g.As,'s')),Sv=L.matrixOuter(g.DvAs),Vs=L.matrixOuter(g.DsAv),Vv=L.matrixOuter(L.matrixDerivative(g.Av,'v'));const c=createContext(2,'oriented +s,+v,-s,-v source loop'),I=N.identity(2),half=a=>N.scale(a,'1/2');const E1=series(c,[I,N.scale(S,-1),half(N.sub(N.mul(S,S),Ss))]);const E2=series(c,[I,N.scale(V,-1),N.sub(half(N.sub(N.mul(V,V),Vv)),Vs)]);const E3=series(c,[I,S,N.add(half(N.add(N.mul(S,S),Ss)),Sv)]);const E4=series(c,[I,V,half(N.add(N.mul(V,V),Vv))]);return Object.freeze({edges:Object.freeze([E1,E2,E3,E4]),holonomy:mul(E4,mul(E3,mul(E2,E1))),curvature:L.matrixOuter(g.Fsv),derivatives:Object.freeze({S,V,Ss,Sv,Vs,Vv})});}
function serialize(x){return {order:x.context.order,height:x.height,width:x.width,coefficients:data(x).map(a=>a.map(r=>r.map(c=>c.toJSON())))};}
module.exports=Object.freeze({createContext,lowerContext,series,coefficients,coefficient,constant,identity,project,add,sub,scale,mul,dagger,equal,derivative,inverse,power,integrateZero,transport,ordinaryExponential,observer,restrictScalar,restrictMatrix,connectionFromMetric,sourcePath,affineBinomialTransport,sourceGeometry,smallLoop,serialize});
