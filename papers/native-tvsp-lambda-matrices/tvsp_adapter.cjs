'use strict';
/** Publication adapter, not a replacement arithmetic engine.
 * All polynomial sums, products, reductions, derivatives and adjoints use the
 * frozen Paninian engine. Rational/complex coefficients are native F/Cut.
 * A depth-N record belongs to C[s,v]/(s,v)^(N+1), with self-adjoint formal s,v.
 * Derivatives map depth N to N-1; no derivation at unchanged depth is asserted.
 * The faithful regular multiplication representation below is NOT a *-rep.
 */
const native = require('../native-common-cone-limit/native_engine/native_operator.cjs');
const panini = require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const {F, Cut, ZERO, ONE} = native;
const {Expr, Presentation, polynomialDerivative, star} = panini;
const contextData = new WeakMap(), valueData = new WeakMap(), matrixData = new WeakMap();

function nonnegative(n, name) {
  if (!Number.isSafeInteger(n) || n < 0) throw new RangeError(name + ' must be a nonnegative safe integer');
  return n;
}
function ctxData(context) {
  const d = contextData.get(context);
  if (!d) throw new TypeError('A context created by this adapter is required');
  return d;
}
function valData(x) {
  const d = valueData.get(x);
  if (!d) throw new TypeError('A scalar record created by this adapter is required');
  return d;
}
function matData(x) {
  const d = matrixData.get(x);
  if (!d) throw new TypeError('A matrix record created by this adapter is required');
  return d;
}
function sameContext(a, b) {
  if (a.context !== b.context) throw new TypeError('Different chart/context or depth: project explicitly within one chart');
}
function word(a, b) {
  nonnegative(a, 's exponent'); nonnegative(b, 'v exponent');
  return Array(a).fill('s').concat(Array(b).fill('v'));
}
function exponents(w) {
  return [w.filter(t => t === 's').length, w.filter(t => t === 'v').length];
}
function makeContext(N, chart) {
  nonnegative(N, 'Order');
  if (chart.contexts.has(N)) return chart.contexts.get(N);
  const rules = [{id:'central_v_s', lhs:['v','s'], rhs:[[['s','v'],ONE]],
    source:'TLM formal commuting local coordinates; v s = s v'}];
  for (let a = 0; a <= N + 1; a++) rules.push({
    id:'truncate_' + a + '_' + (N + 1 - a), lhs:word(a, N + 1 - a), rhs:[],
    source:'TLM finite total-degree quotient (s,v)^(N+1) = 0'
  });
  const system = new Presentation({tokens:['s','v'], rules});
  const context = Object.freeze({kind:'TVSP_LAMBDA_CONTEXT', order:N, label:chart.label, source:chart.source});
  contextData.set(context, {system, chart}); chart.contexts.set(N, context);
  return context;
}
function createContext(N, label='TVSP local chart') {
  if (typeof label !== 'string') throw new TypeError('Context label must be a string');
  return makeContext(N, {label, source:null, contexts:new Map()});
}
function lowerContext(context, N) {
  const d = ctxData(context); nonnegative(N, 'Projected order');
  if (N > context.order) throw new RangeError('Projection cannot invent higher layers');
  return makeContext(N, d.chart);
}
function wrap(context, expr) {
  const {system} = ctxData(context);
  if (!(expr instanceof Expr) || expr.system !== system) throw new TypeError('Native expression belongs to another presentation');
  const normal = system.reduce(expr).normal;
  const x = Object.freeze({kind:'TVSP_LAMBDA_SCALAR', context});
  valueData.set(x, {expr:normal});
  return x;
}
function scalar(context, c=0) {
  return wrap(context, ctxData(context).system.one().scale(Cut.of(c)));
}
function monomial(context, a, b, c=1) {
  return wrap(context, ctxData(context).system.word(word(a,b), Cut.of(c)));
}
function coordinate(context, name) {
  if (name !== 's' && name !== 'v') throw new TypeError('Local coordinate must be s or v');
  return monomial(context, name === 's' ? 1 : 0, name === 'v' ? 1 : 0);
}
function add(a,b) { const x=valData(a),y=valData(b);sameContext(a,b);return wrap(a.context,x.expr.plus(y.expr)); }
function sub(a,b) { const x=valData(a),y=valData(b);sameContext(a,b);return wrap(a.context,x.expr.minus(y.expr)); }
function mul(a,b) { const x=valData(a),y=valData(b);sameContext(a,b);return wrap(a.context,x.expr.times(y.expr)); }
function scale(a,c) { return wrap(a.context,valData(a).expr.scale(Cut.of(c))); }
function coefficient(x,a,b) {
  const {expr}=valData(x);nonnegative(a,'s exponent');nonnegative(b,'v exponent');
  if(a+b>x.context.order)throw new RangeError('Coefficient layer exceeds the retained order');
  return expr.terms.get(JSON.stringify(word(a,b))) || ZERO;
}
function outer(x) { return coefficient(x,0,0); }
function equal(a,b) { const x=valData(a),y=valData(b);sameContext(a,b);return x.expr.equals(y.expr); }
function power(x,n) {
  valData(x);nonnegative(n,'Power');let r=scalar(x.context,1),p=x,k=n;
  while(k){if(k%2)r=mul(r,p);k=Math.floor(k/2);if(k)p=mul(p,p);}
  return r;
}
function inverse(x) {
  valData(x);const a=outer(x);if(a.zero())throw new RangeError('Only a nonzero constant term is invertible in the finite local algebra');
  const ai=a.inv(),q=scale(sub(x,scalar(x.context,a)),ai.neg());
  let sum=scalar(x.context,1),term=sum;
  for(let k=1;k<=x.context.order;k++){term=mul(term,q);sum=add(sum,term);}
  return scale(sum,ai);
}
function div(a,b) { valData(a);valData(b);sameContext(a,b);return mul(a,inverse(b)); }
function dagger(x) {
  const {expr}=valData(x),{system}=ctxData(x.context);
  return wrap(x.context,star(expr,{s:system.word(['s']),v:system.word(['v'])}));
}
function project(x,N) {
  const {expr}=valData(x),context=lowerContext(x.context,N),{system}=ctxData(context);
  if(context===x.context)return x;
  return wrap(context,new Expr(system,'*','*',expr.entries()));
}
function derivative(x,name) {
  const {expr}=valData(x);
  if(name!=='s'&&name!=='v')throw new TypeError('Derivative coordinate must be s or v');
  if(x.context.order===0)throw new RangeError('A depth-zero record has no determined derivative');
  const source=ctxData(x.context).system;
  const raw=polynomialDerivative(expr,{s:name==='s'?source.one():source.zero(),v:name==='v'?source.one():source.zero()});
  const context=lowerContext(x.context,x.context.order-1),{system}=ctxData(context);
  return wrap(context,new Expr(system,'*','*',raw.entries()));
}
function entries(x) {
  return Object.freeze(valData(x).expr.entries().map(([w,c])=>Object.freeze([...exponents(w),c])));
}
function basisIndices(context) {
  ctxData(context);const out=[];
  for(let d=0;d<=context.order;d++)for(let a=d;a>=0;a--)out.push(Object.freeze([a,d-a]));
  return Object.freeze(out);
}
function serialize(x) {
  // Coefficient-only snapshot, not a serialization of the full source context.
  // A complete contextual record also retains its source, H/Lambda and frames.
  valData(x);
  return Object.freeze({kind:x.kind,order:x.context.order,label:x.context.label,
    coefficients:Object.freeze(entries(x).map(([a,b,c])=>Object.freeze({s:a,v:b,value:Object.freeze(c.toJSON())})))});
}
function presentationSpec(context) {
  // Detached JSON snapshot: callers cannot mutate the engine or its rules.
  return JSON.parse(JSON.stringify(ctxData(context).system.spec()));
}
function auditPresentation(context, options={}) { return ctxData(context).system.audit(options); }

function matrix(rows) {
  if(!Array.isArray(rows)||!rows.length||!Array.isArray(rows[0])||!rows[0].length)throw new TypeError('A nonempty rectangular matrix of scalar records is required');
  const n=rows[0].length,context=rows[0][0]&&rows[0][0].context;
  if(rows.some(row=>!Array.isArray(row)||row.length!==n))throw new TypeError('Ragged scalar-record matrix');
  const safe=rows.map(row=>Object.freeze(row.map(x=>{valData(x);if(x.context!==context)throw new TypeError('Matrix entries must share one context and depth');return x;})));
  const x=Object.freeze({kind:'TVSP_LAMBDA_MATRIX',context,rows:Object.freeze(safe),height:rows.length,width:n});
  matrixData.set(x,{rows:x.rows});return x;
}
function matrixAdd(a,b) {
  matData(a);matData(b);sameContext(a,b);
  if(a.height!==b.height||a.width!==b.width)throw new RangeError('Matrix addition shape mismatch');
  return matrix(a.rows.map((r,i)=>r.map((x,j)=>add(x,b.rows[i][j]))));
}
function matrixScale(a,c) { matData(a);return matrix(a.rows.map(r=>r.map(x=>scale(x,c)))); }
function matrixSub(a,b) { return matrixAdd(a,matrixScale(b,-1)); }
function matrixMul(a,b) {
  matData(a);matData(b);sameContext(a,b);if(a.width!==b.height)throw new RangeError('Matrix multiplication shape mismatch');
  return matrix(a.rows.map(r=>Array.from({length:b.width},(_,j)=>r.reduce((sum,x,k)=>add(sum,mul(x,b.rows[k][j])),scalar(a.context)))));
}
function matrixDagger(a) { matData(a);return matrix(Array.from({length:a.width},(_,j)=>a.rows.map(r=>dagger(r[j])))); }
function identity(context,n) {
  ctxData(context);if(!Number.isSafeInteger(n)||n<1)throw new RangeError('Positive matrix dimension required');
  return matrix(Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>scalar(context,i===j?1:0))));
}
function matrixProject(a,N) { matData(a);return matrix(a.rows.map(r=>r.map(x=>project(x,N)))); }
function matrixDerivative(a,name) { matData(a);return matrix(a.rows.map(r=>r.map(x=>derivative(x,name)))); }
function matrixOuter(a) { matData(a);return Object.freeze(a.rows.map(r=>Object.freeze(r.map(outer)))); }
function matrixEqual(a,b) {
  matData(a);matData(b);sameContext(a,b);
  return a.height===b.height&&a.width===b.width&&a.rows.every((r,i)=>r.every((x,j)=>equal(x,b.rows[i][j])));
}
function matrixPower(a,n) {
  matData(a);nonnegative(n,'Matrix power');if(a.height!==a.width)throw new RangeError('Square matrix required');
  let out=identity(a.context,a.height);for(let k=0;k<n;k++)out=matrixMul(out,a);return out;
}
function regularMatrix(a) {
  matData(a);const basis=basisIndices(a.context),d=basis.length;
  const out=native.zeros(d*a.height,d*a.width);
  // Degree-major basis: every diagonal block is the outer matrix. The first
  // monomial column recovers all coefficients, proving algebraic faithfulness.
  for(let k=0;k<d;k++){
    const e=monomial(a.context,...basis[k]);
    for(let i=0;i<a.height;i++)for(let j=0;j<a.width;j++){
      const p=mul(a.rows[i][j],e);
      for(let l=0;l<d;l++)out[l*a.height+i][k*a.width+j]=coefficient(p,...basis[l]);
    }
  }
  return Object.freeze(out.map(r=>Object.freeze(r)));
}
function serializeMatrix(a) { matData(a);return Object.freeze({kind:a.kind,rows:Object.freeze(a.rows.map(r=>Object.freeze(r.map(serialize))))}); }

function buildFixture(mu=0,N=3) {
  nonnegative(N,'Fixture order');const muCut=Cut.of(mu);
  if(!muCut.turn.zero())throw new TypeError('Thermodynamic fixture parameter mu must be real rational');
  const source=Object.freeze({kind:'TLM_SHIFTED_THERMODYNAMIC_SOURCE',mu:muCut,
    coordinates:'S=1+s, V=1+v',potential:'u=10s-10v+s^2+sv+(5/2)v^2+(mu/2)s^2v',
    convention:'T=u_s; P=-u_v; H=[[u_ss,u_sv],[u_sv,u_vv]]'});
  const chart={label:'TVSP source mu='+muCut.rad.toString(),source,contexts:new Map()};
  const precursorContext=makeContext(N+2,chart);
  const uTerms=[[1,0,10],[0,1,-10],[2,0,1],[1,1,1],[0,2,new F(5,2)],[2,1,muCut.mul(new F(1,2))]];
  let u=scalar(precursorContext);for(const [a,b,c] of uTerms)u=add(u,monomial(precursorContext,a,b,c));
  const Tfull=derivative(u,'s'),Pfull=scale(derivative(u,'v'),-1);
  const a=derivative(Tfull,'s'),b=derivative(Tfull,'v'),c=scale(derivative(Pfull,'v'),-1);
  const context=a.context,S=add(scalar(context,1),coordinate(context,'s')),V=add(scalar(context,1),coordinate(context,'v'));
  const T=project(Tfull,N),P=project(Pfull,N),Delta=sub(mul(a,c),mul(b,b));
  const responses=Object.freeze({
    lambdaP:scale(div(mul(S,Delta),c),-1),
    lambdaV:scale(mul(S,a),-1),
    lambdaS:scale(mul(V,c),-1),
    zT:scale(div(mul(V,Delta),a),-1),
    zS:div(P,c),
    lambdaT:div(mul(P,a),Delta)
  });
  const H=matrix([[a,b],[b,c]]),r=scale(responses.lambdaS,-1);
  return Object.freeze({context,ctx:context,source,S,V,T,P,a,b,c,Delta,H,responses,r,
    u:project(u,N),precursor:Object.freeze({context:precursorContext,u,T:Tfull,P:Pfull})});
}

/** Actual NT tensor layers with a flat frame on the six-response fibre E.
 * Returns a frozen array. Layer k has derivativeOrder=k, depth=N-k and
 * components[responseName][word], where a word has k characters in {0,1}.
 * Indices are prepended: component "ijk" means (nabla_i nabla^2 f)_(jk).
 * Gamma_i = (1/2) H^{-1} d_i H is built at depth N-1; every subsequent
 * product is explicitly projected to the output depth. These are formal
 * tensor jets, not claims about continuum limits or physical quantization.
 * Source Hessian integrability and patch positivity are caller hypotheses;
 * the algebraic shape/symmetry checks below do not certify an arbitrary state.
 */
function responseTower(fixture,maxDerivativeOrder=fixture.context.order) {
  if(!fixture||typeof fixture!=='object')throw new TypeError('A TVSP fixture is required');
  const context=fixture.context;ctxData(context);nonnegative(maxDerivativeOrder,'Covariant derivative order');
  const N=context.order;
  if(maxDerivativeOrder>N)throw new RangeError('Covariant derivative order exceeds the retained depth');
  const names=['lambdaP','lambdaV','lambdaS','zT','zS','lambdaT'];
  const initial={};
  for(const name of names){
    const response=fixture.responses&&fixture.responses[name];valData(response);
    if(response.context!==context)throw new TypeError('Responses must share the fixture context');
    initial[name]=Object.freeze({'':response});
  }
  const layers=[Object.freeze({derivativeOrder:0,depth:N,components:Object.freeze(initial)})];
  if(maxDerivativeOrder===0)return Object.freeze(layers);
  const H=fixture.H;matData(H);
  if(H.context!==context||H.height!==2||H.width!==2)throw new TypeError('A two-dimensional Hessian in the fixture context is required');
  if(!equal(H.rows[0][1],H.rows[1][0]))throw new TypeError('The Hessian must be symmetric');
  const [[h00,h01],[h10,h11]]=H.rows;
  const determinant=sub(mul(h00,h11),mul(h01,h10)),detInverse=inverse(determinant);
  const Hinv=matrix([[mul(h11,detInverse),scale(mul(h01,detInverse),-1)],
    [scale(mul(h10,detInverse),-1),mul(h00,detInverse)]]);
  const inverseAtConnectionDepth=matrixProject(Hinv,N-1);
  const gamma=['s','v'].map(name=>matrixScale(
    matrixMul(inverseAtConnectionDepth,matrixDerivative(H,name)),new F(1,2)));
  for(let k=1;k<=maxDerivativeOrder;k++){
    const depth=N-k,previous=layers[k-1].components,components={};
    const connection=gamma.map(g=>matrixProject(g,depth));
    for(const name of names){
      const output={};
      for(const [w,old] of Object.entries(previous[name]))for(let i=0;i<2;i++){
        let value=derivative(old,i===0?'s':'v');
        for(let slot=0;slot<w.length;slot++)for(let r=0;r<2;r++){
          const replaced=w.slice(0,slot)+String(r)+w.slice(slot+1);
          const correction=mul(connection[i].rows[r][Number(w[slot])],project(previous[name][replaced],depth));
          value=sub(value,correction);
        }
        output[String(i)+w]=value;
      }
      components[name]=Object.freeze(output);
    }
    layers.push(Object.freeze({derivativeOrder:k,depth,components:Object.freeze(components)}));
  }
  return Object.freeze(layers);
}

module.exports=Object.freeze({F,Cut,createContext,lowerContext,scalar,monomial,coordinate,
  add,sub,mul,scale,inverse,div,power,dagger,project,derivative,coefficient,equal,outer,
  entries,basisIndices,serialize,presentationSpec,auditPresentation,
  matrix,matrixAdd,matrixSub,matrixMul,matrixScale,matrixDagger,identity,
  matrixProject,matrixDerivative,matrixOuter,matrixEqual,matrixPower,regularMatrix,serializeMatrix,buildFixture,responseTower});
