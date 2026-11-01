'use strict';
/** Typed, proof-carrying finite presentations over the native cut field.
 * Equality rewriting is not priority execution or erasure of physical history.
 * Only degree-lex decreasing monic relations are admitted. No Hilbert premise.
 */
const crypto = require('node:crypto');
const {Cut, ZERO, ONE, matrix, identity, zeros, add: mAdd, scale: mScale,
  mul: mMul, equal: mEqual, rank, flatten} = require('./native_operator.cjs');
function stable(x) {
  if (x instanceof Cut) return x.toJSON();
  if (Array.isArray(x)) return x.map(stable);
  if (x && typeof x === 'object') return Object.fromEntries(Object.keys(x).sort().map(k => [k, stable(x[k])]));
  return x;
}
const digest = x => crypto.createHash('sha256').update(JSON.stringify(stable(x))).digest('hex');
const key = w => JSON.stringify(w);
function bound(n, name, positive = false) {
  if (!Number.isSafeInteger(n) || n < (positive ? 1 : 0)) throw new RangeError(name);
}
function charge(x) {
  if (typeof x === 'bigint') return x;
  if (typeof x === 'number' && Number.isSafeInteger(x)) return BigInt(x);
  if (typeof x === 'string' && /^-?\d+$/.test(x)) return BigInt(x);
  throw new TypeError('Exact integer charge required');
}
class Budget extends Error { constructor(){super('Declared reduction budget exhausted');} }
class Expr {
  constructor(system, from, to, entries = []) {
    this.system = system; this.from = from; this.to = to;
    if (!system.objects.includes(from) || !system.objects.includes(to)) throw new TypeError('Unknown object');
    this.terms = new Map();
    for (const [w0, c0] of entries) {
      const w = [...w0], t = system.wordType(w, from), c = Cut.of(c0);
      if (t.from !== from || t.to !== to) throw new TypeError('Expression has mixed arrow types');
      const k = key(w), v = (this.terms.get(k) || ZERO).add(c);
      if (v.zero()) this.terms.delete(k); else this.terms.set(k, v);
    }
  }
  entries(){return [...this.terms].map(([w,c]) => [JSON.parse(w), c]).sort((a,b) => this.system.compare(a[0],b[0]));}
  plus(b){this.system.same(this,b); return new Expr(this.system,this.from,this.to,this.entries().concat(b.entries()));}
  scale(c){return new Expr(this.system,this.from,this.to,this.entries().map(([w,v]) => [w,v.mul(c)]));}
  minus(b){return this.plus(b.scale(-1));}
  times(b){
    if (this.system !== b.system || this.from !== b.to) throw new TypeError('Unmatched recovered middle object');
    return new Expr(this.system,b.from,this.to,this.entries().flatMap(([u,a]) => b.entries().map(([v,c]) => [u.concat(v),a.mul(c)])));
  }
  equals(b){return this.system === b.system && this.from === b.from && this.to === b.to && this.minus(b).terms.size === 0;}
  toJSON(){return {from:this.from,to:this.to,terms:this.entries().map(([w,c])=>[w,c.toJSON()])};}
}
class Presentation {
  constructor({objects=['*'], tokens, rules=[]}) {
    if (!Array.isArray(objects) || !objects.length || new Set(objects).size !== objects.length || objects.some(x => typeof x !== 'string' || !x)) throw new TypeError('Distinct named objects required');
    if (!Array.isArray(tokens) || !tokens.length) throw new TypeError('Ordered generator catalogue required');
    this.objects = [...objects]; this.tokens = Object.create(null); this.order = new Map(); this.rules = [];
    tokens.forEach((raw,i) => {
      const t = typeof raw === 'string' ? {name:raw,from:objects[0],to:objects[0],charge:'0'} : raw;
      if (typeof t.name !== 'string' || !t.name || this.order.has(t.name) || !objects.includes(t.from) || !objects.includes(t.to)) throw new TypeError('Invalid generator signature');
      this.tokens[t.name] = {...t, charge:charge(t.charge === undefined ? 0 : t.charge).toString()}; this.order.set(t.name,i);
    });
    rules.forEach(r => this.addRule(r));
  }
  compare(a,b){if(a.length!==b.length)return a.length-b.length;for(let i=0;i<a.length;i++){const d=this.order.get(a[i])-this.order.get(b[i]);if(d)return d;}return 0;}
  wordType(w, emptyAt=this.objects[0]) {
    if (!Array.isArray(w)) throw new TypeError('A word is a token array');
    if (!w.length) return {from:emptyAt,to:emptyAt};
    for(const x of w) if(!Object.hasOwn(this.tokens,x)) throw new TypeError('Unknown generator: '+x);
    for(let i=0;i<w.length-1;i++)if(this.tokens[w[i]].from!==this.tokens[w[i+1]].to)throw new TypeError('Noncomposable native word');
    return {from:this.tokens[w.at(-1)].from,to:this.tokens[w[0]].to};
  }
  wordCharge(w){return w.reduce((s,t)=>s+BigInt(this.tokens[t].charge),0n);}
  word(w, c=ONE, at=this.objects[0]){const t=this.wordType(w,at);return new Expr(this,t.from,t.to,[[w,c]]);}
  zero(from=this.objects[0],to=from){return new Expr(this,from,to);}
  one(at=this.objects[0]){return this.word([],ONE,at);}
  fromJSON(x){return new Expr(this,x.from,x.to,x.terms);}
  same(a,b){if(a.system!==this||b.system!==this||a.from!==b.from||a.to!==b.to)throw new TypeError('Typed sum mismatch');}
  spec(){return {objects:this.objects,tokens:[...this.order.keys()].map(k=>this.tokens[k]),rules:this.rules.map(r=>({id:r.id,lhs:r.lhs,rhs:r.rhs.entries().map(([w,c])=>[w,c.toJSON()]),source:r.source}))};}
  hash(){return digest(this.spec());}
  addRule(r){
    if (typeof r.id!=='string'||!r.id||this.rules.some(x=>x.id===r.id)||!Array.isArray(r.lhs)||!r.lhs.length) throw new TypeError('Named nonempty leading word required');
    if (!r.source || typeof r.source !== 'string') throw new TypeError('Every rule needs a declared source/contract');
    const lhs=[...r.lhs],t=this.wordType(lhs),rhs=new Expr(this,t.from,t.to,r.rhs||[]),q=this.wordCharge(lhs);
    for(const [w] of rhs.entries()) {
      if(this.compare(w,lhs)>=0) throw new TypeError('Rule is not strictly degree-lex decreasing');
      if(this.wordCharge(w)!==q)throw new TypeError('Rule would erase a declared integer residue');
    }
    this.rules.push({id:r.id,lhs,rhs,source:r.source}); return this;
  }
  occurrences(w){const out=[];this.rules.forEach((r,ri)=>{for(let p=0;p+r.lhs.length<=w.length;p++)if(r.lhs.every((t,j)=>w[p+j]===t))out.push({ri,p});});return out;}
  replaceWord(w,ri,p,at){
    const r=this.rules[ri]; if(!r || !Number.isSafeInteger(p)||p<0||!r.lhs.every((t,j)=>w[p+j]===t))throw new Error('Invalid claimed rewrite location');
    const t=this.wordType(w,at),a=w.slice(0,p),b=w.slice(p+r.lhs.length);
    return new Expr(this,t.from,t.to,r.rhs.entries().map(([v,c])=>[a.concat(v,b),c]));
  }
  reduce(expr,{strategy='first',maxSteps=100000,maxTerms=100000,witness=false}={}) {
    if(expr.system!==this||!['first','last'].includes(strategy))throw new TypeError('Presentation/strategy mismatch');
    bound(maxSteps,'Invalid step budget');bound(maxTerms,'Invalid term budget',true);
    let e=expr,steps=0;const trace=[];
    while(true){
      if(e.terms.size>maxTerms)throw new Budget();
      let move;
      for(const [w,c] of e.entries().reverse()){const options=this.occurrences(w);if(options.length){move={w,c,...options[strategy==='first'?0:options.length-1]};break;}}
      if(!move)break;if(steps++>=maxSteps)throw new Budget();
      const {w,c,ri,p}=move;
      if(witness)trace.push({word:w,coefficient:c.toJSON(),rule:this.rules[ri].id,position:p});
      e=e.minus(this.word(w,c,expr.from)).plus(this.replaceWord(w,ri,p,expr.from).scale(c));
    }
    return {normal:e,steps,certificate:witness?{protocol:'NATIVE_REWRITE_WITNESS_V1',presentation:this.hash(),input:expr.toJSON(),steps:trace,output:e.toJSON()}:null};
  }
  replay(cert){
    if(cert.protocol!=='NATIVE_REWRITE_WITNESS_V1'||cert.presentation!==this.hash())throw new Error('Wrong proof contract/source hash');
    let e=this.fromJSON(cert.input);
    for(const s of cert.steps){
      const w=s.word,c=Cut.of(s.coefficient),existing=e.terms.get(key(w)),ri=this.rules.findIndex(r=>r.id===s.rule);
      if(!existing||!existing.eq(c)||ri<0)throw new Error('Witness coefficient/rule mismatch');
      e=e.minus(this.word(w,c,e.from)).plus(this.replaceWord(w,ri,s.position,e.from).scale(c));
    }
    if(!e.equals(this.fromJSON(cert.output))||e.entries().some(([w])=>this.occurrences(w).length))throw new Error('Witness does not finish at claimed normal form');
    return true;
  }
  ambiguities(){
    const out=[],seen=new Set();
    this.rules.forEach((a,i)=>this.rules.forEach((b,j)=>{
      for(let shift=1-b.lhs.length;shift<a.lhs.length;shift++){
        const lo=Math.min(0,shift),hi=Math.max(a.lhs.length,shift+b.lhs.length),w=Array(hi-lo),pa=-lo,pb=shift-lo;
        a.lhs.forEach((t,k)=>w[pa+k]=t);let valid=true;
        b.lhs.forEach((t,k)=>{if(w[pb+k]!==undefined&&w[pb+k]!==t)valid=false;w[pb+k]=t;});
        if(!valid||(i===j&&pa===pb))continue;
        try{this.wordType(w);}catch{continue;}
        const positions=[[i,pa],[j,pb]].sort((x,y)=>x[0]-y[0]||x[1]-y[1]);
        const id=key([w,positions]);if(seen.has(id))continue;seen.add(id);
        out.push({word:w,first:positions[0],second:positions[1],kind:(pa<=pb&&pa+a.lhs.length>=pb+b.lhs.length)||(pb<=pa&&pb+b.lhs.length>=pa+a.lhs.length)?'inclusion':'overlap'});
      }
    }));return out;
  }
  audit(options={}){
    const resolutions=[],open=[];let incomplete=false;
    for(const a of this.ambiguities()){
      try{
        const t=this.wordType(a.word),left=this.replaceWord(a.word,...a.first,t.from),right=this.replaceWord(a.word,...a.second,t.from);
        const l=this.reduce(left,{...options,witness:true}),r=this.reduce(right,{...options,witness:true}),closed=l.normal.equals(r.normal);
        const entry={...a,closed,left:l.certificate,right:r.certificate};resolutions.push(entry);if(!closed)open.push(entry);
      }catch(e){if(!(e instanceof Budget))throw e;incomplete=true;resolutions.push({...a,status:'BUDGET_EXHAUSTED'});}
    }
    return {protocol:'NATIVE_CRITICAL_PAIR_AUDIT_V1',presentation:this.hash(),status:incomplete?'INCOMPLETE':open.length?'OPEN_CRITICAL_PAIRS':'CONFLUENT_BY_CHECKED_DIAMONDS',termination:'strict degree-lex decrease',criticalPairs:resolutions.length,resolutions,open};
  }
  complete({maxNewRules=12,...options}={}){
    bound(maxNewRules,'Invalid completion budget');const original=this.hash(),derivations=[];
    for(let k=0;;k++){
      const a=this.audit(options);
      if(a.status!=='OPEN_CRITICAL_PAIRS')return {original,status:a.status,audit:a,derivations};
      if(k>=maxNewRules)return {original,status:'COMPLETION_BUDGET_EXHAUSTED',audit:a,derivations};
      const pair=a.open[0],d=this.fromJSON(pair.left.output).minus(this.fromJSON(pair.right.output)),terms=d.entries();
      const [lead,c]=terms.at(-1);
      if(!lead.length)return {original,status:'LOCAL_UNIT_COLLAPSE',audit:a,derivations};
      const id='derived_'+this.rules.length;
      const rule={id,lhs:lead,rhs:terms.slice(0,-1).map(([w,v])=>[w,v.neg().div(c)]),source:'derived-critical-pair:'+this.hash()};
      derivations.push({before:this.spec(),pair,rule:{...rule,rhs:rule.rhs.map(([w,v])=>[w,v.toJSON()])}});this.addRule(rule);
    }
  }
}
function replayCompletion(original, packet){
  const s=new Presentation(original);if(s.hash()!==packet.original)throw new Error('Completion source mismatch');
  for(const item of packet.derivations){
    if(s.hash()!==digest(item.before))throw new Error('Wrong intermediate rule system');
    const a=item.pair;
    if(!s.ambiguities().some(x=>key([x.word,x.first,x.second])===key([a.word,a.first,a.second])))throw new Error('Fabricated critical pair');
    const t=s.wordType(a.word),l=s.replaceWord(a.word,...a.first,t.from),r=s.replaceWord(a.word,...a.second,t.from);
    if(!l.equals(s.fromJSON(a.left.input))||!r.equals(s.fromJSON(a.right.input)))throw new Error('Wrong critical fork');
    s.replay(a.left);s.replay(a.right);
    const d=s.fromJSON(a.left.output).minus(s.fromJSON(a.right.output)),terms=d.entries(),[w,c]=terms.at(-1)||[[]];
    if(!w.length)throw new Error('Unit-collapse completion is not admitted');
    const expected={id:'derived_'+s.rules.length,lhs:w,rhs:terms.slice(0,-1).map(([v,z])=>[v,z.neg().div(c).toJSON()]),source:'derived-critical-pair:'+s.hash()};
    if(digest(expected)!==digest(item.rule))throw new Error('Completion rule is not derived from its fork');
    s.addRule(item.rule);
  }
  if(s.hash()!==packet.audit.presentation)throw new Error('Wrong completed system hash');
  return s;
}
function bracket(a,b){return a.times(b).minus(b.times(a));}
function polynomialDerivative(e,values){
  const s=e.system;let out=s.zero(e.from,e.to);
  for(const [w,c] of e.entries())for(let i=0;i<w.length;i++){
    const d=values[w[i]];if(!d||d.system!==s)throw new TypeError('Missing typed generator derivative');
    const t=s.tokens[w[i]];if(d.from!==t.from||d.to!==t.to)throw new TypeError('Derivative changes arrow type');
    out=out.plus(new Expr(s,e.from,e.to,d.entries().map(([v,a])=>[w.slice(0,i).concat(v,w.slice(i+1)),a.mul(c)])));
  }return out;
}
function derivationGate(s,values){
  const a=s.audit();if(a.status!=='CONFLUENT_BY_CHECKED_DIAMONDS')return {status:'NEEDS_CONFLUENT_PRESENTATION'};
  for(const t of s.order.keys()){const d=values[t],sig=s.tokens[t];if(!d||d.system!==s||d.from!==sig.from||d.to!==sig.to)throw new TypeError('Missing typed generator derivative');}
  const defects=s.rules.map(r=>({rule:r.id,defect:s.reduce(polynomialDerivative(s.word(r.lhs).minus(r.rhs),values)).normal.toJSON()}));
  return {status:defects.every(x=>x.defect.terms.length===0)?'DERIVATION_DESCENDS':'RELATION_NOT_PRESERVED',defects};
}
function star(e,values){
  const s=e.system;let out=s.zero(e.to,e.from);
  for(const [w,c] of e.entries()){
    let x=s.one(e.from);for(const t of [...w].reverse()){
      const d=values[t],sig=s.tokens[t];if(!d||d.system!==s||d.from!==sig.to||d.to!==sig.from)throw new TypeError('Missing reversed dagger signature');if(d.entries().some(([v])=>s.wordCharge(v)!==-BigInt(sig.charge)))throw new TypeError('Dagger does not reverse declared residue');x=x.times(d);
    }out=out.plus(x.scale(c.dagger()));
  }return out;
}
function daggerGate(s,values){
  if(s.audit().status!=='CONFLUENT_BY_CHECKED_DIAMONDS')return {status:'NEEDS_CONFLUENT_PRESENTATION'};
  for(const t of s.order.keys())if(s.reduce(star(star(s.word([t]),values),values).minus(s.word([t]))).normal.terms.size)return {status:'DAGGER_NOT_INVOLUTIVE'};
  for(const r of s.rules)if(s.reduce(star(s.word(r.lhs).minus(r.rhs),values)).normal.terms.size)return {status:'DAGGER_NOT_RELATION_STABLE'};
  return {status:'DAGGER_DESCENDS'};
}
function seriesMultiply(a,b,degree){bound(degree,'Invalid formal degree');const s=a[0].system,at=a[0].from;if(a.concat(b).some(e=>e.system!==s||e.from!==at||e.to!==at))throw new TypeError('Formal series must share a native corner');return Array.from({length:degree+1},(_,k)=>{let out=s.zero(at);for(let i=0;i<=k;i++)if(a[i]&&b[k-i])out=out.plus(a[i].times(b[k-i]));return s.reduce(out).normal;});}
function formalExp(x,degree){bound(degree,'Invalid formal degree');if(x.from!==x.to)throw new TypeError('Exponential needs a corner');let p=x.system.one(x.from),factor=1n;const out=[p];for(let k=1;k<=degree;k++){p=x.system.reduce(p.times(x)).normal;factor*=BigInt(k);out.push(p.scale(new Cut(new oFraction(1n,factor))));}return out;}
const {F:oFraction} = require('./native_operator.cjs');
function formalLog(a,degree){bound(degree,'Invalid formal degree');const s=a[0].system,at=a[0].from;if(!a[0].equals(s.one(at)))throw new TypeError('Formal logarithm needs constant identity');const z=Array.from({length:degree+1},(_,k)=>k===0?s.zero(at):(a[k]||s.zero(at)));let p=z,out=Array.from({length:degree+1},()=>s.zero(at));for(let j=1;j<=degree;j++){for(let k=0;k<=degree;k++)out[k]=out[k].plus(p[k].scale(new Cut(new oFraction(j%2?1:-1,j))));p=seriesMultiply(p,z,degree);}return out.map(e=>s.reduce(e).normal);}
function proveEquality(a,b){
  const s=a.system;s.same(a,b);const audit=s.audit();
  if(audit.status!=='CONFLUENT_BY_CHECKED_DIAMONDS')return {status:'PRESENTATION_NOT_CERTIFIED',audit};
  const proof=s.reduce(a.minus(b),{witness:true});s.replay(proof.certificate);
  return {status:proof.normal.terms.size?'DISTINCT_IN_DECLARED_QUOTIENT':'EQUAL_IN_DECLARED_QUOTIENT',certificate:proof.certificate};
}
function solveCoefficients(target,basis){
  const s=target.system;if(!Array.isArray(basis)||!basis.length)throw new TypeError('Nonempty typed basis/template required');
  basis.forEach(e=>s.same(target,e));
  if(s.audit().status!=='CONFLUENT_BY_CHECKED_DIAMONDS')return {status:'PRESENTATION_NOT_CERTIFIED'};
  const t=s.reduce(target).normal,b=basis.map(e=>s.reduce(e).normal),keys=[...new Set([t,...b].flatMap(e=>[...e.terms.keys()]))];
  const rows=(keys.length?keys:[null]).map(k=>b.map(e=>e.terms.get(k)||ZERO).concat(t.terms.get(k)||ZERO));
  const rr=require('./native_operator.cjs').rref(rows),n=b.length;
  if(rr.matrix.some(r=>r.slice(0,n).every(x=>x.zero())&&!r[n].zero()))return {status:'NO_SOLUTION_IN_TEMPLATE'};
  const pivots=rr.pivots.filter(i=>i<n),free=Array.from({length:n},(_,i)=>i).filter(i=>!pivots.includes(i)),particular=Array(n).fill(ZERO);
  pivots.forEach((j,i)=>particular[j]=rr.matrix[i][n]);
  const nullspace=free.map(j=>{const v=Array(n).fill(ZERO);v[j]=ONE;pivots.forEach((p,i)=>v[p]=rr.matrix[i][j].neg());return v;});
  let x=s.zero(target.from,target.to);basis.forEach((e,i)=>x=x.plus(e.scale(particular[i])));
  const proof=proveEquality(target,x);if(proof.status!=='EQUAL_IN_DECLARED_QUOTIENT')throw new Error('Independent reconstructed template residual is nonzero');
  return {status:nullspace.length?'AFFINE_SOLUTION_FAMILY':'UNIQUE_SOLUTION',particular:particular.map(x=>x.toJSON()),nullspace:nullspace.map(v=>v.map(x=>x.toJSON())),certificate:proof.certificate};
}
function matrixValue(e,values,n){
  if(e.system.objects.length!==1)throw new TypeError('This backend is for one-object presentations');
  let out=zeros(n);for(const [w,c] of e.entries()){let a=identity(n);for(const t of w)a=mMul(a,values[t]);out=mAdd(out,mScale(a,c));}return out;
}
function representationGate(s,values,n,basis){
  for(const t of s.order.keys()){const m=matrix(values[t]);if(m.length!==n||m[0].length!==n)throw new TypeError('Matrix carrier mismatch');}
  const failed=s.rules.filter(r=>!mEqual(matrixValue(s.word(r.lhs),values,n),matrixValue(r.rhs,values,n))).map(r=>r.id);
  return {status:failed.length?'RELATION_VIOLATION':'REPRESENTATION_RESPECTS_RELATIONS',failed,
    independentProvidedImages:basis?rank(basis.map(w=>flatten(matrixValue(s.word(w),values,n)))):null,
    basisCount:basis?basis.length:null,faithfulnessRequiresProvedSpanningBasis:true};
}
/** A separate operational lane. Priorities and Lopa are not equation proofs. */
function ledgerBeforeLopa(state,markersToErase){
  if(!state||!Array.isArray(state.markers)||!Array.isArray(state.ledger)||!Array.isArray(markersToErase))throw new TypeError('Explicit marker and ledger state required');
  const erase=new Set(markersToErase);if(erase.size!==markersToErase.length||markersToErase.some(m=>!state.markers.includes(m)))throw new TypeError('Cannot erase an absent/duplicate marker');
  const kept=state.markers.filter(m=>!erase.has(m)),ledger=state.ledger.concat(markersToErase.map(m=>({kind:'lopa-memory',marker:m})));
  return {...state,markers:kept,ledger};
}
function resolvePriority(candidates){
  if(!Array.isArray(candidates)||!candidates.length)return {status:'NO_RULE'};
  if(candidates.some(x=>!Number.isSafeInteger(x.priority)))throw new TypeError('Exact declared priority required');
  const max=Math.max(...candidates.map(x=>x.priority)),win=candidates.filter(x=>x.priority===max);
  return win.length===1?{status:'POLICY_SELECTED',rule:win[0].id,provesConfluence:false}:{status:'OPEN_TIE',provesConfluence:false};
}
module.exports={Expr,Presentation,Budget,replayCompletion,digest,stable,bracket,polynomialDerivative,derivationGate,star,daggerGate,seriesMultiply,formalExp,formalLog,proveEquality,solveCoefficients,matrixValue,representationGate,ledgerBeforeLopa,resolvePriority};
