'use strict';
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const f=N.F.of,C=N.Cut.of;
const key=v=>JSON.stringify([...v].sort());
const state=(v=[],c=1)=>new Map([[key(v),C(c)]]),vac=state();
function plus(...xs){const out=new Map();for(const x of xs)for(const[k,c]of x){const z=(out.get(k)||N.ZERO).add(c);if(z.zero())out.delete(k);else out.set(k,z);}return out;}
const scale=(x,c)=>new Map([...x].map(([k,z])=>[k,z.mul(c)]).filter(([,z])=>!z.zero()));
const minus=(x,y)=>plus(x,scale(y,-1)),eq=(x,y)=>minus(x,y).size===0;
function create(x,v){return plus(...[...x].map(([k,c])=>state([...JSON.parse(k),v],c)));}
function deriv(x,v){return plus(...[...x].map(([k,c])=>{const a=JSON.parse(k),r=a.filter(z=>z===v).length;if(!r)return new Map();a.splice(a.indexOf(v),1);return state(a,c.mul(r));}));}
function current(s,n,x,t){if(n===0)return scale(x,t);return n<0?create(x,s+':'+(-n)):scale(deriv(x,s+':'+n),2*n);}
function bound(x){let b=0;for(const k of x.keys())for(const v of JSON.parse(k))b=Math.max(b,Number(v.split(':')[1]));return b;}
function T(s,m,x,t){let out=new Map();const B=bound(x)+Math.abs(m);for(let j=-B;j<=B;j++){
 const a=m-j,z=a>0&&j<0?current(s,j,current(s,a,x,t),t):current(s,a,current(s,j,x,t),t);out=plus(out,z);
}return scale(out,'1/4');}
const A=(s,m,x,t)=>plus(T(s,m,x,t),scale(current(s,m,x,t),new N.Cut(0,f(m).div(2))),m===0?scale(x,'1/4'):new Map());
const L=(m,x,t,u=t)=>minus(A('u',m,x,t),A('v',-m,x,u));
const sourceL=(m,x,t,u=t)=>plus(T('u',m,x,t),scale(T('v',-m,x,u),-1),scale(plus(current('u',m,x,t),current('v',-m,x,u)),new N.Cut(0,f(m).div(2))));
function weight(vars){const counts=new Map();for(const v of vars)counts.set(v,(counts.get(v)||0)+1);let w=f(1);for(const[v,r]of counts){const n=Number(v.split(':')[1]);for(let j=1;j<=r;j++)w=w.mul(2*n).mul(j);}return w;}
function inner(x,y){let z=N.ZERO;for(const[k,c]of x)if(y.has(k))z=z.add(c.dagger().mul(y.get(k)).mul(weight(JSON.parse(k))));return z;}
function tensor(x,y){return plus(...[...x].flatMap(([a,c])=>[...y].map(([b,d])=>state([...JSON.parse(a),...JSON.parse(b)],c.mul(d)))));}
function partitions(n,max=n){if(n===0)return [[]];const out=[];for(let k=Math.min(n,max);k>=1;k--)for(const p of partitions(n-k,k))out.push([k,...p]);return out;}
function descendant(s,p,t){let x=vac;for(const n of [...p].reverse())x=A(s,-n,x,t);return x;}
const grade=(k,s)=>JSON.parse(k).filter(v=>v.startsWith(s+':')).reduce((a,v)=>a+Number(v.split(':')[1]),0);
const project=(x,max)=>new Map([...x].filter(([k])=>grade(k,'u')<=max&&grade(k,'v')<=max));
function paired(level,t){const ps=partitions(level),us=ps.map(p=>descendant('u',p,t)),vs=ps.map(p=>descendant('v',p,t));
 const G=us.map(x=>us.map(y=>inner(x,y))),iv=N.inverse(G);let out=new Map();
 for(let a=0;a<ps.length;a++)for(let b=0;b<ps.length;b++)out=plus(out,scale(tensor(us[a],vs[b]),iv[a][b]));
 return {out,G,iv,ps,us,vs};
}

module.exports={N,f,C,key,state,vac,plus,scale,minus,eq,current,T,A,L,inner,tensor,partitions,descendant,weight,grade,paired};
