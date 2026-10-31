'use strict';
// Frozen predecessor's certified interpolation, adapted as a reusable module.
const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict');
const R=require('./source.cjs');
const {f,iv,point,add,sub,neg,mul,scale,div,mid,absmax,widen,overlap,contains,trig,decimal,PI,floor}=R;

function phaseEngine(directory,check){
 const nodes=new Map();
 for(let block=0;block<9;block++){
  const report=JSON.parse(fs.readFileSync(path.join(directory,`phase_${block}.result.json`)));
  assert.equal(report.status,'completed');assert.equal(report.backend,'cyclotomic64');
  for(const answer of report.answers){
   assert.equal(answer.status,'completed');assert.equal(answer.algebraically_exact,true);
   const [,label,axis]=answer.id.match(/^phase_(\d+)_(X|Y)$/),j=Number(label);
   const enclosure=answer.value.decimal_enclosure,value=iv(decimal(enclosure.lower),decimal(enclosure.upper));
   if(!nodes.has(j))nodes.set(j,{});
   assert(!nodes.get(j)[axis==='X'?'cos':'sin'],'duplicate phase answer');
   nodes.get(j)[axis==='X'?'cos':'sin']=value;
  }
 }
 assert.equal(nodes.size,17);
 for(let j=0;j<=16;j++){
  const node=nodes.get(j),native=trig(scale(PI,f(j).div(32)));
  check('retained_node_cos:'+j,overlap(node.cos,native.cos));
  check('retained_node_sin:'+j,overlap(node.sin,native.sin));
  check('retained_node_unit_circle:'+j,contains(add(mul(node.cos,node.cos),mul(node.sin,node.sin)),1));
 }
 function node(j){
  j=((j%64)+64)%64;
  if(j<=16)return nodes.get(j);
  if(j<32){const a=nodes.get(32-j);return{cos:neg(a.cos),sin:a.sin};}
  if(j<48){const a=nodes.get(j-32);return{cos:neg(a.cos),sin:neg(a.sin)};}
  const a=nodes.get(64-j);return{cos:a.cos,sin:neg(a.sin)};
 }
 const step=scale(PI,'1/32'),error=absmax(mul(mul(step,step),mul(step,step))).div(24);
 function phase(a){
  if(a.lo.zero()&&a.hi.zero())return{cos:point(1),sin:point(0)};
  const z=div(a,step),j=Number(floor(mid(z))),u=sub(z,point(j));
  assert(Number.isSafeInteger(j)&&f('-1/100').le(u.lo)&&u.hi.le('101/100'));
  const um1=sub(u,point(1)),um2=sub(u,point(2)),up1=add(u,point(1));
  const weights=[scale(mul(mul(u,um1),um2),'-1/6'),scale(mul(mul(up1,um1),um2),'1/2'),
   scale(mul(mul(up1,u),um2),'-1/2'),scale(mul(mul(up1,u),um1),'1/6')];
  let cos=point(0),sin=point(0);
  for(let n=0;n<4;n++){const value=node(j+n-1);cos=add(cos,mul(weights[n],value.cos));sin=add(sin,mul(weights[n],value.sin));}
  return{cos:widen(cos,error),sin:widen(sin,error)};
 }
 return{phase,error,nodes};
}
module.exports={phaseEngine};
