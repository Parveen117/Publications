'use strict';
// Local regression witness (not a proof): exact Q+iQ arithmetic with the repo's native engine.
const N=require(process.argv[2]+'/papers/native-common-cone-limit/native_engine/native_operator.cjs');
const C=N.Cut.of, Z=N.ZERO, F=N.F.of;
// psi real trig polynomial, Fourier coeffs psi_m (m != 0), psi_{-m}=conj(psi_m)
const coeffs=new Map([[1,C(['1/2','1/3'])],[2,C(['-2/5','1/7'])],[3,C(['1/4','0'])]]);
for(const [m,v] of [...coeffs]) coeffs.set(-m,v.dagger());
const ps=m=>coeffs.get(m)||Z;
// K_1 e_n = sum_m (m/n) psi_m (e_{m+n} - e_m), projected off the zero mode
const K=(p,n)=>{ if(p===0||n===0) return Z; return C(p-n).mul(ps(p-n)).sub(C(p).mul(ps(p))).div(C(n)); };
const exact=[...coeffs.values()].reduce((s,v)=>s.add(C(v.norm2())),Z).neg(); // -int psi^2
const out={claim:'tr(K1^2) -> -int psi^2 (second-order regularized Liouville identity at sigma=0)', limit:exact.toJSON(), rows:[]};
for(const L of [8,16,32,64]){
  let t=Z; for(let n=-L;n<=L;n++) for(let p=-L;p<=L;p++){ const a=K(n,p); if(a.zero()) continue; const b=K(p,n); if(!b.zero()) t=t.add(a.mul(b)); }
  const err=t.sub(exact); out.rows.push({L, error_times_L: err.mul(C(L)).toJSON(), error: err.toJSON()});
}
// Hida S-transform functional equation exponent identity on random rational vectors
let ok=true; const r=(i)=>F(((i*37)%23-11)+'/'+((i*13)%7+2));
for(let trial=0;trial<50;trial++){ const d=5; const eta=[],k=[]; for(let i=0;i<d;i++){eta.push(r(trial*d+i+1));k.push(r(trial*d+i+101));}
  const dot=(a,b)=>a.reduce((s,x,i)=>s.add(x.mul(b[i])),F(0));
  const em=eta.map((x,i)=>x.sub(k[i].div(2)));
  const lhs=dot(k,eta).div(2).sub(dot(k,k).div(8)).add(dot(em,em).div(2)).sub(dot(eta,eta).div(2));
  ok=ok&&lhs.eq(0); }
out.hida_functional_equation_exponent_zero_50_trials=ok;
console.log(JSON.stringify(out,null,1));
