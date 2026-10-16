'use strict';
// RDE-1 finite exact witnesses. Regression evidence only: no continuum claim
// of ERGODICITY.md is certified here. Arithmetic: unchanged native engine.
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../../native-common-cone-limit/native_engine/native_operator.cjs'),f=N.F.of,C=N.Cut.of;
const checks=[],negative_controls=[],measurements={};
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const Z=f(0),ONE=f(1),sha0=x=>crypto.createHash('sha256').update(x).digest('hex');
const w4=p=>f(1+p*p).pow(4);                      // H^4 Fourier weight (1+p^2)^4
// ---- W1/W2: Fourier-multiplier commutators [M_D, S] -------------------------
const sharp=K=>q=>Math.abs(q)<=K?ONE:Z;
const vp=K=>q=>{const a=Math.abs(q);return a<=K?ONE:(a>=2*K?Z:f(2*K-a).div(K));}; // de la Vallee Poussin
const conv=(D,u)=>{const r=new Map();for(const[l,d]of D)for(const[q,v]of u){const p=l+q;r.set(p,(r.get(p)||Z).add(d.mul(v)));}return r;};
const apply=(m,u)=>new Map([...u].map(([q,v])=>[q,v.mul(m(q))]));
const sub=(a,b)=>{const r=new Map(a);for(const[k,v]of b)r.set(k,(r.get(k)||Z).sub(v));return r;};
const meq=(a,b)=>[...sub(a,b).values()].every(v=>v.zero());
const comm=(D,m,u)=>sub(conv(D,apply(m,u)),apply(m,conv(D,u)));   // M_D S u - S M_D u
const formula=(D,m,u,sgn=1)=>{const r=new Map();for(const[l,d]of D)for(const[q,v]of u){const p=l+q,t=sgn===1?m(q).sub(m(p)):m(p).sub(m(q));r.set(p,(r.get(p)||Z).add(d.mul(v).mul(t)));}return r;};
const trig=rows=>new Map(rows.map(([k,v])=>[k,f(v)]));
const Ds=[trig([[1,'1/2'],[-1,'1/2']]),trig([[2,'1/3'],[-2,'1/3'],[1,'-1/5'],[-1,'-1/5']])];
const us=[trig([[3,1],[-3,1],[7,'2/9'],[-7,'2/9'],[12,'-1/4'],[-12,'-1/4']]),trig([[5,'1/2'],[-5,'1/2'],[9,3],[-9,3]])];
for(const[i,D]of Ds.entries())for(const[j,u]of us.entries())for(const K of [4,6]){
 for(const[name,m]of [['sharp',sharp(K)],['vp',vp(K)]]){
  ck(`commutator_fourier_formula:${name}:${i}${j}K${K}`,meq(comm(D,m,u),formula(D,m,u)));
 }
}
no('commutator_formula_with_reversed_multiplier_difference',!meq(comm(Ds[0],vp(4),us[0]),formula(Ds[0],vp(4),us[0],-1)));
const der=u=>new Map([...u].map(([q,v])=>[q,v.mul(q)]));        // d/dx up to the common factor 2 pi i
ck('multiplier_commutes_with_derivative',meq(apply(vp(4),der(us[0])),der(apply(vp(4),us[0]))));
ck('multiplier_fixes_constants',vp(4)(0).eq(1)&&sharp(4)(0).eq(1));
// Cancellation identity behind Lemma 3: (V_K Phi)' D - V_K(Phi' D) = [M_D,V_K] Phi'
{const phi=us[1],D=Ds[1],m=vp(5);
 const lhs=sub(conv(D,der(apply(m,phi))),apply(m,conv(D,der(phi))));ck('gradient_cancellation_identity',meq(lhs,comm(D,m,der(phi))));
 no('gradient_cancellation_without_commuting_derivative',!meq(sub(conv(D,apply(m,phi)),apply(m,conv(D,der(phi)))),comm(D,m,der(phi))));}
// W2: band norm with the borderline coefficient profile |u_q| = q^-4 (u in H^s for s<7/2 only)
const borderline=Q=>{const u=new Map();for(let q=1;q<=Q;q++){const v=f(1).div(f(q).pow(4));u.set(q,v);u.set(-q,v);}return u;};
const h4norm2=v=>[...v].reduce((s,[p,c])=>s.add(w4(p).mul(c.pow(2))),Z);
const band={sharp:{},vp:{}};const D0=Ds[0];
for(const K of [8,16,32,64,128]){
 const u=borderline(2*K+4);     // exact: both commutators vanish beyond |q| = 2K+1
 band.sharp[K]=h4norm2(comm(D0,sharp(K),u));band.vp[K]=h4norm2(comm(D0,vp(K),u));
}
const Ks=[8,16,32,64,128];
for(let i=1;i<Ks.length;i++){
 ck('vp_band_norm_decreases:'+Ks[i],band.vp[Ks[i]].le(band.vp[Ks[i-1]].mul(f('3/5'))));
 no('sharp_band_norm_decays_like_vp:'+Ks[i],!band.sharp[Ks[i]].mul(Ks[i]).le(f(2)));
}
// Lemma 3 bound, low-frequency case (D has modes |l|<=2<K/2, so the tail term vanishes):
// ||[M_D,V_K]w||_{H4}^2 <= 48^2 ||D||_{A1}^2 ||1_{K/2<|q|<5K/2} w||_{H3}^2
{const A1=[...D0].reduce((s,[l,d])=>s.add(d.mul(Math.abs(l))),Z);
 for(const K of Ks){const u=borderline(2*K+4);let b3=Z;for(const[q,v]of u)if(2*Math.abs(q)>K&&2*Math.abs(q)<5*K)b3=b3.add(f(1+q*q).pow(3).mul(v.pow(2)));
  ck('lemma3_vp_commutator_bound:'+K,band.vp[K].le(f(48*48).mul(A1.pow(2)).mul(b3)));}}
ck('vp_band_norm_times_K_bounded',Ks.every(K=>band.vp[K].mul(K).le(f(2))));
ck('sharp_band_norm_bounded_below_by_one',Ks.every(K=>ONE.le(band.sharp[K])));
// truncation error (I-V_K)phi in H^4 for phi_hat_q = q^-5 (phi in H^s, s<9/2)
const tail=K=>{let s=Z;for(let q=K+1;q<2*K;q++){const t=f(1).sub(vp(K)(q));s=s.add(w4(q).mul(t.pow(2)).div(f(q).pow(10)).mul(2));}for(let q=2*K;q<=40*K;q++)s=s.add(w4(q).div(f(q).pow(10)).mul(2));return s;};
const tl=[8,16,32].map(tail);ck('vp_truncation_error_decreases',tl[1].le(tl[0])&&tl[2].le(tl[1]));
const dec=x=>{const q=x.n*10n**9n/x.d,t=q.toString().padStart(10,'0');return t.slice(0,-9)+'.'+t.slice(-9);}; // truncated decimal, display only
measurements.band_norm_H4_squared={profile:'u_hat(q)=|q|^-4, D=cos(2 pi x), exact rational recorded via sha256, truncated decimal shown',
 sharp:Object.fromEntries(Ks.map(K=>[K,dec(band.sharp[K])])),vallee_poussin:Object.fromEntries(Ks.map(K=>[K,dec(band.vp[K])])),
 exact_sha256:sha0(JSON.stringify(Ks.map(K=>[band.sharp[K].toJSON(),band.vp[K].toJSON()])))};
// ---- W3: finite Gaussian divergence-closure model --------------------------
// E over independent N(0,v1) x N(0,v2); constraint E[F(-d_j psi + x_j psi / v_j)] = 0.
const mom=(v,k)=>{if(k%2)return Z;let r=ONE;for(let i=1;i<k;i+=2)r=r.mul(i);return r.mul(f(v).pow(k/2));};
const monos=d=>{const r=[];for(let i=0;i<=d;i++)for(let j=0;i+j<=d;j++)r.push([i,j]);return r;};
function nullity(vs,dirs,scoreVar){
 const unk=monos(3),rows=[];
 for(const[a,b]of monos(4))for(const j of dirs){
  rows.push(unk.map(([i,k])=>{
   // F-monomial x1^i x2^k times (-d_j psi + x_j psi/scoreVar_j), psi=x1^a x2^b
   let t=Z;const e=[i+a,k+b];
   const dj=j===0?a:b;
   if(dj>0){const ee=[...e];ee[j]-=1;t=t.sub(f(dj).mul(mom(vs[0],ee[0])).mul(mom(vs[1],ee[1])));}
   const ee=[...e];ee[j]+=1;t=t.add(mom(vs[0],ee[0]).mul(mom(vs[1],ee[1])).div(scoreVar[j]));
   return C(t);
  }));
 }
 return unk.length-N.rank(N.matrix(rows));
}
ck('closure_all_directions_leaves_constants_only',nullity([2,2],[0,1],[2,2])===1);
no('closure_one_direction_leaves_constants_only',nullity([2,2],[0],[2,2])!==1);
ck('closure_one_direction_nullity_counts_functions_of_x2',nullity([2,2],[0],[2,2])===4);
no('closure_with_undoubled_score_keeps_constants',nullity([2,2],[0,1],[1,1])!==1);
ck('closure_with_undoubled_score_kills_constants',nullity([2,2],[0,1],[1,1])===0);
// ---- W4: sign sectors, commutant, populations and unitary rescaling --------
const blk=(A,B)=>{const z=[[0,0],[0,0]];return N.matrix([...A.map((r,i)=>[...r,...z[i]]),...B.map((r,i)=>[...z[i],...r])]);};
const sx=[[0,1],[1,0]],sz=[[1,0],[0,-1]],I2=[[1,0],[0,1]],O2=[[0,0],[0,0]];
const obs=[blk(sx,sx),blk(sz,sz),blk(I2,O2)];
const commDim=gs=>{const eqs=[],I4=N.identity(4);for(const g of gs)eqs.push(...N.sub(N.kron(g,I4),N.kron(I4,N.dagger(g))));return 16-N.rank(N.matrix(eqs));}; // vec(gX-Xg), row-major, real g
ck('sector_commutant_dimension_two',commDim(obs)===2);
const mix=N.matrix([[0,0,1,0],[0,0,0,1],[1,0,0,0],[0,1,0,0]]);
ck('sign_mixing_observable_collapses_commutant',commDim([...obs,mix])===1);
const Bw=(cp,cm)=>N.matrix([[cp,0,0,0],[0,cp,0,0],[0,0,cm,0],[0,0,0,cm]]);
const Pp=obs[2],v=[1,0,1,0];
const pop=B=>{const num=N.mul(N.matrix([v]),N.mul(B,N.mul(Pp,N.dagger(N.matrix([v])))))[0][0],den=N.mul(N.matrix([v]),N.mul(B,N.dagger(N.matrix([v]))))[0][0];return num.rad.div(den.rad);};
ck('population_depends_on_sector_weights',pop(Bw(1,1)).eq('1/2')&&pop(Bw(1,4)).eq('1/5'));
const W=Bw(1,2);
ck('rescaling_W_implements_weighted_form',N.equal(N.mul(N.dagger(W),W),Bw(1,4)));
for(const[i,O]of obs.entries())ck('rescaling_W_intertwines_declared_observable:'+i,N.equal(N.mul(W,O),N.mul(O,W)));
no('rescaling_W_intertwines_sign_mixing_observable',!N.equal(N.mul(W,mix),N.mul(mix,W)));
const e3=[0,0,1,0],norm=(B,x)=>N.mul(N.matrix([x]),N.mul(B,N.dagger(N.matrix([x]))))[0][0].rad;
ck('isometric_spectator_contract_holds_for_equal_weights',norm(Bw(1,1),e3).eq(1));
no('isometric_spectator_contract_holds_for_unequal_weights',!norm(Bw(1,4),e3).eq(1));
// ---- certificate ------------------------------------------------------------
const sha=sha0;
const payload={protocol:'RDE1_RESIDUAL_DIFF_ERGODICITY_V1',status:'PASS_FINITE_WITNESSES_FOR_VP_DIVERGENCE_CLOSURE',
 source_pins_sha256:sha(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),checks,negative_controls,measurements,
 results:{continuum_proofs_machine_certified:false,ergodicity_claim:'written proof in ERGODICITY.md, not certified here',
  sharp_dirichlet_D12_route:'band norm does not decay in witness; L1 divergence of sharp sums not decided',
  vallee_poussin_D12_route:'band norm decays in witness; continuum estimate written'},
 summary:{exact_checks:checks.length,negative_controls:negative_controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',hash=sha(bytes),cp=path.join(__dirname,'CERTIFICATE.json'),hp=path.join(__dirname,'EXPECTED.sha256');
if(process.argv.includes('--write')){fs.writeFileSync(cp,bytes);fs.writeFileSync(hp,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cp,'utf8'),bytes);assert.equal(fs.readFileSync(hp,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash,band:measurements.band_norm_H4_squared}));
