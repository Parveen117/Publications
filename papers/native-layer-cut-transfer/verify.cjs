'use strict';
// LCT-1 finite exact witnesses for THEOREM.md. Regression evidence only:
// the theorems are written algebraic/analytic proofs; no finite check here
// certifies a statement about the actual xi source or any all-depth claim.
// Arithmetic: the unchanged native engine over Q + iota Q (Cut). No floats.
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const C=N.Cut.of,Z=N.ZERO,O=N.ONE;
const checks=[],negative_controls=[],records={};
const ck=(name,ok)=>{assert.ok(ok,'CHECK FAILED: '+name);checks.push({name,passed:true});};
const no=(name,rejected)=>{assert.ok(rejected,'CONTROL NOT REJECTED: '+name);negative_controls.push({name,rejected:true});};
const s=x=>JSON.stringify(x.toJSON());

// ---- polynomials in lambda over Q + iota Q (index = degree) ------------------
const trim=p=>{p=p.slice();while(p.length>1&&p[p.length-1].zero())p.pop();return p.length?p:[Z];};
const padd=(a,b)=>trim(Array.from({length:Math.max(a.length,b.length)},(_,i)=>(a[i]||Z).add(b[i]||Z)));
const pscale=(a,c)=>trim(a.map(x=>x.mul(C(c))));
const psub=(a,b)=>padd(a,pscale(b,-1));
const pmul=(a,b)=>{const r=Array(a.length+b.length-1).fill(Z);a.forEach((x,i)=>b.forEach((y,j)=>{r[i+j]=r[i+j].add(x.mul(y));}));return trim(r);};
const isZeroP=p=>trim(p).every(x=>x.zero());
const peq=(a,b)=>isZeroP(psub(a,b));
const peval=(p,x)=>{x=C(x);return p.reduceRight((acc,c)=>acc.mul(x).add(c),Z);};
const mono=k=>trim([...Array(k).fill(Z),O]);
const deg=p=>{p=trim(p);return p.length===1&&p[0].zero()?-1:p.length-1;};
const lowZeros=(p,k)=>{for(let i=0;i<k;i++)if(p[i]&&!p[i].zero())return false;return true;};
const shiftDown=(p,k)=>trim(p.slice(k));
const isReal=x=>x.turn.zero();
const sgn=x=>{assert.ok(isReal(x));return x.rad.zero()?0:(N.F.of(0).le(x.rad)?1:-1);};

// ---- a source is a finite bilinear functional L[g]=sum w g(u), u,w in Q+iQ ---
function source(name,pts,wts){
  const u=pts.map(C),w=wts.map(C);
  const total=w.reduce((a,b)=>a.add(b),Z);assert.ok(total.eq(O),name+' normalized');
  const L=g=>u.reduce((acc,x,i)=>acc.add(w[i].mul(g(x))),Z);
  const mom=n=>L(x=>pow(x,n));
  const D=u.reduce((acc,x)=>pmul(acc,[O,x]),[O]);                 // prod (1+lambda u)
  const Dm=u.map((_,i)=>u.reduce((acc,x,j)=>j===i?acc:pmul(acc,[O,x]),[O])); // D/(1+lambda u_i)
  return {name,u,w,L,mom,D,Dm};
}
const pow=(x,n)=>{let r=O;for(let k=0;k<n;k++)r=r.mul(x);return r;};
const polyY=(p,y)=>p.reduceRight((acc,c)=>acc.mul(y).add(c),Z); // p as polynomial in y
function det(M){M=M.map(r=>r.slice());const n=M.length;let d=O;
  for(let j=0;j<n;j++){let k=-1;for(let i=j;i<n;i++)if(!M[i][j].zero()){k=i;break;}
    if(k<0)return Z;if(k!==j){[M[j],M[k]]=[M[k],M[j]];d=d.mul(-1);}
    const piv=M[j][j];d=d.mul(piv);
    for(let i=j+1;i<n;i++){const q=M[i][j].div(piv);if(!q.zero())M[i]=M[i].map((x,t)=>x.sub(q.mul(M[j][t])));}}
  return d;}
const hankel=(S,n,shift=0)=>Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>S.mom(i+j+shift)));

// ---- the cut at depth N: monic orthogonal P_N, reversed Q_N, Pade numerator R_N
function cut(S,Nd,perturb=null){
  let p;
  if(Nd===0)p=[O];
  else{const H=hankel(S,Nd),rhs=Array.from({length:Nd},(_,i)=>[S.mom(Nd+i).mul(-1)]);
    const x=N.mul(N.inverse(H),rhs).map(r=>r[0]);p=[...x,O];}
  if(perturb)p=p.map((c,i)=>i===perturb.index?c.add(C(perturb.by)):c);
  const Q=trim(Array.from({length:Nd+1},(_,j)=>p[Nd-j].mul(j%2?-1:1)));  // Q_N(lambda)=prod(1+lambda y_j)
  const fser=Array.from({length:Nd},(_,n)=>S.mom(n).mul(n%2?-1:1));       // f series to order N-1
  const R=trim(pmul(Q,fser.length?fser:[Z]).slice(0,Math.max(Nd,1))).slice(0,Math.max(Nd,1));
  const Rt=Nd===0?[Z]:trim(R);
  const ED=S.u.reduce((acc,x,i)=>padd(acc,pscale(S.Dm[i],S.w[i].mul(polyY(p,x).mul(polyY(p,x))))),[Z]); // E_N * D
  const h=S.L(y=>polyY(p,y).mul(polyY(p,y)));
  return {N:Nd,p,Q,R:Rt,ED,h};
}
const fD=S=>S.u.reduce((acc,_,i)=>padd(acc,pscale(S.Dm[i],S.w[i])),[Z]); // f * D

// ---- sources -----------------------------------------------------------------
const S1=source('S1_positive_six_point',['1/2','1/3','1/5','2/7','3/4','1/11'],['1/10','2/10','3/10','1/10','2/10','1/10']);
const S2=source('S2_indefinite_real',['1/2','1/3','1/5','3/4'],['3/5','3/5','-2/5','1/5']);
const S3=source('S3_conjugate_pair_comparison',[['1/2','1/2'],['1/2','-1/2'],'1/4','1/8','1/16'],[['1/8','-1/8'],['1/8','1/8'],'1/4','1/4','1/4']);
const sources=[S1,S2,S3];
records.sources={};

// ---- LCT-T1: the universal generator of the layer tower ----------------------
for(const S of sources){
  const c1=S.mom(1);
  const du=S.u.map(x=>x.mul(O.sub(x)));              // D u = u(1-u)
  const dw=S.w.map((w,i)=>w.mul(c1.sub(S.u[i])));    // D w = (c_1-u) w
  ck(`T1 ${S.name}: generator preserves normalization sum Dw = 0`,dw.reduce((a,b)=>a.add(b),Z).zero());
  for(let n=0;n<=8;n++){
    const lhs=S.u.reduce((acc,x,i)=>acc.add(dw[i].mul(pow(x,n))).add(n?S.w[i].mul(C(n)).mul(pow(x,n-1)).mul(du[i]):Z),Z);
    const rhs=C(n).add(c1).mul(S.mom(n)).sub(C(n+1).mul(S.mom(n+1)));
    ck(`T1 ${S.name}: d c_${n}/d log a = (n+c_1)c_n-(n+1)c_(n+1)`,lhs.eq(rhs));
  }
  if(S===S3)for(let n=0;n<=10;n++)ck(`T1 S3: conjugate cancellation makes c_${n} real`,isReal(S.mom(n)));
}
{ // exact anchor transport of a positive pole source: q(t)=sum nu/(t+mu)
  const mu=['2','5','9','17'].map(C),nu=['1','3','2','1'].map(C);
  const q=t=>mu.reduce((acc,m,i)=>acc.add(nu[i].div(C(t).add(m))),Z);
  const chart=a=>({u:mu.map(m=>C(a).div(C(a).add(m))),w:mu.map((m,i)=>nu[i].div(q(a).mul(C(a).add(m))))});
  const A=chart(3),rho=C('7/3'),B=chart(7);
  const uT=A.u.map(x=>rho.mul(x).div(O.add(rho.sub(O).mul(x))));
  const wraw=A.w.map((w,i)=>w.div(O.add(rho.sub(O).mul(A.u[i]))));const tot=wraw.reduce((a,b)=>a.add(b),Z);
  ck('T1 transport: Moebius chart u -> rho u/(1+(rho-1)u) gives the b-chart spectral values',uT.every((x,i)=>x.eq(B.u[i])));
  ck('T1 transport: weight tilt w/(1+(rho-1)u), renormalized, gives the b-chart weights',wraw.every((w,i)=>w.div(tot).eq(B.w[i])));
  ck('T1 transport: the tilt total is the normalization cocycle q(b)/q(a)',tot.eq(q(7).div(q(3))));
  no('T1 control: the untilted transported weights are not the b-chart weights',!A.w.every((w,i)=>w.eq(B.w[i])));
}

// ---- LCT-T2..T5 on each source --------------------------------------------------
const lam2N=n=>mono(2*n);
for(const S of sources){
  const FD=fD(S),rec=[];const maxN=S===S1?5:(S===S2?3:4);
  for(let n=1;n<=maxN;n++){
    const K=cut(S,n),K1=cut(S,n+1);
    const left=psub(pmul(FD,pmul(K.Q,K.Q)),pmul(pmul(K.R,K.Q),S.D));
    // T2: exact factorization  f - R/Q = lambda^{2N} E_N / Q^2
    ck(`T2 ${S.name} N=${n}: (f-R_N/Q_N) Q_N^2 = lambda^(2N) E_N`,peq(left,pmul(lam2N(n),K.ED)));
    ck(`T2 ${S.name} N=${n}: Q_N(0)=1 and deg Q_N=N, deg R_N<=N-1`,K.Q[0].eq(O)&&deg(K.Q)===n&&deg(K.R)<=n-1);
    ck(`T2 ${S.name} N=${n}: Q_N has real coefficients`,K.Q.every(isReal)&&K.R.every(isReal));
    // T2: the correction at the reference plane is h_N = det H_{N+1}/det H_N = E_N(0)
    const dN=det(hankel(S,n)),dN1=det(hankel(S,n+1));
    ck(`T2 ${S.name} N=${n}: E_N(0)=L[P_N^2]=h_N`,peval(K.ED,0).eq(K.h)); // D(0)=1
    ck(`T2 ${S.name} N=${n}: h_N = det H_(N+1)/det H_N`,K.h.mul(dN).eq(dN1));
    // T3: tower above the cut is re-encoded unipotently by e_(N,k)=L[P_N^2 y^k]
    for(let k=0;k<=3;k++){
      const e=S.L(y=>polyY(K.p,y).mul(polyY(K.p,y)).mul(pow(y,k)));
      let comb=Z;for(let i=0;i<=n;i++)for(let j=0;j<=n;j++)comb=comb.add(K.p[i].mul(K.p[j]).mul(S.mom(i+j+k)));
      ck(`T3 ${S.name} N=${n} k=${k}: e_(N,k) = sum p_i p_j c_(i+j+k)`,e.eq(comb));
      let lower=Z;for(let i=0;i<=n;i++)for(let j=0;j<=n;j++)if(!(i===n&&j===n))lower=lower.add(K.p[i].mul(K.p[j]).mul(S.mom(i+j+k)));
      ck(`T3 ${S.name} N=${n} k=${k}: c_(2N+k) recovered from e_(N,k) and lower layers`,e.sub(lower).eq(S.mom(2*n+k)));
    }
    // T4: one-step depth change with pre-supplied factor g_N = Q_N/Q_(N+1)
    const lhs3=pmul(K.ED,pmul(K1.Q,K1.Q));
    const rhs3=padd(pscale(pmul(pmul(K.Q,K1.Q),S.D),K.h),pmul(mono(2),pmul(pmul(K.Q,K.Q),K1.ED)));
    ck(`T4 ${S.name} N=${n}: E_N Q_(N+1)^2 = h_N Q_N Q_(N+1) + lambda^2 Q_N^2 E_(N+1)`,peq(lhs3,rhs3));
    ck(`T4 ${S.name} N=${n}: R_(N+1)Q_N - R_N Q_(N+1) = h_N lambda^(2N)`,peq(psub(pmul(K1.R,K.Q),pmul(K.R,K1.Q)),pscale(lam2N(n),K.h)));
    no(`T4 control ${S.name} N=${n}: a wrong reference-plane innovation h_N+1/1000 breaks the step`,
       !peq(lhs3,padd(pscale(pmul(pmul(K.Q,K1.Q),S.D),K.h.add(C('1/1000'))),pmul(mono(2),pmul(pmul(K.Q,K.Q),K1.ED)))));
    // T4: two-step composition equals the direct depth-change formula (cocycle)
    if(n+2<=maxN+1){
      const K2=cut(S,n+2);
      const inh=psub(pmul(K2.R,K.Q),pmul(K.R,K2.Q));
      ck(`T4 ${S.name} N=${n}: f_(N+2)-f_N vanishes to order 2N`,lowZeros(inh,2*n));
      const direct=pmul(pmul(K.Q,pmul(K1.Q,K1.Q)),shiftDown(inh,2*n));
      const composed=padd(pscale(pmul(pmul(K.Q,K1.Q),K2.Q),K.h),pscale(pmul(mono(2),pmul(pmul(K.Q,K.Q),K1.Q)),K1.h));
      ck(`T4 ${S.name} N=${n}: composed one-step maps equal the direct N->N+2 transfer`,peq(direct,composed));
    }
    // T5: two sides of the reference plane; sign carried only by E_N
    const sides=['1/2','-1/2','1/3','-1/3','2','-3/4'];
    const sideSigns=[];
    for(const l of sides){
      const Qv=peval(K.Q,l),Dv=peval(S.D,l);if(Qv.zero()||Dv.zero())continue;
      const corr=peval(left,l).div(Qv.mul(Qv).mul(Dv));              // f - f_N
      const fac=pow(C(l),2*n).div(Qv.mul(Qv));                         // K_N(lambda) >= 0 on reals
      const Ev=peval(K.ED,l).div(Dv);
      ck(`T5 ${S.name} N=${n} lambda=${l}: f-f_N = K_N E_N with K_N>=0`,corr.eq(fac.mul(Ev))&&sgn(fac)>=0);
      ck(`T5 ${S.name} N=${n} lambda=${l}: sign(f-f_N)=sign(E_N)`,sgn(corr)===sgn(Ev));
      sideSigns.push({lambda:l,E_sign:sgn(Ev)});
    }
    no(`T5 control ${S.name} N=${n}: the cut alone does not transfer the tower, f-f_N != 0 at lambda=1/2`,!peval(left,'1/2').zero());
    // flatness: K_N = lambda^(2N)(1+O(lambda)), no lower order
    no(`T5 control ${S.name} N=${n}: flatness order 2N-1 instead of 2N fails the factorization`,!peq(left,pmul(mono(2*n-1),K.ED)));
    no(`T5 control ${S.name} N=${n}: an unsquared factor lambda^(2N)/Q_N fails`,!peq(psub(pmul(FD,K.Q),pmul(K.R,S.D)),pmul(lam2N(n),K.ED)));
    // non-orthogonal cut: perturbing P_N destroys exactness
    const Kp=cut(S,n,{index:Math.max(0,n-1),by:'1/7'});
    const leftp=psub(pmul(FD,pmul(Kp.Q,Kp.Q)),pmul(pmul(Kp.R,Kp.Q),S.D));
    no(`T2 control ${S.name} N=${n}: a non-orthogonal cut polynomial breaks the factorization`,!peq(leftp,pmul(lam2N(n),Kp.ED)));
    rec.push({N:n,h_N:K.h.toJSON(),h_sign:sgn(K.h),det_H_N:dN.toJSON(),Q_N:K.Q.map(x=>x.toJSON()),side_signs_of_E_N:sideSigns});
  }
  records.sources[S.name]={points:S.u.map(x=>x.toJSON()),weights:S.w.map(x=>x.toJSON()),cuts:rec};
}

// ---- saturation and existence boundaries --------------------------------------
{ const K6=cut(S1,6);
  ck('T3 S1 saturation: six support points give h_6 = 0',K6.h.zero());
  ck('T3 S1 saturation: E_6 = 0 identically, so f = f_6 exactly',isZeroP(K6.ED));
  let threw=false;try{cut(S1,7);}catch(e){threw=true;}
  no('T2 control S1: no orthogonal P_7 exists (det H_7 = 0); the cut depth cannot exceed the support',threw&&det(hankel(S1,7)).zero());
}
// ---- the sign lives only in E_N: indefinite sources ----------------------------
{ const negS2=records.sources[S2.name].cuts.some(r=>r.h_sign<0||r.side_signs_of_E_N.some(x=>x.E_sign<0));
  const negS3=records.sources[S3.name].cuts.some(r=>r.h_sign<0||r.side_signs_of_E_N.some(x=>x.E_sign<0));
  const posS1=records.sources[S1.name].cuts.every(r=>r.h_sign>0&&r.side_signs_of_E_N.every(x=>x.E_sign>0));
  ck('T5 S1 positive source: every cut correction is positive on both sides of the plane',posS1);
  no('T5 control S2: the universal factor K_N>=0 does not give a sign, the indefinite real source has a negative correction',negS2);
  no('T5 control S3: the conjugate-pair comparison source has a negative correction although every moment is real',negS3);
}
// ---- bilinear, not absolute squares ----------------------------------------
{ const n=2,K=cut(S3,n),FD=fD(S3);
  const EDabs=S3.u.reduce((acc,x,i)=>{const v=polyY(K.p,x);return padd(acc,pscale(S3.Dm[i],S3.w[i].mul(v.mul(v.dagger()))));},[Z]); // w |P(u)|^2
  const left=psub(pmul(FD,pmul(K.Q,K.Q)),pmul(pmul(K.R,K.Q),S3.D));
  no('T2 control S3: replacing the bilinear w P(u)^2 by w |P(u)|^2 breaks the factorization',!peq(left,pmul(lam2N(n),EDabs)));
}

// ---- output -------------------------------------------------------------------
const cert={protocol:'LCT1_LAYER_CUT_TRANSFER_V1',
  scope:'Finite exact regression witnesses for LCT-1. Algebraic identities on finite rational and Gaussian-rational bilinear sources. Not a certificate about the xi source, not an all-depth statement.',
  engine:'papers/native-common-cone-limit/native_engine/native_operator.cjs (Cut = Q + iota Q), unchanged',
  counts:{checks:checks.length,negative_controls:negative_controls.length},checks,negative_controls,records};
const out=JSON.stringify(cert,null,1)+'\n';
const file=path.join(__dirname,'CERTIFICATE.json'),exp=path.join(__dirname,'EXPECTED.sha256');
const digest=crypto.createHash('sha256').update(out).digest('hex');
const mode=process.argv[2];
if(mode==='--write'){fs.writeFileSync(file,out);fs.writeFileSync(exp,digest+'  CERTIFICATE.json\n');}
else if(mode==='--check'){assert.equal(fs.readFileSync(file,'utf8'),out,'CERTIFICATE.json differs');assert.equal(fs.readFileSync(exp,'utf8').split(/\s+/)[0],digest,'digest differs');}
else throw new Error('Usage: node verify.cjs --check|--write');
console.log(`LCT-1: ${checks.length} exact checks, ${negative_controls.length} negative controls passed; sha256 ${digest}`);
