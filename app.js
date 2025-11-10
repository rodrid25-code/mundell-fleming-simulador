// ===== Utilidades numéricas =====
function solve2x2(a11, a12, b1, a21, a22, b2) {
  const det = a11*a22 - a12*a21;
  if (Math.abs(det) < 1e-9) return { ok:false };
  return { ok:true, x:(b1*a22 - a12*b2)/det, y:(a11*b2 - b1*a21)/det };
}
function solve3x3(A,B) {
  const M = A.map(r => r.slice()); const Y = B.slice();
  for (let i=0;i<3;i++) {
    let piv=i; for (let r=i+1;r<3;r++) if (Math.abs(M[r][i])>Math.abs(M[piv][i])) piv=r;
    if (Math.abs(M[piv][i])<1e-10) return { ok:false };
    if (piv!==i){ [M[i],M[piv]]=[M[piv],M[i]]; [Y[i],Y[piv]]=[Y[piv],Y[i]]; }
    const d=M[i][i]; for (let c=i;c<3;c++) M[i][c]/=d; Y[i]/=d;
    for (let r=0;r<3;r++){ if(r===i)continue; const f=M[r][i];
      for (let c=i;c<3;c++) M[r][c]-=f*M[i][c]; Y[r]-=f*Y[i]; }
  }
  return { ok:true, x:Y[0], y:Y[1], z:Y[2] };
}
function makeRange(min, max, steps=200) {
  const arr = new Array(steps); const h=(max-min)/(steps-1);
  for (let i=0;i<steps;i++) arr[i]=min+i*h; return arr;
}
function clamp(v, lo, hi){ return Math.max(lo, Math.min(hi, v)); }

// ===== Estado =====
const state = {
  params: { A0:100, cG:1.0, b:5, phi:2, theta:0.2, tau:0.5,
            L0:0, lY:0.5, li:10,
            n0:0, mY:0.1, me:2, mYf:0.1 },
  policy: { G:100, T:100, M:120, P:1, istar:5, Yf:150, rho:0, k:80 },
  regime: { type:'flex', ebar:100 },
  adas: { gammaAS:1.5, Ypot:150, Pe: null }, // Pe se setea al fijar BASE
  eq: { Y:null, i:null, e:null, M_impl:null },
  base: null // snapshot para curvas punteadas
};
const historyStack = [];

function byId(id){ return document.getElementById(id); }
function setInputAndSpan(id,value){
  const el=byId(id); if(!el) return;
  const min = el.min!==''?parseFloat(el.min):-Infinity;
  const max = el.max!==''?parseFloat(el.max):Infinity;
  const v = clamp(value, min, max); el.value=v;
  const span=byId(id+'_val'); if(span) span.textContent=String(v);
}

// ===== Lectura de UI =====
function readInputs(){
  const p=state.params, pol=state.policy, reg=state.regime, ad=state.adas;
  pol.G=+byId('G').value; byId('G_val').textContent=pol.G;
  pol.T=+byId('T').value; byId('T_val').textContent=pol.T;
  pol.M=+byId('M').value; byId('M_val').textContent=pol.M;
  pol.P=+byId('P').value; byId('P_val').textContent=byId('P').value;
  pol.istar=+byId('istar').value; byId('istar_val').textContent=pol.istar;
  pol.Yf=+byId('Yf').value; byId('Yf_val').textContent=pol.Yf;
  pol.rho=+byId('rho').value; byId('rho_val').textContent=byId('rho').value;
  pol.k=+byId('k').value; byId('k_val').textContent=pol.k;

  p.A0=+byId('A0').value; p.cG=+byId('cG').value; p.b=+byId('b').value;
  p.phi=+byId('phi').value; p.theta=+byId('theta').value; p.tau=+byId('tau').value;
  p.L0=+byId('L0').value; p.lY=+byId('lY').value; p.li=+byId('li').value;
  p.n0=+byId('n0').value; p.mY=+byId('mY').value; p.me=+byId('me').value; p.mYf=+byId('mYf').value;

  reg.type=byId('regimen').value; reg.ebar=+byId('e_bar').value;

  ad.gammaAS=+byId('gammaAS').value; byId('gammaAS_val').textContent=ad.gammaAS;
  ad.Ypot=+byId('Ypot').value; byId('Ypot_val').textContent=ad.Ypot;

  const ebox=byId('ebox'); if(reg.type==='fixed') ebox.classList.remove('hidden'); else ebox.classList.add('hidden');
  byId('M').disabled = (reg.type==='fixed'); byId('M').style.opacity=(reg.type==='fixed')?0.5:1;
}

// ===== Modelo (equilibrios) =====
function equilibriumFlexible(p, pol){
  const RHS1=p.A0 + p.cG*pol.G + p.theta*pol.Yf - p.tau*pol.T;
  const RHS2=(p.L0 - pol.M/pol.P)/p.li;
  const RHS3=pol.k*(pol.istar+pol.rho) - p.mYf*pol.Yf - p.n0;
  const A=[[1,p.b,-p.phi],[-p.lY/p.li,1,0],[-p.mY,pol.k,p.me]];
  const B=[RHS1,RHS2,RHS3];
  const s=solve3x3(A,B); if(!s.ok) return null;
  return { Y:s.x, i:s.y, e:s.z, M_impl:pol.M };
}
function equilibriumFixed(p, pol, ebar){
  const RHS1=p.A0 + p.cG*pol.G + p.phi*ebar + p.theta*pol.Yf - p.tau*pol.T;
  const RHS2=pol.k*(pol.istar+pol.rho) - p.me*ebar - p.mYf*pol.Yf - p.n0;
  const s=solve2x2(1,p.b,RHS1,-p.mY,pol.k,RHS2); if(!s.ok) return null;
  const Y=s.x, i=s.y; const M_impl=pol.P*(p.L0 + p.lY*Y - p.li*i);
  return { Y,i,e:ebar,M_impl };
}
function getEquilibrium(){
  const {params:p, policy:pol, regime:reg}=state;
  return reg.type==='flex' ? equilibriumFlexible(p,pol) : equilibriumFixed(p,pol,reg.ebar);
}

// ===== Curvas para 2×2 =====
// IS, LM en (Y,i)
function curveIS(p,pol,e){ return Y => (p.A0+p.cG*pol.G+p.phi*e+p.theta*pol.Yf-p.tau*pol.T - Y)/p.b; }
function curveLM(p,pol,MP){ return Y => (p.lY*Y + p.L0 - MP)/p.li; }
// BP en (Y,i)
function curveBP(p,pol,e){ return Y => (pol.k*(pol.istar+pol.rho) - p.me*e - p.mYf*pol.Yf - p.n0 + p.mY*Y)/pol.k; }

// Divisas S/D:
// Elegimos una demanda base D_fx(e)=q0 + qY*Y - qE*e (qE>0). Construimos S_fx(e)=D_fx(e)+NX(e)
function buildFX_SD(p,pol,Y){
  const qE = Math.max(0.5, 0.5*p.me);  // pendiente negativa "suave"
  const qY = 0.2;                      // demanda de FC aumenta con Y
  const q0 = 50;                       // nivel base (ajustable)

  const NX = e => (p.n0 - p.mY*Y + p.me*e + p.mYf*pol.Yf);
  const D = e => q0 + qY*Y - qE*e;
  const S = e => D(e) + NX(e); // asegura S-D=NX(e)
  return { D, S };
}

// AD–AS: AD(P) por recomputo del equilibrio para cada P; AS(P)=Y_pot + gamma*(P - Pe)
function computeADcurve(Pgrid){
  const outY=[]; const snapshot = JSON.stringify({state});
  for(const P of Pgrid){
    const savedP = state.policy.P; state.policy.P = P;
    const eq = getEquilibrium(); outY.push(eq ? eq.Y : NaN);
    state.policy.P = savedP; // restore rápido
  }
  return outY;
}
function computeAScurve(Pgrid){
  const { Ypot, gammaAS, Pe } = state.adas;
  return Pgrid.map(P => Ypot + gammaAS * ((Pe ?? P) - (Pe ?? P) + (P - (Pe ?? P)))); // simplifica a Ypot + gamma*(P-Pe)
}

// ===== Plotting =====
function baseColor(line){ return {color:'#aab6c3', width:2, dash:'dot'}; }
function solid(color, width=2.5){ return {color, width}; }

function plotAll(prevEq=null){
  const p=state.params, pol=state.policy, reg=state.regime;
  const animate=byId('animate').checked;
  const eq = state.eq;

  // ---------- Panel 4: IS–LM ----------
  const Ymin = Math.max(0, (eq.Y||100)*0.2), Ymax=(eq.Y||150)*2;
  const Ys = makeRange(Ymin, Math.max(Ymin+50,Ymax));
  const is = curveIS(p,pol,eq.e);
  const lm = curveLM(p,pol, (reg.type==='fixed')? (eq.M_impl/pol.P) : (pol.M/pol.P));
  const dataISLM = [
    { x: Ys, y: Ys.map(is), type:'scatter', mode:'lines', name:'IS', line: solid('#4ea1ff') },
    { x: Ys, y: Ys.map(lm), type:'scatter', mode:'lines', name:'LM', line: solid('#49dcb1') },
    { x:[eq.Y], y:[eq.i], type:'scatter', mode:'markers', name:'Equilibrio', marker:{color:'#ff6fae', size:10} }
  ];
  if (state.base){
    const polB = state.base.policy; const regB = state.base.regime; const pB = state.base.params;
    const eqB = state.base.eq;
    const isB = Y => (pB.A0+pB.cG*polB.G+pB.phi*eqB.e+pB.theta*polB.Yf-pB.tau*polB.T - Y)/pB.b;
    const lmB = Y => (pB.lY*Y + pB.L0 - ((regB.type==='fixed') ? (state.base.eq.M_impl / polB.P) : (polB.M/polB.P)))/pB.li;
    dataISLM.unshift(
      { x: Ys, y: Ys.map(isB), type:'scatter', mode:'lines', name:'IS (base)', line: baseColor() },
      { x: Ys, y: Ys.map(lmB), type:'scatter', mode:'lines', name:'LM (base)', line: baseColor() },
    );
    dataISLM.push({ x:[eqB.Y], y:[eqB.i], type:'scatter', mode:'markers', name:'Eq. base', marker:{color:'#aab6c3', size:8, symbol:'x'} });
  }
  const layoutISLM = {
    paper_bgcolor:'transparent', plot_bgcolor:'transparent',
    margin:{l:50,r:20,t:10,b:40}, xaxis:{title:'Y', color:'#aab6c3', gridcolor:'#1f2937'},
    yaxis:{title:'i (%)', color:'#aab6c3', gridcolor:'#1f2937'},
    legend:{orientation:'h', y:1.13, x:0}
  };

  // ---------- Panel 2: BP ----------
  const bp = curveBP(p,pol,eq.e);
  const dataBP = [
    { x: Ys, y: Ys.map(bp), type:'scatter', mode:'lines', name:'BP', line: solid('#ffb347') },
    { x:[eq.Y], y:[eq.i], type:'scatter', mode:'markers', name:'(Y,i)', marker:{color:'#ff6fae', size:9} }
  ];
  if (state.base){
    const pB=state.base.params, polB=state.base.policy, eqB=state.base.eq;
    const bpB = Y => (polB.k*(polB.istar+polB.rho) - pB.me*eqB.e - pB.mYf*polB.Yf - pB.n0 + pB.mY*Y)/polB.k;
    dataBP.unshift({ x: Ys, y: Ys.map(bpB), type:'scatter', mode:'lines', name:'BP (base)', line: baseColor() });
    dataBP.push({ x:[eqB.Y], y:[eqB.i], type:'scatter', mode:'markers', name:'Punto base', marker:{color:'#aab6c3', size:8, symbol:'x'} });
  }
  const layoutBP = {
    paper_bgcolor:'transparent', plot_bgcolor:'transparent',
    margin:{l:50,r:20,t:10,b:40}, xaxis:{title:'Y', color:'#aab6c3', gridcolor:'#1f2937'},
    yaxis:{title:'i (%)', color:'#aab6c3', gridcolor:'#1f2937'}, legend:{orientation:'h', y:1.13, x:0}
  };

  // ---------- Panel 1: Divisas (S/D) ----------
  const fxSD = buildFX_SD(p,pol,eq.Y);
  const eMin = Math.max(1, eq.e*0.3), eMax = Math.max(eMin+10, eq.e*1.9);
  const Es = makeRange(eMin, eMax);
  const D = Es.map(e => fxSD.D(e));
  const S = Es.map(e => fxSD.S(e));
  const Qeq = fxSD.D(eq.e);
  const dataFX = [
    { x: Es, y: D, type:'scatter', mode:'lines', name:'Demanda FC', line: solid('#4ea1ff') },
    { x: Es, y: S, type:'scatter', mode:'lines', name:'Oferta FC', line: solid('#49dcb1') },
    { x: [eq.e], y:[Qeq], type:'scatter', mode:'markers', name:'Eq FX', marker:{color:'#ff6fae', size:9} }
  ];
  if (state.base){
    const fxB = buildFX_SD(state.base.params, state.base.policy, state.base.eq.Y);
    const DB = Es.map(e => fxB.D(e)); const SB = Es.map(e => fxB.S(e));
    const Qb = fxB.D(state.base.eq.e);
    dataFX.unshift(
      { x: Es, y: DB, type:'scatter', mode:'lines', name:'Dem FC (base)', line: baseColor() },
      { x: Es, y: SB, type:'scatter', mode:'lines', name:'Ofer FC (base)', line: baseColor() }
    );
    dataFX.push({ x:[state.base.eq.e], y:[Qb], type:'scatter', mode:'markers', name:'FX base', marker:{color:'#aab6c3', size:8, symbol:'x'} });
  }
  const layoutFX = {
    paper_bgcolor:'transparent', plot_bgcolor:'transparent',
    margin:{l:50,r:20,t:10,b:40}, xaxis:{title:'e', color:'#aab6c3', gridcolor:'#1f2937'},
    yaxis:{title:'Q_fx', color:'#aab6c3', gridcolor:'#1f2937'}, legend:{orientation:'h', y:1.13, x:0}
  };

  // ---------- Panel 3: AD–AS ----------
  const Pmin = Math.max(0.2, pol.P*0.5), Pmax = pol.P*1.8;
  const Ps = makeRange(Pmin, Pmax);
  const ADy = computeADcurve(Ps);
  const ASy = computeAScurve(Ps);
  const dataADAS = [
    { x: Ps, y: ADy, type:'scatter', mode:'lines', name:'AD', line: solid('#4ea1ff') },
    { x: Ps, y: ASy, type:'scatter', mode:'lines', name:'AS', line: solid('#49dcb1') },
    { x:[pol.P], y:[eq.Y], type:'scatter', mode:'markers', name:'(P,Y)', marker:{color:'#ff6fae', size:9} }
  ];
  if (state.base){
    // AD base: recalculada con parámetros base para cada P
    const polSaved = {...state.policy}; const parSaved = {...state.params}; const regSaved = {...state.regime};
    // swap a base
    state.policy = {...state.base.policy};
    state.params = {...state.base.params};
    state.regime = {...state.base.regime};
    const ADyB = computeADcurve(Ps);
    // restore
    state.policy = polSaved; state.params = parSaved; state.regime = regSaved;

    // AS base: usa Ypot_base=Y_base y Pe=P_base (ya guardado en state.base.adas.Pe)
    const YpotB = state.base.eq.Y; const PeB = state.base.adas.Pe;
    const gammaB = state.base.adas.gammaAS;
    const ASyB = Ps.map(P => YpotB + gammaB*(P - PeB));

    dataADAS.unshift(
      { x: Ps, y: ADyB, type:'scatter', mode:'lines', name:'AD (base)', line: baseColor() },
      { x: Ps, y: ASyB, type:'scatter', mode:'lines', name:'AS (base)', line: baseColor() }
    );
    dataADAS.push({ x:[state.base.policy.P], y:[state.base.eq.Y], type:'scatter', mode:'markers', name:'Punto base', marker:{color:'#aab6c3', size:8, symbol:'x'} });
  }
  const layoutADAS = {
    paper_bgcolor:'transparent', plot_bgcolor:'transparent',
    margin:{l:50,r:20,t:10,b:40}, xaxis:{title:'P', color:'#aab6c3', gridcolor:'#1f2937'},
    yaxis:{title:'Y', color:'#aab6c3', gridcolor:'#1f2937'}, legend:{orientation:'h', y:1.13, x:0}
  };

  const config = { displayModeBar:false, responsive:true };

  // Render (animación sólo en IS–LM para el punto)
  if (prevEq && animate){
    const steps=10, frames=[];
    for (let t=1;t<=steps;t++){
      const w=t/steps, Yt=prevEq.Y+w*(eq.Y-prevEq.Y), it=prevEq.i+w*(eq.i-prevEq.i);
      frames.push({ data:[
        dataISLM[0], dataISLM[1], ...(dataISLM.length>3?[dataISLM[2]]:[]),
        { x:[Yt], y:[it], type:'scatter', mode:'markers', marker:{color:'#ff6fae', size:10} }
      ]});
    }
    Plotly.react('plot_islm', dataISLM, layoutISLM, config).then(()=> {
      Plotly.animate('plot_islm', frames, { frame:{duration:40, redraw:false}, transition:{duration:40, easing:'linear'} });
    });
  } else {
    Plotly.react('plot_islm', dataISLM, layoutISLM, config);
  }
  Plotly.react('plot_bp', dataBP, layoutBP, config);
  Plotly.react('plot_fx_sd', dataFX, layoutFX, config);
  Plotly.react('plot_adas', dataADAS, layoutADAS, config);
}

// ===== Ciclo principal =====
let prevEq=null;
function recomputeAndPlot(){
  readInputs();
  const eq=getEquilibrium(); if(!eq){ console.warn('Sin solución. Revisa parámetros.'); return; }
  state.eq=eq;
  plotAll(prevEq); prevEq=eq;
}

// ===== Vínculos UI =====
function linkRange(id){ const el=byId(id), sp=byId(id+'_val'); if(!el) return;
  el.addEventListener('input', ()=>{ if(sp) sp.textContent=el.value; recomputeAndPlot(); });
}
['G','T','M','P','istar','Yf','rho','k','gammaAS','Ypot'].forEach(linkRange);
['A0','cG','b','phi','theta','tau','L0','lY','li','n0','mY','me','mYf'].forEach(id => {
  byId(id).addEventListener('change', recomputeAndPlot);
});
byId('regimen').addEventListener('change', recomputeAndPlot);
byId('e_bar').addEventListener('change', recomputeAndPlot);
byId('animate').addEventListener('change', recomputeAndPlot);

// Reset
byId('reset').addEventListener('click', ()=>{
  const defaults = {regimen:'flex', e_bar:100, G:100,T:100,M:120,P:1,istar:5,Yf:150,rho:0,k:80, gammaAS:1.5, Ypot:150};
  byId('regimen').value=defaults.regimen; byId('e_bar').value=defaults.e_bar;
  Object.entries({G:100,T:100,M:120,P:1,istar:5,Yf:150,rho:0,k:80,gammaAS:1.5,Ypot:150}).forEach(([k,v])=> setInputAndSpan(k,v));
  const defP={A0:100,cG:1,b:5,phi:2,theta:0.2,tau:0.5,L0:0,lY:0.5,li:10,n0:0,mY:0.1,me:2,mYf:0.1};
  Object.entries(defP).forEach(([k,v])=> setInputAndSpan(k,v));
  state.adas.Pe=null; state.base=null; historyStack.length=0;
  byId('animate').checked=true;
  recomputeAndPlot();
});

// Fijar/Limpiar base
byId('fix_base').addEventListener('click', ()=>{
  // snapshot completo
  state.base = {
    params: {...state.params},
    policy: {...state.policy},
    regime: {...state.regime},
    adas: { gammaAS: state.adas.gammaAS, Pe: state.policy.P },
    eq: {...state.eq}
  };
  // fijar expectativas a P_base
  state.adas.Pe = state.policy.P;
  recomputeAndPlot();
});
byId('clear_base').addEventListener('click', ()=>{
  state.base=null; state.adas.Pe=null; recomputeAndPlot();
});

// Deshacer (último shock)
byId('undo').addEventListener('click', ()=>{
  if (historyStack.length===0) return;
  const s=historyStack.pop();
  byId('regimen').value=s.regimen; setInputAndSpan('e_bar', s.e_bar);
  ['G','T','M','P','istar','Yf','rho','k','gammaAS','Ypot'].forEach(k=> setInputAndSpan(k, s[k]));
  recomputeAndPlot();
});
function snapshotUI(){
  return { regimen:byId('regimen').value, e_bar:+byId('e_bar').value,
    G:+byId('G').value, T:+byId('T').value, M:+byId('M').value, P:+byId('P').value,
    istar:+byId('istar').value, Yf:+byId('Yf').value, rho:+byId('rho').value, k:+byId('k').value,
    gammaAS:+byId('gammaAS').value, Ypot:+byId('Ypot').value };
}

// Shocks
function applyAdd(id,delta){ const el=byId(id); setInputAndSpan(id, +el.value + delta); }
function applyMul(id,f){ const el=byId(id); setInputAndSpan(id, +el.value * f); }
function applyShock(kind){
  historyStack.push(snapshotUI());
  const regType=byId('regimen').value;
  switch(kind){
    case 'fiscal_plus': applyAdd('G',+20); break;
    case 'fiscal_minus': applyAdd('G',-20); break;
    case 'monet_plus': if(regType==='flex') applyAdd('M',+20); break;
    case 'monet_minus': if(regType==='flex') applyAdd('M',-20); break;
    case 'istar_up': applyAdd('istar',+1); break;
    case 'istar_down': applyAdd('istar',-1); break;
    case 'rho_up': applyAdd('rho',+1); break;
    case 'Yf_up': applyAdd('Yf',+10); break;
    case 'deval': if(regType==='fixed') applyMul('e_bar',1.10); break;
    case 'reval': if(regType==='fixed') applyMul('e_bar',0.90); break;
    case 'k_high': setInputAndSpan('k',400); break;
    case 'k_low': setInputAndSpan('k',10); break;
  }
  recomputeAndPlot();
}
function bindShock(id, kind){ const el=byId(id); if(el) el.addEventListener('click', ()=>applyShock(kind)); }
bindShock('btn_fiscal_plus','fiscal_plus');
bindShock('btn_fiscal_minus','fiscal_minus');
bindShock('btn_monet_plus','monet_plus');
bindShock('btn_monet_minus','monet_minus');
bindShock('btn_iStar_up','istar_up');
bindShock('btn_iStar_down','istar_down');
bindShock('btn_rho_up','rho_up');
bindShock('btn_Yf_up','Yf_up');
bindShock('btn_deval','deval');
bindShock('btn_reval','reval');
bindShock('btn_k_high','k_high');
bindShock('btn_k_low','k_low');

// Init
window.addEventListener('DOMContentLoaded', ()=>{ readInputs(); recomputeAndPlot(); });
