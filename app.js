// ---------------------------
// Utilidades de matemáticas
// ---------------------------
function solve2x2(a11, a12, b1, a21, a22, b2) {
  const det = a11 * a22 - a12 * a21;
  if (Math.abs(det) < 1e-9) return { ok: false };
  const x = (b1 * a22 - a12 * b2) / det;
  const y = (a11 * b2 - b1 * a21) / det;
  return { ok: true, x, y };
}

function solve3x3(A, B) {
  const M = A.map(row => row.slice());
  const Y = B.slice();
  for (let i = 0; i < 3; i++) {
    let piv = i;
    for (let r = i + 1; r < 3; r++) {
      if (Math.abs(M[r][i]) > Math.abs(M[piv][i])) piv = r;
    }
    if (Math.abs(M[piv][i]) < 1e-10) return { ok: false };
    if (piv !== i) { [M[i], M[piv]] = [M[piv], M[i]]; [Y[i], Y[piv]] = [Y[piv], Y[i]]; }
    const div = M[i][i];
    for (let c = i; c < 3; c++) M[i][c] /= div;
    Y[i] /= div;
    for (let r = 0; r < 3; r++) {
      if (r === i) continue;
      const factor = M[r][i];
      for (let c = i; c < 3; c++) M[r][c] -= factor * M[i][c];
      Y[r] -= factor * Y[i];
    }
  }
  return { ok: true, x: Y[0], y: Y[1], z: Y[2] };
}

// ---------------------------
// Parámetros y estado
// ---------------------------
const state = {
  params: {
    A0: 100, cG: 1.0, b: 5, phi: 2, theta: 0.2, tau: 0.5,
    L0: 0, lY: 0.5, li: 10,
    n0: 0, mY: 0.1, me: 2, mYf: 0.1
  },
  policy: {
    G: 100, T: 100, M: 120, P: 1, istar: 5, Yf: 150, rho: 0, k: 80
  },
  regime: { type: 'flex', ebar: 100 },
  eq: { Y: null, i: null, e: null, M_impl: null }
};

// Historial para deshacer shocks
const historyStack = [];

// ---------------------------
// Lectura de UI
// ---------------------------
function byId(id) { return document.getElementById(id); }

function setInputAndSpan(id, value) {
  const el = byId(id);
  if (!el) return;
  const min = el.min !== '' ? parseFloat(el.min) : -Infinity;
  const max = el.max !== '' ? parseFloat(el.max) : Infinity;
  let val = value;
  if (!Number.isNaN(min)) val = Math.max(min, val);
  if (!Number.isNaN(max)) val = Math.min(max, val);
  el.value = val;
  const span = byId(id + '_val');
  if (span) span.textContent = el.value;
}

function readInputs() {
  const p = state.params, pol = state.policy, reg = state.regime;

  pol.G = +byId('G').value; byId('G_val').textContent = pol.G;
  pol.T = +byId('T').value; byId('T_val').textContent = pol.T;
  pol.M = +byId('M').value; byId('M_val').textContent = pol.M;
  pol.P = +byId('P').value; byId('P_val').textContent = byId('P').value;
  pol.istar = +byId('istar').value; byId('istar_val').textContent = pol.istar;
  pol.Yf = +byId('Yf').value; byId('Yf_val').textContent = pol.Yf;
  pol.rho = +byId('rho').value; byId('rho_val').textContent = byId('rho').value;
  pol.k = +byId('k').value; byId('k_val').textContent = pol.k;

  p.A0 = +byId('A0').value;
  p.cG = +byId('cG').value;
  p.b = +byId('b').value;
  p.phi = +byId('phi').value;
  p.theta = +byId('theta').value;
  p.tau = +byId('tau').value;

  p.L0 = +byId('L0').value;
  p.lY = +byId('lY').value;
  p.li = +byId('li').value;

  p.n0 = +byId('n0').value;
  p.mY = +byId('mY').value;
  p.me = +byId('me').value;
  p.mYf = +byId('mYf').value;

  reg.type = byId('regimen').value;
  reg.ebar = +byId('e_bar').value;

  const ebox = byId('ebox');
  if (reg.type === 'fixed') ebox.classList.remove('hidden'); else ebox.classList.add('hidden');
  byId('M').disabled = (reg.type === 'fixed');
  byId('M').style.opacity = (reg.type === 'fixed') ? 0.5 : 1;
}

// ---------------------------
// Modelo
// ---------------------------
function equilibriumFlexible(p, pol) {
  const RHS1 = p.A0 + p.cG*pol.G + p.theta*pol.Yf - p.tau*pol.T;
  const RHS2 = (p.L0 - pol.M/pol.P) / p.li;
  const RHS3 = pol.k*(pol.istar + pol.rho) - p.mYf*pol.Yf - p.n0;

  const A = [
    [1, p.b, -p.phi],
    [-p.lY/p.li, 1, 0],
    [-p.mY, pol.k, p.me]
  ];
  const B = [RHS1, RHS2, RHS3];

  const sol = solve3x3(A, B);
  if (!sol.ok) return null;
  return { Y: sol.x, i: sol.y, e: sol.z, M_impl: pol.M };
}

function equilibriumFixed(p, pol, ebar) {
  const RHS1 = p.A0 + p.cG*pol.G + p.phi*ebar + p.theta*pol.Yf - p.tau*pol.T;
  const RHS2 = pol.k*(pol.istar + pol.rho) - p.me*ebar - p.mYf*pol.Yf - p.n0;

  const s = solve2x2(1, p.b, RHS1, -p.mY, pol.k, RHS2);
  if (!s.ok) return null;
  const Y = s.x, i = s.y;
  const M_impl = pol.P * (p.L0 + p.lY*Y - p.li*i);
  return { Y, i, e: ebar, M_impl };
}

function curveIS(p, pol, e) {
  return function(Y) {
    return (p.A0 + p.cG*pol.G + p.phi*e + p.theta*pol.Yf - p.tau*pol.T - Y) / p.b;
  };
}
function curveLM(p, pol, M_over_P) {
  return function(Y) {
    return (p.lY*Y + p.L0 - M_over_P) / p.li;
  };
}
function curveBP(p, pol, e) {
  return function(Y) {
    return (pol.k*(pol.istar + pol.rho) - p.me*e - p.mYf*pol.Yf - p.n0 + p.mY*Y) / pol.k;
  };
}
function fxExcessDemand(p, pol, Y, i) {
  return function(e) {
    return p.n0 - p.mY*Y + p.me*e + p.mYf*pol.Yf + pol.k*(i - pol.istar - pol.rho);
  };
}
function nxAtY(p, pol, Y) {
  return function(e) {
    return p.n0 - p.mY*Y + p.me*e + p.mYf*pol.Yf;
  };
}

// ---------------------------
// Gráficas
// ---------------------------
function makeRange(min, max, steps=200) {
  const arr = new Array(steps);
  const h = (max - min) / (steps - 1);
  for (let i = 0; i < steps; i++) arr[i] = min + i*h;
  return arr;
}

function plotAll(eq, prevEq=null) {
  const p = state.params, pol = state.policy, reg = state.regime;
  const animate = byId('animate').checked;

  const Ymin = Math.max(0, (eq.Y || 0) * 0.2);
  const Ymax = (eq.Y || 200) * 2;
  const Ys = makeRange(Ymin, Math.max(Ymin+50, Ymax));

  const e_use = eq.e;
  const is = curveIS(p, pol, e_use);
  const lm = curveLM(p, pol, (reg.type === 'fixed') ? (eq.M_impl / pol.P) : (pol.M / pol.P));
  const bp = curveBP(p, pol, e_use);

  const IS_y = Ys.map(Y => is(Y));
  const LM_y = Ys.map(Y => lm(Y));
  const BP_y = Ys.map(Y => bp(Y));

  const mainData = [
    { x: Ys, y: IS_y, type: 'scatter', mode: 'lines', name: 'IS', line: { color: '#4ea1ff', width: 2.5 } },
    { x: Ys, y: LM_y, type: 'scatter', mode: 'lines', name: 'LM', line: { color: '#49dcb1', width: 2.5 } },
    { x: Ys, y: BP_y, type: 'scatter', mode: 'lines', name: 'BP', line: { color: '#ffb347', width: 2.5, dash: 'dash' } },
    { x: [eq.Y], y: [eq.i], type: 'scatter', mode: 'markers', name: 'Equilibrio', marker: { color: '#ff6fae', size: 10 } }
  ];

  const mainLayout = {
    paper_bgcolor: 'transparent', plot_bgcolor: 'transparent',
    margin: { l: 50, r: 20, t: 10, b: 40 },
    xaxis: { title: 'Y', color: '#aab6c3', gridcolor: '#1f2937' },
    yaxis: { title: 'i (%)', color: '#aab6c3', gridcolor: '#1f2937' },
    legend: { orientation: 'h', y: 1.1, x: 0 }
  };

  const MP_min = 0, MP_max = Math.max((pol.M/pol.P)*2.2, (eq.M_impl/pol.P)*2.2);
  const MPs = makeRange(MP_min, MP_max, 200);
  const i_MD = MPs.map(MP => (p.lY*eq.Y + p.L0 - MP) / p.li);
  const moneyData = [
    { x: MPs, y: i_MD, type: 'scatter', mode: 'lines', name: 'Demanda de dinero', line: { color: '#4ea1ff' } },
    { x: [ (reg.type === 'fixed') ? (eq.M_impl/pol.P) : (pol.M/pol.P), (reg.type === 'fixed') ? (eq.M_impl/pol.P) : (pol.M/pol.P) ],
      y: [Math.min(...i_MD), Math.max(...i_MD)], type: 'scatter', mode: 'lines',
      name: 'Oferta (M/P)', line: { color: '#49dcb1' } },
    { x: [ (reg.type === 'fixed') ? (eq.M_impl/pol.P) : (pol.M/pol.P) ],
      y: [ (p.lY*eq.Y + p.L0 - ((reg.type === 'fixed') ? (eq.M_impl/pol.P) : (pol.M/pol.P))) / p.li ],
      type: 'scatter', mode: 'markers', name: 'Equilibrio dinero', marker: { color: '#ff6fae', size: 9 } }
  ];
  const moneyLayout = {
    paper_bgcolor: 'transparent', plot_bgcolor: 'transparent',
    margin: { l: 50, r: 20, t: 10, b: 40 },
    xaxis: { title: 'M/P', color: '#aab6c3', gridcolor: '#1f2937' },
    yaxis: { title: 'i (%)', color: '#aab6c3', gridcolor: '#1f2937' },
    showlegend: false
  };

  const eMin = Math.max(1, eq.e * 0.3), eMax = eq.e * 1.9;
  const Es = makeRange(eMin, Math.max(eMin+10, eMax));
  const ED = fxExcessDemand(p, pol, eq.Y, eq.i);
  const ED_y = Es.map(e => ED(e));
  const fxData = [
    { x: Es, y: ED_y, type: 'scatter', mode: 'lines', name: 'Exceso demanda', line: { color: '#ffb347' } },
    { x: [Es[0], Es[Es.length-1]], y: [0, 0], type: 'scatter', mode: 'lines', name: 'BP=0', line: { color: '#49dcb1', dash: 'dot' } },
    { x: [eq.e], y: [0], type: 'scatter', mode: 'markers', name: 'e*', marker: { color: '#ff6fae', size: 9 } }
  ];
  const fxLayout = {
    paper_bgcolor: 'transparent', plot_bgcolor: 'transparent',
    margin: { l: 50, r: 20, t: 10, b: 40 },
    xaxis: { title: 'e', color: '#aab6c3', gridcolor: '#1f2937' },
    yaxis: { title: 'Exceso de demanda de divisas', color: '#aab6c3', gridcolor: '#1f2937', zeroline: true },
    showlegend: false
  };

  const NX = nxAtY(p, pol, eq.Y);
  const NX_y = Es.map(e => NX(e));
  const nxData = [
    { x: Es, y: NX_y, type: 'scatter', mode: 'lines', name: 'NX(e)', line: { color: '#4ea1ff' } },
    { x: [eq.e, eq.e], y: [Math.min(...NX_y), Math.max(...NX_y)], type: 'scatter', mode: 'lines', name: 'e*', line: { color: '#ff6fae', dash: 'dot' } }
  ];
  const nxLayout = {
    paper_bgcolor: 'transparent', plot_bgcolor: 'transparent',
    margin: { l: 50, r: 20, t: 10, b: 40 },
    xaxis: { title: 'e', color: '#aab6c3', gridcolor: '#1f2937' },
    yaxis: { title: 'NX', color: '#aab6c3', gridcolor: '#1f2937' },
    showlegend: false
  };

  const config = { displayModeBar: false, responsive: true };
  if (prevEq && animate) {
    const steps = 10;
    const frames = [];
    for (let t = 1; t <= steps; t++) {
      const w = t / steps;
      const Yt = prevEq.Y + w * (eq.Y - prevEq.Y);
      const it = prevEq.i + w * (eq.i - prevEq.i);
      frames.push({
        data: [
          mainData[0], mainData[1], mainData[2],
          { x: [Yt], y: [it], type: 'scatter', mode: 'markers', marker: { color: '#ff6fae', size: 10 } }
        ]
      });
    }
    Plotly.react('plot_main', mainData, mainLayout, config).then(() => {
      Plotly.animate('plot_main', frames, { frame: { duration: 40, redraw: false }, transition: { duration: 40, easing: 'linear' } });
    });
  } else {
    Plotly.react('plot_main', mainData, mainLayout, config);
  }
  Plotly.react('plot_money', moneyData, moneyLayout, config);
  Plotly.react('plot_fx', fxData, fxLayout, config);
  Plotly.react('plot_nx', nxData, nxLayout, config);
}

// ---------------------------
// Ciclo principal
// ---------------------------
let prevEq = null;

function recomputeAndPlot() {
  readInputs();
  const p = state.params, pol = state.policy, reg = state.regime;

  const eq = (reg.type === 'flex')
    ? equilibriumFlexible(p, pol)
    : equilibriumFixed(p, pol, reg.ebar);

  if (!eq) { console.warn('Sin solución (revisa parámetros).'); return; }
  state.eq = eq;

  plotAll(eq, prevEq);
  prevEq = eq;
}

// vincular sliders a recomputo
function linkRangeToSpan(id) {
  const el = byId(id), span = byId(id + '_val');
  if (!el || !span) return;
  el.addEventListener('input', () => { span.textContent = el.value; recomputeAndPlot(); });
}
['G','T','M','P','istar','Yf','rho','k'].forEach(linkRangeToSpan);
['A0','cG','b','phi','theta','tau','L0','lY','li','n0','mY','me','mYf'].forEach(id => {
  byId(id).addEventListener('change', recomputeAndPlot);
});
byId('regimen').addEventListener('change', () => { recomputeAndPlot(); });
byId('e_bar').addEventListener('change', recomputeAndPlot);
byId('animate').addEventListener('change', recomputeAndPlot);

byId('reset').addEventListener('click', () => {
  document.querySelectorAll('input, select').forEach(inp => { inp.blur(); });
  byId('regimen').value = 'flex';
  byId('e_bar').value = 100;
  ['G','T','M','P','istar','Yf','rho','k'].forEach(id => {
    const defaults = {G:100,T:100,M:120,P:1, istar:5, Yf:150, rho:0, k:80};
    setInputAndSpan(id, defaults[id]);
  });
  const def = {A0:100,cG:1,b:5,phi:2,theta:0.2,tau:0.5,L0:0,lY:0.5,li:10,n0:0,mY:0.1,me:2,mYf:0.1};
  Object.keys(def).forEach(k => setInputAndSpan(k, def[k]));
  byId('animate').checked = true;
  historyStack.length = 0; // limpiar historial
  recomputeAndPlot();
});

// ---------------------------
// Shocks predefinidos + Deshacer
// ---------------------------

// Guardar snapshot anterior a un shock
function snapshotState() {
  return {
    regimen: byId('regimen').value,
    e_bar: +byId('e_bar').value,
    G: +byId('G').value,
    T: +byId('T').value,
    M: +byId('M').value,
    P: +byId('P').value,
    istar: +byId('istar').value,
    Yf: +byId('Yf').value,
    rho: +byId('rho').value,
    k: +byId('k').value
  };
}
function restoreSnapshot(snap) {
  byId('regimen').value = snap.regimen;
  setInputAndSpan('e_bar', snap.e_bar);
  setInputAndSpan('G', snap.G);
  setInputAndSpan('T', snap.T);
  setInputAndSpan('M', snap.M);
  setInputAndSpan('P', snap.P);
  setInputAndSpan('istar', snap.istar);
  setInputAndSpan('Yf', snap.Yf);
  setInputAndSpan('rho', snap.rho);
  setInputAndSpan('k', snap.k);
  recomputeAndPlot();
}

byId('undo').addEventListener('click', () => {
  if (historyStack.length === 0) return;
  const snap = historyStack.pop();
  restoreSnapshot(snap);
});

function applyAdd(id, delta) {
  const el = byId(id);
  const val = parseFloat(el.value) + delta;
  setInputAndSpan(id, val);
}
function applyMul(id, factor) {
  const el = byId(id);
  const val = parseFloat(el.value) * factor;
  setInputAndSpan(id, val);
}

// Handlers de shocks
function applyShock(kind) {
  // Guardar estado previo para poder deshacer
  historyStack.push(snapshotState());

  const regType = byId('regimen').value;
  switch (kind) {
    case 'fiscal_plus':   applyAdd('G', +20); break;
    case 'fiscal_minus':  applyAdd('G', -20); break;

    case 'monet_plus':
      // Bajo fijo, M es endógeno (no tiene efecto directo). Igual actualizamos slider si estuviera en flexible.
      if (regType === 'flex') applyAdd('M', +20);
      break;
    case 'monet_minus':
      if (regType === 'flex') applyAdd('M', -20);
      break;

    case 'istar_up':      applyAdd('istar', +1); break;
    case 'istar_down':    applyAdd('istar', -1); break;

    case 'rho_up':        applyAdd('rho', +1); break;
    case 'Yf_up':         applyAdd('Yf', +10); break;

    case 'deval':
      // Sólo relevante en fijo
      if (regType === 'fixed') applyMul('e_bar', 1.10);
      break;
    case 'reval':
      if (regType === 'fixed') applyMul('e_bar', 0.90);
      break;

    case 'k_high':        setInputAndSpan('k', 400); break;
    case 'k_low':         setInputAndSpan('k', 10); break;
  }
  recomputeAndPlot();
}

// Vincular botones
function bindShockButton(id, kind) {
  const el = byId(id);
  if (!el) return;
  el.addEventListener('click', () => applyShock(kind));
}
bindShockButton('btn_fiscal_plus', 'fiscal_plus');
bindShockButton('btn_fiscal_minus', 'fiscal_minus');
bindShockButton('btn_monet_plus',  'monet_plus');
bindShockButton('btn_monet_minus', 'monet_minus');
bindShockButton('btn_iStar_up',    'istar_up');
bindShockButton('btn_iStar_down',  'istar_down');
bindShockButton('btn_rho_up',      'rho_up');
bindShockButton('btn_Yf_up',       'Yf_up');
bindShockButton('btn_deval',       'deval');
bindShockButton('btn_reval',       'reval');
bindShockButton('btn_k_high',      'k_high');
bindShockButton('btn_k_low',       'k_low');

// Inicializar
window.addEventListener('DOMContentLoaded', () => {
  readInputs();
  recomputeAndPlot();
});