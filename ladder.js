/* Client-side XML → level-ladder (same physics as xmlToLevelLadder.py). */
(function (global) {
  const ENERGY_TOL = 1.0;
  const TINY = -1.17549435e-38;
  const ALPHA = 1.0 / 137.0;
  const PI = 3.141592653589;
  const ELECTRON_MASS = 511.0;
  const GAMMA_FINAL = 0.9999999;

  function num(v, fallback) {
    if (v == null || v === "") return fallback;
    const n = parseFloat(v);
    return Number.isFinite(n) ? n : fallback;
  }

  class FermiDistribution {
    constructor(atomicNumber, qVal, eCharge) {
      if (Math.abs(eCharge) !== 1) throw new Error("eCharge must be ±1");
      if (qVal <= 0) throw new Error("qVal must be positive");
      this.z = atomicNumber;
      this.qVal = qVal;
      this.eCharge = eCharge;
      this.v0 = 1.13 * ALPHA ** 2 * (this.z - this.eCharge) ** (4 / 3);
      this.s = (1.0 - ALPHA ** 2 * this.z * this.z) ** 0.5 - 1.0;
      this.xMin = eCharge > 0 ? 0.0 : (ELECTRON_MASS * this.v0) / this.qVal;
      this.xMax = 1.0;
    }

    eulerGamma(delta0) {
      let gammaValue = 1.0;
      let i = 0.0;
      let denom = (i + 1.0 + this.s) ** 2;
      let step = 1.0 / (1.0 + (delta0 * delta0) / denom);
      while (step <= GAMMA_FINAL) {
        i += 1.0;
        gammaValue *= step;
        denom = (i + 1.0 + this.s) ** 2;
        step = 1.0 / (1.0 + (delta0 * delta0) / denom);
        if (i > 1_000_000) break;
      }
      return gammaValue * step;
    }

    densityX(x) {
      if (x < this.xMin || x > this.xMax) return 0.0;
      const gamma0 = 1.0 + (this.qVal * x) / ELECTRON_MASS - this.eCharge * this.v0;
      if (gamma0 < 1.0) return 0.0;
      const ni0 = Math.sqrt(gamma0 * gamma0 - 1.0);
      if (ni0 <= 0.0) return 0.0;
      const delta0 = (ALPHA * this.z * gamma0) / ni0;
      return (
        gamma0 *
        ni0 *
        (1.0 - x) ** 2 *
        ni0 ** (2.0 * this.s) *
        Math.exp(-this.eCharge * PI * delta0) *
        this.eulerGamma(delta0)
      );
    }

    spectrumKeV(kinetic) {
      return kinetic.map((k) => this.densityX(k / this.qVal));
    }
  }

  function parseNuclide(text, name) {
    const doc = new DOMParser().parseFromString(text, "text/xml");
    const err = doc.querySelector("parsererror");
    if (err) throw new Error("Could not parse XML: " + err.textContent.slice(0, 180));
    const root = doc.documentElement;
    if (!root || root.tagName !== "Nuclide") {
      throw new Error((name || "file") + ": expected <Nuclide> root");
    }
    const z = Math.trunc(num(root.getAttribute("AtomicNumber"), 0));
    const a = Math.trunc(num(root.getAttribute("AtomicMass"), 0));
    const qBeta = num(root.getAttribute("QBeta"), 0) || 0;
    const levels = {};
    const betas = [];
    const gammas = {};
    for (const lev of root.getElementsByTagName("Level")) {
      const e = num(lev.getAttribute("Energy"), 0) || 0;
      levels[e] = true;
      for (const tr of lev.getElementsByTagName("Transition")) {
        const typ = (tr.getAttribute("Type") || "").trim();
        if (typ === "Fake" || typ === "") continue;
        const intensity = num(tr.getAttribute("Intensity"), 0) || 0;
        const qVal = num(tr.getAttribute("TransitionQValue"), 0) || 0;
        const tgt = tr.getElementsByTagName("TargetLevel")[0];
        const eTgt = tgt ? num(tgt.getAttribute("Energy"), 0) || 0 : 0;
        const zTgt = tgt ? Math.trunc(num(tgt.getAttribute("AtomicNumber"), 0)) : 0;
        const aTgt = tgt ? Math.trunc(num(tgt.getAttribute("AtomicMass"), 0)) : 0;
        const rec = { from_e: e, to_e: eTgt, q: qVal, intensity, z_target: zTgt, a_target: aTgt };
        if (typ.toUpperCase().startsWith("B")) {
          rec.charge = typ.includes("-") ? -1 : 1;
          betas.push(rec);
        } else if (typ.toUpperCase().startsWith("G")) {
          if (!gammas[e]) gammas[e] = [];
          gammas[e].push(rec);
        }
      }
    }
    return {
      name: name || "upload.xml",
      z, a, q_beta: qBeta,
      levels: Object.keys(levels).map(Number).sort((x, y) => x - y),
      betas, gammas, text,
    };
  }

  function snapEnergy(e, canonical, tol) {
    if (!canonical.length) return e;
    let nearest = canonical[0];
    let best = Math.abs(nearest - e);
    for (const c of canonical) {
      const d = Math.abs(c - e);
      if (d < best) { nearest = c; best = d; }
    }
    return best <= tol ? nearest : e;
  }

  function snapGraph(nuclides) {
    const canonical = [...new Set(nuclides.flatMap((n) => n.levels))].sort((a, b) => a - b);
    const betas = [];
    const gammas = {};
    let qBeta = Math.max(0, ...nuclides.map((n) => n.q_beta));
    let maxQ = 0;
    let zFermi = null;
    let eCharge = -1;
    for (const n of nuclides) {
      for (const b of n.betas) {
        if (b.intensity <= 0) continue;
        const rec = { ...b };
        rec.to_e = snapEnergy(b.to_e, canonical, ENERGY_TOL);
        rec.from_e = snapEnergy(b.from_e, canonical, ENERGY_TOL);
        betas.push(rec);
        if (rec.z_target) zFermi = rec.z_target;
        eCharge = rec.charge;
        maxQ = Math.max(maxQ, rec.q);
        const implied = rec.q + rec.to_e;
        if (implied > qBeta) qBeta = implied;
      }
      for (const [src, lst] of Object.entries(n.gammas)) {
        const srcS = snapEnergy(Number(src), canonical, ENERGY_TOL);
        for (const g of lst) {
          if (g.intensity <= 0) continue;
          const rec = { ...g };
          rec.from_e = srcS;
          rec.to_e = snapEnergy(g.to_e, canonical, ENERGY_TOL);
          if (!gammas[srcS]) gammas[srcS] = [];
          gammas[srcS].push(rec);
        }
      }
    }
    if (zFermi == null) {
      const parents = nuclides.filter((n) => n.betas.length);
      zFermi = parents.length
        ? parents[0].z + (eCharge < 0 ? 1 : -1)
        : nuclides[0].z;
    }
    return { canonical, betas, gammas, qBeta, maxQ, zFermi: Math.trunc(zFermi), eCharge };
  }

  function normalizeBranchings(betas, gammas) {
    const totB = betas.reduce((s, b) => s + b.intensity, 0);
    if (totB > 0) {
      for (const b of betas) b.intensity = (b.intensity / totB) * 100;
    }
    const outG = {};
    for (const [src, lst] of Object.entries(gammas)) {
      const tot = lst.reduce((s, g) => s + g.intensity, 0);
      outG[src] = lst.map((g) => ({
        ...g,
        intensity: tot > 0 ? (g.intensity / tot) * 100 : g.intensity,
      }));
    }
    return { betas, gammas: outG };
  }

  function lookupGamma(gammas, frm, to) {
    const lst = gammas[frm] || gammas[String(frm)] || [];
    for (const g of lst) {
      if (Math.abs(g.to_e - to) <= ENERGY_TOL) return g.intensity;
    }
    return 100;
  }

  function gammaChains(start, gammas, gs, acc, trail) {
    const key = Math.round(start * 1000) / 1000;
    const seen = trail || [];
    if (seen.indexOf(key) >= 0) return;
    if (Math.abs(start - gs) <= ENERGY_TOL) {
      acc.push([gs]);
      return;
    }
    const branches = gammas[start] || gammas[String(start)] || [];
    if (!branches.length) {
      acc.push([start, gs]);
      return;
    }
    const nxt = seen.concat([key]);
    for (const g of branches) {
      const rest = [];
      gammaChains(g.to_e, gammas, gs, rest, nxt);
      for (const r of rest) acc.push([start, ...r]);
    }
  }

  function enumeratePaths(betas, gammas) {
    let feeds = betas;
    if (!feeds.length) {
      const starts = Object.keys(gammas).map(Number).sort((a, b) => a - b);
      feeds = starts.map((e) => ({ to_e: e, q: null, intensity: 1 }));
      if (!starts.some((e) => Math.abs(e) <= ENERGY_TOL)) {
        feeds.unshift({ to_e: 0, q: null, intensity: 1 });
      }
    }
    const paths = [];
    const seen = new Set();
    for (const b of feeds) {
      const chains = [];
      gammaChains(b.to_e, gammas, 0, chains);
      for (const chain of chains) {
        const key = chain.map((x) => x.toFixed(3)).join(",");
        if (seen.has(key)) continue;
        seen.add(key);
        paths.push({
          feed: b.to_e,
          chain,
          q_beta_branch: b.q,
          beta_intensity: b.intensity,
          gamma_intensities: chain.slice(0, -1).map((frm, i) => lookupGamma(gammas, frm, chain[i + 1])),
        });
      }
    }
    paths.sort((p, q) => p.feed - q.feed || p.chain.length - q.chain.length || p.chain[0] - q.chain[0]);
    return paths;
  }

  function applyPathWeights(paths) {
    const weights = paths.map((p) => {
      let w = p.beta_intensity / 100;
      for (const ig of p.gamma_intensities) w *= ig / 100;
      return w;
    });
    const tot = weights.reduce((s, w) => s + w, 0);
    const normed = tot > 0 ? weights.map((w) => w / tot) : weights;
    paths.forEach((p, i) => { p.weight = normed[i]; });
    return paths;
  }

  function fermiGrid(z, q, eCharge, cache) {
    const key = z + "|" + Math.round(q * 1000) / 1000 + "|" + eCharge;
    if (!cache[key]) {
      const fermi = new FermiDistribution(z, q, eCharge);
      const n = Math.max(Math.ceil(q + 2), 1);
      const grid = new Float64Array(n);
      for (let i = 0; i < n; i++) grid[i] = fermi.densityX((i + 0.5) / q);
      cache[key] = grid;
    }
    return cache[key];
  }

  function fillFermi(row, eX, qVal, eMax, z, eCharge, weight, cache) {
    const nBins = row.length;
    const iX = Math.max(0, Math.min(nBins, Math.round(eX)));
    const amp = -weight;
    if (iX >= nBins) return;
    let q = qVal;
    if (q == null || q <= 0) q = Math.max(eMax - eX, 1);
    const grid = fermiGrid(z, q, eCharge, cache);
    let peak = 0;
    const spec = new Float64Array(nBins - iX);
    for (let i = iX; i < nBins; i++) {
      const k = i + 0.5 - eX;
      let v = 0;
      if (k >= 0 && k <= q + 0.5) {
        const idx = Math.max(0, Math.min(grid.length - 1, Math.floor(k)));
        v = grid[idx];
      }
      spec[i - iX] = v;
      if (v > peak) peak = v;
    }
    if (peak > 0) {
      for (let i = iX; i < nBins; i++) {
        const v = amp * spec[i - iX] / peak;
        row[i] = v < TINY ? v : TINY;
      }
    } else {
      for (let i = iX; i < nBins; i++) row[i] = amp < 0 ? amp : TINY;
    }
  }

  function buildLadder(paths, eMax, z, eCharge, qBeta) {
    const nBins = Math.round(eMax);
    const arr = new Float32Array(paths.length * nBins);
    const cache = {};
    for (let i = 0; i < paths.length; i++) {
      const p = paths[i];
      const w = p.weight;
      const row = arr.subarray(i * nBins, (i + 1) * nBins);
      const chain = p.chain;
      const eX = chain.length ? chain[0] : 0;
      const ints = p.gamma_intensities || [];
      for (let k = 0; k < chain.length - 1; k++) {
        const frm = chain[k];
        const to = chain[k + 1];
        const ig = k < ints.length ? ints[k] : 100;
        const lo = Math.min(frm, to);
        const hi = Math.max(frm, to);
        const i0 = Math.max(0, Math.min(nBins, Math.round(lo)));
        const i1 = Math.max(0, Math.min(nBins, Math.round(hi)));
        const z = w * (ig / 100);
        for (let j = i0; j < i1; j++) row[j] = z;
      }
      let qVal = p.q_beta_branch;
      if (qVal == null) qVal = Math.max(qBeta - eX, 1);
      fillFermi(row, eX, qVal, eMax, z, eCharge, w, cache);
    }
    return { arr, rows: paths.length, cols: nBins };
  }

  function runTranslation(nuclides, eMaxOpt) {
    const snapped = snapGraph(nuclides);
    const normed = normalizeBranchings(snapped.betas, snapped.gammas);
    const paths = applyPathWeights(enumeratePaths(normed.betas, normed.gammas));
    if (!paths.length) throw new Error("no decay paths found");
    let eMax = eMaxOpt;
    if (eMax == null) eMax = snapped.maxQ > 0 ? snapped.maxQ : snapped.qBeta;
    eMax = Math.max(Math.ceil(eMax), 1);
    const ladder = buildLadder(paths, eMax, snapped.zFermi, snapped.eCharge, snapped.qBeta);
    let mn = Infinity;
    let mx = -Infinity;
    for (let i = 0; i < ladder.arr.length; i++) {
      const v = ladder.arr[i];
      if (v < mn) mn = v;
      if (v > mx) mx = v;
    }
    return {
      ...ladder,
      paths: paths.map((p, i) => ({
        index: i + 1,
        weight: p.weight,
        feed: p.feed,
        chain: p.chain,
        beta_intensity: p.beta_intensity,
        gamma_intensities: p.gamma_intensities,
        q_beta_branch: p.q_beta_branch,
      })),
      e_max: eMax,
      q_beta: snapped.qBeta,
      max_q_transition: snapped.maxQ,
      z_fermi: snapped.zFermi,
      e_charge: snapped.eCharge,
      levels: snapped.canonical,
      min: mn,
      max: mx,
      shape: [ladder.rows, ladder.cols],
    };
  }

  function lerpColor(t) {
    const stops = [
      [0.0, [158, 1, 66]],
      [0.1, [213, 62, 79]],
      [0.2, [244, 109, 67]],
      [0.3, [253, 174, 97]],
      [0.4, [254, 224, 139]],
      [0.5, [255, 255, 191]],
      [0.6, [230, 245, 152]],
      [0.7, [171, 221, 164]],
      [0.8, [102, 194, 165]],
      [0.9, [50, 136, 189]],
      [1.0, [94, 79, 162]],
    ];
    const x = Math.max(0, Math.min(1, t));
    let i = 0;
    while (i < stops.length - 2 && x > stops[i + 1][0]) i += 1;
    const [t0, c0] = stops[i];
    const [t1, c1] = stops[i + 1];
    const u = (x - t0) / (t1 - t0 || 1);
    return [
      Math.round(c0[0] + (c1[0] - c0[0]) * u),
      Math.round(c0[1] + (c1[1] - c0[1]) * u),
      Math.round(c0[2] + (c1[2] - c0[2]) * u),
    ];
  }

  function colorFor(v, lim) {
    if (!lim) return [255, 255, 191];
    return lerpColor(0.5 + 0.5 * (v / lim));
  }

  function clamp(v, lo, hi) {
    return Math.max(lo, Math.min(hi, v));
  }

  function normalizeView(result, view) {
    const rows = result.rows;
    const cols = result.cols;
    const v = view || {};
    let p0 = Math.floor(v.p0 != null ? v.p0 : 0);
    let p1 = Math.ceil(v.p1 != null ? v.p1 : rows);
    let e0 = v.e0 != null ? v.e0 : 0;
    let e1 = v.e1 != null ? v.e1 : cols;
    if (p1 < p0) { const t = p0; p0 = p1; p1 = t; }
    if (e1 < e0) { const t = e0; e0 = e1; e1 = t; }
    p0 = clamp(p0, 0, Math.max(0, rows - 1));
    p1 = clamp(p1, p0 + 1, rows);
    e0 = clamp(e0, 0, cols);
    e1 = clamp(e1, 0, cols);
    if (e1 - e0 < 1) e1 = Math.min(cols, e0 + 1);
    if (e1 <= e0) { e0 = Math.max(0, e1 - 1); }
    return { p0, p1, e0, e1 };
  }

  function limForView(result, view) {
    const { arr, cols } = result;
    let mn = Infinity;
    let mx = -Infinity;
    const i0 = Math.max(0, Math.floor(view.e0));
    const i1 = Math.min(cols, Math.ceil(view.e1));
    for (let p = view.p0; p < view.p1; p++) {
      const off = p * cols;
      for (let ie = i0; ie < i1; ie++) {
        const val = arr[off + ie];
        if (val < mn) mn = val;
        if (val > mx) mx = val;
      }
    }
    if (!Number.isFinite(mn)) return 1e-12;
    return Math.max(Math.abs(mn), Math.abs(mx), 1e-12);
  }

  function niceTicks(a, b, n) {
    const span = b - a;
    if (!(span > 0)) return [a];
    const raw = span / Math.max(n, 1);
    const pow = Math.pow(10, Math.floor(Math.log10(raw)));
    const err = raw / pow;
    const step = err >= 7 ? 10 * pow : err >= 3 ? 5 * pow : err >= 1.5 ? 2 * pow : pow;
    const t0 = Math.ceil((a - step * 1e-9) / step) * step;
    const ticks = [];
    for (let t = t0; t <= b + step * 1e-6; t += step) ticks.push(t);
    return ticks;
  }

  function plotLadder(result, title, viewOpt) {
    const { arr, rows, cols } = result;
    const view = normalizeView(result, viewOpt);
    const lim = limForView(result, view);
    const nP = view.p1 - view.p0;
    const nE = view.e1 - view.e0;
    const showProfiles = nP <= 20;
    const nCols = nP <= 16 ? 4 : 6;
    const nRows = showProfiles ? Math.ceil(nP / nCols) : 1;
    const leftPad = 58;
    const barW = 18;
    const top = 36;
    const bot = 48;
    const wfH = 420;
    const wfW = Math.min(1100, Math.max(520, Math.round(36 * nP)));
    const leftW = leftPad + wfW + barW + 70;
    const w = Math.max(16 + leftW + 16, showProfiles ? 780 : 16 + leftW + 16);
    const profH = showProfiles ? Math.max(200, nRows * 96) : 0;
    const h = top + wfH + (showProfiles ? 16 + profH : 0) + bot;
    const canvas = document.createElement("canvas");
    canvas.width = w;
    canvas.height = h;
    const ctx = canvas.getContext("2d");
    ctx.fillStyle = "#ffffff";
    ctx.fillRect(0, 0, w, h);
    ctx.font = "12px ui-sans-serif, system-ui, sans-serif";
    ctx.fillStyle = "#1c1915";
    ctx.fillText(title || "level ladder", 16, 22);

    const x0 = 16 + leftPad;
    const y0 = top;
    const img = ctx.createImageData(wfW, wfH);
    for (let py = 0; py < wfH; py++) {
      const e = view.e0 + ((wfH - 1 - py) / (wfH - 1)) * nE;
      const ie = Math.min(cols - 1, Math.max(0, Math.floor(e)));
      for (let px = 0; px < wfW; px++) {
        const p = view.p0 + ((px + 0.5) / wfW) * nP;
        const ip = Math.min(rows - 1, Math.max(0, Math.floor(p)));
        const v = arr[ip * cols + ie];
        const [r, g, b] = colorFor(v, lim);
        const o = (py * wfW + px) * 4;
        img.data[o] = r;
        img.data[o + 1] = g;
        img.data[o + 2] = b;
        img.data[o + 3] = 255;
      }
    }
    ctx.putImageData(img, x0, y0);
    ctx.save();
    ctx.beginPath();
    ctx.rect(x0, y0, wfW, wfH);
    ctx.clip();
    ctx.strokeStyle = "#111111";
    ctx.lineWidth = 1.2;
    const pathsMeta = result.paths || [];
    for (let i = view.p0; i < view.p1; i++) {
      const xa = x0 + ((i - view.p0) / nP) * wfW + 0.5;
      const xb = x0 + ((i + 1 - view.p0) / nP) * wfW - 0.5;
      let energies = [];
      const chain = (pathsMeta[i] || {}).chain;
      if (chain && chain.length) {
        energies = chain.map((e) => Math.round(e)).filter((e) => e > 0 && e < cols);
      } else {
        const row = arr.subarray(i * cols, (i + 1) * cols);
        for (let ie = 1; ie < cols - 1; ie++) {
          if (row[ie] === 0 && row[ie - 1] > 0 && row[ie + 1] > 0) energies.push(ie);
        }
        let lastPos = -1;
        for (let ie = 0; ie < cols; ie++) if (row[ie] > 0) lastPos = ie;
        if (lastPos >= 0) energies.push(lastPos + 1);
      }
      for (const e of energies) {
        if (e < view.e0 || e > view.e1) continue;
        const y = y0 + wfH * (1 - (e - view.e0) / nE);
        ctx.beginPath();
        ctx.moveTo(xa, y);
        ctx.lineTo(xb, y);
        ctx.stroke();
      }
    }
    ctx.restore();
    ctx.strokeStyle = "#c9c0b0";
    ctx.strokeRect(x0, y0, wfW, wfH);
    ctx.fillStyle = "#1c1915";
    ctx.textAlign = "right";
    for (const tick of niceTicks(view.e0, view.e1, 5)) {
      const y = y0 + wfH * (1 - (tick - view.e0) / nE);
      ctx.fillText(String(Math.round(tick)), x0 - 6, y + 4);
    }
    ctx.save();
    ctx.translate(18, y0 + wfH / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.textAlign = "center";
    ctx.fillText("Energy (keV)", 0, 0);
    ctx.restore();
    ctx.textAlign = "center";
    ctx.fillText("Path", x0 + wfW / 2, y0 + wfH + 28);
    if (nP <= 24) {
      for (let i = view.p0; i < view.p1; i++) {
        const x = x0 + ((i - view.p0 + 0.5) / nP) * wfW;
        ctx.fillText(String(i + 1), x, y0 + wfH + 14);
      }
    }

    const bx = x0 + wfW + 8;
    for (let py = 0; py < wfH; py++) {
      const t = 1 - py / (wfH - 1);
      const v = -lim + 2 * lim * t;
      const [r, g, b] = colorFor(v, lim);
      ctx.fillStyle = `rgb(${r},${g},${b})`;
      ctx.fillRect(bx, y0 + py, barW, 1);
    }
    ctx.strokeStyle = "#c9c0b0";
    ctx.strokeRect(bx, y0, barW, wfH);
    ctx.fillStyle = "#1c1915";
    ctx.textAlign = "left";
    ctx.fillText(lim.toPrecision(3), bx + barW + 4, y0 + 10);
    ctx.fillText("0", bx + barW + 4, y0 + wfH / 2);
    ctx.fillText((-lim).toPrecision(3), bx + barW + 4, y0 + wfH);
    ctx.fillText("W × Iγ/100", bx - 4, y0 - 8);

    if (showProfiles) {
      const rx = 16;
      const rightW = w - 32;
      const cellW = rightW / nCols;
      const cellH = profH / nRows;
      const profTop = y0 + wfH + 28;
      const iE0 = Math.max(0, Math.floor(view.e0));
      const iE1 = Math.min(cols, Math.ceil(view.e1));
      const step = Math.max(1, Math.floor((iE1 - iE0) / 400));
      for (let k = 0; k < nP; k++) {
        const i = view.p0 + k;
        const r = Math.floor(k / nCols);
        const c = k % nCols;
        const cx = rx + c * cellW + 8;
        const cy = profTop + r * cellH + 16;
        const cw = cellW - 16;
        const ch = cellH - 28;
        ctx.strokeStyle = "#ddd4c6";
        ctx.strokeRect(cx, cy, cw, ch);
        ctx.fillStyle = "#6b6258";
        ctx.textAlign = "left";
        ctx.fillText("path " + (i + 1), cx, cy - 4);
        const xMid = cx + cw * 0.5;
        ctx.strokeStyle = "#ccc4b6";
        ctx.beginPath();
        ctx.moveTo(xMid, cy);
        ctx.lineTo(xMid, cy + ch);
        ctx.stroke();
        const chain = ((result.paths || [])[i] || {}).chain || [];
        ctx.strokeStyle = "#888";
        ctx.setLineDash([3, 3]);
        for (const e of chain) {
          if (e < view.e0 || e > view.e1) continue;
          const y = cy + ch * (1 - (e - view.e0) / nE);
          ctx.beginPath();
          ctx.moveTo(cx, y);
          ctx.lineTo(cx + cw, y);
          ctx.stroke();
        }
        ctx.setLineDash([]);
        ctx.beginPath();
        ctx.strokeStyle = "#1c1915";
        ctx.lineWidth = 1;
        let started = false;
        for (let ie = iE0; ie < iE1; ie += step) {
          const v = arr[i * cols + ie];
          const x = cx + cw * ((v - (-lim)) / (2 * lim));
          const y = cy + ch * (1 - (ie - view.e0) / nE);
          if (!started) { ctx.moveTo(x, y); started = true; }
          else ctx.lineTo(x, y);
        }
        ctx.stroke();
        ctx.lineWidth = 1;
      }
    }

    return {
      canvas,
      layout: { x0, y0, wfW, wfH, w, h },
      view,
      lim,
    };
  }

  function mountInteractiveLadder(host, result, title, opts) {
    if (host._ladderAbort) host._ladderAbort.abort();
    const ac = new AbortController();
    const signal = ac.signal;
    host._ladderAbort = ac;
    host.innerHTML = "";
    const plotCanvas = document.createElement("canvas");
    const overlay = document.createElement("canvas");
    overlay.className = "overlay";
    overlay.setAttribute("aria-label", "level-ladder zoom");
    host.appendChild(plotCanvas);
    host.appendChild(overlay);

    let view = normalizeView(result, (opts && opts.view) || null);
    let layout = null;
    let drag = null;

    function redraw() {
      const drawn = plotLadder(result, title, view);
      const src = drawn.canvas;
      plotCanvas.width = src.width;
      plotCanvas.height = src.height;
      plotCanvas.getContext("2d").drawImage(src, 0, 0);
      overlay.width = src.width;
      overlay.height = src.height;
      layout = drawn.layout;
      view = drawn.view;
      if (opts && opts.onRedraw) {
        opts.onRedraw({ view, lim: drawn.lim, canvas: plotCanvas });
      }
    }

    function toCanvas(ev) {
      const r = overlay.getBoundingClientRect();
      return {
        x: (ev.clientX - r.left) * (overlay.width / Math.max(r.width, 1)),
        y: (ev.clientY - r.top) * (overlay.height / Math.max(r.height, 1)),
      };
    }

    function inWf(p) {
      if (!layout) return false;
      return p.x >= layout.x0 && p.x <= layout.x0 + layout.wfW &&
        p.y >= layout.y0 && p.y <= layout.y0 + layout.wfH;
    }

    function xyToData(p) {
      const nP = view.p1 - view.p0;
      const nE = view.e1 - view.e0;
      return {
        path: view.p0 + ((p.x - layout.x0) / layout.wfW) * nP,
        e: view.e0 + (1 - (p.y - layout.y0) / layout.wfH) * nE,
      };
    }

    function clearOverlay() {
      overlay.getContext("2d").clearRect(0, 0, overlay.width, overlay.height);
    }

    overlay.addEventListener("mousedown", (ev) => {
      if (ev.button !== 0) return;
      const p = toCanvas(ev);
      if (!inWf(p)) return;
      drag = { x0: p.x, y0: p.y, x1: p.x, y1: p.y };
      ev.preventDefault();
    }, { signal });

    window.addEventListener("mousemove", (ev) => {
      if (!drag) return;
      const p = toCanvas(ev);
      drag.x1 = p.x;
      drag.y1 = p.y;
      const ctx = overlay.getContext("2d");
      ctx.clearRect(0, 0, overlay.width, overlay.height);
      const x = Math.min(drag.x0, drag.x1);
      const y = Math.min(drag.y0, drag.y1);
      const bw = Math.abs(drag.x1 - drag.x0);
      const bh = Math.abs(drag.y1 - drag.y0);
      ctx.fillStyle = "rgba(196,92,38,0.12)";
      ctx.fillRect(x, y, bw, bh);
      ctx.strokeStyle = "#c45c26";
      ctx.lineWidth = 1.5;
      ctx.setLineDash([5, 4]);
      ctx.strokeRect(x, y, bw, bh);
    }, { signal });

    window.addEventListener("mouseup", () => {
      if (!drag) return;
      const box = drag;
      drag = null;
      clearOverlay();
      if (Math.abs(box.x1 - box.x0) < 8 || Math.abs(box.y1 - box.y0) < 8) return;
      const a = xyToData({ x: Math.min(box.x0, box.x1), y: Math.min(box.y0, box.y1) });
      const b = xyToData({ x: Math.max(box.x0, box.x1), y: Math.max(box.y0, box.y1) });
      view = normalizeView(result, {
        p0: Math.min(a.path, b.path),
        p1: Math.max(a.path, b.path),
        e0: Math.min(a.e, b.e),
        e1: Math.max(a.e, b.e),
      });
      redraw();
    }, { signal });

    overlay.addEventListener("dblclick", () => {
      view = normalizeView(result, null);
      redraw();
    }, { signal });

    overlay.addEventListener("wheel", (ev) => {
      ev.preventDefault();
      const p = toCanvas(ev);
      if (!inWf(p)) return;
      const d = xyToData(p);
      const factor = ev.deltaY > 0 ? 1.25 : 0.8;
      const nP = view.p1 - view.p0;
      const nE = view.e1 - view.e0;
      const newNE = Math.max(1, nE * factor);
      const newNP = Math.max(1, nP * factor);
      const fx = clamp((d.path - view.p0) / nP, 0, 1);
      const fy = clamp((d.e - view.e0) / nE, 0, 1);
      view = normalizeView(result, {
        p0: d.path - fx * newNP,
        p1: d.path + (1 - fx) * newNP,
        e0: d.e - fy * newNE,
        e1: d.e + (1 - fy) * newNE,
      });
      redraw();
    }, { signal, passive: false });

    redraw();
    return {
      reset() { view = normalizeView(result, null); redraw(); },
      getView() { return { ...view }; },
      canvas: plotCanvas,
      abort() { ac.abort(); },
    };
  }

  function encodeNpy(arr, rows, cols) {
    let header = "{'descr': '<f4', 'fortran_order': False, 'shape': (" + rows + ", " + cols + "), }";
    while ((10 + header.length + 1) % 64 !== 0) header += " ";
    header += "\n";
    const headerLen = header.length;
    const buf = new ArrayBuffer(10 + headerLen + arr.byteLength);
    const u8 = new Uint8Array(buf);
    u8.set([0x93, 78, 85, 77, 80, 89, 1, 0, headerLen & 255, (headerLen >> 8) & 255]);
    for (let i = 0; i < headerLen; i++) u8[10 + i] = header.charCodeAt(i);
    new Uint8Array(buf, 10 + headerLen).set(new Uint8Array(arr.buffer, arr.byteOffset, arr.byteLength));
    return new Blob([buf], { type: "application/octet-stream" });
  }

  const SAMPLE_QBETA = `<?xml version="1.0"?>
<Nuclide AtomicNumber="40" AtomicMass="98" QBeta="2238" d_QBeta="10">
	<Level Energy="0.00" Spin="0.00" Parity="+" SpinParity="0+" HalfLifeTime="30.700000000" d_T12="4" TimeUnit="S" Origin="Database">
		<Transition Type="B-" TransitionQValue="2238.00" Intensity="30" d_Intensity="3" Origin="Database">
			<TargetLevel Energy="0.00" AtomicNumber="41" AtomicMass="98" />
		</Transition>
		<Transition Type="B-" TransitionQValue="1888.00" Intensity="20" d_Intensity="3" Origin="Added">
			<TargetLevel Energy="350.00" AtomicNumber="41" AtomicMass="98" />
		</Transition>
		<Transition Type="B-" TransitionQValue="1438.00" Intensity="10" d_Intensity="3" Origin="Added">
			<TargetLevel Energy="800.00" AtomicNumber="41" AtomicMass="98" />
		</Transition>
		<Transition Type="B-" TransitionQValue="738.00" Intensity="30" d_Intensity="1" Origin="Added">
			<TargetLevel Energy="1500.00" AtomicNumber="41" AtomicMass="98" />
		</Transition>
		<Transition Type="B-" TransitionQValue="638.00" Intensity="10" d_Intensity="3" Origin="Added">
			<TargetLevel Energy="1600.00" AtomicNumber="41" AtomicMass="98" />
		</Transition>
	</Level>
</Nuclide>
`;

  const SAMPLE_GAMMA = `<?xml version="1.0"?>
<Nuclide AtomicNumber="41" AtomicMass="98" QBeta="0">
	<Level Energy="0.00" Spin="1.00" Parity="+" SpinParity="1+" HalfLifeTime="0.000000000" TimeUnit="S" Origin="Database">
		<Transition Type="Fake" TransitionQValue="0.00" Intensity="100.0000" Origin="Added">
			<TargetLevel Energy="0.00" AtomicNumber="0" AtomicMass="0" />
		</Transition>
	</Level>
	<Level Energy="350.00" Spin="-1.00" Parity="" HalfLifeTime="0.000000000" TimeUnit="S" Origin="Added">
		<Transition Type="G" TransitionQValue="350.00" Intensity="100.0000" Origin="Added">
			<TargetLevel Energy="0.00" AtomicNumber="41" AtomicMass="98" />
		</Transition>
	</Level>
	<Level Energy="800.00" Spin="-1.00" Parity="" HalfLifeTime="0.000000000" TimeUnit="S" Origin="Added">
		<Transition Type="G" TransitionQValue="800.00" Intensity="84.8912" Origin="Added">
			<TargetLevel Energy="0.00" AtomicNumber="41" AtomicMass="98" />
		</Transition>
		<Transition Type="G" TransitionQValue="450.00" Intensity="15.1088" Origin="Added">
			<TargetLevel Energy="350.00" AtomicNumber="41" AtomicMass="98" />
		</Transition>
	</Level>
	<Level Energy="1500.00" Spin="-1.00" Parity="" HalfLifeTime="0.000000000" TimeUnit="S" Origin="Added">
		<Transition Type="G" TransitionQValue="1500.00" Intensity="90.7746" Origin="Added">
			<TargetLevel Energy="0.00" AtomicNumber="41" AtomicMass="98" />
		</Transition>
		<Transition Type="G" TransitionQValue="700.00" Intensity="9.2254" Origin="Added">
			<TargetLevel Energy="800.00" AtomicNumber="41" AtomicMass="98" />
		</Transition>
	</Level>
	<Level Energy="1600.00" Spin="-1.00" Parity="" HalfLifeTime="0.000000000" TimeUnit="S" Origin="Added">
		<Transition Type="G" TransitionQValue="1600.00" Intensity="53.7533" Origin="Added">
			<TargetLevel Energy="0.00" AtomicNumber="41" AtomicMass="98" />
		</Transition>
		<Transition Type="G" TransitionQValue="1250.00" Intensity="32.8084" Origin="Added">
			<TargetLevel Energy="350.00" AtomicNumber="41" AtomicMass="98" />
		</Transition>
		<Transition Type="G" TransitionQValue="800.00" Intensity="13.4383" Origin="Added">
			<TargetLevel Energy="800.00" AtomicNumber="41" AtomicMass="98" />
		</Transition>
	</Level>
</Nuclide>
`;

  function inspectRecord(n) {
    const gammas = [];
    for (const src of Object.keys(n.gammas).map(Number).sort((a, b) => a - b)) {
      for (const g of n.gammas[src] || []) {
        gammas.push({ from_e: g.from_e, to_e: g.to_e, q: g.q, intensity: g.intensity });
      }
    }
    return {
      name: n.name,
      z: n.z,
      a: n.a,
      q_beta: n.q_beta,
      levels: n.levels,
      betas: n.betas.map((b) => ({ ...b })),
      gammas,
      text: n.text,
      _nuclide: n,
    };
  }

  global.Ladder = {
    parseNuclide,
    runTranslation,
    plotLadder,
    mountInteractiveLadder,
    encodeNpy,
    inspectRecord,
    SAMPLE_QBETA,
    SAMPLE_GAMMA,
  };
})(typeof window !== "undefined" ? window : globalThis);
