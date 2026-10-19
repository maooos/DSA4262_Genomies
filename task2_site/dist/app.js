import {
  parseReadSupport,
  mean,
  finite,
  escapeHTML as esc,
  overviewRows,
  orientPair,
  filterCandidates,
  transcriptRows,
  positionGroups,
  csvText,
} from "./core.mjs";
const $ = (id) => document.getElementById(id);
const fmt = (n) => (finite(n) ? n.toLocaleString("en-US") : "—");
const dec = (n, d = 3) => (finite(n) ? n.toFixed(d) : "—");
const pct = (n) => (finite(n) ? (n * 100).toFixed(2) + "%" : "—");
const compact = (n) => (n >= 1e6 ? (n / 1e6).toFixed(2) + "M" : fmt(n));
const palette = {
  both: "#367c60",
  a: "#b78042",
  b: "#537aa6",
  neither: "#b6c1b2",
};
let study,
  index,
  aliases,
  pair,
  transcript = null,
  filteredCandidates = [];
let state = {
  page: "overview",
  coverage: 20,
  threshold: 0,
  a: "A549",
  b: "HepG2",
  tx: "",
  position: null,
};
let candidatePage = 0,
  positionPage = 0,
  pairRequest = 0,
  overviewRequest = 0,
  txRequest = 0;
let visibleOverviewRows = [];
let scatterPoints = [],
  scatterBox = {},
  toastTimer;
const pairCache = new Map(),
  overviewCache = new Map(),
  shardCache = new Map();
const pageText = {
  overview: [
    "THE SG-NEX COLLECTION",
    "A closer look at RNA modification.",
    "Explore predicted m6A sites across cell lines. Compare patterns, then inspect the evidence behind each site.",
  ],
  compare: [
    "COMPARE & INVESTIGATE",
    "Find the differences worth following.",
    "Compare the same sites in two cell lines, with coverage and replicate support built into every comparison.",
  ],
  transcript: [
    "TRANSCRIPT EXPLORER",
    "Follow the evidence along a transcript.",
    "Inspect individual sequencing runs, identify supported sites and see exactly where observations are missing.",
  ],
  methods: [
    "METHODS & READING GUIDE",
    "Understand the evidence behind the view.",
    "A practical guide to the scores, comparison rules and limits of this exploratory analysis.",
  ],
};
function toast(message) {
  $("toast").textContent = message;
  $("toast").hidden = false;
  clearTimeout(toastTimer);
  toastTimer = setTimeout(() => ($("toast").hidden = true), 4500);
}
async function json(url) {
  const r = await fetch(url);
  if (!r.ok) {
    if (url.startsWith("api/")) {
      const detail = await r.json().catch(() => ({}));
      throw new Error(detail.error || "Custom read filtering needs the local query service. Start it with docker compose up --build -d, then retry.");
    }
    throw new Error(`Could not load ${url} (HTTP ${r.status}).`);
  }
  return r.json();
}
async function compressedJSON(url) {
  const r = await fetch(url);
  if (!r.ok) throw new Error(`Transcript data unavailable (HTTP ${r.status}).`);
  const b = await r.arrayBuffer(),
    bytes = new Uint8Array(b);
  if (bytes[0] === 31 && bytes[1] === 139) {
    if (!("DecompressionStream" in window))
      throw new Error(
        "Please use a current Chrome, Edge, Firefox or Safari browser to open compressed transcript data.",
      );
    return new Response(
      new Blob([b]).stream().pipeThrough(new DecompressionStream("gzip")),
    ).json();
  }
  return JSON.parse(new TextDecoder().decode(b));
}
function download(name, content, type = "text/csv;charset=utf-8") {
  const url = URL.createObjectURL(new Blob([content], { type }));
  const a = document.createElement("a");
  a.href = url;
  a.download = name;
  a.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
  toast("Download ready. Current filters are included in the export.");
}
function csv(name, rows) {
  if (!rows.length) {
    toast("No rows match the current filters.");
    return;
  }
  download(name, csvText(rows, Object.keys(rows[0])));
}
function stat(label, value, note) {
  return `<div class="stat"><span class="stat-label">${esc(label)}</span><strong class="stat-value">${esc(value)}</strong><span class="stat-note">${esc(note)}</span></div>`;
}
function color(value) {
  if (!finite(value)) return "#f2f5ed";
  const t = Math.max(0, Math.min(1, value));
  return `rgb(${Math.round(242 - 211 * t)},${Math.round(247 - 139 * t)},${Math.round(231 - 157 * t)})`;
}
function options(cells) {
  return cells
    .map((c) => `<option value="${esc(c)}">${esc(c)}</option>`)
    .join("");
}
function message(id, text, error = false) {
  $(id).innerHTML = text
    ? `<div class="${error ? "error" : "loading"}">${esc(text)}</div>`
    : "";
}
function setURL() {
  const p = new URLSearchParams();
  if (state.coverage !== 20) p.set("reads", state.coverage);
  if (state.page === "overview" && state.threshold !== study.threshold)
    p.set("threshold", state.threshold);
  if (state.page === "compare") {
    p.set("a", state.a);
    p.set("b", state.b);
  }
  if (state.page === "transcript" && state.tx) {
    p.set("id", state.tx);
    if (state.position !== null) p.set("position", state.position);
  }
  history.replaceState(null, "", `#${state.page}${p.size ? "?" + p : ""}`);
}
function readURL() {
  const [page, params] = location.hash.slice(1).split("?"),
    p = new URLSearchParams(params || "");
  state.page = pageText[page] ? page : "overview";
  if (p.has("reads")) {
    const reads = parseReadSupport(p.get("reads"));
    state.coverage = reads ?? 20;
    if (reads === null) toast("Invalid read support in the link. Showing the default of 20 reads.");
  }
  if (study.cell_lines.includes(p.get("a"))) state.a = p.get("a");
  if (study.cell_lines.includes(p.get("b"))) state.b = p.get("b");
  const threshold = Number(p.get("threshold"));
  if ([0.15, study.threshold, 0.3, 0.5].includes(threshold))
    state.threshold = threshold;
  if (p.get("id")) state.tx = p.get("id");
  const pos = p.get("position");
  if (pos !== null && Number.isInteger(Number(pos)))
    state.position = Number(pos);
}
async function route() {
  ++overviewRequest;
  ++pairRequest;
  ++txRequest;
  readURL();
  window.scrollTo(0, 0);
  const [eyebrow, title, desc] = pageText[state.page];
  $("eyebrow").textContent = eyebrow;
  $("page-title").textContent = title;
  $("page-description").textContent = desc;
  document
    .querySelectorAll(".page")
    .forEach((el) => (el.hidden = el.id !== `page-${state.page}`));
  document.querySelectorAll("nav a").forEach((a) => {
    a.classList.toggle("active", a.dataset.page === state.page);
    if (a.dataset.page === state.page) {
      a.setAttribute("aria-current", "page");
      $("breadcrumb").textContent = a.textContent.slice(1).trim();
    } else a.removeAttribute("aria-current");
  });
  syncControls();
  if (state.page === "overview") renderOverview();
  if (state.page === "compare") await renderCompare();
  if (state.page === "transcript" && state.tx) await loadTranscript(state.tx);
  setURL();
}
function syncControls() {
  $("minimum-reads").value = state.coverage;
  $("minimum-reads").removeAttribute("aria-invalid");
  $("coverage-error").hidden = true;
  document
    .querySelectorAll("[data-coverage]")
    .forEach((b) =>
      b.setAttribute(
        "aria-pressed",
        Number(b.dataset.coverage) === state.coverage,
      ),
    );
  $("cell-a").value = state.a;
  $("cell-b").value = state.b;
  $("overview-threshold").value = state.threshold;
}
async function renderOverview() {
  const request = ++overviewRequest;
  const coverage = state.coverage;
  $("overview-content").hidden = true;
  visibleOverviewRows = [];
  $("download-overview").disabled = true;
  let source = study, scenario = study.scenarios[coverage];
  try {
    if (!scenario) {
      message("overview-message", `Calculating summaries and matched-site comparisons at ≥${coverage} reads… This can take a moment.`);
      if (!overviewCache.has(coverage)) {
        const result = await json(`api/overview?reads=${coverage}`);
        if (overviewCache.size >= 8) overviewCache.delete(overviewCache.keys().next().value);
        overviewCache.set(coverage, result);
      }
      source = overviewCache.get(coverage);
      scenario = source.scenario;
    }
    if (request !== overviewRequest || state.page !== "overview") return;
  } catch (e) {
    if (request === overviewRequest) message("overview-message", e.message, true);
    return;
  }
  message("overview-message", "");
  $("overview-content").hidden = false;
  $("download-overview").disabled = false;
  const rows = overviewRows(source, state.coverage, state.threshold),
    eligible = rows.reduce((s, r) => s + r.eligible_sites, 0);
  visibleOverviewRows = rows;
  $("overview-stats").innerHTML =
    stat(
      "Cell lines",
      fmt(study.cell_lines.length),
      "Across the processed SG-NEx collection",
    ) +
    stat(
      "Sequencing runs",
      fmt(study.samples.length),
      "Technical runs retain their sample identity",
    ) +
    stat(
      "Eligible observations",
      compact(eligible),
      `Run–site observations · ≥${state.coverage} reads`,
    ) +
    stat(
      "Transcript catalogue",
      fmt(Object.keys(index).length),
      "Contains sites exported at ≥10 reads",
    );
  const max =
    Math.max(0.1, ...rows.map((r) => r.positive_fraction || 0)) * 1.15;
  $("fraction-chart").innerHTML =
    study.cell_lines
      .map((c) => {
        const r = rows.filter(
            (r) => r.cell_line === c && finite(r.positive_fraction),
          ),
          v = r.map((r) => r.positive_fraction);
        const min = v.length ? Math.min(...v) : 0,
          high = v.length ? Math.max(...v) : 0;
        return `<div class="fraction-row"><strong>${esc(c)}</strong><div class="fraction-track"><span class="fraction-range" style="left:${(min / max) * 100}%;width:${((high - min) / max) * 100}%"></span>${r.map((s) => `<button class="fraction-dot" style="left:${(s.positive_fraction / max) * 100}%" aria-label="${esc(c)} replicate ${s.replicate} run ${s.run}: ${pct(s.positive_fraction)} predicted positive; ${fmt(s.eligible_sites)} eligible sites" title="${esc(c)} · replicate ${s.replicate}, run ${s.run}\n${pct(s.positive_fraction)} positive (${fmt(s.positive_sites)} / ${fmt(s.eligible_sites)} sites)"></button>`).join("")}</div><small>${r.length} runs</small></div>`;
      })
      .join("") +
    `<div class="axis-labels">${[0, 0.25, 0.5, 0.75, 1].map((x) => `<span>${(x * max * 100).toFixed(1)}%</span>`).join("")}</div>`;
  const values = rows.map((r) => r.positive_fraction).filter(finite);
  $("overview-insight").innerHTML = values.length
    ? `At ≥${state.coverage} reads and threshold ${dec(state.threshold, 4)}, <strong>${pct(Math.min(...values))}–${pct(Math.max(...values))}</strong> of eligible sites per run are predicted positive.`
    : "No eligible observations at these settings.";
  $("overlap-heatmap").innerHTML =
    `<div class="heatmap"><span></span>${study.cell_lines.map((c) => `<span class="heat-label">${esc(c)}</span>`).join("")}${study.cell_lines
      .map(
        (a) =>
          `<span class="heat-label">${esc(a)}</span>` +
          study.cell_lines
            .map((b) => {
              if (a === b) return '<span class="heat-empty">—</span>';
              const p = scenario.pairs.find(
                (p) =>
                  (p.cell_a === a && p.cell_b === b) ||
                  (p.cell_a === b && p.cell_b === a),
              );
              return `<button class="heat-cell" data-pair-a="${esc(a)}" data-pair-b="${esc(b)}" style="background:${color(p?.jaccard)}" title="${esc(a)} / ${esc(b)}: ${dec(p?.jaccard)} overlap on ${fmt(p?.common_sites)} shared sites" aria-label="Compare ${esc(a)} and ${esc(b)}; Jaccard ${dec(p?.jaccard)}">${dec(p?.jaccard, 2)}</button>`;
            })
            .join(""),
      )
      .join("")}</div>`;
  $("agreement").innerHTML =
    `<table class="agreement-table"><thead><tr><th scope="col">Cell line</th><th scope="col">Replicate groups</th><th scope="col">Shared eligible sites</th><th scope="col">Score correlation (Spearman)</th><th scope="col">Positive-site overlap</th></tr></thead><tbody>${scenario.agreement.map((r) => `<tr><td>${esc(r.cell_line)}</td><td>${r.replicate_a} ↔ ${r.replicate_b}</td><td>${fmt(r.common_sites)}</td><td>${dec(r.spearman)} <span class="mini-bar"><i style="width:${Math.max(0, r.spearman || 0) * 100}%"></i></span></td><td>${dec(r.jaccard)}</td></tr>`).join("")}</tbody></table>`;
}
async function renderCompare() {
  const request = ++pairRequest;
  pair = null;
  $("compare-content").hidden = true;
  if (state.a === state.b) {
    message(
      "compare-message",
      "Choose two different cell lines to make a comparison.",
      true,
    );
    return;
  }
  message("compare-message", `Loading matched sites at ≥${state.coverage} reads… Custom comparisons can take a moment.`);
  const ordered = [state.a, state.b].sort(
    (a, b) => study.cell_lines.indexOf(a) - study.cell_lines.indexOf(b),
  );
  const key = `${ordered.join("__")}_min${state.coverage}`;
  try {
    if (!pairCache.has(key)) {
      const url = study.scenarios[state.coverage]
        ? `data/pairs/${key}.json`
        : `api/pair?${new URLSearchParams({reads: state.coverage, a: ordered[0], b: ordered[1]})}`;
      const result = await json(url);
      if (pairCache.size >= 12) pairCache.delete(pairCache.keys().next().value);
      pairCache.set(key, result);
    }
    if (request !== pairRequest) return;
    pair = orientPair(pairCache.get(key), state.a);
    message("compare-message", "");
    $("compare-content").hidden = false;
    const p = pair.summary;
    $("compare-stats").innerHTML =
      stat(
        "Shared eligible sites",
        fmt(p.common_sites),
        "Matched transcript, position and sequence",
      ) +
      stat(
        "Positive in both",
        fmt(p.positive_both),
        "At the saved H04 threshold",
      ) +
      stat(
        "Positive-site overlap",
        dec(p.jaccard),
        "Jaccard · 0 = none, 1 = complete",
      ) +
      stat(
        "Consistent candidates",
        fmt(p.consistent_candidates),
        "Opposite calls across qualifying replicates",
      );
    $("candidate-direction").options[1].textContent = state.a;
    $("candidate-direction").options[2].textContent = state.b;
    const facts = [
      ["both", "Positive in both", p.positive_both],
      ["a", `Positive only in ${state.a}`, p.positive_a_only],
      ["b", `Positive only in ${state.b}`, p.positive_b_only],
      ["neither", "Below threshold in both", p.negative_both],
    ];
    $("pair-explanation").innerHTML =
      facts
        .map(
          ([key, label, value]) =>
            `<div class="pair-fact"><span><i class="key" style="background:${palette[key]}"></i>${esc(label)}</span><b>${fmt(value)}</b></div>`,
        )
        .join("") +
      `<p class="footnote">Score-ranking correlation: <strong>${dec(p.spearman)}</strong> (Spearman). Counts above use all ${fmt(p.common_sites)} shared eligible sites.</p>`;
    $("scatter-caption").textContent =
      `Showing ${fmt(pair.scatter.length)} sampled sites out of ${fmt(pair.scatter_total)} shared eligible sites. Summary counts use the full comparison. Candidate links below provide keyboard-accessible exploration.`;
    candidatePage = 0;
    renderCandidates();
    requestAnimationFrame(drawScatter);
  } catch (e) {
    if (request === pairRequest)
      message("compare-message", e.message + " Refresh to retry.", true);
  }
}
function candidateRows() {
  return filterCandidates(
    pair?.candidates || [],
    $("candidate-search").value,
    $("candidate-direction").value,
    $("candidate-sort").value,
  );
}
function renderCandidates() {
  if (!pair) return;
  filteredCandidates = candidateRows();
  const pages = Math.max(1, Math.ceil(filteredCandidates.length / 20));
  candidatePage = Math.min(candidatePage, pages - 1);
  const rows = filteredCandidates.slice(
    candidatePage * 20,
    (candidatePage + 1) * 20,
  );
  $("candidate-table").innerHTML = rows.length
    ? `<table><thead><tr><th scope="col">Transcript ↗</th><th scope="col">Position</th><th scope="col">7-mer</th><th scope="col">${esc(state.a)} score</th><th scope="col">${esc(state.b)} score</th><th scope="col">Difference A − B</th><th scope="col">Replicates A / B</th></tr></thead><tbody>${rows.map((r) => `<tr><td><button class="table-link" data-transcript="${esc(r.transcript_id)}" data-position="${r.transcript_position}">${esc(r.transcript_id)}</button></td><td>${fmt(r.transcript_position)}</td><td><code>${esc(r.sequence)}</code></td><td>${dec(r.mean_score_a)}</td><td>${dec(r.mean_score_b)}</td><td class="delta">${r.delta_score > 0 ? "+" : ""}${dec(r.delta_score)}</td><td>${r.n_replicates_a} / ${r.n_replicates_b}</td></tr>`).join("")}</tbody></table>`
    : '<div class="empty">No candidates match these filters. Try another transcript or score direction.</div>';
  $("candidate-count").textContent = rows.length
    ? `${candidatePage * 20 + 1}–${Math.min((candidatePage + 1) * 20, filteredCandidates.length)} of ${fmt(filteredCandidates.length)} matching candidates`
    : "0 matching candidates";
  $("candidate-page").textContent = `${candidatePage + 1} / ${pages}`;
  $("candidate-prev").disabled = candidatePage === 0;
  $("candidate-next").disabled = candidatePage >= pages - 1;
  $("download-candidates").disabled = !rows.length;
  $("candidate-limit").textContent =
    `Showing ${fmt(pair.candidates.length)} of ${fmt(pair.summary.consistent_candidates)} consistent candidates (display limit: ${fmt(pair.candidate_limit)}), ranked by absolute score difference. Search, sorting and download apply to this subset; summary counts use every shared eligible site.`;
}
function drawScatter() {
  if (!pair || state.page !== "compare") return;
  const canvas = $("scatter"),
    width = canvas.clientWidth;
  if (!width) return;
  const height = 390,
    dpr = devicePixelRatio || 1;
  canvas.width = width * dpr;
  canvas.height = height * dpr;
  canvas.style.height = height + "px";
  const ctx = canvas.getContext("2d");
  ctx.scale(dpr, dpr);
  const l = 49,
    t = 18,
    w = width - 69,
    h = height - 76;
  scatterBox = { l, t, w, h };
  ctx.font = "10px system-ui";
  ctx.lineWidth = 1;
  for (let i = 0; i <= 5; i++) {
    const v = i / 5,
      x = l + w * v,
      y = t + h * (1 - v);
    ctx.strokeStyle = "#e9eee4";
    ctx.beginPath();
    ctx.moveTo(x, t);
    ctx.lineTo(x, t + h);
    ctx.moveTo(l, y);
    ctx.lineTo(l + w, y);
    ctx.stroke();
    ctx.fillStyle = "#7e8a76";
    ctx.textAlign = "center";
    ctx.fillText(v.toFixed(1), x, t + h + 20);
    ctx.textAlign = "right";
    ctx.fillText(v.toFixed(1), l - 9, y + 3);
  }
  ctx.strokeStyle = "#879c81";
  ctx.setLineDash([4, 4]);
  ctx.beginPath();
  ctx.moveTo(l, t + h);
  ctx.lineTo(l + w, t);
  ctx.stroke();
  ctx.strokeStyle = "#b1bd9d";
  const cutoff = study.threshold;
  ctx.beginPath();
  ctx.moveTo(l + w * cutoff, t);
  ctx.lineTo(l + w * cutoff, t + h);
  ctx.moveTo(l, t + h * (1 - cutoff));
  ctx.lineTo(l + w, t + h * (1 - cutoff));
  ctx.stroke();
  ctx.setLineDash([]);
  scatterPoints = pair.scatter.map((r) => ({
    ...r,
    x: l + r.mean_score_a * w,
    y: t + (1 - r.mean_score_b) * h,
  }));
  for (const r of scatterPoints) {
    const a = r.mean_score_a >= cutoff,
      b = r.mean_score_b >= cutoff;
    ctx.fillStyle = palette[a && b ? "both" : a ? "a" : b ? "b" : "neither"];
    ctx.globalAlpha = a || b ? 0.8 : 0.4;
    ctx.beginPath();
    ctx.arc(r.x, r.y, 2.7, 0, Math.PI * 2);
    ctx.fill();
  }
  ctx.globalAlpha = 1;
  ctx.fillStyle = "#455c3b";
  ctx.textAlign = "center";
  ctx.font = "11px system-ui";
  ctx.fillText(`${state.a} · mean model score`, l + w / 2, height - 10);
  ctx.save();
  ctx.translate(12, t + h / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText(`${state.b} · mean model score`, 0, 0);
  ctx.restore();
}
function nearestPoint(event) {
  const r = $("scatter").getBoundingClientRect(),
    x = event.clientX - r.left,
    y = event.clientY - r.top;
  let found = null,
    d = 100;
  for (const p of scatterPoints) {
    const dist = (p.x - x) ** 2 + (p.y - y) ** 2;
    if (dist < d) {
      d = dist;
      found = p;
    }
  }
  return found;
}
function openTranscript(id, position = null) {
  state.tx = id;
  state.position = position;
  state.page = "transcript";
  const p = new URLSearchParams({ id, reads: state.coverage });
  if (position !== null) p.set("position", position);
  location.hash = "transcript?" + p;
}
async function loadTranscript(input) {
  const request = ++txRequest;
  const raw = input.trim(),
    canonical = aliases[raw] || raw;
  state.tx = raw;
  setURL();
  $("transcript-id").value = raw;
  transcript = null;
  $("transcript-content").hidden = true;
  message("transcript-message", "Loading transcript observations…");
  try {
    if (state.coverage < study.transcript_export_min_reads || !index[canonical]) {
      const result = await json(`api/transcript?${new URLSearchParams({id: raw})}`);
      if (request !== txRequest) return;
      if (!result.rows.length) throw new Error(`“${raw}” has no observations in this study. Check the transcript ID; missing does not mean unmodified.`);
      transcript = result;
    } else {
      const shard = index[canonical];
      if (!shardCache.has(shard)) {
        const d = await compressedJSON(`data/transcripts/${shard}.json.gz`);
        if (shardCache.size >= 4)
          shardCache.delete(shardCache.keys().next().value);
        shardCache.set(shard, d);
      }
      if (request !== txRequest) return;
      const rows = shardCache.get(shard)[canonical];
      if (!rows) throw new Error("The transcript index points to a missing record.");
      transcript = { id: canonical, rows };
    }
    positionPage = 0;
    message(
      "transcript-message",
      raw === transcript.id
        ? ""
        : `Showing ${transcript.id}; ${raw} is a known versioned alias.`,
    );
    $("transcript-content").hidden = false;
    renderTranscript();
    setURL();
  } catch (e) {
    if (request === txRequest)
      message("transcript-message", e.message + " Refresh to retry.", true);
  }
}
function currentTranscriptRows() {
  return transcriptRows(
    transcript.rows,
    study.samples,
    state.coverage,
    $("transcript-cell").value,
  );
}
function renderTranscript() {
  if (!transcript) return;
  const rows = currentTranscriptRows(),
    positions = [...new Set(rows.map((r) => r[1]))].sort((a, b) => a - b);
  const pages = Math.max(1, Math.ceil(positions.length / 24));
  positionPage = Math.min(positionPage, pages - 1);
  if (!positions.includes(state.position))
    state.position = positions[positionPage * 24] ?? null;
  else if (positionPage === 0)
    positionPage = Math.floor(positions.indexOf(state.position) / 24);
  const displayed = positions.slice(positionPage * 24, (positionPage + 1) * 24);
  $("transcript-heading").textContent = transcript.id;
  $("transcript-stats").innerHTML =
    stat(
      "Observed positions",
      fmt(positions.length),
      `In current cell-line and read-support filters`,
    ) +
    stat(
      "Qualifying runs",
      fmt(new Set(rows.map((r) => r[0])).size),
      "With at least one eligible observation",
    ) +
    stat(
      "Eligible observations",
      fmt(rows.length),
      "One record per sequencing run and site",
    ) +
    stat(
      "Positive observations",
      fmt(rows.filter((r) => r[4] >= study.threshold).length),
      "Run-level calls at the saved threshold",
    );
  const sampleIds = study.samples
    .map((s, i) => ({ s, i }))
    .filter(
      ({ s }) =>
        $("transcript-cell").value === "all" ||
        s.cell_line === $("transcript-cell").value,
    );
  const map = new Map(rows.map((r) => [`${r[0]}|${r[1]}`, r]));
  $("transcript-heatmap").innerHTML = positions.length
    ? `<table><thead><tr><th scope="col">Cell line · replicate / run</th>${displayed.map((p) => `<th scope="col"><button class="position-head" data-select-position="${p}" aria-label="Inspect position ${p}">${p}</button></th>`).join("")}</tr></thead><tbody>${sampleIds
        .map(
          ({ s, i }) =>
            `<tr><td>${esc(s.cell_line)} <span class="muted">· ${s.replicate} / ${s.run}</span></td>${displayed
              .map((p) => {
                const r = map.get(`${i}|${p}`);
                return r
                  ? `<td><button class="site-cell ${p === state.position ? "selected" : ""}" data-select-position="${p}" style="background:${color(r[4])};color:${r[4] > 0.65 ? "#fff" : "#254d33"}" aria-label="${esc(s.cell_line)} replicate ${s.replicate} run ${s.run}, position ${p}, score ${dec(r[4], 4)}, ${r[3]} reads" title="${dec(r[4], 4)} score · ${fmt(r[3])} reads">${dec(r[4], 2)}</button></td>`
                  : '<td><span class="site-missing" title="Missing or below selected read support">—</span></td>';
              })
              .join("")}</tr>`,
        )
        .join("")}</tbody></table>`
    : '<div class="empty">No observations meet these filters. Lower the read support or choose another cell line.</div>';
  $("position-count").textContent = positions.length
    ? `Positions ${positionPage * 24 + 1}–${Math.min((positionPage + 1) * 24, positions.length)} of ${fmt(positions.length)} observed positions`
    : "0 eligible positions";
  $("position-prev").disabled = positionPage === 0;
  $("position-next").disabled = positionPage >= pages - 1;
  $("download-transcript").disabled = !rows.length;
  renderPosition(rows);
}
function renderPosition(rows = currentTranscriptRows()) {
  const pos = state.position;
  $("position-heading").textContent =
    pos === null ? "No eligible position selected" : `Position ${fmt(pos)}`;
  const selected = rows.filter((r) => r[1] === pos),
    groups = positionGroups(rows, study.samples, pos, study.threshold);
  $("position-summary").innerHTML = groups.length
    ? `<div class="cell-support">${groups.map((g) => `<div class="support-item">${esc(g.cell)}<b>${dec(g.score)}</b><small>${g.replicates} replicate${g.replicates === 1 ? "" : "s"} · ${g.positiveReplicates} positive<br>${g.replicates >= 2 ? "Eligible replicate support" : "Limited replicate support"}<br><code>${esc(g.sequence)}</code></small></div>`).join("")}</div>`
    : "";
  $("position-table").innerHTML = selected.length
    ? `<table><thead><tr><th scope="col">Cell line</th><th scope="col">Replicate / run</th><th scope="col">Original transcript ID</th><th scope="col">7-mer</th><th scope="col">Reads</th><th scope="col">Model score</th><th scope="col">Prediction</th></tr></thead><tbody>${selected
        .map((r) => {
          const s = study.samples[r[0]],
            positive = r[4] >= study.threshold;
          return `<tr><td>${esc(s.cell_line)}</td><td>${s.replicate} / ${s.run}</td><td>${esc(r[5] || transcript.id)}</td><td><code>${esc(r[2])}</code></td><td>${fmt(r[3])}</td><td>${dec(r[4], 5)}</td><td><span class="${positive ? "positive" : "negative"}">${positive ? "Predicted positive" : "Below threshold"}</span></td></tr>`;
        })
        .join("")}</tbody></table>`
    : "";
  setURL();
}
function renderMethods() {
  const size = study.inventory.reduce((s, r) => s + r.size_bytes, 0);
  $("scope-details").innerHTML =
    `<p>This result bundle covers <strong>${study.samples.length} of ${study.expected_samples} processed direct-RNA runs</strong> across ${study.cell_lines.length} cell lines (${(size / 1e9).toFixed(2)} GB of source JSON). H04 generated ${fmt(study.total_scored_sites)} run–site predictions. This is coverage of the available processed m6Anet collection, not every raw SG-NEx experiment.</p><p>Analysis snapshot: <strong>${esc(study.created_at.slice(0, 10))}</strong>. Enter any positive whole number for minimum reads per run, then select Apply. The default is 20; 10 and 50 are quick sensitivity checks. Custom settings recalculate from the saved predictions. Lower coverage admits less-supported observations; the cutoff is not a guarantee of accuracy.</p>`;
  $("provenance").textContent = JSON.stringify(
    {
      analysis_created_at: study.created_at,
      status: study.status,
      model_sha256: study.model_sha256,
      threshold: study.threshold,
      expected_runs: study.expected_samples,
      missing_samples: study.missing_samples,
      identifier_harmonization: study.identifier_harmonization,
      sequence_conflict_sites: study.sequence_conflict_sites,
    },
    null,
    2,
  );
}
function wire() {
  window.addEventListener("hashchange", route);
  const applyCoverage = (coverage) => {
    state.coverage = coverage;
    syncControls();
    setURL();
    if (state.page === "overview") renderOverview();
    if (state.page === "compare") renderCompare();
    if (state.page === "transcript" && state.tx) {
      positionPage = 0;
      loadTranscript(state.tx);
    }
  };
  $("coverage-form").addEventListener("submit", (event) => {
    event.preventDefault();
    const coverage = parseReadSupport($("minimum-reads").value);
    if (coverage === null) {
      $("minimum-reads").setAttribute("aria-invalid", "true");
      $("coverage-error").textContent = "Enter a positive whole number of reads (1 or more). The current filter has not changed.";
      $("coverage-error").hidden = false;
      $("minimum-reads").focus();
      return;
    }
    applyCoverage(coverage);
  });
  $("minimum-reads").addEventListener("input", () => {
    $("minimum-reads").removeAttribute("aria-invalid");
    $("coverage-error").hidden = true;
  });
  document.querySelectorAll("[data-coverage]").forEach((b) =>
    b.addEventListener("click", () => applyCoverage(Number(b.dataset.coverage))),
  );
  $("overview-threshold").addEventListener("change", () => {
    state.threshold = Number($("overview-threshold").value);
    setURL();
    renderOverview();
  });
  $("reset").addEventListener("click", () => {
    state.coverage = 20;
    state.threshold = study.threshold;
    state.a = "A549";
    state.b = "HepG2";
    state.position = null;
    state.tx = "";
    $("candidate-search").value = "";
    $("candidate-direction").value = "all";
    $("candidate-sort").value = "difference";
    $("transcript-cell").value = "all";
    $("transcript-id").value = "";
    transcript = null;
    $("transcript-content").hidden = true;
    message("transcript-message", "");
    setURL();
    route();
  });
  for (const id of ["cell-a", "cell-b"])
    $(id).addEventListener("change", () => {
      state.a = $("cell-a").value;
      state.b = $("cell-b").value;
      setURL();
      renderCompare();
    });
  $("swap-cells").addEventListener("click", () => {
    [state.a, state.b] = [state.b, state.a];
    syncControls();
    setURL();
    renderCompare();
  });
  document.addEventListener("click", (e) => {
    const dot = e.target.closest(".fraction-dot");
    if (dot) toast(dot.getAttribute("aria-label"));
    const t = e.target.closest("[data-transcript]");
    if (t)
      openTranscript(
        t.dataset.transcript,
        t.dataset.position ? Number(t.dataset.position) : null,
      );
    const p = e.target.closest("[data-pair-a]");
    if (p)
      location.hash = `compare?a=${p.dataset.pairA}&b=${p.dataset.pairB}&reads=${state.coverage}`;
    const pos = e.target.closest("[data-select-position]");
    if (pos) {
      state.position = Number(pos.dataset.selectPosition);
      document
        .querySelectorAll(".site-cell")
        .forEach((b) =>
          b.classList.toggle(
            "selected",
            Number(b.dataset.selectPosition) === state.position,
          ),
        );
      renderPosition();
    }
  });
  for (const id of [
    "candidate-search",
    "candidate-direction",
    "candidate-sort",
  ])
    $(id).addEventListener(
      id === "candidate-search" ? "input" : "change",
      () => {
        candidatePage = 0;
        renderCandidates();
      },
    );
  $("candidate-prev").onclick = () => {
    candidatePage--;
    renderCandidates();
  };
  $("candidate-next").onclick = () => {
    candidatePage++;
    renderCandidates();
  };
  $("position-prev").onclick = () => {
    positionPage--;
    state.position = null;
    renderTranscript();
  };
  $("position-next").onclick = () => {
    positionPage++;
    state.position = null;
    renderTranscript();
  };
  $("transcript-form").addEventListener("submit", (e) => {
    e.preventDefault();
    state.position = null;
    loadTranscript($("transcript-id").value);
  });
  $("transcript-id").addEventListener("input", () => {
    const q = $("transcript-id").value.trim().toUpperCase();
    $("transcript-suggestions").innerHTML =
      q.length >= 5
        ? Object.keys(index)
            .filter((t) => t.startsWith(q))
            .slice(0, 12)
            .map((t) => `<option value="${esc(t)}"></option>`)
            .join("")
        : "";
  });
  $("transcript-cell").onchange = () => {
    positionPage = 0;
    state.position = null;
    renderTranscript();
  };
  $("download-overview").onclick = () =>
    csv(
      `sgnex-runs-min${state.coverage}-threshold${dec(state.threshold, 4)}.csv`,
      visibleOverviewRows.map((r) => ({
        ...r,
        model_sha256: study.model_sha256,
      })),
    );
  $("download-candidates").onclick = () =>
    csv(
      `${state.a}-vs-${state.b}-candidates-min${state.coverage}.csv`,
      filteredCandidates.map((r) => ({
        cell_a: state.a,
        cell_b: state.b,
        min_reads: state.coverage,
        threshold: study.threshold,
        ...r,
        model_sha256: study.model_sha256,
      })),
    );
  $("download-transcript").onclick = () =>
    csv(
      `${transcript.id}-min${state.coverage}.csv`,
      currentTranscriptRows().map((r) => ({
        transcript_id: transcript.id,
        original_transcript_id: r[5] || transcript.id,
        transcript_position: r[1],
        sequence: r[2],
        ...study.samples[r[0]],
        n_reads: r[3],
        score: r[4],
        prediction: Number(r[4] >= study.threshold),
        min_reads: state.coverage,
        threshold: study.threshold,
        model_sha256: study.model_sha256,
      })),
    );
  $("download-provenance").onclick = () =>
    download(
      "sgnex-study-metadata.json",
      JSON.stringify(study, null, 2),
      "application/json",
    );
  $("scatter").onmousemove = (e) => {
    const p = nearestPoint(e),
      tip = $("scatter-tooltip");
    tip.hidden = !p;
    if (p) {
      tip.innerHTML = `<b>${esc(p.transcript_id)} · ${p.transcript_position}</b><br>${esc(state.a)} ${dec(p.mean_score_a)} · ${esc(state.b)} ${dec(p.mean_score_b)}<br>Click to explore this transcript`;
      tip.style.left =
        Math.max(0, Math.min(p.x + 10, $("scatter").clientWidth - 235)) + "px";
      tip.style.top = Math.max(0, p.y - 70) + "px";
    }
  };
  $("scatter").onmouseleave = () => ($("scatter-tooltip").hidden = true);
  $("scatter").onclick = (e) => {
    const p = nearestPoint(e);
    if (p) openTranscript(p.transcript_id, p.transcript_position);
  };
  new ResizeObserver(() => drawScatter()).observe($("scatter"));
}
async function start() {
  try {
    [study, index, aliases] = await Promise.all([
      json("data/summary.json"),
      json("data/transcript_index.json"),
      json("data/transcript_aliases.json"),
    ]);
    state.threshold = study.threshold;
    state.a = study.cell_lines.includes("A549") ? "A549" : study.cell_lines[0];
    state.b = study.cell_lines.includes("HepG2")
      ? "HepG2"
      : study.cell_lines[1];
    $("dataset-badge").textContent =
      `${study.samples.length} / ${study.expected_samples} runs · ${study.missing_samples.length ? "Partial collection" : "Processed collection complete"}`;
    $("cell-a").innerHTML = options(study.cell_lines);
    $("cell-b").innerHTML = options(study.cell_lines);
    $("transcript-cell").innerHTML =
      '<option value="all">All cell lines</option>' + options(study.cell_lines);
    $("overview-threshold").innerHTML = [0.15, study.threshold, 0.3, 0.5]
      .map(
        (t) =>
          `<option value="${t}">${dec(t, 4)}${t === study.threshold ? " · saved model threshold" : ""}</option>`,
      )
      .join("");
    document
      .querySelectorAll(".saved-threshold")
      .forEach((e) => (e.textContent = dec(study.threshold, 4)));
    renderMethods();
    wire();
    $("loading").hidden = true;
    $("app").hidden = false;
    await route();
  } catch (e) {
    $("loading").hidden = true;
    $("load-error").hidden = false;
    $("load-error").textContent =
      `${e.message} Serve this folder through Docker or a local HTTP server, and ensure the bundled data directory is present.`;
  }
}
start();
