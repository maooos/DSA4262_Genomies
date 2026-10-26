import {
  observedExtent,
  clampWindow,
  zoomWindow,
  panWindow,
  nearestRow,
  pageForPosition,
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
} from "./core.mjs?v=dashboard-overview-4";
const $ = (id) => document.getElementById(id);
let transcriptTracks;
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
    selectedTranscriptCells(),
  );
}
function selectedTranscriptCells() {
  return [...document.querySelectorAll(".transcript-cell-choice:checked")].map((input) => input.value);
}
function syncTranscriptCellFilter() {
  const selected = selectedTranscriptCells();
  const all = $("transcript-cell-all");
  all.checked = selected.length === study.cell_lines.length;
  all.indeterminate = selected.length > 0 && !all.checked;
  $("transcript-cell-summary").textContent = all.checked ? "All cell lines" : selected.length ? selected.join(", ") : "No cell lines selected";
  $("transcript-cell-summary").setAttribute("aria-label", `Choose cell lines to show: ${$("transcript-cell-summary").textContent}`);
}
function selectTranscriptPosition(position) {
  if (!transcript) return;
  const positions = [...new Set(currentTranscriptRows().map((r) => r[1]))].sort((a, b) => a - b);
  if (!positions.includes(position)) return;
  state.position = position;
  positionPage = pageForPosition(positions, position);
  renderTranscript();
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
        selectedTranscriptCells().includes(s.cell_line),
    );
  const map = new Map(rows.map((r) => [`${r[0]}|${r[1]}`, r]));
  $("transcript-heatmap").innerHTML = positions.length
    ? `<table><thead><tr><th scope="col">Cell line · replicate / run</th>${displayed.map((p) => `<th scope="col"><button class="position-head ${p === state.position ? "selected" : ""}" data-select-position="${p}" aria-pressed="${p === state.position}" aria-label="Inspect position ${p}">${p}</button></th>`).join("")}</tr></thead><tbody>${sampleIds
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
  try {
    transcriptTracks.update({
      id: transcript.id, rows, allRows: transcript.rows, samples: study.samples,
      cells: selectedTranscriptCells(), position: state.position,
      threshold: study.threshold, coverage: state.coverage,
    });
    $("track-load-message").hidden = true;
  } catch (error) {
    console.error("Transcript track rendering failed:", error);
    $("track-load-message").hidden = false;
    $("track-load-message").className = "error";
    $("track-load-message").textContent = `Transcript tracks could not be drawn: ${error.message}. The heatmap and table are still available. Rebuild the dashboard and hard-refresh to check for mismatched files.`;
  }
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
  transcriptTracks = new TranscriptTracks(selectTranscriptPosition);
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
    document.querySelectorAll(".transcript-cell-choice").forEach((input) => { input.checked = true; });
    syncTranscriptCellFilter();
    $("transcript-id").value = "";
    transcript = null;
    transcriptTracks.id = null;
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
      selectTranscriptPosition(Number(pos.dataset.selectPosition));
      if (e.detail === 0)
        $("transcript-heatmap").querySelector(`.position-head[data-select-position="${state.position}"]`)?.focus({ preventScroll: true });
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
  $("transcript-cell").onchange = (event) => {
    if (event.target.id === "transcript-cell-all")
      document.querySelectorAll(".transcript-cell-choice").forEach((input) => { input.checked = event.target.checked; });
    syncTranscriptCellFilter();
    positionPage = 0;
    state.position = null;
    renderTranscript();
  };
  document.addEventListener("click", (event) => {
    if (!$("transcript-cell-filter").contains(event.target)) $("transcript-cell-filter").open = false;
  });
  $("transcript-cell-filter").addEventListener("keydown", (event) => {
    if (event.key === "Escape") {
      $("transcript-cell-filter").open = false;
      $("transcript-cell-summary").focus();
    }
  });
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
      '<label class="transcript-cell-option"><input id="transcript-cell-all" type="checkbox" checked> All cell lines</label>' +
      study.cell_lines.map((cell) => `<label class="transcript-cell-option"><input class="transcript-cell-choice" type="checkbox" value="${esc(cell)}" checked> ${esc(cell)}</label>`).join("");
    syncTranscriptCellFilter();
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

// Transcript explorer: zoomable run-level score and read-support tracks.
const transcriptTrackColours = ['#17634e','#466fa1','#b87538','#806394','#9b5570','#427d80','#788547'];

class TranscriptTracks {
  constructor(onSelect) {
    this.onSelect = onSelect;
    this.id = null; this.domain = null; this.window = null;
    this.rows = []; this.runs = []; this.position = null; this.drag = null;
    this.canvas = $('transcript-tracks-canvas');
    this.overview = $('transcript-tracks-overview');
    this.tip = $('transcript-tracks-tooltip');
    document.querySelectorAll('[data-track-zoom]').forEach(button => button.addEventListener('click',()=>{
      if(this.domain) this.setWindow(zoomWindow(this.window,Number(button.dataset.trackZoom),this.domain));
    }));
    document.querySelectorAll('[data-track-pan]').forEach(button => button.addEventListener('click',()=>{
      if(this.domain) this.setWindow(panWindow(this.window,Number(button.dataset.trackPan),this.domain));
    }));
    $('track-full-span').addEventListener('click',()=>{if(this.domain) this.setWindow([...this.domain]);});
    $('track-range-form').addEventListener('submit',event=>{
      event.preventDefault();
      const rawStart=$('track-start').value.trim(),rawEnd=$('track-end').value.trim();
      const start=Number(rawStart),end=Number(rawEnd);
      const valid=/^\d+$/.test(rawStart)&&/^\d+$/.test(rawEnd)&&Number.isSafeInteger(start)&&Number.isSafeInteger(end)&&start<end;
      $('track-range-error').hidden=valid;
      $('track-start').setAttribute('aria-invalid',String(!valid));
      $('track-end').setAttribute('aria-invalid',String(!valid));
      if(valid && this.domain) this.setWindow(clampWindow(start,end,this.domain));
    });
    this.canvas.addEventListener('pointerdown',event=>{
      if(!this.window||event.button!==0) return;
      const p=this.local(event);
      if(p.x<this.left||p.x>this.width-this.right) return;
      this.drag={x:p.x,y:p.y,start:[...this.window],moved:false};
      this.canvas.setPointerCapture(event.pointerId);
    });
    this.canvas.addEventListener('pointermove',event=>{
      if(this.drag){
        const p=this.local(event),delta=p.x-this.drag.x;
        if(Math.abs(delta)>5) this.drag.moved=true;
        if(this.drag.moved){
          const shift=-delta/this.plotWidth*(this.drag.start[1]-this.drag.start[0]);
          this.setWindow(clampWindow(this.drag.start[0]+shift,this.drag.start[1]+shift,this.domain),false);
        }
        this.hideTip(); return;
      }
      this.hover(event);
    });
    this.canvas.addEventListener('pointerup',event=>{
      if(!this.drag) return;
      const moved=this.drag.moved;this.drag=null;
      if(this.canvas.hasPointerCapture(event.pointerId)) this.canvas.releasePointerCapture(event.pointerId);
      if(!moved){const row=this.hit(event);if(row) this.onSelect(row[1]);}
      this.describe();
    });
    this.canvas.addEventListener('pointercancel',()=>{this.drag=null;});
    this.canvas.addEventListener('lostpointercapture',()=>{this.drag=null;});
    this.canvas.addEventListener('pointerleave',()=>this.hideTip());
    this.canvas.addEventListener('keydown',event=>{
      if(!this.rows.length) return;
      if(event.key==='+'||event.key==='='){event.preventDefault();this.setWindow(zoomWindow(this.window,2,this.domain));return;}
      if(event.key==='-'){event.preventDefault();this.setWindow(zoomWindow(this.window,.5,this.domain));return;}
      if(!['ArrowLeft','ArrowRight','Home','End'].includes(event.key)) return;
      event.preventDefault();
      const positions=this.positions,index=positions.indexOf(this.position);
      const next=event.key==='Home'?positions[0]:event.key==='End'?positions.at(-1):positions[Math.max(0,Math.min(positions.length-1,index+(event.key==='ArrowLeft'?-1:1)))];
      this.onSelect(next);
    });
    this.canvas.addEventListener('wheel',event=>{
      if(!event.ctrlKey||!this.window) return;
      const p=this.local(event);if(p.x<this.left||p.x>this.width-this.right)return;
      event.preventDefault();
      const anchor=this.toPosition(p.x);
      this.setWindow(zoomWindow(this.window,event.deltaY<0?1.5:1/1.5,this.domain,anchor));
    },{passive:false});
    this.overview.addEventListener('pointerdown',event=>{
      if(event.button!==0||!this.domain||!this.positions.length)return;
      const x=this.overviewX(event);
      if(x<this.left||x>this.width-this.right)return;
      this.overviewDrag={x,start:[...this.window],moved:false,canDrag:this.overviewWindowHit(x)};
      this.overview.setPointerCapture(event.pointerId);
    });
    this.overview.addEventListener('pointermove',event=>{
      const x=this.overviewX(event);
      if(this.overviewDrag){
        const delta=x-this.overviewDrag.x;
        if(Math.abs(delta)>5)this.overviewDrag.moved=true;
        if(this.overviewDrag.moved&&this.overviewDrag.canDrag){
          const shift=delta/this.plotWidth*(this.domain[1]-this.domain[0]);
          this.setWindow(clampWindow(this.overviewDrag.start[0]+shift,this.overviewDrag.start[1]+shift,this.domain));
          this.overview.style.cursor='grabbing';
        }
      } else this.overview.style.cursor=this.overviewWindowHit(x)?'grab':'pointer';
    });
    this.overview.addEventListener('pointerup',event=>{
      if(!this.overviewDrag)return;
      const moved=this.overviewDrag.moved;this.overviewDrag=null;
      if(this.overview.hasPointerCapture(event.pointerId))this.overview.releasePointerCapture(event.pointerId);
      this.overview.style.cursor='pointer';
      if(!moved){
        const position=this.domain[0]+(this.overviewX(event)-this.left)/this.plotWidth*(this.domain[1]-this.domain[0]);
        const row=nearestRow(this.overviewPositions,position);
        if(row)this.onSelect(row[1]);
      }
    });
    for(const type of ['pointercancel','lostpointercapture'])this.overview.addEventListener(type,()=>{
      this.overviewDrag=null;this.overview.style.cursor='pointer';
    });
    this.overview.addEventListener('keydown',event=>{
      if(!this.positions.length||!['ArrowLeft','ArrowRight','Home','End'].includes(event.key))return;
      event.preventDefault();
      const index=this.positions.indexOf(this.position);
      const next=event.key==='Home'?this.positions[0]:event.key==='End'?this.positions.at(-1):this.positions[Math.max(0,Math.min(this.positions.length-1,index+(event.key==='ArrowLeft'?-1:1)))];
      this.onSelect(next);
    });
    $('transcript-run-tracks').addEventListener('toggle',()=>{
      const open=$('transcript-run-tracks').open;
      $('track-toggle-label').textContent=open?'Collapse tracks':'Expand tracks';
      this.overviewDrag=null;this.overview.style.cursor='pointer';
      this.hideTip();
        try { this.draw();this.describe(); }
        catch(error){
          console.error('Transcript track rendering failed:',error);
          $('track-load-message').hidden=false;
          $('track-load-message').className='error';
          $('track-load-message').textContent=`Transcript tracks could not be drawn: ${error.message}. The heatmap and table are still available.`;
        }
    });
    this.observer=new ResizeObserver(()=>this.draw());
    this.observer.observe($('transcript-tracks'));
  }

  update({id,rows,allRows,samples,cells,position,threshold,coverage}) {
    const changed=this.id!==id;
    this.id=id;this.rows=rows;this.threshold=threshold;this.coverage=coverage;
    this.domain=observedExtent(allRows);
    if(changed||!this.window) this.window=this.domain?[...this.domain]:null;
    else if(this.domain) this.window=clampWindow(...this.window,this.domain);
    this.positions=[...new Set(rows.map(r=>r[1]))].sort((a,b)=>a-b);
    this.overviewPositions=this.positions.map(position=>[null,position]);
    this.runs=samples.map((sample,index)=>({sample,index,rows:[]})).filter(run=>cells.includes(run.sample.cell_line));
    const byIndex=new Map(this.runs.map(run=>[run.index,run]));
    for(const row of rows) byIndex.get(row[0])?.rows.push(row);
    for(const run of this.runs) run.rows.sort((a,b)=>a[1]-b[1]);
    const paletteCells=[...new Set(samples.map(s=>s.cell_line))];
    for(const run of this.runs)run.colour=transcriptTrackColours[paletteCells.indexOf(run.sample.cell_line)%transcriptTrackColours.length];
    this.maxReads=rows.reduce((max,r)=>Math.max(max,r[3]),1);
    this.position=position;
    if(!changed&&$('transcript-run-tracks').open) this.ensureVisible(position);
    this.hideTip();this.draw();this.describe();
  }

  setSelection(position) {
    if(this.position===position)return;
    this.position=position;this.ensureVisible(position);this.hideTip();this.draw();this.describe();
  }

  ensureVisible(position) {
    if(position===null||!this.window||!this.domain)return;
    if(position<this.window[0]||position>this.window[1]){
      const half=(this.window[1]-this.window[0])/2;
      this.window=clampWindow(position-half,position+half,this.domain);
    }
  }

  setWindow(window,announce=true) {
    this.window=window;
    $('track-range-error').hidden=true;
    $('track-start').setAttribute('aria-invalid','false');
    $('track-end').setAttribute('aria-invalid','false');
    this.hideTip();this.draw();if(announce)this.describe();
  }

  describe() {
    if(!this.window)return;
    $('track-start').value=this.window[0];$('track-end').value=this.window[1];
    const visible=this.positions.filter(p=>p>=this.window[0]&&p<=this.window[1]).length;
    $('track-window-status').textContent=`Viewing ${fmt(this.window[0])}–${fmt(this.window[1])} · ${fmt(visible)} eligible positions in view${this.position===null?'':` · selected ${fmt(this.position)}`}`;
    $('track-overview-summary').textContent=`${fmt(this.positions.length)} eligible positions · selected position: ${this.position===null?'none':fmt(this.position)}`;
    this.overview.setAttribute('aria-describedby','track-overview-summary track-overview-instruction track-overview-explanation');
    const atFull=this.window[0]===this.domain[0]&&this.window[1]===this.domain[1];
    $('track-full-span').disabled=atFull||!this.rows.length;
    for(const b of document.querySelectorAll('[data-track-zoom]'))b.disabled=!this.rows.length||(Number(b.dataset.trackZoom)>1?this.window[1]-this.window[0]<=1:atFull);
    for(const b of document.querySelectorAll('[data-track-pan]'))b.disabled=!this.rows.length||(Number(b.dataset.trackPan)<0?this.window[0]===this.domain[0]:this.window[1]===this.domain[1]);
    for(const el of [$('track-start'),$('track-end'),$('track-range-apply')])el.disabled=!this.rows.length;
    this.canvas.setAttribute('aria-label',`${this.id}: scores and read support per sequencing run. ${$('track-window-status').textContent}. Arrow keys select sites; plus and minus zoom.`);
  }

  context(canvas,width,height) {
    const dpr=window.devicePixelRatio||1;
    canvas.width=Math.round(width*dpr);canvas.height=Math.round(height*dpr);
    canvas.style.height=`${height}px`;
    const ctx=canvas.getContext('2d');ctx.scale(dpr,dpr);return ctx;
  }

  draw() {
    if(!this.domain||!this.window)return;
    this.width=$('transcript-tracks').clientWidth;
    if(!this.width)return;
    const header='Cell line · replicate / run';
    const measure=this.canvas.getContext('2d');measure.font='600 10px system-ui';
    this.left=Math.ceil(measure.measureText(header).width)+18;
    this.right=16;this.plotWidth=Math.max(1,this.width-this.left-this.right);
    this.rowHeight=108;this.top=40;
    this.drawOverview();
    if(!$('transcript-run-tracks').open)return;
    const height=this.top+this.runs.length*this.rowHeight+28;
    const ctx=this.context(this.canvas,this.width,height);
    ctx.fillStyle='#fff';ctx.fillRect(0,0,this.width,height);
    for(let i=0;i<this.runs.length;i++){
      ctx.fillStyle='#f3f7f2';ctx.fillRect(this.left,this.top+i*this.rowHeight+66,this.plotWidth,34);
    }
    ctx.font='10px system-ui';ctx.fillStyle='#66746e';
    ctx.fillText('Scores 0–1; reads 0–'+fmt(this.maxReads)+' (log scale)',this.left,13);
    ctx.font='600 10px system-ui';ctx.fillStyle='#20372f';
    ctx.fillText(header,1,13);
    ctx.font='10px system-ui';ctx.fillStyle='#66746e';
    const x=p=>this.left+(p-this.window[0])/(this.window[1]-this.window[0])*this.plotWidth;
    const ticks=this.width<500?3:5,seenTicks=new Set();
    for(let i=0;i<=ticks;i++){
      const position=Math.round(this.window[0]+i/ticks*(this.window[1]-this.window[0]));
      if(seenTicks.has(position))continue;seenTicks.add(position);
      ctx.textAlign=i===0?'left':i===ticks?'right':'center';ctx.fillText(fmt(position),x(position),29);
      ctx.strokeStyle='#edf0e9';ctx.beginPath();ctx.moveTo(x(position),this.top);ctx.lineTo(x(position),height-20);ctx.stroke();
    }
    const emptyRuns=[];
    this.runs.forEach((run,i)=>{
      const y=this.top+i*this.rowHeight,scoreTop=y+25,scoreBottom=y+58,readTop=y+71,readBottom=y+96;
      ctx.textAlign='left';ctx.fillStyle=run.colour;ctx.font='600 11px system-ui';
      ctx.fillText(`${run.sample.cell_line} · ${run.sample.replicate} / ${run.sample.run}`,1,y+13);
      ctx.font='600 10px system-ui';ctx.fillStyle='#52645c';ctx.textAlign='right';
      ctx.fillText('Score',this.left-10,y+45);ctx.fillText('Reads',this.left-10,y+87);ctx.textAlign='left';
      ctx.strokeStyle='#dfe5dc';ctx.beginPath();ctx.moveTo(this.left,scoreBottom);ctx.lineTo(this.width-this.right,scoreBottom);ctx.moveTo(this.left,readBottom);ctx.lineTo(this.width-this.right,readBottom);ctx.stroke();
      ctx.save();ctx.beginPath();ctx.rect(this.left,y,this.plotWidth,this.rowHeight);ctx.clip();
      ctx.setLineDash([3,3]);ctx.strokeStyle='#b87538';ctx.beginPath();ctx.moveTo(this.left,scoreBottom-this.threshold*(scoreBottom-scoreTop));ctx.lineTo(this.width-this.right,scoreBottom-this.threshold*(scoreBottom-scoreTop));ctx.stroke();ctx.setLineDash([]);
      let visible=0;
      for(const row of run.rows){
        if(row[1]<this.window[0]||row[1]>this.window[1])continue;visible++;
        const px=x(row[1]),sy=scoreBottom-row[4]*(scoreBottom-scoreTop),ry=readBottom-Math.log1p(row[3])/Math.log1p(this.maxReads)*(readBottom-readTop);
        ctx.strokeStyle=run.colour;ctx.lineWidth=1.5;ctx.beginPath();ctx.moveTo(px,scoreBottom);ctx.lineTo(px,sy);ctx.stroke();
        ctx.fillStyle=run.colour;ctx.beginPath();ctx.arc(px,sy,2.2,0,Math.PI*2);ctx.fill();
        ctx.globalAlpha=.7;ctx.fillRect(px-1.5,ry,3,Math.max(1,readBottom-ry));ctx.globalAlpha=1;
      }
      if(!visible)emptyRuns.push(y);
      ctx.restore();ctx.lineWidth=1;
      ctx.strokeStyle='#a7b8ac';ctx.lineWidth=2;
      ctx.beginPath();ctx.moveTo(0,y+this.rowHeight-3);ctx.lineTo(this.width,y+this.rowHeight-3);ctx.stroke();
      ctx.lineWidth=1;
    });
    if(this.position!==null&&this.position>=this.window[0]&&this.position<=this.window[1]){
      ctx.strokeStyle='#20372f';ctx.lineWidth=1.5;ctx.beginPath();ctx.moveTo(x(this.position),this.top);ctx.lineTo(x(this.position),height-20);ctx.stroke();
    }
    // Empty-state text sits in the run heading, clear of both plot baselines.
    // Draw it last on an opaque backing so selection/grid lines cannot cross it.
    ctx.font='10px system-ui';ctx.textAlign='left';
    for(const y of emptyRuns){
      const lines=this.plotWidth<190?['No qualifying sites','in this window']:['No qualifying sites in this window'];
      const textWidth=Math.max(...lines.map(line=>ctx.measureText(line).width));
      ctx.fillStyle='#fff';ctx.fillRect(this.left+4,y+1,textWidth+10,lines.length*11+5);
      ctx.fillStyle='#66746e';lines.forEach((line,i)=>ctx.fillText(line,this.left+8,y+12+i*11));
    }
    ctx.font='10px system-ui';ctx.fillStyle='#66746e';ctx.textAlign='center';ctx.fillText('Transcript position',this.left+this.plotWidth/2,height-3);
  }

  drawOverview() {
    const axisHeight=22;
    const overview=this.context(this.overview,this.width,axisHeight+58);
    overview.fillStyle='#f5f7f1';overview.fillRect(this.left,axisHeight,this.plotWidth,58);
    const fullX=p=>this.left+(p-this.domain[0])/(this.domain[1]-this.domain[0])*this.plotWidth;
    // Coordinate ticks describe the full observed span, independent of track zoom.
    overview.font='10px system-ui';overview.fillStyle='#66746e';
    const labelWidth=Math.max(overview.measureText(fmt(this.domain[0])).width,overview.measureText(fmt(this.domain[1])).width);
    const intervals=Math.max(1,Math.min(5,Math.floor(this.plotWidth/(labelWidth+24))));
    const seen=new Set();
    for(let i=0;i<=intervals;i++){
      const position=Math.round(this.domain[0]+i/intervals*(this.domain[1]-this.domain[0]));
      if(seen.has(position))continue;seen.add(position);
      const px=fullX(position);
      overview.textAlign=position===this.domain[0]?'left':position===this.domain[1]?'right':'center';
      overview.fillText(fmt(position),px,13);
      overview.strokeStyle='#bbc8ba';overview.beginPath();overview.moveTo(px,17);overview.lineTo(px,axisHeight);overview.stroke();
    }
    overview.textAlign='left';
    // Per-pixel maxima retain high-score sites in a dense full-span overview.
    const peaks=new Map();
    for(const row of this.rows){const pixel=Math.round(fullX(row[1]));peaks.set(pixel,Math.max(peaks.get(pixel)||0,row[4]));}
    overview.strokeStyle='#95b29d';
    for(const [pixel,score]of peaks){overview.beginPath();overview.moveTo(pixel,axisHeight+48);overview.lineTo(pixel,axisHeight+48-score*34);overview.stroke();}
    if($('transcript-run-tracks').open){
      overview.fillStyle='#17634e18';overview.fillRect(fullX(this.window[0]),axisHeight+4,Math.max(2,fullX(this.window[1])-fullX(this.window[0])),48);
      overview.strokeStyle='#17634e';overview.strokeRect(fullX(this.window[0]),axisHeight+4,Math.max(2,fullX(this.window[1])-fullX(this.window[0])),48);
    }
    overview.font='10px system-ui';overview.fillStyle='#66746e';overview.fillText('Observed span',1,axisHeight+22);
    if(this.position!==null){overview.strokeStyle='#20372f';overview.beginPath();overview.moveTo(fullX(this.position),axisHeight+4);overview.lineTo(fullX(this.position),axisHeight+52);overview.stroke();}
    $('track-threshold-value').textContent=this.threshold.toFixed(4);
    $('track-threshold-value').title=`Saved H04 threshold: ${this.threshold}`;
    $('track-empty').hidden=!!this.rows.length;
    $('track-scale-note').textContent=`Both panels use the same coordinates across runs. Score scale: 0–1; dashed line: ${this.threshold.toFixed(4)}. Read support: 0–${fmt(this.maxReads)}, log(1 + reads), shared across all shown runs. Observed span: ${fmt(this.domain[0])}–${fmt(this.domain[1])}; annotated transcript length is unavailable.`;
  }

  overviewX(event){return event.clientX-this.overview.getBoundingClientRect().left;}
  overviewWindowHit(x){
    if(!$('transcript-run-tracks').open||!this.domain||!this.window||!this.positions.length)return false;
    const scale=this.plotWidth/(this.domain[1]-this.domain[0]);
    const selectedX=this.position===null?Infinity:this.left+(this.position-this.domain[0])*scale;
    // The selected-position line is an indicator, never a drag handle.
    if(Math.abs(x-selectedX)<=5)return false;
    return x>=this.left+(this.window[0]-this.domain[0])*scale&&x<=this.left+(this.window[1]-this.domain[0])*scale;
  }
  local(event){const rect=this.canvas.getBoundingClientRect();return{x:event.clientX-rect.left,y:event.clientY-rect.top};}
  toPosition(x){return this.window[0]+(x-this.left)/this.plotWidth*(this.window[1]-this.window[0]);}
  hit(event){
    if(!this.window)return null;
    const p=this.local(event),index=Math.floor((p.y-this.top)/this.rowHeight);
    if(p.x<this.left||p.x>this.width-this.right||index<0||index>=this.runs.length)return null;
    const row=nearestRow(this.runs[index].rows,this.toPosition(p.x),9/this.plotWidth*(this.window[1]-this.window[0]));
    return row&&row[1]>=this.window[0]&&row[1]<=this.window[1]?row:null;
  }
  hover(event){
    const row=this.hit(event);if(!row){this.hideTip();return;}
    const p=this.local(event),run=this.runs.find(r=>r.index===row[0]);
    this.tip.textContent=`${run.sample.cell_line} · replicate ${run.sample.replicate} / run ${run.sample.run}\nPosition ${fmt(row[1])} · ${row[2]}\nScore ${row[4]} · ${fmt(row[3])} reads`;
    this.tip.hidden=false;
    this.tip.style.left=`${Math.max(0,Math.min(p.x+12,this.width-250))}px`;
    this.tip.style.top=`${Math.max(0,p.y-82)}px`;
    this.canvas.style.cursor='pointer';
  }
  hideTip(){this.tip.hidden=true;this.canvas.style.cursor='grab';}
}

start();

