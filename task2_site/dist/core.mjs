// Shared, pure data transformations used by the UI and scientific checks.
export function parseReadSupport(value) {
  const text = String(value).trim();
  const number = Number(text);
  return /^\d+$/.test(text) && Number.isSafeInteger(number) && number >= 1
    ? number
    : null;
}
export const mean = (values) =>
  values.length ? values.reduce((a, b) => a + b, 0) / values.length : null;
export const finite = (value) =>
  typeof value === "number" && Number.isFinite(value);
export const escapeHTML = (value) =>
  String(value).replace(
    /[&<>"']/g,
    (c) =>
      ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[
        c
      ],
  );
export function overviewRows(summary, coverage, threshold) {
  return summary.summaries.filter(
    (r) =>
      r.min_reads === coverage && Math.abs(r.threshold - threshold) < 1e-12,
  );
}
export function orientPair(pair, first) {
  if (pair.summary.cell_a === first) return pair;
  const swap = (row) => {
    const out = { ...row };
    for (const key of Object.keys(row)) {
      if (key.endsWith("_a")) out[key] = row[key.slice(0, -2) + "_b"];
      if (key.endsWith("_b")) out[key] = row[key.slice(0, -2) + "_a"];
    }
    if ("delta_score" in row) out.delta_score = -row.delta_score;
    return out;
  };
  const summary = swap(pair.summary);
  summary.positive_a_only = pair.summary.positive_b_only;
  summary.positive_b_only = pair.summary.positive_a_only;
  return {
    ...pair,
    summary,
    candidates: pair.candidates.map(swap),
    scatter: pair.scatter.map(swap),
  };
}
export function filterCandidates(
  rows,
  query = "",
  direction = "all",
  sort = "difference",
) {
  return rows
    .filter(
      (r) =>
        r.transcript_id.toLowerCase().includes(query.trim().toLowerCase()) &&
        (direction === "all" ||
          (direction === "a" ? r.delta_score > 0 : r.delta_score < 0)),
    )
    .sort((a, b) =>
      sort === "position"
        ? a.transcript_position - b.transcript_position
        : sort === "transcript"
          ? a.transcript_id.localeCompare(b.transcript_id) ||
            a.transcript_position - b.transcript_position
          : Math.abs(b.delta_score) - Math.abs(a.delta_score),
    );
}
export function transcriptRows(rows, samples, coverage, cell = "all") {
  return rows.filter(
    (r) =>
      r[3] >= coverage && (cell === "all" || samples[r[0]].cell_line === cell),
  );
}
export function positionGroups(rows, samples, position, threshold) {
  // Keep sequence contexts separate even at the same transcript coordinate.
  const cells = new Map();
  for (const row of rows.filter((r) => r[1] === position)) {
    const sample = samples[row[0]],
      key = sample.cell_line + "|" + row[2];
    if (!cells.has(key))
      cells.set(key, {
        cell: sample.cell_line,
        sequence: row[2],
        replicates: new Map(),
      });
    const group = cells.get(key);
    if (!group.replicates.has(sample.replicate))
      group.replicates.set(sample.replicate, []);
    group.replicates.get(sample.replicate).push(row[4]);
  }
  return [...cells.values()]
    .map((g) => {
      const scores = [...g.replicates.values()].map(mean);
      return {
        cell: g.cell,
        sequence: g.sequence,
        score: mean(scores),
        replicates: scores.length,
        positiveReplicates: scores.filter((s) => s >= threshold).length,
      };
    })
    .sort((a, b) => a.cell.localeCompare(b.cell));
}
export function csvText(rows, columns) {
  const quote = (value) => {
    let str = value == null ? "" : String(value);
    // Keep spreadsheet applications from interpreting text identifiers as formulas.
    if (typeof value === "string" && /^[=+@\-\t\r]/.test(str)) str = "'" + str;
    return '"' + str.replaceAll('"', '""') + '"';
  };
  return [
    columns.map(quote).join(","),
    ...rows.map((r) => columns.map((c) => quote(r[c])).join(",")),
  ].join("\r\n");
}
