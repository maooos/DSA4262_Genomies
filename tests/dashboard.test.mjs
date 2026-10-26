import test from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { gunzipSync } from "node:zlib";
import {
  observedExtent,
  clampWindow,
  zoomWindow,
  panWindow,
  nearestRow,
  pageForPosition,
  parseReadSupport,
  overviewRows,
  orientPair,
  filterCandidates,
  transcriptRows,
  positionGroups,
  csvText,
} from "../task2_site/dist/core.mjs";
const data = new URL("../task2_site/dist/data/", import.meta.url);
const read = (name) => JSON.parse(readFileSync(new URL(name, data)));

test("custom read support accepts exact positive integers, not invalid or rounded values", () => {
  for (const n of [1, 7, 25, 37, 100, 100000, Number.MAX_SAFE_INTEGER])
    assert.equal(parseReadSupport(String(n)), n);
  for (const value of ["", " ", "0", "-1", "2.5", "NaN", "Infinity", "1e3", "9007199254740992", null])
    assert.equal(parseReadSupport(value), null);
});

test("custom transcript cutoffs filter exact boundaries, including below ten and empty results", () => {
  const samples = [{cell_line: "A", replicate: 1}];
  const rows = [1, 24, 25, 26, 100].map((n) => [0, n, "AAAAAAA", n, 0.5]);
  assert.equal(transcriptRows(rows, samples, 1).length, 5);
  assert.deepEqual(transcriptRows(rows, samples, 25).map((r) => r[3]), [25, 26, 100]);
  assert.deepEqual(transcriptRows(rows, samples, 101), []);
});

test("overview coverage filters match the complete prediction audit", () => {
  const study = read("summary.json");
  const rows = overviewRows(study, 20, study.threshold);
  assert.equal(rows.length, 23);
  assert.equal(
    rows.reduce((n, r) => n + r.eligible_sites, 0),
    2983948,
  );
  assert.equal(
    overviewRows(study, 50, study.threshold).reduce(
      (n, r) => n + r.eligible_sites,
      0,
    ),
    1101143,
  );
});
test("reversing a comparison preserves denominators and swaps all directional quantities", () => {
  const original = read("pairs/A549__HepG2_min20.json");
  const reversed = orientPair(original, "HepG2");
  assert.equal(reversed.summary.cell_a, "HepG2");
  assert.equal(
    reversed.summary.positive_a_only,
    original.summary.positive_b_only,
  );
  assert.equal(reversed.summary.jaccard, original.summary.jaccard);
  assert.equal(
    reversed.candidates[0].mean_score_a,
    original.candidates[0].mean_score_b,
  );
  assert.equal(
    reversed.candidates[0].delta_score,
    -original.candidates[0].delta_score,
  );
  assert.deepEqual(orientPair(reversed, "A549"), original);
});
test("candidate search, direction and ranking use the exported candidate population", () => {
  const rows = read("pairs/A549__HepG2_min20.json").candidates;
  assert.equal(
    filterCandidates(rows, "enst00000451270", "a")[0].transcript_position,
    1212,
  );
  assert.ok(filterCandidates(rows, "", "b").every((r) => r.delta_score < 0));
  assert.equal(filterCandidates(rows, "DOES_NOT_EXIST").length, 0);
});
test("transcript aggregation gives equal weight to replicate groups, not technical runs", () => {
  const samples = [
    { cell_line: "A", replicate: 1 },
    { cell_line: "A", replicate: 1 },
    { cell_line: "A", replicate: 2 },
    { cell_line: "B", replicate: 1 },
  ];
  const rows = [
    [0, 10, "AAAAAAA", 30, 0.9],
    [1, 10, "AAAAAAA", 30, 0.9],
    [2, 10, "AAAAAAA", 25, 0.1],
    [3, 10, "AAAAAAA", 2, 0.8],
  ];
  const eligible = transcriptRows(rows, samples, 20);
  const groups = positionGroups(eligible, samples, 10, 0.5);
  assert.equal(groups.length, 1); // Low coverage B is missing, not a negative.
  assert.equal(groups[0].score, 0.5);
  assert.equal(groups[0].replicates, 2);
  assert.equal(groups[0].positiveReplicates, 1);
});
test("different sequence contexts are never averaged together", () => {
  const groups = positionGroups(
    [
      [0, 1, "AAAAAAA", 20, 0.9],
      [1, 1, "CCCCCCC", 20, 0.1],
    ],
    [
      { cell_line: "A", replicate: 1 },
      { cell_line: "A", replicate: 2 },
    ],
    1,
    0.5,
  );
  assert.equal(groups.length, 2);
  assert.ok(groups.every((g) => g.replicates === 1));
});
test("real candidate transcript agrees with exported pair scores after aggregation", () => {
  const study = read("summary.json");
  const index = read("transcript_index.json"),
    id = "ENST00000451270";
  const rows = JSON.parse(
    gunzipSync(readFileSync(new URL(`transcripts/${index[id]}.json.gz`, data))),
  )[id];
  const groups = positionGroups(
    transcriptRows(rows, study.samples, 20),
    study.samples,
    1212,
    study.threshold,
  );
  const candidate = read("pairs/A549__HepG2_min20.json").candidates.find(
    (r) => r.transcript_id === id && r.transcript_position === 1212,
  );
  assert.ok(
    Math.abs(
      groups.find((g) => g.cell === "A549").score - candidate.mean_score_a,
    ) < 1e-9,
  );
  assert.ok(
    Math.abs(
      groups.find((g) => g.cell === "HepG2").score - candidate.mean_score_b,
    ) < 1e-9,
  );
});
test("CSV exports preserve full score precision and escape labels", () => {
  const study = read("summary.json");
  const text = csvText(
    [
      { id: "=formula", score: study.threshold },
      { id: 'quoted,"id', score: 0 },
    ],
    ["id", "score"],
  );
  assert.ok(text.includes('"\'=formula"'));
  assert.ok(text.includes(String(study.threshold)));
  assert.ok(text.includes('quoted,""id'));
});
test("all pair assets match their complete summary, for every coverage setting", () => {
  const study = read("summary.json");
  for (const [coverage, scenario] of Object.entries(study.scenarios))
    for (const stats of scenario.pairs) {
      const pair = read(
        `pairs/${stats.cell_a}__${stats.cell_b}_min${coverage}.json`,
      );
      assert.deepEqual(pair.summary, stats);
      assert.equal(pair.scatter_total, stats.common_sites);
      assert.equal(
        stats.positive_both +
          stats.positive_a_only +
          stats.positive_b_only +
          stats.negative_both,
        stats.common_sites,
      );
      assert.ok(pair.candidates.length <= pair.candidate_limit);
    }
});


// Transcript explorer: multiselect, coordinate navigation and heatmap paging.
test('multiselect combines chosen cell lines, applies coverage, and supports empty selections',()=>{
  const samples=[{cell_line:'A549'},{cell_line:'HepG2'},{cell_line:'K562'}];
  const rows=[[0,100,'AAAAAAA',30,.5],[1,100,'AAAAAAA',10,.2],[2,100,'AAAAAAA',40,.8]];
  assert.deepEqual(transcriptRows(rows,samples,20,['A549','K562']),[rows[0],rows[2]]);
  assert.deepEqual(transcriptRows(rows,samples,20,['HepG2']),[]);
  assert.deepEqual(transcriptRows(rows,samples,1,[]),[]);
  assert.deepEqual(transcriptRows(rows,samples,20,'all'),[rows[0],rows[2]]);
  assert.deepEqual(transcriptRows(rows,samples,20,'K562'),[rows[2]]);
});

test('observed span handles empty, unsorted and single-site data',()=>{
  assert.equal(observedExtent([]),null);
  assert.deepEqual(observedExtent([[0,500],[0,100],[0,300]]),[100,500]);
  assert.deepEqual(observedExtent([[0,100]]),[100,101]);
});
test('zoom and pan preserve span at boundaries and stay within observed coordinates',()=>{
  const domain=[100,1000];
  assert.deepEqual(zoomWindow([100,1000],2,domain),[325,775]);
  assert.deepEqual(panWindow([325,775],1,domain),[550,1000]);
  assert.deepEqual(panWindow([325,775],-1,domain),[100,550]);
  assert.deepEqual(zoomWindow([100,200],.01,domain),domain);
  assert.deepEqual(clampWindow(0,200,domain),[100,300]);
});
test('pointer-anchored zoom retains its relative anchor and minimum one-base window',()=>{
  assert.deepEqual(zoomWindow([100,500],2,[100,1000],200),[150,350]);
  assert.deepEqual(zoomWindow([100,101],10,[100,1000]),[100,101]);
});
test('nearest-site lookup respects sparse coordinates, tolerance and empty runs',()=>{
  const rows=[[0,100],[0,200],[0,1000]];
  assert.equal(nearestRow(rows,950,60),rows[2]);
  assert.equal(nearestRow(rows,600,60),null);
  assert.equal(nearestRow(rows,150),rows[0]);
  assert.equal(nearestRow([],100),null);
});
test('off-page selection finds the correct 24-column heatmap page',()=>{
  const positions=Array.from({length:60},(_,i)=>100+i*137);
  assert.equal(pageForPosition(positions,positions[59]),2);
  assert.equal(pageForPosition(positions,positions[24]),1);
  assert.equal(pageForPosition(positions,-1),0);
});
