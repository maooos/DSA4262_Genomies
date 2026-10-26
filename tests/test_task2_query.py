import gzip
import json
import os
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory
import threading
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlsplit, parse_qs
from urllib.error import HTTPError
from urllib.request import urlopen

import pandas as pd
import pytest

from scripts.task2_serve import DashboardHandler
from src.task2_analysis import (cell_scores, compare_cells, harmonize_transcripts,
                                replicate_agreement, replicate_scores, sample_summary)
from src.task2_query import CoverageQuery, MAX_READS, parse_reads


@pytest.fixture
def study(tmp_path):
    samples = [{"sample_id": f"run{i}", "cell_line": cell, "replicate": rep, "run": run}
               for i, (cell, rep, run) in enumerate([("A", 1, 1), ("A", 1, 2), ("A", 2, 1), ("B", 1, 1), ("B", 2, 1)])]
    rows = []
    for i, sample in enumerate(samples):
        for position, reads in enumerate([1, 24, 25, 37, 100], 1):
            rows.append({**sample, "transcript_id": "ENST000001.2" if i == 4 else "ENST000001",
                         "transcript_position": position, "sequence": "AAAAAAA" if not (i == 4 and position == 4) else "CCCCCCC",
                         "n_reads": reads if i != 1 else reads - 1 if reads > 1 else 1,
                         "score": [0.9, 0.1, 0.7, 0.2, 0.2][i] if position != 3 else 0.5})
    frame = pd.DataFrame(rows)
    canonical, _ = harmonize_transcripts(frame)
    metadata = {"created_at": "fixture", "model_sha256": "fixture-model", "threshold": 0.5,
                "samples": samples, "cell_lines": ["A", "B"]}
    data, assets = tmp_path / "queries", tmp_path / "assets"
    data.mkdir()
    assets.mkdir()
    canonical[["transcript_id_original", "transcript_id"]].drop_duplicates().rename(
        columns={"transcript_id_original": "original", "transcript_id": "canonical"}).to_parquet(data / "identifiers.parquet")
    bundled = []
    for sample in samples:
        subset = frame[frame.sample_id.eq(sample["sample_id"])]
        filename = sample["sample_id"] + ".parquet"
        subset.to_parquet(data / filename, index=False)
        bundled.append({**sample, "file": filename, "site_count": len(subset)})
    (assets / "summary.json").write_text(json.dumps(metadata))
    (data / "manifest.json").write_text(json.dumps({**metadata, "samples": bundled}))
    engine = CoverageQuery(data, assets)
    yield engine, canonical, samples, assets
    engine.close()


@pytest.mark.parametrize("reads", [1, 7, 25, 37, 51, MAX_READS])
def test_arbitrary_queries_match_original_scientific_analysis(study, reads):
    engine, frame, samples, _ = study
    reps = replicate_scores(frame, reads)
    cells = cell_scores(reps, 0.5)
    _, expected = compare_cells(cells, "A", "B", 0.5)
    result = engine.pair(reads, "A", "B")
    for key, value in expected.items():
        assert result["summary"][key] == pytest.approx(value) if isinstance(value, float) else result["summary"][key] == value
    overview = engine.overview(reads)
    assert overview["scenario"]["pairs"][0] == result["summary"]
    for sample in samples:
        for cutoff in [0.15, 0.3, 0.5]:
            expected = sample_summary(frame[frame.sample_id.eq(sample["sample_id"])], sample, reads, cutoff)
            actual = next(r for r in overview["summaries"] if r["sample_id"] == sample["sample_id"] and r["threshold"] == cutoff)
            assert actual == pytest.approx(expected)
    expected_agreement = json.loads(replicate_agreement(reps, 0.5).to_json(orient="records"))
    for expected, actual in zip(expected_agreement, overview["scenario"]["agreement"]):
        assert actual == pytest.approx(expected)


def test_transcripts_include_low_coverage_preserve_ids_and_do_not_round(study):
    engine, frame, _, _ = study
    result = engine.transcript("ENST000001.2")
    assert result["id"] == "ENST000001"
    assert len(result["rows"]) == len(frame)
    assert min(r[3] for r in result["rows"]) == 1
    assert result["rows"][-1][-1] == "ENST000001.2"
    assert engine.transcript("MISSING")["rows"] == []


@pytest.mark.parametrize("value", ["", "0", "-1", "2.5", "NaN", "1e3", MAX_READS + 1])
def test_invalid_coverage_rejected(value):
    with pytest.raises(ValueError, match="positive whole number"):
        parse_reads(value)


def test_http_validation_and_custom_comparison(study):
    engine, _, _, assets = study
    server = ThreadingHTTPServer(("127.0.0.1", 0), partial(DashboardHandler, directory=str(assets)))
    server.engine, server.query_lock = engine, threading.Lock()
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    base = f"http://127.0.0.1:{server.server_port}"
    try:
        with urlopen(base + "/api/pair?reads=25&a=A&b=B") as response:
            result = json.load(response)
            assert result["summary"]["common_sites"] == 2
        for query in ["reads=0", "reads=25.1", "reads=20&reads=50"]:
            with pytest.raises(HTTPError) as caught:
                urlopen(base + "/api/overview?" + query)
            assert caught.value.code == 400
        with urlopen(base + "/api/overview?reads=9007199254740991") as response:
            assert all(r["eligible_sites"] == 0 for r in json.load(response)["summaries"])
    finally:
        server.shutdown()
        server.server_close()
        thread.join()


# Optional end-to-end checks for the dashboard served alongside the query API.
def _expand_dashboard_tracks(page):
    page.wait_for_selector('#transcript-heatmap table')
    if not page.locator('#transcript-run-tracks').evaluate('(el)=>el.open'):
        page.locator('.track-disclosure').click()
    page.wait_for_selector('#transcript-tracks-canvas:visible')
    page.wait_for_function("document.querySelector('#track-toggle-label').textContent === 'Collapse tracks'")


@pytest.mark.skipif(os.environ.get("RUN_BROWSER_TESTS") != "1",
                    reason="Set RUN_BROWSER_TESTS=1 to run the optional Playwright dashboard test.")
def test_dashboard_browser_interactions():
    sync_playwright = pytest.importorskip("playwright.sync_api").sync_playwright
    root=Path(__file__).resolve().parents[1]
    tx='ENST00000451270'
    positions=[100+i*137 for i in range(60)]
    samples=[dict(sample_id=f's{i}',cell_line=cell,replicate=rep,run=run) for i,(cell,rep,run) in enumerate([
        ('A549',1,1),('A549',1,2),('A549',2,1),('HepG2',1,1),('HepG2',2,1),('K562',1,1)])]
    rows=[]
    for i,s in enumerate(samples):
        for j,pos in enumerate(positions):
            if i==4 and j%5==0:continue
            rows.append([i,pos,'AGAACAT',20+(j%5)*20,((j*7+i*3)%95)/100,f'{tx}.1'])
    # API-only observation checks that low coverage follows the existing loading path.
    all_rows=rows+[[0,9000,'AGAACAT',5,.91,f'{tx}.1']]
    threshold=.20839324867970055
    with TemporaryDirectory() as tmp:
        public=Path(tmp)
        shutil.copytree(root/'task2_site/dist',public,dirs_exist_ok=True)
        data=public/'data';(data/'transcripts').mkdir(parents=True)
        (data/'transcript_index.json').write_text(json.dumps({tx:'fixture'}))
        (data/'transcript_aliases.json').write_text(json.dumps({tx+'.1':tx}))
        with gzip.open(data/'transcripts/fixture.json.gz','wt') as f:json.dump({tx:rows},f)
        stats=dict(cell_a='A549',cell_b='HepG2',common_sites=50,positive_both=20,positive_a_only=10,
                   positive_b_only=5,negative_both=15,positive_a=30,positive_b=25,jaccard=20/35,
                   spearman=.7,consistent_candidates=0)
        study=dict(created_at='2026-10-23T00:00:00Z',threshold=threshold,samples=samples,
            cell_lines=['A549','HepG2','K562'],expected_samples=6,missing_samples=[],inventory=[],
            total_scored_sites=len(all_rows),transcript_export_min_reads=10,model_sha256='SYNTHETIC-TEST-ONLY',
            status='synthetic',sequence_conflict_sites=0,identifier_harmonization={},summaries=[],scenarios={})
        for cutoff in (10,20,50):
            study['scenarios'][str(cutoff)]={'pairs':[stats],'agreement':[]}
            for s in samples:
                eligible=[r for r in rows if samples[r[0]]['sample_id']==s['sample_id'] and r[3]>=cutoff]
                for t in (.15,threshold,.3,.5):
                    n=sum(r[4]>=t for r in eligible)
                    study['summaries'].append({**s,'min_reads':cutoff,'threshold':t,'eligible_sites':len(eligible),
                        'positive_sites':n,'positive_fraction':n/len(eligible) if eligible else None})
        (data/'summary.json').write_text(json.dumps(study))
        (public/'cache-warm.html').write_text('<script type="module" src="app.js"></script>')
        (data/'pairs').mkdir()
        for cutoff in (10,20,50):
            (data/f'pairs/A549__HepG2_min{cutoff}.json').write_text(json.dumps(dict(summary=stats,candidates=[],
                scatter=[],scatter_total=50,candidate_limit=500)))

        class Handler(SimpleHTTPRequestHandler):
            def do_GET(self):
                url=urlsplit(self.path)
                if url.path=='/app.js' and parse_qs(url.query).get('v')!=['dashboard-overview-4']:
                    # Simulate an old app.js already cached before the HTML update.
                    body=b'window.__oldAppLoaded = true;'
                    self.send_response(200);self.send_header('Content-Type','application/javascript')
                    self.send_header('Cache-Control','public, max-age=31536000')
                    self.send_header('Content-Length',str(len(body)));self.end_headers();self.wfile.write(body);return
                if urlsplit(self.path).path=='/api/transcript':
                    body=json.dumps({'id':tx,'rows':all_rows}).encode()
                    self.send_response(200);self.send_header('Content-Type','application/json')
                    self.send_header('Content-Length',str(len(body)));self.end_headers();self.wfile.write(body)
                else:super().do_GET()
            def log_message(self,*args):pass

        server=ThreadingHTTPServer(('127.0.0.1',0),partial(Handler,directory=str(public)))
        thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
        try:
            with sync_playwright() as p:
                executable = os.environ.get("CHROMIUM_EXECUTABLE") or shutil.which("chromium") or shutil.which("chromium-browser")
                launch_options = {"args": ["--no-sandbox"]}
                if executable:
                    launch_options["executable_path"] = executable
                browser=p.chromium.launch(**launch_options)
                page=browser.new_page(viewport={'width':1440,'height':1000},accept_downloads=True)
                page.add_init_script('''
                  window.__trackText=[];
                  window.__overviewRects=[];
                  window.__overviewLabels=[];
                  const originalFillRect=CanvasRenderingContext2D.prototype.fillRect;
                  CanvasRenderingContext2D.prototype.fillRect=function(x,y,w,h){
                    if(this.canvas.id==='transcript-tracks-overview'&&y===22&&h===58){window.__overviewRects=[];window.__overviewLabels=[];}
                    return originalFillRect.call(this,x,y,w,h);
                  };
                  const originalStrokeRect=CanvasRenderingContext2D.prototype.strokeRect;
                  CanvasRenderingContext2D.prototype.strokeRect=function(x,y,w,h){
                    if(this.canvas.id==='transcript-tracks-overview')window.__overviewRects.push({x,y,w,h});
                    return originalStrokeRect.call(this,x,y,w,h);
                  };

                  const original=CanvasRenderingContext2D.prototype.fillText;
                  CanvasRenderingContext2D.prototype.fillText=function(text,x,y,...rest){
                    if(this.canvas.id==='transcript-tracks-overview')window.__overviewLabels.push({text,x,y,width:this.measureText(text).width,align:this.textAlign});
                    if(this.canvas.id==='transcript-tracks-canvas') {
                      window.__trackText.push({text,x,y});
                      if(window.__trackText.length>10000)window.__trackText.shift();
                    }
                    return original.call(this,text,x,y,...rest);
                  };
                ''')
                errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
                obsolete=[];page.on('request',lambda r:obsolete.append(r.url) if any(name in r.url for name in ['transcript-tracks.mjs','transcript-tracks-core.mjs','transcript-tracks.css']) else None)
                page.goto(f'http://127.0.0.1:{server.server_port}/cache-warm.html')
                page.wait_for_function('window.__oldAppLoaded === true')
                page.goto(f'http://127.0.0.1:{server.server_port}/index.html#transcript?id={tx}')
                page.wait_for_selector('#transcript-heatmap table')
                assert page.locator('#transcript-tracks-canvas').is_hidden()
                assert page.locator('#transcript-tracks-overview').is_visible()
                assert page.locator('.track-toolbar').is_hidden()
                assert page.locator('#track-range-form').is_hidden()
                assert page.locator('#track-keyboard-help').is_hidden()
                assert page.locator('#track-window-status').is_hidden()
                assert page.locator('#transcript-tracks-overview').evaluate('(el)=>el.width > 300')
                assert page.locator('#track-overview-summary').evaluate('(el)=>getComputedStyle(el).fontSize') == page.locator('#track-window-status').evaluate('(el)=>getComputedStyle(el).fontSize')
                ticks=page.evaluate("window.__overviewLabels.filter(r=>/^[0-9,]+$/.test(r.text))")
                assert ticks[0]['text']=='100' and ticks[-1]['text']=='8,183'
                assert all(t['y'] < 22 for t in ticks)
                assert len({t['text'] for t in ticks})==len(ticks)

                page.wait_for_function('window.__overviewRects.length===0')
                assert page.locator('#track-overview-summary').inner_text()=='60 eligible positions · selected position: 100'
                # A collapsed overview click selects the nearest eligible site, including off-page.
                overview=page.locator('#transcript-tracks-overview')
                overview.click(position={'x':overview.bounding_box()['width']-18,'y':45})
                assert page.locator('#position-heading').inner_text()==f'Position {positions[-1]:,}'
                assert '49–60' in page.locator('#position-count').inner_text()
                assert f'selected position: {positions[-1]:,}' in page.locator('#track-overview-summary').inner_text()
                # Dragging the black indicator cannot move the selected site or track range.
                overview.scroll_into_view_if_needed();box=overview.bounding_box()
                page.mouse.move(box['x']+box['width']-16,box['y']+45);page.mouse.down()
                page.mouse.move(box['x']+box['width']-66,box['y']+45,steps=6);page.mouse.up()
                assert page.locator('#position-heading').inner_text()==f'Position {positions[-1]:,}'
                assert page.locator('#track-start').input_value()=='100'
                overview.focus();page.keyboard.press('Home')
                assert page.locator('#position-heading').inner_text()=='Position 100'
                assert page.locator('#track-toggle-label').inner_text()=='Expand tracks'
                assert page.locator('link[href*="transcript-tracks.css"]').count()==0
                _expand_dashboard_tracks(page)
                assert page.locator('#track-threshold-legend').is_visible()
                assert page.locator('.track-toolbar').is_visible()
                assert page.evaluate('window.__overviewRects.length')==1
                # Dragging the rectangle pans only the track window and keeps selection fixed.
                page.locator('[data-track-zoom="2"]').click()
                assert page.evaluate("window.__overviewLabels.filter(r=>/^[0-9,]+$/.test(r.text)).map(r=>r.text)")==[t['text'] for t in ticks]
                rect=page.evaluate('window.__overviewRects[0]')
                old_start=int(page.locator('#track-start').input_value())
                old_selection=page.locator('#track-overview-summary').inner_text()
                overview.scroll_into_view_if_needed();box=overview.bounding_box()
                x=box['x']+rect['x']+rect['w']/2
                page.mouse.move(x,box['y']+45);page.mouse.down()
                page.mouse.move(x+45,box['y']+45,steps=8);page.mouse.up()
                assert int(page.locator('#track-start').input_value())>old_start
                assert page.locator('#track-overview-summary').inner_text()==old_selection
                preserved=[page.locator('#track-start').input_value(),page.locator('#track-end').input_value()]
                page.locator('.track-disclosure').click()
                page.wait_for_function('window.__overviewRects.length===0')
                assert preserved==[page.locator('#track-start').input_value(),page.locator('#track-end').input_value()]
                _expand_dashboard_tracks(page)
                assert preserved==[page.locator('#track-start').input_value(),page.locator('#track-end').input_value()]
                # A no-drag click within the rectangle selects a site rather than panning.
                rect=page.evaluate('window.__overviewRects[0]')
                overview.click(position={'x':rect['x']+rect['w']/2,'y':45})
                assert page.locator('#track-overview-summary').inner_text()!=old_selection
                # The selection indicator is not a drag handle, even inside a zoomed rectangle.
                rect=page.evaluate('window.__overviewRects[0]')
                first=int(page.locator('#track-start').input_value());last=int(page.locator('#track-end').input_value())
                selected=int(page.locator('#position-heading').inner_text().replace('Position ','').replace(',',''))
                line_x=rect['x']+(selected-first)/(last-first)*rect['w']
                overview.scroll_into_view_if_needed();box=overview.bounding_box()
                page.mouse.move(box['x']+line_x,box['y']+45);page.mouse.down()
                page.mouse.move(box['x']+line_x+35,box['y']+45,steps=6);page.mouse.up()
                assert [first,last]==[int(page.locator('#track-start').input_value()),int(page.locator('#track-end').input_value())]
                assert page.locator('#position-heading').inner_text()==f'Position {selected:,}'
                page.locator('#track-full-span').click()

                assert 'H04 prediction threshold' in page.locator('#track-threshold-legend').inner_text()
                assert page.locator('#track-threshold-value').inner_text()==f'{threshold:.4f}'
                page.wait_for_function("document.querySelector('#track-start').value === '100'")
                assert page.locator('#track-load-message').is_hidden(), page.locator('#track-load-message').inner_text()
                assert page.locator('#track-end').input_value()==str(positions[-1])
                assert page.locator('a[href="analytics.html"]').count()==0
                canvas=page.locator('#transcript-tracks-canvas')
                # Selecting a distant site jumps to the correct heatmap page (49–60).
                width=canvas.bounding_box()['width']
                canvas.click(position={'x':width-18,'y':55})
                page.wait_for_function(f"document.querySelector('#position-heading').textContent === 'Position {positions[-1]:,}'")
                assert '49–60' in page.locator('#position-count').inner_text()
                assert page.locator(f'.position-head[data-select-position="{positions[-1]}"]').get_attribute('aria-pressed')=='true'
                assert f'position={positions[-1]}' in page.url
                assert page.locator('#position-table tbody tr').count()==6
                # Hover reveals the exact run observation, including read count.
                canvas.hover(position={'x':width-18,'y':110})
                assert page.locator('#transcript-tracks-tooltip').is_visible()
                assert 'reads' in page.locator('#transcript-tracks-tooltip').inner_text()
                # Heatmap selection updates the track selection and detail table.
                target=positions[50]
                page.locator(f'.position-head[data-select-position="{target}"]').click()
                assert f'selected {target:,}' in page.locator('#track-window-status').inner_text()
                assert page.locator('#position-heading').inner_text()==f'Position {target:,}'
                # Collapsing leaves the table/heatmap interactive and reopening stays linked.
                previous_window=[page.locator('#track-start').input_value(),page.locator('#track-end').input_value()]
                page.locator('.track-disclosure').click()
                assert page.locator('#transcript-tracks-canvas').is_hidden()
                assert page.locator('#transcript-tracks-overview').is_visible()
                assert page.locator('#track-threshold-legend').is_hidden()
                page.wait_for_function('window.__overviewRects.length===0')
                assert page.locator('.track-toolbar').is_hidden()
                page.locator(f'.position-head[data-select-position="{positions[51]}"]').click()
                _expand_dashboard_tracks(page)
                assert f'selected {positions[51]:,}' in page.locator('#track-window-status').inner_text()
                assert previous_window==[page.locator('#track-start').input_value(),page.locator('#track-end').input_value()]

                # Zoom, move, full extent, coordinate entry, invalid ranges, and drag.
                page.locator('[data-track-zoom="2"]').click()
                start=int(page.locator('#track-start').input_value());end=int(page.locator('#track-end').input_value())
                assert end-start<positions[-1]-positions[0]
                page.locator('[data-track-pan="0.5"]').click()
                assert int(page.locator('#track-start').input_value())>start
                page.locator('#track-full-span').click()
                canvas.dispatch_event('wheel',{'deltaY':-100,'ctrlKey':True,'clientX':canvas.bounding_box()['x']+250,'clientY':canvas.bounding_box()['y']+55})
                assert int(page.locator('#track-end').input_value())-int(page.locator('#track-start').input_value())<positions[-1]-positions[0]
                page.locator('#track-full-span').click()
                page.locator('#track-start').fill('1000');page.locator('#track-end').fill('3000')
                page.locator('#track-range-apply').click()
                assert page.locator('#track-start').input_value()=='1000'
                page.locator('#track-start').fill('3000');page.locator('#track-end').fill('1000')
                page.locator('#track-range-apply').click()
                assert page.locator('#track-range-error').is_visible()
                # Restore valid range, then drag horizontally to pan.
                page.locator('#track-start').fill('1000');page.locator('#track-end').fill('3000')
                page.locator('#track-range-apply').click()
                assert page.locator('#track-range-error').is_hidden()
                canvas.scroll_into_view_if_needed()
                box=canvas.bounding_box()
                page.mouse.move(box['x']+400,box['y']+55);page.mouse.down()
                page.mouse.move(box['x']+480,box['y']+55,steps=8);page.mouse.up()
                assert int(page.locator('#track-start').input_value())<1000
                # Full-span keyboard selection synchronizes off-page positions.
                page.locator('#track-full-span').click();canvas.focus();page.keyboard.press('End')
                assert page.locator('#position-heading').inner_text()==f'Position {positions[-1]:,}'
                page.keyboard.press('Home')
                assert page.locator('#position-heading').inner_text()=='Position 100'
                assert '1–24' in page.locator('#position-count').inner_text()
                # Existing cell-line and coverage filters update tracks AND tables.
                page.locator('#transcript-cell-summary').click()
                page.locator('.transcript-cell-choice[value="K562"]').uncheck()
                assert 'A549, HepG2' == page.locator('#transcript-cell-summary').inner_text()
                assert page.locator('#transcript-cell-all').evaluate('(input)=>input.indeterminate')
                assert 'K562' not in page.locator('#transcript-heatmap').inner_text()
                with page.expect_download() as multi_download:page.locator('#download-transcript').click()
                multi_csv=Path(multi_download.value.path()).read_text()
                assert 'A549' in multi_csv and 'HepG2' in multi_csv and 'K562' not in multi_csv
                page.locator('#transcript-cell-summary').click()
                page.locator('#transcript-cell-all').check()
                page.locator('#transcript-cell-all').uncheck()
                assert 'No cell lines selected' == page.locator('#transcript-cell-summary').inner_text()
                assert page.locator('#track-empty').is_visible()
                assert page.locator('#position-table').inner_text()==''
                page.locator('#transcript-cell-all').check()
                page.locator('#transcript-cell-summary').click()
                page.locator('#transcript-cell-summary').click()
                page.locator('.transcript-cell-choice[value="A549"]').uncheck()
                page.locator('.transcript-cell-choice[value="K562"]').uncheck()
                page.locator('#transcript-cell-summary').click()
                assert page.locator('#position-table tbody tr').count()==1
                assert 'HepG2' in page.locator('#position-table').inner_text()
                page.locator('[data-coverage="50"]').click()
                page.wait_for_function("document.querySelector('#position-heading').textContent !== 'Position 100'")
                with page.expect_download() as download:page.locator('#download-transcript').click()
                exported=Path(download.value.path()).read_text()
                assert 'SYNTHETIC-TEST-ONLY' in exported
                assert 'A549' not in exported
                # No qualifying observations: clear selection, disable navigation.
                page.locator('#minimum-reads').fill('1000000');page.locator('#coverage-form button').click()
                page.wait_for_selector('#track-empty:visible')
                assert page.locator('#position-heading').inner_text()=='No eligible position selected'
                assert page.locator('[data-track-zoom="2"]').is_disabled()
                assert page.locator('#download-transcript').is_disabled()
                # Existing <10-read query path also supplies the tracks.
                page.locator('#minimum-reads').fill('1');page.locator('#coverage-form button').click()
                page.wait_for_function("document.querySelector('#track-scale-note').textContent.includes('9,000')")
                page.locator('#transcript-cell-summary').click()
                page.locator('#transcript-cell-all').check()
                page.locator('#transcript-cell-summary').click()
                page.locator('#track-full-span').click()
                canvas.focus();page.keyboard.press('End')
                assert page.locator('#position-heading').inner_text()=='Position 9,000'
                assert '5' in page.locator('#position-table').inner_text()
                # Original overview and compare pages still render normally.
                page.locator('a[data-page="overview"]').click()
                # Reset coverage to a saved preset before loading overview.
                page.locator('[data-coverage="20"]').click()
                page.wait_for_selector('#overview-content:visible')
                assert page.locator('.fraction-dot').count()==6
                page.locator('a[data-page="compare"]').click()
                page.wait_for_selector('#compare-content:visible')
                assert page.locator('#compare-stats').is_visible()
                # Return to transcript and capture desktop/mobile for visual review.
                page.goto(f'http://127.0.0.1:{server.server_port}/index.html#transcript?id={tx}')
                _expand_dashboard_tracks(page)
                page.locator('#transcript-tracks').scroll_into_view_if_needed()
                page.evaluate("document.querySelector('#toast').hidden=true;document.activeElement.blur()")
                artifact=root/'preview';artifact.mkdir(exist_ok=True)
                page.locator('#transcript-tracks').screenshot(path=str(artifact/'tracks-desktop-synthetic.png'))
                page.set_viewport_size({'width':400,'height':900})
                page.wait_for_timeout(150)
                assert page.evaluate('document.documentElement.scrollWidth <= innerWidth')
                mobile_ticks=page.evaluate("window.__overviewLabels.filter(r=>/^[0-9,]+$/.test(r.text))")
                assert mobile_ticks[0]['text']=='100' and mobile_ticks[-1]['text']=='8,183'
                bounds=[(t['x'] if t['align']=='left' else t['x']-t['width'] if t['align']=='right' else t['x']-t['width']/2, t['width']) for t in mobile_ticks]
                assert all(a[0]+a[1]<b[0] for a,b in zip(bounds,bounds[1:]))

                page.locator('#transcript-tracks').screenshot(path=str(artifact/'tracks-mobile-synthetic.png'))
                assert page.evaluate("window.__trackText.some(r=>r.text==='Cell line · replicate / run')")
                # An empty run uses a heading-level message, not a plot baseline.
                page.set_viewport_size({'width':1440,'height':1000})
                page.locator('#track-start').fill('100');page.locator('#track-end').fill('101')
                page.locator('#track-range-apply').click()
                assert page.evaluate("window.__trackText.some(r=>r.text==='No qualifying sites in this window' && r.y === 484)")
                page.locator('#transcript-tracks').screenshot(path=str(artifact/'tracks-empty-window-synthetic.png'))
                # Single-site transcript must retain a usable one-base window.
                page.route('**/data/transcripts/fixture.json.gz',lambda route:route.fulfill(
                    content_type='application/json',body=json.dumps({tx:[[0,100,'AGAACAT',50,.6]]})))
                page.reload();_expand_dashboard_tracks(page)
                assert page.locator('#track-end').input_value()=='101'
                assert page.locator('[data-track-zoom="2"]').is_disabled()
                assert not errors,errors
                assert not obsolete,obsolete
                # A renderer failure must surface a message without breaking tables.
                failure=browser.new_page()
                broken_app=(public/'app.js').read_text().replace('\nstart();', '\nTranscriptTracks.prototype.update = function() { throw new Error("Synthetic renderer failure"); };\nstart();')
                failure.route('**/app.js?v=dashboard-overview-4',lambda route:route.fulfill(
                    content_type='application/javascript',body=broken_app))
                failure.goto(f'http://127.0.0.1:{server.server_port}/index.html#transcript?id={tx}')
                failure.wait_for_selector('#transcript-heatmap table')
                failure.locator('.track-disclosure').click()
                failure.wait_for_selector('#track-load-message.error:visible')
                assert 'Synthetic renderer failure' in failure.locator('#track-load-message').inner_text()
                assert failure.locator('#position-table tbody tr').count()>0
                failure.close()
                browser.close()
            print('Browser checks passed: stale-app cache bypass, visible render errors, linked selection, paging, zoom/pan/drag, range validation, keyboard, filters, CSV, empty/single-site data, original pages and mobile layout.')
        finally:
            server.shutdown();server.server_close();thread.join()
