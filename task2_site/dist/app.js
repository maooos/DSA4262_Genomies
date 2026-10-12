const $ = id => document.getElementById(id);
const esc = value => String(value).replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const fmt = value => value == null ? '—' : Number(value).toLocaleString('en-US');
const decimal = (value, n=3) => value == null ? '—' : Number(value).toFixed(n);
const pct = value => value == null ? '—' : `${(100*value).toFixed(1)}%`;
const palette = ['#087d83','#d18c54','#739581','#737ba2','#b06a72','#85a3b0','#c1a44c'];
let data, transcriptIndex, transcriptAliases, currentPair, currentTranscript, transcriptRows=[], visibleRows=200, requestNumber=0;
const coverage = () => +$('coverage').value;
const cutoff = () => +$('threshold').value;
const colour = cell => palette[data.cell_lines.indexOf(cell) % palette.length];
const showError = error => {$('error').textContent=error.message || String(error); $('error').hidden=false;};
const clearError = () => { if(data?.status==='partial')showError('This is an incomplete pilot. Some cell lines or replicates are still missing; results will change when the full inventory is processed.'); else $('error').hidden=true; };
async function getJSON(path) {
  const r=await fetch(path); if(!r.ok) throw Error(`Could not load ${path} (${r.status}).`);
  if(!path.endsWith('.gz'))return r.json();
  const bytes=await r.arrayBuffer(),header=new Uint8Array(bytes);
  // Some hosts transparently decode Content-Encoding; detect the gzip bytes.
  if(header[0]===31&&header[1]===139){
    if(!globalThis.DecompressionStream)throw Error('Transcript lookup requires a browser with gzip decompression support. Please use a current browser.');
    return new Response(new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'))).json();
  }
  return JSON.parse(new TextDecoder().decode(bytes));
}
function metrics(id, entries) { $(id).innerHTML=entries.map(([label,value,note])=>`<div class="metric"><small>${esc(label)}</small><strong>${esc(value)}</strong><p>${esc(note)}</p></div>`).join(''); }
function svg(inner,width=560,height=320,label='Data chart') { return `<svg class="chart" role="img" aria-label="${esc(label)}" viewBox="0 0 ${width} ${height}">${inner}</svg>`; }
function downloadCSV(name, rows) {
  if(!rows?.length) return;
  const cols=Object.keys(rows[0]);
  const quote=v=>'"'+String(v??'').replaceAll('"','""')+'"';
  const text=[cols.map(quote).join(','),...rows.map(r=>cols.map(c=>quote(r[c])).join(','))].join('\r\n');
  const url=URL.createObjectURL(new Blob([text],{type:'text/csv;charset=utf-8'}));
  const a=document.createElement('a'); a.href=url; a.download=name; a.click(); setTimeout(()=>URL.revokeObjectURL(url),1000);
}
function selectedSummaries() { return data.summaries.filter(s=>s.min_reads===coverage() && Math.abs(s.threshold-cutoff())<1e-12); }
function renderHistogram() {
  const sample=$('histogram-sample').value;
  const histogram=data.histograms.find(h=>h.sample_id===sample&&h.min_reads===coverage());
  if(!histogram)return;
  const left=72,top=15,width=950,height=235,plotWidth=850,plotHeight=172,max=Math.max(1,...histogram.counts);
  let inner='';
  for(let i=0;i<=4;i++){const y=top+plotHeight-i*plotHeight/4;inner+=`<line class="grid" x1="${left}" x2="${left+plotWidth}" y1="${y}" y2="${y}"/><text x="${left-10}" y="${y+4}" text-anchor="end">${fmt(Math.round(max*i/4))}</text>`;}
  histogram.counts.forEach((n,i)=>{const h=n/max*plotHeight;inner+=`<rect x="${left+i*plotWidth/20+1}" y="${top+plotHeight-h}" width="${plotWidth/20-2}" height="${h}" fill="#087d83" opacity=".8"><title>${(i*.05).toFixed(2)}–${((i+1)*.05).toFixed(2)}: ${fmt(n)} sites</title></rect>`;});
  for(let i=0;i<=10;i++)inner+=`<text x="${left+i*plotWidth/10}" y="${height-25}" text-anchor="middle">${(i/10).toFixed(1)}</text>`;
  const x=left+cutoff()*plotWidth;inner+=`<line x1="${x}" x2="${x}" y1="${top}" y2="${top+plotHeight}" stroke="#d56b37" stroke-width="2" stroke-dasharray="4 3"/><text x="${left+plotWidth/2}" y="${height-4}" text-anchor="middle" class="axis-label">H04 score</text>`;
  $('histogram-chart').innerHTML=svg(inner,width,height,'Distribution of H04 scores for the selected sequencing run');
}
function horizontalDots(rows, valueKey, max, formatter, label) {
  const width=560,height=42*data.cell_lines.length+50,left=78,right=30,top=18,plot=width-left-right;
  let content='';
  for(let i=0;i<=4;i++){const value=max*i/4,x=left+plot*i/4;content+=`<line class="grid" x1="${x}" x2="${x}" y1="${top-6}" y2="${height-36}"/><text x="${x}" y="${height-13}" text-anchor="middle">${esc(formatter(value))}</text>`;}
  data.cell_lines.forEach((cell,i)=>{
    const cy=top+16+i*42;
    content+=`<text x="${left-13}" y="${cy+4}" text-anchor="end" class="axis-label">${esc(cell)}</text>`;
    const group=rows.filter(r=>r.cell_line===cell&&r[valueKey]!=null);
    group.forEach((r,j)=>{const x=left+plot*r[valueKey]/max,y=cy+(j-(group.length-1)/2)*5;content+=`<circle cx="${x}" cy="${y}" r="5" fill="${colour(cell)}" opacity=".82"><title>${esc(cell)} · ${esc(r.sample_id||`replicates ${r.replicate_a} / ${r.replicate_b}`)}: ${esc(formatter(r[valueKey]))}</title></circle>`;});
    if(!group.length) content+=`<text x="${left+8}" y="${cy+4}">No eligible pair</text>`;
  });
  return svg(content,width,height,label);
}
function renderOverview() {
  const rows=selectedSummaries();
  const reps=new Set(data.samples.map(s=>`${s.cell_line}/${s.replicate}`)).size;
  metrics('metrics',[
    ['Cell lines',data.cell_lines.length,'Available processed direct RNA data'],
    ['Sequencing runs',`${data.samples.length} / ${data.expected_samples}`,`${reps} named replicate groups`],
    ['Eligible observations',fmt(rows.reduce((n,s)=>n+s.eligible_sites,0)),`Sample × site observations · ≥${coverage()} reads`],
    ['Score threshold',decimal(cutoff(),4),Math.abs(cutoff()-data.threshold)<1e-12?'Saved validation threshold':'Sensitivity view · model unchanged']
  ]);
  const max=Math.max(.05,...rows.map(r=>r.positive_fraction||0))*1.12;
  $('fraction-chart').innerHTML=horizontalDots(rows,'positive_fraction',max,pct,'Predicted-positive fraction per sequencing run');
  const agreement=data.scenarios[coverage()].agreement;
  // Spearman values may be negative: use a signed scale rather than clipping them.
  let inner='';const w=520,h=42*data.cell_lines.length+50,left=76,plot=414;
  [-1,-.5,0,.5,1].forEach(v=>{const x=left+(v+1)*plot/2;inner+=`<line class="grid" x1="${x}" x2="${x}" y1="12" y2="${h-36}"/><text x="${x}" y="${h-13}" text-anchor="middle">${v}</text>`;});
  data.cell_lines.forEach((cell,i)=>{let cy=34+i*42;inner+=`<text x="${left-14}" y="${cy+4}" text-anchor="end" class="axis-label">${esc(cell)}</text>`;let rows=agreement.filter(r=>r.cell_line===cell&&r.spearman!=null);rows.forEach((r,j)=>{inner+=`<circle cx="${left+(r.spearman+1)*plot/2}" cy="${cy+5*(j-(rows.length-1)/2)}" r="5" fill="${colour(cell)}"><title>Replicates ${r.replicate_a}/${r.replicate_b}: ρ=${decimal(r.spearman)}; ${fmt(r.common_sites)} shared sites</title></circle>`;});if(!rows.length)inner+=`<text x="${left+10}" y="${cy+4}">No eligible pair</text>`;});
  $('agreement-chart').innerHTML=svg(inner,w,h,'Spearman correlation between replicate pairs');
  $('sample-table').innerHTML=rows.map(s=>`<tr><td class="cell-name">${esc(s.cell_line)}</td><td>${s.replicate} / ${s.run}</td><td>${fmt(s.scored_sites)}</td><td>${fmt(s.eligible_sites)}</td><td>${fmt(s.positive_sites)}</td><td>${pct(s.positive_fraction)}</td><td>${fmt(s.median_reads)}</td></tr>`).join('');
  renderHistogram();
}
async function renderComparison() {
  const a=$('cell-a').value,b=$('cell-b').value,id=++requestNumber;
  currentPair=null;$('download-candidates').disabled=true;$('candidate-table').innerHTML='';
  if(a===b){$('comparison-metrics').innerHTML='';$('scatter-chart').textContent='Choose two different cell lines.';$('overlap-chart').innerHTML='';$('scatter-note').textContent='';return;}
  $('scatter-chart').textContent='Loading shared-site comparison…';
  try {
    const [first,second]=[a,b].sort((x,y)=>data.cell_lines.indexOf(x)-data.cell_lines.indexOf(y));
    const pair=await getJSON(`data/pairs/${first}__${second}_min${coverage()}.json`);
    if(id!==requestNumber)return;
    clearError();
    const reverse=first!==a,s=pair.summary;
    const candidates=pair.candidates.map(r=>({...r,score_a:reverse?r.mean_score_b:r.mean_score_a,score_b:reverse?r.mean_score_a:r.mean_score_b,delta:reverse?-r.delta_score:r.delta_score,reps_a:reverse?r.n_replicates_b:r.n_replicates_a,reps_b:reverse?r.n_replicates_a:r.n_replicates_b}));
    currentPair={...pair,candidates};$('download-candidates').disabled=!candidates.length;
    metrics('comparison-metrics',[
      ['Common eligible sites',fmt(s.common_sites),'≥2 replicates per cell line'],['Score correlation',decimal(s.spearman),'Spearman ρ on shared eligible sites'],
      ['Positive-site overlap',decimal(s.jaccard),'Jaccard · shared eligible universe'],['Consistent candidates',fmt(s.consistent_candidates),'Opposite calls in all observed replicates']
    ]);
    const size=380,left=62,top=18,plot=290;
    let inner='';
    for(let i=0;i<=5;i++){const v=i/5,x=left+v*plot,y=top+(1-v)*plot;inner+=`<line class="grid" x1="${x}" x2="${x}" y1="${top}" y2="${top+plot}"/><line class="grid" x1="${left}" x2="${left+plot}" y1="${y}" y2="${y}"/><text x="${x}" y="${top+plot+19}" text-anchor="middle">${v}</text><text x="${left-12}" y="${y+4}" text-anchor="end">${v}</text>`;}
    inner+=`<line x1="${left}" y1="${top+plot}" x2="${left+plot}" y2="${top}" stroke="#b5c7ba" stroke-dasharray="4 4"/>`;
    const tx=left+data.threshold*plot,ty=top+(1-data.threshold)*plot;
    inner+=`<path d="M${tx} ${top}V${top+plot}M${left} ${ty}H${left+plot}" stroke="#d8b090" fill="none" stroke-dasharray="3 4"/>`;
    pair.scatter.forEach(r=>{const x=reverse?r.mean_score_b:r.mean_score_a,y=reverse?r.mean_score_a:r.mean_score_b;inner+=`<circle cx="${left+x*plot}" cy="${top+(1-y)*plot}" r="2.2" fill="#087d83" opacity=".35"><title>${esc(r.transcript_id)}:${r.transcript_position} · ${decimal(x)} / ${decimal(y)}</title></circle>`;});
    inner+=`<text class="axis-label" x="${left+plot/2}" y="${size-10}" text-anchor="middle">${esc(a)} · mean replicate score</text><text class="axis-label" transform="translate(14 ${top+plot/2}) rotate(-90)" text-anchor="middle">${esc(b)} · mean replicate score</text>`;
    $('scatter-chart').innerHTML=svg(inner,size,size,`${a} versus ${b} mean scores`);
    $('scatter-note').textContent=`${fmt(pair.scatter.length)} deterministically sampled sites plotted out of ${fmt(pair.scatter_total)}. Statistics use all shared sites. Dashed orange lines mark the saved threshold.`;
    const groups=[['Both positive',s.positive_both,'#087d83'],[`${a} only`,reverse?s.positive_b_only:s.positive_a_only,colour(a)],[`${b} only`,reverse?s.positive_a_only:s.positive_b_only,colour(b)],['Both negative',s.negative_both,'#bdcabe']];
    const max=Math.max(1,...groups.map(r=>r[1]));let bars='';
    groups.forEach(([name,count,c],i)=>{const y=40+i*68;bars+=`<text x="15" y="${y}" class="axis-label">${esc(name)}</text><text x="480" y="${y}" text-anchor="end">${fmt(count)}</text><rect x="15" y="${y+11}" width="465" height="13" rx="4" fill="#f0f4ed"/><rect x="15" y="${y+11}" width="${465*count/max}" height="13" rx="4" fill="${c}"/>`;});
    $('overlap-chart').innerHTML=svg(bars,500,340,'Counts of positive and negative calls on shared eligible sites');
    $('candidate-table').innerHTML=candidates.length?candidates.map(r=>`<tr><td><button class="link-button" data-transcript="${esc(r.transcript_id)}">${esc(r.transcript_id)} : ${r.transcript_position}</button></td><td>${esc(r.sequence)}</td><td>${decimal(r.score_a)}</td><td>${decimal(r.score_b)}</td><td>${r.delta>0?'+':''}${decimal(r.delta)}</td><td>${r.reps_a} / ${r.reps_b}</td></tr>`).join(''):'<tr><td colspan="6">No consistent candidates under this coverage setting.</td></tr>';
  } catch(error) {if(id===requestNumber){$('scatter-chart').textContent='Comparison unavailable.';showError(error);}}
}
async function loadTranscript(tx) {
  if(!transcriptIndex)transcriptIndex=await getJSON('data/transcript_index.json');
  if(!transcriptAliases)transcriptAliases=await getJSON('data/transcript_aliases.json');
  tx=transcriptAliases[tx]||tx;
  if(!Object.hasOwn(transcriptIndex,tx))throw Error('Transcript not found among exported observations with at least 10 reads. Use the exact ID from the suggestions.');
  $('transcript-title').textContent=`${tx} · loading…`;
  const shard=await getJSON(`data/transcripts/${transcriptIndex[tx]}.json.gz`);
  currentTranscript=tx;transcriptRows=shard[tx];visibleRows=200;$('transcript-input').value=tx;clearError();renderTranscript();
}
function transcriptObservations(){return transcriptRows.filter(r=>r[3]>=coverage()).map(([sampleIndex,position,sequence,n_reads,score,originalId])=>({...data.samples[sampleIndex],transcript_id:currentTranscript,transcript_id_original:originalId||currentTranscript,position,sequence,n_reads,score,prediction:Number(score>=data.threshold)})).sort((a,b)=>a.position-b.position||a.cell_line.localeCompare(b.cell_line)||a.replicate-b.replicate||a.run-b.run);}
function renderTranscript(){
  if(!currentTranscript)return;
  const rows=transcriptObservations();$('transcript-title').textContent=`${currentTranscript} · ${fmt(rows.length)} observations`;$('download-transcript').disabled=!rows.length;
  $('transcript-table').innerHTML=rows.slice(0,visibleRows).map(r=>`<tr><td>${r.position}</td><td>${esc(r.sequence)}</td><td class="cell-name">${esc(r.cell_line)}</td><td>${r.replicate} / ${r.run}</td><td>${fmt(r.n_reads)}</td><td>${decimal(r.score,4)}</td><td><span class="pill ${r.prediction?'':'negative'}">${r.prediction?'Positive':'Negative'}</span></td></tr>`).join('');
  $('more-observations').hidden=rows.length<=visibleRows;
  $('transcript-chart').classList.remove('empty-state');
  if(!rows.length){$('transcript-chart').textContent='No observations meet the selected coverage.';return;}
  const positions=rows.map(r=>r.position),min=Math.min(...positions),max=Math.max(...positions),left=80,plot=870,h=72*data.cell_lines.length+45;
  let inner='';
  for(let i=0;i<=5;i++){const p=Math.round(min+(max-min)*i/5),x=left+plot*i/5;inner+=`<line class="grid" x1="${x}" x2="${x}" y1="12" y2="${h-35}"/><text x="${x}" y="${h-12}" text-anchor="middle">${fmt(p)}</text>`;}
  data.cell_lines.forEach((cell,i)=>{const top=12+i*72;inner+=`<text x="${left-14}" y="${top+32}" text-anchor="end" class="axis-label">${esc(cell)}</text><line x1="${left}" x2="${left+plot}" y1="${top+54-data.threshold*48}" y2="${top+54-data.threshold*48}" stroke="#dfc6ab" stroke-dasharray="3 4"/>`;rows.filter(r=>r.cell_line===cell).forEach(r=>{const x=left+plot*(r.position-min)/(max-min||1),y=top+54-r.score*48;inner+=`<circle cx="${x}" cy="${y}" r="3.5" fill="${colour(cell)}" opacity=".7"><title>Position ${r.position} · ${esc(cell)} replicate ${r.replicate}, run ${r.run} · score ${decimal(r.score,4)} · ${r.n_reads} reads</title></circle>`;});});
  $('transcript-chart').innerHTML=svg(inner,980,h,'Transcript positions and per-run scores, faceted by cell line; vertical scale 0 to 1 within each cell line');
}
const views={overview:['THE BIG PICTURE','Signals across cell lines.','Explore predicted m6A sites, compare shared signals, and inspect the evidence across replicates.'],compare:['SHARED SITES, DIFFERENT SIGNALS','Compare cell lines.','Focus on the same covered sites, then inspect differences supported across replicates.'],explorer:['FROM ATLAS TO EVIDENCE','Follow a transcript.','Inspect individual sites, sequencing runs and read coverage in one place.'],methods:['ANALYSIS NOTES','Understand the evidence.','The model, comparison rules, data provenance and limits behind this atlas.']};
function route(){if(!data)return;const view=location.hash.slice(1).split('?')[0]||'overview',v=Object.hasOwn(views,view)?view:'overview';document.querySelectorAll('.view').forEach(s=>s.hidden=s.id!==v);document.querySelectorAll('nav a').forEach(a=>a.classList.toggle('active',a.dataset.view===v));[$('eyebrow').textContent,$('page-title').textContent,$('page-description').textContent]=views[v];$('filters').hidden=v==='methods';$('threshold-control').hidden=v!=='overview';$('filter-note').textContent=v==='overview'?'Coverage filters determine which observations are eligible.':`Comparisons and calls use the saved threshold: ${decimal(data.threshold,6)}.`;if(v==='compare')renderComparison();if(v==='explorer')renderTranscript();}
async function init(){
  data=await getJSON('data/summary.json');
  $('dataset-status').textContent=data.status==='partial'?`PILOT · ${data.samples.length}/${data.expected_samples} runs`:`${data.samples.length} processed runs · ${data.cell_lines.length} cell lines`;
  const thresholds=[...new Set(data.summaries.map(s=>s.threshold))].sort((a,b)=>a-b);
  $('threshold').innerHTML=thresholds.map(t=>`<option value="${t}" ${t===data.threshold?'selected':''}>${decimal(t,4)}${t===data.threshold?' · saved':''}</option>`).join('');
  for(const id of ['cell-a','cell-b'])$(id).innerHTML=data.cell_lines.map(c=>`<option>${esc(c)}</option>`).join('');
  $('cell-b').selectedIndex=Math.min(1,data.cell_lines.length-1);
  $('histogram-sample').innerHTML=data.samples.map(s=>`<option value="${esc(s.sample_id)}">${esc(s.cell_line)} · replicate ${s.replicate} / run ${s.run}</option>`).join('');
  $('threshold-method').textContent=`The saved validation threshold is ${decimal(data.threshold,9)}. Scores equal to the threshold are positive. The overview offers threshold sensitivity settings; cell-line comparisons use the saved threshold consistently.`;
  $('scope-method').textContent=`This build contains ${data.samples.length} of ${data.expected_samples} inventoried processed runs and ${fmt(data.total_scored_sites)} sample–site predictions. Transcript lookup retains every observation with at least 10 reads. ${data.missing_samples.length?'Missing runs: '+data.missing_samples.join(', '):'All runs in the processed inventory have been scored.'}`;
  $('provenance').textContent=`Built ${data.created_at} · H04 SHA256: ${data.model_sha256} · ${data.sequence_conflict_sites} transcript-position keys have multiple sequences; sequence is included in the matching key. ${data.identifier_harmonization?.changed_identifiers||0} versioned identifiers harmonized; original IDs are retained in observation downloads.`;
  if(data.status==='partial')showError('This is an incomplete pilot. Some cell lines or replicates are still missing; results will change when the full inventory is processed.');
  renderOverview();route();
  transcriptIndex=await getJSON('data/transcript_index.json');
  $('transcript-suggestions').innerHTML=Object.keys(transcriptIndex).sort().slice(0,80).map(tx=>`<option value="${esc(tx)}"></option>`).join('');
}
$('coverage').addEventListener('change',()=>{renderOverview();route();});
$('threshold').addEventListener('change',renderOverview);
$('histogram-sample').addEventListener('change',renderHistogram);
for(const id of ['cell-a','cell-b'])$(id).addEventListener('change',renderComparison);
$('download-summary').addEventListener('click',()=>downloadCSV(`sgnex_summary_min${coverage()}.csv`,selectedSummaries()));
$('download-candidates').addEventListener('click',()=>downloadCSV(`candidates_${$('cell-a').value}_${$('cell-b').value}_min${coverage()}.csv`,currentPair?.candidates));
$('download-transcript').addEventListener('click',()=>downloadCSV(`${currentTranscript}_min${coverage()}.csv`,transcriptObservations()));
$('more-observations').addEventListener('click',()=>{visibleRows+=200;renderTranscript();});
$('transcript-form').addEventListener('submit',async event=>{event.preventDefault();try{await loadTranscript($('transcript-input').value.trim());}catch(error){showError(error);}});
$('transcript-input').addEventListener('input',()=>{if(!transcriptIndex)return;const value=$('transcript-input').value.trim().toUpperCase();$('transcript-suggestions').innerHTML=Object.keys(transcriptIndex).filter(tx=>tx.toUpperCase().includes(value)).sort().slice(0,80).map(tx=>`<option value="${esc(tx)}"></option>`).join('');});
$('candidate-table').addEventListener('click',async event=>{const button=event.target.closest('[data-transcript]');if(button){location.hash='explorer';try{await loadTranscript(button.dataset.transcript);}catch(error){showError(error);}}});
window.addEventListener('hashchange',route);
init().catch(showError);
