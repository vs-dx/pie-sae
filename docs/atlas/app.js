/* PIE-SAE feature explorer for the bundled JSON export. */
'use strict';
const GH = 'https://github.com/vs-dx/pie-sae';
const NOTEBOOKS = {
  '01': '01_pie_architecture', '02': '02_extract_activations',
  '03': '03_train_layerwise_sae', '04': '04_sae_feature_analysis',
  '05': '05_cross_layer_analysis', '06a': '06a_statistical_significance',
  '06b': '06b_full_statistical_validation', '07': '07_causal_validation'
};
const names = ['Overview','Layer explorer','Feature explorer','Gene search','Pathway search','Cross-layer flow','Interventions','Comparisons','Methods & sources'];
const pages = ['overview','layers','features','genes','pathways','cross-layer','interventions','comparisons','methods'];
const layerNames = [...Array.from({length:4}, (_,i)=>`encoder.${i}`), ...Array.from({length:6},(_,i)=>`processor.${i}`), ...Array.from({length:4},(_,i)=>`decoder.${i}`)];
const decoderNames=layerNames.filter(x=>x.startsWith('decoder.'));
const cache=new Map();
const ui={layer:'decoder.0',feature:189,sort:'frequency',query:'',gene:'PPCS',geneFilter:'all',termQuery:'',termId:0,termFilter:'all',cross:'decoder_0_to_decoder_1',crossQuery:'',annotationLimit:12,featureFilter:'all'};
let activePage='overview', renderSequence=0, lastAppliedHash=null;
const main=document.querySelector('#main-content');
const esc=x=>String(x??'').replace(/[&<>"']/g,m=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[m]));
const fmt=(v,n=3)=>Number.isFinite(Number(v))?Number(v).toLocaleString('en-US',{minimumFractionDigits:n,maximumFractionDigits:n}):'—';
const int=v=>Number.isFinite(Number(v))?Math.round(Number(v)).toLocaleString('en-US'):'—';
const pct=(v,n=2)=>`${fmt(Number(v)*100,n)}%`;
const short=(v)=>{let x=Number(v);return !Number.isFinite(x)?'—':x===0?'0':(x<0.0001?x.toExponential(2):fmt(x,x<.01?5:4))};
const fileSource=k=>`${GH}/blob/main/notebooks/${NOTEBOOKS[k]}.ipynb`;
const link=(url,label,cls='')=>`<a href="${esc(url)}"${cls?` class="${cls}"`:''} target="_blank" rel="noopener noreferrer">${esc(label)}</a>`;
const linkNotebook=k=>link(fileSource(k),`Notebook ${k} ↗`);
const disclaimer=`<div class="callout"><strong>Inference boundary.</strong> Feature–gene associations and historical functional annotations are exploratory. Persistent gene identity, dependence within conditions, and prior feature selection limit mechanistic interpretation. No biological causality is established.</div>`;
async function load(path){
 if(!cache.has(path)){
   const promise=fetch(`./data/${path}`).then(r=>{if(!r.ok)throw Error(`Data file ${path} is missing (${r.status}). Serve this folder using a local HTTP server or publish all data/ files.`);return r.json()}).catch(err=>{cache.delete(path);throw err});
   cache.set(path,promise);
 }
 return cache.get(path);
}
function route(){
 const hash=decodeURI((location.hash||'#overview').slice(1));
 const [raw,query='']=hash.split('?');const page=pages.includes(raw)?raw:'overview';const p=new URLSearchParams(query);
 return {page,p};
}
function setHash(page,p={}){
 const q=new URLSearchParams();for(const [k,v] of Object.entries(p)){if(v!==undefined&&v!==null&&v!=='')q.set(k,v)}
 const h=`#${page}${q.toString()?'?'+q.toString():''}`;
 if(location.hash===h){render();return}location.hash=h;
}
function pageHead(k,title,sub,tools=''){
 return `<div class="page-head"><div class="section-kicker">${esc(k)}</div><div class="head-row"><div><h1 class="page-title">${title}</h1><p class="page-subtitle">${sub}</p></div>${tools}</div></div>`;
}
const metric=(label,value,meta='')=>`<div class="metric-card">${meta?`<div class="meta">${esc(meta)}</div>`:''}<strong>${value}</strong><small>${esc(label)}</small></div>`;
const panel=(title,sub,html,extra='')=>`<section class="panel ${extra}"><h2 class="block-title">${title}</h2>${sub?`<p class="block-sub">${sub}</p>`:''}${html}</section>`;
const field=(name,label,options,selected)=>`<div class="field"><label for="${esc(name)}">${esc(label)}</label><select class="select" id="${esc(name)}" data-control="${esc(name)}">${options.map(v=>`<option value="${esc(typeof v==='string'?v:v[0])}" ${String(typeof v==='string'?v:v[0])===String(selected)?'selected':''}>${esc(typeof v==='string'?v:v[1])}</option>`).join('')}</select></div>`;
const textField=(name,label,value,placeholder='Search...',cls='')=>`<div class="field ${cls}"><label for="${esc(name)}">${esc(label)}</label><input class="input" id="${esc(name)}" data-input="${esc(name)}" value="${esc(value)}" placeholder="${esc(placeholder)}" autocomplete="off" /></div>`;
const statsLine=(name,value)=>`<div class="statline"><span>${name}</span><b>${value}</b></div>`;
function barrows(data,{max=1,limit=null}={}){
 const rows=limit?data.slice(0,limit):data;
 return `<div class="bar-list">${rows.map(d=>`<div class="barline"><span class="barname">${esc(d.name)}</span><div class="bar-track" title="${esc(d.name)}: ${esc(d.value)}"><div class="bar-inner ${esc(d.group||'')}" style="width:${Math.min(100,Math.max(0,100*d.value/max))}%"></div></div><span class="barval">${esc(d.label??fmt(d.value,4))}</span></div>`).join('')}</div>`;
}
function histogram(values,{bins=12,color='#71d9cd',log=false}={}){
 const W=760,H=223,pad={left:35,right:12,top:14,bottom:30};let counts=Array(bins).fill(0);
 let max= Math.max(...values,1e-9);
 for(let val of values){const scaled=log?Math.log1p(Math.max(0,val)*1000)/Math.log1p(max*1000):Math.max(0,val)/max;const i=Math.min(bins-1,Math.floor(scaled*bins));counts[i]++}
 let top=Math.max(...counts,1),dx=(W-pad.left-pad.right)/bins;
 let bars=counts.map((count,i)=>{let bh=(H-pad.top-pad.bottom)*count/top;return `<rect x="${(pad.left+i*dx+2).toFixed(2)}" y="${(H-pad.bottom-bh).toFixed(2)}" width="${Math.max(1,dx-4).toFixed(2)}" height="${Math.max(0,bh).toFixed(2)}" rx="2" fill="${color}"><title>Bin ${i+1}: ${count.toLocaleString()} features</title></rect>`}).join('');
 let grid=[0,.25,.5,.75,1].map(t=>`<line x1="${pad.left}" y1="${pad.top+(H-pad.top-pad.bottom)*(1-t)}" x2="${W-pad.right}" y2="${pad.top+(H-pad.top-pad.bottom)*(1-t)}" stroke="#284353" stroke-dasharray="3,5"/><text x="${pad.left-8}" y="${pad.top+(H-pad.top-pad.bottom)*(1-t)+4}" font-size="10" fill="#8ba5b4" text-anchor="end">${Math.round(top*t)}</text>`).join('');
 return `<div class="chart"><svg viewBox="0 0 ${W} ${H}" width="100%" role="img" aria-label="Histogram of ${int(values.length)} features across ${bins} bins">${grid}${bars}<text x="${pad.left}" y="${H-8}" class="chart-label">0</text><text x="${W-pad.right}" y="${H-8}" class="chart-label" text-anchor="end">${short(max)}</text></svg></div><div class="hist-caption"><span>${log?'Logarithmic-like bin spacing (log1p transform)':'Equal-width bins'} · ${int(values.length)} observations</span><span>Maximum ${short(max)}</span></div>`;
}
function table(columns,rows){return `<div class="tablewrap"><table><thead><tr>${columns.map(c=>`<th>${c}</th>`).join('')}</tr></thead><tbody>${rows.join('')}</tbody></table></div>`}
function row(cells){return `<tr>${cells.map(c=>`<td>${c}</td>`).join('')}</tr>`}
function featureHref(layer,feature){return `#features?layer=${encodeURIComponent(layer)}&feature=${int(feature).replace(/,/g,'')}`}
const featureLink=(layer,id)=>`<a href="${featureHref(layer,id)}">${esc(layer)} / ${int(id)}</a>`;
function geneHref(gene){return `#genes?gene=${encodeURIComponent(gene)}`}
async function overviewPage(){
 const o=await load('overview.json');const data=o.layers;const [genes,meta]=await Promise.all([load('genes.json'),load('dataset.json')]);
 const bar=barrows(data.map(x=>({name:x.layer,value:x.r2,label:fmt(x.r2,4),group:x.layer.split('.')[0]})),{max:1});
 const clusters=o.cross_layer_transitions.map(x=>({name:`${x.layer_a} → ${x.layer_b}`,n:x.n_mutual_pairs}));
 return `<div class="hero"><div class="hero-text"><h1 class="page-title">PIE <span class="title-accent">feature atlas</span></h1><p class="page-subtitle">Explore layer-wise sparse features in the pretrained Perturbation Is Everything model. Browse ${int(meta.n_decoder_features)} decoder features, their training and validation gene rankings, historical annotations, and signed-cosine cross-layer matches.</p><div class="btnrow"><a class="action-btn" href="#features?layer=decoder.0&feature=189">Explore feature 189 →</a><a class="muted-btn" href="#cross-layer">Cross-layer relationships ↗</a></div></div>
 <div class="panel hero-figure"><div class="section-kicker">STUDY PIPELINE</div>${[['01–02','Pretrained PIE · 14 blocks'],['03–04','TopK SAE · 3,072 features/block'],['05–06','Correspondence · statistical checks'],['07','Ablations · patching · downstream response']].map((v,i)=>`<div class="flow-node"><b>${v[1]}</b><span>${v[0]}</span></div>${i<3?'<div class="flow-arrow">↓</div>':''}`).join('')}<div class="footnote">Diagram shows the research workflow, not a discovered biological mechanism.</div></div></div>
 <div class="metrics">${metric('Block-specific SAE dictionaries','14','ARCHITECTURE')}${metric('All decoder SAE features','12,288','FEATURE SPACE')}${metric('Distinct PIE output genes',int(genes.length),'GENE VOCABULARY')}${metric('Within-component adjacent transitions','11','CORRESPONDENCE')}</div>
 <div class="panel"><h2 class="block-title">Decoder enrichment · 04 and corrected 06b</h2><p class="block-sub">${esc(meta.annotations)}</p><p><a href="https://github.com/vs-dx/pie-sae/blob/main/reports/Case_Normalization_Correction.md" target="_blank" rel="noopener noreferrer">Correction record ↗</a></p></div>
 <div class="twocol">${panel('Validation reconstruction R²','The recorded R² for each independently initialized block-specific SAE.',bar+'<div class="legend"><span><i></i>Encoder</span><span><i class="processor"></i>Processor</span><span><i class="decoder"></i>Decoder</span></div>')}
 ${panel('Explore the data','',`<div class="liststack">${[['Layer explorer','3,072 validation-feature statistics in every block','#layers'],['Feature explorer','Training/validation top genes and historical annotations','#features'],['Gene search','6,642 searchable PIE output genes','#genes'],['Cross-layer flow','12,411 mutual matches over 11 transitions','#cross-layer'],['Interventions','Evidence and negative controls from notebook 07','#interventions']].map(v=>`<div class="statline"><span>${link(v[2],v[0])}<br/><small class="tiny">${v[1]}</small></span><b>↗</b></div>`).join('')}</div>`)}
 </div><div class="section-heading"><h2>Scientific interpretation</h2><small>Association ≠ causation</small></div><div class="threecol">${[
 ['01 / REPRESENTATION','Sparse reconstruction','R² = 0.8119–0.9821 across fixed block-specific dictionaries. Reconstruction does not guarantee semantic decomposition.'],
 ['02 / ASSOCIATIONS','Gene and pathway enrichment','Feature cards and pathway search use the per-layer analysis in notebook 04. Notebook 06b applies a global adjustment across all four decoder layers.'],
 ['03 / INTERVENTIONS','Computational influence','Feature 189 affects PIE outputs, but donor patching was inconsistent and unique mediation remains unproven.']
 ].map(v=>`<div class="panel"><div class="section-kicker">${v[0]}</div><h3 class="block-title" style="margin-top:14px">${v[1]}</h3><p class="block-sub">${v[2]}</p></div>`).join('')}</div>${disclaimer}
 <p class="data-credit">Sources: ${linkNotebook('01')}, ${linkNotebook('03')}, ${linkNotebook('04')}, ${linkNotebook('06b')}, ${linkNotebook('07')}. Statistical and annotation caveats apply.</p>`;
}
async function layerPage(){
 const cur=ui.layer;
 const d=await load('layers/'+cur.replace('.','_')+'.json');const f=d.features;
 const frequencies=f.map(x=>x.activation_frequency);const active=f.filter(x=>x.fire_count>0).length;const avg=frequencies.reduce((a,b)=>a+b,0)/f.length;
 const ranked=[...f].sort((a,b)=> ui.sort==='frequency'?b.activation_frequency-a.activation_frequency:ui.sort==='amplitude'?b.mean_active_amplitude-a.mean_active_amplitude:a.feature_id-b.feature_id).slice(0,28);
 const options=layerNames.map(s=>[s,s]);
 return pageHead('LAYER EXPLORER',`Inside <span class="title-accent">${esc(cur)}</span>`,'Held-out validation feature activity across the 3,072-dimensional TopK dictionary. Encoder/processor features represent latent slots; decoder features are associated with output genes.')+
 `<div class="controls">${field('layer','PIE block',options,cur)}${field('sort','Sort feature table',[['frequency','Highest firing frequency'],['amplitude','Highest active amplitude'],['id','Feature ID']],ui.sort)}<span class="chip green">VALIDATION SPLIT</span></div>
 <div class="metrics">${metric('Held-out reconstruction R²',fmt(d.r2,4))}${metric('Active on selected validation data',int(active)+' / 3,072')}${metric('Total representation vectors',int(d.n_vectors))}${metric('Mean activation frequency',pct(avg,3))}</div>
 <div class="twocol">${panel('Feature activation frequency','Binned distribution over all 3,072 features; each bar counts features, not biological conditions.',histogram(frequencies,{bins:15}))}${panel('Layer diagnostics','Statistics are computed from the saved validation feature arrays. Values refer to representation positions, not independently replicated conditions.',statsLine('Component',esc(cur.split('.')[0]))+statsLine('Features with zero firings',int(f.length-active))+statsLine('Maximum firing frequency',pct(Math.max(...frequencies),3))+statsLine('Median firing frequency',pct([...frequencies].sort((a,b)=>a-b)[1535],4))+statsLine('Mean L0 (expected TopK)',fmt(frequencies.reduce((a,b)=>a+b,0),3))+`<div class="callout info">For decoder features, open Feature Explorer to inspect the Top-20 training/validation gene rankings. The encoder and processor have no direct output-gene semantics.</div>`)}
 </div><div class="section-heading"><h2>Feature statistics</h2><small>First ${ranked.length} rows according to selected order</small></div>
 ${table(['Feature ID','Firing frequency','Mean active amplitude','Firing count','Inspect'],ranked.map(x=>row([`<span class="signal-value mono">${x.feature_id}</span>`,pct(x.activation_frequency,4),fmt(x.mean_active_amplitude,4),int(x.fire_count),`<a href="${featureHref(cur,x.feature_id)}">Open ↗</a>`])))}
 <p class="data-credit">Source: local 14-layer validation feature statistics exported from the fixed SAE series; ${linkNotebook('04')}. Original arrays were not modified.</p>`;
}
function subsetView(slice,genes,layer){
 if(!slice.length)return `<p class="footnote">No genes ranked in this split for this feature.</p>`;
 const max=Math.max(...slice.map(v=>v[1]),1);
 return `<div class="liststack">${slice.map((t,i)=>`<div style="display:grid;grid-template-columns:27px 1fr 90px;align-items:center;gap:12px"><span class="tiny mono">${String(i+1).padStart(2,'0')}</span><div><a href="${geneHref(genes[t[0]])}" style="font-size:12px;text-decoration:none;font-weight:700">${esc(genes[t[0]])}</a><div class="bar-track" style="margin-top:5px;height:5px"><div class="bar-inner" style="width:${Math.min(100,t[1]/max*100).toFixed(3)}%"></div></div></div><span style="text-align:right;color:var(--muted);font-size:11px" title="${int(t[2])} supporting conditions">${fmt(t[1],2)}<br/><small>${int(t[2])} conditions</small></span></div>`).join('')}</div>`;
}
async function featurePage(){
 const layer=ui.layer,feature=Number(ui.feature);
 const l=await load('layers/'+layer.replace('.','_')+'.json');const f=l.features[feature];
 if(!f)return `<div class="error-banner">Feature ID must be between 0 and 3,071.</div>`;
 const isDecoder=layer.startsWith('decoder.');let detail=null,genes=null,terms=null;
 if(isDecoder){[detail,genes,terms]=await Promise.all([load(`decoder/${layer.replace('.','_')}.json`),load('genes.json'),load('terms.json')]);}
 const fr=detail?.features[feature];const topTrain=fr?fr[5]:[],topValidation=fr?fr[6]:[],an=fr?fr[7]:[];
 const labels=layerNames.map(n=>[n,n]);
 return pageHead('FEATURE EXPLORER',`Feature <span class="title-accent">${int(feature)}</span>`,'Inspect a single sparse dictionary direction. Decoder features include historical functional annotations and separate training and validation Top-20 gene rankings; other blocks offer activation statistics only.')+
 `<div class="controls">${field('layer','PIE block',labels,layer)}${textField('featureId','Feature ID',feature,'0–3071')}<button class="action-btn" data-action="go-feature">Open feature</button><span class="chip ${isDecoder?'green':'purple'}">${isDecoder?'DECODER · GENE-ALIGNED':'LATENT-SLOT REPRESENTATIONS'}</span></div>
 <div class="metrics">${metric('Validation firing frequency',pct(f.activation_frequency,4))}${metric('Validation positive firings',int(f.fire_count))}${metric('Mean amplitude when firing',fmt(f.mean_active_amplitude,4))}${metric('Historical annotations',isDecoder?int(an.length):'N/A')}</div>
 ${isDecoder?`<div class="callout"><strong>Historical annotations (04).</strong> Card q-values use per-layer BH across 4,422 terms including TRRUST. Corrected 06b uses global BH across four layers and 4,421 ontology terms. <a href="https://github.com/vs-dx/pie-sae/blob/main/reports/Case_Normalization_Correction.md" target="_blank" rel="noopener noreferrer">Correction record ↗</a></div>`:''}
 ${isDecoder?`<div class="twocol" style="margin-top:18px">${panel('Training Top-20 associated genes','Ranked by mean positive feature activation conditional on firing; support counts are out of 128 selected training conditions.',subsetView(topTrain,genes,layer))}${panel('Validation Top-20 associated genes','The same fixed SAE weights evaluated on 32 selected validation conditions. Train and validation ranks may differ.',subsetView(topValidation,genes,layer))}</div>
 <div class="section-heading"><h2>Historical functional annotations</h2><small>${int(an.length)} records · fixed 04 analysis</small></div>
 <div class="controls">${field('annotationLimit','Display annotations',[['12','First 12'],['40','First 40'],['100','First 100'],['9999','All']],String(ui.annotationLimit))}</div>
 ${an.length?table(['Source','Term','q-value (historical)','Overlap','Inspect'],[...an].sort((a,b)=>a[2]-b[2]).slice(0,ui.annotationLimit).map(a=>row([`<span class="chip amber">${esc(terms[a[0]][0])}</span>`,esc(terms[a[0]][1]),short(a[2]),int(a[3]),`<a href="#pathways?term=${a[0]}">Term ↗</a>`]))):`<div class="empty">No historically significant annotation was saved for this feature under the original 04 protocol. This does not imply the feature has no biological role.</div>`}
 <div class="section-heading"><h2>Related analysis</h2></div><div class="twocol">${panel('Feature provenance','Dictionary is fixed; gene lists are based on the originally selected conditions.',statsLine('Source layer',esc(layer))+statsLine('Feature index',int(feature))+statsLine('Training firing frequency',pct(fr[1],4))+statsLine('Validation firing frequency',pct(fr[2],4))+statsLine('Training firing count',int(fr[3]))+statsLine('Validation firing count',int(fr[4])))}${panel('Interpretation limits','A rank is not a causal attribution.',`<p class="footnote">Gene-specific feature firing can reflect persistent gene identity. Apparent functional enrichment depends on gene set definitions and tested hypothesis families.</p><div class="btnrow"><a class="muted-btn" href="#cross-layer">Explore layer matching →</a>${[0,189,2120].includes(feature)&&layer==='decoder.0'?'<a class="muted-btn" href="#interventions">Intervention synthesis →</a>':''}</div>`)}
 </div>`:panel('Feature-level validation statistics','This encoder/processor block does not have direct output-gene annotations; its 512 positions per condition are latent slots.',statsLine('Representation type','Latent-slot vectors')+statsLine('Firing count',int(f.fire_count))+statsLine('Frequency',pct(f.activation_frequency,4))+statsLine('Mean positive amplitude',fmt(f.mean_active_amplitude,4)))}
 <p class="data-credit">Sources: ${linkNotebook('03')}, ${linkNotebook('04')}, ${linkNotebook('06b')}. No independent re-training or biological replication.</p>`;
}
async function genePage(){
 const [idx,genes]=await Promise.all([load('gene_index.json'),load('genes.json')]);
 const q=ui.query.trim().toUpperCase();const suggestions=genes.filter(x=>!q||x.toUpperCase().includes(q)).slice(0,70);
 const key=genes.find(x=>x.toUpperCase()===ui.gene.toUpperCase())||'';
 const hits=key?(idx.genes[key]||[]):[];
 const filtered=hits.filter(x=>ui.geneFilter==='all'||`decoder.${x[0]}`===ui.geneFilter).slice().sort((a,b)=>a[2]-b[2]||a[0]-b[0]);
 return pageHead('GENE SEARCH',`Search <span class="title-accent">output genes</span>`,'Locate which decoder SAE features include a gene among their training Top-20 gene associations. This is a gene-to-feature index, not a perturbation response assay.')+
 `<div class="horizontal-grid"><div class="panel">${textField('geneQuery','Gene symbol',ui.query,'e.g. PPCS, TP53','grow')}<p class="block-sub" style="margin-top:13px">${int(genes.length)} searchable output genes</p><div class="search-list" id="gene-suggestions">${suggestions.map(g=>`<button class="search-item ${g===key?'active':''}" data-action="gene-select" data-gene="${esc(g)}"><b>${esc(g)}</b><span class="count-pill">${int((idx.genes[g]||[]).length)} matches</span></button>`).join('')}</div></div>
 <div class="stack"><div class="panel"><div class="feature-jump"><div><span class="section-kicker">SELECTED GENE</span><h2 class="block-title" style="font-size:27px;margin-top:6px">${key?esc(key):'Select a gene'}</h2><p class="block-sub">${int(hits.length)} training Top-20 feature memberships across 4 decoder layers</p></div>${key?`<span class="chip green">GENE / FEATURE LOOKUP</span>`:''}</div>${field('geneFilter','Decoder block',[['all','All decoder blocks'],...decoderNames.map(n=>[n,n])],ui.geneFilter)}</div>
 ${key?table(['SAE block','Feature','Top-20 rank','Inspect'],filtered.slice(0,160).map(x=>row([`decoder.${x[0]}`,featureLink(`decoder.${x[0]}`,x[1]),`#${x[2]}`,`<a href="${featureHref(`decoder.${x[0]}`,x[1])}">Open ↗</a>`]))):'<div class="empty">Select a gene on the left to inspect its feature associations.</div>'}
 ${filtered.length>160?`<p class="footnote">Showing first 160 of ${int(filtered.length)} matches by rank.</p>`:''}
 </div></div><div class="callout info"><strong>Scope.</strong> This index covers the <em>training Top-20 gene lists</em> of decoder features only. A feature not returned for a gene might still respond at that position outside its 20 highest-ranked genes. ${linkNotebook('04')}.</div>`;
}
async function pathwayPage(){
 const [terms,index]=await Promise.all([load('terms.json'),load('pathway_index.json')]);
 const q=ui.termQuery.trim().toLowerCase();const found=terms.map((x,i)=>({x,i})).filter(z=>(z.x[1].toLowerCase().includes(q)||z.x[0].toLowerCase().includes(q))).sort((a,b)=>b.x[2]-a.x[2]).slice(0,65);
 const tid=ui.termId,term=terms[tid]||null;
 const entries=term?index[tid]||[]:[];
 const filtered=entries.filter(x=>ui.termFilter==='all'||`decoder.${x[0]}`===ui.termFilter).sort((a,b)=>a[2]-b[2]);
 return pageHead('PATHWAY SEARCH',`Browse <span class="title-accent">historical annotations</span>`,'Search GO Biological Process, KEGG, Reactome and curated TF terms from the saved notebook-04 annotation outputs. These are historical annotations, not independent 06b-confirmed biological functions.')+
 `<div class="horizontal-grid"><div class="panel">${textField('termQuery','Functional term',ui.termQuery,'Search pathway or GO term','grow')}<p class="block-sub" style="margin-top:12px">${int(terms.length)} terms represented among saved historical significant associations</p><div class="search-list" id="term-suggestions">${found.map(z=>`<button class="search-item ${z.i===tid?'active':''}" data-action="term-select" data-term="${z.i}"><span><b>${esc(z.x[1])}</b><br/><small>${esc(z.x[0])}</small></span><span class="count-pill">${int(z.x[2])}</span></button>`).join('')}</div></div>
 <div class="stack"><div class="panel">${term?`<span class="chip amber">${esc(term[0])} · HISTORICAL</span><h2 class="block-title" style="font-size:21px;margin-top:12px">${esc(term[1])}</h2><p class="block-sub">${int(entries.length)} feature–term associations in the original annotation outputs</p>`:'<div class="empty">Select a term to examine its associated SAE features.</div>'}${field('termFilter','Decoder block',[['all','All decoder blocks'],...decoderNames.map(n=>[n,n])],ui.termFilter)}</div>
 ${term?table(['Layer','Feature','Historical q-value','Overlap genes'],filtered.slice(0,180).map(e=>row([`decoder.${e[0]}`,featureLink(`decoder.${e[0]}`,e[1]),short(e[2]),int(e[3])]))):''}
 ${filtered.length>180?`<p class="footnote">Showing first 180 of ${int(filtered.length)} associated features.</p>`:''}</div></div>
 <div class="callout"><strong>Historical pathway search (04).</strong> The table displays original per-layer q-values for GO BP, KEGG, Reactome and TRRUST. The corrected 06b analysis has a separate global family. <a href="https://github.com/vs-dx/pie-sae/blob/main/reports/Case_Normalization_Correction.md" target="_blank" rel="noopener noreferrer">Correction record ↗</a></div>`;
}
async function crossLayerPage(){
 const o=await load('overview.json');const file=ui.cross;
 const d=await load(`cross_layer/${file}.json`);
 let pairs=d.pairs;const query=ui.crossQuery.trim();
 if(query){const q=Number(query);pairs=pairs.filter(x=>Number.isFinite(q)&&(x.feature_a===q||x.feature_b===q))}
 pairs=[...pairs].sort((a,b)=>b.signed_cosine-a.signed_cosine);
 const all=d.pairs,cos=all.map(x=>x.signed_cosine);const cosmean=cos.reduce((a,b)=>a+b,0)/cos.length;
 const transitionNames=o.cross_layer_transitions.map(x=>{const n=x.layer_a.replace('.','_')+'_to_'+x.layer_b.replace('.','_');return [n,`${x.layer_a} → ${x.layer_b}`]});
 const notable= d.layer_a==='decoder.0'&&d.layer_b==='decoder.1' ? all.find(x=>x.feature_a===189&&x.feature_b===413):null;
 return pageHead('CROSS-LAYER FLOW',`Adjacent <span class="title-accent">dictionary matches</span>`,'Explore signed-cosine mutual nearest-neighbor correspondences between fixed SAE dictionaries. These are geometric matches, not experimentally established biological pathways.')+
 `<div class="controls">${field('cross','Adjacent transition',transitionNames,file)}${textField('crossQuery','Feature ID (either side)',ui.crossQuery,'Filter by feature ID')}<button class="muted-btn" data-action="clear-cross">Clear filter</button></div>
 <div class="metrics">${metric('Mutual nearest-neighbor matches',int(all.length))}${metric('Mean signed cosine',fmt(cosmean,4))}${metric('Highest signed cosine',fmt(Math.max(...cos),4))}${metric('Source + target feature spaces','3,072 × 3,072')}</div>
 <div class="twocol">${panel('Signed-cosine distribution','Binned cosine values for all saved mutual matches in this transition.',histogram(cos,{bins:16,color:'#9b96f8'}))}${panel('Transition methodology','The original matching identifies mutual nearest neighbors in signed decoder-direction cosine space.',statsLine('Source block',esc(d.layer_a))+statsLine('Target block',esc(d.layer_b))+statsLine('Number of pairs',int(all.length))+statsLine('Pair selection','Signed-cosine MNN')+`<div class="callout info">The co-activation and gene-overlap analyses use additional null models. Geometric similarity alone cannot demonstrate information flow or a causal feature chain.</div>`)}
 </div>
 ${notable?`<div class="callout" style="margin-bottom:18px"><strong>Featured pair 189 → 413.</strong> Decoder.0 feature 189 is paired with decoder.1 feature 413, with signed cosine ${fmt(notable.signed_cosine,4)}. Notebook 07 found a pronounced downstream response under upstream ablation, but stronger geometry-matched controls would be needed to claim unique mediation. <a href="#interventions">Review the intervention findings →</a></div>`:''}
 <div class="section-heading"><h2>Matched pairs</h2><small>${int(pairs.length)} displayed matches</small></div>
 ${pairs.length?table(['Source feature','Target feature','Signed cosine','Open source','Open target'],pairs.slice(0,180).map(x=>row([int(x.feature_a),int(x.feature_b),fmt(x.signed_cosine,5),`<a href="${featureHref(d.layer_a,x.feature_a)}">${esc(d.layer_a)} / ${x.feature_a} ↗</a>`,`<a href="${featureHref(d.layer_b,x.feature_b)}">${esc(d.layer_b)} / ${x.feature_b} ↗</a>`]))):'<div class="empty">No matching features have the entered ID in this transition.</div>'}
 ${pairs.length>180?'<p class="footnote">Showing the first 180 matches sorted by signed cosine. Enter a feature ID to filter.</p>':''}
 <p class="data-credit">Sources: ${linkNotebook('05')}, ${linkNotebook('06b')}, and saved layer-matching arrays. No downstream causal inference from this table alone.</p>`;
}
function interventionPage(){
 const heads=[['DE logits',3],['LFC',2],['Delta-p logits',3]];
 return pageHead('NOTEBOOK 07 · INTERVENTIONS',`Computational <span class="title-accent">effects</span>`,'Effects of subtracting selected decoded SAE contributions or patching selected feature coefficients in the fixed pretrained PIE model. This view summarizes saved findings; per-condition intervention tables were not part of the exported atlas dataset.')+
 `<div class="metrics">${metric('Validation conditions analyzed','32')}${metric('Primary upstream feature','189')}${metric('Donor–recipient patch pairs','9')}${metric('Chosen downstream feature','413')}</div>
 <div class="twocol">${panel('Feature 189 · ablation','An intervention at decoder.0 produced nonzero output changes across 32 selected conditions.',statsLine('DE logits: ratio to fixed random controls','≈ 1.21×')+statsLine('LFC: ratio to fixed random controls','≈ 1.30×')+statsLine('Delta-p logits: ratio to fixed random controls','≈ 1.36×')+statsLine('Advantage across conditions','17–19 / 32')+`<div class="callout info">Ratios compare means across the same three output heads, not between heads. The feature was not consistently better than the learned-direction controls.</div>`)}
 ${panel('Activation patching · negative evidence','How often patching a feature-189 coefficient made a recipient prediction closer to a donor, for nine directed pairs.',heads.map(x=>`<div class="statline"><span>${x[0]}</span><b class="text-red">${x[1]} / 9 improved</b></div><div class="bar-track" style="margin:1px 0 11px"><div class="bar-inner" style="background:#ee9897;width:${x[1]/9*100}%"></div></div>`).join('')+`<p class="footnote">Patching did not reliably transfer donor predictions. Alternative donor and reverse-patching controls further limited the interpretation.</p>`)}
 </div>
 <div class="section-heading"><h2>Downstream decoder response</h2><span class="chip amber">NOT UNIQUE MEDIATION</span></div>
 <div class="panel"><div class="twocol"><div><div class="section-kicker">DECODER.0 / 189 → DECODER.1 / 413</div><div class="big-number" style="margin:14px 0">Rank #1 of 3,072</div><p class="block-sub">After feature-189 ablation, feature 413 ranked first by mean absolute activation response among all decoder.1 SAE features in six selected validation conditions. The statistic was computed on source-feature-active gene positions.</p><div class="tagrow"><span class="chip green">6 SELECTED CONDITIONS</span><span class="chip purple">COSINE ≈ 0.3023</span></div></div><div><div class="callout"><strong>Why this is not proof of a causal feature pathway.</strong> The chosen 189 direction was much more aligned geometrically with feature 413 than the sixteen learned-direction controls. No cosine-matched control or direct downstream mediation experiment was performed. TopK competition can also affect activation responses.</div><p class="footnote">The intervention subtracts a decoded feature component from the original hidden state; it does not enforce an exact zero latent coefficient after re-encoding.</p></div></div></div>
 ${disclaimer}<p class="data-credit">All findings on this page summarize the saved conclusions of ${linkNotebook('07')}. No per-condition response table is included in the current atlas export.</p>`;
}
async function comparisonsPage(){
 const o=await load('overview.json'),arr=o.layers;const group={encoder:[],processor:[],decoder:[]};for(let d of arr)group[d.layer.split('.')[0]].push(d.r2);
 const mean=vals=>vals.reduce((a,b)=>a+b,0)/vals.length;
 const vgroup=Object.entries(group).map(([name,vals])=>({name,value:mean(vals),label:fmt(mean(vals),4),group:name}));
 const transition=o.cross_layer_transitions;
 const rec=barrows(arr.map(x=>({name:x.layer,value:x.r2,label:fmt(x.r2,4),group:x.layer.split('.')[0]})),{max:1});
 const transitions=barrows(transition.map(x=>({name:x.layer_a+'→'+x.layer_b.split('.')[1],value:x.n_mutual_pairs,label:int(x.n_mutual_pairs),group:x.layer_a.split('.')[0]})),{max:3072});
 return pageHead('COMPARISONS',`Cross-block <span class="title-accent">measurements</span>`,'Aggregate comparisons use the original held-out reconstruction table and saved mutual-nearest-neighbor counts. Higher reconstruction or more matches do not independently establish interpretability.')+
 `<div class="threecol">${Object.entries(group).map(([name,values])=>`<div class="metric-card"><div class="meta">${name.toUpperCase()} BLOCKS</div><strong>${fmt(mean(values),4)}</strong><small>Mean validation reconstruction R² · ${values.length} blocks</small></div>`).join('')}</div>
 <div class="twocol" style="margin-top:18px">${panel('SAE reconstruction by block','Held-out R²; original numerical values rounded for display.',rec+'<div class="legend"><span><i></i>Encoder</span><span><i class="processor"></i>Processor</span><span><i class="decoder"></i>Decoder</span></div>')}
 ${panel('Adjacent-layer mutual matches','Count of mutual-nearest-neighbor matched feature pairs; maximum theoretically 3,072 per transition.',transitions)}
 </div><div class="callout info"><strong>Interpreting differences.</strong> Each block was fit once with its own fixed initialization seed. A greater number of mutual matches need not imply greater cosine alignment; decoder and latent-slot features also differ in their biological indexing.</div>
 <p class="data-credit">Sources: ${linkNotebook('03')}, ${linkNotebook('05')}, and compact source-verified tables. Group means calculated from the 14 published block R² values.</p>`;
}
function methodsPage(){
 return pageHead('METHODS & SOURCES',`Methodology, <span class="title-accent">with boundaries</span>`,'The atlas is a read-only interface to saved experimental summaries. It is not a new inference engine or an independent reproduction of the pretrained PIE model.')+
 `<div class="metrics">${metric('Experimental notebook stages','01–07')}${metric('Original representation dimension','768')}${metric('Dictionary features per block','3,072')}${metric('TopK sparsity','32')}</div>
 <div class="twocol">${panel('Data provenance','Each visualization refers to source-saved arrays, tables or experiment conclusions.',`<div class="liststack">${[
 ['Full research report',`${GH}/blob/main/reports/Research_Report.md`],
 ['Read the repository README',`${GH}/blob/main/README.md`],
 ['Published source notebooks',`${GH}/tree/main/notebooks`],
 ['CSV tables used for this Atlas',`${GH}/tree/c70dccff03d617187812db72ddf48a9e5062876f/reports/tables`],
 ['Archived research snapshot (c70dccf)','https://github.com/vs-dx/pie-sae/tree/c70dccff03d617187812db72ddf48a9e5062876f'],
 ['Pinned upstream PIE','https://github.com/ArcInstitute/pie/tree/21a526add9ff3a81b4f962d9b16bbae570052c63'],
 ['Pinned upstream bio-sae','https://github.com/Biodyn-AI/bio-sae/tree/3ac41b1e7725c224dc1225a3311e8802eb96da85']
 ].map(([t,u])=>`<div class="statline"><span>${link(u,t+' ↗')}</span><b>↗</b></div>`).join('')}</div>`)}
 ${panel('Export scope & limitations','Scientific boundaries are part of the data interface, not hidden in a footnote.',`<ul class="liststack" style="list-style:none;padding:0;font-size:12px;color:var(--sub)"><li>• 14 layer feature-statistics summaries on selected validation conditions.</li><li>• Training/validation Top-20 gene catalogs for four decoder blocks.</li><li>• Historical notebook-04 annotations with per-layer BH and the curated TRRUST term.</li><li>• 11 saved signed-cosine mutual-match lists.</li><li>• Search index over training Top-20 decoder-gene associations only.</li><li>• No large HDF5 arrays, checkpoint weights, or per-condition 06b/07 intervention tables.</li><li>• No independent clean-environment reproduction or biological validation.</li></ul>`)}
 </div><div class="section-heading"><h2>Original experimental notebooks</h2><small>Direct links to code and saved outputs</small></div>
 <div class="notebook-links">${Object.entries(NOTEBOOKS).map(([k,name])=>`<a href="${fileSource(k)}" target="_blank" rel="noopener noreferrer"><b>${k} · ${esc(name.replaceAll('_',' '))} ↗</b><small>Notebook code, methods, and stored outputs</small></a>`).join('')}</div>
 <div class="callout"><strong>Annotation reconciliation.</strong> Matching gene symbols in uppercase restored 340 explicit tests in 06b. The curated TRRUST term adds 3,439 explicit tests in 04. Together these changes reconcile the counts in all four decoder layers. <a href="https://github.com/vs-dx/pie-sae/blob/main/reports/Case_Normalization_Correction.md" target="_blank" rel="noopener noreferrer">Correction record ↗</a></div>
 <div class="callout info"><strong>Methods distinction.</strong> Association in a fixed dataset, sensitivity of a computational model to edited hidden states, and biological causality are three different levels of evidence. This atlas supports exploration of the first two, within the documented controls and limitations.</div>`;
}
async function render(){
 const seq=++renderSequence; const r=route();activePage=r.page;
 document.querySelector('#breadcrumb-current').textContent=names[pages.indexOf(activePage)];
 document.querySelectorAll('[data-nav]').forEach(x=>{const on=x.getAttribute('data-nav')===activePage;x.classList.toggle('active',on);if(on)x.setAttribute('aria-current','page');else x.removeAttribute('aria-current')});
 const p=r.p;
 if(lastAppliedHash!==location.hash){
   if(p.has('layer') && layerNames.includes(p.get('layer')))ui.layer=p.get('layer');
   if(p.has('feature')){const f=Number(p.get('feature'));if(Number.isInteger(f)&&f>=0&&f<3072)ui.feature=f}
   if(p.has('gene')){ui.gene=p.get('gene').trim().toUpperCase();ui.query=ui.gene}
   if(p.has('term')){let t=Number(p.get('term'));if(Number.isInteger(t)&&t>=0)ui.termId=t}
   if(p.has('cross'))ui.cross=p.get('cross');
   lastAppliedHash=location.hash;
 }
 main.innerHTML='<div class="loader">Loading research data <span class="dots">•••</span></div>';
 try{
 let html;
 switch(activePage){
  case 'overview':html=await overviewPage();break;
  case 'layers':html=await layerPage();break;
  case 'features':html=await featurePage();break;
  case 'genes':html=await genePage();break;
  case 'pathways':html=await pathwayPage();break;
  case 'cross-layer':html=await crossLayerPage();break;
  case 'interventions':html=interventionPage();break;
  case 'comparisons':html=await comparisonsPage();break;
  default:html=methodsPage();
 }
 if(seq!==renderSequence)return;
 main.innerHTML=html;
 document.querySelector('.sidebar').classList.remove('open');
 }catch(err){if(seq!==renderSequence)return;main.innerHTML=`<div class="error-banner"><b>Unable to display this research page.</b><p>${esc(err.message)}</p><p>Make sure the complete atlas <code>data/</code> directory is present and serve the folder over HTTP, not by opening the file directly from Explorer.</p></div>`;console.error(err)}
}
function debounce(fn,delay){let t;return (...args)=>{clearTimeout(t);t=setTimeout(()=>fn(...args),delay)}}
const onSearch=debounce(async(key,cursor)=>{
 await render();
 const input=document.querySelector(`[data-input="${key}"]`);
 if(input){input.focus();if(typeof cursor==='number')input.setSelectionRange(Math.min(cursor,input.value.length),Math.min(cursor,input.value.length));}
},240);
main.addEventListener('click',event=>{
 const btn=event.target.closest('[data-action]');if(!btn)return;
 const action=btn.dataset.action;
 if(action==='go-feature'){
  const x=document.querySelector('[data-input="featureId"]');let n=Number(x?.value);
  if(!x?.value?.trim()||!Number.isInteger(n)||n<0||n>=3072){x?.setCustomValidity('Enter a whole feature ID between 0 and 3071');x?.reportValidity();return}
  ui.feature=n;setHash('features',{layer:ui.layer,feature:n});return;
 }
 if(action==='gene-select'){ui.gene=btn.dataset.gene;ui.query=ui.gene;setHash('genes',{gene:ui.gene});return}
 if(action==='term-select'){ui.termId=Number(btn.dataset.term);setHash('pathways',{term:ui.termId});return}
 if(action==='clear-cross'){ui.crossQuery='';render();return}
});
main.addEventListener('change',e=>{
 const el=e.target;if(!(el instanceof HTMLElement))return;
 const key=el.dataset.control;if(!key)return;
 switch(key){
 case 'layer':ui.layer=el.value;setHash(activePage,activePage==='features'?{layer:ui.layer,feature:ui.feature}:{layer:ui.layer});break;
 case 'sort':ui.sort=el.value;render();break;
 case 'geneFilter':ui.geneFilter=el.value;render();break;
 case 'termFilter':ui.termFilter=el.value;render();break;
 case 'cross':ui.cross=el.value;setHash('cross-layer',{cross:ui.cross});break;
 case 'annotationLimit':ui.annotationLimit=Number(el.value);render();break;
 }
});
main.addEventListener('input',e=>{
 const el=e.target;const key=el?.dataset?.input;if(!key)return;
 if(key==='geneQuery'){ui.query=el.value;ui.gene='';onSearch(key,el.selectionStart)}
 if(key==='termQuery'){ui.termQuery=el.value;onSearch(key,el.selectionStart)}
 if(key==='crossQuery'){ui.crossQuery=el.value;onSearch(key,el.selectionStart)}
});
main.addEventListener('keydown',e=>{if(e.key==='Enter'&&e.target?.dataset?.input==='featureId'){e.preventDefault();document.querySelector('[data-action="go-feature"]')?.click()}});
document.querySelector('#mobile-toggle').addEventListener('click',()=>document.querySelector('.sidebar').classList.toggle('open'));
window.addEventListener('hashchange',()=>{window.scrollTo({top:0,behavior:'instant'});render()});
render();
