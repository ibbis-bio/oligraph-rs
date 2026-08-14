mod graph_view;
mod layout;

use gloo_timers::future::TimeoutFuture;
use graph_view::{EdgeKind, EdgeView, GraphData, GraphView};
use leptos::prelude::*;
use leptos::task::spawn_local;
use oligraph_rs::{
    assemble_contigs, build_overlap_graph, parse_fasta_str, AssemblyMethod, Strand, Topology, LIMBS,
};
use wasm_bindgen::JsCast;
use web_sys::HtmlInputElement;

const MAX_GRAPH_NODES: usize = 10000;

#[derive(Clone)]
struct Analysis {
    input_count: usize,
    skipped: u32,
    len_min: usize,
    len_max: usize,
    edge_count: usize,
    isolated_count: usize,
    contigs: Vec<ContigSummary>,
}

#[derive(Clone)]
struct ContigSummary {
    index: usize,
    component: usize,
    oligos: usize,
    length: usize,
    topology: Topology,
    branches: u32,
    fasta: String,
}

struct GraphLayoutInput {
    n: usize,
    seq_lengths: Vec<usize>,
    edge_views: Vec<EdgeView>,
    node_components: Vec<Option<usize>>,
    n_components: usize,
    edge_pairs: Vec<(u32, u32)>,
    contig_paths: Vec<Vec<u32>>,
    has_edge: Vec<bool>,
}

fn analyse_stats(
    content: &str,
    l_min: u32,
    method: AssemblyMethod,
) -> Result<(Analysis, Result<GraphLayoutInput, String>), String> {
    let (seqs, skipped) = parse_fasta_str(content);
    if seqs.is_empty() {
        return Err("No valid sequences found in FASTA".into());
    }
    let max_len = seqs.iter().map(|s| s.len()).max().unwrap();
    let min_len = seqs.iter().map(|s| s.len()).min().unwrap();
    if (l_min as usize) > min_len {
        return Err(format!(
            "min_overlap ({}) exceeds shortest sequence length ({}).",
            l_min, min_len
        ));
    }

    let seq_refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
    let edges =
        build_overlap_graph::<LIMBS>(&seq_refs, l_min, method).map_err(|e| e.to_string())?;
    let contigs = assemble_contigs(&seq_refs, &edges).map_err(|e| e.to_string())?;

    let mut has_edge = vec![false; seqs.len()];
    for e in &edges {
        if e.from_id != e.to_id {
            has_edge[e.from_id as usize] = true;
            has_edge[e.to_id as usize] = true;
        }
    }
    let isolated_count = has_edge.iter().filter(|&&v| !v).count();

    let summaries: Vec<ContigSummary> = contigs
        .iter()
        .enumerate()
        .map(|(i, c)| {
            let header = format!(
                ">contig_{} component={} oligos={} length={}",
                i,
                c.component,
                c.path.len(),
                c.sequence.len()
            );
            let seq_str = String::from_utf8_lossy(&c.sequence).into_owned();
            ContigSummary {
                index: i,
                component: c.component,
                oligos: c.path.len(),
                length: c.sequence.len(),
                topology: c.topology,
                branches: c.branches,
                fasta: format!("{}\n{}\n", header, seq_str),
            }
        })
        .collect();

    let n = seqs.len();
    let mut node_components: Vec<Option<usize>> = vec![None; n];
    let mut max_comp = 0usize;
    for c in &contigs {
        if c.component > max_comp {
            max_comp = c.component;
        }
        for &(id, _) in &c.path {
            let i = id as usize;
            if i < n {
                node_components[i] = Some(c.component);
            }
        }
    }
    let n_components = if contigs.is_empty() { 0 } else { max_comp + 1 };

    let graph_result = if n > MAX_GRAPH_NODES {
        Err(format!(
            "Graph view skipped: {} sequences exceeds limit of {}.",
            n, MAX_GRAPH_NODES
        ))
    } else {
        let edge_views: Vec<EdgeView> = edges
            .iter()
            .filter_map(|e| {
                let kind = edge_kind_from(e.from_strand, e.to_strand);
                // PCA mode must never show Fwd→Fwd edges; filter defensively in
                // case any slip through the backend filter.
                if method == AssemblyMethod::Pca && kind == EdgeKind::FwdFwd {
                    return None;
                }
                Some(EdgeView {
                    from: e.from_id,
                    to: e.to_id,
                    overlap: e.overlap_len,
                    kind,
                })
            })
            .collect();
        let edge_pairs: Vec<(u32, u32)> = edges.iter().map(|e| (e.from_id, e.to_id)).collect();
        let contig_paths: Vec<Vec<u32>> = contigs
            .iter()
            .map(|c| c.path.iter().map(|&(id, _)| id).collect())
            .collect();
        let seq_lengths: Vec<usize> = seqs.iter().map(|s| s.len()).collect();
        Ok(GraphLayoutInput {
            n,
            seq_lengths,
            edge_views,
            node_components,
            n_components,
            edge_pairs,
            contig_paths,
            has_edge,
        })
    };

    Ok((
        Analysis {
            input_count: n,
            skipped,
            len_min: min_len,
            len_max: max_len,
            edge_count: edges.len(),
            isolated_count,
            contigs: summaries,
        },
        graph_result,
    ))
}

fn edge_kind_from(from: Strand, to: Strand) -> EdgeKind {
    match (from, to) {
        (Strand::Fwd, Strand::Fwd) => EdgeKind::FwdFwd,
        (Strand::Fwd, Strand::Rev) => EdgeKind::FwdRev,
        (Strand::Rev, Strand::Fwd) => EdgeKind::RevFwd,
        (Strand::Rev, Strand::Rev) => EdgeKind::RevRev,
    }
}

#[component]
fn App() -> impl IntoView {
    let (l_min, set_l_min) = signal(20u32);
    let (method, set_method) = signal("all".to_string());
    let (busy, set_busy) = signal(false);
    let (error, set_error) = signal::<Option<String>>(None);
    let (analysis, set_analysis) = signal::<Option<Analysis>>(None);
    let (graph, set_graph) = signal::<Option<GraphData>>(None);
    let (graph_skip, set_graph_skip) = signal::<Option<String>>(None);
    let (layout_busy, set_layout_busy) = signal(false);
    let (filename, set_filename) = signal::<Option<String>>(None);
    let (highlight_comp, set_highlight_comp) = signal::<Option<usize>>(None);

    let file_input_ref: NodeRef<leptos::html::Input> = NodeRef::new();

    let on_analyse = move |_| {
        let Some(input) = file_input_ref.get() else {
            set_error.set(Some("File input not ready".into()));
            return;
        };
        let Some(files) = input.files() else {
            set_error.set(Some("Choose a FASTA file first".into()));
            return;
        };
        let Some(file) = files.get(0) else {
            set_error.set(Some("Choose a FASTA file first".into()));
            return;
        };

        set_filename.set(Some(file.name()));
        set_busy.set(true);
        set_error.set(None);
        set_analysis.set(None);
        set_graph.set(None);
        set_graph_skip.set(None);

        let l = l_min.get_untracked();
        let m = match method.get_untracked().as_str() {
            "pca" => AssemblyMethod::Pca,
            _ => AssemblyMethod::All,
        };

        let blob: gloo_file::Blob = file.into();

        spawn_local(async move {
            match gloo_file::futures::read_as_text(&blob).await {
                Ok(text) => match analyse_stats(&text, l, m) {
                    Ok((a, graph_result)) => {
                        set_analysis.set(Some(a));
                        match graph_result {
                            Ok(li) => {
                                set_layout_busy.set(true);
                                TimeoutFuture::new(0).await;

                                // Compact to connected-only nodes so isolated nodes
                                // don't participate in the O(n²) F-R repulsion.
                                let connected_indices: Vec<usize> =
                                    (0..li.n).filter(|&i| li.has_edge[i]).collect();
                                let connected_n = connected_indices.len();
                                let mut old_to_new = vec![None::<usize>; li.n];
                                for (ni, &oi) in connected_indices.iter().enumerate() {
                                    old_to_new[oi] = Some(ni);
                                }
                                let compact_edges: Vec<(u32, u32)> = li
                                    .edge_pairs
                                    .iter()
                                    .filter_map(|&(u, v)| {
                                        let nu = old_to_new[u as usize]?;
                                        let nv = old_to_new[v as usize]?;
                                        Some((nu as u32, nv as u32))
                                    })
                                    .collect();
                                let compact_paths: Vec<Vec<u32>> = li
                                    .contig_paths
                                    .iter()
                                    .map(|p| {
                                        p.iter()
                                            .filter_map(|&id| {
                                                old_to_new[id as usize].map(|ni| ni as u32)
                                            })
                                            .collect::<Vec<_>>()
                                    })
                                    .filter(|p| !p.is_empty())
                                    .collect();

                                let seed =
                                    layout::path_seeded_positions(connected_n, &compact_paths);
                                let compact_pos = layout::fruchterman_reingold(
                                    connected_n,
                                    &compact_edges,
                                    60,
                                    Some(seed),
                                );

                                // Expand back: place isolated nodes in a grid below
                                // the connected graph so they don't pile at the origin.
                                let mut positions = vec![(0.0_f32, 0.0_f32); li.n];
                                for (ni, &oi) in connected_indices.iter().enumerate() {
                                    positions[oi] = compact_pos[ni];
                                }
                                let max_y = compact_pos
                                    .iter()
                                    .map(|p| p.1)
                                    .fold(f32::NEG_INFINITY, f32::max);
                                let iso_base_y = if max_y.is_finite() {
                                    max_y + 100.0
                                } else {
                                    0.0
                                };
                                let mut iso_k = 0usize;
                                for i in 0..li.n {
                                    if !li.has_edge[i] {
                                        positions[i] = (
                                            (iso_k % 10) as f32 * 60.0,
                                            iso_base_y + (iso_k / 10) as f32 * 60.0,
                                        );
                                        iso_k += 1;
                                    }
                                }

                                let initial_viewbox = compact_paths
                                    .iter()
                                    .max_by_key(|p| p.len())
                                    .filter(|p| !p.is_empty())
                                    .map(|path| {
                                        let pts: Vec<(f32, f32)> = path
                                            .iter()
                                            .filter_map(|&id| compact_pos.get(id as usize).copied())
                                            .collect();
                                        layout::bounding_box(&pts)
                                    })
                                    .unwrap_or_else(|| layout::bounding_box(&compact_pos));

                                set_graph.set(Some(GraphData {
                                    n_nodes: li.n,
                                    seq_lengths: li.seq_lengths,
                                    edges: li.edge_views,
                                    components: li.node_components,
                                    n_components: li.n_components,
                                    initial_positions: positions,
                                    initial_viewbox,
                                    connected: li.has_edge,
                                }));
                                set_layout_busy.set(false);
                            }
                            Err(skip_reason) => {
                                set_graph_skip.set(Some(skip_reason));
                            }
                        }
                        set_busy.set(false);
                    }
                    Err(e) => {
                        set_error.set(Some(e));
                        set_busy.set(false);
                    }
                },
                Err(e) => {
                    set_error.set(Some(format!("File read failed: {}", e)));
                    set_busy.set(false);
                }
            }
        });
    };

    let on_l_min_change = move |ev: leptos::ev::Event| {
        let target: HtmlInputElement = ev.target().unwrap().dyn_into().unwrap();
        if let Ok(v) = target.value().parse::<u32>() {
            set_l_min.set(v.clamp(1, 64));
        }
    };

    let on_method_change = move |ev: leptos::ev::Event| {
        use wasm_bindgen::JsCast;
        if let Some(target) = ev.target() {
            if let Ok(sel) = target.dyn_into::<web_sys::HtmlSelectElement>() {
                set_method.set(sel.value());
            }
        }
    };

    view! {
        <main>
            <Logo />

            <section class="panel">
                <div class="row">
                    <div class="field">
                        <label>"FASTA file"</label>
                        <input
                            type="file"
                            accept=".fasta,.fa,.fna,.txt"
                            node_ref=file_input_ref
                        />
                    </div>
                    <div class="field">
                        <label>"min overlap (bp)"</label>
                        <input
                            type="number"
                            min="1"
                            max="64"
                            prop:value=move || l_min.get().to_string()
                            on:change=on_l_min_change
                        />
                    </div>
                    <div class="field">
                        <label>"assembly method"</label>
                        <select on:change=on_method_change prop:value=move || method.get()>
                            <option value="all">"all"</option>
                            <option value="pca">"pca"</option>
                        </select>
                    </div>
                    <button on:click=on_analyse prop:disabled=move || busy.get()>
                        {move || if busy.get() { "analysing…" } else { "analyse" }}
                    </button>
                </div>
                {move || error.get().map(|e| view! { <div class="err">{e}</div> })}
                {move || filename.get().map(|n| view! {
                    <div style="color:var(--muted); margin-top: 1rem; font-size: 0.8rem; display: flex; align-items: center; gap: 0.5rem;">
                        <span style="opacity: 0.6;">"ACTIVE FILE:"</span> <code>{n}</code>
                    </div>
                })}
            </section>

            {move || analysis.get().map(|a| view! { <Results analysis=a set_highlight_comp /> })}

            {move || layout_busy.get().then(|| view! {
                <section class="panel" style="display: flex; align-items: center; justify-content: center; padding: 2rem;">
                    <LoadingDots label="Computing graph layout" />
                </section>
            })}

            {move || graph.get().map(|g| view! { <GraphView graph=g highlight_comp /> })}

            {move || graph_skip.get().map(|r| view! {
                <section class="panel" style="border-left: 4px solid var(--warn);">
                    <p style="color:var(--muted); margin:0; font-size: 0.9rem;">
                        <strong>"Notice:"</strong> " " {r}
                    </p>
                </section>
            })}

            <footer style="margin-top: auto; padding: 4rem 0 2rem; border-top: 1px solid var(--border); text-align: center;">
                <div style="display: flex; align-items: center; justify-content: center; gap: 0.5rem; color: var(--muted); font-size: 0.85rem;">
                    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" style="opacity: 0.8;">
                        <rect x="3" y="11" width="18" height="11" rx="2" ry="2"></rect>
                        <path d="M7 11V7a5 5 0 0 1 10 0v4"></path>
                    </svg>
                    <span>"Privacy First: All sequence analysis and graph layout is performed locally in your browser. "</span>
                    <strong style="color: var(--text); font-weight: 600;">"No data is sent to any server."</strong>
                </div>
            </footer>
        </main>
    }
}

#[component]
fn Logo() -> impl IntoView {
    view! {
        <header style="margin-bottom: 3.5rem; display: flex; flex-direction: column; align-items: center; text-align: center;">
            <div style="display:flex; align-items:center; gap:1.25rem; margin-bottom:1.25rem;">
                <div class="logo-icon-container" style="position: relative; width: 48px; height: 48px; background: rgba(255,255,255,0.03); border-radius: 12px; display: flex; align-items: center; justify-content: center; border: 1px solid var(--border);">
                    <svg width="32" height="32" viewBox="0 0 44 44" style="overflow:visible;">
                        <g class="logo-helix">
                            <line x1="14" y1="10" x2="30" y2="10" stroke="currentColor" stroke-width="1.2" opacity="0.25"/>
                            <line x1="14" y1="22" x2="30" y2="22" stroke="currentColor" stroke-width="1.2" opacity="0.25"/>
                            <line x1="14" y1="34" x2="30" y2="34" stroke="currentColor" stroke-width="1.2" opacity="0.25"/>

                            <path d="M 16 4 C 34 4, 34 16, 22 16 S 10 28, 22 28 S 34 40, 16 40" fill="none" stroke="#4ade80" stroke-width="2.8" stroke-linecap="round"/>
                            <path d="M 28 4 C 10 4, 10 16, 22 16 S 34 28, 22 28 S 10 40, 28 40" fill="none" stroke="#60a5fa" stroke-width="2.8" stroke-linecap="round"/>
                        </g>
                        <g class="logo-graph">
                            <line class="logo-edge" x1="22" y1="5"  x2="5"  y2="22" stroke="currentColor" stroke-width="1.5" opacity="0.15"/>
                            <line class="logo-edge" x1="39" y1="22" x2="22" y2="39" stroke="currentColor" stroke-width="1.5" opacity="0.15"/>
                            <line class="logo-edge-flow" x1="22" y1="5"  x2="5"  y2="22" stroke="var(--accent)" stroke-width="1.5" stroke-linecap="round"/>
                            <line class="logo-edge-flow" x1="39" y1="22" x2="22" y2="39" stroke="var(--accent)" stroke-width="1.5" stroke-linecap="round"/>

                            <circle cx="22" cy="5"  r="5.5" fill="#4ade80" stroke="#0b0d11" stroke-width="2"/>
                            <circle cx="5"  cy="22" r="5.5" fill="#60a5fa" stroke="#0b0d11" stroke-width="2"/>
                            <circle cx="39" cy="22" r="5.5" fill="#fbbf24" stroke="#0b0d11" stroke-width="2"/>
                            <circle cx="22" cy="39" r="5.5" fill="#a78bfa" stroke="#0b0d11" stroke-width="2"/>
                        </g>
                    </svg>
                </div>
                <div>
                    <h1 class="logo-wordmark">"OliGraph"</h1>
                    <div style="font-size:0.95rem; color:var(--muted); margin-top:0.15rem; font-weight: 450;">
                        "Oligonucleotide overlap graph explorer"
                    </div>
                </div>
            </div>
        </header>
    }
}

#[component]
fn LoadingDots(label: &'static str) -> impl IntoView {
    view! {
        <div style="display: flex; align-items: center; gap: 0.75rem; color: var(--text); font-weight: 500;">
            {label}
            <div style="display: flex; gap: 0.25rem;">
                <span class="dot">"●"</span>
                <span class="dot" style="animation-delay: 0.2s">"●"</span>
                <span class="dot" style="animation-delay: 0.4s">"●"</span>
            </div>
        </div>
    }
}

#[component]
fn Results(analysis: Analysis, set_highlight_comp: WriteSignal<Option<usize>>) -> impl IntoView {
    let stat_grid = view! {
        <div class="stats">
            <div class="stat">
                <div class="k">"sequences"</div>
                <div class="v">{analysis.input_count}</div>
            </div>
            <div class="stat">
                <div class="k">"skipped"</div>
                <div class="v" style={move || if analysis.skipped > 0 { "color: var(--warn)" } else { "" }}>
                    {analysis.skipped}
                </div>
            </div>
            <div class="stat">
                <div class="k">"length range"</div>
                <div class="v" style="font-size: 1.1rem;">
                    {format!("{}–{} ", analysis.len_min, analysis.len_max)}
                    <span style="font-size: 0.8rem; font-weight: 500; opacity: 0.6;">"BP"</span>
                </div>
            </div>
            <div class="stat">
                <div class="k">"edges"</div>
                <div class="v">{analysis.edge_count}</div>
            </div>
            <div class="stat">
                <div class="k">"contigs"</div>
                <div class="v" style="color: var(--pass)">{analysis.contigs.len()}</div>
            </div>
        </div>
    };

    let rows: Vec<_> = analysis
        .contigs
        .iter()
        .map(|c| {
            let topo = match c.topology {
                Topology::Linear => "linear",
                Topology::Cyclic => "cyclic",
            };
            let href = format!(
                "data:text/plain;charset=utf-8,{}",
                js_sys::encode_uri_component(&c.fasta)
            );
            let dl_name = format!("contig_{}.fasta", c.index);
            let comp = c.component;
            view! {
                <tr
                    on:mouseenter=move |_| set_highlight_comp.set(Some(comp))
                    on:mouseleave=move |_| set_highlight_comp.set(None)
                    style="cursor: default;"
                >
                    <td class="num" style="color: var(--accent); font-weight: 600;">{format!("c_{}", c.index)}</td>
                    <td class="num">{c.component}</td>
                    <td class="num">{c.oligos}</td>
                    <td class="num">{c.length}</td>
                    <td style="text-transform: capitalize;">{topo}</td>
                    <td class="num">{c.branches}</td>
                    <td>
                        <a
                            href=href
                            download=dl_name
                            style="color:var(--accent); font-size:0.8rem; font-weight: 600; text-decoration: none; display: flex; align-items: center; gap: 0.25rem;"
                        >
                            <span>"↓"</span>
                            <span>"FASTA"</span>
                        </a>
                    </td>
                </tr>
            }
        })
        .collect();

    view! {
        <div style="display: grid; grid-template-columns: 1fr; gap: 1.5rem; margin-bottom: 1.5rem;">
            <section class="panel">
                <h2 style="margin-bottom: 1.25rem; display: flex; align-items: center; gap: 0.5rem;">
                    <span style="width: 8px; height: 8px; border-radius: 2px; background: var(--accent);"></span>
                    "Analysis Summary"
                </h2>
                {stat_grid}
            </section>

            <section class="panel" style="padding: 0; overflow: hidden;">
                <div style="padding: 1.25rem 1.5rem; border-bottom: 1px solid var(--border);">
                    <h2 style="display: flex; align-items: center; gap: 0.5rem;">
                        <span style="width: 8px; height: 8px; border-radius: 2px; background: var(--pass);"></span>
                        "Assembled Contigs"
                    </h2>
                </div>
                <div style="overflow-x: auto;">
                    {if analysis.contigs.is_empty() {
                        view! { <p style="color:var(--muted); padding: 2rem; margin:0; text-align: center;">"No contigs assembled."</p> }.into_any()
                    } else {
                        view! {
                            <table>
                                <thead>
                                    <tr>
                                        <th>"ID"</th>
                                        <th>"Comp"</th>
                                        <th>"Oligos"</th>
                                        <th>"Len (bp)"</th>
                                        <th>"Topology"</th>
                                        <th>"Br."</th>
                                        <th>"Download"</th>
                                    </tr>
                                </thead>
                                <tbody>{rows}</tbody>
                            </table>
                        }.into_any()
                    }}
                </div>
            </section>
        </div>
    }
}

fn main() {
    console_error_panic_hook::set_once();
    leptos::mount::mount_to_body(App);
}
