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

fn analyze_stats(
    content: &str,
    l_min: u32,
    method: AssemblyMethod,
) -> Result<(Analysis, Result<GraphLayoutInput, String>), String> {
    let (seqs, skipped) = parse_fasta_str(content);
    if seqs.is_empty() {
        return Err("No valid sequences found in FASTA".into());
    }
    let max_len = seqs.iter().map(|s| s.len()).max().unwrap();
    if max_len > LIMBS * 32 {
        return Err(format!(
            "Sequence length {} exceeds capacity (max {}bp).",
            max_len,
            LIMBS * 32
        ));
    }
    let min_len = seqs.iter().map(|s| s.len()).min().unwrap();
    if (l_min as usize) > min_len {
        return Err(format!(
            "min_overlap ({}) exceeds shortest sequence length ({}).",
            l_min, min_len
        ));
    }

    let seq_refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
    let edges = build_overlap_graph::<LIMBS>(&seq_refs, l_min, method);
    let contigs = assemble_contigs(&seq_refs, &edges);

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
            .map(|e| EdgeView {
                from: e.from_id,
                to: e.to_id,
                overlap: e.overlap_len,
                kind: edge_kind_from(e.from_strand, e.to_strand),
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

    let file_input_ref: NodeRef<leptos::html::Input> = NodeRef::new();

    let on_analyze = move |_| {
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
                Ok(text) => match analyze_stats(&text, l, m) {
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

                                let seed = layout::path_seeded_positions(
                                    connected_n,
                                    &compact_paths,
                                );
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
                                let iso_base_y = if max_y.is_finite() { max_y + 100.0 } else { 0.0 };
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
        let target: HtmlInputElement = ev.target().unwrap().dyn_into().unwrap();
        set_method.set(target.value());
    };

    view! {
        <main>
            <h1>"OliGraph"</h1>
            <p class="lede">"Upload a FASTA of oligos. Get overlap graph stats and assembled contigs."</p>

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
                    <button on:click=on_analyze prop:disabled=move || busy.get()>
                        {move || if busy.get() { "Analyzing…" } else { "Analyze" }}
                    </button>
                </div>
                {move || error.get().map(|e| view! { <p class="err">{e}</p> })}
                {move || filename.get().map(|n| view! {
                    <p style="color:var(--muted); margin: 0.5rem 0 0; font-size: 0.85rem;">
                        "File: " <code>{n}</code>
                    </p>
                })}
            </section>

            {move || analysis.get().map(|a| view! { <Results analysis=a /> })}

            {move || layout_busy.get().then(|| view! {
                <section class="panel">
                    <p style="color:var(--muted); margin:0;">"Computing graph layout…"</p>
                </section>
            })}

            {move || graph.get().map(|g| view! { <GraphView graph=g /> })}

            {move || graph_skip.get().map(|r| view! {
                <section class="panel">
                    <p style="color:var(--muted); margin:0;">{r}</p>
                </section>
            })}
        </main>
    }
}

#[component]
fn Results(analysis: Analysis) -> impl IntoView {
    let stat_grid = view! {
        <div class="stats">
            <div class="stat"><div class="k">"input sequences"</div><div class="v">{analysis.input_count}</div></div>
            <div class="stat"><div class="k">"skipped (non-ACGT)"</div><div class="v">{analysis.skipped}</div></div>
            <div class="stat"><div class="k">"length range"</div><div class="v">{format!("{}–{} bp", analysis.len_min, analysis.len_max)}</div></div>
            <div class="stat"><div class="k">"overlap edges"</div><div class="v">{analysis.edge_count}</div></div>
            <div class="stat"><div class="k">"isolated oligos"</div><div class="v">{analysis.isolated_count}</div></div>
            <div class="stat"><div class="k">"contigs"</div><div class="v">{analysis.contigs.len()}</div></div>
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
            view! {
                <tr>
                    <td class="num">{format!("contig_{}", c.index)}</td>
                    <td class="num">{c.component}</td>
                    <td class="num">{c.oligos}</td>
                    <td class="num">{c.length}</td>
                    <td>{topo}</td>
                    <td class="num">{c.branches}</td>
                    <td>
                        <a
                            href=href
                            download=dl_name
                            style="color:var(--accent); font-size:0.85rem;"
                        >"FASTA"</a>
                    </td>
                </tr>
            }
        })
        .collect();

    view! {
        <section class="panel">
            <h2 style="margin-top:0; font-size:1.1rem;">"Summary"</h2>
            {stat_grid}
        </section>
        <section class="panel">
            <h2 style="margin-top:0; font-size:1.1rem;">"Contigs"</h2>
            {if analysis.contigs.is_empty() {
                view! { <p style="color:var(--muted); margin:0;">"No contigs assembled."</p> }.into_any()
            } else {
                view! {
                    <table>
                        <thead>
                            <tr>
                                <th>"id"</th>
                                <th>"component"</th>
                                <th>"oligos"</th>
                                <th>"length (bp)"</th>
                                <th>"topology"</th>
                                <th>"branches"</th>
                                <th>"download"</th>
                            </tr>
                        </thead>
                        <tbody>{rows}</tbody>
                    </table>
                }.into_any()
            }}
        </section>
    }
}

fn main() {
    console_error_panic_hook::set_once();
    leptos::mount::mount_to_body(App);
}
