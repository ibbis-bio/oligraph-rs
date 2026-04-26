mod graph_view;
mod layout;

use graph_view::{EdgeKind, EdgeView, GraphData, GraphView};
use leptos::prelude::*;
use leptos::task::spawn_local;
use oligraph_rs::{
    AssemblyMethod, Contig, LIMBS, Strand, Topology, assemble_contigs, build_overlap_graph,
    parse_fasta_str,
};
use wasm_bindgen::JsCast;
use web_sys::HtmlInputElement;

const MAX_GRAPH_NODES: usize = 5000;

#[derive(Clone)]
struct Analysis {
    input_count: usize,
    skipped: u32,
    len_min: usize,
    len_max: usize,
    edge_count: usize,
    isolated_count: usize,
    contigs: Vec<ContigSummary>,
    graph: Option<GraphData>,
    graph_skipped_reason: Option<String>,
}

#[derive(Clone)]
struct ContigSummary {
    index: usize,
    component: usize,
    oligos: usize,
    length: usize,
    topology: Topology,
    branches: u32,
}

fn analyze(content: &str, l_min: u32, method: AssemblyMethod) -> Result<Analysis, String> {
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
        .map(|(i, c): (usize, &Contig)| ContigSummary {
            index: i,
            component: c.component,
            oligos: c.path.len(),
            length: c.sequence.len(),
            topology: c.topology,
            branches: c.branches,
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

    let (graph, graph_skipped_reason) = if n > MAX_GRAPH_NODES {
        (
            None,
            Some(format!(
                "Graph view skipped: {} sequences exceeds limit of {}.",
                n, MAX_GRAPH_NODES
            )),
        )
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
        let seed = layout::path_seeded_positions(n, &contig_paths);
        let positions = layout::fruchterman_reingold(n, &edge_pairs, 60, Some(seed));

        let seq_lengths: Vec<usize> = seqs.iter().map(|s| s.len()).collect();

        (
            Some(GraphData {
                n_nodes: n,
                seq_lengths,
                edges: edge_views,
                components: node_components,
                n_components,
                initial_positions: positions,
            }),
            None,
        )
    };

    Ok(Analysis {
        input_count: n,
        skipped,
        len_min: min_len,
        len_max: max_len,
        edge_count: edges.len(),
        isolated_count,
        contigs: summaries,
        graph,
        graph_skipped_reason,
    })
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

        let l = l_min.get_untracked();
        let m = match method.get_untracked().as_str() {
            "pca" => AssemblyMethod::Pca,
            _ => AssemblyMethod::All,
        };

        let blob: gloo_file::Blob = file.into();

        spawn_local(async move {
            match gloo_file::futures::read_as_text(&blob).await {
                Ok(text) => match analyze(&text, l, m) {
                    Ok(a) => set_analysis.set(Some(a)),
                    Err(e) => set_error.set(Some(e)),
                },
                Err(e) => set_error.set(Some(format!("File read failed: {}", e))),
            }
            set_busy.set(false);
        });
    };

    let on_l_min_change = move |ev: leptos::ev::Event| {
        let target: HtmlInputElement = ev.target().unwrap().dyn_into().unwrap();
        if let Ok(v) = target.value().parse::<u32>() {
            set_l_min.set(v.clamp(1, 32));
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
                            max="32"
                            prop:value=move || l_min.get().to_string()
                            on:change=on_l_min_change
                        />
                    </div>
                    <div class="field">
                        <label>"assembly method"</label>
                        <select on:change=on_method_change prop:value=move || method.get()>
                            <option value="all">"all"</option>
                            <option value="pca">"pca (3'-3' only)"</option>
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

            {move || analysis.get().map(|a| {
                let graph = a.graph.clone();
                let skipped_reason = a.graph_skipped_reason.clone();
                view! {
                    <Results analysis=a />
                    {graph.map(|g| view! { <GraphView graph=g /> })}
                    {skipped_reason.map(|r| view! {
                        <section class="panel"><p style="color:var(--muted); margin:0;">{r}</p></section>
                    })}
                }
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
            view! {
                <tr>
                    <td class="num">{format!("contig_{}", c.index)}</td>
                    <td class="num">{c.component}</td>
                    <td class="num">{c.oligos}</td>
                    <td class="num">{c.length}</td>
                    <td>{topo}</td>
                    <td class="num">{c.branches}</td>
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
