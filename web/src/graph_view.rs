use leptos::prelude::*;

use crate::layout;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum EdgeKind {
    FwdFwd,
    FwdRev,
    RevFwd,
    RevRev,
}

impl EdgeKind {
    fn color(self) -> &'static str {
        match self {
            EdgeKind::FwdFwd => "#4ade80",
            EdgeKind::FwdRev => "#6aa9ff",
            EdgeKind::RevFwd => "#fbbf24",
            EdgeKind::RevRev => "#a78bfa",
        }
    }
    pub fn label(self) -> &'static str {
        match self {
            EdgeKind::FwdFwd => "Fwd→Fwd",
            EdgeKind::FwdRev => "Fwd→Rev",
            EdgeKind::RevFwd => "Rev→Fwd",
            EdgeKind::RevRev => "Rev→Rev",
        }
    }
}

#[derive(Clone, Copy)]
pub struct EdgeView {
    pub from: u32,
    pub to: u32,
    pub overlap: u32,
    pub kind: EdgeKind,
}

#[derive(Clone)]
pub struct GraphData {
    pub n_nodes: usize,
    pub seq_lengths: Vec<usize>,
    pub edges: Vec<EdgeView>,
    pub components: Vec<Option<usize>>,
    pub n_components: usize,
    pub initial_positions: Vec<(f32, f32)>,
}

const COMPONENT_COLORS: &[&str] = &[
    "#5470c6", "#91cc75", "#fac858", "#ee6666", "#73c0de", "#3ba272", "#fc8452", "#9a60b4",
    "#ea7ccc", "#2f4b7c", "#a0d911", "#fa541c",
];

fn component_color(comp: Option<usize>) -> &'static str {
    match comp {
        None => "#6b7280",
        Some(i) => COMPONENT_COLORS[i % COMPONENT_COLORS.len()],
    }
}

#[derive(Clone, Copy)]
enum DragState {
    None,
    Node {
        id: usize,
        start_x: f32,
        start_y: f32,
        start_client_x: f32,
        start_client_y: f32,
    },
    Pan {
        start_vb_x: f32,
        start_vb_y: f32,
        start_client_x: f32,
        start_client_y: f32,
    },
}

#[component]
pub fn GraphView(graph: GraphData) -> impl IntoView {
    let n = graph.n_nodes;
    let n_components = graph.n_components;
    let positions = RwSignal::new(graph.initial_positions.clone());
    let bb = layout::bounding_box(&graph.initial_positions);
    let viewbox = RwSignal::new(bb);
    let drag = StoredValue::new(DragState::None);
    let svg_ref: NodeRef<leptos::svg::Svg> = NodeRef::new();

    let initial_positions = StoredValue::new(graph.initial_positions.clone());
    let edges = StoredValue::new(graph.edges);
    let components = StoredValue::new(graph.components);
    let seq_lengths = StoredValue::new(graph.seq_lengths);

    let svg_rect = move || -> Option<(f32, f32, f32, f32)> {
        svg_ref.get_untracked().map(|el| {
            let r = el.get_bounding_client_rect();
            (r.x() as f32, r.y() as f32, r.width() as f32, r.height() as f32)
        })
    };

    let on_bg_down = move |ev: web_sys::PointerEvent| {
        let (vx, vy, _, _) = viewbox.get_untracked();
        drag.set_value(DragState::Pan {
            start_vb_x: vx,
            start_vb_y: vy,
            start_client_x: ev.client_x() as f32,
            start_client_y: ev.client_y() as f32,
        });
    };

    let on_move = move |ev: web_sys::PointerEvent| {
        let cur = drag.get_value();
        if matches!(cur, DragState::None) {
            return;
        }
        let Some((_, _, sw, _sh)) = svg_rect() else {
            return;
        };
        let (_, _, vw, _) = viewbox.get_untracked();
        let scale = if sw > 1.0 { vw / sw } else { 1.0 };
        match cur {
            DragState::None => {}
            DragState::Pan {
                start_vb_x,
                start_vb_y,
                start_client_x,
                start_client_y,
            } => {
                let dx = (ev.client_x() as f32 - start_client_x) * scale;
                let dy = (ev.client_y() as f32 - start_client_y) * scale;
                viewbox.update(|vb| {
                    vb.0 = start_vb_x - dx;
                    vb.1 = start_vb_y - dy;
                });
            }
            DragState::Node {
                id,
                start_x,
                start_y,
                start_client_x,
                start_client_y,
            } => {
                let dx = (ev.client_x() as f32 - start_client_x) * scale;
                let dy = (ev.client_y() as f32 - start_client_y) * scale;
                positions.update(|p| {
                    if id < p.len() {
                        p[id] = (start_x + dx, start_y + dy);
                    }
                });
            }
        }
    };

    let on_up = move |_ev: web_sys::PointerEvent| {
        drag.set_value(DragState::None);
    };

    let on_wheel = move |ev: web_sys::WheelEvent| {
        ev.prevent_default();
        let Some((rx, ry, sw, sh)) = svg_rect() else {
            return;
        };
        if sw < 1.0 || sh < 1.0 {
            return;
        }
        let factor = if ev.delta_y() > 0.0 { 1.1 } else { 1.0 / 1.1 };
        let cx = ev.client_x() as f32 - rx;
        let cy = ev.client_y() as f32 - ry;
        viewbox.update(|vb| {
            let new_w = vb.2 * factor;
            let new_h = vb.3 * factor;
            let svg_cx = vb.0 + (cx / sw) * vb.2;
            let svg_cy = vb.1 + (cy / sh) * vb.3;
            vb.0 = svg_cx - (cx / sw) * new_w;
            vb.1 = svg_cy - (cy / sh) * new_h;
            vb.2 = new_w;
            vb.3 = new_h;
        });
    };

    let reset_view = move |_| {
        let init = initial_positions.with_value(|p| p.clone());
        let new_bb = layout::bounding_box(&init);
        viewbox.set(new_bb);
        positions.set(init);
    };

    let viewbox_attr = move || {
        let (x, y, w, h) = viewbox.get();
        format!("{} {} {} {}", x, y, w, h)
    };

    let edge_lines = move || {
        let p = positions.get();
        edges.with_value(|edges_vec| {
            edges_vec
                .iter()
                .map(|e| {
                    let from = e.from as usize;
                    let to = e.to as usize;
                    let (x1, y1) = p.get(from).copied().unwrap_or((0.0, 0.0));
                    let (x2, y2) = p.get(to).copied().unwrap_or((0.0, 0.0));
                    let color = e.kind.color();
                    let sw = (e.overlap as f32 / 8.0).clamp(0.6, 3.5);
                    let title_text = format!(
                        "{} → {}  ({}, overlap {} bp)",
                        from,
                        to,
                        e.kind.label(),
                        e.overlap
                    );
                    if from == to {
                        let r = 8.0;
                        let path = format!(
                            "M {} {} a {r} {r} 0 1 1 {} {}",
                            x1 + r,
                            y1,
                            -2.0 * r,
                            0.001
                        );
                        view! {
                            <path
                                d=path
                                fill="none"
                                stroke=color
                                stroke-width=sw
                                opacity="0.7"
                            >
                                <title>{title_text}</title>
                            </path>
                        }
                        .into_any()
                    } else {
                        view! {
                            <line
                                x1=x1
                                y1=y1
                                x2=x2
                                y2=y2
                                stroke=color
                                stroke-width=sw
                                opacity="0.55"
                            >
                                <title>{title_text}</title>
                            </line>
                        }
                        .into_any()
                    }
                })
                .collect::<Vec<_>>()
        })
    };

    let node_circles = move || {
        let p = positions.get();
        components.with_value(|comps| {
            seq_lengths.with_value(|lens| {
                (0..n)
                    .map(|i| {
                        let (x, y) = p.get(i).copied().unwrap_or((0.0, 0.0));
                        let color = component_color(comps[i]);
                        let length = lens[i];
                        let title_text = format!("oligo {} ({} bp)", i, length);
                        let on_node_down = move |ev: web_sys::PointerEvent| {
                            ev.stop_propagation();
                            let pos_i = positions.get_untracked().get(i).copied().unwrap_or((0.0, 0.0));
                            drag.set_value(DragState::Node {
                                id: i,
                                start_x: pos_i.0,
                                start_y: pos_i.1,
                                start_client_x: ev.client_x() as f32,
                                start_client_y: ev.client_y() as f32,
                            });
                        };
                        view! {
                            <circle
                                cx=x
                                cy=y
                                r="5"
                                fill=color
                                stroke="#0a0c10"
                                stroke-width="1.2"
                                style="cursor:grab;"
                                on:pointerdown=on_node_down
                            >
                                <title>{title_text}</title>
                            </circle>
                        }
                    })
                    .collect::<Vec<_>>()
            })
        })
    };

    view! {
        <section class="panel">
            <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:0.5rem;">
                <h2 style="margin:0; font-size:1.1rem;">"Overlap graph"</h2>
                <button on:click=reset_view style="font-size:0.85rem; padding:0.3rem 0.7rem;">"Reset layout"</button>
            </div>
            <div style="font-size:0.8rem; color:var(--muted); margin-bottom:0.5rem;">
                "Drag nodes · drag background to pan · scroll to zoom · "
                {format!("{} nodes, {} edges, {} component(s)", n, edges.with_value(|e| e.len()), n_components)}
            </div>
            <svg
                node_ref=svg_ref
                viewBox=viewbox_attr
                preserveAspectRatio="xMidYMid meet"
                on:pointerdown=on_bg_down
                on:pointermove=on_move
                on:pointerup=on_up
                on:pointerleave=on_up
                on:wheel=on_wheel
                style="display:block; width:100%; height:600px; background:#0a0c10; border-radius:6px; cursor:grab; touch-action:none;"
            >
                <g>{edge_lines}</g>
                <g>{node_circles}</g>
            </svg>
            <Legend />
        </section>
    }
}

#[component]
fn Legend() -> impl IntoView {
    view! {
        <div style="display:flex; flex-wrap:wrap; gap:1rem; margin-top:0.6rem; font-size:0.8rem;">
            <span style="color:var(--muted);">"Edges:"</span>
            <LegendSwatch color=EdgeKind::FwdFwd.color() label=EdgeKind::FwdFwd.label() />
            <LegendSwatch color=EdgeKind::FwdRev.color() label=EdgeKind::FwdRev.label() />
            <LegendSwatch color=EdgeKind::RevFwd.color() label=EdgeKind::RevFwd.label() />
            <LegendSwatch color=EdgeKind::RevRev.color() label=EdgeKind::RevRev.label() />
            <span style="color:var(--muted); margin-left:1rem;">"Nodes colored by connected component"</span>
        </div>
    }
}

#[component]
fn LegendSwatch(color: &'static str, label: &'static str) -> impl IntoView {
    view! {
        <span style="display:inline-flex; align-items:center; gap:0.35rem;">
            <span style=move || format!("display:inline-block; width:14px; height:3px; background:{};", color)></span>
            <span>{label}</span>
        </span>
    }
}
