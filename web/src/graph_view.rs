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
    fn colour(self) -> &'static str {
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
    pub initial_viewbox: (f32, f32, f32, f32),
    pub connected: Vec<bool>,
}

const COMPONENT_COLOURS: &[&str] = &[
    "#60a5fa", "#4ade80", "#fbbf24", "#f87171", "#67e8f9", "#34d399", "#fb923c", "#c084fc",
    "#f472b6", "#818cf8", "#a3e635", "#fb7185",
];

fn component_colour(comp: Option<usize>) -> &'static str {
    match comp {
        None => "#9ca3af",
        Some(i) => COMPONENT_COLOURS[i % COMPONENT_COLOURS.len()],
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
pub fn GraphView(graph: GraphData, highlight_comp: ReadSignal<Option<usize>>) -> impl IntoView {
    let n = graph.n_nodes;
    let n_components = graph.n_components;
    let positions = RwSignal::new(graph.initial_positions.clone());
    let viewbox = RwSignal::new(graph.initial_viewbox);
    let drag = StoredValue::new(DragState::None);
    let svg_ref: NodeRef<leptos::svg::Svg> = NodeRef::new();

    let initial_positions = StoredValue::new(graph.initial_positions.clone());
    let initial_viewbox = StoredValue::new(graph.initial_viewbox);
    let edges = StoredValue::new(graph.edges);
    let components = StoredValue::new(graph.components);
    let seq_lengths = StoredValue::new(graph.seq_lengths);
    let connected = StoredValue::new(graph.connected);
    let hide_isolated = RwSignal::new(true);

    let n_isolated =
        connected.with_value(|c| c.iter().filter(|&&v| !v).count());

    let svg_rect = move || -> Option<(f32, f32, f32, f32)> {
        svg_ref.get_untracked().map(|el| {
            let r = el.get_bounding_client_rect();
            (
                r.x() as f32,
                r.y() as f32,
                r.width() as f32,
                r.height() as f32,
            )
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

    // Restore seeded positions and zoom to the largest contig.
    let reset_view = move |_| {
        let init = initial_positions.with_value(|p| p.clone());
        positions.set(init);
        viewbox.set(initial_viewbox.get_value());
    };

    // Fit visible node positions in view without moving nodes.
    let view_all = move |_| {
        let p = positions.get_untracked();
        let hide = hide_isolated.get_untracked();
        let pts: Vec<(f32, f32)> = connected.with_value(|c| {
            p.iter()
                .enumerate()
                .filter(|&(i, _)| !hide || c.get(i).copied().unwrap_or(false))
                .map(|(_, &pos)| pos)
                .collect()
        });
        if !pts.is_empty() {
            viewbox.set(layout::bounding_box(&pts));
        }
    };

    let viewbox_attr = move || {
        let (x, y, w, h) = viewbox.get();
        format!("{} {} {} {}", x, y, w, h)
    };

    let edge_lines = move || {
        let p = positions.get();
        let hide = hide_isolated.get();
        let highlighted = highlight_comp.get();
        edges.with_value(|edges_vec| {
            connected.with_value(|conn| {
                components.with_value(|comps| {
                    edges_vec
                        .iter()
                        .filter_map(|e| {
                            let from = e.from as usize;
                            let to = e.to as usize;
                            if hide && !conn.get(from).copied().unwrap_or(false) {
                                return None;
                            }
                            let (x1, y1) = p.get(from).copied().unwrap_or((0.0, 0.0));
                            let (x2, y2) = p.get(to).copied().unwrap_or((0.0, 0.0));
                            let colour = e.kind.colour();
                            let sw = (e.overlap as f32 / 8.0).clamp(0.6, 3.5);
                            
                            let comp_from = comps[from];
                            let is_highlighted = highlighted.is_some() && highlighted == comp_from;
                            let opacity = if highlighted.is_none() || is_highlighted { "0.85" } else { "0.1" };
                            let stroke_width = if is_highlighted { sw * 1.5 } else { sw };

                            let title_text = format!(
                                "{} → {}  ({}, overlap {} bp)",
                                from,
                                to,
                                e.kind.label(),
                                e.overlap
                            );
                            Some(if from == to {
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
                                        stroke=colour
                                        stroke-width=stroke_width
                                        opacity=opacity
                                        style="transition: opacity 0.2s, stroke-width 0.2s;"
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
                                        stroke=colour
                                        stroke-width=stroke_width
                                        opacity=opacity
                                        style="transition: opacity 0.2s, stroke-width 0.2s;"
                                    >
                                        <title>{title_text}</title>
                                    </line>
                                }
                                .into_any()
                            })
                        })
                        .collect::<Vec<_>>()
                })
            })
        })
    };

    let node_circles = move || {
        let p = positions.get();
        let hide = hide_isolated.get();
        let highlighted = highlight_comp.get();
        connected.with_value(|conn| {
            components.with_value(|comps| {
            seq_lengths.with_value(|lens| {
                (0..n)
                    .filter_map(|i| {
                        if hide && !conn.get(i).copied().unwrap_or(false) {
                            return None;
                        }
                        let (x, y) = p.get(i).copied().unwrap_or((0.0, 0.0));
                        let colour = component_colour(comps[i]);
                        let length = lens[i];
                        
                        let comp = comps[i];
                        let is_highlighted = highlighted.is_some() && highlighted == comp;
                        let opacity = if highlighted.is_none() || is_highlighted { "1.0" } else { "0.15" };
                        let radius = if is_highlighted { "10" } else { "7" };

                        let title_text = format!("oligo {} ({} bp)", i, length);
                        let on_node_down = move |ev: web_sys::PointerEvent| {
                            ev.stop_propagation();
                            let pos_i = positions
                                .get_untracked()
                                .get(i)
                                .copied()
                                .unwrap_or((0.0, 0.0));
                            drag.set_value(DragState::Node {
                                id: i,
                                start_x: pos_i.0,
                                start_y: pos_i.1,
                                start_client_x: ev.client_x() as f32,
                                start_client_y: ev.client_y() as f32,
                            });
                        };
                        Some(view! {
                            <circle
                                cx=x
                                cy=y
                                r=radius
                                fill=colour
                                stroke="#1e2330"
                                stroke-width="1.2"
                                opacity=opacity
                                style="cursor:grab; transition: opacity 0.2s, r 0.2s;"
                                on:pointerdown=on_node_down
                            >
                                <title>{title_text}</title>
                            </circle>
                        })
                    })
                    .collect::<Vec<_>>()
            })
            })
        })
    };

    view! {
        <section class="panel" style="padding: 1.25rem;">
            <div style="display:flex; justify-content:space-between; align-items:flex-start; margin-bottom:1.25rem;">
                <div>
                    <h2 style="margin-bottom: 0.25rem; display: flex; align-items: center; gap: 0.5rem;">
                        <span style="width: 8px; height: 8px; border-radius: 2px; background: var(--warn);"></span>
                        "Overlap Graph"
                    </h2>
                    <div style="font-size:0.8rem; color:var(--muted); font-family: var(--font-mono);">
                        {format!("{} nodes • {} edges • {} components", n, edges.with_value(|e| e.len()), n_components)}
                    </div>
                </div>
                <div style="display:flex; gap:0.6rem;">
                    <button
                        class="btn-secondary"
                        on:click=move |_| hide_isolated.update(|v| *v = !*v)
                    >
                        {move || if hide_isolated.get() {
                            format!("Show isolated ({})", n_isolated)
                        } else {
                            format!("Hide isolated ({})", n_isolated)
                        }}
                    </button>
                    <button class="btn-secondary" on:click=view_all>"View all"</button>
                    <button class="btn-secondary" on:click=reset_view>"Reset view"</button>
                </div>
            </div>

            <div class="graph-container">
                <svg
                    node_ref=svg_ref
                    viewBox=viewbox_attr
                    preserveAspectRatio="xMidYMid meet"
                    on:pointerdown=on_bg_down
                    on:pointermove=on_move
                    on:pointerup=on_up
                    on:pointerleave=on_up
                    on:wheel=on_wheel
                    style="display:block; width:100%; height:600px; cursor:grab; touch-action:none;"
                >
                    <defs>
                        <pattern id="grid" width="100" height="100" patternUnits="userSpaceOnUse">
                            <path d="M 100 0 L 0 0 0 100" fill="none" stroke="rgba(255,255,255,0.03)" stroke-width="1"/>
                        </pattern>
                    </defs>
                    <rect width="100000" height="100000" x="-50000" y="-50000" fill="url(#grid)" />
                    
                    <g>{edge_lines}</g>
                    <g>{node_circles}</g>
                </svg>
                
                <div style="position: absolute; bottom: 1rem; left: 1rem; pointer-events: none; background: rgba(0,0,0,0.4); backdrop-filter: blur(4px); padding: 0.5rem 0.75rem; border-radius: 6px; border: 1px solid rgba(255,255,255,0.05); font-size: 0.75rem; color: var(--muted);">
                    "DRAG TO PAN • SCROLL TO ZOOM • DRAG NODES"
                </div>
            </div>

            <Legend />
        </section>
    }
}

#[component]
fn Legend() -> impl IntoView {
    view! {
        <div style="display:flex; flex-wrap:wrap; align-items: center; gap:1.25rem; margin-top:1.25rem; padding: 0.75rem 1rem; background: rgba(0,0,0,0.15); border-radius: 6px; font-size: 0.75rem;">
            <div style="display: flex; align-items: center; gap: 0.75rem;">
                <span style="color:var(--muted); text-transform: uppercase; letter-spacing: 0.05em; font-weight: 600;">"Edges"</span>
                <LegendSwatch colour=EdgeKind::FwdFwd.colour() label=EdgeKind::FwdFwd.label() />
                <LegendSwatch colour=EdgeKind::FwdRev.colour() label=EdgeKind::FwdRev.label() />
                <LegendSwatch colour=EdgeKind::RevFwd.colour() label=EdgeKind::RevFwd.label() />
                <LegendSwatch colour=EdgeKind::RevRev.colour() label=EdgeKind::RevRev.label() />
            </div>
            <div style="height: 12px; width: 1px; background: var(--border);"></div>
            <div style="display: flex; align-items: center; gap: 0.5rem;">
                <span style="color:var(--muted); text-transform: uppercase; letter-spacing: 0.05em; font-weight: 600;">"Nodes"</span>
                <span style="color: var(--text);">"Coloured by connected component"</span>
            </div>
        </div>
    }
}

#[component]
fn LegendSwatch(colour: &'static str, label: &'static str) -> impl IntoView {
    view! {
        <span style="display:inline-flex; align-items:center; gap:0.4rem; background: rgba(255,255,255,0.03); padding: 0.2rem 0.5rem; border-radius: 4px; border: 1px solid rgba(255,255,255,0.03);">
            <span style=move || format!("display:inline-block; width:10px; height:10px; border-radius: 2px; background:{};", colour)></span>
            <span style="color: var(--text); font-family: var(--font-mono); font-size: 0.7rem;">{label}</span>
        </span>
    }
}



