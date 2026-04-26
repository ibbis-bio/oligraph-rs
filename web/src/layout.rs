use std::f32::consts::PI;

pub fn path_seeded_positions(n_nodes: usize, contig_paths: &[Vec<u32>]) -> Vec<(f32, f32)> {
    const X_STEP: f32 = 60.0;
    const Y_STEP: f32 = 120.0;

    let mut pos = vec![(0.0_f32, 0.0_f32); n_nodes];
    let mut placed = vec![false; n_nodes];

    for (row, path) in contig_paths.iter().enumerate() {
        let y = row as f32 * Y_STEP;
        for (col, &id) in path.iter().enumerate() {
            let i = id as usize;
            if i < n_nodes {
                pos[i] = (col as f32 * X_STEP, y);
                placed[i] = true;
            }
        }
    }

    let n_contigs = contig_paths.len();
    let base_y = n_contigs as f32 * Y_STEP + Y_STEP;
    let mut iso_col = 0usize;
    let mut iso_row = 0usize;
    for i in 0..n_nodes {
        if !placed[i] {
            pos[i] = (iso_col as f32 * X_STEP, base_y + iso_row as f32 * Y_STEP);
            iso_col += 1;
            if iso_col >= 10 {
                iso_col = 0;
                iso_row += 1;
            }
        }
    }

    pos
}

pub fn fruchterman_reingold(
    n_nodes: usize,
    edges: &[(u32, u32)],
    iterations: usize,
    initial_positions: Option<Vec<(f32, f32)>>,
) -> Vec<(f32, f32)> {
    if n_nodes == 0 {
        return Vec::new();
    }
    let area = 500.0_f32 * 500.0_f32;
    let k = (area / n_nodes as f32).sqrt();

    let mut pos: Vec<(f32, f32)> = match initial_positions {
        Some(p) if p.len() == n_nodes => p,
        _ => (0..n_nodes)
            .map(|i| {
                let theta = i as f32 * 137.5 * PI / 180.0;
                let r = (i as f32 + 1.0).sqrt() * 25.0;
                (r * theta.cos(), r * theta.sin())
            })
            .collect(),
    };
    let mut disp: Vec<(f32, f32)> = vec![(0.0, 0.0); n_nodes];

    let mut t: f32 = (n_nodes as f32).sqrt() * 10.0;
    let cooling = (0.001_f32).powf(1.0 / iterations.max(1) as f32);

    for _ in 0..iterations {
        for d in disp.iter_mut() {
            *d = (0.0, 0.0);
        }

        for i in 0..n_nodes {
            for j in (i + 1)..n_nodes {
                let dx = pos[i].0 - pos[j].0;
                let dy = pos[i].1 - pos[j].1;
                let d2 = dx * dx + dy * dy;
                let d = d2.sqrt().max(0.01);
                let force = k * k / d;
                let fx = dx / d * force;
                let fy = dy / d * force;
                disp[i].0 += fx;
                disp[i].1 += fy;
                disp[j].0 -= fx;
                disp[j].1 -= fy;
            }
        }

        for &(u, v) in edges {
            let u = u as usize;
            let v = v as usize;
            if u == v || u >= n_nodes || v >= n_nodes {
                continue;
            }
            let dx = pos[u].0 - pos[v].0;
            let dy = pos[u].1 - pos[v].1;
            let d = (dx * dx + dy * dy).sqrt().max(0.01);
            let force = d * d / k;
            let fx = dx / d * force;
            let fy = dy / d * force;
            disp[u].0 -= fx;
            disp[u].1 -= fy;
            disp[v].0 += fx;
            disp[v].1 += fy;
        }

        for i in 0..n_nodes {
            let (dx, dy) = disp[i];
            let dlen = (dx * dx + dy * dy).sqrt().max(0.01);
            let limited = dlen.min(t);
            pos[i].0 += dx / dlen * limited;
            pos[i].1 += dy / dlen * limited;
        }

        t *= cooling;
    }

    pos
}

pub fn bounding_box(positions: &[(f32, f32)]) -> (f32, f32, f32, f32) {
    if positions.is_empty() {
        return (-50.0, -50.0, 100.0, 100.0);
    }
    let mut min_x = f32::INFINITY;
    let mut min_y = f32::INFINITY;
    let mut max_x = f32::NEG_INFINITY;
    let mut max_y = f32::NEG_INFINITY;
    for &(x, y) in positions {
        if x < min_x {
            min_x = x;
        }
        if y < min_y {
            min_y = y;
        }
        if x > max_x {
            max_x = x;
        }
        if y > max_y {
            max_y = y;
        }
    }
    let w = (max_x - min_x).max(50.0);
    let h = (max_y - min_y).max(50.0);
    let pad_x = w * 0.1;
    let pad_y = h * 0.1;
    (
        min_x - pad_x,
        min_y - pad_y,
        w + 2.0 * pad_x,
        h + 2.0 * pad_y,
    )
}
