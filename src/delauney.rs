use delaunator::*;

fn delauney_triangulation() {
    let points: Vec<delaunator::Point> = (0..100)
        .map(|_| delaunator::Point {
            x: rng.gen_range(0.0..x as f64),
            y: rng.gen_range(0.0..y as f64),
        })
        .collect();
}
