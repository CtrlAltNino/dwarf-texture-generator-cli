use std::usize;

use rand::Rng;

#[derive(Debug, Copy, Clone)]
struct Vector2 {
    x: f64,
    y: f64,
}

impl Vector2 {
    fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    fn dot(&self, other: Vector2) -> f64 {
        self.x * other.x + self.y * other.y
    }
}

fn shuffle(array_to_shuffle: &mut [usize]) {
    let mut rng = rand::thread_rng();
    for e in (1..array_to_shuffle.len()).rev() {
        let index = rng.gen_range(0..e);
        array_to_shuffle.swap(e, index);
    }
}

fn get_constant_vector(v: usize) -> Vector2 {
    let h = v & 3;
    match h {
        0 => Vector2::new(1.0, 1.0),
        1 => Vector2::new(-1.0, 1.0),
        2 => Vector2::new(-1.0, -1.0),
        _ => Vector2::new(1.0, -1.0),
    }
}

fn fade(t: f64) -> f64 {
    ((6.0 * t - 15.0) * t + 10.0) * t * t * t
}

fn lerp(t: f64, a1: f64, a2: f64) -> f64 {
    a1 + t * (a2 - a1)
}

fn make_permutation() -> Vec<usize> {
    let mut permutation = Vec::with_capacity(256);
    for i in 0..256 {
        permutation.push(i);
    }

    shuffle(&mut permutation);

    for i in 0..256 {
        permutation.push(permutation[i]);
    }

    permutation
}

pub struct Noisegenerator {
    permutation: Vec<usize>,
}

impl Noisegenerator {
    pub fn new() -> Self {
        Self {
            permutation: make_permutation(),
        }
    }

    pub fn noise_2d(&self, x: f64, y: f64) -> f64 {
        let x_floor = x.floor() as usize & 255;
        let y_floor = y.floor() as usize & 255;

        let xf = x - x.floor();
        let yf = y - y.floor();

        let top_right = Vector2::new(xf - 1.0, yf - 1.0);
        let top_left = Vector2::new(xf, yf - 1.0);
        let bottom_right = Vector2::new(xf - 1.0, yf);
        let bottom_left = Vector2::new(xf, yf);

        let value_top_right = self.permutation[self.permutation[x_floor + 1] + y_floor + 1];
        let value_top_left = self.permutation[self.permutation[x_floor] + y_floor + 1];
        let value_bottom_right = self.permutation[self.permutation[x_floor + 1] + y_floor];
        let value_bottom_left = self.permutation[self.permutation[x_floor] + y_floor];

        let dot_top_right = top_right.dot(get_constant_vector(value_top_right));
        let dot_top_left = top_left.dot(get_constant_vector(value_top_left));
        let dot_bottom_right = bottom_right.dot(get_constant_vector(value_bottom_right));
        let dot_bottom_left = bottom_left.dot(get_constant_vector(value_bottom_left));

        let u = fade(xf);
        let v = fade(yf);

        lerp(
            u,
            lerp(v, dot_bottom_left, dot_top_left),
            lerp(v, dot_bottom_right, dot_top_right),
        )
    }
}
