use std::ops::Mul;

use crate::output;

use image::{Pixel, Rgb};

fn comp<T, U, V>(g: impl Fn(U) -> V, h: impl Fn(T) -> U) -> impl Fn(T) -> V {
    move |x| g(h(x))
}

pub trait Interpolation {
    fn interpolate(&self, y: &Self, a: f32) -> Self;
}

impl Interpolation for f32 {
    fn interpolate(&self, y: &Self, a: f32) -> f32 {
        a * y + (1.0 - a) * self
    }
}

impl<T> Interpolation for Vec<T>
where
    T: Interpolation,
{
    fn interpolate(&self, y: &Vec<T>, a: f32) -> Vec<T> {
        self.iter()
            .zip(y.iter())
            .map(|(u, v)| u.interpolate(v, a))
            .collect()
    }
}

impl Interpolation for image::Rgb<f32> {
    fn interpolate(&self, y: &Self, a: f32) -> Self {
        let image::Rgb([r1, g1, b1]) = self;
        let image::Rgb([r2, g2, b2]) = y;
        image::Rgb([
            r1.interpolate(r2, a),
            g1.interpolate(g2, a),
            b1.interpolate(b2, a),
        ])
    }
}

pub fn concat_paths<U>(v: Vec<impl Fn(f32) -> U>) -> impl Fn(f32) -> U {
    let l = v.len() as f32;
    move |x| {
        let i = x * l;
        let a = i.floor();
        v[a.min(l - 1.0) as usize](i - a)
    }
}

pub fn build_gradient<'a>(
    colors: &'a Vec<image::Rgb<f32>>,
) -> impl Fn(f32) -> image::Rgb<f32> + 'a {
    let fs: Vec<_> = colors
        .windows(2)
        .map(|window| |x| window[0].interpolate(&window[1], x))
        .collect();
    concat_paths(fs)
}
