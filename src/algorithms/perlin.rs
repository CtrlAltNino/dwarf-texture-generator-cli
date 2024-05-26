use clap::Args;

use crate::algorithms::Noise;
use crate::chatpgt;

#[derive(Debug, Args)]
pub struct Perlin {
    // #[arg(short, long)]
    // seed: u32,
    #[arg(short, long)]
    octaves: u32,

    #[arg(short, long)]
    persistence: f64,

    #[arg(short, long)]
    scale: f64,
}

impl Perlin {
    fn perlin_normed(&self, x: usize, y: usize) -> Vec<Vec<f64>> {
        let mut array: Vec<Vec<f64>> = vec![vec![0.0; y]; x];
        let mut ng = chatpgt::Noisegenerator::new();
        let max_val = (0..self.octaves).fold(0.0, |acc, x| acc + self.persistence.powf(x as f64));

        // let amplitude = self.persistence.powf(o as f64);
        for (px, row) in array.iter_mut().enumerate() {
            for (py, element) in row.iter_mut().enumerate() {
                let mut ofactor = 1.0;
                let mut amplitude = 1.0;
                for _ in 0..self.octaves {
                    // let ofactor = 2_i64.pow(o);
                    let freq = self.scale * (ofactor as f64);

                    *element += amplitude * ng.noise_2d(px as f64 * freq, py as f64 * freq);

                    //max_val += amplitude;
                    amplitude *= self.persistence;
                    ofactor *= 2.0;
                }
                *element /= max_val;
            }
        }
        array
    }
}

fn comp<T, U, V>(g: impl Fn(U) -> V, h: impl Fn(T) -> U) -> impl Fn(T) -> V {
    move |x| g(h(x))
}

impl Noise for Perlin {
    fn generate(&self, x: u32, y: u32) -> image::ImageBuffer<image::Rgb<f32>, Vec<f32>> {
        let mut imgbuf = image::ImageBuffer::new(x, y);
        let generated_noise = self.perlin_normed(x as usize, y as usize);

        for (px, py, pixel) in imgbuf.enumerate_pixels_mut() {
            let a = ((generated_noise[px as usize][py as usize] + 1.0) * 0.5) as f32;
            *pixel = image::Rgb([a, a, a]);
        }
        imgbuf
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic_case() {
        let p = Perlin {
            octaves: 13,
            persistence: 1.0,
            scale: 0.001,
        };

        let r = p.generate(10, 10);

        assert_eq!((10, 10), r.dimensions());
    }
}
