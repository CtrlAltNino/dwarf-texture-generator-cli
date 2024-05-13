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
    persistence: f32,

    #[arg(short, long)]
    scale: f32,
}

impl Perlin {
    fn perlin_normed(&self, x: usize, y: usize) -> Vec<Vec<f64>> {
        let mut array: Vec<Vec<f64>> = vec![vec![0.0; y]; x];
        let mut ng = chatpgt::Noisegenerator::new();
        for o in 0..self.octaves {
            let ofactor = 2_i64.pow(o);
            let freq = self.scale * (ofactor as f32);
            let amplitude = 1.0 / ofactor as f64;
            for (px, row) in array.iter_mut().enumerate() {
                for (py, element) in row.iter_mut().enumerate() {
                    let a = amplitude
                        * ng.noise_2d((px as f32 * freq) as f64, (py as f32 * freq) as f64);
                    *element += a;
                }
            }
        }
        array
    }
}

impl Noise for Perlin {
    fn generate(&self, x: u32, y: u32) -> image::ImageBuffer<image::Rgb<u8>, Vec<u8>> {
        let mut imgbuf = image::ImageBuffer::new(x, y);
        let generated_noise = self.perlin_normed(x as usize, y as usize);

        for (px, py, pixel) in imgbuf.enumerate_pixels_mut() {
            let a = ((generated_noise[px as usize][py as usize] + 1.0) * 0.5 * 255.0) as u8;
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
