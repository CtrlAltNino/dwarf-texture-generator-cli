use clap::Args;

use crate::algorithms::Noise;

#[derive(Debug, Args)]
pub struct Constant {
    // #[arg(short, long)]
    // seed: u32,
}

impl Noise for Constant {
    fn generate(&self, x: u32, y: u32) -> image::ImageBuffer<image::Rgb<f32>, Vec<f32>> {
        let mut imgbuf = image::ImageBuffer::new(x, y);

        for (_px, _py, pixel) in imgbuf.enumerate_pixels_mut() {
            *pixel = image::Rgb([0.1, 0.3, 0.8]);
        }
        imgbuf
    }
}

// #[cfg(test)]
// mod tests {
//     use super::*;
//
//     #[test]
//     fn basic_case() {
//         let p = Perlin {
//             octaves: 13,
//             persistence: 1.0,
//             scale: 0.001,
//         };
//
//         let r = p.generate(10, 10);
//
//         assert_eq!((10, 10), r.dimensions());
//     }
// }
