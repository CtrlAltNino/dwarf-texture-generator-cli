#[derive(Debug)]
pub struct PerlinParameters {
    pub seed: u32,
    pub octaves: u32,
    pub persistence: f32,
    pub scale: f32,
}

pub fn perlin_2d(
    x: u32,
    y: u32,
    params: PerlinParameters,
) -> image::ImageBuffer<image::Rgb<u8>, Vec<u8>> {
    let mut imgbuf = image::ImageBuffer::new(x, y);

    for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
        let r = (0.8 * x as f32) as u8;
        let b = (0.3 * y as f32) as u8;
        *pixel = image::Rgb([r, 0, b]);
    }

    imgbuf
}
