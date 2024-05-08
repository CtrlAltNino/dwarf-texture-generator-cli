#[derive(Debug)]
pub struct PerlinParameters {
    pub seed: u32,
    pub octaves: u32,
    pub persistence: f32,
    pub scale: f32,
}

pub fn perlin_2d(x: usize, y: usize, params: PerlinParameters) -> f32 {
    1.0
}
