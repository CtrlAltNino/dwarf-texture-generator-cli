use clap::Parser;
use clap_num::number_range;

mod algorithms {
    pub mod perlin;
}

use crate::algorithms::perlin;

fn less_than_3(s: &str) -> Result<u8, String> {
    number_range(s, 0, 3)
}

#[derive(Parser, Debug)]
struct Cli {
    /// The pattern to look for
    #[arg(short, long)]
    noise_type: String,

    #[arg(short, long, value_parser=less_than_3, default_value_t=2)]
    dimension: u8,
    // The path to the outputfile
    #[arg(short, long)]
    path: std::path::PathBuf,
}

fn main() {
    let args = Cli::parse();

    println!("Hello, world! {:?}", args);
    let x = perlin::PerlinParameters {
        seed: 1,
        scale: 1.0,
        octaves: 1,
        persistence: 1.0,
    };

    let out = perlin::perlin_2d(800, 800, x);

    out.save(args.path.as_path()).unwrap();
}
