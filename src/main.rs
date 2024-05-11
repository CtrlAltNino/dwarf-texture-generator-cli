use clap::Parser;
use clap_num::number_range;

mod luatryout;
mod algorithms {
    pub mod perlin;
}

mod chatpgt;

use crate::algorithms::perlin::Noise;
use crate::algorithms::perlin::Perlin;

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

    // luatryout::main().unwrap();

    println!("Hello, world! {:?}", args);
    let x = Perlin {
        seed: 1,
        scale: 0.001,
        octaves: 13,
        persistence: 1.0,
    };

    let out = x.generate(1000, 1000);

    out.save(args.path.as_path()).unwrap();
}
