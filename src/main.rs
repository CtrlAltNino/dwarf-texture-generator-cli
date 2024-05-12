use clap::Args;
use clap::Parser;
use clap::Subcommand;
use clap::ValueEnum;
use clap_num::number_range;
use image::ImageFormat;
use std::io::Cursor;
use std::io::{self, Write};

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
    #[arg(short, long, value_parser=less_than_3, default_value_t=2)]
    dimension: u8,
    /// Path to the outputfile, image format is inffered from the filename (*.png, *.tff, ...)
    // #[arg(short, long)]
    // path: std::path::PathBuf,

    /// help texto
    #[command(subcommand)]
    noise_type: NoiseType,

    #[command(flatten)]
    out: Output,

    image_format: Option<Imageformat>,
}

#[derive(Args, Debug)]
#[group(required = true, multiple = false)]
struct Output {
    #[arg(short, long)]
    path: Option<std::path::PathBuf>,

    #[arg(short, long)]
    stdout: bool,
}

#[derive(Args, Debug)]
#[group(multiple = true, required = true)]
struct StdoutRequirements {
    #[arg(short, long)]
    stdout: bool,

    #[arg(short, long)]
    imageformat: Imageformat,
}

#[derive(ValueEnum, Debug, Clone)]
enum Imageformat {
    png,
    tff,
}

#[derive(Subcommand, Debug)]
#[command(subcommand_value_name = "NOISETYPE")]
#[command(subcommand_help_heading("Noise types"))]
#[command(disable_help_subcommand = true)]
enum NoiseType {
    Perlin(Perlin),
    Voronoi,
}

fn main() {
    let args = Cli::parse();

    let x = match args.noise_type {
        NoiseType::Perlin(args) => Some(args.generate(1000, 1000)),
        _ => None,
    };

    let mut cursor = Cursor::new(Vec::new());

    x.unwrap().write_to(&mut cursor, ImageFormat::Png).unwrap();

    io::stdout().write_all(&cursor.into_inner()).unwrap();
    io::stdout().flush().unwrap();
}
