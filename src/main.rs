use clap::builder::ArgPredicate;
use clap::Args;
use clap::CommandFactory;
use clap::Parser;
use clap::Subcommand;
use clap::ValueEnum;
use clap_num::number_range;
use image::ImageFormat;
use std::fmt::format;
use std::io::Cursor;
use std::io::ErrorKind;
use std::io::{self, Write};
use std::path::PathBuf;

mod luatryout;
mod output;
mod algorithms {
    pub mod constant;
    pub mod perlin;

    pub trait Noise {
        fn generate(&self, x: u32, y: u32) -> image::ImageBuffer<image::Rgb<u8>, Vec<u8>>;
    }
}

mod chatpgt;

use crate::algorithms::constant::Constant;
use crate::algorithms::perlin::Perlin;
use crate::algorithms::*;

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

    /// png, tif, ...
    #[arg(short, long, value_parser = crate::output::parse_file_format, default_value = "png", default_value_if("path" , ArgPredicate::IsPresent,  None))]
    image_format: Option<ImageFormat>,

    #[command(flatten)]
    out: Output,
}

#[derive(Args, Debug)]
#[group(required = true, multiple = false)]
struct Output {
    #[arg(short, long)]
    path: Option<std::path::PathBuf>,

    #[arg(short, long)]
    stdout: bool,
}

#[derive(Subcommand, Debug)]
#[command(subcommand_value_name = "NOISETYPE")]
#[command(subcommand_help_heading("Noise types"))]
#[command(disable_help_subcommand = true)]
enum NoiseType {
    Perlin(Perlin),
    Constant(Constant),
}

fn post_process_cli_args(args: &mut crate::Cli) -> () {
    output::infer_file_format(args);
}

fn main() {
    let mut args = Cli::parse();
    post_process_cli_args(&mut args);

    // I dont like this repitition
    let image = match args.noise_type {
        NoiseType::Perlin(args) => args.generate(1000, 1000),
        NoiseType::Constant(args) => args.generate(1000, 1000),
    };

    if args.out.stdout {
        output::spit_to_stdout(image, args.image_format.unwrap());
    } else {
        output::save_to_file(
            image,
            args.image_format.unwrap(),
            args.out.path.unwrap().as_path(),
        );
    }
}

#[cfg(test)]
mod tests {

    use std::error::Error;

    use super::*;

    trait Testable {
        fn run(cmd: &[&str]) -> Result<Self, clap::Error>
        where
            Self: Sized;
    }

    impl Testable for Cli {
        fn run(cmd: &[&str]) -> Result<Self, clap::Error>
        where
            Self: Sized,
        {
            let mut r = Self::try_parse_from(cmd)?;
            post_process_cli_args(&mut r);
            Ok(r)
        }
    }

    #[test]
    fn inferred_image_format() -> Result<(), clap::Error> {
        let r = Cli::run(&[
            "",
            "--path",
            "test.png",
            "perlin",
            "--octaves",
            "1",
            "--scale",
            "1",
            "--persistence",
            "1",
        ])?;
        assert_eq!(r.image_format, Some(ImageFormat::Png));

        Ok(())
    }
}
