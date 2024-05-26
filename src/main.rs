use std::time::Instant;

use algorithms::voronoi::Voronoi;
use clap::builder::ArgPredicate;
use clap::ArgAction;
use clap::Args;
use clap::Parser;
use clap::Subcommand;
use clap_num::number_range;
use image::buffer::ConvertBuffer;
use image::ImageError;
use image::ImageFormat;

mod composition;
mod luatryout;
mod output;
mod algorithms {
    pub mod constant;
    pub mod perlin;
    pub mod voronoi;

    pub trait Noise {
        fn generate(&self, x: u32, y: u32) -> image::ImageBuffer<image::Rgb<f32>, Vec<f32>>;
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
#[command(disable_help_flag = true)]
struct Cli {
    // #[arg(short, long, value_parser=less_than_3, default_value_t=2)]
    // dimension: u8,
    // #[arg(short, long)]
    // path: std::path::PathBuf,
    #[arg(short, long, default_value_t = 1024)]
    width: u32,

    #[arg(short, long, default_value_t = 1024)]
    height: u32,

    #[arg(short = '?', long = "help", action = ArgAction::Help)]
    holp: Option<bool>,

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
    /// Path to the outputfile, image format is inffered from the filename (*.png, *.tff, ...)
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
    Voronoi(Voronoi),
}

fn post_process_cli_args(args: &mut crate::Cli) -> Result<(), ImageError> {
    output::infer_file_format(args)?;
    Ok(())
}

fn main() {
    let start = Instant::now();
    let mut args = Cli::parse();
    post_process_cli_args(&mut args).unwrap();

    // I dont like this repitition
    let mut image = match args.noise_type {
        NoiseType::Perlin(nargs) => nargs.generate(args.width, args.height),
        NoiseType::Constant(nargs) => nargs.generate(args.width, args.height),
        NoiseType::Voronoi(nargs) => nargs.generate(args.width, args.height),
    };

    // let t = composition::color_gradient3();

    let tttt = vec![
        image::Rgb([0.0, 0.0, 0.0]),
        image::Rgb([0.0, 0.0, 0.0]),
        image::Rgb([1.0, 0.0, 0.0]),
        image::Rgb([1.0, 1.0, 0.0]),
        image::Rgb([1.0, 1.0, 1.0]),
        image::Rgb([0.647, 0.0, 0.969]),
    ];

    let t = composition::build_gradient(&tttt);

    for (px, py, pixel) in image.enumerate_pixels_mut() {
        let image::Rgb([d, _, _]) = *pixel;
        *pixel = t(d);
    }

    let u8image: image::RgbImage = image.convert();

    if args.out.stdout {
        output::spit_to_stdout(u8image, args.image_format.unwrap());
    } else {
        output::save_to_file(
            u8image,
            args.image_format.unwrap(),
            args.out.path.unwrap().as_path(),
        );
    }
    let duration = start.elapsed();
    eprintln!("Time elapsed in expensive_function() is: {:?}", duration);
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
