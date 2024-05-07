use clap::Parser;
use clap_num::number_range;

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
    // path: std::path::PathBuf,
}

fn main() {
    let args = Cli::parse();
    println!("Hello, world! {:?}", args);
}
