//! Quantise points on all six grids at one resolution, from 0 to 19, the range they share.
//!
//! ```text
//! cargo run --example quantise -- <resolution> <lat> <lon> [<lat> <lon> ...]
//! cargo run --example quantise -- 10 38.7223 -9.1393 51.5074 -0.1278
//! ```

use std::process::ExitCode;

use rs4dggs::{GRID_NAMES, get_grid};

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() < 3 || args.len() % 2 == 0 {
        eprintln!("usage: quantise <resolution> <lat> <lon> [<lat> <lon> ...]");
        return ExitCode::FAILURE;
    }
    match run(&args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run(args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    // The whole table is computed before anything is printed, so that a bad argument, one that
    // does not parse or one a grid refuses, is reported alone rather than after a header and
    // part of the table.
    let res: u8 = parse("resolution", &args[0])?;
    let mut rows = Vec::new();
    for (i, pair) in args[1..].chunks(2).enumerate() {
        let lat: f64 = parse(&format!("latitude {}", i + 1), &pair[0])?;
        let lon: f64 = parse(&format!("longitude {}", i + 1), &pair[1])?;
        let mut row = format!("{lat}\t{lon}");
        for name in GRID_NAMES {
            let grid = get_grid(name)?;
            let id = grid.zone_from_geo(lat, lon, res)?;
            row.push('\t');
            row.push_str(&grid.text_id(id));
        }
        rows.push(row);
    }
    println!("lat\tlon\t{}", GRID_NAMES.join("\t"));
    for row in rows {
        println!("{row}");
    }
    Ok(())
}

/// Reads one argument, naming it in the error if it does not parse.
fn parse<T>(what: &str, text: &str) -> Result<T, String>
where
    T: std::str::FromStr,
    T::Err: std::fmt::Display,
{
    text.parse().map_err(|e| format!("{what} {text:?}: {e}"))
}
