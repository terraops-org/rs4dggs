//! A grid chosen by name at run time, as a web service receives it, and the disk around a
//! zone served one ring at a time.
//!
//! ```text
//! cargo run --example service -- <grid> <zone id> <radius>
//! cargo run --example service -- rtea7h 006415650342 2
//! ```

use std::process::ExitCode;

use rs4dggs::{GRID_NAMES, get_grid};

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let [name, text, k] = args.as_slice() else {
        eprintln!(
            "usage: service <grid> <zone id> <radius>   (grids: {})",
            GRID_NAMES.join(", ")
        );
        return ExitCode::FAILURE;
    };
    match run(name, text, k) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run(name: &str, text: &str, k: &str) -> Result<(), Box<dyn std::error::Error>> {
    let grid = get_grid(name)?; // any casing; an unknown name lists the valid ones
    let centre = grid.zone_from_text(text)?;
    let k: u32 = k.parse().map_err(|e| format!("radius {k:?}: {e}"))?;
    // Each ring is one page of the answer. The walk computes a ring only when it is asked
    // for, so a service can stop after any page, however large the radius requested.
    for (n, ring) in grid.disk(centre, k).rings().enumerate() {
        let ids: Vec<String> = ring.iter().map(|&id| grid.text_id(id)).collect();
        println!("ring {n}: {}", ids.join(" "));
    }
    Ok(())
}
