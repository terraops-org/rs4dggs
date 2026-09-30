//! The `rs4dggs` command: DGGAL's discrete global grids from the shell.
#![forbid(unsafe_code)]

use std::io::{self, BufWriter};
use std::process::ExitCode;

fn main() -> ExitCode {
    // Lossily: an argument that is not UTF-8 is answered as a malformed one, never a panic,
    // which is what `std::env::args` would do.
    let args: Vec<String> = std::env::args_os()
        .skip(1)
        .map(|a| a.to_string_lossy().into_owned())
        .collect();
    let mut input = io::stdin().lock();
    let mut out = BufWriter::new(io::stdout().lock());
    let mut err = io::stderr().lock();
    ExitCode::from(rs4dggs_cli::run(&args, &mut input, &mut out, &mut err))
}
