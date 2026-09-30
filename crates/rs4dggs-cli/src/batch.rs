//! Standard input, one point or zone per line.

use std::io::{self, BufRead, Write};

use crate::Failure;

/// Hands each input line to `each`, and answers how many lines were bad. Blank lines and lines
/// beginning with `#` are skipped; the first other line, if `each` refuses it and it holds a
/// letter, is a header and skipped silently. A bad line is reported on `err` with its number, and
/// reading continues. Lines are read as bytes and converted lossily, so no byte stops the run, and
/// a trailing `\r` is dropped.
pub fn lines(
    input: &mut dyn BufRead,
    err: &mut dyn Write,
    each: &mut dyn FnMut(&str) -> Result<(), Failure>,
) -> Result<u64, Failure> {
    let (mut number, mut bad, mut seen_data) = (0u64, 0u64, false);
    let mut buf = Vec::new();
    loop {
        buf.clear();
        let read = input.read_until(b'\n', &mut buf).map_err(|e| {
            Failure::Io(io::Error::new(
                e.kind(),
                format!("cannot read standard input: {e}"),
            ))
        })?;
        if read == 0 {
            return Ok(bad);
        }
        number += 1;
        let text = String::from_utf8_lossy(&buf);
        let line = text.trim_end_matches(['\n', '\r']).trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let first = !seen_data;
        seen_data = true;
        match each(line) {
            Ok(()) => {}
            Err(Failure::Input(_)) if first && line.chars().any(|c| c.is_ascii_alphabetic()) => {}
            Err(Failure::Input(m) | Failure::Usage(m)) => {
                bad += 1;
                writeln!(err, "line {number}: {m}")?;
            }
            Err(e @ Failure::Io(_)) => return Err(e),
        }
    }
}
