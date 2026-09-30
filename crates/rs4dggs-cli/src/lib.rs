//! The `rs4dggs` command-line tool, kept as a library so that its tests run it in process. Its
//! one entry point is [`run`]; nothing else here is meant for other programs.
#![forbid(unsafe_code)]

mod args;
mod batch;
mod commands;
mod csv;
mod geojson;
mod rings;
mod text;

use std::io::{self, BufRead, Write};

/// What went wrong, and so the exit code a run ends with.
#[derive(Debug)]
pub(crate) enum Failure {
    /// A malformed invocation: exit code 2.
    Usage(String),
    /// Input that names something that is not so: exit code 1.
    Input(String),
    /// The answer could not be written.
    Io(io::Error),
}

impl From<io::Error> for Failure {
    fn from(e: io::Error) -> Self {
        Failure::Io(e)
    }
}

/// Runs the tool on `args` (the arguments after the program's name), reading `input` where a
/// command is given `-`, writing its answer to `out` (or to the file `-o` names) and its
/// complaints to `err`. Answers the exit code: 0, 1 or 2, as the README states.
pub fn run(
    args: &[String],
    input: &mut dyn BufRead,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> u8 {
    let result = args::parse(args).and_then(|inv| match &inv.options.output {
        Some(path) => {
            let mut file = LazyFile { path, file: None };
            let answered = commands::execute(&inv, input, &mut file, err);
            file.flush()?;
            answered
        }
        None => {
            let answered = commands::execute(&inv, input, out, err);
            out.flush()?;
            answered
        }
    });
    match result {
        Ok(code) => code,
        Err(Failure::Usage(m)) => {
            let _ = writeln!(err, "error: {m}");
            2
        }
        Err(Failure::Input(m)) => {
            let _ = writeln!(err, "error: {m}");
            1
        }
        // The reader has gone (`rs4dggs ... | head`): nothing more is wanted, and nothing failed.
        Err(Failure::Io(e)) if e.kind() == io::ErrorKind::BrokenPipe => 0,
        Err(Failure::Io(e)) => {
            let _ = writeln!(err, "error: {e}");
            1
        }
    }
}

/// The file `-o` names, created and truncated when the first byte is written, so that a command
/// that fails before it has an answer leaves an existing file as it was.
struct LazyFile<'a> {
    path: &'a std::path::Path,
    file: Option<io::BufWriter<std::fs::File>>,
}

impl Write for LazyFile<'_> {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        if self.file.is_none() {
            let file = std::fs::File::create(self.path).map_err(|e| {
                io::Error::new(
                    e.kind(),
                    format!("cannot write to {}: {e}", self.path.display()),
                )
            })?;
            self.file = Some(io::BufWriter::new(file));
        }
        self.file.as_mut().expect("just created").write(bytes)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.file.as_mut().map_or(Ok(()), Write::flush)
    }
}

/// The parsed invocation in a short form the tests compare: `zone <point> res <res>`.
#[doc(hidden)]
pub fn parse_for_tests(args: &[&str]) -> Result<String, String> {
    let args: Vec<String> = args.iter().map(|a| a.to_string()).collect();
    match args::parse(&args) {
        Ok(inv) => Ok(match inv.command {
            args::Command::Zone { point, res } => format!("zone {point} res {res:?}"),
            other => format!("{other:?}"),
        }),
        Err(Failure::Usage(m) | Failure::Input(m)) => Err(m),
        Err(Failure::Io(e)) => Err(e.to_string()),
    }
}
