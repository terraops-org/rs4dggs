//! The tool must never panic, whatever it is given: every command, on every grid, with hostile
//! arguments in every arity the command takes, and each hostile value as a line of batch input.
//! Each run must end with the exit code 0, 1 or 2.

use std::cell::Cell;
use std::io::Cursor;
use std::panic::{AssertUnwindSafe, catch_unwind};

/// Hostile values, of every kind the arguments read: numbers, points, identifiers, options.
fn hostile() -> Vec<String> {
    let mut v: Vec<String> = [
        "",
        "-",
        "nan,nan",
        "inf,0",
        "1e308,0",
        "90,180",
        "-90,-180",
        "(null)",
        "0x",
        "0xFFFFFFFFFFFFFFFF",
        "18446744073709551615",
        "99999999999999999999",
        "00",
        "11",
        "12",
        "A0-0-A",
        "L0-0-A",
        "A4-0-A",
        "é",
        "-depth",
        "255",
        "4294967295",
        "38.72,-9.14",
        "0",
        "5",
        "-1",
        "1e-320,1e-320",
        " ",
        "0x0",
        "0xffffffffffffffff",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    v.push("1".repeat(21));
    v.push("7".repeat(400));
    v
}

thread_local! {
    /// The runs of the test on this thread.
    static RUNS: Cell<usize> = const { Cell::new(0) };
    /// The grid the test on this thread is sweeping, named in its failures.
    static GRID: Cell<&'static str> = const { Cell::new("no grid") };
}

/// Runs the tool on `args` and `stdin`, which must neither panic nor answer any code but 0, 1
/// or 2.
fn once(args: &[&str], stdin: &[u8]) {
    RUNS.with(|r| r.set(r.get() + 1));
    let args: Vec<String> = args.iter().map(|a| a.to_string()).collect();
    let outcome = catch_unwind(AssertUnwindSafe(|| {
        let mut input = Cursor::new(stdin.to_vec());
        let (mut out, mut err) = (Vec::new(), Vec::new());
        rs4dggs_cli::run(&args, &mut input, &mut out, &mut err)
    }));
    match outcome {
        Ok(0..=2) => {}
        Ok(code) => panic!(
            "{}: exit code {code} for {args:?} on {:?}",
            GRID.with(Cell::get),
            String::from_utf8_lossy(stdin)
        ),
        Err(_) => panic!(
            "{}: panic for {args:?} on {:?}",
            GRID.with(Cell::get),
            String::from_utf8_lossy(stdin)
        ),
    }
}

/// Every command of `grid`, with the hostile values in every arity it takes.
fn sweep(grid: &'static str, valid_zone: &str) {
    GRID.with(|g| g.set(grid));
    let before = RUNS.with(Cell::get);
    let values = hostile();
    let none: &[u8] = b"";

    once(&[grid], none);
    once(&[grid, "info"], none);
    once(&[grid, "help"], none);
    for command in [
        "zone",
        "info",
        "geom",
        "neighbours",
        "neighbors",
        "disk",
        "rel",
        "sub",
        "index",
        "help",
    ] {
        // Too few arguments, and the command alone.
        once(&[grid, command], none);
    }
    for a in &values {
        for command in ["zone", "info", "geom", "neighbours", "sub", "help"] {
            once(&[grid, command, a], none);
        }
        // Output formats and options on their own, and as the argument of a command.
        once(&[grid, "info", valid_zone, "-f", a], none);
        once(&[grid, "info", valid_zone, "-precision", a], none);
        once(&[grid, "sub", valid_zone, "-depth", a], none);
        once(&[grid, "info", valid_zone, a], none);
        for b in &values {
            once(&[grid, "zone", a, b], none);
            // A radius of four thousand million about a fine zone asks for a walk over the
            // whole sphere, as it should; the coarsest zone below ends it at once.
            if b != "4294967295" {
                once(&[grid, "disk", a, b], none);
            }
            once(&[grid, "rel", a, b], none);
            once(&[grid, "index", a, b], none);
            once(&[grid, "sub", a, b], none);
        }
        // A zone against the hostile value, both ways round.
        once(&[grid, "rel", valid_zone, a], none);
        once(&[grid, "rel", a, valid_zone], none);
        once(&[grid, "index", valid_zone, a], none);
        once(&[grid, "index", a, valid_zone], none);
        once(&[grid, "disk", valid_zone, a], none);
        once(&[grid, "sub", valid_zone, a, "-depth", "2"], none);
        // A point as two arguments, and with a resolution.
        for r in ["0", "5", "19", "33", "34", "255"] {
            once(&[grid, "zone", a, r], none);
            once(&[grid, "zone", a, a, r], none);
        }
        once(&[grid, "zone", "10", a, "5"], none);
    }

    // Each hostile value as a batch line, in every format, for every command that reads a line.
    for a in &values {
        let line = format!("{a}\n");
        let pair = format!("{a},{a}\n");
        for format in ["text", "csv", "geojson"] {
            for (args, stdin) in [
                (vec!["zone", "-", "5"], &line),
                (vec!["zone", "-"], &line),
                (vec!["zone", "-", "5"], &pair),
                (vec!["info", "-"], &line),
                (vec!["geom", "-"], &line),
                (vec!["neighbours", "-"], &line),
                (vec!["disk", "-", "2"], &line),
                (vec!["rel", "-", valid_zone], &line),
                (vec!["sub", "-"], &line),
                (vec!["index", valid_zone, "-"], &line),
            ] {
                let mut full = vec![grid];
                full.extend(args);
                full.extend(["-f", format]);
                once(&full, stdin.as_bytes());
            }
        }
    }
    // Line endings, a last line without a newline, a header, comments, blank lines and bytes
    // that are not UTF-8.
    let mixed: Vec<&[u8]> = vec![
        b"38.72,-9.14\r\n-90,0\r\n",
        b"38.72,-9.14",
        b"lat,lon\n38.72,-9.14\n# a comment\n\n   \n1,2,3\n",
        b"\xff\xfe\n\xc3\x28\n",
        b"\0\0\0\n",
        b"\n\n\n",
        b"",
    ];
    for stdin in mixed {
        for command in [
            &["zone", "-", "5"][..],
            &["zone", "-"],
            &["info", "-"],
            &["geom", "-"],
            &["neighbours", "-"],
            &["disk", "-", "1"],
        ] {
            for format in ["text", "csv", "geojson"] {
                let mut full = vec![grid];
                full.extend(command);
                full.extend(["-f", format]);
                once(&full, stdin);
            }
        }
    }
    println!("{grid}: {} runs", RUNS.with(Cell::get) - before);
}

#[test]
fn igeo7_survives_hostile_input() {
    sweep("igeo7", "00");
}

#[test]
fn ivea7h_survives_hostile_input() {
    sweep("ivea7h", "00");
}

#[test]
fn rtea7h_survives_hostile_input() {
    sweep("rtea7h", "00");
}

#[test]
fn isea3h_survives_hostile_input() {
    sweep("isea3h", "A4-0-A");
}

#[test]
fn ivea3h_survives_hostile_input() {
    sweep("ivea3h", "A4-0-A");
}

#[test]
fn rtea3h_survives_hostile_input() {
    sweep("rtea3h", "A4-0-A");
}

/// What has no grid: the grid names, the options and the commands before any grid.
#[test]
fn hostile_invocations_without_a_grid() {
    let before = RUNS.with(Cell::get);
    let none: &[u8] = b"";
    // Output into a directory that does not exist, and a missing value for `-o`.
    once(&["igeo7", "info", "-o", "/no/such/directory/x"], none);
    once(&["igeo7", "info", "-o"], none);
    for a in hostile() {
        once(&[&a], none);
        once(&[&a, "info"], none);
        once(&["grids", &a], none);
        once(&["help", &a], none);
        once(&["-f", &a], none);
        once(&["--precision", &a], none);
        once(
            &["igeo7", "zone", "38.72,-9.14", "5", "-precision", &a],
            none,
        );
        once(&["igeo7", "sub", "00", "-depth", &a], none);
        once(&["igeo7", "info", "--centroids", &a], none);
    }
    println!("no grid: {} runs", RUNS.with(Cell::get) - before);
}
