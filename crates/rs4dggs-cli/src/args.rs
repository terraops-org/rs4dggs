//! The command line: which grid, which command, which options; and the help texts.

use std::path::PathBuf;

use rs4dggs::AnyGrid;

use crate::Failure;
use crate::commands::{facts, grid_by_name};

/// The output formats.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    Text,
    Csv,
    GeoJson,
}

/// The options, each already checked.
#[derive(Debug, Default)]
pub struct Options {
    pub format: Option<Format>,
    pub output: Option<PathBuf>,
    pub precision: Option<usize>,
    pub depth: Option<u8>,
    pub centroids: bool,
}

/// A command and its arguments, as typed: points and zones are read against the grid later.
#[derive(Debug, PartialEq, Eq)]
pub enum Command {
    Help(Option<String>),
    Version,
    Grids,
    GridInfo,
    /// `help <grid>`: the grid's card, which is help and is never refused for its options.
    HelpGrid,
    /// `point` is `lat,lon` or `-`; `res` absent means every resolution.
    Zone {
        point: String,
        res: Option<u8>,
    },
    Info(String),
    Geom(String),
    Neighbours(String),
    Disk {
        zone: String,
        k: u32,
    },
    Rel(String, String),
    Sub {
        zone: String,
        index: Option<u64>,
    },
    Index {
        parent: String,
        sub: String,
    },
}

#[derive(Debug)]
pub struct Invocation {
    pub grid: Option<AnyGrid>,
    pub command: Command,
    pub options: Options,
}

/// Whether `a` is an option. `-` alone is standard input, and a dash before a digit or a point
/// begins a negative number, so neither is an option.
fn is_option(a: &str) -> bool {
    a.len() > 1
        && a.starts_with('-')
        && !a[1..].starts_with(|c: char| c.is_ascii_digit() || c == '.')
}

/// A whole number from 0 to 255, as a resolution or a depth is.
fn small(v: &str) -> Option<u8> {
    v.parse().ok()
}

/// Two zones that exist on `grid`, for the examples: the library's own on aperture 7, and
/// DGGAL's on aperture 3, where two neighbouring zones are also of one resolution.
fn sample(grid: &AnyGrid) -> (&'static str, &'static str) {
    if facts(grid).aperture == 7 {
        ("006415654636", "006415654634")
    } else {
        ("C2-23-C", "C4-6-B")
    }
}

/// A zone of `grid`, a depth, an index and the sub-zone at that index of the zone's order at
/// that depth, for the examples of `sub` and `index`.
fn sub_sample(grid: &AnyGrid) -> (&'static str, u8, u64, &'static str) {
    if facts(grid).aperture == 7 {
        ("0064156", 2, 27, "006415600")
    } else {
        ("A4-0-A", 3, 8, "B2-5-C")
    }
}

pub fn parse(args: &[String]) -> Result<Invocation, Failure> {
    let mut options = Options::default();
    let mut words: Vec<&str> = Vec::new();
    let mut help = false;
    let mut version = false;
    let mut it = args.iter();
    while let Some(a) = it.next() {
        if !is_option(a) {
            words.push(a);
            continue;
        }
        let name = a.trim_start_matches('-');
        let mut value = |what: &str| -> Result<&String, Failure> {
            match it.next() {
                None => Err(Failure::Usage(format!("{a} needs {what} after it"))),
                Some(v) if is_option(v) => Err(Failure::Usage(format!(
                    "{a} needs {what} after it, not the option {v}"
                ))),
                Some(v) => Ok(v),
            }
        };
        match name {
            "h" | "help" => help = true,
            "V" | "version" => version = true,
            "centroids" => options.centroids = true,
            "f" | "format" => {
                options.format = Some(match value("a format")?.as_str() {
                    "text" => Format::Text,
                    "csv" => Format::Csv,
                    "geojson" | "json" => Format::GeoJson,
                    other => {
                        return Err(Failure::Usage(format!(
                            "{other:?} is not a format; use text, csv or geojson"
                        )));
                    }
                })
            }
            "o" | "output" => options.output = Some(PathBuf::from(value("a file name")?)),
            "precision" => {
                let v = value("a number of decimal places")?;
                options.precision =
                    Some(v.parse().ok().filter(|p: &usize| *p <= 17).ok_or_else(|| {
                        Failure::Usage(format!(
                            "-precision takes 0 to 17 decimal places, not {v:?}"
                        ))
                    })?)
            }
            "depth" => {
                let v = value("a depth")?;
                options.depth = Some(small(v).ok_or_else(|| {
                    Failure::Usage(format!(
                        "-depth takes a whole number of resolutions from 0 to 255, not {v:?}"
                    ))
                })?)
            }
            _ => {
                return Err(Failure::Usage(format!(
                    "unknown option {a}; `rs4dggs help` lists the options"
                )));
            }
        }
    }
    let usage = |m: String| Err(Failure::Usage(m));
    if version {
        return Ok(Invocation {
            grid: None,
            command: Command::Version,
            options,
        });
    }
    let command_of = |words: &[&str]| words.get(1).map(|w| w.to_string());
    if help || words.is_empty() || words[0] == "help" {
        let mut topic = if words.first() == Some(&"help") {
            command_of(&words)
        } else {
            words.get(1).map(|w| w.to_string())
        };
        if words.first() == Some(&"help") {
            // `help igeo7` is that grid's info; `help grids` and `help help` are the general help.
            if let Some(grid) = topic.as_deref().and_then(grid_by_name) {
                // `help igeo7 zone` is that command's help, on that grid.
                let command = match command_of(&words[1..]) {
                    Some(c) => Command::Help(Some(c)),
                    None => Command::HelpGrid,
                };
                return Ok(Invocation {
                    grid: Some(grid),
                    command,
                    options,
                });
            }
            if matches!(topic.as_deref(), Some("grids" | "help")) {
                topic = None;
            }
        }
        let grid = words.first().and_then(|w| grid_by_name(w));
        return Ok(Invocation {
            grid,
            command: Command::Help(topic),
            options,
        });
    }
    if words[0] == "grids" {
        if words.len() > 1 {
            return usage("`grids` takes no arguments".into());
        }
        return Ok(Invocation {
            grid: None,
            command: Command::Grids,
            options,
        });
    }
    let grid = grid_by_name(words[0]).ok_or_else(|| {
        Failure::Usage(format!(
            "unknown grid {:?}; the grids are igeo7, ivea7h, rtea7h, isea3h, ivea3h and rtea3h \
             (`rs4dggs grids` describes them)",
            words[0]
        ))
    })?;
    let rest = &words[2.min(words.len())..];
    let (z1, z2) = sample(&grid);
    let (sz, sdepth, _, ssub) = sub_sample(&grid);
    let g = words[0];
    let arg = |i: usize, example: String| -> Result<String, Failure> {
        rest.get(i)
            .map(|w| w.to_string())
            .ok_or_else(|| Failure::Usage(format!("missing argument; for example: {example}")))
    };
    let at_most = |n: usize, example: String| -> Result<(), Failure> {
        if rest.len() > n {
            Err(Failure::Usage(format!(
                "too many arguments (unwanted: {}); for example: {example}",
                rest[n..].join(" ")
            )))
        } else {
            Ok(())
        }
    };
    let command = match words.get(1).copied() {
        None => Command::GridInfo,
        Some("info") if rest.is_empty() => Command::GridInfo,
        Some("help") => Command::Help(rest.first().map(|w| w.to_string())),
        Some("zone") => {
            let (point, after) = match rest {
                // A point pasted from a web map, `38.72, -9.14`, reaches us as two words; a
                // complete point with a trailing comma, `38.72,-9.14,`, is one.
                [p, lon, ..] if p.ends_with(',') && p.matches(',').count() == 1 => {
                    (format!("{p}{lon}"), 2)
                }
                [p, ..] if p.contains(',') || *p == "-" => (p.to_string(), 1),
                [lat, lon, ..] => (format!("{lat},{lon}"), 2),
                _ => {
                    return usage(format!(
                        "zone needs a point, as 38.72,-9.14 or 38.72 -9.14; for example: \
                         rs4dggs {g} zone 38.72,-9.14 10"
                    ));
                }
            };
            at_most(after + 1, format!("rs4dggs {g} zone 38.72,-9.14 10"))?;
            let res = match rest.get(after) {
                None => None,
                Some(r) => Some(small(r).ok_or_else(|| {
                    Failure::Usage(format!(
                        "{r:?} is not a resolution; a resolution is a whole number from 0 to 255 \
                         (`rs4dggs grids` gives each grid's finest)"
                    ))
                })?),
            };
            Command::Zone { point, res }
        }
        Some("info") => {
            let e = format!("rs4dggs {g} info {z1}");
            at_most(1, e.clone())?;
            Command::Info(arg(0, e)?)
        }
        Some("geom") => {
            let e = format!("rs4dggs {g} geom {z1}");
            at_most(1, e.clone())?;
            Command::Geom(arg(0, e)?)
        }
        Some("neighbours" | "neighbors") => {
            let e = format!("rs4dggs {g} neighbours {z1}");
            at_most(1, e.clone())?;
            Command::Neighbours(arg(0, e)?)
        }
        Some("disk") => {
            let e = format!("rs4dggs {g} disk {z1} 2");
            at_most(2, e.clone())?;
            let zone = arg(0, e.clone())?;
            let k = arg(1, e)?;
            let k = k.parse().map_err(|_| {
                Failure::Usage(format!(
                    "{k:?} is not a radius; a radius is a whole number from 0"
                ))
            })?;
            Command::Disk { zone, k }
        }
        Some("rel") => {
            let e = format!("rs4dggs {g} rel {z1} {z2}");
            at_most(2, e.clone())?;
            Command::Rel(arg(0, e.clone())?, arg(1, e)?)
        }
        Some("sub") => {
            let e = format!("rs4dggs {g} sub {sz} -depth {sdepth}");
            at_most(2, e.clone())?;
            let zone = arg(0, e)?;
            let index = match rest.get(1) {
                None => None,
                Some(i) => Some(i.parse().map_err(|_| {
                    Failure::Usage(format!(
                        "{i:?} is not an index; an index is a whole number from 0"
                    ))
                })?),
            };
            Command::Sub { zone, index }
        }
        Some("index") => {
            if options.depth.is_some() {
                return usage(
                    "index takes its depth from the two zones' resolutions, so -depth is refused"
                        .into(),
                );
            }
            let e = format!("rs4dggs {g} index {sz} {ssub}");
            at_most(2, e.clone())?;
            Command::Index {
                parent: arg(0, e.clone())?,
                sub: arg(1, e)?,
            }
        }
        Some(other) => {
            return usage(format!(
                "unknown command {other:?}; `rs4dggs help` lists the commands"
            ));
        }
    };
    Ok(Invocation {
        grid: Some(grid),
        command,
        options,
    })
}

/// The general help.
pub const HELP: &str = "\
rs4dggs: DGGAL's discrete global grids from the shell

Usage:
  rs4dggs <grid> <command> [options] [arguments]
  rs4dggs grids
  rs4dggs help
  rs4dggs <grid> help <command>

Grids:
  igeo7, ivea7h, rtea7h (aperture 7), isea3h, ivea3h, rtea3h (aperture 3), in any case;
  DGGAL's names isea7h_z7, ivea7h_z7 and rtea7h_z7 are accepted too.

Commands:
  info                   the grid itself
  zone <lat,lon> [res]   the zone containing a point; without a resolution, at every one
  info <zone>            a zone's card: identifiers, centroid, parents, children, neighbours,
                         vertices
  geom <zone>            a zone's boundary, as GeoJSON
  neighbours <zone>      a zone's neighbours
  disk <zone> <k>        every zone within k steps, ring by ring
  rel <zone> <zone>      how two zones are related
  sub <zone> [index]     the sub-zones at -depth, in DGGAL's order, or the one at that index
  index <zone> <zone>    the index of a sub-zone within a zone

Options:
  -f, --format text|csv|geojson   the output format
  -o, --output <file>             write to a file instead of the terminal
  -precision <n>                  decimal places for coordinates (default: every digit)
  -depth <n>                      the relative depth for sub (default: 1)
  -centroids                      GeoJSON points at centroids instead of polygons

A point is latitude and longitude in degrees, as 38.72,-9.14 or 38.72 -9.14.
A zone is its text identifier, its 64-bit identifier in decimal, or in hexadecimal after 0x
(a decimal identifier that also reads as text is read as text; write 0x to be certain).
Give - in place of a point or a zone to read one per line from standard input; a first line
that does not parse and holds a letter is taken as a header and skipped.

Examples:
  rs4dggs igeo7 zone 38.72,-9.14 10
  rs4dggs igeo7 info 006415654636
  rs4dggs isea3h geom C2-23-C > cell.geojson
  rs4dggs igeo7 zone - 12 < points.csv > zones.csv
";

/// The help on one command, with examples on `grid`; `None` for a word that is no command.
pub fn command_help(command: &str, grid: &AnyGrid) -> Option<String> {
    let g = grid.name().to_ascii_lowercase();
    let (z1, z2) = sample(grid);
    let (sz, sdepth, sindex, ssub) = sub_sample(grid);
    // What DGGAL's hierarchy is on the grids of this aperture.
    let (parents, children) = if facts(grid).aperture == 7 {
        ("one parent or two", "thirteen children (eleven")
    } else {
        ("one parent or three", "seven children (six")
    };
    // A head and its description, and a second line of the description under the first.
    let two =
        |head: String, l1: &str, l2: &str| format!("{head}{l1}\n{}{l2}\n", " ".repeat(head.len()));
    let text = match command {
        "info" => format!(
            "rs4dggs {g} info             the grid: aperture, resolutions, projection, indexing\n\
             rs4dggs {g} info <zone>      a zone's card\n\n\
             The card gives the zone's text and 64-bit identifiers, its resolution and shape, its\n\
             centroid, its parents and children, its neighbours and its vertices.\n\n\
             The parents and children are DGGAL's own: on this grid a zone has {parents}, and\n\
             {children} under a pentagon, none at the finest resolution). The parent that\n\
             is itself a centroid child is marked, and so is the child at the zone's centre.\n\n\
             With `-` for the zone, one is read per line from standard input; a first line that\n\
             does not parse and holds a letter is taken as a header and skipped.\n\n\
             Examples:\n  rs4dggs {g} info {z1}\n  rs4dggs {g} info {z1} -f geojson\n  \
             rs4dggs {g} info - < zones.txt > cards.csv\n"
        ),
        "zone" => format!(
            "rs4dggs {g} zone <lat,lon> <res>   the zone containing the point, as a card\n\
             rs4dggs {g} zone <lat,lon>         that point's zone at every resolution\n\n\
             With `-` for the point, one is read per line from standard input; a first line that\n\
             does not parse and holds a letter is taken as a header and skipped.\n\n\
             Examples:\n  rs4dggs {g} zone 38.72,-9.14 10\n  rs4dggs {g} zone 38.72 -9.14\n  \
             rs4dggs {g} zone - 12 < points.csv > zones.csv\n"
        ),
        "geom" => format!(
            "{}\nExamples:\n  rs4dggs {g} geom {z1} > cell.geojson\n  \
             rs4dggs {g} geom - < zones.txt > cells.geojson\n",
            two(
                format!("rs4dggs {g} geom <zone>   "),
                "the boundary as GeoJSON, drawn correctly across the antimeridian",
                "and about the poles"
            )
        ),
        "neighbours" | "neighbors" => format!(
            "rs4dggs {g} neighbours <zone>   the zones sharing an edge with it, in DGGAL's \
             order\n\n\
             Examples:\n  rs4dggs {g} neighbours {z1}\n  \
             rs4dggs {g} neighbours - -f csv < zones.txt\n"
        ),
        "disk" => format!(
            "rs4dggs {g} disk <zone> <k>   every zone within k neighbour steps, ring by ring, the \
             zone first\n\n\
             Examples:\n  rs4dggs {g} disk {z1} 2\n  \
             rs4dggs {g} disk {z1} 3 -f geojson > disk.geojson\n"
        ),
        "rel" => format!(
            "{}\nExample:\n  rs4dggs {g} rel {z1} {z2}\n",
            two(
                format!("rs4dggs {g} rel <zone> <zone>   "),
                "how the two zones are related: resolutions, neighbours,",
                "parent and child, ancestor and descendant, siblings"
            )
        ),
        "sub" => format!(
            "rs4dggs {g} sub <zone> [-depth n]     the sub-zones n resolutions below (default 1), \
             in order\n\
             rs4dggs {g} sub <zone> <index> [-depth n]   the one at that index\n\n\
             The order is DGGAL's own. A position printed as (no cell) is one at which the grid\n\
             names no zone. A list longer than the library's limit is refused; one sub-zone by\n\
             its index is answered at any depth.\n\n\
             Examples:\n  rs4dggs {g} sub {sz} -depth {sdepth}\n  \
             rs4dggs {g} sub {sz} {sindex} -depth {sdepth}\n"
        ),
        "index" => format!(
            "{}\nExample:\n  rs4dggs {g} index {sz} {ssub}\n",
            two(
                format!("rs4dggs {g} index <zone> <sub-zone>   "),
                "the sub-zone's index within the zone, at the depth",
                "their resolutions set"
            )
        ),
        _ => return None,
    };
    Some(text)
}
