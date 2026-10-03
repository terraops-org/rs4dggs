//! Each command: asks the library, and writes the answer.

use std::io::{BufRead, Write};

use rs4dggs::{AnyGrid, Error, GeoPoint, NULL_TEXT, ZoneId, get_grid};

use crate::args::{Command, Format, HELP, Invocation, command_help};
use crate::batch;
use crate::csv::Csv;
use crate::geojson::Collection;
use crate::{Failure, text};

/// What the tool says about a grid beyond what the library reports.
pub struct GridFacts {
    pub name: &'static str,
    pub aperture: u8,
    pub projection: &'static str,
    pub projection_name: &'static str,
    pub indexing: &'static str,
}

pub const FACTS: [GridFacts; 6] = [
    GridFacts {
        name: "IGEO7",
        aperture: 7,
        projection: "ISEA",
        projection_name: "icosahedral Snyder equal-area",
        indexing: "Z7",
    },
    GridFacts {
        name: "IVEA7H",
        aperture: 7,
        projection: "IVEA",
        projection_name: "icosahedral vertex-oriented great-circle equal-area",
        indexing: "Z7",
    },
    GridFacts {
        name: "RTEA7H",
        aperture: 7,
        projection: "RTEA",
        projection_name: "rhombic triacontahedral equal-area",
        indexing: "Z7",
    },
    GridFacts {
        name: "ISEA3H",
        aperture: 3,
        projection: "ISEA",
        projection_name: "icosahedral Snyder equal-area",
        indexing: "I3H",
    },
    GridFacts {
        name: "IVEA3H",
        aperture: 3,
        projection: "IVEA",
        projection_name: "icosahedral vertex-oriented great-circle equal-area",
        indexing: "I3H",
    },
    GridFacts {
        name: "RTEA3H",
        aperture: 3,
        projection: "RTEA",
        projection_name: "rhombic triacontahedral equal-area",
        indexing: "I3H",
    },
];

pub fn facts(grid: &AnyGrid) -> &'static GridFacts {
    FACTS
        .iter()
        .find(|f| f.name == grid.name())
        .expect("every grid of the library is in FACTS")
}

/// The grid a name gives, the library's names in any case, and DGGAL's for the aperture-7 grids.
pub fn grid_by_name(name: &str) -> Option<AnyGrid> {
    let name = match name.to_ascii_lowercase().as_str() {
        "isea7h_z7" => "IGEO7".to_string(),
        "ivea7h_z7" => "IVEA7H".to_string(),
        "rtea7h_z7" => "RTEA7H".to_string(),
        _ => name.to_string(),
    };
    get_grid(&name).ok()
}

/// A coordinate as written: its shortest exact form (the fewest digits that read back to the
/// same double, which is what `Display` writes), or `precision` decimal places.
pub fn coord(x: f64, precision: Option<usize>) -> String {
    match precision {
        Some(p) => {
            let s = format!("{x:.p$}");
            // A small negative rounds to `-0.000`, which is no better than `0.000`.
            match s.strip_prefix('-') {
                Some(t) if t.bytes().all(|b| b == b'0' || b == b'.') => t.to_string(),
                _ => s,
            }
        }
        None => format!("{x}"),
    }
}

/// A point, `lat,lon` or `lat lon` (a comma, spaces or a tab between), in range.
pub fn parse_point(s: &str) -> Result<(f64, f64), String> {
    let fields: Vec<&str> = s
        .split(|c: char| c == ',' || c.is_whitespace())
        .filter(|f| !f.is_empty())
        .collect();
    let [lat, lon] = fields.as_slice() else {
        return Err(format!(
            "{s:?} is not a point; write latitude and longitude, as 38.72,-9.14"
        ));
    };
    let number = |f: &str, what: &str| -> Result<f64, String> {
        f.parse::<f64>()
            .ok()
            .filter(|v| v.is_finite())
            .ok_or_else(|| format!("{f:?} is not a {what}"))
    };
    let (lat, lon) = (number(lat, "latitude")?, number(lon, "longitude")?);
    if !(-90.0..=90.0).contains(&lat) {
        return Err(format!("latitude {lat} is outside -90 to 90"));
    }
    if !(-180.0..=180.0).contains(&lon) {
        return Err(format!("longitude {lon} is outside -180 to 180"));
    }
    Ok((lat, lon))
}

/// A zone of `grid`, from its text, its decimal identifier, or its hexadecimal one after `0x`.
/// An identifier is accepted only if the grid reads its text back to it.
pub fn parse_zone(grid: &AnyGrid, s: &str) -> Result<ZoneId, String> {
    let s = s.trim();
    let name = grid.name();
    if s == NULL_TEXT {
        return Err(format!("{NULL_TEXT} is the null zone, which names no cell"));
    }
    let checked = |id: ZoneId| -> Result<ZoneId, String> {
        if id == ZoneId::NULL {
            return Err(format!("{s} is the null zone, which names no cell"));
        }
        match grid.zone_from_text(&grid.text_id(id)) {
            Ok(back) if back == id => Ok(id),
            _ => Err(format!("{s} is not a zone of {name}")),
        }
    };
    if let Some(hex) = s.strip_prefix("0x").or_else(|| s.strip_prefix("0X")) {
        let id = u64::from_str_radix(hex, 16)
            .map_err(|_| format!("{s:?} is not a hexadecimal identifier"))?;
        return checked(ZoneId(id));
    }
    match grid.zone_from_text(s) {
        Ok(id) if id != ZoneId::NULL => Ok(id),
        Ok(_) => Err(format!("{s} is the null zone, which names no cell")),
        Err(e) => match s.parse::<u64>() {
            Ok(id) if s.bytes().all(|b| b.is_ascii_digit()) => checked(ZoneId(id)),
            _ => Err(format!("{s:?} is not a zone of {name}: {e}")),
        },
    }
}

pub struct Named {
    pub text: String,
    pub mark: bool,
}

pub struct Card {
    pub grid: &'static str,
    pub id: ZoneId,
    pub text: String,
    pub resolution: u8,
    pub pentagon: bool,
    pub centroid: GeoPoint,
    /// The centroid parent marked.
    pub parents: Vec<Named>,
    /// The centroid child marked.
    pub children: Vec<Named>,
    pub neighbours: Vec<String>,
    pub vertices: Vec<GeoPoint>,
}

pub fn card(grid: &AnyGrid, id: ZoneId) -> Card {
    let centroid_parent = grid.centroid_parent(id);
    let parents = grid.parents(id);
    Card {
        grid: grid.name(),
        id,
        text: grid.text_id(id),
        resolution: grid.resolution(id),
        pentagon: grid.is_pentagon(id),
        centroid: grid.centroid(id),
        parents: parents
            .into_iter()
            .map(|p| Named {
                text: grid.text_id(p),
                // The parent that is itself a centroid child, a lone parent included. On the
                // aperture-3 grids `dgg` marks the same one; on the aperture-7 grids it marks
                // none, since DGGAL answers no centroid parent there, and the library's is marked.
                mark: Some(p) == centroid_parent,
            })
            .collect(),
        // The child at the zone's centre: a centroid child whose parent is the zone. Along the
        // two broken seams of the aperture-7 grids a child may be the centroid child of another
        // zone, and is not marked, as `dgg` does not mark it.
        children: grid
            .children(id)
            .into_iter()
            .map(|c| Named {
                text: grid.text_id(c),
                mark: grid.is_centroid_child(c) && grid.parent(c) == Some(id),
            })
            .collect(),
        neighbours: grid
            .neighbors(id)
            .into_iter()
            .map(|n| grid.text_id(n))
            .collect(),
        vertices: grid.vertices(id),
    }
}

/// A library refusal, as a message for the person who typed the arguments.
fn refusal(grid: &AnyGrid, e: Error) -> Failure {
    let name = grid.name();
    Failure::Input(match e {
        Error::ResolutionOutOfRange { res, max } => {
            format!("resolution {res} is beyond {name}'s finest, {max}")
        }
        // No grid of the library answers this today: all six order their sub-zones.
        Error::NoSubZoneOrder => format!("{name} has no sub-zone order"),
        Error::TooManySubZones { count, limit } => format!(
            "{count} sub-zones are more than the limit of {limit}; ask for a smaller -depth, or \
             for one sub-zone by its index"
        ),
        e => e.to_string(),
    })
}

/// What the text output prints in place of a zone where the grid has none: the null zone.
const NO_CELL: &str = "(no cell)";

/// The text of a zone, or [`NO_CELL`] for the null zone, which names no cell.
fn text_or_no_cell(grid: &AnyGrid, id: ZoneId) -> String {
    if id == ZoneId::NULL {
        NO_CELL.to_string()
    } else {
        grid.text_id(id)
    }
}

fn zone_arg(grid: &AnyGrid, s: &str) -> Result<ZoneId, Failure> {
    parse_zone(grid, s).map_err(Failure::Input)
}

/// The name of a command, for the messages that concern it.
fn command_name(c: &Command) -> &'static str {
    match c {
        Command::Zone { .. } => "zone",
        Command::Info(_) => "info",
        Command::Geom(_) => "geom",
        Command::Neighbours(_) => "neighbours",
        Command::Disk { .. } => "disk",
        Command::Rel(..) => "rel",
        Command::Sub { .. } => "sub",
        Command::Index { .. } => "index",
        _ => "this command",
    }
}

fn reads_standard_input(c: &Command) -> bool {
    match c {
        Command::Zone { point: a, .. }
        | Command::Info(a)
        | Command::Geom(a)
        | Command::Neighbours(a)
        | Command::Disk { zone: a, .. }
        | Command::Sub { zone: a, .. } => a == "-",
        Command::Rel(a, b) => a == "-" || b == "-",
        Command::Index { parent, sub } => parent == "-" || sub == "-",
        _ => false,
    }
}

pub fn execute(
    inv: &Invocation,
    input: &mut dyn BufRead,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> Result<u8, Failure> {
    match &inv.command {
        Command::Help(None) => write!(out, "{HELP}")?,
        Command::Help(Some(c)) => {
            let grid = inv
                .grid
                .unwrap_or_else(|| grid_by_name("igeo7").expect("a registered grid"));
            match command_help(c, &grid) {
                Some(text) => write!(out, "{text}")?,
                None => {
                    return Err(Failure::Usage(format!(
                        "unknown command {c:?}; `rs4dggs help` lists the commands"
                    )));
                }
            }
        }
        // The library's version cannot be had without adding to its surface, which is frozen, so
        // the tool's own is given; the two share the workspace's version.
        Command::Version => writeln!(out, "rs4dggs-cli {}", env!("CARGO_PKG_VERSION"))?,
        Command::Grids => {
            text_only(inv, "grids writes text")?;
            no_coordinate(inv, "grids")?;
            writeln!(out, "grid    aperture  resolutions  projection  indexing")?;
            for f in &FACTS {
                let g = grid_by_name(f.name).expect("a registered grid");
                writeln!(
                    out,
                    "{:<7} {:<9} {:<12} {:<11} {}",
                    f.name,
                    f.aperture,
                    format!("0 to {}", g.max_resolution()),
                    f.projection,
                    f.indexing
                )?;
            }
        }
        Command::GridInfo | Command::HelpGrid => {
            if inv.command == Command::GridInfo {
                text_only(
                    inv,
                    "the grid's info writes text; `info <zone>` writes CSV and GeoJSON too",
                )?;
                no_coordinate(inv, "the grid's info")?;
            }
            let grid = inv.grid.expect("a grid command has a grid");
            let f = facts(&grid);
            writeln!(out, "{}", f.name)?;
            writeln!(
                out,
                "  aperture     {} (each zone has about {} times the area of its children)",
                f.aperture, f.aperture
            )?;
            writeln!(out, "  resolutions  0 to {}", grid.max_resolution())?;
            writeln!(
                out,
                "  projection   {} ({})",
                f.projection, f.projection_name
            )?;
            writeln!(out, "  indexing     {}", f.indexing)?;
        }
        command => {
            let grid = inv.grid.as_ref().expect("a grid command has a grid");
            let name = command_name(command);
            let batch = reads_standard_input(command);
            let (allowed, says) = match command {
                Command::Zone { .. }
                | Command::Info(_)
                | Command::Neighbours(_)
                | Command::Disk { .. } => (
                    &[Format::Text, Format::Csv, Format::GeoJson][..],
                    "text, csv and geojson",
                ),
                Command::Geom(_) => (&[Format::GeoJson][..], "geojson"),
                _ => (&[Format::Text][..], "text"),
            };
            let format = inv.options.format.unwrap_or(match command {
                Command::Geom(_) => Format::GeoJson,
                _ if batch => Format::Csv,
                _ => Format::Text,
            });
            if batch && allowed == [Format::Text] {
                return Err(Failure::Usage(format!(
                    "{name} does not read standard input"
                )));
            }
            if allowed == [Format::Text] {
                no_coordinate(inv, name)?;
            }
            if inv.options.centroids && format != Format::GeoJson {
                return Err(Failure::Usage(
                    "-centroids applies to GeoJSON output only; add -f geojson".into(),
                ));
            }
            if !allowed.contains(&format) {
                return Err(Failure::Usage(format!("{name} writes {says}")));
            }
            let mut sink = match format {
                Format::Text => Sink::Text {
                    out,
                    written: false,
                },
                Format::Csv => Sink::Csv(Csv::new(
                    out,
                    match command {
                        Command::Zone { res: Some(_), .. } => "lat,lon,zone",
                        Command::Zone { .. } => "lat,lon,resolution,zone",
                        Command::Info(_) => {
                            "zone,resolution,shape,centroid_lat,centroid_lon,parent"
                        }
                        Command::Neighbours(_) => "zone,neighbour",
                        _ => "zone,ring,member",
                    },
                )),
                Format::GeoJson => Sink::Geo(Collection::new(
                    out,
                    inv.options.precision,
                    inv.options.centroids,
                )),
            };
            let mut bad = 0;
            if batch {
                let mut each = |line: &str| answer(grid, inv, line, &mut sink, true);
                bad = batch::lines(input, err, &mut each)?;
            } else {
                answer(grid, inv, argument(command), &mut sink, false)?;
            }
            match &mut sink {
                Sink::Text { .. } => {}
                Sink::Csv(c) => c.begin()?,
                Sink::Geo(g) => g.end(err)?,
            }
            return Ok(u8::from(bad > 0));
        }
    }
    Ok(0)
}

/// Refuses `-f csv` and `-f geojson` for a command that writes text alone.
fn text_only(inv: &Invocation, message: &str) -> Result<(), Failure> {
    match inv.options.format {
        None | Some(Format::Text) => Ok(()),
        Some(_) => Err(Failure::Usage(message.to_string())),
    }
}

/// Refuses `-centroids` and `-precision` for a command that prints no coordinate.
fn no_coordinate(inv: &Invocation, what: &str) -> Result<(), Failure> {
    if inv.options.centroids || inv.options.precision.is_some() {
        return Err(Failure::Usage(format!(
            "{what} prints no coordinate, so -centroids and -precision do not apply to it"
        )));
    }
    Ok(())
}

/// `n` resolutions, or `1 resolution`.
fn resolutions(n: u8) -> String {
    format!("{n} resolution{}", if n == 1 { "" } else { "s" })
}

/// The point or zone a one-input command is given.
fn argument(c: &Command) -> &str {
    match c {
        Command::Zone { point: a, .. }
        | Command::Info(a)
        | Command::Geom(a)
        | Command::Neighbours(a)
        | Command::Disk { zone: a, .. }
        | Command::Sub { zone: a, .. } => a,
        Command::Rel(a, _) => a,
        Command::Index { parent, .. } => parent,
        _ => "",
    }
}

/// Where the answers go.
enum Sink<'a> {
    /// In batch, answers are separated by a blank line.
    Text {
        out: &'a mut dyn Write,
        written: bool,
    },
    Csv(Csv<'a>),
    Geo(Collection<'a>),
}

/// Writes to `out`, putting a blank line before the first byte if one is due.
struct Separated<'a> {
    out: &'a mut dyn Write,
    blank_line_first: bool,
    started: bool,
}

impl Write for Separated<'_> {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        if !self.started && !bytes.is_empty() {
            self.started = true;
            if self.blank_line_first {
                self.out.write_all(b"\n")?;
            }
        }
        self.out.write(bytes)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.out.flush()
    }
}

/// Answers one command for one input, `arg` being its point or zone.
fn answer(
    grid: &AnyGrid,
    inv: &Invocation,
    arg: &str,
    sink: &mut Sink,
    batch: bool,
) -> Result<(), Failure> {
    if let Sink::Text { out, written } = sink {
        // The answer streams; only its separator waits, for the first byte, so that an input
        // refused writes nothing at all.
        let mut first = Separated {
            out: &mut **out,
            blank_line_first: batch && *written,
            started: false,
        };
        let mut alone = Sink::Text {
            out: &mut first,
            written: false,
        };
        let answered = write_answer(grid, inv, arg, &mut alone);
        *written |= first.started;
        return answered;
    }
    write_answer(grid, inv, arg, sink)
}

/// Writes the answer to `sink` as it is made.
fn write_answer(
    grid: &AnyGrid,
    inv: &Invocation,
    arg: &str,
    sink: &mut Sink,
) -> Result<(), Failure> {
    let precision = inv.options.precision;
    let name = grid.name();
    let text_id = |id: ZoneId| {
        if id == ZoneId::NULL {
            String::new()
        } else {
            grid.text_id(id)
        }
    };
    match &inv.command {
        Command::Zone { res, .. } => {
            let (lat, lon) = parse_point(arg).map_err(Failure::Input)?;
            let rows = match res {
                Some(r) => vec![(
                    *r,
                    grid.zone_from_geo(lat, lon, *r)
                        .map_err(|e| refusal(grid, e))?,
                )],
                None => (0..=grid.max_resolution())
                    .map(|r| grid.zone_from_geo(lat, lon, r).map(|id| (r, id)))
                    .collect::<Result<Vec<_>, _>>()
                    .map_err(|e| refusal(grid, e))?,
            };
            let (slat, slon) = (coord(lat, precision), coord(lon, precision));
            match sink {
                Sink::Text { out, .. } => match res {
                    Some(r) => {
                        let id = rows[0].1;
                        if id == ZoneId::NULL {
                            writeln!(
                                out,
                                "{name}: no zone contains {slat}, {slon} at resolution {r}; the \
                                 grid assigns no cell there (see the rs4dggs library's README, \
                                 \"What a caller must know\", \"The null zone is an answer\")"
                            )?;
                        } else {
                            text::card(out, &card(grid, id), precision)?;
                        }
                    }
                    None => {
                        let rows: Vec<(u8, String)> = rows
                            .into_iter()
                            .map(|(r, id)| (r, text_or_no_cell(grid, id)))
                            .collect();
                        text::zones_at(out, name, (lat, lon), &rows, precision)?;
                    }
                },
                Sink::Csv(c) => {
                    for (r, id) in rows {
                        let z = text_id(id);
                        if res.is_some() {
                            c.row(&[&slat, &slon, &z])?;
                        } else {
                            c.row(&[&slat, &slon, &r.to_string(), &z])?;
                        }
                    }
                }
                Sink::Geo(g) => {
                    let extra = format!(r#","lat":{slat},"lon":{slon}"#);
                    for (_, id) in rows {
                        g.zone(grid, id, &extra)?;
                    }
                }
            }
        }
        Command::Info(_) => {
            let id = zone_arg(grid, arg)?;
            match sink {
                Sink::Text { out, .. } => text::card(out, &card(grid, id), precision)?,
                Sink::Csv(c) => {
                    let (centroid, parents) = (grid.centroid(id), grid.parents(id));
                    let parents: Vec<String> =
                        parents.into_iter().map(|p| grid.text_id(p)).collect();
                    c.row(&[
                        &grid.text_id(id),
                        &grid.resolution(id).to_string(),
                        if grid.is_pentagon(id) {
                            "pentagon"
                        } else {
                            "hexagon"
                        },
                        &coord(centroid.lat, precision),
                        &coord(centroid.lon, precision),
                        &parents.join(" "),
                    ])?;
                }
                Sink::Geo(g) => g.zone(grid, id, "")?,
            }
        }
        Command::Geom(_) => {
            let id = zone_arg(grid, arg)?;
            if let Sink::Geo(g) = sink {
                g.zone(grid, id, "")?;
            }
        }
        Command::Neighbours(_) => {
            let id = zone_arg(grid, arg)?;
            let list = grid.neighbors(id);
            match sink {
                Sink::Text { out, .. } => {
                    let list: Vec<String> = list.into_iter().map(|n| grid.text_id(n)).collect();
                    text::neighbours(out, name, &grid.text_id(id), &list)?;
                }
                Sink::Csv(c) => {
                    let z = grid.text_id(id);
                    for n in list {
                        c.row(&[&z, &grid.text_id(n)])?;
                    }
                }
                Sink::Geo(g) => {
                    let extra = format!(r#","neighbour_of":"{}""#, grid.text_id(id));
                    for n in list {
                        g.zone(grid, n, &extra)?;
                    }
                }
            }
        }
        Command::Disk { k, .. } => {
            let id = zone_arg(grid, arg)?;
            let z = grid.text_id(id);
            if let Sink::Text { out, .. } = sink {
                text::disk_head(out, name, *k, &z)?;
            }
            for (i, ring) in grid.disk(id, *k).rings().enumerate() {
                match sink {
                    Sink::Text { out, .. } => {
                        let ids: Vec<String> = ring.into_iter().map(|m| grid.text_id(m)).collect();
                        text::ring(out, i, &ids)?;
                    }
                    Sink::Csv(c) => {
                        for m in ring {
                            c.row(&[&z, &i.to_string(), &grid.text_id(m)])?;
                        }
                    }
                    Sink::Geo(g) => {
                        let extra = format!(r#","ring":{i},"centre":"{z}""#);
                        for m in ring {
                            g.zone(grid, m, &extra)?;
                        }
                    }
                }
            }
        }
        _ => {
            let Sink::Text { out, .. } = sink else {
                unreachable!("the other commands write text alone");
            };
            other(grid, inv, out)?;
        }
    }
    Ok(())
}

/// Answers `rel`, `sub` and `index`, which write text alone.
fn other(grid: &AnyGrid, inv: &Invocation, out: &mut dyn Write) -> Result<(), Failure> {
    let name = grid.name();
    match &inv.command {
        Command::Rel(a, b) => {
            let (ia, ib) = (zone_arg(grid, a)?, zone_arg(grid, b)?);
            let (ta, tb) = (grid.text_id(ia), grid.text_id(ib));
            let mut lines = Vec::new();
            if ia == ib {
                lines.push(format!("{ta} and {tb} are the same zone"));
            } else {
                let (ra, rb) = (grid.resolution(ia), grid.resolution(ib));
                lines.push(match ra.cmp(&rb) {
                    std::cmp::Ordering::Less => {
                        format!("{ta} is coarser than {tb} by {}", resolutions(rb - ra))
                    }
                    std::cmp::Ordering::Greater => {
                        format!("{ta} is finer than {tb} by {}", resolutions(ra - rb))
                    }
                    std::cmp::Ordering::Equal => format!("{ta} and {tb} have the same resolution"),
                });
                lines.push(if grid.neighbors(ia).contains(&ib) {
                    format!("{ta} and {tb} are neighbours")
                } else {
                    format!("{ta} and {tb} are not neighbours")
                });
                if grid.is_immediate_child_of(ib, ia) {
                    lines.push(format!("{ta} is an immediate parent of {tb}"));
                }
                if grid.is_immediate_child_of(ia, ib) {
                    lines.push(format!("{ta} is an immediate child of {tb}"));
                }
                let ancestor = grid.is_ancestor_of(ia, ib);
                if ancestor {
                    lines.push(format!("{ta} is an ancestor of {tb}"));
                }
                let descendant = grid.is_ancestor_of(ib, ia);
                if descendant {
                    lines.push(format!("{ta} is a descendant of {tb}"));
                }
                if grid.is_sibling_of(ia, ib) {
                    lines.push(format!("{ta} and {tb} are siblings"));
                }
                // Where the library refuses the index, the order being longer than it lists,
                // nothing is said.
                if ancestor {
                    if let Ok(Some(i)) = grid.sub_zone_index(ia, ib) {
                        lines.push(format!(
                            "{tb} is sub-zone {i} of {ta}, at depth {}",
                            rb - ra
                        ));
                    }
                }
                if descendant {
                    if let Ok(Some(i)) = grid.sub_zone_index(ib, ia) {
                        lines.push(format!(
                            "{ta} is sub-zone {i} of {tb}, at depth {}",
                            ra - rb
                        ));
                    }
                }
            }
            for l in lines {
                writeln!(out, "{l}")?;
            }
        }
        Command::Sub { zone, index } => {
            let id = zone_arg(grid, zone)?;
            let depth = inv.options.depth.unwrap_or(1);
            let t = grid.text_id(id);
            match index {
                Some(i) => {
                    let s = grid
                        .sub_zone_at_index(id, depth, *i)
                        .map_err(|e| refusal(grid, e))?;
                    // A position of an order may hold the null zone, as DGGAL's own does.
                    writeln!(
                        out,
                        "sub-zone {i} of {t} at depth {depth}: {}",
                        text_or_no_cell(grid, s)
                    )?;
                }
                None => {
                    let subs = grid.sub_zones(id, depth).map_err(|e| refusal(grid, e))?;
                    writeln!(
                        out,
                        "{name} zone {t}: {} at depth {depth}",
                        match subs.len() {
                            1 => "1 sub-zone".to_string(),
                            n => format!("{n} sub-zones"),
                        }
                    )?;
                    for (i, s) in subs.into_iter().enumerate() {
                        writeln!(out, "{i:>7}  {}", text_or_no_cell(grid, s))?;
                    }
                }
            }
        }
        Command::Index { parent, sub } => {
            let (ip, is) = (zone_arg(grid, parent)?, zone_arg(grid, sub)?);
            let (tp, ts) = (grid.text_id(ip), grid.text_id(is));
            match grid.sub_zone_index(ip, is).map_err(|e| refusal(grid, e))? {
                Some(i) => {
                    let depth = grid.resolution(is).saturating_sub(grid.resolution(ip));
                    writeln!(out, "{ts} is sub-zone {i} of {tp}, at depth {depth}")?;
                }
                // No more is said than the library knows: along the two broken seams of the
                // aperture-7 grids an order may hold a zone to which it gives no index.
                None => writeln!(out, "{ts} has no index among the sub-zones of {tp}")?,
            }
        }
        _ => unreachable!("answered by `answer`"),
    }
    Ok(())
}
