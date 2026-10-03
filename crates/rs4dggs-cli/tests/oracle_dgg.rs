//! The tool against DGGAL's own `dgg`, on the same points and zones, for the six grids.
//!
//! The library's suites already compare every value in process; these tests check what the tool
//! adds: its reading of arguments, its choice of what to print and its formatting. `dgg` is found
//! under `DGGAL_SITE_PACKAGES`, and its absence fails the tests.

use std::io::Cursor;
use std::path::PathBuf;
use std::process::Command;

use rs4dggs::{AnyGrid, ZoneId, get_grid};

/// The tool, run in process: exit code, standard output, standard error.
fn tool(args: &[&str]) -> (u8, String, String) {
    let args: Vec<String> = args.iter().map(|a| a.to_string()).collect();
    let mut input = Cursor::new(Vec::new());
    let (mut out, mut err) = (Vec::new(), Vec::new());
    let code = rs4dggs_cli::run(&args, &mut input, &mut out, &mut err);
    (
        code,
        String::from_utf8(out).unwrap(),
        String::from_utf8(err).unwrap(),
    )
}

/// The tool's answer, which must be a success.
fn tool_ok(args: &[&str]) -> String {
    let (code, out, err) = tool(args);
    assert_eq!(code, 0, "{args:?}: {err}");
    out
}

/// The seconds after which a `dgg` still running is stopped, so that no call can hold a test for
/// ever. No call made here takes a second by itself; the bound is generous on purpose, and is
/// no measure of anything.
const DGG_SECONDS: &str = "600";

/// `dgg` run with the arguments, under the `timeout` of the coreutils: whether it answered (exit
/// code 0) or refused (exit code 1, as for an index at which it finds no zone), and its standard
/// output. Any other end fails the test, so that it is never taken for an answer: DGGAL ends the
/// process at some zones, which shows here as a signal, and a `dgg` stopped for running too long
/// shows as exit code 124.
fn dgg_raw(args: &[&str]) -> (bool, String) {
    let site = std::env::var("DGGAL_SITE_PACKAGES")
        .expect("DGGAL_SITE_PACKAGES must name a site-packages with dggal installed");
    let exe: PathBuf = [site.as_str(), "dggal", "bin", "dgg"].iter().collect();
    assert!(exe.is_file(), "{} does not exist", exe.display());
    let libs = format!("{site}/dggal/lib:{site}/ecrt/lib");
    let ld = match std::env::var("LD_LIBRARY_PATH") {
        Ok(old) if !old.is_empty() => format!("{libs}:{old}"),
        _ => libs,
    };
    let output = Command::new("timeout")
        .arg(DGG_SECONDS)
        .arg(&exe)
        .args(args)
        .env("LD_LIBRARY_PATH", ld)
        .output()
        .expect("timeout and dgg run");
    let code = output.status.code();
    assert!(
        matches!(code, Some(0 | 1)),
        "dgg {args:?} neither answered nor refused: {} (a signal is the engine ending the \
         process, and 124 the time allowed running out)",
        output.status
    );
    (
        code == Some(0),
        String::from_utf8(output.stdout).expect("dgg writes UTF-8"),
    )
}

/// `dgg`'s standard output, which must be a success.
fn dgg(args: &[&str]) -> String {
    let (ok, out) = dgg_raw(args);
    assert!(ok, "dgg {args:?} failed: {out}");
    out
}

/// A 64-bit linear congruential generator, so that both tools get the same sample everywhere.
struct Lcg(u64);

impl Lcg {
    fn unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// The points, as text with six decimal places: forty from the generator, both poles, and the
/// antimeridian on both sides.
fn sample_points() -> Vec<String> {
    let mut rng = Lcg(0x5EED_0FD6_65EE_D001);
    let mut points: Vec<String> = (0..40)
        .map(|_| {
            let lat = (2.0 * rng.unit() - 1.0).asin().to_degrees();
            let lon = 360.0 * rng.unit() - 180.0;
            format!("{lat:.6},{lon:.6}")
        })
        .collect();
    for fixed in ["90,0", "-90,0", "0,180", "0,-180"] {
        points.push(fixed.to_string());
    }
    points
}

/// The eC runtime's `cleanFinalDigits`, which trims what `%.nlf` leaves at the end of a number.
/// `number` is the C string with room for its terminator; a zero byte ends it.
fn clean_final_digits(number: &mut [u8], num_digits: i64) {
    let len = number.iter().position(|&b| b == 0).unwrap();
    let dot = number[..len].iter().position(|&b| b == b'.');
    let (mut last, mut check_for_1, mut check_for_9, mut first9) = (0usize, true, true, 0usize);
    for c in (0..len).rev() {
        let ch = number[c];
        if ch != b'0' {
            if let Some(d) = dot {
                let after_dot = c as i64 - d as i64;
                if ch == b'1' && after_dot >= num_digits - 1 && c == len - 1 && check_for_1 {
                    check_for_1 = false;
                } else if ch == b'9' && after_dot >= num_digits - 1 && c == len - 1 && check_for_9 {
                    first9 = c;
                } else {
                    last = last.max(c);
                    check_for_9 = false;
                    check_for_1 = false;
                }
            }
        }
        if ch == b'.' {
            if last == c {
                number[c] = 0;
            } else {
                number[last + 1] = 0;
                if first9 != 0 {
                    loop {
                        first9 -= 1;
                        if first9 == 0 {
                            break;
                        }
                        if first9 != c {
                            let ch = number[first9];
                            if ch == 0 || ch == b'.' {
                            } else if ch < b'9' {
                                number[first9] += 1;
                                let mut j = first9 + 1;
                                while j < dot.unwrap() {
                                    number[j] = b'0';
                                    j += 1;
                                }
                                number[j] = 0;
                                break;
                            }
                        }
                    }
                }
            }
            break;
        }
    }
}

/// A double as the eC runtime prints it (`Double_OnGetString`), which is how `dgg` writes a
/// coordinate: fixed notation with 17 decimals less one for each power of ten of its magnitude
/// above 0.01, and its final digits cleaned. That is dgg's printed precision, about 15
/// significant digits (14 decimals for |x| in [1, 10), 13 in [10, 100), 12 in [100, 180)), so
/// a comparison of two values through it holds them equal at that precision and no finer; the
/// comparison of the tool's own digits with the library's is a separate one, bit for bit.
/// Valid for |x| from 1e-20 to 1e20 (and zero), which covers every coordinate; outside that
/// range the eC prints `%.15e`, which is not modelled.
fn ec_double(x: f64) -> String {
    let mut digits = 17i64;
    let mut num = 0.01f64;
    let magnitude = x.abs();
    while digits > 0 && num < magnitude {
        digits -= 1;
        num *= 10.0;
    }
    let mut buffer = format!("{x:.*}", digits as usize).into_bytes();
    buffer.push(0);
    clean_final_digits(&mut buffer, digits);
    let end = buffer.iter().position(|&b| b == 0).unwrap();
    String::from_utf8(buffer[..end].to_vec()).unwrap()
}

/// A coordinate pair as `dgg` writes it.
fn ec_pair(p: (f64, f64)) -> (String, String) {
    (ec_double(p.0), ec_double(p.1))
}

/// The list items of a card, each with whether it carries the mark.
type Marked = Vec<(String, bool)>;

#[derive(Debug, Default)]
struct ToolCard {
    id: String,
    centroid: (f64, f64),
    parents: Marked,
    children: Marked,
    neighbours: Vec<String>,
    vertices: Vec<(f64, f64)>,
}

fn pair(s: &str) -> (f64, f64) {
    let (a, b) = s
        .split_once(',')
        .unwrap_or_else(|| panic!("no pair in {s:?}"));
    (a.trim().parse().unwrap(), b.trim().parse().unwrap())
}

/// The items after `n: ` in `a, b (mark), c`.
fn marked_list(s: &str, mark: &str) -> Marked {
    let s = s.trim();
    if s == "none" {
        return Vec::new();
    }
    let (_, items) = s.split_once(": ").unwrap_or(("", ""));
    let suffix = format!(" ({mark})");
    items
        .split(", ")
        .filter(|i| !i.is_empty())
        .map(|i| match i.strip_suffix(suffix.as_str()) {
            Some(bare) => (bare.to_string(), true),
            None => (i.to_string(), false),
        })
        .collect()
}

fn parse_tool_card(text: &str) -> ToolCard {
    let mut card = ToolCard::default();
    let mut lines = text.lines();
    card.id = lines
        .next()
        .and_then(|l| l.split_whitespace().nth(2))
        .expect("a card starts with the zone")
        .to_string();
    let mut in_vertices = false;
    for line in lines {
        if in_vertices {
            card.vertices.push(pair(line));
        } else if let Some(v) = line.strip_prefix("  centroid    ") {
            card.centroid = pair(v.split("   (").next().unwrap());
        } else if let Some(v) = line.strip_prefix("  parents     ") {
            card.parents = marked_list(v, "centroid parent");
        } else if let Some(v) = line.strip_prefix("  children    ") {
            card.children = marked_list(v, "centre");
        } else if let Some(v) = line.strip_prefix("  neighbours  ") {
            card.neighbours = marked_list(v, "").into_iter().map(|(n, _)| n).collect();
        } else if line.starts_with("  vertices    ") {
            in_vertices = true;
        }
    }
    card
}

#[derive(Debug, Default)]
struct DggCard {
    id: String,
    centroid: (String, String),
    parents: Marked,
    children: Marked,
    neighbours: Vec<String>,
    vertices: Vec<(String, String)>,
}

/// The two coordinates of a `dgg` line, as written.
fn raw_pair(s: &str) -> (String, String) {
    let (a, b) = s
        .split_once(',')
        .unwrap_or_else(|| panic!("no pair in {s:?}"));
    (a.trim().to_string(), b.trim().to_string())
}

fn parse_dgg_card(text: &str) -> DggCard {
    #[derive(PartialEq)]
    enum In {
        None,
        Parents,
        Children,
        Neighbours,
        Vertices,
    }
    let mut card = DggCard::default();
    let mut section = In::None;
    for line in text.lines() {
        if let Some(v) = line.strip_prefix("Textual Zone ID: ") {
            card.id = v.trim().to_string();
        } else if let Some(v) = line.strip_prefix("WGS84 Centroid (lat, lon): ") {
            card.centroid = raw_pair(v);
        } else if line.starts_with("Parent") {
            section = In::Parents;
        } else if line.starts_with("Children") {
            section = In::Children;
        } else if line.starts_with("Neighbors") {
            section = In::Neighbours;
        } else if line.starts_with("[EPSG:4326] Vertices") {
            section = In::Vertices;
        } else if line.starts_with("   ") {
            let item = line.trim();
            match section {
                In::Parents => card
                    .parents
                    .push(match item.strip_suffix(" (centroid child)") {
                        Some(bare) => (bare.to_string(), true),
                        None => (item.to_string(), false),
                    }),
                In::Children => card.children.push(match item.strip_suffix(" (centroid)") {
                    Some(bare) => (bare.to_string(), true),
                    None => (item.to_string(), false),
                }),
                In::Neighbours => card
                    .neighbours
                    .push(item.split_once("): ").expect("a direction").1.to_string()),
                In::Vertices => card.vertices.push(raw_pair(item)),
                In::None => {}
            }
        } else {
            section = In::None;
        }
    }
    card
}

/// `dgg`'s neighbours less the four kinds the library drops: the null zone, the zone itself,
/// repeats, and identifiers the engine cannot read back.
fn engine_neighbours_less_four(grid: &AnyGrid, own: &str, list: &[String]) -> Vec<String> {
    let mut kept: Vec<String> = Vec::new();
    for n in list {
        let readable = matches!(grid.zone_from_text(n), Ok(id) if id != ZoneId::NULL);
        if n != "(null)" && n != own && readable && !kept.contains(n) {
            kept.push(n.clone());
        }
    }
    kept
}

/// Whether the grid reads `text` as a zone. `dgg` prints the null zone as `(null)`, and along the
/// two broken seams of the aperture-7 grids some identifiers that the engine itself cannot read.
fn readable(grid: &AnyGrid, text: &str) -> bool {
    matches!(grid.zone_from_text(text), Ok(id) if id != ZoneId::NULL)
}

/// `dgg`'s parents or children less the three kinds the library never hands out: the null zone,
/// an identifier the engine cannot read back, and a repeat. The marks stay with what is kept.
fn engine_list_less_three(grid: &AnyGrid, list: &Marked) -> Marked {
    let mut kept: Marked = Vec::new();
    for (z, mark) in list {
        if readable(grid, z) && !kept.iter().any(|(k, _)| k == z) {
            kept.push((z.clone(), *mark));
        }
    }
    kept
}

/// The parents and children of a card of an aperture-7 grid against `dgg`'s: the same zones in
/// the same order, less the three kinds dropped; the same child marked; and the parent marked
/// that the library names as the centroid parent, where `dgg` marks none.
fn compare_hierarchy_on_aperture_7(
    grid: &AnyGrid,
    ours: &ToolCard,
    theirs: &DggCard,
    context: &str,
    counts: &mut Counts,
) {
    let id = grid.zone_from_text(&ours.id).unwrap();
    let parents = engine_list_less_three(grid, &theirs.parents);
    let children = engine_list_less_three(grid, &theirs.children);
    let dropped = theirs.parents.len() - parents.len() + theirs.children.len() - children.len();
    counts.hierarchy_cards += 1;
    counts.hierarchy_dropped += dropped;
    counts.hierarchy_cards_with_drops += usize::from(dropped > 0);
    // DGGAL answers no centroid parent for any zone of these grids, so `dgg` marks no parent.
    // A `dgg` that marked one would fail here, and so be noticed.
    assert!(
        theirs.parents.iter().all(|p| !p.1),
        "dgg marks a parent: {context}"
    );
    // The tool marks the library's: the first parent that is itself a centroid child.
    let centroid_parent = grid.centroid_parent(id).map(|z| grid.text_id(z));
    let expected: Marked = parents
        .iter()
        .map(|(p, _)| (p.clone(), Some(p) == centroid_parent.as_ref()))
        .collect();
    assert_eq!(ours.parents, expected, "parents: {context}");
    counts.marked_parents += ours.parents.iter().filter(|p| p.1).count();
    assert_eq!(ours.children, children, "children: {context}");
}

/// Identifiers on one of the two broken seams of the aperture-7 grids, which their text alone
/// reaches, and at which `dgg`'s lists hold what the library drops: an identifier the engine
/// cannot read and a repeat among the children of the first; a repeated parent, and those two
/// kinds among the children, at the second; and at the third no parent at all, above
/// resolution 0, and the null zone among the children.
const SEAM_CARDS: [&str; 3] = [
    "0000000000000001644",
    "0000000000000000136",
    "00000000000000001311",
];

/// The leading integer of `s` after ` by `.
fn after_by(s: &str) -> i64 {
    let rest = s.split(" by ").nth(1).expect("a difference by");
    rest.split_whitespace().next().unwrap().parse().unwrap()
}

/// The signed difference of resolutions, second minus first, from either tool's line.
fn resolution_difference(text: &str) -> i64 {
    for line in text.lines() {
        if line.contains("coarser than") {
            return after_by(line);
        }
        if line.contains("finer than") {
            return -after_by(line);
        }
        if line.contains("same resolution") || line.contains("same refinement level") {
            return 0;
        }
    }
    panic!("no resolution line in {text}");
}

/// The six relations, each true where the tool's line, or `dgg`'s positive line, is present.
const RELATIONS: [(&str, &str, &str); 6] = [
    (
        "immediate child",
        " is an immediate child of ",
        "Zone A is an immediate child of",
    ),
    (
        "immediate parent",
        " is an immediate parent of ",
        "Zone A is an immediate parent of",
    ),
    (
        "descendant",
        " is a descendant of ",
        "Zone A is a descendant of",
    ),
    (
        "ancestor",
        " is an ancestor of ",
        "Zone A is an ancestor of",
    ),
    ("neighbours", " are neighbours", "These zones are neighbors"),
    ("siblings", " are siblings", "These zones are siblings"),
];

#[derive(Default)]
struct Counts {
    zones: usize,
    nulls: usize,
    relations: usize,
    relation_lines: usize,
    neighbour_lists: usize,
    dropped_neighbours: usize,
    hierarchy_cards: usize,
    hierarchy_cards_with_drops: usize,
    hierarchy_dropped: usize,
    marked_parents: usize,
    no_cell_entries: usize,
    unreadable_entries: usize,
    firsts: usize,
    firsts_not_asked: usize,
    depth_0_firsts_differing: usize,
    coarse_zones: usize,
    sub_lists: usize,
    sub_corrected: usize,
    index_not_found: usize,
    indices: usize,
}

fn compare_rel(tool_name: &str, dgg_name: &str, a: &str, b: &str, counts: &mut Counts) {
    let t = tool_ok(&[tool_name, "rel", a, b]);
    let d = dgg(&[dgg_name, "rel", a, b]);
    let context = format!("{tool_name} rel {a} {b}\ntool:\n{t}\ndgg:\n{d}");
    counts.relations += 1;
    if d.contains("These two zones are the same") {
        assert!(t.contains("are the same zone"), "{context}");
        return;
    }
    assert_eq!(
        resolution_difference(&t),
        resolution_difference(&d),
        "{context}"
    );
    counts.relation_lines += 1;
    for (name, ours, theirs) in RELATIONS {
        // The tool's neighbours line must be the positive one, not "are not neighbours".
        let ours_says = t
            .lines()
            .any(|l| l.contains(ours) && !l.contains("are not"));
        let theirs_says = d.lines().any(|l| l.contains(theirs));
        assert_eq!(ours_says, theirs_says, "{name}: {context}");
        counts.relation_lines += 1;
    }
}

/// The quoted strings of a JSON array of strings.
fn json_strings(text: &str) -> Vec<String> {
    text.split('"')
        .skip(1)
        .step_by(2)
        .map(str::to_string)
        .collect()
}

fn run_grid(tool_name: &str, dgg_name: &str, aperture3: bool) {
    let grid = get_grid(&tool_name.to_ascii_uppercase()).expect("a grid");
    let resolutions: &[u8] = if aperture3 {
        &[0, 7, 15, 25]
    } else {
        &[0, 5, 10, 15]
    };
    let points = sample_points();
    let mut counts = Counts::default();
    for res in resolutions {
        let res_text = res.to_string();
        // The zone of every point, from the tool, kept for the relations to a zone elsewhere.
        let mut others: Vec<String> = Vec::new();
        for p in &points {
            let card_text = tool_ok(&[tool_name, "zone", p, &res_text]);
            let ours = parse_tool_card(&card_text);
            if card_text.contains("no zone contains") {
                counts.nulls += 1;
                continue;
            }
            others.push(ours.id.clone());
        }
        for (k, p) in points.iter().enumerate() {
            let card_text = tool_ok(&[tool_name, "zone", p, &res_text]);
            let ours = parse_tool_card(&card_text);
            if card_text.contains("no zone contains") {
                continue;
            }
            let theirs = parse_dgg_card(&dgg(&[dgg_name, "zone", p, &res_text]));
            let context =
                format!("{tool_name} zone {p} {res}\ntool:\n{card_text}\ndgg: {theirs:?}");
            counts.zones += 1;
            assert_eq!(ours.id, theirs.id, "{context}");
            assert_eq!(
                ec_pair(ours.centroid),
                theirs.centroid,
                "centroid: {context}"
            );
            // `dgg` prints the engine's order; the library turns an aperture-3 ring about its
            // first vertex, save at the two polar pentagons (roots 10 and 11), so that every
            // ring runs anticlockwise.
            let id = grid.zone_from_text(&ours.id).unwrap();
            let mut dgg_ring = theirs.vertices.clone();
            if aperture3 && (id.0 >> 53) & 0xF < 10 {
                dgg_ring[1..].reverse();
            }
            assert_eq!(ours.vertices.len(), dgg_ring.len(), "{context}");
            for (o, t) in ours.vertices.iter().zip(&dgg_ring) {
                assert_eq!(ec_pair(*o), *t, "vertex: {context}");
            }
            // The digits the tool prints read back as the library's own doubles, bit for bit.
            let bits = |p: (f64, f64)| (p.0.to_bits(), p.1.to_bits());
            let centroid = grid.centroid(id);
            assert_eq!(
                bits(ours.centroid),
                bits((centroid.lat, centroid.lon)),
                "centroid bits: {context}"
            );
            let vertices = grid.vertices(id);
            assert_eq!(ours.vertices.len(), vertices.len(), "{context}");
            for (o, v) in ours.vertices.iter().zip(&vertices) {
                assert_eq!(bits(*o), bits((v.lat, v.lon)), "vertex bits: {context}");
            }
            let csv = tool_ok(&[tool_name, "info", &ours.id, "-f", "csv"]);
            let row: Vec<&str> = csv.lines().nth(1).unwrap().split(',').collect();
            assert_eq!(
                (
                    row[3].parse::<f64>().unwrap().to_bits(),
                    row[4].parse::<f64>().unwrap().to_bits()
                ),
                bits((centroid.lat, centroid.lon)),
                "csv centroid bits: {context}: {csv}"
            );
            let all = engine_neighbours_less_four(&grid, &theirs.id, &theirs.neighbours);
            counts.dropped_neighbours += theirs.neighbours.len() - all.len();
            counts.neighbour_lists += 1;
            assert_eq!(ours.neighbours, all, "neighbours: {context}");
            if aperture3 {
                assert_eq!(ours.parents, theirs.parents, "parents: {context}");
                assert_eq!(ours.children, theirs.children, "children: {context}");
            } else {
                compare_hierarchy_on_aperture_7(&grid, &ours, &theirs, &context, &mut counts);
            }

            // The relations of this zone.
            let mut pairs: Vec<(String, String)> = Vec::new();
            if let Some(n) = ours.neighbours.first() {
                pairs.push((ours.id.clone(), n.clone()));
            }
            if let Some((child, _)) = ours.children.first() {
                pairs.push((ours.id.clone(), child.clone()));
                pairs.push((child.clone(), ours.id.clone()));
                let cid = grid.zone_from_text(child).unwrap();
                if let Some(g) = grid.children(cid).first() {
                    pairs.push((ours.id.clone(), grid.text_id(*g)));
                }
            }
            // The last child, which on aperture 7 lies under a neighbouring zone as well; and
            // the first and the last, which are siblings.
            if let (Some((first, _)), Some((last, _))) =
                (ours.children.first(), ours.children.last())
            {
                pairs.push((ours.id.clone(), last.clone()));
                pairs.push((first.clone(), last.clone()));
            }
            if !others.is_empty() {
                pairs.push((ours.id.clone(), others[(k + 1) % others.len()].clone()));
            }
            for (a, b) in pairs {
                compare_rel(tool_name, dgg_name, &a, &b, &mut counts);
            }
        }
    }

    if aperture3 {
        sub_and_index(&grid, tool_name, dgg_name, &points, &mut counts);
    } else {
        for text in SEAM_CARDS {
            let ours = parse_tool_card(&tool_ok(&[tool_name, "info", text]));
            let theirs = parse_dgg_card(&dgg(&[dgg_name, "info", text]));
            let context = format!("{tool_name} info {text}\ntool: {ours:?}\ndgg: {theirs:?}");
            assert_eq!(ours.id, theirs.id, "{context}");
            compare_hierarchy_on_aperture_7(&grid, &ours, &theirs, &context, &mut counts);
        }
        sub_and_index_on_aperture_7(&grid, tool_name, dgg_name, &points, &mut counts);
    }

    println!(
        "{tool_name}: {} zone cards ({} points answered as the null zone, left out), \
         {} neighbour lists ({} entries of the engine's dropped), {} relations \
         ({} lines compared)",
        counts.zones,
        counts.nulls,
        counts.neighbour_lists,
        counts.dropped_neighbours,
        counts.relations,
        counts.relation_lines,
    );
    assert!(counts.zones >= 150, "zone cards: {}", counts.zones);
    assert!(counts.relations >= 150, "relations: {}", counts.relations);
    if aperture3 {
        println!(
            "{tool_name}: {} sub-zone lists (three depths of each of {} zones at resolutions 0 \
             and 1, 30 others, 6 at the south pole and one rhombus-edge zone; {} where dgg's \
             holds a repeat or a non-descendant, with {} entries of those that `dgg index` does \
             not find), {} indices",
            counts.sub_lists,
            counts.coarse_zones,
            counts.sub_corrected,
            counts.index_not_found,
            counts.indices,
        );
        assert!(
            counts.sub_corrected >= 1,
            "no rhombus-edge order was corrected"
        );
        assert!(
            counts.sub_lists >= 20,
            "sub-zone lists: {}",
            counts.sub_lists
        );
    } else {
        println!(
            "{tool_name}: {} cards whose parents and children are dgg's less the three kinds \
             dropped ({} of them with a drop, {} entries dropped in all; {} parents marked, \
             where dgg marks none), {} sub-zone lists (three depths of each of {} zones at \
             resolutions 0 and 1, 30 others, and one on a seam at depth 5; {} positions with no \
             cell, {} of them an identifier of dgg's that cannot be read), {} first sub-zones \
             ({} not asked of dgg, which the engine would end), {} of {} first sub-zones at \
             depth 0 where dgg by index is not the zone itself, {} indices",
            counts.hierarchy_cards,
            counts.hierarchy_cards_with_drops,
            counts.hierarchy_dropped,
            counts.marked_parents,
            counts.sub_lists,
            counts.coarse_zones,
            counts.no_cell_entries,
            counts.unreadable_entries,
            counts.firsts,
            counts.firsts_not_asked,
            counts.depth_0_firsts_differing,
            counts.coarse_zones,
            counts.indices,
        );
        // Each as measured. Of the cards, one of the sample has a drop (the zone at the north
        // pole at resolution 15, on a seam, among whose children dgg names two twice) and the
        // three seam cards have 25 between them: 9, 5 and 11.
        assert_eq!(counts.hierarchy_cards, counts.zones + SEAM_CARDS.len());
        assert_eq!(
            (counts.hierarchy_cards_with_drops, counts.hierarchy_dropped),
            (4, 27),
            "{tool_name}"
        );
        assert!(
            counts.marked_parents >= 30,
            "parents marked: {}",
            counts.marked_parents
        );
        assert_eq!(counts.coarse_zones, 84, "{tool_name}");
        assert_eq!(
            (counts.sub_lists, counts.indices),
            (283, 283),
            "{tool_name}"
        );
        // The 99 positions are those of the seam zone's order at depth 5, each the null zone
        // in dgg's own list.
        assert_eq!(
            (counts.no_cell_entries, counts.unreadable_entries),
            (99, 0),
            "{tool_name}"
        );
        // The one first sub-zone not asked is that of `110` at depth 2.
        assert_eq!(
            (counts.firsts, counts.firsts_not_asked),
            (282, 1),
            "{tool_name}"
        );
        assert_eq!(counts.depth_0_firsts_differing, 84, "{tool_name}");
    }
}

/// What the tool prints in a sub-zone order at a position that holds no zone.
const NO_CELL: &str = "(no cell)";

/// Whether `dgg` is not to be asked for sub-zone 0 of `zone` at `depth`: at the pentagon of base
/// cell 11, the southern polar one, at an odd resolution and at depth 2, DGGAL's first sub-zone
/// ends the process. The list of that same order is answered, and is asked.
fn the_first_sub_zone_ends_dgg(zone: &str, depth: u8) -> bool {
    depth == 2
        && zone.len() % 2 == 1
        && zone.starts_with("11")
        && zone[2..].bytes().all(|b| b == b'0')
}

/// The entries of the tool's list of sub-zones, in order.
fn tool_sub_zones(tool_name: &str, zone: &str, depth: u8) -> Vec<String> {
    tool_ok(&[tool_name, "sub", zone, "-depth", &depth.to_string()])
        .lines()
        .skip(1)
        .map(|l| l.trim_start().split_once("  ").unwrap().1.to_string())
        .collect()
}

/// The tool's sub-zone at an index, which must be the entry given.
fn assert_tool_sub_zone_at(tool_name: &str, zone: &str, depth: u8, index: usize, entry: &str) {
    let args = [
        tool_name,
        "sub",
        zone,
        &index.to_string(),
        "-depth",
        &depth.to_string(),
    ];
    let expected = format!("sub-zone {index} of {zone} at depth {depth}: {entry}");
    assert_eq!(tool_ok(&args).trim_end(), expected, "{args:?}");
}

/// `dgg`'s sub-zone at an index: the zone, or `None` where it says that the index is invalid,
/// as it does where the zone at the index is the null zone.
fn dgg_sub_zone_at(dgg_name: &str, zone: &str, depth: u8, index: usize) -> Option<String> {
    let args = [
        dgg_name,
        "sub",
        zone,
        "-depth",
        &depth.to_string(),
        &index.to_string(),
    ];
    let (answered, out) = dgg_raw(&args);
    if !answered {
        assert!(out.contains("Invalid zone index"), "dgg {args:?}: {out}");
        return None;
    }
    let zones = json_strings(&out);
    assert_eq!(zones.len(), 1, "dgg {args:?}: {out}");
    zones.into_iter().next()
}

/// The sub-zone orders of an aperture-7 grid against `dgg`'s: the list entry by entry, the first
/// sub-zone by its index, and the index of an entry half-way along.
fn sub_and_index_on_aperture_7(
    grid: &AnyGrid,
    tool_name: &str,
    dgg_name: &str,
    points: &[String],
    counts: &mut Counts,
) {
    // Thirty zones at resolutions 0 to 16, each at one of the depths 1 to 3.
    let mut rng = Lcg(0xC0FF_EE5E_ED12_3400);
    let mut cases: Vec<(ZoneId, u8)> = Vec::new();
    for k in 0..30 {
        let res = (rng.unit() * 17.0) as u8;
        let p = &points[(rng.unit() * points.len() as f64) as usize];
        let (lat, lon) = p.split_once(',').unwrap();
        let id = grid
            .zone_from_geo(lat.parse().unwrap(), lon.parse().unwrap(), res)
            .unwrap();
        assert_ne!(id, ZoneId::NULL, "{tool_name}: {p} at resolution {res}");
        cases.push((id, (k % 3 + 1) as u8));
    }
    // Every zone of resolutions 0 and 1, the twelve pentagons of each among them, at depths 1 to
    // 3; and at depth 0, where the order is the zone alone.
    let mut coarse: Vec<ZoneId> = (0..12)
        .map(|base| grid.zone_from_text(&format!("{base:02}")).unwrap())
        .collect();
    for id in coarse.clone() {
        coarse.extend(grid.sub_zones(id, 1).unwrap());
    }
    coarse.sort_by_key(|z| z.0);
    coarse.dedup();
    counts.coarse_zones = coarse.len();
    for id in coarse {
        for depth in 0..=3 {
            cases.push((id, depth));
        }
    }
    // A zone on one of the two broken seams, whose order at depth 5 has positions with no zone.
    cases.push((grid.zone_from_text("010004000400").unwrap(), 5));

    for (id, depth) in cases {
        let text = grid.text_id(id);
        let depth_text = depth.to_string();
        let context = format!("{tool_name} sub {text} depth {depth}");

        // The list: dgg's own, with the tool's words wherever dgg's entry is no zone.
        let ours = tool_sub_zones(tool_name, &text, depth);
        let theirs = json_strings(&dgg(&[dgg_name, "sub", &text, "-depth", &depth_text]));
        assert_eq!(ours.len(), theirs.len(), "{context}");
        for (i, (o, t)) in ours.iter().zip(&theirs).enumerate() {
            if readable(grid, t) {
                assert_eq!(o, t, "{context}: position {i}");
            } else {
                assert_eq!(o, NO_CELL, "{context}: position {i}, where dgg has {t}");
                counts.no_cell_entries += 1;
                counts.unreadable_entries += usize::from(t != "(null)");
            }
        }

        // The first sub-zone, by its index: the first entry of the list.
        assert_tool_sub_zone_at(tool_name, &text, depth, 0, &ours[0]);
        if depth == 0 {
            // The order is the zone alone, and the tool answers it; dgg by index answers
            // another zone.
            assert_eq!(theirs, [text.as_str()], "{context}");
            let by_index = dgg_sub_zone_at(dgg_name, &text, depth, 0);
            counts.depth_0_firsts_differing += usize::from(by_index.as_ref() != Some(&text));
            continue;
        }
        counts.sub_lists += 1;
        if the_first_sub_zone_ends_dgg(&text, depth) {
            counts.firsts_not_asked += 1;
        } else {
            let by_index = dgg_sub_zone_at(dgg_name, &text, depth, 0);
            let expected = readable(grid, &theirs[0]).then(|| theirs[0].clone());
            assert_eq!(by_index, expected, "{context}: dgg's first sub-zone");
            counts.firsts += 1;
        }

        // The index of the entry half-way along, or of the first after it that is a zone and
        // that the order holds once.
        let middle = (ours.len() / 2..ours.len())
            .find(|&i| ours[i] != NO_CELL && ours.iter().filter(|o| **o == ours[i]).count() == 1)
            .unwrap_or_else(|| panic!("{context}: no entry to ask the index of"));
        let sub = &ours[middle];
        let index_text = tool_ok(&[tool_name, "index", &text, sub]);
        let expected = format!("{sub} is sub-zone {middle} of {text}, at depth {depth}");
        assert_eq!(index_text.trim(), expected, "{context}");
        let d = dgg(&[dgg_name, "index", &text, sub]);
        let expected = format!("{sub} is at index {middle} of {text} at depth {depth}");
        assert_eq!(d.lines().last().unwrap().trim(), expected, "{context}: {d}");
        counts.indices += 1;
    }

    // Along the seams the engine's index of a sub-zone may be a position at which its own order
    // holds another zone. The tool states an index only where the order holds the sub-zone at
    // it, and here says that there is none: the order holds the sub-zone at position 3, and
    // `dgg` says 6.
    let (zone, sub) = ("00000000000000000", "000000000000000001");
    let list = json_strings(&dgg(&[dgg_name, "sub", zone, "-depth", "1"]));
    assert_eq!(tool_sub_zones(tool_name, zone, 1), list, "{tool_name}");
    assert_eq!(
        (list[3].as_str(), list[6].as_str()),
        (sub, "000000000000000000"),
        "{tool_name}"
    );
    let d = dgg(&[dgg_name, "index", zone, sub]);
    assert_eq!(
        d.lines().last().unwrap().trim(),
        format!("{sub} is at index 6 of {zone} at depth 1"),
        "{tool_name}"
    );
    assert_eq!(
        tool_ok(&[tool_name, "index", zone, sub]),
        format!("{sub} has no index among the sub-zones of {zone}\n"),
        "{tool_name}"
    );
}

fn sub_and_index(
    grid: &AnyGrid,
    tool_name: &str,
    dgg_name: &str,
    points: &[String],
    counts: &mut Counts,
) {
    // Thirty zones at resolutions 0 to 20, and every zone of resolutions 0 and 1, where the
    // rhombus edges are most of the grid.
    let mut rng = Lcg(0xC0FF_EE5E_ED12_3400);
    let mut cases: Vec<(ZoneId, u8)> = Vec::new();
    for k in 0..30 {
        let res = (rng.unit() * 21.0) as u8;
        let p = &points[(rng.unit() * points.len() as f64) as usize];
        let (lat, lon) = p.split_once(',').unwrap();
        let id = grid
            .zone_from_geo(lat.parse().unwrap(), lon.parse().unwrap(), res)
            .unwrap();
        cases.push((id, (k % 3 + 1) as u8));
    }
    // The southern polar pentagon at resolutions where the engine's order is known to fail.
    for res in [11, 25] {
        let pole = grid.zone_from_geo(-90.0, 0.0, res).unwrap();
        for depth in [3, 5, 7] {
            cases.push((pole, depth));
        }
    }
    // A rhombus-edge zone whose engine order is unsound: entry 4 is no descendant.
    cases.push((grid.zone_from_text("K5-E6A9-A").unwrap(), 2));
    let mut coarse: Vec<ZoneId> = Vec::new();
    for _ in 0..2000 {
        let (lat, lon) = (
            (2.0 * rng.unit() - 1.0).asin().to_degrees(),
            360.0 * rng.unit() - 180.0,
        );
        let id = grid.zone_from_geo(lat, lon, 0).unwrap();
        if !coarse.contains(&id) {
            coarse.push(id);
        }
    }
    for id in coarse.clone() {
        coarse.extend(grid.sub_zones(id, 1).unwrap());
    }
    coarse.sort_by_key(|z| z.0);
    coarse.dedup();
    counts.coarse_zones = coarse.len();
    for id in coarse {
        for depth in 1..=3 {
            cases.push((id, depth));
        }
    }
    for (id, depth) in cases {
        if id == ZoneId::NULL {
            continue;
        }
        let text = grid.text_id(id);
        let depth_text = depth.to_string();
        let context = format!("{tool_name} sub {text} depth {depth}");

        let ours: Vec<String> = tool_ok(&[tool_name, "sub", &text, "-depth", &depth_text])
            .lines()
            .skip(1)
            .map(|l| l.split_whitespace().nth(1).unwrap().to_string())
            .collect();
        let library: Vec<String> = grid
            .sub_zones(id, depth)
            .unwrap()
            .into_iter()
            .map(|z| grid.text_id(z))
            .collect();
        assert_eq!(ours, library, "{context}");
        let theirs = json_strings(&dgg(&[dgg_name, "sub", &text, "-depth", &depth_text]));
        counts.sub_lists += 1;
        // Whether each entry of dgg's list is a descendant of the zone, and whether it occurs once.
        let descendant =
            |z: &String| matches!(grid.zone_from_text(z), Ok(zid) if grid.is_ancestor_of(id, zid));
        let entry_sound: Vec<bool> = theirs
            .iter()
            .map(|z| descendant(z) && theirs.iter().filter(|t| *t == z).count() == 1)
            .collect();
        let sound = entry_sound.iter().all(|s| *s);
        if sound {
            assert_eq!(ours, theirs, "{context}");
        } else {
            counts.sub_corrected += 1;
            // The tool's list: no repeat, every entry a descendant, dgg's length, and dgg's own
            // entry wherever that is sound.
            for (i, z) in ours.iter().enumerate() {
                assert!(descendant(z), "{context}: entry {i} {z} is no descendant");
                assert!(
                    ours.iter().filter(|o| *o == z).count() == 1,
                    "{context}: {z} repeats"
                );
            }
            assert_eq!(ours.len(), theirs.len(), "{context}");
            for (i, sound_entry) in entry_sound.iter().enumerate() {
                if *sound_entry {
                    assert_eq!(ours[i], theirs[i], "{context}: position {i}");
                }
            }
            // Where `dgg index` does not find an entry the tool places, say so.
            for z in &ours {
                if dgg_raw(&[dgg_name, "index", &text, z])
                    .1
                    .contains("not found")
                {
                    counts.index_not_found += 1;
                }
            }
        }

        // The index of the sub-zone half-way along.
        let middle = ours.len() / 2;
        let sub = &ours[middle];
        let index_text = tool_ok(&[tool_name, "index", &text, sub]);
        let expected = format!("{sub} is sub-zone {middle} of {text}, at depth {depth}");
        assert_eq!(index_text.trim(), expected, "{context}");
        counts.indices += 1;
        if sound {
            let d = dgg(&[dgg_name, "index", &text, sub]);
            let expected = format!("{sub} is at index {middle} of {text} at depth {depth}");
            assert_eq!(d.lines().last().unwrap().trim(), expected, "{context}: {d}");
        }
    }
}

/// `dgg`'s sub-zone by its index against its own list, and the tool against that list, at every
/// index of the orders of one zone at depths 1 to 3: how many indices were compared, at how
/// many of each depth `dgg` by index is not the entry of its list, and how many were not asked.
fn by_index_at(tool_name: &str, dgg_name: &str, zone: &str) -> (usize, [usize; 3], usize) {
    let (mut compared, mut differing, mut not_asked) = (0, [0; 3], 0);
    for depth in 1..=3u8 {
        let list = json_strings(&dgg(&[dgg_name, "sub", zone, "-depth", &depth.to_string()]));
        for (i, entry) in list.iter().enumerate() {
            // The tool answers the entry of the order at the index, at every index.
            assert_tool_sub_zone_at(tool_name, zone, depth, i, entry);
            if i == 0 && the_first_sub_zone_ends_dgg(zone, depth) {
                not_asked += 1;
                continue;
            }
            compared += 1;
            if dgg_sub_zone_at(dgg_name, zone, depth, i).as_ref() != Some(entry) {
                differing[usize::from(depth) - 1] += 1;
            }
        }
    }
    (compared, differing, not_asked)
}

/// DGGAL's search for a sub-zone by its index is at fault at the pentagons: there `dgg` by index
/// answers, at some indices, a zone that is not the entry of its own list, where the tool
/// answers the entry. Measured here at the twelve pentagons of resolutions 0 and 1 and at three
/// hexagons, at depths 1 to 3 and at every index, each zone in a thread of its own.
fn sub_zone_by_index_against_dgg(tool_name: &str, dgg_name: &str) {
    // Each zone with the indices at which `dgg` by index is not the entry of its own list, at
    // depths 1, 2 and 3, exactly as measured, so that an engine corrected, or at fault
    // elsewhere, is noticed. The fault lies at the odd depths, and at depth 1 only at a
    // pentagon of an even resolution; the five pentagons of base cells 6 to 10 have about half
    // as many faulty indices as the other seven. At depth 2, and at a hexagon, there is none.
    let mut zones: Vec<(String, [usize; 3])> = Vec::new();
    // The pentagons are the base cells and, below them, their centre children: every digit 0.
    for level in 0..=1usize {
        for base in 0..12 {
            let expected = match ((6..=10).contains(&base), level) {
                (false, 0) => [2, 0, 115],
                (true, 0) => [1, 0, 54],
                (false, _) => [0, 0, 97],
                (true, _) => [0, 0, 53],
            };
            zones.push((format!("{base:02}{}", "0".repeat(level)), expected));
        }
    }
    let pentagons = zones.len();
    for hexagon in ["006", "0064", "00641"] {
        zones.push((hexagon.to_string(), [0, 0, 0]));
    }
    let results: Vec<(usize, [usize; 3], usize)> = std::thread::scope(|scope| {
        let threads: Vec<_> = zones
            .iter()
            .map(|(zone, _)| scope.spawn(move || by_index_at(tool_name, dgg_name, zone)))
            .collect();
        threads.into_iter().map(|t| t.join().unwrap()).collect()
    });

    let (mut compared, mut differing, mut not_asked) = (0, 0, 0);
    for (k, ((zone, expected), (c, d, n))) in zones.iter().zip(&results).enumerate() {
        compared += c;
        differing += d.iter().sum::<usize>();
        not_asked += n;
        assert_eq!(d, expected, "{tool_name} {zone}: indices where dgg differs");
        // Whatever the table above becomes, every pentagon shows the fault at an odd depth.
        if k < pentagons {
            assert!(d[0] + d[2] >= 1, "{tool_name} {zone}: no index differs");
        }
    }
    // One of them by name: sub-zone 9 of `00` at depth 1 is `006`, and `dgg` by index says `030`,
    // which is not in its list at all.
    let list = json_strings(&dgg(&[dgg_name, "sub", "00", "-depth", "1"]));
    assert_eq!(list[9], "006", "{tool_name}");
    assert!(list.iter().all(|z| z != "030"), "{tool_name}: {list:?}");
    assert_eq!(
        dgg_sub_zone_at(dgg_name, "00", 1, 9).as_deref(),
        Some("030"),
        "{tool_name}"
    );
    assert_tool_sub_zone_at(tool_name, "00", 1, 9, "006");

    println!(
        "{tool_name}: {compared} indices of {} orders asked of dgg by index and of the tool \
         ({pentagons} pentagons and {} hexagons at depths 1 to 3); at {differing} of them, all \
         at pentagons at depths 1 and 3, dgg by index is not the entry of its own list, which \
         the tool answers at every one; {not_asked} not asked of dgg, which the engine would end",
        3 * zones.len(),
        zones.len() - pentagons,
    );
    // 373 indices at each pentagon, 447 at each hexagon, less the one not asked: the first
    // sub-zone of `110` at depth 2.
    assert_eq!(
        (compared, not_asked),
        (24 * 373 + 3 * 447 - 1, 1),
        "{tool_name}"
    );
    assert_eq!(
        differing,
        7 * (2 + 115) + 5 * (1 + 54) + 7 * 97 + 5 * 53,
        "{tool_name}"
    );
}

#[test]
fn igeo7_sub_zone_by_index_against_dgg() {
    sub_zone_by_index_against_dgg("igeo7", "isea7h_z7");
}

#[test]
fn ivea7h_sub_zone_by_index_against_dgg() {
    sub_zone_by_index_against_dgg("ivea7h", "ivea7h_z7");
}

#[test]
fn rtea7h_sub_zone_by_index_against_dgg() {
    sub_zone_by_index_against_dgg("rtea7h", "rtea7h_z7");
}

#[test]
fn igeo7_against_dgg() {
    run_grid("igeo7", "isea7h_z7", false);
}

#[test]
fn ivea7h_against_dgg() {
    run_grid("ivea7h", "ivea7h_z7", false);
}

#[test]
fn rtea7h_against_dgg() {
    run_grid("rtea7h", "rtea7h_z7", false);
}

#[test]
fn isea3h_against_dgg() {
    run_grid("isea3h", "isea3h", true);
}

#[test]
fn ivea3h_against_dgg() {
    run_grid("ivea3h", "ivea3h", true);
}

#[test]
fn rtea3h_against_dgg() {
    run_grid("rtea3h", "rtea3h", true);
}
