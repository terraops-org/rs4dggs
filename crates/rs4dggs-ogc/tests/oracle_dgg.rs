//! The encodings against DGGAL's own `dgg`. `dgg` is found under `DGGAL_SITE_PACKAGES`, and its
//! absence fails the tests.

use std::path::PathBuf;
use std::process::Command;

/// The seconds after which a `dgg` still running is stopped, so that no call can hold a test for
/// ever. The bound is generous on purpose, and is no measure of anything.
const DGG_SECONDS: &str = "600";

/// `dgg` run with the arguments, under the `timeout` of the coreutils: its standard output.
///
/// Exit codes 0 and 1 are both accepted, since `dgg togeo` ends with 1 even when it has written
/// its answer, so that an answer is judged by what was written and never by the code. Any other
/// end fails the test: DGGAL ends the process at some zones, which shows here as a signal, and a
/// `dgg` stopped for running too long shows as exit code 124.
fn dgg(args: &[&str]) -> String {
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
    assert!(
        matches!(output.status.code(), Some(0 | 1)),
        "dgg {args:?} ended neither with 0 nor with 1: {} (a signal is the engine ending the \
         process, and 124 the time allowed running out)",
        output.status
    );
    String::from_utf8(output.stdout).expect("dgg writes UTF-8")
}

#[test]
fn dgg_answers_for_a_grid_it_knows_and_not_for_one_it_does_not() {
    assert!(dgg(&["isea3h", "info"]).contains("Refinement Ratio: 3"));
    // `IGEO7` is the library's name, not the engine's: `dgg` prints its usage and ends with 1.
    assert!(!dgg(&["IGEO7", "info"]).contains("Refinement Ratio"));
}

// The zone lists.

/// The library's name of each grid and the engine's, which `dgg` takes in lower case.
const GRIDS: [(&str, &str); 6] = [
    ("ISEA3H", "isea3h"),
    ("IVEA3H", "ivea3h"),
    ("RTEA3H", "rtea3h"),
    ("IGEO7", "isea7h_z7"),
    ("IVEA7H", "ivea7h_z7"),
    ("RTEA7H", "rtea7h_z7"),
];

/// SplitMix64, so that the sample is the same everywhere.
struct SplitMix64(u64);

impl SplitMix64 {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    /// A number in [0, 1).
    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// A hundred zones of `grid`, the same at every run: the twelve zones of level 0, two pentagons
/// under each of them at levels taken in turn, and 64 zones about random points at every level
/// in turn.
fn sample(grid: rs4dggs::AnyGrid, seed: u64) -> Vec<rs4dggs::ZoneId> {
    let mut rng = SplitMix64(seed);
    let max = grid.max_resolution();
    let mut bases = Vec::new();
    for _ in 0..100_000 {
        let z = grid
            .zone_from_geo(rng.unit() * 180.0 - 90.0, rng.unit() * 360.0 - 180.0, 0)
            .unwrap();
        if !bases.contains(&z) {
            bases.push(z);
        }
        if bases.len() == 12 {
            break;
        }
    }
    assert_eq!(bases.len(), 12, "{}: the zones of level 0", grid.name());
    let mut zones = bases.clone();
    for i in 0..24 {
        let c = grid.centroid(bases[i % 12]);
        let level = 1 + (i as u8 % max);
        zones.push(grid.zone_from_geo(c.lat, c.lon, level).unwrap());
    }
    for i in 0..64u8 {
        let level = i % (max + 1);
        let (lat, lon) = (rng.unit() * 180.0 - 90.0, rng.unit() * 360.0 - 180.0);
        zones.push(grid.zone_from_geo(lat, lon, level).unwrap());
    }
    zones
}

#[test]
fn every_identifier_of_the_binary_list_is_the_engines_64_bit_integer() {
    for (seed, (library, engine)) in GRIDS.iter().enumerate() {
        let grid = rs4dggs::get_grid(library).unwrap();
        let zones = sample(grid, 0x0067_c000 + seed as u64);
        assert_eq!(zones.len(), 100);
        // The sample holds every level, pentagons below level 0, and, on aperture 7, zones of
        // levels 1 and finer under the zones 08 to 11 of level 0, whose identifiers set bit 63.
        for level in 0..=grid.max_resolution() {
            assert!(
                zones.iter().any(|&z| grid.resolution(z) == level),
                "{library}: no zone of level {level} in the sample"
            );
        }
        let pentagons = zones
            .iter()
            .filter(|&&z| grid.is_pentagon(z) && grid.resolution(z) > 0)
            .count();
        assert!(
            pentagons >= 20,
            "{library}: {pentagons} pentagons below level 0"
        );
        if grid.refinement_ratio() == 7 {
            let high = zones
                .iter()
                .filter(|&&z| z.0 >> 63 == 1 && grid.resolution(z) > 0)
                .count();
            assert!(
                high >= 10,
                "{library}: {high} zones with bit 63 below level 0"
            );
        }
        for &zone in &zones {
            let text = grid.text_id(zone);
            let info = dgg(&[engine, "info", &text]);
            let line = info
                .lines()
                .find_map(|l| l.strip_prefix("64-bit integer ID: "))
                .unwrap_or_else(|| panic!("dgg {engine} info {text} prints no 64-bit ID:\n{info}"));
            assert_eq!(
                line,
                format!("{} (0x{:X})", zone.0, zone.0),
                "{library} {text}"
            );
        }
    }
}

#[test]
fn the_engines_list_of_level_1_is_the_golden_of_level_1() {
    let listed = dgg(&["isea3h", "list", "1"]);
    let engine: Vec<&str> = listed
        .trim()
        .trim_start_matches('[')
        .trim_end_matches(']')
        .split(',')
        .map(|s| s.trim().trim_matches('"'))
        .collect();
    let golden = std::fs::read_to_string(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/goldens/zones-ISEA3H-level1.json"),
    )
    .unwrap();
    let zones = golden
        .strip_prefix("{\"zones\":[")
        .and_then(|rest| rest.split_once(']'))
        .expect("the golden begins with its zones")
        .0;
    let golden: Vec<&str> = zones.split(',').map(|s| s.trim_matches('"')).collect();
    assert_eq!(golden.len(), 32);
    assert_eq!(engine, golden);
}

// The six grids.

#[test]
fn the_engine_knows_each_grid_by_the_last_segment_of_its_uri_and_agrees_on_its_facts() {
    for dggrs in rs4dggs_ogc::Dggrs::all() {
        let name = dggrs.uri().rsplit('/').next().unwrap();
        let grid = dggrs.grid();
        let info = dgg(&[name, "info"]);
        for line in [
            format!("Refinement Ratio: {}", grid.refinement_ratio()),
            format!(
                "Maximum level for 64-bit global identifiers (DGGAL DGGRSZone): {}",
                grid.max_resolution()
            ),
            format!(
                "Default ~64K sub-zones relative depth: {}",
                grid.depth_64k()
            ),
        ] {
            assert!(
                info.lines().any(|l| l == line),
                "dgg {name} info does not print {line:?}:\n{info}"
            );
        }
    }
    // The library's second name of ISEA7H_Z7 is unknown to the engine, so no URI may end in it.
    let igeo7 = dgg(&["IGEO7", "info"]);
    for fact in [
        "Refinement Ratio",
        "Maximum level for 64-bit global identifiers",
        "Default ~64K sub-zones relative depth",
    ] {
        assert!(!igeo7.contains(fact), "dgg IGEO7 info prints {fact:?}");
    }
}

// Zone data.

/// The zones whose data the engine reads back, on each grid: a hexagon and a pentagon at two
/// levels, and, on the two grids of the standard's examples, the zones of those examples and the
/// pentagons of level 0.
const DATA_ZONES: [(&str, &[&str]); 6] = [
    (
        "ISEA3H",
        &[
            "C4-6-A",
            "CA-0-A",
            "E7-12C9-C",
            "EB-0-B",
            "C2-23-C",
            "AA-0-A",
        ],
    ),
    ("IVEA3H", &["C4-6-A", "CA-0-A", "E7-131B-B", "EB-0-B"]),
    ("RTEA3H", &["C4-6-A", "CA-0-A", "E7-12C9-C", "EB-0-B"]),
    (
        "ISEA7H_Z7",
        &["0064", "0000", "1015505", "1100000", "00", "11"],
    ),
    ("IVEA7H_Z7", &["0064", "0000", "1015505", "1100000"]),
    ("RTEA7H_Z7", &["0064", "0000", "1015505", "1100000"]),
];

/// The zones on which `dgg togeo` places each value of a DGGS-JSON document of `zone` at
/// `depth`, whose value at each position is the position: the engine's zone at each position, in
/// order. The document is written under the tests' temporary directory.
fn engine_order(
    dggrs: &'static rs4dggs_ogc::Dggrs,
    zone: rs4dggs::ZoneId,
    depth: u8,
) -> Vec<String> {
    use rs4dggs_ogc::{Property, Values, ZoneData, write_dggs_json};
    let grid = dggrs.grid();
    let count = grid.count_sub_zones(zone, depth).unwrap();
    let positions: Vec<Option<u32>> = (0..count as u32).map(Some).collect();
    let values = [Values::U32(&positions)];
    let properties = [Property::new("i", &values)];
    let depths = [depth];
    let data = ZoneData::new(dggrs, zone, &depths, &properties);
    // A number of its own for each document, since tests that run at once may ask for one zone.
    static DOCUMENTS: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let n = DOCUMENTS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let path = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!(
        "togeo-{}-{}-{depth}-{n}.json",
        dggrs.id(),
        grid.text_id(zone)
    ));
    let file = std::io::BufWriter::new(std::fs::File::create(&path).unwrap());
    write_dggs_json(file, &data).unwrap();
    let geo = dgg(&["togeo", path.to_str().unwrap()]);
    // Each feature's properties are two lines: `"zoneID" : "<zone>",` and `"i" : <position>`.
    let mut zones = Vec::new();
    let mut pending = None;
    for line in geo.lines().map(str::trim) {
        if let Some(rest) = line.strip_prefix("\"zoneID\" : \"") {
            assert!(pending.is_none(), "two zones without a value:\n{geo}");
            pending = Some(rest.trim_end_matches(',').trim_end_matches('"').to_string());
        } else if let Some(rest) = line.strip_prefix("\"i\" : ") {
            let position: usize = rest.trim_end_matches(',').parse().unwrap();
            assert_eq!(position, zones.len(), "the values out of order:\n{geo}");
            zones.push(pending.take().expect("a value without a zone"));
        }
    }
    assert!(pending.is_none(), "a zone without a value:\n{geo}");
    zones
}

#[test]
fn the_engine_reads_each_value_of_dggs_json_on_the_librarys_sub_zone_at_its_position() {
    let mut documents = 0;
    let mut values = 0;
    for (id, zones) in DATA_ZONES {
        let dggrs = rs4dggs_ogc::Dggrs::from_id(id).unwrap();
        let grid = dggrs.grid();
        for text in zones {
            let zone = grid.zone_from_text(text).unwrap();
            for depth in 1..=3 {
                let library: Vec<String> = grid
                    .sub_zones(zone, depth)
                    .unwrap()
                    .into_iter()
                    .map(|z| grid.text_id(z))
                    .collect();
                assert_eq!(
                    engine_order(dggrs, zone, depth),
                    library,
                    "{id} {text} at depth {depth}"
                );
                documents += 1;
                values += library.len();
            }
        }
    }
    assert_eq!(documents, 84);
    assert!(values > 5_000, "{values}");
}

#[test]
fn at_a_broken_seam_the_engine_names_a_zone_where_the_library_has_none() {
    let dggrs = rs4dggs_ogc::Dggrs::from_id("ISEA7H_Z7").unwrap();
    let grid = dggrs.grid();
    let zone = grid.zone_from_text("00055353260226021").unwrap();
    let library = grid.sub_zones(zone, 1).unwrap();
    let engine = engine_order(dggrs, zone, 1);
    assert_eq!(engine.len(), 13);
    assert_eq!(library.len(), 13);
    for (i, (e, &l)) in engine.iter().zip(&library).enumerate() {
        if i == 8 || i == 12 {
            // The one difference: the null zone in the library, a zone far away in the engine.
            assert_eq!(l, rs4dggs::ZoneId::NULL, "{i}");
            assert_eq!(e, "012222222222222220", "{i}");
        } else {
            assert_eq!(*e, grid.text_id(l), "{i}");
        }
    }
}

/// Every number in `text`, in order: the characters of numbers between any others.
fn numbers(text: &str) -> Vec<f64> {
    text.split(|c: char| !(c.is_ascii_digit() || matches!(c, '-' | '+' | '.' | 'e' | 'E')))
        .filter(|t| !t.is_empty())
        .map(|t| t.parse().unwrap_or_else(|e| panic!("{t}: {e}")))
        .collect()
}

/// The positions of the coordinates of a GeoJSON geometry in `text`, the first that follows
/// `"coordinates"`, up to `end`.
fn positions(text: &str, end: &str) -> Vec<[f64; 2]> {
    let start = text.find("\"coordinates\"").expect("coordinates") + "\"coordinates\"".len();
    let rest = &text[start..];
    let n = numbers(&rest[..rest.find(end).expect("the end of the coordinates")]);
    assert_eq!(n.len() % 2, 0, "an odd count of coordinates");
    n.chunks_exact(2).map(|p| [p[0], p[1]]).collect()
}

/// For a seeded sample of zones of each grid, at every level, whose extent neither crosses the
/// antimeridian nor reaches a pole, the polygon written is the ring of `dgg geom`, position for
/// position, within 1e-9 degrees (`dgg` writes fifteen significant digits).
#[test]
fn the_region_of_a_zone_is_the_ring_dgg_draws() {
    use rs4dggs_ogc::{Coordinates, ZoneGeometry, write_zone_geometry};
    const PER_GRID: usize = 50;
    let mut rng = SplitMix64(0x5eed_0005);
    let mut unit = || rng.unit();
    for (name, engine) in GRIDS {
        let grid = rs4dggs::get_grid(name).unwrap();
        let mut seen = std::collections::HashSet::new();
        let (mut positions_compared, mut largest, mut levels) = (0, 0.0_f64, (u8::MAX, 0));
        while seen.len() < PER_GRID {
            let lat = unit() * 180.0 - 90.0;
            let lon = unit() * 360.0 - 180.0;
            let level = (unit() * f64::from(grid.max_resolution() + 1)) as u8;
            let zone = grid.zone_from_geo(lat, lon, level).unwrap();
            let Some(e) = grid.extent(zone) else { continue };
            let away = e.ll.lon <= e.ur.lon
                && e.ll.lon > -180.0
                && e.ur.lon < 180.0
                && e.ll.lat > -90.0
                && e.ur.lat < 90.0;
            if zone == rs4dggs::ZoneId::NULL || !away || !seen.insert(zone) {
                continue;
            }
            let text = grid.text_id(zone);
            let mut out = Vec::new();
            let region = ZoneGeometry::Region { edge_refinement: 0 };
            // The shortest exact form, so that the comparison is of the library's own doubles.
            let shortest = Coordinates::shortest();
            let drawn = write_zone_geometry(&mut out, grid, zone, region, shortest).unwrap();
            assert!(drawn, "{name} {text}");
            let ours = String::from_utf8(out).unwrap();
            assert!(
                ours.starts_with(r#"{"type":"Polygon","#),
                "{name} {text}: {ours}"
            );
            let ours = positions(&ours, "}");
            let theirs = positions(&dgg(&[engine, "geom", &text]), "\"properties\"");
            assert_eq!(
                ours.len(),
                theirs.len(),
                "{name} {text}: the count of positions"
            );
            assert_eq!(
                ours.first(),
                ours.last(),
                "{name} {text}: the ring is closed"
            );
            for (i, (a, b)) in ours.iter().zip(&theirs).enumerate() {
                let d = (a[0] - b[0]).abs().max((a[1] - b[1]).abs());
                largest = largest.max(d);
                assert!(d <= 1e-9, "{name} {text}, position {i}: {a:?} and {b:?}");
            }
            positions_compared += ours.len();
            levels = (levels.0.min(level), levels.1.max(level));
        }
        assert!(
            positions_compared >= 50 * 6,
            "{name}: {positions_compared} positions"
        );
        eprintln!(
            "{name}: {PER_GRID} zones of levels {} to {}, {positions_compared} positions, the \
             largest difference {largest:e} degrees",
            levels.0, levels.1
        );
    }
}

// Zone data as GeoJSON.

/// The features of the data of `zone` at `depth` as GeoJSON without geometry, the value of each
/// entry of the order being its position: the identifier and the value of each feature, in order.
fn geojson_features(
    dggrs: &'static rs4dggs_ogc::Dggrs,
    zone: rs4dggs::ZoneId,
    depth: u8,
) -> Vec<(String, usize)> {
    use rs4dggs_ogc::{
        Coordinates, Property, Values, ZoneData, ZoneGeometry, write_zone_data_geojson,
    };
    let count = dggrs.grid().count_sub_zones(zone, depth).unwrap();
    let positions: Vec<Option<u32>> = (0..count as u32).map(Some).collect();
    let values = [Values::U32(&positions)];
    let properties = [Property::new("i", &values)];
    let depths = [depth];
    let data = ZoneData::new(dggrs, zone, &depths, &properties);
    let mut out = Vec::new();
    write_zone_data_geojson(&mut out, &data, ZoneGeometry::None, Coordinates::default()).unwrap();
    let doc = String::from_utf8(out).unwrap();
    doc.split(r#""properties":{"i":"#)
        .skip(1)
        .map(|f| {
            let (value, rest) = f.split_once(r#","zoneID":""#).unwrap();
            (
                rest[..rest.find('"').unwrap()].to_string(),
                value.parse().unwrap(),
            )
        })
        .collect()
}

/// For every document the engine reads above, the features written as GeoJSON are those that
/// `dgg togeo` makes of the DGGS-JSON: one zone and one value each, in the same order.
#[test]
fn the_features_of_zone_data_are_those_the_engine_makes_of_dggs_json() {
    let mut documents = 0;
    let mut features = 0;
    for (id, zones) in DATA_ZONES {
        let dggrs = rs4dggs_ogc::Dggrs::from_id(id).unwrap();
        for text in zones {
            let zone = dggrs.grid().zone_from_text(text).unwrap();
            for depth in 1..=3 {
                let engine: Vec<(String, usize)> = engine_order(dggrs, zone, depth)
                    .into_iter()
                    .zip(0..)
                    .collect();
                let ours = geojson_features(dggrs, zone, depth);
                assert_eq!(ours, engine, "{id} {text} at depth {depth}");
                documents += 1;
                features += ours.len();
            }
        }
    }
    assert_eq!(documents, 84);
    assert!(features > 5_000, "{features}");
}

/// At a broken seam, the two entries where the engine names a zone far away and the library has
/// none are left out, and every other feature is the engine's.
#[test]
fn at_a_broken_seam_the_features_are_the_engines_but_the_two_it_places_far_away() {
    let dggrs = rs4dggs_ogc::Dggrs::from_id("ISEA7H_Z7").unwrap();
    let zone = dggrs.grid().zone_from_text("00055353260226021").unwrap();
    let engine = engine_order(dggrs, zone, 1);
    assert_eq!(engine.len(), 13);
    assert_eq!(engine[8], "012222222222222220");
    assert_eq!(engine[12], "012222222222222220");
    let expected: Vec<(String, usize)> = engine
        .into_iter()
        .zip(0..)
        .filter(|&(_, i)| i != 8 && i != 12)
        .collect();
    let ours = geojson_features(dggrs, zone, 1);
    assert_eq!(ours.len(), 11);
    assert_eq!(ours, expected);
}
