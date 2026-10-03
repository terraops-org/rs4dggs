//! DGGAL's compaction of a set of zones of one level. On the aperture-3 grids a zone two
//! levels up stands for those of its sub-zones that it alone holds, wherever all of them are
//! in the set, level after level, and the answer comes in ascending order of the identifiers.
//! On the aperture-7 grids a zone one level up stands for the one child at its centre, wherever
//! all its children are in the set, level after level, and the answer comes in the order of the
//! engine's own identifiers, which is that of the lattice of a level and not that of the Z7
//! identifiers.
//!
//! What needs no engine is tested in every build: the answers that are known by name, the
//! properties of an answer, and the refusals. The comparison with the engine's own
//! `compactZones`, sequence for sequence, is in the module `engine`, under the `oracle` feature.
use crate::interfaces::sealed::{Private, TopologyPlumbing};
use rs4dggs::indexings::Z7;
use rs4dggs::topologies::HexA7;
use rs4dggs::{AnyGrid, Error, Extent, GeoPoint, Indexing, ZoneId};
use std::collections::BTreeSet;

// The boxes of the comparison with the engine, which the oracle suites share.
#[cfg(feature = "oracle")]
use crate::common::boxes;

const APERTURE_3: [&str; 3] = ["ISEA3H", "IVEA3H", "RTEA3H"];
const APERTURE_7: [&str; 3] = ["IGEO7", "IVEA7H", "RTEA7H"];

fn grid(name: &str) -> AnyGrid {
    rs4dggs::get_grid(name).unwrap()
}

fn zone(grid: AnyGrid, text: &str) -> ZoneId {
    grid.zone_from_text(text).unwrap()
}

fn zones(grid: AnyGrid, texts: &str) -> Vec<ZoneId> {
    texts.split_whitespace().map(|t| zone(grid, t)).collect()
}

fn texts(grid: AnyGrid, zones: &[ZoneId]) -> String {
    let texts: Vec<String> = zones.iter().map(|&z| grid.text_id(z)).collect();
    texts.join(" ")
}

/// A box `[south, west, north, east]`, in degrees, as the extent a caller hands over.
fn degrees(b: &[f64; 4]) -> Extent {
    Extent {
        ll: GeoPoint {
            lat: b[0],
            lon: b[1],
        },
        ur: GeoPoint {
            lat: b[2],
            lon: b[3],
        },
    }
}

/// The zones of `level` whose extent meets the box, the first `most` of them.
fn in_box(grid: AnyGrid, level: u8, b: &[f64; 4], most: usize) -> Vec<ZoneId> {
    grid.zones_in_box(level, &degrees(b))
        .unwrap()
        .take(most)
        .collect()
}

/// A sequence of numbers that repeats from its seed: enough to shuffle a set and to draw
/// identifiers, and no more.
struct Draw(u64);

impl Draw {
    fn next(&mut self) -> u64 {
        // The multiplier and the increment of Knuth's MMIX.
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0 ^ (self.0 >> 29)
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }

    fn shuffled(&mut self, zones: &[ZoneId]) -> Vec<ZoneId> {
        let mut zones = zones.to_vec();
        for i in (1..zones.len()).rev() {
            zones.swap(i, self.below(i + 1));
        }
        zones
    }
}

/// Whether the grid is one of aperture 7.
fn aperture_7(grid: AnyGrid) -> bool {
    grid.refinement_ratio() == 7
}

/// The place of a zone in the engine's order of an answer. On the aperture-3 grids it is the
/// identifier. On the aperture-7 grids it is the key of the zone in the lattice of its level,
/// the value of the engine's own identifier, which runs through the levels two by two (a
/// level and the odd level above it share their cells) and is not the order of the Z7
/// identifiers; the zone must be one that the lattice holds.
fn order_key(grid: AnyGrid, z: ZoneId) -> u64 {
    if aperture_7(grid) {
        let place = HexA7::locate(Private, &Z7::decode(z));
        place.expect("a zone that the lattice holds").3
    } else {
        z.0
    }
}

/// The compaction of a set of zones of `level`, with what every answer must show: it is in
/// the ascending order of [`order_key`], none twice; each of its zones is of the level of the
/// set or above it, on the aperture-3 grids an even number of levels above it (or of level 0,
/// which stands for the whole of level 1); the sub-zones, at the level of the set, of the
/// zones of the answer hold every zone of the set; and the same set in other orders gives the
/// same answer.
///
/// With the answer comes the number of zones that those sub-zones hold beyond the set. It is
/// nought for a set without a gap, which `whole` says the set is: the zones of a box, the
/// sub-zones of a zone, a whole level. On the aperture-3 grids a set with a gap may be covered
/// beyond itself, the engine not looking for the zone at the centre of a coarser one; on the
/// aperture-7 grids no set is, a zone being taken only where the set has all its children.
fn compacted(
    grid: AnyGrid,
    level: u8,
    set: &[ZoneId],
    whole: bool,
    draw: &mut Draw,
    what: &str,
) -> (Vec<ZoneId>, usize) {
    let answer = grid.compact_zones(set).unwrap();
    let keys: Vec<u64> = answer.iter().map(|&z| order_key(grid, z)).collect();
    assert!(
        keys.windows(2).all(|w| w[0] < w[1]),
        "{what}: the answer is not in ascending order"
    );
    let mut expanded = BTreeSet::new();
    for &z in &answer {
        let of = grid.resolution(z);
        assert!(
            of <= level && (aperture_7(grid) || (level - of) % 2 == 0 || of == 0),
            "{what}: {} is of level {of}",
            grid.text_id(z)
        );
        expanded.extend(grid.sub_zones(z, level - of).unwrap());
    }
    let given: BTreeSet<ZoneId> = set.iter().copied().collect();
    assert!(
        expanded.is_superset(&given),
        "{what}: a zone of the set is not among the sub-zones of the answer"
    );
    let beyond = expanded.len() - given.len();
    assert!(
        !(whole || aperture_7(grid)) || beyond == 0,
        "{what}: the answer expands to {} zones, the set has {}",
        expanded.len(),
        given.len()
    );
    assert!(answer.len() <= given.len(), "{what}");
    // Three other orders of a set of up to a hundred and fifty zones, and one of a larger set.
    for _ in 0..if set.len() <= 150 { 3 } else { 1 } {
        assert!(
            grid.compact_zones(&draw.shuffled(set)).unwrap() == answer,
            "{what}: the answer depends on the order of the set"
        );
    }
    (answer, beyond)
}

/// What the thirteen sub-zones of `B4-2-A`, two levels down, compact to: the zone and the six
/// on its vertices.
const THE_SEVEN: &str = "B4-2-A C2-11-A C2-19-A C2-2C-A C4-E-A C4-11-A C4-19-A";

/// The sub-zones of one zone of ISEA3H, `B4-2-A`, the zone of level 2 that holds Lisbon. At
/// depth 2 they are thirteen: the seven that the zone alone holds give way to it, and the six
/// on its vertices, each shared with two other zones of level 2, stay. At depth 1 nothing can
/// give way, the compaction going two levels at a time. At depth 4 the zone stands for all
/// but thirty-six of the ninety-one.
#[test]
fn the_sub_zones_of_one_zone_give_way_to_it_and_those_on_its_vertices_stay() {
    let g = grid("ISEA3H");
    let parent = zone(g, "B4-2-A");
    let mut draw = Draw(0x51);

    let thirteen = g.sub_zones(parent, 2).unwrap();
    assert_eq!(thirteen.len(), 13);
    let (answer, _) = compacted(g, 4, &thirteen, true, &mut draw, "thirteen");
    assert_eq!(texts(g, &answer), THE_SEVEN);

    let seven = g.sub_zones(parent, 1).unwrap();
    assert_eq!(seven.len(), 7);
    let mut sorted = seven.clone();
    sorted.sort();
    assert_eq!(compacted(g, 3, &seven, true, &mut draw, "seven").0, sorted);

    let ninety_one = g.sub_zones(parent, 4).unwrap();
    assert_eq!(ninety_one.len(), 91);
    let (answer, _) = compacted(g, 6, &ninety_one, true, &mut draw, "ninety-one");
    assert_eq!(answer.len(), 37);
    assert_eq!(answer[0], parent);
    assert!(answer[1..].iter().all(|&z| g.resolution(z) == 6));
}

/// The engine takes a coarser zone when the six neighbours of the sub-zone at its centre are in
/// the set, and does not look for that sub-zone itself. So the thirteen sub-zones of `B4-2-A`
/// without the one at its centre compact to the same seven zones as the thirteen, which cover
/// the zone that the set lacked; without one of the six about the centre, nothing gives way.
#[test]
fn a_set_that_lacks_the_zone_at_the_centre_is_covered_beyond_itself() {
    let g = grid("ISEA3H");
    let thirteen = g.sub_zones(zone(g, "B4-2-A"), 2).unwrap();
    let within = |z: ZoneId| g.neighbors(z).iter().all(|n| thirteen.contains(n));
    let centre = *thirteen.iter().find(|&&z| within(z)).unwrap();
    assert_eq!(g.centroid(centre), g.centroid(zone(g, "B4-2-A")));
    let mut draw = Draw(0x52);

    let twelve: Vec<ZoneId> = thirteen.iter().copied().filter(|&z| z != centre).collect();
    let (answer, beyond) = compacted(g, 4, &twelve, false, &mut draw, "without the centre");
    assert_eq!(texts(g, &answer), THE_SEVEN);
    assert_eq!(beyond, 1);

    let beside = g.neighbors(centre)[0];
    let twelve: Vec<ZoneId> = thirteen.iter().copied().filter(|&z| z != beside).collect();
    let mut sorted = twelve.clone();
    sorted.sort();
    let (answer, beyond) = compacted(g, 4, &twelve, false, &mut draw, "without one beside it");
    assert_eq!(answer, sorted);
    assert_eq!(beyond, 0);
}

/// Every zone of a level compacts to the twelve zones of level 0, at every level but level 0
/// itself, which stays as it is.
#[test]
fn the_whole_world_compacts_to_the_twelve_zones_of_level_0() {
    for name in APERTURE_3 {
        let g = grid(name);
        let mut twelve: Vec<ZoneId> = g.zones(0).unwrap().collect();
        twelve.sort();
        assert_eq!(twelve.len(), 12);
        let mut draw = Draw(0x77);
        for level in 0..=7 {
            let world: Vec<ZoneId> = g.zones(level).unwrap().collect();
            let what = format!("{name}, the whole of level {level}");
            let (answer, _) = compacted(g, level, &world, true, &mut draw, &what);
            assert_eq!(answer, twelve, "{what}");
        }
        // One zone fewer at level 1, and nothing gives way.
        let mut level_1: Vec<ZoneId> = g.zones(1).unwrap().collect();
        level_1.sort();
        level_1.pop();
        assert_eq!(g.compact_zones(&level_1).unwrap(), level_1);
    }
}

/// What the zones of boxes compacted to, over a run of [`boxes_are_compacted`].
struct OfBoxes {
    sets: usize,
    smaller: usize,
    covered_beyond: usize,
    zones_in: usize,
    zones_out: usize,
}

/// The zones of five boxes, each some eight zones of the level across, at each of `levels` on
/// each of `names`, and the same with one zone in twenty taken out: every answer has the
/// properties of [`compacted`]. The boxes are about Lisbon, on either pole, over the
/// antimeridian and about a vertex of the icosahedron.
fn boxes_are_compacted(names: [&str; 3], levels: &[u8]) -> OfBoxes {
    let boxes: [[f64; 4]; 5] = [
        [38.0, -10.0, 39.0, -9.0],
        [89.0, -180.0, 90.0, 180.0],
        [-90.0, -180.0, -89.0, 180.0],
        [-1.0, 179.0, 1.0, -179.0],
        [57.9, 10.7, 58.9, 11.7],
    ];
    let mut met = OfBoxes {
        sets: 0,
        smaller: 0,
        covered_beyond: 0,
        zones_in: 0,
        zones_out: 0,
    };
    for name in names {
        let g = grid(name);
        let mut draw = Draw(0x1234);
        for &level in levels {
            // The box about each place, eight zones of the level across.
            let width = (4.0 * std::f64::consts::PI / g.count_zones(level).unwrap() as f64)
                .sqrt()
                .to_degrees();
            for b in &boxes {
                let (lat, lon) = ((b[0] + b[2]) / 2.0, (b[1] + b[3]) / 2.0);
                let scaled = if b[1] == -180.0 {
                    let cap = (4.0 * width).min(89.0);
                    if b[2] == 90.0 {
                        [90.0 - cap, -180.0, 90.0, 180.0]
                    } else {
                        [-90.0, -180.0, cap - 90.0, 180.0]
                    }
                } else if b[1] > b[3] {
                    [
                        -4.0 * width,
                        180.0 - 4.0 * width,
                        4.0 * width,
                        -180.0 + 4.0 * width,
                    ]
                } else {
                    let half = (4.0 * width).min(20.0);
                    [lat - half, lon - half, lat + half, lon + half]
                };
                let whole = in_box(g, level, &scaled, 2_000);
                let mut thinned = whole.clone();
                let mut i = draw.below(20);
                while i < thinned.len() {
                    thinned.remove(i);
                    i += 19;
                }
                for (set, whole) in [(whole, true), (thinned, false)] {
                    let what = format!("{name} level {level}, {} zones of {scaled:?}", set.len());
                    let (answer, beyond) = compacted(g, level, &set, whole, &mut draw, &what);
                    met.covered_beyond += usize::from(beyond > 0);
                    met.sets += 1;
                    met.smaller += usize::from(answer.len() < set.len());
                    met.zones_in += set.len();
                    met.zones_out += answer.len();
                }
            }
        }
    }
    eprintln!(
        "{} sets, {} made smaller, {} covered beyond themselves; {} zones in, {} out",
        met.sets, met.smaller, met.covered_beyond, met.zones_in, met.zones_out
    );
    met
}

/// The zones of boxes at every sort of level, the finest among them, and the same with one
/// zone in twenty taken out: every answer has the properties of [`compacted`], the sets
/// together are made smaller, and some of those with zones taken out are covered beyond
/// themselves.
#[test]
fn the_zones_of_a_box_compact_to_zones_that_expand_to_them() {
    let levels: Vec<u8> = (2..=12).chain([20, 26, 33]).collect();
    let met = boxes_are_compacted(APERTURE_3, &levels);
    assert!(met.covered_beyond >= 1, "{}", met.covered_beyond);
    assert!(met.sets >= 420, "{}", met.sets);
    assert!(met.smaller >= 200, "{}", met.smaller);
    assert!(
        met.zones_in >= 25_000 && met.zones_out < met.zones_in,
        "{} {}",
        met.zones_in,
        met.zones_out
    );
}

/// The same on the aperture-7 grids, at a coarse level and at the finest at which the zones
/// of a box are listed: the sets together are made smaller, and none is covered beyond
/// itself, with or without zones taken out.
#[test]
fn on_aperture_7_the_zones_of_a_box_compact_to_zones_that_expand_to_them() {
    let met = boxes_are_compacted(APERTURE_7, &[4, 14]);
    assert_eq!(met.covered_beyond, 0);
    assert_eq!(met.sets, 60);
    assert!(met.smaller >= 30, "{}", met.smaller);
    assert!(
        met.zones_in >= 3_000 && met.zones_out < met.zones_in,
        "{} {}",
        met.zones_in,
        met.zones_out
    );
}

/// A zone given twice counts once; a slice of no zone gives none; and the answer is a list
/// of its own, which leaves the slice as it was.
#[test]
fn a_zone_given_twice_counts_once_and_no_zone_gives_none() {
    let g = grid("ISEA3H");
    let thirteen = g.sub_zones(zone(g, "B4-2-A"), 2).unwrap();
    let mut twice = thirteen.clone();
    twice.extend(&thirteen);
    twice.push(thirteen[4]);
    assert_eq!(g.compact_zones(&twice), g.compact_zones(&thirteen));
    assert_eq!(g.compact_zones(&thirteen).unwrap().len(), 7);
    for name in APERTURE_3 {
        assert_eq!(grid(name).compact_zones(&[]), Ok(Vec::new()));
    }
    // One zone is its own compaction, at any level.
    for text in ["A0-0-A", "A5-0-B", "B4-2-A", "Q3-6954FE0D5D2E0-A", "QA-0-B"] {
        assert_eq!(
            g.compact_zones(&[zone(g, text)]),
            Ok(vec![zone(g, text)]),
            "{text}"
        );
    }
    // The grid taken by its type answers as the handle on it.
    assert_eq!(
        rs4dggs::isea3h().compact_zones(&thirteen),
        g.compact_zones(&thirteen)
    );
}

/// Zones of more than one level are refused, with the two levels first met: the engine's
/// compaction drops a zone that has no zone two levels up to give way to as soon as the set
/// holds a finer one, so that `A2-0-C B4-2-A` would lose `A2-0-C`, and a caller would lose a
/// zone in silence.
#[test]
fn zones_of_more_than_one_level_are_refused() {
    let g = grid("ISEA3H");
    for (set, level, other) in [
        ("A2-0-C B4-2-A", 1, 2),
        ("A2-0-C C4-22-B", 1, 5),
        ("B5-0-A A9-0-C A0-0-B A3-0-A", 2, 1),
        ("A0-0-A B4-2-A", 0, 2),
        ("C4-E-A C4-11-A B4-2-A", 4, 2),
    ] {
        assert_eq!(
            g.compact_zones(&zones(g, set)),
            Err(Error::MixedLevels { level, other }),
            "{set}"
        );
    }
    // An answer of the compaction is of more than one level, and is no set to compact again.
    let answer = g
        .compact_zones(&g.sub_zones(zone(g, "B4-2-A"), 2).unwrap())
        .unwrap();
    assert_eq!(
        g.compact_zones(&answer),
        Err(Error::MixedLevels { level: 2, other: 4 })
    );
    assert_eq!(
        Error::MixedLevels { level: 2, other: 4 }.to_string(),
        "the zones are of more than one level: 2 and 4"
    );
}

/// An identifier that is no zone of the grid is refused, wherever it stands in the slice and
/// before the levels are compared: the null zone, which the engine leaves out of the set, the
/// identifiers that the grid reads no cell for, and any pattern of bits. Nothing is refused
/// by a panic, and nothing costs more than the slice is long.
#[test]
fn an_identifier_that_is_no_zone_is_refused() {
    for name in APERTURE_3 {
        let g = grid(name);
        let thirteen = g.sub_zones(zone(g, "B4-2-A"), 2).unwrap();
        let hostile = [
            ZoneId::NULL,
            ZoneId(u64::MAX - 3),
            ZoneId(1 << 63),
            ZoneId(1 << 62),
            ZoneId(0x7fff_ffff_ffff_ffff),
            // A polar root with an index, and the sub-hexes of a pole that no zone has.
            zone(g, "BA-0-A").0.checked_add(4).map(ZoneId).unwrap(),
            ZoneId(zone(g, "AA-0-B").0 + 1),
            ZoneId(zone(g, "AA-0-B").0 + 2),
        ];
        for bad in hostile {
            for at in [0, 6, 13] {
                let mut set = thirteen.clone();
                set.insert(at, bad);
                assert!(
                    matches!(g.compact_zones(&set), Err(Error::InvalidZone(_))),
                    "{name}: {:#x} at {at}",
                    bad.0
                );
            }
            assert!(matches!(
                g.compact_zones(&[bad]),
                Err(Error::InvalidZone(_))
            ));
            // Refused before the levels are compared.
            let mixed = [zone(g, "A0-0-A"), zone(g, "B4-2-A"), bad];
            assert!(matches!(
                g.compact_zones(&mixed),
                Err(Error::InvalidZone(_))
            ));
        }
        // Any pattern of bits: a zone, compacted to itself, or a refusal.
        let mut draw = Draw(0xbad);
        let (mut read, mut refused) = (0, 0);
        for i in 0..20_000 {
            let bits = match i % 4 {
                0 => draw.next(),
                1 => draw.next() >> 2,
                2 => draw.next() >> (draw.below(60) + 2),
                _ => thirteen[draw.below(13)].0 ^ (1 << draw.below(64)),
            };
            match g.compact_zones(&[ZoneId(bits)]) {
                Ok(answer) => {
                    assert_eq!(answer, [ZoneId(bits)]);
                    read += 1;
                }
                Err(Error::InvalidZone(_)) => refused += 1,
                Err(other) => panic!("{name}: {bits:#x} gives {other:?}"),
            }
        }
        assert!(
            read >= 1_000 && refused >= 2_000,
            "{name}: {read} {refused}"
        );
        // A slice of any zones the grid reads, of one level, far from one another.
        for level in [0, 1, 2, 5, 12, 21, 32, 33] {
            let mut set = Vec::new();
            while set.len() < 300 {
                let lat = (draw.below(1_800_001) as f64 - 900_000.0) / 10_000.0;
                let lon = (draw.below(3_600_001) as f64 - 1_800_000.0) / 10_000.0;
                let z = g.zone_from_geo(lat, lon, level).unwrap();
                if z != ZoneId::NULL {
                    set.push(z);
                }
            }
            let what = format!("{name}, zones of level {level} far from one another");
            compacted(g, level, &set, false, &mut draw, &what);
        }
    }
}

/// What the thirteen sub-zones of `0064`, one level down, compact to on ISEA7H_Z7: the zone in
/// the place of the child at its centre, `00640`, and the twelve others, of which six are the
/// zones that the identifier names with one digit more and six lie across its vertices. The
/// zone of level 2 stands among those of level 3: the order is not that of the levels.
const THE_THIRTEEN: &str =
    "00422 03332 03323 0064 00641 00645 00644 00646 00642 00643 00656 00665 00604";

/// On the aperture-7 grids a zone gives way to a coarser one where the set has all the
/// children of each of its parents, and a zone that is not at the centre of its parent has
/// two: so the thirteen children of one zone give way at its centre alone, thirteen zones for
/// thirteen, the coarser one overlapping the twelve that stay. The seven zones that the
/// identifier names with one digit more are not all the children, and nothing gives way; they
/// come back in the engine's order, which is not that of the digits. A pentagon has eleven
/// children, and stands for the one at its centre likewise.
#[test]
fn on_aperture_7_a_zone_stands_for_the_child_at_its_centre_alone() {
    let g = grid("IGEO7");
    let mut draw = Draw(0x71);

    let thirteen = g.sub_zones(zone(g, "0064"), 1).unwrap();
    assert_eq!(thirteen.len(), 13);
    let (answer, _) = compacted(g, 3, &thirteen, true, &mut draw, "thirteen");
    assert_eq!(texts(g, &answer), THE_THIRTEEN);

    let seven = zones(g, "00640 00641 00642 00643 00644 00645 00646");
    let (answer, _) = compacted(g, 3, &seven, false, &mut draw, "seven");
    assert_eq!(
        texts(g, &answer),
        "00640 00641 00645 00644 00646 00642 00643"
    );

    let eleven = g.sub_zones(zone(g, "00"), 1).unwrap();
    assert_eq!(eleven.len(), 11);
    let (answer, _) = compacted(g, 1, &eleven, true, &mut draw, "eleven");
    assert_eq!(
        texts(g, &answer),
        "013 023 033 043 053 00 005 004 006 003 001"
    );

    // Two levels down the zone of level 1 at the centre stands for one zone of forty-six.
    let forty_six = g.sub_zones(zone(g, "00"), 2).unwrap();
    assert_eq!(forty_six.len(), 46);
    let (answer, _) = compacted(g, 2, &forty_six, true, &mut draw, "forty-six");
    assert_eq!(answer.len(), 46);
    assert_eq!(texts(g, &answer[..1]), "000");
    assert!(answer[1..].iter().all(|&z| g.resolution(z) == 2));
    assert!(!answer.contains(&zone(g, "0000")));
}

/// On the aperture-7 grids a coarser zone is taken only where the set has every one of its
/// children, the one at its centre among them: the thirteen children of `0064` without any
/// one of them come back as they are, in the engine's order, and no set is covered beyond
/// itself.
#[test]
fn on_aperture_7_a_set_with_a_gap_is_not_covered_beyond_itself() {
    let g = grid("IGEO7");
    let thirteen = g.sub_zones(zone(g, "0064"), 1).unwrap();
    let mut draw = Draw(0x72);
    for missing in ["00640", "00641", "00422"] {
        let mut twelve: Vec<ZoneId> = thirteen
            .iter()
            .copied()
            .filter(|&z| z != zone(g, missing))
            .collect();
        assert_eq!(twelve.len(), 12);
        let what = format!("without {missing}");
        let (answer, beyond) = compacted(g, 3, &twelve, false, &mut draw, &what);
        twelve.sort_by_key(|&z| order_key(g, z));
        assert_eq!(answer, twelve, "{what}");
        assert_eq!(beyond, 0);
    }
}

/// Every zone of a level of an aperture-7 grid compacts to the twelve zones of level 0, in
/// the order of the level; one zone fewer at level 1, and the two zones of level 0 that are
/// its parents are not taken.
#[test]
fn on_aperture_7_the_whole_world_compacts_to_the_twelve_zones_of_level_0() {
    for name in APERTURE_7 {
        let g = grid(name);
        let twelve: Vec<ZoneId> = g.zones(0).unwrap().collect();
        assert_eq!(twelve.len(), 12);
        let mut draw = Draw(0x77);
        for level in 0..=WORLD_7 {
            let world: Vec<ZoneId> = g.zones(level).unwrap().collect();
            let what = format!("{name}, the whole of level {level}");
            let (answer, _) = compacted(g, level, &world, true, &mut draw, &what);
            assert_eq!(answer, twelve, "{what}");
        }
        let mut level_1: Vec<ZoneId> = g.zones(1).unwrap().collect();
        assert_eq!(level_1.len(), 72);
        let last = level_1.pop().unwrap();
        let parents = g.parents(last);
        assert_eq!(parents.len(), 2);
        let answer = g.compact_zones(&level_1).unwrap();
        assert!(parents.iter().all(|p| !answer.contains(p)), "{name}");
        let (level_0, stay): (Vec<ZoneId>, Vec<ZoneId>) =
            answer.iter().partition(|&&z| g.resolution(z) == 0);
        assert_eq!(level_0.len(), 10, "{name}");
        // A zone of level 1 stays where one of its parents is one of the two.
        assert!(!stay.is_empty(), "{name}");
        for z in stay {
            assert!(g.parents(z).iter().any(|p| parents.contains(p)), "{name}");
        }
    }
}

/// The finest level of an aperture-7 grid whose every zone is compacted here in every build.
const WORLD_7: u8 = 3;

/// The Z7 identifier of a base cell and its digits, whatever they name: also the child a
/// pentagon does not have, and twenty digits, which no text may name.
fn z7(base: u64, digits: &[u8]) -> ZoneId {
    Z7::encode(&rs4dggs::Address::new(base, digits))
}

/// On the aperture-7 grids as on the aperture-3 ones: no zone gives no zone, one zone is its
/// own compaction at any level, a zone given twice counts once, and zones of two levels are
/// refused. An identifier that is no zone of the grid is refused wherever it stands and
/// before the levels are compared: the null zone, a base cell beyond the twelve, and a Z7
/// identifier of twenty digits, which the engine reads as its null zone. Any other pattern of
/// bits is read as a zone and answered, without a panic and at the cost of one zone.
#[test]
fn on_aperture_7_the_refusals_are_those_of_aperture_3() {
    for name in APERTURE_7 {
        let g = grid(name);
        assert_eq!(g.compact_zones(&[]), Ok(Vec::new()), "{name}");
        for text in [
            "00",
            "11",
            "060",
            "0064",
            "006415654636",
            "006415654636011111111",
        ] {
            assert_eq!(
                g.compact_zones(&[zone(g, text)]),
                Ok(vec![zone(g, text)]),
                "{name} {text}"
            );
        }
        let thirteen = g.sub_zones(zone(g, "0064"), 1).unwrap();
        let mut twice = thirteen.clone();
        twice.extend(&thirteen);
        twice.push(thirteen[4]);
        assert_eq!(g.compact_zones(&twice), g.compact_zones(&thirteen));
        assert_eq!(texts(g, &g.compact_zones(&twice).unwrap()), THE_THIRTEEN);

        let hostile = [
            ZoneId::NULL,
            ZoneId(0xc3ff_ffff_ffff_ffff),
            ZoneId(0xf000_0000_0000_0000),
            // Twenty digits: nought throughout, and under a zone of level 3.
            z7(0, &[0; 20]),
            z7(3, &[1; 20]),
        ];
        for bad in hostile {
            for at in [0, 6, 13] {
                let mut set = thirteen.clone();
                set.insert(at, bad);
                assert!(
                    matches!(g.compact_zones(&set), Err(Error::InvalidZone(_))),
                    "{name}: {:#x} at {at}",
                    bad.0
                );
            }
            assert!(
                matches!(g.compact_zones(&[bad]), Err(Error::InvalidZone(_))),
                "{name}: {:#x}",
                bad.0
            );
            // Refused before the levels are compared.
            let mixed = [zone(g, "00"), zone(g, "0064"), bad];
            assert!(matches!(
                g.compact_zones(&mixed),
                Err(Error::InvalidZone(_))
            ));
        }
        for (set, level, other) in [("00 0064", 0, 2), ("05 011", 0, 1), ("0064 00641 00", 2, 3)] {
            assert_eq!(
                g.compact_zones(&zones(g, set)),
                Err(Error::MixedLevels { level, other }),
                "{name} {set}"
            );
        }
        // An answer of the compaction is of two levels, and is no set to compact again.
        let answer = g.compact_zones(&thirteen).unwrap();
        assert_eq!(
            g.compact_zones(&answer),
            Err(Error::MixedLevels { level: 3, other: 2 })
        );
        // The grid taken by its type answers as the handle on it.
        assert_eq!(
            rs4dggs::igeo7().compact_zones(&thirteen).unwrap(),
            grid("IGEO7").compact_zones(&thirteen).unwrap()
        );

        // Any pattern of bits: refused, or read as one zone and answered with one zone of its
        // level, which is the zone the engine reads the bits as. That is the zone the bits
        // name, save where they name the child a pentagon does not have.
        let mut draw = Draw(0xbad7);
        let (mut itself, mut another, mut refused) = (0, 0, 0);
        for i in 0..3_000 {
            let bits = match i % 4 {
                0 => draw.next(),
                1 => draw.next() | (0x0000_0fff_ffff_ffff >> (3 * draw.below(14))),
                2 => draw.next() | u64::MAX >> (4 + 3 * draw.below(20)),
                _ => thirteen[draw.below(13)].0 ^ (1 << draw.below(64)),
            };
            match g.compact_zones(&[ZoneId(bits)]) {
                Ok(answer) => {
                    assert_eq!(answer.len(), 1, "{name}: {bits:#x}");
                    assert_eq!(g.resolution(answer[0]), g.resolution(ZoneId(bits)));
                    let named = Z7::encode(&Z7::decode(ZoneId(bits)));
                    if answer[0] == named {
                        itself += 1;
                    } else {
                        assert!(HexA7::locate(Private, &Z7::decode(ZoneId(bits))).is_none());
                        another += 1;
                    }
                }
                Err(Error::InvalidZone(_)) => refused += 1,
                Err(other) => panic!("{name}: {bits:#x} gives {other:?}"),
            }
        }
        eprintln!(
            "{name}: {itself} zones answered as themselves, {another} as another, {refused} refused"
        );
        assert!(
            itself >= 2_000 && another >= 200 && refused >= 550,
            "{name}: {itself} {another} {refused}"
        );
        // The child a pentagon does not have is read as one that it has.
        assert_eq!(g.compact_zones(&[z7(0, &[2])]), Ok(vec![zone(g, "003")]));
        assert_eq!(g.compact_zones(&[z7(6, &[5])]), Ok(vec![zone(g, "064")]));

        // A slice of any zones the grid reads, of one level, far from one another.
        for level in [0, 1, 5, 14, 19] {
            let mut set = Vec::new();
            while set.len() < 120 {
                let lat = (draw.below(1_800_001) as f64 - 900_000.0) / 10_000.0;
                let lon = (draw.below(3_600_001) as f64 - 1_800_000.0) / 10_000.0;
                let z = g.zone_from_geo(lat, lon, level).unwrap();
                if z != ZoneId::NULL {
                    set.push(z);
                }
            }
            let what = format!("{name}, zones of level {level} far from one another");
            compacted(g, level, &set, false, &mut draw, &what);
        }
    }
}

/// The children of `00000000000000005`, a zone of level 15 beside the pole of base cell 0,
/// where the engine's lists of parents and children are not consistent with themselves:
/// eleven zones of level 16.
const SEAM_ELEVEN: &str = "000000000000000050 000000000000000053 000000000000000051 \
     000000000000000055 000000000000000054 000000000000000056 000000000000000052 \
     000000000000000004 000000000000000001 000000000000000010 000000000000000522";

/// What the engine compacts [`SEAM_ELEVEN`] to: the zone of level 15 in the place of three of
/// the eleven, although the set has not all that the engine lists as its children.
const SEAM_NINE: &str = "00000000000000005 000000000000000051 000000000000000055 \
     000000000000000052 000000000000000056 000000000000000004 000000000000000001 \
     000000000000000010 000000000000000522";

/// Two sets of level 15 in the same seam, each the children of a zone of level 14 whose list
/// in the engine holds its null zone besides, four times for the hexagon and once for the
/// pentagon; and what the engine compacts each to. The null entries are passed over, and the
/// zone of level 14 is taken.
const SEAM_WITH_NULL_CHILDREN: [(&str, &str); 2] = [
    (
        "00000000000000050 00000000000000052 00000000000000053 00000000000000051 \
         00000000000000055 00000000000000054 00000000000000056 00000000000000005 \
         00000000000000532",
        "00000000000000532 0000000000000005 00000000000000052 00000000000000053 \
         00000000000000051 00000000000000056 00000000000000005",
    ),
    (
        "00000000000000000 00000000000000005 00000000000000004 00000000000000006 \
         00000000000000003 00000000000000001 00000000000000016 00000000000000061 \
         00000000000000052 00000000000000034",
        "00000000000000052 00000000000000061 00000000000000034 00000000000000016 \
         0000000000000000 00000000000000004 00000000000000006 00000000000000003",
    ),
];

/// In the broken seam over the pole of base cell 0 the answer is the engine's, whatever it
/// is, each of these read from the engine's `compactZones` on the three grids: a set of
/// level 16 in which a zone stands for three; two sets of level 15 whose parent is taken
/// although the engine lists its null zone among the children; a zone of level 18 for which
/// the engine finds no parent, which leaves the answer with nothing in its place; and a zone
/// of level 17 that the engine reads as another, which comes back as that one.
///
/// And the one departure: for `0100000000000001623` the engine answers an identifier that it
/// cannot itself read back, `0122222222222222203`, which this crate never hands out; the
/// answer here is empty.
#[test]
fn in_a_broken_seam_the_answer_is_the_engines() {
    for name in APERTURE_7 {
        let g = grid(name);
        let answer = g.compact_zones(&zones(g, SEAM_ELEVEN)).unwrap();
        assert_eq!(texts(g, &answer), SEAM_NINE, "{name}");
        for (set, compact) in SEAM_WITH_NULL_CHILDREN {
            let answer = g.compact_zones(&zones(g, set)).unwrap();
            assert_eq!(texts(g, &answer), compact, "{name}");
        }
        assert_eq!(
            g.compact_zones(&[zone(g, "00000000000000001311")]),
            Ok(Vec::new()),
            "{name}"
        );
        assert_eq!(
            g.compact_zones(&[zone(g, "0000000000000001644")]),
            Ok(vec![zone(g, "0000000000000000056")]),
            "{name}"
        );
        assert_eq!(
            g.compact_zones(&[zone(g, "0100000000000001623")]),
            Ok(Vec::new()),
            "{name}"
        );
    }
}

#[cfg(feature = "oracle")]
mod engine {
    use super::boxes::{Kind, Rng, make_box, sites, zone_width_deg};
    use super::*;
    use dggal_oracle::NULL_ZONE;

    fn ids(zones: &[ZoneId]) -> Vec<u64> {
        zones.iter().map(|z| z.0).collect()
    }

    /// The engine's `compactZones` of a set, as zones.
    fn the_engines(engine: &str, set: &[ZoneId]) -> Vec<ZoneId> {
        dggal_oracle::compact_zones(engine, &ids(set))
            .into_iter()
            .map(ZoneId)
            .collect()
    }

    /// What a run over the sets of a grid met, at one level.
    #[derive(Default, Debug)]
    struct Tally {
        /// Sets compared with the engine: the zones of boxes, the same thinned, the sub-zones
        /// of one zone, the whole world.
        of_boxes: usize,
        /// Of those, the sets that are the first zones of a larger answer.
        pages: usize,
        thinned: usize,
        of_one_zone: usize,
        world: usize,
        /// Sets that the compaction made smaller.
        smaller: usize,
        /// Sets, of those with zones taken out or cut to a page, that the answer covers beyond
        /// themselves.
        covered_beyond: usize,
        /// Zones in the sets, and in the answers.
        zones_in: usize,
        zones_out: usize,
    }

    /// The crate's compaction of a set against the engine's, sequence for sequence, with the
    /// properties of [`compacted`].
    fn is_the_engines(
        g: AnyGrid,
        engine: &str,
        level: u8,
        set: &[ZoneId],
        whole: bool,
        draw: &mut Draw,
        tally: &mut Tally,
        what: &str,
    ) {
        let (ours, beyond) = compacted(g, level, set, whole, draw, what);
        tally.covered_beyond += usize::from(beyond > 0);
        let theirs = the_engines(engine, set);
        assert!(
            ours == theirs,
            "{what}: {} zones where the engine has {}; the first to differ {:?}",
            ours.len(),
            theirs.len(),
            ours.iter().zip(&theirs).find(|(a, b)| a != b)
        );
        tally.smaller += usize::from(ours.len() < set.len());
        tally.zones_in += set.len();
        tally.zones_out += ours.len();
    }

    /// The aperture-7 grids, each with the name of the engine's class.
    const ENGINES_7: [(&str, &str); 3] = [
        ("IGEO7", dggal_oracle::IGEO7),
        ("IVEA7H", dggal_oracle::IVEA7H),
        ("RTEA7H", dggal_oracle::RTEA7H),
    ];

    /// The kinds of box whose zones are compacted: boxes anywhere, over the antimeridian, on
    /// a pole, about a vertex of the icosahedron and about the vertices and centroids of
    /// zones, each some zones of the level across.
    const KINDS: [Kind; 5] = [
        Kind::Rand,
        Kind::Anti,
        Kind::Pole,
        Kind::Vertex,
        Kind::Corner,
    ];

    /// The most zones of a box that are compacted on aperture 3. A box that holds more gives
    /// its first, in the order of the level, which is a page of its answer and no longer a set
    /// without a gap: the zone on a pole comes last in that order, after its neighbours, so
    /// that a page may hold the neighbours without it.
    const MOST: usize = 600;

    /// The sets of a run on the grids of one aperture.
    struct Sets {
        /// The levels at which sets are compared with the engine.
        levels: &'static [u8],
        /// The finest of them at which the whole world is a set too.
        world_to: u8,
        /// The depths at which the sub-zones of a zone are a set.
        depths: [u8; 3],
        /// The boxes of each of [`KINDS`] at a level, and the most zones of a box that are
        /// compacted: fewer on aperture 7, where a zone costs several times what it costs on
        /// aperture 3, the parents and the children of a zone being found by geometry.
        boxes: usize,
        most: usize,
        /// The fewest sets made smaller, and the fewest zones in the sets, at each level.
        smaller: usize,
        zones_in: usize,
    }

    /// On aperture 3: the coarse levels and three finer ones, the finest among them; the
    /// sub-zones of a zone from two levels down, the compaction going two levels at a time.
    const SETS_3: Sets = Sets {
        levels: &[2, 3, 4, 5, 6, 7, 8, 12, 20, 33],
        world_to: 8,
        depths: [2, 3, 4],
        boxes: 30,
        most: MOST,
        smaller: 80,
        zones_in: 6_000,
    };

    /// On aperture 7: the coarse levels, at which the engine's cost is that of the set; the
    /// sub-zones of a zone from one level down, its thirteen children.
    const SETS_7: Sets = Sets {
        levels: &[2, 3, 4, 5, 6],
        world_to: 3,
        depths: [1, 2, 3],
        boxes: 12,
        most: 300,
        smaller: 20,
        zones_in: 6_000,
    };

    /// At each level of `sets`, the engine's sequence for some two hundred sets of zones of
    /// the level, or some hundred on aperture 7: the zones of boxes of each of [`KINDS`];
    /// every other of those sets of twenty zones or more with one zone in twenty taken out;
    /// the sub-zones of zones, at three depths; and, at the coarser levels, the whole world,
    /// which the engine too compacts to the twelve zones of level 0.
    fn the_compaction_is_the_engines(
        name: &str,
        engine: &str,
        seed: u64,
        sets: &Sets,
    ) -> Vec<Tally> {
        let g = grid(name);
        let mut tallies = Vec::new();
        for &level in sets.levels {
            let mut tally = Tally::default();
            let mut rng = Rng::new(seed + u64::from(level));
            let mut draw = Draw(seed ^ u64::from(level));
            let sites = sites(engine, i32::from(level), &mut rng);
            // Boxes four times the size of those the generator makes for a zone of the level,
            // so that most of them hold dozens of zones or hundreds.
            let size = 4.0 * zone_width_deg(engine, i32::from(level));
            for kind in KINDS {
                for i in 0..sets.boxes {
                    let b = make_box(kind, &mut rng, size, &sites);
                    let set = in_box(g, level, &b, sets.most);
                    let what = format!(
                        "{name} level {level}, the {} zones of the {} box {b:?}",
                        set.len(),
                        kind.name()
                    );
                    let whole = set.len() < sets.most;
                    is_the_engines(g, engine, level, &set, whole, &mut draw, &mut tally, &what);
                    tally.of_boxes += 1;
                    tally.pages += usize::from(!whole);
                    if i % 2 == 0 && set.len() >= 20 {
                        let mut thinned = set.clone();
                        let mut at = draw.below(20);
                        while at < thinned.len() {
                            thinned.remove(at);
                            at += 19;
                        }
                        let what = format!("{what}, one in twenty taken out");
                        is_the_engines(
                            g, engine, level, &thinned, false, &mut draw, &mut tally, &what,
                        );
                        tally.thinned += 1;
                    }
                }
            }
            for depth in sets.depths {
                let mut made = 0;
                while made < 7 {
                    let (lat, lon) = rng.point();
                    let parent = g.zone_from_geo(lat, lon, level - depth.min(level)).unwrap();
                    if parent == ZoneId::NULL || depth > level {
                        made += usize::from(depth > level);
                        continue;
                    }
                    let set = g.sub_zones(parent, depth).unwrap();
                    let what = format!(
                        "{name} level {level}, the {} sub-zones of {} at depth {depth}",
                        set.len(),
                        g.text_id(parent)
                    );
                    is_the_engines(g, engine, level, &set, true, &mut draw, &mut tally, &what);
                    tally.of_one_zone += 1;
                    made += 1;
                }
            }
            if level <= sets.world_to {
                let world: Vec<ZoneId> = g.zones(level).unwrap().collect();
                let what = format!("{name}, the whole of level {level}");
                is_the_engines(g, engine, level, &world, true, &mut draw, &mut tally, &what);
                assert_eq!(the_engines(engine, &world).len(), 12, "{what}");
                tally.world += 1;
            }
            eprintln!("{name} level {level}: {tally:?}");
            tallies.push(tally);
        }
        tallies
    }

    /// The floors of a run: at every level the sets of each sort, many of them made smaller,
    /// and thousands of zones. On aperture 7 no set is covered beyond itself.
    fn has_at_least(tallies: &[Tally], sets: &Sets, aperture_7: bool) {
        assert_eq!(tallies.len(), sets.levels.len());
        for (tally, &level) in tallies.iter().zip(sets.levels) {
            assert_eq!(tally.of_boxes, 5 * sets.boxes, "level {level}");
            assert!(tally.thinned >= sets.boxes, "level {level}: {tally:?}");
            // The sub-zones of seven zones at each of the depths the level allows.
            let depths = sets.depths.iter().filter(|&&depth| depth <= level).count();
            assert_eq!(tally.of_one_zone, 7 * depths, "level {level}");
            assert_eq!(
                tally.world,
                usize::from(level <= sets.world_to),
                "level {level}"
            );
            assert!(tally.smaller >= sets.smaller, "level {level}: {tally:?}");
            assert!(tally.zones_in >= sets.zones_in, "level {level}: {tally:?}");
            assert!(tally.zones_out < tally.zones_in, "level {level}: {tally:?}");
            assert!(
                !aperture_7 || tally.covered_beyond == 0,
                "level {level}: {tally:?}"
            );
        }
    }

    #[test]
    fn the_compaction_is_the_engines_on_isea3h() {
        let tallies =
            the_compaction_is_the_engines("ISEA3H", dggal_oracle::ISEA3H, 0xc0_0000, &SETS_3);
        has_at_least(&tallies, &SETS_3, false);
    }

    #[test]
    fn the_compaction_is_the_engines_on_ivea3h() {
        let tallies =
            the_compaction_is_the_engines("IVEA3H", dggal_oracle::IVEA3H, 0xc1_0000, &SETS_3);
        has_at_least(&tallies, &SETS_3, false);
    }

    #[test]
    fn the_compaction_is_the_engines_on_rtea3h() {
        let tallies =
            the_compaction_is_the_engines("RTEA3H", dggal_oracle::RTEA3H, 0xc2_0000, &SETS_3);
        has_at_least(&tallies, &SETS_3, false);
    }

    #[test]
    fn the_compaction_is_the_engines_on_igeo7() {
        let tallies =
            the_compaction_is_the_engines("IGEO7", dggal_oracle::IGEO7, 0xc7_0000, &SETS_7);
        has_at_least(&tallies, &SETS_7, true);
    }

    #[test]
    fn the_compaction_is_the_engines_on_ivea7h() {
        let tallies =
            the_compaction_is_the_engines("IVEA7H", dggal_oracle::IVEA7H, 0xc8_0000, &SETS_7);
        has_at_least(&tallies, &SETS_7, true);
    }

    #[test]
    fn the_compaction_is_the_engines_on_rtea7h() {
        let tallies =
            the_compaction_is_the_engines("RTEA7H", dggal_oracle::RTEA7H, 0xc9_0000, &SETS_7);
        has_at_least(&tallies, &SETS_7, true);
    }

    /// A set of the finest level, by itself: the sub-zones of a zone of level 29 and the zones
    /// of a box about Lisbon, as the engine compacts them.
    #[test]
    fn a_set_of_the_finest_level_is_compacted_as_the_engine_does() {
        let g = grid("ISEA3H");
        let mut draw = Draw(0x33);
        let mut tally = Tally::default();
        let parent = g.zone_from_geo(38.7223, -9.1393, 29).unwrap();
        let set = g.sub_zones(parent, 4).unwrap();
        assert_eq!(set.len(), 91);
        is_the_engines(
            g,
            dggal_oracle::ISEA3H,
            33,
            &set,
            true,
            &mut draw,
            &mut tally,
            "level 33",
        );
        let b = [38.7223, -9.1393, 38.7223 + 1.5e-5, -9.1393 + 1.5e-5];
        let set = in_box(g, 33, &b, MOST);
        assert!(set.len() >= 100, "{}", set.len());
        is_the_engines(
            g,
            dggal_oracle::ISEA3H,
            33,
            &set,
            true,
            &mut draw,
            &mut tally,
            "level 33",
        );
        assert_eq!(tally.smaller, 2);
    }

    /// The same on the aperture-7 grids, far from the two broken seams: the fifty-five
    /// sub-zones, two levels down, of the zone of level 17 that holds Lisbon, which are of
    /// level 19, and the zones of a box about it at level 14, the finest at which the zones of
    /// a box are listed.
    #[test]
    fn on_aperture_7_a_set_of_the_finest_level_is_compacted_as_the_engine_does() {
        for (name, engine) in ENGINES_7 {
            let g = grid(name);
            let mut draw = Draw(0x37);
            let mut tally = Tally::default();
            let parent = g.zone_from_geo(38.7223, -9.1393, 17).unwrap();
            let set = g.sub_zones(parent, 2).unwrap();
            assert_eq!(set.len(), 55);
            let what = format!("{name} level 19");
            is_the_engines(g, engine, 19, &set, true, &mut draw, &mut tally, &what);
            let b = [38.7223, -9.1393, 38.7223 + 2e-3, -9.1393 + 2e-3];
            let set = in_box(g, 14, &b, MOST);
            assert!(set.len() >= 100, "{}", set.len());
            let what = format!("{name} level 14");
            is_the_engines(g, engine, 14, &set, true, &mut draw, &mut tally, &what);
            eprintln!("{name}: {tally:?}");
            assert!(tally.smaller >= 1, "{name}");
        }
    }

    /// The whole world, by itself: every zone of level 8 of an aperture-3 grid, and of level 4
    /// of an aperture-7 grid, compacts to the twelve of level 0, in the engine as here.
    #[test]
    fn the_whole_world_is_compacted_as_the_engine_does() {
        for (name, engine) in [
            ("ISEA3H", dggal_oracle::ISEA3H),
            ("IVEA3H", dggal_oracle::IVEA3H),
            ("RTEA3H", dggal_oracle::RTEA3H),
        ] {
            let g = grid(name);
            let world: Vec<ZoneId> = g.zones(8).unwrap().collect();
            assert_eq!(world.len(), 65_612);
            let theirs = the_engines(engine, &world);
            assert_eq!(theirs.len(), 12);
            assert_eq!(g.compact_zones(&world).unwrap(), theirs, "{name}");
        }
        // On the aperture-7 grids, every zone of level 4.
        for (name, engine) in ENGINES_7 {
            let g = grid(name);
            let world: Vec<ZoneId> = g.zones(4).unwrap().collect();
            assert_eq!(world.len(), 24_012);
            let theirs = the_engines(engine, &world);
            assert_eq!(theirs.len(), 12);
            assert_eq!(g.compact_zones(&world).unwrap(), theirs, "{name}");
        }
    }

    /// The engine leaves the null zone out of the set and counts a zone given twice once: its
    /// answer for a set with either is its answer for the set without. This crate refuses the
    /// null zone, as any identifier that is no zone, and counts a zone given twice once.
    #[test]
    fn the_engine_leaves_out_the_null_zone_and_a_zone_given_twice() {
        let g = grid("ISEA3H");
        let thirteen = g.sub_zones(zone(g, "B4-2-A"), 2).unwrap();
        let plain = the_engines(dggal_oracle::ISEA3H, &thirteen);
        assert_eq!(texts(g, &plain), THE_SEVEN);
        let mut with_null = thirteen.clone();
        with_null.insert(5, ZoneId(NULL_ZONE));
        with_null.push(ZoneId(NULL_ZONE));
        assert_eq!(the_engines(dggal_oracle::ISEA3H, &with_null), plain);
        assert!(matches!(
            g.compact_zones(&with_null),
            Err(Error::InvalidZone(_))
        ));
        assert_eq!(the_engines(dggal_oracle::ISEA3H, &[ZoneId(NULL_ZONE)]), []);
        let mut twice = thirteen.clone();
        twice.extend(&thirteen);
        assert_eq!(the_engines(dggal_oracle::ISEA3H, &twice), plain);
        assert_eq!(g.compact_zones(&twice).unwrap(), plain);
    }

    /// The engine does not look for the zone at the centre of a coarser one: the thirteen
    /// sub-zones of `B4-2-A` without the one at its centre compact, in the engine as here, to
    /// the seven zones that the thirteen compact to.
    #[test]
    fn the_engine_covers_a_set_that_lacks_the_zone_at_the_centre() {
        let g = grid("ISEA3H");
        let thirteen = g.sub_zones(zone(g, "B4-2-A"), 2).unwrap();
        let within = |z: ZoneId| g.neighbors(z).iter().all(|n| thirteen.contains(n));
        let twelve: Vec<ZoneId> = thirteen.iter().copied().filter(|&z| !within(z)).collect();
        assert_eq!(twelve.len(), 12);
        let theirs = the_engines(dggal_oracle::ISEA3H, &twelve);
        assert_eq!(texts(g, &theirs), THE_SEVEN);
        assert_eq!(g.compact_zones(&twelve).unwrap(), theirs);
    }

    /// The departure met live: given zones of two levels the engine drops the coarser zone
    /// that has no zone two levels up to give way to, or keeps both, where this crate refuses
    /// the set; and the engine's own answer, being of two levels, is not a set that it
    /// compacts to itself.
    #[test]
    fn the_engine_loses_zones_of_a_set_of_more_than_one_level() {
        let g = grid("ISEA3H");
        for (set, answer) in [
            ("A2-0-C C4-22-B", "C4-22-B"),
            ("A2-0-C B4-2-A", "B4-2-A"),
            ("B5-0-A A9-0-C A0-0-B A3-0-A", "A0-0-B A3-0-A B5-0-A"),
            ("A0-0-A B4-2-A", "A0-0-A B4-2-A"),
            // The engine goes on after a pass that takes nothing: the first pass, from level
            // 6, takes none, and the second takes `B4-2-A` for the six zones of level 4 about
            // its centre, on the evidence of the six of level 6 that the first pass kept.
            (
                "C2-1A-A C2-23-A C4-5-A C4-7-A C4-F-A C4-10-A D2-69-A D2-B7-A D2-15C-A D4-60-A \
                 D4-69-A D4-B7-A",
                "B4-2-A C2-1A-A C2-23-A C4-5-A C4-7-A C4-F-A C4-10-A D2-69-A D2-B7-A D2-15C-A \
                 D4-60-A D4-69-A D4-B7-A",
            ),
        ] {
            let set = zones(g, set);
            assert_eq!(texts(g, &the_engines(dggal_oracle::ISEA3H, &set)), answer);
            assert!(matches!(
                g.compact_zones(&set),
                Err(Error::MixedLevels { .. })
            ));
        }
        // Every zone of level 2 with every zone of level 0: the twelve of level 0.
        let mut twelve: Vec<ZoneId> = g.zones(0).unwrap().collect();
        twelve.sort();
        let mut set: Vec<ZoneId> = g.zones(2).unwrap().collect();
        set.extend(&twelve);
        assert_eq!(the_engines(dggal_oracle::ISEA3H, &set), twelve);
    }

    /// On the aperture-7 grids too the engine leaves its null zone out of the set and counts
    /// a zone given twice once; and it leaves out a Z7 identifier of twenty digits, which it
    /// reads as its null zone. This crate refuses the null zone and the twenty digits.
    #[test]
    fn on_aperture_7_the_engine_leaves_out_the_null_zone_and_twenty_digits() {
        for (name, engine) in ENGINES_7 {
            let g = grid(name);
            let thirteen = g.sub_zones(zone(g, "0064"), 1).unwrap();
            let plain = the_engines(engine, &thirteen);
            assert_eq!(texts(g, &plain), THE_THIRTEEN, "{name}");
            for left_out in [ZoneId(NULL_ZONE), z7(0, &[0; 20]), z7(3, &[1; 20])] {
                let mut with = thirteen.clone();
                with.insert(5, left_out);
                with.push(left_out);
                assert_eq!(the_engines(engine, &with), plain, "{name}");
                assert!(matches!(g.compact_zones(&with), Err(Error::InvalidZone(_))));
                assert_eq!(the_engines(engine, &[left_out]), [], "{name}");
            }
            let mut twice = thirteen.clone();
            twice.extend(&thirteen);
            assert_eq!(the_engines(engine, &twice), plain);
            assert_eq!(g.compact_zones(&twice).unwrap(), plain);
            // The child a pentagon does not have is read as one that it has.
            for (given, read_as) in [(z7(0, &[2]), "003"), (z7(6, &[5]), "064")] {
                assert_eq!(the_engines(engine, &[given]), [zone(g, read_as)], "{name}");
                assert_eq!(g.compact_zones(&[given]).unwrap(), [zone(g, read_as)]);
            }
        }
    }

    /// On the aperture-7 grids the engine takes a zone only where the set has every one of
    /// its children: the thirteen children of `0064` without any one of them come back as
    /// they are, in the engine as here.
    #[test]
    fn on_aperture_7_the_engine_does_not_cover_a_set_with_a_gap() {
        for (name, engine) in ENGINES_7 {
            let g = grid(name);
            let thirteen = g.sub_zones(zone(g, "0064"), 1).unwrap();
            for missing in &thirteen {
                let mut twelve: Vec<ZoneId> =
                    thirteen.iter().copied().filter(|z| z != missing).collect();
                let theirs = the_engines(engine, &twelve);
                assert_eq!(g.compact_zones(&twelve).unwrap(), theirs, "{name}");
                twelve.sort_by_key(|&z| order_key(g, z));
                assert_eq!(theirs, twelve, "{name}");
            }
        }
    }

    /// The departure met live on the aperture-7 grids: given zones of two levels the engine
    /// drops the zone of level 0, which has no parent to give way to, where this crate
    /// refuses the set.
    ///
    /// And a set of two levels reaches the short cut that the engine takes when its answer
    /// holds the seventy-two zones of level 1: of those zones, the ones whose digit is 1, 3
    /// or 4 as they are, and the children of the others. The engine answers ten zones of ten
    /// even levels under base cell 1 and its null zone twice, and not the twelve zones of
    /// level 0; the compaction of this crate, asked beneath the refusal, answers the ten, which
    /// the topology's own unit test of the compaction holds to the same ten identifiers.
    #[test]
    fn on_aperture_7_the_engine_loses_zones_of_a_set_of_more_than_one_level() {
        for (name, engine) in ENGINES_7 {
            let g = grid(name);
            for (set, answer) in [
                ("00 0064", "0064"),
                ("05 011", "011"),
                ("0064 00641 00", "0064 00641"),
                // The engine goes on after a pass that takes nothing: the first pass takes
                // `00640` for its thirteen children, the second takes none, and the third
                // takes `0064`, every child of which is then among the zones kept.
                (
                    "00641 00645 00644 00646 00642 00643 00665 00604 00656 00422 03332 03323 \
                     006400 006401 006405 006404 006406 006402 006403 006461 006425 006434 \
                     006416 006452 006443",
                    "00422 03332 03323 0064 00640 00641 00645 00644 00646 00642 00643 00656 \
                     00665 00604 006416 006452 006405 006404 006443 006406 006401 006461 \
                     006402 006403 006434 006425",
                ),
            ] {
                let set = zones(g, set);
                assert_eq!(texts(g, &the_engines(engine, &set)), answer, "{name}");
                assert!(matches!(
                    g.compact_zones(&set),
                    Err(Error::MixedLevels { .. })
                ));
            }

            let mut set: Vec<ZoneId> = Vec::new();
            for z in g.zones(1).unwrap() {
                if matches!(g.text_id(z).as_bytes()[2], b'1' | b'3' | b'4') {
                    set.push(z);
                } else {
                    set.extend(g.children(z));
                }
            }
            assert!(matches!(
                g.compact_zones(&set),
                Err(Error::MixedLevels { level: 2, other: 1 })
            ));
            let theirs = the_engines(engine, &set);
            assert_eq!(theirs.len(), 12, "{name}");
            assert_eq!(
                texts(g, &theirs[..10]),
                "01 0100 010000 01000000 0100000000 010000000000 01000000000000 \
                 0100000000000000 010000000000000000 01000000000000000000",
                "{name}"
            );
            assert_eq!(theirs[10..], [ZoneId(NULL_ZONE); 2], "{name}");
        }
    }

    /// The literal answers of the broken seam, as the engine gives them on the three grids;
    /// and the one identifier, among those of the test without the engine, that the engine
    /// answers and cannot itself read back.
    #[test]
    fn in_a_broken_seam_the_engine_answers_what_is_pinned() {
        for (name, engine) in ENGINES_7 {
            let g = grid(name);
            let theirs = the_engines(engine, &zones(g, SEAM_ELEVEN));
            assert_eq!(texts(g, &theirs), SEAM_NINE, "{name}");
            for (set, compact) in SEAM_WITH_NULL_CHILDREN {
                let theirs = the_engines(engine, &zones(g, set));
                assert_eq!(texts(g, &theirs), compact, "{name}");
            }
            assert_eq!(
                the_engines(engine, &[zone(g, "00000000000000001311")]),
                [],
                "{name}"
            );
            assert_eq!(
                the_engines(engine, &[zone(g, "0000000000000001644")]),
                [zone(g, "0000000000000000056")],
                "{name}"
            );
            let theirs = the_engines(engine, &[zone(g, "0100000000000001623")]);
            assert_eq!(theirs.len(), 1, "{name}");
            assert_eq!(
                dggal_oracle::text_id(engine, theirs[0].0),
                "0122222222222222203"
            );
            assert!(
                !dggal_oracle::engine_can_read(engine, theirs[0].0),
                "{name}"
            );
        }
    }

    /// What a run over the sets of a broken seam met.
    #[derive(Default, Debug, PartialEq)]
    struct SeamTally {
        /// Sets compared with the engine, and those of them that its answer makes smaller.
        sets: usize,
        smaller: usize,
        /// Sets whose answer is the engine's, zone for zone.
        equal: usize,
        /// Sets whose answer is the engine's less the identifiers in it that the engine
        /// cannot itself read back, which this crate never hands out.
        less_unreadable: usize,
    }

    /// This crate's compaction of a set of a broken seam against the engine's: the same
    /// sequence, or the same less the identifiers of the engine's answer that the engine
    /// cannot read back; anything else fails.
    fn is_the_engines_in_a_seam(
        g: AnyGrid,
        engine: &str,
        set: &[ZoneId],
        tally: &mut SeamTally,
        what: &str,
    ) {
        let ours = g.compact_zones(set).unwrap();
        let theirs = the_engines(engine, set);
        tally.sets += 1;
        tally.smaller += usize::from(theirs.len() < set.len());
        if ours == theirs {
            tally.equal += 1;
            return;
        }
        let readable: Vec<ZoneId> = theirs
            .iter()
            .copied()
            .filter(|z| dggal_oracle::engine_can_read(engine, z.0))
            .collect();
        assert!(
            ours == readable,
            "{what}: {} where the engine has {}",
            texts(g, &ours),
            texts(g, &theirs)
        );
        tally.less_unreadable += 1;
    }

    /// The Z7 identifier of `level` digits under `base`: noughts, and then the four digits of
    /// `tail` in base seven. Beside the pole of base cell 0 and along the edges from base cell
    /// 0 to 1 and from 1 to 6 such identifiers lie in the two broken seams from level 15.
    fn seam_zone(base: u64, level: usize, tail: usize) -> ZoneId {
        let mut digits = vec![0u8; level];
        let mut rest = tail;
        for d in digits.iter_mut().rev().take(4) {
            *d = (rest % 7) as u8;
            rest /= 7;
        }
        z7(base, &digits)
    }

    /// Sets of level 16 on a broken seam, where the engine's lists of children hold its null
    /// zone and identifiers it cannot read back, a zone may have no parent and the first
    /// parent need not be the one that the digits name. Each set is the children of a zone
    /// of level 15 and of its neighbours, some eighty zones: about zones named by their
    /// digits under base cell 0, and about the zones found at the centroids of those. Every
    /// answer is the engine's, zone for zone: no identifier that the engine cannot read back
    /// was met in the answer for a set made so.
    #[test]
    fn in_a_broken_seam_sets_of_level_16_are_compacted_as_the_engine_does() {
        for (name, engine) in ENGINES_7 {
            let g = grid(name);
            let mut tally = SeamTally::default();
            // One identifier in forty-one of the 2,401.
            for tail in (0..2_401).step_by(41) {
                let named = seam_zone(0, 15, tail);
                let c = g.centroid(named);
                let found = g.zone_from_geo(c.lat, c.lon, 15).unwrap();
                for about in [named, found] {
                    if about == ZoneId::NULL {
                        continue;
                    }
                    let mut set: Vec<ZoneId> = Vec::new();
                    for z in g.disk(about, 1) {
                        set.extend(g.children(z));
                    }
                    assert!(set.iter().all(|&z| g.resolution(z) == 16));
                    let what = format!("{name}, the children about {}", g.text_id(about));
                    is_the_engines_in_a_seam(g, engine, &set, &mut tally, &what);
                }
            }
            eprintln!("{name}: {tally:?}");
            assert!(tally.sets >= 100, "{name}: {tally:?}");
            assert_eq!(tally.equal, tally.sets, "{name}: {tally:?}");
            assert!(tally.smaller >= 100, "{name}: {tally:?}");
        }
    }

    /// The one departure from the engine's answer, counted where it is met: single zones of
    /// levels 17 and 19 named by their digits under base cell 1, beside the edge from base
    /// cell 0 to 1. For some of them the engine's answer is an identifier that the engine
    /// cannot itself read back, and this crate's answer is empty; for every other the answer
    /// is the engine's: the zone, another that the engine reads it as, or none.
    #[test]
    fn in_a_broken_seam_an_identifier_the_engine_cannot_read_back_is_left_out() {
        for (name, engine) in ENGINES_7 {
            let g = grid(name);
            let mut tally = SeamTally::default();
            for level in [17, 19] {
                for tail in (0..2_401).step_by(13) {
                    let named = seam_zone(1, level, tail);
                    let what = format!("{name}, {}", g.text_id(named));
                    is_the_engines_in_a_seam(g, engine, &[named], &mut tally, &what);
                }
            }
            eprintln!("{name}: {tally:?}");
            assert_eq!(tally.sets, 2 * 185, "{name}");
            assert_eq!(tally.equal + tally.less_unreadable, tally.sets, "{name}");
            assert!(tally.less_unreadable >= 10, "{name}: {tally:?}");
        }
    }
}
