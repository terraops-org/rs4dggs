//! The zones of a bounding box, their order, the estimate of their number and the compaction of
//! a set of zones, against the live engine, for the grid under test.
//!
//! The tests of each of these (`tests/zones.rs`, `tests/zones_in_box.rs`,
//! `tests/estimate_zones_in_box.rs`, `tests/compact_zones.rs`) hold the answers themselves, on
//! thousands of boxes and sets. What is held here is what they leave: the same kinds of box at
//! EVERY level of every grid, each put to the engine or counted as kept from it; the one place
//! in which the answer to a box departs from the engine's, classified against the live engine
//! at every instance and counted by kind of box; the levels next beyond those that the tests
//! of the answer take whole; the margin of the test of two extents, shown on the engine itself;
//! the compaction at every level, and on the suites' own sample of zones; and what hostile
//! input does. Every count is asserted exactly, from a table for each grid: the boxes are
//! seeded, the engine answers the same every time, and a grid without a table fails.
//!
//! **The departure.** `zones_in_box` yields every zone of the level whose extent meets the box
//! by the engine's own test of two extents, where the engine's `listZones` leaves some of them
//! out. Wherever the two answers differ, `expected_divergences::engine_omits_zones_of_a_box`
//! must hold: the engine's answer is this crate's less some zones, in the same order, and each
//! zone it leaves out meets the box by the engine's own extent of it. Never the reverse: a zone
//! of the engine's that this crate lacks fails the test.
//!
//! **The boxes kept from the engine.** The engine does not survive every box, and for some it
//! costs more than a test may spend (`boxes::engine_survives`). Such a box is not put to it,
//! which is the one exception here to the rule that no oracle call is skipped; it is counted,
//! by kind, and this crate's answer is held there by what can be asked of the engine zone by
//! zone: every zone returned meets the box by the engine's own extent of it, no neighbour of a
//! returned zone that meets the box is missing, the zone at the centre of the box is there,
//! and the answer taken in pages is the same answer, with nothing after its last zone.
use std::collections::{BTreeSet, HashSet};

use dggal_oracle as o;
use rs4dggs::{Error, Extent, GeoPoint, Indexing, Topology, ZoneId, get_grid};

use super::boxes::{self, KINDS, Kind};
use super::expected_divergences as div;
use super::geometry::hostile_identifiers;
use super::subject::Subject;
use super::suite::{SEED, zones as sample};

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

/// The first `most` zones of `level` whose extent meets the box.
fn first_in_box<S: Subject>(level: u8, b: &[f64; 4], most: usize) -> Vec<ZoneId> {
    S::grid()
        .zones_in_box(level, &degrees(b))
        .unwrap()
        .take(most)
        .collect()
}

/// The answer taken in pages of `size`, each entered after the last zone of the one before,
/// until a page comes back empty or `most` zones are had; cut to `most`.
fn paged<S: Subject>(level: u8, b: &[f64; 4], size: usize, most: usize) -> Vec<ZoneId> {
    let mut out = first_in_box::<S>(level, b, size);
    while let Some(&last) = out.last() {
        if out.len() >= most {
            break;
        }
        let page: Vec<ZoneId> = S::grid()
            .zones_in_box(level, &degrees(b))
            .unwrap()
            .after(last)
            .unwrap()
            .take(size)
            .collect();
        if page.is_empty() {
            break;
        }
        out.extend(page);
    }
    out.truncate(most);
    out
}

/// The zones are of `level`, in the order of the level, none twice: their keys in the lattice
/// ascend.
fn are_in_the_order_of_the_level<S: Subject>(level: u8, zones: &[ZoneId], what: &str) {
    let g = S::grid();
    let key = |z: ZoneId| -> u64 {
        assert_eq!(g.resolution(z), level, "{what}: {}", g.text_id(z));
        S::T::locate(&S::I::decode(z))
            .unwrap_or_else(|| panic!("{what}: {} is no zone of the lattice", g.text_id(z)))
            .3
    };
    let keys: Vec<u64> = zones.iter().map(|&z| key(z)).collect();
    assert!(
        keys.windows(2).all(|w| w[0] < w[1]),
        "{what}: the keys do not ascend"
    );
}

/// The engine's extent of a zone, in radians.
fn engines_extent<S: Subject>(zone: ZoneId) -> [f64; 4] {
    o::extent(S::ORACLE, zone.0)
}

/// The zone at the centre of a box, in degrees, whose height and whose width both exceed 1e-9
/// radian: a zone that every answer to the box must hold, its extent holding the centre, about
/// which the box reaches half of that on every side. `None` for a box without such a height
/// and width (a point on a line at which the extents of the zones on either side end yields no
/// zone), and where the centre has no zone.
fn the_zone_at_the_centre<S: Subject>(level: u8, b: &[f64; 4]) -> Option<ZoneId> {
    let [s, w, n, e] = boxes::radians(b);
    let width = if e < w {
        e - w + 2.0 * std::f64::consts::PI
    } else {
        e - w
    };
    if n - s <= 1e-9 || width <= 1e-9 {
        return None;
    }
    let span = if b[3] < b[1] {
        b[3] - b[1] + 360.0
    } else {
        b[3] - b[1]
    };
    let lon = b[1] + span / 2.0;
    let lon = if lon > 180.0 { lon - 360.0 } else { lon };
    let at = S::grid()
        .zone_from_geo((b[0] + b[2]) / 2.0, lon, level)
        .unwrap()
        .id();
    (at != ZoneId::NULL).then_some(at)
}

/// The most zones taken of the answer to a box that is not put to the engine and is too large
/// to take whole: one of more than [`WHOLE`] zones' area.
const MOST: usize = 300;

/// The greatest area of a box that is not put to the engine, in zones of the level, at which
/// its answer is taken whole all the same. A box that is put to the engine has its answer
/// taken whole whatever its area, which `boxes::engine_survives` keeps to
/// `boxes::MAX_ZONES_IN_BOX`.
const WHOLE: f64 = 5_000.0;

/// The most zones of an answer that is also compacted, whole and with one zone in twenty
/// taken out, and the compaction compared with the engine's.
const COMPACTED: usize = 150;

/// What the run over the boxes of a grid met. Every count is exact.
#[derive(Debug, Default, PartialEq, Eq)]
struct Boxes {
    /// Boxes made, of each kind, in the order of `boxes::KINDS`: `rand`, `tiny`, `big`,
    /// `anti`, `pole`, `vertex`, `point`, `line`, `corner`.
    made: [usize; 9],
    /// Boxes put to the engine's `listZones`.
    asked: [usize; 9],
    /// Boxes kept from it: those it does not survive, and those that would cost it more than
    /// a test may spend.
    kept: [usize; 9],
    /// Boxes, of those put to the engine, in which this crate's answer holds zones that the
    /// engine's lacks: the departure.
    larger: [usize; 9],
    /// Of those, the boxes for which the engine returns no zone at all.
    engines_empty: [usize; 9],
    /// The zones that the engine left out, over those boxes.
    omitted: usize,
    /// The zones this crate returned, over every box, and the zones the engine returned.
    zones: usize,
    engines: usize,
    /// Answers taken whole, each of which was held against its estimate.
    whole: usize,
    /// Of the boxes kept from the engine, those whose answer was held to the zone at the
    /// centre of the box.
    centred: usize,
    /// Answers that were compacted, and those of twenty zones or more compacted again with one
    /// zone in twenty taken out.
    compacted: usize,
    thinned: usize,
    /// Of both, the sets that the compaction made smaller; the zones in the sets and in the
    /// answers; and the thinned sets that the answer covers beyond themselves.
    smaller: usize,
    zones_in: usize,
    zones_out: usize,
    covered_beyond: usize,
}

/// What the boxes of each grid meet, measured.
fn boxes_recorded<S: Subject>() -> Boxes {
    match S::NAME {
        "IGEO7" => Boxes {
            made: [15, 15, 15, 15, 15, 15, 15, 15, 15],
            asked: [15, 15, 5, 15, 15, 15, 15, 15, 15],
            kept: [0, 0, 10, 0, 0, 0, 0, 0, 0],
            larger: [0, 0, 0, 0, 0, 0, 0, 0, 0],
            engines_empty: [0, 0, 0, 0, 0, 0, 0, 0, 0],
            omitted: 0,
            zones: 9288,
            engines: 6288,
            whole: 125,
            centred: 0,
            compacted: 97,
            thinned: 23,
            smaller: 10,
            zones_in: 3035,
            zones_out: 2953,
            covered_beyond: 0,
        },
        "IVEA7H" => Boxes {
            made: [15, 15, 15, 15, 15, 15, 15, 15, 15],
            asked: [15, 15, 5, 15, 15, 15, 15, 15, 15],
            kept: [0, 0, 10, 0, 0, 0, 0, 0, 0],
            larger: [0, 0, 0, 0, 0, 0, 0, 0, 0],
            engines_empty: [0, 0, 0, 0, 0, 0, 0, 0, 0],
            omitted: 0,
            zones: 9297,
            engines: 6297,
            whole: 125,
            centred: 0,
            compacted: 94,
            thinned: 25,
            smaller: 9,
            zones_in: 3096,
            zones_out: 3022,
            covered_beyond: 0,
        },
        "RTEA7H" => Boxes {
            made: [15, 15, 15, 15, 15, 15, 15, 15, 15],
            asked: [15, 15, 5, 15, 15, 15, 15, 15, 15],
            kept: [0, 0, 10, 0, 0, 0, 0, 0, 0],
            larger: [0, 0, 0, 0, 0, 0, 0, 0, 0],
            engines_empty: [0, 0, 0, 0, 0, 0, 0, 0, 0],
            omitted: 0,
            zones: 9304,
            engines: 6304,
            whole: 125,
            centred: 0,
            compacted: 92,
            thinned: 24,
            smaller: 11,
            zones_in: 3079,
            zones_out: 2993,
            covered_beyond: 0,
        },
        "ISEA3H" => Boxes {
            made: [68, 68, 68, 68, 68, 68, 68, 68, 68],
            asked: [67, 60, 14, 68, 26, 55, 68, 68, 68],
            kept: [1, 8, 54, 0, 42, 13, 0, 0, 0],
            larger: [4, 1, 0, 1, 18, 5, 27, 5, 6],
            engines_empty: [0, 0, 0, 0, 13, 0, 0, 0, 0],
            omitted: 101,
            zones: 32411,
            engines: 12693,
            whole: 559,
            centred: 57,
            compacted: 224,
            thinned: 53,
            smaller: 38,
            zones_in: 6076,
            zones_out: 5549,
            covered_beyond: 1,
        },
        "IVEA3H" => Boxes {
            made: [68, 68, 68, 68, 68, 68, 68, 68, 68],
            asked: [67, 60, 14, 68, 26, 55, 68, 68, 68],
            kept: [1, 8, 54, 0, 42, 13, 0, 0, 0],
            larger: [5, 2, 0, 1, 18, 2, 31, 2, 5],
            engines_empty: [0, 0, 0, 0, 13, 0, 0, 0, 0],
            omitted: 101,
            zones: 32392,
            engines: 12694,
            whole: 559,
            centred: 57,
            compacted: 225,
            thinned: 54,
            smaller: 40,
            zones_in: 6099,
            zones_out: 5542,
            covered_beyond: 1,
        },
        "RTEA3H" => Boxes {
            made: [68, 68, 68, 68, 68, 68, 68, 68, 68],
            asked: [67, 60, 14, 68, 26, 55, 68, 68, 68],
            kept: [1, 8, 54, 0, 42, 13, 0, 0, 0],
            larger: [4, 1, 0, 2, 18, 1, 25, 5, 3],
            engines_empty: [0, 0, 0, 0, 13, 0, 0, 0, 0],
            omitted: 90,
            zones: 32350,
            engines: 12653,
            whole: 559,
            centred: 57,
            compacted: 215,
            thinned: 54,
            smaller: 38,
            zones_in: 6018,
            zones_out: 5485,
            covered_beyond: 3,
        },
        other => panic!("nothing has been measured of the boxes of {other}"),
    }
}

/// This crate's compaction of a set of zones of `level` against the engine's, sequence for
/// sequence; and that the sub-zones, at the level of the set, of the zones of the answer hold
/// every zone of the set. Returns the length of the answer and the number of zones that those
/// sub-zones hold beyond the set.
fn is_compacted_as_the_engine_compacts<S: Subject>(
    level: u8,
    set: &[ZoneId],
    what: &str,
) -> (usize, usize) {
    let g = S::grid();
    let ours = g
        .compact_zones(set)
        .unwrap_or_else(|e| panic!("{what}: the set is refused: {e}"));
    let ids: Vec<u64> = set.iter().map(|z| z.0).collect();
    let theirs: Vec<ZoneId> = o::compact_zones(S::ORACLE, &ids)
        .into_iter()
        .map(ZoneId)
        .collect();
    assert!(
        ours == theirs,
        "{what}: the set of {} zones compacts to {} where the engine has {}; the first to \
         differ {:?}",
        set.len(),
        ours.len(),
        theirs.len(),
        ours.iter()
            .zip(&theirs)
            .find(|(a, b)| a != b)
            .map(|(a, b)| (g.text_id(*a), g.text_id(*b)))
    );
    let mut expanded = BTreeSet::new();
    for &z in &ours {
        let of = g.resolution(z);
        assert!(of <= level, "{what}: {} is of level {of}", g.text_id(z));
        expanded.extend(g.sub_zones(z, level - of).unwrap());
    }
    let given: BTreeSet<ZoneId> = set.iter().copied().collect();
    assert!(
        expanded.is_superset(&given),
        "{what}: a zone of the set is not among the sub-zones of the answer"
    );
    (ours.len(), expanded.len() - given.len())
}

/// Boxes of the nine kinds at every level at which the grid answers a box, two of each kind
/// at each level on aperture 3 and one on aperture 7, each put to the engine where it survives
/// it.
///
/// - A box put to the engine: the two answers are equal, or differ as
///   `expected_divergences::engine_omits_zones_of_a_box` allows, which re-checks its evidence
///   against the engine there and then.
/// - A box kept from the engine: every zone returned meets the box by the engine's own extent
///   of it; no neighbour of a returned zone that meets the box is missing; the zone at the
///   centre of the box is among them; and the answer in pages of 17 is the same answer, with
///   nothing after its last zone. An answer too large to take whole is held to its first
///   [`MOST`] zones.
/// - Every answer is in the order of the level. Every whole answer is at most its estimate.
/// - Every whole answer of two to [`COMPACTED`] zones, of the first box of each kind, is
///   compacted as the engine compacts it, and so is a set of twenty zones or more with one
///   zone in twenty taken out. A whole answer is a set without a gap, and its compaction
///   covers it and no more; on aperture 7 the same holds of the thinned set, and on aperture 3
///   the thinned sets covered beyond themselves are counted.
///
/// The counts, by kind of box, are those recorded for the grid, exactly.
pub fn the_zones_of_a_box_are_the_engines_and_those_it_leaves_out<S: Subject>() {
    let g = S::grid();
    let seven = S::T::APERTURE == 7;
    let a_kind = if seven { 1 } else { 2 };
    let mut met = Boxes::default();
    for level in 0..=g.max_box_level() {
        let mut rng = boxes::Rng::new(SEED ^ (0xb0c5 << 8) ^ u64::from(level));
        let sites = boxes::sites(S::ORACLE, i32::from(level), &mut rng);
        let size = boxes::zone_width_deg(S::ORACLE, i32::from(level));
        for (k, kind) in KINDS.into_iter().enumerate() {
            for i in 0..a_kind {
                let b = boxes::make_box(kind, &mut rng, size, &sites);
                let r = boxes::radians(&b);
                let what = format!(
                    "{} level {level}, the {} box {b:?} (seed {SEED:#x})",
                    S::NAME,
                    kind.name()
                );
                met.made[k] += 1;
                let asked = boxes::engine_survives(S::ORACLE, i32::from(level), &b);
                let whole = asked || boxes::area_in_zones(S::ORACLE, i32::from(level), &b) <= WHOLE;
                let most = if whole { usize::MAX } else { MOST };
                let ours = first_in_box::<S>(level, &b, most);
                met.zones += ours.len();
                are_in_the_order_of_the_level::<S>(level, &ours, &what);
                if whole {
                    let estimate = g.estimate_zones_in_box(level, &degrees(&b)).unwrap();
                    assert!(
                        estimate >= ours.len() as u64,
                        "{what}: {} zones, estimated at {estimate}",
                        ours.len()
                    );
                    met.whole += 1;
                }

                if asked {
                    let theirs = o::list_zones(S::ORACLE, i32::from(level), Some(r));
                    met.asked[k] += 1;
                    met.engines += theirs.len();
                    if theirs.len() != ours.len()
                        || ours.iter().zip(&theirs).any(|(a, b)| a.0 != *b)
                    {
                        let ids: Vec<u64> = ours.iter().map(|z| z.0).collect();
                        match div::engine_omits_zones_of_a_box::<S>(level, &b, &ids, &theirs) {
                            Ok(left_out) => {
                                met.larger[k] += 1;
                                met.omitted += left_out;
                                met.engines_empty[k] += usize::from(theirs.is_empty());
                            }
                            Err(why) => panic!("{what}: {why}"),
                        }
                    }
                } else {
                    met.kept[k] += 1;
                    for &z in &ours {
                        assert!(
                            boxes::meets(&engines_extent::<S>(z), &r),
                            "{what}: {} is returned and its extent does not meet the box",
                            g.text_id(z)
                        );
                    }
                    let most = if whole { ours.len() + 1 } else { MOST };
                    assert!(
                        paged::<S>(level, &b, 17, most) == ours,
                        "{what}: the pages of 17 are not the answer"
                    );
                    if whole {
                        if let Some(at) = the_zone_at_the_centre::<S>(level, &b) {
                            assert!(
                                ours.contains(&at),
                                "{what}: {}, the zone at the centre of the box, is missing",
                                g.text_id(at)
                            );
                            met.centred += 1;
                        }
                        let held: HashSet<u64> = ours.iter().map(|z| z.0).collect();
                        for &z in &ours {
                            for n in o::neighbors(S::ORACLE, z.0) {
                                assert!(
                                    held.contains(&n)
                                        || !boxes::meets(&engines_extent::<S>(ZoneId(n)), &r),
                                    "{what}: {}, a neighbour of {}, meets the box and is missing",
                                    g.text_id(ZoneId(n)),
                                    g.text_id(z)
                                );
                            }
                        }
                    } else {
                        assert_eq!(ours.len(), MOST, "{what}");
                    }
                }

                if whole && (2..=COMPACTED).contains(&ours.len()) && i == 0 {
                    let (out, beyond) =
                        is_compacted_as_the_engine_compacts::<S>(level, &ours, &what);
                    assert_eq!(
                        beyond, 0,
                        "{what}: the whole answer is covered beyond itself"
                    );
                    met.compacted += 1;
                    met.smaller += usize::from(out < ours.len());
                    met.zones_in += ours.len();
                    met.zones_out += out;
                    if ours.len() >= 20 {
                        let thinned: Vec<ZoneId> = ours
                            .iter()
                            .enumerate()
                            .filter(|(i, _)| i % 20 != 7)
                            .map(|(_, &z)| z)
                            .collect();
                        let what = format!("{what}, one zone in twenty taken out");
                        let (out, beyond) =
                            is_compacted_as_the_engine_compacts::<S>(level, &thinned, &what);
                        assert!(
                            !seven || beyond == 0,
                            "{what}: the set is covered beyond itself"
                        );
                        met.thinned += 1;
                        met.smaller += usize::from(out < thinned.len());
                        met.zones_in += thinned.len();
                        met.zones_out += out;
                        met.covered_beyond += usize::from(beyond > 0);
                    }
                }
            }
        }
    }
    eprintln!("the zones of a box on {}: {met:?}", S::NAME);
    assert_eq!(
        met,
        boxes_recorded::<S>(),
        "the boxes of {}, by kind ({:?}), are not as recorded",
        S::NAME,
        KINDS.map(Kind::name)
    );
}

/// The box about Lisbon, and a small cap about each pole.
const LISBON: [f64; 4] = [38.0, -10.0, 39.0, -9.0];
const NORTHERN_CAP: [f64; 4] = [89.5, -180.0, 90.0, 180.0];
const SOUTHERN_CAP: [f64; 4] = [-90.0, -180.0, -89.5, 180.0];

/// The boxes known by name, each with its level: the box about Lisbon at levels 2, 3 and 4,
/// where it lies across the side of two rhombi, and the two polar caps at level 6, where each
/// is smaller than a zone.
const NAMED: [(&str, u8, [f64; 4]); 5] = [
    ("the box about Lisbon", 2, LISBON),
    ("the box about Lisbon", 3, LISBON),
    ("the box about Lisbon", 4, LISBON),
    ("the cap about the north pole", 6, NORTHERN_CAP),
    ("the cap about the south pole", 6, SOUTHERN_CAP),
];

/// How many zones this crate returns for each of the boxes of [`NAMED`], and how many the
/// engine, measured.
fn named_recorded<S: Subject>() -> [(usize, usize); 5] {
    match S::NAME {
        "IGEO7" => [(3, 3), (3, 3), (3, 3), (38, 38), (38, 38)],
        "IVEA7H" => [(3, 3), (3, 3), (3, 3), (38, 38), (38, 38)],
        "RTEA7H" => [(3, 3), (3, 3), (3, 3), (38, 38), (38, 38)],
        "ISEA3H" => [(1, 0), (2, 1), (2, 1), (2, 0), (2, 0)],
        "IVEA3H" => [(1, 0), (2, 1), (2, 1), (2, 0), (2, 0)],
        "RTEA3H" => [(1, 0), (2, 1), (2, 1), (2, 0), (2, 0)],
        other => panic!("nothing has been measured of the named boxes of {other}"),
    }
}

/// The boxes known by name, on the grid under test: the box about Lisbon at levels 2, 3 and 4
/// and a cap of half a degree about each pole at level 6. On the aperture-3 grids each is a
/// box for which the engine returns fewer zones than have an extent that meets it, and none
/// at all for the two caps; on the aperture-7 grids the two answers are one. The numbers of
/// zones are those recorded, exactly, and each difference is one that
/// `expected_divergences::engine_omits_zones_of_a_box` allows.
pub fn the_boxes_known_by_name_depart_from_the_engine_as_recorded<S: Subject>() {
    let mut met = [(0, 0); 5];
    for ((name, level, b), met) in NAMED.into_iter().zip(&mut met) {
        let what = format!("{}, {name} at level {level}", S::NAME);
        assert!(
            boxes::engine_survives(S::ORACLE, i32::from(level), &b),
            "{what}: not a box to put to the engine"
        );
        let ours = first_in_box::<S>(level, &b, usize::MAX);
        are_in_the_order_of_the_level::<S>(level, &ours, &what);
        let ids: Vec<u64> = ours.iter().map(|z| z.0).collect();
        let theirs = o::list_zones(S::ORACLE, i32::from(level), Some(boxes::radians(&b)));
        *met = (ours.len(), theirs.len());
        if ids != theirs {
            let left_out = div::engine_omits_zones_of_a_box::<S>(level, &b, &ids, &theirs)
                .unwrap_or_else(|why| panic!("{what}: {why}"));
            assert_eq!(left_out, ours.len() - theirs.len(), "{what}");
        }
    }
    eprintln!("the boxes known by name on {}: {met:?}", S::NAME);
    assert_eq!(
        met,
        named_recorded::<S>(),
        "{}: the zones of this crate and of the engine for the boxes known by name",
        S::NAME
    );
}

/// The levels at which every zone of the grid is put to the rule, those next finer than the
/// levels at which the tests of the answer do so: levels 7 and 8 on aperture 3, of 21,872 and
/// 65,612 zones, and level 5 on aperture 7, of 168,072, each zone with the engine's own extent.
fn levels_taken_whole<S: Subject>() -> &'static [u8] {
    if S::T::APERTURE == 7 { &[5] } else { &[7, 8] }
}

/// How many boxes of the levels taken whole are put to the rule, and how many zones their
/// answers hold, measured.
fn whole_levels_recorded<S: Subject>() -> (usize, usize) {
    match S::NAME {
        "IGEO7" => (50, 31091),
        "IVEA7H" => (50, 31091),
        "RTEA7H" => (50, 31085),
        "ISEA3H" => (100, 30120),
        "IVEA3H" => (100, 30116),
        "RTEA3H" => (100, 30113),
        other => panic!("nothing has been measured of the levels taken whole of {other}"),
    }
}

/// At the levels of [`levels_taken_whole`], every zone of the grid is put to the rule. The
/// engine lists its whole world there, which is this crate's `zones` of the level, sequence
/// for sequence; each zone is given the engine's own extent; and for boxes of every kind, six
/// of each and two of the large ones, the answer of `zones_in_box` is the zones of that list
/// whose extent meets the box, the same zones in the same order. Nothing is left to the search
/// here: a zone wrongly left out is seen, as it is not where the engine itself leaves zones
/// out.
pub fn at_a_level_taken_whole_the_zones_of_a_box_are_those_of_the_rule<S: Subject>() {
    let g = S::grid();
    let (mut made, mut zones) = (0, 0);
    for &level in levels_taken_whole::<S>() {
        let listed = o::list_zones(S::ORACLE, i32::from(level), None);
        assert_eq!(
            listed.len() as u64,
            o::count_zones(S::ORACLE, i32::from(level))
        );
        let every: Vec<u64> = g.zones(level).unwrap().map(|z| z.0).collect();
        assert!(
            every == listed,
            "{} level {level}: every zone of the level is not the engine's list of the world",
            S::NAME
        );
        let world: Vec<(u64, [f64; 4])> = listed
            .into_iter()
            .map(|z| (z, o::extent(S::ORACLE, z)))
            .collect();

        let mut rng = boxes::Rng::new(SEED ^ (0x3a11 << 8) ^ u64::from(level));
        let sites = boxes::sites(S::ORACLE, i32::from(level), &mut rng);
        let size = boxes::zone_width_deg(S::ORACLE, i32::from(level));
        for kind in KINDS {
            for _ in 0..if kind == Kind::Big { 2 } else { 6 } {
                let b = boxes::make_box(kind, &mut rng, size, &sites);
                let r = boxes::radians(&b);
                let expected: Vec<u64> = world
                    .iter()
                    .filter(|(_, e)| boxes::meets(e, &r))
                    .map(|(z, _)| *z)
                    .collect();
                let ours: Vec<u64> = first_in_box::<S>(level, &b, usize::MAX)
                    .into_iter()
                    .map(|z| z.0)
                    .collect();
                assert!(
                    ours == expected,
                    "{} level {level}, the {} box {b:?} (seed {SEED:#x}): {} zones where the \
                     rule has {}; the first to differ {:?}",
                    S::NAME,
                    kind.name(),
                    ours.len(),
                    expected.len(),
                    ours.iter().zip(&expected).find(|(a, b)| a != b)
                );
                made += 1;
                zones += ours.len();
            }
        }
    }
    eprintln!(
        "{}, the levels taken whole: {made} boxes, {zones} zones",
        S::NAME
    );
    assert_eq!(
        (made, zones),
        whole_levels_recorded::<S>(),
        "{}: the boxes put to the rule, and the zones of their answers",
        S::NAME
    );
}

/// The cap above 70 degrees north, a box of tens of thousands of zones at the level of
/// [`large_box_level`].
const LARGE_BOX: [f64; 4] = [70.0, -180.0, 90.0, 180.0];

/// The level at which the cap above 70 degrees north is counted.
fn large_box_level<S: Subject>() -> u8 {
    if S::T::APERTURE == 7 { 6 } else { 11 }
}

/// The number of zones of the cap above 70 degrees north at that level, measured.
fn large_box_recorded<S: Subject>() -> u64 {
    match S::NAME {
        "IGEO7" => 36180,
        "IVEA7H" => 36186,
        "RTEA7H" => 36170,
        "ISEA3H" => 54362,
        "IVEA3H" => 54356,
        "RTEA3H" => 54330,
        other => panic!("nothing has been measured of the large box of {other}"),
    }
}

/// The estimate of a large box is not under its count, and is little above it: the cap above
/// 70 degrees north, tens of thousands of zones, on the grid under test. The engine is not
/// asked, the box being far beyond what a test may spend on it; the count is this crate's
/// own, and is the number recorded. On a box of this size the estimate stands some two parts
/// in a hundred above the count at the most.
pub fn the_estimate_of_a_large_box_is_not_under_its_count<S: Subject>() {
    let g = S::grid();
    let level = large_box_level::<S>();
    assert!(!boxes::engine_survives(
        S::ORACLE,
        i32::from(level),
        &LARGE_BOX
    ));
    let cap = degrees(&LARGE_BOX);
    let counted = g.zones_in_box(level, &cap).unwrap().count() as u64;
    let estimated = g.estimate_zones_in_box(level, &cap).unwrap();
    eprintln!(
        "{} level {level}, the cap above 70 degrees north: {counted} zones, estimated at \
         {estimated}",
        S::NAME
    );
    assert_eq!(
        counted,
        large_box_recorded::<S>(),
        "{} level {level}",
        S::NAME
    );
    assert!(
        estimated >= counted && estimated as f64 <= 1.025 * counted as f64,
        "{} level {level}: {counted} zones, estimated at {estimated}",
        S::NAME
    );
}

/// How many zones of the suites' sample give a set to compact: every `SETS_STRIDE`-th zone of
/// it that has zones below it.
const SETS_STRIDE: usize = 12;

/// The sets made of the suites' sample, those in whose answer a zone coarser than the set
/// stands, and the zones in the sets and in the answers, measured.
fn sample_sets_recorded<S: Subject>() -> [usize; 4] {
    match S::NAME {
        "IGEO7" => [433, 433, 5591, 5586],
        "IVEA7H" => [439, 439, 5675, 5675],
        "RTEA7H" => [427, 427, 5523, 5521],
        "ISEA3H" => [417, 417, 5387, 2902],
        "IVEA3H" => [417, 417, 5383, 2900],
        "RTEA3H" => [417, 417, 5387, 2902],
        other => panic!("nothing has been measured of the compaction on the sample of {other}"),
    }
}

/// The compaction on the suites' own sample of zones, which holds zones of every level, on
/// the poles, beside the seams of the layout and at the pentagons: for every twelfth zone of
/// it, the zones just below it are a set of one level, and this crate compacts that set as the
/// engine does, sequence for sequence. On aperture 7 the set is the children of the zone,
/// thirteen for a hexagon, of which the one at the centre gives way to the zone. On aperture
/// 3 it is the sub-zones two levels down, thirteen for a hexagon, of which the seven within
/// the zone give way to it.
pub fn the_compaction_on_the_sample_is_the_engines<S: Subject>() {
    let g = S::grid();
    let [mut sets, mut taken, mut zones_in, mut zones_out] = [0usize; 4];
    for zone in sample::<S>().into_iter().step_by(SETS_STRIDE) {
        let set = if S::T::APERTURE == 7 {
            g.children(zone)
        } else {
            g.sub_zones(zone, 2).unwrap_or_default()
        };
        if set.is_empty() {
            continue;
        }
        let what = format!(
            "{}, the zones below {} ({:#018x}, seed {SEED:#x})",
            S::NAME,
            g.text_id(zone),
            zone.0
        );
        let ours = g
            .compact_zones(&set)
            .unwrap_or_else(|e| panic!("{what}: the set is refused: {e}"));
        let ids: Vec<u64> = set.iter().map(|z| z.0).collect();
        let theirs: Vec<ZoneId> = o::compact_zones(S::ORACLE, &ids)
            .into_iter()
            .map(ZoneId)
            .collect();
        assert!(
            ours == theirs,
            "{what}: the set of {} zones compacts to {:?} where the engine has {:?}",
            set.len(),
            ours.iter().map(|&z| g.text_id(z)).collect::<Vec<_>>(),
            theirs.iter().map(|&z| g.text_id(z)).collect::<Vec<_>>()
        );
        sets += 1;
        let level = g.resolution(set[0]);
        taken += usize::from(ours.iter().any(|&z| g.resolution(z) < level));
        zones_in += set.len();
        zones_out += ours.len();
    }
    eprintln!(
        "the compaction on the sample of {}: {sets} sets, a coarser zone taken in {taken}, \
         {zones_in} zones in, {zones_out} out",
        S::NAME
    );
    assert_eq!(
        [sets, taken, zones_in, zones_out],
        sample_sets_recorded::<S>(),
        "{}: the sets of the sample, those with a coarser zone taken, the zones in and out",
        S::NAME
    );
}

/// The margin that the engine's source gives for the test of two extents, ten times the
/// spacing of the doubles at one, 2.22e-15: the library holds another, `boxes::ENGINE_MARGIN`,
/// 2e-15, and it is the library's that decides.
const SOURCES_MARGIN: f64 = 10.0 * f64::EPSILON;

/// How a coordinate of a box, in radians, lies below the `edge` of an extent that it must
/// overlap: by no more than the library's margin, so that the two do not meet; by more than
/// the library's margin and no more than the source's, where the two margins part; or by more
/// than both.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Overlap {
    Short,
    Between,
    Clear,
}

fn overlap(edge: f64, coordinate: f64) -> Overlap {
    match (
        edge - boxes::ENGINE_MARGIN > coordinate,
        edge - SOURCES_MARGIN > coordinate,
    ) {
        (false, _) => Overlap::Short,
        (true, false) => Overlap::Between,
        (true, true) => Overlap::Clear,
    }
}

/// A coordinate in degrees whose value in radians, as the crate converts it, lies below
/// `edge` as `wanted` says: sought among the doubles about the one that lies `by` radians
/// below it. `None` where the doubles there are too far apart for one to fall so.
fn degrees_below(edge: f64, by: f64, wanted: Overlap) -> Option<f64> {
    let about = (edge - by).to_degrees();
    (-12i64..=12)
        .map(|step| f64::from_bits(about.to_bits().wrapping_add_signed(step)))
        .find(|&deg| overlap(edge, boxes::radians(&[deg; 4])[0]) == wanted)
}

/// The levels at which the margin is tried, two zones at each.
fn margin_levels<S: Subject>() -> [u8; 6] {
    if S::T::APERTURE == 7 {
        [1, 3, 5, 8, 11, 14]
    } else {
        [1, 4, 8, 14, 20, 33]
    }
}

/// The boxes that the two margins decide differently: those made, those put to the engine,
/// and those of them in which the engine's answer holds the zone, measured.
fn margin_recorded<S: Subject>() -> [usize; 3] {
    match S::NAME {
        "IGEO7" => [19, 19, 19],
        "IVEA7H" => [20, 20, 20],
        "RTEA7H" => [20, 20, 20],
        "ISEA3H" => [21, 21, 19],
        "IVEA3H" => [20, 20, 18],
        "RTEA3H" => [21, 21, 20],
        other => panic!("nothing has been measured of the margin on {other}"),
    }
}

/// The margin by which the extent of a zone must overlap a box is the one the engine's
/// library holds, 2e-15 radian, and not the 2.22e-15 of its source. For two zones at each of
/// six levels, boxes are laid against the northern edge of the zone's extent, and against its
/// eastern edge, so as to overlap it by less than the library's margin, by more than the
/// library's and no more than the source's, and by more than both. The zone is in the answer
/// of `zones_in_box` exactly where the engine's own extent of it meets the box by the
/// library's margin. Where the two margins part, the box is put to the engine as well: its
/// answer holds the zone in as many boxes as are recorded, which the source's margin would
/// leave out of every one; where it does not, on aperture 3, the engine's answer is one that
/// `expected_divergences::engine_omits_zones_of_a_box` allows, its search leaving zones out
/// there as elsewhere. And the engine's answer never holds the zone where the overlap is
/// short.
pub fn the_margin_of_the_test_of_two_extents_is_the_librarys<S: Subject>() {
    use std::f64::consts::PI;
    let g = S::grid();
    let mut rng = boxes::Rng::new(SEED ^ (0x3a26 << 8));
    let [mut between, mut asked, mut held] = [0usize; 3];
    let mut tried = 0;
    for level in margin_levels::<S>() {
        let mut zones = 0;
        while zones < 2 {
            let (lat, lon) = rng.point();
            let zone = g.zone_from_geo(lat, lon, level).unwrap().id();
            if zone == ZoneId::NULL {
                continue;
            }
            let e = engines_extent::<S>(zone);
            let (height, width) = (e[2] - e[0], e[3] - e[1]);
            // A zone whose extent lies clear of the poles and of the antimeridian, with room
            // for a box of its own size to its north and to its east.
            if width <= 0.0 || e[2] + height >= PI / 2.0 || e[3] + width >= PI || e[1] <= -PI {
                continue;
            }
            zones += 1;
            let mid = |a: f64, b: f64, part: f64| (a + (b - a) * part).to_degrees();
            for (by, wanted) in [
                (1.5e-15, Overlap::Short),
                (2.1e-15, Overlap::Between),
                (3.0e-15, Overlap::Clear),
            ] {
                // Against the northern edge, and against the eastern edge.
                let north = degrees_below(e[2], by, wanted).map(|south| {
                    [
                        south,
                        mid(e[1], e[3], 0.25),
                        south + height.to_degrees(),
                        mid(e[1], e[3], 0.75),
                    ]
                });
                let east = degrees_below(e[3], by, wanted).map(|west| {
                    [
                        mid(e[0], e[2], 0.25),
                        west,
                        mid(e[0], e[2], 0.75),
                        west + width.to_degrees(),
                    ]
                });
                for b in [north, east].into_iter().flatten() {
                    let r = boxes::radians(&b);
                    let what = format!(
                        "{} level {level}, {} ({:#018x}), the box {b:?} (seed {SEED:#x})",
                        S::NAME,
                        g.text_id(zone),
                        zone.0
                    );
                    let meets = boxes::meets(&e, &r);
                    assert_eq!(meets, wanted != Overlap::Short, "{what}: {e:?}");
                    let ours = first_in_box::<S>(level, &b, usize::MAX);
                    assert_eq!(
                        ours.contains(&zone),
                        meets,
                        "{what}: the extent {e:?} overlaps the box {wanted:?}"
                    );
                    tried += 1;
                    if wanted == Overlap::Clear {
                        continue;
                    }
                    between += usize::from(wanted == Overlap::Between);
                    if !boxes::engine_survives(S::ORACLE, i32::from(level), &b) {
                        continue;
                    }
                    let theirs = o::list_zones(S::ORACLE, i32::from(level), Some(r));
                    let holds = theirs.contains(&zone.0);
                    if wanted == Overlap::Short {
                        assert!(!holds, "{what}: the engine's answer holds the zone");
                        continue;
                    }
                    asked += 1;
                    held += usize::from(holds);
                    if !holds {
                        // The engine's answer lacks the zone for another reason than the
                        // margin: its search leaves zones out, as it does for other boxes.
                        let ids: Vec<u64> = ours.iter().map(|z| z.0).collect();
                        div::engine_omits_zones_of_a_box::<S>(level, &b, &ids, &theirs)
                            .unwrap_or_else(|why| panic!("{what}: {why}"));
                    }
                }
            }
        }
    }
    eprintln!(
        "the margin on {}: {tried} boxes; {between} where the two margins part, {asked} of \
         them put to the engine, whose answer holds the zone in {held}",
        S::NAME
    );
    assert!(tried >= 60, "{tried} boxes tried");
    assert_eq!(
        [between, asked, held],
        margin_recorded::<S>(),
        "{}: the boxes that the two margins decide differently, those put to the engine, and \
         those in which its answer holds the zone",
        S::NAME
    );
}

/// Coordinates that a caller may hand over as a box: the limits, values just within and just
/// beyond them, noughts of either sign, the smallest and the greatest doubles, and what is no
/// number.
const HOSTILE_COORDINATES: [f64; 24] = [
    f64::NAN,
    f64::INFINITY,
    f64::NEG_INFINITY,
    f64::MAX,
    f64::MIN,
    f64::MIN_POSITIVE,
    5e-324,
    0.0,
    -0.0,
    90.0,
    -90.0,
    180.0,
    -180.0,
    90.000_000_000_000_01,
    -90.000_000_000_000_01,
    180.000_000_000_000_03,
    -180.000_000_000_000_03,
    89.999_999_999_999_99,
    -89.999_999_999_999_99,
    179.999_999_999_999_97,
    -179.999_999_999_999_97,
    38.7223,
    -9.1393,
    1e-9,
];

/// The enumeration of a level, the zones of a box, the estimate and the compaction, asked with
/// hostile input: none panics, none works without end, and each refuses what it refuses with
/// the error that names the cause.
///
/// - **Levels.** Beyond `max_box_level`, the three that take a level refuse it with that
///   limit, whatever the box: on the aperture-7 grids levels 15 and 19 among them, which the
///   grids have and the enumeration does not reach.
/// - **Boxes.** Two thousand boxes drawn from [`HOSTILE_COORDINATES`], at three levels: the
///   zones of the box and the estimate refuse the same boxes with the same error, and a box
///   that is not refused yields its first zones and an estimate that is not nought.
/// - **Cursors.** Every hostile identifier, as the zone to enter after: the enumeration of the
///   level and that of a box refuse the same identifiers, as no zones of the level, and go on
///   after the others.
/// - **Slices.** Every hostile identifier alone, and slices of them of every length to six:
///   an answer, or one of the two refusals of the compaction.
pub fn the_enumerations_the_estimate_and_the_compaction_withstand_hostile_input<S: Subject>() {
    let g = S::grid();
    let any = get_grid(S::NAME).unwrap();
    let max = g.max_box_level();
    let lisbon = degrees(&LISBON);

    // Levels.
    assert_eq!(max, if S::T::APERTURE == 7 { 14 } else { 33 });
    assert_eq!(any.max_box_level(), max);
    for level in [max + 1, 15, 19, 20, 33, 34, 35, 100, 255] {
        let beyond = Error::ResolutionOutOfRange { res: level, max };
        if level <= max {
            assert!(g.zones(level).is_ok() && g.zones_in_box(level, &lisbon).is_ok());
            assert!(g.estimate_zones_in_box(level, &lisbon).is_ok());
            continue;
        }
        let what = format!("{} level {level}", S::NAME);
        assert_eq!(g.zones(level).err(), Some(beyond.clone()), "{what}");
        assert_eq!(any.zones(level).err(), Some(beyond.clone()), "{what}");
        assert_eq!(
            g.zones_in_box(level, &lisbon).err(),
            Some(beyond.clone()),
            "{what}"
        );
        assert_eq!(
            any.zones_in_box(level, &lisbon).err(),
            Some(beyond.clone()),
            "{what}"
        );
        assert_eq!(
            g.estimate_zones_in_box(level, &lisbon).err(),
            Some(beyond.clone()),
            "{what}"
        );
        // The level is refused before the box is looked at.
        let no_box = degrees(&[f64::NAN, 200.0, -100.0, f64::INFINITY]);
        assert_eq!(
            g.zones_in_box(level, &no_box).err(),
            Some(beyond.clone()),
            "{what}"
        );
        assert_eq!(
            g.estimate_zones_in_box(level, &no_box).err(),
            Some(beyond),
            "{what}"
        );
    }

    // Boxes.
    let mut rng = boxes::Rng::new(SEED ^ (0x0b5e << 8));
    let (mut refused, mut answered) = (0usize, 0usize);
    for i in 0..2_000 {
        let b: [f64; 4] = std::array::from_fn(|_| {
            HOSTILE_COORDINATES[(rng.next_u64() % HOSTILE_COORDINATES.len() as u64) as usize]
        });
        let level = [0, max / 2, max][i % 3];
        let what = format!("{} level {level}, the box {b:?} (seed {SEED:#x})", S::NAME);
        let estimate = g.estimate_zones_in_box(level, &degrees(&b));
        match g.zones_in_box(level, &degrees(&b)) {
            Err(e) => {
                assert!(
                    matches!(
                        e,
                        Error::NonFinite { .. }
                            | Error::LatitudeOutOfRange { .. }
                            | Error::LongitudeOutOfRange { .. }
                            | Error::InvertedBox { .. }
                    ),
                    "{what}: {e:?}"
                );
                assert_eq!(estimate.err(), Some(e.clone()), "{what}");
                assert_eq!(
                    any.zones_in_box(level, &degrees(&b)).err(),
                    Some(e),
                    "{what}"
                );
                refused += 1;
            }
            Ok(zones) => {
                let first: Vec<ZoneId> = zones.take(3).collect();
                let estimate = estimate.unwrap_or_else(|e| panic!("{what}: {e:?}"));
                assert!(estimate >= first.len().max(1) as u64, "{what}: {estimate}");
                are_in_the_order_of_the_level::<S>(level, &first, &what);
                answered += 1;
            }
        }
    }
    assert!(
        refused >= 1_500 && answered >= 80,
        "{}: {refused} boxes refused, {answered} answered",
        S::NAME
    );

    // Cursors.
    let hostile: Vec<ZoneId> = hostile_identifiers::<S>().into_iter().map(ZoneId).collect();
    let (mut entered, mut not_of_the_level) = (0usize, 0usize);
    for level in [0, max / 2, max] {
        for &id in &hostile {
            let what = format!("{} level {level}, after {:#018x}", S::NAME, id.0);
            let of_the_level = g.zones(level).unwrap().after(id);
            let of_the_box = g.zones_in_box(level, &lisbon).unwrap().after(id);
            match (of_the_level, of_the_box) {
                (Ok(mut every), Ok(mut some)) => {
                    assert_eq!(g.resolution(id), level, "{what}");
                    let _ = (every.next(), some.next());
                    entered += 1;
                }
                (Err(a), Err(b)) => {
                    assert!(matches!(a, Error::InvalidZone(_)), "{what}: {a:?}");
                    assert_eq!(a, b, "{what}");
                    not_of_the_level += 1;
                }
                (a, b) => panic!(
                    "{what}: the level {} it and the box {} it",
                    if a.is_ok() { "accepts" } else { "refuses" },
                    if b.is_ok() { "accepts" } else { "refuses" }
                ),
            }
        }
    }
    assert!(
        entered >= 3 && not_of_the_level >= 8_000,
        "{}: {entered} cursors entered, {not_of_the_level} refused",
        S::NAME
    );

    // Slices.
    let (mut compacted, mut invalid, mut mixed) = (0usize, 0usize, 0usize);
    let mut tried = |slice: &[ZoneId]| match g.compact_zones(slice) {
        Ok(answer) => {
            assert!(answer.len() <= slice.len(), "{slice:x?}");
            assert_eq!(any.compact_zones(slice).ok(), Some(answer), "{slice:x?}");
            compacted += 1;
        }
        Err(Error::InvalidZone(_)) => invalid += 1,
        Err(Error::MixedLevels { level, other }) => {
            assert_ne!(level, other, "{slice:x?}");
            mixed += 1;
        }
        Err(e) => panic!("{}: {slice:x?} is refused with {e:?}", S::NAME),
    };
    tried(&[]);
    for (i, &id) in hostile.iter().enumerate() {
        tried(&[id]);
        tried(&[id, id]);
        let to = (i + 2 + i % 5).min(hostile.len());
        tried(&hostile[i..to]);
    }
    // The zones of the sample, which are zones: by twos they are of one level or of two.
    let zones = sample::<S>();
    for pair in zones.chunks(2).step_by(7) {
        tried(pair);
    }
    assert!(
        compacted >= 300 && invalid >= 3_000 && mixed >= 100,
        "{}: {compacted} slices compacted, {invalid} refused as holding what is no zone, \
         {mixed} as of more than one level",
        S::NAME
    );
    eprintln!(
        "hostile input on {}: {refused} boxes refused and {answered} answered; {entered} \
         cursors entered and {not_of_the_level} refused; {compacted} slices compacted, \
         {invalid} refused as holding what is no zone, {mixed} as of more than one level",
        S::NAME
    );
}
