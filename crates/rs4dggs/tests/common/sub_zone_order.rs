//! The sub-zone order of the aperture-7 grids against the live engine's, entry by entry.
//!
//! One comparison, [`Orders::compare`], is made of a zone at a depth: the count against the
//! engine's `countSubZones`; the order against its `getSubZones`, entry by entry; the first
//! sub-zone against its `getFirstSubZone`; the zone at every index against the entry that
//! the engine's order holds at that index; and, at the orders and entries of
//! [`Orders::index_is_asked`], the index of an entry against the engine's `getSubZoneIndex`
//! and against the order. Four tests run it, each over a sample of its own:
//! the hexagons that positions give, at depths 1 to 3; the twelve pentagons of every level,
//! at depths 1 to 3, and those of levels 0 and 1 to depth 5; a few hexagons at depths 4 and
//! 5; and identifiers of a broken seam that are built by text. A fifth counts what the engine
//! answers for the first sub-zone at depth 0, and a sixth puts hostile identifiers, depths and
//! indices to the five methods, without the engine.
//!
//! Nothing is admitted here: no class rule and no unresolved region applies to an order, and
//! no tie is allowed. An entry may differ from the engine's in one way alone, which is this
//! crate's stated answer and is re-checked live at every instance: where the engine's entry
//! is an identifier that it cannot read back, the entry here is the null zone, at the same
//! place. The length of an order and every other entry are the engine's own, its null zones
//! and its repeats included.
//!
//! The index of a sub-zone is the engine's wherever it is answered, and it is answered only
//! where the order holds the zone at that index. It parts from the engine's in one way
//! alone, which is counted: where the engine answers an index at which its own order holds
//! another zone, there is none here. Below level 15 the index of every entry asked is its
//! place, here and in the engine.
//!
//! What the broken seams do to an order is counted by class and by the level of the
//! sub-zones, and held, exactly, to what is recorded of each grid (see [`Recorded`]): the
//! entries that are the engine's null zone, those that it cannot read back, those that the
//! order holds already, the orders that hold any of the three, which are not sound, the
//! orders at depth 1 that hold none of them and are all the same not the set of the zone's
//! children, and the entries that have no index, or one that is not their place. Every such
//! order has its sub-zones at level 15 or finer, and its zone reaches
//! the band about the seams; below level 15 every order is asserted sound, and at depth 1
//! the set of the children.
//!
//! The engine is not asked three things (see `common::expected_divergences`): its
//! `getSubZoneAtIndex`, whose answers at a pentagon at an odd depth are not the entries of
//! its own order, and which the vendored binding declares without its depth, so that the
//! call cannot be made from here; its `getFirstSubZone` where that call ends the process,
//! counted by [`div::engine_first_sub_zone_ends_the_process`]; and an order longer than
//! [`div::LONGEST_ORDER_ASKED`].

use std::collections::{BTreeMap, BTreeSet};

use super::expected_divergences as div;
use super::subject::Subject;
use super::suite::{
    SEED, broken_seam_distance_deg, engine_may_be_asked, engine_neighbours, seam_texts, zones,
};
use super::*;
use dggal_oracle as o;
use rs4dggs::indexings::Z7;
use rs4dggs::topologies::HexA7;
use rs4dggs::{Error, ZoneId};

/// The coarsest level of the sub-zones at which a broken seam was measured to spoil an
/// order. Below it every order compared is asserted sound.
const SEAMS_FROM_SUB_ZONE_LEVEL: u8 = 15;

/// The classes that the comparison counts, each by the level of the sub-zones. Every one of
/// them was met in the broken seams alone.
///
/// An entry that is the engine's null zone, and the null zone here.
const ENGINES_NULL_ZONE: &str = "an entry that is the engine's null zone";
/// An entry that is the null zone here where the engine's is an identifier that it cannot
/// read back.
const UNREADABLE: &str = "an entry that the engine cannot read back, the null zone here";
/// An entry that names a zone which the order holds at an earlier place, as the engine's
/// does.
const REPEAT: &str = "an entry that the order holds already";
/// An order that holds an entry of any of the three kinds above.
const NOT_SOUND: &str = "an order that is not sound";
/// An order at depth 1 with none of the three, which all the same is not the set of the
/// zone's children: the two lists are made by two functions of the engine, which part at a
/// broken seam.
const NOT_THE_CHILDREN: &str = "a sound order at depth 1 that is not the set of the children";

/// An entry of an order for which the engine's `getSubZoneIndex` answers -1: it has no
/// index, in the engine and here.
const NO_INDEX: &str = "an entry that has no index, as in the engine";
/// An entry for which the engine answers an index at which its own order holds another
/// zone: it has no index here.
const INDEX_OF_ANOTHER_ZONE: &str =
    "an entry whose index in the engine is the place of another zone, and that has none here";
/// An entry whose index, in the engine and here, is another place of the order, which holds
/// the same zone there.
const INDEX_OF_ANOTHER_PLACE: &str = "an entry whose index is another place of the same zone";
/// A zone of a neighbour's order that this order does not hold, for which the engine
/// answers an index all the same: it has none here.
const INDEX_OF_NO_ENTRY: &str =
    "a zone that the order does not hold, with an index in the engine and none here";

/// Every class counted, with its count at each level of the sub-zones at which it was met.
type ClassCounts = &'static [(&'static str, &'static [(u8, usize)])];

/// What is recorded of one sample of one grid. The comparison is of the whole of it, so
/// that an order more or fewer, an entry more or fewer, or a class that appears, vanishes
/// or moves by one at one level fails it.
struct Recorded {
    /// At each depth: the orders compared with the engine's, and the entries in them.
    by_depth: &'static [(u8, usize, usize)],
    /// The orders of which the engine's `getFirstSubZone` was not asked.
    first_not_asked: usize,
    /// The orders of pentagons at odd depths, at which the engine's search by index is at
    /// fault and the zone at each index is held to the engine's order.
    pentagons_at_odd_depths: usize,
    /// The entries whose index was asked, here and of the engine; the zones of a neighbour's
    /// order asked likewise; and those of the latter that have an index here, being
    /// sub-zones of both.
    indices: [usize; 3],
    classes: ClassCounts,
}

/// Tallies of the sub-zone orders of one sample, each compared with the engine's.
#[derive(Default)]
struct Orders {
    /// By depth: orders compared, and entries in them.
    by_depth: BTreeMap<u8, (usize, usize)>,
    first_not_asked: usize,
    pentagons_at_odd_depths: usize,
    /// Indices asked of entries, of the zones of a neighbour's order, and those of the
    /// latter that have one.
    indices: [usize; 3],
    /// The longest order asked of the engine.
    longest: u64,
    /// Every class met, by the level of the sub-zones.
    classes: BTreeMap<&'static str, BTreeMap<u8, usize>>,
    /// The orders counted in a class, and the farthest that the centroid of the zone of one
    /// lies from a broken seam, in degrees.
    in_a_class: usize,
    reach_deg: f64,
}

impl Orders {
    fn count(&mut self, class: &'static str, level: u8, n: usize) {
        if n > 0 {
            *self
                .classes
                .entry(class)
                .or_default()
                .entry(level)
                .or_default() += n;
        }
    }

    /// Holds an order that is counted in a class to the broken seams: its sub-zones lie at
    /// level 15 or finer, and its zone reaches the band about a seam, which is to say that
    /// the centroid of the zone lies within [`div::BROKEN_SEAM_BAND_DEG`] of one, or farther
    /// by no more than the distance from the centroid to the farthest vertex.
    fn held_to_the_seams<S: Subject<T = HexA7, I = Z7>>(
        &mut self,
        id: ZoneId,
        level: u8,
        at: &str,
    ) {
        let g = S::grid();
        assert!(
            level >= SEAMS_FROM_SUB_ZONE_LEVEL,
            "{at}: an order with its sub-zones at level {level}, coarser than \
             {SEAMS_FROM_SUB_ZONE_LEVEL}, is counted in a class of the broken seams"
        );
        let c = g.centroid(id);
        let radius = g
            .vertices(id)
            .iter()
            .map(|v| arc_deg((c.lat, c.lon), (v.lat, v.lon)))
            .fold(0.0, f64::max);
        let distance = broken_seam_distance_deg::<S>(c.lat, c.lon);
        assert!(
            distance <= div::BROKEN_SEAM_BAND_DEG + radius,
            "{at}: the order is counted in a class of the broken seams, and its zone, of \
             radius {radius:e} degrees, lies {distance:e} degrees from them"
        );
        self.in_a_class += 1;
        self.reach_deg = self.reach_deg.max(distance);
    }

    /// Whether the index of entry `i` of an order of `count` entries at `depth` is asked, of
    /// this crate and of the engine. The walk that finds an index costs one centroid for
    /// each entry before it, so that asking every entry of an order costs as the square of
    /// its length: every entry is asked at depths 1 and 2, and at the deeper orders one in
    /// 8, 64 and 512, with the last.
    fn index_is_asked(depth: u8, i: usize, count: u64) -> bool {
        let stride = match depth {
            0..=2 => 1,
            3 => 8,
            4 => 64,
            _ => 512,
        };
        i % stride == 0 || i as u64 + 1 == count
    }

    /// The index of `sub` among the sub-zones `ours` of `id`, here and in the engine, where
    /// `place` is the place at which the order holds it, if it does. The index answered
    /// here is the engine's, and the order holds the zone there; where none is answered,
    /// the engine answers -1, or an index at which the order holds another zone. Answers
    /// whether the pair fell in a class of the broken seams, which is counted by the level
    /// of the sub-zones; and whether the zone has an index here.
    fn index_of<S: Subject<T = HexA7, I = Z7>>(
        &mut self,
        id: ZoneId,
        ours: &[ZoneId],
        sub: ZoneId,
        place: Option<usize>,
        level: u8,
        at: &str,
    ) -> (bool, bool) {
        let g = S::grid();
        let at = format!("{at}: the index of {}", g.text_id(sub));
        let engine = o::sub_zone_index(S::ORACLE, id.0, sub.0);
        let index = g.sub_zone_index(id, sub);
        let class = match index {
            Ok(Some(j)) => {
                assert_eq!(j as i64, engine, "{at}, against the engine's");
                assert_eq!(
                    ours[j as usize], sub,
                    "{at} is {j}, and the order holds another zone there"
                );
                let place = place.unwrap_or_else(|| panic!("{at} is {j}, and no place was found"));
                (place != j as usize).then_some(INDEX_OF_ANOTHER_PLACE)
            }
            Ok(None) if engine == -1 => place.map(|_| NO_INDEX),
            Ok(None) => {
                let there = usize::try_from(engine).ok().and_then(|j| ours.get(j));
                assert!(
                    there.is_some_and(|&zone| zone != sub),
                    "{at} is none here and {engine} in the engine, where the order holds {:?}",
                    there.map(|&zone| g.text_id(zone))
                );
                Some(if place.is_some() {
                    INDEX_OF_ANOTHER_ZONE
                } else {
                    INDEX_OF_NO_ENTRY
                })
            }
            Err(e) => panic!("{at} is refused: {e}"),
        };
        if let Some(class) = class {
            self.count(class, level, 1);
        }
        (class.is_some(), index != Ok(None))
    }

    /// The sub-zones of `id` at `depth`, compared with the engine's:
    /// - the count against `countSubZones`, and the length of both orders against it;
    /// - the order against `getSubZones`, entry by entry: the same identifier, the null zone
    ///   where the engine's is the null zone, and the null zone where the engine's is an
    ///   identifier that it cannot read back, which is re-checked live;
    /// - the zone at every index against that entry, and the index past the last refused;
    /// - the first sub-zone against `getFirstSubZone`, which must itself be the first entry
    ///   of the engine's order; where that call would end the process it is not made, the
    ///   order is counted, and the first sub-zone here is held to the first entry of the
    ///   engine's order;
    /// - where the order holds a null zone or a repeat: the classes counted, the order held
    ///   to the broken seams, and every other entry of it put to the engine's reader;
    /// - from level 15, every entry that is a zone is one that this crate reads back from
    ///   its own text;
    /// - at depth 1, where the order is sound, the order as a set against the children, and
    ///   an order that is not that set counted and held to the broken seams;
    /// - with `ask_index`: the index of each entry that [`Orders::index_is_asked`] names,
    ///   against the engine's `getSubZoneIndex` and against the order, as
    ///   [`Orders::index_of`] does; the null zone refused as a sub-zone; and, at depth 1, the
    ///   index of every zone of the order of the zone's first neighbour, which is a
    ///   sub-zone of this zone too where the two share it across their boundary and no
    ///   sub-zone of it otherwise. An entry whose index is not its place, or a zone with an
    ///   index in the engine and none here, is counted, and its order held to the broken
    ///   seams.
    ///
    /// An order is held to the seams once, whatever the number of its classes.
    fn compare<S: Subject<T = HexA7, I = Z7>>(&mut self, id: ZoneId, depth: u8, ask_index: bool) {
        let g = S::grid();
        let text = g.text_id(id);
        let at = format!("{text} at depth {depth} on {} (seed {SEED:#x})", S::NAME);
        let level = g.resolution(id) + depth;
        assert!(
            engine_may_be_asked::<S>(id.0),
            "{at}: the engine cannot read the zone back"
        );
        let count = g.count_sub_zones(id, depth).unwrap();
        assert_eq!(
            count,
            o::count_sub_zones(S::ORACLE, id.0, i32::from(depth)),
            "{at}: the count, against the engine's"
        );
        assert!(
            count <= div::LONGEST_ORDER_ASKED,
            "{at}: an order of {count} sub-zones is longer than the engine is asked for"
        );
        let theirs = o::sub_zones(S::ORACLE, id.0, i32::from(depth));
        let ours = g.sub_zones(id, depth).unwrap();
        assert_eq!(
            (ours.len() as u64, theirs.len() as u64),
            (count, count),
            "{at}: the lengths of this crate's order and of the engine's, against the count"
        );

        let mut seen = BTreeSet::new();
        let (mut nulls, mut unreadable, mut repeats) = (0usize, 0usize, 0usize);
        for (i, (&sub, &engine)) in ours.iter().zip(&theirs).enumerate() {
            if sub == ZoneId::NULL {
                if engine == o::NULL_ZONE {
                    nulls += 1;
                } else {
                    assert!(
                        !o::engine_can_read(S::ORACLE, engine),
                        "{at}: entry {i} is the null zone here and {} in the engine's order, \
                         which the engine reads back",
                        o::text_id(S::ORACLE, engine)
                    );
                    unreadable += 1;
                }
            } else {
                assert_eq!(
                    sub.0,
                    engine,
                    "{at}: entry {i} is {} here and {} in the engine's order",
                    g.text_id(sub),
                    o::text_id(S::ORACLE, engine)
                );
                if !seen.insert(sub) {
                    repeats += 1;
                }
                if level >= SEAMS_FROM_SUB_ZONE_LEVEL {
                    assert_eq!(
                        g.zone_from_text(&g.text_id(sub)).map(|z| z.id()),
                        Ok(sub),
                        "{at}: entry {i} is an identifier that this crate does not read back"
                    );
                }
            }
            assert_eq!(
                g.sub_zone_at_index(id, depth, i as u64),
                Ok(sub),
                "{at}: the zone at index {i}, against entry {i} of the engine's order"
            );
        }
        assert_eq!(
            g.sub_zone_at_index(id, depth, count),
            Err(Error::IndexOutOfRange {
                index: count,
                count
            }),
            "{at}: the index past the last"
        );

        let first = g.first_sub_zone(id, depth).unwrap();
        assert_eq!(
            first, ours[0],
            "{at}: the first sub-zone, against the first entry of the order"
        );
        if div::engine_first_sub_zone_ends_the_process::<S>(id, depth) {
            assert_eq!(
                first.0, theirs[0],
                "{at}: the first sub-zone, where the engine's own call is not made, against \
                 the first entry of the engine's order"
            );
            self.first_not_asked += 1;
        } else {
            assert_eq!(
                o::first_sub_zone(S::ORACLE, id.0, i32::from(depth)),
                theirs[0],
                "{at}: the engine's first sub-zone, against the first entry of its own order"
            );
        }

        let sound = nulls + unreadable + repeats == 0;
        let mut in_a_class = !sound;
        if !sound {
            for &sub in ours.iter().filter(|&&sub| sub != ZoneId::NULL) {
                assert!(
                    o::engine_can_read(S::ORACLE, sub.0),
                    "{at}: the order holds {}, which the engine cannot read back",
                    g.text_id(sub)
                );
            }
            self.count(ENGINES_NULL_ZONE, level, nulls);
            self.count(UNREADABLE, level, unreadable);
            self.count(REPEAT, level, repeats);
            self.count(NOT_SOUND, level, 1);
        } else if depth == 1 {
            let children: BTreeSet<ZoneId> = g.children(id).into_iter().collect();
            if children != seen {
                self.count(NOT_THE_CHILDREN, level, 1);
                in_a_class = true;
            }
        }

        if ask_index {
            for (i, &sub) in ours.iter().enumerate() {
                if sub == ZoneId::NULL {
                    assert!(
                        matches!(g.sub_zone_index(id, sub), Err(Error::InvalidZone(_))),
                        "{at}: the index of the null zone, entry {i}"
                    );
                } else if Self::index_is_asked(depth, i, count) {
                    in_a_class |= self.index_of::<S>(id, &ours, sub, Some(i), level, &at).0;
                    self.indices[0] += 1;
                }
            }
            if let (1, Some(&neighbour)) = (depth, g.neighbors(id).first()) {
                for sub in g.sub_zones(neighbour, depth).unwrap() {
                    if sub != ZoneId::NULL {
                        let place = ours.iter().position(|&entry| entry == sub);
                        let (class, found) = self.index_of::<S>(id, &ours, sub, place, level, &at);
                        in_a_class |= class;
                        self.indices[1] += 1;
                        self.indices[2] += usize::from(found);
                    }
                }
            }
        }
        if in_a_class {
            self.held_to_the_seams::<S>(id, level, &at);
        }

        let tally = self.by_depth.entry(depth).or_default();
        tally.0 += 1;
        tally.1 += ours.len();
        self.longest = self.longest.max(count);
        if g.is_pentagon(id) && depth % 2 == 1 {
            self.pentagons_at_odd_depths += 1;
        }
    }

    /// Prints the tallies of one sample, and holds them to what is recorded of it on the
    /// grid under test, exactly. A grid of which nothing is recorded fails, with its tallies
    /// printed, and is never handed another grid's.
    fn finish<S: Subject>(&self, what: &str, recorded: Option<Recorded>) {
        let by_depth: Vec<(u8, usize, usize)> = self
            .by_depth
            .iter()
            .map(|(&depth, &(orders, entries))| (depth, orders, entries))
            .collect();
        eprintln!(
            "{what}: {} orders as the engine's and {} entries in them, by depth (depth, \
             orders, entries) {by_depth:?}; the longest of {} entries; {} orders of pentagons \
             at odd depths; the engine's first sub-zone not asked at {}; the index asked of \
             entries, of the zones of a neighbour's order, and those of the latter that have \
             one {:?}; {} orders counted in a class of the broken seams, the farthest zone \
             {:e} degrees from them",
            by_depth.iter().map(|d| d.1).sum::<usize>(),
            by_depth.iter().map(|d| d.2).sum::<usize>(),
            self.longest,
            self.pentagons_at_odd_depths,
            self.first_not_asked,
            self.indices,
            self.in_a_class,
            self.reach_deg
        );
        for (class, by_level) in &self.classes {
            eprintln!(
                "  {class}: {} {:?}",
                by_level.values().sum::<usize>(),
                by_level.iter().collect::<Vec<_>>()
            );
        }
        let Some(recorded) = recorded else {
            panic!(
                "nothing is recorded of {what} on {}: the tallies printed above are what this \
                 run met, and the table is to be filled from them only once each class has \
                 been characterised, and never from another grid",
                S::NAME
            )
        };
        assert_eq!(
            (
                by_depth.as_slice(),
                self.first_not_asked,
                self.pentagons_at_odd_depths,
                self.indices
            ),
            (
                recorded.by_depth,
                recorded.first_not_asked,
                recorded.pentagons_at_odd_depths,
                recorded.indices
            ),
            "{what} on {}: the orders and entries compared at each depth, the orders at which \
             the engine's first sub-zone was not asked, the orders of pentagons at odd \
             depths, and the indices asked of entries, of the zones of a neighbour's order, \
             and those of the latter that have one, against what is recorded of the grid \
             (seed {SEED:#x})",
            S::NAME
        );
        let classes: BTreeMap<&str, BTreeMap<u8, usize>> = recorded
            .classes
            .iter()
            .map(|(class, by_level)| (*class, by_level.iter().copied().collect()))
            .collect();
        assert_eq!(
            self.classes,
            classes,
            "{what} on {}: the classes counted by the level of the sub-zones, against what \
             is recorded of the grid (seed {SEED:#x})",
            S::NAME
        );
    }
}

/// One hexagon in this many of the sample is compared at each depth: every one at depth 1,
/// one in six at depth 2 and one in forty at depth 3. An order grows sevenfold with each
/// depth, and with it the comparison of the zone at every index.
const HEXAGON_STRIDES: [(u8, usize); 3] = [(1, 1), (2, 6), (3, 40)];

/// Of the hexagons compared at a depth, one in this many has the index of its entries
/// asked, of this crate and of the engine: at the pentagons, in the deep orders and at the
/// identifiers of the seam every order has.
const HEXAGON_INDEX_STRIDE: usize = 4;

/// What differs by grid in what is recorded of the hexagons: the orders and entries at
/// each depth, the indices asked, and the classes.
type RecordedOfTheHexagons = (&'static [(u8, usize, usize)], [usize; 3], ClassCounts);

/// What is recorded of the orders of the hexagons that positions give. The zones are those
/// that each grid's projection makes of the points, and so the counts differ by grid.
fn hexagon_orders_recorded<S: Subject>() -> Option<Recorded> {
    let (by_depth, indices, classes): RecordedOfTheHexagons = match S::NAME {
        "IGEO7" => (
            &[(1, 4_993, 64_909), (2, 746, 41_030), (3, 106, 40_174)],
            [27_529, 16_166, 2_480],
            &[
                (ENGINES_NULL_ZONE, &[(15, 116), (17, 66), (19, 24)]),
                (UNREADABLE, &[(16, 14), (18, 6)]),
                (NOT_SOUND, &[(15, 14), (16, 5), (17, 11), (18, 4), (19, 4)]),
                (NOT_THE_CHILDREN, &[(16, 3), (18, 2)]),
                (NO_INDEX, &[(15, 6), (16, 8), (17, 8), (18, 6)]),
                (INDEX_OF_ANOTHER_ZONE, &[(16, 19), (18, 5)]),
                (INDEX_OF_NO_ENTRY, &[(16, 2)]),
            ],
        ),
        "IVEA7H" => (
            &[(1, 4_965, 64_545), (2, 740, 40_700), (3, 101, 38_279)],
            [27_113, 15_961, 2_465],
            &[
                (ENGINES_NULL_ZONE, &[(15, 44), (17, 42), (19, 24)]),
                (UNREADABLE, &[(16, 1), (18, 4)]),
                (NOT_SOUND, &[(15, 7), (16, 1), (17, 7), (18, 3), (19, 4)]),
                (NOT_THE_CHILDREN, &[(16, 1), (18, 2)]),
                (NO_INDEX, &[(15, 7), (17, 12), (18, 53), (19, 8)]),
            ],
        ),
        "RTEA7H" => (
            &[(1, 4_960, 64_480), (2, 762, 41_910), (3, 108, 40_932)],
            [27_888, 16_148, 2_492],
            &[
                (ENGINES_NULL_ZONE, &[(15, 54), (17, 122), (19, 24)]),
                (UNREADABLE, &[(16, 1), (18, 4)]),
                (NOT_SOUND, &[(15, 8), (16, 1), (17, 9), (18, 3), (19, 4)]),
                (NOT_THE_CHILDREN, &[(16, 1), (18, 2)]),
                (NO_INDEX, &[(15, 2), (16, 13), (17, 19), (19, 4)]),
                (INDEX_OF_ANOTHER_ZONE, &[(16, 32)]),
                (INDEX_OF_NO_ENTRY, &[(16, 6)]),
            ],
        ),
        _ => return None,
    };
    Some(Recorded {
        by_depth,
        first_not_asked: 0,
        pentagons_at_odd_depths: 0,
        indices,
        classes,
    })
}

/// The hexagons that positions give (those of `suite::zones`), each compared with the
/// engine as [`Orders::compare`] does, at the depths and strides of [`HEXAGON_STRIDES`],
/// wherever the sub-zones stay within the finest level. The engine's `getFirstSubZone` is
/// asked at every one of these orders.
pub fn sub_zone_orders_of_the_hexagons_as_the_engine_lists_them<S: Subject<T = HexA7, I = Z7>>() {
    the_limit_on_a_list_of_sub_zones_is_its_ceiling();
    let g = S::grid();
    let mut orders = Orders::default();
    let hexagons = zones::<S>()
        .into_iter()
        .filter(|&id| id != ZoneId::NULL && !g.is_pentagon(id));
    for (n, id) in hexagons.enumerate() {
        for (depth, stride) in HEXAGON_STRIDES {
            if n % stride == 0 && g.resolution(id) + depth <= g.max_resolution() {
                orders.compare::<S>(id, depth, n % (stride * HEXAGON_INDEX_STRIDE) == 0);
            }
        }
    }
    orders.finish::<S>(
        "sub-zone orders of the hexagons",
        hexagon_orders_recorded::<S>(),
    );
}

/// The levels whose twelve pentagons are compared to depth 5: one even and one odd, since
/// the scanlines of an order run in another direction at the sub-zones of an odd level.
const PENTAGON_LEVELS_TO_DEPTH_5: [u8; 2] = [0, 1];

/// What is recorded of the orders of the pentagons: the same zones on every grid, and, as
/// measured on IGEO7, IVEA7H and RTEA7H, the same counts.
fn pentagon_orders_recorded<S: Subject>() -> Option<Recorded> {
    match S::NAME {
        "IGEO7" | "IVEA7H" | "RTEA7H" => Some(Recorded {
            by_depth: &[
                (1, 228, 2_508),
                (2, 216, 9_936),
                (3, 204, 64_464),
                (4, 24, 48_984),
                (5, 24, 341_064),
            ],
            // The pentagon of the south pole at the levels 1, 3, 5 and so on to 17, at
            // depth 2.
            first_not_asked: 9,
            pentagons_at_odd_depths: 228 + 204 + 24,
            indices: [21_973, 2_900, 441],
            classes: &[
                (ENGINES_NULL_ZONE, &[(15, 93), (17, 568), (19, 568)]),
                (UNREADABLE, &[(16, 31), (18, 115)]),
                (REPEAT, &[(16, 6), (17, 7), (18, 12), (19, 7)]),
                (NOT_SOUND, &[(15, 9), (16, 4), (17, 9), (18, 5), (19, 9)]),
                (
                    NO_INDEX,
                    &[(15, 12), (16, 16), (17, 18), (18, 18), (19, 18)],
                ),
                (INDEX_OF_ANOTHER_ZONE, &[(16, 110), (18, 130)]),
                (
                    INDEX_OF_ANOTHER_PLACE,
                    &[(16, 2), (17, 1), (18, 2), (19, 1)],
                ),
                (INDEX_OF_NO_ENTRY, &[(16, 8), (18, 2)]),
            ],
        }),
        _ => None,
    }
}

/// The twelve pentagons of every level, built by text, each compared with the engine as
/// [`Orders::compare`] does at depths 1 to 3, and those of the levels of
/// [`PENTAGON_LEVELS_TO_DEPTH_5`] at depths 4 and 5 as well, wherever the sub-zones stay
/// within the finest level.
///
/// At a pentagon at an odd depth the engine's own search by index answers another zone
/// than its order holds at that index, from the scanline after the pentagon's own; here the
/// zone at every index is the entry of the engine's order, and these are the orders that
/// tell the two apart.
///
/// The pentagon of the south pole at an odd level, at depth 2, is where the engine's
/// `getFirstSubZone` ends the process: nine orders, of which the call is not made, and at
/// which the first sub-zone here is the first entry of the engine's order.
pub fn sub_zone_orders_of_the_pentagons_as_the_engine_lists_them<S: Subject<T = HexA7, I = Z7>>() {
    the_limit_on_a_list_of_sub_zones_is_its_ceiling();
    let g = S::grid();
    let mut orders = Orders::default();
    for level in 0..g.max_resolution() {
        for base in 0..12 {
            let text = format!("{base:02}{}", "0".repeat(usize::from(level)));
            let id = g.zone_from_text(&text).unwrap().id();
            assert!(g.is_pentagon(id), "{text} is no pentagon");
            let deepest = if PENTAGON_LEVELS_TO_DEPTH_5.contains(&level) {
                5
            } else {
                3
            };
            for depth in 1..=deepest {
                if level + depth <= g.max_resolution() {
                    orders.compare::<S>(id, depth, true);
                }
            }
        }
    }
    orders.finish::<S>(
        "sub-zone orders of the pentagons",
        pentagon_orders_recorded::<S>(),
    );
    // The same answer as literals at four of the nine, read from the engine's order in a
    // process of its own.
    for (zone, first) in [
        ("110", "11022"),
        ("11000", "1100022"),
        ("1100000", "110000022"),
        ("110000000", "11000000022"),
    ] {
        let id = g.zone_from_text(zone).unwrap().id();
        assert!(
            div::engine_first_sub_zone_ends_the_process::<S>(id, 2),
            "{zone} at depth 2 is not kept from the engine's first sub-zone"
        );
        assert_eq!(
            g.first_sub_zone(id, 2).map(|sub| g.text_id(sub)),
            Ok(first.to_string()),
            "the first sub-zone of {zone} at depth 2"
        );
    }
}

/// One hexagon in this many of the sample, of those whose sub-zones stay within the finest
/// level, is compared at depth 4, and one in this many at depth 5.
const DEEP_STRIDES: [(u8, usize); 2] = [(4, 400), (5, 800)];

/// A zone of level 10 on a broken seam, whose order at depth 5 holds the engine's null zone
/// at 99 of its 17,053 places.
const DEEP_SEAM_ZONE: &str = "010004000400";

/// What is recorded of the deep orders. The hexagons are those that each grid's projection
/// makes of the points; as measured on IGEO7, IVEA7H and RTEA7H they are as many on each,
/// and the one order among them that is not sound is that of [`DEEP_SEAM_ZONE`].
fn deep_orders_recorded<S: Subject>() -> Option<Recorded> {
    match S::NAME {
        "IGEO7" | "IVEA7H" | "RTEA7H" => Some(Recorded {
            by_depth: &[(4, 10, 24_490), (5, 6, 102_318)],
            first_not_asked: 0,
            pentagons_at_odd_depths: 0,
            // One entry in 64 at depth 4 and one in 512 at depth 5, with the last of each.
            indices: [610, 0, 0],
            classes: &[(ENGINES_NULL_ZONE, &[(15, 99)]), (NOT_SOUND, &[(15, 1)])],
        }),
        _ => None,
    }
}

/// Orders of hexagons at depths 4 and 5, of 2,449 and 17,053 entries, each compared with
/// the engine as [`Orders::compare`] does, the zone at every index included: the hexagons
/// that positions give, at the strides of [`DEEP_STRIDES`], and [`DEEP_SEAM_ZONE`] at depth
/// 5. The pentagons at these depths are in
/// `sub_zone_orders_of_the_pentagons_as_the_engine_lists_them`.
pub fn deep_sub_zone_orders_as_the_engine_lists_them<S: Subject<T = HexA7, I = Z7>>() {
    the_limit_on_a_list_of_sub_zones_is_its_ceiling();
    let g = S::grid();
    let mut orders = Orders::default();
    let hexagons: Vec<ZoneId> = zones::<S>()
        .into_iter()
        .filter(|&id| id != ZoneId::NULL && !g.is_pentagon(id))
        .collect();
    for (depth, stride) in DEEP_STRIDES {
        let within = hexagons
            .iter()
            .filter(|&&id| g.resolution(id) + depth <= g.max_resolution());
        for &id in within.step_by(stride) {
            orders.compare::<S>(id, depth, true);
        }
    }
    orders.compare::<S>(g.zone_from_text(DEEP_SEAM_ZONE).unwrap().id(), 5, true);
    orders.finish::<S>("deep sub-zone orders", deep_orders_recorded::<S>());
}

/// One identifier in this many of the sample is compared at each depth: every one at depth
/// 1, one in four at depth 2 and one in sixteen at depth 3.
const SEAM_STRIDES: [(u8, usize); 3] = [(1, 1), (2, 4), (3, 16)];

/// What is recorded of the orders at the identifiers of `suite::seam_texts`: the same texts
/// on every grid, and, as measured on IGEO7, IVEA7H and RTEA7H, the same counts.
fn seam_orders_recorded<S: Subject>() -> Option<Recorded> {
    match S::NAME {
        "IGEO7" | "IVEA7H" | "RTEA7H" => Some(Recorded {
            by_depth: &[(1, 1_898, 24_672), (2, 379, 20_845), (3, 71, 26_909)],
            first_not_asked: 0,
            // The one pentagon of the sample, `0000000000000000`, at depth 1.
            pentagons_at_odd_depths: 1,
            indices: [47_441, 23_694, 3_401],
            classes: &[
                (ENGINES_NULL_ZONE, &[(15, 24), (17, 227), (19, 1_293)]),
                (UNREADABLE, &[(16, 21), (18, 78), (19, 10)]),
                (REPEAT, &[(18, 3)]),
                (
                    NOT_SOUND,
                    &[(15, 7), (16, 4), (17, 32), (18, 11), (19, 136)],
                ),
                (NOT_THE_CHILDREN, &[(16, 5), (18, 52)]),
                (
                    NO_INDEX,
                    &[(15, 5), (16, 80), (17, 236), (18, 1_110), (19, 1_570)],
                ),
                (INDEX_OF_ANOTHER_ZONE, &[(16, 80), (18, 59)]),
                (INDEX_OF_NO_ENTRY, &[(16, 4), (18, 2)]),
            ],
        }),
        _ => None,
    }
}

/// The orders at identifiers of the broken seam that are built by text (those of
/// `suite::seam_texts`, of levels 14 to 19), each compared with the engine as
/// [`Orders::compare`] does, at the depths and strides of [`SEAM_STRIDES`], wherever the
/// sub-zones stay within the finest level. Here the classes of the seams are met in
/// number: entries that are the engine's null zone, entries that it cannot read back,
/// repeats, and orders at depth 1 that are not the set of the children.
///
/// A text that this crate refuses passes through a deleted child of the pentagon, and the
/// engine must answer its null zone for it.
pub fn sub_zone_orders_at_seam_identifiers_built_by_text<S: Subject<T = HexA7, I = Z7>>() {
    the_limit_on_a_list_of_sub_zones_is_its_ceiling();
    let g = S::grid();
    let mut orders = Orders::default();
    let mut read = 0usize;
    for text in seam_texts() {
        let engine_id = o::zone_from_text(S::ORACLE, &text);
        let Ok(zone) = g.zone_from_text(&text) else {
            assert_eq!(
                engine_id,
                o::NULL_ZONE,
                "{text}: the engine reads a zone this crate refuses"
            );
            continue;
        };
        assert_eq!(
            zone.id().0,
            engine_id,
            "{text}: the engine reads another zone"
        );
        for (depth, stride) in SEAM_STRIDES {
            if read % stride == 0 && zone.resolution() + depth <= g.max_resolution() {
                orders.compare::<S>(zone.id(), depth, true);
            }
        }
        read += 1;
    }
    orders.finish::<S>(
        "sub-zone orders at seam identifiers built by text",
        seam_orders_recorded::<S>(),
    );
}

/// What is recorded of the engine at depth 0 over the zones of `suite::zones`: the zones
/// asked; those at which its first sub-zone is a neighbour of the zone, its null zone, the
/// zone itself, and another identifier; and those at which its own order of one entry holds
/// another identifier than the zone.
fn depth_0_recorded<S: Subject>() -> Option<[usize; 6]> {
    match S::NAME {
        "IGEO7" => Some([5_607, 5_594, 6, 3, 4, 10]),
        "IVEA7H" => Some([5_579, 5_570, 5, 1, 3, 5]),
        "RTEA7H" => Some([5_574, 5_566, 4, 1, 3, 5]),
        _ => None,
    }
}

/// The sub-zones at depth 0, a characterised divergence, re-checked live at every zone that
/// positions give: this crate answers the zone itself, which is its only sub-zone at depth
/// 0, as the engine's own `countSubZones` has it. The engine's `getFirstSubZone` quantises
/// a vertex of the zone at the zone's own level, and so names a neighbour of the zone at
/// most zones, re-checked to be in the engine's own list of the zone's neighbours, and at
/// some its null zone. Its `getSubZones` quantises the centroid of the zone, and holds the
/// zone itself. At a broken seam either may name another identifier still, and every zone
/// at which one does is held to the band about the seams. The counts are held exactly, per
/// grid.
pub fn first_sub_zone_at_depth_0_departs_from_the_engine_as_recorded<
    S: Subject<T = HexA7, I = Z7>,
>() {
    let g = S::grid();
    let [
        mut asked,
        mut neighbour,
        mut null,
        mut itself,
        mut another,
        mut another_order,
    ] = [0usize; 6];
    for id in zones::<S>() {
        if id == ZoneId::NULL {
            continue;
        }
        let text = g.text_id(id);
        assert_eq!(
            g.count_sub_zones(id, 0),
            Ok(1),
            "count of {text} at depth 0"
        );
        assert_eq!(
            g.first_sub_zone(id, 0),
            Ok(id),
            "first sub-zone of {text} at depth 0"
        );
        assert_eq!(
            g.sub_zone_at_index(id, 0, 0),
            Ok(id),
            "sub-zone 0 of {text} at depth 0"
        );
        assert_eq!(
            g.sub_zones(id, 0),
            Ok(vec![id]),
            "sub-zones of {text} at depth 0"
        );
        assert!(
            engine_may_be_asked::<S>(id.0),
            "the engine cannot read {text} back (seed {SEED:#x})"
        );
        assert_eq!(
            o::count_sub_zones(S::ORACLE, id.0, 0),
            1,
            "the engine's count of {text} at depth 0"
        );
        let held_to_the_seams = |what: &str, answer: u64| {
            let c = g.centroid(id);
            assert!(
                div::near_broken_seam::<S>(c.lat, c.lon),
                "{what} of {text} at depth 0 is {}, and the zone lies {:e} degrees from the \
                 broken seams (seed {SEED:#x})",
                o::text_id(S::ORACLE, answer),
                broken_seam_distance_deg::<S>(c.lat, c.lon)
            );
        };
        let theirs = o::sub_zones(S::ORACLE, id.0, 0);
        assert_eq!(
            theirs.len(),
            1,
            "the engine's sub-zones of {text} at depth 0"
        );
        if theirs[0] != id.0 {
            held_to_the_seams("the one entry of the engine's order", theirs[0]);
            another_order += 1;
        }
        let first = o::first_sub_zone(S::ORACLE, id.0, 0);
        if first == id.0 {
            itself += 1;
        } else if first == o::NULL_ZONE {
            null += 1;
        } else if engine_neighbours::<S>(id.0).contains(&first) {
            neighbour += 1;
        } else {
            held_to_the_seams(
                "the engine's first sub-zone, no neighbour of the zone,",
                first,
            );
            another += 1;
        }
        asked += 1;
    }
    eprintln!(
        "sub-zones at depth 0: of {asked} zones the engine's first sub-zone is a neighbour at \
         {neighbour}, its null zone at {null}, the zone itself at {itself} and another \
         identifier at {another}; its order holds another identifier than the zone at \
         {another_order}"
    );
    let Some(recorded) = depth_0_recorded::<S>() else {
        panic!(
            "nothing is recorded of the engine's sub-zones at depth 0 on {}: the counts \
             printed above are what this run met",
            S::NAME
        )
    };
    assert_eq!(
        [asked, neighbour, null, itself, another, another_order],
        recorded,
        "the zones asked; those at which the engine's first sub-zone at depth 0 is a \
         neighbour, its null zone, the zone itself and another identifier; and those at \
         which its order holds another identifier than the zone (seed {SEED:#x})"
    );
}

/// The five sub-zone methods, asked of hostile identifiers, depths and indices: the null
/// zone, the identifiers at the ends of the range, those with no cell behind them, and 3,000
/// drawn at random from every level; depths from 0 to 255; indices from 0 to the greatest.
/// None panics or fails to end; the typed grid, `AnyGrid` and `Zone` give one answer; an
/// identifier with no cell behind it is refused by every method, and a depth past the
/// finest level by the four that take one; the count, the first sub-zone, the zone at an
/// index and the order agree with one another; an entry is the null zone or a zone of the
/// level asked for; an index from the count upwards is refused by its name; and a list
/// above the limit is refused, never built. The engine is not asked: these are not
/// identifiers it may be handed.
///
/// What is left out is a matter of cost alone, and nothing is asked of a clock. The order
/// is built at depth 1 for every zone, at depth 2 for one in eight and at depth 3 for one
/// in thirty-two, and above the limit, where it is refused. The last index is asked to
/// depth 8, since the search reckons every scanline before the one it wants. The index of
/// another identifier is asked at every pair, whatever the depth between the two, and is
/// refused where the order is longer than a list may be; that of every entry of an order
/// at depth 1 is asked at one zone in eight, and is a place at which the order holds the
/// entry: the first, at every level coarser than 15 below an identifier that is the zone
/// found at its own centroid, and otherwise the first, a later one, or none, counted
/// exactly. `AnyGrid` and `Zone` are set beside the typed grid at one
/// identifier in four, and at every one that is refused.
pub fn the_sub_zone_methods_do_not_panic_on_hostile_input<S: Subject<T = HexA7, I = Z7>>() {
    the_limit_on_a_list_of_sub_zones_is_its_ceiling();
    let g = S::grid();
    let any = rs4dggs::get_grid(S::NAME).unwrap();
    let ids = super::geometry::hostile_identifiers::<S>();
    let finest = g.max_resolution();
    let limit = rs4dggs::grid::max_materialised_sub_zones();
    let invalid = |r: Result<(), Error>| matches!(r, Err(Error::InvalidZone(_)));
    // An identifier with a cell behind it, which `sub_zone_index` takes; and one that is a
    // zone of a level up to the finest besides, which the other four methods take: an
    // identifier of twenty digits is the first and not the second.
    let has_a_cell = |id: ZoneId| g.sub_zone_index(id, id).is_ok();
    let is_zone = |id: ZoneId| g.count_sub_zones(id, 0).is_ok();
    let real = g.zone_from_text("0064156").unwrap().id();
    let (mut zones_asked, mut refused, mut orders, mut entries, mut indices) =
        (0usize, 0usize, 0usize, 0usize, 0usize);
    // The indices refused because the order is longer than a list may be.
    let mut too_long = 0usize;
    // The entries whose index was asked; those with none, below a zone and below an
    // identifier that is not the zone found at its own centroid; those with the index of a
    // later place.
    let (mut entries_asked, mut without_index, mut without_index_of_no_zone, mut at_a_later_place) =
        (0usize, 0usize, 0usize, 0usize);
    // An identifier that the grid finds again at its own centroid, at its own level. One
    // that it does not is taken by the sub-zone methods as though it were a zone, as the
    // engine takes it: the child that a pentagon lacks is one.
    let is_its_own_zone = |id: ZoneId| {
        let c = g.centroid(id);
        g.zone_from_geo(c.lat, c.lon, g.resolution(id))
            .is_ok_and(|zone| zone.id() == id)
    };
    let mut previous = ZoneId::NULL;
    for (n, id) in ids.into_iter().enumerate() {
        let id = ZoneId(id);
        let z = g.zone(id);
        let at = format!("{:#018x}", id.0);
        let every_surface = n % 4 == 0 || !is_zone(id);

        // The index of a sub-zone: of the zone in itself, of the null zone, of a real zone
        // and of the identifier drawn before, at whatever depth the one lies below the
        // other: an order longer than a list may be is refused, and no other is built.
        for sub in [id, ZoneId::NULL, real, previous] {
            let at = format!("{at} and {:#018x}", sub.0);
            let index = g.sub_zone_index(id, sub);
            if every_surface {
                assert_eq!(any.sub_zone_index(id, sub), index, "{at}");
                assert_eq!(z.sub_zone_index(&g.zone(sub)), index, "{at}");
            }
            if !has_a_cell(id) || !has_a_cell(sub) {
                assert!(invalid(index.map(drop)), "{at}");
            } else if sub == id {
                assert_eq!(index, Ok(Some(0)), "{at}");
            } else if let Ok(Some(i)) = index {
                assert!(is_zone(id) && is_zone(sub), "{at}");
                let depth = g.resolution(sub) - g.resolution(id);
                assert_eq!(
                    g.sub_zone_at_index(id, depth, i),
                    Ok(sub),
                    "{at}: the zone at the index answered"
                );
                indices += 1;
            } else if let Err(e) = index {
                // Refused for the length of the order alone, where both are zones.
                let depth = g.resolution(sub).checked_sub(g.resolution(id));
                let count = depth.and_then(|depth| g.count_sub_zones(id, depth).ok());
                assert!(
                    count.is_some_and(
                        |count| count > limit && e == Error::TooManySubZones { count, limit }
                    ) || !is_zone(id)
                        || !is_zone(sub),
                    "{at}: {e}"
                );
                too_long += usize::from(count.is_some_and(|count| count > limit));
            }
        }
        previous = id;

        if !is_zone(id) {
            for depth in [0u8, 1, 2, 19, 20, 255] {
                assert!(invalid(g.count_sub_zones(id, depth).map(drop)), "{at}");
                assert!(invalid(any.count_sub_zones(id, depth).map(drop)), "{at}");
                assert!(invalid(z.count_sub_zones(depth).map(drop)), "{at}");
                assert!(invalid(g.first_sub_zone(id, depth).map(drop)), "{at}");
                assert!(invalid(any.first_sub_zone(id, depth).map(drop)), "{at}");
                assert!(invalid(z.first_sub_zone(depth).map(drop)), "{at}");
                assert!(invalid(g.sub_zones(id, depth).map(drop)), "{at}");
                assert!(invalid(any.sub_zones(id, depth).map(drop)), "{at}");
                assert!(invalid(z.sub_zones(depth).map(drop)), "{at}");
                for index in [0, 1, u64::MAX] {
                    assert!(
                        invalid(g.sub_zone_at_index(id, depth, index).map(drop)),
                        "{at}"
                    );
                    assert!(
                        invalid(any.sub_zone_at_index(id, depth, index).map(drop)),
                        "{at}"
                    );
                    assert!(invalid(z.sub_zone_at_index(depth, index).map(drop)), "{at}");
                }
            }
            refused += 1;
            continue;
        }

        let level = g.resolution(id);
        assert!(level <= finest, "{at}: a zone of level {level}");
        let deepest = finest - level;
        // An entry is the null zone or a zone of the level asked for.
        let is_entry = |sub: ZoneId, depth: u8| {
            sub == ZoneId::NULL || (is_zone(sub) && g.resolution(sub) == level + depth)
        };
        for depth in [0u8, 1, 2, 3, 8, deepest, deepest.saturating_add(1), 255] {
            let at = format!("{at} at depth {depth}");
            let count = g.count_sub_zones(id, depth);
            let first = g.first_sub_zone(id, depth);
            if every_surface {
                assert_eq!(any.count_sub_zones(id, depth), count, "{at}");
                assert_eq!(z.count_sub_zones(depth), count, "{at}");
                assert_eq!(any.first_sub_zone(id, depth), first, "{at}");
                assert_eq!(z.first_sub_zone(depth).map(|s| s.id()), first, "{at}");
            }
            if depth > deepest {
                assert!(invalid(count.map(drop)), "{at}");
                assert!(invalid(first.map(drop)), "{at}");
                assert!(invalid(g.sub_zones(id, depth).map(drop)), "{at}");
                assert!(invalid(g.sub_zone_at_index(id, depth, 0).map(drop)), "{at}");
                continue;
            }
            let (count, first) = (count.unwrap(), first.unwrap());
            assert!(count >= 1 && is_entry(first, depth), "{at}");
            assert_eq!(depth == 0, count == 1, "{at}");

            let mut asked = vec![0, 1, count, u64::MAX, 1 << 63];
            if depth <= 8 {
                asked.push(count - 1);
            }
            for index in asked {
                let sub = g.sub_zone_at_index(id, depth, index);
                if every_surface {
                    assert_eq!(any.sub_zone_at_index(id, depth, index), sub, "{at}");
                    assert_eq!(
                        z.sub_zone_at_index(depth, index).map(|s| s.id()),
                        sub,
                        "{at}"
                    );
                }
                if index >= count {
                    assert_eq!(sub, Err(Error::IndexOutOfRange { index, count }), "{at}");
                } else {
                    let sub = sub.unwrap();
                    assert!(is_entry(sub, depth), "{at}, index {index}");
                    assert!(index != 0 || sub == first, "{at}, index 0");
                }
            }

            let built = match depth {
                0 | 1 => true,
                2 => n % 8 == 0,
                3 => n % 32 == 0,
                _ => count > limit,
            };
            if !built {
                continue;
            }
            let order = g.sub_zones(id, depth);
            if every_surface {
                assert_eq!(any.sub_zones(id, depth), order, "{at}");
                assert_eq!(
                    z.sub_zones(depth)
                        .map(|subs| subs.iter().map(|s| s.id()).collect::<Vec<_>>()),
                    order,
                    "{at}"
                );
            }
            if count > limit {
                assert_eq!(order, Err(Error::TooManySubZones { count, limit }), "{at}");
                continue;
            }
            let order = order.unwrap();
            assert_eq!(order.len() as u64, count, "{at}");
            assert_eq!(order[0], first, "{at}");
            assert_eq!(
                g.sub_zone_at_index(id, depth, count - 1),
                Ok(order[order.len() - 1]),
                "{at}"
            );
            if depth == 0 {
                assert_eq!(order, vec![id], "{at}");
            }
            for (i, &sub) in order.iter().enumerate() {
                assert!(is_entry(sub, depth), "{at}, entry {i}");
                if depth == 1 && n % 8 == 0 && sub != ZoneId::NULL {
                    // The index of an entry is a place at which the order holds it: the
                    // first, and so the entry's own where the order names it once. From
                    // level 15 an entry may have no index, or that of a later place; and
                    // so may an entry of the order of an identifier that is not the zone
                    // found at its own centroid, whose order is that of no zone.
                    let first = order.iter().position(|&s| s == sub).map(|p| p as u64);
                    let index = g.sub_zone_index(id, sub).unwrap();
                    entries_asked += 1;
                    if let Some(j) = index {
                        assert_eq!(order[j as usize], sub, "{at}, entry {i}");
                    }
                    if index != first {
                        assert!(
                            level + depth >= 15 || !is_its_own_zone(id),
                            "{at}, entry {i}: the index is {index:?} and the first place \
                             {first:?}, at a level coarser than 15, of a zone that is found \
                             at its own centroid"
                        );
                        match (index, is_its_own_zone(id)) {
                            (None, true) => without_index += 1,
                            (None, false) => without_index_of_no_zone += 1,
                            (Some(_), _) => at_a_later_place += 1,
                        }
                    }
                }
            }
            orders += 1;
            entries += order.len();
        }
        zones_asked += 1;
    }
    eprintln!(
        "hostile sub-zone requests on {}: {zones_asked} identifiers with sub-zones and \
         {refused} refused; {orders} orders built, of {entries} entries; {indices} indices \
         answered for another identifier and {too_long} refused for the length of the \
         order; of {entries_asked} entries asked for their index \
         {without_index} have none, {without_index_of_no_zone} have none below an identifier \
         that is not the zone at its own centroid, and {at_a_later_place} have that of a \
         later place",
        S::NAME
    );
    // The identifiers are the same on every grid, and the order of a zone is made before
    // any projection: the counts are the same on the three, and held exactly.
    assert_eq!(
        [zones_asked, refused, orders, entries, indices],
        [2_166, 844, 4_682, 69_470, 46],
        "the identifiers with sub-zones that the sweep met; those refused; the orders built; \
         the entries in them; and the indices answered for another identifier (seed \
         {SEED:#x})"
    );
    assert_eq!(
        [
            entries_asked,
            without_index,
            without_index_of_no_zone,
            at_a_later_place
        ],
        [3_307, 0, 26, 0],
        "the entries asked for their index; those that have none; those that have none \
         below an identifier that is not the zone found at its own centroid; and those that \
         have the index of a later place of the same zone (seed {SEED:#x})"
    );
}
