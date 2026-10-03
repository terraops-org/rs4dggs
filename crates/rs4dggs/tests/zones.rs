//! Every zone of a level, in the order of DGGAL's `listZones`, and the same sequence entered
//! after any of its zones.
//!
//! What needs no engine is tested in every build; the comparison with the engine's own list of
//! the whole world is in the module `engine`, under the `oracle` feature.
use rs4dggs::indexings::{I3h, Z7};
use rs4dggs::topologies::{HexA3, HexA7};
use rs4dggs::{Address, AnyGrid, Error, Grid, Indexing, Projection, Topology, ZoneId};
use std::collections::HashSet;

const APERTURE_3: [&str; 3] = ["ISEA3H", "IVEA3H", "RTEA3H"];
const APERTURE_7: [&str; 3] = ["IGEO7", "IVEA7H", "RTEA7H"];

fn grid(name: &str) -> AnyGrid {
    rs4dggs::get_grid(name).unwrap()
}

/// Every zone of the level, whole.
fn whole(grid: AnyGrid, level: u8) -> Vec<ZoneId> {
    grid.zones(level).unwrap().collect()
}

/// The zones of the level taken in pages of `size`, each page entered after the last zone of
/// the one before, until a page comes back empty. `count` is the number of zones the level has:
/// a walk that yields more is stopped there, so that a cursor that does not advance fails the
/// comparison instead of running without end.
fn paged(grid: AnyGrid, level: u8, size: usize, count: usize) -> Vec<ZoneId> {
    let mut out: Vec<ZoneId> = grid.zones(level).unwrap().take(size).collect();
    while let Some(&last) = out.last() {
        let page: Vec<ZoneId> = grid
            .zones(level)
            .unwrap()
            .after(last)
            .unwrap()
            .take(size)
            .collect();
        if page.is_empty() || out.len() > count {
            break;
        }
        out.extend(page);
    }
    out
}

/// The zones of the lattice of `level` from the cell at `row` and `col` of the rhombus `root`
/// to the end of the level, read cell by cell from the topology, each with its key: what the
/// iterator must give, taken without it. `cells` bounds how many cells of rhombi are read.
fn from_the_lattice<T: Topology, I: Indexing>(
    level: u8,
    (root, row, col): (u8, u64, u64),
    cells: usize,
) -> Vec<(ZoneId, u64)> {
    let edge = T::lattice_edge(level);
    let mut out = Vec::new();
    let (mut root, mut row, mut col) = (root, row, col);
    for _ in 0..cells {
        if root > 9 {
            break;
        }
        out.extend(T::cell_zones(level, root, row, col));
        col += 1;
        if col == edge {
            (row, col) = (row + 1, 0);
        }
        if row == edge {
            (root, row) = (root + 1, 0);
        }
    }
    if root > 9 {
        out.extend(T::polar_zones(level, 10));
        out.extend(T::polar_zones(level, 11));
    }
    out.into_iter()
        .map(|(a, key)| (I::encode(&a), key))
        .collect()
}

/// The key of a zone in the order of its level, as the lattice gives it.
fn key<T: Topology, I: Indexing>(zone: ZoneId) -> u64 {
    T::locate(&I::decode(zone))
        .unwrap_or_else(|| panic!("{:#018x} is no zone of the lattice", zone.0))
        .3
}

/// Each of `zones` is a zone of `level` that the grid reads back from its own text, their keys
/// ascend strictly, and none comes twice.
fn are_zones_of_the_level_in_order<T: Topology, I: Indexing>(
    grid: AnyGrid,
    level: u8,
    zones: &[ZoneId],
) {
    for &z in zones {
        assert_eq!(grid.resolution(z), level, "{}", grid.text_id(z));
        assert_eq!(grid.zone_from_text(&grid.text_id(z)), Ok(z));
        assert!(grid.extent(z).is_some(), "{}", grid.text_id(z));
    }
    let keys: Vec<u64> = zones.iter().map(|&z| key::<T, I>(z)).collect();
    assert!(
        keys.windows(2).all(|w| w[0] < w[1]),
        "{} level {level}: the keys do not ascend",
        grid.name()
    );
    let distinct: HashSet<ZoneId> = zones.iter().copied().collect();
    assert_eq!(distinct.len(), zones.len(), "{} level {level}", grid.name());
}

/// The enumeration is refused beyond 14 on the aperture-7 grids, whose last level is 19, and
/// beyond the last level, 33, on the aperture-3 grids; the error carries the limit.
#[test]
fn a_level_beyond_the_limit_of_the_enumeration_is_refused() {
    for name in APERTURE_7 {
        let g = grid(name);
        assert_eq!(g.max_box_level(), 14, "{name}");
        assert_eq!(g.max_resolution(), 19, "{name}");
        assert!(g.zones(14).is_ok(), "{name}");
        for level in [15, 16, 19, 20, 255] {
            assert_eq!(
                g.zones(level).err(),
                Some(Error::ResolutionOutOfRange {
                    res: level,
                    max: 14
                }),
                "{name} level {level}"
            );
        }
    }
    for name in APERTURE_3 {
        let g = grid(name);
        assert_eq!(g.max_box_level(), 33, "{name}");
        assert_eq!(g.max_resolution(), 33, "{name}");
        assert!(g.zones(33).is_ok(), "{name}");
        for level in [34, 255] {
            assert_eq!(
                g.zones(level).err(),
                Some(Error::ResolutionOutOfRange {
                    res: level,
                    max: 33
                }),
                "{name} level {level}"
            );
        }
    }
    // The grid taken by its type answers as the handle does.
    assert_eq!(rs4dggs::igeo7().max_box_level(), 14);
    assert_eq!(rs4dggs::rtea3h().max_box_level(), 33);
    assert_eq!(
        rs4dggs::ivea7h().zones(15).err(),
        Some(Error::ResolutionOutOfRange { res: 15, max: 14 })
    );
    assert_eq!(
        rs4dggs::isea3h().zones(34).err(),
        Some(Error::ResolutionOutOfRange { res: 34, max: 33 })
    );
    assert_eq!(
        Error::ResolutionOutOfRange { res: 15, max: 14 }.to_string(),
        "resolution 15 out of range 0 to 14"
    );
}

/// The grid taken by its type and the handle on it give one sequence, whole and entered after
/// a zone, and it holds every zone of the level once.
fn the_grid_and_its_handle_agree<P: Projection, T: Topology, I: Indexing>(
    typed: &Grid<P, T, I>,
    levels: std::ops::RangeInclusive<u8>,
) {
    let any = grid(typed.name());
    for level in levels {
        let by_type: Vec<ZoneId> = typed.zones(level).unwrap().collect();
        let by_handle = whole(any, level);
        assert_eq!(by_type, by_handle, "{} level {level}", typed.name());
        assert_eq!(
            by_type.len() as u64,
            typed.count_zones(level).unwrap(),
            "{} level {level}",
            typed.name()
        );
        are_zones_of_the_level_in_order::<T, I>(any, level, &by_type);
        let middle = by_type.len() / 2;
        let after_by_type: Vec<ZoneId> = typed
            .zones(level)
            .unwrap()
            .after(by_type[middle])
            .unwrap()
            .collect();
        let after_by_handle: Vec<ZoneId> = any
            .zones(level)
            .unwrap()
            .after(by_type[middle])
            .unwrap()
            .collect();
        assert_eq!(after_by_type, by_type[middle + 1..]);
        assert_eq!(after_by_handle, by_type[middle + 1..]);
    }
}

#[test]
fn the_grid_and_its_handle_give_one_sequence_of_every_zone_of_a_level() {
    the_grid_and_its_handle_agree(rs4dggs::isea3h(), 0..=5);
    the_grid_and_its_handle_agree(rs4dggs::ivea3h(), 0..=5);
    the_grid_and_its_handle_agree(rs4dggs::rtea3h(), 0..=5);
    the_grid_and_its_handle_agree(rs4dggs::igeo7(), 0..=3);
    the_grid_and_its_handle_agree(rs4dggs::ivea7h(), 0..=3);
    the_grid_and_its_handle_agree(rs4dggs::rtea7h(), 0..=3);
}

/// On the aperture-3 grids the order is the ascending order of the identifiers; on the
/// aperture-7 grids it is not, from level 1, and level 0 comes in the order of the rhombi.
#[test]
fn the_order_is_that_of_the_identifiers_on_aperture_3_alone() {
    for name in APERTURE_3 {
        for level in 0..=6 {
            let zones = whole(grid(name), level);
            assert!(
                zones.windows(2).all(|w| w[0] < w[1]),
                "{name} level {level}"
            );
        }
    }
    for name in APERTURE_7 {
        let g = grid(name);
        let level_0: Vec<String> = whole(g, 0).iter().map(|&z| g.text_id(z)).collect();
        assert_eq!(
            level_0,
            [
                "01", "06", "02", "07", "03", "08", "04", "09", "05", "10", "00", "11"
            ],
            "{name}"
        );
        for level in 1..=3 {
            let zones = whole(g, level);
            assert!(
                zones.windows(2).any(|w| w[0] > w[1]),
                "{name} level {level}"
            );
        }
    }
}

/// The sequence cut into pages of 1, 7 and 97, each entered after the last zone of the one
/// before, is the whole sequence; the page after the last zone is empty.
#[test]
fn pages_entered_after_the_last_zone_of_the_one_before_make_the_whole() {
    let mut compared = 0;
    for (names, levels) in [(APERTURE_3, 0..=5u8), (APERTURE_7, 0..=3u8)] {
        for name in names {
            let g = grid(name);
            for level in levels.clone() {
                let all = whole(g, level);
                for size in [1, 7, 97] {
                    let cut = paged(g, level, size, all.len());
                    assert!(
                        cut == all,
                        "{name} level {level}, pages of {size}: {} zones for {}, the first \
                         departure at entry {:?}",
                        cut.len(),
                        all.len(),
                        cut.iter().zip(&all).position(|(a, b)| a != b)
                    );
                    compared += all.len();
                }
            }
        }
    }
    // Levels 0 to 5 of aperture 3 hold 3,652 zones and levels 0 to 3 of aperture 7 hold 4,008:
    // three grids of each, and three sizes of page.
    assert_eq!(compared, 3 * 3 * (3_652 + 4_008));
}

/// After each zone of a level come the zones that follow it, and after the last none; and an
/// iterator that has ended stays ended.
#[test]
fn after_any_zone_come_the_zones_that_follow_it() {
    for (name, level) in [
        ("ISEA3H", 0),
        ("IVEA3H", 2),
        ("RTEA3H", 3),
        ("IGEO7", 0),
        ("IVEA7H", 1),
        ("RTEA7H", 2),
    ] {
        let g = grid(name);
        let all = whole(g, level);
        for (i, &z) in all.iter().enumerate() {
            let rest: Vec<ZoneId> = g.zones(level).unwrap().after(z).unwrap().collect();
            assert_eq!(rest, all[i + 1..], "{name} level {level}, after entry {i}");
        }
        let mut ended = g.zones(level).unwrap().after(*all.last().unwrap()).unwrap();
        for _ in 0..3 {
            assert_eq!(ended.next(), None, "{name} level {level}");
        }
        let mut consumed = g.zones(level).unwrap();
        assert_eq!(consumed.by_ref().count(), all.len());
        for _ in 0..3 {
            assert_eq!(consumed.next(), None, "{name} level {level}");
        }
        // An iterator partly consumed, or ended, is entered anew: what it gave is not counted.
        let again: Vec<ZoneId> = consumed.after(all[0]).unwrap().collect();
        assert_eq!(again, all[1..], "{name} level {level}");
    }
}

/// A zone to enter after is accepted exactly where the enumeration of the level yields it.
/// `candidates` are identifiers of every provenance; the count of those refused is returned.
fn is_accepted_exactly_where_it_is_yielded(
    g: AnyGrid,
    level: u8,
    candidates: impl IntoIterator<Item = ZoneId>,
) -> usize {
    let all: HashSet<ZoneId> = whole(g, level).into_iter().collect();
    let mut refused = 0;
    for z in candidates {
        let entered = g.zones(level).unwrap().after(z);
        assert_eq!(
            entered.is_ok(),
            all.contains(&z),
            "{} level {level}: {:#018x}",
            g.name(),
            z.0
        );
        if let Err(e) = entered {
            assert!(
                matches!(&e, Error::InvalidZone(message) if message.len() < 200),
                "{e:?}"
            );
            refused += 1;
        }
    }
    refused
}

/// Identifiers that no grid gave out: the null zone, both ends of the range, and patterns of
/// bits, with each of the sixteen values of the highest four bits among them.
fn hostile() -> Vec<ZoneId> {
    let mut out = vec![
        ZoneId::NULL,
        ZoneId(0),
        ZoneId(1),
        ZoneId(u64::MAX - 1),
        ZoneId(0xdead_beef_dead_beef),
        ZoneId(0x5555_5555_5555_5555),
        ZoneId(0xaaaa_aaaa_aaaa_aaaa),
        ZoneId(0x8000_0000_0000_0000),
    ];
    for top in 0..16u64 {
        out.push(ZoneId(top << 60));
        out.push(ZoneId((top << 60) | 0x0123_4567_89ab_cdef));
        out.push(ZoneId((top << 60) | 0x0fff_ffff_ffff_ffff));
    }
    out
}

/// A zone of another level, the null zone, any other pattern of bits and the identifiers of
/// another grid's zones are refused as a zone to enter after, unless the level holds that very
/// identifier.
#[test]
fn a_zone_to_enter_after_that_the_level_does_not_hold_is_refused() {
    for (names, others, levels) in [
        (APERTURE_3, APERTURE_7, 0..=3u8),
        (APERTURE_7, APERTURE_3, 0..=2u8),
    ] {
        for name in names {
            let g = grid(name);
            for level in levels.clone() {
                // The null zone, by itself, and the zones of the levels beside this one.
                assert!(matches!(
                    g.zones(level).unwrap().after(ZoneId::NULL),
                    Err(Error::InvalidZone(_))
                ));
                let mut beside = whole(g, level + 1);
                if level > 0 {
                    beside.extend(whole(g, level - 1));
                }
                let count = beside.len();
                assert_eq!(
                    is_accepted_exactly_where_it_is_yielded(g, level, beside),
                    count,
                    "{name} level {level}: a zone of another level was accepted"
                );
                // Patterns of bits: none makes the call panic, and all are refused but those
                // that are a zone of the level, as twelve of them are at level 0 of aperture 7.
                let patterns = hostile();
                let count = patterns.len();
                assert!(
                    is_accepted_exactly_where_it_is_yielded(g, level, patterns) >= count - 12,
                    "{name} level {level}"
                );
                // The identifiers of the other aperture's zones, of this level and of others.
                let other = grid(others[0]);
                let foreign: Vec<ZoneId> = (0..=3).flat_map(|l| whole(other, l)).collect();
                let count = foreign.len();
                assert!(
                    is_accepted_exactly_where_it_is_yielded(g, level, foreign) >= count * 9 / 10,
                    "{name} level {level}"
                );
            }
        }
    }
}

/// Identifiers that decode to something and that the lattice does not hold are refused: on
/// aperture 7 the child a pentagon does not have, and an identifier with a digit written after
/// its last; on aperture 3 the sub-hexagons C and D of a polar root, which the engine reads
/// from text and never lists.
#[test]
fn an_identifier_that_the_lattice_does_not_hold_is_refused() {
    for name in APERTURE_7 {
        let g = grid(name);
        for (base, digits) in [
            (0, &[2][..]),
            (5, &[2]),
            (6, &[5]),
            (11, &[5]),
            (0, &[0, 2]),
            (11, &[0, 0, 5]),
        ] {
            let level = digits.len() as u8;
            let deleted = Z7::encode(&Address::new(base, digits));
            assert_eq!(g.resolution(deleted), level);
            assert!(
                matches!(
                    g.zones(level).unwrap().after(deleted),
                    Err(Error::InvalidZone(_))
                ),
                "{name}: base cell {base}, digits {digits:?}"
            );
        }
        // A zone of level 1 with the lowest digit of the padding overwritten: it decodes to
        // the same address, and it is not that zone's identifier.
        for &z in &whole(g, 1) {
            let stray = ZoneId(z.0 & !0b111);
            assert_ne!(stray, z);
            assert_eq!(Z7::decode(stray), Z7::decode(z));
            assert!(g.zones(1).unwrap().after(z).is_ok());
            assert!(
                matches!(g.zones(1).unwrap().after(stray), Err(Error::InvalidZone(_))),
                "{name}: {:#018x}",
                stray.0
            );
        }
    }
    for name in APERTURE_3 {
        let g = grid(name);
        for level in [1, 3, 33] {
            // The last two zones of an odd level are the sub-hexagons B of the two polar roots;
            // C and D are the two identifiers after each.
            let polar = from_the_lattice::<HexA3, I3h>(level, (10, 0, 0), 0);
            assert_eq!(polar.len(), 2);
            for (b, _) in polar {
                assert!(g.text_id(b).ends_with("-B"), "{}", g.text_id(b));
                assert!(g.zones(level).unwrap().after(b).is_ok());
                for (step, letter) in [(1, "-C"), (2, "-D")] {
                    let other = ZoneId(b.0 + step);
                    assert!(g.text_id(other).ends_with(letter), "{}", g.text_id(other));
                    assert!(
                        matches!(
                            g.zones(level).unwrap().after(other),
                            Err(Error::InvalidZone(_))
                        ),
                        "{name}: {}",
                        g.text_id(other)
                    );
                }
            }
        }
    }
}

/// At the finest level enumerated, where no list of the level could be held: the first thousand
/// zones; the thousand after a zone in the middle of a rhombus; the passage from the end of a
/// row to the next row, from the end of a rhombus to the next rhombus, and from the last rhombus
/// to the polar roots and the end. Each is what the lattice gives cell by cell, and each zone is
/// a zone of the level, the keys ascending, none twice.
fn the_finest_level_is_walked_from_any_place<T: Topology, I: Indexing>(name: &str, level: u8) {
    let g = grid(name);
    assert_eq!(level, g.max_box_level());
    let edge = T::lattice_edge(level);
    assert!(edge > 2_000);

    let first: Vec<ZoneId> = g.zones(level).unwrap().take(1_000).collect();
    assert_eq!(first.len(), 1_000);
    are_zones_of_the_level_in_order::<T, I>(g, level, &first);
    let expected = from_the_lattice::<T, I>(level, (0, 0, 0), 1_000);
    assert!(first.iter().eq(expected.iter().map(|(z, _)| z).take(1_000)));

    // After the first zone of a cell, and after the last, for each place.
    let last = edge - 1;
    for (place, cells, to_the_end) in [
        ((4, edge / 2, edge / 2), 1_001, false),
        ((2, 5, last - 2), 1_001, false),
        ((3, last, last - 2), 1_001, false),
        ((9, last, last - 100), 101, true),
    ] {
        let hosted = from_the_lattice::<T, I>(level, place, cells);
        for at in [0, T::cell_zones(level, place.0, place.1, place.2).len() - 1] {
            let (cursor, cursor_key) = hosted[at];
            let page: Vec<ZoneId> = g
                .zones(level)
                .unwrap()
                .after(cursor)
                .unwrap()
                .take(1_000)
                .collect();
            let expected: Vec<ZoneId> = hosted[at + 1..]
                .iter()
                .map(|(z, _)| *z)
                .take(1_000)
                .collect();
            assert_eq!(page, expected, "{name} level {level}, after {place:?}");
            if to_the_end {
                // Fewer than a thousand remain: the walk passed the polar roots and ended.
                assert_eq!(page.len(), hosted.len() - at - 1);
                assert!(page.len() >= 100);
                let polar = from_the_lattice::<T, I>(level, (10, 0, 0), 0);
                assert!(!polar.is_empty());
                assert!(page.ends_with(&polar.iter().map(|(z, _)| *z).collect::<Vec<_>>()));
            } else {
                assert_eq!(page.len(), 1_000);
            }
            are_zones_of_the_level_in_order::<T, I>(g, level, &page);
            assert!(page.iter().all(|&z| key::<T, I>(z) > cursor_key));
        }
    }
}

#[test]
fn level_33_of_aperture_3_is_walked_from_any_place() {
    for name in APERTURE_3 {
        the_finest_level_is_walked_from_any_place::<HexA3, I3h>(name, 33);
    }
}

#[test]
fn level_14_of_aperture_7_is_walked_from_any_place() {
    for name in APERTURE_7 {
        the_finest_level_is_walked_from_any_place::<HexA7, Z7>(name, 14);
    }
}

/// The finest odd level of aperture 7, where a cell hosts seven zones and a pentagon's six, and
/// level 32 of aperture 3, where a cell hosts one.
#[test]
fn the_level_above_the_finest_is_walked_from_the_middle_of_a_rhombus() {
    for (name, level) in [("ISEA3H", 32), ("RTEA3H", 31)] {
        the_middle_of_a_rhombus::<HexA3, I3h>(name, level);
    }
    for (name, level) in [("IGEO7", 13), ("IVEA7H", 12)] {
        the_middle_of_a_rhombus::<HexA7, Z7>(name, level);
    }
}

fn the_middle_of_a_rhombus<T: Topology, I: Indexing>(name: &str, level: u8) {
    let g = grid(name);
    let edge = T::lattice_edge(level);
    let hosted = from_the_lattice::<T, I>(level, (7, edge / 2, edge / 3), 1_003);
    for at in [0, 1, 2] {
        let page: Vec<ZoneId> = g
            .zones(level)
            .unwrap()
            .after(hosted[at].0)
            .unwrap()
            .take(1_000)
            .collect();
        let expected: Vec<ZoneId> = hosted[at + 1..]
            .iter()
            .map(|(z, _)| *z)
            .take(1_000)
            .collect();
        assert_eq!(page, expected, "{name} level {level}");
        are_zones_of_the_level_in_order::<T, I>(g, level, &page);
    }
}

/// The enumeration against the live engine.
#[cfg(feature = "oracle")]
mod engine {
    use super::*;
    use dggal_oracle::{count_zones, list_zones};

    /// The crate's name of a grid, and the engine's.
    const GRIDS_3: [(&str, &str); 3] = [
        ("ISEA3H", dggal_oracle::ISEA3H),
        ("IVEA3H", dggal_oracle::IVEA3H),
        ("RTEA3H", dggal_oracle::RTEA3H),
    ];
    const GRIDS_7: [(&str, &str); 3] = [
        ("IGEO7", dggal_oracle::IGEO7),
        ("IVEA7H", dggal_oracle::IVEA7H),
        ("RTEA7H", dggal_oracle::RTEA7H),
    ];

    /// The engine's list of the whole world at `level`, which holds as many zones as the engine
    /// counts there.
    fn the_engines(engine: &str, level: u8) -> Vec<ZoneId> {
        let listed = list_zones(engine, i32::from(level), None);
        assert_eq!(
            listed.len() as u64,
            count_zones(engine, i32::from(level)),
            "{engine} level {level}"
        );
        listed.into_iter().map(ZoneId).collect()
    }

    fn is_the_engines(what: &str, ours: &[ZoneId], theirs: &[ZoneId]) {
        assert!(
            ours == theirs,
            "{what}: {} zones depart from the engine's list of {} at entry {:?}",
            ours.len(),
            theirs.len(),
            ours.iter().zip(theirs).position(|(a, b)| a != b)
        );
    }

    /// Every zone of a level is the engine's list of the whole world, sequence for sequence
    /// and never sorted: levels 0 to 8 on the three grids of aperture 3 and 0 to 5 on the three
    /// of aperture 7, by the grid's type through its handle. The same sequence cut into pages,
    /// each entered after the last zone of the one before, is the engine's list too: pages of
    /// 97 at every level, of 7 to the level before the last, of 1 two levels before it.
    #[test]
    fn every_zone_of_a_level_is_the_engines_whole_world_in_its_order() {
        let mut compared = 0;
        for (grids, last) in [(GRIDS_3, 8u8), (GRIDS_7, 5u8)] {
            for (name, engine) in grids {
                let g = grid(name);
                for level in 0..=last {
                    let theirs = the_engines(engine, level);
                    let what = format!("{name} level {level}");
                    is_the_engines(&what, &whole(g, level), &theirs);
                    compared += theirs.len();
                    for (size, to) in [(97, last), (7, last - 1), (1, last - 2)] {
                        if level <= to {
                            let cut = paged(g, level, size, theirs.len());
                            is_the_engines(&format!("{what}, pages of {size}"), &cut, &theirs);
                        }
                    }
                }
            }
        }
        let zones = |aperture: usize, last: u32| {
            (0..=last).map(|l| 10 * aperture.pow(l) + 2).sum::<usize>()
        };
        assert_eq!(compared, 3 * (zones(3, 8) + zones(7, 5)));
    }
}
