//! The zones of a level as a lattice: ten root rhombi of cells and two polar roots. Walked
//! rhombus by rhombus, row by row and cell by cell, and then the two polar roots, the zones the
//! cells host are every zone of the level, each once, in the order of DGGAL's `listZones`.
//!
//! What needs no engine is tested in every build; the comparison with the engine's own list is
//! in the module `engine`, under the `oracle` feature.
use crate::interfaces::sealed::{Private, TopologyPlumbing};
use rs4dggs::indexings::{I3h, Z7};
use rs4dggs::topologies::{HexA3, HexA7};
use rs4dggs::{Address, Indexing, Topology};

/// A cell of the lattice: the root, the row and the column, both nought under a polar root.
type Cell = (u8, u64, u64);
/// The zones a cell hosts, each with its key.
type Hosted = Vec<(Address, u64)>;

/// The cells of the lattice of `level` in the order of the walk, each with the zones it hosts.
fn cells<T: Topology>(level: u8) -> Vec<(Cell, Hosted)> {
    let edge = T::lattice_edge(Private, level);
    let mut out = Vec::new();
    for root in 0..10 {
        for row in 0..edge {
            for col in 0..edge {
                out.push((
                    (root, row, col),
                    T::cell_zones(Private, level, root, row, col),
                ));
            }
        }
    }
    for root in 10..12 {
        out.push(((root, 0, 0), T::polar_zones(Private, level, root)));
    }
    out
}

/// The zones of `level` in the order of the walk, each with its key.
fn lattice<T: Topology>(level: u8) -> Hosted {
    cells::<T>(level)
        .into_iter()
        .flat_map(|(_, zones)| zones)
        .collect()
}

/// Every zone of the walk is located in the cell that gave it, with the key that cell gave it,
/// and that cell holds it once.
fn every_zone_is_located_in_its_own_cell<T: Topology>(level: u8) -> usize {
    let mut zones = 0;
    for (cell, hosted) in cells::<T>(level) {
        for (a, key) in &hosted {
            assert_eq!(
                T::locate(Private, a),
                Some((cell.0, cell.1, cell.2, *key)),
                "level {level}, {a:?} of the cell {cell:?}"
            );
            assert_eq!(
                hosted.iter().filter(|(b, _)| b == a).count(),
                1,
                "level {level}, {a:?} in the cell {cell:?}"
            );
            zones += 1;
        }
    }
    zones
}

/// How many times the cell in which `a` is located holds `a`: nought where it is located
/// nowhere, and one for a zone of the lattice.
fn times_hosted<T: Topology>(level: u8, a: &Address) -> usize {
    let Some((root, row, col, key)) = T::locate(Private, a) else {
        return 0;
    };
    let hosted = if root < 10 {
        T::cell_zones(Private, level, root, row, col)
    } else {
        T::polar_zones(Private, level, root)
    };
    hosted.iter().filter(|(b, k)| b == a && *k == key).count()
}

/// On aperture 3 the walk gives the zones of a level in the ascending order of their
/// identifiers, each its own key, and as many as the grid counts at the level.
#[test]
fn the_aperture_3_lattice_ascends_and_holds_as_many_zones_as_the_level_has() {
    let grid = rs4dggs::isea3h();
    for level in 0..=8 {
        let zones = lattice::<HexA3>(level);
        assert_eq!(
            zones.len() as u64,
            grid.count_zones(level).unwrap(),
            "level {level}"
        );
        for (a, key) in &zones {
            assert_eq!(I3h::encode(a).0, *key, "level {level}");
            assert_eq!(I3h::resolution(I3h::encode(a)), level);
        }
        assert!(zones.windows(2).all(|w| w[0].1 < w[1].1), "level {level}");
    }
}

/// At every level of aperture 3, the cells counted and one cell of each kind asked: ten rhombi
/// of one zone a cell at an even level and three at an odd one, and one zone a polar root.
#[test]
fn the_aperture_3_lattice_counts_the_zones_of_every_level() {
    let grid = rs4dggs::isea3h();
    for level in 0..=I3h::MAX_RESOLUTION {
        let edge = HexA3::lattice_edge(Private, level);
        let per_cell = if level % 2 == 0 { 1 } else { 3 };
        for (row, col) in [(0, 0), (edge - 1, 0), (0, edge - 1), (edge / 2, edge / 3)] {
            for root in 0..10 {
                assert_eq!(
                    HexA3::cell_zones(Private, level, root, row, col).len(),
                    per_cell,
                    "level {level}, root {root}, row {row}, column {col}"
                );
            }
        }
        let polar = HexA3::polar_zones(Private, level, 10).len()
            + HexA3::polar_zones(Private, level, 11).len();
        assert_eq!(polar, 2, "level {level}");
        assert_eq!(
            10 * edge * edge * per_cell as u64 + polar as u64,
            grid.count_zones(level).unwrap(),
            "level {level}"
        );
    }
}

/// `locate` is the inverse of `cell_zones` and `polar_zones` at every zone of levels 0 to 6 of
/// aperture 3, and at the cells of the corners at the finest two levels.
#[test]
fn an_aperture_3_zone_is_located_in_the_cell_that_hosts_it() {
    let mut zones = 0;
    for level in 0..=6 {
        zones += every_zone_is_located_in_its_own_cell::<HexA3>(level);
    }
    assert_eq!(
        zones,
        (0..=6).map(|l| 10 * 3usize.pow(l) + 2).sum::<usize>()
    );
    for level in [32, 33] {
        let last = HexA3::lattice_edge(Private, level) - 1;
        for (root, row, col) in [(0, 0, 0), (9, last, last), (4, last, 0), (5, 1, last)] {
            let hosted = HexA3::cell_zones(Private, level, root, row, col);
            assert!(!hosted.is_empty());
            for (a, key) in hosted {
                assert_eq!(HexA3::locate(Private, &a), Some((root, row, col, key)));
            }
        }
    }
}

/// Outside the lattice there is no cell and no zone: a level beyond the finest, a root that is
/// not of the kind asked, a row or a column beyond the side, an identifier that names no zone.
#[test]
fn the_aperture_3_lattice_answers_nothing_outside_itself() {
    for level in [34, 35, 66, 255] {
        assert_eq!(HexA3::lattice_edge(Private, level), 0, "level {level}");
        assert!(HexA3::cell_zones(Private, level, 0, 0, 0).is_empty());
        assert!(HexA3::polar_zones(Private, level, 10).is_empty());
    }
    let edge = HexA3::lattice_edge(Private, 5);
    assert_eq!(edge, 9);
    for (root, row, col) in [
        (10, 0, 0),
        (11, 0, 0),
        (12, 0, 0),
        (255, 0, 0),
        (0, edge, 0),
        (0, 0, edge),
        (0, u64::MAX, 0),
        (0, 0, u64::MAX),
    ] {
        assert!(
            HexA3::cell_zones(Private, 5, root, row, col).is_empty(),
            "{root} {row} {col}"
        );
    }
    for root in [0, 9, 12, 255] {
        assert!(
            HexA3::polar_zones(Private, 5, root).is_empty(),
            "root {root}"
        );
    }
    // The null zone, the sub-hexagon C of a polar root, a polar root with an index, an index
    // beyond the rhombus, a level beyond the finest, and an identifier with its top bit set.
    for id in [
        u64::MAX,
        (10 << 53) | 2,
        (1 << 57) | (10 << 53) | (1 << 2),
        (1 << 57) | (9 << 2),
        17 << 57,
        1 << 63,
    ] {
        assert_eq!(
            HexA3::locate(Private, &Address::new(id, &[])),
            None,
            "{id:#x}"
        );
    }
}

/// The address a Z7 text names.
fn z7(text: &str) -> Address {
    Z7::decode(Z7::from_text(text).unwrap())
}

/// The digit a pentagon's children never bear, by the base cell: 2 under the cells 0 to 5, 5
/// under 6 to 11. An address whose first digit that is not 0 is that one names the child a
/// pentagon does not have.
fn names_a_deleted_child(a: &Address) -> bool {
    let deleted = if a.base <= 5 { 2 } else { 5 };
    a.digits().iter().find(|&&d| d != 0) == Some(&deleted)
}

/// On aperture 7 the order is not that of the Z7 identifiers: level 0 comes rhombus by rhombus
/// and the polar roots last, which in base cells is this.
#[test]
fn the_aperture_7_base_cells_come_in_the_order_of_the_rhombi() {
    let texts: Vec<String> = lattice::<HexA7>(0)
        .iter()
        .map(|(a, _)| Z7::to_text(Z7::encode(a)))
        .collect();
    assert_eq!(texts.join(" "), "01 06 02 07 03 08 04 09 05 10 00 11");
}

/// At levels 0 to 4 of aperture 7 the walk gives as many zones as the grid counts, no zone
/// twice, each of the level and readable from its own text, and the keys ascend.
#[test]
fn the_aperture_7_lattice_ascends_and_holds_every_zone_of_the_level_once() {
    let grid = rs4dggs::igeo7();
    for level in 0..=4 {
        let zones = lattice::<HexA7>(level);
        assert_eq!(
            zones.len() as u64,
            grid.count_zones(level).unwrap(),
            "level {level}"
        );
        assert!(zones.windows(2).all(|w| w[0].1 < w[1].1), "level {level}");
        let mut ids = std::collections::HashSet::new();
        for (a, _) in &zones {
            let id = Z7::encode(a);
            assert_eq!(a.len(), usize::from(level), "{a:?}");
            assert_eq!(Z7::from_text(&Z7::to_text(id)), Ok(id), "{a:?}");
            assert!(ids.insert(id), "level {level}: {a:?} twice");
        }
    }
}

/// At every level of aperture 7 the cells counted, and cells of each kind asked: ten rhombi of
/// one zone a cell at an even level and seven at an odd one, save the ten cells of the
/// pentagons, which host six at an odd level; one zone a polar root, six at an odd level. The
/// sum is the count of the level, with nothing enumerated.
///
/// This is an identity on the cells asked, which lie inside the rhombi and at their pentagons,
/// and no census of the lattice. From level 16 some cells beside two edges of the icosahedron
/// host nothing, so that the lattice itself holds fewer zones than this sum there:
/// `from_level_16_cells_beside_two_edges_host_nothing`.
#[test]
fn the_regular_cells_of_the_aperture_7_lattice_host_what_the_count_of_a_level_assumes() {
    let grid = rs4dggs::igeo7();
    for level in 0..=Z7::MAX_RESOLUTION {
        let edge = HexA7::lattice_edge(Private, level);
        assert_eq!(edge, 7u64.pow(u32::from(level / 2)), "level {level}");
        let (per_cell, per_pentagon) = if level % 2 == 1 { (7, 6) } else { (1, 1) };
        for root in 0..10 {
            assert_eq!(
                HexA7::cell_zones(Private, level, root, 0, 0).len(),
                per_pentagon,
                "level {level}, the pentagon's cell of root {root}"
            );
            if edge == 1 {
                continue;
            }
            for (row, col) in [
                (edge / 2, edge / 3),
                (edge / 3, edge / 2),
                (edge / 5, 3 * edge / 5),
                (4 * edge / 7, 2 * edge / 7),
                (5 * edge / 7, 5 * edge / 7),
            ] {
                assert_eq!(
                    HexA7::cell_zones(Private, level, root, row, col).len(),
                    per_cell,
                    "level {level}, root {root}, row {row}, column {col}"
                );
            }
        }
        let polar = HexA7::polar_zones(Private, level, 10).len()
            + HexA7::polar_zones(Private, level, 11).len();
        assert_eq!(polar, 2 * per_pentagon, "level {level}");
        assert_eq!(
            (10 * edge * edge - 10) * per_cell as u64 + 10 * per_pentagon as u64 + polar as u64,
            grid.count_zones(level).unwrap(),
            "level {level}"
        );
    }
}

/// `locate` is the inverse of `cell_zones` and `polar_zones` at every zone of levels 0 to 4 of
/// aperture 7.
#[test]
fn an_aperture_7_zone_is_located_in_the_cell_that_hosts_it() {
    let mut zones = 0;
    for level in 0..=4 {
        zones += every_zone_is_located_in_its_own_cell::<HexA7>(level);
    }
    assert_eq!(
        zones,
        (0..=4).map(|l| 10 * 7usize.pow(l) + 2).sum::<usize>()
    );
}

/// The child a pentagon does not have is in no cell: at every level from 1 to 19, the cells of
/// the ten pentagons and the two polar roots host no address that names it, and at an odd
/// level, where the pentagon's cell would hold it, they host the six others. Nor is it located
/// anywhere.
#[test]
fn no_cell_hosts_the_child_a_pentagon_does_not_have() {
    for level in 1..=Z7::MAX_RESOLUTION {
        let mut hosted: Vec<Hosted> = (0..10)
            .map(|root| HexA7::cell_zones(Private, level, root, 0, 0))
            .collect();
        hosted.push(HexA7::polar_zones(Private, level, 10));
        hosted.push(HexA7::polar_zones(Private, level, 11));
        for zones in hosted {
            assert_eq!(zones.len(), if level % 2 == 1 { 6 } else { 1 });
            for (a, _) in zones {
                assert!(!names_a_deleted_child(&a), "level {level}: {a:?}");
                let id = Z7::encode(&a);
                assert_eq!(Z7::from_text(&Z7::to_text(id)), Ok(id), "{a:?}");
            }
        }
    }
    for level in 1..=4 {
        for (a, _) in lattice::<HexA7>(level) {
            assert!(!names_a_deleted_child(&a), "level {level}: {a:?}");
        }
    }
    for (base, digits) in [
        (0, [2].as_slice()),
        (5, &[0, 2]),
        (6, &[5]),
        (11, &[0, 0, 5, 1]),
        (3, &[0, 0, 0, 0, 2, 6, 6]),
    ] {
        assert_eq!(HexA7::locate(Private, &Address::new(base, digits)), None);
    }
}

/// Outside the lattice there is no cell and no zone: a level beyond the finest, a root that is
/// not of the kind asked, a row or a column beyond the side, an address that names no zone.
#[test]
fn the_aperture_7_lattice_answers_nothing_outside_itself() {
    for level in [20, 21, 40, 255] {
        assert_eq!(HexA7::lattice_edge(Private, level), 0, "level {level}");
        assert!(HexA7::cell_zones(Private, level, 0, 0, 0).is_empty());
        assert!(HexA7::polar_zones(Private, level, 10).is_empty());
    }
    let edge = HexA7::lattice_edge(Private, 5);
    assert_eq!(edge, 49);
    for (root, row, col) in [
        (10, 0, 0),
        (11, 0, 0),
        (12, 0, 0),
        (255, 0, 0),
        (0, edge, 0),
        (0, 0, edge),
        (0, u64::MAX, 0),
        (0, 0, u64::MAX),
    ] {
        assert!(
            HexA7::cell_zones(Private, 5, root, row, col).is_empty(),
            "{root} {row} {col}"
        );
    }
    for root in [0, 9, 12, 255] {
        assert!(
            HexA7::polar_zones(Private, 5, root).is_empty(),
            "root {root}"
        );
    }
    // A base cell the icosahedron does not have, a digit beyond the aperture, the terminator
    // as a digit, and the twentieth level, which the packing holds and the engine does not draw.
    for a in [
        Address::new(12, &[]),
        Address::new(15, &[1, 2]),
        Address::new(3, &[1, 8]),
        Address::new(3, &[1, 7, 1]),
        Address::new(3, &[1; 20]),
    ] {
        assert_eq!(HexA7::locate(Private, &a), None, "{a:?}");
    }
}

/// From level 15 the engine puts some zones along two edges of the icosahedron (from base cell
/// 0 to 1, over the north pole, and from 1 to 6) far from the cell their identifiers give: the
/// centroid of the first two below lies two cells from the cell that hosts them, that of the
/// third four cells from it, on the edge of the rhombus. The lattice hosts them where their
/// identifiers put them, each once; so it does the zones of those edges at which the engine's
/// hierarchy and sub-zone orders are known to fail.
#[test]
fn the_lattice_hosts_the_zones_that_the_engine_displaces_along_two_edges() {
    for (text, cell) in [
        ("0132323232323232316", (8, 2_882_400, 5_764_799)),
        ("0132323232323232351", (8, 2_882_401, 5_764_799)),
        ("0132323232323232620", (0, 4, 2_882_399)),
    ] {
        let a = z7(text);
        let located = HexA7::locate(Private, &a).map(|(root, row, col, _)| (root, row, col));
        assert_eq!(located, Some(cell), "{text}");
        assert_eq!(times_hosted::<HexA7>(17, &a), 1, "{text}");
    }
    for text in [
        "010004000400",
        "0100000000000004",
        "0100040004000400",
        "0000000000000000",
        "0000000000000005",
        "01000000000000003",
        "00000000000000005",
        "00052626050026015",
        "010000000000000042",
        "010000000000000065",
    ] {
        let a = z7(text);
        assert_eq!(times_hosted::<HexA7>(a.len() as u8, &a), 1, "{text}");
    }
}

/// Deep under the pentagon of base cell 0 some addresses read to a zone that the engine itself
/// names by another identifier: the first and the third below it draws where it draws the
/// pentagon `0000000000000000000`. The lattice hosts a zone under the one address its cell
/// gives it, and such an address is located nowhere.
const ADDRESSES_OF_ZONES_NAMED_OTHERWISE: [&str; 4] = [
    "0000000000000000120",
    "0000000000000000121",
    "0000000000000000130",
    "0000000000000000131",
];

#[test]
fn an_address_of_a_zone_that_is_named_otherwise_is_located_nowhere() {
    for text in ADDRESSES_OF_ZONES_NAMED_OTHERWISE {
        assert_eq!(HexA7::locate(Private, &z7(text)), None, "{text}");
    }
    assert_eq!(times_hosted::<HexA7>(17, &z7("0000000000000000000")), 1);
}

/// From level 16, beside the two edges of the icosahedron from base cell 0 to 1 and from 1 to
/// 6, there are cells whose zone has no identifier of its own: its conversion gives none, where
/// the engine prints one that it cannot read back, or gives the identifier of another zone.
/// The lattice hosts nothing there, and so from level 16 it holds fewer zones than the level
/// counts. Here, four cells of the first row of rhombus 8 at levels 16 and 17, and twenty-nine
/// cells of one column of rhombus 0 at levels 18 and 19; the cells about them host what a cell
/// hosts, and so do the same cells at the levels below.
#[test]
fn from_level_16_cells_beside_two_edges_host_nothing() {
    let full = |level: u8| if level % 2 == 1 { 7 } else { 1 };
    for level in 14..=17 {
        let edge = HexA7::lattice_edge(Private, level);
        for back in 1..=5 {
            let hosted = HexA7::cell_zones(Private, level, 8, 0, edge - back).len();
            let expected = if level >= 16 && back >= 2 {
                0
            } else {
                full(level)
            };
            assert_eq!(hosted, expected, "level {level}, column {}", edge - back);
        }
    }
    for level in 16..=19 {
        for row in 0..=45 {
            let hosted = HexA7::cell_zones(Private, level, 0, row, 100).len();
            let expected = if level >= 18 && (12..=40).contains(&row) {
                0
            } else {
                full(level)
            };
            assert_eq!(hosted, expected, "level {level}, row {row}");
        }
    }
}

/// The lattice against the live engine.
#[cfg(feature = "oracle")]
mod engine {
    use super::*;
    use dggal_oracle::{
        IGEO7, ISEA3H, IVEA3H, IVEA7H, RTEA3H, RTEA7H, centroid, count_zones, engine_can_read,
        list_zones, zone_from_text,
    };
    use rs4dggs::ZoneId;
    use std::collections::BTreeSet;
    use std::f64::consts::PI;

    /// The walk of the aperture-3 lattice is the engine's list of the whole world, sequence for
    /// sequence, at levels 0 to 8 on the three grids.
    #[test]
    fn the_aperture_3_lattice_is_the_engines_whole_world_in_its_order() {
        let mut compared = 0;
        for level in 0..=8u8 {
            let ours: Vec<u64> = lattice::<HexA3>(level)
                .iter()
                .map(|(a, _)| I3h::encode(a).0)
                .collect();
            for grid in [ISEA3H, IVEA3H, RTEA3H] {
                let theirs = list_zones(grid, i32::from(level), None);
                assert_eq!(
                    theirs.len() as u64,
                    count_zones(grid, i32::from(level)),
                    "{grid} level {level}"
                );
                assert!(
                    ours == theirs,
                    "{grid} level {level}: the walk of {} zones departs from the engine's list of \
                     {} at entry {:?}",
                    ours.len(),
                    theirs.len(),
                    ours.iter().zip(&theirs).position(|(a, b)| a != b)
                );
                compared += theirs.len();
            }
        }
        assert_eq!(
            compared,
            3 * (0..=8).map(|l| 10 * 3usize.pow(l) + 2).sum::<usize>()
        );
    }

    /// The walk of the aperture-7 lattice is the engine's list of the whole world, sequence
    /// for sequence, at levels 0 to 5 on the three grids. Neither list is sorted: the engine's
    /// order is not that of the Z7 identifiers, which the walk leaves from level 1.
    #[test]
    fn the_aperture_7_lattice_is_the_engines_whole_world_in_its_order() {
        let mut compared = 0;
        for level in 0..=5u8 {
            let ours: Vec<u64> = lattice::<HexA7>(level)
                .iter()
                .map(|(a, _)| Z7::encode(a).0)
                .collect();
            for grid in [IGEO7, IVEA7H, RTEA7H] {
                let theirs = list_zones(grid, i32::from(level), None);
                assert_eq!(
                    theirs.len() as u64,
                    count_zones(grid, i32::from(level)),
                    "{grid} level {level}"
                );
                assert!(
                    ours == theirs,
                    "{grid} level {level}: the walk of {} zones departs from the engine's list of \
                     {} at entry {:?}",
                    ours.len(),
                    theirs.len(),
                    ours.iter().zip(&theirs).position(|(a, b)| a != b)
                );
                compared += theirs.len();
            }
            assert!(
                ours.windows(2).any(|w| w[0] > w[1]),
                "level {level}: the order would be that of the Z7 identifiers"
            );
        }
        assert_eq!(
            compared,
            3 * (0..=5).map(|l| 10 * 7usize.pow(l) + 2).sum::<usize>()
        );
    }

    /// A box of six zone widths about a point, in degrees: south, west, north, east.
    fn box_about(level: u8, lat: f64, lon: f64) -> [f64; 4] {
        let width = (4.0 * PI / (10.0 * 7f64.powi(i32::from(level)) + 2.0)).sqrt();
        let half = 3.0 * width.to_degrees();
        let half_lon = half / lat.to_radians().cos();
        [
            lat - half,
            lon - half_lon,
            (lat + half).min(90.0),
            lon + half_lon,
        ]
    }

    /// Points of the edge of the icosahedron from base cell 1 to base cell 6, in degrees.
    const ON_THE_EDGE_FROM_1_TO_6: [(f64, f64); 5] = [
        (57.361162809, -167.416422128),
        (41.799471943, -153.893098119),
        (30.057436593, -147.871189351),
        (12.079280893, -141.051897277),
        (1.208873949, -137.473618746),
    ];

    /// Small boxes along the two edges where the engine displaces zones. The edge from base
    /// cell 0 to 1 runs over the north pole, along the meridians of 11.2 degrees east and 168.8
    /// degrees west from the latitude of the two cells: eight boxes on each meridian, the
    /// first about the cell itself, and a cap about the pole. On the edge from 1 to 6, five
    /// boxes, to level 17 alone.
    ///
    /// Every box here was first put to the engine in a process of its own, on the three grids:
    /// from level 15 the engine answers a box in which it finds no zone with a null array,
    /// which ends the process that asked. It does so for a box of this size on the edge from 1
    /// to 6 at levels 18 and 19, `[30.057431818, -147.871194868, 30.057441368, -147.871183834]`
    /// at level 18 among them, and such a box is not asked.
    fn boxes_along_the_two_edges(level: u8) -> Vec<[f64; 4]> {
        let mut boxes = Vec::new();
        for lon in [11.2, -168.8] {
            for lat in [
                58.39714590743112,
                58.5,
                60.0,
                67.5,
                77.0,
                86.0,
                89.9,
                89.999,
            ] {
                boxes.push(box_about(level, lat, lon));
            }
        }
        let cap = box_about(level, 90.0, 0.0);
        boxes.push([cap[0], -180.0, 90.0, 180.0]);
        if level <= 17 {
            for (lat, lon) in ON_THE_EDGE_FROM_1_TO_6 {
                boxes.push(box_about(level, lat, lon));
            }
        }
        boxes
    }

    /// The engine's answer for a box given in degrees, each zone once, in two parts: those the
    /// engine reads back from their own text, and those it does not.
    fn listed(grid: &str, level: u8, b: &[f64; 4]) -> (BTreeSet<u64>, BTreeSet<u64>) {
        let radians = b.map(|x| x * (PI / 180.0));
        let zones: BTreeSet<u64> = list_zones(grid, i32::from(level), Some(radians))
            .into_iter()
            .collect();
        assert!(!zones.is_empty(), "{grid} level {level} {b:?}");
        zones.into_iter().partition(|&z| engine_can_read(grid, z))
    }

    /// Every zone that the engine lists there and can read back is hosted once by the cell in
    /// which it is located; an identifier that the engine lists and cannot read back is located
    /// nowhere. Returns how many of each.
    fn hosts_what_the_engine_lists(grid: &str, level: u8, b: &[f64; 4]) -> (usize, usize) {
        let (readable, unreadable) = listed(grid, level, b);
        for &z in &readable {
            let a = Z7::decode(ZoneId(z));
            assert_eq!(
                times_hosted::<HexA7>(level, &a),
                1,
                "{grid} level {level} {b:?}: {}",
                Z7::to_text(ZoneId(z))
            );
        }
        for &z in &unreadable {
            assert_eq!(
                HexA7::locate(Private, &Z7::decode(ZoneId(z))),
                None,
                "{grid} level {level} {b:?}: {z:#x}"
            );
        }
        (readable.len(), unreadable.len())
    }

    /// The lattice hosts the zones that the engine displaces along two edges of the
    /// icosahedron at levels 15 to 19: whatever the engine lists in a small box there and can
    /// itself read back, the lattice hosts once. Among what it lists at levels 16, 18 and 19
    /// are identifiers it cannot read back, which name no zone.
    #[test]
    fn the_lattice_hosts_what_the_engine_lists_along_the_two_edges() {
        // Floors a tenth under what was met on the three grids together.
        for (level, readable_floor, unreadable_floor) in [
            (15u8, 3_000, 0),
            (16, 3_500, 43),
            (17, 6_300, 0),
            (18, 4_300, 45),
            (19, 10_600, 318),
        ] {
            let (mut readable, mut unreadable) = (0, 0);
            for grid in [IGEO7, IVEA7H, RTEA7H] {
                for b in boxes_along_the_two_edges(level) {
                    let (r, u) = hosts_what_the_engine_lists(grid, level, &b);
                    readable += r;
                    unreadable += u;
                }
            }
            eprintln!(
                "level {level}: {readable} zones hosted, {unreadable} identifiers unreadable"
            );
            assert!(readable >= readable_floor, "level {level}: {readable}");
            assert!(
                unreadable >= unreadable_floor,
                "level {level}: {unreadable}"
            );
        }
    }

    /// Two boxes at the north pole at level 17 for which the engine lists zones that lie two
    /// and four cells from the cell that hosts them: the lattice hosts them, and all it lists
    /// there.
    #[test]
    fn the_lattice_hosts_the_displaced_zones_that_the_engine_lists_at_the_pole() {
        for grid in [IGEO7, IVEA7H, RTEA7H] {
            for (b, texts) in [
                (
                    [89.999999940, -180.0, 90.0, 180.0],
                    ["0132323232323232316", "0132323232323232351"].as_slice(),
                ),
                (
                    [89.999985393, -41.492241524, 90.0, -129.455779006],
                    &["0132323232323232620"],
                ),
            ] {
                let (readable, unreadable) = listed(grid, 17, &b);
                assert!(unreadable.is_empty(), "{grid} {b:?}");
                for text in texts {
                    assert!(
                        readable.contains(&zone_from_text(grid, text)),
                        "{grid} {b:?}: {text}"
                    );
                }
                let (hosted, _) = hosts_what_the_engine_lists(grid, 17, &b);
                assert_eq!(hosted, readable.len());
            }
        }
    }

    /// The addresses that the lattice locates nowhere deep under the pentagon of base cell 0
    /// name no zone of the engine's own: it draws two of them where it draws the pentagon, the
    /// other two at one place, and lists none of the four in a box that holds their centroids.
    #[test]
    fn the_engine_names_otherwise_the_zones_that_are_located_nowhere() {
        for grid in [IGEO7, IVEA7H, RTEA7H] {
            let at = |text: &str| centroid(grid, zone_from_text(grid, text));
            let pentagon = at("0000000000000000000");
            assert_eq!(at("0000000000000000120"), pentagon, "{grid}");
            assert_eq!(at("0000000000000000130"), pentagon, "{grid}");
            assert_eq!(
                at("0000000000000000121"),
                at("0000000000000000131"),
                "{grid}"
            );
            let (readable, _) = listed(grid, 17, &box_about(17, 58.39714590743112, 11.2));
            assert!(readable.contains(&zone_from_text(grid, "0000000000000000000")));
            for text in ADDRESSES_OF_ZONES_NAMED_OTHERWISE {
                assert!(
                    !readable.contains(&zone_from_text(grid, text)),
                    "{grid} {text}"
                );
            }
        }
    }
}
