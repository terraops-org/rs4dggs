//! The zones of a level whose extent meets a bounding box, in the order of DGGAL's `listZones`,
//! and the same answer entered after any zone of the level.
//!
//! What needs no engine is tested in every build: the answers for boxes whose zones are known,
//! the pages, the entry after a zone and the refusals. The comparison with the engine, by its
//! own extents and by its own answers, is in the module `engine`, under the `oracle` feature.
use crate::interfaces::sealed::{Private, TopologyPlumbing};
use rs4dggs::indexings::{I3h, Z7};
use rs4dggs::topologies::{HexA3, HexA7};
use rs4dggs::{AnyGrid, Error, Extent, GeoPoint, Indexing, ZoneId};

// The boxes of the comparison with the engine, which the oracle suites share.
#[cfg(feature = "oracle")]
use crate::common::boxes;

const APERTURE_3: [&str; 3] = ["ISEA3H", "IVEA3H", "RTEA3H"];
const APERTURE_7: [&str; 3] = ["IGEO7", "IVEA7H", "RTEA7H"];

/// The box about Lisbon, `[south, west, north, east]` in degrees, as every box here.
const LISBON: [f64; 4] = [38.0, -10.0, 39.0, -9.0];

fn grid(name: &str) -> AnyGrid {
    rs4dggs::get_grid(name).unwrap()
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

/// The zones of `level` whose extent meets the box, whole. The estimate of the box is held
/// against every answer taken so: it is never under the number of zones.
fn in_box(grid: AnyGrid, level: u8, b: &[f64; 4]) -> Vec<ZoneId> {
    let zones: Vec<ZoneId> = grid.zones_in_box(level, &degrees(b)).unwrap().collect();
    let estimate = grid.estimate_zones_in_box(level, &degrees(b)).unwrap();
    assert!(
        estimate >= zones.len() as u64,
        "{} level {level}, the box {b:?}: {} zones, estimated at {estimate}",
        grid.name(),
        zones.len()
    );
    zones
}

/// The first `most` zones of that answer.
fn first_in_box(grid: AnyGrid, level: u8, b: &[f64; 4], most: usize) -> Vec<ZoneId> {
    grid.zones_in_box(level, &degrees(b))
        .unwrap()
        .take(most)
        .collect()
}

fn texts(grid: AnyGrid, zones: &[ZoneId]) -> Vec<String> {
    zones.iter().map(|&z| grid.text_id(z)).collect()
}

/// The answer taken in pages of `size`, each page entered after the last zone of the one
/// before, until a page comes back empty or `most` zones are had; cut to `most`. A cursor that
/// did not advance would yield beyond any answer, and is stopped at `most` by the same test.
fn paged(grid: AnyGrid, level: u8, b: &[f64; 4], size: usize, most: usize) -> Vec<ZoneId> {
    let mut out = first_in_box(grid, level, b, size);
    while let Some(&last) = out.last() {
        if out.len() >= most {
            break;
        }
        let page: Vec<ZoneId> = grid
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

/// The key of a zone in the order of its level, as the lattice gives it.
fn key(grid: AnyGrid, zone: ZoneId) -> u64 {
    let place = if APERTURE_3.contains(&grid.name()) {
        HexA3::locate(Private, &I3h::decode(zone))
    } else {
        HexA7::locate(Private, &Z7::decode(zone))
    };
    place
        .unwrap_or_else(|| panic!("{:#018x} is no zone of the lattice", zone.0))
        .3
}

/// The zones are of `level`, in the order of the level, none twice.
fn are_in_the_order_of_the_level(grid: AnyGrid, level: u8, zones: &[ZoneId], what: &str) {
    for &z in zones {
        assert_eq!(grid.resolution(z), level, "{what}: {}", grid.text_id(z));
    }
    let keys: Vec<u64> = zones.iter().map(|&z| key(grid, z)).collect();
    assert!(
        keys.windows(2).all(|w| w[0] < w[1]),
        "{what}: the keys do not ascend"
    );
}

/// The box about Lisbon on ISEA3H. At level 2 one zone holds the whole box, and it belongs to
/// another rhombus than the one the box lies in; at levels 3 and 4 two zones share the box, one
/// from each rhombus.
#[test]
fn the_box_about_lisbon_has_the_zones_of_both_rhombi() {
    let g = grid("ISEA3H");
    assert_eq!(texts(g, &in_box(g, 2, &LISBON)), ["B4-2-A"]);
    assert_eq!(texts(g, &in_box(g, 3, &LISBON)), ["B2-5-C", "B4-2-B"]);
    assert_eq!(texts(g, &in_box(g, 4, &LISBON)), ["C2-23-A", "C4-6-A"]);
    assert_eq!(texts(g, &in_box(g, 5, &LISBON)), ["C2-23-C"]);
    assert_eq!(
        texts(g, &in_box(g, 6, &LISBON)),
        ["D2-10C-A", "D2-10D-A", "D2-128-A"]
    );
    for (level, count) in [(7, 3), (8, 6), (10, 20)] {
        assert_eq!(in_box(g, level, &LISBON).len(), count, "level {level}");
    }
}

/// A cap about a pole that is smaller than a zone has the zones that cover it, and so has a
/// part of such a cap.
#[test]
fn a_small_cap_about_a_pole_has_the_zones_that_cover_it() {
    let g = grid("ISEA3H");
    let north = [89.5, -180.0, 90.0, 180.0];
    assert_eq!(texts(g, &in_box(g, 6, &north)), ["D0-D-A", "D0-E-A"]);
    let south = [-90.0, -180.0, -89.5, 180.0];
    assert_eq!(texts(g, &in_box(g, 6, &south)), ["D5-15F-A", "D5-17A-A"]);
    let quarter = [89.9, 0.0, 90.0, 90.0];
    assert_eq!(texts(g, &in_box(g, 8, &quarter)), ["E0-29-A"]);
    assert_eq!(in_box(g, 6, &[85.0, -180.0, 90.0, 180.0]).len(), 24);
    // On aperture 7 as well, at both poles.
    let g = grid("IGEO7");
    assert_eq!(in_box(g, 6, &north).len(), 38);
    assert_eq!(in_box(g, 6, &south).len(), 38);
    assert_eq!(in_box(g, 8, &quarter).len(), 21);
}

/// A small box upon a pole at a deep level, where the engine's own search does not return:
/// eight zones, of two rhombi, at once.
#[test]
fn a_polar_box_at_a_deep_level_has_its_eight_zones() {
    let g = grid("ISEA3H");
    let cap = [89.997928885, 11.131869680, 90.0, 108.081826856];
    assert_eq!(
        texts(g, &in_box(g, 20, &cap)),
        [
            "K0-7354-A",
            "K0-7355-A",
            "K0-7356-A",
            "K8-67E8B3CA-A",
            "K8-67E99A72-A",
            "K8-67E99A73-A",
            "K8-67EA811B-A",
            "K8-67EA811C-A"
        ]
    );
    for (level, count) in [(8, 2), (10, 2), (12, 2)] {
        assert_eq!(in_box(g, level, &cap).len(), count, "level {level}");
    }
    // At the finest level the box holds millions of zones, and the first of them cost what the
    // first of any answer cost.
    let first = first_in_box(g, 33, &cap, 500);
    assert_eq!(first.len(), 500);
    are_in_the_order_of_the_level(g, 33, &first, "the polar box at level 33");
}

/// A west above the east is a box that runs eastwards over the antimeridian: from 10 degrees
/// east to 5 degrees east it is 355 degrees wide, and has every zone of its band of latitude
/// but those that lie wholly between the two meridians.
#[test]
fn a_west_above_the_east_is_a_box_over_the_antimeridian() {
    let wide = [30.0, 10.0, 40.0, 5.0];
    let narrow = [30.0, 5.0, 40.0, 10.0];
    let band = [30.0, -180.0, 40.0, 180.0];
    for (name, level, count) in [("ISEA3H", 6, 663), ("IGEO7", 4, 1972)] {
        let g = grid(name);
        let zones = in_box(g, level, &wide);
        assert_eq!(zones.len(), count, "{name}");
        are_in_the_order_of_the_level(g, level, &zones, name);
        // The band is the wide box and the narrow one together, and they share the zones that
        // straddle either meridian.
        let all = in_box(g, level, &band);
        let between = in_box(g, level, &narrow);
        assert!(zones.iter().all(|z| all.contains(z)), "{name}");
        assert!(between.iter().all(|z| all.contains(z)), "{name}");
        assert!(
            all.iter().all(|z| zones.contains(z) || between.contains(z)),
            "{name}"
        );
        assert!(zones.len() < all.len(), "{name}");
    }
    // A box about the antimeridian itself.
    let about = [-5.0, 175.0, 5.0, -175.0];
    assert_eq!(in_box(grid("ISEA3H"), 6, &about).len(), 25);
    assert_eq!(in_box(grid("IGEO7"), 4, &about).len(), 77);
    for name in ["ISEA3H", "IGEO7"] {
        let g = grid(name);
        for z in in_box(g, 4, &about) {
            let lon = g.centroid(z).lon;
            assert!(lon.abs() > 165.0, "{name}: {} at {lon}", g.text_id(z));
        }
    }
}

/// A box may have no area. A point has the zones whose extents strictly hold it, from one to
/// three; a segment those whose extents it strictly crosses; a point exactly on a pole none,
/// since no extent reaches beyond a pole. A point of a line at which the extents of the zones
/// on either side end has none either: of the equator, and of a meridian along which rhombi
/// are joined at an odd level of the aperture-3 grids.
#[test]
fn a_point_and_a_segment_are_boxes() {
    let point = [38.7223, -9.1393, 38.7223, -9.1393];
    let (g3, g7) = (grid("ISEA3H"), grid("IGEO7"));
    assert_eq!(texts(g3, &in_box(g3, 4, &point)), ["C4-6-A"]);
    assert_eq!(texts(g3, &in_box(g3, 6, &point)), ["D2-10D-A", "D2-128-A"]);
    assert_eq!(texts(g7, &in_box(g7, 4, &point)), ["006415"]);
    assert_eq!(texts(g7, &in_box(g7, 6, &point)), ["00641565"]);
    // The zone the point quantises to is among them at every level.
    for (g, levels) in [(g3, 0..=33u8), (g7, 0..=14u8)] {
        for level in levels {
            let zones = in_box(g, level, &point);
            let at = g.zone_from_geo(point[0], point[1], level).unwrap();
            assert!(zones.contains(&at), "{} level {level}", g.name());
            assert!((1..=3).contains(&zones.len()), "{} level {level}", g.name());
        }
    }
    // A segment of a meridian and a segment of a parallel.
    let meridian = [30.0, 10.0, 40.0, 10.0];
    let parallel = [30.0, 10.0, 30.0, 20.0];
    assert_eq!(in_box(g3, 6, &meridian).len(), 6);
    assert_eq!(in_box(g3, 6, &parallel).len(), 4);
    assert_eq!(in_box(g7, 4, &meridian).len(), 10);
    assert_eq!(in_box(g7, 4, &parallel).len(), 7);
    // A point on a pole.
    for name in APERTURE_3.into_iter().chain(APERTURE_7) {
        for pole in [90.0, -90.0] {
            for level in [0, 4, 6] {
                let zones = in_box(grid(name), level, &[pole, 0.0, pole, 0.0]);
                assert!(zones.is_empty(), "{name} level {level}");
            }
        }
    }
}

/// The whole world, as a box, has every zone of the level, in the order of the level.
#[test]
fn the_whole_world_as_a_box_is_every_zone_of_the_level() {
    let world = [-90.0, -180.0, 90.0, 180.0];
    for (names, levels) in [(APERTURE_3, 0..=5u8), (APERTURE_7, 0..=3u8)] {
        for name in names {
            let g = grid(name);
            for level in levels.clone() {
                let every: Vec<ZoneId> = g.zones(level).unwrap().collect();
                assert_eq!(in_box(g, level, &world), every, "{name} level {level}");
            }
        }
    }
}

/// An answer taken in pages of any size, each entered after the last zone of the one before,
/// is the whole answer, each zone once and in order; after the last zone the page is empty, and
/// an ended iterator stays ended.
#[test]
fn pages_entered_after_the_last_zone_of_the_one_before_make_the_whole_answer() {
    let iberia = [36.0, -10.0, 44.0, 3.0];
    let over = [-20.0, 170.0, 10.0, -160.0];
    let cap = [80.0, -180.0, 90.0, 180.0];
    let mut pages = 0;
    for (names, levels) in [(APERTURE_3, [3u8, 6, 8]), (APERTURE_7, [2u8, 3, 5])] {
        for name in names {
            let g = grid(name);
            for level in levels {
                for b in [&iberia, &over, &cap, &LISBON] {
                    let whole = in_box(g, level, b);
                    assert!(!whole.is_empty(), "{name} level {level} {b:?}");
                    are_in_the_order_of_the_level(g, level, &whole, name);
                    for size in [1, 7, 13, 97] {
                        if whole.len() > 600 && size == 1 {
                            continue;
                        }
                        let cut = paged(g, level, b, size, whole.len() + 1);
                        assert_eq!(cut, whole, "{name} level {level} {b:?}, pages of {size}");
                        pages += whole.len() / size + 1;
                    }
                    // After the last zone there is nothing, and there stays nothing.
                    let last = *whole.last().unwrap();
                    let mut after = g
                        .zones_in_box(level, &degrees(b))
                        .unwrap()
                        .after(last)
                        .unwrap();
                    for _ in 0..3 {
                        assert_eq!(after.next(), None, "{name} level {level} {b:?}");
                    }
                }
            }
        }
    }
    assert!(pages >= 10_000, "{pages} pages");
}

/// The answer may be entered after any zone of the level, whether or not the zone is of the
/// answer: what follows is the zones of the answer that come after it in the order of the
/// level.
#[test]
fn after_any_zone_of_the_level_come_the_zones_of_the_answer_that_follow_it() {
    let iberia = [36.0, -10.0, 44.0, 3.0];
    let over = [-20.0, 170.0, 10.0, -160.0];
    let mut entered = 0;
    for (name, level) in [("ISEA3H", 4u8), ("RTEA3H", 3), ("IGEO7", 2), ("IVEA7H", 1)] {
        let g = grid(name);
        let every: Vec<ZoneId> = g.zones(level).unwrap().collect();
        for b in [&iberia, &over, &LISBON] {
            let whole = in_box(g, level, b);
            let keys: Vec<u64> = whole.iter().map(|&z| key(g, z)).collect();
            for &z in &every {
                let k = key(g, z);
                let expected: Vec<ZoneId> = whole
                    .iter()
                    .zip(&keys)
                    .filter(|&(_, &other)| other > k)
                    .map(|(&w, _)| w)
                    .collect();
                let got: Vec<ZoneId> = g
                    .zones_in_box(level, &degrees(b))
                    .unwrap()
                    .after(z)
                    .unwrap()
                    .collect();
                assert_eq!(
                    got,
                    expected,
                    "{name} level {level} {b:?}, after {}",
                    g.text_id(z)
                );
                entered += 1;
            }
        }
    }
    assert_eq!(entered, 3 * (812 + 272 + 492 + 72));
    // An iterator partly consumed, or ended, is entered anew: what it gave before is not
    // counted.
    let g = grid("ISEA3H");
    let whole = in_box(g, 5, &iberia);
    let mut it = g.zones_in_box(5, &degrees(&iberia)).unwrap();
    assert_eq!(it.by_ref().count(), whole.len());
    assert_eq!(it.next(), None);
    let again: Vec<ZoneId> = it.after(whole[2]).unwrap().collect();
    assert_eq!(again, whole[3..]);
}

/// The grid taken by its type and the handle on it give one answer, whole and entered after a
/// zone; an iterator copied in the middle of an answer goes on as the one it was copied from.
#[test]
fn the_grid_and_its_handle_give_one_answer() {
    let iberia = degrees(&[36.0, -10.0, 44.0, 3.0]);
    let by_type: Vec<ZoneId> = rs4dggs::isea3h()
        .zones_in_box(7, &iberia)
        .unwrap()
        .collect();
    let by_handle: Vec<ZoneId> = grid("ISEA3H").zones_in_box(7, &iberia).unwrap().collect();
    assert_eq!(by_type, by_handle);
    assert_eq!(by_type.len(), 60);
    let after_by_type: Vec<ZoneId> = rs4dggs::isea3h()
        .zones_in_box(7, &iberia)
        .unwrap()
        .after(by_type[30])
        .unwrap()
        .collect();
    assert_eq!(after_by_type, by_type[31..]);
    let by_type: Vec<ZoneId> = rs4dggs::rtea7h()
        .zones_in_box(5, &iberia)
        .unwrap()
        .collect();
    let by_handle: Vec<ZoneId> = grid("RTEA7H").zones_in_box(5, &iberia).unwrap().collect();
    assert_eq!(by_type, by_handle);
    let mut it = grid("RTEA7H").zones_in_box(5, &iberia).unwrap();
    let head: Vec<ZoneId> = it.by_ref().take(40).collect();
    let copy = it.clone();
    assert_eq!(head, by_type[..40]);
    assert_eq!(it.collect::<Vec<_>>(), by_type[40..]);
    assert_eq!(copy.collect::<Vec<_>>(), by_type[40..]);
    assert!(
        format!("{:?}", grid("RTEA7H").zones_in_box(5, &iberia).unwrap()).contains("ZonesInBox")
    );
}

/// A level beyond the limit is refused, with the limit; a box that is none is refused, with
/// the value that makes it none; and the level is refused before the box. Nothing is refused by
/// a panic, and nothing is computed first.
#[test]
fn a_level_beyond_the_limit_and_a_box_that_is_none_are_refused() {
    let lisbon = degrees(&LISBON);
    let refused = |name: &str, level: u8, b: &[f64; 4]| {
        let g = grid(name);
        g.zones_in_box(level, &degrees(b)).err()
    };
    for name in APERTURE_7 {
        assert!(grid(name).zones_in_box(14, &lisbon).is_ok(), "{name}");
        for level in [15, 16, 19, 20, 255] {
            assert_eq!(
                refused(name, level, &LISBON),
                Some(Error::ResolutionOutOfRange {
                    res: level,
                    max: 14
                }),
                "{name} level {level}"
            );
        }
    }
    for name in APERTURE_3 {
        assert!(grid(name).zones_in_box(33, &lisbon).is_ok(), "{name}");
        for level in [34, 255] {
            assert_eq!(
                refused(name, level, &LISBON),
                Some(Error::ResolutionOutOfRange {
                    res: level,
                    max: 33
                }),
                "{name} level {level}"
            );
        }
    }
    assert_eq!(
        rs4dggs::igeo7().zones_in_box(15, &lisbon).err(),
        Some(Error::ResolutionOutOfRange { res: 15, max: 14 })
    );
    for name in APERTURE_3.into_iter().chain(APERTURE_7) {
        for (i, field) in ["south", "west", "north", "east"].into_iter().enumerate() {
            for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
                let mut b = LISBON;
                b[i] = bad;
                assert_eq!(
                    refused(name, 3, &b),
                    Some(Error::NonFinite { name: field }),
                    "{name}"
                );
            }
        }
        assert_eq!(
            refused(name, 3, &[80.0, 0.0, 90.000001, 10.0]),
            Some(Error::LatitudeOutOfRange { name: "north" })
        );
        assert_eq!(
            refused(name, 3, &[-91.0, 0.0, 10.0, 10.0]),
            Some(Error::LatitudeOutOfRange { name: "south" })
        );
        assert_eq!(
            refused(name, 3, &[10.0, 170.0, 20.0, -181.0]),
            Some(Error::LongitudeOutOfRange { name: "east" })
        );
        assert_eq!(
            refused(name, 3, &[10.0, 180.5, 20.0, 170.0]),
            Some(Error::LongitudeOutOfRange { name: "west" })
        );
        assert_eq!(
            refused(name, 3, &[40.0, 10.0, 30.0, 20.0]),
            Some(Error::InvertedBox)
        );
        // The level is refused first.
        assert!(matches!(
            refused(name, 40, &[40.0, 10.0, 30.0, f64::NAN]),
            Some(Error::ResolutionOutOfRange { res: 40, .. })
        ));
        // The limits themselves are boxes.
        assert!(refused(name, 3, &[-90.0, -180.0, 90.0, 180.0]).is_none());
        assert!(refused(name, 3, &[90.0, 180.0, 90.0, -180.0]).is_none());
    }
}

/// A zone to enter after that the level does not hold is refused, as the enumeration of the
/// whole level refuses it: a zone of another level, the null zone, an identifier of the other
/// aperture, any pattern of bits that is no zone of the level.
#[test]
fn a_zone_to_enter_after_that_the_level_does_not_hold_is_refused() {
    let lisbon = degrees(&LISBON);
    let mut refused = 0;
    for (names, level, finest) in [(APERTURE_3, 4u8, 33u8), (APERTURE_7, 3u8, 19u8)] {
        for name in names {
            let g = grid(name);
            let after = |z: ZoneId| g.zones_in_box(level, &lisbon).unwrap().after(z);
            let mut patterns = vec![ZoneId::NULL, ZoneId(0), ZoneId(1), ZoneId(u64::MAX - 1)];
            patterns.extend((0..64).map(|bit| ZoneId(1 << bit)));
            // Zones of the levels above and below, and of the finest level of the grid, which
            // on aperture 7 lies beyond the limit of the enumeration.
            for other in [level - 1, level + 1, 14, 15, finest] {
                patterns.push(g.zone_from_geo(38.7, -9.1, other).unwrap());
                patterns.push(g.zone_from_geo(-89.9, 120.0, other).unwrap());
            }
            // Zones of the other aperture.
            let other = grid(if finest == 33 { "IGEO7" } else { "ISEA3H" });
            patterns.extend(other.zones(level).unwrap().take(40));
            for z in patterns {
                let held = g.resolution(z) == level && g.zones(level).unwrap().any(|w| w == z);
                match after(z) {
                    Ok(_) => assert!(held, "{name}: {:#018x} accepted", z.0),
                    Err(e) => {
                        assert!(!held, "{name}: {:#018x} refused", z.0);
                        assert!(matches!(e, Error::InvalidZone(_)), "{name}: {e:?}");
                        assert_eq!(
                            g.zones(level).unwrap().after(z).err(),
                            Some(e),
                            "{name}: {:#018x}",
                            z.0
                        );
                        refused += 1;
                    }
                }
            }
            // A zone of the level is accepted wherever it lies, in the answer or not.
            for z in g.zones(level).unwrap().step_by(37) {
                assert!(after(z).is_ok(), "{name}: {}", g.text_id(z));
            }
        }
    }
    assert!(refused >= 6 * 70, "{refused} refused");
}

#[cfg(feature = "oracle")]
mod engine {
    use super::boxes::{
        KINDS, Kind, Rng, area_in_zones, engine_survives, make_box, meets, radians, sites,
        zone_width_deg,
    };
    use super::*;
    use dggal_oracle::{count_zones, extent, list_zones, neighbors};
    use std::collections::HashSet;
    use std::f64::consts::PI;
    use std::ops::RangeInclusive;

    /// The engine's extent of a zone, which is one of a zone with geometry.
    fn engines_extent(engine: &str, zone: ZoneId) -> [f64; 4] {
        let e = extent(engine, zone.0);
        assert!(e[0].abs() <= PI / 2.0 && e[2].abs() <= PI / 2.0, "{e:?}");
        e
    }

    /// Every zone of a level as the engine lists the whole world, in its order, each with the
    /// engine's own extent.
    fn the_engines_world(engine: &str, level: u8) -> Vec<(ZoneId, [f64; 4])> {
        let listed = list_zones(engine, i32::from(level), None);
        assert_eq!(listed.len() as u64, count_zones(engine, i32::from(level)));
        listed
            .into_iter()
            .map(|z| (ZoneId(z), engines_extent(engine, ZoneId(z))))
            .collect()
    }

    /// The rule, taken over every zone of the engine's world with the engine's own extents:
    /// the zones whose extent meets the box `r`, in radians, in the engine's order.
    fn by_the_rule(world: &[(ZoneId, [f64; 4])], r: &[f64; 4]) -> Vec<ZoneId> {
        world
            .iter()
            .filter(|(_, e)| meets(e, r))
            .map(|(z, _)| *z)
            .collect()
    }

    /// The first zone of the engine's answer that the crate's answer does not hold at its
    /// place: the engine's zones, in its order, must all be found in the crate's, in order.
    fn not_within(theirs: &[u64], ours: &[ZoneId]) -> Option<u64> {
        let mut ours = ours.iter();
        theirs.iter().copied().find(|&t| !ours.any(|o| o.0 == t))
    }

    /// What a run over the boxes of a grid met, by kind of box.
    #[derive(Default)]
    struct Tally {
        /// Boxes made, of each kind.
        made: [usize; 9],
        /// Boxes put to the engine's `listZones`.
        asked: [usize; 9],
        /// Boxes that the engine does not survive, or that would cost it more than a test may
        /// spend, and that it was therefore not asked.
        kept: [usize; 9],
        /// Zones the crate returned.
        zones: usize,
        /// Zones the engine returned.
        engines: usize,
        /// Boxes in which the crate returned more zones than the engine.
        larger: usize,
        /// Boxes in which the engine returned no zone and the crate some.
        empty: usize,
        /// The estimate of a box over the number of its zones, the greatest met, by the size
        /// of the answer: 1 to 9 zones, 10 to 99, 100 to 999, a thousand and more.
        over: [f64; 4],
        /// How many whole answers of each of those sizes were held against their estimate.
        estimated: [usize; 4],
        /// Whole answers, of boxes with a height and a width, that were found to hold the
        /// zone at the centre of the box.
        centred: usize,
    }

    impl Tally {
        /// Records the estimate of a box over the number of zones of its whole answer, which
        /// [`in_box`] has already shown to be at least one.
        fn estimate(&mut self, g: AnyGrid, level: u8, b: &[f64; 4], zones: usize) {
            let class = match zones {
                0 => return,
                1..=9 => 0,
                10..=99 => 1,
                100..=999 => 2,
                _ => 3,
            };
            let estimate = g.estimate_zones_in_box(level, &degrees(b)).unwrap();
            self.over[class] = self.over[class].max(estimate as f64 / zones as f64);
            self.estimated[class] += 1;
        }

        fn report(&self, what: &str) {
            let by_kind = |counts: &[usize; 9]| -> String {
                KINDS
                    .iter()
                    .zip(counts)
                    .map(|(k, c)| format!("{} {c}", k.name()))
                    .collect::<Vec<_>>()
                    .join(", ")
            };
            eprintln!(
                "{what}: {} boxes ({}); asked of the engine {} ({}); kept from it {} ({}); \
                 {} zones, the engine's {}; the answer larger than the engine's in {} boxes, \
                 the engine's empty in {}; the estimate over the count at most {:?} in {:?} \
                 answers of 1 to 9, 10 to 99, 100 to 999 and more zones; the zone at the \
                 centre of the box in {} answers",
                self.made.iter().sum::<usize>(),
                by_kind(&self.made),
                self.asked.iter().sum::<usize>(),
                by_kind(&self.asked),
                self.kept.iter().sum::<usize>(),
                by_kind(&self.kept),
                self.zones,
                self.engines,
                self.larger,
                self.empty,
                self.over,
                self.estimated,
                self.centred
            );
        }
    }

    fn index_of(kind: Kind) -> usize {
        KINDS.iter().position(|&k| k == kind).unwrap()
    }

    /// Puts the box to the engine, where it survives it, and holds the crate's answer against
    /// the engine's: every zone of the engine's is in the crate's, in the engine's order.
    /// `ours` is the crate's whole answer.
    fn the_engines_answer_lies_within(
        engine: &str,
        level: u8,
        kind: Kind,
        b: &[f64; 4],
        ours: &[ZoneId],
        tally: &mut Tally,
    ) {
        let k = index_of(kind);
        if !engine_survives(engine, i32::from(level), b) {
            tally.kept[k] += 1;
            return;
        }
        let theirs = list_zones(engine, i32::from(level), Some(radians(b)));
        tally.asked[k] += 1;
        tally.engines += theirs.len();
        assert_eq!(
            not_within(&theirs, ours),
            None,
            "{engine} level {level}, the box {b:?}: a zone of the engine's {} is not among \
             the crate's {}",
            theirs.len(),
            ours.len()
        );
        if ours.len() > theirs.len() {
            tally.larger += 1;
            if theirs.is_empty() {
                tally.empty += 1;
            }
        }
    }

    /// How many boxes of a kind are made at a level of `count` zones, where every zone of the
    /// level is put to the rule: a hundred, and thirty of the large ones where the level has
    /// more than five thousand zones, an answer there being thousands of zones.
    fn boxes_of(kind: Kind, count: usize) -> usize {
        if kind == Kind::Big && count > 5_000 {
            30
        } else {
            100
        }
    }

    /// How many of them are also put to the engine's `listZones`: the first `asked`, and the
    /// first three of the large ones where the level has more than two thousand zones, each
    /// of which costs the engine from a tenth to half of a second.
    fn asked_of(kind: Kind, count: usize, asked: usize) -> usize {
        if kind == Kind::Big && count > 2_000 {
            3
        } else {
            asked
        }
    }

    /// At the levels at which every zone can be put to the rule, the crate's answer is the
    /// rule's: the zones of the engine's whole world whose extent, the engine's own, meets the
    /// box, the same set in the same order. Every fifth box is also taken in pages of 13. The
    /// first boxes of each kind at each level, `asked` of them or as [`asked_of`] has it, are
    /// put to the engine as well, where it survives them, and its answer lies within the
    /// crate's.
    fn the_answer_is_the_rule_over_the_engines_world(
        name: &str,
        engine: &str,
        levels: RangeInclusive<u8>,
        asked: usize,
        seed: u64,
    ) -> Tally {
        let g = grid(name);
        let mut tally = Tally::default();
        for level in levels {
            let world = the_engines_world(engine, level);
            let mut rng = Rng::new(seed + u64::from(level));
            let sites = sites(engine, i32::from(level), &mut rng);
            let size = zone_width_deg(engine, i32::from(level));
            for kind in KINDS {
                for i in 0..boxes_of(kind, world.len()) {
                    let b = make_box(kind, &mut rng, size, &sites);
                    let what = format!("{name} level {level}, the {} box {b:?}", kind.name());
                    let expected = by_the_rule(&world, &radians(&b));
                    let ours = in_box(g, level, &b);
                    assert!(
                        ours == expected,
                        "{what}: {} zones where the rule has {}; the first to differ {:?}",
                        ours.len(),
                        expected.len(),
                        ours.iter().zip(&expected).find(|(a, b)| a != b)
                    );
                    if i % 5 == 0 {
                        let cut = paged(g, level, &b, 13, expected.len() + 1);
                        assert!(
                            cut == expected,
                            "{what}: the pages of 13 are not the answer"
                        );
                    }
                    tally.made[index_of(kind)] += 1;
                    tally.zones += ours.len();
                    tally.estimate(g, level, &b, ours.len());
                    if i < asked_of(kind, world.len(), asked) {
                        the_engines_answer_lies_within(engine, level, kind, &b, &ours, &mut tally);
                    }
                }
            }
        }
        tally.report(name);
        tally
    }

    /// The most zones taken of an answer too large to take whole.
    const MOST: usize = 300;

    /// The zone at the centre of a box, in degrees, whose height and whose width both exceed
    /// 1e-9 radian: a zone that every answer to the box must hold, its extent holding the
    /// centre, about which the box reaches half of that on every side. `None` for a box
    /// without such a height and width, a point or a segment among them: a point on a
    /// meridian along which rhombi are joined may yield no zone, the extents of the zones on
    /// either side stopping up to 2.6e-10 radian short of it. `None` too where the centre has
    /// no zone.
    fn the_zone_at_the_centre(g: AnyGrid, level: u8, b: &[f64; 4]) -> Option<ZoneId> {
        let [s, w, n, e] = radians(b);
        let width = if e < w { e - w + 2.0 * PI } else { e - w };
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
        let at = g.zone_from_geo((b[0] + b[2]) / 2.0, lon, level).unwrap();
        (at != ZoneId::NULL).then_some(at)
    }

    /// Beyond the levels at which every zone can be put to the rule, on boxes of every kind:
    /// every zone the crate returns passes the rule on the engine's own extent; they are in the
    /// order of the level, none twice; pages of 13 give the same sequence, and nothing after
    /// its last zone; no neighbour of a returned zone that passes the rule is missing; and
    /// where the engine survives the box, its answer lies within the crate's. An answer of a
    /// box of more than five thousand zones' area is taken to its first [`MOST`] zones, and
    /// is held to all but the last two.
    ///
    /// On the aperture-3 grids the engine's own answer omits zones, and may be empty, so that
    /// an answer that were wrongly empty would pass every test above. A whole answer of a box
    /// with a height and a width is therefore held to the zone at the centre of the box
    /// ([`the_zone_at_the_centre`]).
    fn the_answer_holds_at_depth(
        name: &str,
        engine: &str,
        levels: RangeInclusive<u8>,
        boxes_a_kind: usize,
        seed: u64,
    ) -> Tally {
        let g = grid(name);
        let mut tally = Tally::default();
        for level in levels {
            let mut rng = Rng::new(seed + u64::from(level));
            let sites = sites(engine, i32::from(level), &mut rng);
            let size = zone_width_deg(engine, i32::from(level));
            for kind in KINDS {
                for _ in 0..boxes_a_kind {
                    let b = make_box(kind, &mut rng, size, &sites);
                    let r = radians(&b);
                    let what = format!("{name} level {level}, the {} box {b:?}", kind.name());
                    let whole = area_in_zones(engine, i32::from(level), &b) <= 5_000.0;
                    let ours = if whole {
                        in_box(g, level, &b)
                    } else {
                        first_in_box(g, level, &b, MOST)
                    };
                    tally.made[index_of(kind)] += 1;
                    tally.zones += ours.len();
                    are_in_the_order_of_the_level(g, level, &ours, &what);
                    for &z in &ours {
                        assert!(
                            meets(&engines_extent(engine, z), &r),
                            "{what}: {} is returned and its extent does not meet the box",
                            g.text_id(z)
                        );
                    }
                    // A whole answer is paged to one zone beyond its length, so that a zone
                    // yielded after the last is seen.
                    let most = if whole { ours.len() + 1 } else { MOST };
                    let cut = paged(g, level, &b, 13, most);
                    assert!(cut == ours, "{what}: the pages of 13 are not the answer");
                    if !whole {
                        assert_eq!(ours.len(), MOST, "{what}");
                        tally.kept[index_of(kind)] += 1;
                        continue;
                    }
                    tally.estimate(g, level, &b, ours.len());
                    if let Some(at) = the_zone_at_the_centre(g, level, &b) {
                        assert!(
                            ours.contains(&at),
                            "{what}: {}, the zone at the centre of the box, is missing",
                            g.text_id(at)
                        );
                        tally.centred += 1;
                    }
                    let held: HashSet<u64> = ours.iter().map(|z| z.0).collect();
                    for &z in &ours {
                        for n in neighbors(engine, z.0) {
                            assert!(
                                held.contains(&n) || !meets(&engines_extent(engine, ZoneId(n)), &r),
                                "{what}: {}, a neighbour of {}, meets the box and is missing",
                                g.text_id(ZoneId(n)),
                                g.text_id(z)
                            );
                        }
                    }
                    the_engines_answer_lies_within(engine, level, kind, &b, &ours, &mut tally);
                }
            }
        }
        tally.report(name);
        tally
    }

    /// The floors of a run: so many boxes made, so many of them put to the engine, and every
    /// kind of box among those put to it.
    fn has_at_least(tally: &Tally, made: usize, asked: usize, zones: usize, kinds_asked: &[Kind]) {
        assert!(tally.made.iter().sum::<usize>() >= made);
        assert!(tally.asked.iter().sum::<usize>() >= asked);
        assert!(tally.zones >= zones, "{} zones", tally.zones);
        for &kind in kinds_asked {
            assert!(
                tally.asked[index_of(kind)] > 0,
                "no {} box asked",
                kind.name()
            );
        }
        // Every box is put to the engine, or counted as kept from it, or was beyond the number
        // of its kind that are put to it at the level.
        for k in 0..9 {
            assert!(tally.asked[k] + tally.kept[k] <= tally.made[k]);
        }
        // Three boxes in four at the least had their whole answer held against its estimate,
        // which was under none; and it was no farther above the answer than the greatest
        // measured on boxes of every kind, with a little room: some twenty times an answer of
        // under ten zones, three tenths above one of a thousand zones or more.
        assert!(tally.estimated.iter().sum::<usize>() * 4 >= tally.made.iter().sum::<usize>() * 3);
        for (class, ceiling) in [22.0, 6.0, 2.3, 1.3].into_iter().enumerate() {
            assert!(
                tally.over[class] <= ceiling,
                "the estimate over the count, by size: {:?}",
                tally.over
            );
        }
    }

    /// On the aperture-3 grids the engine's own search omits zones, so that in some boxes the
    /// answer is larger than the engine's, and in some of those the engine's is empty: the
    /// departure is met as the tests run, in hundreds of boxes and in dozens.
    fn larger_than_the_engines(tally: &Tally) {
        assert!(tally.larger >= 80, "larger in {} boxes", tally.larger);
        assert!(
            tally.empty >= 10,
            "the engine's empty in {} boxes",
            tally.empty
        );
    }

    #[test]
    fn the_answer_is_the_rule_on_isea3h() {
        let t = the_answer_is_the_rule_over_the_engines_world(
            "ISEA3H",
            dggal_oracle::ISEA3H,
            0..=6,
            25,
            0x3a_0000,
        );
        has_at_least(&t, 6_230, 1_500, 220_000, &KINDS);
        larger_than_the_engines(&t);
    }

    #[test]
    fn the_answer_is_the_rule_on_ivea3h() {
        let t = the_answer_is_the_rule_over_the_engines_world(
            "IVEA3H",
            dggal_oracle::IVEA3H,
            0..=6,
            25,
            0x3b_0000,
        );
        has_at_least(&t, 6_230, 1_500, 220_000, &KINDS);
        larger_than_the_engines(&t);
    }

    #[test]
    fn the_answer_is_the_rule_on_rtea3h() {
        let t = the_answer_is_the_rule_over_the_engines_world(
            "RTEA3H",
            dggal_oracle::RTEA3H,
            0..=6,
            25,
            0x3c_0000,
        );
        has_at_least(&t, 6_230, 1_500, 220_000, &KINDS);
        larger_than_the_engines(&t);
    }

    #[test]
    fn the_answer_is_the_rule_on_igeo7() {
        let t = the_answer_is_the_rule_over_the_engines_world(
            "IGEO7",
            dggal_oracle::IGEO7,
            0..=4,
            10,
            0x7a_0000,
        );
        has_at_least(&t, 4_430, 430, 330_000, &KINDS);
    }

    #[test]
    fn the_answer_is_the_rule_on_ivea7h() {
        let t = the_answer_is_the_rule_over_the_engines_world(
            "IVEA7H",
            dggal_oracle::IVEA7H,
            0..=4,
            10,
            0x7b_0000,
        );
        has_at_least(&t, 4_430, 430, 330_000, &KINDS);
    }

    #[test]
    fn the_answer_is_the_rule_on_rtea7h() {
        let t = the_answer_is_the_rule_over_the_engines_world(
            "RTEA7H",
            dggal_oracle::RTEA7H,
            0..=4,
            10,
            0x7c_0000,
        );
        has_at_least(&t, 4_430, 430, 330_000, &KINDS);
    }

    /// The kinds of box that the engine is asked at depth: the large ones hold more zones than
    /// a test may ask it for.
    const ASKED_AT_DEPTH: [Kind; 8] = [
        Kind::Rand,
        Kind::Tiny,
        Kind::Anti,
        Kind::Pole,
        Kind::Vertex,
        Kind::Point,
        Kind::Line,
        Kind::Corner,
    ];

    #[test]
    fn the_answer_holds_at_depth_on_isea3h() {
        let t = the_answer_holds_at_depth("ISEA3H", dggal_oracle::ISEA3H, 7..=33, 5, 0x3d_0000);
        has_at_least(&t, 1_215, 890, 64_000, &ASKED_AT_DEPTH);
        larger_than_the_engines(&t);
        assert!(
            t.centred >= 730,
            "{} answers held to their centre",
            t.centred
        );
    }

    #[test]
    fn the_answer_holds_at_depth_on_ivea3h() {
        let t = the_answer_holds_at_depth("IVEA3H", dggal_oracle::IVEA3H, 7..=33, 5, 0x3e_0000);
        has_at_least(&t, 1_215, 890, 64_000, &ASKED_AT_DEPTH);
        larger_than_the_engines(&t);
        assert!(
            t.centred >= 730,
            "{} answers held to their centre",
            t.centred
        );
    }

    #[test]
    fn the_answer_holds_at_depth_on_rtea3h() {
        let t = the_answer_holds_at_depth("RTEA3H", dggal_oracle::RTEA3H, 7..=33, 5, 0x3f_0000);
        has_at_least(&t, 1_215, 890, 64_000, &ASKED_AT_DEPTH);
        larger_than_the_engines(&t);
        assert!(
            t.centred >= 730,
            "{} answers held to their centre",
            t.centred
        );
    }

    #[test]
    fn the_answer_holds_at_depth_on_igeo7() {
        let t = the_answer_holds_at_depth("IGEO7", dggal_oracle::IGEO7, 5..=14, 5, 0x7d_0000);
        has_at_least(&t, 450, 395, 20_500, &ASKED_AT_DEPTH);
        assert!(
            t.centred >= 280,
            "{} answers held to their centre",
            t.centred
        );
    }

    #[test]
    fn the_answer_holds_at_depth_on_ivea7h() {
        let t = the_answer_holds_at_depth("IVEA7H", dggal_oracle::IVEA7H, 5..=14, 5, 0x7e_0000);
        has_at_least(&t, 450, 395, 20_500, &ASKED_AT_DEPTH);
        assert!(
            t.centred >= 280,
            "{} answers held to their centre",
            t.centred
        );
    }

    #[test]
    fn the_answer_holds_at_depth_on_rtea7h() {
        let t = the_answer_holds_at_depth("RTEA7H", dggal_oracle::RTEA7H, 5..=14, 5, 0x7f_0000);
        has_at_least(&t, 450, 395, 20_500, &ASKED_AT_DEPTH);
        assert!(
            t.centred >= 280,
            "{} answers held to their centre",
            t.centred
        );
    }

    /// Where the answer departs from the engine's, shown against the engine as it runs: on the
    /// aperture-3 grids the engine returns no zone for the box about Lisbon at level 2 and one
    /// at levels 3 and 4, and none for a small cap about the north pole, where the rule has
    /// the zones that hold them; from level 6 the two agree on the box about Lisbon. On the
    /// aperture-7 grids the two agree on all of these.
    #[test]
    fn the_answer_has_the_zones_that_the_engines_search_does_not_reach() {
        let g = grid("ISEA3H");
        let engines = |level: u8, b: &[f64; 4]| -> Vec<String> {
            assert!(
                engine_survives(dggal_oracle::ISEA3H, i32::from(level), b),
                "ISEA3H level {level} {b:?}"
            );
            let listed = list_zones(dggal_oracle::ISEA3H, i32::from(level), Some(radians(b)));
            texts(g, &listed.into_iter().map(ZoneId).collect::<Vec<_>>())
        };
        assert!(engines(2, &LISBON).is_empty());
        assert_eq!(engines(3, &LISBON), ["B2-5-C"]);
        assert_eq!(engines(4, &LISBON), ["C2-23-A"]);
        assert_eq!(texts(g, &in_box(g, 2, &LISBON)), ["B4-2-A"]);
        for level in [6, 8, 10] {
            assert_eq!(
                engines(level, &LISBON),
                texts(g, &in_box(g, level, &LISBON))
            );
        }
        let north = [89.5, -180.0, 90.0, 180.0];
        assert!(engines(6, &north).is_empty());
        assert_eq!(in_box(g, 6, &north).len(), 2);
        // The box of 355 degrees: the same zones, at level 5, the finest at which a box so
        // long is put to the engine.
        let wide = [30.0, 10.0, 40.0, 5.0];
        assert_eq!(engines(5, &wide), texts(g, &in_box(g, 5, &wide)));
        assert_eq!(in_box(g, 5, &wide).len(), 258);
        let g = grid("IGEO7");
        for (level, b) in [
            (2u8, &LISBON),
            (3, &LISBON),
            (4, &LISBON),
            (4, &wide),
            (6, &north),
        ] {
            assert!(
                engine_survives(dggal_oracle::IGEO7, i32::from(level), b),
                "IGEO7 level {level} {b:?}"
            );
            let listed = list_zones(dggal_oracle::IGEO7, i32::from(level), Some(radians(b)));
            let listed: Vec<ZoneId> = listed.into_iter().map(ZoneId).collect();
            assert_eq!(listed, in_box(g, level, b), "IGEO7 level {level} {b:?}");
        }
    }
}
