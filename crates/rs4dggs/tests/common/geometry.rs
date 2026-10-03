//! The refined boundary and the extent of a zone against the live engine: `getZoneRefinedWGS84Vertices`
//! at the edge refinements 0 (the engine's own choice), 1, 2 and 3, and `getZoneWGS84Extent`, at every
//! zone of the suites' sample and at a fixed list the suites always visit (the zones that
//! pin each compiled-order site, and on aperture 3 the sub-hexes C and D
//! of the polar roots); and, at one zone of each kind, the refinements 2,000, 20,000 and 100,000.
//!
//! Every comparison is exact: the engine answers in radians, and this crate in degrees by the
//! product with `180 / Pi`, so that each double compared is the engine's own, converted as the
//! crate converts it (`RAD2DEG` below is the crate's own expression); a ring compares as a
//! sequence, its length and every coordinate of every point; an extent as four doubles. Every
//! count has a floor, or is asserted exactly where it is a finding about the engine.
//!
//! The ring of the engine is NOT put through `in_crate_order`: that turn is for the plain
//! vertices. The engine's refined ring runs anticlockwise on all six grids, and this
//! crate gives it as it is.
//!
//! An identifier this crate has no geometry for answers the empty ring, no extent and no area. Each is classified against the live engine by
//! `expected_divergences::engine_answers_without_geometry`, which says how the engine answers it
//! and re-checks that at every instance; the kinds are counted exactly. On aperture 3 the engine
//! gives a null array as the ring of an identifier it cannot read, and the vendored binding
//! crashes on it: its ring is not asked (the one exception to the rule that no oracle call is
//! skipped), and the number of them is counted. Their extent and area are asked.
use dggal_oracle as o;
use rs4dggs::{Error, Topology, ZoneId, get_grid};

use super::expected_divergences::{self as div, NoGeometry};
use super::subject::Subject;
use super::suite::{SEED, zones};
use super::{Floor, Winding, cross, dot, normalised, unit_vector, winding};

/// The crate's own conversion of a radian to a degree (`math::RAD2DEG`).
const RAD2DEG: f64 = 180.0 / std::f64::consts::PI;

/// The pinning zones on aperture 7, the pinning zones of each compiled-order site and one
/// zone of each kind, among them the zone centred on a pole at level 19, the identifier the
/// engine writes for a neighbour and cannot read back (`01222222222222222220`), the two
/// identifiers of twenty digits and the null zone.
const PINNED_OF_APERTURE_7: [u64; 50] = [
    0x203f_ffff_ffff_ffff,
    0x713f_ffff_ffff_ffff,
    0x787f_ffff_ffff_ffff,
    0xafff_ffff_ffff_ffff,
    0x269f_ffff_ffff_ffff,
    0x6107_ffff_ffff_ffff,
    0x095f_ffff_ffff_ffff,
    0x1687_ffff_ffff_ffff,
    0x6fff_ffff_ffff_ffff,
    0x11ff_ffff_ffff_ffff,
    0x19ff_ffff_ffff_ffff,
    0x1dff_ffff_ffff_ffff,
    0x51ff_ffff_ffff_ffff,
    0x00ff_ffff_ffff_ffff,
    0x023f_ffff_ffff_ffff,
    0x02ff_ffff_ffff_ffff,
    0x0fff_ffff_ffff_ffff,
    0x0000_ffff_ffff_ffff,
    0x0000_0000_3fff_ffff,
    0x0000_03ff_ffff_ffff,
    0x0000_0000_00ff_ffff,
    0x0000_0000_0000_003f,
    0xb03f_ffff_ffff_ffff,
    0xb000_0000_3fff_ffff,
    0xb000_0000_00ff_ffff,
    0xb000_ffff_ffff_ffff,
    0x0007_ffff_ffff_ffff,
    0x0000_0001_ffff_ffff,
    0x0000_0000_0000_7fff,
    0x01ff_ffff_ffff_ffff,
    0x0000_0000_07ff_ffff,
    0x0000_0000_001f_ffff,
    0x0000_0000_0000_0007,
    0xb007_ffff_ffff_ffff,
    0xb000_1fff_ffff_ffff,
    0xb000_0001_ffff_ffff,
    0x32ea_5fff_ffff_ffff,
    0x1000_1fff_ffff_ffff,
    0x0000_1fff_ffff_ffff,
    0x1003_7fff_ffff_ffff,
    0x169b_dfff_ffff_ffff,
    0x169a_7fff_ffff_ffff,
    0x169a_69a6_9a69_a69f,
    0x0a69_a69a_69a6_9a6f,
    0x0000_0000_0000_b3ff,
    0x169a_69a6_9a69_7fff,
    0x1492_4924_9249_243f,
    0x0000_0000_0000_0000,
    0x329c_b814_e5c0_a72e,
    0xffff_ffff_ffff_ffff,
];

/// The pinning zones on aperture 3: the pinning zones of each site, the eight arms of the
/// `crs84` base vertices, one zone of each kind, and the thirty identifiers of a polar root with
/// a non-zero index whose ring the engine answers with a null array.
const PINNED_OF_APERTURE_3: [u64; 118] = [
    0x0360_0000_0000_0000,
    0x0800_0000_0000_0010,
    0x20ad_2a9f_c95e_2b84,
    0x20ad_2a9f_d3a1_8888,
    0x20ad_2a9f_e828_4290,
    0x0a20_0000_0000_5b20,
    0x0260_0000_0000_0018,
    0x02c0_0000_0000_0004,
    0x0460_0000_0000_00fc,
    0x0800_0000_0000_008c,
    0x0400_0000_0000_000c,
    0x0460_0000_0000_00b4,
    0x0020_0000_0000_0000,
    0x0000_0000_0000_0001,
    0x0000_0000_0000_0003,
    0x0020_0000_0000_0001,
    0x0280_0000_0000_0000,
    0x0280_0000_0000_0010,
    0x0280_0000_0000_0014,
    0x0e00_0000_0000_32d4,
    0x0e00_0000_0000_32d6,
    0x0c00_0000_0000_129d,
    0x0c00_0000_0000_1119,
    0x1000_0000_0000_9958,
    0x0c00_0000_0000_12a1,
    0x0c20_0000_0006_fc59,
    0x0c20_0000_001d_8001,
    0x0c00_0000_0000_1118,
    0x0c00_0000_001f_2578,
    0x0e00_0000_0000_3344,
    0x1000_0000_0000_0000,
    0x1400_0000_0000_0000,
    0x0c20_0000_0007_1e82,
    0x0c00_0000_0000_0000,
    0x0c00_0000_0000_058f,
    0x0e00_0000_0000_0000,
    0x0c40_0000_0000_0959,
    0x0c40_0000_0000_0a81,
    0x0c40_0000_0000_0acd,
    0x0c40_0000_0000_0961,
    0x0c40_0000_0000_09a9,
    0x0c40_0000_0000_0ad1,
    0x0c40_0000_0000_0891,
    0x0c40_0000_0000_09b5,
    0x0c40_0000_0000_09fd,
    0x0c60_0000_0018_f625,
    0x0c60_0000_001d_e681,
    0x0c60_0000_001f_9759,
    0x0c60_0000_0019_0ced,
    0x0c60_0000_001e_ca51,
    0x0c60_0000_0018_6211,
    0x0c60_0000_001a_d48d,
    0x0c60_0000_001b_a195,
    0x1240_0000_0001_2b7c,
    0x1e40_0000_0354_d5d0,
    0x1240_0000_0001_03fc,
    0x1240_0000_0000_ed8c,
    0x1240_0000_0000_f54c,
    0x1240_0000_0001_0c98,
    0x1240_0000_0000_e74c,
    0x1240_0000_0000_f6b8,
    0x1240_0000_0000_fe68,
    0x0c60_0000_0019_ce90,
    0x0c60_0000_001a_9b98,
    0x0c60_0000_001b_68a0,
    0x0a60_0000_0001_cb6c,
    0x0c60_0000_001a_b260,
    0x0c60_0000_001b_7f68,
    0x0c60_0000_0019_3a7c,
    0x0c60_0000_001a_0784,
    0x0540_0000_0000_0000,
    0x0560_0000_0000_0000,
    0x0500_0000_0000_00b0,
    0x0400_0000_0000_0000,
    0x0540_0000_0000_0001,
    0x0560_0000_0000_0001,
    0x0500_0000_0000_00b1,
    0x0400_0000_0000_0001,
    0x0400_0000_0000_000e,
    0x04a0_0000_0000_006f,
    0x0480_0000_0000_0125,
    0x0400_0000_0000_0018,
    0x0800_0000_0000_01e8,
    0x0400_0000_0000_0012,
    0x0a00_0000_0000_01e4,
    0x1800_0000_0010_37e0,
    0x2000_0000_0521_ae82,
    0x210d_2a9f_c95e_2b82,
    0x0340_0000_0000_0008,
    0x0360_0000_0000_000c,
    0x0360_0000_0000_0018,
    0x0340_0000_0000_0009,
    0x0360_0000_0000_000d,
    0x0360_0000_0000_0019,
    0x0540_0000_0000_0014,
    0x0560_0000_0000_0048,
    0x0560_0000_0000_0120,
    0x0560_0000_0000_0024,
    0x0540_0000_0000_0015,
    0x0560_0000_0000_0049,
    0x0560_0000_0000_0121,
    0x0560_0000_0000_0025,
    0x0740_0000_0000_0038,
    0x0760_0000_0000_006c,
    0x0760_0000_0000_0af8,
    0x0760_0000_0000_00d8,
    0x0740_0000_0000_0039,
    0x0760_0000_0000_006d,
    0x0760_0000_0000_0af9,
    0x0760_0000_0000_00d9,
    0x0940_0000_0000_00a4,
    0x0940_0000_0000_00a5,
    0x0b40_0000_0000_01e8,
    0x0b40_0000_0000_01e9,
    0x0d40_0000_0000_05b4,
    0x0d40_0000_0000_05b5,
    0x0f40_0000_0000_1118,
    0x0f40_0000_0000_1119,
];

/// The ring at a refinement between 400 and 100,000: one zone of each kind of ring the engine
/// gives. The identifiers are the same on the three grids of an aperture.
const KINDS_OF_APERTURE_7: [(&str, u64); 11] = [
    ("a hexagon", 0x32ea_5fff_ffff_ffff),
    ("a pentagon", 0x1000_1fff_ffff_ffff),
    (
        "a pentagon at the north pole of the icosahedron",
        0x0000_1fff_ffff_ffff,
    ),
    (
        "a pentagon at the south pole of the icosahedron",
        0xb000_ffff_ffff_ffff,
    ),
    ("a hexagon across an interruption", 0x1003_7fff_ffff_ffff),
    (
        "a hexagon with an edge beside a pole",
        0x169b_dfff_ffff_ffff,
    ),
    ("a hexagon that touches a pole", 0x169a_7fff_ffff_ffff),
    (
        "the zone centred on a pole, a ring of one point",
        0x169a_69a6_9a69_a69f,
    ),
    ("a zone beside it", 0x0a69_a69a_69a6_9a6f),
    ("a zone that loses a vertex", 0x0000_0000_0000_b3ff),
    ("a ring that crosses itself", 0x169a_69a6_9a69_7fff),
];
const KINDS_OF_APERTURE_3: [(&str, u64); 11] = [
    ("a hexagon", 0x0480_0000_0000_0125),
    ("a pentagon", 0x0400_0000_0000_0001),
    (
        "a pentagon at the north pole of the icosahedron",
        0x0540_0000_0000_0001,
    ),
    (
        "a pentagon at the south pole of the icosahedron",
        0x0560_0000_0000_0000,
    ),
    ("a hexagon across an interruption", 0x0400_0000_0000_0018),
    (
        "a hexagon with an edge beside a pole",
        0x0800_0000_0000_01e8,
    ),
    ("a hexagon that touches a pole", 0x0400_0000_0000_0012),
    ("the pole named twice", 0x0a00_0000_0000_01e4),
    (
        "a ring that names the pole at five longitudes",
        0x1800_0000_0010_37e0,
    ),
    (
        "the zone centred on a pole at level 33",
        0x2000_0000_0521_ae82,
    ),
    ("a zone beside it", 0x210d_2a9f_c95e_2b82),
];

/// The identifiers of the fixed list, of the grid under test, which the suites always visit.
pub fn fixed_zones<S: Subject>() -> Vec<ZoneId> {
    let mut ids: Vec<u64> = match S::T::APERTURE {
        7 => PINNED_OF_APERTURE_7.to_vec(),
        3 => {
            let mut v = PINNED_OF_APERTURE_3.to_vec();
            v.extend(div::polar_sub_hexes::<S>());
            v
        }
        a => panic!("no fixed zones are recorded for aperture {a}"),
    };
    if S::T::APERTURE == 7 {
        // The third zone centred on a pole at level 19, of the three known.
        ids.push(
            S::grid()
                .zone_from_text("005151515151515151551")
                .unwrap()
                .id()
                .0,
        );
    }
    ids.extend(div::FABRICATED);
    ids.sort_unstable();
    ids.dedup();
    ids.into_iter().map(ZoneId).collect()
}

/// The zones the boundary, the extent and the area are compared at: the suites' sample and the
/// fixed list.
pub fn geometry_zones<S: Subject>() -> Vec<ZoneId> {
    let mut ids = zones::<S>();
    ids.extend(fixed_zones::<S>());
    ids.sort();
    ids.dedup();
    ids
}

/// The ring in the bits of its doubles, to compare what `PartialEq` of `f64` would let pass
/// (the sign of a zero).
fn bits(ring: &[rs4dggs::GeoPoint]) -> Vec<(u64, u64)> {
    ring.iter()
        .map(|p| (p.lat.to_bits(), p.lon.to_bits()))
        .collect()
}

/// A point of the engine in radians, in degrees as this crate gives it.
fn degrees(p: (f64, f64)) -> (f64, f64) {
    (p.0 * RAD2DEG, p.1 * RAD2DEG)
}

/// How many pairs of edges of a ring, not adjacent to one another, cross properly on the
/// sphere, seen in the orthographic plane about the ring's mean: the edges with an end of one
/// strictly on either side of the other and the reverse. A vertex on a pole is the pole itself,
/// whatever its longitude, as in [`winding`]: the edges that meet at a pole touch there and
/// never cross. A point within the rounding of the line through an edge, some 4e-15 times the
/// ring's own radius in the unit of an area, is taken to lie on it, so that the points a refined
/// edge puts along a straight line never cross each other by rounding. (Read otherwise, with
/// `cos(Pi / 2) = 6.1e-17` the several longitudes of a pole would be distinct points.)
fn proper_crossings(ring: &[(f64, f64)]) -> usize {
    let n = ring.len();
    if n < 4 {
        return 0;
    }
    let at = |p: &(f64, f64)| {
        if p.0.abs() == 90.0 {
            [0.0, 0.0, p.0.signum()]
        } else {
            unit_vector(p.0, p.1)
        }
    };
    let v: Vec<[f64; 3]> = ring.iter().map(at).collect();
    let m = normalised([0, 1, 2].map(|k| v.iter().map(|p| p[k]).sum::<f64>()));
    let axis = if m[0].abs() <= m[1].abs() && m[0].abs() <= m[2].abs() {
        [1.0, 0.0, 0.0]
    } else if m[1].abs() <= m[2].abs() {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    let e1 = normalised(cross(m, axis));
    let e2 = cross(m, e1);
    let p: Vec<(f64, f64)> = v.iter().map(|q| (dot(*q, e1), dot(*q, e2))).collect();
    let (cx, cy) = (
        p.iter().map(|q| q.0).sum::<f64>() / n as f64,
        p.iter().map(|q| q.1).sum::<f64>() / n as f64,
    );
    let r2 = p
        .iter()
        .map(|q| (q.0 - cx).powi(2) + (q.1 - cy).powi(2))
        .fold(0.0, f64::max);
    let eps = 4e-15 * r2.sqrt();
    // The side of c of the line from a to b: 1, -1, or 0 within the tolerance.
    let side = |a: (f64, f64), b: (f64, f64), c: (f64, f64)| {
        let s = (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0);
        if s > eps {
            1
        } else if s < -eps {
            -1
        } else {
            0
        }
    };
    let mut crossings = 0;
    for i in 0..n {
        let (a, b) = (p[i], p[(i + 1) % n]);
        for j in (i + 2)..n {
            if i == 0 && j == n - 1 {
                continue;
            }
            let (c, d) = (p[j], p[(j + 1) % n]);
            if side(a, b, c) * side(a, b, d) == -1 && side(c, d, a) * side(c, d, b) == -1 {
                crossings += 1;
            }
        }
    }
    crossings
}

/// How many pairs of edges of a ring, not adjacent to one another, cross properly in the plane
/// of the longitudes and latitudes as they are, with strict signs.
fn crossings_in_degrees(ring: &[rs4dggs::GeoPoint]) -> usize {
    let n = ring.len();
    let side = |a: &rs4dggs::GeoPoint, b: &rs4dggs::GeoPoint, c: &rs4dggs::GeoPoint| {
        let s = (b.lon - a.lon) * (c.lat - a.lat) - (b.lat - a.lat) * (c.lon - a.lon);
        (s > 0.0) as i32 - (s < 0.0) as i32
    };
    let mut crossings = 0;
    for i in 0..n {
        let (a, b) = (&ring[i], &ring[(i + 1) % n]);
        for j in (i + 2)..n {
            if i == 0 && j == n - 1 {
                continue;
            }
            let (c, d) = (&ring[j], &ring[(j + 1) % n]);
            if side(a, b, c) * side(a, b, d) == -1 && side(c, d, a) * side(c, d, b) == -1 {
                crossings += 1;
            }
        }
    }
    crossings
}

/// What the comparison of the boundary and the extent met. The count of coordinates is a
/// floor; every other is a finding about the engine, asserted exactly per grid.
#[derive(Debug, Default, PartialEq, Eq)]
struct Kinds {
    /// Zones with a geometry, compared at the automatic refinement.
    rings: usize,
    /// The coordinates compared (a point is two), over every refinement and the extents.
    coordinates: usize,
    /// Zones without geometry, by the way the engine answers them.
    null_zone: usize,
    undrawn: usize,
    undrawn_without_area: usize,
    drawn: usize,
    /// Of those, identifiers on aperture 3 whose ring was not asked of the engine.
    ring_not_asked: usize,
    /// Rings of one point.
    one_point: usize,
    /// Rings of more than one point that do not run anticlockwise (the crate mirrors them):
    /// four zones of level 33 on aperture 3, centred on or beside a pole, whose every point but
    /// two is on the pole, and which so enclose no area by `winding`'s measure.
    not_anticlockwise: usize,
    /// Rings that cross themselves on the sphere, beyond rounding, the pole one point
    /// (`proper_crossings`): three on each aperture-7 grid, none on aperture 3.
    crossing: usize,
    /// Rings that cross themselves in the plane of the longitudes and latitudes the crate
    /// returns, as a consumer who reads them as planar coordinates sees them: proper crossings
    /// of non-adjacent edges with strict signs, no unwrapping of a longitude and no collapsing
    /// of a pole. On aperture 3 from level 10 (`F0-79-A`, level 10, whose ring names the pole
    /// twice, and `M0-40DF8-A`, level 24, which names it at five longitudes); on aperture 7 the
    /// rings beside the poles at levels 15 to 19.
    crossing_in_degrees: usize,
    /// Rings that hold two points on a pole.
    pole_named_twice: usize,
    /// Rings that hold the same point at two places.
    repeated_point: usize,
    /// Rings at a refinement of 1 with fewer points than the zone has edges: points dropped
    /// where the inverse projection fails, and, on aperture 7, the rings of one point, whose
    /// points the deduplication loses.
    dropped: usize,
    /// Rings with a longitude beyond a half turn, in the engine's radians.
    beyond_half_turn: usize,
    /// Extents with `ll.lon > ur.lon`; the same zones, one for one.
    western_end_the_greater: usize,
    /// Extents that reach a pole exactly; the rings that hold a pole, one for one.
    reach_a_pole: usize,
}

/// What each grid's sample and fixed list meets, measured, asserted exactly (the count of
/// coordinates, which is a floor, is left at zero).
fn expected_kinds<S: Subject>() -> Kinds {
    match S::NAME {
        "IGEO7" => Kinds {
            rings: 5624,
            null_zone: 1,
            undrawn: 5,
            undrawn_without_area: 0,
            drawn: 0,
            ring_not_asked: 0,
            one_point: 3,
            not_anticlockwise: 0,
            crossing: 3,
            pole_named_twice: 0,
            repeated_point: 0,
            dropped: 4,
            beyond_half_turn: 94,
            western_end_the_greater: 94,
            reach_a_pole: 60,
            crossing_in_degrees: 3,
            ..Kinds::default()
        },
        "IVEA7H" => Kinds {
            rings: 5596,
            null_zone: 1,
            undrawn: 5,
            undrawn_without_area: 0,
            drawn: 0,
            ring_not_asked: 0,
            one_point: 3,
            not_anticlockwise: 0,
            crossing: 3,
            pole_named_twice: 0,
            repeated_point: 0,
            dropped: 4,
            beyond_half_turn: 95,
            western_end_the_greater: 95,
            reach_a_pole: 61,
            crossing_in_degrees: 3,
            ..Kinds::default()
        },
        "RTEA7H" => Kinds {
            rings: 5591,
            null_zone: 1,
            undrawn: 5,
            undrawn_without_area: 0,
            drawn: 0,
            ring_not_asked: 0,
            one_point: 3,
            not_anticlockwise: 0,
            crossing: 3,
            pole_named_twice: 0,
            repeated_point: 0,
            dropped: 4,
            beyond_half_turn: 95,
            western_end_the_greater: 95,
            reach_a_pole: 60,
            crossing_in_degrees: 6,
            ..Kinds::default()
        },
        "ISEA3H" => Kinds {
            rings: 5937,
            null_zone: 1,
            undrawn: 32,
            undrawn_without_area: 1,
            drawn: 68,
            ring_not_asked: 34,
            one_point: 0,
            not_anticlockwise: 4,
            crossing: 0,
            pole_named_twice: 38,
            repeated_point: 14,
            dropped: 0,
            beyond_half_turn: 119,
            western_end_the_greater: 119,
            reach_a_pole: 67,
            crossing_in_degrees: 17,
            ..Kinds::default()
        },
        "IVEA3H" => Kinds {
            rings: 5936,
            null_zone: 1,
            undrawn: 32,
            undrawn_without_area: 1,
            drawn: 68,
            ring_not_asked: 34,
            one_point: 0,
            not_anticlockwise: 4,
            crossing: 0,
            pole_named_twice: 36,
            repeated_point: 12,
            dropped: 1,
            beyond_half_turn: 117,
            western_end_the_greater: 117,
            reach_a_pole: 67,
            crossing_in_degrees: 12,
            ..Kinds::default()
        },
        "RTEA3H" => Kinds {
            rings: 5932,
            null_zone: 1,
            undrawn: 32,
            undrawn_without_area: 1,
            drawn: 68,
            ring_not_asked: 34,
            one_point: 0,
            not_anticlockwise: 4,
            crossing: 0,
            pole_named_twice: 38,
            repeated_point: 11,
            dropped: 0,
            beyond_half_turn: 117,
            western_end_the_greater: 117,
            reach_a_pole: 67,
            crossing_in_degrees: 19,
            ..Kinds::default()
        },
        other => panic!("nothing has been measured of the boundary of {other}"),
    }
}

/// The zones without geometry that the sample and the fixed list meet, by the way the engine
/// answers them: the null zone, the identifiers it does not draw, those among them with no
/// area, and the ones it draws.
pub fn expected_without_geometry<S: Subject>() -> [usize; 4] {
    let k = expected_kinds::<S>();
    [k.null_zone, k.undrawn, k.undrawn_without_area, k.drawn]
}

/// The ring of every zone of the sample and of the fixed list at the refinements 0, 1, 2 and 3,
/// and its extent, against the engine's, exactly; the same through `AnyGrid` and `Zone`; the
/// counts of the kinds of ring and of extent, exactly; and every identifier without geometry
/// classified against the live engine.
pub fn the_refined_boundary_and_the_extent_are_the_engines<S: Subject>() {
    assert_eq!(
        rs4dggs::grid::max_edge_refinement(),
        rs4dggs::grid::MAX_EDGE_REFINEMENT,
        "the limit of the edge refinement is lowered: unset RS4DGGS_MAX_EDGE_REFINEMENT to run \
         the suites"
    );
    let g = S::grid();
    let any = get_grid(S::NAME).unwrap();
    let mut floor = Floor::new("boundary and extent", 5000);
    let mut k = Kinds::default();
    for id in geometry_zones::<S>() {
        let text = g.text_id(id);
        let at = format!("{text} ({:#018x}, seed {SEED:#x})", id.0);
        let ours_extent = g.extent(id);
        // Every public path gives the same answer, bit for bit.
        let ext_bits = |e: Option<rs4dggs::Extent>| {
            e.map(|e| [e.ll.lat, e.ll.lon, e.ur.lat, e.ur.lon].map(f64::to_bits))
        };
        assert_eq!(
            ext_bits(g.zone(id).extent()),
            ext_bits(ours_extent),
            "Zone::extent of {at}"
        );
        assert_eq!(
            ext_bits(any.extent(id)),
            ext_bits(ours_extent),
            "AnyGrid::extent of {at}"
        );

        let Some(extent) = ours_extent else {
            // No geometry: the ring, the extent and the area are the empty ones, the six zero
            // points are the vertices, and the engine's side is classified live.
            for r in [0, 1, 2, 3] {
                assert_eq!(
                    g.refined_vertices(id, r),
                    Ok(Vec::new()),
                    "ring of {at} at {r}"
                );
            }
            assert_eq!(
                any.refined_vertices(id, 0),
                Ok(Vec::new()),
                "AnyGrid ring of {at}"
            );
            assert_eq!(
                g.zone(id).refined_vertices(0),
                Ok(Vec::new()),
                "Zone ring of {at}"
            );
            assert_eq!(g.area(id), None, "area of {at}");
            let zeros = g.vertices(id);
            assert!(
                zeros.len() == 6 && zeros.iter().all(|p| (p.lat, p.lon) == (0.0, 0.0)),
                "vertices of {at} are not the six zero points: {zeros:?}"
            );
            match div::engine_answers_without_geometry::<S>(id.0) {
                Ok(NoGeometry::NullZone) => k.null_zone += 1,
                Ok(NoGeometry::Undrawn) => k.undrawn += 1,
                Ok(NoGeometry::UndrawnWithoutArea) => k.undrawn_without_area += 1,
                Ok(NoGeometry::Drawn) => k.drawn += 1,
                Err(why) => panic!("{at}: {why}"),
            }
            k.ring_not_asked +=
                usize::from(S::T::APERTURE == 3 && !o::engine_can_read(S::ORACLE, id.0));
            floor.hit();
            continue;
        };
        assert!(
            S::T::APERTURE != 3 || o::engine_can_read(S::ORACLE, id.0),
            "{at}: this crate has a ring for an identifier the engine cannot read"
        );

        // The extent.
        let theirs = o::extent(S::ORACLE, id.0);
        let ours = [extent.ll.lat, extent.ll.lon, extent.ur.lat, extent.ur.lon];
        assert_eq!(
            ours.map(f64::to_bits),
            theirs.map(|m| (m * RAD2DEG).to_bits()),
            "extent of {at}: ours {ours:?}, the engine's {theirs:?} radians"
        );
        k.coordinates += 4;
        k.western_end_the_greater += usize::from(extent.ll.lon > extent.ur.lon);

        // The ring, at each refinement.
        let edges = if g.is_pentagon(id) { 5 } else { 6 };
        for r in 0u32..=3 {
            let engine = o::refined_vertices(S::ORACLE, id.0, r as i32);
            let want: Vec<(f64, f64)> = engine.iter().copied().map(degrees).collect();
            let ring = g.refined_vertices(id, r).unwrap();
            assert_eq!(
                ring.len(),
                want.len(),
                "points of the ring of {at} at a refinement of {r}"
            );
            for (i, (a, b)) in ring.iter().zip(&want).enumerate() {
                if (a.lat.to_bits(), a.lon.to_bits()) != (b.0.to_bits(), b.1.to_bits()) {
                    let mut worst = 0.0;
                    super::assert_bit_identical(
                        &format!("point {i} of the ring of {at} at a refinement of {r}"),
                        (a.lat, a.lon),
                        *b,
                        &mut worst,
                    );
                }
            }
            k.coordinates += 2 * want.len();
            if r < 2 {
                // Through the other two paths, at the refinements of 0 and 1.
                assert_eq!(
                    bits(&any.refined_vertices(id, r).unwrap()),
                    bits(&ring),
                    "AnyGrid ring of {at} at {r}"
                );
                assert_eq!(
                    bits(&g.zone(id).refined_vertices(r).unwrap()),
                    bits(&ring),
                    "Zone ring of {at} at {r}"
                );
            }
            if r == 1 {
                k.dropped += usize::from(want.len() < edges);
            }
            if r != 0 {
                continue;
            }
            // The kinds of ring, at the automatic refinement.
            k.rings += 1;
            k.one_point += usize::from(want.len() == 1);
            let poles = want.iter().filter(|p| p.0.abs() == 90.0).count();
            match winding(&want) {
                Winding::Anticlockwise => {}
                Winding::NoArea if want.len() == 1 => {}
                Winding::NoArea => {
                    // The four rings of level 33 on aperture 3 that every point of which but two
                    // is on the pole: no ring that runs clockwise passes here.
                    assert!(
                        S::T::APERTURE == 3 && g.resolution(id) == 33 && poles + 2 == want.len(),
                        "the ring of {at} encloses no area and has {} points, {poles} on a pole",
                        want.len()
                    );
                    k.not_anticlockwise += 1;
                }
                w => panic!("the ring of {at} runs {w:?}"),
            }
            k.crossing_in_degrees += usize::from(crossings_in_degrees(&ring) > 0);
            k.crossing += usize::from(proper_crossings(&want) > 0);
            k.pole_named_twice += usize::from(poles >= 2);
            let mut sorted: Vec<(u64, u64)> = want
                .iter()
                .map(|p| (p.0.to_bits(), p.1.to_bits()))
                .collect();
            sorted.sort_unstable();
            sorted.dedup();
            k.repeated_point += usize::from(sorted.len() < want.len());
            let beyond = engine.iter().any(|p| p.1.abs() > std::f64::consts::PI);
            k.beyond_half_turn += usize::from(beyond);
            // A ring beyond a half turn is a zone over the antimeridian, one for one.
            assert_eq!(
                beyond,
                extent.ll.lon > extent.ur.lon,
                "the ring of {at} holds a longitude beyond 180 degrees and its extent is over \
                 the antimeridian, or neither"
            );
            // A ring that holds a pole is an extent that reaches it, one for one.
            let reaches = extent.ur.lat == 90.0 || extent.ll.lat == -90.0;
            assert_eq!(
                poles > 0,
                reaches,
                "the ring of {at} holds a pole, or the extent reaches one"
            );
            k.reach_a_pole += usize::from(reaches);
        }
        floor.hit();
    }
    eprintln!("boundary and extent on {}: {k:?}", S::NAME);
    assert!(
        k.coordinates >= 900_000,
        "only {} coordinates compared",
        k.coordinates
    );
    k.coordinates = 0;
    assert_eq!(k, expected_kinds::<S>(), "the kinds met on {}", S::NAME);
}

/// The ring at a refinement of 2,000, 20,000 and 100,000, at one zone of each kind of ring
/// the engine gives, against the engine's, exactly: the edge is divided into that many
/// parts, with no ceiling in the engine, and this crate is compared with it up to its own.
pub fn the_finest_refinements_are_the_engines<S: Subject>() {
    assert_eq!(
        rs4dggs::grid::max_edge_refinement(),
        rs4dggs::grid::MAX_EDGE_REFINEMENT,
        "the limit of the edge refinement is lowered: unset RS4DGGS_MAX_EDGE_REFINEMENT to run \
         the suites"
    );
    let g = S::grid();
    let any = get_grid(S::NAME).unwrap();
    let kinds = match S::T::APERTURE {
        7 => &KINDS_OF_APERTURE_7,
        3 => &KINDS_OF_APERTURE_3,
        a => panic!("no zones of each kind are recorded for aperture {a}"),
    };
    let refinements = [2_000u32, 20_000, 100_000];
    let mut floor = Floor::new("fine refinements", kinds.len() * refinements.len());
    let mut points = 0usize;
    for &(kind, id) in kinds {
        let id = ZoneId(id);
        let at = format!("{kind}, {} ({:#018x})", g.text_id(id), id.0);
        for r in refinements {
            let theirs = o::refined_vertices(S::ORACLE, id.0, r as i32);
            let ours = g.refined_vertices(id, r).unwrap();
            assert_eq!(
                ours.len(),
                theirs.len(),
                "points of {at} at a refinement of {r}"
            );
            for (i, (a, b)) in ours.iter().zip(&theirs).enumerate() {
                assert_eq!(
                    (a.lat.to_bits(), a.lon.to_bits()),
                    ((b.0 * RAD2DEG).to_bits(), (b.1 * RAD2DEG).to_bits()),
                    "point {i} of {at} at a refinement of {r}"
                );
            }
            assert!(
                !theirs.is_empty(),
                "the engine's ring of {at} at {r} is empty"
            );
            points += theirs.len();
            if r == 100_000 {
                assert_eq!(
                    bits(&any.refined_vertices(id, r).unwrap()),
                    bits(&ours),
                    "AnyGrid ring of {at}"
                );
            }
            floor.hit();
        }
    }
    // The limit is the ceiling, and a finer refinement is refused with both numbers.
    let zone = ZoneId(kinds[0].1);
    for r in [100_001u32, 1 << 20, u32::MAX] {
        let refused = Err(Error::EdgeRefinementOutOfRange {
            refinement: r,
            max: rs4dggs::grid::MAX_EDGE_REFINEMENT,
        });
        assert_eq!(g.refined_vertices(zone, r), refused);
        assert_eq!(any.refined_vertices(zone, r), refused);
        assert_eq!(g.zone(zone).refined_vertices(r), refused);
    }
    eprintln!("fine refinements on {}: {points} points compared", S::NAME);
}

/// Hostile identifiers for the sweeps that must not panic: the null zone, the identifiers at
/// the ends of the range, those with no cell behind them, and 3,000 drawn at random from
/// every level, real zones among them as often as junk.
pub(super) fn hostile_identifiers<S: Subject>() -> Vec<u64> {
    let mut ids = vec![
        0u64,
        1,
        u64::MAX,
        1 << 63,
        (1 << 63) - 1,
        0x5555_5555_5555_5555,
    ];
    ids.extend(div::FABRICATED);
    let mut state = SEED;
    for _ in 0..3_000 {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        if S::T::APERTURE == 7 {
            // Part of the draws keep the bits below a level, as an identifier does that has used
            // few of its digit slots, so that real zones are met as often as junk.
            let level = (state >> 59) % 34;
            let kept = u64::MAX.checked_shr(4 + 2 * level as u32).unwrap_or(0);
            ids.push(if state & 1 == 0 {
                state ^ 0x9e37
            } else {
                (state << 3) | kept
            });
        } else {
            // The packing of aperture 3: a letter, a root, an index and a sub-hex, the index
            // below the count of the zones of a rhombus, and sometimes beyond it.
            let (letter, root, sub_hex) =
                ((state >> 59) % 18, (state >> 55) % 13, (state >> 53) & 3);
            let index = (state >> 3)
                % 9u64.pow(letter as u32 % 17).max(2)
                % (1u64 << ((state % 5 + 1) * 6));
            ids.push((letter << 57) | (root << 53) | (index << 2) | sub_hex);
        }
    }
    ids
}

/// The three new answers of a zone, asked of hostile identifiers and arguments: the null zone,
/// the identifiers at the ends of the range, those with no cell behind them, and 3,000 drawn at
/// random from every level; at refinements from 0 to far beyond the ceiling. None panics, an
/// identifier has a ring, an extent and an area together or none of the three, and a refinement
/// above the ceiling is refused whatever the identifier.
pub fn the_new_methods_do_not_panic_on_hostile_input<S: Subject>() {
    assert_eq!(
        rs4dggs::grid::max_edge_refinement(),
        rs4dggs::grid::MAX_EDGE_REFINEMENT,
        "the limit of the edge refinement is lowered: unset RS4DGGS_MAX_EDGE_REFINEMENT to run \
         the suites"
    );
    let g = S::grid();
    let any = get_grid(S::NAME).unwrap();
    let ids = hostile_identifiers::<S>();
    let floor_with = ids.len() / 20;

    let (mut with_geometry, mut without) = (0usize, 0usize);
    for id in ids {
        let id = ZoneId(id);
        let at = format!("{:#018x}", id.0);
        let extent = g.extent(id);
        let area = g.area(id);
        assert_eq!(any.extent(id).is_some(), extent.is_some(), "{at}");
        assert_eq!(
            any.area(id).map(f64::to_bits),
            area.map(f64::to_bits),
            "{at}"
        );
        assert_eq!(g.zone(id).area().is_some(), area.is_some(), "{at}");
        assert_eq!(
            extent.is_some(),
            area.is_some(),
            "{at}: an extent and an area go together"
        );
        for r in [0u32, 1, 2, 7, 100] {
            let ring = g.refined_vertices(id, r).unwrap();
            assert_eq!(ring.is_empty(), extent.is_none(), "{at} at {r}");
            assert_eq!(
                any.refined_vertices(id, r).unwrap().len(),
                ring.len(),
                "{at} at {r}"
            );
            assert!(
                ring.iter().all(|p| !p.lat.is_nan() && !p.lon.is_nan()),
                "{at} at {r}: a NaN"
            );
        }
        for r in [100_001u32, u32::MAX] {
            assert!(
                matches!(
                    g.refined_vertices(id, r),
                    Err(Error::EdgeRefinementOutOfRange { refinement, .. }) if refinement == r
                ),
                "{at} at {r}"
            );
        }
        if extent.is_some() {
            with_geometry += 1;
        } else {
            without += 1;
        }
    }
    assert!(
        with_geometry >= floor_with && without >= floor_with,
        "{with_geometry} identifiers with geometry, {without} without"
    );
    eprintln!(
        "hostile input on {}: {with_geometry} with geometry, {without} without",
        S::NAME
    );
}
