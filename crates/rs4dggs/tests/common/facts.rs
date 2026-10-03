//! The area of a zone and the facts of a grid against the live engine: `getZoneArea` for every
//! zone of the suites' sample, `countZones`, `getRefZoneArea` and `getMetersPerSubZoneFromLevel`
//! at every level a grid has, the grid's constants, and the two functions that give a level from
//! an area, `getLevelFromRefZoneArea` and `getLevelFromMetersPerSubZone`, which the crate
//! mirrors for every input and which a census of areas tests with every kind: spread
//! logarithmically over thirty-two decades, each level's own area and its neighbours to the bit,
//! and the special values.
//!
//! Every comparison is exact, and every count has a floor. The engine is asked for its own
//! answer at levels and areas the crate refuses (its tables run out at 25 and 41 levels, and a
//! `getRefZoneArea` beyond the maximum counts a hexagon's sub-zones), so that the sweep reaches
//! the arithmetic that the crate keeps behind its bounds. Nothing here is skipped.
use std::collections::BTreeSet;

use dggal_oracle as o;
use rs4dggs::indexings::Z7;
use rs4dggs::topologies::HexA7;
use rs4dggs::{Address, Error, Indexing, Topology, ZoneId, get_grid};

use super::Floor;
use super::expected_divergences::{self as div, NoGeometry};
use super::geometry::geometry_zones;
use super::subject::Subject;
use super::suite::SEED;

/// `wholeWorld.geodeticArea`, as the crate pins it: the census builds a hexagon's area from it.
const EARTH_AREA: f64 = f64::from_bits(0x42fc_fe6e_831e_bb8c);

/// `getZoneArea`, of every zone of the main sample and of the fixed list of
/// [`super::geometry::fixed_zones`]: this crate's area is the engine's double, to the bit, for a
/// pentagon as for a hexagon. A zone with no area is exactly a zone with no geometry (by the
/// crate's rule for a zone without geometry: its vertices are the six zero points, its ring the empty one and its extent none), and each is
/// a departure of this crate, classified against the live engine and counted exactly:
/// the null zone, for which the engine answers `+inf`; an identifier it does not draw, for which
/// it answers the area of its level's formula, or `+inf` where that formula's divisor is zero;
/// and, on aperture 3, the sub-hexes C and D of a polar root, which it draws and gives an area.
pub fn the_area_of_every_sampled_zone_is_the_engines<S: Subject>() {
    let mut floor = Floor::new("area", 5000);
    let any = get_grid(S::NAME).unwrap();
    let mut pentagons = 0usize;
    let mut without = [0usize; 4];
    for id in geometry_zones::<S>() {
        let text = S::grid().text_id(id);
        let ours = S::grid().area(id);
        let theirs = o::area(S::ORACLE, id.0);
        assert_eq!(
            S::grid().zone(id).area().map(f64::to_bits),
            ours.map(f64::to_bits),
            "Zone::area of {text}"
        );
        assert_eq!(
            any.area(id).map(f64::to_bits),
            ours.map(f64::to_bits),
            "AnyGrid::area of {text}"
        );
        match ours {
            Some(ours) => {
                assert_eq!(
                    ours.to_bits(),
                    theirs.to_bits(),
                    "area of {text} (seed {SEED:#x}): ours {ours:e}, the engine's {theirs:e}"
                );
                assert!(ours.is_finite() && ours > 0.0, "area of {text}");
                assert!(
                    S::grid()
                        .vertices(id)
                        .iter()
                        .any(|p| (p.lat, p.lon) != (0.0, 0.0)),
                    "{text} has an area and its vertices are the six zero points"
                );
                pentagons += usize::from(S::grid().is_pentagon(id));
            }
            None => {
                // No area exactly where there is no geometry.
                let vertices = S::grid().vertices(id);
                assert!(
                    vertices.len() == 6 && vertices.iter().all(|p| (p.lat, p.lon) == (0.0, 0.0)),
                    "{text} has no area and its vertices are not the six zero points"
                );
                assert_eq!(
                    S::grid().extent(id),
                    None,
                    "{text} has no area and has an extent"
                );
                match div::engine_answers_without_geometry::<S>(id.0) {
                    Ok(NoGeometry::NullZone) => {
                        assert_eq!(theirs, f64::INFINITY, "the engine's area of the null zone");
                        without[0] += 1;
                    }
                    Ok(kind @ (NoGeometry::Undrawn | NoGeometry::Drawn)) => {
                        assert!(theirs.is_finite() && theirs > 0.0, "{text}: {theirs:e}");
                        without[if kind == NoGeometry::Undrawn { 1 } else { 3 }] += 1;
                    }
                    Ok(NoGeometry::UndrawnWithoutArea) => {
                        assert_eq!(theirs, f64::INFINITY, "{text}");
                        without[2] += 1;
                    }
                    Err(why) => panic!("{text}: {why}"),
                }
            }
        }
        floor.hit();
    }
    assert!(pentagons >= 24, "only {pentagons} pentagons compared");
    // The null zone, the undrawn, the undrawn with no area, the drawn: the counts that the
    // boundary's test holds.
    assert_eq!(
        without,
        super::geometry::expected_without_geometry::<S>(),
        "{without:?}"
    );
    eprintln!(
        "area: {pentagons} pentagons; no area (null zone, undrawn, undrawn without area, drawn): {without:?}"
    );
}

/// An identifier with no cell behind it has no area here, and the engine's answers for the four
/// of [`div::FABRICATED`] are characterised exactly. The engine reads the level out of the bits of
/// such an identifier as it would of any other: on aperture 7 a base cell of 12 to 14 is a
/// pentagon of level 0 for it, and on aperture 3 the fields give levels 15 and 31, hexagons, and
/// then a level whose count settles at 2, the divisor 0. Its null zone is `+inf` on both.
pub fn an_identifier_with_no_cell_has_no_area<S: Subject>() {
    let mut floor = Floor::new("area of an identifier with no cell", div::FABRICATED.len());
    let infinity = f64::INFINITY.to_bits();
    let engine = match S::T::APERTURE {
        7 => [
            0x42c3_5449_acbf_27b3,
            0x42c3_5449_acbf_27b3,
            0x42c3_5449_acbf_27b3,
            infinity,
        ],
        3 => [
            0x414b_1ed7_8235_ceda,
            0x3fb5_23dd_dd42_e91b,
            infinity,
            infinity,
        ],
        a => panic!("no answers of the engine are recorded for aperture {a}"),
    };
    for (id, bits) in div::FABRICATED.into_iter().zip(engine) {
        assert_eq!(S::grid().area(ZoneId(id)), None, "{id:#018x}");
        assert_eq!(S::grid().zone(ZoneId(id)).area(), None, "{id:#018x}");
        assert_eq!(
            o::area(S::ORACLE, id).to_bits(),
            bits,
            "the engine's area of {id:#018x} on {}",
            S::NAME
        );
        floor.hit();
    }
}

/// On aperture 7 two kinds of identifier are odd to the engine. One it cannot read back as text,
/// `01222222222222222220` among them, a neighbour it lists and then refuses: it has its ordinary
/// geometry, and so its ordinary area, which this crate keeps. The other has twenty
/// digits and no geometry, so no area here, where the engine answers the area of its level's
/// formula, a finite number for a zone it cannot draw: `earthArea` over `countZones(20) - 2`, and
/// five sixths of it for the pentagon.
pub fn identifiers_the_engine_cannot_read_or_cannot_draw_have_the_areas_decided<
    S: Subject<T = HexA7, I = Z7>,
>() {
    // The unreadable identifier is a hexagon at level 18.
    let unreadable = ZoneId(0x1492_4924_9249_243f);
    let theirs = o::area(S::ORACLE, unreadable.0);
    assert_eq!(theirs.to_bits(), 0x3fa0_098c_9747_bc25);
    assert_eq!(
        S::grid().area(unreadable).map(f64::to_bits),
        Some(theirs.to_bits()),
        "the area of an identifier the engine cannot read back"
    );

    // `countZones(20)` is `10 * 7^20 + 2`, which the crate does not offer past level 19: the
    // hexagons at level 20 number `10 * 7^20`.
    let hexagons = 797_922_662_976_120_010u64 as f64;
    let one_digits = Z7::encode(&Address::new(1, &[1; 20]));
    for (id, pentagon, bits) in [
        (
            ZoneId(0x329c_b814_e5c0_a72e),
            false,
            0x3f44_f259_9691_f084u64,
        ),
        (one_digits, false, 0x3f44_f259_9691_f084),
        (ZoneId(0), true, 0x3f41_749f_fd79_9dc4),
    ] {
        assert_eq!(S::grid().area(id), None, "{:#018x}", id.0);
        let formula = EARTH_AREA / hexagons * if pentagon { 5.0 / 6.0 } else { 1.0 };
        let theirs = o::area(S::ORACLE, id.0);
        assert_eq!(
            theirs.to_bits(),
            bits,
            "the engine's area of {:#018x}",
            id.0
        );
        assert_eq!(theirs.to_bits(), formula.to_bits(), "{:#018x}", id.0);
    }
}

/// The grid's constants, and `countZones`, `getRefZoneArea` and `getMetersPerSubZoneFromLevel`
/// at every level the grid has, and at every pair of a level and a depth within it. A level
/// beyond the maximum is refused with the error that says so.
pub fn the_grid_facts_are_the_engines<S: Subject>() {
    let g = S::grid();
    let any = get_grid(S::NAME).unwrap();
    let max = g.max_resolution();
    // Six facts; the count and the reference area at each level; the size at each level and depth.
    let levels = usize::from(max) + 1;
    let mut floor = Floor::new("grid facts", 6 + 2 * levels + levels * (levels + 1) / 2);

    for (what, ours, theirs) in [
        (
            "refinement ratio",
            g.refinement_ratio(),
            o::refinement_ratio(S::ORACLE),
        ),
        ("max parents", g.max_parents(), o::max_parents(S::ORACLE)),
        ("max children", g.max_children(), o::max_children(S::ORACLE)),
        (
            "max neighbours",
            g.max_neighbors(),
            o::max_neighbors(S::ORACLE),
        ),
        ("64K depth", g.depth_64k(), o::depth_64k(S::ORACLE)),
        ("max depth", g.max_depth(), o::max_depth(S::ORACLE)),
    ] {
        assert_eq!(i32::from(ours), theirs, "{what} of {}", S::NAME);
        floor.hit();
    }
    assert_eq!(
        (
            any.refinement_ratio(),
            any.max_parents(),
            any.max_children()
        ),
        (g.refinement_ratio(), g.max_parents(), g.max_children())
    );
    assert_eq!(
        (any.max_neighbors(), any.depth_64k(), any.max_depth()),
        (g.max_neighbors(), g.depth_64k(), g.max_depth())
    );
    assert_eq!(i32::from(S::T::APERTURE), o::refinement_ratio(S::ORACLE));

    for level in 0..=max {
        let l = i32::from(level);
        assert_eq!(
            g.count_zones(level).unwrap(),
            o::count_zones(S::ORACLE, l),
            "countZones({level}) on {}",
            S::NAME
        );
        assert_eq!(
            g.ref_zone_area(level).unwrap().to_bits(),
            o::ref_zone_area(S::ORACLE, l).to_bits(),
            "getRefZoneArea({level}) on {}",
            S::NAME
        );
        assert_eq!(any.count_zones(level), g.count_zones(level));
        assert_eq!(
            any.ref_zone_area(level).unwrap().to_bits(),
            g.ref_zone_area(level).unwrap().to_bits()
        );
        floor.hit();
        floor.hit();
        for depth in 0..=(max - level) {
            let theirs = o::meters_per_sub_zone(S::ORACLE, l, i32::from(depth));
            assert_eq!(
                g.meters_per_sub_zone(level, depth).unwrap().to_bits(),
                theirs.to_bits(),
                "getMetersPerSubZoneFromLevel({level}, {depth}) on {}",
                S::NAME
            );
            assert_eq!(
                any.meters_per_sub_zone(level, depth).unwrap().to_bits(),
                theirs.to_bits()
            );
            floor.hit();
        }
    }

    // The tables of the engine run out past the maximum; this crate stops at it.
    let refused = |res: u8| Error::ResolutionOutOfRange { res, max };
    assert_eq!(g.count_zones(max + 1), Err(refused(max + 1)));
    assert_eq!(g.ref_zone_area(max + 1), Err(refused(max + 1)));
    assert_eq!(g.meters_per_sub_zone(max, 1), Err(refused(max + 1)));
    assert_eq!(g.meters_per_sub_zone(255, 255), Err(refused(255)));
    assert_eq!(any.count_zones(255), Err(refused(255)));
    // ... where the engine answers, with counts that wrap and then settle at 2.
    assert_eq!(o::count_zones(S::ORACLE, i32::from(max) + 30), 2);
    assert_eq!(o::count_zones(S::ORACLE, -1), 2);
}

/// The areas of the census: spread logarithmically from 1e-14 to 1e18 square metres, each level's
/// own reference area and a hexagon's area at it three units in the last place either side, from
/// below level 0 to above the last level the function answers, and the special values.
fn census_areas<S: Subject>() -> Vec<f64> {
    let max = i32::from(S::grid().max_resolution());
    let d64 = i32::from(S::grid().depth_64k());
    let mut areas: Vec<f64> = (0..200_000)
        .map(|i| 10f64.powf(-14.0 + 32.0 * f64::from(i) / 200_000.0))
        .collect();
    let beside = |x: f64| (-3i64..=3).map(move |k| f64::from_bits((x.to_bits() as i64 + k) as u64));
    for level in -3..=(max + d64 + 3) {
        areas.extend(beside(o::ref_zone_area(S::ORACLE, level)));
        // The area of a hexagon at the level, from the engine's own count of zones.
        let hexagons = o::count_zones(S::ORACLE, level).wrapping_sub(2);
        areas.extend(beside(EARTH_AREA / hexagons as f64));
    }
    areas.extend([
        0.0,
        -0.0,
        -1.0,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NAN,
        5e-324,
        f64::MAX,
        f64::MIN_POSITIVE,
        EARTH_AREA,
        1e-300,
    ]);
    areas
}

/// `getLevelFromRefZoneArea` and `getLevelFromMetersPerSubZone`, over [`census_areas`]: the
/// first at each area, the second at its square root (so NaN for a negative area) and three
/// depths.
pub fn the_levels_from_an_area_are_the_engines<S: Subject>() {
    let g = S::grid();
    let any = get_grid(S::NAME).unwrap();
    let top = g.max_resolution() + g.depth_64k();
    let areas = census_areas::<S>();
    let mut floor = Floor::new("levels from an area", 4 * 200_000);
    let mut answered = BTreeSet::new();
    for a in areas {
        let theirs = o::level_from_ref_zone_area(S::ORACLE, a);
        let ours = i32::from(g.level_from_ref_zone_area(a));
        assert_eq!(
            ours,
            theirs,
            "getLevelFromRefZoneArea({a:e}) on {}",
            S::NAME
        );
        assert_eq!(i32::from(any.level_from_ref_zone_area(a)), ours);
        answered.insert(theirs);
        floor.hit();
        let metres = a.sqrt();
        for depth in [0u8, 2, 40] {
            let theirs = o::level_from_meters_per_sub_zone(S::ORACLE, metres, i32::from(depth));
            let ours = i32::from(g.level_from_meters_per_sub_zone(metres, depth));
            assert_eq!(
                i32::from(any.level_from_meters_per_sub_zone(metres, depth)),
                ours,
                "AnyGrid::level_from_meters_per_sub_zone({metres:e}, {depth}) on {}",
                S::NAME
            );
            assert_eq!(
                ours,
                theirs,
                "getLevelFromMetersPerSubZone({metres:e}, {depth}) on {}",
                S::NAME
            );
            floor.hit();
        }
    }
    // The special values: 0 and NaN reach the top, a negative or infinite area is level 0.
    for (a, level) in [
        (0.0, top),
        (-0.0, 0),
        (f64::NAN, top),
        (-1.0, 0),
        (f64::INFINITY, 0),
        (f64::NEG_INFINITY, 0),
    ] {
        assert_eq!(g.level_from_ref_zone_area(a), level, "{a:e}");
        assert_eq!(
            o::level_from_ref_zone_area(S::ORACLE, a),
            i32::from(level),
            "{a:e}"
        );
    }
    // The answers span the range, but for the levels whose wrapped count of zones lies below that
    // of a level before them, which the loop can never stop at: 22 on aperture 7, and 39 to 42 on
    // aperture 3, where the counts of 41 and 42 are 2. The top, 25 or 43, answers an area too small
    // for any count, and the areas of 0 and NaN.
    let expected: BTreeSet<i32> = match S::T::APERTURE {
        7 => (0..=25).filter(|&l| l != 22).collect(),
        3 => (0..=38).chain([43]).collect(),
        a => panic!("no answers of the engine are recorded for aperture {a}"),
    };
    assert_eq!(
        answered,
        expected,
        "the levels the engine answers on {}",
        S::NAME
    );
    eprintln!(
        "levels from an area on {}: the engine answered {answered:?}",
        S::NAME
    );
}
