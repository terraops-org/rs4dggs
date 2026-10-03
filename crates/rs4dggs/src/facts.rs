//! The area of a zone and the facts of a grid, as DGGAL's `DGGRS` computes them: the count of
//! zones at a level, the reference area of a level, the level that matches an area, and the
//! grid's constants.
//!
//! Each function keys on the aperture alone, 7 for `RhombicIcosahedral7H` (and its Z7 class) and
//! 3 for `RhombicIcosahedral3H`, which is all that separates the six grids here: the three
//! projections of an aperture give the same doubles, to the bit (the oracle suites assert it).
//! The arithmetic is that of the library `libdggal.so`, BuildID
//! `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, in the order gcc compiled it; a `// faithful:`
//! comment names the eC and the address wherever that order is not the source's, or where the
//! source holds something this port leaves out.
//!
//! Source: `dggrs.ec`, `RI7H.ec`, `RI7H_Z7.ec`, `RI3H.ec` and `GeoExtent.ec` of DGGAL v0.0.6.

use crate::math;
use crate::topologies::pow3;

/// `wholeWorld.geodeticArea`, the area of the whole ellipsoid in square metres, which every
/// function here divides. It is pinned as bits because neither of its two source forms gives it:
/// the eC's order, evaluated with glibc, gives `0x42fcfe6e831ebb8b`, one unit in the last place
/// below, and `4 * Pi * R^2` with the authalic radius gives `...bb87`. gcc rewrote the
/// unguarded `geodeticArea` (`GeoExtent.ec:108-128`, `0x429e0`) into five operations on six
/// folded constants, which the tests of this module rebuild; the double is the one the engine
/// holds in its `earthArea` statics (`0x56ac0` and its siblings).
pub(crate) const EARTH_AREA: f64 = f64::from_bits(0x42fc_fe6e_831e_bb8c);

/// `5/6.0`, by which `getZoneArea` multiplies a pentagon's area (`0x350bd`).
const FIVE_SIXTHS: f64 = f64::from_bits(0x3fea_aaaa_aaaa_aaab);

/// `log(65536)`, the numerator of `get64KDepth` (`dggrs.ec:224-227`), which gcc folded at
/// `0x448d0`. Inert: `16 * ln 2` gives the same double, so the folding moves no answer.
const LN_65536: f64 = f64::from_bits(0x4026_2e42_fefa_39ef);

/// The eC's `powersOf7` (`RI7H.ec:837-841`, `0x562a0`): `7^0` to `7^24`, the last two modulo
/// `2^64`.
const POWERS_OF_7: [u64; 25] = [
    1,
    7,
    49,
    343,
    2_401,
    16_807,
    117_649,
    823_543,
    5_764_801,
    40_353_607,
    282_475_249,
    1_977_326_743,
    13_841_287_201,
    96_889_010_407,
    678_223_072_849,
    4_747_561_509_943,
    33_232_930_569_601,
    232_630_513_987_207,
    1_628_413_597_910_449,
    11_398_895_185_373_143,
    79_792_266_297_612_001,
    558_545_864_083_284_007,
    3_909_821_048_582_988_049,
    8_922_003_266_371_364_727,
    7_113_790_643_470_898_241,
];

/// The engine's `countZones` on aperture 3 at levels 34 to 42, which lie beyond its table of
/// powers of 3. There the engine computes `10 * (uint64)(exp(level * ln 3) + 0.1) + 2`, and the
/// rounding of `exp` makes the count at 34 differ from `10 * 3^34 + 2` (that is
/// `166771816996665692`); the counts at 39 and 40 wrap modulo `2^64`; and from 41 the double is
/// beyond the range of `uint64`, which gcc's conversion answers with 0, so that the count is 2.
/// They are pinned as the integers the engine answers read from it, so that no
/// `exp` of this machine's libm enters; `the_pinned_counts_and_every_level_are_the_engines`
/// compares each with the live engine.
const A3_COUNTS_34_TO_42: [u64; 9] = [
    166_771_816_996_666_582,
    500_315_450_989_998_562,
    1_500_946_352_969_992_002,
    4_502_839_058_909_997_442,
    13_508_517_176_729_960_962,
    3_632_063_382_770_684_930,
    10_896_190_148_312_633_346,
    2,
    2,
];

/// Whether the aperture is 7 and not 3. The trait that carries the aperture is sealed, and these
/// are its two.
///
/// No eC counterpart: the eC tells the two apertures apart by their classes, `RhombicIcosahedral7H`
/// and `RhombicIcosahedral3H`, where this module takes the number.
fn is_a7(aperture: u8) -> bool {
    debug_assert!(aperture == 3 || aperture == 7);
    aperture == 7
}

/// `getMaxDGGRSZoneLevel`: 19 on aperture 7 (`RI7H.ec:47`), 33 on aperture 3 (`RI3H.ec:42`).
pub(crate) fn max_level(aperture: u8) -> u8 {
    if is_a7(aperture) { 19 } else { 33 }
}

/// `getMaxParents`: 2 on aperture 7 (`RI7H.ec:49`), 3 on aperture 3 (`RI3H.ec:44`).
pub(crate) fn max_parents(aperture: u8) -> u8 {
    if is_a7(aperture) { 2 } else { 3 }
}

/// `getMaxChildren`: 13 on aperture 7 (`RI7H.ec:51`), 7 on aperture 3 (`RI3H.ec:46`).
pub(crate) fn max_children(aperture: u8) -> u8 {
    if is_a7(aperture) { 13 } else { 7 }
}

/// `getMaxNeighbors`: 6 on both apertures (`RI7H.ec:50`, `RI3H.ec:45`).
pub(crate) fn max_neighbors(_aperture: u8) -> u8 {
    6
}

/// `get64KDepth` (`dggrs.ec:224-227`): the depth at which a zone has about 65,536 sub-zones,
/// `(int)(log(65536) / log(ratio) + 0.5)`: 6 on aperture 7, 10 on aperture 3.
pub(crate) fn depth_64k(aperture: u8) -> u8 {
    // faithful: gcc folded `log(65536)` into the constant `LN_65536` (`0x448d0`), and divides it
    // at run time by `log` of the refinement ratio, which is the aperture.
    (LN_65536 / math::ln(f64::from(aperture)) + 0.5) as u8
}

/// The count of sub-zones of a hexagon `depth` levels below it, `getSubZonesCount` of the
/// engine's zone classes (`RI7H_Z7.ec:492-500` on aperture 7, `RI3H.ec:2198-2202` on aperture 3)
/// with six points, for which their `(nHexSubZones * nPoints + 5) / 6` is `nHexSubZones`. Only
/// `max_depth` asks, at depths to 19, within both tables of powers.
fn count_sub_zones_hexagon(aperture: u8, depth: u8) -> u64 {
    let d = usize::from(depth);
    if depth == 0 {
        1
    } else if is_a7(aperture) {
        POWERS_OF_7[d]
            + if depth & 1 != 0 {
                5 * POWERS_OF_7[(d - 1) / 2] + 1
            } else {
                POWERS_OF_7[d / 2] - 1
            }
    } else {
        pow3(d as u64) + pow3(d.div_ceil(2) as u64) + 1
    }
}

/// `getMaxDepth` (`dggrs.ec:230-241`): the deepest relative depth at which a zone's sub-zones
/// number at most `2^28`, found from twice the 64K depth less one, counting a hexagon's
/// sub-zones: 9 on aperture 7, 17 on aperture 3.
pub(crate) fn max_depth(aperture: u8) -> u8 {
    // The eC counts the sub-zones of the zone at (0, 10) at level 2, chosen "avoiding ISEA3H
    // pentagon at 0,0 at level 1 and 2": a hexagon on all six grids, as a test of this module
    // checks.
    let mut depth = 2 * depth_64k(aperture) - 1;
    while count_sub_zones_hexagon(aperture, depth) << 4 > 1 << 32 {
        depth -= 1;
    }
    depth
}

/// `countZones` (`RI7H.ec:42-45`, `0x1e0e0`; `RI3H.ec:37-40`, `0xb4f0`): `10 * POW(level) + 2`,
/// with the wrap of `uint64`. The engine answers the same for any `int`, and this answers for any
/// `u8`: past its tables the count is the pinned one or 2.
pub(crate) fn count_zones(aperture: u8, level: u8) -> u64 {
    // faithful: the eC's `POW7` and `POW3` read a table of powers below its length (`powersOf7`
    // at `0x562a0`, 25 entries; `powersOf3` at `0x56180`, 34), and beyond it gcc turned
    // `pow(k, level)` into `exp(level * ln k) + 0.1` converted to `uint64`. This port has no
    // `exp`. Aperture 7 beyond its table is 2 at every level, `7^25` being past `2^64`, where the
    // conversion answers 0; aperture 3 beyond its table is the pinned counts, and 2 from level 43.
    let power = if is_a7(aperture) {
        match POWERS_OF_7.get(usize::from(level)) {
            Some(&p) => p,
            None => return 2,
        }
    } else if level <= 33 {
        pow3(u64::from(level))
    } else {
        return A3_COUNTS_34_TO_42
            .get(usize::from(level) - 34)
            .copied()
            .unwrap_or(2);
    };
    // faithful: `10 * p + 2` as two `lea`, so that the product wraps modulo `2^64`.
    power.wrapping_mul(10).wrapping_add(2)
}

/// `getZoneArea` (`RI7H_Z7.ec:514-525`, `0x34f90`; `RI3H.ec:66-88`, `0xb580`): the area of a zone
/// of the given level in square metres, `earthArea / (countZones(level) - 2)`, and five sixths of
/// that for a pentagon. The Z7 override that the three aperture-7 grids run is unguarded and
/// nothing in it is re-associated; the aperture-3 function is guarded and compiled as written.
/// Both perform the same three operations in the same order.
pub(crate) fn zone_area(aperture: u8, level: u8, pentagon: bool) -> f64 {
    let hexagons = count_zones(aperture, level).wrapping_sub(2);
    let area = EARTH_AREA / hexagons as f64;
    if pentagon { area * FIVE_SIXTHS } else { area }
}

/// `getRefZoneArea` (`dggrs.ec:266-284`, `0x44a30` to `0x44c8f`), for the levels 0 to the grid's
/// maximum: `earthArea / countZones(level)`, which is neither a hexagon's area nor a pentagon's.
pub(crate) fn ref_zone_area(aperture: u8, level: u8) -> f64 {
    // faithful: the eC takes the count at `min(level, max)` and divides by it where `level` is
    // below `max`. At `level >= max` it divides by `zc * countSubZones(testZone, level - max)`
    // (the branch at `0x44bc6`). On this domain the only such level is `max`, the count of
    // sub-zones at depth 0 is 1, and `earthArea / (zc * 1.0)` is `earthArea / zc` to the bit, so
    // the test zone and its longitude (`0x49eb0`) are not needed. The eC's test `zc == 0` is never
    // true here, the counts being at least 12.
    EARTH_AREA / count_zones(aperture, level) as f64
}

/// `getMetersPerSubZoneFromLevel` (`dggrs.ec:305-308`, `0x44cf0`) at the level the parent's and
/// the depth make: the square root of that level's reference area.
pub(crate) fn meters_per_sub_zone(aperture: u8, level: u8) -> f64 {
    math::sqrt(ref_zone_area(aperture, level))
}

/// `getLevelFromRefZoneArea` (`dggrs.ec:250-265`, `0x44d10` to `0x44e4f`): the first level whose
/// count of zones reaches `earthArea` divided by the area, and `max + get64KDepth` (25 on aperture
/// 7, 43 on aperture 3) where none does, which is the answer for an area of 0 or NaN. A negative
/// or infinite area answers 0. It mirrors the engine for every `f64`, and is not the inverse of
/// [`ref_zone_area`]: the reference area of level 9 on ISEA3H answers 10.
pub(crate) fn level_from_ref_zone_area(aperture: u8, square_metres: f64) -> u8 {
    let max_level = max_level(aperture) + depth_64k(aperture);
    let target = EARTH_AREA / square_metres;
    // The loop is `level < max_level`, so that 25 and 43 are answers and never counted.
    (0..max_level)
        .find(|&level| count_zones(aperture, level) as f64 >= target)
        .unwrap_or(max_level)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math;
    use crate::{Error, get_grid, igeo7, isea3h, rtea3h};

    /// `geodeticArea` of `wholeWorld` (`GeoExtent.ec:108-128`, `0x429e0`), rebuilt from the six
    /// constants gcc folded into it, in the five operations of its tail (`0x42c1b` to `0x42c48`),
    /// and from the eC's own order.
    #[test]
    fn earth_area_is_the_compiled_double_and_not_the_source_order() {
        let asinh_term = f64::from_bits(0x3fb4_fe1d_b8b2_8e4d);
        let sin_term = 6378137.0;
        let (k1, k2, k3) = (
            f64::from_bits(0x442b_d973_89aa_9363),
            f64::from_bits(0x4288_22f2_51dc_fbf3),
            f64::from_bits(0x3eb0_20f2_a503_db0b),
        );
        let pi = std::f64::consts::PI;
        let x = (asinh_term - -asinh_term) * k1; // one
        let y = (sin_term - -sin_term) * k2; // two
        let sum = x + y; // three
        let compiled = k3 * sum * (pi - -pi); // four and five
        assert_eq!(compiled.to_bits(), 0x42fc_fe6e_831e_bb8c);
        assert_eq!(EARTH_AREA.to_bits(), 0x42fc_fe6e_831e_bb8c);

        // The eC's order, at the two poles, where the pole branch makes `sin` exactly +-1 and
        // `asinh` the folded term, and no `atan` or `tan` is reached.
        let major = 6378137.0f64;
        let minor = major - (major / 298.257223563);
        let (ooa, a2) = (1.0 / minor, minor * minor);
        let bmabpa = (major - minor) * (major + minor);
        let sr = math::sqrt(bmabpa);
        let (north, south) = (1.0f64, -1.0f64);
        let source = (pi - -pi)
            * (major
                * (a2 * minor * (asinh_term - -asinh_term)
                    + minor.abs()
                        * sr
                        * (north * math::sqrt(bmabpa * north * north + a2)
                            - south * math::sqrt(bmabpa * south * south + a2))))
            * 0.5
            * ooa
            / sr;
        assert_eq!(source.to_bits(), 0x42fc_fe6e_831e_bb8b);
        assert_ne!(source, EARTH_AREA);
    }

    #[test]
    fn the_folded_constants_are_pinned_as_bits() {
        assert_eq!(FIVE_SIXTHS.to_bits(), 0x3fea_aaaa_aaaa_aaab);
        assert_eq!(FIVE_SIXTHS, 5.0 / 6.0);
        assert_eq!(LN_65536.to_bits(), 0x4026_2e42_fefa_39ef);
        assert_eq!(LN_65536, 16.0 * std::f64::consts::LN_2);
    }

    /// The facts' maximum level is the indexing's own maximum resolution.
    #[test]
    fn the_maximum_level_is_the_indexings_maximum_resolution() {
        use crate::interfaces::Indexing;
        assert_eq!(max_level(7), crate::indexings::Z7::MAX_RESOLUTION);
        assert_eq!(max_level(3), crate::indexings::I3h::MAX_RESOLUTION);
    }

    /// The areas read from the engine, as bits, at four levels of each topology.
    #[test]
    fn zone_areas_are_the_engines_to_the_bit() {
        // (aperture, level, hexagon, pentagon)
        for (aperture, level, hexagon, pentagon) in [
            (7, 1, 0x429a_8239_276c_8e37, 0x4296_172f_a0da_7683),
            (7, 7, 0x418d_887c_eec2_05e0, 0x4188_9c68_1c4c_5a3b),
            (7, 13, 0x4080_738b_b47f_9a45, 0x407b_6b3e_2cd4_abc9),
            (7, 19, 0x3f72_540e_63bf_b274, 0x3f6e_8c17_fb94_d417),
            (3, 1, 0x42ae_ed42_adfe_a5eb, 0x42a9_c5b7_90fe_df99),
            (3, 10, 0x41c9_be46_9099_1359, 0x41c5_73e5_787f_9020),
            (3, 22, 0x4099_6592_adac_3e23, 0x4095_29fa_3b64_de73),
            (3, 33, 0x3f82_ca8c_52e6_248a, 0x3f7f_5194_8a2a_3ce6),
        ] {
            assert_eq!(
                zone_area(aperture, level, false).to_bits(),
                hexagon,
                "hexagon, aperture {aperture}, level {level}"
            );
            assert_eq!(
                zone_area(aperture, level, true).to_bits(),
                pentagon,
                "pentagon, aperture {aperture}, level {level}"
            );
        }
        // Level 0 has twelve pentagons and no hexagon.
        for aperture in [7, 3] {
            assert_eq!(
                zone_area(aperture, 0, true).to_bits(),
                0x42c3_5449_acbf_27b3
            );
        }
    }

    #[test]
    fn the_grid_facts_are_the_engines() {
        // (aperture, max level, max parents, max children, max neighbours, 64K depth, max depth)
        for (aperture, max, parents, children, neighbours, d64, depth) in
            [(7, 19, 2, 13, 6, 6, 9), (3, 33, 3, 7, 6, 10, 17)]
        {
            assert_eq!(max_level(aperture), max);
            assert_eq!(max_parents(aperture), parents);
            assert_eq!(max_children(aperture), children);
            assert_eq!(max_neighbors(aperture), neighbours);
            assert_eq!(depth_64k(aperture), d64);
            assert_eq!(max_depth(aperture), depth);
        }
        assert_eq!(count_zones(7, 19), 113_988_951_853_731_432);
        assert_eq!(count_zones(3, 33), 55_590_605_665_555_232);
        assert_eq!(ref_zone_area(7, 0), 4.25054684770074e13);
        assert_eq!(ref_zone_area(3, 0), 4.25054684770074e13);
        assert_eq!(ref_zone_area(7, 19), 4.474693498178629e-3);
        assert_eq!(ref_zone_area(3, 33), 9.175392417789991e-3);
    }

    /// The hexagon that `getMaxDepth` (`dggrs.ec:230-241`) counts the sub-zones of is the one
    /// at (0, 10) at level 2, which the eC chose to avoid ISEA3H's pentagon at (0, 0): it must be
    /// a hexagon on every grid, for the count in [`max_depth`] to be a hexagon's.
    #[test]
    fn the_test_zone_of_max_depth_is_a_hexagon_on_every_grid() {
        for name in crate::GRID_NAMES {
            let grid = get_grid(name).unwrap();
            let id = grid.zone_from_geo(0.0, 10.0, 2).unwrap();
            assert!(!grid.is_pentagon(id), "{name}");
        }
    }

    /// The counts of zones past the end of the engine's tables. On aperture 3 the engine
    /// computes `exp(level * ln 3) + 0.1` there, and converts it to `uint64`: the rounding of
    /// `exp` makes the count at 34 differ from `10 * 3^34 + 2`, the counts at 39 and 40 wrap, and
    /// from 41 the double is out of the range of `uint64` and converts to 0. They are pinned as
    /// integers, and the oracle test below checks each against the engine.
    #[test]
    fn counts_past_the_tables_are_pinned() {
        assert_eq!(A3_COUNTS_34_TO_42[0], 166_771_816_996_666_582);
        assert_eq!(count_zones(3, 34), 166_771_816_996_666_582);
        assert_ne!(count_zones(3, 34), 10 * 3u64.pow(34) + 2);
        assert_eq!(&A3_COUNTS_34_TO_42[7..], &[2, 2]);
        // Aperture 7: the last two entries of `powersOf7` are 7^23 and 7^24 modulo 2^64, and the
        // count from level 25 is 2.
        assert_eq!(POWERS_OF_7[23], 8_922_003_266_371_364_727);
        assert_eq!(POWERS_OF_7[24], 7_113_790_643_470_898_241);
        for n in 0..25u32 {
            assert_eq!(POWERS_OF_7[n as usize], 7u64.wrapping_pow(n), "7^{n}");
        }
        assert_eq!(
            count_zones(7, 23),
            8_922_003_266_371_364_727u64.wrapping_mul(10) + 2
        );
        assert_eq!((count_zones(7, 25), count_zones(7, 200)), (2, 2));
        assert_eq!((count_zones(3, 43), count_zones(3, 255)), (2, 2));
    }

    #[test]
    fn the_levels_from_an_area_at_the_special_values() {
        for (aperture, top) in [(7, 25), (3, 43)] {
            assert_eq!(level_from_ref_zone_area(aperture, 0.0), top);
            assert_eq!(level_from_ref_zone_area(aperture, -0.0), 0);
            assert_eq!(level_from_ref_zone_area(aperture, f64::NAN), top);
            assert_eq!(level_from_ref_zone_area(aperture, -1.0), 0);
            assert_eq!(level_from_ref_zone_area(aperture, f64::INFINITY), 0);
            assert_eq!(level_from_ref_zone_area(aperture, f64::NEG_INFINITY), 0);
            assert_eq!(level_from_ref_zone_area(aperture, f64::MAX), 0);
        }
    }

    /// The engine is not the inverse of itself: `earthArea / area` for the reference area of
    /// level 9 on ISEA3H rounds above the count of zones at level 9, so the answer is 10.
    #[test]
    fn the_reference_area_of_level_9_on_isea3h_answers_level_10() {
        assert_eq!(level_from_ref_zone_area(3, ref_zone_area(3, 9)), 10);
        let g = isea3h();
        assert_eq!(g.level_from_ref_zone_area(g.ref_zone_area(9).unwrap()), 10);
    }

    #[test]
    fn the_public_methods_agree_with_the_arithmetic() {
        let g = igeo7();
        let pentagon = g.zone_from_text("00").unwrap();
        let hexagon = g.zone_from_text("0064156").unwrap();
        assert_eq!(pentagon.area().unwrap().to_bits(), 0x42c3_5449_acbf_27b3);
        assert_eq!(hexagon.area().unwrap().to_bits(), 0x41e6_9c7f_a6cc_8c7f);
        assert_eq!(g.area(hexagon.id()), hexagon.area());
        assert_eq!(
            get_grid("igeo7").unwrap().area(hexagon.id()),
            hexagon.area()
        );

        assert_eq!(g.count_zones(5).unwrap(), 168_072);
        assert_eq!(g.ref_zone_area(0).unwrap(), 4.25054684770074e13);
        assert_eq!(
            g.meters_per_sub_zone(3, 2).unwrap(),
            math::sqrt(g.ref_zone_area(5).unwrap())
        );
        assert_eq!(g.level_from_meters_per_sub_zone(1000.0, 0), 10);
        assert_eq!(g.level_from_meters_per_sub_zone(1000.0, 3), 7);
        assert_eq!(g.level_from_meters_per_sub_zone(1e9, 40), 0);

        let any = get_grid("RTEA3H").unwrap();
        let rtea = rtea3h();
        assert_eq!(any.count_zones(33).unwrap(), rtea.count_zones(33).unwrap());
        assert_eq!(
            (
                any.refinement_ratio(),
                any.max_parents(),
                any.max_children()
            ),
            (3, 3, 7)
        );
        assert_eq!(
            (any.max_neighbors(), any.depth_64k(), any.max_depth()),
            (6, 10, 17)
        );
    }

    /// A zone without geometry has no area, by the crate's rule for a zone without geometry: the engine answers `+inf` for its null
    /// zone and, for a Z7 identifier of twenty digits, a finite area for a zone it cannot draw.
    #[test]
    fn a_zone_without_geometry_has_no_area() {
        use crate::{Address, Indexing, ZoneId};
        let level_20 = crate::indexings::Z7::encode(&Address::new(1, &[1; 20]));
        for name in crate::GRID_NAMES {
            let g = get_grid(name).unwrap();
            assert_eq!(g.area(ZoneId::NULL), None, "{name}");
            // Base cells 12 to 15 name no cell.
            assert_eq!(g.area(ZoneId(0xCFFF_FFFF_FFFF_FFFF)), None, "{name}");
        }
        for name in ["IGEO7", "IVEA7H", "RTEA7H"] {
            assert_eq!(get_grid(name).unwrap().area(level_20), None, "{name}");
        }
        assert_eq!(igeo7().zone(level_20).area(), None);
        assert_eq!(isea3h().zone(ZoneId::NULL).area(), None);
    }

    #[test]
    fn a_level_beyond_the_grids_maximum_is_refused() {
        for (name, max) in [
            ("IGEO7", 19u8),
            ("IVEA7H", 19),
            ("ISEA3H", 33),
            ("RTEA3H", 33),
        ] {
            let g = get_grid(name).unwrap();
            let refused = Error::ResolutionOutOfRange { res: max + 1, max };
            assert_eq!(g.count_zones(max + 1).unwrap_err(), refused);
            assert_eq!(g.ref_zone_area(max + 1).unwrap_err(), refused);
            assert_eq!(g.meters_per_sub_zone(max + 1, 0).unwrap_err(), refused);
            assert_eq!(g.meters_per_sub_zone(max, 1).unwrap_err(), refused);
            assert_eq!(g.meters_per_sub_zone(0, max + 1).unwrap_err(), refused);
            assert!(g.meters_per_sub_zone(max, 0).is_ok());
            assert!(g.meters_per_sub_zone(0, max).is_ok());
            // The sum is saturated at 255 in the error, as a `u8` cannot say more.
            assert_eq!(
                g.meters_per_sub_zone(255, 255).unwrap_err(),
                Error::ResolutionOutOfRange { res: 255, max }
            );
        }
    }

    #[cfg(feature = "oracle")]
    #[test]
    fn the_literals_of_the_public_tests_are_the_engines() {
        let o = dggal_oracle::IGEO7;
        assert_eq!(
            dggal_oracle::level_from_meters_per_sub_zone(o, 1000.0, 0),
            10
        );
        assert_eq!(
            dggal_oracle::level_from_meters_per_sub_zone(o, 1000.0, 3),
            7
        );
        assert_eq!(
            dggal_oracle::level_from_meters_per_sub_zone(o, 100.0, 0),
            12
        );
        assert_eq!(dggal_oracle::level_from_meters_per_sub_zone(o, 1e9, 40), 0);
        assert_eq!(dggal_oracle::count_zones(o, 5), 168_072);
        assert_eq!(
            dggal_oracle::level_from_ref_zone_area(
                dggal_oracle::ISEA3H,
                dggal_oracle::ref_zone_area(dggal_oracle::ISEA3H, 9)
            ),
            10
        );
    }

    #[cfg(feature = "oracle")]
    #[test]
    fn the_pinned_counts_and_every_level_are_the_engines() {
        for (name, oracle, aperture) in [
            ("ISEA3H", dggal_oracle::ISEA3H, 3u8),
            ("IGEO7", dggal_oracle::IGEO7, 7),
        ] {
            for level in 0..=255u8 {
                assert_eq!(
                    count_zones(aperture, level),
                    dggal_oracle::count_zones(oracle, i32::from(level)),
                    "countZones({level}) on {name}"
                );
            }
        }
    }
}
