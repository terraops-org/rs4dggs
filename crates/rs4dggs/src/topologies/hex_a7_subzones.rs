//! The sub-zone order of the aperture-7 topology: how many sub-zones a zone has at a relative
//! depth, where their order begins, and the order itself.
//!
//! The engine generates the sub-zones of a zone as centroids, row by row across the zone
//! ("scanlines"), beginning beside the vertex of the zone that stands uppermost on the
//! icosahedral net, and quantises each at the sub-zone level. This module holds that generator:
//! the net, the uppermost vertex, the centroid of the first sub-zone with the direction of the
//! scanlines turned as the engine turns it, the closed form of the count, the scanlines
//! themselves, the search of the entry at an index, the walk that finds the index of a
//! sub-zone, with the test of whether a zone is a sub-zone at all, and the quantisation of a
//! centroid; and, at its foot, the five answers that the topology gives through them: the
//! count, the first sub-zone, the order, the entry at an index, and the index of an entry.
//!
//! Provenance: of `dggrs/RI7H.ec` of DGGAL v0.0.6, `getSubZonesCount` (`:3063-3116`),
//! `getTopIcoVertex` (`:3126-3180`), `getFirstSubZoneCentroid` (`:3182-3258`),
//! `iterateI7HSubZones` (`:3260-3705`), `getVerticesDirections` (`:1897-1931`),
//! `zoneHasSubZone` (`:258-308`), `findSubZone` (`:212-222`) and `getSubZoneIndex`
//! (`:224-239`); of `projections/ri5x6.ec`, `toIcosahedronNet` (`:639-647`). None of them carries
//! `__attribute__ ((optimize("-fno-unsafe-math-optimizations")))`, and DGGAL is built with
//! `-O2 -ffast-math` (`Makefile.dggal:133`), so wherever `libdggal.so`, BuildID
//! `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, performs another order of operations than the eC
//! writes, the port performs the library's, and the comments marked `faithful:` record both.
//!
//! One site moves an answer: FS-1, the offset from the uppermost vertex to the first centroid,
//! where gcc replaced each division by three times a power of seven by a folded numerator
//! divided by the power. The others give the source's own doubles (TI-1, TI-2, IT-2), the
//! source's own test (TI-3, FS-2), the source's own integers (IT-1) or, at every direction that
//! the generator holds, the source's own numbers (IT-3, IT-4, IT-5), and are written as compiled
//! all the same. So are the three of the walk that finds an index: the directions to the
//! vertices (VD-1) and the step along each (ZH-1) give the source's own doubles, and the
//! comparison of two centroids (FZ-1) the source's own answer for every pair of numbers.
//!
//! One departure from the engine, in the search of the entry at an index: see
//! [`iterate_sub_zones`].
use super::hex_a7::{
    EVEN_HEX_A, EVEN_HEX_B, HexA7, NULL_GEOMETRY_LEVEL, ODD_HEX_A, ODD_HEX_B, ODD_HEX_C, Z,
    fold_centroid, from_centroid, get_vertices, pow7, rotate5x6_offset, zone_centroid,
    zone_from_steps, zone_level, zone_npoints,
};
use crate::fivebysix::{canonicalize5x6, move5x6};
use crate::{Address, PlanarPoint, Topology};

/// `triWidthOver2` (`ri5x6.ec:84`), half the width of a triangle of the icosahedral net: the
/// `.rodata` word at `0x49c90`.
const TRI_WIDTH_OVER_2: f64 = f64::from_bits(0x414d_4696_7cb2_2755);
/// faithful (TI-1): `triWidthOver2 * sqrt3` (`ri5x6.ec:644`), which gcc folded into the `.rodata`
/// word at `0x49c98`, the second half of the pair that `toIcosahedronNet` multiplies by at
/// `0x360b5`. Inert: the product formed at run time is the same double.
const TRI_WIDTH_OVER_2_SQRT3: f64 = f64::from_bits(0x4159_5a80_e2c9_cd20);
/// faithful (TI-2): `5*triWidthOver2` (`RI7H.ec:3159-3160`), which gcc folded into the `.rodata`
/// word at `0x4a1d0`, compared at `0x25ea4` and `0x25f3e`. Inert: the product formed at run time
/// is the same double.
const FIVE_TRI_WIDTH_OVER_2: f64 = f64::from_bits(0x4172_4c1e_0def_5895);

/// The six numerators of site FS-1, as `libdggal.so` holds them in `.rodata`: each is
/// `fl(k * fl(1/3))` for the `k` of its name, which gcc folded where the eC writes
/// `k / (3.0 * szp)` (`RI7H.ec:3199-3219`). For `k = -5` and `k = -10` that is one unit in the
/// last place away from the double nearest to `k/3`, so none of the six is written as a quotient.
const FS1_NUMERATOR_2: f64 = f64::from_bits(0x3fe5_5555_5555_5555); // 0x4a1e8
const FS1_NUMERATOR_1: f64 = f64::from_bits(0x3fd5_5555_5555_5555); // 0x49ba0
const FS1_NUMERATOR_MINUS_2: f64 = f64::from_bits(0xbfe5_5555_5555_5555); // 0x4a1d8
const FS1_NUMERATOR_MINUS_4: f64 = f64::from_bits(0xbff5_5555_5555_5555); // 0x4a1f0
const FS1_NUMERATOR_MINUS_5: f64 = f64::from_bits(0xbffa_aaaa_aaaa_aaaa); // 0x4a1f8
const FS1_NUMERATOR_MINUS_10: f64 = f64::from_bits(0xc00a_aaaa_aaaa_aaaa); // 0x4a1e0

/// The part of the way from a zone's centroid to one of its vertices at which
/// `zoneHasSubZone` takes its point, the `0.99` of `RI7H.ec:297`: the `.rodata` word at
/// `0x4a120`.
const VERTEX_STEP_FACTOR: f64 = 0.99;

/// A point of the 5x6 layout on the icosahedral net, in metres:
/// `RI5x6Projection::toIcosahedronNet` (`ri5x6.ec:639-647`), compiled at `0x36090`.
fn to_icosahedron_net(v: (f64, f64)) -> (f64, f64) {
    // faithful (TI-1): the sum and the difference of the two coordinates are formed as one packed
    // pair (`0x360a9`, `0x360ad`) and multiplied at `0x360b5` by the pair of `.rodata` words at
    // `0x49c90`: `triWidthOver2` and the folded `triWidthOver2 * sqrt3`. Inert: the eC groups
    // the second product in the same way, and the first is commutative.
    (
        (v.0 + v.1) * TRI_WIDTH_OVER_2,
        (v.0 - v.1) * TRI_WIDTH_OVER_2_SQRT3,
    )
}

/// The vertex of `z` that stands uppermost on the icosahedral net, and how many vertices share
/// that height: 2 where two do, or where the special case of an odd level applies, and 1
/// otherwise. Two vertices share a height where their ordinates on the net lie within 1e-6 of
/// each other; of the two, the one further left by more than 1e-8 is taken, unless it lies five
/// half-widths of a triangle or more away, across the cut of the net. Ports
/// `I7HZone::getTopIcoVertex` (`RI7H.ec:3126-3180`), compiled at `0x25d20` to `0x26200` in the
/// eC's own order: each comparison is taken as written, the two constants TI-1 and TI-2 are
/// folded, and the test TI-3 is rewritten.
///
/// At the pentagon of the north pole (root 0xA) the net is read upside down: the vertex of the
/// least ordinate is the uppermost, and of two at one height the one further right is taken,
/// where the two abscissae sum to more than five half-widths.
fn get_top_ico_vertex(z: &Z) -> ((f64, f64), i32) {
    // The centroid property, unfolded, as `c = this.centroid` (RI7H.ec:3129); `getVertices`
    // folds its own copy (RI7H.ec:1806-1809).
    let c = zone_centroid(z);
    let (fx, fy) = fold_centroid(c.0, c.1);
    let n_points = zone_npoints(z);
    let level = zone_level(z);
    // Five or six vertices, never none.
    let vertices = get_vertices(z.root, z.sub_hex, n_points, fx, fy, level);
    let n = vertices.len();
    // faithful (TI-3): `c.x - c.y - 1E-11 > 0` (RI7H.ec:3134) is compiled as
    // `c.x - c.y > 1e-11`: the difference is taken at `0x25d8f` and compared at `0x25dd9`, and
    // at `0x261bd` for a zone whose sub-hex is 2 or more, with the `.rodata` word at `0x49ec0`.
    // Inert: the same test in IEEE arithmetic, since a difference of two doubles is zero only
    // when they are equal.
    let north = c.0 - c.1 > 1e-11 || (n_points == 5 && z.root & 1 == 0);
    let odd_level = level & 1 != 0;
    let north_pole = n_points == 5 && z.root == 0xA;

    let mut top = 0;
    let mut equal_top = false;
    let mut top_ico = to_icosahedron_net(vertices[0]);
    for (i, &vertex) in vertices.iter().enumerate().skip(1) {
        let ico = to_icosahedron_net(vertex);
        let higher = if north_pole {
            ico.1 < top_ico.1 - 1e-6
        } else {
            ico.1 > top_ico.1 + 1e-6
        };
        let level_with = if north_pole {
            ico.1 <= top_ico.1 + 1e-6
        } else {
            ico.1 >= top_ico.1 - 1e-6
        };
        if higher {
            equal_top = false;
            top_ico = ico;
            top = i;
        } else if level_with {
            if !odd_level || north {
                equal_top = true;
            }
            // faithful (TI-2): `5*triWidthOver2` (RI7H.ec:3159-3160) is the folded word at
            // `0x4a1d0`, compared at `0x25ea4` with the sum and at `0x25f3e` with the
            // difference. Inert: the same double.
            let further = if north_pole {
                ico.0 > top_ico.0 + 1e-8 && top_ico.0 + ico.0 > FIVE_TRI_WIDTH_OVER_2
            } else {
                ico.0 < top_ico.0 - 1e-8 && top_ico.0 - ico.0 < FIVE_TRI_WIDTH_OVER_2
            };
            if further {
                top_ico = ico;
                top = i;
            }
        }
    }

    // The eC's remark: the first vertex cannot be to the right of the northern interruption
    // (RI7H.ec:3168-3175). The two floors are expanded in line (`0x25fc1` to `0x25fe8` and
    // `0x26012` to `0x26039`) as a truncation, less one where it exceeds its operand: that is
    // `f64::floor` save in the sign of a zero, which the comparison at `0x2603d` does not read.
    let next = (top + 1) % n;
    let special_odd_case = odd_level
        && north
        && (vertices[top].1 + 1e-11).floor() > (vertices[next].1 + 1e-11).floor()
        && vertices[top].1 - vertices[next].1 < 3.0;
    if special_odd_case {
        top = next;
    }
    (
        vertices[top],
        if equal_top || special_odd_case { 2 } else { 1 },
    )
}

/// The centroid of the first sub-zone `r_depth` levels below `z`: one of four offsets from the
/// uppermost vertex, chosen by the parities of the depth and of the zone's level, turned by sixty
/// degrees where the zone has two uppermost vertices astride an interruption or is a pentagon,
/// and then walked across the layout by [`move5x6`]. Ports `I7HZone::getFirstSubZoneCentroid`
/// (`RI7H.ec:3182-3258`), compiled at `0x26200` to `0x265a0`.
///
/// `s` is the eC's pair `sx` and `sy`, which it heeds only where neither is null: the direction
/// of the scanlines, which the generator passes so that it is turned in place as the offset is
/// turned, and which `getFirstSubZone` does not pass. It is turned before the walk, so it is the
/// engine's whether or not the walk ends.
///
/// `None` is the answer of [`move5x6`] where its walk does not end within its bound, and the
/// engine does not return; no call for a valid zone is known to reach it.
///
/// `r_depth` is not negative, and the sub-zone level, `r_depth` below the zone's, is at most 19
/// for every caller: there the power of seven that scales the offset is one of the eleven that
/// `pow7` holds as integers, as the engine reads it from its own table (`0x2625d`). The four
/// answers at the foot of this module, through which alone the topology reaches this function,
/// refuse a deeper request before they ask.
///
/// At the pentagon of the south pole at an odd level and depth 2 the eC turns `*sx` and `*sy`
/// five times without asking whether it was given them (RI7H.ec:3246-3255; the loads at
/// `0x26568` and `0x26576`), and `getFirstSubZone`, which gives none, ends the process there.
/// Here a caller that gives no direction has none turned, and the centroid is the one that the
/// eC's arithmetic gives, which is the one the engine answers when it is given a direction.
pub(super) fn get_first_sub_zone_centroid(
    z: &Z,
    r_depth: i64,
    mut s: Option<&mut (f64, f64)>,
) -> Option<(f64, f64)> {
    let (v, n_top) = get_top_ico_vertex(z);
    let level = zone_level(z);
    let odd_level = level & 1 != 0;
    let odd_depth = r_depth & 1 != 0;
    let szp = pow7((level + 1 + r_depth) / 2) as f64;
    let n_points = zone_npoints(z);
    let north_pole = z.root == 0xA && n_points == 5;
    let pgon_tweak = n_points == 5 && !odd_depth && odd_level;
    let north = z.root & 1 == 0;

    // faithful (FS-1): the eC writes each offset as `k / (3.0 * szp)` (RI7H.ec:3199-3219). In
    // `libdggal.so` gcc folded `k / 3.0` into the numerator `fl(k * fl(1/3))` and divides it by
    // `(double)szp`, a `cvtsi2sd` of the 64-bit integer and two `divsd`: at `0x26500` to
    // `0x26514` for an odd depth below an odd level, `0x262c3` to `0x262d7` for an odd depth
    // below an even level, `0x26385` to `0x263a7` and, for a pentagon, `0x26531` to `0x26552`
    // for an even depth below an odd level, and `0x264a0` to `0x264bd` for an even depth below
    // an even level. Active: the library divides once where the eC multiplies and then
    // divides, and for -5 and -10 its numerator is not the double nearest to `k/3`, so the
    // offset differs by a unit in the last place at some zones, and the first centroid with it.
    // The first sub-zone centroid is one sub-zone edge away from the vertex at an odd depth, and
    // two sub-zone edges towards the next vertex at an even one (the eC's remarks).
    let (mut dx, mut dy) = match (odd_depth, odd_level) {
        (true, true) => (FS1_NUMERATOR_2 / szp, FS1_NUMERATOR_1 / szp),
        (true, false) => (FS1_NUMERATOR_MINUS_4 / szp, FS1_NUMERATOR_MINUS_5 / szp),
        (false, true) => (FS1_NUMERATOR_MINUS_10 / szp, FS1_NUMERATOR_MINUS_2 / szp),
        (false, false) => (FS1_NUMERATOR_MINUS_4 / szp, FS1_NUMERATOR_MINUS_2 / szp),
    };

    // faithful (FS-2): `v.x - v.y - 1E-11 > 0` (RI7H.ec:3223) is compiled as
    // `v.x - v.y > 1e-11`, the difference at `0x264d0` and the comparison at `0x264d6`. Inert:
    // the same test in IEEE arithmetic.
    if (pgon_tweak && north) || (!pgon_tweak && n_top == 2 && (odd_level || v.0 - v.1 > 1e-11)) {
        // A hexagon astride the interruption between its two uppermost vertices: the offset is
        // turned by sixty degrees anticlockwise, three times at the north pole.
        let n_rotation = if north_pole { 3 } else { 1 };
        for i in 0..n_rotation {
            (dx, dy) = (dy, dy - dx);
            // The north pole at depth 2 skips the last turn of the direction (the eC's remark,
            // RI7H.ec:3236-3237).
            let skipped = north_pole && pgon_tweak && r_depth == 2 && i == n_rotation - 1;
            if let Some(s) = s.as_deref_mut().filter(|_| !skipped) {
                *s = (s.1, s.1 - s.0);
            }
        }
    } else if pgon_tweak && z.root == 0xB && r_depth == 2 {
        // The pentagon of the south pole: five turns of the direction, and none of the offset.
        // The eC makes them through its two pointers whether or not they are null.
        if let Some(s) = s {
            for _ in 0..5 {
                *s = (s.1, s.1 - s.0);
            }
        }
    }
    move5x6(v, dx, dy, 1, None, false)
}

/// The number of sub-zones `r_depth` levels below a zone of `n_points` sides, and 1 for a depth
/// of 0 or less: the engine's closed form, `getSubZonesCount` (`RI7H.ec:3063-3116`), which
/// `countSubZones` of the Z7 class repeats (`RI7H_Z7.ec:492-500`). A hexagon has
/// `7^d + 5 * 7^((d-1)/2) + 1` at an odd depth `d` and `7^d + 7^(d/2) - 1` at an even one; a
/// zone of `n` sides has that number times `n`, plus 5, over 6.
///
/// The powers are integers: `pow7` holds eleven and goes through a double beyond them, and
/// `7^19`, which a pentagon of level 0 asks at depth 19, is no double. The arithmetic is the
/// engine's own, 64-bit and signed (the product, the sum and the division at `0x1e311` to
/// `0x1e331`), and from depth 22 it wraps as the engine's does, whose table `powersOf7`
/// (`RI7H.ec:837-841`) holds `7^0` to `7^24` modulo `2^64`. Beyond depth 24 the engine takes its
/// powers through `exp`, which is not followed: the function answers, and its answer there is
/// its own. No zone within level 19 asks beyond depth 19, and no public call reaches a depth
/// beyond it: [`count`], through which the topology answers, refuses a request whose sub-zones
/// would lie beyond level 19 before it asks, as the grid does before the topology.
pub(super) fn sub_zones_count(n_points: i64, r_depth: i64) -> u64 {
    if r_depth <= 0 {
        return 1;
    }
    let pow = |n: i64| 7_i64.wrapping_pow(n as u32);
    let n_hex_sub_zones = if r_depth & 1 != 0 {
        pow(r_depth)
            .wrapping_add(pow((r_depth - 1) / 2).wrapping_mul(5))
            .wrapping_add(1)
    } else {
        pow(r_depth).wrapping_add(pow(r_depth / 2)).wrapping_sub(1)
    };
    (n_hex_sub_zones.wrapping_mul(n_points).wrapping_add(5) / 6) as u64
}

/// The callback through which the generator delivers each sub-zone centroid with its index, and
/// which answers whether to go on, as the eC's `centroidCallback`. The centroid is `None` where
/// a walk of [`move5x6`] on the way to it did not end within its bound.
type Callback<'a> = &'a mut dyn FnMut(i64, Option<(f64, f64)>) -> bool;

/// At an odd depth, the numbers of scanlines of a cap and of the middle band of the order,
/// `nCapSL` and `nMidSL`, from the number of scanlines between two vertices, `nInterSL`, which
/// is a power of seven: the ceilings of a third of it and of two thirds of it.
fn odd_depth_scanline_counts(n_inter_sl: i64) -> (i64, i64) {
    // faithful (IT-1): the eC writes `(int64)(ceil(nInterSL / 3.0) + 0.5)` and
    // `(int64)(ceil(nInterSL * 2 / 3.0) + 0.5)` (RI7H.ec:3299-3300). In `libdggal.so` each
    // division is a product with the double nearest to a third, the `.rodata` word at
    // `0x49ba0` loaded at `0x26894`: the count is converted and multiplied at `0x268a4` and
    // `0x268b5`, and its double, formed at `0x268ee`, at `0x268ff` and `0x26909`; each ceiling
    // is expanded in line and truncated at `0x268f6` and `0x2694b`. Inert: for every power of
    // seven to `7^18` both are the integers below, and no depth within level 19 asks beyond
    // `7^9`; the two forms part at `7^19` alone, which is no double.
    ((n_inter_sl + 2) / 3, (2 * n_inter_sl + 2) / 3)
}

/// The eC's post-increment, `counter++`: the value before, and the counter one more.
///
/// No eC counterpart: the eC writes the operator in line, in the counters of the scanlines.
fn post_increment(counter: &mut i64) -> i64 {
    let before = *counter;
    *counter += 1;
    before
}

/// A direction turned by sixty degrees: the eC's `rotate5x6Offset` (`RI7H.ec:1933-1947`) with
/// its result assigned back, as the generator always uses it.
fn turned(t: (f64, f64), clockwise: bool) -> (f64, f64) {
    rotate5x6_offset(t.0, t.1, clockwise)
}

/// A direction turned by 120 degrees clockwise, as the library forms it in one step wherever
/// the eC turns a direction twice clockwise, or four times anticlockwise, past a polar pentagon.
fn turned_120_degrees_clockwise((x, y): (f64, f64)) -> (f64, f64) {
    // faithful (IT-3): two calls of `rotate5x6Offset` clockwise give `(x - y, x)` and then
    // `((x - y) - x, x - y)` (RI7H.ec:3464-3469 and :3475-3480 at an odd depth, :3626-3630 at an
    // even one). gcc folded the second abscissa into `-y`: the pair is `x - y` and the negation
    // of `y`, at `0x2798c`, `0x27958` and `0x27962` for the odd depths and at `0x27d3c` and
    // `0x27d40` for the even ones.
    // faithful (IT-4): four calls anticlockwise (RI7H.ec:3634-3646, and :3664-3671) give the
    // same direction, and gcc formed it in the same way: at `0x2856e` and `0x28576` at the
    // north pole, at `0x28440` and `0x27f00` at the south.
    // Both inert as numbers: on the six directions that the generator holds every difference
    // of the eC's turns is exact, and the pair is the same. A zero may differ in its sign (the
    // eC's `(x - y) - x` is `+0` where `y` is `+0`, and the library's `-y` is `-0`), so the
    // library's form is the one written.
    (-y, x - y)
}

/// A direction turned five times anticlockwise, as the library forms it in one step at the
/// pentagon of the south pole.
fn turned_five_times_anticlockwise((x, y): (f64, f64)) -> (f64, f64) {
    // faithful (IT-4): five calls of `rotate5x6Offset` anticlockwise (RI7H.ec:3658-3671) come
    // to one turn clockwise, `(x - y, x)`. gcc formed the abscissa as the negation of `y - x`,
    // at `0x27efc` and `0x27f00`, and left the ordinate as `x`. Inert as numbers, for the
    // reason given at IT-3.
    (-(y - x), x)
}

/// Generates the centroids of the sub-zones `r_depth` levels below `z` in the engine's order,
/// handing each to `callback` with its index until it answers `false`; answers the index at
/// which it was stopped, or -1 where the callback let it run to the end. Ports
/// `I7HZone::iterateI7HSubZones` (`RI7H.ec:3260-3705`), compiled at `0x26640` to `0x28a30`.
///
/// The order is made of scanlines across the zone, from the vertex that stands uppermost on the
/// icosahedral net ([`get_first_sub_zone_centroid`]): the integers decide how many centroids
/// each scanline has and where it starts, and every centroid is a point reached by
/// [`move5x6`], from the first centroid to the start of its scanline and from there along it.
/// Past a pentagon a scanline is bent at its middle, and its second half is reached by a
/// second walk in a direction turned for the pole or the hemisphere.
///
/// With `search_index` other than -1 only the scanline that holds that index is generated, and
/// of it one centroid, so the search costs the integers of every scanline before it and a
/// handful of walks: it grows with the number of scanlines, and never with the number of
/// centroids before the index.
///
/// One departure from the engine, in that search. The second half of a scanline past a
/// pentagon, at an odd depth, begins at an offset reckoned from the half-width of the
/// pentagon's own scanline, `pSLZones`, which the eC records inside the block that generates a
/// scanline (`RI7H.ec:3401-3402`). A search enters that block for one scanline alone, so for
/// any later scanline it reads the half-width as zero, and the engine's `getSubZoneAtIndex`
/// answers there another zone than its `getSubZones` lists at that index: at every pentagon, at
/// odd depths. Here the half-width is recorded on the way past, and the entry found by index is
/// the entry of the sequence.
///
/// A centroid is `None` where a walk of [`move5x6`] on the way to it did not end within its
/// bound: the first centroid, the start of the scanline, or the centroid's own walk. This is
/// the one place where the crate gives an answer and the engine would not return; no call for a
/// valid zone is known to reach it (of a thousand million walks of real orders the longest
/// took 4 turns of a bound of 64). The length of the order and the place of every other entry
/// stay the engine's, and the direction that an unfinished walk leaves behind is not read.
///
/// `r_depth` is not negative, and the sub-zone level, `r_depth` below the zone's, is at most 19
/// for every caller, as for [`get_first_sub_zone_centroid`]: there each power of seven asked is
/// one that `pow7` holds as an integer. [`sub_zones`] and [`sub_zone_at_index`] refuse a deeper
/// request before they ask.
///
/// It is compiled in the eC's own order save at five sites, each recorded where it stands: the
/// two integer ceilings (IT-1), the direction of the scanlines (IT-2), the turns past a polar
/// pentagon (IT-3, IT-4) and one sum past the pentagon of the south pole (IT-5). Every product
/// of a count with a direction is one conversion and one multiplication, and each difference
/// of two products is taken as written (`0x26ca8` to `0x26ccb`, `0x273ad` to `0x273de`,
/// `0x278b0` to `0x278d9`).
pub(super) fn iterate_sub_zones(
    z: &Z,
    r_depth: i64,
    search_index: i64,
    callback: Callback<'_>,
) -> i64 {
    if r_depth == 0 {
        return if callback(0, Some(zone_centroid(z))) {
            -1
        } else {
            0
        };
    }
    let mut keep_going = true;
    let level = zone_level(z);
    let sz_level = level + r_depth;
    let odd_ancestor = level & 1 != 0;
    let odd_depth = r_depth & 1 != 0;
    let odd_level_sz = sz_level & 1 != 0;
    let szp = pow7((sz_level + i64::from(odd_level_sz)) / 2);
    // Centroid to centroid distance between sub-zones along 5x6 x and y axes (the eC's remark).
    let c2c = 1.0 / szp as f64;
    let mut c_start: i64 = 0;
    let mut index: i64 = 0;
    let mut zones_per_sl: i64;
    // Direction along scanlines (the eC's remark).
    // faithful (IT-2): the eC writes `c2c * (oddLevelSZ ? 3 : 1)` and `c2c * (oddLevelSZ ? 2 :
    // 1)` (RI7H.ec:3281-3282). At an odd sub-zone level the library multiplies the first by 3.0
    // at `0x26729` and forms the second as the sum `c2c + c2c` at `0x26725`; at an even one it
    // multiplies neither. Inert: the same doubles.
    let mut direction = if odd_level_sz {
        (c2c * 3.0, c2c + c2c)
    } else {
        (c2c, c2c)
    };
    let n_points = zone_npoints(z);
    let root = z.root;
    let south = root & 1 != 0;

    let first = get_first_sub_zone_centroid(z, r_depth, Some(&mut direction));
    let (sx, sy) = direction;

    // Rotate scanline direction to get direction to next scanline (hexagon immediately to the
    // left: the eC's remarks).
    let (nsx, nsy) = if !odd_ancestor && odd_depth {
        (sx - sy, sx) // 60 degrees clockwise
    } else {
        (-sy, sx - sy) // 120 degrees clockwise
    };

    // The start of scanline `s`, `start` centroids along the direction of the scanlines from
    // the foot of the walk to the next scanline, with the direction as the walk leaves it: the
    // eC's `move5x6(sc, first, s * nsx - cStart * tsx, s * nsy - cStart * tsy, 1, &tsx, &tsy,
    // true)`. `None` where the first centroid is none or the walk does not end.
    let scanline_start = |s: i64, start: i64| -> Option<((f64, f64), (f64, f64))> {
        let mut ts = (sx, sy);
        let dx = s as f64 * nsx - start as f64 * sx;
        let dy = s as f64 * nsy - start as f64 * sy;
        let sc = move5x6(first?, dx, dy, 1, Some(&mut ts), true)?;
        Some((sc, ts))
    };

    if odd_depth {
        let n_inter_sl = pow7((r_depth - 1) / 2);
        let (n_cap_sl, n_mid_sl) = odd_depth_scanline_counts(n_inter_sl);
        let b = n_cap_sl;
        let c = b + n_inter_sl;
        let d = c + n_mid_sl;
        let e = d + n_inter_sl;
        let (mut left_ac_counter, mut left_ce_counter) = (0_i64, 0_i64);
        let (mut right_bd_counter, mut right_df_counter) = (0_i64, 0_i64);
        let mut n_until_pentagon = i64::MAX;
        let mut past_pentagon = false;
        let (mut p_sl_zones, mut p_right_bd_counter_mod5) = (0_i64, 0_i64);
        let big_s = 2 * (n_inter_sl - 1);

        let mut n_scanlines = 2 * n_cap_sl + 2 * n_inter_sl + n_mid_sl;
        if n_points == 5 {
            n_until_pentagon = n_inter_sl + n_mid_sl;
            n_scanlines -= (4 * n_inter_sl + 11) / 15;
        }
        // For South, number of scanlines on left side of interruption (the eC's remark).
        let n = n_scanlines - i64::from(r_depth > 1) - big_s / 5;

        zones_per_sl = 1;

        for s in 0..n_scanlines {
            if !keep_going {
                break;
            }
            if s == n_until_pentagon {
                past_pentagon = true;
            }

            let (mut left, mut right): (i64, i64);
            if s < b {
                left = i64::from(s > 0 && post_increment(&mut left_ac_counter) % 4 == 0);
                right = if s > 0 { 5 } else { 0 };
            } else if s < c {
                left = i64::from(post_increment(&mut left_ac_counter) % 4 == 0);
                right = if s == b {
                    2
                } else {
                    i64::from(post_increment(&mut right_bd_counter) % 5 != 0)
                };
            } else if n_points == 5 {
                if s < d {
                    left = if s == c || post_increment(&mut left_ce_counter) % 5 != 0 {
                        -1
                    } else {
                        0
                    };
                    right = i64::from(post_increment(&mut right_bd_counter) % 5 != 0);
                } else {
                    left = if post_increment(&mut left_ce_counter) % 5 != 0 {
                        -1
                    } else {
                        0
                    };
                    right = if s == d {
                        1
                    } else if post_increment(&mut right_df_counter) % 4 == 0 {
                        -1
                    } else {
                        0
                    };
                }

                if s >= c + b {
                    if odd_ancestor {
                        right -= 1;
                    } else {
                        left -= 1;
                    }
                }
                if s >= e {
                    left -= 1 + i64::from(s > e) * 3 + i64::from(left_ce_counter % 5 == 1);
                }

                if s == n_until_pentagon {
                    p_right_bd_counter_mod5 = (right_bd_counter + 4) % 5;
                }
            } else if s < d {
                left = if s == c || post_increment(&mut left_ce_counter) % 5 != 0 {
                    -1
                } else {
                    0
                };
                right = i64::from(post_increment(&mut right_bd_counter) % 5 != 0);
            } else if s < e {
                left = if post_increment(&mut left_ce_counter) % 5 != 0 {
                    -1
                } else {
                    0
                };
                right = if s == d {
                    1
                } else if post_increment(&mut right_df_counter) % 4 == 0 {
                    -1
                } else {
                    0
                };
            } else {
                left = if s == e { -2 } else { -5 };
                right = if post_increment(&mut right_df_counter) % 4 == 0 {
                    -1
                } else {
                    0
                };
            }

            c_start += if odd_ancestor { left } else { right };
            zones_per_sl += left + right;

            // The departure recorded above: the eC assigns `pSLZones` inside the block below
            // (RI7H.ec:3401-3402), which a search skips for every scanline but its own.
            if s == n_until_pentagon {
                p_sl_zones = zones_per_sl / 2;
            }

            if search_index == -1 || (search_index >= index && search_index < index + zones_per_sl)
            {
                let pgon_redir = past_pentagon && (!south || s > n_until_pentagon || root == 0xB);
                let mut h = zones_per_sl / 2;

                if pgon_redir && s > n_until_pentagon {
                    if odd_ancestor {
                        let m = s.min(e) - n_until_pentagon;
                        let shift = (5 - p_right_bd_counter_mod5) % 5;
                        let n5 = (m + shift - 1) / 5;

                        h = (p_sl_zones - (m - n5)).max(1);
                        if s >= e {
                            h -= (s - e) * 5 + 1;
                        }
                    } else {
                        let dms = s.min(d) - n_until_pentagon;

                        h = p_sl_zones;
                        if dms != 0 {
                            h -= (dms + p_right_bd_counter_mod5) / 5;
                        }
                        if s > d {
                            let smd = s - d;
                            h -= smd + (smd + 3) / 4;
                        }
                    }
                }

                let mut i = if search_index != -1 {
                    // `i = (int)(searchIndex - index)` (RI7H.ec:3435): a subtraction of 32
                    // bits, extended (`0x27a0a` to `0x27a12`). No scanline within level 19
                    // holds `2^31` centroids.
                    let i = i64::from((search_index - index) as i32);
                    index = search_index;
                    i
                } else {
                    0
                };

                // Start of scanline (the eC's remark).
                let start = scanline_start(s, c_start);

                while i < zones_per_sl {
                    let centroid = start.and_then(|(sc, ts)| {
                        if pgon_redir && i > h {
                            let cc = move5x6(sc, h as f64 * ts.0, h as f64 * ts.1, 1, None, true)?;
                            let tt = if root == 0xA && s == n_until_pentagon && !odd_ancestor {
                                turned(ts, false)
                            } else if root == 0xA || root == 0xB {
                                // faithful (IT-3): the second of the two turns clockwise is
                                // folded into the first; see `turned_120_degrees_clockwise`.
                                if s < n || odd_ancestor {
                                    turned_120_degrees_clockwise(ts)
                                } else {
                                    turned(ts, true)
                                }
                            } else if !south || (!odd_ancestor && s >= n) {
                                turned(ts, true)
                            } else {
                                ts
                            };
                            let k = (i - h) as f64;
                            move5x6(cc, k * tt.0, k * tt.1, 1, None, true)
                        } else {
                            move5x6(sc, i as f64 * ts.0, i as f64 * ts.1, 1, None, true)
                        }
                    });
                    keep_going = callback(index, centroid);
                    if search_index != -1 || !keep_going {
                        break;
                    }
                    index += 1;
                    i += 1;
                }
            } else {
                index += zones_per_sl;
            }
        }
    } else {
        // Even depths
        let p = pow7(r_depth / 2);
        let n_cap_sl = (p - 1) / 3;
        let n_mid_sl = (2 * p + 1) / 3;
        let b = n_cap_sl;
        let c = b + n_mid_sl;
        let mut n_until_pentagon = i64::MAX;
        let mut past_pentagon = false;
        let odd_south_pentagon = n_points == 5 && odd_ancestor && south;
        // n is the number of scanlines starting on the left side of the interruption; for more
        // complicated odd parent / south pentagons (the eC's remarks).
        let (mut n, mut flip_start) = (0_i64, 0_i64);

        let mut n_scanlines = 2 * n_cap_sl + n_mid_sl;
        if n_points == 5 {
            n_scanlines -= n_cap_sl / 2;
            n_until_pentagon = n_mid_sl - 1;
        }

        zones_per_sl = 3;

        if n_points == 5 && odd_ancestor {
            let big_s = 2 * (p - 1);
            n = n_until_pentagon + big_s / 5 + i64::from(big_s % 5 == 4);
        }

        for s in 0..n_scanlines {
            if !keep_going {
                break;
            }
            if s == n_until_pentagon {
                past_pentagon = true;
            }

            let (left, mut right): (i64, i64);
            if s < b {
                left = i64::from(s > 0);
                right = if s > 0 { 2 } else { 0 };
            } else if s < c {
                let odd = s == b || (s - b) & 1 != 0;
                left = if odd { 0 } else { -1 };
                right = i64::from(odd);
            } else {
                left = if s == c { -1 } else { -2 };
                right = if s == c { 0 } else { -1 };
            }

            if past_pentagon && s > n_until_pentagon {
                right -= 1;
            }

            c_start += left;
            zones_per_sl += left + right;

            if s == n {
                flip_start = c_start;
            }

            if search_index == -1 || (search_index >= index && search_index < index + zones_per_sl)
            {
                let pgon_redir = past_pentagon
                    && (!south || s > n_until_pentagon || (odd_ancestor && root == 0xB));

                // Start of scanline (the eC's remark).
                let start = if odd_south_pentagon && s > n {
                    // This scenario is rather complicated as rotations not based on
                    // interruption are needed (the eC's remark).
                    scanline_start(n, flip_start).and_then(|(sc, mut ts)| {
                        // This needs a special case for depth = 2 (the eC's remark).
                        let tsc = if n_scanlines == n + 2 {
                            sc
                        } else {
                            let step = ts;
                            move5x6(sc, step.0, step.1, 1, Some(&mut ts), true)?
                        };

                        ts = turned(ts, true);

                        let step = ts;
                        let mut sc = move5x6(tsc, step.0, step.1, 1, Some(&mut ts), true)?;
                        let before = ts.0;
                        ts = turned(ts, false);

                        if s > n + 1 {
                            // faithful (IT-5): the eC turns the direction clockwise into `dd`
                            // and adds, `tsx + dd.x` and `tsy + dd.y` (RI7H.ec:3592-3596),
                            // where `dd.x` is `tsx - tsy` of a direction that the line above
                            // has just turned anticlockwise. gcc took for `dd.x` the abscissa
                            // from before that turn, loaded at `0x274aa` and added at
                            // `0x283d4`; the other sum is as written, at `0x283d0`, and each
                            // is multiplied once, at `0x2840f` and `0x28413`. Inert: on the
                            // six directions that the generator holds a turn anticlockwise
                            // and one clockwise give the abscissa back exactly.
                            let k = (s - n - 1) as f64;
                            sc =
                                move5x6(sc, k * (ts.0 + before), k * (ts.1 + ts.0), 1, None, true)?;
                        }
                        Some((sc, ts))
                    })
                } else {
                    scanline_start(s, c_start)
                };

                let mut i = if search_index != -1 {
                    let i = search_index - index;
                    index = search_index;
                    i
                } else {
                    0
                };

                let h = zones_per_sl / 2;
                while i < zones_per_sl {
                    let centroid = start.and_then(|(sc, ts)| {
                        if pgon_redir && i > h {
                            let cc = move5x6(sc, h as f64 * ts.0, h as f64 * ts.1, 1, None, true)?;
                            let tt = if root == 0xA {
                                if !odd_ancestor {
                                    // faithful (IT-3): two turns clockwise, as one.
                                    turned_120_degrees_clockwise(ts)
                                } else {
                                    // faithful (IT-4): the first turn anticlockwise and the
                                    // three that follow it past the pentagon, as one; the
                                    // turn of a scanline beyond `n` is made as written, at
                                    // `0x281b2` to `0x281be`.
                                    let t = if s > n_until_pentagon {
                                        turned_120_degrees_clockwise(ts)
                                    } else {
                                        turned(ts, false)
                                    };
                                    if s > n { turned(t, false) } else { t }
                                }
                            } else if root == 0xB {
                                // faithful (IT-4): four turns anticlockwise, or five beyond
                                // `n` at an odd level, as one.
                                if odd_ancestor && s > n {
                                    turned_five_times_anticlockwise(ts)
                                } else {
                                    turned_120_degrees_clockwise(ts)
                                }
                            } else if south {
                                if odd_ancestor && s > n {
                                    turned(ts, true)
                                } else {
                                    ts
                                }
                            } else {
                                turned(ts, true)
                            };
                            let k = (i - h) as f64;
                            move5x6(cc, k * tt.0, k * tt.1, 1, None, true)
                        } else {
                            let n_rotations = if !odd_ancestor && past_pentagon { 2 } else { 1 };
                            move5x6(
                                sc,
                                i as f64 * ts.0,
                                i as f64 * ts.1,
                                n_rotations,
                                None,
                                true,
                            )
                        }
                    });
                    keep_going = callback(index, centroid);
                    if search_index != -1 || !keep_going {
                        break;
                    }
                    index += 1;
                    i += 1;
                }
            } else {
                index += zones_per_sl;
            }
        }
    }
    if keep_going { -1 } else { index }
}

/// The centroid at `index` of the order of [`iterate_sub_zones`], found by its search, which
/// generates one scanline and one centroid of it: the walk of `getSubZoneAtIndex`
/// (`RI7H.ec:241-256`) with its callback, less the fault recorded at [`iterate_sub_zones`], so
/// that the centroid found is the centroid of the sequence at every index. The cost is the
/// integers of the scanlines before the one that holds the index, at most some
/// `(10 / 3) * 7^(depth / 2)` of them, and a handful of walks.
///
/// `None` where the order holds no such index, and where a walk on the way to the centroid did
/// not end (see [`iterate_sub_zones`]); a caller that must tell the two apart compares the
/// index with [`sub_zones_count`] first.
pub(super) fn sub_zone_centroid_at_index(z: &Z, r_depth: i64, index: u64) -> Option<(f64, f64)> {
    // An index that is no `int64` is in no order, and -1 would ask for the whole of it.
    let search_index = i64::try_from(index).ok()?;
    let mut found = None;
    iterate_sub_zones(z, r_depth, search_index, &mut |at, centroid| {
        // At depth 0 the generator hands over the zone's own centroid whatever the index asked.
        if at == search_index {
            found = centroid;
        }
        false
    });
    found
}

/// The six directions from the centroid of `z` to its vertices, as offsets in the 5x6 layout:
/// `I7HZone::getVerticesDirections` (`RI7H.ec:1897-1931`), compiled at `0x1e990` to `0x1ea9e`.
/// They are the offsets of a hexagon of the zone's level, six of them for a pentagon too, as
/// the eC's own remark asks; none is carried across an interruption.
fn get_vertices_directions(z: &Z) -> [(f64, f64); 6] {
    // faithful (VD-1): the eC forms `oonp = 1.0 / (7 * p)` and writes each coordinate as
    // `- oonp * A` or `+ oonp * A` (RI7H.ec:1901-1928). The library negates the reciprocal
    // once, at `0x1e9ee`, multiplies it by each offset (`0x1ea00` and `0x1ea08` at an even
    // level; `0x1ea70`, `0x1ea78` and `0x1ea7c` at an odd one), and takes each coordinate that
    // the eC writes with a plus as the negation of such a product (`0x1ea14`, `0x1ea18`;
    // `0x1ea88` to `0x1ea94`). Inert: a negation is exact, so each is the eC's own double.
    let m = -(1.0 / (7 * pow7(z.l49r)) as f64);
    if z.sub_hex == 0 {
        // Even level (the eC's remark).
        let (a, b) = (m * EVEN_HEX_A, EVEN_HEX_B * m);
        [(a, b), (b, a), (a, -a), (-a, -b), (-b, -a), (-a, a)]
    } else {
        // Odd level (the eC's remark).
        let (a, b, c) = (m * ODD_HEX_A, ODD_HEX_B * m, m * ODD_HEX_C);
        [(a, b), (b, c), (c, -a), (-a, -b), (-b, -c), (-c, a)]
    }
}

/// Whether `needle` is a sub-zone of `hay_stack`, as the engine decides it:
/// `RhombicIcosahedral7H::zoneHasSubZone` (`RI7H.ec:258-308`), compiled at `0x1eaa0` to
/// `0x1ebd7`. A zone is a sub-zone of a coarser one where a point 99 hundredths of the way from
/// its centroid to one of its six vertices lies in the coarser zone: its centroid need not, so
/// that the zones across the boundary of the coarser one are its sub-zones too.
///
/// A walk of [`move5x6`] that does not end within its bound finds nothing in its direction,
/// as a point that no zone holds finds nothing.
fn zone_has_sub_zone(hay_stack: &Z, needle: &Z) -> bool {
    let z_level = zone_level(hay_stack);
    if zone_level(needle) <= z_level {
        return false;
    }
    // `getZoneCRSCentroid(needle, 0, c)` (RI7H.ec:267) is the centroid property, unfolded.
    let c = zone_centroid(needle);
    get_vertices_directions(needle).iter().any(|v| {
        // faithful (ZH-1): the eC writes `v[i].x * 0.99` and `v[i].y * 0.99` (RI7H.ec:297);
        // the library loads the `.rodata` word at `0x4a120` and multiplies it by each
        // coordinate, at `0x1eb51` and `0x1eb69`. Inert: a product commutes.
        let m = move5x6(
            c,
            VERTEX_STEP_FACTOR * v.0,
            VERTEX_STEP_FACTOR * v.1,
            1,
            None,
            false,
        );
        // `getZoneFromCRSCentroid(zLevel, 0, m)` (RI7H.ec:299) is `fromCentroid` for a level
        // of 19 or less (RI7H.ec:99-106), which the level of a zone with sub-zones is.
        m.and_then(|m| from_centroid(z_level, m.0, m.1)) == Some(*hay_stack)
    })
}

/// Whether the centroid `c` of an order is another than `sz_centroid`, the canonical centroid
/// of the sub-zone sought: the callback of the engine's walk, `findSubZone`
/// (`RI7H.ec:212-222`), compiled at `0x1e020` to `0x1e078`. `c` is made canonical first
/// (`canonicalize5x6`, `ri5x6.ec:1562-1606`), and the two are one centroid where each
/// coordinate of the one lies within 1e-11 of the other's.
///
/// `true` is the engine's "go on", `false` its "this is the one".
fn find_sub_zone(sz_centroid: (f64, f64), c: (f64, f64)) -> bool {
    let centroid = canonicalize5x6(c.0, c.1);
    // faithful (FZ-1): the eC answers `false` where `fabs(centroid.x - szCentroid.x) < 1E-11
    // && fabs(centroid.y - szCentroid.y) < 1E-11`, and `true` otherwise (RI7H.ec:217-220).
    // The library takes the first test as written (`comisd` at `0x1e055`, `jbe` at
    // `0x1e059`) and the second by its complement: it answers whether the second difference
    // is 1e-11 or more (`comisd` at `0x1e06c`, `setae` at `0x1e070`). Inert: the two forms
    // part only where that difference is no number, which the eC's takes for another
    // centroid and the library's for the one sought, and no centroid of an order is such a
    // point.
    if (centroid.0 - sz_centroid.0).abs() < 1e-11 {
        (centroid.1 - sz_centroid.1).abs() >= 1e-11
    } else {
        true
    }
}

/// The index of `sub_zone` in the order of `parent`, as the engine's walk finds it:
/// `RhombicIcosahedral7H::getSubZoneIndex` (`RI7H.ec:224-239`), compiled at `0x28a30` to
/// `0x28b64` in the eC's own order. It answers 0 for two zones of one level, whatever they are,
/// and -1 for a coarser `sub_zone` and for one that [`zone_has_sub_zone`] refuses; otherwise it
/// generates the order from its head and answers the index of the first centroid that
/// [`find_sub_zone`] takes for the sub-zone's own, or -1 where none is.
///
/// The cost grows with the index, one centroid of the generator for each entry before it, so
/// the caller bounds the length of the order: see [`sub_zone_index`]. The sub-zone level is
/// at most 19, as for [`iterate_sub_zones`].
///
/// A centroid that the generator does not reach (see [`iterate_sub_zones`]) is no centroid,
/// and so not the one sought.
fn index_by_the_walk(parent: &Z, sub_zone: &Z) -> i64 {
    let (level, sz_level) = (zone_level(parent), zone_level(sub_zone));
    if sz_level == level {
        return 0;
    }
    if sz_level < level || !zone_has_sub_zone(parent, sub_zone) {
        return -1;
    }
    // `canonicalize5x6(subZone.centroid, zCentroid)` (RI7H.ec:235): the centroid property at
    // `0x28ae3`, and the fold of `ri5x6.ec:1562-1606` at `0x28af3`.
    let zc = zone_centroid(sub_zone);
    let zc = canonicalize5x6(zc.0, zc.1);
    iterate_sub_zones(parent, sz_level - level, -1, &mut |_, centroid| {
        centroid.is_none_or(|c| find_sub_zone(zc, c))
    })
}

/// The address of no zone: base cell 15 and no digit, which the Z7 packing writes as its null
/// identifier, every bit set. It is no valid address, and names nothing else.
///
/// No eC counterpart: the eC holds its null zone as an identifier, never as an address.
pub(super) fn null_address() -> Address {
    Address::new(0xF, &[])
}

/// The address at level `res` of the zone that holds the centroid `c`, as the engine's
/// `getZoneFromCRSCentroid` answers in the 5x6 CRS (`RI7H_Z7.ec:527-530`: `fromCentroid`, and
/// then `from7H`), through [`HexA7::quantize`]; and [`null_address`] where the centroid is
/// none, where the engine answers its null zone, and where it answers an identifier that it
/// cannot read back, for which `from7H` answers nothing here.
pub(super) fn address_of_centroid(c: Option<(f64, f64)>, res: u8) -> Address {
    c.and_then(|(x, y)| HexA7::quantize(PlanarPoint { face: -1, x, y }, res))
        .unwrap_or_else(null_address)
}

/// The zone that the address names and the level of its sub-zones `depth` levels below, or
/// `None` where the address names no zone (it fails [`HexA7::is_valid_address`]) or the
/// sub-zones would lie beyond level 19, the finest that has geometry. Every answer below
/// begins here, so that none of them reaches the generator, the first centroid or the count
/// with a level that their powers of seven are not held for, whatever the depth asked.
///
/// No eC counterpart: the engine checks neither. Its `countSubZones` answers its formula for
/// a level beyond 19, and its `getSubZones` answers nothing there (`dggrs.ec:195-222`).
fn sub_zone_request(a: &Address, depth: u8) -> Option<(Z, u8)> {
    let sz_level = a.len() + usize::from(depth);
    (HexA7::is_valid_address(a) && sz_level < NULL_GEOMETRY_LEVEL)
        .then(|| (zone_from_steps(a.base, a.digits()), sz_level as u8))
}

/// The number of sub-zones `depth` levels below the zone of the address, as the engine's
/// `countSubZones` answers (`RI7H_Z7.ec:492-500`), through [`sub_zones_count`]; `None` where
/// [`sub_zone_request`] refuses.
pub(super) fn count(a: &Address, depth: u8) -> Option<u64> {
    let (z, _) = sub_zone_request(a, depth)?;
    Some(sub_zones_count(zone_npoints(&z), i64::from(depth)))
}

/// The first sub-zone `depth` levels below the zone of the address: the first centroid,
/// quantised at the sub-zone level, as `I7HZone::getFirstSubZone` answers
/// (`RI7H.ec:3118-3124`, which the Z7 class reaches at `RI7H_Z7.ec:576-579`); `None` where
/// [`sub_zone_request`] refuses, and [`null_address`] where no zone holds the centroid.
///
/// Two departures from the engine. At depth 0 this answers the zone itself, as the grid does
/// before it asks: the engine has no case for that depth, quantises a point beside a vertex of
/// the zone at the zone's own level, and so names a neighbour at most zones, while its
/// `getSubZones` names the zone itself. And at the pentagon of the south pole at an odd level,
/// at depth 2, the engine's call ends the process (see [`get_first_sub_zone_centroid`]); here
/// the answer is the zone of the centroid that its arithmetic gives, which is the first entry
/// of the order.
pub(super) fn first_sub_zone(a: &Address, depth: u8) -> Option<Address> {
    let (z, sz_level) = sub_zone_request(a, depth)?;
    if depth == 0 {
        return Some(*a);
    }
    let centroid = get_first_sub_zone_centroid(&z, i64::from(depth), None);
    Some(address_of_centroid(centroid, sz_level))
}

/// Every sub-zone `depth` levels below the zone of the address, in the engine's order: the
/// centroids of [`iterate_sub_zones`], each stored at its index in a list of the count's
/// length, as `getSubZoneCentroids` stores them (`RI7H.ec:3707-3726`), and each quantised at
/// the sub-zone level, as `getSubZones` quantises them (`dggrs.ec:195-222`).
///
/// `None` where [`sub_zone_request`] refuses, and above
/// [`crate::grid::max_materialised_sub_zones`], far below the `2^28` centroids at which the
/// engine's own list ends the process and the `2^32` from which it answers nothing; the grid
/// refuses both before it asks.
///
/// An entry is [`null_address`] where no zone holds its centroid, where the engine's entry is
/// its null zone or an identifier that it cannot read back: the length of the order and the
/// place of every other entry are the engine's. Whatever else the engine's order holds in the
/// broken seams is kept as it is: a zone named twice, a zone that is no descendant.
pub(super) fn sub_zones(a: &Address, depth: u8) -> Option<Vec<Address>> {
    let (z, sz_level) = sub_zone_request(a, depth)?;
    let n_sub_zones = sub_zones_count(zone_npoints(&z), i64::from(depth));
    if n_sub_zones > crate::grid::max_materialised_sub_zones() {
        return None;
    }
    if depth == 0 {
        return Some(vec![*a]);
    }
    let mut order = vec![null_address(); n_sub_zones as usize];
    iterate_sub_zones(&z, i64::from(depth), -1, &mut |index, centroid| {
        // `centroids[(uint)index] = centroid` (RI7H.ec:3709): an index beyond the count would
        // write outside the engine's array; here it is dropped.
        if let Some(entry) = order.get_mut(index as usize) {
            *entry = address_of_centroid(centroid, sz_level);
        }
        true
    });
    Some(order)
}

/// The sub-zone at `index` in the order of [`sub_zones`], found by the generator's search,
/// without building the order: `getSubZoneAtIndex` (`RI7H.ec:241-256`), less its fault at the
/// pentagons (see [`iterate_sub_zones`]), and with index 0 found by the same search and not
/// by `getFirstSubZone`, which gives the same zone wherever it gives one. `None` where
/// [`sub_zone_request`] refuses and where `index` is not below the count; [`null_address`]
/// where the entry at `index` is no zone. At depth 0 this answers the zone itself, as
/// [`first_sub_zone`] does.
pub(super) fn sub_zone_at_index(a: &Address, depth: u8, index: u64) -> Option<Address> {
    let (z, sz_level) = sub_zone_request(a, depth)?;
    if index >= sub_zones_count(zone_npoints(&z), i64::from(depth)) {
        return None;
    }
    if depth == 0 {
        return Some(*a);
    }
    let centroid = sub_zone_centroid_at_index(&z, i64::from(depth), index);
    Some(address_of_centroid(centroid, sz_level))
}

/// The index of the zone of `sub` in the order of the sub-zones of the zone of `parent`, by
/// the engine's own walk ([`index_by_the_walk`]), as its `getSubZoneIndex` answers
/// (`RI7H_Z7.ec:599-602`): `Some(Some(i))` where the walk finds it, `Some(None)` where it
/// finds none. The index is the walk's and no more: whether the entry of
/// the order at `i` is the zone of `sub` is for the caller to confirm, as the grid does, since
/// the walk compares centroids and the order holds what the quantiser makes of them, which in
/// the broken seams need not be the same zone.
///
/// The crate's own rule stands where the two share a level, which the grid settles before it
/// asks: a zone is its own sub-zone at index 0, and no other zone of its level is, where the
/// engine answers 0 for any two zones of one level. A coarser `sub` is no sub-zone.
///
/// `None` where either address names no zone, where [`sub_zone_request`] refuses the depth,
/// and where the order is longer than [`crate::grid::max_materialised_sub_zones`]: the walk
/// costs one centroid for each entry before the one sought, so it is bounded as the list is,
/// and the grid refuses such a request before it asks.
pub(super) fn sub_zone_index(parent: &Address, sub: &Address) -> Option<Option<u64>> {
    if !HexA7::is_valid_address(parent) || !HexA7::is_valid_address(sub) {
        return None;
    }
    if sub.len() <= parent.len() {
        return Some((sub == parent).then_some(0));
    }
    // An address holds twenty digits at most, so the difference is a `u8`.
    let depth = (sub.len() - parent.len()) as u8;
    let (z, _) = sub_zone_request(parent, depth)?;
    if sub_zones_count(zone_npoints(&z), i64::from(depth))
        > crate::grid::max_materialised_sub_zones()
    {
        return None;
    }
    let sub_zone = zone_from_steps(sub.base, sub.digits());
    Some(u64::try_from(index_by_the_walk(&z, &sub_zone)).ok())
}

#[cfg(test)]
mod tests {
    //! Every expected value below is the shipped engine's own answer, read from `libdggal.so`,
    //! BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, and held as bits where it is a double;
    //! the comparisons are bit for bit. The three aperture-7 grids share these answers, which are
    //! taken in the 5x6 layout, before any projection.
    use super::*;

    /// The zone of a Z7 identifier in text: two digits of base cell, then the digits.
    fn zone(text: &str) -> Z {
        let digits: Vec<u8> = text.bytes().skip(2).map(|b| b - b'0').collect();
        zone_from_steps(text[..2].parse().unwrap(), &digits)
    }

    /// The bits of a point, or of a direction.
    type Bits = (u64, u64);

    fn bits(p: (f64, f64)) -> Bits {
        (p.0.to_bits(), p.1.to_bits())
    }

    /// Asserts the first centroid of each `(zone, depth)` against the engine's, and names every
    /// row that differs.
    fn assert_first_centroids(rows: &[(&str, i64, Bits)]) {
        let wrong: Vec<String> = rows
            .iter()
            .filter_map(|&(text, depth, engine)| {
                let got = get_first_sub_zone_centroid(&zone(text), depth, None).map(bits);
                (got != Some(engine)).then(|| {
                    format!("{text} at depth {depth}: {got:x?}, where the engine has {engine:x?}")
                })
            })
            .collect();
        assert!(wrong.is_empty(), "{wrong:#?}");
    }

    #[test]
    fn the_folded_constants_match_the_engines_rodata() {
        // Site FS-1: the six numerators.
        assert_eq!(FS1_NUMERATOR_MINUS_2.to_bits(), 0xbfe5_5555_5555_5555); // 0x4a1d8
        assert_eq!(FS1_NUMERATOR_MINUS_10.to_bits(), 0xc00a_aaaa_aaaa_aaaa); // 0x4a1e0
        assert_eq!(FS1_NUMERATOR_2.to_bits(), 0x3fe5_5555_5555_5555); // 0x4a1e8
        assert_eq!(FS1_NUMERATOR_MINUS_4.to_bits(), 0xbff5_5555_5555_5555); // 0x4a1f0
        assert_eq!(FS1_NUMERATOR_MINUS_5.to_bits(), 0xbffa_aaaa_aaaa_aaaa); // 0x4a1f8
        assert_eq!(FS1_NUMERATOR_1.to_bits(), 0x3fd5_5555_5555_5555); // 0x49ba0
        // Each is the product of its `k` and the double nearest to a third, which is how gcc
        // folded it; for -5 and -10 the double nearest to `k/3` is another.
        let third = 1.0_f64 / 3.0;
        for (k, numerator) in [
            (2.0, FS1_NUMERATOR_2),
            (1.0, FS1_NUMERATOR_1),
            (-2.0, FS1_NUMERATOR_MINUS_2),
            (-4.0, FS1_NUMERATOR_MINUS_4),
            (-5.0, FS1_NUMERATOR_MINUS_5),
            (-10.0, FS1_NUMERATOR_MINUS_10),
        ] {
            assert_eq!((k * third).to_bits(), numerator.to_bits(), "k = {k}");
            assert_eq!(
                (k / 3.0).to_bits() != numerator.to_bits(),
                k == -5.0 || k == -10.0,
                "k = {k}"
            );
        }
        assert_eq!((-5.0_f64 / 3.0).to_bits(), 0xbffa_aaaa_aaaa_aaab);
        assert_eq!((-10.0_f64 / 3.0).to_bits(), 0xc00a_aaaa_aaaa_aaab);

        // The net: `triWidthOver2`, and sites TI-1 and TI-2, each the double that the eC's own
        // product gives.
        assert_eq!(TRI_WIDTH_OVER_2.to_bits(), 0x414d_4696_7cb2_2755); // 0x49c90
        assert_eq!(TRI_WIDTH_OVER_2, 3_837_228.974_186_817_6);
        assert_eq!(TRI_WIDTH_OVER_2_SQRT3.to_bits(), 0x4159_5a80_e2c9_cd20); // 0x49c98
        assert_eq!(
            TRI_WIDTH_OVER_2 * 1.732_050_807_568_877_2,
            TRI_WIDTH_OVER_2_SQRT3
        );
        assert_eq!(FIVE_TRI_WIDTH_OVER_2.to_bits(), 0x4172_4c1e_0def_5895); // 0x4a1d0
        assert_eq!(5.0 * TRI_WIDTH_OVER_2, FIVE_TRI_WIDTH_OVER_2);

        // The tolerances of the uppermost vertex and the epsilon of the two hemisphere tests.
        assert_eq!(1e-6_f64.to_bits(), 0x3eb0_c6f7_a0b5_ed8d); // 0x49f50
        assert_eq!(1e-8_f64.to_bits(), 0x3e45_798e_e230_8c3a); // 0x4a1c8
        assert_eq!(1e-11_f64.to_bits(), 0x3da5_fd7f_e179_6495); // 0x49ec0

        // Site VD-1: the five offsets of a hexagon's vertices that the directions are made of.
        assert_eq!(EVEN_HEX_A.to_bits(), 0x4002_aaaa_aaaa_aaab); // 0x4a110
        assert_eq!(EVEN_HEX_B.to_bits(), 0x4012_aaaa_aaaa_aaab); // 0x4a118
        assert_eq!(ODD_HEX_A.to_bits(), 0x3ff5_5555_5555_5555); // 0x49ba8
        assert_eq!(ODD_HEX_B.to_bits(), 0x3ffa_aaaa_aaaa_aaab); // 0x49c00
        assert_eq!(ODD_HEX_C.to_bits(), 0x3fd5_5555_5555_5555); // 0x49ba0
        // Site ZH-1: the 99 hundredths of the step towards a vertex. Site FZ-1 compares with
        // the 1e-11 pinned above.
        assert_eq!(VERTEX_STEP_FACTOR.to_bits(), 0x3fef_ae14_7ae1_47ae); // 0x4a120
    }

    // The four tests that follow pin site FS-1, one for each pair of numerators. At every row the
    // eC's own order, `k / (3.0 * szp)`, gives another centroid, a unit or more in the last place
    // away in one coordinate or in both: put back to that order, the pair fails its test.

    #[test]
    fn the_first_centroid_at_an_odd_depth_below_an_odd_level_is_the_engines() {
        // The numerators of 2 and 1.
        assert_first_centroids(&[(
            "01000000000000003",
            1,
            (0x3e8f_0b04_ec9e_dd69, 0x3e7f_0b04_ee9e_dd69),
        )]);
    }

    #[test]
    fn the_first_centroid_at_an_odd_depth_below_an_even_level_is_the_engines() {
        // The numerators of -4 and -5.
        assert_first_centroids(&[
            ("0004", 1, (0x3fed_6343_eb1a_1f59, 0x3fa4_e5e0_a72f_0552)),
            ("0050", 1, (0x3fe7_829c_bc14_e5de, 0x3faf_58d0_fac6_87d6)),
            (
                "0100000000000004",
                1,
                (0x3e87_4843_b1ee_4c20, 0x3ea1_7632_c572_b916),
            ),
        ]);
    }

    #[test]
    fn the_first_centroid_at_an_even_depth_below_an_odd_level_is_the_engines() {
        // The numerators of -10 and -2. The third row is the pentagon of the north pole, where
        // the offset is turned three times.
        assert_first_centroids(&[
            ("014", 2, (0x3fbf_58d0_fac6_87d6, 0x3fcc_bc14_e5e0_a733)),
            ("016", 2, (0x3fda_1f58_d0fa_c68a, 0x3fb4_e5e0_a72f_0535)),
            ("000", 2, (0x3fea_1f58_d0fa_c687, 0x3f94_e5e0_a72f_0568)),
            (
                "01000000000000003",
                2,
                (0x3e61_bd27_6083_3e42, 0x3e66_2c71_3e80_a60d),
            ),
        ]);
    }

    #[test]
    fn the_first_centroid_at_an_even_depth_below_an_even_level_is_the_engines() {
        // The numerators of -4 and -2.
        assert_first_centroids(&[
            (
                "0100000000000004",
                2,
                (0x3e87_4843_b1ee_4c20, 0x3ea7_4843_b1ee_4c1d),
            ),
            (
                "0100040004000400",
                2,
                (0x3e87_4843_b1ee_4c20, 0x3f95_5526_86b7_e79e),
            ),
        ]);
    }

    #[test]
    fn the_count_is_the_engines_closed_form() {
        // Depths 1 to 10, a hexagon and a pentagon.
        const HEXAGON: [u64; 10] = [
            13,
            55,
            379,
            2_449,
            17_053,
            117_991,
            825_259,
            5_767_201,
            40_365_613,
            282_492_055,
        ];
        const PENTAGON: [u64; 10] = [
            11,
            46,
            316,
            2_041,
            14_211,
            98_326,
            687_716,
            4_806_001,
            33_638_011,
            235_410_046,
        ];
        for (i, (&hexagon, &pentagon)) in HEXAGON.iter().zip(&PENTAGON).enumerate() {
            let depth = i as i64 + 1;
            assert_eq!(sub_zones_count(6, depth), hexagon, "depth {depth}");
            assert_eq!(sub_zones_count(5, depth), pentagon, "depth {depth}");
        }
        // A zone is its own one sub-zone at depth 0, and the engine answers 1 below it too.
        for depth in [0, -1, i64::MIN] {
            assert_eq!(sub_zones_count(6, depth), 1, "depth {depth}");
            assert_eq!(sub_zones_count(5, depth), 1, "depth {depth}");
        }
    }

    #[test]
    fn the_count_is_exact_where_a_power_of_seven_is_no_double() {
        // The two greatest counts a zone can have within level 19: a pentagon of level 0 at depth
        // 19 and a hexagon of level 1 at depth 18. `7^19` is no double, and a count taken through
        // one answers a unit more for the first.
        let (pentagon, hexagon) = (zone("00"), zone("001"));
        assert_eq!((zone_level(&pentagon), zone_npoints(&pentagon)), (0, 5));
        assert_eq!((zone_level(&hexagon), zone_npoints(&hexagon)), (1, 6));
        assert_eq!(
            sub_zones_count(zone_npoints(&pentagon), 19),
            9_499_079_489_284_316
        );
        assert_eq!(
            sub_zones_count(zone_npoints(&hexagon), 18),
            1_628_413_638_264_055
        );
        // A hexagon at depth 19 lies outside what a grid can ask, since every zone of level 0 is
        // a pentagon; the closed form is pinned there on its own.
        assert_eq!(sub_zones_count(6, 19), 11_398_895_387_141_179);
    }

    #[test]
    fn the_count_wraps_as_the_engines_integers_and_never_panics() {
        // (depth, hexagon, pentagon). Depths 20 and 21 are exact still; from depth 22 the engine's
        // 64-bit integers wrap, and it divides the wrapped sum as a signed one.
        for (depth, hexagon, pentagon) in [
            (20, 79_792_266_580_087_249_u64, 66_493_555_483_406_041_u64),
            (21, 558_545_865_495_660_253, 465_454_887_913_050_211),
            (22, 835_363_704_942_056_189, 183_726_863_182_003_724),
            (23, 18_145_375_313_112_774_252, 1_286_088_038_978_481_498),
            (24, 964_875_966_075_668_236, 18_225_988_263_566_522_279),
        ] {
            assert_eq!(sub_zones_count(6, depth), hexagon, "depth {depth}");
            assert_eq!(sub_zones_count(5, depth), pentagon, "depth {depth}");
        }
        // Beyond its table of powers the engine goes through `exp`, which is not followed; the
        // function answers all the same.
        for depth in [25, 64, 1 << 40, i64::MAX] {
            let _ = sub_zones_count(6, depth);
            let _ = sub_zones_count(5, depth);
        }
    }

    #[test]
    fn the_uppermost_vertex_is_the_engines() {
        // (zone, the vertex, how many share its height). The pentagons of the two poles, at levels
        // 0 and 1, and of base cells 1 and 6; hexagons with one uppermost vertex and with two.
        const ROWS: [(&str, Bits, i32); 18] = [
            ("0004", (0x3fee_79e7_9e79_e79f, 0x3fa8_6186_1861_8630), 2),
            ("0050", (0x3fe8_6186_1861_8616, 0x3fb8_6186_1861_8618), 1),
            ("014", (0x3fc8_6186_1861_8618, 0x3fce_79e7_9e79_e7a2), 1),
            ("016", (0x3fde_79e7_9e79_e7a0, 0x3fb8_6186_1861_8614), 1),
            ("000", (0x4013_cf3c_f3cf_3cf4, 0x4010_c30c_30c3_0c31), 2),
            (
                "01000000000000003",
                (0x4014_0000_07c2_c13b, 0x4014_0000_03e1_609e),
                1,
            ),
            (
                "0100000000000004",
                (0x3e9b_29a4_4f96_0379, 0x3eab_29a4_4f96_0377),
                1,
            ),
            (
                "0100040004000400",
                (0x3e9b_29a4_4f96_0379, 0x3f95_552e_4979_22ed),
                1,
            ),
            ("00", (0x4012_aaaa_aaaa_aaab, 0x4011_5555_5555_5555), 2),
            ("11", (0x3fd5_5555_5555_5555, 0x3ffa_aaaa_aaaa_aaab), 2),
            ("110", (0x3fa8_6186_1861_8618, 0x3ffc_f3cf_3cf3_cf3d), 1),
            ("11000", (0x3f7b_dd2b_8994_06f6, 0x3fff_908b_51d9_afe4), 1),
            ("01", (0xbfd5_5555_5555_5555, 0xbfe5_5555_5555_5555), 2),
            ("010", (0xbfc8_6186_1861_8618, 0xbfce_79e7_9e79_e79e), 2),
            ("06", (0x3fd5_5555_5555_5557, 0x3fe5_5555_5555_5557), 1),
            ("060", (0x3fa8_6186_1861_8618, 0x3fe9_e79e_79e7_9e7a), 1),
            ("0064156", (0x3fff_68bd_13a7_6eb6, 0x3ff5_a0f6_cb81_9dfb), 1),
            (
                "006415654636",
                (0x3fff_55f9_1ae7_832f, 0x3ff5_a846_78e3_a9f3),
                1,
            ),
        ];
        for (text, vertex, n_top) in ROWS {
            let (v, n) = get_top_ico_vertex(&zone(text));
            assert_eq!((bits(v), n), (vertex, n_top), "{text}");
        }
    }

    /// The direction `(0.3, 0.2)` as the engine leaves it: untouched, and turned by sixty degrees
    /// once, twice, three times and five times.
    const UNTURNED: Bits = (0x3fd3_3333_3333_3333, 0x3fc9_9999_9999_999a);
    const TURNED_ONCE: Bits = (0x3fc9_9999_9999_999a, 0xbfb9_9999_9999_9998);
    const TURNED_TWICE: Bits = (0xbfb9_9999_9999_9998, 0xbfd3_3333_3333_3333);
    const TURNED_THRICE: Bits = (0xbfd3_3333_3333_3333, 0xbfc9_9999_9999_999a);
    const TURNED_FIVE_TIMES: Bits = (0x3fb9_9999_9999_9998, 0x3fd3_3333_3333_3333);

    #[test]
    fn the_first_centroid_and_the_direction_it_turns_are_the_engines() {
        // (zone, depth, the first centroid, the direction (0.3, 0.2) as the engine leaves it).
        // The pentagon of the north pole turns the direction three times, save at an odd level and
        // depth 2, where it turns the offset three times and the direction twice; that of the
        // south pole at an odd level and depth 2 turns it five times; a northern pentagon and a
        // hexagon with two uppermost vertices north of the interruption turn it once.
        const ROWS: [(&str, i64, Bits, Bits); 34] = [
            (
                "0004",
                1,
                (0x3fed_6343_eb1a_1f59, 0x3fa4_e5e0_a72f_0552),
                TURNED_ONCE,
            ),
            (
                "0004",
                2,
                (0x3fee_0a72_f053_9783, 0x3faf_58d0_fac6_87ee),
                TURNED_ONCE,
            ),
            (
                "0004",
                3,
                (0x3fee_5219_6047_a670, 0x3fa7_e225_515a_4f35),
                TURNED_ONCE,
            ),
            (
                "0050",
                1,
                (0x3fe7_829c_bc14_e5de, 0x3faf_58d0_fac6_87d6),
                UNTURNED,
            ),
            (
                "000",
                1,
                (0x4013_6db6_db6d_b6dc, 0x4010_9249_2492_4925),
                TURNED_THRICE,
            ),
            (
                "000",
                2,
                (0x3fea_1f58_d0fa_c687, 0x3f94_e5e0_a72f_0568),
                TURNED_TWICE,
            ),
            (
                "000",
                3,
                (0x4013_c14e_5e0a_72f1, 0x4010_bc14_e5e0_a72f),
                TURNED_THRICE,
            ),
            (
                "00",
                1,
                (0x4013_6db6_db6d_b6dc, 0x4012_4924_9249_2492),
                TURNED_THRICE,
            ),
            (
                "00",
                2,
                (0x4013_6db6_db6d_b6dc, 0x4011_b6db_6db6_db6d),
                TURNED_THRICE,
            ),
            (
                "00",
                3,
                (0x4012_c687_d634_3eb2, 0x4011_7829_cbc1_4e5e),
                TURNED_THRICE,
            ),
            (
                "11",
                1,
                (0x3fc2_4924_9249_2492, 0x3ff6_db6d_b6db_6db7),
                UNTURNED,
            ),
            (
                "11",
                2,
                (0x3fc2_4924_9249_2492, 0x3ff9_2492_4924_924a),
                UNTURNED,
            ),
            (
                "11",
                3,
                (0x3fd3_9782_9cbc_14e6, 0x3ffa_1f58_d0fa_c688),
                UNTURNED,
            ),
            (
                "110",
                1,
                (0x3fc2_4924_9249_2492, 0x3ffd_b6db_6db6_db6e),
                UNTURNED,
            ),
            (
                "110",
                2,
                (0x4010_bc14_e5e0_a72f, 0x4017_eb1a_1f58_d0fb),
                TURNED_FIVE_TIMES,
            ),
            (
                "110",
                3,
                (0x3faf_58d0_fac6_87d6, 0x3ffd_0fac_687d_6344),
                UNTURNED,
            ),
            (
                "11000",
                2,
                (0x4010_1ade_69fb_8599, 0x4017_fd03_bb55_d4b6),
                TURNED_FIVE_TIMES,
            ),
            (
                "01",
                1,
                (0x4011_b6db_6db6_db6e, 0x4011_2492_4924_9249),
                TURNED_ONCE,
            ),
            (
                "01",
                2,
                (0x4012_4924_9249_2493, 0x4011_b6db_6db6_db6d),
                TURNED_ONCE,
            ),
            (
                "01",
                3,
                (0x4012_87d6_343e_b1a2, 0x4011_4e5e_0a72_f053),
                TURNED_ONCE,
            ),
            (
                "010",
                1,
                (0x4013_6db6_db6d_b6db, 0x4012_db6d_b6db_6db7),
                TURNED_ONCE,
            ),
            (
                "010",
                2,
                (0x4013_2f05_3978_29cc, 0x4013_43eb_1a1f_58d1),
                TURNED_ONCE,
            ),
            (
                "010",
                3,
                (0x4013_43eb_1a1f_58d1, 0x4013_0539_7829_cbc1),
                TURNED_ONCE,
            ),
            (
                "06",
                1,
                (0x3fc2_4924_9249_2496, 0x3fdb_6db6_db6d_b6df),
                UNTURNED,
            ),
            (
                "06",
                2,
                (0x3fc2_4924_9249_2496, 0x3fe2_4924_9249_2494),
                UNTURNED,
            ),
            (
                "06",
                3,
                (0x3fd3_9782_9cbc_14e8, 0x3fe4_3eb1_a1f5_8d11),
                UNTURNED,
            ),
            (
                "060",
                1,
                (0x3fc2_4924_9249_2492, 0x3feb_6db6_db6d_b6dc),
                UNTURNED,
            ),
            (
                "060",
                2,
                (0x4013_eb1a_1f58_d0fb, 0x4017_2f05_3978_29cc),
                UNTURNED,
            ),
            (
                "060",
                3,
                (0x3faf_58d0_fac6_87d6, 0x3fea_1f58_d0fa_c688),
                UNTURNED,
            ),
            (
                "0064156",
                1,
                (0x3fff_70b3_2017_e226, 0x3ff5_a4f1_d1b9_d7b3),
                UNTURNED,
            ),
            (
                "0064156",
                2,
                (0x3fff_630d_53e9_6566, 0x3ff5_9fd3_a528_68eb),
                UNTURNED,
            ),
            (
                "0064156",
                3,
                (0x3fff_69e0_3a00_a3c6, 0x3ff5_a188_5eae_3883),
                UNTURNED,
            ),
            (
                "006415654636",
                1,
                (0x3fff_55ed_38af_a05d, 0x3ff5_a837_9e1d_ce6d),
                UNTURNED,
            ),
            (
                "006415654636",
                2,
                (0x3fff_55ed_38af_a05d, 0x3ff5_a840_87c7_b88a),
                UNTURNED,
            ),
        ];
        for (text, depth, centroid, direction) in ROWS {
            let mut s = (0.3, 0.2);
            let c = get_first_sub_zone_centroid(&zone(text), depth, Some(&mut s));
            assert_eq!(
                (c.map(bits), bits(s)),
                (Some(centroid), direction),
                "{text} at depth {depth}"
            );
            // The centroid does not depend on a direction being given.
            assert_eq!(
                get_first_sub_zone_centroid(&zone(text), depth, None).map(bits),
                Some(centroid),
                "{text} at depth {depth}, no direction given"
            );
        }
    }

    #[test]
    fn the_south_polar_pentagon_answers_where_the_engine_ends_the_process() {
        // At an odd level and depth 2 the engine turns a direction it was not given, through a
        // null pointer. With none given there is none to turn, and the centroid is the one its
        // arithmetic gives, which is the one it answers when it is given a direction.
        for (text, centroid) in [
            ("110", (0x4010_bc14_e5e0_a72f, 0x4017_eb1a_1f58_d0fb)),
            ("11000", (0x4010_1ade_69fb_8599, 0x4017_fd03_bb55_d4b6)),
        ] {
            let z = zone(text);
            assert_eq!((z.root, zone_npoints(&z), zone_level(&z) & 1), (0xB, 5, 1));
            assert_eq!(
                get_first_sub_zone_centroid(&z, 2, None).map(bits),
                Some(centroid),
                "{text}"
            );
        }
    }

    /// The order of `z` at `depth`, as the generator delivers it: every centroid, the indices
    /// checked to run from 0 without a gap.
    fn order(z: &Z, depth: i64) -> Vec<Option<(f64, f64)>> {
        let mut centroids = vec![];
        let stop = iterate_sub_zones(z, depth, -1, &mut |index, centroid| {
            assert_eq!(index, centroids.len() as i64, "the order of the indices");
            centroids.push(centroid);
            true
        });
        assert_eq!(stop, -1);
        centroids
    }

    /// The Z7 text of an address.
    fn text_of(a: &Address) -> String {
        use crate::Indexing;
        use crate::indexings::Z7;
        Z7::to_text(Z7::encode(a))
    }

    #[test]
    fn the_scanline_counts_of_an_odd_depth_are_the_engines_ceilings() {
        // Site IT-1. The library takes each ceiling in doubles: the count of the scanlines
        // between two vertices, or twice it, converted and multiplied by the double nearest to a
        // third (the `.rodata` word at `0x49ba0`), raised to the next integer, a half added, and
        // truncated. Over every power of seven that a depth within level 19 asks, and nine
        // powers beyond, that is the integer the port computes.
        let third = f64::from_bits(0x3fd5_5555_5555_5555);
        let as_compiled = |n: i64| (((n as f64) * third).ceil() + 0.5) as i64;
        let mut n_inter_sl: i64 = 1;
        for exponent in 0..=18 {
            assert_eq!(
                odd_depth_scanline_counts(n_inter_sl),
                (
                    as_compiled(n_inter_sl),
                    as_compiled(n_inter_sl + n_inter_sl)
                ),
                "7^{exponent}"
            );
            // A power of seven is no multiple of three, so each ceiling is the quotient raised.
            assert_eq!(
                odd_depth_scanline_counts(n_inter_sl),
                (n_inter_sl / 3 + 1, 2 * n_inter_sl / 3 + 1),
                "7^{exponent}"
            );
            n_inter_sl *= 7;
        }
        // The two part at `7^19`, which is no double; the greatest power asked is `7^9`, at
        // depth 19.
        assert_eq!(n_inter_sl, 11_398_895_185_373_143);
        assert_ne!(
            odd_depth_scanline_counts(n_inter_sl),
            (
                as_compiled(n_inter_sl),
                as_compiled(n_inter_sl + n_inter_sl)
            )
        );
        assert_eq!(pow7((19 - 1) / 2), 40_353_607);
    }

    #[test]
    fn the_folded_turns_are_the_turns_made_one_after_another() {
        // Sites IT-3, IT-4 and IT-5. The directions that the generator holds are the six turns
        // of the direction of the scanlines, which is `(c, c)` at an even sub-zone level and
        // `(3c, c + c)` at an odd one (site IT-2), `c` being one over a power of seven. On each
        // of them, at every sub-zone level, the pair that the library forms in one step is, as
        // numbers, the pair that the eC's turns give one after another; a zero may differ in its
        // sign, which is why the comparison is of numbers and why the port takes the library's
        // form.
        let mut directions = 0;
        for sz_level in 1..=19_i64 {
            let odd = sz_level & 1 != 0;
            let c2c = 1.0 / pow7((sz_level + i64::from(odd)) / 2) as f64;
            let mut t = if odd {
                (c2c * 3.0, c2c + c2c)
            } else {
                (c2c, c2c)
            };
            let first = t;
            for _ in 0..6 {
                let twice_clockwise = turned(turned(t, true), true);
                assert_eq!(turned_120_degrees_clockwise(t), twice_clockwise, "{t:?}");
                let mut anticlockwise = t;
                for _ in 0..4 {
                    anticlockwise = turned(anticlockwise, false);
                }
                assert_eq!(turned_120_degrees_clockwise(t), anticlockwise, "{t:?}");
                assert_eq!(
                    turned_five_times_anticlockwise(t),
                    turned(anticlockwise, false),
                    "{t:?}"
                );
                // A turn anticlockwise and then one clockwise give the abscissa back.
                assert_eq!(turned(turned(t, false), true).0, t.0, "{t:?}");
                t = turned(t, false);
                directions += 1;
            }
            // Six turns come back to the first direction: the six are closed under the turns.
            assert_eq!(t, first, "level {sz_level}");
        }
        assert_eq!(directions, 19 * 6);
        // One turn is the eC's `rotate5x6Offset`.
        assert_eq!(turned((3.0, 2.0), true), (1.0, 3.0));
        assert_eq!(turned((3.0, 2.0), false), (2.0, -1.0));
    }

    #[test]
    fn the_order_below_the_south_polar_pentagon_at_depth_2_is_the_engines() {
        // The pentagon of the south pole at level 1, where the scanlines past the pentagon start
        // from a point reached by turns that no interruption gives, and where depth 2 has an arm
        // of its own. The engine's `getFirstSubZone` ends the process at this zone and depth;
        // its `getSubZones` and its centroids are answered, and are these, read in a process of
        // their own. The comparison with the live engine holds the same centroids; these
        // literals hold them where no engine is at hand.
        const ENGINE: [(&str, Bits); 46] = [
            ("11022", (0x4010_bc14_e5e0_a72f, 0x4017_eb1a_1f58_d0fb)),
            ("11026", (0x3fa4_e5e0_a72f_0556, 0x3ffd_6343_eb1a_1f58)),
            ("11611", (0x3fba_1f58_d0fa_c696, 0x3ffe_0a72_f053_9782)),
            ("11255", (0x4010_687d_6343_eb1a, 0x4017_829c_bc14_e5e1)),
            ("11023", (0x4010_7d63_43eb_1a1f, 0x4017_c14e_5e0a_72f1)),
            ("11020", (0x0000_0000_0000_0000, 0x3ffd_b6db_6db6_db70)),
            ("11024", (0x3faf_58d0_fac6_8800, 0x3ffe_5e0a_72f0_539b)),
            ("11062", (0x3fbf_58d0_fac6_87ec, 0x3fff_0539_7829_cbc4)),
            ("11066", (0x3fc7_829c_bc14_e5ec, 0x3fff_ac68_7d63_43ee)),
            ("11032", (0x4010_29cb_c14e_5e0a, 0x4017_58d0_fac6_87d6)),
            ("11036", (0x4010_3eb1_a1f5_8d0f, 0x4017_9782_9cbc_14e6)),
            ("11021", (0x4010_5397_829c_bc14, 0x4017_d634_3eb1_a1f5)),
            ("11025", (0x3f94_e5e0_a72f_0500, 0x3ffe_b1a1_f58d_0fad)),
            ("11063", (0x3fb4_e5e0_a72f_052c, 0x3fff_58d0_fac6_87d7)),
            ("11060", (0x3ff0_0000_0000_0000, 0x4006_db6d_b6db_6db8)),
            ("11064", (0x3ff0_a72f_0539_782b, 0x4006_b1a1_f58d_0fad)),
            ("11033", (0x4009_7829_cbc1_4e5d, 0x4013_eb1a_1f58_d0fb)),
            ("11030", (0x4010_0000_0000_0000, 0x4017_6db6_db6d_b6dc)),
            ("11034", (0x4010_14e5_e0a7_2f05, 0x4017_ac68_7d63_43ec)),
            ("11002", (0x4010_29cb_c14e_5e0b, 0x4017_eb1a_1f58_d0fc)),
            ("11006", (0x3fa4_e5e0_a72f_05d8, 0x3fff_ac68_7d63_43ef)),
            ("11061", (0x3ff0_5397_829c_bc19, 0x4007_58d0_fac6_87d6)),
            ("11065", (0x3ff0_fac6_87d6_3443, 0x4007_2f05_3978_29cb)),
            ("11433", (0x3ff1_a1f5_8d0f_ac6d, 0x4007_0539_7829_cbc1)),
            ("11031", (0x4008_fac6_87d6_343e, 0x4013_c14e_5e0a_72f0)),
            ("11035", (0x4008_a72f_0539_7829, 0x4013_d634_3eb1_a1f5)),
            ("11003", (0x4008_5397_829c_bc14, 0x4013_eb1a_1f58_d0fa)),
            ("11000", (0x4000_0000_0000_0000, 0x4010_0000_0000_0000)),
            ("11004", (0x3ff0_a72f_0539_782a, 0x4007_d634_3eb1_a1f6)),
            ("11042", (0x3ff1_4e5e_0a72_f054, 0x4007_ac68_7d63_43eb)),
            ("11046", (0x3ff1_f58d_0fac_687d, 0x4007_829c_bc14_e5e1)),
            ("11344", (0x4008_d0fa_c687_d633, 0x4013_829c_bc14_e5e1)),
            ("11012", (0x4008_7d63_43eb_1a1e, 0x4013_9782_9cbc_14e6)),
            ("11016", (0x4008_29cb_c14e_5e09, 0x4013_ac68_7d63_43eb)),
            ("11001", (0x4000_5397_829c_bc12, 0x400f_d634_3eb1_a1f4)),
            ("11043", (0x4000_29cb_c14e_5e08, 0x400f_58d0_fac6_87d5)),
            ("11040", (0x3ff2_4924_9249_2484, 0x4008_0000_0000_0000)),
            ("11044", (0x3ff2_f053_9782_9cb9, 0x4007_d634_3eb1_a1f3)),
            ("11013", (0x4008_5397_829c_bc14, 0x4013_58d0_fac6_87d6)),
            ("11010", (0x4001_2492_4924_924a, 0x4010_0000_0000_0000)),
            ("11014", (0x4000_a72f_0539_7829, 0x400f_ac68_7d63_43ea)),
            ("11041", (0x4000_7d63_43eb_1a1f, 0x400f_2f05_3978_29cb)),
            ("11045", (0x4000_5397_829c_bc14, 0x400e_b1a1_f58d_0fab)),
            ("11011", (0x4001_7829_cbc1_4e5e, 0x400f_d634_3eb1_a1f5)),
            ("11015", (0x4000_fac6_87d6_343f, 0x400f_829c_bc14_e5e0)),
            ("11166", (0x4000_d0fa_c687_d635, 0x400f_0539_7829_cbc1)),
        ];
        let z = zone("110");
        let ours = order(&z, 2);
        assert_eq!(ours.len(), ENGINE.len());
        for (index, (centroid, (text, engine))) in ours.iter().zip(ENGINE).enumerate() {
            assert_eq!(centroid.map(bits), Some(engine), "the centroid at {index}");
            assert_eq!(
                text_of(&address_of_centroid(*centroid, 3)),
                text,
                "the zone at {index}"
            );
            assert_eq!(
                sub_zone_centroid_at_index(&z, 2, index as u64).map(bits),
                Some(engine),
                "the centroid found by the index {index}"
            );
        }
    }

    #[test]
    fn the_entry_found_by_index_is_the_entry_of_the_sequence() {
        // The search skips every scanline but the one that holds the index. At a pentagon, at an
        // odd depth, the engine's own search answers another entry than its sequence holds from
        // the scanline after the pentagon's, since it records the half-width of that scanline
        // only where it generates it; here the entry found is the entry of the sequence, at
        // every index, the pentagons of both poles and of both hemispheres among the zones.
        let mut entries = 0;
        for text in [
            "00", "000", "0000", "01", "010", "0100", "06", "060", "10", "100", "1000", "11",
            "110", "1100", "0064156", "00641565", "013", "0642",
        ] {
            let z = zone(text);
            for depth in 0..=4 {
                let sequence = order(&z, depth);
                assert_eq!(
                    sequence.len() as u64,
                    sub_zones_count(zone_npoints(&z), depth),
                    "{text} at depth {depth}"
                );
                for (index, centroid) in sequence.iter().enumerate() {
                    assert!(centroid.is_some(), "{text} at depth {depth}, index {index}");
                    assert_eq!(
                        sub_zone_centroid_at_index(&z, depth, index as u64).map(bits),
                        centroid.map(bits),
                        "{text} at depth {depth}, index {index}"
                    );
                    entries += 1;
                }
                // An index that the order does not hold finds nothing.
                assert_eq!(
                    sub_zone_centroid_at_index(&z, depth, sequence.len() as u64),
                    None
                );
                assert_eq!(sub_zone_centroid_at_index(&z, depth, u64::MAX), None);
            }
        }
        assert_eq!(
            entries,
            14 * (1 + 11 + 46 + 316 + 2_041) + 4 * (1 + 13 + 55 + 379 + 2_449)
        );
    }

    #[test]
    fn a_centroid_that_no_zone_holds_is_the_null_address() {
        use crate::indexings::Z7;
        use crate::{Indexing, ZoneId};
        // The address of no zone is the one that the Z7 packing writes with every bit set.
        assert_eq!(Z7::encode(&null_address()), ZoneId::NULL);
        assert!(!HexA7::is_valid_address(&null_address()));
        // No centroid: a walk that did not end.
        assert_eq!(address_of_centroid(None, 5), null_address());
        // A centroid in the strip of the seam from base cell 0 to base cell 1 where the
        // quantiser finds no zone at an odd level of 15 or more.
        let z = zone("010004000400");
        let nulls = order(&z, 5)
            .into_iter()
            .filter(|&c| address_of_centroid(c, 15) == null_address())
            .count();
        assert_eq!(nulls, 99);
        // A centroid that a zone holds is that zone's address: the centroid descendant stands in
        // the middle of the order.
        let z = zone("0064156");
        assert_eq!(
            text_of(&address_of_centroid(order(&z, 1)[6], 6)),
            "00641560"
        );
        assert_eq!(
            text_of(&address_of_centroid(order(&z, 2)[27], 7)),
            "006415600"
        );
    }

    /// The address of a Z7 identifier in text.
    fn address(text: &str) -> Address {
        let digits: Vec<u8> = text.bytes().skip(2).map(|b| b - b'0').collect();
        Address::new(text[..2].parse().unwrap(), &digits)
    }

    #[test]
    fn the_four_answers_of_the_topology_agree_with_one_another() {
        // The count is the length of the order; the first sub-zone is its first entry; the
        // zone at an index is its entry there, at every index; and at depth 0 each is the zone
        // itself.
        for text in [
            "00", "000", "01", "010", "06", "060", "10", "100", "11", "110", "0064156", "013",
            "0642",
        ] {
            let a = address(text);
            for depth in 0..=4_u8 {
                let at = format!("{text} at depth {depth}");
                let count = HexA7::count_sub_zones(&a, depth).unwrap();
                let order = HexA7::sub_zones(&a, depth).unwrap();
                assert_eq!(order.len() as u64, count, "{at}");
                assert_eq!(HexA7::first_sub_zone(&a, depth), Some(order[0]), "{at}");
                for (index, entry) in order.iter().enumerate() {
                    assert_eq!(
                        HexA7::sub_zone_at_index(&a, depth, index as u64),
                        Some(*entry),
                        "{at}, index {index}"
                    );
                    assert_eq!(entry.len(), a.len() + usize::from(depth), "{at}");
                }
                assert_eq!(HexA7::sub_zone_at_index(&a, depth, count), None, "{at}");
                assert_eq!(HexA7::sub_zone_at_index(&a, depth, u64::MAX), None, "{at}");
                if depth == 0 {
                    assert_eq!(order, [a], "{at}");
                }
            }
        }
    }

    #[test]
    fn a_depth_beyond_the_finest_level_is_refused_and_never_panics() {
        // Asked directly, the topology answers every address at every depth that a `u8` holds:
        // nothing where the sub-zones would lie beyond level 19, or where the address names no
        // zone. The count beyond depth 24, which is not the engine's, is never reached.
        let addresses = [
            address("00"),
            address("0064156"),
            Address::new(1, &[1; 19]),
            Address::new(11, &[0; 19]),
            Address::new(1, &[1; 20]),
            Address::new(12, &[]),
            Address::new(0, &[7]),
            null_address(),
        ];
        for a in &addresses {
            for depth in 0..=u8::MAX {
                let within = HexA7::is_valid_address(a) && a.len() + usize::from(depth) <= 19;
                let count = HexA7::count_sub_zones(a, depth);
                assert_eq!(count.is_some(), within, "{a:?} at depth {depth}");
                let first = HexA7::first_sub_zone(a, depth);
                assert_eq!(first.is_some(), within, "{a:?} at depth {depth}");
                for index in [0, 1, u64::MAX] {
                    let found = HexA7::sub_zone_at_index(a, depth, index);
                    assert_eq!(
                        found.is_some(),
                        within && index < count.unwrap(),
                        "{a:?} at depth {depth}, index {index}"
                    );
                }
                // The list is built only where it is short.
                if !within || count.unwrap() > crate::grid::max_materialised_sub_zones() {
                    assert_eq!(HexA7::sub_zones(a, depth), None, "{a:?} at depth {depth}");
                } else if depth <= 3 {
                    let order = HexA7::sub_zones(a, depth).unwrap();
                    assert_eq!(order.len() as u64, count.unwrap(), "{a:?} at depth {depth}");
                }
            }
        }
        // The greatest depths that are answered: 19 below level 0, and 0 below level 19.
        assert_eq!(
            HexA7::count_sub_zones(&address("00"), 19),
            Some(9_499_079_489_284_316)
        );
        assert_eq!(HexA7::count_sub_zones(&address("00"), 20), None);
        assert_eq!(
            HexA7::count_sub_zones(&Address::new(1, &[1; 19]), 0),
            Some(1)
        );
        assert_eq!(HexA7::count_sub_zones(&Address::new(1, &[1; 19]), 1), None);
    }

    #[test]
    fn a_list_above_the_limit_is_not_built() {
        // A hexagon at depth 8 and a pentagon at depth 8 have more sub-zones than a list may
        // hold; the first of them and the last are answered all the same.
        for (text, count) in [("0064156", 5_767_201_u64), ("00", 4_806_001)] {
            let a = address(text);
            assert_eq!(HexA7::count_sub_zones(&a, 8), Some(count));
            assert!(count > crate::grid::MAX_MATERIALISED_SUB_ZONES);
            assert_eq!(HexA7::sub_zones(&a, 8), None);
            assert!(HexA7::first_sub_zone(&a, 8).is_some());
            assert!(HexA7::sub_zone_at_index(&a, 8, count - 1).is_some());
        }
        // The topology reads the limit in force, which is the one that the grid reads: a
        // list is given exactly where it is no longer than that limit, whatever the
        // environment of the run has made of it.
        let limit = crate::grid::max_materialised_sub_zones();
        let a = address("0064156");
        for depth in 1..=4 {
            let count = HexA7::count_sub_zones(&a, depth).unwrap();
            assert_eq!(
                HexA7::sub_zones(&a, depth).is_some(),
                count <= limit,
                "depth {depth}"
            );
        }
    }

    #[test]
    fn an_entry_that_is_no_zone_is_the_null_address_and_not_an_absence() {
        // In the seam from base cell 0 to base cell 1: the order of `010004000400` at depth 5
        // holds 99 entries that are no zone. Each is `Some` of the null address, in the list
        // and by its index, so that the grid writes its null identifier there and reports no
        // error.
        let a = address("010004000400");
        let order = HexA7::sub_zones(&a, 5).unwrap();
        assert_eq!(order.len(), 17_053);
        let nulls: Vec<usize> = (0..order.len())
            .filter(|&i| order[i] == null_address())
            .collect();
        assert_eq!(nulls.len(), 99);
        assert_eq!((nulls[0], nulls[98]), (204, 16_893));
        for &index in &nulls {
            assert_eq!(
                HexA7::sub_zone_at_index(&a, 5, index as u64),
                Some(null_address()),
                "index {index}"
            );
        }
    }

    // The tests that follow are of the walk that finds the index of a sub-zone, and of its
    // three sites, VD-1, ZH-1 and FZ-1, each of which is inert and pinned all the same.

    #[test]
    fn the_directions_to_the_vertices_are_the_sources_doubles() {
        // Site VD-1. The six directions at a pentagon of level 0 (an even level) and at a
        // hexagon of level 5 (an odd one), the first three as bits; the last three are their
        // negations.
        for (text, first_three) in [
            (
                "00",
                [
                    (0xbfd5_5555_5555_5555_u64, 0xbfe5_5555_5555_5555_u64),
                    (0xbfe5_5555_5555_5555, 0xbfd5_5555_5555_5555),
                    (0xbfd5_5555_5555_5555, 0x3fd5_5555_5555_5555),
                ],
            ),
            (
                "0064156",
                [
                    (0xbf6f_d831_c1cd_bed1, 0xbf73_e71f_1920_9743),
                    (0xbf73_e71f_1920_9743, 0xbf4f_d831_c1cd_bed1),
                    (0xbf4f_d831_c1cd_bed1, 0x3f6f_d831_c1cd_bed1),
                ],
            ),
        ] {
            let v = get_vertices_directions(&zone(text));
            for i in 0..3 {
                assert_eq!(bits(v[i]), first_three[i], "{text}, direction {i}");
                assert_eq!(
                    bits(v[i + 3]),
                    bits((-v[i].0, -v[i].1)),
                    "{text}, direction {}",
                    i + 3
                );
            }
        }
        // Inert: at every level the library's products with the negated reciprocal are the
        // doubles that the eC's own expressions give, each offset multiplied by the reciprocal
        // and negated where the eC negates it.
        for level in 0..=19_usize {
            let z = zone_from_steps(1, &vec![3; level]);
            let oonp = 1.0 / (7 * pow7(level as i64 / 2)) as f64;
            let source = if level % 2 == 0 {
                let (a, b) = (EVEN_HEX_A, EVEN_HEX_B);
                [
                    (-oonp * a, -oonp * b),
                    (-oonp * b, -oonp * a),
                    (-oonp * a, oonp * a),
                    (oonp * a, oonp * b),
                    (oonp * b, oonp * a),
                    (oonp * a, -oonp * a),
                ]
            } else {
                let (a, b, c) = (ODD_HEX_A, ODD_HEX_B, ODD_HEX_C);
                [
                    (-oonp * a, -oonp * b),
                    (-oonp * b, -oonp * c),
                    (-oonp * c, oonp * a),
                    (oonp * a, oonp * b),
                    (oonp * b, oonp * c),
                    (oonp * c, -oonp * a),
                ]
            };
            assert_eq!(
                get_vertices_directions(&z).map(bits),
                source.map(bits),
                "level {level}"
            );
        }
    }

    #[test]
    fn a_zone_is_a_sub_zone_where_the_engine_says_that_it_is() {
        // Each answer is the engine's `zoneHasSubZone`, the same on the three grids.
        let wrong: Vec<String> = [
            // A hexagon: a sub-zone within it, two across its boundary, its centroid
            // descendant seven levels down, and a sub-zone of its neighbour alone.
            ("0064156", "00641560", true),
            ("0064156", "00641565", true),
            ("0064156", "00641524", true),
            ("0064156", "00641560000000", true),
            ("0064156", "00641545", false),
            // Never a zone of its own level, itself included, nor a coarser one.
            ("0064156", "0064156", false),
            ("0064156", "0064155", false),
            ("0064156", "006415", false),
            // A pentagon at an odd depth: two entries of its order, and the zone that the
            // engine's search by index answers for the first of them.
            ("00", "006", true),
            ("00", "043", true),
            ("00", "030", false),
            // The pentagon of the south pole.
            ("110", "11022", true),
            ("110", "11000", true),
            ("11000", "1100022", true),
            // Nineteen levels down: the last entry of the order of `00`, and the zone that
            // the engine's search by index answers for it.
            ("00", "006620115066201156311", true),
            ("00", "034503320445033204454", false),
            // In the broken seams the test parts from the order: of the entries of four
            // orders at depth 1, some are refused.
            ("010004000400", "01000400040000000", true),
            ("0000000000000005", "00000000000000053", true),
            ("0000000000000005", "00000000000000055", true),
            ("0000000000000005", "00000000000000005", false),
            ("0000000000000000136", "00000000000000003555", true),
            ("0000000000000000136", "00000000000000000005", false),
            ("00052626050026015", "000526260500261626", true),
            ("00052626050026015", "000526260500260533", false),
            ("01000000000000003", "010000000000000006", false),
            ("0000000000000001644", "00000000000000000010", false),
        ]
        .iter()
        .filter(|&&(hay_stack, needle, engine)| {
            zone_has_sub_zone(&zone(hay_stack), &zone(needle)) != engine
        })
        .map(|&(hay_stack, needle, engine)| format!("{needle} in {hay_stack}: {engine}"))
        .collect();
        assert!(wrong.is_empty(), "the engine's answers missed: {wrong:#?}");
    }

    #[test]
    fn two_centroids_are_one_within_1e_11_in_each_coordinate() {
        // Site FZ-1. About a canonical centroid, every pair of differences on either side of
        // the tolerance: the answer is the eC's `!(fabs(dx) < 1E-11 && fabs(dy) < 1E-11)`.
        let zc = zone_centroid(&zone("0064156"));
        assert_eq!(bits(canonicalize5x6(zc.0, zc.1)), bits(zc));
        let steps = [
            0.0, 5e-12, -5e-12, 9.9e-12, -9.9e-12, 1.01e-11, -1.01e-11, 3e-11, -3e-11, 1e-3,
        ];
        let (mut same, mut other) = (0, 0);
        for dx in steps {
            for dy in steps {
                let c = (zc.0 + dx, zc.1 + dy);
                assert_eq!(bits(canonicalize5x6(c.0, c.1)), bits(c));
                let source = !((c.0 - zc.0).abs() < 1e-11 && (c.1 - zc.1).abs() < 1e-11);
                assert_eq!(find_sub_zone(zc, c), source, "{dx:e}, {dy:e}");
                if source {
                    other += 1;
                } else {
                    same += 1;
                }
            }
        }
        assert_eq!((same, other), (25, 75));
        // The one input at which the library parts from the eC's text: an ordinate that is
        // no number, beside an abscissa within the tolerance, is taken for the centroid
        // sought. No centroid of an order is such a point.
        assert!(!find_sub_zone((zc.0, f64::NAN), zc));
        assert!(find_sub_zone((f64::NAN, zc.1), zc));
    }

    #[test]
    fn the_walk_finds_the_index_that_the_engine_finds() {
        // Each answer is the engine's `getSubZoneIndex`, the same on the three grids.
        let wrong: Vec<String> = [
            ("0064156", "00641524", 0_i64),
            ("0064156", "00641565", 5),
            ("0064156", "00641560", 6),
            ("0064156", "00641562", 7),
            ("0064156", "00641545", -1),
            // 0 for two zones of one level, whatever they are; -1 for a coarser one.
            ("0064156", "0064156", 0),
            ("0064156", "0064155", 0),
            ("0064156", "006415", -1),
            // A pentagon at an odd depth, where the engine's search by index is at fault
            // and its walk is not.
            ("00", "006", 9),
            ("00", "043", 10),
            ("00", "030", -1),
            // The pentagon of the south pole, at the depth where the engine's first
            // sub-zone ends the process, and beside it.
            ("110", "11022", 0),
            ("110", "11000", 27),
            ("11000", "1100022", 0),
            ("11000", "1100000", 27),
            ("11", "110", 6),
            // The broken seams: a zone that the order names and the test of a sub-zone
            // refuses, and one that the test accepts and whose centroid no centroid of
            // the order meets, have no index.
            ("010004000400", "01000400040000000", 8_526),
            ("0000000000000000", "00000000000000000", 6),
            ("0000000000000000", "00000000000000001", -1),
            // An index at which the order holds another zone, `000000000000000000`.
            ("00000000000000000", "000000000000000001", 6),
            ("0000000000000005", "00000000000000053", 2),
            ("0000000000000005", "00000000000000052", 3),
            ("0000000000000005", "00000000000000005", -1),
            ("0000000000000005", "00000000000000055", -1),
            ("0000000000000000136", "00000000000000003555", 8),
            ("0000000000000000136", "00000000000000000005", -1),
            ("00052626050026015", "000526260500260105", 1),
            ("00052626050026015", "000526260500260533", -1),
            ("00052626050026015", "000526260500261626", 12),
            // Refused by the test of a sub-zone alone: a centroid of the order lies within
            // the tolerance of this zone's own, and the walk is not begun.
            ("0000000000000001644", "00000000000000000010", -1),
        ]
        .iter()
        .filter_map(|&(parent, sub, engine)| {
            let got = index_by_the_walk(&zone(parent), &zone(sub));
            (got != engine).then(|| format!("{sub} in {parent}: {got}, the engine's {engine}"))
        })
        .collect();
        assert!(wrong.is_empty(), "{wrong:#?}");
    }

    #[test]
    fn the_walk_is_the_inverse_of_the_order_where_the_order_is_sound() {
        // Every entry of an order found again at its own index: hexagons and pentagons of
        // both parities of the level, both poles, at depths 1 to 3.
        let mut entries = 0;
        for text in [
            "00", "10", "11", "000", "110", "0064156", "00641565", "0313522",
        ] {
            let a = address(text);
            let z = zone(text);
            for depth in 1..=3_u8 {
                let order = HexA7::sub_zones(&a, depth).unwrap();
                for (i, sub) in order.iter().enumerate() {
                    assert_eq!(
                        index_by_the_walk(&z, &zone_from_steps(sub.base, sub.digits())),
                        i as i64,
                        "{} in {text} at depth {depth}",
                        text_of(sub)
                    );
                    assert_eq!(
                        HexA7::sub_zone_index(&a, sub),
                        Some(Some(i as u64)),
                        "{} in {text} at depth {depth}",
                        text_of(sub)
                    );
                    entries += 1;
                }
            }
        }
        assert_eq!(entries, 3 * (13 + 55 + 379) + 5 * (11 + 46 + 316));
    }

    #[test]
    fn the_index_of_a_sub_zone_is_refused_where_the_walk_is_not_bounded() {
        let parent = address("0064156");
        // An order of 13, which the walk is made over unless the environment of the run
        // has lowered the limit below it.
        let within = crate::grid::max_materialised_sub_zones() >= 13;
        assert_eq!(
            HexA7::sub_zone_index(&parent, &address("00641565")),
            within.then_some(Some(5))
        );
        assert_eq!(
            HexA7::sub_zone_index(&parent, &address("00641545")),
            within.then_some(None)
        );
        // The walk is made exactly where the order is no longer than the limit in force,
        // which is the one that the grid reads: over the 55 sub-zones two levels down as
        // over the 13.
        assert_eq!(
            HexA7::sub_zone_index(&parent, &address("006415600")),
            (crate::grid::max_materialised_sub_zones() >= 55).then_some(Some(27))
        );
        // Of two zones of one level a zone is its own sub-zone at index 0 and no other is,
        // where the engine's walk answers 0 for both; a coarser zone is none.
        assert_eq!(HexA7::sub_zone_index(&parent, &parent), Some(Some(0)));
        assert_eq!(
            HexA7::sub_zone_index(&parent, &address("0064155")),
            Some(None)
        );
        assert_eq!(
            HexA7::sub_zone_index(&parent, &address("006415")),
            Some(None)
        );
        // An address that names no zone, on either side.
        for bad in [null_address(), Address::new(12, &[]), Address::new(0, &[7])] {
            assert!(!HexA7::is_valid_address(&bad), "{bad:?}");
            assert_eq!(HexA7::sub_zone_index(&parent, &bad), None, "{bad:?}");
            assert_eq!(HexA7::sub_zone_index(&bad, &parent), None, "{bad:?}");
            assert_eq!(HexA7::sub_zone_index(&bad, &bad), None, "{bad:?}");
        }
        // A sub-zone of the twentieth level, which has no geometry.
        assert_eq!(
            HexA7::sub_zone_index(&Address::new(1, &[1; 5]), &Address::new(1, &[1; 20])),
            None
        );
        // An order longer than a list may be: the walk is not begun. Eight levels below a
        // hexagon and a pentagon, and nineteen below a pentagon, where it would be of
        // 9,499,079,489,284,316 centroids.
        for (text, depth) in [("0064156", 8), ("00", 8), ("00", 19)] {
            let a = address(text);
            let mut sub = a;
            for _ in 0..depth {
                sub.push(0);
            }
            assert!(
                HexA7::count_sub_zones(&a, depth).unwrap()
                    > crate::grid::max_materialised_sub_zones()
            );
            assert_eq!(
                HexA7::sub_zone_index(&a, &sub),
                None,
                "{text} at depth {depth}"
            );
        }
    }

    /// The generator's centroids against the engine's, bit for bit and index by index, through
    /// its `getSubZoneCRSCentroids` in the 5x6 CRS, which hands back its generator's output
    /// untouched (`RI7H.ec:615-644`); and, for each, the centroid that the search by index
    /// finds against the centroid of the sequence. The identifiers cannot see the order of the
    /// arithmetic, since a sub-zone centroid is a cell centre, far from any boundary of its
    /// quantisation; the centroids can.
    ///
    /// Sample: every zone of levels 0 to 2 and the twelve pentagons of each level from 3 to 8,
    /// at depths 1 to 5: both poles, both hemispheres, both parities of the level and of the
    /// depth. No pair is left out: at the pentagon of the south pole at an odd level, at depth
    /// 2, where the engine's `getFirstSubZone` ends the process, the call made here is answered
    /// (it was asked at each of the four in a process of its own), and the order of one of them
    /// is pinned as literals besides, by
    /// `the_order_below_the_south_polar_pentagon_at_depth_2_is_the_engines`.
    ///
    /// Gated behind the `oracle` feature, meaningful only within this repository.
    #[test]
    #[cfg(feature = "oracle")]
    fn the_generators_centroids_match_the_engines_bit_for_bit() {
        use crate::Indexing;
        use crate::indexings::Z7;
        use dggal_oracle as o;

        // Every zone of a level, by its text: the addresses that are the zone found at their own
        // centroid, which leaves out the child that a pentagon lacks.
        let mut zones: Vec<Address> = vec![];
        for (level, count) in [(0_u32, 12), (1, 72), (2, 492)] {
            let before = zones.len();
            for base in 0..12 {
                for path in 0..7_u32.pow(level) {
                    let digits: Vec<u8> = (0..level)
                        .map(|k| (path / 7_u32.pow(k) % 7) as u8)
                        .collect();
                    let a = Address::new(base, &digits);
                    if HexA7::quantize(HexA7::planar_centroid(&a), level as u8) == Some(a) {
                        zones.push(a);
                    }
                }
            }
            assert_eq!(zones.len() - before, count, "the zones of level {level}");
        }
        for level in 3..=8 {
            zones.extend((0..12).map(|base| Address::new(base, &vec![0; level])));
        }

        let (mut pairs, mut centroids) = (0_usize, 0_usize);
        for a in &zones {
            let z = zone_from_steps(a.base, a.digits());
            let id = Z7::encode(a).0;
            for depth in 1..=5 {
                let at = format!("{} at depth {depth}", text_of(a));
                let ours = order(&z, depth);
                let theirs = o::sub_zone_crs_centroids_5x6(o::IGEO7, id, depth as i32);
                assert_eq!(ours.len(), theirs.len(), "{at}");
                assert_eq!(
                    ours.len() as u64,
                    sub_zones_count(zone_npoints(&z), depth),
                    "{at}"
                );
                pairs += 1;
                for (index, (&centroid, &engine)) in ours.iter().zip(&theirs).enumerate() {
                    assert_eq!(
                        centroid.map(bits),
                        Some(bits(engine)),
                        "{at}: the centroid at {index}"
                    );
                    assert_eq!(
                        sub_zone_centroid_at_index(&z, depth, index as u64).map(bits),
                        centroid.map(bits),
                        "{at}: the centroid found by the index {index}"
                    );
                    centroids += 1;
                }
            }
        }
        // Measured, and asserted exactly, so that the test cannot quietly compare less than it
        // claims: 648 zones at five depths, 540 of them hexagons and 108 pentagons.
        assert_eq!(zones.len(), 12 + 72 + 492 + 6 * 12);
        assert_eq!(pairs, 648 * 5);
        assert_eq!(
            centroids,
            540 * (13 + 55 + 379 + 2_449 + 17_053) + 108 * (11 + 46 + 316 + 2_041 + 14_211)
        );
    }

    /// The test of a sub-zone and the walk against the engine's `zoneHasSubZone` and
    /// `getSubZoneIndex`, answer for answer, over pairs of both kinds: every entry of an order
    /// against its zone, and every entry of the order of the zone before it in the sample,
    /// most of which are no sub-zones of this one and some of which, across a shared
    /// boundary, are.
    ///
    /// Sample: every zone of levels 0 to 2 at depths 1 and 2, and those of levels 0 and 1 at
    /// depth 3; the twelve pentagons of each level from 3 to 17 at depths 1 and 2, the deepest
    /// of which have their sub-zones in the broken seams; six zones of a broken seam at depth
    /// 1; and one entry in 64 of the 17,053 of a seam order at depth 5. An entry that is no
    /// zone is not asked.
    ///
    /// What the tallies record besides, all of it in the seams: the entries of an order that
    /// have no index, and the indices that the walk answers, as the engine does, at which the
    /// order holds another zone, which is why the grid confirms an index before it answers
    /// it.
    ///
    /// Gated behind the `oracle` feature, meaningful only within this repository.
    #[test]
    #[cfg(feature = "oracle")]
    fn the_test_of_a_sub_zone_and_the_walk_match_the_engines() {
        use crate::Indexing;
        use crate::indexings::Z7;
        use dggal_oracle as o;

        // (zone, depth, one entry in this many).
        let mut sample: Vec<(Address, u8, usize)> = vec![];
        for level in 0..=2_u32 {
            for base in 0..12 {
                for path in 0..7_u32.pow(level) {
                    let digits: Vec<u8> = (0..level)
                        .map(|k| (path / 7_u32.pow(k) % 7) as u8)
                        .collect();
                    let a = Address::new(base, &digits);
                    if HexA7::quantize(HexA7::planar_centroid(&a), level as u8) == Some(a) {
                        let deepest = if level < 2 { 3 } else { 2 };
                        sample.extend((1..=deepest).map(|depth| (a, depth, 1)));
                    }
                }
            }
        }
        for level in 3..=17 {
            for base in 0..12 {
                let a = Address::new(base, &vec![0; level]);
                sample.extend([(a, 1, 1), (a, 2, 1)]);
            }
        }
        for text in [
            "0000000000000000",
            "0000000000000005",
            "0000000000000000136",
            "0000000000000001644",
            "00052626050026015",
            "01000000000000003",
        ] {
            sample.push((address(text), 1, 1));
        }
        sample.push((address("010004000400"), 5, 64));

        // [pairs, accepted by the test, with an index, entries of an order without one,
        // indices at which the order holds another zone].
        let mut tally = [0_usize; 5];
        // By depth: the level of the zone asked last at that depth, and its order.
        let mut before: [Option<(usize, Vec<Address>)>; 6] = Default::default();
        for &(a, depth, stride) in &sample {
            let z = zone_from_steps(a.base, a.digits());
            let id = Z7::encode(&a).0;
            let order = HexA7::sub_zones(&a, depth).unwrap();
            let own = order.iter().step_by(stride).map(|sub| (sub, true));
            // The entries of the order of the zone before this one, where it is of this level.
            let foreign = before[usize::from(depth)]
                .iter()
                .filter(|(level, _)| *level == a.len())
                .flat_map(|(_, entries)| entries.iter().step_by(stride))
                .map(|sub| (sub, false));
            for (sub, entry) in own.chain(foreign) {
                if *sub == null_address() {
                    continue;
                }
                let at = format!("{} in {} at depth {depth}", text_of(sub), text_of(&a));
                let sz = zone_from_steps(sub.base, sub.digits());
                let sub_id = Z7::encode(sub).0;
                let has = zone_has_sub_zone(&z, &sz);
                assert_eq!(has, o::zone_has_sub_zone(o::IGEO7, id, sub_id), "{at}");
                let index = index_by_the_walk(&z, &sz);
                assert_eq!(index, o::sub_zone_index(o::IGEO7, id, sub_id), "{at}");
                tally[0] += 1;
                tally[1] += usize::from(has);
                tally[2] += usize::from(index >= 0);
                tally[3] += usize::from(entry && index < 0);
                if index >= 0 {
                    assert!(has, "{at}: an index for a zone that the test refuses");
                    tally[4] += usize::from(order[index as usize] != *sub);
                }
            }
            before[usize::from(depth)] = Some((a.len(), order));
        }
        // Measured, and asserted exactly, so that the test cannot quietly compare less than it
        // claims. The 98 entries without an index and the 146 indices at which the order holds
        // another zone all have their sub-zones at level 15 or finer.
        assert_eq!(tally, [157_016, 80_278, 80_240, 98, 146]);
    }
}
