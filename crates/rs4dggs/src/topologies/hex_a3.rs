//! Aperture-3 hexagonal topology on the 5x6 planar grid: the topology half of ISEA3H, IVEA3H
//! and RTEA3H.
//!
//! It turns a projection's planar point `(x, y)` in DGGAL's 5x6 layout into an I3H zone,
//! rebuilds a zone's planar centroid and vertices, and lists its neighbours, parents and
//! children; the sibling module `hex_a3_subzones` orders its sub-zones. The zone travels whole
//! in [`Address::base`], packed as the eC packs an `I3HZone` (`RI3H.ec:855-861`), since the
//! aperture-3 hierarchy is not a digit path; [`crate::indexings::I3h`] reads and prints it.
//!
//! Provenance: `dggrs/RI3H.ec` of DGGAL v0.0.6: `I3HZone::fromI9R` (`:951-963`), the property
//! `parent0` (`:977-1003`), `getNeighbor` (`:1005-1161`), `getNeighbors` (`:1163-1201`), the
//! property `centroidParent` (`:1203-1223`), `getParents` (`:1247-1328`), `fromCentroid`
//! (`:1330-1672`), `getVertices` (`:1674-1790`), the property `centroidChild` (`:2012-2042`),
//! `getChildren` (`:2044-2085`) and the properties `centroid` (`:2138-2168`) and
//! `isCentroidChild` (`:2170-2196`); `iLRCFromLRtI` of `dggrs/RI9R.ec` (`:710-729`); and the 5x6
//! helpers of `projections/ri5x6.ec` from [`crate::fivebysix`]. None of these members carries
//! `__attribute__ ((optimize("-fno-unsafe-math-optimizations")))`, and DGGAL is built with
//! `-O2 -ffast-math` (`Makefile.dggal:133`), so gcc was free to re-associate their arithmetic,
//! and it did. At every such site on this path this module performs the order read
//! out of the shipped `libdggal.so`, BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, and a
//! `// faithful:` comment names the address and the eC line it replaces. Eight of those sites can
//! change an answer (a zone, or a coordinate compared bit for bit); eight cannot, and follow the
//! library all the same, so that the port rests on no argument about its inputs. Everywhere
//! else the eC's text decides, integer arithmetic included: `uint64` arithmetic wraps as the
//! eC's does, an `int` that gcc computes in 32 bits is truncated to them, and a conversion from
//! a double to an integer is performed as gcc emitted it.
//!
//! The order of a zone's vertices departs from the engine's on purpose: the engine's ring runs
//! clockwise as seen from outside the sphere, save at the two polar pentagons, and
//! [`HexA3::planar_vertices`] turns it about its first vertex, so that every ring runs
//! anticlockwise, as on the aperture-7 grids. The coordinates are the engine's own, and the
//! sub-zones are generated from the engine's own order.
use super::hex_a3_subzones;
use crate::fivebysix::{
    canonicalize5x6, cvtt_i32, floor_i32, move5x6_vertex, move5x6_vertex2, refine5x6,
};
use crate::indexings::{fields, is_readable, pack};
use crate::{Address, PlanarPoint, Topology};

/// DGGAL's aperture-3 hexagonal topology: the I9R rhombic grid of every even level, refined at
/// the odd levels into the sub-hexagons B, C and D (`RI3H.ec:855-861`).
///
/// A zone's vertices run anticlockwise as seen from outside the sphere, which departs from the
/// engine's order everywhere but at the two polar pentagons; see [`HexA3::planar_vertices`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct HexA3;

impl crate::interfaces::sealed::Sealed for HexA3 {}

/// `fl(1/3)`, which gcc substituted for divisions by three; `.rodata` at `0x49a98`.
const T13: f64 = f64::from_bits(0x3fd5_5555_5555_5555);
/// `fl(2/3)`, substituted likewise; `.rodata` at `0x49a90`.
const T23: f64 = f64::from_bits(0x3fe5_5555_5555_5555);
/// `fl(-1/3)`, substituted for the negated division in `getVertices`; `.rodata` at `0x49f60`.
const MINUS_T13: f64 = f64::from_bits(0xbfd5_5555_5555_5555);
/// `fl(10 - 1e-11)`, folded into the corner test of `fromCentroid`; `.rodata` at `0x49f58`.
const TEN_MINUS_1E11: f64 = 10.0 - 1e-11;
/// `fl(5 - 1e-11)`, folded into the corner wrap of `fromCentroid`; `.rodata` at `0x49f48`.
const FIVE_MINUS_1E11: f64 = 5.0 - 1e-11;
/// `fl(1 + 1e-11)`, folded into the southern-half test of `getNeighbor` and into a crossing test
/// of `getParents`; `.rodata` at `0x4a010`.
const ONE_PLUS_1E11: f64 = 1.0 + 1e-11;
/// `fl(1 - 1e-11)`, folded into the other crossing test of `getParents`; `.rodata` at `0x4a018`.
const ONE_MINUS_1E11: f64 = 1.0 - 1e-11;
/// `2^63`, the threshold of gcc's conversion of a double to `uint64`; `.rodata` at `0x49ef0`.
const TWO_POW_63: f64 = 9_223_372_036_854_775_808.0;
/// The finest level `getZoneFromCRSCentroid` quantises at (`RI3H.ec:96`).
pub(super) const MAX_LEVEL: u8 = 33;

/// `cvttsd2si` into a 64-bit register, as [`cvtt_i32`] into a 32-bit one.
///
/// No eC counterpart: a machine conversion, which gcc emits within the eC's `(uint64)` casts
/// and on which [`to_uint64`] is built.
fn cvtt_i64(v: f64) -> i64 {
    if (-TWO_POW_63..TWO_POW_63).contains(&v) {
        v as i64
    } else {
        i64::MIN
    }
}

/// The eC's `(uint64)` of a double, as gcc emitted it at `0xc8c2` to `0xc90f`: below `2^63` a
/// signed truncation reinterpreted, and from `2^63` the value less `2^63`, truncated, with bit 63
/// flipped back. From `2^64` that gives 0, where Rust's `as u64` saturates, and the two answer
/// different zones for a point far outside the layout. Below `-1` it wraps where Rust's gives 0,
/// which no answer shows: the wrapped column is out of range and the saturated one leaves an
/// offset below `-1`, and `fromCentroid` answers the null zone for either. The casts it serves
/// are those of `fromCentroid` (`RI3H.ec:1367-1368`).
fn to_uint64(v: f64) -> u64 {
    if v >= TWO_POW_63 {
        (cvtt_i64(v - TWO_POW_63) as u64) ^ (1 << 63)
    } else {
        cvtt_i64(v) as u64
    }
}

/// The eC's `POW3` (`RI3H.ec:29`) for the levels an identifier's five-bit field can hold, all
/// within `u64`.
pub(crate) fn pow3(level: u64) -> u64 {
    3u64.pow(level as u32)
}

/// The row and column of a zone's I9R cell in the whole layout, from its root rhombus and its
/// index within it, with `p = 3^levelI9R`, as the `centroid` property, `getVertices`,
/// `iLRCFromLRtI` and `centroidChild` all form them (`RI3H.ec:2153-2156`, `:1680-1683`,
/// `RI9R.ec:719-722`, `RI3H.ec:2033-2034`). The eC's `uint64` arithmetic wraps: `colOP - ixOP`
/// is negative for every cell below the first row of a southern rhombus, and the product wraps
/// back. Where the eC holds the two in an `int`, the caller truncates them, which keeps the low
/// 32 bits that gcc's 32-bit arithmetic computes there.
fn row_col(p: u64, root: u64, rix: u64) -> (u64, u64) {
    let row_op = (root + 1) >> 1;
    let col_op = root >> 1;
    let ix_op = rix / p;
    (
        row_op.wrapping_mul(p).wrapping_add(ix_op),
        col_op.wrapping_sub(ix_op).wrapping_mul(p).wrapping_add(rix),
    )
}

/// The top-left corner of a zone's I9R cell and `d = 1 / p`, as the `centroid` property and
/// `getVertices` both compute them (`RI3H.ec:2150-2158`, `:1679-1685`). The polar rows and
/// columns are `getVertices`'s; the property returns the poles before it reaches them.
fn top_left(level_i9r: u64, root: u64, rix: u64) -> ((f64, f64), f64) {
    let p = pow3(level_i9r);
    let (row, col) = match root {
        10 => (0, p - 1),
        11 => (6 * p - 1, 4 * p),
        _ => row_col(p, root, rix),
    };
    let d = 1.0 / p as f64;
    ((col as f64 * d, row as f64 * d), d)
}

/// The zone at `level` whose cell contains the planar point `(cx, cy)`, packed, or `None` for the
/// eC's null zone. Ports `I3HZone::fromCentroid` (`RI3H.ec:1330-1672`), compiled at `0xc760` to
/// `0xd4e4`; the eC's local names are kept.
#[allow(
    clippy::if_same_then_else,
    clippy::collapsible_else_if,
    reason = "the eC's branches are kept apart, as written, for the line-by-line correspondence"
)]
fn from_centroid(level: u32, cx_in: f64, cy_in: f64) -> Option<u64> {
    let (mut x, mut y) = (cx_in, cy_in);
    let l9r = level / 2;
    let p = pow3(u64::from(l9r));
    let pf = p as f64;
    let d = 1.0 / pf;
    let mut is_north_pole = false;
    let mut is_south_pole = false;
    // faithful: the north-pole test is compiled in the eC's order (`0xc7f7` to `0xc81a`), but
    // the south-pole test `fabs(c.y - c.x - 2)` (RI3H.ec:1339) is compiled as
    // `fabs((c.y - 2) - c.x)`, in DGGAL's `-O2 -ffast-math` build, `libdggal.so`, BuildID
    // `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, at `0xc820` to `0xc838`. Within a unit of the
    // threshold the two groupings decide differently whether the point is the south pole.
    if ((x - y) - 1.0).abs() < 1e-10 {
        is_north_pole = true;
    } else if ((y - 2.0) - x).abs() < 1e-10 {
        is_south_pole = true;
    } else if y < -1e-11 && x > -1e-11 {
        x -= y;
        y = 0.0;
    } else if floor_i32(x + 1e-11) > floor_i32(y + 1e-11) {
        // Over top dent to the right
        let cy = 5.min(floor_i32(y + 1e-11));
        // faithful: `c.x += (cy+1 - c.y)` (RI3H.ec:1351) is compiled as `(c.x - c.y) + (cy + 1)`,
        // reusing the `c.x - c.y` of the pole test, in the same library at `0xcad1` to `0xcaf0`.
        x = (x - y) + f64::from(cy.wrapping_add(1));
        y = f64::from(cy.wrapping_add(1));
    } else if floor_i32(y + 1e-11).wrapping_sub(floor_i32(x + 1e-11)) > 1 {
        // Over bottom dent to the right
        let cx = 4.min(floor_i32(x + 1e-11));
        // faithful: `c.y += (cx+1 - c.x)` (RI3H.ec:1357) is compiled as `(c.y - c.x) + (cx + 1)`,
        // in the same library at `0xcecb` to `0xcee9`.
        y = (y - x) + f64::from(cx.wrapping_add(1));
        x = f64::from(cx.wrapping_add(1));
    } else if x < -1e-11 || y < -1e-11 {
        (x, y) = move5x6_vertex(5.0, 5.0, x, y);
    }
    // faithful: the corner test `c.x + c.y > 5.0 + 5.0 - d - 1E-11` (RI3H.ec:1361-1362) is
    // compiled as `c.x + c.y > fl(10 - 1e-11) - d`, in the same library at `0xc862` to `0xc883`.
    // It cannot change an answer: the first two conjuncts force `c.x + c.y > 10 - 2e-11`, while
    // `d` is at least `3^-16`, so neither threshold rises above `10 - 2.3e-8`.
    //
    // faithful: on the wrap, the ordinate the truncation is given, `(c.y - 5) + 1e-11`
    // (RI3H.ec:1363, :1365), is compiled as `c.y - fl(5 - 1e-11)` of the unwrapped `c.y`, at
    // `0xcdb0` to `0xcdd4`; the abscissa keeps `(c.x - 5) + 1e-11`. The two forms flip at the same
    // double at every crossing the `Min(5, .)` below admits, so this cannot change an answer.
    let mut wrapped_y = None;
    if x > 5.0 - 1e-11 && y > 5.0 - 1e-11 && x + y > TEN_MINUS_1E11 - d {
        wrapped_y = Some(y - FIVE_MINUS_1E11);
        x -= 5.0;
        y -= 5.0;
    }
    // The coordinate of the root rhombus, then the row and column within it.
    let mut cx = 4.min(cvtt_i32(x + 1e-11));
    let mut cy = 5.min(cvtt_i32(wrapped_y.unwrap_or(y + 1e-11)));
    let mut root = cx.wrapping_add(cy) as u32;
    let mut xi = to_uint64((x - f64::from(cx)) * pf + 1e-6);
    let mut yi = to_uint64((y - f64::from(cy)) * pf + 1e-6);
    // The re-anchoring branches; with `Some(true)` the abscissa was re-anchored, with
    // `Some(false)` the ordinate.
    let mut reanchored = None;

    // REVIEW: Valid scenarios where x == p or y == p are currently possible (the eC's note).
    if yi == p {
        cy = cy.wrapping_add(1);
        root = root.wrapping_add(1);
        yi = yi.wrapping_sub(p);
        y = f64::from(cy) + yi as f64 / pf;
    }
    if xi == p {
        cx = cx.wrapping_add(1);
        root = root.wrapping_add(1);
        if root == 10 {
            cx = cx.wrapping_sub(5);
            cy = cy.wrapping_sub(5);
            root = 0;
            y = f64::from(cy) + yi as f64 / pf;
        }
        xi = xi.wrapping_sub(p);
        x = f64::from(cx) + xi as f64 / pf;
    }
    if cy.wrapping_sub(cx) > 1 && yi == 0 {
        cx = cx.wrapping_add(1);
        yi = p.wrapping_sub(xi);
        xi = 0;
        root = root.wrapping_add(1);
        y = f64::from(cy) + yi as f64 / pf;
        x = f64::from(cx);
        reanchored = Some(false);
    } else if cy < cx && xi == 0 {
        cy = cy.wrapping_add(1);
        xi = p.wrapping_sub(yi);
        yi = 0;
        root = root.wrapping_add(1);
        x = f64::from(cx) + xi as f64 / pf;
        y = f64::from(cy);
        reanchored = Some(true);
    }

    // Index within root rhombus
    let mut rix = yi.wrapping_mul(p).wrapping_add(xi);
    // faithful: `xd = (c.x - cx) * p - x` and `yd = (c.y - cy) * p - y` (RI3H.ec:1416-1417) are
    // kept apart as `X - x` and `Y - y`, since the diagonals below read `X` and `Y` themselves.
    // On the two re-anchoring branches (RI3H.ec:1396-1413) gcc cancelled `cy` (resp. `cx`) and
    // forms the re-anchored offset as `fl(y / p) * p` (resp. `fl(x / p) * p`), at `0xcc1b` to
    // `0xcc5f` and `0xc953` to `0xc9ac`; the other offset is zero either way. The two values
    // part at dozens of cell corners and decide every comparison below alike, so this cannot
    // change an answer.
    let big_x = match reanchored {
        Some(true) => (xi as f64 / pf) * pf,
        _ => (x - f64::from(cx)) * pf,
    };
    let big_y = match reanchored {
        Some(false) => (yi as f64 / pf) * pf,
        _ => (y - f64::from(cy)) * pf,
    };
    let (xf, yf) = (xi as f64, yi as f64);
    let xd = big_x - xf;
    let yd = big_y - yf;

    let sh: u64;
    if is_north_pole {
        sh = u64::from(level & 1);
        root = 10;
        rix = 0;
    } else if is_south_pole {
        sh = u64::from(level & 1);
        root = 11;
        rix = 0;
    } else {
        let right_sr = xi == p - 1;
        let top_sr = yi == 0;
        let left_sr = xi == 0;
        let bottom_sr = yi == p - 1;
        let np_sub_rhombus = right_sr && top_sr && root & 1 == 0;
        let sp_sub_rhombus = bottom_sr && left_sr && root & 1 != 0;

        if cy < cx || xd < -1.0 || yd < -1.0 || xi >= p || yi >= p || rix >= p.wrapping_mul(p) {
            return None; // y cannot be smaller than x
        }

        // From here `xi` and `yi` are below `p`, so `p - 1 - xi` and its like cannot wrap. Beyond
        // the I9R level 16 of any cell the products with `p` can: from level 21 `p * p` exceeds
        // 64 bits, and `rix` is only below its wrapped value. They wrap as the eC's `uint64`
        // does; the sums with `rix` wrap likewise, although that bound keeps them in range at
        // every level the field can hold. The root can wrap too, for a point far outside the
        // layout.
        if level & 1 != 0 {
            // Odd level
            //
            // faithful: none of the products `3*xd`, `3*yd` (RI3H.ec:1436, :1442) is formed:
            // `3*xd < 1` is compiled as `xd < fl(1/3)` and `3*xd > 2` as `xd > fl(2/3)`, and the
            // same for `yd`, in the same library at `0xcf70`, `0xcf7e`, `0xd018`, `0xd026`,
            // `0xd138` and `0xd1b8`. The first pair decides alike at every double; the second
            // does not at `next(fl(2/3))`, where `3*xd` rounds to 2.
            let left_third = xd < T13;
            let top_third = yd < T13;

            if left_third && top_third {
                sh = 1; // B
            } else {
                let right_third = xd > T23;
                let bottom_third = yd > T23;
                if right_third && bottom_third {
                    if bottom_sr && right_sr {
                        // Non-polar pentagon
                        root = root.wrapping_add(2) % 10;
                        rix = 0;
                    } else if bottom_sr {
                        // Indexed to another root rhombus
                        if root & 1 != 0 {
                            // Crossing South interruption to the right
                            root = root.wrapping_add(2) % 10;
                            rix = (p - 1 - xi).wrapping_mul(p);
                        } else {
                            root = root.wrapping_add(1);
                            rix = xi + 1;
                        }
                    } else if right_sr {
                        // Indexed to another root rhombus
                        if root & 1 == 0 {
                            // Crossing North interruption to the right
                            root = root.wrapping_add(2) % 10;
                            rix = p - 1 - yi;
                        } else {
                            root = root.wrapping_add(1) % 10;
                            rix = p.wrapping_mul(yi + 1);
                        }
                    } else {
                        rix = rix.wrapping_add(p + 1);
                    }
                    sh = 1; // B
                } else if bottom_third {
                    // faithful: `3 * (yd - xd) > 2` (RI3H.ec:1487) is compiled as
                    // `(x + yd) - X > fl(2/3)`, at `0xd13e` to `0xd14a`.
                    if (xf + yd) - big_x > T23 {
                        if sp_sub_rhombus {
                            // "South" pole B
                            root = 11;
                            rix = 0;
                            sh = 1;
                        } else {
                            if bottom_sr {
                                // Indexed to another root rhombus
                                if root & 1 != 0 {
                                    // Crossing South interruption to the right
                                    root = root.wrapping_add(2) % 10;
                                    rix = (p - xi).wrapping_mul(p);
                                } else {
                                    rix = xi;
                                    root = root.wrapping_add(1);
                                }
                            } else {
                                rix = rix.wrapping_add(p);
                            }
                            sh = 1; // B
                        }
                    } else {
                        sh = 3; // D
                    }
                } else if right_third {
                    // faithful: `3 * (xd - yd) > 2` (RI3H.ec:1518) is compiled as
                    // `(y + xd) - Y > fl(2/3)`, at `0xd1be` to `0xd1ca`.
                    if (yf + xd) - big_y > T23 {
                        if np_sub_rhombus {
                            // "North" pole B
                            root = 10;
                            rix = 0;
                            sh = 1;
                        } else {
                            if right_sr {
                                if root & 1 == 0 {
                                    // Crossing North interruption to the right
                                    root = root.wrapping_add(2) % 10;
                                    rix = p - yi;
                                } else {
                                    root = root.wrapping_add(1) % 10;
                                    rix = p.wrapping_mul(yi);
                                }
                            } else {
                                rix = rix.wrapping_add(1);
                            }
                            sh = 1; // B
                        }
                    } else {
                        sh = 2; // C
                    }
                } else if xd > yd {
                    sh = 2; // C
                } else {
                    sh = 3; // D
                }
            }
        } else {
            // Even level
            let (mut top_right, mut bottom_left, mut bottom_right) = (false, false, false);
            if xd - 1.0 > -yd {
                // Bottom-Right portion
                if xd > yd * 2.0 {
                    top_right = true;
                } else if 2.0 * xd < yd {
                    bottom_left = true;
                } else {
                    bottom_right = true;
                }
            } else {
                // Top-Left portion
                if 2.0 * xd > yd + 1.0 {
                    top_right = true;
                } else if xd + 1.0 < 2.0 * yd {
                    bottom_left = true;
                }
            }

            sh = 0; // A
            if top_right {
                if np_sub_rhombus {
                    // "North" pole A
                    root = 10;
                    rix = 0;
                } else if right_sr {
                    if root & 1 == 0 {
                        // Crossing North interruption to the right
                        root = root.wrapping_add(2) % 10;
                        rix = p - yi;
                    } else {
                        root = root.wrapping_add(1) % 10;
                        rix = p.wrapping_mul(yi);
                    }
                } else {
                    rix = rix.wrapping_add(1);
                }
            } else if bottom_left {
                if sp_sub_rhombus {
                    // "South" pole A
                    root = 11;
                    rix = 0;
                } else if bottom_sr {
                    // Indexed to another root rhombus
                    if root & 1 != 0 {
                        // Crossing South interruption to the right
                        root = root.wrapping_add(2) % 10;
                        rix = (p - xi).wrapping_mul(p);
                    } else {
                        rix = xi;
                        root = root.wrapping_add(1);
                    }
                } else {
                    rix = rix.wrapping_add(p);
                }
            } else if bottom_right {
                if bottom_sr && right_sr {
                    // Non-polar pentagon
                    root = root.wrapping_add(2) % 10;
                    rix = 0;
                } else if bottom_sr {
                    // Indexed to another root rhombus
                    if root & 1 != 0 {
                        // Crossing South interruption to the right
                        root = root.wrapping_add(2) % 10;
                        rix = (p - 1 - xi).wrapping_mul(p);
                    } else {
                        root = root.wrapping_add(1);
                        rix = xi + 1;
                    }
                } else if right_sr {
                    // Indexed to another root rhombus
                    if root & 1 == 0 {
                        // Crossing North interruption to the right
                        root = root.wrapping_add(2) % 10;
                        rix = p - 1 - yi;
                    } else {
                        root = root.wrapping_add(1) % 10;
                        rix = p.wrapping_mul(yi + 1);
                    }
                } else {
                    rix = rix.wrapping_add(p + 1);
                }
            }
        }
    }
    Some(pack(u64::from(l9r), u64::from(root), rix, sh))
}

/// The planar centroid of a zone. Ports the `I3HZone` `centroid` property (`RI3H.ec:2138-2168`),
/// compiled at `0xa9c0` to `0xaaee`.
fn centroid(id: u64) -> (f64, f64) {
    let (level_i9r, root, rix, sh) = fields(id);
    match root {
        10 => return (1.0, 0.0), // "North" pole (Even level A and Odd level B)
        11 => return (4.0, 6.0), // "South" pole (Even level A and Odd level B)
        _ => {}
    }
    let ((tx, ty), d) = top_left(level_i9r, root, rix);
    // faithful: `2*d/3` and `d/3` (RI3H.ec:2163, :2165) are compiled as `d * fl(2/3)` and
    // `d * fl(1/3)`, both lanes in one `mulpd`, in DGGAL's `-O2 -ffast-math` build,
    // `libdggal.so`, BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, at `0xaabf` (sub-hex D)
    // and `0xaae0` (sub-hex C). The two forms differ by a unit in the last place at several
    // levels, which the sum with the corner rarely but measurably carries into the centroid.
    match sh {
        0 | 1 => (tx, ty),                 // Even level A or Odd level B hex
        2 => (tx + d * T23, ty + d * T13), // Odd level C hex
        _ => (tx + d * T13, ty + d * T23), // Odd level D hex
    }
}

/// The planar vertices of a zone, five or six, before any canonicalisation. Ports
/// `I3HZone::getVertices` (`RI3H.ec:1674-1790`), compiled at `0xd8a0` to `0xe03f`.
pub(super) fn vertices(id: u64) -> Vec<(f64, f64)> {
    vertices_from(id, false)
}

/// The vertices from which the sub-zones of zone `id` are generated: [`vertices`], except for a
/// hexagon one of whose vertices [`vertices`] carries across an interruption that the vertex does
/// not reach, which takes them from its exact corner.
///
/// A departure from the engine, which generates the sub-zones of such a zone from a point
/// outside it. A corner on a rhombus edge is a whole number, but `col * fl(1/p)` falls a unit in
/// the last place short of it at I9R levels 6, 9, 10 and 16 for the edges at 1, 2 and 4. A vertex
/// that does not move along that axis then floors into the column (or row) before it, and
/// `move5x6Vertex` gives its image across the interruption, which `getI3HFirstSubZoneCentroid`
/// compares with the others as though it lay in the zone's frame. The scanlines that start from
/// it leave the zone, and the engine's order names zones twice, names zones that are not
/// descendants and omits descendants. From the exact corner no vertex is carried, and the order
/// is the zone's descendants, in the scanline order the generator gives elsewhere. Every other
/// zone keeps the engine's vertices, and with them its order, down to the last bit.
///
/// It stands where `getI3HFirstSubZoneCentroid` calls `getVertices` (`I3HSubZones.ec:105`).
pub(super) fn sub_zone_vertices(id: u64) -> Vec<(f64, f64)> {
    let faithful = vertices(id);
    if faithful.len() != 6 {
        return faithful;
    }
    let exact = vertices_from(id, true);
    // An image across an interruption lies at least a third of a cell (1e-8 at level 33) from its
    // vertex; the two corners differ by a unit in the last place.
    let carried = faithful
        .iter()
        .zip(&exact)
        .any(|(a, b)| (a.0 - b.0).abs() > 1e-11 || (a.1 - b.1).abs() > 1e-11);
    if carried { exact } else { faithful }
}

/// [`vertices`], from the exact corner where `exact_corner` is set: a corner within 1e-11 below a
/// whole number is that number, which only a corner on a rhombus edge can be, since an interior
/// one lies at least `1/p >= 3^-16` from one.
fn vertices_from(id: u64, exact_corner: bool) -> Vec<(f64, f64)> {
    let (level_i9r, root, rix, sh) = fields(id);
    let ((mut tx, mut ty), d) = top_left(level_i9r, root, rix);
    if exact_corner {
        let exact = |v: f64| if v.ceil() - v < 1e-11 { v.ceil() } else { v };
        (tx, ty) = (exact(tx), exact(ty));
    }
    let south = root & 1 != 0;
    let mv = |dx: f64, dy: f64| move5x6_vertex(tx, ty, dx, dy);
    // faithful: every `d/3` of `getVertices` (RI3H.ec:1699-1786) is compiled as `d * fl(1/3)`
    // and every `2*d/3` as `d * fl(2/3)`; `-d/3` is `d * fl(-1/3)` in the even arm (`0xd9e3`) and
    // `-(d * fl(1/3))` in the odd one (`0xdbca`), the same double, and `-2*d/3` is
    // `-(d * fl(2/3))` (`0xda83`). In DGGAL's `-O2 -ffast-math` build, `libdggal.so`, BuildID
    // `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, at `0xd9df`, `0xda04`, `0xdb41`, `0xdc38`,
    // `0xdc80`, `0xdd10`, `0xdd46` and `0xdf9d`. The guarded `move5x6Vertex` then adds each
    // offset to the corner in the eC's order, and the polar fans add their whole numbers as the
    // eC writes them (`0xdec7` to `0xe018`).
    let d3 = d * T13;
    let d23 = d * T23;
    let mut v = Vec::with_capacity(6);
    match sh {
        // Even level
        0 => {
            let m3 = d * MINUS_T13;
            if root == 10 {
                // "North" pole
                let (x, y) = mv(d3, m3);
                v.extend((1..=5).map(|k| (x + f64::from(k), y + f64::from(k))));
            } else if root == 11 {
                // "South" pole
                let (x, y) = mv(m3, d3);
                v.extend([
                    (x - 0.0, y - 0.0),
                    (x - 1.0, y - 1.0),
                    (x - 2.0, y - 2.0),
                    (x - 3.0, y - 3.0),
                    (x + 1.0, y + 1.0),
                ]);
            } else {
                // Regular case
                v.push(mv(d23, d3));
                v.push(mv(d3, d23));
                if !south || rix != 0 {
                    // 0 rhombusIndex are pentagons
                    v.push(mv(m3, d3));
                }
                v.push(mv(-d23, m3));
                if south || rix != 0 {
                    v.push(mv(m3, -d23));
                }
                v.push(mv(d3, m3));
            }
        }
        // Odd level, type B
        1 => {
            if root == 10 {
                // "North" pole
                let (x, y) = mv(d23, 0.0);
                v.extend((1..=5).map(|k| (x + f64::from(k), y + f64::from(k))));
            } else if root == 11 {
                // "South" pole
                let (x, y) = mv(d3, d);
                v.extend((0..=4).map(|k| (x - f64::from(k), y - f64::from(k))));
            } else {
                if south || rix != 0 {
                    // 0 rhombusIndex are pentagons
                    v.push(mv(d3, 0.0));
                }
                v.push(mv(d3, d3));
                v.push(mv(0.0, d3));
                if !south || rix != 0 {
                    v.push(mv(-d3, 0.0));
                }
                v.push(mv(-d3, -d3));
                v.push(mv(0.0, -d3));
            }
        }
        // Odd level, type C
        2 => v.extend([
            mv(d3, 0.0),
            mv(d23, 0.0),
            mv(d, d3),
            mv(d, d23),
            mv(d23, d23),
            mv(d3, d3),
        ]),
        // Odd level, type D
        _ => v.extend([
            mv(0.0, d3),
            mv(d3, d3),
            mv(d23, d23),
            mv(d23, d),
            mv(d3, d),
            mv(0.0, d23),
        ]),
    }
    v
}

/// The corners in the 5x6 plane from which the refined boundary of a zone is drawn, for a ring
/// bound for WGS84: five for a polar pentagon and six for every other zone, or `None` where the
/// eC answers no point. Ports `I3HZone::getBaseRefinedVertices(crs84, vertices)`
/// (`RI3H.ec:1792-2010`) with `crs84` true, compiled at `0xe540` to `0xf3bf`; the arms it takes
/// only towards the 5x6 CRS, the caps of the polar pentagons and the doubled corners of the odd
/// level, are not on the path of a ring bound for WGS84 and are not ported.
///
/// It is not [`vertices`], which ports `getVertices`, another function. This one runs
/// anticlockwise about every zone, where `getVertices` does so about the two polar pentagons
/// alone, and it begins at another corner; it gives a pentagon off the poles six corners, of
/// which [`refine5x6`] makes five; on a sub-hex C or D it draws two corners 2e-11 within the
/// cell's edge, where `getVertices` has them on it; and it does not leave a corner across an
/// interruption to be canonicalised later, but brings a corner that falls before the layout
/// forward by a period as it goes.
///
/// `None` is the eC's count of zero, for an identifier of a polar root with a non-zero index,
/// of sub-hex A or B (`RI3H.ec:1846`, `:1875`, `:1921`, `:1944`), from which
/// `getRefinedVertices` answers a null array; and a root that is none of the twelve, for which the eC reads a row and a column it
/// never set. `Grid` passes neither.
fn base_refined_vertices(id: u64) -> Option<Vec<(f64, f64)>> {
    let (level_i9r, root, rix, sh) = fields(id);
    if root > 11 || (root >= 10 && rix != 0) {
        return None;
    }
    // RI3H.ec:1796-1807, at `0xe567` to `0xe5b5`: the eC takes the row and the column from
    // `iLRCFromLRtI`, or from the two polar cells, and the corner from `I9RZone::ri5x6Extent`
    // (`RI9R.ec:526`, called at `0xe5aa`), which forms the same two products as `top_left`.
    let ((tx, ty), d) = top_left(level_i9r, root, rix);
    let mv = |dx: f64, dy: f64| move5x6_vertex(tx, ty, dx, dy);
    // faithful: sites GB3-1 and GB3-2, in DGGAL v0.0.6 as built with `-O2 -ffast-math`
    // (`libdggal.so`, BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`). Every `d/3` of this
    // function is compiled as `d * fl(1/3)` (GB3-1: `0xe610`, `0xe861`, `0xea0a`, `0xeaad`,
    // `0xec19`, `0xee1e` and `0xf014`) and every `2*d/3` as `d * fl(2/3)` (GB3-2: `0xe9e5`,
    // `0xead6`, `0xed43` and `0xee2e`): eleven products where the eC divides. The negated
    // offsets are those products with the sign turned (`xorpd` at `0xe61c`, `0xe87f`, `0xee3f`
    // and `0xee43`), save in the north polar arm of the even level, which multiplies by
    // `fl(-1/3)` (`0xec0c`), the same double. The two forms differ by a unit in the last
    // place at several levels, and either site put back in the eC's order moves rings of the
    // engine: few, some ten of the 88,000 sampled on a grid, and each is pinned by a ring read
    // from the engine. The guarded `move5x6Vertex` then adds each offset
    // to the corner in the eC's order, and the polar fans add their whole numbers as the eC
    // writes them.
    let d3 = d * T13;
    let d23 = d * T23;
    // RI3H.ec:1882-1883 and its fellows: a corner that falls before the layout, in either
    // coordinate, is brought forward by a period before the next corner is formed (`0xe8d4`
    // to `0xe8f2` and `0xf2c8`, and so on for each).
    let forward = |v: &mut Vec<(f64, f64)>| {
        if let Some(last) = v.last_mut() {
            if last.1 < 0.0 || last.0 < 0.0 {
                *last = (last.0 + 5.0, last.1 + 5.0);
            }
        }
    };
    // The eC opens both regular arms with the same test on `vertices[numPoints-1]` while
    // `numPoints` is still zero (RI3H.ec:1879-1880, :1948-1949), and the library compiles it as
    // written: it loads `vertices[-1]` at `0xe895` and `0xee61` and, where a member is
    // negative, stores it back plus five at `0xf258` and `0xf350`. That slot lies outside the
    // array, in the caller's frame (`rsp+0xc0` of `getRefinedVertices`), a temporary of
    // `ap.Add` that is dead at that moment and written afresh before each of its uses: the
    // read is of stack residue and the write reaches nothing. It is left out here.
    let mut v: Vec<(f64, f64)> = Vec::with_capacity(6);
    match sh {
        // Even level
        0 => {
            if root == 10 {
                // "North" pole (RI3H.ec:1816-1847, `0xebe8` to `0xec82` and `0xf390`): the
                // five corners are the first, one to each rhombus of the northern row.
                let (x, y) = mv(d3, d * MINUS_T13);
                v.extend([
                    (x + 5.0, y + 5.0),
                    (x + 1.0, y + 1.0),
                    (x + 2.0, y + 2.0),
                    (x + 3.0, y + 3.0),
                    (x + 4.0, y + 4.0),
                ]);
            } else if root == 11 {
                // "South" pole (RI3H.ec:1848-1876, `0xe5ec` to `0xe680` and `0xf3b0`). The
                // eC's `v.x - 0` is `v.x`.
                let (x, y) = mv(-d3, d3);
                v.extend([
                    (x, y),
                    (x - 1.0, y - 1.0),
                    (x - 2.0, y - 2.0),
                    (x - 3.0, y - 3.0),
                    (x + 1.0, y + 1.0),
                ]);
            } else {
                // Regular A (RI3H.ec:1877-1897, `0xee08` to `0xefc3` and `0xf225`): every
                // corner but the last is brought forward.
                v.push(mv(d3, -d3));
                forward(&mut v);
                v.push(mv(-d3, -d23));
                forward(&mut v);
                v.push(mv(-d23, -d3));
                forward(&mut v);
                v.push(mv(-d3, d3));
                forward(&mut v);
                v.push(mv(d3, d23));
                forward(&mut v);
                v.push(mv(d23, d3));
            }
        }
        // Odd level, type B
        1 => {
            if root == 10 {
                // "North" pole (RI3H.ec:1900-1922, `0xed20` to `0xedb0`). The eC's `v.x + 0`
                // is `v.x`.
                let (x, y) = mv(d23, 0.0);
                v.extend([
                    (x, y),
                    (x + 1.0, y + 1.0),
                    (x + 2.0, y + 2.0),
                    (x + 3.0, y + 3.0),
                    (x + 4.0, y + 4.0),
                ]);
            } else if root == 11 {
                // "South" pole (RI3H.ec:1923-1945, `0xeff0` to `0xf08a`).
                let (x, y) = mv(d3, d);
                v.extend([
                    (x, y),
                    (x - 1.0, y - 1.0),
                    (x - 2.0, y - 2.0),
                    (x - 3.0, y - 3.0),
                    (x - 4.0, y - 4.0),
                ]);
            } else {
                // Regular B (RI3H.ec:1946-1988, `0xe854` to `0xe9c6`). Neither the fourth
                // corner nor the last is brought forward: the eC has no test after the
                // fourth (RI3H.ec:1967-1975), and the library none either (`0xe979` to
                // `0xe98d`). With that test put in, no ring of some 87,000 zones a grid
                // moves: the fourth corner did not fall before the layout at any of them.
                v.push(mv(0.0, -d3));
                forward(&mut v);
                v.push(mv(-d3, -d3));
                forward(&mut v);
                v.push(mv(-d3, 0.0));
                forward(&mut v);
                v.push(mv(0.0, d3));
                v.push(mv(d3, d3));
                forward(&mut v);
                v.push(mv(d3, 0.0));
            }
        }
        // Odd level, type C (RI3H.ec:1990-1997, `0xeaa0` to `0xeb49`)
        2 => v.extend([
            mv(d3, d3),
            mv(d23, d23),
            mv(d, d23),
            mv(d, d3),
            mv(d23, 2e-11),
            mv(d3, 2e-11),
        ]),
        // Odd level, type D (RI3H.ec:1998-2005, `0xe9d0` to `0xea7b`)
        _ => v.extend([
            mv(2e-11, d23),
            mv(d3, d),
            mv(d23, d),
            mv(d23, d23),
            mv(d3, d3),
            mv(2e-11, d3),
        ]),
    }
    Some(v)
}

/// The eight directions of `I3HNeighbor` (`RI3H.ec:841-852`), in the order in which
/// `getNeighbors` asks `getNeighbor` for them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Direction {
    Top,
    Bottom,
    Left,
    Right,
    TopLeft,
    TopRight,
    BottomLeft,
    BottomRight,
}

impl Direction {
    const ALL: [Self; 8] = [
        Self::Top,
        Self::Bottom,
        Self::Left,
        Self::Right,
        Self::TopLeft,
        Self::TopRight,
        Self::BottomLeft,
        Self::BottomRight,
    ];
}

/// The zone across the side of zone `id` that `which` names, or `None` for the eC's null zone,
/// which it answers where the direction names no side of the zone. Ports
/// `I3HZone::getNeighbor` (`RI3H.ec:1005-1161`), compiled at `0x10ca0` to `0x11557`; the eC's
/// local names are kept.
fn neighbor(id: u64, which: Direction) -> Option<u64> {
    use Direction::{Bottom, BottomLeft, BottomRight, Left, Right, Top, TopLeft, TopRight};
    let (level_i9r, _, _, sh) = fields(id);
    let (c_x, c_y) = centroid(id);
    let cx = floor_i32(c_x + 1e-11);
    let cy = floor_i32(c_y + 1e-11);
    // faithful: `centroid.y - centroid.x - 1E-11 > 1` (RI3H.ec:1010) is compiled as
    // `c.y - c.x > fl(1 + 1e-11)`, and so is its re-test for `top` (:1029) and for `right`
    // (:1137), in DGGAL's `-O2 -ffast-math` build, `libdggal.so`, BuildID
    // `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, at `0x10d50`, `0x110b0` and `0x113b0`.
    // `c.y - c.x` lies within about `1e-15` of a multiple of `1 / (3p)`, never within a unit
    // of either threshold, so this cannot change an answer.
    let south = c_y - c_x > ONE_PLUS_1E11; // Not counting pentagons as south or north
    // faithful: `centroid.x - centroid.y - 1E-11 > 0` (RI3H.ec:1011) is compiled as
    // `c.x - c.y > 1e-11`, and so is its re-test for `right` (:1134), in DGGAL's
    // `-O2 -ffast-math` build, `libdggal.so`, BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`,
    // at `0x10d22` and `0x10fc0`; `fl(a - e) > 0` exactly when `a > e`, so the two are
    // equivalent.
    let north = c_x - c_y > 1e-11;
    let north_pole = north && (c_x - c_y - 1.0).abs() < 1e-11;
    let south_pole = south && (c_y - c_x - 2.0).abs() < 1e-11;
    let d = 1.0 / pow3(level_i9r) as f64;
    let (mut x, mut y) = (0.0, 0.0);
    let mut cross_early = true;

    if sh == 0 {
        // Even level
        //
        // NOTE: See getNeighbors() for special interruption cases (the eC's note)
        match which {
            Top => {
                if south && c_x - f64::from(cx) < 1e-11 {
                    cross_early = false;
                    if south_pole {
                        (x, y) = (-3.0, -3.0 - d);
                    } else {
                        // Extra top neighbor at south interruptions
                        y = -d;
                    }
                }
            }
            Bottom => {
                if north && c_y - f64::from(cy) < 1e-11 {
                    cross_early = false;
                    if north_pole {
                        (x, y) = (2.0 - d, 2.0);
                    } else {
                        // Extra bottom neighbor at north interruptions
                        x = -d;
                    }
                }
            }
            Left => (x, y) = (-d, -d),
            Right => (x, y) = (d, d),
            TopLeft => {
                if north_pole {
                    cross_early = false;
                    (x, y) = (3.0 - d, 3.0);
                } else if south_pole {
                    cross_early = false;
                    y = -d;
                } else {
                    y = -d;
                }
            }
            BottomLeft => {
                if south_pole {
                    cross_early = false;
                    (x, y) = (-2.0, -2.0 - d);
                } else {
                    x = -d;
                }
            }
            TopRight => {
                if north_pole {
                    cross_early = false;
                    (x, y) = (4.0 - d, 4.0);
                } else if south_pole {
                    cross_early = false;
                    (x, y) = (-4.0, -d - 4.0);
                } else {
                    x = d;
                }
            }
            BottomRight => {
                if south_pole {
                    cross_early = false;
                    (x, y) = (-1.0, -1.0 - d);
                } else {
                    y = d;
                }
            }
        }
    } else {
        // Odd level
        //
        // faithful: `do3 = d/3` (RI3H.ec:1083) is compiled as `d * fl(1/3)`, in DGGAL's
        // `-O2 -ffast-math` build, `libdggal.so`, BuildID
        // `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, at `0x10e00`, so that every offset below
        // is formed from that double: `2*do3` as `do3 + do3` (`0x11448`, `0x114bc`) and
        // `-2*do3` as `-2.0 * do3` (`0x11418`, `0x11488`), each the same double as the product
        // written here. The offset only places the target on a neighbour's centre, which the
        // tolerances of `move5x6Vertex2` and the quantiser absorb, so this cannot change an
        // answer.
        let do3 = d * T13;

        // NOTE: See getNeighbors() for special interruption cases (the eC's note)
        match which {
            Top => {
                if south_pole {
                    (x, y) = (do3 - 5.0, -do3 - 5.0);
                    cross_early = false;
                } else if !north_pole {
                    (x, y) = (do3, -do3);
                }
            }
            Bottom => {
                if north_pole {
                    (x, y) = (1.0 - do3, 1.0 + do3);
                    cross_early = false;
                } else if !south_pole {
                    (x, y) = (-do3, do3);
                }
            }
            TopLeft => {
                if north_pole {
                    (x, y) = (2.0 - do3, 2.0 + do3);
                    cross_early = false;
                } else if south_pole {
                    (x, y) = (do3, -do3);
                } else {
                    (x, y) = (-do3, -2.0 * do3);
                }
            }
            BottomLeft => {
                if north_pole {
                    (x, y) = (4.0 - do3, 4.0 + do3);
                    cross_early = false;
                } else if south_pole {
                    (x, y) = (do3 - 2.0, -do3 - 2.0);
                } else {
                    (x, y) = (-2.0 * do3, -do3);
                }
            }
            TopRight => {
                if north_pole {
                    (x, y) = (3.0 - do3, 3.0 + do3);
                    cross_early = false;
                } else if south_pole {
                    (x, y) = (do3 - 4.0, -do3 - 4.0);
                } else {
                    (x, y) = (2.0 * do3, do3);
                }
            }
            BottomRight => {
                if north_pole {
                    (x, y) = (5.0 - do3, 5.0 + do3);
                    cross_early = false;
                } else if south_pole {
                    (x, y) = (do3 - 1.0, -do3 - 1.0);
                } else {
                    (x, y) = (do3, 2.0 * do3);
                }
            }
            // Currently stand-in for second bottom / top neighbor (the eC's note)
            Right => {
                if north && !north_pole && c_y - f64::from(cy) < 1e-11 {
                    // Extra bottom neighbor at north interruptions
                    cross_early = false;
                    (x, y) = (-do3, do3);
                } else if south && !south_pole && c_x - f64::from(cx) < 1e-11 {
                    // Extra bottom neighbor at south interruptions (the eC's words)
                    cross_early = false;
                    (x, y) = (do3, -do3);
                }
            }
            Left => {}
        }
    }
    if x != 0.0 || y != 0.0 {
        // REVIEW: This is the only place we use moveISEAVertex2() (the eC's note)
        let (vx, vy) = move5x6_vertex2(c_x, c_y, x, y, cross_early);
        let level = 2 * level_i9r as u32 + u32::from(sh > 0);
        let result = from_centroid(level, vx, vy)?;
        // This should not happen (the eC's note)
        (result != id).then_some(result)
    } else {
        None
    }
}

/// The neighbours of zone `id` in the engine's order. Ports `I3HZone::getNeighbors`
/// (`RI3H.ec:1163-1201`): the directions in the order of `I3HNeighbor`, each kept where it
/// names a zone, except that a `topRight` naming the zone that the `topLeft` just kept names is
/// merged into it, and a `bottomRight` into a `bottomLeft` alike. No other repeat is removed.
fn neighbors(id: u64) -> Vec<u64> {
    use Direction::{Bottom, BottomLeft, BottomRight, Top, TopLeft, TopRight};
    let mut kept: Vec<(u64, Direction)> = Vec::with_capacity(6);
    for n in Direction::ALL {
        let Some(nb) = neighbor(id, n) else {
            continue;
        };
        let mut which = n;
        if let Some(last) = kept.last_mut() {
            // Handle special cases here so that getNeighbor() can still return same neighbor
            // for multiple directions (the eC's note)
            if n == TopRight && last.1 == TopLeft && last.0 == nb {
                last.1 = Top;
                continue;
            } else if n == BottomRight && last.1 == BottomLeft && last.0 == nb {
                last.1 = Bottom;
                continue;
            } else if n == TopRight && last.1 != TopLeft {
                which = Top;
            } else if n == BottomRight && last.1 != BottomLeft {
                which = Bottom;
            }
        }
        kept.push((nb, which));
    }
    kept.into_iter().map(|(nb, _)| nb).collect()
}

/// The eC's `nullZone`, which the hierarchy passes along as a value, as the eC does.
const NULL_ZONE: u64 = crate::ZoneId::NULL.0;

/// The row and column of a cell of rhombus 0 to 9 at I9R level `level`, or `None` for the eC's
/// `-1`, where the level exceeds 16 or the index lies beyond the rhombus. Ports
/// `iLRCFromLRtI` (`RI9R.ec:710-729`), compiled at `0x12dd0` to `0x12e41`, as the I3H members
/// call it, with the level letter of `levelI9R`; the row and column are `int`s.
fn i_lrc_from_lrt_i(level: u64, root: u64, ix: u64) -> Option<(i32, i32)> {
    if level > 16 || root > 9 {
        return None;
    }
    let p = pow3(level);
    if ix >= p * p {
        return None;
    }
    let (row, col) = row_col(p, root, ix);
    Some((row as i32, col as i32))
}

/// The zone at I9R level `level` whose cell lies at `row` and `col` in the whole layout, with
/// the sub-hex `sub_hex` (0 to 3 for A to D), or the pole `pole` where it is not 0; the null
/// zone where the cell is out of range. Ports `I3HZone::fromI9R` (`RI3H.ec:951-963`),
/// compiled at `0xb910` to `0xba43`; the eC's local names are kept, and the sub-hex, which
/// every caller passes in range, is not re-checked.
///
/// Only `parent0` passes a negative level, for an identifier with no cell. `POW3` is then
/// `(uint64)(pow(3, level) + 0.1)`, which is 0, and the eC divides by it and dies with
/// `SIGFPE` (`0xb940`); this answers the null zone instead.
fn from_i9r(level: i32, row: u32, col: u32, sub_hex: u64, pole: u32) -> u64 {
    let Ok(level) = u64::try_from(level) else {
        return NULL_ZONE;
    };
    let p = pow3(level);
    // The quotients are `uint`s and the differences `int`s, which gcc forms in 32 bits: the low
    // half of the eC's `uint64` arithmetic, as the truncations here keep it.
    let row_op = (u64::from(row) / p) as u32;
    let col_op = (u64::from(col) / p) as u32;
    let root = (if pole != 0 {
        pole
    } else {
        row_op.wrapping_add(col_op)
    }) as i32;
    let y = row.wrapping_sub(row_op.wrapping_mul(p as u32)) as i32;
    let x = col.wrapping_sub(col_op.wrapping_mul(p as u32)) as i32;
    // Each `int` widens to `uint64` with its sign, for the index and the tests alike.
    let (y, x) = (i64::from(y) as u64, i64::from(x) as u64);
    let ix = if pole != 0 {
        0
    } else {
        y.wrapping_mul(p).wrapping_add(x)
    };
    // Avoid returning bad key (the eC's note)
    if root > 11 || (root < 10 && (row_op < col_op || row_op - col_op > 1 || y >= p || x >= p)) {
        return NULL_ZONE;
    }
    pack(level, u64::from(root as u32), ix, sub_hex)
}

/// A zone's primary parent, or the null zone at level 0. Ports the `I3HZone` property
/// `parent0` (`RI3H.ec:977-1003`), compiled at `0xba50` to `0xbc48`: at an odd level, the same
/// zone with its sub-hex cleared, the even-level cell whose centre it shares; at an even level,
/// the odd-level cell one I9R level up at `(row / 3, col / 3)`, with the sub-hex B, C or D that
/// the residues of its index modulo 3 decide, and B for a pole.
///
/// The eC divides by zero and dies with `SIGFPE` on two kinds of identifier with no cell: a
/// rhombus of 0 to 9 whose level or index `iLRCFromLRtI` refuses, in its own division by
/// `POW3(-1)` (`0xbbed`); and a root of 13 to 15, or of 10 to 12 with a non-zero index, in
/// `fromI9R`'s division by `POW3(-2)`, gcc having dropped this function's own division there,
/// whose quotient no sub-hex of a pole reads. This answers the null zone for both.
fn parent0(id: u64) -> u64 {
    let (level_i9r, root, rix, sh) = fields(id);
    if level_i9r == 0 && sh == 0 {
        return NULL_ZONE;
    }
    if sh > 0 {
        return id & !3; // key = this; key.subHex = 0
    }
    if root > 9 {
        // The eC leaves the row and column unset here (`0xbaad`); `fromI9R` reads neither for a
        // pole, and answers the null zone for a root of 12.
        let level = if root <= 12 && rix == 0 {
            level_i9r as i32
        } else {
            -1
        };
        return from_i9r(level - 1, 0, 0, 1, root as u32);
    }
    let Some((row, col)) = i_lrc_from_lrt_i(level_i9r, root, rix) else {
        return NULL_ZONE;
    };
    let p = pow3(level_i9r);
    let (r, c) = (rix / p, rix % p);
    let (rm3, cm3) = (r % 3, c % 3);
    let sub_hex = if cm3 > 1 {
        2
    } else if rm3 > 1 {
        3
    } else {
        1
    };
    from_i9r(
        level_i9r as i32 - 1,
        (row / 3) as u32,
        (col / 3) as u32,
        sub_hex,
        0,
    )
}

/// Whether a zone is a centroid child: every odd-level B, both poles, and an even-level A whose
/// row and column within its rhombus are both multiples of 3, or whose sum is. Ports the
/// `I3HZone` property `isCentroidChild` (`RI3H.ec:2170-2196`), compiled at `0xb430` to
/// `0xb4ec`.
fn is_centroid_child(id: u64) -> bool {
    if id == NULL_ZONE {
        return false;
    }
    let (level_i9r, root, rix, sh) = fields(id);
    if sh > 0 {
        return sh == 1; // B are centroid children
    }
    if root == 10 || root == 11 {
        return true; // Polar A are centroid children
    }
    let p = pow3(level_i9r);
    let (r, c) = (rix / p, rix % p);
    (r % 3 == 0 && c % 3 == 0) || (r + c) % 3 == 0
}

/// A zone's parents, in the engine's order: none at level 0; `parent0` alone for a centroid
/// child; otherwise `parent0` and two of its neighbours, three entries whether or not those
/// neighbours exist, each the null zone where `getNeighbor` names none. Ports
/// `I3HZone::getParents` (`RI3H.ec:1247-1328`), compiled at `0x11710` to `0x11a54`, the one
/// copy, which `centroidParent` and `getZoneParents` call; the eC's local names are kept.
fn parents(id: u64) -> Vec<u64> {
    use Direction::{Bottom, BottomLeft, BottomRight, Right, Top, TopLeft, TopRight};
    let parent0 = parent0(id);
    if is_centroid_child(id) {
        return if parent0 == NULL_ZONE {
            Vec::new()
        } else {
            vec![parent0]
        };
    }
    let sh = id & 3;
    let (first, second) = if sh > 0 {
        // Odd level
        (Right, if sh == 2 { TopRight } else { BottomRight })
    } else {
        // Even level
        let (c_x, c_y) = centroid(id);
        let (p0_x, p0_y) = centroid(parent0);
        let dx = c_x - p0_x;
        let dy = c_y - p0_y;
        // Both floors are inlined as `floor_i32` describes, at `0x117ac` and `0x1189f`.
        let p0cx = floor_i32(p0_x + 1e-11);
        let p0cy = floor_i32(p0_y + 1e-11);
        // faithful: `p0Centroid.y - p0Centroid.x + 1E-11 > 1` (RI3H.ec:1273) is compiled as
        // `p0.y - p0.x > fl(1 - 1e-11)`, in DGGAL's `-O2 -ffast-math` build, `libdggal.so`,
        // BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, at `0x117c8`. A difference of
        // two centroids lies within about `1e-15` of a multiple of `1 / (3p)`, never within a
        // unit of either threshold, so this cannot change an answer.
        let on_bottom_crossing_left =
            p0_y - p0_x > ONE_MINUS_1E11 && p0_x - f64::from(p0cx) < 1e-11;
        // faithful: `p0Centroid.x - p0Centroid.y + 1E-11 > 0` (RI3H.ec:1274) is compiled as
        // `p0.x - p0.y > -1e-11`, in the same library at `0x11800`; `fl(a + e) > 0` exactly when
        // `a > -e`, so the two are equivalent.
        let on_top_crossing_right = p0_x - p0_y > -1e-11 && p0_y - f64::from(p0cy) < 1e-11;

        if dx.abs() < 1e-11 {
            if dy > 0.0 {
                // faithful: `p0Centroid.x - p0Centroid.y - 1E-11 > 0` (RI3H.ec:1280) is compiled
                // as `p0.x - p0.y > 1e-11`, in the same library at `0x1194e`; `fl(a - e) > 0`
                // exactly when `a > e`, so the two are equivalent.
                let on_top_crossing_right_neg_epsilon =
                    p0_x - p0_y > 1e-11 && p0_y - f64::from(p0cy) < 1e-11;
                // Bottom-Right vertex child of p0
                let second = if on_bottom_crossing_left {
                    BottomLeft
                } else if on_top_crossing_right_neg_epsilon {
                    Right
                } else {
                    Bottom
                };
                (BottomRight, second)
            } else {
                // Top-Left vertex child of p0
                (TopLeft, Top)
            }
        } else if dy.abs() < 1e-11 {
            if dx > 0.0 {
                // faithful: `p0Centroid.y - p0Centroid.x - 1E-11 > 1` (RI3H.ec:1297) is compiled
                // as `p0.y - p0.x > fl(1 + 1e-11)`, in the same library at `0x118f4`; it cannot
                // change an answer, by the argument given at `0x117c8`.
                let on_bottom_crossing_left_neg_epsilon =
                    p0_y - p0_x > ONE_PLUS_1E11 && p0_x - f64::from(p0cx) < 1e-11;
                // Top-Right vertex child of p0
                let second = if on_top_crossing_right {
                    TopLeft
                } else if on_bottom_crossing_left_neg_epsilon {
                    Right
                } else {
                    Top
                };
                (TopRight, second)
            } else {
                // Bottom-Left vertex child of p0
                (BottomLeft, Bottom)
            }
        } else if dx > 0.0 {
            // Right vertex child of p0
            (TopRight, BottomRight)
        } else {
            // Left vertex child of p0
            (TopLeft, BottomLeft)
        }
    };
    let nb = |which| neighbor(parent0, which).unwrap_or(NULL_ZONE);
    vec![parent0, nb(first), nb(second)]
}

/// The first of a zone's parents that is itself a centroid child: `parent0` if it is one,
/// otherwise the first of the two others that is one, and otherwise the null zone. Ports the
/// `I3HZone` property `centroidParent` (`RI3H.ec:1203-1223`), compiled at `0x11a60` to
/// `0x11acd`.
fn centroid_parent(id: u64) -> u64 {
    let c_parent = parent0(id);
    if c_parent != NULL_ZONE && is_centroid_child(c_parent) {
        return c_parent;
    }
    parents(id)
        .into_iter()
        .skip(1)
        .find(|&p| is_centroid_child(p))
        .unwrap_or(NULL_ZONE)
}

/// The child that shares a zone's centre, or the null zone for the null zone: at an even level
/// the same cell's sub-hex B; at an odd level the even-level cell one I9R level down at
/// `(3 row + a, 3 col + b)`, `a` and `b` by the sub-hex; for a pole the pole of the next level.
/// Ports the `I3HZone` property `centroidChild` (`RI3H.ec:2012-2042`), compiled at `0xbc50` to
/// `0xbd35`, whose `int` arithmetic gcc forms in 32 bits, as the truncations below keep it.
pub(super) fn centroid_child(id: u64) -> u64 {
    if id == NULL_ZONE {
        return NULL_ZONE;
    }
    let (l9r, root, rix, sh) = fields(id);
    if sh == 0 {
        // Centroid child for Even level (including poles)
        return pack(l9r, root, rix, 1);
    }
    if root > 9 {
        // Odd level "North" and "South" Poles
        return pack(l9r + 1, root, 0, 0);
    }
    // Centroid child for Odd level
    let (row, col) = row_col(pow3(l9r), root, rix);
    let r = (row as u32).wrapping_mul(3).wrapping_add(match sh {
        3 => 2,
        2 => 1,
        _ => 0,
    });
    let c = (col as u32).wrapping_mul(3).wrapping_add(match sh {
        2 => 2,
        3 => 1,
        _ => 0,
    });
    from_i9r(l9r as i32 + 1, r, c, 0, 0)
}

/// A zone's children, in the engine's order: its centroid child, then one child on each of its
/// vertices, seven in all for a hexagon and six for a pentagon; none at resolution 33. Ports
/// `I3HZone::getChildren` (`RI3H.ec:2044-2085`), compiled at `0x10a50` to `0x10c8b`: a pole's
/// vertex children by closed formulas, and every other zone's by quantising each of its
/// vertices at the next level, as `getVertices` gives them, before any canonicalisation.
///
/// One deliberate departure. At resolution 33, the finest, the engine answers seven level-34
/// identifiers, whose indices spill into the root and level fields of the identifier, whose
/// text it cannot read back, and whose parents it cannot list without dividing by zero; this
/// answers none, as the aperture-7 grids do at their finest resolution.
fn children(id: u64) -> Vec<u64> {
    let (l9r, root, _, sh) = fields(id);
    if 2 * l9r + u64::from(sh > 0) >= u64::from(MAX_LEVEL) {
        return Vec::new();
    }
    let mut children = vec![centroid_child(id)];
    if root > 9 {
        // Special cases for the poles; every root above 10 takes the southern formula, as in
        // the eC.
        let p = pow3(l9r);
        let north = root == 10;
        children.extend((0..5).map(|i: u64| match (sh == 0, north) {
            // Even level
            (true, true) => pack(l9r, i * 2, p - 1, 2),
            (true, false) => pack(l9r, i * 2 + 1, p * (p - 1), 3),
            // Odd level
            (false, true) => pack(l9r + 1, i * 2, 3 * p - 1, 0),
            (false, false) => pack(l9r + 1, i * 2 + 1, 3 * p * (3 * p - 1), 0),
        }));
    } else {
        let next_level = 2 * l9r as u32 + 1 + u32::from(sh > 0);
        children.extend(
            vertices(id)
                .into_iter()
                .map(|(x, y)| from_centroid(next_level, x, y).unwrap_or(NULL_ZONE)),
        );
    }
    children
}

/// A zone's level, `levelI9R * 2 + (subHex > 0)`, read from the fields of the identifier as
/// they stand. Ports the `I3HZone` property `level` (`RI3H.ec:878-881`).
fn level(id: u64) -> u64 {
    let (level_i9r, _, _, sub_hex) = fields(id);
    2 * level_i9r + u64::from(sub_hex > 0)
}

/// The zones two levels up that hold the whole of this one: one for a zone that lies wholly
/// within a grandparent, three for one that three grandparents share. The entries are the
/// eC's as they come, the null zone among them where a parent does not exist. Ports
/// `I3HZone::getContainingGrandParents` (`RI3H.ec:1225-1245`).
fn containing_grand_parents(id: u64) -> Vec<u64> {
    let c_parent = centroid_parent(id);
    if is_centroid_child(id) {
        if c_parent != NULL_ZONE {
            vec![parent0(c_parent)]
        } else {
            parents(parent0(id))
        }
    } else {
        parents(c_parent)
    }
}

/// The compaction of a set of zones: every grandparent whose grandchildren are all in the set
/// stands for those of them that it alone holds, two levels a pass from the finest level of
/// the set upwards, in ascending order of the identifiers. Ports
/// `RhombicIcosahedral3H::compactZones` (`RI3H.ec:161-186`) with `compactI3HZones`
/// (`RI3H.ec:2250-2386`), statement for statement; the eC's local names are kept.
///
/// The eC reads the array into a tree ordered by the identifier's value, so that a zone given
/// twice counts once, and leaves the null zone out; `BTreeSet` is that tree. A grandparent is
/// taken when the grandchild at its centre has every neighbour in the set, and when each of
/// the six grandchildren on its vertices is in the set too or, failing that, has the zone two
/// levels finer at its own centre already kept in the output. A zone leaves the set when
/// every grandparent that holds it was taken, which is one for the seven about the centre and
/// three for one on a vertex: so the six on the vertices of a lone grandparent stay beside
/// it, and overlap it.
///
/// Three things in the eC look like oversights and are the engine's answer, so all are kept.
/// The grandchild at the centre is not itself looked for, so that a set which lacks it and
/// has its neighbours is compacted as if it had it, and the answer then covers a zone that
/// the set did not hold. `output` is not emptied between passes. And the loop goes on after a
/// pass that found nothing to take (the eC's `break` is commented out).
///
/// The answer covers a zone that the set did not hold by a second route as well, which is
/// the rule as written and no oversight: from the second pass on, a grandchild on a vertex
/// that the pass lacks is taken as present on the evidence of the one zone at its own centre,
/// kept in the output, and the grandparent then covers the part of itself about that vertex
/// whatever of it the set held.
///
/// The zones must be ones the engine's own parser reads back, which `Grid` checks first, and
/// the function is meant for zones of one level, which `Grid` checks too: given zones of
/// several levels the engine drops those with no grandparent to be dropped for, as this
/// does. The cost is that of the set: at most seventeen passes, each a bounded amount of
/// work for every zone it holds, and a search of a tree for each.
fn compact(input: &[u64]) -> Vec<u64> {
    use std::collections::BTreeSet;
    let mut max_level = 0;
    let mut zones: BTreeSet<u64> = BTreeSet::new();
    for &zone in input {
        if zone != NULL_ZONE {
            max_level = max_level.max(level(zone) as i64);
            zones.insert(zone);
        }
    }

    let mut output: BTreeSet<u64> = BTreeSet::new();
    let mut next: BTreeSet<u64> = BTreeSet::new();
    let mut l = max_level - 2;
    while l >= 0 {
        for &zone in &zones {
            for g_parent in containing_grand_parents(zone) {
                // The eC tests the grandparent against the null zone alone. One value passes
                // that test and is no zone: `parent0` of the null zone is the null zone with
                // its sub-hex cleared, `0xffff_ffff_ffff_fffc`, which the grandparents of a
                // zone of level 0 come to. The eC carries on with it, reads relations for it
                // that mean nothing, and never keeps it; this crate's relations for it would
                // be empty and would let it through, so it is refused here with every other
                // identifier that is no zone. A set of one level never reaches it: a zone of
                // level 0 is in the set of a pass only where the input holds a finer zone
                // beside it.
                if g_parent == NULL_ZONE || !is_readable(g_parent) || next.contains(&g_parent) {
                    continue;
                }
                let c_zone = centroid_child(centroid_child(g_parent));
                let mut parent_all_in = neighbors(c_zone)
                    .into_iter()
                    .all(|nb| nb == NULL_ZONE || zones.contains(&nb));
                if parent_all_in {
                    // Grandparent vertex children's centroid children are partially within it
                    // and must be present to perform replacement (the eC's own note).
                    for ch in children(g_parent).into_iter().skip(1) {
                        if ch == NULL_ZONE {
                            continue;
                        }
                        let c_child = centroid_child(ch);
                        if !zones.contains(&c_child) {
                            let (cx, cy) = centroid(c_child);
                            let sub = from_centroid(level(c_child) as u32 + 2, cx, cy)
                                .unwrap_or(NULL_ZONE);
                            if !output.contains(&sub) {
                                parent_all_in = false;
                            }
                        }
                    }
                    if parent_all_in {
                        next.insert(g_parent);
                    }
                }
            }
        }

        for &zone in &zones {
            let all_in = containing_grand_parents(zone)
                .iter()
                .all(|g_parent| next.contains(g_parent));
            if !all_in {
                output.insert(zone);
            }
        }

        if l - 2 >= 0 && !next.is_empty() {
            // Not done: the next level becomes the zones to compact.
            zones = std::mem::take(&mut next);
        } else {
            // Done: next is combined with output into the final zones.
            zones = output.clone();
            zones.extend(next.iter().copied());
        }
        l -= 2;
    }

    if zones.len() >= 32 && zones.first().is_some_and(|&zone| level(zone) == 1) {
        let mut n_l1 = 0;
        for &zone in &zones {
            match level(zone) {
                1 => n_l1 += 1,
                0 => {}
                _ => break,
            }
        }
        if n_l1 == 32 {
            // Simplifying the full globe to the zones of level 0.
            zones = (0..12).map(|root| pack(0, root, 0, 0)).collect();
        }
    }
    zones.into_iter().collect()
}

impl Topology for HexA3 {
    const APERTURE: u8 = 3;

    /// Exactly the identifiers [`crate::indexings::I3h`] reads back from its own text, as the eC's
    /// `fromZoneID` reads back only the text it prints (`RI3H.ec:908-931`).
    ///
    /// Every other identifier answers as the null zone through `Grid`: the zero point as its
    /// centroid, six zero vertices and no relations. That covers the null zone itself, a
    /// departure from the engine, which gives it a meaningless centroid and vertices and a
    /// parent chain on which it divides by zero; the sub-hexes C and D of a polar pentagon,
    /// which the engine reads from text although no relation of any zone reaches them; the
    /// identifiers of a polar root with a non-zero index, which the engine's quantiser answers
    /// along two seams of the layout and its own parser refuses (see [`HexA3::quantize`]); and
    /// every identifier with a level, root or index out of range, bits 62 and 63 set, or a level
    /// of 34 or more.
    fn is_valid_address(a: &Address) -> bool {
        is_readable(a.base)
    }

    /// Quantises the planar point `p` to the zone at resolution `res` that contains it, as
    /// `getZoneFromCRSCentroid` does in the 5x6 CRS (`RI3H.ec:94-117`), through
    /// `I3HZone::fromCentroid`. `p.face` is not read, and a resolution above 33 gives `None`.
    ///
    /// `None` is the eC's null zone, where `fromCentroid` finds no cell (`RI3H.ec:1430-1431`).
    /// A geographic point reaches it near both poles at resolutions 0 and 1, and at the coarsest
    /// resolutions in a narrow band beside each of the ten icosahedron edges that end at the
    /// vertices of the two polar pentagons: a probe of the engine met it there at resolutions 0
    /// to 7, never further than 2e-6 radians from the edge.
    ///
    /// One deliberate departure, at two seams of the layout. A point within
    /// `1e-6 / 3^(res / 2)` of the seam `x = 5` in rhombus 8, or of the seam `y = 6` in rhombus
    /// 9, is re-anchored twice (`RI3H.ec:1374-1413`), first across the seam and then into the
    /// next rhombus, and leaves with a root of 10, or of 11, without the modulo the later
    /// branches apply. The eC's quantiser then answers an identifier of a polar root with a
    /// non-zero index. The engine's own text parser refuses it (`RI3H.ec:919-920`); the engine
    /// gives it the geometry of the polar pentagon, up to some 63 degrees away, lists seven
    /// neighbours for it into arrays of six, and divides by zero among its parents. This answers
    /// `None` wherever the root is 10 or more and the index is not zero, as the aperture-7 grids
    /// answer the null zone where the engine names a zone it cannot read back.
    ///
    /// Neither class is a single identifier near a pole. Measured on ISEA3H, and met again by the
    /// oracle suites on all three aperture-3 grids, both are thin bands met at every resolution
    /// from 2 to 21, each with many identifiers. Root 10 (`BA-2-A`, `BA-1-A`, `CA-3-A` and more)
    /// lies along the whole icosahedron edge that runs over the geographic north pole, on the
    /// meridians 11.2 E and 168.8 W from the pole down to 58.4 N.
    /// Root 11 (`BB-6-A`, `CB-12-A` and more) lies in the southern hemisphere, from the equator
    /// to 58.4 S between the meridians 137 W and 169 W, along the edge that ends at the southern
    /// polar vertex, 58.4 S 168.8 W.
    fn quantize(p: PlanarPoint, res: u8) -> Option<Address> {
        if res > MAX_LEVEL {
            return None;
        }
        let id = from_centroid(u32::from(res), p.x, p.y)?;
        let (_, root, rix, _) = fields(id);
        if root >= 10 && rix != 0 {
            return None;
        }
        Some(Address::new(id, &[]))
    }

    /// The zone's planar centroid in the 5x6 layout, exactly as the `centroid` property gives it
    /// and as `getZoneCRSCentroid` returns it; `face` carries the `-1` sentinel.
    ///
    /// Unlike the aperture-7 topology's, which folds the centroid as the eC's `getVertices` does,
    /// this is not canonicalised: the engine passes the property to the inverse projection as it
    /// stands (`RI3H.ec:268`), and so does `Grid`. The poles are `(1, 0)` and `(4, 6)`.
    ///
    /// `a` must satisfy [`HexA3::is_valid_address`]; on any other address the answer is
    /// meaningless, though the arithmetic wraps as the eC's does and does not panic.
    fn planar_centroid(a: &Address) -> PlanarPoint {
        let (x, y) = centroid(a.base);
        PlanarPoint { face: -1, x, y }
    }

    /// The zone's planar vertices, five for a pentagon and six for a hexagon, each passed through
    /// `canonicalize5x6` (`ri5x6.ec:1562-1606`) as the engine passes every vertex before the
    /// inverse projection (`RI3H.ec:380-381`), so that a vertex on an interruption or beyond the
    /// layout is given in the engine's form. They are the vertices `getZoneCRSVertices` returns
    /// in the 5x6 CRS, which canonicalises them alike (`RI3H.ec:342-343`), to the last bit.
    ///
    /// One deliberate departure, in their order: the ring runs anticlockwise as seen from
    /// outside the sphere, as it does on every grid of this crate (with the aperture-7
    /// exceptions that [`Grid::vertices`](crate::Grid::vertices) sets out: rings with every
    /// vertex on a pole, and rings that hold the point (0, 0)). The engine's ring runs clockwise
    /// for every zone but the two polar pentagons, roots 10 and 11, whose ring runs anticlockwise
    /// at every resolution. For every other zone the
    /// engine's `[v0, v1, ..., v(n-1)]` is given as `[v0, v(n-1), ..., v1]`: it starts at the
    /// engine's own first vertex and turns the other way. The polar pentagons keep the engine's
    /// order. The sub-zone orders are generated from the engine's own order, which this does
    /// not touch.
    ///
    /// `a` must satisfy [`HexA3::is_valid_address`], as for [`HexA3::planar_centroid`].
    fn planar_vertices(a: &Address) -> Vec<PlanarPoint> {
        let mut ring: Vec<PlanarPoint> = vertices(a.base)
            .into_iter()
            .map(|(vx, vy)| {
                let (x, y) = canonicalize5x6(vx, vy);
                PlanarPoint { face: -1, x, y }
            })
            .collect();
        // departure: the engine's ring, clockwise but at the polar pentagons, is turned about
        // its first vertex so that every ring runs anticlockwise, as on the aperture-7 grids;
        // every coordinate stays the engine's own.
        let (_, root, _, _) = fields(a.base);
        if root < 10 {
            ring[1..].reverse();
        }
        ring
    }

    /// The zone's refined boundary in the 5x6 plane, before any projection, as
    /// `getRefinedVertices` builds it for WGS84 (`RI3H.ec:456`, `:507`): the corners of
    /// `getBaseRefinedVertices` with `crs84` true, every side then divided into `n_divisions`
    /// parts by `refine5x6` with `wrap` true. The points are in the engine's own sequence,
    /// which runs anticlockwise; nothing is turned here, as [`HexA3::planar_vertices`] turns
    /// the plain ring.
    ///
    /// `a` must satisfy [`HexA3::is_valid_address`]. `None`, which `Grid` gives as the empty
    /// ring, is the eC's null array, for an identifier of a polar root with a non-zero index;
    /// no valid address is one.
    fn planar_refined_vertices(a: &Address, n_divisions: i32) -> Option<Vec<(f64, f64)>> {
        Some(refine5x6(&base_refined_vertices(a.base)?, n_divisions))
    }

    /// The zone's neighbours, as the engine's `getZoneNeighbors` (`RI3H.ec:119-122`) lists them:
    /// six for a hexagon and five for a pentagon, in the engine's order, which is the order in
    /// which [`crate::Disk`] walks them. They are found in the planar layout, by stepping from the
    /// centroid across each side and quantising where the step lands, so no projection is
    /// involved.
    ///
    /// `a` must satisfy [`HexA3::is_valid_address`], which `Grid` checks first, so that the null
    /// zone and every other identifier with no cell have no neighbours. The engine lists seven
    /// for an identifier of root 10 with a non-zero index, into arrays of six; that class has no
    /// cell here. On any other address the answer is meaningless, though the arithmetic wraps
    /// as the eC's does and does not panic.
    fn neighbors(a: &Address) -> Option<Vec<Address>> {
        Some(
            neighbors(a.base)
                .into_iter()
                .map(|id| Address::new(id, &[]))
                .collect(),
        )
    }

    /// The zone's parents, as the engine's `getZoneParents` (`RI3H.ec:134-137`) lists them and in
    /// its order: none at resolution 0; one for a centroid child, the cell one level up that
    /// shares its centre or contains it, DGGAL's `parent0`; and three for every other zone, since
    /// aperture 3 is not congruent and a vertex child straddles three cells of the level above:
    /// `parent0` and two of its neighbours. The engine lists three entries whether or not those
    /// neighbours exist, and so does this, with the null zone for one that does not.
    ///
    /// `a` must satisfy [`HexA3::is_valid_address`], which `Grid` checks first, so that the
    /// null zone and every other identifier with no cell have no parents. On any other address
    /// the answer is meaningless, though it does not panic: the arithmetic wraps as the eC's
    /// does, and where the eC divides by `POW3(-1)`, which is 0, this answers the null zone.
    fn parents(a: &Address) -> Option<Vec<Address>> {
        Some(
            parents(a.base)
                .into_iter()
                .map(|id| Address::new(id, &[]))
                .collect(),
        )
    }

    /// The zone's children, as the engine's `getZoneChildren` (`RI3H.ec:139-142`) lists them and
    /// in its order: the centroid child, whose centre is the zone's, then one child on each
    /// vertex, seven in all for a hexagon and six for a pentagon. A vertex child straddles three
    /// cells of the zone's level and lists all three among its parents.
    ///
    /// One deliberate departure: a zone at resolution 33, the finest, has none here. The
    /// engine answers seven identifiers of level 34, whose indices spill into the root and
    /// level fields, which its own parser refuses and among whose parents it divides by zero.
    ///
    /// `a` must satisfy [`HexA3::is_valid_address`], as for [`HexA3::parents`].
    fn children(a: &Address) -> Option<Vec<Address>> {
        Some(
            children(a.base)
                .into_iter()
                .map(|id| Address::new(id, &[]))
                .collect(),
        )
    }

    /// The engine's centroid parent (`getZoneCentroidParent`, `RI3H.ec:124-127`): the first of the
    /// zone's parents that is itself a centroid child, `parent0` if it is one and otherwise the
    /// first of the two others that is, or `None` where none is. See
    /// [`crate::Grid::centroid_parent`] for how far that is from the parent of which the zone is
    /// the centroid child.
    ///
    /// `a` must satisfy [`HexA3::is_valid_address`], as for [`HexA3::parents`].
    fn centroid_parent(a: &Address) -> Option<Option<Address>> {
        let p = centroid_parent(a.base);
        Some((p != NULL_ZONE).then(|| Address::new(p, &[])))
    }

    /// Whether the zone is a centroid child, as the engine's `isZoneCentroidChild`
    /// (`RI3H.ec:60-63`) answers: every odd-level sub-hex B, both poles, and the even-level zones
    /// whose row and column within their rhombus are both multiples of 3, or whose sum is, every
    /// zone of resolution 0 among them. Above resolution 0, each is the centroid child of
    /// `parent0`, its only parent.
    ///
    /// `a` must satisfy [`HexA3::is_valid_address`], as for [`HexA3::parents`].
    fn is_centroid_child(a: &Address) -> Option<bool> {
        Some(is_centroid_child(a.base))
    }

    /// The aperture-3 grids order their sub-zones, as `I3HSubZones.ec` generates them: row by
    /// row across the zone, from a vertex of it at an even depth and from an edge at an odd one.
    const HAS_SUB_ZONE_ORDER: bool = true;

    /// The number of sub-zones `depth` levels below the zone, as the engine's `countSubZones`
    /// (`RI3H.ec:48-51`, through `getSubZonesCount` at `RI3H.ec:2198-2202`) answers:
    /// `((3^d + 3^((d + 1) / 2) + 1) * n + 5) / 6` for a zone of `n` sides, 7 for a hexagon and 6
    /// for a pentagon at depth 1, and 1 at depth 0. `None` where the sub-zones would lie beyond
    /// resolution 33, which `Grid` refuses first.
    ///
    /// `a` must satisfy [`HexA3::is_valid_address`], which `Grid` checks first, here and in the
    /// two methods below.
    fn count_sub_zones(a: &Address, depth: u8) -> Option<u64> {
        hex_a3_subzones::count(a.base, depth)
    }

    /// The first sub-zone `depth` levels below the zone, in the engine's order, as its
    /// `getFirstSubZone` (`RI3H.ec:156-159`, `RI3H.ec:2204-2210`) answers at every depth but 0;
    /// `None` beyond resolution 33, as for [`HexA3::count_sub_zones`]. Where
    /// [`HexA3::sub_zones`] departs from the engine's order, this is the first of its own.
    ///
    /// One deliberate departure: at depth 0 this answers the zone itself, as `Grid` does before
    /// it asks. The engine quantises a vertex of the zone at the zone's own level there, and so
    /// names a neighbour for two zones in three, while its `getSubZones` names the zone itself.
    fn first_sub_zone(a: &Address, depth: u8) -> Option<Address> {
        hex_a3_subzones::first_sub_zone(a.base, depth).map(|id| Address::new(id, &[]))
    }

    /// Every sub-zone `depth` levels below the zone, in the engine's order, as its `getSubZones`
    /// (`dggrs.ec:195-222`, through `getSubZoneCRSCentroids` at `RI3H.ec:676-705`) answers in the
    /// 5x6 CRS: the sub-zone centroids, generated as the engine generates them, each quantised at
    /// the sub-zone level.
    ///
    /// `None` beyond resolution 33, as for [`HexA3::count_sub_zones`], and from `2^28`
    /// sub-zones, where the engine answers nothing either; `Grid` refuses both before it asks.
    ///
    /// One departure, at some zones on the edges of the rhombi, where a unit in the last place
    /// carries a point across an interruption of the layout that it should only reach, and the
    /// engine's order names zones twice, names zones that are not descendants of this one and
    /// omits descendants: this order is the zone's descendants, each once, in the generator's
    /// own scanline order. Every other order is the engine's.
    ///
    /// It also answers `None` above [`crate::grid::max_materialised_sub_zones`], well below the
    /// engine's own `2^28`: an `Address` costs about 32 bytes against a centroid's 16, and `Grid`
    /// already refuses a request above that cap before it ever calls here, so this only guards a
    /// caller that reaches this trait method directly, which would otherwise be free to
    /// materialise gigabytes of addresses.
    fn sub_zones(a: &Address, depth: u8) -> Option<Vec<Address>> {
        let count = hex_a3_subzones::count(a.base, depth)?;
        if count > crate::grid::max_materialised_sub_zones() {
            return None;
        }
        hex_a3_subzones::sub_zones(a.base, depth)
            .map(|ids| ids.into_iter().map(|id| Address::new(id, &[])).collect())
    }

    /// The sub-zone at `index`, `depth` levels below the zone, found through the generators' own
    /// index skip, without materialising [`HexA3::sub_zones`]'s whole sequence: a cost of the
    /// order of the number of scanlines, `3^(level / 2)` at the sub-zone's own level, rather than
    /// of `index`. Ports `getSubZoneAtIndex` (`RI3H.ec:212-227`).
    ///
    /// `None` beyond resolution 33 or at an `index` at or beyond the count, as
    /// [`HexA3::count_sub_zones`] answers `None` for the first and [`HexA3::sub_zones`] would
    /// answer nothing past the last for the second; a direct caller of this trait method must
    /// tell those apart from "no override" itself, as [`crate::Topology::sub_zone_at_index`]
    /// documents. `Grid` never meets either case, since it checks both before it asks.
    ///
    /// One deliberate departure, at depth 0, where this answers the zone itself: see
    /// [`HexA3::first_sub_zone`].
    fn sub_zone_at_index(a: &Address, depth: u8, index: u64) -> Option<Address> {
        hex_a3_subzones::sub_zone_at_index(a.base, depth, index).map(|id| Address::new(id, &[]))
    }

    /// `3^(level / 2)` cells along a rhombus: the I9R grid of the even level at or below
    /// `level`, whose index within a rhombus runs left to right and top to bottom
    /// (`RI3H.ec:858-861`). Nought beyond resolution 33.
    fn lattice_edge(level: u8) -> u64 {
        if level > MAX_LEVEL {
            return 0;
        }
        pow3(u64::from(level / 2))
    }

    /// The zones of one I9R cell: the cell's own zone at an even level, its sub-hexagons B, C
    /// and D at an odd one, as `listZones` names them (`RI3H.ec:779-786`).
    ///
    /// The identifier is its own key. `listZones` sorts its answer (`RI3H.ec:824`) by
    /// `I3HZone::OnCompare` (`RI3H.ec:865-876`), which compares the levels and then the values;
    /// at one level the value reads, from its highest field down, the root rhombus, the index
    /// `row * 3^(level / 2) + col` and the sub-hexagon, which is the walk rhombus by rhombus,
    /// row by row and cell by cell, and the polar roots 10 and 11 after the ten rhombi.
    fn cell_zones(level: u8, root: u8, row: u64, col: u64) -> Vec<(Address, u64)> {
        let p = Self::lattice_edge(level);
        if root > 9 || row >= p || col >= p {
            return Vec::new();
        }
        let sub_hexes: &[u64] = if level % 2 == 0 { &[0] } else { &[1, 2, 3] };
        sub_hexes
            .iter()
            .map(|&sub_hex| {
                let id = pack(
                    u64::from(level / 2),
                    u64::from(root),
                    row * p + col,
                    sub_hex,
                );
                (Address::new(id, &[]), id)
            })
            .collect()
    }

    /// The one zone of a polar root: the polar pentagon itself at an even level, its
    /// sub-hexagon B at an odd one, as `listZones` adds the two poles to every answer
    /// (`RI3H.ec:795-799`). The engine reads the sub-hexagons C and D of a polar root from text
    /// and lists neither; see [`HexA3::is_valid_address`].
    fn polar_zones(level: u8, root: u8) -> Vec<(Address, u64)> {
        if level > MAX_LEVEL || !(10..=11).contains(&root) {
            return Vec::new();
        }
        let id = pack(
            u64::from(level / 2),
            u64::from(root),
            0,
            u64::from(level % 2),
        );
        vec![(Address::new(id, &[]), id)]
    }

    /// The root, and the row and the column that the index within the root holds as
    /// `row * 3^levelI9R + col`, left to right and top to bottom (`RI3H.ec:858-861`). `None`
    /// for an address that [`HexA3::is_valid_address`] refuses.
    fn locate(a: &Address) -> Option<(u8, u64, u64, u64)> {
        if !Self::is_valid_address(a) {
            return None;
        }
        let (level_i9r, root, ix, _) = fields(a.base);
        let p = pow3(level_i9r);
        Some((root as u8, ix / p, ix % p, a.base))
    }

    /// DGGAL's compaction of a set of zones (`compactZones`, `RI3H.ec:161-186`), in ascending
    /// order of the identifiers: see [`compact`]. Every address must satisfy
    /// [`HexA3::is_valid_address`], which `Grid` checks first.
    fn compact(zones: &[Address]) -> Option<Vec<Address>> {
        let ids: Vec<u64> = zones.iter().map(|a| a.base).collect();
        Some(
            compact(&ids)
                .into_iter()
                .map(|id| Address::new(id, &[]))
                .collect(),
        )
    }
}

#[cfg(test)]
mod tests {
    //! Every expected coordinate and identifier below is the shipped engine's own answer, DGGAL
    //! v0.0.6 on ISEA3H, read through `getZoneFromCRSCentroid`, `getZoneCRSCentroid` and
    //! `getZoneCRSVertices` in the 5x6 CRS (`CRS(ogc, 153456)`), `getZoneFromWGS84Centroid` and
    //! `getZoneNeighbors`, and written out as bit patterns; the comparisons are bit for bit, and
    //! neighbour lists are compared in the engine's order.
    use super::*;
    use crate::indexings::I3h;
    use crate::projections::{Isea, Ivea, Rtea};
    use crate::{Error, GeoPoint, Grid, GridConfig, Indexing, Projection, ZoneId};

    fn bits(p: &PlanarPoint) -> (u64, u64) {
        (p.x.to_bits(), p.y.to_bits())
    }

    fn quantize(res: u8, x: u64, y: u64) -> Option<u64> {
        let p = PlanarPoint {
            face: -1,
            x: f64::from_bits(x),
            y: f64::from_bits(y),
        };
        HexA3::quantize(p, res).map(|a| a.base)
    }

    #[test]
    fn the_folded_constants_match_the_engines_rodata() {
        assert_eq!(T13.to_bits(), 0x3fd5_5555_5555_5555); // 0x49a98
        assert_eq!(T23.to_bits(), 0x3fe5_5555_5555_5555); // 0x49a90
        assert_eq!(MINUS_T13.to_bits(), 0xbfd5_5555_5555_5555); // 0x49f60
        assert_eq!(TEN_MINUS_1E11.to_bits(), 0x4023_ffff_ffff_ea03); // 0x49f58
        assert_eq!(FIVE_MINUS_1E11.to_bits(), 0x4013_ffff_ffff_d405); // 0x49f48
        assert_eq!(ONE_PLUS_1E11.to_bits(), 0x3ff0_0000_0000_afec); // 0x4a010
        assert_eq!(ONE_MINUS_1E11.to_bits(), 0x3fef_ffff_fffe_a028); // 0x4a018
        assert_eq!(TWO_POW_63.to_bits(), 0x43e0_0000_0000_0000); // 0x49ef0
    }

    /// The first five points are ones at which the eC's order and the library's part company,
    /// so that each answer-changing site of `fromCentroid` has a point that fails without it.
    #[test]
    fn quantisation_matches_the_engine() {
        for (what, res, x, y, want) in [
            // The south-pole test, `fabs(c.y - c.x - 2)` (RI3H.ec:1339), compiled as
            // `fabs((c.y - 2) - c.x)` at `0xc820` to `0xc838`: the eC's grouping answers the
            // null zone.
            (
                "the south-pole test, RI3H.ec:1339, at 0xc820",
                0,
                0x3fef_ffff_ffff_dcd1,
                0x4008_0000_0003_66d0,
                Some(0x0160_0000_0000_0000),
            ),
            // `next(fl(2/3))`, the one double at which `3*xd > 2` (RI3H.ec:1436) and its compiled
            // form `xd > fl(2/3)` (`0xd138`) differ, and the diagonal with it, `3 * (xd - yd) > 2`
            // (RI3H.ec:1518), compiled as `(y + xd) - Y > fl(2/3)` (`0xd1be` to `0xd1ca`); the
            // eC's order answers `A0-0-C`.
            (
                "the right third and diagonal, RI3H.ec:1436 and :1518, at 0xd138 and 0xd1be",
                1,
                0x3fe5_5555_5555_5556,
                0x8000_0000_0000_0001,
                Some(0x0140_0000_0000_0001),
            ),
            // The same double on the other axis: `3 * (yd - xd) > 2` (RI3H.ec:1487), compiled as
            // `(x + yd) - X > fl(2/3)` (`0xd13e` to `0xd14a`); the eC's order answers `A0-0-D`.
            (
                "the bottom third and diagonal, RI3H.ec:1487, at 0xd13e",
                1,
                0x8000_0000_0000_0001,
                0x3fe5_5555_5555_5556,
                Some(0x0020_0000_0000_0001),
            ),
            // The top-dent fold, `c.x += (cy+1 - c.y)` (RI3H.ec:1351), compiled as
            // `(c.x - c.y) + (cy + 1)` at `0xcad1` to `0xcaf0`.
            (
                "the top-dent fold, RI3H.ec:1351, at 0xcad1",
                3,
                0x3fef_ffff_ffff_ffff,
                0x3fdc_71c7_1c71_c71b,
                Some(0x0240_0000_0000_0009),
            ),
            // The bottom-dent fold, `c.y += (cx+1 - c.x)` (RI3H.ec:1357), compiled as
            // `(c.y - c.x) + (cx + 1)` at `0xcecb` to `0xcee9`.
            (
                "the bottom-dent fold, RI3H.ec:1357, at 0xcecb",
                13,
                0x3fac_8fbe_72f2_c5d2,
                0x4000_0000_0000_0000,
                Some(0x0c60_0000_001e_9cc1),
            ),
            (
                "resolution 0",
                0,
                0x4002_6666_6666_6666,
                0x4005_9999_9999_999a,
                Some(0x00a0_0000_0000_0000),
            ),
            (
                "resolution 33",
                33,
                0x4002_6666_6666_6666,
                0x4005_9999_9999_999a,
                Some(0x2092_6edf_a839_bd03),
            ),
            (
                "resolution 9",
                9,
                0x4013_ffff_ffff_fb9a,
                0x4011_0000_0000_0000,
                Some(0x0800_0000_0000_00f5),
            ),
            // The engine answers `BA-2-A` and `BB-6-A`, which its parser refuses: the departure.
            (
                "root 10",
                2,
                0x4013_ffff_ffee_d1f4,
                0x4012_0000_0000_0000,
                None,
            ),
            (
                "root 11",
                2,
                0x4011_5555_5555_5555,
                0x4017_ffff_f94a_035a,
                None,
            ),
        ] {
            assert_eq!(quantize(res, x, y), want, "{what}");
        }
        // Far outside the layout: `(int)floor(-inf)` is `i32::MAX` as gcc inlined it, and
        // `move5x6Vertex`'s integers are 32 bits wide, so that `-1e15` overflows them.
        let far = (-5.0f64).to_bits();
        assert_eq!(
            quantize(0, far, f64::NEG_INFINITY.to_bits()),
            Some(0x0140_0000_0000_0000)
        );
        assert_eq!(
            quantize(20, far, (-1e15f64).to_bits()),
            Some(0x1400_0000_0000_0004)
        );
        // The bottom-dent fold leaves this point more than a cell left of its rhombus.
        assert_eq!(quantize(2, (-3.0f64).to_bits(), 0.5f64.to_bits()), None);
        // `getZoneFromCRSCentroid` quantises at no level above 33.
        assert_eq!(quantize(34, 2.3f64.to_bits(), 2.7f64.to_bits()), None);
    }

    /// A pentagon of each parity, both polar pentagons, a cell of each kind across an
    /// interruption, whose vertices `canonicalize5x6` moves, the finest resolution, and three
    /// zones at which the eC's order and the library's part company: the centroid at the first
    /// two, on a sub-hex D and a C, where the `centroid` property's `2*d/3` and `d/3`
    /// (`RI3H.ec:2163`, `:2165`) are compiled as products with `fl(2/3)` and `fl(1/3)` (`0xaabf`,
    /// `0xaae0`); the vertices at the third, where `getVertices` (`RI3H.ec:1699-1786`) forms its
    /// thirds of `d` in the same way (`0xd9df` and the sites beside it). Each table lists the
    /// vertices in the engine's order, which the faithful `getVertices` gives; `planar_vertices`
    /// gives them anticlockwise, reversed from the engine's first vertex save at the two polar
    /// pentagons, and a pentagon of roots 0 to 9, a hexagon and a zone across an interruption
    /// are among them.
    #[test]
    #[allow(clippy::type_complexity, reason = "the table of engine answers")]
    fn centroids_and_vertices_match_the_engine() {
        let zones: [(u64, (u64, u64), &[(u64, u64)]); 11] = [
            // BA-0-A, the north polar pentagon, resolution 2.
            (
                0x0340_0000_0000_0000,
                (0x3ff0_0000_0000_0000, 0x0000_0000_0000_0000),
                &[
                    (0x3fec_71c7_1c71_c71c, 0x3fbc_71c7_1c71_c720),
                    (0x3ffe_38e3_8e38_e38e, 0x3ff1_c71c_71c7_1c72),
                    (0x4007_1c71_c71c_71c7, 0x4000_e38e_38e3_8e39),
                    (0x400f_1c71_c71c_71c7, 0x4008_e38e_38e3_8e39),
                    (0x4013_8e38_e38e_38e4, 0x4010_71c7_1c71_c71c),
                ],
            ),
            // CB-0-B, the south polar pentagon, resolution 5.
            (
                0x0560_0000_0000_0001,
                (0x4010_0000_0000_0000, 0x4018_0000_0000_0000),
                &[
                    (0xbcd0_0000_0000_0000, 0x3fff_684b_da12_f684),
                    (0x400f_ffff_ffff_fffe, 0x4017_da12_f684_bda1),
                    (0x4007_ffff_ffff_fffe, 0x4013_da12_f684_bda1),
                    (0x3fff_ffff_ffff_fffc, 0x400f_b425_ed09_7b42),
                    (0x3fef_ffff_ffff_fff8, 0x4007_b425_ed09_7b42),
                ],
            ),
            // C4-0-A, a pentagon, resolution 4.
            (
                0x0480_0000_0000_0000,
                (0x4000_0000_0000_0000, 0x4000_0000_0000_0000),
                &[
                    (0x4000_97b4_25ed_097b, 0x4000_4bda_12f6_84be),
                    (0x4000_4bda_12f6_84be, 0x4000_97b4_25ed_097b),
                    (0x3fff_684b_da12_f685, 0x4000_4bda_12f6_84be),
                    (0x3ffe_d097_b425_ed0a, 0x3fff_684b_da12_f685),
                    (0x3fff_684b_da12_f685, 0x3ffe_d097_b425_ed0a),
                ],
            ),
            // D7-0-B, a pentagon of a southern rhombus, resolution 7.
            (
                0x06e0_0000_0000_0001,
                (0x4008_0000_0000_0000, 0x4010_0000_0000_0000),
                &[
                    (0x4008_1948_b0fc_d6ea, 0x4010_0000_0000_0000),
                    (0x4008_1948_b0fc_d6ea, 0x4010_0ca4_587e_6b75),
                    (0x4008_0000_0000_0000, 0x4010_0ca4_587e_6b75),
                    (0x4007_e6b7_4f03_2916, 0x400f_e6b7_4f03_2916),
                    (0x4008_0000_0000_0000, 0x400f_e6b7_4f03_2916),
                ],
            ),
            // C1-3F-A, across an interruption, resolution 4.
            (
                0x0420_0000_0000_00fc,
                (0x0000_0000_0000_0000, 0x3ffc_71c7_1c71_c71c),
                &[
                    (0x3fb2_f684_bda1_2f68, 0x3ffd_097b_425e_d097),
                    (0x3fa2_f684_bda1_2f68, 0x3ffd_a12f_684b_da12),
                    (0x4010_97b4_25ed_097c, 0x4017_da12_f684_bda1),
                    (0x4010_bda1_2f68_4bda, 0x4017_b425_ed09_7b42),
                    (0x4011_097b_425e_d098, 0x4017_da12_f684_bda1),
                    (0x3fa2_f684_bda1_2f68, 0x3ffb_da12_f684_bda1),
                ],
            ),
            // C5-4D-D, across an interruption, resolution 5.
            (
                0x04a0_0000_0000_0137,
                (0x4004_bda1_2f68_4bda, 0x400f_b425_ed09_7b42),
                &[
                    (0x4004_71c7_1c71_c71c, 0x400f_684b_da12_f685),
                    (0x4004_bda1_2f68_4bda, 0x400f_684b_da12_f685),
                    (0x4005_097b_425e_d097, 0x400f_b425_ed09_7b42),
                    (0x4008_0000_0000_0000, 0x4011_7b42_5ed0_97b4),
                    (0x4008_0000_0000_0000, 0x4011_a12f_684b_da13),
                    (0x4004_71c7_1c71_c71c, 0x400f_b425_ed09_7b42),
                ],
            ),
            // D0-BC-C, across an interruption, resolution 7.
            (
                0x0600_0000_0000_02f2,
                (0x3fef_9add_3c0c_a458, 0x3fce_0652_2c3f_35ba),
                &[
                    (0x3fef_35ba_7819_48b0, 0x3fcc_71c7_1c71_c71c),
                    (0x3fef_9add_3c0c_a458, 0x3fcc_71c7_1c71_c71c),
                    (0x3ffc_3f35_ba78_1949, 0x3ff0_0000_0000_0000),
                    (0x3ffc_0ca4_587e_6b75, 0x3ff0_0000_0000_0000),
                    (0x3fef_9add_3c0c_a458, 0x3fcf_9add_3c0c_a458),
                    (0x3fef_35ba_7819_48b0, 0x3fce_0652_2c3f_35ba),
                ],
            ),
            // G0-63226-D, resolution 13: the sub-hex D centroid, `RI3H.ec:2165`, at `0xaabf`.
            (
                0x0c00_0000_0018_c89b,
                (0x3f5d_f756_80fe_b65e, 0x3fe8_7aac_8a20_12ba),
                &[
                    (0x3f56_7980_e0bf_08c7, 0x3fe8_76ed_9f4f_f2e3),
                    (0x3f5d_f756_80fe_b65e, 0x3fe8_76ed_9f4f_f2e3),
                    (0x3f62_ba96_109f_31fb, 0x3fe8_7aac_8a20_12ba),
                    (0x3f62_ba96_109f_31fb, 0x3fe8_7e6b_74f0_3290),
                    (0x3f5d_f756_80fe_b65e, 0x3fe8_7e6b_74f0_3290),
                    (0x3f56_7980_e0bf_08c7, 0x3fe8_7aac_8a20_12ba),
                ],
            ),
            // G0-92-C, resolution 13: the sub-hex C centroid, `RI3H.ec:2163`, at `0xaae0`.
            (
                0x0c00_0000_0000_024a,
                (0x3fc9_c08e_56da_e4b9, 0x3f3d_f756_80fe_b65e),
                &[
                    (0x3fc9_b192_ab9a_655e, 0x0000_0000_0000_0000),
                    (0x3fc9_c08e_56da_e4b9, 0x0000_0000_0000_0000),
                    (0x3fc9_cf8a_021b_6415, 0x3f3d_f756_80fe_b65e),
                    (0x3fc9_cf8a_021b_6415, 0x3f4d_f756_80fe_b65e),
                    (0x3fc9_c08e_56da_e4b9, 0x3f4d_f756_80fe_b65e),
                    (0x3fc9_b192_ab9a_655e, 0x3f3d_f756_80fe_b65e),
                ],
            ),
            // G2-2D8-A, resolution 12: the vertices, `RI3H.ec:1699-1786`, from `0xd9df`.
            (
                0x0c40_0000_0000_0b60,
                (0x3fff_fa61_9fc7_d03d, 0x3fef_ffff_ffff_ffff),
                &[
                    (0x3fff_fe20_8a97_f014, 0x3ff0_01df_7568_0feb),
                    (0x3fff_fc41_152f_e028, 0x3ff0_03be_ead0_1fd6),
                    (0x3fff_f882_2a5f_c052, 0x3ff0_01df_7568_0feb),
                    (0x3fef_fc41_152f_e029, 0x3f5d_f756_80fe_b998),
                    (0x3fef_f882_2a5f_c052, 0x3f4d_f756_80fe_bcd2),
                    (0x3fef_fc41_152f_e029, 0x3f3d_f756_80fe_c344),
                ],
            ),
            // Q4-49BB7EA0E6F40-D, resolution 33.
            (
                0x2092_6edf_a839_bd03,
                (0x4002_6666_6681_01ad, 0x4005_9999_997e_fe52),
                &[
                    (0x4002_6666_6576_f0ea, 0x4005_9999_9874_ed90),
                    (0x4002_6666_6681_01ad, 0x4005_9999_9874_ed90),
                    (0x4002_6666_678b_126f, 0x4005_9999_997e_fe52),
                    (0x4002_6666_678b_126f, 0x4005_9999_9a89_0f15),
                    (0x4002_6666_6681_01ad, 0x4005_9999_9a89_0f15),
                    (0x4002_6666_6576_f0ea, 0x4005_9999_997e_fe52),
                ],
            ),
        ];
        for (id, centroid, engine) in zones {
            let a = Address::new(id, &[]);
            assert!(HexA3::is_valid_address(&a), "{id:#018x}");
            assert_eq!(bits(&HexA3::planar_centroid(&a)), centroid, "{id:#018x}");
            // `getVertices` in the engine's order, canonicalised as the engine passes them on.
            let faithful: Vec<_> = vertices(id)
                .into_iter()
                .map(|(x, y)| canonicalize5x6(x, y))
                .map(|(x, y)| (x.to_bits(), y.to_bits()))
                .collect();
            assert_eq!(faithful, engine, "{id:#018x}");
            // `planar_vertices` turns every ring but a polar pentagon's the other way from the
            // engine's first vertex: `[v0, v1, ..., v(n-1)]` becomes `[v0, v(n-1), ..., v1]`.
            let mut anticlockwise = engine.to_vec();
            if (id >> 53) & 0xF < 10 {
                anticlockwise[1..].reverse();
            }
            let ours: Vec<_> = HexA3::planar_vertices(&a).iter().map(bits).collect();
            assert_eq!(ours, anticlockwise, "{id:#018x}");
        }
    }

    fn isea3h() -> Grid<Isea, HexA3, I3h> {
        Grid::new(GridConfig::default(), "ISEA3H").unwrap()
    }

    /// A hexagon of each parity, the finest resolution among them; a pentagon of each parity,
    /// the even one merging `topRight` into `topLeft`; an even pentagon of an odd root, which
    /// merges `bottomRight` into `bottomLeft`; both polar pentagons; three zones on an
    /// interruption, which step across it and, at the south one, take the extra `top`
    /// neighbour; the odd `right` stand-in at a north interruption; and a hexagon on an
    /// interruption whose `topRight` is merged into its `topLeft`. Each list is the engine's,
    /// in its order.
    #[test]
    fn neighbours_match_the_engine_in_its_order() {
        let grid = isea3h();
        let zones: [(u64, &[u64]); 11] = [
            // G2-2D8-A, resolution 12: G2-2D7-A G0-2D8-A G4-2D8-A GA-0-A G0-5B1-A G2-5B1-A.
            (
                0x0c40_0000_0000_0b60,
                &[
                    0x0c40_0000_0000_0b5c,
                    0x0c00_0000_0000_0b60,
                    0x0c80_0000_0000_0b60,
                    0x0d40_0000_0000_0000,
                    0x0c00_0000_0000_16c4,
                    0x0c40_0000_0000_16c4,
                ],
            ),
            // Q4-49BB7EA0E6F40-D, resolution 33: Q4-49BB7EA0E6F40-C Q4-49BB7EC9F4681-B
            // Q4-49BB7EA0E6F40-B Q4-49BB7EC9F4682-B Q4-49BB7EA0E6F3F-C Q4-49BB7EC9F4681-C.
            (
                0x2092_6edf_a839_bd03,
                &[
                    0x2092_6edf_a839_bd02,
                    0x2092_6edf_b27d_1a05,
                    0x2092_6edf_a839_bd01,
                    0x2092_6edf_b27d_1a09,
                    0x2092_6edf_a839_bcfe,
                    0x2092_6edf_b27d_1a06,
                ],
            ),
            // C4-0-A, a pentagon, resolution 4: C2-50-A C4-A-A C4-1-A C3-8-A C4-9-A.
            (
                0x0480_0000_0000_0000,
                &[
                    0x0440_0000_0000_0140,
                    0x0480_0000_0000_0028,
                    0x0480_0000_0000_0004,
                    0x0460_0000_0000_0020,
                    0x0480_0000_0000_0024,
                ],
            ),
            // C1-0-A, a pentagon of a southern rhombus, resolution 4: C9-50-A C1-A-A C0-48-A
            // C1-1-A C1-9-A.
            (
                0x0420_0000_0000_0000,
                &[
                    0x0520_0000_0000_0140,
                    0x0420_0000_0000_0028,
                    0x0400_0000_0000_0120,
                    0x0420_0000_0000_0004,
                    0x0420_0000_0000_0024,
                ],
            ),
            // D7-0-B, a pentagon of a southern rhombus, resolution 7: D6-2BE-D D5-2D8-C D7-0-C
            // D5-2D8-D D7-0-D.
            (
                0x06e0_0000_0000_0001,
                &[
                    0x06c0_0000_0000_0afb,
                    0x06a0_0000_0000_0b62,
                    0x06e0_0000_0000_0002,
                    0x06a0_0000_0000_0b63,
                    0x06e0_0000_0000_0003,
                ],
            ),
            // BA-0-A, the north polar pentagon, resolution 2: B4-2-A B6-2-A B8-2-A B0-2-A B2-2-A.
            (
                0x0340_0000_0000_0000,
                &[
                    0x0280_0000_0000_0008,
                    0x02c0_0000_0000_0008,
                    0x0300_0000_0000_0008,
                    0x0200_0000_0000_0008,
                    0x0240_0000_0000_0008,
                ],
            ),
            // CB-0-B, the south polar pentagon, resolution 5: C9-48-D C1-48-D C3-48-D C5-48-D
            // C7-48-D.
            (
                0x0560_0000_0000_0001,
                &[
                    0x0520_0000_0000_0123,
                    0x0420_0000_0000_0123,
                    0x0460_0000_0000_0123,
                    0x04a0_0000_0000_0123,
                    0x04e0_0000_0000_0123,
                ],
            ),
            // C1-3F-A, on a south interruption, resolution 4: C1-36-A C9-49-A C1-49-A C9-4A-A
            // C1-40-A C1-48-A.
            (
                0x0420_0000_0000_00fc,
                &[
                    0x0420_0000_0000_00d8,
                    0x0520_0000_0000_0124,
                    0x0420_0000_0000_0124,
                    0x0520_0000_0000_0128,
                    0x0420_0000_0000_0100,
                    0x0420_0000_0000_0120,
                ],
            ),
            // C5-4D-D, on an interruption, resolution 5: C5-4D-C C7-24-B C5-4D-B C7-1B-B
            // C5-4C-C C7-1B-D.
            (
                0x04a0_0000_0000_0137,
                &[
                    0x04a0_0000_0000_0136,
                    0x04e0_0000_0000_0091,
                    0x04a0_0000_0000_0135,
                    0x04e0_0000_0000_006d,
                    0x04a0_0000_0000_0132,
                    0x04e0_0000_0000_006f,
                ],
            ),
            // B0-1-B, on a north interruption, resolution 3: B8-8-C B0-0-C B8-5-C B0-1-C B8-5-D
            // B0-1-D.
            (
                0x0200_0000_0000_0005,
                &[
                    0x0300_0000_0000_0022,
                    0x0200_0000_0000_0002,
                    0x0300_0000_0000_0016,
                    0x0200_0000_0000_0006,
                    0x0300_0000_0000_0017,
                    0x0200_0000_0000_0007,
                ],
            ),
            // B0-2-A, on a north interruption, resolution 2: B0-1-A B8-2-A B2-2-A BA-0-A B8-5-A
            // B0-5-A.
            (
                0x0200_0000_0000_0008,
                &[
                    0x0200_0000_0000_0004,
                    0x0300_0000_0000_0008,
                    0x0240_0000_0000_0008,
                    0x0340_0000_0000_0000,
                    0x0300_0000_0000_0014,
                    0x0200_0000_0000_0014,
                ],
            ),
        ];
        for (id, want) in zones {
            let ours: Vec<u64> = grid.neighbors(ZoneId(id)).iter().map(|n| n.0).collect();
            assert_eq!(ours, want, "{id:#018x}");
        }
    }

    /// Whether `id` is the null zone or an identifier the indexing reads back from its text.
    fn null_or_readable(id: ZoneId) -> bool {
        id == ZoneId::NULL || I3h::from_text(&I3h::to_text(id)) == Ok(id)
    }

    /// Near both geographic poles, at resolutions 2 to 21, every answer is the null zone or an
    /// identifier that reads back. The sweep crosses the band of the departure, so some answers
    /// are the null zone, and two points found on it, on either meridian of the edge that runs
    /// over the north pole, are the null zone where the engine answers `BA-2-A` and `BA-1-A`.
    #[test]
    fn the_unreadable_north_pole_class_is_the_null_zone() {
        fn sweep<P: Projection>() -> usize {
            let grid: Grid<P, HexA3, I3h> = Grid::new(GridConfig::default(), "aperture 3").unwrap();
            let mut nulls = 0;
            for res in 2..=21 {
                for colatitude in [1e-4, 1e-7] {
                    for i in 0..360 {
                        let lon = f64::from(i) - 179.5;
                        for lat in [90.0 - colatitude, colatitude - 90.0] {
                            let id = grid.zone_from_geo(lat, lon, res).unwrap().id();
                            assert!(null_or_readable(id), "({lat}, {lon}) at {res}: {id:?}");
                            nulls += usize::from(id == ZoneId::NULL);
                        }
                    }
                }
            }
            nulls
        }
        assert!(sweep::<Isea>() > 0);
        assert!(sweep::<Ivea>() > 0);
        assert!(sweep::<Rtea>() > 0);
        let grid = isea3h();
        for (lat, lon) in [
            (84.822_824_872_311_86, 11.200_000_070_111_09),
            (70.745_192_359_828_7, -168.800_000_019_691_62),
        ] {
            assert_eq!(grid.zone_from_geo(lat, lon, 2).unwrap().id(), ZoneId::NULL);
        }
    }

    /// The 347 points of a sample of 680,000, planar points `1e-12` to `1e-3` from the lines of
    /// the layout lifted by ISEA's inverse at every resolution, at which the engine answers a
    /// root of 11 with a non-zero index: resolution, latitude and longitude, the shortest decimal
    /// of each double. The first is the engine's `BB-6-A`.
    const ROOT_11_SAMPLE: &str = "
        2 -27.777308938516747 -146.89252234152625   2 -12.080285213029388 -141.0522302271511
        2 -16.243130536790687 -142.4901088818944   2 -14.762602851426433 -141.97224327809596
        2 -6.616068330613575 -139.23267517938848   2 -0.881106367570185 -137.3675583692437
        2 -19.191769361322883 -143.54723308472714   2 -28.281645415154 -147.10507214908844
        2 -8.699165894485912 -139.9194115539424   2 -17.12314786036956 -142.8018090677669
        2 -31.198934994504143 -148.37928731423904   2 -19.088203425241055 -143.5094591044675
        2 -32.08068915740835 -148.78087699057457   2 -10.032067129039568 -140.3629374359994
        2 -28.518755870656864 -147.20574578929256   2 -7.089438871111131 -139.388128432186
        2 -35.337203410117326 -150.34096304777668   2 -36.65404054111939 -151.01051075442172
        2 -36.21959778165924 -150.78692872619064   2 -20.643962622895746 -144.082377559302
        2 -38.04409641541246 -151.74485432985367   2 -27.53479331715992 -146.79106466613183
        2 -35.58482641429763 -150.46501423782044   2 -3.1654875451167563 -138.10755517123408
        3 -2.4201768812405087 -137.8658641088481   3 -18.42490391481147 -143.26870048828718
        3 -23.771934940533914 -145.27372706526688   3 -29.81462636650685 -147.7647232999311
        3 -5.348838947077956 -138.81801028351322   3 -11.312461764994291 -140.79263421227623
        3 -19.631121195020008 -143.7080452037545   3 -37.51880186780795 -151.4638509714011
        3 -17.8825751791318 -143.0733109719046   3 -15.045570966591958 -142.07061449871458
        3 -0.8042578531832558 -137.34269530767557   3 -1.620664373221203 -137.60690475186993
        3 -1.2730156354664106 -137.4943831872062   3 -26.12655271052739 -146.21113134583507
        3 -0.675093311948068 -137.30090802471878   3 -23.972975830081687 -145.3523037685614
        3 -18.90949539442374 -143.4444062032475   3 -31.189903875034634 -148.37521618699583
        3 -11.431457553022755 -140.83276974551094   3 -30.842567968971167 -148.21928016220608
        3 -5.2007115516974665 -138.76967195595668   4 -18.125172425711053 -143.16055465960125
        4 -48.682734459437285 -158.67312958600385   4 -29.02724044069413 -147.42328556716762
        4 -21.5731786851937 -144.43047936168034   4 -49.979323262351734 -159.73940546909174
        4 -36.81112202367475 -151.09202178744994   4 -32.906461326541525 -149.16458214396013
        4 -40.66473672323775 -153.21624091214207   4 -44.08359927721403 -155.3396957961701
        4 -17.001958545716096 -142.7587027356167   4 -42.894073259436176 -154.5715120275598
        4 -49.107861749643426 -159.0157103446817   4 -16.210096226513194 -142.47847021749837
        4 -30.61974911855004 -148.11986742155293   4 -3.5404605262660707 -138.22928535455125
        4 -27.462777213104786 -146.76102885522178   4 -17.200400840891533 -142.82931876656386
        4 -17.558076942435378 -142.95700552988265   4 -22.0863336042455 -144.62474105567665
        4 -12.244299869084324 -141.1078831487883   4 -49.421932170721206 -159.2731308848764
        4 -51.855133407313325 -161.40485576419866   4 -16.252222366227237 -142.4933137542634
        4 -34.39990706482167 -149.8786272456542   4 -34.70401760148747 -150.02738237709255
        4 -11.688627934618342 -140.91962247241597   4 -27.162643480159414 -146.6363057989103
        4 -5.333689304592607 -138.8130639788036   5 -25.120961935540628 -145.80610392210485
        5 -49.628475833108816 -159.44448494515763   5 -16.424773427760368 -142.55419938316237
        5 -42.15970285701023 -154.11348412970762   5 -0.21012017462017168 -137.1504953565719
        5 -1.6603235533340244 -137.61974427573682   5 -43.77459786037566 -155.13689253953711
        5 -25.252191740285216 -145.85855165633274   5 -13.770674403287147 -141.6294589100536
        5 -18.0496596386593 -143.13337367003058   5 -46.58204179329902 -157.07172372141304
        5 -41.7331929713067 -153.85283065923028   5 -29.23309697314883 -147.51200682152853
        5 -23.288942937082567 -145.08599781072618   5 -27.575021156227486 -146.80786082077762
        5 -48.137066333726096 -158.24294611246515   5 -8.964957964070809 -140.0075688315893
        5 -48.51765863637428 -158.54188190276935   5 -52.449569960736994 -161.9666389863454
        5 -7.711898566492905 -139.59305071172153   5 -45.827382292862715 -156.53029958031166
        5 -23.711200016923566 -145.25004069787855   5 -9.159369861623931 -140.07214297330404
        5 -15.045277058489052 -142.07051403410645   6 -37.134682773616014 -151.26108756537172
        6 -56.36456410547509 -166.17403139432102   6 -44.894384803212695 -155.88324145610127
        6 -3.8426520641087425 -138.32746306188642   6 -32.85332691118151 -149.13966213547909
        6 -4.83244875795975 -138.64957791637838   6 -3.265117272726723 -138.13988918266713
        6 -37.26439389829997 -151.32930560589355   6 -9.559924782793965 -140.20541482324884
        6 -24.18755484407214 -145.43646024915589   6 -34.50193333252348 -149.9284020649387
        6 -49.95799225562589 -159.72133357658618   6 -22.137407617328797 -144.64415750386453
        6 -29.56662667517534 -147.65656830890822   6 -14.797357147178412 -141.984301844576
        6 -52.34770035779056 -161.8691088643057   6 -34.57466028446431 -149.96396338899575
        6 -47.35319363406426 -157.64279109487427   6 -3.5839415322415245 -138.24340739351663
        6 -49.34766589880095 -159.2119211006649   6 -30.587019703564103 -148.10530774511957
        6 -9.82287621968261 -140.29308579219307   6 -23.142835202547236 -145.0294950368223
        6 -45.41814837425578 -156.24354058491397   6 -12.704862627906952 -141.2645538400971
        6 -50.04802095581923 -159.79773279827106   7 -53.164342066304 -162.66615096401037
        7 -53.64016331184543 -163.14722668721004   7 -46.304380919974356 -156.87056983816782
        7 -49.438700908593425 -159.28698081687128   7 -12.354155573733246 -141.14520139192058
        7 -26.356872816128526 -146.30493810731429   7 -12.781250585600672 -141.29059639708643
        7 -50.83844118146436 -160.4829438523217   7 -51.69990794498755 -161.26098531084185
        7 -52.2037495672698 -161.7321864630761   7 -37.07964078401734 -151.2322164481956
        7 -49.5268299966686 -159.35994987243276   7 -25.53346125541689 -145.97137085774878
        7 -24.093611515011865 -145.39957942750692   7 -24.253797704796522 -145.4625010483957
        7 -36.27001307668915 -150.81271862131118   7 -45.87871025211096 -156.5665994448691
        7 -23.974347752521474 -145.3528409079973   7 -40.302128051745356 -153.00528632401125
        7 -51.98412159581232 -161.52526689165688   7 -13.982537580407154 -141.70240558129737
        7 -23.583596450341105 -145.20034914461488   7 -4.550718369612074 -138.55780054324543
        7 -43.69507791054715 -155.08507921977724   7 -6.76829012042182 -139.28262825410442
        8 -24.469839459519385 -145.54763086495703   8 -5.155676418805736 -138.75497179467325
        8 -8.013549717375604 -139.6925858249284   8 -49.32596513510222 -159.19407535382732
        8 -27.10947787769016 -146.61428212940422   8 -43.10915321509375 -154.70794555667266
        8 -42.38485109685742 -154.25264835637986   8 -55.32131004665728 -164.95586500416357
        8 -47.501226552170735 -157.75457109767692   8 -45.06243272349785 -155.9980365385631
        8 -15.483421026361045 -142.2233897527763   8 -4.211390584150817 -138.4473606696288
        8 -24.134415203175006 -145.4155912709395   8 -43.88204435252174 -155.20714536637365
        8 -48.91854909993002 -158.8623359714405   8 -55.683636620160776 -165.36999036420679
        8 -48.006068633903325 -158.14121164428317   8 -6.531839707734886 -139.20504888970464
        8 -27.61290661083796 -146.82369094410188   8 -30.320096578250517 -147.98695892875344
        8 -37.21853328531144 -151.30515719839715   9 -34.38193867426581 -149.86987539074164
        9 -17.482814705351917 -142.9300936281378   9 -5.536975983801087 -138.87944744417203
        9 -29.375244253122315 -147.57349325661252   9 -27.06798170113734 -146.59711078721344
        9 -51.00683283727417 -160.63238994374387   9 -12.801209447715614 -141.2974035581314
        9 -57.489154107302014 -167.58203977049814   9 -1.287075477243471 -137.49892619967542
        9 -3.7444496446724598 -138.2955506771238   9 -49.09435766102771 -159.00472494970057
        9 -37.97812115129874 -151.7093158171485   9 -9.870313794665273 -140.30891757414193
        9 -18.96751764740096 -143.46550348505943   9 -24.27770243057538 -145.4719053415501
        9 -52.527337233563884 -162.04144583449997   9 -34.8290308994375 -150.08887471055792
        9 -39.89064221665071 -152.76888550819228   9 -18.591046321055696 -143.32881787001676
        9 -6.500577230463197 -139.1947975607967   9 -50.1441238620757 -159.87964648901053
        9 -43.528272764053455 -154.97688735342837   9 -45.07761309201219 -156.00844369756652
        9 -28.55894853942546 -147.2228582942737   9 -14.407220277390255 -141.84905832533676
        9 -17.01183760553381 -142.76221419017813   9 -30.99296403045641 -148.28664496428755
        9 -47.885481821660385 -158.04808326597194   9 -50.71013493201859 -160.36990742139776
        10 -56.20817214185752 -165.9862869957093   10 -49.961853331491746 -159.72460326779043
        10 -4.547639704379901 -138.55679811146524   10 -24.07038171983388 -145.390468578964
        10 -31.20992556729815 -148.38424285155364   10 -11.292606170462815 -140.7859385538828
        10 -29.63810867888837 -147.68768366003408   10 -24.798654219560884 -145.67779753887808
        10 -2.5187824036599533 -137.897822193314   10 -4.371807474118398 -138.499557837692
        10 -17.503956940071564 -142.9376511201122   10 -49.725592247740174 -159.5256296904324
        10 -55.81033637924359 -165.5170038918907   10 -29.117077842916718 -147.46195692447995
        10 -45.5772488685416 -156.35446929943154   10 -35.26205720261352 -150.30345110087237
        10 -34.45517536297258 -149.90557446348456   10 -14.844546047025043 -142.000694807126
        10 -51.803244503503976 -161.3566255942401   10 -34.02547781012661 -149.69707281872886
        10 -39.30792842562556 -152.43937867879637   10 -15.826318482367467 -142.34351138523488
        10 -23.43959977620799 -145.14439827178447   10 -40.273255512615165 -152.98859651076305
        11 -35.28045076084947 -150.31262036650315   11 -9.309060507143565 -140.12190695762533
        11 -26.76244435266748 -146.47109704935482   11 -14.792946908710038 -141.98277014482042
        11 -11.602047471438269 -140.89036154288254   11 -6.284879480249349 -139.1241040899179
        11 -20.552152611240444 -144.04823099762325   11 -57.354915558659634 -167.40837336760907
        11 -12.453258059256772 -141.17889192160882   11 -56.63192555318788 -166.49939770277126
        11 -23.562325421694993 -145.19207592478818   11 -31.326759508776195 -148.4369990230212
        11 -40.172202367832405 -152.93030400103135   11 -51.4203336063458 -161.004759548723
        11 -50.81756354093646 -160.46450193320794   11 -41.956335235747645 -153.98871936607372
        11 -13.265839982969377 -141.45619672313364   11 -49.05322078867172 -158.97130458028187
        11 -31.264705002894484 -148.4089609211875   11 -46.228680657475266 -156.81612594221045
        11 -3.091847036623139 -138.08366012513278   11 -13.122796148842493 -141.40724255349932
        11 -42.62265322388713 -154.40083652629662   11 -8.740672743104783 -139.9331629514374
        11 -41.33398665806173 -153.61229646493328   11 -5.208712155485857 -138.77227567003234
        11 -49.400911452545216 -159.2557841897075   12 -56.248830296572876 -166.03491570415247
        12 -42.40891626111808 -154.26758804978735   12 -10.785955456359709 -140.61547489557336
        12 -0.34426549814547597 -137.19388566440634   12 -13.736442529558547 -141.61768571136213
        12 -27.948193654622532 -146.96430176003938   12 -6.9439823189492795 -139.34032381393146
        12 -16.868606788578955 -142.7113377684148   12 -1.0178781048114454 -137.41181238889098
        12 -46.981961963725674 -157.36556088522792   12 -1.907377543450194 -137.6997398748424
        12 -53.429678449381214 -162.93284604335983   12 -41.3864049699148 -153.64369413889992
        12 -19.45050819143759 -143.64182551474835   12 -28.172393800770116 -147.05884689320996
        12 -5.175798343065759 -138.7615365649812   12 -25.834738778658654 -146.09284380938382
        12 -23.970673702023934 -145.35140255810762   12 -11.139798439095976 -140.734460517414
        13 -36.374109561751695 -150.86612091333936   13 -36.81462613297243 -151.09384317655812
        13 -25.781164584238727 -146.07119507672155   13 -21.825496215757052 -144.52581181086242
        13 -13.955607262504966 -141.6931253422146   13 -31.422550413991523 -148.4803578045089
        13 -34.44622150945894 -149.90120600632025   13 -32.25442999698655 -148.86097353308497
        13 -18.97178051559491 -143.46705492894412   13 -47.18648186782412 -157.51775172689756
        13 -3.8320286210982046 -138.3240104188097   13 -45.074568768996386 -156.00635608278137
        13 -36.18962804452329 -150.77158579131216   13 -23.438523556206086 -145.14398042137594
        13 -24.526490455649963 -145.5700051705231   14 -52.01066639341805 -161.55014972020055
        14 -29.283112819694527 -147.533620212243   14 -13.364821742215653 -141.49010704294537
        14 -2.828856313527209 -137.9983514586373   14 -55.76654204357243 -165.4660568892826
        14 -51.16283346129372 -160.77196387363253   14 -9.810185529481142 -140.28885122599993
        15 -17.16839198515434 -142.8179174706039   15 -42.26874798386079 -154.18074778720217
        15 -33.158291313645584 -149.28313199521244   15 -54.360059991761865 -163.90005027644705
        15 -58.10297815680756 -168.39689213067365   15 -23.529888269785637 -145.17946518689834
        15 -23.944393047057837 -145.341116346468   15 -26.679051959469994 -146.4368291815967
        15 -48.29747856196324 -158.36831989921657   15 -9.924248232354143 -140.3269233820328
        15 -35.35838909316294 -150.3515231534977   15 -28.58820835673691 -147.23532457402672
        15 -57.51164408724341 -167.6112907369589   15 -48.291691006029716 -158.3637808008047
        16 -50.71327508513398 -160.37266529496824   16 -31.135142018518415 -148.35054810490726
        16 -10.141884422516869 -140.39964444247224   16 -53.16168391565476 -162.66349857055442
        16 -51.902630284857665 -161.44909222801732   16 -26.82657219442564 -146.4974852099962
        16 -30.16623240760651 -147.91905255468592   16 -29.28296490704945 -147.53355623948573
        16 -37.38114497779341 -151.39092680836833   16 -37.59063456003472 -151.50202090256928
        16 -6.309610414433272 -139.13220623406417   16 -32.66049783017761 -149.04949324841908
        16 -33.234840488674855 -149.3193131867565   16 -45.15883113197017 -156.06422898601116
        16 -43.2470867574418 -154.79600121777744   16 -41.40265477853297 -153.65343879883494
        16 -18.638013218033652 -143.3458349713197   16 -13.232825004036913 -141.44489255603088
        17 -29.806772210630566 -147.76128924155918   17 -28.312326783225764 -147.11807140283952
        17 -49.403915153269274 -159.2582618437504   17 -43.60854023569678 -155.02886663823583
        17 -56.67423583482331 -166.5514061634829   17 -9.91567991639461 -140.32406245456613
        17 -53.559647328759475 -163.0649206050478   17 -45.18870444050869 -156.08479245169465
        17 -38.57722670051069 -152.03456743158384   17 -38.09037011733147 -151.76981143026748
        18 -50.698403650610615 -160.35960792231535   18 -2.976751305604898 -138.04632014109708
        18 -17.834887800863914 -143.05619048889793   18 -40.13151697441302 -152.90688851672826
        18 -24.40561838988588 -145.52229260725625   19 -7.555096466058241 -139.54137044179322
        19 -28.98140291030597 -147.40358123946112   19 -3.6529719380559404 -138.2658304655989
        19 -48.594320007866266 -158.60270996041396   19 -50.890824332327895 -160.52929943477218
        19 -23.37383955663761 -145.11888959733014   19 -0.6215549172122506 -137.2835866019327
        20 -8.377470965032675 -139.8128747054684   21 -8.286382617793492 -139.7827447968157
        21 -20.983301861953734 -144.20896468485407
    ";

    /// At every point of the root-11 sample the engine answers an identifier it cannot read
    /// back; this answers the null zone at each, so that everything `zone_from_geo` answers
    /// across the sample is the null zone or an identifier that reads back.
    #[test]
    fn the_unreadable_root_11_answers_are_the_null_zone() {
        let grid = isea3h();
        let values: Vec<f64> = ROOT_11_SAMPLE
            .split_whitespace()
            .map(|v| v.parse().unwrap())
            .collect();
        assert_eq!(values.len(), 3 * 347);
        for p in values.chunks_exact(3) {
            let (res, lat, lon) = (p[0] as u8, p[1], p[2]);
            let id = grid.zone_from_geo(lat, lon, res).unwrap().id();
            assert!(null_or_readable(id), "({lat}, {lon}) at {res}: {id:?}");
            assert_eq!(id, ZoneId::NULL, "({lat}, {lon}) at {res}");
        }
        // The engine's `BB-6-A`, which names no cell and has no geometry.
        let id = grid
            .zone_from_geo(-27.777_308_938_516_747, -146.892_522_341_526_25, 2)
            .unwrap()
            .id();
        assert_eq!(id, ZoneId::NULL);
        assert!(!HexA3::is_valid_address(&I3h::decode(ZoneId(
            0x0360_0000_0000_0018
        ))));
    }

    /// Points `1e-7 / 3^(res / 2)` inside the seam `x = 5` of rhombus 8 and the seam `y = 6` of
    /// rhombus 9, lifted by each projection's inverse, at resolutions 2 to 21: every answer is the
    /// null zone or an identifier that reads back. The engine answers a polar root with a
    /// non-zero index at 256 to 288 of the 380 points of each seam on each grid, at every
    /// resolution from 2 to 17, so the departure answers the null zone at each of those.
    #[test]
    fn both_seams_answer_the_null_zone_or_a_readable_zone() {
        fn sweep<P: Projection>() {
            let grid: Grid<P, HexA3, I3h> = Grid::new(GridConfig::default(), "aperture 3").unwrap();
            let geom = P::build_geometry(&GridConfig::default());
            for res in 2..=21u8 {
                let delta = 1e-7 / 3u64.pow(u32::from(res / 2)) as f64;
                let mut nulls = [0; 2];
                for i in 1..20 {
                    let t = f64::from(i) / 20.0;
                    for (seam, (x, y)) in [(5.0 - delta, 4.0 + t), (4.0 + t, 6.0 - delta)]
                        .into_iter()
                        .enumerate()
                    {
                        let p = PlanarPoint { face: -1, x, y };
                        let g = P::inverse(&geom, p, res % 2 == 1);
                        let id = grid.zone_from_geo(g.lat, g.lon, res).unwrap().id();
                        assert!(null_or_readable(id), "({x}, {y}) at {res}: {id:?}");
                        nulls[seam] += usize::from(id == ZoneId::NULL);
                    }
                }
                if res <= 17 {
                    assert!(nulls[0] > 0 && nulls[1] > 0, "resolution {res}: {nulls:?}");
                }
            }
        }
        sweep::<Isea>();
        sweep::<Ivea>();
        sweep::<Rtea>();
    }

    /// The null zone, a phantom polar sub-hex (`AA-0-C`), the two polar-root classes (`BA-2-A`,
    /// `BB-6-A`), a level-34 identifier and 10,000 random values: every one that names no cell
    /// has no geometry and no relations, and none panics.
    #[test]
    fn identifiers_with_no_cell_have_no_geometry() {
        let grid = isea3h();
        let null_point = GeoPoint { lat: 0.0, lon: 0.0 };
        let mut ids = vec![
            ZoneId::NULL.0,
            pack(0, 10, 0, 2),
            0x0340_0000_0000_0008,
            0x0360_0000_0000_0018,
            pack(17, 0, 0, 0),
        ];
        assert!(ids.iter().all(|&raw| !is_readable(raw)));
        let mut state: u64 = 0xA3A3_A3A3;
        for _ in 0..10_000 {
            state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = state;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            ids.push(z ^ (z >> 31));
        }
        let mut readable = 0;
        for &raw in &ids {
            let id = ZoneId(raw);
            if is_readable(raw) {
                readable += 1;
                assert!(grid.centroid(id).lat.is_finite(), "{raw:#018x}");
                let n = grid.vertices(id).len();
                assert!(n == 5 || n == 6, "{raw:#018x}");
                continue;
            }
            assert_eq!(grid.centroid(id), null_point, "{raw:#018x}");
            assert_eq!(grid.vertices(id), vec![null_point; 6], "{raw:#018x}");
            assert!(grid.neighbors(id).is_empty());
            assert!(grid.parents(id).is_empty());
            assert!(grid.children(id).is_empty());
            assert!(matches!(
                grid.count_sub_zones(id, 1),
                Err(Error::InvalidZone(_))
            ));
        }
        assert!(
            readable > 0 && readable < 100,
            "{readable} random identifiers name a cell"
        );
    }

    /// `HexA3::neighbors` called directly, without `Grid`'s check, on the identifiers of
    /// `arbitrary_ids`, most of which name no cell. The answers are meaningless, but none
    /// panics: the index arithmetic wraps, as the eC's `uint64` does, where
    /// `3^level * 3^level` exceeds 64 bits from level 21.
    #[test]
    fn neighbours_of_any_address_do_not_panic() {
        for id in arbitrary_ids() {
            let ns = HexA3::neighbors(&Address::new(id, &[])).unwrap();
            assert!(ns.len() <= 8, "{id:#018x}");
        }
        // The quantiser at every level `neighbor` can pass it, up to 63, at points in the first
        // and last two rows and columns of cells of each rhombus, which reach every branch of
        // the index arithmetic.
        for level in 0..64 {
            let p = pow3(u64::from(level / 2)) as f64;
            for root in 0..10 {
                let (x0, y0) = (f64::from(root / 2), f64::from((root + 1) / 2));
                let at = |i: f64, t: f64| if i < 2.0 { i + t } else { p - 4.0 + i + t } / p;
                for (i, j) in (0..4).flat_map(|i| (0..4).map(move |j| (i, j))) {
                    for (s, t) in [
                        (0.1, 0.6),
                        (0.6, 0.1),
                        (0.05, 0.8),
                        (0.8, 0.05),
                        (0.5, 0.5),
                        (0.7, 0.7),
                        (0.1, 0.9),
                        (0.9, 0.1),
                    ] {
                        let (x, y) = (x0 + at(f64::from(i), s), y0 + at(f64::from(j), t));
                        let _ = from_centroid(level, x, y);
                    }
                }
            }
        }
    }

    /// A hexagon of each parity, both sub-hexes C and D among the odd ones; a pentagon of each
    /// parity; each polar pentagon at each parity, which reach the four closed formulas of
    /// `getChildren`; four zones on an interruption, one for each crossing test of
    /// `getParents`, `onBottomCrossingLeft` (`C1-9-A`), `onTopCrossingRightNegEpsilon`
    /// (`C0-C-A`), `onBottomCrossingLeftNegEpsilon` (`C1-1C-A`) and `onTopCrossingRight`
    /// (`C0-1-A`); a zone at resolution 32, whose children are at 33; a zone at resolution 33,
    /// which has none (the departure); and a zone at resolution 1. Each list is the engine's,
    /// in its order, read through `getZoneParents`, `getZoneChildren`, `getZoneCentroidParent`
    /// and `isZoneCentroidChild`.
    #[test]
    #[allow(clippy::type_complexity, reason = "the table of engine answers")]
    fn the_hierarchy_matches_the_engine() {
        let grid = isea3h();
        let zones: [(u64, bool, &[u64], &[u64], Option<u64>); 16] = [
            // G2-2D8-A, resolution 12.
            (
                0x0c40_0000_0000_0b60,
                false,
                // F2-F2-C F0-F2-C FA-0-B
                &[
                    0x0a40_0000_0000_03ca,
                    0x0a00_0000_0000_03ca,
                    0x0b40_0000_0000_0001,
                ],
                // G2-2D8-B G2-2D8-C G2-2D8-D G2-2D7-C G0-5B1-C G0-2D8-D G0-2D8-C
                &[
                    0x0c40_0000_0000_0b61,
                    0x0c40_0000_0000_0b62,
                    0x0c40_0000_0000_0b63,
                    0x0c40_0000_0000_0b5e,
                    0x0c00_0000_0000_16c6,
                    0x0c00_0000_0000_0b63,
                    0x0c00_0000_0000_0b62,
                ],
                Some(0x0b40_0000_0000_0001), // FA-0-B
            ),
            // G0-92-C, resolution 13.
            (
                0x0c00_0000_0000_024a,
                false,
                // G0-92-A G0-36C-A G0-93-A
                &[
                    0x0c00_0000_0000_0248,
                    0x0c00_0000_0000_0db0,
                    0x0c00_0000_0000_024c,
                ],
                // H0-A43-A H0-1B7-A H0-1B8-A H0-A44-A H0-12CF-A H0-12CE-A H0-A42-A
                &[
                    0x0e00_0000_0000_290c,
                    0x0e00_0000_0000_06dc,
                    0x0e00_0000_0000_06e0,
                    0x0e00_0000_0000_2910,
                    0x0e00_0000_0000_4b3c,
                    0x0e00_0000_0000_4b38,
                    0x0e00_0000_0000_2908,
                ],
                Some(0x0c00_0000_0000_024c), // G0-93-A
            ),
            // G0-63226-D, resolution 13.
            (
                0x0c00_0000_0018_c89b,
                false,
                // G0-63226-A G0-63500-A G0-634FF-A
                &[
                    0x0c00_0000_0018_c898,
                    0x0c00_0000_0018_d400,
                    0x0c00_0000_0018_d3fc,
                ],
                // H0-37D467-A H0-37CBDB-A H0-37CBDC-A H0-37D468-A H0-37DCF3-A H0-37DCF2-A
                // H0-37D466-A
                &[
                    0x0e00_0000_00df_519c,
                    0x0e00_0000_00df_2f6c,
                    0x0e00_0000_00df_2f70,
                    0x0e00_0000_00df_51a0,
                    0x0e00_0000_00df_73cc,
                    0x0e00_0000_00df_73c8,
                    0x0e00_0000_00df_5198,
                ],
                Some(0x0c00_0000_0018_c898), // G0-63226-A
            ),
            // C4-0-A, resolution 4.
            (
                0x0480_0000_0000_0000,
                true,
                // B4-0-B
                &[0x0280_0000_0000_0001],
                // C4-0-B C4-0-C C4-0-D C3-8-C C2-50-D C2-50-C
                &[
                    0x0480_0000_0000_0001,
                    0x0480_0000_0000_0002,
                    0x0480_0000_0000_0003,
                    0x0460_0000_0000_0022,
                    0x0440_0000_0000_0143,
                    0x0440_0000_0000_0142,
                ],
                Some(0x0280_0000_0000_0001), // B4-0-B
            ),
            // D7-0-B, resolution 7.
            (
                0x06e0_0000_0000_0001,
                true,
                // D7-0-A
                &[0x06e0_0000_0000_0000],
                // E7-0-A E7-1-A E7-52-A E7-51-A E5-19A0-A E6-1950-A
                &[
                    0x08e0_0000_0000_0000,
                    0x08e0_0000_0000_0004,
                    0x08e0_0000_0000_0148,
                    0x08e0_0000_0000_0144,
                    0x08a0_0000_0000_6680,
                    0x08c0_0000_0000_6540,
                ],
                Some(0x06e0_0000_0000_0000), // D7-0-A
            ),
            // BA-0-A, resolution 2.
            (
                0x0340_0000_0000_0000,
                true,
                // AA-0-B
                &[0x0140_0000_0000_0001],
                // BA-0-B B0-2-C B2-2-C B4-2-C B6-2-C B8-2-C
                &[
                    0x0340_0000_0000_0001,
                    0x0200_0000_0000_000a,
                    0x0240_0000_0000_000a,
                    0x0280_0000_0000_000a,
                    0x02c0_0000_0000_000a,
                    0x0300_0000_0000_000a,
                ],
                Some(0x0140_0000_0000_0001), // AA-0-B
            ),
            // CB-0-B, resolution 5.
            (
                0x0560_0000_0000_0001,
                true,
                // CB-0-A
                &[0x0560_0000_0000_0000],
                // DB-0-A D1-2BE-A D3-2BE-A D5-2BE-A D7-2BE-A D9-2BE-A
                &[
                    0x0760_0000_0000_0000,
                    0x0620_0000_0000_0af8,
                    0x0660_0000_0000_0af8,
                    0x06a0_0000_0000_0af8,
                    0x06e0_0000_0000_0af8,
                    0x0720_0000_0000_0af8,
                ],
                Some(0x0560_0000_0000_0000), // CB-0-A
            ),
            // CB-0-A, resolution 4.
            (
                0x0560_0000_0000_0000,
                true,
                // BB-0-B
                &[0x0360_0000_0000_0001],
                // CB-0-B C1-48-D C3-48-D C5-48-D C7-48-D C9-48-D
                &[
                    0x0560_0000_0000_0001,
                    0x0420_0000_0000_0123,
                    0x0460_0000_0000_0123,
                    0x04a0_0000_0000_0123,
                    0x04e0_0000_0000_0123,
                    0x0520_0000_0000_0123,
                ],
                Some(0x0360_0000_0000_0001), // BB-0-B
            ),
            // CA-0-B, resolution 5.
            (
                0x0540_0000_0000_0001,
                true,
                // CA-0-A
                &[0x0540_0000_0000_0000],
                // DA-0-A D0-1A-A D2-1A-A D4-1A-A D6-1A-A D8-1A-A
                &[
                    0x0740_0000_0000_0000,
                    0x0600_0000_0000_0068,
                    0x0640_0000_0000_0068,
                    0x0680_0000_0000_0068,
                    0x06c0_0000_0000_0068,
                    0x0700_0000_0000_0068,
                ],
                Some(0x0540_0000_0000_0000), // CA-0-A
            ),
            // C1-9-A, resolution 4.
            (
                0x0420_0000_0000_0024,
                false,
                // B1-0-B B1-0-D B9-8-D
                &[
                    0x0220_0000_0000_0001,
                    0x0220_0000_0000_0003,
                    0x0320_0000_0000_0023,
                ],
                // C1-9-B C1-9-C C1-9-D C9-4F-D C9-4F-C C9-50-D C1-0-D
                &[
                    0x0420_0000_0000_0025,
                    0x0420_0000_0000_0026,
                    0x0420_0000_0000_0027,
                    0x0520_0000_0000_013f,
                    0x0520_0000_0000_013e,
                    0x0520_0000_0000_0143,
                    0x0420_0000_0000_0003,
                ],
                Some(0x0220_0000_0000_0001), // B1-0-B
            ),
            // C0-C-A, resolution 4.
            (
                0x0400_0000_0000_0030,
                false,
                // B0-1-B B0-1-D B0-0-C
                &[
                    0x0200_0000_0000_0005,
                    0x0200_0000_0000_0007,
                    0x0200_0000_0000_0002,
                ],
                // C0-C-B C0-C-C C0-C-D C0-B-C C0-2-D C0-2-C C0-3-D
                &[
                    0x0400_0000_0000_0031,
                    0x0400_0000_0000_0032,
                    0x0400_0000_0000_0033,
                    0x0400_0000_0000_002e,
                    0x0400_0000_0000_000b,
                    0x0400_0000_0000_000a,
                    0x0400_0000_0000_000f,
                ],
                Some(0x0200_0000_0000_0005), // B0-1-B
            ),
            // C1-1C-A, resolution 4.
            (
                0x0420_0000_0000_0070,
                false,
                // B1-3-B B1-3-C B1-0-D
                &[
                    0x0220_0000_0000_000d,
                    0x0220_0000_0000_000e,
                    0x0220_0000_0000_0003,
                ],
                // C1-1C-B C1-1C-C C1-1C-D C1-1B-C C1-12-D C1-12-C C1-13-D
                &[
                    0x0420_0000_0000_0071,
                    0x0420_0000_0000_0072,
                    0x0420_0000_0000_0073,
                    0x0420_0000_0000_006e,
                    0x0420_0000_0000_004b,
                    0x0420_0000_0000_004a,
                    0x0420_0000_0000_004f,
                ],
                Some(0x0220_0000_0000_000d), // B1-3-B
            ),
            // C0-1-A, resolution 4.
            (
                0x0400_0000_0000_0004,
                false,
                // B0-0-B B0-0-C B8-8-C
                &[
                    0x0200_0000_0000_0001,
                    0x0200_0000_0000_0002,
                    0x0300_0000_0000_0022,
                ],
                // C0-1-B C0-1-C C0-1-D C0-0-C C8-50-C C8-47-D C8-47-C
                &[
                    0x0400_0000_0000_0005,
                    0x0400_0000_0000_0006,
                    0x0400_0000_0000_0007,
                    0x0400_0000_0000_0002,
                    0x0500_0000_0000_0142,
                    0x0500_0000_0000_011f,
                    0x0500_0000_0000_011e,
                ],
                Some(0x0200_0000_0000_0001), // B0-0-B
            ),
            // Q4-49BB7EA0E6F40-D, resolution 33.
            (
                0x2092_6edf_a839_bd03,
                false,
                // Q4-49BB7EA0E6F40-A Q4-49BB7EC9F4682-A Q4-49BB7EC9F4681-A
                &[
                    0x2092_6edf_a839_bd00,
                    0x2092_6edf_b27d_1a08,
                    0x2092_6edf_b27d_1a04,
                ],
                // The engine answers seven level-34 identifiers, which do not read back:
                // R5-179774548A487-A R5-179773D961EC3-A R5-179773D961EC4-A R5-179774548A488-A
                // R5-179774CFB2A4B-A R5-179774CFB2A4A-A R5-179774548A486-A
                &[],
                Some(0x2092_6edf_b27d_1a04), // Q4-49BB7EC9F4681-A
            ),
            // Q4-49BB7EA0E6F40-A, resolution 32.
            (
                0x2092_6edf_a839_bd00,
                false,
                // P4-83146EF0C30E-D P4-83146EF0C30D-C P4-83146FCBB579-B
                &[
                    0x1e82_0c51_bbc3_0c3b,
                    0x1e82_0c51_bbc3_0c36,
                    0x1e82_0c51_bf2e_d5e5,
                ],
                // Q4-49BB7EA0E6F40-B Q4-49BB7EA0E6F40-C Q4-49BB7EA0E6F40-D Q4-49BB7EA0E6F3F-C
                // Q4-49BB7E77D97FE-D Q4-49BB7E77D97FE-C Q4-49BB7E77D97FF-D
                &[
                    0x2092_6edf_a839_bd01,
                    0x2092_6edf_a839_bd02,
                    0x2092_6edf_a839_bd03,
                    0x2092_6edf_a839_bcfe,
                    0x2092_6edf_9df6_5ffb,
                    0x2092_6edf_9df6_5ffa,
                    0x2092_6edf_9df6_5fff,
                ],
                Some(0x1e82_0c51_bf2e_d5e5), // P4-83146FCBB579-B
            ),
            // A0-0-C, resolution 1.
            (
                0x0000_0000_0000_0002,
                false,
                // A0-0-A A2-0-A AA-0-A
                &[
                    0x0000_0000_0000_0000,
                    0x0040_0000_0000_0000,
                    0x0140_0000_0000_0000,
                ],
                // B0-5-A B0-1-A B0-2-A B2-2-A B2-1-A B0-8-A B0-4-A
                &[
                    0x0200_0000_0000_0014,
                    0x0200_0000_0000_0004,
                    0x0200_0000_0000_0008,
                    0x0240_0000_0000_0008,
                    0x0240_0000_0000_0004,
                    0x0200_0000_0000_0020,
                    0x0200_0000_0000_0010,
                ],
                Some(0x0000_0000_0000_0000), // A0-0-A
            ),
        ];
        for (id, is_centroid_child, parents, children, centroid_parent) in zones {
            let z = ZoneId(id);
            let ids = |zs: Vec<ZoneId>| zs.into_iter().map(|p| p.0).collect::<Vec<_>>();
            assert_eq!(grid.is_centroid_child(z), is_centroid_child, "{id:#018x}");
            assert_eq!(ids(grid.parents(z)), parents, "{id:#018x}");
            assert_eq!(ids(grid.children(z)), children, "{id:#018x}");
            assert_eq!(
                grid.centroid_parent(z).map(|p| p.0),
                centroid_parent,
                "{id:#018x}"
            );
        }
    }

    /// Identifiers the hierarchy can meet: every I9R level the five-bit field can hold, every
    /// root the four-bit field can hold, indices at and beyond the edges of the rhombus, every
    /// sub-hex, and 2,000 random values.
    fn arbitrary_ids() -> Vec<u64> {
        let mut ids = Vec::new();
        for level_i9r in 0..32 {
            let p = pow3(level_i9r);
            for root in 0..16 {
                for ix in [
                    0,
                    1,
                    p - 1,
                    p,
                    p.wrapping_mul(p).wrapping_sub(1),
                    (1 << 51) - 1,
                ] {
                    for sh in 0..4 {
                        ids.push(pack(level_i9r, root, ix & ((1 << 51) - 1), sh));
                    }
                }
            }
        }
        let mut state: u64 = 0x6E16_6E16;
        for _ in 0..2_000 {
            state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = state;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            ids.push(z ^ (z >> 31));
        }
        ids
    }

    /// The null zone, a phantom polar sub-hex (`AA-0-C`), the two polar-root classes (`BA-2-A`,
    /// `BB-6-A`), a level-34 identifier and the engine's own level-34 child of a zone at 33,
    /// 10,000 random values, and the identifiers of `arbitrary_ids`, through every hierarchy
    /// method and predicate: `Grid` answers no relation for any that names no cell, and none
    /// panics, whether through `Grid` and `Zone` or through `HexA3`'s own methods called
    /// directly, where the answers are meaningless. Run in release as well as in debug, since
    /// only debug builds check integer overflow.
    #[test]
    fn the_hierarchy_of_any_identifier_does_not_panic() {
        let grid = isea3h();
        let mut ids = vec![
            ZoneId::NULL.0,
            pack(0, 10, 0, 2),
            0x0340_0000_0000_0008,
            0x0360_0000_0000_0018,
            pack(17, 0, 0, 0),
            0x22a5_e5dd_1522_921c,
        ];
        let mut state: u64 = 0x7E57_0007;
        for _ in 0..10_000 {
            state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = state;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            ids.push(z ^ (z >> 31));
        }
        ids.extend(arbitrary_ids());
        // A zone with three parents, two of them across an interruption, and its relatives.
        let known = grid.zone(ZoneId(0x0420_0000_0000_0024));
        let mut readable = 0;
        for &raw in &ids {
            let a = Address::new(raw, &[]);
            let parents = HexA3::parents(&a).unwrap();
            let children = HexA3::children(&a).unwrap();
            assert!(parents.len() <= 3 && children.len() <= 7, "{raw:#018x}");
            let _ = HexA3::centroid_parent(&a).unwrap();
            let _ = HexA3::is_centroid_child(&a).unwrap();
            let z = grid.zone(ZoneId(raw));
            for other in [&known, &z] {
                let _ = z.is_ancestor_of(other);
                let _ = other.is_ancestor_of(&z);
                let _ = z.is_sibling_of(other);
                let _ = z.is_immediate_child_of(other);
                let _ = other.is_immediate_child_of(&z);
            }
            if is_readable(raw) {
                readable += 1;
                continue;
            }
            assert!(z.parents().is_empty(), "{raw:#018x}");
            assert!(z.parent().is_none(), "{raw:#018x}");
            assert!(z.children().is_empty(), "{raw:#018x}");
            assert!(z.centroid_parent().is_none(), "{raw:#018x}");
            assert!(!z.is_centroid_child(), "{raw:#018x}");
            assert!(!z.is_ancestor_of(&known) && !known.is_ancestor_of(&z));
            assert!(!z.is_sibling_of(&known) && !known.is_sibling_of(&z));
            assert!(!z.is_immediate_child_of(&known) && !known.is_immediate_child_of(&z));
        }
        assert!(readable > 0, "no identifier names a cell");
    }

    /// No zone at resolution 33 has children, the engine's being level-34 identifiers that it
    /// cannot read back (the departure): the twelve pentagons, the corners of every rhombus in
    /// each of the sub-hexes B, C and D, and 3,000 random zones, through `Grid` and through
    /// `HexA3::children` directly.
    #[test]
    fn zones_at_resolution_33_have_no_children() {
        let grid = isea3h();
        let p = pow3(16);
        let mut ids: Vec<u64> = (0..12).map(|root| pack(16, root, 0, 1)).collect();
        for root in 0..10 {
            for ix in [0, 1, p - 1, p, p * p - p, p * p - 1] {
                for sh in 1..4 {
                    ids.push(pack(16, root, ix, sh));
                }
            }
        }
        let mut state: u64 = 0x3333_0033;
        for _ in 0..3_000 {
            state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = state;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^= z >> 31;
            ids.push(pack(16, z % 10, (z >> 4) % (p * p), 1 + (z >> 60) % 3));
        }
        for id in ids {
            assert!(is_readable(id), "{id:#018x}");
            assert_eq!(grid.resolution(ZoneId(id)), 33);
            assert!(grid.children(ZoneId(id)).is_empty(), "{id:#018x}");
            assert_eq!(HexA3::children(&Address::new(id, &[])), Some(Vec::new()));
        }
    }

    /// The depth-0 branches of `HexA3::first_sub_zone`, `HexA3::sub_zones` and
    /// `HexA3::sub_zone_at_index` are unreachable through `Grid`, which settles depth 0 itself
    /// before it ever asks the topology (`check_sub_zone_depth`); called directly, each still
    /// answers the zone itself, as its own documentation says.
    #[test]
    fn sub_zone_methods_at_depth_0_answer_the_zone_itself_for_a_direct_caller() {
        let a = Address::new(pack(0, 6, 1, 0), &[]);
        assert_eq!(HexA3::first_sub_zone(&a, 0), Some(a));
        assert_eq!(HexA3::sub_zones(&a, 0), Some(vec![a]));
        assert_eq!(HexA3::sub_zone_at_index(&a, 0, 0), Some(a));
    }

    /// `HexA3::sub_zones` called directly, without `Grid`'s check, at a depth whose count already
    /// exceeds `MAX_MATERIALISED_SUB_ZONES`: a hexagon at depth 14 counts 4,785,157
    /// sub-zones, above the four-million cap and still well below the engine's own `2^28`. It
    /// answers `None` promptly, without materialising a single centroid.
    #[test]
    fn sub_zones_refuses_above_the_materialisation_cap() {
        let a = Address::new(pack(0, 6, 1, 0), &[]);
        let count = hex_a3_subzones::count(a.base, 14).unwrap();
        assert_eq!(count, 4_785_157);
        assert!(count > crate::grid::MAX_MATERIALISED_SUB_ZONES);
        let start = std::time::Instant::now();
        assert_eq!(HexA3::sub_zones(&a, 14), None);
        assert!(start.elapsed() < std::time::Duration::from_millis(200));
        // The topology reads the limit in force, which is the one that the grid reads: a
        // list is given exactly where it is no longer than that limit, whatever the
        // environment of the run has made of it.
        let limit = crate::grid::max_materialised_sub_zones();
        for depth in 1..=6 {
            let count = hex_a3_subzones::count(a.base, depth).unwrap();
            assert_eq!(
                HexA3::sub_zones(&a, depth).is_some(),
                count <= limit,
                "depth {depth}"
            );
        }
    }

    // --- the refined boundary in the plane: the base vertices and `refine5x6` ---
    //
    // The engine has no planar accessor that honours the refinement on aperture 3
    // (`RI3H.ec:656-663`), so the planar code is held to the engine through the ring: this
    // crate's refined boundary in radians, as `Grid::refined_ring_radians` gives it, against
    // the engine's `getZoneRefinedWGS84Vertices`.

    /// This crate's refined boundary of the zone `id` in radians, a latitude and a longitude to
    /// each point, with every edge divided into `n` parts, 0 for the engine's own choice.
    fn ring<P: Projection>(grid: &Grid<P, HexA3, I3h>, id: u64, n: i32) -> Vec<(f64, f64)> {
        let zone = ZoneId(id);
        let a = grid.geometry_address(zone).expect("a zone with geometry");
        grid.refined_ring_radians(zone, &a, n)
    }

    /// Where the ring `ours` leaves the engine's, as a sequence of doubles compared by their
    /// bits: `None` where the two are the same.
    fn ring_difference(ours: &[(f64, f64)], theirs: &[(f64, f64)]) -> Option<String> {
        if ours.len() != theirs.len() {
            return Some(format!(
                "{} points, the engine has {}",
                ours.len(),
                theirs.len()
            ));
        }
        ours.iter()
            .zip(theirs)
            .position(|(a, b)| bits2(*a) != bits2(*b))
            .map(|i| {
                format!(
                    "point {i} of {} is {:?}, the engine has {:?}",
                    ours.len(),
                    ours[i],
                    theirs[i]
                )
            })
    }

    fn bits2(p: (f64, f64)) -> (u64, u64) {
        (p.0.to_bits(), p.1.to_bits())
    }

    /// One ring of the engine: `getZoneRefinedWGS84Vertices` of the zone on the grid named, at
    /// the edge refinement given, read from DGGAL v0.0.6, BuildID
    /// `e75d6ab18f8b460713ebcd21df6f07922fb909c3`: the bits of a latitude and of a longitude,
    /// in radians, to each point.
    struct EngineRing {
        /// The arm of the base vertices, or the site of the compiled order, that the ring holds.
        what: &'static str,
        grid: &'static str,
        text: &'static str,
        id: u64,
        refinement: i32,
        ring: &'static [(u64, u64)],
    }

    /// One zone for each arm of `getBaseRefinedVertices`, with every edge in one part, so that
    /// the ring is the arm's own corners and nothing else, save that a pentagon off the poles
    /// has five of the arm's six; and one ring for each site of the compiled order that moves
    /// a ring of the engine, at a zone and a refinement at which it does.
    const ENGINE_RINGS: [EngineRing; 17] = [
        EngineRing {
            what: "the even level, the north polar pentagon",
            grid: "ISEA3H",
            text: "CA-0-A",
            id: 0x0540_0000_0000_0000,
            refinement: 1,
            ring: &[
                (0x3ff1_3e9c_b6b5_1637, 0x3fd2_6e77_0b7a_3a87),
                (0x3ff1_3e9c_b6b5_1635, 0x3fba_5baf_5f18_eceb),
                (0x3fef_c0c0_aeab_98df, 0x3fb0_d0aa_1998_1e44),
                (0x3fee_3a73_d8f5_0744, 0x3fc9_0562_e340_75b8),
                (0x3fef_c0c0_aeab_98dd, 0x3fd4_d138_5cda_6e31),
            ],
        },
        EngineRing {
            what: "the even level, the south polar pentagon",
            grid: "ISEA3H",
            text: "CB-0-A",
            id: 0x0560_0000_0000_0000,
            refinement: 1,
            ring: &[
                (0xbfef_c0c0_aeab_98de, 0xc008_9b76_0377_6c28),
                (0xbff1_3e9c_b6b5_1637, 0xc008_4f1d_d94b_65b2),
                (0xbff1_3e9c_b6b5_1637, 0xc006_d42c_72d4_e5c7),
                (0xbfef_c0c0_aeab_98dd, 0xc006_87d4_48a8_df52),
                (0xbfee_3a73_d8f5_0744, 0xc007_91a5_2610_25bd),
            ],
        },
        EngineRing {
            what: "the even level, a hexagon of the regular arm",
            grid: "ISEA3H",
            text: "C8-2C-A",
            id: 0x0500_0000_0000_00b0,
            refinement: 1,
            ring: &[
                (0x3ff7_7b67_4a6f_9f04, 0x3ff1_b3e1_5e5f_519c),
                (0x3ff6_5655_4405_5fc0, 0x3ff6_497f_8070_ac89),
                (0x3ff5_e45f_d6bc_f93a, 0x3ffc_42a7_b0ac_3bae),
                (0x3ff6_5655_4405_5fc4, 0x4001_1de7_f073_e572),
                (0x3ff7_7b67_4a6f_9f0b, 0x4003_68b7_017c_92d7),
                (0x3ff8_7b1d_1ab4_ece5, 0x3ffc_42a7_b0ac_3b28),
            ],
        },
        EngineRing {
            what: "the even level, a pentagon of the regular arm",
            grid: "ISEA3H",
            text: "C0-0-A",
            id: 0x0400_0000_0000_0000,
            refinement: 1,
            ring: &[
                (0x3ff1_3e9c_b6b5_1644, 0xc008_4f1d_d94b_65b3),
                (0x3fef_c0c0_aeab_98f8, 0xc008_9b76_0377_6c28),
                (0x3fee_3a73_d8f5_0760, 0xc007_91a5_2610_25bc),
                (0x3fef_c0c0_aeab_98fa, 0xc006_87d4_48a8_df51),
                (0x3ff1_3e9c_b6b5_1643, 0xc006_d42c_72d4_e5c7),
            ],
        },
        EngineRing {
            what: "the odd level, the north polar pentagon",
            grid: "ISEA3H",
            text: "CA-0-B",
            id: 0x0540_0000_0000_0001,
            refinement: 1,
            ring: &[
                (0x3ff0_ef07_32f4_e64e, 0x3fc9_0562_e340_75b6),
                (0x3ff0_7b93_03b5_e3f3, 0x3fbf_89ef_34ad_96be),
                (0x3fef_96b7_3b98_1826, 0x3fc3_aedc_79c5_1bc9),
                (0x3fef_96b7_3b98_1826, 0x3fce_5be9_4cbb_cfa4),
                (0x3ff0_7b93_03b5_e3f3, 0x3fd1_22e7_1615_1012),
            ],
        },
        EngineRing {
            what: "the odd level, the south polar pentagon",
            grid: "ISEA3H",
            text: "CB-0-B",
            id: 0x0560_0000_0000_0001,
            refinement: 1,
            ring: &[
                (0xbfef_96b7_3b98_1817, 0xc007_3c3c_bf78_701f),
                (0xbfef_96b7_3b98_1822, 0xc007_e70d_8ca7_db62),
                (0xbff0_7b93_03b5_e3f8, 0xc008_25ab_da9e_c062),
                (0xbff0_ef07_32f4_e650, 0xc007_91a5_2610_25b8),
                (0xbff0_7b93_03b5_e3ef, 0xc006_fd9e_7181_8b14),
            ],
        },
        EngineRing {
            what: "the odd level, a hexagon of the regular arm of sub-hex B",
            grid: "ISEA3H",
            text: "C8-2C-B",
            id: 0x0500_0000_0000_00b1,
            refinement: 1,
            ring: &[
                (0x3ff7_112c_5bcf_3468, 0x3ff6_e9e5_e83e_f64a),
                (0x3ff6_839c_6e4c_20ad, 0x3ffa_2ce7_035c_b43e),
                (0x3ff6_839c_6e4c_20ae, 0x3ffe_5868_5dfb_c324),
                (0x3ff7_112c_5bcf_346b, 0x4000_cdb4_bc8c_c083),
                (0x3ff7_c9b4_4b0e_777f, 0x4000_28d0_992d_7d03),
                (0x3ff7_c9b4_4b0e_777c, 0x3ff8_33ae_2efd_7d24),
            ],
        },
        EngineRing {
            what: "the odd level, a pentagon of the regular arm of sub-hex B",
            grid: "ISEA3H",
            text: "C0-0-B",
            id: 0x0400_0000_0000_0001,
            refinement: 1,
            ring: &[
                (0x3ff0_ef07_32f4_e65e, 0xc007_91a5_2610_25bc),
                (0x3ff0_7b93_03b5_e400, 0xc008_25ab_da9e_c05e),
                (0x3fef_96b7_3b98_1841, 0xc007_e70d_8ca7_db5c),
                (0x3fef_96b7_3b98_1841, 0xc007_3c3c_bf78_701d),
                (0x3ff0_7b93_03b5_e400, 0xc006_fd9e_7181_8b18),
            ],
        },
        EngineRing {
            what: "the odd level, sub-hex C",
            grid: "ISEA3H",
            text: "C0-3-C",
            id: 0x0400_0000_0000_000e,
            refinement: 1,
            ring: &[
                (0x3ff6_5f1d_1465_e233, 0xc005_a31f_8db4_50cd),
                (0x3ff6_6f6c_c10c_e110, 0xc003_7dd4_36cc_39be),
                (0x3ff6_ffd5_e013_a8a4, 0xc002_48d5_918f_c5e5),
                (0x3ff7_a36f_f83a_2508, 0xc003_f202_5958_2138),
                (0x3ff7_75b5_09ff_b676, 0xc007_91a5_2608_df40),
                (0x3ff6_cae3_8238_c55b, 0xc007_91a5_260a_ec81),
            ],
        },
        EngineRing {
            what: "the odd level, sub-hex D",
            grid: "ISEA3H",
            text: "C5-1B-D",
            id: 0x04a0_0000_0000_006f,
            refinement: 1,
            ring: &[
                (0xbff7_75b5_09ff_b674, 0x3fc9_0562_e3b4_dd68),
                (0xbff7_a36f_f83a_2507, 0x3fe4_bfe3_ebb0_2f7c),
                (0xbff6_ffd5_e013_a8a3, 0x3feb_6497_0ad1_9cb8),
                (0xbff6_6f6c_c10c_e10f, 0x3fe6_909c_75df_cd55),
                (0xbff6_5f1d_1465_e232, 0x3fdb_f6de_347e_e25e),
                (0xbff6_cae3_8238_c55b, 0x3fc9_0562_e394_0960),
            ],
        },
        EngineRing {
            what: "GB3-1",
            grid: "ISEA3H",
            text: "H0-CB5-A",
            id: 0x0e00_0000_0000_32d4,
            refinement: 0,
            ring: &[
                (0x3ff8_e7c1_0a87_e675, 0xc007_613e_813a_dad0),
                (0x3ff8_e78b_b3eb_4f4b, 0xc007_663d_822e_7b95),
                (0x3ff8_e756_46c5_0757, 0xc007_6b33_6bf2_6913),
                (0x3ff8_e720_c352_9e3e, 0xc007_7020_5383_6df5),
                (0x3ff8_e6eb_29d0_ddf3, 0xc007_7504_4dad_8f65),
                (0x3ff8_e6b5_7a7b_cd4e, 0xc007_79df_6f0c_4291),
                (0x3ff8_e67e_ccba_c9e0, 0xc007_7538_72d5_7965),
                (0x3ff8_e648_0aeb_d1eb, 0xc007_7099_fe68_4e6c),
                (0x3ff8_e611_3545_b174, 0xc007_6c03_fd75_ba78),
                (0x3ff8_e5da_4bfe_805c, 0xc007_6776_5be2_7387),
                (0x3ff8_e5a3_4f4b_a4e9, 0xc007_62f1_05c6_8c96),
                (0x3ff8_e5a1_7261_8660, 0xc007_599b_0661_5af5),
                (0x3ff8_e59f_434c_93dd, 0xc007_5045_ab02_f73a),
                (0x3ff8_e59c_c216_1c2c, 0xc007_46f1_0d00_1613),
                (0x3ff8_e599_eec8_bbc7, 0xc007_3d9d_45a1_5c34),
                (0x3ff8_e596_c970_5bfd, 0xc007_344a_6e21_c538),
                (0x3ff8_e5cb_3897_cda3, 0xc007_2f49_31d1_5275),
                (0x3ff8_e5ff_901f_56c3, 0xc007_2a3f_38ad_3946),
                (0x3ff8_e633_cfc8_baad, 0xc007_252c_6fc0_353d),
                (0x3ff8_e667_f754_fed7, 0xc007_2010_c3ee_b4b5),
                (0x3ff8_e69c_0684_6891, 0xc007_1aec_21f6_d828),
                (0x3ff8_e6d4_2db5_b6c6, 0xc007_1f40_ceb3_b6f0),
                (0x3ff8_e70c_435b_9dd0, 0xc007_239d_bb99_57a8),
                (0x3ff8_e744_4744_0889, 0xc007_2802_fdac_2a82),
                (0x3ff8_e77c_393c_31fa, 0xc007_2c70_aa2c_fa80),
                (0x3ff8_e7b4_1910_a2a0, 0xc007_30e6_d699_86ec),
                (0x3ff8_e7b7_59df_b8c1, 0xc007_3a90_1ccb_f214),
                (0x3ff8_e7ba_45ac_767c, 0xc007_443a_6d0d_c888),
                (0x3ff8_e7bc_dc69_a9ef, 0xc007_4de5_ab4e_7c42),
                (0x3ff8_e7bf_1e0b_936b, 0xc007_5791_bb6d_3bd5),
            ],
        },
        EngineRing {
            what: "GB3-2",
            grid: "IVEA3H",
            text: "I0-0-A",
            id: 0x1000_0000_0000_0000,
            refinement: 1,
            ring: &[
                (0x3ff0_4f0e_e947_864e, 0xc007_91de_692f_c837),
                (0x3ff0_4e9c_a08b_4e49, 0xc007_9201_c8d5_0d5a),
                (0x3ff0_4e56_0175_86f4, 0xc007_91a5_2610_25bb),
                (0x3ff0_4e9c_a08b_4e44, 0xc007_9148_834b_3e18),
                (0x3ff0_4f0e_e947_8652, 0xc007_916b_e2f0_833f),
            ],
        },
        EngineRing {
            what: "R56-1 and R56-5",
            grid: "ISEA3H",
            text: "G2-256-B",
            id: 0x0c40_0000_0000_0959,
            refinement: 2,
            ring: &[
                (0x3ff0_cb4e_caf7_c3dc, 0xbfc6_e7c1_ac76_959f),
                (0x3ff0_cb4c_8862_89ab, 0xbfc6_f840_9df1_c21d),
                (0x3ff0_cb4a_286a_0939, 0xbfc7_08bf_b707_d074),
                (0x3ff0_ca37_b8a9_1f15, 0xbfc7_10e6_9456_cfe2),
                (0x3ff0_c925_4761_3910, 0xbfc7_190b_b78c_1edb),
                (0x3ff0_c815_590d_b238, 0xbfc7_1094_8bc9_4b29),
                (0x3ff0_c705_5d04_639a, 0xbfc7_081f_7ca9_2fae),
                (0x3ff0_c707_af89_49df, 0xbfc6_f7a8_03b9_8ac8),
                (0x3ff0_c709_e4b6_2b76, 0xbfc6_e730_b286_8f86),
                (0x3ff0_c81c_2f14_c032, 0xbfc6_df28_291a_f48f),
                (0x3ff0_c92e_72b5_3806, 0xbfc6_d71d_989d_672f),
                (0x3ff0_ca3e_a2d0_8e98, 0xbfc6_df6e_c147_4696),
            ],
        },
        EngineRing {
            what: "R56-4 and R56-8",
            grid: "ISEA3H",
            text: "G3-63D89-B",
            id: 0x0c60_0000_0018_f625,
            refinement: 3,
            ring: &[
                (0xbff0_bbd9_b328_2e36, 0xc003_aae5_1504_2ad8),
                (0xbff0_bb2f_cd3e_f6cf, 0xc003_ab51_aae3_acd4),
                (0xbff0_ba85_e076_d5c1, 0xc003_abbe_3134_0b2e),
                (0xbff0_b9db_eccf_77ee, 0xc003_ac2a_a7f6_fdf5),
                (0xbff0_b9ee_928b_84ea, 0xc003_acd9_f45e_3445),
                (0xbff0_ba01_2b21_a85e, 0xc003_ad89_4263_9ab4),
                (0xbff0_ba13_b691_c7db, 0xc003_ae38_9205_24a0),
                (0xbff0_bad0_4dd0_b2a9, 0xc003_ae7b_9f81_1b06),
                (0xbff0_bb8c_e3af_ecda, 0xc003_aebe_b894_bf42),
                (0xbff0_bc49_782f_15b1, 0xc003_af01_dd43_29e7),
                (0xbff0_bcf3_f0fd_4fbe, 0xc003_ae99_1dad_c719),
                (0xbff0_bd9e_64b4_5f10, 0xc003_ae30_501f_b1ca),
                (0xbff0_be48_d353_79e7, 0xc003_adc7_7495_daa0),
                (0xbff0_be36_3899_b9a7, 0xc003_ad17_d70f_d7ba),
                (0xbff0_be23_90b5_48c5, 0xc003_ac68_3b2a_4b63),
                (0xbff0_be10_dba6_4224, 0xc003_abb8_a0e7_44ab),
                (0xbff0_bd53_ce8c_9ec3, 0xc003_ab72_1227_09cc),
                (0xbff0_bc96_c10d_7dbc, 0xc003_ab2b_8e31_1dd7),
            ],
        },
        EngineRing {
            what: "R56-15",
            grid: "RTEA3H",
            text: "J2-4ADF-A",
            id: 0x1240_0000_0001_2b7c,
            refinement: 2,
            ring: &[
                (0x3ff0_71a3_69ef_54f9, 0x3fc2_0c1c_e738_217d),
                (0x3ff0_71b0_7477_44f8, 0x3fc2_0b7c_d8c5_520c),
                (0x3ff0_71bd_7ef7_36aa, 0x3fc2_0adc_c872_d702),
                (0x3ff0_71b8_6913_2c0d, 0x3fc2_09bc_428b_07be),
                (0x3ff0_71b3_5309_4e2a, 0x3fc2_089b_bd99_cc11),
                (0x3ff0_71a1_329a_686b, 0x3fc2_081b_4c3b_da6c),
                (0x3ff0_718d_5e84_ca9b, 0x3fc2_07fa_d629_6e41),
                (0x3ff0_717f_7a6e_202c, 0x3fc2_08ca_e23d_d14c),
                (0x3ff0_7171_9648_1ee0, 0x3fc2_099a_ebdb_1807),
                (0x3ff0_7177_860f_287a, 0x3fc2_0a8b_6f2d_7303),
                (0x3ff0_717d_75bc_107b, 0x3fc2_0b7b_f364_3434),
                (0x3ff0_7191_49bc_1880, 0x3fc2_0b9c_702d_d455),
            ],
        },
        EngineRing {
            what: "R56-17",
            grid: "IVEA3H",
            text: "G3-673A4-A",
            id: 0x0c60_0000_0019_ce90,
            refinement: 2,
            ring: &[
                (0xbff0_c6d2_868f_486d, 0xc004_1352_97a0_a0fc),
                (0xbff0_c4b1_10f3_f8c7, 0xc004_1374_733e_22b1),
                (0xbff0_c28f_9dea_903d, 0xc004_1396_796c_8179),
                (0xbff0_c197_388c_53e1, 0xc004_1531_0f56_af1d),
                (0xbff0_c09e_89ac_cd96, 0xc004_16cb_5bef_b2a1),
                (0xbff0_c1c7_19ec_e3c6, 0xc004_1844_26c0_5a53),
                (0xbff0_c2ef_71a5_ef42, 0xc004_19bd_5c10_4a0e),
                (0xbff0_c511_0699_575a, 0xc004_199c_bebe_3f53),
                (0xbff0_c732_9e0a_67f8, 0xc004_197c_4d34_1c9f),
                (0xbff0_c82b_c34d_e41b, 0xc004_17e1_8152_ebc2),
                (0xbff0_c924_9eef_4ef9, 0xc004_1646_6b52_322f),
                (0xbff0_c7fb_af1d_b0dd, 0xc004_14cc_4bb5_0654),
            ],
        },
        EngineRing {
            what: "AIP-3, the points of a plain edge, at a pentagon at three parts to an edge",
            grid: "IVEA3H",
            text: "A1-0-B",
            id: 0x0020_0000_0000_0001,
            refinement: 3,
            ring: &[
                (0x3fd6_580c_33af_b5a0, 0xc004_15a5_3e7f_3b0f),
                (0x3fcf_4727_39e5_831a, 0xc004_e292_a35d_2443),
                (0x3fc0_25db_d921_d615, 0xc005_8c1a_ad13_3672),
                (0x3cdb_1666_a542_2886, 0xc006_1250_a972_e75a),
                (0xbfc0_25db_d921_d594, 0xc005_8c1a_ad13_3672),
                (0xbfcf_4727_39e5_82b6, 0xc004_e292_a35d_2443),
                (0xbfd6_580c_33af_b56b, 0xc004_15a5_3e7f_3b10),
                (0xbfd4_d229_7274_ecd9, 0xc002_e6d1_9fa2_6ffc),
                (0xbfd1_d984_48fa_1a14, 0xc001_c473_693e_59b5),
                (0xbfcb_4691_d7dc_6051, 0xc000_bb74_c86b_ae1f),
                (0xbfb2_6523_4fb6_6b8e, 0xc000_9977_f0d8_0da8),
                (0x3fb2_6523_4fb6_6c20, 0xc000_9977_f0d8_0da8),
                (0x3fcb_4691_d7dc_6097, 0xc000_bb74_c86b_ae1d),
                (0x3fd1_d984_48fa_1a3c, 0xc001_c473_693e_59b3),
                (0x3fd4_d229_7274_ed06, 0xc002_e6d1_9fa2_6ffa),
            ],
        },
    ];

    /// The base vertices and `refine5x6` against rings read from the engine, without DGGAL.
    ///
    /// The first ten rings hold the eight arms of `getBaseRefinedVertices`, the two regular
    /// arms each at a hexagon and at a pentagon. The next six pin the eight sites of the
    /// planar code whose order moves a ring: put back alone in the eC's order, GB3-1 moves the
    /// ring named for it here, and so do GB3-2, R56-15 and R56-17 theirs; R56-1 and R56-5 each
    /// move the one ring named for both, and R56-4 and R56-8 likewise. Each was so put back,
    /// and this test seen to fail at its ring. The last, a pentagon of IVEA3H at three parts to
    /// an edge, gives site AIP-3 of `addIntermediatePoints` a literal ring of the engine on
    /// aperture 3: with that site put back in the eC's order alone, it fails here.
    #[test]
    fn the_base_vertices_and_refine5x6_follow_the_compiled_order() {
        let mut differences: Vec<String> = Vec::new();
        for engine in &ENGINE_RINGS {
            assert_eq!(I3h::to_text(ZoneId(engine.id)), engine.text);
            let ours = match engine.grid {
                "ISEA3H" => ring(crate::registry::isea3h(), engine.id, engine.refinement),
                "IVEA3H" => ring(crate::registry::ivea3h(), engine.id, engine.refinement),
                "RTEA3H" => ring(crate::registry::rtea3h(), engine.id, engine.refinement),
                other => panic!("no grid {other}"),
            };
            let theirs: Vec<(f64, f64)> = engine
                .ring
                .iter()
                .map(|&(lat, lon)| (f64::from_bits(lat), f64::from_bits(lon)))
                .collect();
            if let Some(difference) = ring_difference(&ours, &theirs) {
                differences.push(format!(
                    "{}: {} {} at a refinement of {}: {difference}",
                    engine.what, engine.grid, engine.text, engine.refinement
                ));
            }
        }
        assert!(
            differences.is_empty(),
            "{} rings leave the engine's:\n{}",
            differences.len(),
            differences.join("\n")
        );
    }

    /// The refinements at which the ring is compared with the engine's: the automatic one, and
    /// those at which each site was found to move a ring.
    #[cfg(feature = "oracle")]
    const REFINEMENTS: [i32; 5] = [0, 1, 2, 3, 7];

    /// The pinning zones on aperture 3, by the site of the planar code each
    /// pins: reverted alone to the eC's order, the site moves the ring of each of its zones,
    /// on one of the three grids at least, at some refinement of `REFINEMENTS`. A zone that
    /// pins a site on one grid need not pin it on another, since the projection decides the
    /// last bits; every zone is compared on the three. The last entry holds one zone for each
    /// arm of the base vertices.
    #[cfg(feature = "oracle")]
    const PINNED_ZONES: [(&str, &[u64]); 7] = [
        // `d/3` as `d * fl(1/3)`: `H0-CB5-A`, `H0-CB5-C`, `G0-4A7-B`, `G0-446-B`, `I0-2656-A`,
        // `G0-4A8-B`, `G1-1BF16-B`, `G1-76000-B` and `G0-446-A`.
        (
            "GB3-1",
            &[
                0x0e00_0000_0000_32d4,
                0x0e00_0000_0000_32d6,
                0x0c00_0000_0000_129d,
                0x0c00_0000_0000_1119,
                0x1000_0000_0000_9958,
                0x0c00_0000_0000_12a1,
                0x0c20_0000_0006_fc59,
                0x0c20_0000_001d_8001,
                0x0c00_0000_0000_1118,
            ],
        ),
        // `2*d/3` as `d * fl(2/3)`: `G0-7C95E-A`, `H0-CD1-A`, `I0-0-A`, `K0-0-A`, `G1-1C7A0-C`,
        // `G0-0-A`, `G0-163-D` and `H0-0-A`.
        (
            "GB3-2",
            &[
                0x0c00_0000_001f_2578,
                0x0e00_0000_0000_3344,
                0x1000_0000_0000_0000,
                0x1400_0000_0000_0000,
                0x0c20_0000_0007_1e82,
                0x0c00_0000_0000_0000,
                0x0c00_0000_0000_058f,
                0x0e00_0000_0000_0000,
            ],
        ),
        // The early crossing at a top dent, of `p` to the right and of `next`: `G2-256-B`,
        // `G2-2A0-B`, `G2-2B3-B`, `G2-258-B`, `G2-26A-B`, `G2-2B4-B`, `G2-224-B`, `G2-26D-B`
        // and `G2-27F-B`. Each zone pins both sites.
        (
            "R56-1 and R56-5",
            &[
                0x0c40_0000_0000_0959,
                0x0c40_0000_0000_0a81,
                0x0c40_0000_0000_0acd,
                0x0c40_0000_0000_0961,
                0x0c40_0000_0000_09a9,
                0x0c40_0000_0000_0ad1,
                0x0c40_0000_0000_0891,
                0x0c40_0000_0000_09b5,
                0x0c40_0000_0000_09fd,
            ],
        ),
        // The early crossing at a bottom dent, of `p` to the right and of `next`: `G3-63D89-B`,
        // `G3-779A0-B`, `G3-7E5D6-B`, `G3-6433B-B`, `G3-7B294-B`, `G3-61884-B`, `G3-6B523-B`
        // and `G3-6E865-B`. Each zone pins both sites.
        (
            "R56-4 and R56-8",
            &[
                0x0c60_0000_0018_f625,
                0x0c60_0000_001d_e681,
                0x0c60_0000_001f_9759,
                0x0c60_0000_0019_0ced,
                0x0c60_0000_001e_ca51,
                0x0c60_0000_0018_6211,
                0x0c60_0000_001a_d48d,
                0x0c60_0000_001b_a195,
            ],
        ),
        // The second point of an interruption at a top dent crossed to the right: `J2-4ADF-A`,
        // `P2-D53574-A`, `J2-40FF-A`, `J2-3B63-A`, `J2-3D53-A`, `J2-4326-A`, `J2-39D3-A`,
        // `J2-3DAE-A` and `J2-3F9A-A`.
        (
            "R56-15",
            &[
                0x1240_0000_0001_2b7c,
                0x1e40_0000_0354_d5d0,
                0x1240_0000_0001_03fc,
                0x1240_0000_0000_ed8c,
                0x1240_0000_0000_f54c,
                0x1240_0000_0001_0c98,
                0x1240_0000_0000_e74c,
                0x1240_0000_0000_f6b8,
                0x1240_0000_0000_fe68,
            ],
        ),
        // The second point of an interruption at a bottom dent crossed to the right:
        // `G3-673A4-A`, `G3-6A6E6-A`, `G3-6DA28-A`, `F3-72DB-A`, `G3-6AC98-A`, `G3-6DFDA-A`,
        // `G3-64E9F-A` and `G3-681E1-A`.
        (
            "R56-17",
            &[
                0x0c60_0000_0019_ce90,
                0x0c60_0000_001a_9b98,
                0x0c60_0000_001b_68a0,
                0x0a60_0000_0001_cb6c,
                0x0c60_0000_001a_b260,
                0x0c60_0000_001b_7f68,
                0x0c60_0000_0019_3a7c,
                0x0c60_0000_001a_0784,
            ],
        ),
        // The even level: `CA-0-A` and `CB-0-A`, the polar pentagons; `C8-2C-A`, a hexagon, and
        // `C0-0-A`, a pentagon, both of the regular arm. The odd level: `CA-0-B` and `CB-0-B`;
        // `C8-2C-B` and `C0-0-B`, of the regular arm of sub-hex B; `C0-3-C` and `C5-1B-D`.
        (
            "an arm of the base vertices",
            &[
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
            ],
        ),
    ];

    /// Each compiled-order site of the planar code, at the zones that pin it, and each arm of
    /// the base vertices, against the engine's ring on the three grids. Every difference is
    /// gathered before the test fails, so that a site put back in the eC's order names its
    /// zones in the failure.
    #[test]
    #[cfg(feature = "oracle")]
    fn the_compiled_order_sites_hold_at_the_zones_that_pin_them() {
        fn on<P: Projection>(
            oracle: &str,
            grid: &Grid<P, HexA3, I3h>,
            differences: &mut Vec<String>,
        ) {
            for (site, zones) in PINNED_ZONES {
                for &id in zones {
                    let text = dggal_oracle::text_id(oracle, id);
                    assert_eq!(grid.text_id(ZoneId(id)), text);
                    for n in REFINEMENTS {
                        let theirs = dggal_oracle::refined_vertices(oracle, id, n);
                        if let Some(difference) = ring_difference(&ring(grid, id, n), &theirs) {
                            differences.push(format!(
                                "{site}: {oracle} {text} ({id:#018x}) at a refinement of {n}: \
                                 {difference}"
                            ));
                        }
                    }
                }
            }
        }

        let mut differences: Vec<String> = Vec::new();
        on(
            dggal_oracle::ISEA3H,
            crate::registry::isea3h(),
            &mut differences,
        );
        on(
            dggal_oracle::IVEA3H,
            crate::registry::ivea3h(),
            &mut differences,
        );
        on(
            dggal_oracle::RTEA3H,
            crate::registry::rtea3h(),
            &mut differences,
        );
        assert!(
            differences.is_empty(),
            "{} rings leave the engine's:\n{}",
            differences.len(),
            differences.join("\n")
        );
    }

    fn ids(texts: &str) -> Vec<u64> {
        texts
            .split_whitespace()
            .map(|text| I3h::from_text(text).unwrap().0)
            .collect()
    }

    fn texts(ids: &[u64]) -> String {
        let texts: Vec<String> = ids.iter().map(|&id| I3h::to_text(ZoneId(id))).collect();
        texts.join(" ")
    }

    /// The zones two levels up that hold the whole of a zone, on ISEA3H as the engine
    /// answers: one for each of the seven sub-zones about the centre of `B4-2-A`, that zone
    /// itself, and three for each of the six on its vertices, that zone among them. A zone of
    /// level 1 on a vertex has none, and a zone of level 0 has the parents of the null zone,
    /// whose first is the null zone with its sub-hex cleared and is no zone.
    #[test]
    fn the_grand_parents_that_hold_a_zone_are_one_or_three() {
        let parent = ids("B4-2-A")[0];
        let (mut within, mut on_a_vertex) = (0, 0);
        for sub in hex_a3_subzones::sub_zones(parent, 2).unwrap() {
            let holders = containing_grand_parents(sub);
            match holders.len() {
                1 => {
                    assert_eq!(holders, [parent]);
                    within += 1;
                }
                3 => {
                    assert!(holders.contains(&parent), "{}", texts(&holders));
                    assert!(holders.iter().all(|&g| is_readable(g) && level(g) == 2));
                    on_a_vertex += 1;
                }
                n => panic!("{n} zones hold {}", texts(&[sub])),
            }
        }
        assert_eq!((within, on_a_vertex), (7, 6));
        let of_level_2 = containing_grand_parents(parent);
        assert!(of_level_2.len() == 1 && is_readable(of_level_2[0]) && level(of_level_2[0]) == 0);
        assert!(containing_grand_parents(ids("A2-0-C")[0]).is_empty());
        let of_level_0 = containing_grand_parents(ids("A0-0-A")[0]);
        assert_eq!(of_level_0[0], 0xffff_ffff_ffff_fffc);
        assert!(of_level_0.iter().all(|&g| !is_readable(g)));
    }

    /// The compaction of sets of more than one level, which `Grid` refuses and the engine
    /// answers so: a zone with no zone two levels up to give way to is lost as soon as the
    /// set holds a finer zone, and a zone of level 0 beside a finer one is kept, the
    /// identifier that its grandparents come to, which is no zone, never entering the answer.
    /// The null zone is left out, and a zone given twice counts once.
    ///
    /// The last set shows that the loop goes on after a pass that takes nothing: the six
    /// zones of level 4 about the centre of `B4-2-A`, with the zone of level 6 at the centre
    /// of each of the six of level 4 on its vertices. The first pass, from level 6, takes
    /// nothing; the second takes `B4-2-A`, on the evidence of those six zones of level 6, kept
    /// by the first. Were the loop left at the first pass, the answer would be the twelve.
    #[test]
    fn the_compaction_of_several_levels_is_the_engines() {
        for (set, answer) in [
            ("A2-0-C C4-22-B", "C4-22-B"),
            ("A2-0-C B4-2-A", "B4-2-A"),
            ("B5-0-A A9-0-C A0-0-B A3-0-A", "A0-0-B A3-0-A B5-0-A"),
            ("A0-0-A B4-2-A", "A0-0-A B4-2-A"),
            ("A5-0-A A0-0-B", "A0-0-B A5-0-A"),
            (
                "C2-1A-A C2-23-A C4-5-A C4-7-A C4-F-A C4-10-A D2-69-A D2-B7-A D2-15C-A D4-60-A \
                 D4-69-A D4-B7-A",
                "B4-2-A C2-1A-A C2-23-A C4-5-A C4-7-A C4-F-A C4-10-A D2-69-A D2-B7-A D2-15C-A \
                 D4-60-A D4-69-A D4-B7-A",
            ),
        ] {
            assert_eq!(texts(&compact(&ids(set))), answer, "{set}");
        }
        // Every zone of level 2 with every zone of level 0: the twelve of level 0, as the
        // engine answers, and not a thirteenth that is no zone, the null zone with its
        // sub-hex cleared, which the grandparents of a zone of level 0 come to.
        let twelve: Vec<u64> = (0..12).map(|root| pack(0, root, 0, 0)).collect();
        let mut set: Vec<u64> = crate::isea3h().zones(2).unwrap().map(|z| z.0).collect();
        assert_eq!(set.len(), 92);
        set.extend(&twelve);
        assert_eq!(compact(&set), twelve);

        let mut set = ids("A0-0-A B4-2-A A0-0-A");
        set.insert(1, NULL_ZONE);
        assert_eq!(texts(&compact(&set)), "A0-0-A B4-2-A");
        assert!(compact(&[NULL_ZONE]).is_empty());
        assert!(compact(&[]).is_empty());
    }
}
