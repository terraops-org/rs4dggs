//! The sub-zone order of the aperture-3 topology: how many sub-zones a zone has at a relative
//! depth, which comes first, and the whole sequence in the engine's order, which the index of a
//! sub-zone refers to.
//!
//! The engine generates the sub-zones of a zone as centroids, row by row across the zone
//! ("scanlines"), starting on a vertex of the zone for an even relative depth and on an edge for
//! an odd one, and quantises each at the sub-zone level (`I3HSubZones.ec:1775-1847`,
//! `dggrs.ec:195-222`). Which of four generators runs depends on the parities of the zone's level
//! and of the depth; within each, the northern or southern rhombus, an edge hexagon, a pentagon
//! and a polar pentagon each take their own path, 28 cases in all (`I3HSubZones.ec:11-98`).
//!
//! Provenance: `dggrs/I3HSubZones.ec` of DGGAL v0.0.6: `getI3HFirstSubZoneCentroid`
//! (`:102-247`), `generateOddParentOddDepth` (`:250-432`), `generateEvenParentOddDepth`
//! (`:435-821`), `generateEvenParentEvenDepth` (`:824-1363`), `generateOddParentEvenDepth`
//! (`:1366-1767`), `getI3HSubZoneCentroids` (`:1775-1800`) and `iterateI3HSubZones`
//! (`:1802-1847`); and of `dggrs/RI3H.ec`, `getSubZoneAtIndex` (`:212-227`), `getSubZonesCount`
//! (`:2198-2202`) and `getFirstSubZone` (`:2204-2210`). None of them carries
//! `__attribute__ ((optimize("-fno-unsafe-math-optimizations")))`, and DGGAL is built with
//! `-O2 -ffast-math` (`Makefile.dggal:133`). In `generateEvenParentOddDepth` and
//! `generateOddParentEvenDepth` gcc re-associated every division by three and regrouped several
//! sums; in `generateEvenParentEvenDepth`, which divides by nothing, it regrouped four sums; in
//! the others it kept the eC's order. The port performs the compiled order, and the comments
//! marked `faithful:` record both. The re-associated sites these functions reach through their
//! calls lie in `getVertices`, `move5x6Vertex3` and `cross5x6Interruption`, whose ports perform
//! the compiled order.
//!
//! One departure from the engine's order, at some zones on the edges of the rhombi, where a unit
//! in the last place carries a point across an interruption of the layout that it should only
//! reach: a vertex of the zone, from which the first centroid is read (`sub_zone_vertices`, in
//! `hex_a3.rs`), or a step of a generator that is then carried across a second time
//! ([`onto_interruption`]). The engine's order there names zones twice, names zones that are not
//! descendants of the zone and omits descendants. This order is the zone's descendants, each
//! once, as many as the count, in the generator's own scanline order, and the first sub-zone and
//! the sub-zone at each index follow it. Every other order is the engine's, down to the bits of
//! the generator's centroids.
use super::hex_a3::{HexA3, MAX_LEVEL, pow3, sub_zone_vertices};
use crate::fivebysix::{cross5x6_interruption, cvtt_i32, move5x6_vertex, move5x6_vertex3};
use crate::indexings::fields;
use crate::{PlanarPoint, Topology, ZoneId};

/// The eC's `nullZone`.
const NULL_ZONE: u64 = ZoneId::NULL.0;

/// `fl(1/3)`, which gcc substituted for the divisions by three of `generateEvenParentOddDepth`
/// and `generateOddParentEvenDepth`; the copy both read, at `.rodata` `0x49ef8`.
const T13: f64 = f64::from_bits(0x3fd5_5555_5555_5555);
/// `fl(1/6)`, which gcc substituted for `/ 2.0 * u / 3` there (`I3HSubZones.ec:688`, `:698`);
/// `.rodata` at `0x49f18`.
const T16: f64 = f64::from_bits(0x3fc5_5555_5555_5555);

/// A zone's level, `levelI9R * 2 + (subHex > 0)` (`RI3H.ec:878-881`).
fn level(id: u64) -> u32 {
    let (level_i9r, _, _, sub_hex) = fields(id);
    2 * level_i9r as u32 + u32::from(sub_hex > 0)
}

/// A zone's number of sides, five for a pentagon, the polar ones included (`RI3H.ec:883-892`).
fn n_points(id: u64) -> i32 {
    let (_, _, rix, sub_hex) = fields(id);
    if rix == 0 && sub_hex <= 1 { 5 } else { 6 }
}

/// The level of the sub-zones `depth` levels below zone `id`, or `None` beyond 33, where the
/// engine's `getSubZones` answers nothing (`dggrs.ec:200`) and `Grid` refuses the request first.
fn sub_zone_level(id: u64, depth: u8) -> Option<u32> {
    let sz_level = level(id) + u32::from(depth);
    (sz_level <= u32::from(MAX_LEVEL)).then_some(sz_level)
}

/// The zone at `level` whose cell contains the planar point `(x, y)`, or the null zone, as
/// `getZoneFromCRSCentroid` answers in the 5x6 CRS (`RI3H.ec:94-117`), through
/// [`HexA3::quantize`] and its departure at two seams of the layout.
fn quantize(level: u32, (x, y): (f64, f64)) -> u64 {
    HexA3::quantize(PlanarPoint { face: -1, x, y }, level as u8).map_or(NULL_ZONE, |a| a.base)
}

/// The number of sub-zones `depth` levels below zone `id`, or `None` where they would lie beyond
/// level 33. Ports `I3HZone::getSubZonesCount` (`RI3H.ec:2198-2202`), which the engine's
/// `countSubZones` answers: `((3^d + 3^((d + 1) / 2) + 1) * nPoints + 5) / 6`, and 1 at depth 0.
#[allow(
    clippy::manual_div_ceil,
    reason = "the eC's two expressions are kept as it writes them"
)]
pub(super) fn count(id: u64, depth: u8) -> Option<u64> {
    sub_zone_level(id, depth)?;
    let r_depth = u64::from(depth);
    let n_hex_sub_zones = if r_depth > 0 {
        pow3(r_depth) + pow3((r_depth + 1) / 2) + 1
    } else {
        1
    };
    Some((n_hex_sub_zones * n_points(id) as u64 + 5) / 6)
}

/// The planar point from which the sub-zones of zone `id` at relative depth `r_depth` are
/// generated: a vertex of the zone, as `getVertices` gives it, folded back into the layout. Ports
/// `getI3HFirstSubZoneCentroid` (`I3HSubZones.ec:102-247`), compiled at `0x3880` to `0x3cc8`,
/// whose own arithmetic is comparisons in the eC's form. The eC's local names are kept, and each
/// search below keeps its first candidate, as the eC's `-1` start does.
fn first_centroid(id: u64, r_depth: u32) -> (f64, f64) {
    let vertices = sub_zone_vertices(id);
    let nv = vertices.len();
    let (level_i9r, root, rhombus_ix, sub_hex) = fields(id);
    let odd_depth = r_depth & 1 != 0;
    let odd_parent = sub_hex > 0;
    let south_rhombus = root & 1 != 0;
    let (mut x, mut y) = if odd_parent {
        // Odd Level parents (the eC's note)
        if odd_depth {
            // Left-to-Right Scanlines; Top to Bottom Scanline order; start from hexagon /
            // pentagon edge (the eC's note)
            let mut tl = 0;
            for i in 1..nv {
                let (v, t) = (vertices[i], vertices[tl]);
                if v.1 < t.1 || ((v.1 - t.1).abs() < 1e-11 && v.0 > t.0) {
                    tl = i;
                }
            }
            if root == 10 && sub_hex == 1 {
                tl = 0;
            } else if nv == 5 {
                tl = if south_rhombus {
                    if root == 11 && sub_hex == 1 { 1 } else { 4 }
                } else {
                    3
                };
            }
            // This is the vertex immediately to the right of the odd parent (the eC's note)
            vertices[tl]
        } else {
            // Bottom-to-Top Scanlines; Left to Right Scanline order; start from hexagon /
            // pentagon vertex, the leftmost in planar ISEA (the eC's note)
            let mut tl = 0;
            for i in 1..nv {
                let (v, t) = (vertices[i], vertices[tl]);
                if v.1 < t.1 || ((v.1 - t.1).abs() < 1e-11 && v.0 < t.0) {
                    tl = i;
                }
            }
            if root == 10 && sub_hex == 1 {
                tl = 4;
            } else if nv == 5 {
                tl = if root == 11 && sub_hex == 1 { 0 } else { 3 };
            }
            vertices[tl]
        }
    } else if odd_depth {
        // Even Level parents; Bottom-to-Top Scanlines; Left to Right Scanline order; start from
        // hexagon / pentagon edge (the eC's note)
        let mut left = 0;
        for i in 1..nv {
            if vertices[i].0 < vertices[left].0 {
                left = i;
            }
        }
        if root == 10 && sub_hex == 0 {
            left = 0;
        } else if nv == 5 && south_rhombus {
            left = if root == 11 && sub_hex == 0 { 4 } else { 2 };
        }
        vertices[left]
    } else if nv == 5 && (root < 10 || sub_hex != 0) && !south_rhombus {
        // Left-to-Right Scanlines; Bottom to Top Scanline order; start from hexagon / pentagon
        // vertex (the eC's note)
        vertices[4]
    } else {
        // The first vertex is between the topmost and rightmost vertices in rotated ISEA 5x6
        // space (the eC's note)
        let mut ix = if root == 10 && sub_hex == 0 {
            0
        } else if nv == 6 {
            5
        } else {
            4
        };
        let divs = pow3(level_i9r);
        let edge_hex = nv == 6
            && rhombus_ix != 0
            && if south_rhombus {
                rhombus_ix % divs == 0
            } else {
                rhombus_ix / divs == 0
            };
        if edge_hex && south_rhombus {
            ix = 4;
        }
        vertices[ix]
    };
    // getVertices() currently return negative vertices... (the eC's note)
    if x > 5.0 && y > 5.0 {
        x -= 5.0;
        y -= 5.0;
    } else if x < 0.0 && y < 1.0 {
        x += 5.0;
        y += 5.0;
    }
    (x, y)
}

/// The first sub-zone `depth` levels below zone `id`, or `None` where it would lie beyond level
/// 33. Ports `I3HZone::getFirstSubZone` (`RI3H.ec:2204-2210`): the first centroid, quantised at
/// the sub-zone level, which is the first of [`sub_zones`].
///
/// One deliberate departure, at depth 0, where this answers the zone itself, as `Grid` does
/// before it asks. `getI3HFirstSubZoneCentroid` has no case for depth 0: the engine quantises a
/// vertex of the zone at the zone's own level, and so names a neighbour, for two zones in three
/// (`B4-5-D` for `B6-1-B`), while its `getSubZones` answers the zone itself
/// (`I3HSubZones.ec:1789-1790`).
pub(super) fn first_sub_zone(id: u64, depth: u8) -> Option<u64> {
    let sz_level = sub_zone_level(id, depth)?;
    if depth == 0 {
        return Some(id);
    }
    Some(quantize(sz_level, first_centroid(id, u32::from(depth))))
}

/// `move5x6Vertex3` (`ri5x6.ec:1438-1463`) without its crossing: the step to a point on an
/// interruption of the layout, which the caller then carries across itself with
/// `cross5x6Interruption`.
///
/// A departure from the engine, which takes that step with `move5x6Vertex3`. Where the point
/// should lie exactly on the interruption, the sum can land a unit in the last place beyond it;
/// `move5x6Vertex3` then carries the point across, and the crossing that follows carries it a
/// second time, so that the engine's order names a zone that is not a descendant in place of one
/// that is, and at greater depths names zones twice. Here the point is carried across once.
/// Wherever `move5x6Vertex3` does not cross, this is its own double, and the order the engine's.
fn onto_interruption(c: (f64, f64), dx: f64, dy: f64) -> (f64, f64) {
    let (mut vx, mut vy) = (dx + c.0, c.1 + dy);
    if vx > 5.0 && vy > 5.0 {
        vx -= 5.0;
        vy -= 5.0;
    } else if vx < 0.0 && vy < 1.0 {
        vx += 5.0;
        vy += 5.0;
    }
    (vx, vy)
}

/// The callback through which a generator delivers each sub-zone centroid with its index, and
/// which answers whether to go on, as the eC's `centroidCallback`.
type Callback<'a> = &'a mut dyn FnMut(i64, (f64, f64)) -> bool;

/// Generates the sub-zone centroids of zone `id` at relative depth `r_depth`, in the engine's
/// order but for the departure of the module note, and hands each to `callback`; with an
/// `index` other than `-1`, it skips to that index and delivers it first. Answers `-1` if the
/// callback never stopped it, and otherwise the index at which it did. Ports `iterateI3HSubZones`
/// (`I3HSubZones.ec:1802-1847`), compiled at `0x9670` to `0xa7fd` with `generateOddParentOddDepth`
/// inlined.
fn iterate(id: u64, r_depth: u32, callback: Callback<'_>, index: i64) -> i64 {
    let (level_i9r, root, rhombus_ix, sub_hex) = fields(id);
    let nv = n_points(id);
    let odd_depth = r_depth & 1 != 0;
    let odd_parent = sub_hex > 0;
    let sz_level = level(id) + r_depth;
    // faithful: `u = 1.0 / POW3((szLevel / 2))` (I3HSubZones.ec:1812) is a division in DGGAL's
    // `-O2 -ffast-math` build, `libdggal.so`, BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`,
    // at `0x971c`, as the eC writes it.
    let u = 1.0 / pow3(u64::from(sz_level / 2)) as f64;
    let south_rhombus = root & 1 != 0;
    let polar_pentagon = root > 9;
    let divs = pow3(level_i9r);
    // Edge Hexagons are either -A or -D (the eC's note, in an older text form; in this one, the
    // sub-hexes A and B)
    let edge_hex = nv == 6
        && rhombus_ix != 0
        && (sub_hex == 0 || sub_hex == 1)
        && if south_rhombus {
            rhombus_ix % divs == 0
        } else {
            rhombus_ix / divs == 0
        };
    let first = first_centroid(id, r_depth);
    let r_depth = r_depth as i32;
    let generate = match (odd_parent, odd_depth) {
        (true, true) => generate_odd_parent_odd_depth,
        (true, false) => generate_odd_parent_even_depth,
        (false, true) => generate_even_parent_odd_depth,
        (false, false) => generate_even_parent_even_depth,
    };
    generate(
        callback,
        first,
        r_depth,
        u,
        nv,
        polar_pentagon,
        south_rhombus,
        edge_hex,
        index,
    )
}

/// The generator of an odd-level zone at an odd relative depth: the cases 1 (a hexagon), 5 and 6
/// (a northern and a southern edge hexagon), 13 and 14 (a northern and a southern pentagon) and
/// 21 and 22 (the two polar pentagons) of `I3HSubZones.ec:11-98`, which run left-to-right
/// scanlines, starting from an edge of the zone, in two halves of `nHalfRows + 1` and
/// `nHalfRows` rows. Ports `generateOddParentOddDepth` (`I3HSubZones.ec:250-432`), with its
/// index skip: whole scanlines at once, then a jump within the one that holds the index, and in
/// the second half, over the columns that a southern pentagon leaves out at the interruption.
/// The eC's local names are kept.
///
/// faithful: gcc inlined this function into `iterateI3HSubZones` and kept the eC's order: in
/// DGGAL's `-O2 -ffast-math` build, `libdggal.so`, BuildID
/// `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, every offset `k * u` whose `k` is an `int`
/// expression, the sign factors folded into it (I3HSubZones.ec:291-404), is one product
/// `(double)k * u` at the 24 call sites `0x9bb5` to `0xa770`, and the start of a polar scanline
/// of the second half, `1 + r * u` or `-1 - r * u` (I3HSubZones.ec:372), is formed as written
/// (`0xa6d4` to `0xa6dc`, `0xa798`, passed at `0xa6f5`). So is every offset here, `ku(k)`.
#[allow(
    clippy::too_many_arguments,
    reason = "the eC's signature, which all four generators share, is kept"
)]
fn generate_odd_parent_odd_depth(
    centroid_callback: Callback<'_>,
    first_centroid: (f64, f64),
    r_depth: i32,
    u: f64,
    nv: i32,
    polar_pentagon: bool,
    south_rhombus: bool,
    edge_hex: bool,
    index: i64,
) -> i64 {
    let ku = |k: i32| f64::from(k) * u;
    let mv = |c: (f64, f64), dx: f64, dy: f64| move5x6_vertex(c.0, c.1, dx, dy);
    let mv3 = |c: (f64, f64), dx: f64, dy: f64| move5x6_vertex3(c.0, c.1, dx, dy);
    let cross =
        |s: (f64, f64), south: bool, left: bool| cross5x6_interruption(s.0, s.1, south, left);
    let mut keep_going = true;
    let north_pentagon = nv == 5 && !south_rhombus;
    let south_pentagon = nv == 5 && south_rhombus;
    // Start from hexagon / pentagon edge; Left-To-Right Scanlines (the eC's note)
    let n_half_rows = pow3(((r_depth - 1) / 2) as u64) as i32;
    let n_rows = 2 * n_half_rows + 1;
    let max_cols = n_rows;
    let min_cols = max_cols - (max_cols - 1) / 2;
    let mut rc = (0.0, 0.0);
    let mut i: i64 = 0;

    // First half
    let (mut r, mut n_cols) = (0, min_cols);
    while keep_going && r <= n_half_rows {
        if index != -1 && i + i64::from(n_cols) <= index {
            i += i64::from(n_cols);
        } else {
            // Computing start of scanline
            rc = if polar_pentagon {
                let k = if south_pentagon { -r } else { r };
                mv(first_centroid, ku(k), ku(k))
            } else if north_pentagon || (edge_hex && !south_rhombus) {
                mv(first_centroid, 0.0, ku(r))
            } else {
                mv(first_centroid, ku(-r), 0.0)
            };

            // Iterating through scanline
            let mut col = 0;
            while keep_going && col < n_cols {
                if index != -1 {
                    // Jump to index
                    col = (index - i) as i32;
                    i = index;
                }

                let centroid = if polar_pentagon {
                    let a = col.min(r);
                    let mut b = col - a;
                    // REVIEW: Using move5x6Vertex() breaks things here (the eC's note). Where
                    // `cross` follows, the step stops on the interruption, which `cross` then
                    // crosses once: a departure from the engine (see `onto_interruption`).
                    let k = ku(if south_pentagon { a } else { -a });
                    let t = if b != 0 {
                        onto_interruption(rc, 0.0, k)
                    } else {
                        mv3(rc, 0.0, k)
                    };
                    if b != 0 {
                        let max_b = n_cols - 1 - 2 * r;
                        let mut c = 0;
                        if b > max_b {
                            c = b - max_b;
                            b = max_b;
                        }
                        let left = north_pentagon || r >= n_half_rows;
                        let i2 = cross(t, south_pentagon, left);
                        let k = if south_pentagon { b } else { -b };
                        let centroid = if c != 0 {
                            mv(i2, ku(k), ku(k))
                        } else {
                            mv3(i2, ku(k), ku(k))
                        };
                        if c != 0 {
                            let i2 = cross(centroid, south_pentagon, left);
                            mv3(
                                i2,
                                if south_pentagon {
                                    if r >= n_half_rows { 0.0 } else { ku(c) }
                                } else {
                                    ku(-c)
                                },
                                if south_pentagon && r >= n_half_rows {
                                    ku(-c)
                                } else {
                                    0.0
                                },
                            )
                        } else {
                            centroid
                        }
                    } else {
                        t
                    }
                } else if edge_hex || north_pentagon {
                    let a = col.min(n_half_rows);
                    let b = col - a;
                    let dy = if south_rhombus { ku(a) } else { 0.0 };
                    let t = if b != 0 {
                        mv(rc, ku(a), dy)
                    } else {
                        mv3(rc, ku(a), dy)
                    };
                    if b != 0 {
                        let i2 = cross(t, south_rhombus, false);
                        mv3(i2, ku(b), if south_rhombus { 0.0 } else { ku(b) })
                    } else {
                        t
                    }
                } else {
                    mv(rc, ku(col), ku(col))
                };

                keep_going = centroid_callback(i, centroid);
                if keep_going {
                    i += 1;
                }
                col += 1;
            }
        }
        r += 1;
        n_cols += 1;
    }

    // Second half
    let (mut r, mut n_cols) = (1, max_cols - 1);
    while keep_going && r <= n_half_rows {
        let mut n = -1;

        // Computing start of scanline
        if north_pentagon || polar_pentagon {
            n_cols -= 1;
        } else if south_pentagon {
            n = (n_cols - r) / 2;
        }

        if index != -1 && i + i64::from(n_cols) <= index {
            // nCols is not always the actual number of indices (the eC's note)
        } else if polar_pentagon {
            let d = if south_pentagon {
                -1.0 - ku(r)
            } else {
                1.0 + ku(r)
            };
            rc = mv(first_centroid, d, d);
        } else if north_pentagon || (edge_hex && !south_rhombus) {
            rc = mv(first_centroid, ku(r), ku(r + n_half_rows));
        } else {
            rc = mv(first_centroid, ku(-n_half_rows), ku(r));
        }

        // Iterating through scanline
        let mut col = 0;
        while keep_going && col < n_cols {
            if index == -1 || index == i {
                let centroid = if north_pentagon || edge_hex || polar_pentagon {
                    let a = col.min(n_half_rows - r);
                    let b = col - a;
                    let dx = ku(if polar_pentagon && south_pentagon {
                        -a
                    } else {
                        a
                    });
                    let dy = if south_rhombus && !polar_pentagon {
                        ku(a)
                    } else {
                        0.0
                    };
                    if b != 0 {
                        let t = if polar_pentagon || edge_hex {
                            let t = mv(rc, dx, dy);
                            cross(t, south_rhombus, south_pentagon && polar_pentagon)
                        } else {
                            mv3(rc, dx, dy)
                        };
                        mv(
                            t,
                            if polar_pentagon { 0.0 } else { ku(b) },
                            if south_rhombus {
                                if polar_pentagon { ku(-b) } else { 0.0 }
                            } else {
                                ku(b)
                            },
                        )
                    } else {
                        mv3(rc, dx, dy)
                    }
                } else {
                    mv(rc, ku(col), ku(col))
                };

                keep_going = centroid_callback(i, centroid);
                if keep_going {
                    i += 1;
                }
            } else if nv == 5 && n > col && i64::from(n - col) < index - i {
                i += i64::from(n) - i64::from(col) + 1;
                col = n;
            } else {
                let ff = (index - i).min(i64::from(n_cols - col)) as i32;
                col += ff - 1;
                i += i64::from(ff);
            }

            if col == n {
                col = n_cols - 1 - n; // Skip interruption
            }
            col += 1;
        }
        r += 1;
        n_cols -= 1;
    }

    if keep_going { -1 } else { i }
}

/// The corner of the interruption from which the scanlines of a southern edge hexagon or a
/// southern pentagon start beyond the middle of the zone: `i1` of `I3HSubZones.ec:505-508` and
/// `:682-685`, built from the first centroid through `ix` and `iy`, its coordinates truncated.
///
/// faithful: in DGGAL's `-O2 -ffast-math` build, `libdggal.so`, BuildID
/// `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, gcc never forms `i0`: it simplified
/// `i1 = { i0.y - 1, iy + 1 + (ix + 1 - i0.x) }`, with `i0 = { FC.x + ((iy + 1) - FC.y), iy + 1 }`
/// (I3HSubZones.ec:507-508, :684-685), to `i1 = { (iy + 1) - 1.0, (FC.y - FC.x) + (ix + 1) }`, at
/// `0x69df` to `0x6a6e` and again at `0x7064` to `0x7116` (`subsd` at `0x6a45` and `0x70c9`,
/// `addpd` at `0x6a5d` and `0x70f7`), and so rounds `FC.y - FC.x` once where the eC rounds `i0.x`
/// and then the inner difference. `ix` and `iy` truncate `FC - 1e-11` as the eC writes it
/// (I3HSubZones.ec:505-506, :682-683; a packed `cvttpd2dq`, then a `paddd` of 1). Both members
/// are computed from the first centroid alone.
fn interruption_corner((x, y): (f64, f64)) -> (f64, f64) {
    let ix_1 = cvtt_i32(x - 1e-11).wrapping_add(1);
    let iy_1 = cvtt_i32(y - 1e-11).wrapping_add(1);
    (f64::from(iy_1) - 1.0, (y - x) + f64::from(ix_1))
}

/// The generator of an even-level zone at an odd relative depth: the cases 2 (a hexagon), 7 and
/// 8 (a northern and a southern edge hexagon), 15 and 16 (a northern and a southern pentagon)
/// and 23 and 24 (the two polar pentagons) of `I3HSubZones.ec:11-98`, which run bottom-to-top
/// scanlines, starting from an edge of the zone, in two halves of `nHalfRows + 1` and
/// `nHalfRows` rows. Ports `generateEvenParentOddDepth` (`I3HSubZones.ec:435-821`), with its
/// index skip: whole scanlines at once, in the second half without the `skip` columns that a
/// pentagon lacks, then a jump within the one that holds the index. The eC's local names are
/// kept.
///
/// faithful: gcc compiled this function on its own, at `0x5960`, and re-associated every
/// division by three in it. In DGGAL's `-O2 -ffast-math` build, `libdggal.so`, BuildID
/// `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, it computes `u3 = fl(T13 * u)` once (`0x5a18` to
/// `0x5a23`) and forms every offset `k * u / 3` whose `k` is an `int` expression, the factors
/// `2`, `-2` and `(southPentagon ? -1 : 1)` folded into it, as the one product `(double)k * u3`,
/// at the 42 call sites `0x5c48` to `0x705a` (I3HSubZones.ec:496-808), where the eC rounds
/// `k * u` and then its third. The coordinate updates follow the same rule: `t -= s * u/3` and
/// `i2 -= s * u/3` (I3HSubZones.ec:580, :584; `t` at `0x622b` and `0x6827`, `i2` at `0x6299` and
/// after `0x6845`), `rc += s * u/3` (I3HSubZones.ec:675, `0x5cd3`) and the point `ii`
/// (I3HSubZones.ec:693, `0x6f76` to `0x6f7e`). That product is `ku3(k)` here. The offsets gcc
/// regrouped further are marked where they occur (I3HSubZones.ec:688, :698-699, :706-707), as is
/// [`interruption_corner`]. The two orders part only in the last bits of a centroid, which the
/// quantisation of a cell centre does not see, so the sub-zone identifiers agree either way; the
/// port follows the library all the same, so that its centroids are the engine's bit for bit.
#[allow(
    clippy::too_many_arguments,
    reason = "the eC's signature, which all four generators share, is kept"
)]
fn generate_even_parent_odd_depth(
    centroid_callback: Callback<'_>,
    first_centroid: (f64, f64),
    r_depth: i32,
    u: f64,
    nv: i32,
    polar_pentagon: bool,
    south_rhombus: bool,
    edge_hex: bool,
    index: i64,
) -> i64 {
    let u3 = T13 * u;
    let ku3 = |k: i32| f64::from(k) * u3;
    let mv = |c: (f64, f64), dx: f64, dy: f64| move5x6_vertex(c.0, c.1, dx, dy);
    let mv3 = |c: (f64, f64), dx: f64, dy: f64| move5x6_vertex3(c.0, c.1, dx, dy);
    // The eC's `(k ? move5x6Vertex : move5x6Vertex3)(...)`.
    let mv_if = |k: i32, c: (f64, f64), dx: f64, dy: f64| {
        if k != 0 {
            mv(c, dx, dy)
        } else {
            mv3(c, dx, dy)
        }
    };
    let cross =
        |s: (f64, f64), south: bool, left: bool| cross5x6_interruption(s.0, s.1, south, left);
    let mut keep_going = true;
    let north_pentagon = nv == 5 && !south_rhombus;
    let south_pentagon = nv == 5 && south_rhombus;
    // The eC's `(southPentagon ? -1 : 1)`.
    let s = if south_pentagon { -1 } else { 1 };
    // Start from hexagon / pentagon edge; Bottom-To-Top Scanlines ("rows" are the vertical
    // scanlines) (the eC's notes)
    let n_half_rows = pow3(((r_depth - 1) / 2) as u64) as i32;
    let n_rows = 2 * n_half_rows + 1;
    let max_cols = n_rows;
    let min_cols = max_cols - (max_cols - 1) / 2;
    let mut i: i64 = 0;

    // First half
    let (mut r, mut n_cols) = (0, min_cols);
    while keep_going && r <= n_half_rows {
        if index != -1 && i + i64::from(n_cols) <= index {
            i += i64::from(n_cols);
        } else {
            // Computing start of scanline
            let rc = if polar_pentagon {
                if r > n_half_rows / 2 {
                    let a = n_half_rows / 2;
                    let b = r - (a + 1);
                    // REVIEW: Avoid calling move5x6Vertex3() before cross5x6Interruption() (the
                    // eC's note)
                    let rc = mv(first_centroid, ku3((a * 2 + 1) * s), ku3((a + 1) * s));
                    let i1 = cross(rc, south_pentagon, south_pentagon);
                    mv3(i1, ku3((1 + b) * s), ku3((1 + b * 2) * s))
                } else {
                    mv(first_centroid, ku3(r * 2 * s), ku3(r * s))
                }
            } else if south_rhombus && (edge_hex || south_pentagon) && r > n_half_rows / 2 {
                let i1 = interruption_corner(first_centroid);
                let b = r - n_half_rows / 2 - 1;
                mv3(i1, ku3(b * 2 + 1), ku3((n_half_rows + 1) / 2 + b))
            } else {
                mv(first_centroid, ku3(r), ku3(r * 2))
            };

            // Iterating through scanline
            let mut col = 0;
            while keep_going && col < n_cols {
                if index != -1 {
                    // Jump to index
                    col = (index - i) as i32;
                    i = index;
                }

                let centroid = if polar_pentagon {
                    let mut col_rem = col;
                    let crossing_left = if south_pentagon {
                        r >= n_half_rows
                    } else {
                        r < n_half_rows
                    };
                    let (start, n) = if r > n_half_rows / 2 {
                        let n = (r - n_half_rows / 2) * 2 - 1;
                        let a = col.min(n);
                        col_rem = col - a;
                        // REVIEW: Avoid calling move5x6Vertex3() before cross5x6Interruption()
                        // (the eC's note)
                        let start = mv_if(col_rem, rc, ku3(a * s), ku3(-a * s));
                        if col_rem != 0 {
                            let n = (n_half_rows - r) / 2 + (n_half_rows - r);
                            (cross(start, south_pentagon, crossing_left), n)
                        } else {
                            (start, n)
                        }
                    } else {
                        (rc, (n_half_rows + r) / 2)
                    };

                    if col_rem != 0 {
                        let a = col_rem.min(n);
                        let mut b = col_rem - a;
                        let mut t = mv_if(b, start, ku3(-a * s), ku3(-2 * a * s));
                        if b != 0 {
                            let odd_row = r & 1 != 0;
                            if !odd_row {
                                b -= 1;
                                t = (t.0 - ku3(s), t.1 - ku3(s));
                            }
                            let mut i2 = cross(t, south_pentagon, crossing_left);
                            if !odd_row {
                                i2 = (i2.0 - ku3(s), i2.1 - ku3(s));
                            }

                            if b > n {
                                let a = b.min(n);
                                b -= a;
                                let t = mv_if(a, i2, ku3(a * s * -2), ku3(-a * s));
                                let i2 = if a != 0 {
                                    cross(t, south_pentagon, crossing_left)
                                } else {
                                    t
                                };
                                mv3(i2, ku3(-b * s), ku3(b * s))
                            } else {
                                mv3(i2, ku3(b * -2 * s), ku3(-b * s))
                            }
                        } else {
                            t
                        }
                    } else {
                        start
                    }
                } else if south_rhombus && (edge_hex || south_pentagon) && r > n_half_rows / 2 {
                    let n = (r - n_half_rows / 2) * 2 - 1;
                    let a = col.min(n);
                    let b = col - a;
                    if b != 0 {
                        // REVIEW: crossing when at intersection corner (the eC's note)
                        let i2 = if !south_pentagon || r < n_half_rows {
                            let t = mv(rc, ku3(-a), ku3(-2 * a));
                            cross(t, true, true)
                        } else {
                            mv3(rc, ku3(-a), ku3(-2 * a))
                        };
                        mv3(i2, ku3(b), ku3(-b))
                    } else {
                        mv3(rc, ku3(-a), ku3(-2 * a))
                    }
                } else if north_pentagon && r > n_half_rows / 2 {
                    let n = 2 * n_half_rows - r;
                    let a = col.min(n);
                    let b = col - a;
                    let t = mv_if(b, rc, ku3(a), ku3(-a));
                    if b != 0 {
                        let i2 = cross(t, false, false);
                        mv3(i2, ku3(b * 2), ku3(b))
                    } else {
                        t
                    }
                } else {
                    mv3(rc, ku3(col), ku3(-col))
                };

                keep_going = centroid_callback(i, centroid);
                if keep_going {
                    i += 1;
                }
                col += 1;
            }
        }
        r += 1;
        n_cols += 1;
    }

    // Second half
    let (mut r, mut n_cols) = (1, max_cols - 1);
    while keep_going && r <= n_half_rows {
        let skip = if nv == 5 { r } else { 0 };
        let mut crosses = false;

        if index != -1 && i + i64::from(n_cols - skip) <= index {
            i += i64::from(n_cols - skip);
        } else {
            // Computing start of scanline
            let rc = if polar_pentagon {
                let a = n_half_rows / 2;
                let b = n_half_rows - (a + 1);
                let aa = r.min(n_half_rows / 2);
                let rc = mv(first_centroid, ku3((a * 2 + 1) * s), ku3((a + 1) * s));
                let i1 = cross(rc, south_pentagon, south_pentagon);
                let rc = mv3(i1, ku3((1 + b + aa * 2) * s), ku3((1 + b * 2 + aa) * s));
                if r > n_half_rows / 2 {
                    // The eC's `bb = (r > nHalfRows/2) ? r - (aa + 1) : 0`, whose test holds
                    // here.
                    let bb = r - (aa + 1);
                    let rc = (rc.0 + ku3(s), rc.1 + ku3(s));
                    let i1 = cross(rc, south_pentagon, south_pentagon);
                    mv3(i1, ku3((1 + bb) * s), ku3((1 + bb * 2) * s))
                } else {
                    rc
                }
            } else if south_rhombus && (edge_hex || south_pentagon) {
                let i1 = interruption_corner(first_centroid);
                let b = n_half_rows - n_half_rows / 2 - 1;
                // faithful: in DGGAL's `-O2 -ffast-math` build, `libdggal.so`, BuildID
                // `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, the offset
                // `u / 3 + b * 2 * u / 3 + r * u / 3` (I3HSubZones.ec:688) is one product
                // `((2b + 1.0) + r) * u3` (`0x7147` to `0x7157`), and
                // `(nHalfRows + 1)/2.0 * u / 3 + b * u / 3 - r * u / 3` is
                // `((nHalfRows + 1) * fl(1/6) + b * fl(1/3)) * u - u3 * r`, the `/2.0` and the
                // `/3` folded into `T16` (`0x70ff` to `0x714f`, `T16` read at `0x7126`), both
                // passed at `0x715b`.
                let dx = (f64::from(2 * b) + 1.0 + f64::from(r)) * u3;
                let dy =
                    (f64::from(n_half_rows + 1) * T16 + f64::from(b) * T13) * u - u3 * f64::from(r);
                mv3(i1, dx, dy)
            } else if !south_rhombus && edge_hex && r > n_half_rows / 2 {
                // The eC's `i`, a name the index already bears here.
                let row = r - n_half_rows / 2 - 1;
                let ii = (
                    first_centroid.0 + ku3(n_half_rows * 2),
                    first_centroid.1 + ku3(n_half_rows * 2),
                );
                let iy = cvtt_i32(ii.1 - 1e-11);
                // faithful: in DGGAL's `-O2 -ffast-math` build, `libdggal.so`, BuildID
                // `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, the start point
                // `t = { iy + 2 - (ii.y - iy), ii.x }` (I3HSubZones.ec:695) is formed as
                // `(iy + 2) + (iy - ii.y)` (`0x6f9a` to `0x6fb4`), the same double; the offset
                // `-(nHalfRows - 1)/2.0 * u / 3 + i * u / 3` (I3HSubZones.ec:698) is
                // `((1 - nHalfRows) * fl(1/6) + i * fl(1/3)) * u` (`0x6fe0` to `0x6ffb`), and
                // `u / 3 + i * 2 * u / 3` (I3HSubZones.ec:699) is `(2i) * u3 + u3` (`0x7004`),
                // both passed at `0x7008`.
                let t = (f64::from(iy + 2) + (f64::from(iy) - ii.1), ii.0);
                crosses = true;
                mv(
                    t,
                    (f64::from(1 - n_half_rows) * T16 + f64::from(row) * T13) * u,
                    ku3(2 * row) + u3,
                )
            } else {
                // faithful: in DGGAL's `-O2 -ffast-math` build, `libdggal.so`, BuildID
                // `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, the two sums
                // `nHalfRows * u / 3 + r * 2 * u / 3` and `nHalfRows * 2 * u / 3 + r * u / 3`
                // (I3HSubZones.ec:706-707) are each one product, `(2r + nHalfRows) * u3` and
                // `(2 nHalfRows + r) * u3` (`0x6da0` to `0x6df5`, passed at `0x6df9`).
                let dx = (f64::from(2 * r) + f64::from(n_half_rows)) * u3;
                let dy = (f64::from(2 * n_half_rows) + f64::from(r)) * u3;
                mv(first_centroid, dx, dy)
            };

            // Iterating through scanline
            let mut col = 0;
            while keep_going && col < n_cols - skip {
                if index != -1 {
                    col = (index - i) as i32;
                    i = index;
                }

                let centroid = if polar_pentagon {
                    let mut b = col;
                    let mut n = n_half_rows - 2 * r;
                    let t = if r <= n_half_rows / 2 {
                        let a = col.min(n);
                        b = col - a;
                        let t = mv_if(b, rc, ku3(s * a), ku3(s * -a));
                        if b != 0 {
                            n = r;
                            cross(t, south_pentagon, south_pentagon)
                        } else {
                            t
                        }
                    } else {
                        n = n_half_rows - r;
                        rc
                    };

                    if b != 0 {
                        if b > n {
                            let a = b.min(n);
                            b -= a;
                            let i2 = mv3(t, ku3(s * a * 2), ku3(s * a));
                            if b > n {
                                let a = b.min(n);
                                b -= a;
                                // REVIEW: move5x6Vertex3() currently crosses on arrival, so
                                // should not be used followed by cross5x6Interruption() ? (the
                                // eC's note)
                                let t = mv(i2, ku3(s * a), ku3(s * a * 2));
                                let i2 = cross(t, south_pentagon, south_pentagon);
                                mv3(i2, ku3(s * -b), ku3(s * b))
                            } else {
                                mv3(i2, ku3(s * b), ku3(s * b * 2))
                            }
                        } else {
                            mv3(t, ku3(s * b * 2), ku3(s * b))
                        }
                    } else {
                        t
                    }
                } else if crosses {
                    mv(rc, ku3(col * 2), ku3(col))
                } else if south_rhombus && (edge_hex || south_pentagon) {
                    let n = if south_pentagon {
                        n_half_rows - r
                    } else {
                        n_cols - (n_half_rows - 2 * r) - 1
                    };
                    let a = col.min(n);
                    let b = (col - n).max(0);
                    if b != 0 {
                        let i2 = if south_pentagon {
                            mv3(rc, ku3(-a), ku3(-2 * a))
                        } else {
                            let t = mv(rc, ku3(-a), ku3(-2 * a));
                            cross(t, true, true)
                        };
                        mv3(i2, ku3(b), ku3(-b))
                    } else {
                        mv3(rc, ku3(-a), ku3(-2 * a))
                    }
                } else if north_pentagon {
                    let n = n_half_rows - r;
                    let a = col.min(n);
                    let b = col - a;
                    let t = mv3(rc, ku3(a), ku3(-a));
                    if b != 0 {
                        mv3(t, ku3(b * 2), ku3(b))
                    } else {
                        t
                    }
                } else {
                    mv3(rc, ku3(col), ku3(-col))
                };

                keep_going = centroid_callback(i, centroid);
                if keep_going {
                    i += 1;
                }
                col += 1;
            }
        }
        r += 1;
        n_cols -= 1;
    }

    if keep_going { -1 } else { i }
}

/// The generator of an odd-level zone at an even relative depth: the cases 4 (a hexagon), 11 and
/// 12 (a northern and a southern edge hexagon), 19 and 20 (a northern and a southern pentagon)
/// and 27 and 28 (the two polar pentagons) of `I3HSubZones.ec:11-98`, which run bottom-to-top
/// scanlines, starting from a vertex of the zone, in three sections: a first cap of `nCapRows`
/// rows, each three columns longer than the last; a main section of `nMidRows` rows; and a
/// second cap, each row three columns shorter, which for a pentagon ends `endCapSkip` rows
/// early. Ports `generateOddParentEvenDepth` (`I3HSubZones.ec:1366-1767`), with its index skip:
/// whole scanlines at once, in the main section without the `skip` columns that a pentagon
/// lacks, then a jump within the one that holds the index, and in the second cap over the
/// columns that a pentagon leaves out. The eC's local names are kept, and `nCols` runs on from
/// one section to the next, as in the eC.
///
/// faithful: gcc compiled this function on its own, at `0x3cd0`, and re-associated every
/// division by three in it. In DGGAL's `-O2 -ffast-math` build, `libdggal.so`, BuildID
/// `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, it computes `u3 = fl(T13 * u)` once (`0x3db8` to
/// `0x3dc8`) and forms every offset `k * u / 3` whose `k` is an `int` expression, the factors
/// `2`, `-2` and `(southPentagon ? -1 : 1)` folded into it, as the one product `(double)k * u3`,
/// at the call sites from `0x3f99` to `0x58ac` (I3HSubZones.ec:1453-1739), where the eC rounds
/// `k * u` and then its third. That product is `ku3(k)` here. Every offset `k * u` is one
/// product `(double)k * u`, as the eC writes it (`ku(k)`), and the sums of
/// I3HSubZones.ec:1497-1498, :1509 and :1696, whose terms follow the two rules, are added in the
/// eC's order (`0x5139`, `0x5018`, `0x5774`). The offsets gcc regrouped further are marked where
/// they occur (I3HSubZones.ec:1517-1518, :1527-1528), as are the two whose factor is a double
/// (I3HSubZones.ec:1460, :1574) and the term gcc dropped (I3HSubZones.ec:1671-1673). The two
/// orders part only in the last bits of a centroid, which the quantisation of a cell centre does
/// not see, so the sub-zone identifiers agree either way; the port follows the library all the
/// same, so that its centroids are the engine's bit for bit.
#[allow(
    clippy::too_many_arguments,
    reason = "the eC's signature, which all four generators share, is kept"
)]
fn generate_odd_parent_even_depth(
    centroid_callback: Callback<'_>,
    first_centroid: (f64, f64),
    r_depth: i32,
    u: f64,
    nv: i32,
    polar_pentagon: bool,
    south_rhombus: bool,
    edge_hex: bool,
    index: i64,
) -> i64 {
    let u3 = T13 * u;
    let ku = |k: i32| f64::from(k) * u;
    let ku3 = |k: i32| f64::from(k) * u3;
    let mv = |c: (f64, f64), dx: f64, dy: f64| move5x6_vertex(c.0, c.1, dx, dy);
    let cross =
        |s: (f64, f64), south: bool, left: bool| cross5x6_interruption(s.0, s.1, south, left);
    let mut keep_going = true;
    let north_pentagon = nv == 5 && !south_rhombus;
    let south_pentagon = nv == 5 && south_rhombus;
    // The eC's `(southPentagon ? -1 : 1)`.
    let s = if south_pentagon { -1 } else { 1 };
    // Start from hexagon / pentagon vertex; Bottom-To-Top Scanlines ("rows" are the vertical
    // scanlines) (the eC's notes)
    let n_cap_rows = pow3(((r_depth - 2) / 2) as u64) as i32;
    let n_mid_rows = 2 * n_cap_rows + 1;
    let end_cap_skip = if nv == 5 { (n_cap_rows + 1) / 2 } else { 0 };
    let min_cols = 1;
    let mut i: i64 = 0;

    // First cap
    let (mut r, mut n_cols) = (0, min_cols);
    while keep_going && r < n_cap_rows {
        if index != -1 && i + i64::from(n_cols) <= index {
            i += i64::from(n_cols);
        } else {
            // Compute start of scanline
            let rc = if polar_pentagon {
                mv(first_centroid, ku(s * r), ku(s * r))
            } else {
                mv(first_centroid, 0.0, ku(r))
            };

            // Iterate through scanline
            let mut col = 0;
            while keep_going && col < n_cols {
                if index != -1 {
                    col = (index - i) as i32;
                    i = index;
                }

                let centroid = if polar_pentagon {
                    let n = r / 2 + r;
                    let a = col.min(n);
                    let mut b = col - a;
                    let mut t = mv(rc, ku3(s * -a), ku3(s * -2 * a));
                    if b != 0 {
                        if r & 1 != 0 {
                            b -= 1;
                            // faithful: the offset `(southPentagon ? -1 : 1) * -u / 3`
                            // (I3HSubZones.ec:1460) multiplies a double, not an `int`: in
                            // DGGAL's `-O2 -ffast-math` build, `libdggal.so`, BuildID
                            // `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, it is
                            // `-((southPentagon ? -1.0 : 1.0) * u3)` (passed at `0x3fd0`), the
                            // same double as `ku3(-s)`, since `u3` is not zero.
                            t = mv(t, ku3(-s), ku3(-s));
                        }
                        let i2 = cross(t, south_pentagon, !south_pentagon);
                        mv(i2, ku3(s * (-b * 2 - (r & 1))), ku3(s * (-b - (r & 1))))
                    } else {
                        t
                    }
                } else {
                    mv(rc, ku3(col), ku3(-col))
                };

                keep_going = centroid_callback(i, centroid);
                if keep_going {
                    i += 1;
                }
                col += 1;
            }
        }
        r += 1;
        n_cols += 3;
    }

    // Main section
    let mut r = 0;
    while keep_going && r < n_mid_rows {
        let skip = if nv == 5 && r > n_cap_rows {
            r - n_cap_rows
        } else {
            0
        };

        if index != -1 && i + i64::from(n_cols - skip) <= index {
            i += i64::from(n_cols - skip);
        } else {
            // Compute start of scanline
            let rc = if polar_pentagon {
                let t = mv(first_centroid, ku(s * n_cap_rows), ku(s * n_cap_rows));
                let i2 = cross(t, south_pentagon, south_pentagon);
                if r != 0 {
                    let (a, b) = (r >> 1, r & 1);
                    mv(i2, ku(s * a) + ku3(s * b * 2), ku(s * a) + ku3(s * b))
                } else {
                    i2
                }
            } else if south_pentagon && r > n_cap_rows {
                let (a, b) = (r / 2, r & 1);
                // This is the pentagon centroid (the eC's note)
                let t = mv(first_centroid, ku(n_cap_rows), ku(n_cap_rows));
                mv(t, ku(a) + ku3(b), ku(n_cap_rows) - ku3(b))
            } else if (edge_hex || south_pentagon) && south_rhombus && r != 0 {
                let n = r + r / 2;
                // faithful: in DGGAL's `-O2 -ffast-math` build, `libdggal.so`, BuildID
                // `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, the two sums
                // `(nCapRows + r) * 2 * u / 3 - (2 * nCapRows - r) * u / 3` and
                // `(nCapRows + r) * u / 3 + (2 * nCapRows - r) * u / 3`
                // (I3HSubZones.ec:1517-1518) are each one product, `(2(nCapRows + r) -
                // (2 nCapRows - r)) * u3` and `((nCapRows + r) + (2 nCapRows - r)) * u3`
                // (`0x505f` to `0x50c0`, passed at `0x50c4`).
                let t = mv(
                    first_centroid,
                    (f64::from(2 * (n_cap_rows + r)) - f64::from(2 * n_cap_rows - r)) * u3,
                    (f64::from(n_cap_rows + r) + f64::from(2 * n_cap_rows - r)) * u3,
                );
                let i2 = if nv == 5 && r > n_cap_rows {
                    t
                } else {
                    cross(t, true, false)
                };
                mv(i2, ku3(n), ku3(n * 2))
            } else {
                // faithful: in DGGAL's `-O2 -ffast-math` build, `libdggal.so`, BuildID
                // `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, the offset
                // `(r >> 1) * u + (r & 1) * 2 * u / 3` (I3HSubZones.ec:1527) is
                // `(2(r & 1) * fl(1/3) + (r >> 1)) * u`, and
                // `(nCapRows + (r >> 1)) * u + (r & 1) * u / 3` (I3HSubZones.ec:1528) is
                // `((nCapRows + (r >> 1)) + (r & 1) * fl(1/3)) * u`, both with `T13` itself
                // rather than `u3` (`0x51d0` to `0x523e`, `T13` read at `0x5202` and
                // `0x522e`, passed at `0x5242`).
                mv(
                    first_centroid,
                    (f64::from(2 * (r & 1)) * T13 + f64::from(r >> 1)) * u,
                    (f64::from(n_cap_rows + (r >> 1)) + f64::from(r & 1) * T13) * u,
                )
            };

            // Iterate through scanline
            let mut col = 0;
            while keep_going && col < n_cols - skip {
                if index != -1 {
                    col = (index - i) as i32;
                    i = index;
                }

                let centroid = if polar_pentagon {
                    let r2 = if r > n_cap_rows {
                        2 * n_cap_rows - r
                    } else {
                        r
                    };
                    let n = r2 + r2 / 2;
                    let a = col.min(n);
                    let mut b = col - a;
                    let t = if a != 0 {
                        mv(rc, ku3(s * a), ku3(s * -a))
                    } else {
                        rc
                    };
                    if b != 0 {
                        let max_b = if r < n_cap_rows {
                            n_cap_rows - r + (n_cap_rows - r) / 2
                        } else {
                            r - n_cap_rows
                        };
                        let mut c = 0;
                        let crossing_left = if south_pentagon {
                            r >= n_cap_rows
                        } else {
                            r < n_cap_rows
                        };

                        if b > max_b {
                            c = b - max_b;
                            b = max_b;
                        }

                        let i2 = cross(t, south_pentagon, crossing_left);
                        let mut centroid = if r >= n_cap_rows {
                            mv(i2, ku3(s * b * 2), ku3(s * b))
                        } else {
                            mv(i2, ku3(s * -b), ku3(s * -b * 2))
                        };

                        if r < n_cap_rows {
                            if c != 0 {
                                let odd_r = r & 1;

                                if odd_r == 0 {
                                    c -= 1;
                                    // faithful: as at I3HSubZones.ec:1460, the offset
                                    // `(southPentagon ? -1 : 1) * -u / 3`
                                    // (I3HSubZones.ec:1574) is, in DGGAL's `-O2 -ffast-math`
                                    // build, `libdggal.so`, BuildID
                                    // `e75d6ab18f8b460713ebcd21df6f07922fb909c3`,
                                    // `-((southPentagon ? -1.0 : 1.0) * u3)` (passed at
                                    // `0x4959`), the same double as `ku3(-s)`.
                                    centroid = mv(centroid, ku3(-s), ku3(-s));
                                }
                                let i2 = cross(centroid, south_pentagon, north_pentagon);

                                centroid = if c > max_b {
                                    let a = c.min(max_b);
                                    let b = c - a;
                                    let t = mv(
                                        i2,
                                        ku3(s * (-a * 2 - (1 - odd_r))),
                                        ku3(s * (-a - (1 - odd_r))),
                                    );
                                    let i2 = cross(t, south_pentagon, north_pentagon);
                                    mv(i2, ku3(s * -b), ku3(s * b))
                                } else {
                                    mv(
                                        i2,
                                        ku3(s * (-c * 2 - (1 - odd_r))),
                                        ku3(s * (-c - (1 - odd_r))),
                                    )
                                };
                            }
                        } else {
                            let i2 = centroid;
                            centroid = if c > max_b {
                                let a = c.min(max_b);
                                let b = c - a;
                                let t = mv(i2, ku3(s * a), ku3(s * a * 2));
                                let i2 = cross(t, south_pentagon, crossing_left);
                                mv(i2, ku3(s * -b), ku3(s * b))
                            } else {
                                mv(i2, ku3(s * c), ku3(s * c * 2))
                            };
                        }
                        centroid
                    } else {
                        t
                    }
                } else if north_pentagon || (edge_hex && !south_rhombus) {
                    let n = if north_pentagon && r > n_cap_rows {
                        n_cap_rows + n_cap_rows / 2 - (r - n_cap_rows) / 2
                    } else {
                        3 * n_cap_rows - r - (r + 1) / 2
                    };
                    let a = col.min(n);
                    let b = col - a;
                    let t = mv(rc, ku3(a), ku3(-a));
                    if b != 0 {
                        let i2 = if north_pentagon && r > n_cap_rows {
                            t
                        } else {
                            cross(t, false, false)
                        };
                        mv(i2, ku3(b * 2), ku3(b))
                    } else {
                        t
                    }
                } else if (edge_hex || south_pentagon) && south_rhombus && r != 0 {
                    let n = if south_pentagon && r > n_cap_rows {
                        n_cap_rows + n_cap_rows / 2 - (r - n_cap_rows) / 2
                    } else {
                        r + r / 2
                    };
                    let a = col.min(n);
                    let b = col - a;
                    let t = mv(rc, ku3(-a), ku3(-a * 2));
                    if b != 0 {
                        let i2 = if south_pentagon && r >= n_cap_rows {
                            t
                        } else {
                            cross(t, true, true)
                        };
                        mv(i2, ku3(b), ku3(-b))
                    } else {
                        t
                    }
                } else {
                    mv(rc, ku3(col), ku3(-col))
                };

                keep_going = centroid_callback(i, centroid);
                if keep_going {
                    i += 1;
                }
                col += 1;
            }
        }
        r += 1;
        n_cols += if n_cols & 1 != 0 { 1 } else { -1 };
    }

    // Second cap
    let mut r = 0;
    n_cols -= 2;
    // Left as it is where a scanline is skipped, as the eC leaves it; never read there.
    let mut rc = (0.0, 0.0);
    while keep_going && r < n_cap_rows - end_cap_skip {
        let n = if nv == 5 {
            if r >= n_cap_rows / 2 {
                0
            } else {
                n_cap_rows - 2 * (r + 1)
            }
        } else {
            -1
        };

        if index != -1 && i + i64::from(n_cols) <= index {
            // The eC's empty statement: the scanline is skipped by the jumps below.
        } else if polar_pentagon {
            // The eC's `b = r2 & 1`, with `r2 = nMidRows - 1`, is 0, since `nMidRows` is odd.
            let a = (n_mid_rows - 1) >> 1;
            let t = mv(first_centroid, ku(s * n_cap_rows), ku(s * n_cap_rows));
            let i2 = cross(t, south_pentagon, south_pentagon);
            // faithful: the offset `(southPentagon ? -1 : 1) * a * u + (southPentagon ? -1 : 1)
            // * b * 2 * u / 3` and its sibling (I3HSubZones.ec:1671-1673) are each `(s * a) * u`
            // alone: in DGGAL's `-O2 -ffast-math` build, `libdggal.so`, BuildID
            // `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, gcc dropped the term of `b`, which is
            // zero (passed at `0x449a`); adding it would not change the offset either.
            let t = mv(i2, ku(s * a), ku(s * a));
            let i2 = cross(t, south_pentagon, south_pentagon);
            rc = mv(i2, ku(s * (r + 1)), ku(s * (r + 1)));
        } else if south_pentagon {
            // This is the pentagon centroid (the eC's note)
            let t = mv(first_centroid, ku(n_cap_rows), ku(n_cap_rows));
            rc = mv(t, ku(n_mid_rows / 2), ku(n_cap_rows - r - 1));
        } else if edge_hex && !south_rhombus {
            let t = mv(first_centroid, ku3(3 * n_cap_rows), ku3(3 * n_cap_rows * 2));
            let i2 = cross(t, false, false);
            rc = mv(i2, ku(r + 1), ku(r + 1));
        } else if edge_hex && south_rhombus {
            let t = mv(first_centroid, ku3(3 * n_cap_rows * 2), ku3(3 * n_cap_rows));
            let i2 = cross(t, true, false);
            rc = mv(i2, ku3(3 * n_cap_rows), ku3(3 * n_cap_rows * 2) - ku(r + 1));
        } else {
            rc = mv(
                first_centroid,
                ku((n_mid_rows + 1) / 2 + r),
                ku(n_cap_rows + ((n_mid_rows - 1) >> 1)),
            );
        }

        // Iterate through scanline
        let mut col = 0;
        while keep_going && col < n_cols {
            let (mut a, mut b) = (0, 0);

            if nv == 5 {
                a = col.min(n);
                b = col - a;
                if b != 0 {
                    b -= (n_cols / 2 - n / 2) + r + r / 2 + 1;
                }
            }

            if index == -1 || index == i {
                let centroid = if polar_pentagon {
                    let t = mv(rc, ku3(s * a * 2), ku3(s * a));
                    mv(t, ku3(s * b), ku3(s * b * 2))
                } else if north_pentagon {
                    let t = mv(rc, ku3(a), ku3(-a));
                    mv(t, ku3(b * 2), ku3(b))
                } else if south_pentagon {
                    let t = mv(rc, ku3(-a), ku3(-2 * a));
                    mv(t, ku3(b), ku3(-b))
                } else if edge_hex && !south_rhombus {
                    mv(rc, ku3(col * 2), ku3(col))
                } else if edge_hex && south_rhombus {
                    mv(rc, ku3(-col), ku3(-2 * col))
                } else {
                    mv(rc, ku3(col), ku3(-col))
                };

                keep_going = centroid_callback(i, centroid);
                if keep_going {
                    i += 1;
                }
            } else if nv == 5 && n > col && i64::from(n - col) < index - i {
                i += i64::from(n) - i64::from(col) + 1;
                col = n;
            } else {
                let ff = (index - i).min(i64::from(n_cols - col)) as i32;
                col += ff - 1;
                i += i64::from(ff);
            }

            if nv == 5 && col == n {
                col = n_cols - n - 1;
            }
            col += 1;
        }
        r += 1;
        n_cols -= 3;
    }

    if keep_going { -1 } else { i }
}

/// The generator of an even-level zone at an even relative depth: the cases 3 (a hexagon), 9 and
/// 10 (a northern and a southern edge hexagon), 17 and 18 (a northern and a southern pentagon)
/// and 25 and 26 (the two polar pentagons) of `I3HSubZones.ec:11-98`, which run left-to-right
/// scanlines, starting from a vertex of the zone, in three sections: a first cap of `nCapRows`
/// rows, each three columns longer than the last; a main portion of `nMidRows` rows, in which a
/// pentagon lacks `colSkip` columns beyond the middle row; and a second cap, each row three
/// columns shorter, which for a pentagon ends `endCapSkip` rows early. Ports
/// `generateEvenParentEvenDepth` (`I3HSubZones.ec:824-1363`), with its index skip: whole
/// scanlines at once, in the main portion without the columns that a pentagon lacks, then a jump
/// within the one that holds the index, and in the second cap over the columns that a pentagon
/// leaves out. The eC's local names are kept, and `nCols` runs on from one section to the next,
/// as in the eC.
///
/// faithful: gcc compiled this function on its own, at `0x7170`, and kept the eC's order in it
/// but for four sums. It has no division by three, and the library no `u3`. In DGGAL's
/// `-O2 -ffast-math` build, `libdggal.so`, BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`,
/// every offset `k * u` whose `k` is an `int` expression, a sign or a factor such as the `-2` of
/// `r * -2 * u` included, is one product `(double)k * u`, as the eC writes it, at the call sites
/// from `0x7372` to `0x9658` (I3HSubZones.ec:867-1335); that product is `ku(k)` here. The two
/// corners `{ 2 + nCapRows * u, 4 - nCapRows * u }` and `{ 3 - nCapRows * u, 2 + nCapRows * u }`
/// (I3HSubZones.ec:1206, :1211) are formed as written, once for the whole second cap (`0x8a54` to
/// `0x8adc`, the constants at `0x49f20` to `0x49f30`). The four sums of such products that gcc
/// gathered into one product are marked where they occur (I3HSubZones.ec:1041-1042, :1046-1047,
/// :1223-1224, :1255-1256). The two orders part only in the last bits of a centroid, which the
/// quantisation of a cell centre does not see, so the sub-zone identifiers agree either way; the
/// port follows the library all the same, so that its centroids are the engine's bit for bit.
#[allow(
    clippy::too_many_arguments,
    reason = "the eC's signature, which all four generators share, is kept"
)]
fn generate_even_parent_even_depth(
    centroid_callback: Callback<'_>,
    first_centroid: (f64, f64),
    r_depth: i32,
    u: f64,
    nv: i32,
    polar_pentagon: bool,
    south_rhombus: bool,
    edge_hex: bool,
    index: i64,
) -> i64 {
    let ku = |k: i32| f64::from(k) * u;
    let mv = |c: (f64, f64), dx: f64, dy: f64| move5x6_vertex(c.0, c.1, dx, dy);
    let mv3 = |c: (f64, f64), dx: f64, dy: f64| move5x6_vertex3(c.0, c.1, dx, dy);
    // The eC's `(k ? move5x6Vertex : move5x6Vertex3)(...)`.
    let mv_if = |k: i32, c: (f64, f64), dx: f64, dy: f64| {
        if k != 0 {
            mv(c, dx, dy)
        } else {
            mv3(c, dx, dy)
        }
    };
    let cross =
        |s: (f64, f64), south: bool, left: bool| cross5x6_interruption(s.0, s.1, south, left);
    let mut keep_going = true;
    let south_pentagon = nv == 5 && south_rhombus;
    // Start from hexagon / pentagon vertex; Left-To-Right Scanlines (the eC's notes)
    let n_cap_rows = pow3(((r_depth - 2) / 2) as u64) as i32;
    let n_mid_rows = 2 * n_cap_rows + 1;
    let end_cap_skip = if nv == 5 { (n_cap_rows + 1) / 2 } else { 0 };
    let min_cols = 1;
    let mut i: i64 = 0;

    // First cap
    let (mut r, mut n_cols) = (0, min_cols);
    while keep_going && r < n_cap_rows {
        if index != -1 && i + i64::from(n_cols) <= index {
            i += i64::from(n_cols);
        } else {
            // Computing start of scanline
            let rc = if polar_pentagon && south_rhombus {
                if r > n_cap_rows / 2 {
                    let t = mv(first_centroid, ku(-n_cap_rows), ku(-(n_cap_rows - r)));
                    let i2 = cross(t, true, true);
                    mv(i2, 0.0, ku(-(2 * (r - n_cap_rows / 2) - 1)))
                } else {
                    mv3(first_centroid, ku(-2 * r), ku(-r))
                }
            } else if (edge_hex || nv == 5) && !south_rhombus {
                if polar_pentagon {
                    if r > n_cap_rows / 2 {
                        let t = mv(first_centroid, ku(n_cap_rows), ku(n_cap_rows - r));
                        let i2 = cross(t, false, false);
                        mv(i2, 0.0, ku(2 * (r - n_cap_rows / 2) - 1))
                    } else {
                        mv(first_centroid, ku(2 * r), ku(r))
                    }
                } else {
                    mv(first_centroid, ku(-r), ku(r))
                }
            } else {
                // The eC's `r * -1 * u`, the same `int` as `-r`.
                mv(first_centroid, ku(r * -2), ku(-r))
            };

            // Iterating through scanline
            let mut col = 0;
            while keep_going && col < n_cols {
                if index != -1 {
                    // Jump to index
                    col = (index - i) as i32;
                    i = index;
                }

                let centroid = if (edge_hex || nv == 5) && !south_rhombus {
                    if r > n_cap_rows / 2 {
                        let n = if polar_pentagon {
                            2 * r - n_cap_rows
                        } else {
                            n_cap_rows + r
                        };
                        let a = col.min(n);
                        let mut b = col - a;
                        // REVIEW: Avoid calling move5x6Vertex3() before cross5x6Interruption()
                        // (the eC's note)
                        let t = if polar_pentagon {
                            // move5x6Vertex() breaks things here (the eC's note)
                            mv_if(b, rc, 0.0, ku(-a))
                        } else {
                            mv_if(b, rc, ku(a), 0.0)
                        };
                        if b != 0 {
                            let i2 = cross(t, false, polar_pentagon);
                            if polar_pentagon {
                                // The eC's `c = col - (a + (b = Min(b, nb)))`.
                                let nb = 2 * n_cap_rows - r;
                                b = b.min(nb);
                                let c = col - (a + b);
                                let t2 = mv_if(c, i2, ku(-b), ku(-b));
                                if c != 0 {
                                    let i2 = cross(t2, false, true);
                                    mv3(i2, ku(-c), 0.0)
                                } else {
                                    t2
                                }
                            } else {
                                mv3(i2, ku(b), ku(b))
                            }
                        } else {
                            t
                        }
                    } else if polar_pentagon {
                        mv(rc, ku(-col), ku(-col))
                    } else {
                        mv(rc, ku(col), 0.0)
                    }
                } else if (edge_hex || nv == 5) && south_rhombus {
                    if r > n_cap_rows / 2 {
                        let n = if polar_pentagon {
                            2 * r - n_cap_rows
                        } else {
                            n_cap_rows + r
                        };
                        let a = col.min(n);
                        let mut b = col - a;
                        let t = if polar_pentagon {
                            mv_if(b, rc, 0.0, ku(a))
                        } else {
                            mv_if(b, rc, ku(a), ku(a))
                        };
                        if b != 0 {
                            let i2 = cross(t, true, false);
                            if polar_pentagon {
                                // The eC's `c = col - (a + (b = Min(b, nb)))`.
                                let nb = 2 * n_cap_rows - r;
                                b = b.min(nb);
                                let c = col - (a + b);
                                // REVIEW: Avoid calling move5x6Vertex3() before
                                // cross5x6Interruption() (the eC's note)
                                let t2 = mv_if(c, i2, ku(b), ku(b));
                                if c != 0 {
                                    let i2 = cross(t2, true, false);
                                    mv3(i2, ku(c), 0.0)
                                } else {
                                    t2
                                }
                            } else if nv == 5 {
                                // The southern pentagon moves on from `t`, not from the crossing.
                                mv3(t, ku(b), ku(b))
                            } else {
                                mv3(i2, ku(b), 0.0)
                            }
                        } else {
                            t
                        }
                    } else {
                        mv(rc, ku(col), ku(col))
                    }
                } else {
                    mv(rc, ku(col), ku(col))
                };

                keep_going = centroid_callback(i, centroid);
                if keep_going {
                    i += 1;
                }
                col += 1;
            }
        }
        r += 1;
        n_cols += 3;
    }

    // Main portion
    let mut r = 0;
    while keep_going && r < n_mid_rows {
        let col_skip = if nv == 5 && r > n_mid_rows / 2 {
            r - n_mid_rows / 2
        } else {
            0
        };

        if index != -1 && i + i64::from(n_cols - col_skip) <= index {
            i += i64::from(n_cols - col_skip);
        } else {
            // Computing start of scanline
            let rc = if polar_pentagon && south_rhombus {
                let a = r.min(n_cap_rows);
                let b = r - a;
                let t = mv(first_centroid, ku(-n_cap_rows), 0.0);
                let i2 = cross(t, true, true);
                let t2 = mv(i2, ku(-a), ku(-(a / 2 + n_cap_rows)));
                if b != 0 {
                    let i2 = cross(t2, true, true);
                    mv(i2, ku(-(b / 2)), ku(-b))
                } else {
                    t2
                }
            } else if (edge_hex || nv == 5) && !south_rhombus {
                if polar_pentagon {
                    let a = r.min(n_cap_rows);
                    let b = r - a;
                    let t = mv(first_centroid, ku(n_cap_rows), 0.0);
                    let i2 = cross(t, false, false);
                    let t2 = mv(i2, ku(a), ku(a / 2 + n_cap_rows));
                    if b != 0 {
                        let i2 = cross(t2, false, false);
                        mv(i2, ku(b / 2), ku(b))
                    } else {
                        t2
                    }
                } else {
                    // faithful: in DGGAL's `-O2 -ffast-math` build, `libdggal.so`, BuildID
                    // `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, the offsets
                    // `nCapRows * -u + ((r + 1) >> 1) * u` and `nCapRows * u + r * u`
                    // (I3HSubZones.ec:1041-1042) are each one product,
                    // `(((r + 1) >> 1) - nCapRows) * u` and `(r + nCapRows) * u` (`0x880e` to
                    // `0x8852`, passed at `0x8856`).
                    mv(
                        first_centroid,
                        (f64::from((r + 1) >> 1) - f64::from(n_cap_rows)) * u,
                        (f64::from(r) + f64::from(n_cap_rows)) * u,
                    )
                }
            } else {
                // faithful: in DGGAL's `-O2 -ffast-math` build, `libdggal.so`, BuildID
                // `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, the offsets
                // `nCapRows * -2 * u + (r >> 1) * -u` and
                // `nCapRows * -1 * u + (r >> 1) * u + (r & 1) * u` (I3HSubZones.ec:1046-1047)
                // are each one product, `(-2 nCapRows - (r >> 1)) * u` and
                // `(((r & 1) + (r >> 1)) + (-nCapRows)) * u` (`0x7cbc` to `0x7d17`, the double
                // `-nCapRows` converted once at `0x7c0c`, passed at `0x7d1b`).
                mv(
                    first_centroid,
                    (f64::from(-2 * n_cap_rows) - f64::from(r >> 1)) * u,
                    ((f64::from(r & 1) + f64::from(r >> 1)) + f64::from(-n_cap_rows)) * u,
                )
            };

            // Iterating through scanline
            let mut col = 0;
            while keep_going && col < n_cols - col_skip {
                if index != -1 {
                    col = (index - i) as i32;
                    i = index;
                }

                let centroid = if col_skip != 0 && south_pentagon {
                    let jump_col = (n_cols - col_skip) >> 1;
                    if col <= jump_col {
                        if polar_pentagon {
                            mv(rc, ku(-col), 0.0)
                        } else {
                            mv(rc, ku(col), ku(col))
                        }
                    } else if polar_pentagon {
                        // Avoid interruption (the eC's note)
                        let t = mv(rc, ku(-jump_col), 0.0);
                        let i2 = cross(t, true, true);
                        mv(i2, 0.0, ku(-(col - jump_col)))
                    } else {
                        let start_shift = col_skip + jump_col + 1;
                        let extra_cols = col - jump_col - 1;
                        let t = mv(rc, ku(start_shift), ku(start_shift));
                        mv(t, ku(extra_cols), ku(extra_cols))
                    }
                } else if (edge_hex || nv == 5) && !south_rhombus {
                    let n = if polar_pentagon {
                        n_cap_rows + r / 2
                    } else {
                        n_cols - 1 - n_cap_rows - r / 2
                    };
                    let a = col.min(n);
                    let mut b = col - a;
                    let (dx, dy) = if polar_pentagon && r <= n_cap_rows {
                        (0.0, ku(-a))
                    } else {
                        (ku(a), 0.0)
                    };
                    // REVIEW: Issues with move5x6Vertex3() followed by cross5x6Interruption()
                    // (the eC's note)
                    let t = if b == 0 || (nv == 5 && r > n_mid_rows / 2) {
                        mv3(rc, dx, dy)
                    } else {
                        mv(rc, dx, dy)
                    };
                    if b != 0 {
                        if nv == 5 && r > n_mid_rows / 2 {
                            mv(t, if polar_pentagon { 0.0 } else { ku(b) }, ku(b))
                        } else {
                            let i2 = cross(t, false, polar_pentagon);
                            if polar_pentagon {
                                // The eC's `c = col - (a + (b = Min(b, nb)))`.
                                let nb = n_cap_rows - r;
                                b = b.min(nb);
                                let c = col - (a + b);
                                let t2 = mv(i2, ku(-b), ku(-b));
                                if c != 0 {
                                    let i2 = cross(t2, false, true);
                                    mv3(i2, ku(-c), 0.0)
                                } else {
                                    t2
                                }
                            } else {
                                mv(i2, ku(b), ku(b))
                            }
                        }
                    } else {
                        t
                    }
                } else if (edge_hex || nv == 5) && south_rhombus {
                    let n = if polar_pentagon {
                        n_cap_rows + r / 2
                    } else {
                        n_cols - 1 - n_cap_rows - r / 2
                    };
                    let a = col.min(n);
                    let mut b = col - a;
                    let dx = if polar_pentagon && r <= n_cap_rows {
                        0.0
                    } else {
                        ku(a)
                    };
                    if b != 0 {
                        if nv == 5 && r <= n_mid_rows / 2 && !polar_pentagon {
                            let t = mv3(rc, dx, ku(a));
                            mv(t, ku(b), ku(b))
                        } else {
                            // REVIEW: move5x6Vertex3() is required here (the eC's note). The
                            // step stops on the interruption, which `cross` then crosses once: a
                            // departure from the engine (see `onto_interruption`).
                            let t = onto_interruption(rc, dx, ku(a));
                            let i2 = cross(t, true, polar_pentagon && r >= n_mid_rows / 2);
                            if polar_pentagon {
                                // The eC's `c = col - (a + (b = Min(b, nb)))`.
                                let nb = n_cap_rows - r;
                                b = b.min(nb);
                                let c = col - (a + b);
                                let t2 = mv(i2, ku(b), ku(b));
                                if c != 0 {
                                    let i2 = cross(t2, true, r >= n_mid_rows / 2);
                                    if r >= n_mid_rows / 2 {
                                        mv3(i2, 0.0, ku(-c))
                                    } else {
                                        mv3(i2, ku(c), 0.0)
                                    }
                                } else {
                                    t2
                                }
                            } else {
                                mv(i2, ku(b), if nv == 5 { ku(b) } else { 0.0 })
                            }
                        }
                    } else {
                        mv3(rc, dx, ku(a))
                    }
                } else {
                    mv(rc, ku(col), ku(col))
                };

                keep_going = centroid_callback(i, centroid);
                if keep_going {
                    i += 1;
                }
                col += 1;
            }
        }
        r += 1;
        n_cols += if n_cols & 1 != 0 { 1 } else { -1 };
    }

    // Second cap
    let mut r = 0;
    n_cols -= 2;
    // Left as it is where a scanline is skipped, as the eC leaves it; never read there.
    let mut rc = (0.0, 0.0);
    while keep_going && r < n_cap_rows - end_cap_skip {
        let n = if (edge_hex || nv == 5) && r < n_cap_rows / 2 {
            n_cap_rows - 2 * (r + 1)
        } else {
            0
        };

        if index != -1 && i + i64::from(n_cols) <= index {
            // The eC's empty statement: the scanline is skipped by the jumps below.
        } else if polar_pentagon {
            rc = if south_rhombus {
                let i2 = (2.0 + ku(n_cap_rows), 4.0 - ku(n_cap_rows));
                mv(i2, ku(-(r + 1) * 2), ku(-(r + 1)))
            } else {
                let i2 = (3.0 - ku(n_cap_rows), 2.0 + ku(n_cap_rows));
                mv(i2, ku((r + 1) * 2), ku(r + 1))
            };
        } else if (edge_hex || nv == 5) && !south_rhombus {
            rc = if r < n_cap_rows / 2 {
                let a = (r + 1).min(n_cap_rows / 2);
                let b = r + 1 - a;
                // faithful: in DGGAL's `-O2 -ffast-math` build, `libdggal.so`, BuildID
                // `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, the offsets
                // `nCapRows * -u + (nMidRows >> 1) * u + a * 2 * u` and
                // `nCapRows * u + (nMidRows - 1) * u + a * u` (I3HSubZones.ec:1223-1224) are
                // each one product, `(((nMidRows >> 1) - nCapRows) + 2a) * u` and
                // `(((nMidRows - 1) + nCapRows) + a) * u` (`0x9194` to `0x9204`, passed at
                // `0x9208`). The first is the same double in either order, since
                // `nMidRows >> 1` is `nCapRows`.
                let t = mv(
                    first_centroid,
                    ((f64::from(n_mid_rows >> 1) - f64::from(n_cap_rows)) + f64::from(a * 2)) * u,
                    ((f64::from(n_mid_rows - 1) + f64::from(n_cap_rows)) + f64::from(a)) * u,
                );
                if b != 0 {
                    let i2 = cross(t, false, false);
                    mv(i2, ku(b * 2), ku(b))
                } else {
                    t
                }
            } else {
                let a = n_cap_rows;
                let ay = 3 * n_cap_rows - 1 - r;
                let b = 1 + 2 * (r - n_cap_rows / 2);
                let t = mv(first_centroid, ku(a), ku(a));
                let i2 = cross(t, false, false);
                mv3(i2, ku(-ay), ku(b))
            };
        } else if (edge_hex || nv == 5) && south_rhombus && r >= n_cap_rows / 2 {
            let a = n_cap_rows;
            let ay = 2 * n_cap_rows + r + 1;
            let b = 1 + 2 * (r - n_cap_rows / 2);
            let t = mv(first_centroid, 0.0, ku(a));
            let i2 = cross(t, true, false);
            rc = mv3(i2, ku(b), ku(ay));
        } else {
            // faithful: in DGGAL's `-O2 -ffast-math` build, `libdggal.so`, BuildID
            // `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, the offsets
            // `nCapRows * -2 * u + ((nMidRows - 1) >> 1) * -u + (r + 1) * u` and
            // `nCapRows * -1 * u + ((nMidRows - 1) >> 1) * u + (r + 1) * 2 * u`
            // (I3HSubZones.ec:1255-1256) are each one product,
            // `((-2 nCapRows - ((nMidRows - 1) >> 1)) + (r + 1)) * u` and
            // `((-nCapRows + ((nMidRows - 1) >> 1)) + 2(r + 1)) * u` (`0x8b47` to `0x8bb9`,
            // passed at `0x8bbd`). The second is the same double in either order, since
            // `(nMidRows - 1) >> 1` is `nCapRows`.
            rc = mv3(
                first_centroid,
                ((f64::from(-2 * n_cap_rows) - f64::from((n_mid_rows - 1) >> 1))
                    + f64::from(r + 1))
                    * u,
                ((f64::from(-n_cap_rows) + f64::from((n_mid_rows - 1) >> 1))
                    + f64::from((r + 1) * 2))
                    * u,
            );
        }

        // Iterate through scanline
        let mut col = 0;
        while keep_going && col < n_cols {
            if index == -1 || i == index {
                let centroid = if (edge_hex || nv == 5) && !south_rhombus {
                    let a = col.min(n);
                    let mut b = col - a;
                    if a != 0 {
                        if b != 0 {
                            let i2 = if nv == 5 {
                                b -= (n_cols / 2 - n / 2) + r + r / 2 + 1;
                                if polar_pentagon {
                                    let t = mv(rc, ku(a), 0.0);
                                    cross(t, false, false)
                                } else {
                                    mv3(rc, ku(a), 0.0)
                                }
                            } else {
                                let t = mv(rc, ku(a), 0.0);
                                cross(t, false, false)
                            };
                            if polar_pentagon {
                                mv3(i2, 0.0, ku(b))
                            } else {
                                mv3(i2, ku(b), ku(b))
                            }
                        } else {
                            mv3(rc, ku(a), 0.0)
                        }
                    } else {
                        mv3(rc, ku(b), ku(b))
                    }
                } else if (edge_hex || nv == 5) && south_rhombus {
                    let a = col.min(n);
                    let mut b = col - a;
                    if a != 0 {
                        let t = if polar_pentagon {
                            // REVIEW: move5x6Vertex3() breaks things here (the eC's note)
                            mv(rc, ku(-a), 0.0)
                        } else {
                            // REVIEW: Avoid using move5x6Vertex3() before
                            // cross5x6Interruption() (the eC's note)
                            mv_if(b, rc, ku(a), ku(a))
                        };
                        if b != 0 {
                            let i2 = cross(t, true, polar_pentagon);
                            if nv == 5 {
                                b -= (n_cols / 2 - n / 2) + r + r / 2 + 1;
                                if polar_pentagon {
                                    mv3(i2, 0.0, ku(-b))
                                } else {
                                    mv3(i2, ku(b), ku(b))
                                }
                            } else {
                                mv3(i2, ku(b), 0.0)
                            }
                        } else {
                            t
                        }
                    } else {
                        mv3(rc, ku(b), 0.0)
                    }
                } else {
                    mv(rc, ku(col), ku(col))
                };

                keep_going = centroid_callback(i, centroid);
                if keep_going {
                    i += 1;
                }
            } else if nv == 5 && n > col && i64::from(n - col) < index - i {
                i += i64::from(n) - i64::from(col) + 1;
                col = n;
            } else {
                let ff = (index - i).min(i64::from(n_cols - col)) as i32;
                col += ff - 1;
                i += i64::from(ff);
            }

            if nv == 5 && col == n {
                col = n_cols - n - 1;
            }
            col += 1;
        }
        r += 1;
        n_cols -= 3;
    }

    if keep_going { -1 } else { i }
}

/// Every sub-zone `depth` levels below zone `id`, in the engine's order but for the departure
/// of the module note, or `None` where they would lie beyond level 33 or number `2^28` or more,
/// where the engine's `getSubZones` answers nothing either. Ports `getI3HSubZoneCentroids`
/// (`I3HSubZones.ec:1775-1800`), which stores each centroid at its index in an array of the
/// count's length, and `getSubZones` (`dggrs.ec:195-222`), which quantises each at the sub-zone
/// level. `Grid` settles depth 0 before it asks.
pub(super) fn sub_zones(id: u64, depth: u8) -> Option<Vec<u64>> {
    let sz_level = sub_zone_level(id, depth)?;
    let n_sub_zones = count(id, depth)?;
    // Each centroid is 16 bytes and array memory allocation currently does not support more
    // than 4G (the eC's note)
    if n_sub_zones >= 1 << (32 - 4) {
        return None;
    }
    if depth == 0 {
        return Some(vec![id]);
    }
    let mut centroids = vec![(0.0, 0.0); n_sub_zones as usize];
    iterate(
        id,
        u32::from(depth),
        &mut |index, centroid| {
            if let Some(slot) = centroids.get_mut(index as usize) {
                *slot = centroid;
            }
            true
        },
        -1,
    );
    Some(
        centroids
            .into_iter()
            .map(|c| quantize(sz_level, c))
            .collect(),
    )
}

/// The sub-zone at `index` in the order of [`sub_zones`], found through the generator's index
/// skip, without materialising the order, in steps of the order of the number of scanlines,
/// `3^(depth / 2)`; `None` where `index` is not below the count or beyond level 33. Ports
/// `getSubZoneAtIndex` (`RI3H.ec:212-227`) with its callback
/// `findByIndex` (`RI3H.ec:2388-2392`), except that index 0 is found by the same walk rather than
/// by `getFirstSubZone`, which gives the same zone at every depth but 0; at depth 0 this answers
/// the zone itself, as [`first_sub_zone`] does. [`super::HexA3::sub_zone_at_index`] answers
/// through it.
pub(super) fn sub_zone_at_index(id: u64, depth: u8, index: u64) -> Option<u64> {
    let sz_level = sub_zone_level(id, depth)?;
    if index >= count(id, depth)? {
        return None;
    }
    if depth == 0 {
        return Some(id);
    }
    let mut centroid = (0.0, 0.0);
    iterate(
        id,
        u32::from(depth),
        &mut |_, c| {
            centroid = c;
            false
        },
        index as i64,
    );
    Some(quantize(sz_level, centroid))
}

#[cfg(test)]
mod tests {
    //! Every expected identifier and count below is the shipped engine's own answer, DGGAL
    //! v0.0.6 on ISEA3H, read through `getSubZones`, `getFirstSubZone` and `countSubZones`;
    //! sub-zone lists are compared in the engine's order.
    use super::*;
    use crate::indexings::I3h;
    use crate::projections::Isea;
    use crate::topologies::HexA3;
    use crate::{Grid, GridConfig, ZoneId};

    fn isea3h() -> Grid<Isea, HexA3, I3h> {
        Grid::new(GridConfig::default(), "ISEA3H").unwrap()
    }

    fn id(grid: &Grid<Isea, HexA3, I3h>, text: &str) -> ZoneId {
        grid.zone_from_text(text).unwrap().id()
    }

    fn texts(grid: &Grid<Isea, HexA3, I3h>, ids: &[ZoneId]) -> Vec<String> {
        ids.iter().map(|&z| grid.zone(z).text_id()).collect()
    }

    #[test]
    fn counts_match_the_engine() {
        let grid = isea3h();
        for (zone, want) in [
            // A pentagon at resolution 0, and a hexagon and a pentagon at resolution 1.
            (
                "A0-0-A",
                [
                    (1, 6),
                    (2, 11),
                    (3, 31),
                    (6, 631),
                    (13, 1_330_426),
                    (14, 3_987_631),
                    (17, 107_633_206),
                    (18, 322_866_811),
                    (33, 4_632_550_579_746_406),
                ],
            ),
            (
                "A6-0-C",
                [
                    (1, 7),
                    (2, 13),
                    (3, 37),
                    (6, 757),
                    (13, 1_596_511),
                    (14, 4_785_157),
                    (17, 129_159_847),
                    (18, 387_440_173),
                    (32, 1_853_020_231_898_563),
                ],
            ),
            (
                "A6-0-B",
                [
                    (1, 6),
                    (2, 11),
                    (3, 31),
                    (6, 631),
                    (13, 1_330_426),
                    (14, 3_987_631),
                    (17, 107_633_206),
                    (18, 322_866_811),
                    (32, 1_544_183_526_582_136),
                ],
            ),
        ] {
            let z = id(&grid, zone);
            for (depth, count) in want {
                assert_eq!(
                    grid.count_sub_zones(z, depth),
                    Ok(count),
                    "{zone} at {depth}"
                );
            }
        }
    }

    /// The seven cases of the odd parent at an odd depth (`I3HSubZones.ec:16`, `:29`, `:31`,
    /// `:51`, `:53`, `:73` and `:76`): a hexagon, a northern and a southern edge hexagon, a
    /// northern and a southern pentagon, and the two polar pentagons.
    #[test]
    fn the_order_at_depth_1_matches_the_engine() {
        let grid = isea3h();
        for (zone, want) in [
            (
                "A6-0-C",
                &[
                    "B6-2-A", "B8-2-A", "B6-1-A", "B6-5-A", "B8-1-A", "B6-4-A", "B6-8-A",
                ][..],
            ),
            (
                "B6-1-B",
                &[
                    "C4-35-A", "C6-4-A", "C4-3E-A", "C6-3-A", "C6-D-A", "C6-2-A", "C6-C-A",
                ],
            ),
            (
                "B5-3-B",
                &[
                    "C3-4E-A", "C5-12-A", "C3-4D-A", "C5-1B-A", "C5-1C-A", "C5-24-A", "C5-25-A",
                ],
            ),
            (
                "A6-0-B",
                &["B4-8-A", "B6-1-A", "B5-2-A", "B6-0-A", "B6-4-A", "B6-3-A"],
            ),
            (
                "A5-0-B",
                &["B4-6-A", "B5-1-A", "B3-8-A", "B5-0-A", "B5-4-A", "B5-3-A"],
            ),
            (
                "AA-0-B",
                &["B2-2-A", "B0-2-A", "B4-2-A", "BA-0-A", "B8-2-A", "B6-2-A"],
            ),
            (
                "AB-0-B",
                &["B9-6-A", "B1-6-A", "B7-6-A", "BB-0-A", "B3-6-A", "B5-6-A"],
            ),
        ] {
            let subs = grid.sub_zones(id(&grid, zone), 1).unwrap();
            assert_eq!(texts(&grid, &subs), want, "{zone}");
        }
    }

    /// Depth 3 reaches both halves of the generator, the interruption crossings and, for the
    /// southern pentagon, the columns it skips.
    #[test]
    fn the_order_at_depth_3_matches_the_engine() {
        let grid = isea3h();
        for (zone, want) in [
            (
                "B6-4-C",
                &[
                    "D6-102-A", "D6-11E-A", "D6-13A-A", "D6-156-A", "D6-101-A", "D6-11D-A",
                    "D6-139-A", "D6-155-A", "D6-171-A", "D6-100-A", "D6-11C-A", "D6-138-A",
                    "D6-154-A", "D6-170-A", "D6-18C-A", "D6-FF-A", "D6-11B-A", "D6-137-A",
                    "D6-153-A", "D6-16F-A", "D6-18B-A", "D6-1A7-A", "D6-11A-A", "D6-136-A",
                    "D6-152-A", "D6-16E-A", "D6-18A-A", "D6-1A6-A", "D6-135-A", "D6-151-A",
                    "D6-16D-A", "D6-189-A", "D6-1A5-A", "D6-150-A", "D6-16C-A", "D6-188-A",
                    "D6-1A4-A",
                ][..],
            ),
            (
                "B6-1-B",
                &[
                    "D4-1AD-A", "D4-1AE-A", "D4-1AF-A", "D6-C-A", "D4-1C8-A", "D4-1C9-A",
                    "D4-1CA-A", "D6-B-A", "D6-27-A", "D4-1E3-A", "D4-1E4-A", "D4-1E5-A", "D6-A-A",
                    "D6-26-A", "D6-42-A", "D4-1FE-A", "D4-1FF-A", "D4-200-A", "D6-9-A", "D6-25-A",
                    "D6-41-A", "D6-5D-A", "D4-21A-A", "D4-21B-A", "D6-8-A", "D6-24-A", "D6-40-A",
                    "D6-5C-A", "D4-236-A", "D6-7-A", "D6-23-A", "D6-3F-A", "D6-5B-A", "D6-6-A",
                    "D6-22-A", "D6-3E-A", "D6-5A-A",
                ],
            ),
            (
                "A5-0-B",
                &[
                    "C4-36-A", "C4-40-A", "C4-4A-A", "C5-3-A", "C3-3E-A", "C4-3F-A", "C4-49-A",
                    "C5-2-A", "C5-C-A", "C3-3D-A", "C3-47-A", "C4-48-A", "C5-1-A", "C5-B-A",
                    "C5-15-A", "C3-3C-A", "C3-46-A", "C3-50-A", "C5-0-A", "C5-A-A", "C5-14-A",
                    "C5-1E-A", "C3-45-A", "C3-4F-A", "C5-9-A", "C5-13-A", "C5-1D-A", "C3-4E-A",
                    "C5-12-A", "C5-1C-A", "C5-1B-A",
                ],
            ),
            (
                "AA-0-B",
                &[
                    "C2-6-A", "C0-1A-A", "C0-10-A", "C0-6-A", "C2-10-A", "C2-7-A", "C0-11-A",
                    "C0-7-A", "C8-1A-A", "C2-1A-A", "C2-11-A", "C2-8-A", "C0-8-A", "C8-11-A",
                    "C8-10-A", "C4-6-A", "C4-7-A", "C4-8-A", "CA-0-A", "C8-8-A", "C8-7-A",
                    "C8-6-A", "C4-10-A", "C4-11-A", "C6-8-A", "C6-11-A", "C6-1A-A", "C4-1A-A",
                    "C6-7-A", "C6-10-A", "C6-6-A",
                ],
            ),
            (
                "AB-0-B",
                &[
                    "C9-36-A", "C9-40-A", "C9-4A-A", "C1-36-A", "C7-4A-A", "C9-3F-A", "C9-49-A",
                    "C1-3F-A", "C1-40-A", "C7-40-A", "C7-49-A", "C9-48-A", "C1-48-A", "C1-49-A",
                    "C1-4A-A", "C7-36-A", "C7-3F-A", "C7-48-A", "CB-0-A", "C3-48-A", "C3-3F-A",
                    "C3-36-A", "C5-4A-A", "C5-49-A", "C5-48-A", "C3-49-A", "C3-40-A", "C5-40-A",
                    "C5-3F-A", "C3-4A-A", "C5-36-A",
                ],
            ),
        ] {
            let subs = grid.sub_zones(id(&grid, zone), 3).unwrap();
            assert_eq!(texts(&grid, &subs), want, "{zone}");
        }
    }

    /// The seven cases of the even parent at an odd depth (`I3HSubZones.ec:18`, `:33`, `:37`,
    /// `:55`, `:57`, `:78` and `:80`): a hexagon, a northern and a southern edge hexagon, a
    /// northern and a southern pentagon, and the two polar pentagons.
    #[test]
    fn the_even_parent_order_at_depth_1_matches_the_engine() {
        let grid = isea3h();
        for (zone, want) in [
            (
                "B6-5-A",
                &[
                    "B6-1-D", "B6-1-C", "B6-4-C", "B6-5-B", "B6-2-D", "B6-5-D", "B6-5-C",
                ][..],
            ),
            (
                "B6-1-A",
                &[
                    "B4-5-D", "B4-5-C", "B4-8-C", "B6-1-B", "B6-1-C", "B6-0-C", "B6-1-D",
                ],
            ),
            (
                "B5-3-A",
                &[
                    "B3-7-D", "B3-7-C", "B5-3-D", "B5-3-B", "B3-8-D", "B5-3-C", "B5-0-D",
                ],
            ),
            (
                "A2-0-A",
                &["A0-0-D", "A0-0-C", "A1-0-C", "A2-0-B", "A2-0-C", "A2-0-D"],
            ),
            (
                "A3-0-A",
                &["A1-0-D", "A1-0-C", "A3-0-D", "A3-0-B", "A2-0-D", "A3-0-C"],
            ),
            (
                "AA-0-A",
                &["A0-0-C", "A8-0-C", "A2-0-C", "AA-0-B", "A6-0-C", "A4-0-C"],
            ),
            (
                "AB-0-A",
                &["A9-0-D", "A1-0-D", "A7-0-D", "AB-0-B", "A3-0-D", "A5-0-D"],
            ),
        ] {
            let subs = grid.sub_zones(id(&grid, zone), 1).unwrap();
            assert_eq!(texts(&grid, &subs), want, "{zone}");
        }
    }

    /// Depth 3 of an even parent reaches both halves of the generator: for the edge hexagons
    /// and the southern pentagon, the scanlines that start from a corner of the interruption,
    /// and for the polar pentagons, the scanlines that cross it.
    #[test]
    fn the_even_parent_order_at_depth_3_matches_the_engine() {
        let grid = isea3h();
        for (zone, want) in [
            (
                "B6-5-A",
                &[
                    "C6-16-B", "C6-D-D", "C6-D-C", "C6-E-B", "C6-16-D", "C6-16-C", "C6-17-B",
                    "C6-E-D", "C6-E-C", "C6-1F-C", "C6-20-B", "C6-17-D", "C6-17-C", "C6-18-B",
                    "C6-F-D", "C6-29-B", "C6-20-D", "C6-20-C", "C6-21-B", "C6-18-D", "C6-18-C",
                    "C6-19-B", "C6-29-C", "C6-2A-B", "C6-21-D", "C6-21-C", "C6-22-B", "C6-19-D",
                    "C6-2A-D", "C6-2A-C", "C6-2B-B", "C6-22-D", "C6-22-C", "C6-34-B", "C6-2B-D",
                    "C6-2B-C", "C6-2C-B",
                ][..],
            ),
            (
                "B6-1-A",
                &[
                    "C4-34-B", "C4-2B-D", "C4-2B-C", "C4-2C-B", "C4-34-D", "C4-34-C", "C4-35-B",
                    "C4-2C-D", "C4-2C-C", "C4-3D-C", "C4-3E-B", "C4-35-D", "C4-35-C", "C6-4-B",
                    "C6-4-C", "C4-47-B", "C4-3E-D", "C4-3E-C", "C6-3-B", "C6-3-C", "C6-4-D",
                    "C6-E-B", "C4-47-C", "C6-2-B", "C6-2-C", "C6-3-D", "C6-D-B", "C6-D-C",
                    "C6-1-C", "C6-2-D", "C6-C-B", "C6-C-C", "C6-D-D", "C6-B-B", "C6-B-C", "C6-C-D",
                    "C6-16-B",
                ],
            ),
            (
                "B5-3-A",
                &[
                    "C3-4C-B", "C3-43-D", "C3-43-C", "C3-44-B", "C3-4C-D", "C3-4C-C", "C3-4D-B",
                    "C3-44-D", "C3-44-C", "C5-24-D", "C5-24-B", "C3-4D-D", "C3-4D-C", "C3-4E-B",
                    "C3-45-D", "C5-2E-B", "C5-24-C", "C5-1B-D", "C5-1B-B", "C3-4E-D", "C3-4E-C",
                    "C3-4F-B", "C5-25-D", "C5-25-B", "C5-1B-C", "C5-12-D", "C5-12-B", "C3-4F-D",
                    "C5-25-C", "C5-1C-D", "C5-1C-B", "C5-12-C", "C5-9-D", "C5-26-B", "C5-1C-C",
                    "C5-13-D", "C5-13-B",
                ],
            ),
            (
                "A3-0-A",
                &[
                    "B1-7-B", "B1-4-D", "B1-4-C", "B1-5-B", "B1-7-D", "B1-7-C", "B1-8-B", "B1-5-D",
                    "B1-5-C", "B3-3-D", "B3-3-B", "B1-8-D", "B1-8-C", "B2-6-B", "B2-3-D", "B3-7-B",
                    "B3-3-C", "B3-0-D", "B3-0-B", "B2-6-D", "B2-6-C", "B2-7-B", "B3-4-D", "B3-4-B",
                    "B3-0-C", "B3-1-B", "B2-7-D", "B3-4-C", "B3-1-D", "B3-1-C", "B3-5-B",
                ],
            ),
            (
                "A2-0-A",
                &[
                    "B0-7-B", "B0-4-D", "B0-4-C", "B0-5-B", "B0-7-D", "B0-7-C", "B0-8-B", "B0-5-D",
                    "B0-5-C", "B1-1-C", "B1-2-B", "B0-8-D", "B0-8-C", "B2-1-B", "B2-1-C", "B1-5-B",
                    "B1-2-D", "B1-2-C", "B2-0-B", "B2-0-C", "B2-1-D", "B2-5-B", "B1-5-C", "B2-3-B",
                    "B2-0-D", "B2-4-B", "B2-4-C", "B2-3-D", "B2-3-C", "B2-4-D", "B2-7-B",
                ],
            ),
            (
                "AA-0-A",
                &[
                    "B0-5-B", "B0-1-C", "B8-5-C", "B8-5-B", "B0-5-C", "B0-2-D", "B0-2-B", "B8-2-D",
                    "B8-1-C", "B2-1-C", "B2-2-B", "B0-2-C", "B8-2-C", "B8-2-B", "B6-5-C", "B2-5-B",
                    "B2-2-D", "B2-2-C", "BA-0-B", "B6-2-C", "B6-2-D", "B6-5-B", "B2-5-C", "B4-2-B",
                    "B4-2-C", "B6-2-B", "B6-1-C", "B4-1-C", "B4-2-D", "B4-5-C", "B4-5-B",
                ],
            ),
            (
                "AB-0-A",
                &[
                    "B9-7-B", "B9-7-D", "B1-3-D", "B1-7-B", "B9-3-D", "B9-6-C", "B1-6-B", "B1-6-C",
                    "B1-7-D", "B7-7-D", "B9-6-B", "B9-6-D", "B1-6-D", "B3-6-B", "B3-3-D", "B7-7-B",
                    "B7-6-C", "B7-6-D", "BB-0-B", "B3-6-D", "B3-6-C", "B3-7-B", "B7-3-D", "B7-6-B",
                    "B5-6-D", "B5-6-B", "B3-7-D", "B5-7-D", "B5-6-C", "B5-3-D", "B5-7-B",
                ],
            ),
            // The polar pentagons at levels 16 and 22, where calling `move5x6Vertex` rather
            // than `move5x6Vertex3` at I3HSubZones.ec:590 and :758 decides a sub-zone.
            (
                "IA-0-A",
                &[
                    "J0-99C5-B",
                    "J0-4CE1-C",
                    "J8-99C5-C",
                    "J8-99C5-B",
                    "J0-99C5-C",
                    "J0-4CE2-D",
                    "J0-4CE2-B",
                    "J8-4CE2-D",
                    "J8-4CE1-C",
                    "J2-4CE1-C",
                    "J2-4CE2-B",
                    "J0-4CE2-C",
                    "J8-4CE2-C",
                    "J8-4CE2-B",
                    "J6-99C5-C",
                    "J2-99C5-B",
                    "J2-4CE2-D",
                    "J2-4CE2-C",
                    "JA-0-B",
                    "J6-4CE2-C",
                    "J6-4CE2-D",
                    "J6-99C5-B",
                    "J2-99C5-C",
                    "J4-4CE2-B",
                    "J4-4CE2-C",
                    "J6-4CE2-B",
                    "J6-4CE1-C",
                    "J4-4CE1-C",
                    "J4-4CE2-D",
                    "J4-99C5-C",
                    "J4-99C5-B",
                ],
            ),
            (
                "LB-0-A",
                &[
                    "M9-41C2149CF1-B",
                    "M9-41C2149CF1-D",
                    "M1-41C20C80FF-D",
                    "M1-41C2149CF1-B",
                    "M9-41C20C80FF-D",
                    "M9-41C2149CF0-C",
                    "M1-41C2149CF0-B",
                    "M1-41C2149CF0-C",
                    "M1-41C2149CF1-D",
                    "M7-41C2149CF1-D",
                    "M9-41C2149CF0-B",
                    "M9-41C2149CF0-D",
                    "M1-41C2149CF0-D",
                    "M3-41C2149CF0-B",
                    "M3-41C20C80FF-D",
                    "M7-41C2149CF1-B",
                    "M7-41C2149CF0-C",
                    "M7-41C2149CF0-D",
                    "MB-0-B",
                    "M3-41C2149CF0-D",
                    "M3-41C2149CF0-C",
                    "M3-41C2149CF1-B",
                    "M7-41C20C80FF-D",
                    "M7-41C2149CF0-B",
                    "M5-41C2149CF0-D",
                    "M5-41C2149CF0-B",
                    "M3-41C2149CF1-D",
                    "M5-41C2149CF1-D",
                    "M5-41C2149CF0-C",
                    "M5-41C20C80FF-D",
                    "M5-41C2149CF1-B",
                ],
            ),
        ] {
            let subs = grid.sub_zones(id(&grid, zone), 3).unwrap();
            assert_eq!(texts(&grid, &subs), want, "{zone}");
        }
    }

    /// The seven cases of the odd parent at an even depth (`I3HSubZones.ec:24`, `:45`, `:47`,
    /// `:66`, `:69`, `:91` and `:93`): a hexagon, a northern and a southern edge hexagon, a
    /// northern and a southern pentagon, and the two polar pentagons.
    #[test]
    fn the_odd_parent_order_at_depth_2_matches_the_engine() {
        let grid = isea3h();
        for (zone, want) in [
            (
                "A6-0-C",
                &[
                    "B6-1-B", "B6-4-B", "B6-1-D", "B6-1-C", "B6-2-B", "B6-4-C", "B6-5-B", "B6-2-D",
                    "B6-8-B", "B6-5-D", "B6-5-C", "B8-2-B", "B8-1-B",
                ][..],
            ),
            (
                "B6-1-B",
                &[
                    "C4-35-B", "C4-3E-B", "C4-35-D", "C4-35-C", "C6-4-B", "C4-3E-C", "C6-3-B",
                    "C6-3-C", "C6-2-B", "C6-2-C", "C6-3-D", "C6-D-B", "C6-C-B",
                ],
            ),
            (
                "B5-3-B",
                &[
                    "C3-4D-B", "C5-24-B", "C3-4D-D", "C3-4D-C", "C3-4E-B", "C5-1B-D", "C5-1B-B",
                    "C3-4E-D", "C5-25-B", "C5-1B-C", "C5-12-D", "C5-12-B", "C5-1C-B",
                ],
            ),
            (
                "A6-0-B",
                &[
                    "B4-8-B", "B5-2-B", "B4-8-D", "B4-8-C", "B6-1-B", "B5-2-C", "B6-0-B", "B6-0-C",
                    "B6-3-B", "B6-0-D", "B6-4-B",
                ],
            ),
            (
                "A5-0-B",
                &[
                    "B3-8-B", "B5-3-B", "B3-8-D", "B3-8-C", "B4-6-B", "B5-0-D", "B5-0-B", "B4-6-D",
                    "B5-4-B", "B5-0-C", "B5-1-B",
                ],
            ),
            (
                "AA-0-B",
                &[
                    "B0-2-B", "B2-2-B", "B0-2-C", "B8-2-C", "B8-2-B", "B2-2-C", "BA-0-B", "B6-2-C",
                    "B4-2-B", "B4-2-C", "B6-2-B",
                ],
            ),
            (
                "AB-0-B",
                &[
                    "B1-6-B", "B9-6-B", "B9-6-D", "B1-6-D", "B3-6-B", "B7-6-D", "BB-0-B", "B3-6-D",
                    "B7-6-B", "B5-6-D", "B5-6-B",
                ],
            ),
        ] {
            let subs = grid.sub_zones(id(&grid, zone), 2).unwrap();
            assert_eq!(texts(&grid, &subs), want, "{zone}");
        }
    }

    /// Depth 4 of an odd parent reaches the three sections of the generator, the first cap,
    /// the main section and the second cap, with the columns a pentagon lacks and, for the
    /// polar pentagon, the crossings of the interruptions.
    #[test]
    fn the_odd_parent_order_at_depth_4_matches_the_engine() {
        let grid = isea3h();
        for (zone, want) in [
            (
                "A6-0-C",
                &[
                    "C6-3-B", "C6-C-B", "C6-3-D", "C6-3-C", "C6-4-B", "C6-15-B", "C6-C-D",
                    "C6-C-C", "C6-D-B", "C6-4-D", "C6-4-C", "C6-5-B", "C6-1E-B", "C6-15-D",
                    "C6-15-C", "C6-16-B", "C6-D-D", "C6-D-C", "C6-E-B", "C6-5-D", "C6-5-C",
                    "C6-6-B", "C6-1E-C", "C6-1F-B", "C6-16-D", "C6-16-C", "C6-17-B", "C6-E-D",
                    "C6-E-C", "C6-F-B", "C6-6-D", "C6-28-B", "C6-1F-D", "C6-1F-C", "C6-20-B",
                    "C6-17-D", "C6-17-C", "C6-18-B", "C6-F-D", "C6-F-C", "C6-10-B", "C6-28-C",
                    "C6-29-B", "C6-20-D", "C6-20-C", "C6-21-B", "C6-18-D", "C6-18-C", "C6-19-B",
                    "C6-10-D", "C6-32-B", "C6-29-D", "C6-29-C", "C6-2A-B", "C6-21-D", "C6-21-C",
                    "C6-22-B", "C6-19-D", "C6-19-C", "C6-1A-B", "C6-32-C", "C6-33-B", "C6-2A-D",
                    "C6-2A-C", "C6-2B-B", "C6-22-D", "C6-22-C", "C6-23-B", "C6-1A-D", "C6-3C-B",
                    "C6-33-D", "C6-33-C", "C6-34-B", "C6-2B-D", "C6-2B-C", "C6-2C-B", "C6-23-D",
                    "C6-23-C", "C8-6-B", "C6-3D-B", "C6-34-D", "C6-34-C", "C6-35-B", "C6-2C-D",
                    "C6-2C-C", "C8-5-B", "C6-3E-B", "C6-35-D", "C6-35-C", "C8-4-B", "C8-3-B",
                ][..],
            ),
            (
                "A5-0-B",
                &[
                    "C3-3C-B", "C3-45-B", "C3-3C-D", "C3-3C-C", "C3-3D-B", "C3-4E-B", "C3-45-D",
                    "C3-45-C", "C3-46-B", "C3-3D-D", "C3-3D-C", "C3-3E-B", "C5-1B-B", "C3-4E-D",
                    "C3-4E-C", "C3-4F-B", "C3-46-D", "C3-46-C", "C3-47-B", "C3-3E-D", "C3-3E-C",
                    "C4-36-B", "C5-12-D", "C5-12-B", "C3-4F-D", "C3-4F-C", "C3-50-B", "C3-47-D",
                    "C3-47-C", "C4-3F-B", "C4-36-D", "C5-1C-B", "C5-12-C", "C5-9-D", "C5-9-B",
                    "C3-50-D", "C3-50-C", "C4-48-B", "C4-3F-D", "C4-3F-C", "C4-40-B", "C5-13-D",
                    "C5-13-B", "C5-9-C", "C5-0-D", "C5-0-B", "C4-48-D", "C4-48-C", "C4-49-B",
                    "C4-40-D", "C5-1D-B", "C5-13-C", "C5-A-D", "C5-A-B", "C5-0-C", "C5-1-B",
                    "C4-49-D", "C4-49-C", "C4-4A-B", "C5-14-D", "C5-14-B", "C5-A-C", "C5-1-D",
                    "C5-1-C", "C5-2-B", "C4-4A-D", "C5-1E-B", "C5-14-C", "C5-B-D", "C5-B-B",
                    "C5-2-D", "C5-2-C", "C5-3-B", "C5-15-B", "C5-B-C", "C5-C-B",
                ],
            ),
            (
                "AB-0-B",
                &[
                    "C1-36-B", "C9-4A-B", "C9-4A-D", "C1-36-D", "C1-40-B", "C9-40-B", "C9-40-D",
                    "C9-49-C", "C1-3F-B", "C1-3F-C", "C1-40-D", "C1-4A-B", "C9-36-B", "C9-36-D",
                    "C9-3F-C", "C9-49-B", "C9-49-D", "C1-3F-D", "C1-49-B", "C1-49-C", "C1-4A-D",
                    "C3-36-B", "C7-4A-D", "C9-3F-B", "C9-3F-D", "C9-48-C", "C1-48-B", "C1-48-C",
                    "C1-49-D", "C3-3F-B", "C3-36-D", "C7-4A-B", "C7-49-C", "C7-49-D", "C9-48-B",
                    "C9-48-D", "C1-48-D", "C3-48-B", "C3-3F-D", "C3-3F-C", "C3-40-B", "C7-40-D",
                    "C7-49-B", "C7-48-C", "C7-48-D", "CB-0-B", "C3-48-D", "C3-48-C", "C3-49-B",
                    "C3-40-D", "C7-40-B", "C7-3F-C", "C7-3F-D", "C7-48-B", "C5-48-D", "C5-48-B",
                    "C3-49-D", "C3-49-C", "C3-4A-B", "C7-36-D", "C7-3F-B", "C5-49-D", "C5-48-C",
                    "C5-3F-D", "C5-3F-B", "C3-4A-D", "C7-36-B", "C5-4A-D", "C5-49-C", "C5-49-B",
                    "C5-3F-C", "C5-36-D", "C5-36-B", "C5-4A-B", "C5-40-D", "C5-40-B",
                ],
            ),
        ] {
            let subs = grid.sub_zones(id(&grid, zone), 4).unwrap();
            assert_eq!(texts(&grid, &subs), want, "{zone}");
        }
    }

    /// Depth 6 is the first depth at which the second cap of a pentagon has more than one row,
    /// each with the columns it lacks; its last twelve sub-zones, for the four kinds of
    /// pentagon, pin them.
    #[test]
    fn the_pentagons_odd_parent_order_at_depth_6_ends_as_the_engines() {
        let grid = isea3h();
        for (zone, want) in [
            (
                "A6-0-B",
                &[
                    "D6-DE-D", "D6-FA-B", "D6-F6-B", "D6-DB-D", "D6-DB-C", "D6-DC-B", "D6-DC-C",
                    "D6-DD-D", "D6-F9-B", "D6-F7-B", "D6-DC-D", "D6-F8-B",
                ][..],
            ),
            (
                "A5-0-B",
                &[
                    "D5-3E-C", "D5-3F-B", "D5-AB-B", "D5-8F-C", "D5-74-D", "D5-74-B", "D5-59-D",
                    "D5-59-C", "D5-5A-B", "D5-90-B", "D5-74-C", "D5-75-B",
                ],
            ),
            (
                "AA-0-B",
                &[
                    "D4-BA-C", "D4-D6-B", "D4-66-B", "D4-66-C", "D4-67-D", "D4-83-B", "D4-83-D",
                    "D4-9E-C", "D4-BA-B", "D4-82-B", "D4-82-C", "D4-9E-B",
                ],
            ),
            (
                "AB-0-B",
                &[
                    "D5-21E-D", "D5-21E-B", "D5-28E-B", "D5-272-D", "D5-271-C", "D5-271-B",
                    "D5-255-C", "D5-23A-D", "D5-23A-B", "D5-272-B", "D5-256-D", "D5-256-B",
                ],
            ),
        ] {
            let subs = grid.sub_zones(id(&grid, zone), 6).unwrap();
            assert_eq!(subs.len(), 631, "{zone}");
            assert_eq!(texts(&grid, &subs[619..]), want, "{zone}");
        }
    }

    /// The seven cases of the even parent at an even depth (`I3HSubZones.ec:20`, `:41`, `:43`,
    /// `:59`, `:64`, `:82` and `:89`): a hexagon, a northern and a southern edge hexagon, a
    /// northern and a southern pentagon, and the two polar pentagons.
    #[test]
    fn the_even_parent_order_at_depth_2_matches_the_engine() {
        let grid = isea3h();
        for (zone, want) in [
            (
                "B6-5-A",
                &[
                    "C6-19-A", "C6-E-A", "C6-18-A", "C6-22-A", "C6-2C-A", "C6-17-A", "C6-21-A",
                    "C6-2B-A", "C6-16-A", "C6-20-A", "C6-2A-A", "C6-34-A", "C6-29-A",
                ][..],
            ),
            (
                "B6-1-A",
                &[
                    "C4-2C-A", "C4-34-A", "C4-35-A", "C6-4-A", "C6-E-A", "C4-3E-A", "C6-3-A",
                    "C6-D-A", "C4-47-A", "C6-2-A", "C6-C-A", "C6-16-A", "C6-B-A",
                ],
            ),
            (
                "B5-3-A",
                &[
                    "C3-4F-A", "C3-44-A", "C3-4E-A", "C5-12-A", "C5-13-A", "C3-4D-A", "C5-1B-A",
                    "C5-1C-A", "C3-4C-A", "C5-24-A", "C5-25-A", "C5-26-A", "C5-2E-A",
                ],
            ),
            (
                "A2-0-A",
                &[
                    "B0-5-A", "B0-7-A", "B0-8-A", "B2-1-A", "B2-5-A", "B1-2-A", "B2-0-A", "B2-4-A",
                    "B1-5-A", "B2-3-A", "B2-7-A",
                ],
            ),
            (
                "A3-0-A",
                &[
                    "B2-7-A", "B1-5-A", "B2-6-A", "B3-1-A", "B3-5-A", "B1-8-A", "B3-0-A", "B3-4-A",
                    "B1-7-A", "B3-3-A", "B3-7-A",
                ],
            ),
            (
                "AA-0-A",
                &[
                    "B0-5-A", "B2-5-A", "B2-2-A", "B0-2-A", "B8-5-A", "B4-2-A", "BA-0-A", "B8-2-A",
                    "B4-5-A", "B6-2-A", "B6-5-A",
                ],
            ),
            (
                "AB-0-A",
                &[
                    "B9-7-A", "B7-7-A", "B9-6-A", "B1-6-A", "B1-7-A", "B7-6-A", "BB-0-A", "B3-6-A",
                    "B5-7-A", "B5-6-A", "B3-7-A",
                ],
            ),
        ] {
            let subs = grid.sub_zones(id(&grid, zone), 2).unwrap();
            assert_eq!(texts(&grid, &subs), want, "{zone}");
        }
    }

    /// Depth 4 of an even parent reaches the three sections of the generator, each cap of a
    /// hexagon with three rows: for the northern edge hexagon, the start of a scanline of its
    /// second cap that gcc regrouped (`I3HSubZones.ec:1223-1224`); for the southern edge
    /// hexagon, the scanlines that cross the interruption; for the southern pentagons, the
    /// columns they skip in the main portion; for the polar pentagons, the crossings of the
    /// interruptions.
    #[test]
    fn the_even_parent_order_at_depth_4_matches_the_engine() {
        let grid = isea3h();
        for (zone, want) in [
            (
                "B6-5-A",
                &[
                    "D6-B7-A", "D6-9A-A", "D6-B6-A", "D6-D2-A", "D6-EE-A", "D6-7D-A", "D6-99-A",
                    "D6-B5-A", "D6-D1-A", "D6-ED-A", "D6-109-A", "D6-125-A", "D6-60-A", "D6-7C-A",
                    "D6-98-A", "D6-B4-A", "D6-D0-A", "D6-EC-A", "D6-108-A", "D6-124-A", "D6-140-A",
                    "D6-15C-A", "D6-7B-A", "D6-97-A", "D6-B3-A", "D6-CF-A", "D6-EB-A", "D6-107-A",
                    "D6-123-A", "D6-13F-A", "D6-15B-A", "D6-7A-A", "D6-96-A", "D6-B2-A", "D6-CE-A",
                    "D6-EA-A", "D6-106-A", "D6-122-A", "D6-13E-A", "D6-15A-A", "D6-176-A",
                    "D6-95-A", "D6-B1-A", "D6-CD-A", "D6-E9-A", "D6-105-A", "D6-121-A", "D6-13D-A",
                    "D6-159-A", "D6-175-A", "D6-94-A", "D6-B0-A", "D6-CC-A", "D6-E8-A", "D6-104-A",
                    "D6-120-A", "D6-13C-A", "D6-158-A", "D6-174-A", "D6-190-A", "D6-AF-A",
                    "D6-CB-A", "D6-E7-A", "D6-103-A", "D6-11F-A", "D6-13B-A", "D6-157-A",
                    "D6-173-A", "D6-18F-A", "D6-AE-A", "D6-CA-A", "D6-E6-A", "D6-102-A",
                    "D6-11E-A", "D6-13A-A", "D6-156-A", "D6-172-A", "D6-18E-A", "D6-1AA-A",
                    "D6-E5-A", "D6-101-A", "D6-11D-A", "D6-139-A", "D6-155-A", "D6-171-A",
                    "D6-18D-A", "D6-11C-A", "D6-138-A", "D6-154-A", "D6-170-A", "D6-153-A",
                ][..],
            ),
            (
                "B6-1-A",
                &[
                    "D4-15C-A", "D4-176-A", "D4-177-A", "D4-178-A", "D4-179-A", "D4-190-A",
                    "D4-191-A", "D4-192-A", "D4-193-A", "D4-194-A", "D6-D-A", "D6-29-A",
                    "D4-1AA-A", "D4-1AB-A", "D4-1AC-A", "D4-1AD-A", "D4-1AE-A", "D4-1AF-A",
                    "D6-C-A", "D6-28-A", "D6-44-A", "D6-60-A", "D4-1C6-A", "D4-1C7-A", "D4-1C8-A",
                    "D4-1C9-A", "D4-1CA-A", "D6-B-A", "D6-27-A", "D6-43-A", "D6-5F-A", "D4-1E1-A",
                    "D4-1E2-A", "D4-1E3-A", "D4-1E4-A", "D4-1E5-A", "D6-A-A", "D6-26-A", "D6-42-A",
                    "D6-5E-A", "D6-7A-A", "D4-1FD-A", "D4-1FE-A", "D4-1FF-A", "D4-200-A", "D6-9-A",
                    "D6-25-A", "D6-41-A", "D6-5D-A", "D6-79-A", "D4-218-A", "D4-219-A", "D4-21A-A",
                    "D4-21B-A", "D6-8-A", "D6-24-A", "D6-40-A", "D6-5C-A", "D6-78-A", "D6-94-A",
                    "D4-234-A", "D4-235-A", "D4-236-A", "D6-7-A", "D6-23-A", "D6-3F-A", "D6-5B-A",
                    "D6-77-A", "D6-93-A", "D4-24F-A", "D4-250-A", "D4-251-A", "D6-6-A", "D6-22-A",
                    "D6-3E-A", "D6-5A-A", "D6-76-A", "D6-92-A", "D6-AE-A", "D4-26C-A", "D6-5-A",
                    "D6-21-A", "D6-3D-A", "D6-59-A", "D6-75-A", "D6-91-A", "D6-20-A", "D6-3C-A",
                    "D6-58-A", "D6-74-A", "D6-57-A",
                ],
            ),
            (
                "B1-3-A",
                &[
                    "D9-29D-A", "D9-280-A", "D9-29C-A", "D9-2B8-A", "D9-2D4-A", "D9-263-A",
                    "D9-27F-A", "D9-29B-A", "D9-2B7-A", "D9-2D3-A", "D1-87-A", "D1-88-A",
                    "D9-246-A", "D9-262-A", "D9-27E-A", "D9-29A-A", "D9-2B6-A", "D9-2D2-A",
                    "D1-A2-A", "D1-A3-A", "D1-A4-A", "D1-A5-A", "D9-261-A", "D9-27D-A", "D9-299-A",
                    "D9-2B5-A", "D9-2D1-A", "D1-BD-A", "D1-BE-A", "D1-BF-A", "D1-C0-A", "D9-260-A",
                    "D9-27C-A", "D9-298-A", "D9-2B4-A", "D9-2D0-A", "D1-D8-A", "D1-D9-A",
                    "D1-DA-A", "D1-DB-A", "D1-DC-A", "D9-27B-A", "D9-297-A", "D9-2B3-A",
                    "D9-2CF-A", "D1-F3-A", "D1-F4-A", "D1-F5-A", "D1-F6-A", "D1-F7-A", "D9-27A-A",
                    "D9-296-A", "D9-2B2-A", "D9-2CE-A", "D1-10E-A", "D1-10F-A", "D1-110-A",
                    "D1-111-A", "D1-112-A", "D1-113-A", "D9-295-A", "D9-2B1-A", "D9-2CD-A",
                    "D1-129-A", "D1-12A-A", "D1-12B-A", "D1-12C-A", "D1-12D-A", "D1-12E-A",
                    "D9-294-A", "D9-2B0-A", "D9-2CC-A", "D1-144-A", "D1-145-A", "D1-146-A",
                    "D1-147-A", "D1-148-A", "D1-149-A", "D1-14A-A", "D9-2CB-A", "D1-15F-A",
                    "D1-160-A", "D1-161-A", "D1-162-A", "D1-163-A", "D1-164-A", "D1-17B-A",
                    "D1-17C-A", "D1-17D-A", "D1-17E-A", "D1-198-A",
                ],
            ),
            (
                "A3-0-A",
                &[
                    "C2-39-A", "C2-2E-A", "C2-38-A", "C2-42-A", "C2-4C-A", "C1-2C-A", "C2-2D-A",
                    "C2-37-A", "C2-41-A", "C2-4B-A", "C3-4-A", "C3-E-A", "C1-21-A", "C1-2B-A",
                    "C1-35-A", "C2-36-A", "C2-40-A", "C2-4A-A", "C3-3-A", "C3-D-A", "C3-17-A",
                    "C3-21-A", "C1-2A-A", "C1-34-A", "C1-3E-A", "C2-3F-A", "C2-49-A", "C3-2-A",
                    "C3-C-A", "C3-16-A", "C3-20-A", "C1-29-A", "C1-33-A", "C1-3D-A", "C1-47-A",
                    "C2-48-A", "C3-1-A", "C3-B-A", "C3-15-A", "C3-1F-A", "C3-29-A", "C1-32-A",
                    "C1-3C-A", "C1-46-A", "C1-50-A", "C3-0-A", "C3-A-A", "C3-14-A", "C3-1E-A",
                    "C3-28-A", "C1-31-A", "C1-3B-A", "C1-45-A", "C1-4F-A", "C3-9-A", "C3-13-A",
                    "C3-1D-A", "C3-27-A", "C3-31-A", "C1-3A-A", "C1-44-A", "C1-4E-A", "C3-12-A",
                    "C3-1C-A", "C3-26-A", "C3-30-A", "C1-39-A", "C1-43-A", "C1-4D-A", "C3-1B-A",
                    "C3-25-A", "C3-2F-A", "C3-39-A", "C1-4C-A", "C3-24-A", "C3-2E-A",
                ],
            ),
            (
                "A1-0-A",
                &[
                    "C0-39-A", "C0-2E-A", "C0-38-A", "C0-42-A", "C0-4C-A", "C9-2C-A", "C0-2D-A",
                    "C0-37-A", "C0-41-A", "C0-4B-A", "C1-4-A", "C1-E-A", "C9-21-A", "C9-2B-A",
                    "C9-35-A", "C0-36-A", "C0-40-A", "C0-4A-A", "C1-3-A", "C1-D-A", "C1-17-A",
                    "C1-21-A", "C9-2A-A", "C9-34-A", "C9-3E-A", "C0-3F-A", "C0-49-A", "C1-2-A",
                    "C1-C-A", "C1-16-A", "C1-20-A", "C9-29-A", "C9-33-A", "C9-3D-A", "C9-47-A",
                    "C0-48-A", "C1-1-A", "C1-B-A", "C1-15-A", "C1-1F-A", "C1-29-A", "C9-32-A",
                    "C9-3C-A", "C9-46-A", "C9-50-A", "C1-0-A", "C1-A-A", "C1-14-A", "C1-1E-A",
                    "C1-28-A", "C9-31-A", "C9-3B-A", "C9-45-A", "C9-4F-A", "C1-9-A", "C1-13-A",
                    "C1-1D-A", "C1-27-A", "C1-31-A", "C9-3A-A", "C9-44-A", "C9-4E-A", "C1-12-A",
                    "C1-1C-A", "C1-26-A", "C1-30-A", "C9-39-A", "C9-43-A", "C9-4D-A", "C1-1B-A",
                    "C1-25-A", "C1-2F-A", "C1-39-A", "C9-4C-A", "C1-24-A", "C1-2E-A",
                ],
            ),
            (
                "AA-0-A",
                &[
                    "C0-21-A", "C0-2C-A", "C0-22-A", "C0-18-A", "C0-E-A", "C2-E-A", "C2-5-A",
                    "C0-23-A", "C0-19-A", "C0-F-A", "C0-5-A", "C8-2C-A", "C2-21-A", "C2-18-A",
                    "C2-F-A", "C2-6-A", "C0-1A-A", "C0-10-A", "C0-6-A", "C8-23-A", "C8-22-A",
                    "C8-21-A", "C2-22-A", "C2-19-A", "C2-10-A", "C2-7-A", "C0-11-A", "C0-7-A",
                    "C8-1A-A", "C8-19-A", "C8-18-A", "C2-2C-A", "C2-23-A", "C2-1A-A", "C2-11-A",
                    "C2-8-A", "C0-8-A", "C8-11-A", "C8-10-A", "C8-F-A", "C8-E-A", "C4-5-A",
                    "C4-6-A", "C4-7-A", "C4-8-A", "CA-0-A", "C8-8-A", "C8-7-A", "C8-6-A", "C8-5-A",
                    "C4-E-A", "C4-F-A", "C4-10-A", "C4-11-A", "C6-8-A", "C6-11-A", "C6-1A-A",
                    "C6-23-A", "C6-2C-A", "C4-18-A", "C4-19-A", "C4-1A-A", "C6-7-A", "C6-10-A",
                    "C6-19-A", "C6-22-A", "C4-21-A", "C4-22-A", "C4-23-A", "C6-6-A", "C6-F-A",
                    "C6-18-A", "C6-21-A", "C4-2C-A", "C6-5-A", "C6-E-A",
                ],
            ),
            (
                "AB-0-A",
                &[
                    "C9-39-A", "C9-2E-A", "C9-38-A", "C9-42-A", "C9-4C-A", "C7-4C-A", "C9-2D-A",
                    "C9-37-A", "C9-41-A", "C9-4B-A", "C1-2D-A", "C1-2E-A", "C7-39-A", "C7-42-A",
                    "C7-4B-A", "C9-36-A", "C9-40-A", "C9-4A-A", "C1-36-A", "C1-37-A", "C1-38-A",
                    "C1-39-A", "C7-38-A", "C7-41-A", "C7-4A-A", "C9-3F-A", "C9-49-A", "C1-3F-A",
                    "C1-40-A", "C1-41-A", "C1-42-A", "C7-2E-A", "C7-37-A", "C7-40-A", "C7-49-A",
                    "C9-48-A", "C1-48-A", "C1-49-A", "C1-4A-A", "C1-4B-A", "C1-4C-A", "C7-2D-A",
                    "C7-36-A", "C7-3F-A", "C7-48-A", "CB-0-A", "C3-48-A", "C3-3F-A", "C3-36-A",
                    "C3-2D-A", "C5-4C-A", "C5-4B-A", "C5-4A-A", "C5-49-A", "C5-48-A", "C3-49-A",
                    "C3-40-A", "C3-37-A", "C3-2E-A", "C5-42-A", "C5-41-A", "C5-40-A", "C5-3F-A",
                    "C3-4A-A", "C3-41-A", "C3-38-A", "C5-39-A", "C5-38-A", "C5-37-A", "C5-36-A",
                    "C3-4B-A", "C3-42-A", "C3-39-A", "C5-2E-A", "C5-2D-A", "C3-4C-A",
                ],
            ),
            // The southern polar pentagon at level 16, where calling `move5x6Vertex3` rather
            // than `move5x6Vertex` at I3HSubZones.ec:1309 decides its last sub-zone.
            (
                "IB-0-A",
                &[
                    "K9-CFD16799-A",
                    "K9-CFD080EE-A",
                    "K9-CFD16798-A",
                    "K9-CFD24E42-A",
                    "K9-CFD334EC-A",
                    "K7-CFD334EC-A",
                    "K9-CFD080ED-A",
                    "K9-CFD16797-A",
                    "K9-CFD24E41-A",
                    "K9-CFD334EB-A",
                    "K1-CFD080ED-A",
                    "K1-CFD080EE-A",
                    "K7-CFD16799-A",
                    "K7-CFD24E42-A",
                    "K7-CFD334EB-A",
                    "K9-CFD16796-A",
                    "K9-CFD24E40-A",
                    "K9-CFD334EA-A",
                    "K1-CFD16796-A",
                    "K1-CFD16797-A",
                    "K1-CFD16798-A",
                    "K1-CFD16799-A",
                    "K7-CFD16798-A",
                    "K7-CFD24E41-A",
                    "K7-CFD334EA-A",
                    "K9-CFD24E3F-A",
                    "K9-CFD334E9-A",
                    "K1-CFD24E3F-A",
                    "K1-CFD24E40-A",
                    "K1-CFD24E41-A",
                    "K1-CFD24E42-A",
                    "K7-CFD080EE-A",
                    "K7-CFD16797-A",
                    "K7-CFD24E40-A",
                    "K7-CFD334E9-A",
                    "K9-CFD334E8-A",
                    "K1-CFD334E8-A",
                    "K1-CFD334E9-A",
                    "K1-CFD334EA-A",
                    "K1-CFD334EB-A",
                    "K1-CFD334EC-A",
                    "K7-CFD080ED-A",
                    "K7-CFD16796-A",
                    "K7-CFD24E3F-A",
                    "K7-CFD334E8-A",
                    "KB-0-A",
                    "K3-CFD334E8-A",
                    "K3-CFD24E3F-A",
                    "K3-CFD16796-A",
                    "K3-CFD080ED-A",
                    "K5-CFD334EC-A",
                    "K5-CFD334EB-A",
                    "K5-CFD334EA-A",
                    "K5-CFD334E9-A",
                    "K5-CFD334E8-A",
                    "K3-CFD334E9-A",
                    "K3-CFD24E40-A",
                    "K3-CFD16797-A",
                    "K3-CFD080EE-A",
                    "K5-CFD24E42-A",
                    "K5-CFD24E41-A",
                    "K5-CFD24E40-A",
                    "K5-CFD24E3F-A",
                    "K3-CFD334EA-A",
                    "K3-CFD24E41-A",
                    "K3-CFD16798-A",
                    "K5-CFD16799-A",
                    "K5-CFD16798-A",
                    "K5-CFD16797-A",
                    "K5-CFD16796-A",
                    "K3-CFD334EB-A",
                    "K3-CFD24E42-A",
                    "K3-CFD16799-A",
                    "K5-CFD080EE-A",
                    "K5-CFD080ED-A",
                    "K3-CFD334EC-A",
                ],
            ),
        ] {
            let subs = grid.sub_zones(id(&grid, zone), 4).unwrap();
            assert_eq!(texts(&grid, &subs), want, "{zone}");
        }
    }

    /// Depth 6 is the first depth at which the second cap of an even pentagon has more than one
    /// row, each with the columns it lacks; its last twelve sub-zones, for the four kinds of
    /// pentagon, pin them.
    #[test]
    fn the_pentagons_even_parent_order_at_depth_6_ends_as_the_engines() {
        let grid = isea3h();
        for (zone, want) in [
            (
                "A2-0-A",
                &[
                    "D2-199-A", "D2-1B5-A", "D1-15C-A", "D1-15D-A", "D1-15E-A", "D2-144-A",
                    "D2-160-A", "D2-17C-A", "D2-198-A", "D1-179-A", "D2-15F-A", "D2-17B-A",
                ][..],
            ),
            (
                "A3-0-A",
                &[
                    "D3-199-A", "D3-1B5-A", "D1-294-A", "D1-2B0-A", "D1-2CC-A", "D3-144-A",
                    "D3-160-A", "D3-17C-A", "D3-198-A", "D1-2CB-A", "D3-15F-A", "D3-17B-A",
                ],
            ),
            (
                "AA-0-A",
                &[
                    "D6-7C-A", "D6-97-A", "D4-15C-A", "D4-15D-A", "D4-15E-A", "D6-F-A", "D6-2A-A",
                    "D6-45-A", "D6-60-A", "D4-179-A", "D6-E-A", "D6-29-A",
                ],
            ),
            (
                "AB-0-A",
                &[
                    "D3-278-A", "D3-25D-A", "D5-198-A", "D5-197-A", "D5-196-A", "D5-195-A",
                    "D3-2CA-A", "D3-2AF-A", "D3-294-A", "D5-17B-A", "D5-17A-A", "D3-2CB-A",
                ],
            ),
        ] {
            let subs = grid.sub_zones(id(&grid, zone), 6).unwrap();
            assert_eq!(subs.len(), 631, "{zone}");
            assert_eq!(texts(&grid, &subs[619..]), want, "{zone}");
        }
    }

    /// The count, the first sub-zone and the order answer nothing where the sub-zones would lie
    /// beyond level 33, whatever the parities, and answer at level 33 itself.
    #[test]
    fn nothing_answers_beyond_level_33() {
        use crate::indexings::pack;
        // A zone at level 33 (odd), one at level 32 (even) and one at level 0.
        for (z, deepest) in [
            (pack(16, 6, 5, 1), 0),
            (pack(16, 6, 5, 0), 1),
            (pack(0, 6, 0, 0), 33),
        ] {
            assert!(count(z, deepest).is_some(), "{z:#018x} at {deepest}");
            assert!(
                first_sub_zone(z, deepest).is_some(),
                "{z:#018x} at {deepest}"
            );
            for depth in [deepest + 1, deepest + 2, u8::MAX] {
                assert_eq!(count(z, depth), None, "{z:#018x} at {depth}");
                assert_eq!(first_sub_zone(z, depth), None, "{z:#018x} at {depth}");
                assert_eq!(sub_zones(z, depth), None, "{z:#018x} at {depth}");
                assert_eq!(sub_zone_at_index(z, depth, 0), None, "{z:#018x} at {depth}");
            }
        }
    }

    /// The first sub-zone in all four combinations of parent and depth parity, and for every
    /// kind of zone, through the branches of `getI3HFirstSubZoneCentroid`
    /// (`I3HSubZones.ec:102-247`).
    #[test]
    fn first_sub_zones_match_the_engine() {
        let grid = isea3h();
        for (zone, depth, want) in [
            // Even parent, even depth.
            ("B6-5-A", 2, "C6-19-A"),
            ("A2-0-A", 2, "B0-5-A"),
            ("A3-0-A", 2, "B2-7-A"),
            ("AA-0-A", 2, "B0-5-A"),
            ("AB-0-A", 2, "B9-7-A"),
            ("B6-1-A", 2, "C4-2C-A"),
            ("B5-3-A", 2, "C3-4F-A"),
            // Even parent, odd depth.
            ("B6-5-A", 3, "C6-16-B"),
            ("A2-0-A", 3, "B0-7-B"),
            ("A3-0-A", 3, "B1-7-B"),
            ("AA-0-A", 3, "B0-5-B"),
            ("AB-0-A", 3, "B9-7-B"),
            // Odd parent, even depth.
            ("A6-0-C", 2, "B6-1-B"),
            ("A6-0-B", 2, "B4-8-B"),
            ("A5-0-B", 2, "B3-8-B"),
            ("AA-0-B", 2, "B0-2-B"),
            ("AB-0-B", 2, "B1-6-B"),
            ("B6-1-B", 2, "C4-35-B"),
            ("B5-3-B", 2, "C3-4D-B"),
            // Odd parent, odd depth.
            ("B6-4-C", 3, "D6-102-A"),
            ("A5-0-B", 3, "C4-36-A"),
            ("AB-0-B", 3, "C9-36-A"),
        ] {
            let first = grid.first_sub_zone(id(&grid, zone), depth).unwrap();
            assert_eq!(grid.zone(first).text_id(), want, "{zone} at {depth}");
        }
        // At depth 0 the engine's `getFirstSubZone` names `B4-5-D`, a neighbour, where its
        // `getSubZones` names the zone itself; the crate answers the zone.
        let z = id(&grid, "B6-1-B");
        assert_eq!(grid.first_sub_zone(z, 0), Ok(z));
    }

    /// Checks that the index skip finds each of the sub-zones of `z` at `depth` at its index in
    /// the order, and nothing at the count; that the walk hands the callback the index it was
    /// asked for and stops there; and that a walk the callback never stops hands it every index
    /// below the count once, in order, and answers `-1`, as the eC's generators account for
    /// them.
    fn check_index_skip(z: u64, depth: u8) {
        let subs = sub_zones(z, depth).unwrap();
        assert_eq!(
            Some(subs.len() as u64),
            count(z, depth),
            "{z:#018x} at {depth}"
        );
        let r_depth = u32::from(depth);
        let mut delivered = vec![];
        let stop = iterate(
            z,
            r_depth,
            &mut |i, _| {
                delivered.push(i);
                true
            },
            -1,
        );
        assert_eq!(stop, -1, "{z:#018x} at {depth}");
        assert!(
            delivered.iter().copied().eq(0..subs.len() as i64),
            "{z:#018x} at {depth}: the indices delivered"
        );
        for (i, &sub) in subs.iter().enumerate() {
            let at = sub_zone_at_index(z, depth, i as u64);
            assert_eq!(at, Some(sub), "{z:#018x} at {depth}, index {i}");
            let mut handed = None;
            let stop = iterate(
                z,
                r_depth,
                &mut |j, _| {
                    handed = Some(j);
                    false
                },
                i as i64,
            );
            assert_eq!((stop, handed), (i as i64, Some(i as i64)));
        }
        assert_eq!(sub_zone_at_index(z, depth, subs.len() as u64), None);
    }

    /// The index skip over the seven cases of each of the four generators: at levels 0 to 5, the
    /// polar pentagons, every pentagon, the edge hexagons at both ends of the edge of every
    /// rhombus, and the zones A, B, C and D at the first, the middle and the last index of every
    /// rhombus, at depths 1 to 5; and one zone of each case at depths 6 and 7.
    #[test]
    fn the_index_skip_finds_every_sub_zone() {
        use crate::indexings::pack;
        for level_i9r in 0..3 {
            let p = 3u64.pow(level_i9r as u32);
            let mut zones = vec![];
            for sub_hex in [0, 1] {
                zones.extend([10, 11].map(|root| pack(level_i9r, root, 0, sub_hex)));
            }
            for root in 0..10 {
                for rix in [0, p * p / 2, p * p - 1] {
                    zones.extend((0..=3).map(|sub_hex| pack(level_i9r, root, rix, sub_hex)));
                }
                if p > 1 {
                    let edge = if root & 1 == 0 {
                        [1, p - 1]
                    } else {
                        [p, p * (p - 1)]
                    };
                    for sub_hex in [0, 1] {
                        zones.extend(edge.map(|rix| pack(level_i9r, root, rix, sub_hex)));
                    }
                }
            }
            for z in zones {
                for depth in 1..=5 {
                    check_index_skip(z, depth);
                }
            }
        }
        let grid = isea3h();
        for zone in [
            "B6-4-C", "B6-1-B", "B5-3-B", "B6-0-B", "B5-0-B", "BA-0-B", "BB-0-B", "B6-5-A",
            "B6-1-A", "B5-3-A", "B6-0-A", "B5-0-A", "BA-0-A", "BB-0-A",
        ] {
            check_index_skip(id(&grid, zone).0, 6);
            check_index_skip(id(&grid, zone).0, 7);
        }
    }

    /// The engine's own test of whether `needle` lies under `hay`, its `zoneHasSubZone`
    /// (`dggrs.ec:151-189`), ported for these tests: some vertex of the needle's centroid child,
    /// canonicalised and quantised at the level of `hay`, is `hay` itself. The engine's test
    /// fails for a needle at level 33, whose centroid child would lie at 34, and is not asked
    /// of one here.
    fn has_sub_zone(hay: u64, needle: u64) -> bool {
        use super::super::hex_a3::{centroid_child, vertices};
        use crate::fivebysix::canonicalize5x6;
        assert!(
            level(needle) < u32::from(MAX_LEVEL),
            "{needle:#018x} lies at level 33, where the descendant test cannot be asked"
        );
        level(needle) > level(hay)
            && vertices(centroid_child(needle))
                .into_iter()
                .any(|(x, y)| quantize(level(hay), canonicalize5x6(x, y)) == hay)
    }

    /// What is wrong with the sub-zones of `z` at `depth`, if anything, as the zone's
    /// descendants and as paged: a zone named twice; a zone that [`has_sub_zone`] rejects; a
    /// number other than the count; at depth 1, a set other than the children; a first sub-zone
    /// other than the first entry; the index skip giving another zone at some place; and the
    /// index search giving another place for some entry. The index search builds the order at
    /// each call, so beyond depth 5 it is asked of some seventeen places spread over it.
    fn faults(grid: &Grid<Isea, HexA3, I3h>, z: u64, depth: u8) -> Vec<String> {
        let mut faults = vec![];
        let subs = sub_zones(z, depth).unwrap();
        let distinct: std::collections::BTreeSet<u64> = subs.iter().copied().collect();
        if distinct.len() != subs.len() {
            faults.push(format!("{} repeated", subs.len() - distinct.len()));
        }
        let strays = subs.iter().filter(|&&s| !has_sub_zone(z, s)).count();
        if strays != 0 {
            faults.push(format!("{strays} not descendants"));
        }
        if Some(subs.len() as u64) != count(z, depth) {
            faults.push(format!("{} of {:?}", subs.len(), count(z, depth)));
        }
        if depth == 1 {
            let children: std::collections::BTreeSet<u64> =
                grid.children(ZoneId(z)).iter().map(|c| c.0).collect();
            if children != distinct {
                faults.push("not the children".into());
            }
        }
        if first_sub_zone(z, depth) != subs.first().copied() {
            faults.push("the first sub-zone is not the first entry".into());
        }
        if let Some(i) =
            (0..subs.len()).find(|&i| sub_zone_at_index(z, depth, i as u64) != Some(subs[i]))
        {
            faults.push(format!("the index skip differs at {i}"));
        }
        let step = if depth <= 5 {
            1
        } else {
            subs.len().div_ceil(16)
        };
        let mut places = (0..subs.len()).step_by(step).chain([subs.len() - 1]);
        if let Some(i) =
            places.find(|&i| grid.sub_zone_index(ZoneId(z), ZoneId(subs[i])) != Ok(Some(i as u64)))
        {
            faults.push(format!("the index search differs at {i}"));
        }
        faults
    }

    /// Zones on the edges of the rhombi at which the engine's order is not the zone's
    /// descendants, by each of the three ways in which a unit in the last place carries a point
    /// across an interruption of the layout, at depths 1 to 3 and at deeper depths at which each
    /// way occurs: here each order is the zone's descendants, each once, as many as the count,
    /// the children at depth 1, and paged consistently (see [`faults`]).
    ///
    /// - A vertex of the zone carried across the interruption it lies one unit short of, the
    ///   first centroid then read from it: `K5-67E99A74-D` and `K5-89336A36-D` (sub-hex D on
    ///   the left edge of rhombus 5 at level 21), `G2-5-B` and `G2-5-C` (sub-hexes B and C on
    ///   the top edge of rhombus 2 at level 13).
    /// - A point stepped a unit beyond an interruption and carried across it twice, by the
    ///   generator of an even-level zone at an even depth: the hexagons of level 16 and 20 on
    ///   the left edge of rhombus 5, `K5-8855C3D0-A` among them.
    /// - The same, by the generator of an odd-level zone at an odd depth, at the southern polar
    ///   pentagon of levels 11 and 25, `FB-0-B` and `MB-0-B`.
    #[test]
    fn the_orders_at_the_rhombus_edges_are_the_zones_descendants() {
        use crate::indexings::pack;
        let grid = isea3h();
        let p8 = pow3(8);
        let witnesses: [(u64, &[u8]); 10] = [
            (id(&grid, "K5-67E99A74-D").0, &[1, 2, 3, 4, 7]),
            (id(&grid, "K5-89336A36-D").0, &[1, 2, 3, 4, 6]),
            (id(&grid, "G2-5-B").0, &[1, 2, 3, 5]),
            (id(&grid, "G2-5-C").0, &[1, 2, 3, 8]),
            (pack(8, 5, p8, 0), &[1, 2, 3, 4, 6, 8]),
            (pack(8, 5, p8 * (p8 / 2), 0), &[2, 4]),
            (id(&grid, "K5-8855C3D0-A").0, &[1, 2, 3, 4, 6]),
            (id(&grid, "K5-E6A9-A").0, &[2, 8]),
            (id(&grid, "FB-0-B").0, &[1, 2, 3, 5, 7]),
            (id(&grid, "MB-0-B").0, &[1, 2, 3, 5]),
        ];
        let mut failed = vec![];
        let mut pairs = 0;
        for (z, depths) in witnesses {
            for &depth in depths {
                let found = faults(&grid, z, depth);
                if !found.is_empty() {
                    let text = grid.zone(ZoneId(z)).text_id();
                    failed.push(format!("{text} at {depth}: {}", found.join(", ")));
                }
                pairs += 1;
            }
        }
        assert!(failed.is_empty(), "{}", failed.join("\n"));
        assert_eq!(pairs, 42);
    }

    /// The corrected orders keep the generator's own scanline order: where the engine's order
    /// goes wrong at a single place, only that place changes, and where its first centroid is
    /// read from a vertex carried across an interruption, the whole order moves to the one the
    /// generator gives from the vertex in the zone's own frame.
    #[test]
    fn the_corrected_orders_keep_the_generators_scanline_order() {
        let grid = isea3h();
        for (zone, depth, want) in [
            (
                "K5-67E99A74-D",
                1,
                &[
                    "L5-3A7392210-A",
                    "L5-3A73BD60C-A",
                    "L5-3A739220F-A",
                    "L5-3A73BD60B-A",
                    "L5-3A73E8A07-A",
                    "L5-3A73BD60A-A",
                    "L5-3A73E8A06-A",
                ][..],
            ),
            // The engine names `L5-3-A` at place 4.
            (
                "K5-E6A9-A",
                2,
                &[
                    "L3-74E74F817-A",
                    "L3-74E72441A-A",
                    "L3-74E74F816-A",
                    "L5-567F6-A",
                    "L5-567F7-A",
                    "L3-74E74F815-A",
                    "L5-81BF1-A",
                    "L5-81BF2-A",
                    "L3-74E74F814-A",
                    "L5-ACFEC-A",
                    "L5-ACFED-A",
                    "L5-ACFEE-A",
                    "L5-D83E8-A",
                ],
            ),
            // The engine names `H9-1115-A`, `H0-1116-A` and `H8-3-A` at places 6 to 8.
            (
                "FB-0-B",
                3,
                &[
                    "H9-48E1D8-A",
                    "H9-48EA64-A",
                    "H9-48F2F0-A",
                    "H1-48E1D8-A",
                    "H7-48F2F0-A",
                    "H9-48EA63-A",
                    "H9-48F2EF-A",
                    "H1-48EA63-A",
                    "H1-48EA64-A",
                    "H7-48EA64-A",
                    "H7-48F2EF-A",
                    "H9-48F2EE-A",
                    "H1-48F2EE-A",
                    "H1-48F2EF-A",
                    "H1-48F2F0-A",
                    "H7-48E1D8-A",
                    "H7-48EA63-A",
                    "H7-48F2EE-A",
                    "HB-0-A",
                    "H3-48F2EE-A",
                    "H3-48EA63-A",
                    "H3-48E1D8-A",
                    "H5-48F2F0-A",
                    "H5-48F2EF-A",
                    "H5-48F2EE-A",
                    "H3-48F2EF-A",
                    "H3-48EA64-A",
                    "H5-48EA64-A",
                    "H5-48EA63-A",
                    "H3-48F2F0-A",
                    "H5-48E1D8-A",
                ],
            ),
        ] {
            let subs = grid.sub_zones(id(&grid, zone), depth).unwrap();
            assert_eq!(texts(&grid, &subs), want, "{zone} at {depth}");
        }
    }

    /// The orders the engine gets right are kept, bit for bit: at every zone of the index-skip
    /// test's sample (every case of the four generators, in both hemispheres, the pentagons and
    /// the polar pentagons among them), and at zones beside those at which the engine's order
    /// goes wrong (the other sub-hexes and rhombi at the same place, the depths at which it
    /// does not, the southern pentagons of sub-hex B, which carry a vertex across an
    /// interruption and whose first vertex the generator takes by its index, and the edges of
    /// levels 12 and 18, whose corners fall a unit short of a whole number with no vertex
    /// carried), at depths 1 to 4, each order is the zone's descendants (see [`faults`]), and
    /// its raw centroids, its sub-zones and its first sub-zone hash to the value pinned, which
    /// the engine's own order gives.
    #[test]
    fn the_orders_the_engine_gets_right_are_kept_bit_for_bit() {
        use crate::indexings::pack;
        let grid = isea3h();
        let mut zones = vec![];
        for level_i9r in 0..3 {
            let p = 3u64.pow(level_i9r as u32);
            for sub_hex in [0, 1] {
                zones.extend([10, 11].map(|root| pack(level_i9r, root, 0, sub_hex)));
            }
            for root in 0..10 {
                for rix in [0, p * p / 2, p * p - 1] {
                    zones.extend((0..=3).map(|sub_hex| pack(level_i9r, root, rix, sub_hex)));
                }
                if p > 1 {
                    let edge = if root & 1 == 0 {
                        [1, p - 1]
                    } else {
                        [p, p * (p - 1)]
                    };
                    for sub_hex in [0, 1] {
                        zones.extend(edge.map(|rix| pack(level_i9r, root, rix, sub_hex)));
                    }
                }
            }
        }
        let p6 = pow3(6);
        let p8 = pow3(8);
        let p9 = pow3(9);
        zones.extend(
            [
                "K5-67E99A74-B",
                "K5-67E99A74-C",
                "K4-67E99A74-D",
                "K7-67E99A74-D",
                "G2-5-D",
                "G0-5-B",
                "G6-5-C",
                "G3-0-B",
                "G5-0-B",
                "G9-0-B",
                "J3-0-B",
                "K5-0-B",
                "K9-0-B",
                "FA-0-B",
                "MA-0-B",
                "K7-E6A9-A",
                "K1-E6A9-A",
            ]
            .map(|t| id(&grid, t).0),
        );
        zones.extend([
            // Level 12, left edge of rhombus 5 and top edge of rhombus 2.
            pack(6, 5, p6 * (p6 / 2), 0),
            pack(6, 2, p6 / 2, 0),
            // Level 16, left edges of rhombi 1 and 9.
            pack(8, 1, p8 * (p8 / 2), 0),
            pack(8, 9, p8 * (p8 / 2), 0),
            // Level 18, left edge of rhombus 5.
            pack(9, 5, p9 * (p9 / 2), 0),
        ]);
        // Zones at which the engine's order goes wrong at some depths, at the others: the
        // edge hexagons of rhombus 5 at levels 16 and 20 at odd depths, sub-hex B on the top
        // edge of rhombus 2 at level 13 at even ones, and the southern polar pentagon of
        // levels 11 and 25 at depths 1, 2 and 4.
        let (b, c) = (&[1, 3][..], &[1, 2, 4][..]);
        let some_depths = [
            (pack(8, 5, p8, 0), b),
            (id(&grid, "K5-8855C3D0-A").0, b),
            (id(&grid, "G2-5-B").0, &[2, 4]),
            (id(&grid, "FB-0-B").0, c),
            (id(&grid, "MB-0-B").0, c),
        ];
        let mut digest: u64 = 0xcbf2_9ce4_8422_2325;
        let mut eat = |w: u64| {
            for b in w.to_le_bytes() {
                digest = (digest ^ u64::from(b)).wrapping_mul(0x0100_0000_01b3);
            }
        };
        let (mut pairs, mut failed) = (0, vec![]);
        let pairs_of = zones
            .iter()
            .flat_map(|&z| (1..=4).map(move |d| (z, d)))
            .chain(
                some_depths
                    .iter()
                    .flat_map(|&(z, ds)| ds.iter().map(move |&d| (z, d))),
            );
        for (z, depth) in pairs_of {
            let found = faults(&grid, z, depth);
            if !found.is_empty() {
                let text = grid.zone(ZoneId(z)).text_id();
                failed.push(format!("{text} at {depth}: {}", found.join(", ")));
            }
            eat(z);
            eat(u64::from(depth));
            iterate(
                z,
                u32::from(depth),
                &mut |i, (x, y)| {
                    eat(i as u64);
                    eat(x.to_bits());
                    eat(y.to_bits());
                    true
                },
                -1,
            );
            sub_zones(z, depth).unwrap().into_iter().for_each(&mut eat);
            eat(first_sub_zone(z, depth).unwrap());
            pairs += 1;
        }
        assert!(failed.is_empty(), "{}", failed.join("\n"));
        assert_eq!((pairs, digest), (1_908, 0x9676_9d7d_94c6_f06d));
    }

    /// The case of `I3HSubZones.ec:11-98` that the sub-zones of zone `id` at relative depth
    /// `depth` fall in, as `iterateI3HSubZones` tells them apart (`I3HSubZones.ec:1807-1817`).
    #[cfg(feature = "oracle")]
    fn case(id: u64, depth: u32) -> u32 {
        let (level_i9r, root, rhombus_ix, sub_hex) = fields(id);
        let south = u32::from(root & 1 != 0);
        let divs = pow3(level_i9r);
        let generator = match (sub_hex > 0, depth & 1 != 0) {
            (true, true) => 0,
            (false, true) => 1,
            (false, false) => 2,
            (true, false) => 3,
        };
        let edge_hex = n_points(id) == 6
            && rhombus_ix != 0
            && sub_hex <= 1
            && if south == 1 {
                rhombus_ix % divs == 0
            } else {
                rhombus_ix / divs == 0
            };
        if root > 9 {
            21 + 2 * generator + south
        } else if n_points(id) == 5 {
            13 + 2 * generator + south
        } else if edge_hex {
            5 + 2 * generator + south
        } else {
            1 + generator
        }
    }

    /// Identifiers of every I9R level the five-bit field holds and every root the four-bit field
    /// holds, with indices at and beyond the edges of a rhombus and every sub-hex, the null zone
    /// and 2,000 random values, at every depth to 34 and at 255, and at indices about the count
    /// and at the end of `u64`: the vertices the order is generated from, the count, the first
    /// sub-zone, the order where it is short and the index skip answer without a panic, called
    /// directly with no check of the identifier, and so do the same questions through `Grid`.
    /// Run in release as well as in debug, since only debug builds check integer overflow. The
    /// index skip walks the scanlines before the index, some `3^(depth / 2)`, so the indices far
    /// from the first are asked to depth 12.
    #[test]
    fn the_sub_zones_of_any_identifier_do_not_panic() {
        use super::super::hex_a3::sub_zone_vertices;
        use crate::indexings::pack;
        let grid = isea3h();
        let mut ids = vec![NULL_ZONE];
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
        let mut state: u64 = 0x5B20_0003;
        for _ in 0..2_000 {
            state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = state;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            ids.push(z ^ (z >> 31));
        }
        let (mut answered, mut carried) = (0, 0);
        for &z in &ids {
            let from = sub_zone_vertices(z);
            assert!(from.len() <= 6, "{z:#018x}");
            carried += usize::from(from != super::super::hex_a3::vertices(z));
            for depth in (0..=34).chain([u8::MAX]) {
                let Some(n) = count(z, depth) else {
                    assert_eq!(first_sub_zone(z, depth), None, "{z:#018x} at {depth}");
                    assert_eq!(sub_zone_at_index(z, depth, 0), None, "{z:#018x} at {depth}");
                    continue;
                };
                answered += 1;
                let _ = first_sub_zone(z, depth);
                if n <= 400 {
                    assert_eq!(sub_zones(z, depth).map(|s| s.len() as u64), Some(n));
                }
                let mut indices = vec![0, n, n.saturating_add(1), u64::MAX];
                if depth <= 12 {
                    indices.extend([1, n / 2, n - 1]);
                }
                for index in indices {
                    let at = sub_zone_at_index(z, depth, index);
                    assert_eq!(
                        at.is_some(),
                        index < n,
                        "{z:#018x} at {depth}, index {index}"
                    );
                }
                if depth <= 12 {
                    let zone = ZoneId(z);
                    let _ = grid.first_sub_zone(zone, depth);
                    let _ = grid.sub_zone_at_index(zone, depth, n / 2);
                    if n <= 400 {
                        let _ = grid.sub_zones(zone, depth);
                    }
                }
            }
        }
        assert!(answered > 0, "no identifier has a sub-zone order");
        assert!(carried > 0, "no identifier reaches the exact corner");
    }

    /// Asserts against the live engine that its sub-zone order of `z` at `depth` is not the
    /// zone's descendants, and that this crate's is. The engine's names a zone twice, names one
    /// that its own `zoneHasSubZone` rejects, or holds fewer distinct zones than its own count.
    /// This crate's names each zone once, each accepted by the engine's `zoneHasSubZone`, as
    /// many as the engine's count, at depth 1 the engine's children as a set, and has none of
    /// the [`faults`]. The engine is asked about no zone it cannot read back, and about no
    /// sub-zone at level 33, whose descendant test it cannot answer.
    #[cfg(feature = "oracle")]
    fn assert_corrected(grid: &Grid<Isea, HexA3, I3h>, z: u64, depth: u8) {
        use dggal_oracle as o;
        use std::collections::BTreeSet;
        let at = format!("{} at {depth}", grid.zone(ZoneId(z)).text_id());
        assert!(
            level(z) + u32::from(depth) < u32::from(MAX_LEVEL),
            "{at}: the sub-zones lie at level 33, where the engine's descendant test fails"
        );
        let has = |s: u64| {
            assert!(
                o::engine_can_read(o::ISEA3H, s),
                "{at}: the engine cannot read {s:#018x} back"
            );
            o::zone_has_sub_zone(o::ISEA3H, z, s)
        };
        let count = o::count_sub_zones(o::ISEA3H, z, i32::from(depth));
        let theirs = o::sub_zones(o::ISEA3H, z, i32::from(depth));
        let distinct: BTreeSet<u64> = theirs.iter().copied().collect();
        assert!(
            distinct.len() < theirs.len()
                || (distinct.len() as u64) < count
                || !theirs.iter().all(|&s| has(s)),
            "{at}: the engine's order is the zone's descendants, and this crate's departs from it"
        );
        let ours = sub_zones(z, depth).unwrap();
        let ours_distinct: BTreeSet<u64> = ours.iter().copied().collect();
        assert_eq!(ours_distinct.len(), ours.len(), "{at}: a zone named twice");
        assert_eq!(ours.len() as u64, count, "{at}: against the engine's count");
        if let Some(s) = ours.iter().find(|&&s| !has(s)) {
            panic!("{at}: the engine's zoneHasSubZone rejects {s:#018x}");
        }
        if depth == 1 {
            let children: BTreeSet<u64> = o::children(o::ISEA3H, z).into_iter().collect();
            assert_eq!(
                ours_distinct, children,
                "{at}: against the engine's children"
            );
        }
        assert_eq!(faults(grid, z, depth), Vec::<String>::new(), "{at}");
    }

    /// The raw centroids of the four generators, before their quantisation, against the
    /// engine's, bit for bit and index by index, through its `getSubZoneCRSCentroids` in the
    /// 5x6 CRS, which hands back its generators' output untouched (`RI3H.ec:676-705`). The
    /// sub-zone identifiers cannot see the order of the arithmetic, since a sub-zone centroid
    /// is a cell centre, far from any boundary of its quantisation; the centroids can. So this
    /// is the check that the compiled orders ported in `generateEvenParentOddDepth` and
    /// `generateOddParentEvenDepth` (the `u3` rule and the sums gcc regrouped) and in
    /// `generateEvenParentEvenDepth` (the sums gcc regrouped) are the engine's, and that nothing
    /// was lost in the transcription of any generator.
    ///
    /// Sample: at every level from 0 to 32, the two polar pentagons, the ten other pentagons,
    /// every edge hexagon to level 8 and above it those at both ends of the edge of every
    /// rhombus, and the zones in the middle of every rhombus (at an odd level, its sub-hexes B,
    /// C and D); at the depths 1 to 6 that stay within level 33, at depth 7 to level 4 and, for
    /// an even parent, at depth 8 to level 2. All twenty-eight cases of the four generators are
    /// reached at each of their depths. Reverting to the eC's order the `u3` rule of either
    /// generator that has it, any of the sums gcc regrouped at `I3HSubZones.ec:688`, `:698`,
    /// `:706-707`, `:1041-1042`, `:1046-1047`, `:1223-1224`, `:1255-1256`, `:1517-1518` and
    /// `:1527-1528`, or the crossing formulae of `move5x6Vertex3`, fails this test; depth 8 is
    /// there for the first sum of `:1041` and the second of `:1047`, which first part from the
    /// eC's order there. It cannot see the corner of [`interruption_corner`], whose two orders
    /// agree on every input that reaches it, nor the sums of `:1223` and `:1256`, which are the
    /// same double in either order, since `nMidRows >> 1` and `(nMidRows - 1) >> 1` are both
    /// `nCapRows`.
    ///
    /// Where the engine's order is not the zone's descendants, this crate's departs from it (see
    /// the module note), and the centroids differ: at each order in which any does, the
    /// engine's order is checked to be faulty and this crate's to be sound, by the engine's own
    /// tests ([`assert_corrected`]), and the orders and the centroids that differ are counted
    /// exactly. Every other order is the engine's, bit for bit.
    ///
    /// Gated behind the `oracle` feature, meaningful only within this repository, as the oracle
    /// test of `hex_a7.rs` is.
    #[test]
    #[cfg(feature = "oracle")]
    fn the_generators_centroids_match_the_engines_bit_for_bit() {
        use crate::indexings::pack;
        use dggal_oracle as o;
        let mut zones = std::collections::BTreeSet::new();
        for level in 0..=32u32 {
            let level_i9r = u64::from(level / 2);
            let odd = u64::from(level & 1);
            let p = pow3(level_i9r);
            zones.extend([10, 11].map(|root| (level, pack(level_i9r, root, 0, odd))));
            for root in 0..10 {
                zones.insert((level, pack(level_i9r, root, 0, odd)));
                // Every edge hexagon to level 8, then those at both ends of the edge.
                let edge: Vec<u64> = if level <= 8 {
                    (1..p).collect()
                } else {
                    vec![1, p - 1]
                };
                for j in edge {
                    let rix = if root & 1 == 0 { j } else { p * j };
                    zones.insert((level, pack(level_i9r, root, rix, odd)));
                }
                let sub_hexes: &[u64] = if odd == 1 { &[1, 2, 3] } else { &[0] };
                for &sub_hex in sub_hexes {
                    zones.insert((level, pack(level_i9r, root, p * p / 2, sub_hex)));
                }
            }
        }
        let grid = isea3h();
        let (mut pairs, mut centroids, mut unreadable) = (0usize, 0usize, 0usize);
        let mut reached = std::collections::BTreeSet::new();
        // The orders corrected, and the centroids that differ in them.
        let mut corrected = [0usize; 2];
        for (level, z) in zones {
            if !o::engine_can_read(o::ISEA3H, z) {
                unreadable += 1;
                continue;
            }
            for depth in 1..=8 {
                if level + depth > 33
                    || (depth == 7 && level > 4)
                    || (depth == 8 && (level > 2 || level & 1 != 0))
                {
                    continue;
                }
                let mut ours = vec![];
                let stop = iterate(
                    z,
                    depth,
                    &mut |i, c| {
                        ours.push((i, c));
                        true
                    },
                    -1,
                );
                assert_eq!(stop, -1, "{z:#018x} at {depth}");
                let theirs = o::sub_zone_crs_centroids_5x6(o::ISEA3H, z, depth as i32);
                assert_eq!(ours.len(), theirs.len(), "{z:#018x} at {depth}");
                pairs += 1;
                reached.insert((case(z, depth), depth));
                let mut differ = 0;
                for (k, (&(i, a), &b)) in ours.iter().zip(&theirs).enumerate() {
                    assert_eq!(
                        i, k as i64,
                        "{z:#018x} at {depth}: the order of the indices"
                    );
                    centroids += 1;
                    if (a.0.to_bits(), a.1.to_bits()) != (b.0.to_bits(), b.1.to_bits()) {
                        differ += 1;
                    }
                }
                if differ != 0 {
                    assert_corrected(&grid, z, depth as u8);
                    corrected[0] += 1;
                    corrected[1] += differ;
                }
            }
        }
        // Measured, and asserted exactly, so that the test cannot quietly compare less than it
        // claims, nor correct more or fewer orders than it did when measured.
        assert_eq!((pairs, centroids, unreadable), (17_710, 4_105_210, 0));
        assert_eq!(
            corrected,
            [98, 6_884],
            "the orders corrected, and their centroids"
        );
        let odd_depths = [1, 5, 6, 13, 14, 21, 22, 2, 7, 8, 15, 16, 23, 24];
        let odd_parents = [4, 11, 12, 19, 20, 27, 28];
        let even_parents = [3, 9, 10, 17, 18, 25, 26];
        let want: std::collections::BTreeSet<_> = odd_depths
            .iter()
            .flat_map(|&c| [1, 3, 5, 7].map(|d| (c, d)))
            .chain(odd_parents.iter().flat_map(|&c| [2, 4, 6].map(|d| (c, d))))
            .chain(
                even_parents
                    .iter()
                    .flat_map(|&c| [2, 4, 6, 8].map(|d| (c, d))),
            )
            .collect();
        assert_eq!(reached, want, "the cases and depths reached");
    }
}
