//! The helpers of the 5x6 planar layout, shared by the aperture-7 and aperture-3 topologies.
//!
//! Each ports a function of `projections/ri5x6.ec` in DGGAL v0.0.6. Three of them,
//! `move5x6Vertex`, `move5x6Vertex2` and `crosses5x6Interruption`, carry the eC's
//! `__attribute__ ((optimize("-fno-unsafe-math-optimizations")))`, and the library performs them
//! in the order of the eC's text, which is the order ported here. The other four,
//! `intersects5x6Interruption`, `cross5x6Interruption`, `move5x6Vertex3` and `canonicalize5x6`,
//! carry no such mark, and DGGAL is built with `-O2 -ffast-math` (`Makefile.dggal:133`), so gcc
//! was free to re-associate them, and in places it did. At every such site this module performs
//! the order read out of the shipped `libdggal.so`, BuildID
//! `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, and a `// faithful:` comment names the address and
//! the eC line it replaces. Several of those sites cannot change an answer for the arguments this
//! crate passes; they follow the library all the same, so that the port does not rest on an
//! argument about its callers which the next caller may invalidate. A guarded helper is still
//! reached from unguarded callers, which compute its arguments under unsafe mathematics; that is
//! the callers' concern, and each caller ports its own order.
//!
//! The aperture-7 topology uses [`move5x6_vertex`], [`sgn`] and [`Crossing`]; the others serve the
//! aperture-3 topology, whose vertices, neighbours and sub-zones call them. The aperture-7
//! quantiser's domain fold, which the eC performs inline in `I7HZone::fromCentroid`, stays in
//! `topologies/hex_a7.rs`: it bears the name of the eC function ported here as
//! [`canonicalize5x6`], but it is a different routine.

/// `fl(1 - 1e-11)`, a constant gcc folded into `cross5x6Interruption`; `.rodata` at `0x4a018`.
const ONE_MINUS_1E11: f64 = 1.0 - 1e-11;
/// `fl(1 + 1e-11)`, folded likewise; `.rodata` at `0x4a010`.
const ONE_PLUS_1E11: f64 = 1.0 + 1e-11;
/// `fl(2 - 1e-11)`, folded likewise; `.rodata` at `0x4a258`.
const TWO_MINUS_1E11: f64 = 2.0 - 1e-11;
/// `fl(1 - 1e-12)`, folded into `intersects5x6Interruption`; `.rodata` at `0x4a398`.
const ONE_MINUS_1E12: f64 = 1.0 - 1e-12;
/// The tolerance `e` of `intersects5x6Interruption` (ri5x6.ec:1130); `.rodata` at `0x4a390`.
const INTERSECTION_EPSILON: f64 = 1e-12;

/// A segment of the 5x6 plane, by its two end points.
type Segment = [(f64, f64); 2];

/// The twenty interruption segments of the 5x6 layout, as `crosses5x6Interruption` lists them
/// (ri5x6.ec:1307-1323): by hemisphere, north first, then by root rhombus, then the left and the
/// right segment of each, each segment by its two end points.
const INTERRUPTIONS: [[[Segment; 2]; 5]; 2] = [
    [
        [[(0.0, 0.0), (1.0, 0.0)], [(1.0, 0.0), (1.0, 1.0)]],
        [[(1.0, 1.0), (2.0, 1.0)], [(2.0, 1.0), (2.0, 2.0)]],
        [[(2.0, 2.0), (3.0, 2.0)], [(3.0, 2.0), (3.0, 3.0)]],
        [[(3.0, 3.0), (4.0, 3.0)], [(4.0, 3.0), (4.0, 4.0)]],
        [[(4.0, 4.0), (5.0, 4.0)], [(5.0, 4.0), (5.0, 5.0)]],
    ],
    [
        [[(0.0, 1.0), (0.0, 2.0)], [(0.0, 2.0), (1.0, 2.0)]],
        [[(1.0, 2.0), (1.0, 3.0)], [(1.0, 3.0), (2.0, 3.0)]],
        [[(2.0, 3.0), (2.0, 4.0)], [(2.0, 4.0), (3.0, 4.0)]],
        [[(3.0, 4.0), (3.0, 5.0)], [(3.0, 5.0), (4.0, 5.0)]],
        [[(4.0, 5.0), (4.0, 6.0)], [(4.0, 6.0), (5.0, 6.0)]],
    ],
];

/// The eC's `Sgn`, which is zero at zero.
///
/// `Sgn` is not a function of the eC runtime but an intrinsic of the eC compiler, which rewrites
/// `Sgn(a)` as `((!(a))?(0):(((a)<0)?(-1):(1)))` (the eC compiler's `ectp/src/pass15.ec:9807-9815`,
/// in `dgbuild`). The shipped
/// library bears this out in both guarded functions that use it: `move5x6Vertex` compares its
/// offset with zero and adds no epsilon when they are equal (`0x379ed` and `0x379fd` for `dx`,
/// `0x37a72` and `0x37a7e` for `dy`), and `crosses5x6InterruptionV2` does the same (`0x38946`
/// and `0x3894a` for `dx`, `0x389be` and `0x389c2` for `dy`). A negative zero counts as zero,
/// as `!(a)` has it.
pub(crate) fn sgn(a: f64) -> f64 {
    if a == 0.0 {
        0.0
    } else if a < 0.0 {
        -1.0
    } else {
        1.0
    }
}

/// The eC's `(int)` of a double, a 32-bit `cvttsd2si` in the library: a truncation, and the
/// integer indefinite `i32::MIN` for NaN and for any value outside the range of `int`, where
/// Rust's `as` would saturate. The eC writes it inline, as in `move5x6Vertex` (`ri5x6.ec:1390`).
pub(crate) fn cvtt_i32(v: f64) -> i32 {
    if v > -2_147_483_649.0 && v < 2_147_483_648.0 {
        v as i32
    } else {
        i32::MIN
    }
}

/// A crossing of a rhombus interruption: the point at which a step leaves one rhombus, the point
/// at which it re-enters the layout on the far side of the interruption, and whether the
/// interruption is a northern one. [`crosses5x6_interruption`] returns it, and so does the
/// aperture-7 walker's `crosses5x6InterruptionV2`.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Crossing {
    pub(crate) i_src: (f64, f64),
    pub(crate) i_dst: (f64, f64),
    pub(crate) in_north: bool,
}

/// Offsets a vertex from the centroid `(cx, cy)` by `(dx, dy)`, re-expressing the offset across
/// a rhombus interruption when it crosses one. Ports `move5x6Vertex` (ri5x6.ec:1388-1434), and
/// not `move5x6Vertex2` (ri5x6.ec:1470), which is [`move5x6_vertex2`], a genuinely different
/// routine that the eC itself warns "does not generate correct geometry at level 2"
/// (RI7H.ec:2954). The function is guarded, and the library performs it in the eC's order
/// (`0x379b0` to `0x37bdc`).
///
/// The four `ivx`/`ivy` branches mirror the four crossing directions in the eC and deliberately
/// share their bodies; they are left apart rather than merged, for faithfulness.
#[allow(
    clippy::if_same_then_else,
    reason = "the eC keeps these four crossing directions as separate branches with shared bodies; merging them would break the line-by-line correspondence"
)]
pub(crate) fn move5x6_vertex(cx: f64, cy: f64, dx: f64, dy: f64) -> (f64, f64) {
    let mut vx = cx + dx;
    let mut vy = cy + dy;
    // faithful: `int cx = (int)(c.x + 1E-11)` (ri5x6.ec:1390) truncates, and so does the library,
    // DGGAL's `-O2 -ffast-math` build with this function under its guard (`libdggal.so`, BuildID
    // `e75d6ab18f8b460713ebcd21df6f07922fb909c3`), with a bare `cvttsd2si` at `0x379e2` and
    // `0x379f1`. py4dggs's aperture-7 copy floors
    // instead, as this crate's did; on aperture 7 the two agree, since its one caller passes
    // coordinates that are not negative, but a point left of the origin tells them apart:
    // `(-0.5, 0.2)` moved by `(0, 1)` is `(-1.2, 0.0)` in the engine, and `(0.0, 0.5)` with a
    // floor. The `ivx`/`ivy` pair below does floor, as the eC does (ri5x6.ec:1394-1395), with
    // the floor inlined at `0x37a15` to `0x37a61`. All four are the eC's 32-bit `int`, converted
    // by a 32-bit `cvttsd2si` (`0x37a6e` for `ivx`), so a coordinate beyond the range of `int`,
    // far outside the layout, gives the integer indefinite, as in the engine.
    let icx = cvtt_i32(cx + 1e-11);
    let icy = cvtt_i32(cy + 1e-11);
    // faithful: `Sgn(0)` is 0 (see [`sgn`]), so an offset of zero takes no epsilon
    // (ri5x6.ec:1394-1395, `0x379fd` and `0x37a7e` in the same library). Aperture 3
    // passes such offsets (`RI3H.ec:1738`, `:1763-1769`), and on an integer coordinate the
    // epsilon would floor to the cell below; aperture 7 never passes one.
    let ivx = cvtt_i32((cx + dx - sgn(dx) * 1e-11).floor());
    let ivy = cvtt_i32((cy + dy - sgn(dy) * 1e-11).floor());

    if ((ivx != icx && (vy - f64::from(ivy)).abs() > 1e-11)
        || (ivy != icy && (vx - f64::from(ivx)).abs() > 1e-11))
        && (ivy.wrapping_sub(ivx) > 1 || ivy < ivx)
    {
        let fcx = cx - f64::from(icx);
        let fcy = cy - f64::from(icy);
        if ivx < icx {
            vx = f64::from(icx) - fcy + dx - dy;
            vy = f64::from(icy) + dx;
        } else if ivx > icx {
            vx = f64::from(icx) - fcy + dx - dy;
            vy = f64::from(icy) + dx;
        } else if ivy < icy {
            vx = f64::from(icx) + dy;
            vy = f64::from(icy) - fcx - dx + dy;
        } else if ivy > icy {
            vx = f64::from(icx) + dy;
            vy = f64::from(icy) - fcx - dx + dy;
        }
    }
    (vx, vy)
}

/// Offsets the point `(sx, sy)` by `(dx, dy)`, re-expressing a step that crosses a rhombus
/// interruption through the point where it meets it, snapping an answer on a polar diagonal to
/// its pole and, with `cross_early`, first carrying a point that lies on an interruption across
/// it. Ports `move5x6Vertex2` (ri5x6.ec:1470-1560), which this crate reaches from the ports of
/// `I3HZone::getNeighbor` (RI3H.ec:1146) and `I7HZone::getNeighbors` (RI7H.ec:1206). The
/// function is guarded, and is ported in the eC's order.
///
/// The eC assigns the new centre of an early crossing as a compound literal, such as
/// `c = { cx, cy - (c.x - cx) }`, which reads `c.x` in the member it does not assign first. The
/// engine computes both members from the old centre before storing either, and the tuple
/// assignments below do the same.
pub(crate) fn move5x6_vertex2(sx: f64, sy: f64, dx: f64, dy: f64, cross_early: bool) -> (f64, f64) {
    let e = 1e-11;
    let (mut x, mut y) = (sx, sy);
    let mut cx = (x + e).floor() as i64;
    let mut cy = (y + e).floor() as i64;
    let cx2 = (x - e).floor() as i64;
    let cy2 = (y - e).floor() as i64;
    let (mut nx, mut ny) = (x + dx, y + dy);

    if nx < 0.0 {
        nx += 5.0;
    } else if nx > 5.0 {
        nx -= 5.0;
    }
    if ny < 0.0 {
        ny += 5.0;
    } else if ny > 5.0 && y > 6.0 - e {
        ny -= 5.0;
    }

    let n_top_right_of_p = (nx > x + e && nx - x < 3.0) || x - nx > 3.0;
    let n_top_left_of_p = (nx < x - e && x - nx < 3.0) || nx - x > 3.0;
    let n_bottom_right_of_p = (ny > y && ny - y < 3.0) || y - ny > 3.0;
    let n_bottom_left_of_p = (ny < y - e && y - ny < 3.0) || ny - y > 3.0;
    let at_top_dent_crossing_right = cx2 != cx && x > y + e && n_top_right_of_p;
    let at_top_dent_crossing_left = cy2 != cy && x > y + e && n_top_left_of_p;
    let at_bottom_dent_crossing_left = cx2 != cx && y > x + 1.0 + e && n_bottom_left_of_p;
    let at_bottom_dent_crossing_right = cy2 != cy && y > x + 1.0 + e && n_bottom_right_of_p;

    // Cross already for cases where crossing does not happen mid-edge.
    if cross_early {
        if at_top_dent_crossing_right {
            (x, y) = (cx as f64 + 1.0 - (y - cy as f64), (cy + 1) as f64);
        } else if at_top_dent_crossing_left {
            (x, y) = (cx as f64, cy as f64 - (x - cx as f64));
        } else if at_bottom_dent_crossing_left {
            (x, y) = (cx as f64 - (y - cy as f64), cy as f64);
        } else if at_bottom_dent_crossing_right {
            (x, y) = ((cx + 1) as f64, cy as f64 + ((cx + 1) as f64 - x));
        }

        if x > 5.0 || y > 6.0 + e {
            x -= 5.0;
            y -= 5.0;
        } else if x < 0.0 || y < -e {
            x += 5.0;
            y += 5.0;
        }
        cx = (x + e).floor() as i64;
        cy = (y + e).floor() as i64;
    }

    let (mut vx, mut vy) = (x + dx, y + dy);
    let ivx = (x + dx + 1e-11).floor() as i64;
    let ivy = (y + dy + 1e-11).floor() as i64;

    if ((ivx != cx && (vy - ivy as f64).abs() > 1e-11)
        || (ivy != cy && (vx - ivx as f64).abs() > 1e-11))
        && (ivy - ivx > 1 || ivy < ivx)
    {
        // Assuming the crossing point is at half the distance.
        if (vx - vy - 1.0).abs() < 1e-10 {
            (vx, vy) = (1.0, 0.0); // "North" pole
        } else if (vy - vx - 2.0).abs() < 1e-10 {
            (vx, vy) = (4.0, 6.0); // "South" pole
        } else if ivx < cx && vx - (ivx as f64) < 1.0 - e {
            // Stepping over bottom dent to the left
            let pi1 = (cx as f64, y + 0.5 * (cx as f64 - x));
            let pi2 = (cx as f64 - (pi1.1 - cy as f64), cy as f64);
            vx = pi2.0 + pi1.1 - y;
            vy = pi2.1 + pi1.0 - x;
        } else if ivx > cx && vx - (ivx as f64) > e {
            // Stepping over top dent to the right
            let pi1 = ((cx + 1) as f64, y + 0.5 * ((cx + 1) as f64 - x));
            let pi2 = ((cx + 1) as f64 + ((cy + 1) as f64 - pi1.1), (cy + 1) as f64);
            vx = pi2.0 + pi1.1 - y;
            vy = pi2.1 + pi1.0 - x;
        } else if ivy < cy && vy - (ivy as f64) < 1.0 - e {
            // Stepping over top dent to the left
            let pi1 = (x + 0.5 * (cy as f64 - y), cy as f64);
            let pi2 = (cx as f64, cy as f64 - (pi1.0 - cx as f64));
            vx = pi2.0 + pi1.1 - y;
            vy = pi2.1 + pi1.0 - x;
        } else if ivy > cy && vy - (ivy as f64) > e {
            // Stepping over bottom dent to the right. The eC starts `pi1` from `v.x` here where
            // its three siblings start from the centre; faithful, and kept.
            let pi1 = (vx + 0.5 * ((cy + 1) as f64 - y), (cy + 1) as f64);
            let pi2 = ((cx + 1) as f64, (cy + 1) as f64 + ((cx + 1) as f64 - pi1.0));
            vx = pi2.0 + pi1.1 - y;
            vy = pi2.1 + pi1.0 - x;
        }
    }
    if vx > 5.0 {
        vx -= 5.0;
        vy -= 5.0;
    } else if vx < 0.0 {
        vx += 5.0;
        vy += 5.0;
    }
    (vx, vy)
}

/// Offsets the point `(cx, cy)` by `(dx, dy)` and, when the step crosses a rhombus interruption,
/// re-expresses the remainder of the step on the far side of it, assuming a crossing to the
/// right. Ports `move5x6Vertex3` (ri5x6.ec:1438-1463), which the sub-zone generators call; the
/// eC's own comment notes that it does not have safe optimisations disabled (ri5x6.ec:1436-1437),
/// and it is compiled at `0x3bef0` to `0x3c046`.
pub(crate) fn move5x6_vertex3(cx: f64, cy: f64, dx: f64, dy: f64) -> (f64, f64) {
    let (mut vx, mut vy);
    if let Some(crossing) = crosses5x6_interruption(cx, cy, dx, dy) {
        let (i1, i2) = (crossing.i_src, crossing.i_dst);
        // faithful: the function is unguarded, and in DGGAL's `-O2 -ffast-math` build gcc
        // re-associated all four assignments of the crossing branch (ri5x6.ec:1447-1453). In
        // `libdggal.so`, BuildID
        // `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, the northern pair is
        //
        // ```text
        //   3bf70  subsd  xmm1,xmm15   ; i1.y - dy
        //   3bf75  subsd  xmm5,xmm6    ; i2.y - i1.x
        //   3bf79  subsd  xmm1,xmm4    ; (i1.y - dy) - c.y
        //   3bf7d  addsd  xmm1,xmm1    ; doubled
        //   3bf81  addsd  xmm0,xmm1    ; v.x = i2.x + that
        //   3bf8a  addsd  xmm3,[rbx]   ; dx + c.x
        //   3bf92  addsd  xmm2,xmm5    ; v.y = (dx + c.x) + (i2.y - i1.x)
        // ```
        //
        // for the eC's `i2.x - 2 * (dy - (i1.y - c.y))` and `i2.y + dx - (i1.x - c.x)`, and the
        // southern pair is
        //
        // ```text
        //   3bf45  addsd  xmm2,xmm15   ; c.y + dy, hoisted above the call's test
        //   3c020  subsd  xmm0,xmm1    ; i2.x - i1.y
        //   3c024  subsd  xmm3,xmm6    ; dx - i1.x
        //   3c028  addsd  xmm0,xmm2    ; v.x = (i2.x - i1.y) + (c.y + dy)
        //   3c035  addsd  xmm2,xmm3    ; c.x + (dx - i1.x)
        //   3c039  addsd  xmm2,xmm2    ; doubled
        //   3c03d  addsd  xmm2,xmm5    ; v.y = that + i2.y
        // ```
        //
        // for `i2.x + dy - (i1.y - c.y)` and `i2.y + 2 * (dx - (i1.x - c.x))`. A doubling by
        // `addsd x,x` is the same double as a product by two. On 300,000 calls near the
        // interruptions, 16,995 of them crossing, the eC's order moved the answer at 8,557.
        if crossing.in_north {
            vx = i2.0 + 2.0 * ((i1.1 - dy) - cy);
            vy = (dx + cx) + (i2.1 - i1.0);
        } else {
            vx = (i2.0 - i1.1) + (cy + dy);
            vy = 2.0 * (cx + (dx - i1.0)) + i2.1;
        }
    } else {
        // `{ c.x + dx, c.y + dy }`, formed as `dx + c.x` at `0x3bfa0`: the same double.
        (vx, vy) = (dx + cx, cy + dy);
    }

    if vx > 5.0 && vy > 5.0 {
        vx -= 5.0;
        vy -= 5.0;
    } else if vx < 0.0 && vy < 1.0 {
        vx += 5.0;
        vy += 5.0;
    }
    (vx, vy)
}

/// Whether the segment from `a0` to `a1` meets the interruption segment from `b0` to `b1`, and if
/// so where, with the distance along the first segment as a fraction of its length. Ports
/// `intersects5x6Interruption` (ri5x6.ec:1128-1155), which is not inlined into its one caller,
/// [`crosses5x6_interruption`], but compiled on its own at `0x3ba30` to `0x3bb1b`, its eight
/// coordinates arriving in registers.
pub(crate) fn intersects5x6_interruption(
    a0: (f64, f64),
    a1: (f64, f64),
    b0: (f64, f64),
    b1: (f64, f64),
) -> Option<((f64, f64), f64)> {
    let (s1y, s1x) = (a1.1 - a0.1, a1.0 - a0.0);
    let (s2y, s2x) = (b1.1 - b0.1, b1.0 - b0.0);
    let (dy, dx) = (a0.1 - b0.1, a0.0 - b0.0);
    let d = s1x * s2y - s2x * s1y;

    // Parallel segments do not meet.
    if d.abs() > 1e-13 {
        let factor = 1.0 / d;
        let s = (s1x * dy - s1y * dx) * factor;
        // faithful: the eC tests `s - e >= 0 && s + e <= 1` (ri5x6.ec:1140), and the library,
        // which compiles this unguarded function under `-O2 -ffast-math`, compares `s` with `e`
        // and with `fl(1 - e)`: `comisd` against `.rodata` `0x4a390`
        // (`1e-12`) at `0x3bab4` and `0x4a398` (`fl(1 - 1e-12)`) at `0x3babe`, in
        // `libdggal.so`, BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`. The first pair is
        // exactly equivalent; the second switches at the same double as the eC's, which the test
        // below pins, so the rewriting changes no answer, and it is ported as compiled.
        if (INTERSECTION_EPSILON..=ONE_MINUS_1E12).contains(&s) {
            let t = (s2x * dy - s2y * dx) * factor;
            // faithful: the same for `t` (ri5x6.ec:1143), at `0x3bad9` and `0x3bae3`.
            if (INTERSECTION_EPSILON..=ONE_MINUS_1E12).contains(&t) {
                return Some(((a0.0 + t * s1x, a0.1 + t * s1y), t));
            }
        }
    }
    None
}

/// Whether a step of `(dx, dy)` from the point `(cx, cy)` crosses one of the twenty rhombus
/// interruptions and, if it crosses several, the one it meets first, with its far-side image.
/// Ports `crosses5x6Interruption` (ri5x6.ec:1304-1374), which is guarded and ported in the eC's
/// order; the two helpers it calls are not guarded, and follow the library.
pub(crate) fn crosses5x6_interruption(cx: f64, cy: f64, dx: f64, dy: f64) -> Option<Crossing> {
    let (mut x, mut y) = (cx, cy);
    // A start within 1e-12 of a whole number is moved onto it; the eC's chains of `else if`
    // (ri5x6.ec:1327-1339) take the first whole number that matches, as the loops do.
    for k in 0..=5 {
        if (x - k as f64).abs() < 1e-12 {
            x = k as f64;
            break;
        }
    }
    for k in 0..=6 {
        if (y - k as f64).abs() < 1e-12 {
            y = k as f64;
            break;
        }
    }

    // +1 avoids rounding the wrong way for < 0
    if (x + 1.0) as i64 == (x + dx + 1.0) as i64 && (y + 1.0) as i64 == (y + dy + 1.0) as i64 {
        return None;
    }

    let mut min_d = f64::MAX;
    let mut nearest: Option<((f64, f64), usize, usize)> = None;
    for (h, hemisphere) in INTERRUPTIONS.iter().enumerate() {
        for rhombus in hemisphere {
            for (s, segment) in rhombus.iter().enumerate() {
                // NOTE: This will currently not detect intersections from edge to edge across
                // the interruption (the eC's own note, ri5x6.ec:1353).
                if let Some((i_cur, d)) =
                    intersects5x6_interruption((x, y), (x + dx, y + dy), segment[0], segment[1])
                {
                    if d < min_d {
                        nearest = Some((i_cur, h, s));
                        min_d = d;
                    }
                }
            }
        }
    }
    let (i_src, cross_h, cross_s) = nearest?;
    Some(Crossing {
        i_src,
        i_dst: cross5x6_interruption(i_src.0, i_src.1, cross_h == 1, cross_s == 0),
        in_north: cross_h == 0,
    })
}

/// The image, on the far side of a rhombus interruption, of the point `(sx, sy)` lying on it:
/// `south` names the hemisphere and `left` the side of the rhombus crossed. Ports
/// `cross5x6Interruption` (ri5x6.ec:1247-1299), compiled at `0x360d0` to `0x36288`.
///
/// faithful: the function is unguarded, and in DGGAL's `-O2 -ffast-math` build gcc re-associated
/// each of its four arms; the arms below perform the library's order, in `libdggal.so`, BuildID
/// `e75d6ab18f8b460713ebcd21df6f07922fb909c3`. For a point on
/// an interruption segment, which is what this crate's callers pass, each gives the eC's answer;
/// for an arbitrary point the two orders part at hundreds of points in a hundred thousand, and
/// the tests below pin a few of them. The three forms `k - (s - i)`, compiled as `(i - s) + k`,
/// are the same double either way, since IEEE subtraction is exactly antisymmetric; they are
/// written as compiled.
pub(crate) fn cross5x6_interruption(sx: f64, sy: f64, south: bool, left: bool) -> (f64, f64) {
    let (mut x, mut y) = match (south, left) {
        // Crossing the northern hemisphere to the left.
        (false, true) => {
            let ix = (sy + 1e-11) as i64 as f64;
            // faithful: `ix - (iSrc.x - ix)` (ri5x6.ec:1266) as `(ix + ix) - iSrc.x`
            // (`0x361e6 addsd xmm1,xmm1`, `0x361ea subsd xmm1,[rdi]`).
            (sy, (ix + ix) - sx)
        }
        // Crossing the northern hemisphere to the right.
        (false, false) => {
            // faithful: `(int)(iSrc.x - 1 + 1E-11)` (ri5x6.ec:1271) as
            // `(int)(iSrc.x - fl(1 - 1e-11))`, the constant folded to `.rodata` `0x4a018`
            // (`0x3620c`); then `(iy - iSrc.y) + (iy + 2)` (`0x36223`, `0x3622c`).
            let iy = (sx - ONE_MINUS_1E11) as i64;
            ((iy as f64 - sy) + (iy + 2) as f64, sx)
        }
        // Crossing the southern hemisphere to the left.
        (true, true) => {
            // faithful: `(int)(iSrc.x + 1 + 1E-11)` (ri5x6.ec:1282) as
            // `(int)(fl(1 + 1e-11) + iSrc.x)`, the constant at `.rodata` `0x4a010` (`0x3616d`
            // to `0x36185`); then `(iy - iSrc.y) + (iy - 1)` (`0x36194`, `0x3619d`).
            let iy = (ONE_PLUS_1E11 + sx) as i64;
            ((iy as f64 - sy) + (iy - 1) as f64, sx + 1.0)
        }
        // Crossing the southern hemisphere to the right.
        (true, false) => {
            // faithful: `(int)(iSrc.y - 2 + 1E-11)` (ri5x6.ec:1288) as
            // `(int)(iSrc.y - fl(2 - 1e-11))`, the constant at `.rodata` `0x4a258` (`0x3624d`);
            // then `(ix - iSrc.x) + (ix + 3)` (`0x3626c`, `0x36274`).
            let ix = (sy - TWO_MINUS_1E11) as i64;
            (sy - 1.0, (ix as f64 - sx) + (ix + 3) as f64)
        }
    };
    if x > 5.0 - 1e-11 && y > 5.0 - 1e-11 {
        x -= 5.0;
        y -= 5.0;
    } else if x < 1e-11 && y < 1.0 - 1e-11 {
        x += 5.0;
        y += 5.0;
    }
    (x, y)
}

/// Folds a 5x6 point into the canonical form the engine gives every vertex before the inverse
/// projection: back into the layout from beyond its two ends, onto a pole from any of its images,
/// and from an interruption to its image on the far side. Ports the eC's
/// `canonicalize5x6` (ri5x6.ec:1562-1606), compiled at `0x3a080` to `0x3a494`; the
/// aperture-7 quantiser's inline fold of the same name, in `topologies/hex_a7.rs`, is another
/// routine.
pub(crate) fn canonicalize5x6(sx: f64, sy: f64) -> (f64, f64) {
    let (mut x, mut y) = (sx, sy);
    // faithful: the eC tests the second fold after the first (ri5x6.ec:1569-1572), and the
    // library, which compiles this unguarded function under `-O2 -ffast-math`, does not: once
    // the first fold fires (`0x3a0ac`, `0x3a0b6`, `0x3a178`) it goes straight to the floors, in
    // `libdggal.so`, BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`. It cannot matter: a
    // coordinate the first fold moves exceeded `fl(5 - 1e-11)`, so it is at least a unit in the
    // last place above it, and after the fold it is at least `-9.9991e-12`, which the second
    // test does not catch. It is ported as compiled all the same.
    if x > 5.0 - 1e-11 && y > 5.0 - 1e-11 {
        x -= 5.0;
        y -= 5.0;
    } else if x < -1e-11 || y < -1e-11 {
        x += 5.0;
        y += 5.0;
    }

    let cx = (x + 1e-11).floor() as i64;
    let cy = (y + 1e-11).floor() as i64;
    let (mut south, mut cross, mut np, mut sp) = (false, false, false, false);

    // The eC's two `switch` statements, whose cases 0 to 4 each test the line through that
    // cell: `fabs(src.x - (cy + 1)) < 1E-11` with the north pole at `src.y == cy`, and
    // `fabs(src.y - (cx + 2)) < 1E-11` with the south pole at `src.x == cx`.
    if (0..=4).contains(&cy) {
        cross = (x - (cy + 1) as f64).abs() < 1e-11;
        np = cross && (y - cy as f64).abs() < 1e-11;
    }
    if (0..=4).contains(&cx) && (y - (cx + 2) as f64).abs() < 1e-11 {
        cross = true;
        south = true;
        sp = (x - cx as f64).abs() < 1e-11;
    }

    if sp {
        (4.0, 6.0)
    } else if np {
        (1.0, 0.0)
    } else if cross {
        cross5x6_interruption(x, y, south, false)
    } else if (x - 5.0).abs() < 1e-11 {
        (0.0, y - 5.0)
    } else {
        (x, y)
    }
}

#[cfg(test)]
#[allow(clippy::type_complexity, reason = "the tables of test vectors")]
mod tests {
    //! Every expected value below is the shipped engine's own answer, read by calling the
    //! function at its address in `libdggal.so`, BuildID
    //! `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, and written out as the shortest decimal that
    //! reads back to the same double; the comparisons are bit for bit.
    use super::*;

    fn bits(p: (f64, f64)) -> (u64, u64) {
        (p.0.to_bits(), p.1.to_bits())
    }

    #[test]
    fn the_folded_constants_match_the_engines_rodata() {
        assert_eq!(ONE_MINUS_1E11.to_bits(), 0x3fef_ffff_fffe_a028); // 0x4a018
        assert_eq!(ONE_PLUS_1E11.to_bits(), 0x3ff0_0000_0000_afec); // 0x4a010
        assert_eq!(TWO_MINUS_1E11.to_bits(), 0x3fff_ffff_ffff_5014); // 0x4a258
        assert_eq!(ONE_MINUS_1E12.to_bits(), 0x3fef_ffff_ffff_dcd1); // 0x4a398
        assert_eq!(INTERSECTION_EPSILON.to_bits(), 0x3d71_9799_812d_ea11); // 0x4a390
    }

    #[test]
    fn sgn_is_zero_at_zero_as_the_ec_compiler_expands_it() {
        assert_eq!(sgn(0.0), 0.0);
        assert_eq!(sgn(-0.0), 0.0);
        assert_eq!(sgn(1e-300), 1.0);
        assert_eq!(sgn(-1e-300), -1.0);
    }

    #[test]
    fn cross5x6_interruption_reproduces_the_ecs_own_crossing_vectors() {
        // `cross5x6InterruptionTest` (ri5x6.ec:1166-1245), compiled only under
        // `TEST_CROSSING`, accepts an answer within 1e-8 of its expected point; the
        // engine's answer is also required here, to the last bit.
        // (src, south, left, the eC's expected dst, the engine's dst)
        const CROSS_EC_VECTORS: [((f64, f64), bool, bool, (f64, f64), (f64, f64)); 48] = [
            ((1.0, 0.3), false, false, (1.7, 1.0), (1.7, 1.0)),
            ((2.0, 1.3), false, false, (2.7, 2.0), (2.7, 2.0)),
            ((3.0, 2.3), false, false, (3.7, 3.0), (3.7, 3.0)),
            ((4.0, 3.3), false, false, (4.7, 4.0), (4.7, 4.0)),
            (
                (1.7, 1.0),
                false,
                true,
                (1.0, 0.3),
                (1.0, 0.30000000000000004),
            ),
            (
                (2.7, 2.0),
                false,
                true,
                (2.0, 1.3),
                (2.0, 1.2999999999999998),
            ),
            ((3.7, 3.0), false, true, (3.0, 2.3), (3.0, 2.3)),
            ((4.7, 4.0), false, true, (4.0, 3.3), (4.0, 3.3)),
            ((0.3, 2.0), true, false, (1.0, 2.7), (1.0, 2.7)),
            ((1.3, 3.0), true, false, (2.0, 3.7), (2.0, 3.7)),
            ((2.3, 4.0), true, false, (3.0, 4.7), (3.0, 4.7)),
            ((3.3, 5.0), true, false, (4.0, 5.7), (4.0, 5.7)),
            (
                (1.0, 2.7),
                true,
                true,
                (0.3, 2.0),
                (0.2999999999999998, 2.0),
            ),
            (
                (2.0, 3.7),
                true,
                true,
                (1.3, 3.0),
                (1.2999999999999998, 3.0),
            ),
            ((3.0, 4.7), true, true, (2.3, 4.0), (2.3, 4.0)),
            ((4.0, 5.7), true, true, (3.3, 5.0), (3.3, 5.0)),
            ((1.0, 0.0), false, false, (2.0, 1.0), (2.0, 1.0)),
            ((2.0, 1.0), false, false, (3.0, 2.0), (3.0, 2.0)),
            ((3.0, 2.0), false, false, (4.0, 3.0), (4.0, 3.0)),
            ((4.0, 3.0), false, false, (5.0, 4.0), (5.0, 4.0)),
            ((1.0, 1.0), false, false, (1.0, 1.0), (1.0, 1.0)),
            ((2.0, 2.0), false, false, (2.0, 2.0), (2.0, 2.0)),
            ((3.0, 3.0), false, false, (3.0, 3.0), (3.0, 3.0)),
            ((4.0, 4.0), false, false, (4.0, 4.0), (4.0, 4.0)),
            ((2.0, 1.0), false, true, (1.0, 0.0), (1.0, 0.0)),
            ((3.0, 2.0), false, true, (2.0, 1.0), (2.0, 1.0)),
            ((4.0, 3.0), false, true, (3.0, 2.0), (3.0, 2.0)),
            ((5.0, 4.0), false, true, (4.0, 3.0), (4.0, 3.0)),
            ((1.0, 1.0), false, true, (1.0, 1.0), (1.0, 1.0)),
            ((2.0, 2.0), false, true, (2.0, 2.0), (2.0, 2.0)),
            ((3.0, 3.0), false, true, (3.0, 3.0), (3.0, 3.0)),
            ((4.0, 4.0), false, true, (4.0, 4.0), (4.0, 4.0)),
            ((0.0, 2.0), true, false, (1.0, 3.0), (1.0, 3.0)),
            ((1.0, 3.0), true, false, (2.0, 4.0), (2.0, 4.0)),
            ((2.0, 4.0), true, false, (3.0, 5.0), (3.0, 5.0)),
            ((3.0, 5.0), true, false, (4.0, 6.0), (4.0, 6.0)),
            ((1.0, 2.0), true, false, (1.0, 2.0), (1.0, 2.0)),
            ((2.0, 3.0), true, false, (2.0, 3.0), (2.0, 3.0)),
            ((3.0, 4.0), true, false, (3.0, 4.0), (3.0, 4.0)),
            ((4.0, 5.0), true, false, (4.0, 5.0), (4.0, 5.0)),
            ((1.0, 3.0), true, true, (0.0, 2.0), (0.0, 2.0)),
            ((2.0, 4.0), true, true, (1.0, 3.0), (1.0, 3.0)),
            ((3.0, 5.0), true, true, (2.0, 4.0), (2.0, 4.0)),
            ((4.0, 6.0), true, true, (3.0, 5.0), (3.0, 5.0)),
            ((1.0, 2.0), true, true, (1.0, 2.0), (1.0, 2.0)),
            ((2.0, 3.0), true, true, (2.0, 3.0), (2.0, 3.0)),
            ((3.0, 4.0), true, true, (3.0, 4.0), (3.0, 4.0)),
            ((4.0, 5.0), true, true, (4.0, 5.0), (4.0, 5.0)),
        ];
        for &(src, south, left, expected, engine) in &CROSS_EC_VECTORS {
            let got = cross5x6_interruption(src.0, src.1, south, left);
            assert!(
                (got.0 - expected.0).abs() <= 1e-8 && (got.1 - expected.1).abs() <= 1e-8,
                "{src:?} south {south} left {left}: {got:?}, the eC expects {expected:?}"
            );
            assert_eq!(bits(got), bits(engine), "{src:?} south {south} left {left}");
        }
    }

    #[test]
    fn cross5x6_interruption_follows_the_compiled_order() {
        // One or two points for each of the four arms at which the eC's own grouping
        // gives another double than the engine; the comment on each is that other answer.
        // (src, south, left, the engine's dst): the eC's own order answers otherwise at each
        const CROSS_COMPILED_VECTORS: &[((f64, f64), bool, bool, (f64, f64))] = &[
            (
                (-1.0000000000099998, -1.0000000000100002),
                false,
                false,
                (1.0000000000100002, -1.0000000000099998),
            ), // eC order: (4.00000000001, 3.99999999999)
            (
                (-1.0000000000099998, 1.3967260141776432),
                false,
                false,
                (3.603273985822357, 3.99999999999),
            ), // eC order: (1.6032739858223568, 3.99999999999)
            (
                (-0.7932577379529034, 2.3787405476839245),
                false,
                true,
                (2.3787405476839245, 4.793257737952904),
            ), // eC order: (2.3787405476839245, 4.793257737952903)
            (
                (-0.12335382376281345, 4.155512865853154),
                false,
                true,
                (4.155512865853154, 8.123353823762814),
            ), // eC order: (4.155512865853154, 8.123353823762812)
            (
                (0.99999999999, 0.9999999999900001),
                true,
                false,
                (-9.999889805101247e-12, 2.00000000001),
            ), // eC order: (4.99999999999, 5.00000000001)
            (
                (0.9999999999900003, 0.9999999999900001),
                true,
                false,
                (-9.999889805101247e-12, 2.0000000000099996),
            ), // eC order: (4.99999999999, 5.000000000009999)
            (
                (7.999999999989999, 7.999999999990003),
                true,
                true,
                (4.000000000009997, 3.999999999989999),
            ), // eC order: (2.0000000000099973, 3.999999999989999)
            (
                (3.9999999999899996, 5.036883343107645),
                true,
                true,
                (3.963116656892355, 4.999999999989999),
            ), // eC order: (1.963116656892355, 4.999999999989999)
        ];
        for &(src, south, left, engine) in CROSS_COMPILED_VECTORS {
            let got = cross5x6_interruption(src.0, src.1, south, left);
            assert_eq!(bits(got), bits(engine), "{src:?} south {south} left {left}");
        }
    }

    #[test]
    fn canonicalize5x6_matches_the_engine() {
        // (label, point, the engine's canonicalize5x6)
        const CANON_VECTORS: &[(&str, (f64, f64), (f64, f64))] = &[
            ("plain", (2.3, 3.1), (2.3, 3.1)),
            (
                "first fold",
                (6.2, 5.4),
                (1.2000000000000002, 0.40000000000000036),
            ),
            ("second fold", (-0.25, 0.6), (4.75, 5.6)),
            ("second fold on y", (0.6, -0.25), (5.6, 4.75)),
            ("north pole", (1.0, 0.0), (1.0, 0.0)),
            ("north pole image", (3.0, 2.0), (1.0, 0.0)),
            (
                "north pole within 1e-11",
                (4.000000000005, 2.999999999995),
                (1.0, 0.0),
            ),
            ("south pole", (4.0, 6.0), (4.0, 6.0)),
            ("south pole image", (2.0, 4.0), (4.0, 6.0)),
            ("north interruption", (2.0, 1.4), (2.6, 2.0)),
            (
                "north interruption, x just below",
                (2.9999999999999987, 2.7),
                (3.3, 2.9999999999999987),
            ),
            ("north interruption at x = 5", (5.0, 4.5), (0.5, 0.0)),
            ("south interruption", (1.6, 3.0), (2.0, 3.4)),
            (
                "south interruption, y just above",
                (0.3, 2.000000000000001),
                (1.0000000000000009, 2.7),
            ),
            (
                "east edge",
                (4.999999999999999, 3.2),
                (0.0, -1.7999999999999998),
            ),
            ("first fold then interruption", (7.0, 6.3), (2.7, 2.0)),
            (
                "second fold then south interruption",
                (-0.4, 1.0),
                (0.0, 1.4000000000000004),
            ),
        ];
        for &(label, p, engine) in CANON_VECTORS {
            assert_eq!(
                bits(canonicalize5x6(p.0, p.1)),
                bits(engine),
                "{label}: {p:?}"
            );
        }
    }

    #[test]
    fn intersects5x6_interruption_matches_the_engine() {
        // (label, a0, a1, b0, b1, the engine's answer: crossing, point, distance)
        const INTERSECTS_VECTORS: &[(&str, [f64; 8], Option<((f64, f64), f64)>)] = &[
            (
                "crossing",
                [0.5, 0.2, 1.5, 0.2, 1.0, 0.0, 1.0, 1.0],
                Some(((1.0, 0.2), 0.5)),
            ),
            (
                "crossing, oblique",
                [0.25, 0.1, 1.75, 0.9, 1.0, 0.0, 1.0, 1.0],
                Some(((1.0, 0.5), 0.5)),
            ),
            ("parallel", [0.5, 0.2, 0.5, 0.9, 1.0, 0.0, 1.0, 1.0], None),
            (
                "short of the segment",
                [0.5, 0.2, 0.9, 0.2, 1.0, 0.0, 1.0, 1.0],
                None,
            ),
            (
                "beyond the segment",
                [0.5, 1.2, 1.5, 1.2, 1.0, 0.0, 1.0, 1.0],
                None,
            ),
            (
                "south interruption",
                [0.3, 1.5, 0.9, 2.5, 0.0, 2.0, 1.0, 2.0],
                Some(((0.6000000000000001, 2.0), 0.5)),
            ),
            (
                "s at fl(1 - 1e-12)",
                [0.5, 0.999999999999, 1.5, 0.999999999999, 1.0, 0.0, 1.0, 1.0],
                Some(((1.0, 0.999999999999), 0.5)),
            ),
            (
                "s one double above fl(1 - 1e-12)",
                [
                    0.5,
                    0.9999999999990001,
                    1.5,
                    0.9999999999990001,
                    1.0,
                    0.0,
                    1.0,
                    1.0,
                ],
                None,
            ),
            (
                "s at 1e-12",
                [0.5, 1e-12, 1.5, 1e-12, 1.0, 0.0, 1.0, 1.0],
                Some(((1.0, 1e-12), 0.5)),
            ),
            (
                "s one double below 1e-12",
                [
                    0.5,
                    9.999999999999998e-13,
                    1.5,
                    9.999999999999998e-13,
                    1.0,
                    0.0,
                    1.0,
                    1.0,
                ],
                None,
            ),
        ];
        for &(label, a, engine) in INTERSECTS_VECTORS {
            let got =
                intersects5x6_interruption((a[0], a[1]), (a[2], a[3]), (a[4], a[5]), (a[6], a[7]));
            match (got, engine) {
                (None, None) => {}
                (Some((i, t)), Some((ei, et))) => {
                    assert_eq!(bits(i), bits(ei), "{label}");
                    assert_eq!(t.to_bits(), et.to_bits(), "{label}");
                }
                _ => panic!("{label}: {got:?}, the engine answers {engine:?}"),
            }
        }
    }

    #[test]
    fn crosses5x6_interruption_matches_the_engine() {
        // (label, c, dx, dy, the engine's answer: i_src, i_dst, north)
        const CROSSES_VECTORS: &[(&str, [f64; 4], Option<((f64, f64), (f64, f64), bool)>)] = &[
            (
                "within one cell: the early return",
                [2.3, 3.1, 0.1, 0.1],
                None,
            ),
            (
                "north, crossing right",
                [0.5, 0.2, 1.0, 0.0],
                Some(((1.0, 0.2), (1.8, 1.0), true)),
            ),
            (
                "north, crossing left",
                [1.6, 1.3, 0.0, -0.6],
                Some(((1.6, 1.0), (1.0, 0.3999999999999999), true)),
            ),
            (
                "south, crossing right",
                [0.3, 1.5, 0.6, 1.0],
                Some(((0.5999999999999999, 2.0), (1.0, 2.4000000000000004), false)),
            ),
            (
                "south, crossing left",
                [1.4, 2.6, -0.8, -0.2],
                Some(((1.0, 2.5), (0.5, 2.0), false)),
            ),
            (
                "the start snapped to an integer first",
                [1.9999999999999996, 1.5, -0.5, -0.8],
                Some(((1.6875, 1.0), (1.0, 0.3125), true)),
            ),
            (
                "two interruptions, the nearer taken",
                [0.5, 0.5, 2.0, 1.0],
                Some(((1.0, 0.75), (1.25, 1.0), true)),
            ),
            (
                "a cell boundary but no interruption",
                [1.5, 1.5, 0.0, 1.0],
                None,
            ),
            (
                "from the far side, wrapping",
                [4.5, 4.2, 1.0, 0.1],
                Some(((5.0, 4.25), (0.75, 0.0), true)),
            ),
        ];
        for &(label, a, engine) in CROSSES_VECTORS {
            let got = crosses5x6_interruption(a[0], a[1], a[2], a[3]);
            match (got, engine) {
                (None, None) => {}
                (Some(c), Some((i_src, i_dst, north))) => {
                    assert_eq!(bits(c.i_src), bits(i_src), "{label}");
                    assert_eq!(bits(c.i_dst), bits(i_dst), "{label}");
                    assert_eq!(c.in_north, north, "{label}");
                }
                _ => panic!("{label}: {got:?}, the engine answers {engine:?}"),
            }
        }
    }

    #[test]
    fn move5x6_vertex_matches_the_engine() {
        // The first three take a zero offset, as aperture 3 does and aperture 7 never
        // does; the comment on each gives the answer of this crate's earlier aperture-7
        // copy, which took `Sgn(0)` as 1 and floored where the eC truncates.
        // (label, c, dx, dy, the engine's move5x6Vertex)
        const MOVE_VECTORS: &[(&str, [f64; 4], (f64, f64))] = &[
            ("zero dy: Sgn(0) is 0", [1.0, 1.0, 0.2, 0.0], (1.2, 1.0)), // before: (1.0, 0.8)
            (
                "zero dx at a negative abscissa: truncation",
                [-0.5, 0.2, 0.0, 1.0],
                (-1.2, 0.0),
            ), // before: (0.0, 0.5)
            (
                "zero dy on an integer ordinate, the aperture-3 vertex",
                [2.3333333333333335, 2.0, 0.2222222222222222, 0.0],
                (2.555555555555556, 2.0),
            ), // before: (2.0, 1.4444444444444442)
            (
                "an aperture-7 sub-hexagon offset",
                [1.3, 1.2, -0.02040816326530612, -0.061224489795918366],
                (1.279591836734694, 1.1387755102040815),
            ),
            (
                "crossing: bottom-left",
                [
                    4.0,
                    5.9198302990490905,
                    -0.037037037037037035,
                    0.012345679012345678,
                ],
                (3.0307869849015265, 4.962962962962963),
            ),
            (
                "crossing: top-right",
                [0.7161147962813589, 0.0, 1.0, 0.3333333333333333],
                (0.6666666666666667, 1.0),
            ),
            (
                "crossing: top-left",
                [
                    1.5736969044751,
                    1.0,
                    -0.3333333333333333,
                    -0.6666666666666666,
                ],
                (0.33333333333333337, 0.09296976219156672),
            ),
            (
                "crossing: bottom-right",
                [
                    4.0,
                    5.714462689120146,
                    0.3333333333333333,
                    0.3333333333333333,
                ],
                (4.333333333333333, 5.0),
            ),
        ];
        for &(label, a, engine) in MOVE_VECTORS {
            assert_eq!(
                bits(move5x6_vertex(a[0], a[1], a[2], a[3])),
                bits(engine),
                "{label}"
            );
        }
    }

    #[test]
    fn move5x6_vertex2_matches_the_engine() {
        // One case for each arm, with and without the early crossing. The first pins how
        // the eC assigns a compound literal to the point it reads: the new centre of the
        // early crossing, `c = { cx, cy - (c.x - cx) }`, is computed wholly from the old
        // one, since assigning `c.x` before reading it would give (1.9, 1.05) here.
        // (label, c, dx, dy, cross_early, the engine's move5x6Vertex2)
        const MOVE2_VECTORS: &[(&str, [f64; 4], bool, (f64, f64))] = &[
            (
                "early crossing of the top dent to the left, all of the new centre read from the old",
                [2.25, 1.0, -0.1, 0.05],
                true,
                (1.25, 0.0),
            ),
            (
                "bottom-left",
                [4.0, 5.275866396957509, -0.07407407407407407, 0.0],
                false,
                (3.7241336030424907, 5.0),
            ),
            (
                "bottom-right",
                [
                    0.0,
                    1.5402500200404046,
                    0.6666666666666666,
                    0.6666666666666666,
                ],
                false,
                (1.4597499799595954, 3.0),
            ),
            (
                "cross-none",
                [1.0, 0.5361699555039579, 1.0, -0.3333333333333333],
                false,
                (2.0, 0.20283662217062454),
            ),
            (
                "north-pole",
                [1.9346699542210173, 0.9346699542210173, -1.0, -1.0],
                false,
                (1.0, 0.0),
            ),
            (
                "plain",
                [
                    4.0,
                    3.5774956669561826,
                    0.037037037037037035,
                    0.07407407407407407,
                ],
                false,
                (4.037037037037037, 3.6515697410302566),
            ),
            (
                "south-pole",
                [2.6057015974562416, 4.605701597456242, -2.0, -2.0],
                false,
                (4.0, 6.0),
            ),
            (
                "top-left",
                [
                    1.599854830593052,
                    0.599854830593052,
                    0.3333333333333333,
                    -0.6666666666666666,
                ],
                false,
                (0.400145169406948, -0.599854830593052),
            ),
            (
                "top-right",
                [1.9719175970658638, 3.0, 3.0, -1.0],
                false,
                (3.0, 4.028082402934136),
            ),
            (
                "early-bl bottom-left",
                [2.0, 3.6271958759386225, -2.0, -1.0],
                true,
                (1.0, 2.6271958759386225),
            ),
            (
                "early-bl plain",
                [
                    1.0,
                    2.0426756586388746,
                    0.037037037037037035,
                    -0.07407407407407407,
                ],
                true,
                (0.9943613783981624, 1.925925925925926),
            ),
            (
                "early-bl wrap+ bottom-left",
                [0.0, 1.5458156032271617, -1.0, -0.6666666666666666],
                true,
                (4.0, 5.545815603227162),
            ),
            (
                "early-br bottom-right",
                [
                    2.0835120187573044,
                    4.0,
                    0.6666666666666666,
                    0.6666666666666666,
                ],
                true,
                (4.083512018757304, 6.0),
            ),
            (
                "early-br plain",
                [2.8170705454377507, 4.0, 2.0, 1.0],
                true,
                (5.0, 5.182929454562249),
            ),
            (
                "early-tl plain",
                [
                    2.4582421390105247,
                    2.0,
                    -0.3333333333333333,
                    -0.1111111111111111,
                ],
                true,
                (1.6666666666666667, 1.4306467498783642),
            ),
            (
                "early-tl top-left",
                [2.803022762706295, 2.0, 3.0, -2.0],
                true,
                (1.803022762706295, 1.0),
            ),
            (
                "early-tr plain",
                [
                    3.0,
                    2.305241871895376,
                    0.1111111111111111,
                    0.3333333333333333,
                ],
                true,
                (3.805869239215735, 3.3333333333333335),
            ),
            (
                "early-tr top-right",
                [3.0, 2.382428138617534, 1.0, 0.6666666666666666],
                true,
                (5.0, 4.382428138617534),
            ),
            (
                "early-tr wrap- top-right",
                [
                    5.0,
                    4.018558129785283,
                    0.07407407407407407,
                    -0.1111111111111111,
                ],
                true,
                (2.0, 1.0185581297852826),
            ),
            (
                "wrap+ top-left",
                [
                    0.19665275671036564,
                    -0.8033472432896344,
                    0.3333333333333333,
                    -1.0,
                ],
                true,
                (4.803347243289634, 3.8033472432896342),
            ),
            (
                "wrap- bottom-right",
                [
                    5.1622450229351236,
                    -0.007919236379444694,
                    0.1111111111111111,
                    0.037037037037037035,
                ],
                true,
                (1.0079192363794451, -4.162245022935123),
            ),
        ];
        for &(label, a, early, engine) in MOVE2_VECTORS {
            let got = move5x6_vertex2(a[0], a[1], a[2], a[3], early);
            assert_eq!(bits(got), bits(engine), "{label}, cross_early {early}");
        }
    }

    #[test]
    fn move5x6_vertex3_follows_the_compiled_order() {
        // The comment on each crossing case gives the eC's own order's answer where it
        // differs from the engine's.
        // (label, c, dx, dy, the engine's move5x6Vertex3)
        const MOVE3_VECTORS: &[(&str, [f64; 4], (f64, f64))] = &[
            (
                "no crossing",
                [
                    0.9306445648387689,
                    0.9736714625397489,
                    6.774035123372114e-05,
                    1.6935087808430286e-05,
                ],
                (0.9307123051900026, 0.9736883976275573),
            ),
            (
                "no crossing, wrapped",
                [
                    5.0,
                    6.732837617648027,
                    1.6935087808430286e-05,
                    -5.080526342529086e-05,
                ],
                (1.6935087808533922e-05, 1.7327868123846013),
            ),
            (
                "north",
                [3.0, 2.0357656412939136, 1.0, 1.6666666666666667],
                (1.5951353840788394, 2.8429187695526963),
            ),
            (
                "north, the eC order differs",
                [4.0, 4.167828519619697, 2.0, 0.3333333333333333],
                (0.3321714803803033, 1.0),
            ), // eC order: (0.33217148038030275, 1.0)
            (
                "north, the eC order differs",
                [
                    3.8019405415132965,
                    3.0466418341183297,
                    0.0,
                    -0.06172839506172839,
                ],
                (3.030173121886797, 2.1980594584867035),
            ), // eC order: (3.0301731218867975, 2.1980594584867035)
            (
                "north, the eC order differs",
                [
                    3.0,
                    2.736955713020523,
                    -0.6666666666666666,
                    -1.3333333333333333,
                ],
                (3.19275524062562, 1.0702890463538561),
            ), // eC order: (3.1927552406256203, 1.0702890463538561)
            (
                "south",
                [
                    1.0,
                    2.0803531357654292,
                    -0.3333333333333333,
                    -0.5555555555555556,
                ],
                (0.5247975802098737, 1.4779689777111065),
            ),
            (
                "south, the eC order differs",
                [
                    1.9918212236384822,
                    3.246883404828425,
                    0.024691358024691357,
                    0.00411522633744856,
                ],
                (1.7545055627218513, 3.0330251633263474),
            ), // eC order: (1.7545055627218515, 3.0330251633263474)
            (
                "south, the eC order differs",
                [
                    1.8024132824512558,
                    1.888344775676165,
                    -1.6666666666666667,
                    1.6666666666666667,
                ],
                (1.1734953260879903, 0.2714932315691776),
            ), // eC order: (1.1734953260879901, 0.27149323156917804)
            (
                "south, the eC order differs",
                [
                    4.0,
                    5.954158674408842,
                    0.04938271604938271,
                    0.06172839506172839,
                ],
                (0.015887069470570836, 1.988746250679986),
            ), // eC order: (0.01588706947057056, 1.9887462506799856)
        ];
        for &(label, a, engine) in MOVE3_VECTORS {
            assert_eq!(
                bits(move5x6_vertex3(a[0], a[1], a[2], a[3])),
                bits(engine),
                "{label}"
            );
        }
    }
}
