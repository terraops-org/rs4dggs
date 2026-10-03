//! The helpers of the 5x6 planar layout, shared by the aperture-7 and aperture-3 topologies.
//!
//! Each ports a function of `projections/ri5x6.ec` in DGGAL v0.0.6. Three of them,
//! `move5x6Vertex`, `move5x6Vertex2` and `crosses5x6Interruption`, carry the eC's
//! `__attribute__ ((optimize("-fno-unsafe-math-optimizations")))`, and the library performs them
//! in the order of the eC's text, which is the order ported here. The other seven,
//! `intersects5x6Interruption`, `cross5x6Interruption`, `move5x6Vertex3`, `canonicalize5x6`,
//! `addIntermediatePoints`, `refine5x6` and `move5x6`,
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
//! The aperture-7 topology uses [`move5x6_vertex`], [`sgn`] and [`Crossing`], and both
//! topologies divide the edges of a zone's refined boundary with [`add_intermediate_points`];
//! [`move5x6`] is the walk on which the aperture-7 sub-zone order is built; the
//! others serve the
//! aperture-3 topology, whose vertices, neighbours, sub-zones and refined boundary call them.
//! The aperture-7
//! quantiser's domain fold, which the eC performs inline in `I7HZone::fromCentroid`, stays in
//! `topologies/hex_a7.rs`: it bears the name of the eC function ported here as
//! [`canonicalize5x6`], but it is a different routine.

/// `fl(1 - 1e-11)`, a constant gcc folded into `cross5x6Interruption` and into the crossings of
/// `move5x6`, as it did the next two; `.rodata` at `0x4a018`.
const ONE_MINUS_1E11: f64 = 1.0 - 1e-11;
/// `fl(1 + 1e-11)`, folded likewise; `.rodata` at `0x4a010`.
const ONE_PLUS_1E11: f64 = 1.0 + 1e-11;
/// `fl(2 - 1e-11)`, folded likewise; `.rodata` at `0x4a258`.
const TWO_MINUS_1E11: f64 = 2.0 - 1e-11;
/// `fl(1 - 1e-12)`, folded into `intersects5x6Interruption`; `.rodata` at `0x4a398`.
const ONE_MINUS_1E12: f64 = 1.0 - 1e-12;
/// The tolerance `e` of `intersects5x6Interruption` (ri5x6.ec:1130); `.rodata` at `0x4a390`.
const INTERSECTION_EPSILON: f64 = 1e-12;
/// `fl(5 + 1e-11)`, folded into `refine5x6`: the bound of the layout along `x`, and the term of
/// the floors of a `next` brought forward by a period; `.rodata` at `0x4a188`.
const FIVE_PLUS_1E11: f64 = 5.0 + 1e-11;
/// `fl(6 + 1e-11)`, folded into `refine5x6` and into `move5x6`: the bound of the layout along
/// `y`; `.rodata` at `0x4a2d0`.
const SIX_PLUS_1E11: f64 = 6.0 + 1e-11;
/// `fl(5 - 1e-11)`, folded into `refine5x6`, as the term of the floors of a `next` brought back
/// by a period, and into `move5x6`, as a bound of the layout; `.rodata` at `0x49f48`.
const FIVE_MINUS_1E11: f64 = 5.0 - 1e-11;
/// `fl(2e-11)`, to which gcc folded the `2 * EPS` of `2*Sgn(d) * EPS` in `move5x6`
/// (ri5x6.ec:1752, :1759); `.rodata` at `0x49fa8`.
const TWICE_1E11: f64 = 2e-11;

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

/// The eC's `(int)floor(v)`, as gcc inlines it: the truncation, less one where the truncation
/// exceeds `v` (in `fromCentroid` at `0xca99` to `0xcac6`; in `refine5x6` at `0x3aabf` to
/// `0x3aad6` and at each of its other floors). That is the floor for every `v` within the range
/// of `int`;
/// below it the subtraction wraps to `i32::MAX`. The eC writes it inline, as in `refine5x6`
/// (`ri5x6.ec:808-813`, `:886-887`) and, for the aperture-3 topology, in `fromCentroid`
/// (`RI3H.ec:1347-1356`) and in `getNeighbor` (`RI3H.ec:1008-1009`).
pub(crate) fn floor_i32(v: f64) -> i32 {
    let t = cvtt_i32(v);
    if f64::from(t) > v {
        t.wrapping_sub(1)
    } else {
        t
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

/// Adds to `points` the points of one edge of a zone's refined boundary in the 5x6 plane: the
/// vertex `p`, then the edge's divisions towards the next vertex `n`, which is itself left to the
/// next edge. Ports `addIntermediatePoints` (ri5x6.ec:651-726), compiled at `0x36e70` to
/// `0x379a2`, on which the refined boundaries of both topologies are built.
///
/// `n_divisions` is the eC's `nDivisions`, an `int`: zero counts as one, and a negative number
/// divides nothing. `interruption` is the eC's pair `i1` and `i2`, which are null together or
/// not at all: where the edge crosses an interruption of the layout, the point at which it
/// leaves one rhombus and the point at which it enters the next, in the frame of `p` and of `n`
/// respectively. The edge is then divided on either side of the interruption in proportion to
/// the two lengths, and the two points are added between the sides. `crs84` is the eC's flag of
/// the same name, true where the points are bound for WGS84 and false where they stay in the 5x6
/// CRS; it decides two things alone: an edge of one division shows its interruption only in the
/// 5x6 CRS, and an interruption at a pole's image is looked for only towards WGS84.
///
/// The pole rule (ri5x6.ec:655-676): an edge whose interruption lies within 1e-6 of one of the
/// four images of a pole, or whose two ends lie on either side of one along the layout's edge,
/// is divided twenty times as finely, since it is there that a straight edge of the plane is
/// furthest from its image on the sphere; and the walk then takes those divisions twenty at a
/// time, save within twenty of the interruption, where it takes each one.
pub(crate) fn add_intermediate_points(
    points: &mut Vec<(f64, f64)>,
    p: (f64, f64),
    n: (f64, f64),
    n_divisions: i32,
    interruption: Option<((f64, f64), (f64, f64))>,
    crs84: bool,
) {
    let (mut dx, mut dy) = (n.0 - p.0, n.1 - p.1);
    // The eC sets the flag from the interruption and then, in an `if` and three `else if`, from the
    // two ends (ri5x6.ec:655-671); every arm of the chain can only set it, so the chain is a
    // disjunction. The library compares as the eC writes, `fabs(i.x - 0.5) < 1E-6` and its
    // fellows by a subtraction, a mask and a `comisd` each (`0x36f2c` to `0x36f5b`, `0x375a0` to
    // `0x3795e`), and the ends against 0.5, 0.01, 4.99, 4.5, 1.5, 2.99, 2.01 and 3.5 from
    // `.rodata` (`0x36f61` to `0x36fae`, `0x372e5` to `0x37378`, `0x37470` to `0x37493`,
    // `0x3774d` to `0x377fb`).
    let at_a_pole = |i: (f64, f64)| {
        let near = |x: f64, y: f64| (i.0 - x).abs() < 1e-6 && (i.1 - y).abs() < 1e-6;
        near(0.5, 0.0) || near(5.0, 4.5) || near(2.0, 3.5) || near(1.5, 3.0)
    };
    let interruption_near_pole = (crs84
        && interruption.is_some_and(|(i1, i2)| at_a_pole(i1) || at_a_pole(i2)))
        || (n.0 < 0.5 && p.0 > 0.5 && p.1 < 0.01 && n.1 < 0.01)
        || (n.0 > 4.99 && p.0 > 4.99 && p.1 > 4.5 && n.1 < 4.5)
        || (n.0 > 1.5 && p.0 < 1.5 && p.1 > 2.99 && n.1 > 2.99)
        || (n.0 < 2.01 && p.0 < 2.01 && p.1 < 3.5 && n.1 > 3.5);

    // faithful: `nDivisions *= 20` (ri5x6.ec:676) is a product of the eC's 32-bit `int`, which
    // the library takes with a `lea` and a `shl` on `ebx` (`0x37398`, `0x3739b`), so it wraps;
    // the zero made one before it is the `cmove` at `0x36fc7`.
    let mut n_divisions = if n_divisions == 0 { 1 } else { n_divisions };
    if interruption_near_pole {
        n_divisions = n_divisions.wrapping_mul(20);
    }
    let divisions = f64::from(n_divisions);

    if dx < -3.0 {
        dx += 5.0;
        dy += 5.0;
    }
    if dx > 3.0 {
        dx -= 5.0;
        dy -= 5.0;
    }

    let Some((i1, i2)) = interruption else {
        if n_divisions == 1 {
            points.push(p);
            return;
        }
        // faithful: site AIP-3. The eC's `p.x + j * dx / nDivisions` (ri5x6.ec:724) is, in
        // DGGAL v0.0.6 as built with `-O2 -ffast-math` (`libdggal.so`, BuildID
        // `e75d6ab18f8b460713ebcd21df6f07922fb909c3`), `p.x + j * (dx * (1 / nDivisions))`: the
        // reciprocal is taken once (`divsd` at `0x373f3`), both offsets multiply it before the
        // loop (`mulpd` at `0x373fb`), and each turn multiplies the pair by `j` and adds `p`
        // (`0x37436`, `0x3743b`). The eC's order, a product by `j` and then a division, gives
        // another double at some point of 17,553 of the 49,695 rings sampled on IGEO7.
        let per_division = 1.0 / divisions;
        let step = (dx * per_division, dy * per_division);
        for j in 0..n_divisions {
            let j = f64::from(j);
            points.push((p.0 + j * step.0, p.1 + j * step.1));
        }
        return;
    };

    let pi1 = (i1.0 - p.0, i1.1 - p.1);
    let i2n = (n.0 - i2.0, n.1 - i2.1);
    let l1 = crate::math::sqrt(pi1.0 * pi1.0 + pi1.1 * pi1.1);
    let l2 = crate::math::sqrt(i2n.0 * i2n.0 + i2n.1 * i2n.1);
    let length = l1 + l2;
    // As written (ri5x6.ec:690): the product, then the quotient (`0x37068` to `0x37070`).
    let t = divisions * l1 / length;

    points.push(p);
    if n_divisions == 1 && !crs84 {
        points.push(i1);
        points.push(i2);
    }

    // faithful: `(int)t` (ri5x6.ec:706) is a 32-bit `cvttsd2si` (`0x370ee`), taken once before
    // the loop: where the two sides of the interruption have no length `t` is NaN, the
    // conversion gives the integer indefinite, which no `j` equals, and neither point of the
    // interruption is added. Rust's `as` would give 0 there, and both points at `j == 1`.
    let whole_t = cvtt_i32(t);
    // faithful: sites AIP-1 and AIP-2, in the same library. The eC's `p.x + j * pi1.x / t`
    // (ri5x6.ec:702-703) is `p.x + j * (pi1.x * (1 / t))`, and its
    // `i2.x + (j - t) * i2n.x / (nDivisions - t)` (ri5x6.ec:714-715) is
    // `i2.x + (j - t) * (i2n.x * (1 / (nDivisions - t)))`: both reciprocals are taken before
    // the loop (`divsd` at `0x370f6` and, on the difference formed at `0x370e4`, at `0x37105`),
    // the offsets multiply them there (`mulpd` at `0x3710d` and `0x3711a`), and each turn
    // multiplies the stored pair by `j` or by `j - t` and adds `p` or `i2` (`0x37223` and
    // `0x37229`, `0x37185` and `0x3718b`). Either put back in the eC's order moves some three
    // hundred of the 49,695 rings sampled on IGEO7.
    let before = {
        let per_t = 1.0 / t;
        (pi1.0 * per_t, pi1.1 * per_t)
    };
    let after = {
        let per_rest = 1.0 / (divisions - t);
        (i2n.0 * per_rest, i2n.1 * per_rest)
    };
    let mut j: i32 = 1;
    while j < n_divisions {
        let at = f64::from(j);
        if at < t {
            points.push((p.0 + at * before.0, p.1 + at * before.1));
        }
        if j == whole_t || (j == 1 && whole_t == 0) {
            points.push(i1);
            points.push(i2);
        }
        if at > t {
            points.push((i2.0 + (at - t) * after.0, i2.1 + (at - t) * after.1));
        }
        // The eC's stride (ri5x6.ec:698), on the eC's `int`, which the library adds without a
        // check (`0x371ae` to `0x371ca`): twenty where the pole rule holds and `j` is twenty or
        // more from `t`, one otherwise, a NaN `t` among the otherwise.
        let stride = if interruption_near_pole && (at - t).abs() >= 20.0 {
            20
        } else {
            1
        };
        j = j.wrapping_add(stride);
    }
}

/// Refines the outline `src` of a zone in the 5x6 plane: for each side, from a vertex to the
/// next and from the last back to the first, the vertex and the divisions of the side, the next
/// vertex being left to the next side. Ports `refine5x6` (ri5x6.ec:796-1027; its lines 1028 to
/// 1121 lie under `#if 0`), compiled at `0x3a7f0` to `0x3ba2f`, on which the refined boundary of
/// the aperture-3 topology is built.
///
/// It is ported with the eC's `wrap` true, the one form in which a ring bound for WGS84 asks it
/// (`RI3H.ec:507`): a vertex that lies on an interruption of the layout, with its neighbour
/// beyond, is carried across before the side is divided; a side that leaves the layout by one
/// end and re-enters by the other is divided as one; and a side of no length is left out.
/// `n_divisions` is the eC's `nDivisions`, an `int`, handed to [`add_intermediate_points`] as it
/// stands. A side that crosses an interruption between its ends is given the two points of the
/// crossing, in the frame of each end, which [`add_intermediate_points`] places between its
/// divisions.
///
/// The eC warns that this logic "is likely not as correct as `crosses5x6Interruption()` and
/// `rotate5x6Offset()`" (ri5x6.ec:804); it is the engine's all the same, and it is ported as it
/// is compiled.
#[allow(
    clippy::too_many_lines,
    reason = "one eC function, kept whole so that it reads against its source line by line"
)]
pub(crate) fn refine5x6(src: &[(f64, f64)], n_divisions: i32) -> Vec<(f64, f64)> {
    // faithful: DGGAL v0.0.6 as built with `-O2 -ffast-math` (`libdggal.so`, BuildID
    // `e75d6ab18f8b460713ebcd21df6f07922fb909c3`). `refine5x6` carries no attribute, so the
    // library decides its order, and eighteen of its expressions are re-associated there: the
    // sites R56-1 to R56-18 below, each with its address. Six of them move a ring of the
    // engine when put back in the eC's order (R56-1, -4, -5, -8, -15 and -17) and are pinned
    // by rings read from the engine; the other twelve moved none of the 264,375 rings
    // compared with the engine, and follow the library all the same. The comparisons and the floors
    // are the eC's own.
    const E: f64 = 1e-11;
    let count = src.len();
    // The eC's `minAllocSize` (ri5x6.ec:798-799), which is a hint and no more.
    let mut points: Vec<(f64, f64)> =
        Vec::with_capacity((count + 2) * usize::try_from(n_divisions).unwrap_or(0));
    // The eC's `int` arithmetic, which the library performs in 32 bits without a check.
    let whole = f64::from;
    let above = |c: i32| f64::from(c.wrapping_add(1));

    for i in 0..count {
        let mut p = src[i];
        let mut next = src[if i < count - 1 { i + 1 } else { 0 }];
        // ri5x6.ec:808-813, the eight floors at `0x3aabf` to `0x3abc9`.
        let (cpx1, cpy1) = (floor_i32(p.0 + E), floor_i32(p.1 + E));
        let (cnx1, cny1) = (floor_i32(next.0 + E), floor_i32(next.1 + E));
        let (cpx2, cpy2) = (floor_i32(p.0 - E), floor_i32(p.1 - E));
        let (cnx2, cny2) = (floor_i32(next.0 - E), floor_i32(next.1 - E));
        // ri5x6.ec:814-817, at `0x3abbf` to `0x3ac17` and `0x3a870` to `0x3a924`. The third of
        // the four has no epsilon, in the eC as in the library (`0x3a8c2`); given one, it
        // moves no ring of some 87,000 zones a grid.
        let n_top_right_of_p = (next.0 > p.0 + E && next.0 - p.0 < 3.0) || p.0 - next.0 > 3.0;
        let n_top_left_of_p = (next.0 < p.0 - E && p.0 - next.0 < 3.0) || next.0 - p.0 > 3.0;
        let n_bottom_right_of_p = (next.1 > p.1 && next.1 - p.1 < 3.0) || p.1 - next.1 > 3.0;
        let n_bottom_left_of_p = (next.1 < p.1 - E && p.1 - next.1 < 3.0) || next.1 - p.1 > 3.0;
        // ri5x6.ec:818-825, at `0x3a927` to `0x3aa01` and the blocks it branches to.
        let at_top_dent_crossing_right = cpx2 != cpx1 && p.0 > p.1 && n_top_right_of_p;
        let at_top_dent_crossing_left = cpy2 != cpy1 && p.0 > p.1 && n_top_left_of_p;
        let at_bottom_dent_crossing_left = cpx2 != cpx1 && p.1 > p.0 + 1.0 && n_bottom_left_of_p;
        let at_bottom_dent_crossing_right = cpy2 != cpy1 && p.1 > p.0 + 1.0 && n_bottom_right_of_p;
        let next_at_top_dent_crossing_right = cnx2 != cnx1 && next.0 > next.1 && n_top_left_of_p;
        let next_at_top_dent_crossing_left = cny2 != cny1 && next.0 > next.1 && n_top_right_of_p;
        let next_at_bottom_dent_crossing_left =
            cnx2 != cnx1 && next.1 > next.0 + 1.0 && n_bottom_right_of_p;
        let next_at_bottom_dent_crossing_right =
            cny2 != cny1 && next.1 > next.0 + 1.0 && n_bottom_left_of_p;

        if (at_top_dent_crossing_right && next_at_top_dent_crossing_left)
            || (at_top_dent_crossing_left && next_at_top_dent_crossing_right)
            || (at_bottom_dent_crossing_right && next_at_bottom_dent_crossing_left)
            || (at_bottom_dent_crossing_left && next_at_bottom_dent_crossing_right)
        {
            // A point on each side of an interruption: both are added as they stand
            // (ri5x6.ec:834-843, `0x3aa07` to `0x3aa51`).
            points.push(p);
            points.push(next);
            continue;
        }

        // ri5x6.ec:852, at `0x3ac5d` to `0x3ac88`: a side of no length.
        if (next.0 - p.0).abs() < E && (next.1 - p.1).abs() < E {
            continue;
        }

        // The eC crosses at once where the crossing is not within the side (ri5x6.ec:858-865).
        // faithful: sites R56-1 to R56-4. The library gathers the whole numbers of each
        // expression into one sum before the coordinate enters it:
        //   R56-1, `cpx1 + 1.0 - (p.y - cpy1)` (:859), is `(cpx1 + cpy1) + (1.0 - p.y)`
        //          (`0x3acff` to `0x3ad22`);
        //   R56-2, `cpy1 - (p.x - cpx1)` (:861), is `(cpy1 + cpx1) - p.x` (`0x3b350` to
        //          `0x3b364`);
        //   R56-3, `cpx1 - (p.y - cpy1)` (:863), is `(cpx1 + cpy1) - p.y` (`0x3b308` to
        //          `0x3b31c`);
        //   R56-4, `cpy1 + (cpx1 + 1 - p.x)` (:865), is `(cpy1 - p.x) + (cpx1 + 1)` (`0x3b6c7`
        //          to `0x3b6de`).
        // R56-1 and R56-4 move rings of the engine; R56-2 and R56-3 moved none.
        if at_top_dent_crossing_right {
            p = ((whole(cpx1) + whole(cpy1)) + (1.0 - p.1), above(cpy1));
        } else if at_top_dent_crossing_left {
            p = (whole(cpx1), (whole(cpy1) + whole(cpx1)) - p.0);
        } else if at_bottom_dent_crossing_left {
            p = ((whole(cpx1) + whole(cpy1)) - p.1, whole(cpy1));
        } else if at_bottom_dent_crossing_right {
            p = (above(cpx1), (whole(cpy1) - p.0) + above(cpx1));
        }

        // faithful: sites R56-5 to R56-8, the same four forms on `next` (ri5x6.ec:868-874):
        // R56-5 at `0x3ad43` to `0x3ad6e`, R56-6 at `0x3b4c0` to `0x3b4d8`, R56-7 at `0x3b387`
        // to `0x3b39b`, R56-8 at `0x3b691` to `0x3b6ac`. R56-5 and R56-8 move rings of the
        // engine; R56-6 and R56-7 moved none.
        if next_at_top_dent_crossing_right {
            next = ((whole(cnx1) + whole(cny1)) + (1.0 - next.1), above(cny1));
        } else if next_at_top_dent_crossing_left {
            next = (whole(cnx1), (whole(cny1) + whole(cnx1)) - next.0);
        } else if next_at_bottom_dent_crossing_left {
            next = ((whole(cnx1) + whole(cny1)) - next.1, whole(cny1));
        } else if next_at_bottom_dent_crossing_right {
            next = (above(cnx1), (whole(cny1) - next.0) + above(cnx1));
        }

        // ri5x6.ec:876-879, at `0x3ad7f` to `0x3adc4` and `0x3b4f0` to `0x3b51c`: a `p` beyond
        // the layout is brought back by a period. The bounds are the folded `fl(5 + 1e-11)`
        // and `fl(6 + 1e-11)`.
        if p.0 > FIVE_PLUS_1E11 || p.1 > SIX_PLUS_1E11 {
            p = (p.0 - 5.0, p.1 - 5.0);
        } else if p.0 < -E || p.1 < -E {
            p = (p.0 + 5.0, p.1 + 5.0);
        }
        // ri5x6.ec:886, at `0x3ae10` to `0x3ae82`, `0x3b07e` to `0x3b0f1` and `0x3b862` to
        // `0x3b8e5`, one copy to each arm below.
        let (cpx, cpy) = (floor_i32(p.0 + E), floor_i32(p.1 + E));

        // ri5x6.ec:880-883 and :887, at `0x3adca` to `0x3adf4`: `next` likewise, and its floors.
        // faithful: sites R56-9 and R56-10. Where `next` is moved by a period, the eC floors
        // the moved coordinate plus the epsilon, `floor((next.x - 5) + e)`; the library floors
        // the coordinate as it was, less the folded `fl(5 - 1e-11)` (`0x3b0a2` to `0x3b0b2`,
        // floored at `0x3b0b9` and `0x3b0f5`) or plus the folded `fl(5 + 1e-11)` (`0x3ae14` and
        // `0x3ae21`, floored at `0x3ae4a` and `0x3b0f5`). Neither moved a ring of the sample;
        // both constants are pinned against the library's `.rodata`.
        let (cnx, cny);
        if next.0 > FIVE_PLUS_1E11 || next.1 > SIX_PLUS_1E11 {
            (cnx, cny) = (
                floor_i32(next.0 - FIVE_MINUS_1E11),
                floor_i32(next.1 - FIVE_MINUS_1E11),
            );
            next = (next.0 - 5.0, next.1 - 5.0);
        } else if next.0 < -E || next.1 < -E {
            (cnx, cny) = (
                floor_i32(next.0 + FIVE_PLUS_1E11),
                floor_i32(next.1 + FIVE_PLUS_1E11),
            );
            next = (next.0 + 5.0, next.1 + 5.0);
        } else {
            (cnx, cny) = (floor_i32(next.0 + E), floor_i32(next.1 + E));
        }

        // ri5x6.ec:890-892, at `0x3b0f5` to `0x3b106`, `0x3b866` to `0x3b89c` and `0x3b3b8` to
        // `0x3b3dd`.
        let (mut dx, mut dy) = (next.0 - p.0, next.1 - p.1);
        let n_south = cnx.wrapping_add(cny) & 1 != 0;
        let p_south = cpx.wrapping_add(cpy) & 1 != 0;
        // The eC's `pi1` and `pi2` where it sets `interrupted`: the point at which the side
        // leaves the rhombus of `p`, and the point at which it enters that of `next`.
        let mut interruption: Option<((f64, f64), (f64, f64))> = None;

        // faithful: sites R56-11 to R56-18, the second member of each interruption
        // (ri5x6.ec:913, :921, :951, :960, :979, :987, :995, :1003). The eC writes each as a
        // whole number less, or plus, a difference; the library gathers the whole numbers
        // first, and in the four arms that wrap, where the cells are known, folds them:
        //   R56-11, `cny + 1 - (pi1.x - cpx)`, is `(cpx + 5.0) - pi1.x` (`0x3b55c` to `0x3b575`);
        //   R56-12, `cnx + 1 - (pi1.y - cpy)`, is `6.0 - pi1.y` (`0x3b5e0` to `0x3b5f2`);
        //   R56-13, `cnx + 1 - (pi1.y - cpy)`, is `5.0 - pi1.y` (`0x3b66b`);
        //   R56-14, `cny + (cpx + 1 - pi1.x)`, is `((cpx + 1) + 1.0) - pi1.x` (`0x3b7aa` to
        //           `0x3b7b6`);
        //   R56-15, `cnx + (cpy + 1 - pi1.y)`, is `((cpy + 1) + cnx) - pi1.y` (`0x3b947` to
        //           `0x3b94f`);
        //   R56-16, `cny + 1 - (pi1.x - cpx)`, is `((cny + 1) + cpx) - pi1.x` (`0x3b845`,
        //           `0x3b46f` to `0x3b478`);
        //   R56-17, `cny + (cpx + 1 - pi1.x)`, is `((cpx + 1) + cny) - pi1.x` (`0x3b465` to
        //           `0x3b478`);
        //   R56-18, `cnx + 1 - (pi1.y - cpy)`, is `((cnx + 1) + cpy) - pi1.y` (`0x3b9dd` to
        //           `0x3b9f1`).
        // R56-15 and R56-17 move rings of the engine; the other six moved none. The first
        // member of each, `pi1`, is performed as the eC writes it.
        if cpx != cnx && cpy != cny {
            if cnx.wrapping_sub(cpx) > 3 {
                // ri5x6.ec:896-931, at `0x3b14d` to `0x3b1c0` and `0x3b528` to `0x3b602`.
                dx -= 5.0;
                dy -= 5.0;
                if cpy == 0 && cny == 5 && cnx == 4 {
                    // An uninterrupted wrap to the left. The eC computes `pi1` and `pi2` here
                    // (ri5x6.ec:904-905), and so does the library (`0x3b18a` to `0x3b1c0`), but
                    // `interrupted` stays false and neither is read again.
                } else if cpy == 0 && cny == 4 && cnx == 4 {
                    // Crossing the top dent to the left while wrapping. With `cpy` known to be
                    // zero the library takes `0.5 * (cpy - p.y)` as `p.y * -0.5` (`0x3b53e`),
                    // the same double.
                    let pi1 = (p.0 + p.1 * -0.5, 0.0);
                    interruption = Some((pi1, (5.0, (whole(cpx) + 5.0) - pi1.0)));
                } else if cpy == 1 && cny == 5 && cnx == 4 {
                    // Crossing the bottom dent to the left while wrapping.
                    let pi1 = (whole(cpx), p.1 + 0.5 * (whole(cpx) - p.0));
                    interruption = Some((pi1, (6.0 - pi1.1, 6.0)));
                } else {
                    // "This happens on exact same point": `wrapped` is set, so the side is left
                    // out (ri5x6.ec:929, `0x3b2b0`).
                    continue;
                }
            } else if cnx.wrapping_sub(cpx) < -3 {
                // ri5x6.ec:932-970, at `0x3b269` to `0x3b2a5`, `0x3b610` to `0x3b684`,
                // `0x3b6f0` to `0x3b738` and `0x3b74a` to `0x3b7cc`.
                dx += 5.0;
                dy += 5.0;
                if cpy == 5 && cny == 0 && cnx == 0 {
                    // An uninterrupted wrap to the right. The eC then takes from `dx` and `dy`
                    // the differences `pi2.x + 5 - pi1.x` and `pi2.y + 5 - pi1.y`
                    // (ri5x6.ec:942-943), of a `pi2` that is `pi1` less five: the library folds
                    // both to zero and subtracts nothing (`0x3b6f0` to `0x3b738`, where `dx`
                    // and `dy` reach the test below as they were).
                } else if cpy == 4 && cny == 0 && cnx == 0 {
                    // Crossing the top dent to the right while wrapping.
                    let pi1 = (above(cpx), p.1 + 0.5 * (above(cpx) - p.0));
                    interruption = Some((pi1, (5.0 - pi1.1, 0.0)));
                } else if cpy == 5 && cny == 1 && cnx == 0 {
                    // Crossing the bottom dent to the right while wrapping.
                    let pi1 = (p.0 + 0.5 * (6.0 - p.1), 6.0);
                    interruption = Some((pi1, (0.0, (above(cpx) + 1.0) - pi1.0)));
                } else {
                    continue;
                }
            } else if n_south == p_south {
                // ri5x6.ec:971-1013, at `0x3b3e5` to `0x3b489`, `0x3b7d1` to `0x3b849` and
                // `0x3b904` to `0x3ba29`. Where none of the four holds `wrapped` is not set,
                // and the side goes on as a plain one.
                if cnx > cpx && !n_south && !p_south {
                    // Crossing the top dent to the right.
                    let pi1 = (above(cpx), p.1 + 0.5 * (above(cpx) - p.0));
                    let pi2 = ((above(cpy) + whole(cnx)) - pi1.1, whole(cny));
                    interruption = Some((pi1, pi2));
                } else if cnx < cpx && !n_south && !p_south {
                    // Crossing the top dent to the left.
                    let pi1 = (p.0 + 0.5 * (whole(cpy) - p.1), whole(cpy));
                    let pi2 = (above(cnx), (above(cny) + whole(cpx)) - pi1.0);
                    interruption = Some((pi1, pi2));
                } else if cnx > cpx && n_south && p_south {
                    // Crossing the bottom dent to the right.
                    let pi1 = (p.0 + 0.5 * (above(cpy) - p.1), above(cpy));
                    let pi2 = (whole(cnx), (above(cpx) + whole(cny)) - pi1.0);
                    interruption = Some((pi1, pi2));
                } else if cnx < cpx && n_south && p_south {
                    // Crossing the bottom dent to the left.
                    let pi1 = (whole(cpx), p.1 + 0.5 * (whole(cpx) - p.0));
                    let pi2 = ((above(cnx) + whole(cpy)) - pi1.1, above(cny));
                    interruption = Some((pi1, pi2));
                }
            }
        }

        // ri5x6.ec:1016-1017, at `0x3b1d0` to `0x3b1f6` and its copies: no length is left.
        if dx.abs() < E && dy.abs() < E {
            continue;
        }
        // The eC then sets `dx` and `dy` anew where the side is interrupted (ri5x6.ec:1019-1024)
        // and reads neither again: `addIntermediatePoints` takes the two ends and forms its
        // own. The library forms nothing there (`0x3b4b0`, `0x3acb6` to `0x3acd0`).
        add_intermediate_points(&mut points, p, next, n_divisions, interruption, true);
    }
    points
}

/// The most turns of its loop that [`move5x6`] makes before it gives the walk up.
///
/// The eC's loop (ri5x6.ec:1736) has no bound, and the engine's function does not return from
/// some arguments that its own callers never give it. Measured over the calls that the
/// aperture-7 sub-zone generator makes on IGEO7 (for the order, the first sub-zone, the entry at
/// an index, the index of a sub-zone and the test for one; at every zone of levels 0 to 3 and a
/// sample of each level to 18, the pentagons among them; at depths 1 to 3 and, for the pentagons
/// and every fortieth zone, to 5: 1,022,074,348 calls), no walk takes more than four turns: four
/// in five take one, and one in a thousand takes four. The bound is sixteen times that maximum.
const MOVE5X6_MAX_TURNS: u32 = 64;

/// Walks from `o` by the offset `(dx, dy)` across the 5x6 layout and returns the point reached:
/// the walk goes rhombus by rhombus, and at each interruption it crosses it carries the point to
/// the far side and turns what is left of the offset by sixty degrees. Ports `move5x6`
/// (ri5x6.ec:1723-1931), compiled at `0x38d70` to `0x39880`, on which the sub-zone order of the
/// aperture-7 topology is built.
///
/// `None` is this crate's own answer, where the engine has none: the walk did not end within
/// [`MOVE5X6_MAX_TURNS`] turns of its loop. The engine's function does not return from some
/// arguments that its own callers never give it (a point 1e-12 off an edge with a long offset, a
/// corner of the layout, two turns a crossing from an arbitrary point, as in the call that the
/// eC's own test sets aside as endless, ri5x6.ec:1940-1941), and the crate does not mirror
/// that. Every call that the sub-zone generator makes for a valid zone ends well within the
/// bound, and returns `Some` of the engine's point. Where the answer is `None`, `adjust` is
/// left as the turns made until then have it, and means nothing. The bound counts the turns of
/// the loop, not the rotations made within each: a call costs at most [`MOVE5X6_MAX_TURNS`]
/// times `n_rotations` iterations, and the callers pass 1 or 2.
///
/// `n_rotations` is the eC's `nRotations`, an `int`: how many turns of sixty degrees a crossing
/// gives the offset, one everywhere save on the scanlines past the centre of a pentagon of an
/// even level, where the generator asks two (`RI7H.ec:3692`); none where it is zero or negative.
/// `adjust` is the eC's pair `adjX` and `adjY`, which it heeds only where neither is null and
/// which every caller passes together or not at all: a direction of the caller's own, turned in
/// place as the offset is turned, so that the caller goes on in the frame of the rhombus
/// reached. `final_cross` is the eC's `finalCross`: where it is true an offset that ends on an
/// interruption is carried across it, and where it is false the walk stops as soon as the offset
/// is spent. An offset shorter than 1e-11 in both coordinates moves nothing: the origin is
/// returned, brought forward by a period of five where it lies before the layout.
///
/// faithful: the function is unguarded, and DGGAL is built with `-O2 -ffast-math`. In
/// `libdggal.so`, BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, the stepping is compiled as
/// the eC writes it, with four constants folded (`1 + 1e-11`, `5 - 1e-11`, `6 + 1e-11` and
/// `2e-11`); the four crossings are not, and are performed here as compiled: sites MV-1 to MV-7,
/// marked in [`move5x6_within`]. None of the seven moves an answer for the arguments the sub-zone
/// generator passes, whose origins lie inside the layout: there the differences they regroup
/// are exact on a crossing, and the folded constants were not seen to part from the eC's sums.
/// That does not hold for every input: from the origin `(0.0, 4.714285714285714)`, outside the
/// layout, with the offset `(1/7, -5/7)`, one turn a crossing and the last crossing asked, the
/// compiled order of MV-1 gives `y` = `0x401f6db6db6db6db`, as the engine does, and the eC's
/// order one unit in the last place more. The sites follow the library all the same. So does
/// site MV-8, the bounds of the nudged point, which parts from the eC at a NaN alone.
pub(crate) fn move5x6(
    o: (f64, f64),
    dx: f64,
    dy: f64,
    n_rotations: i32,
    adjust: Option<&mut (f64, f64)>,
    final_cross: bool,
) -> Option<(f64, f64)> {
    move5x6_within(
        o,
        dx,
        dy,
        n_rotations,
        adjust,
        final_cross,
        MOVE5X6_MAX_TURNS,
    )
}

/// The walk of [`move5x6`], given up after `max_turns` turns of its loop: the body of the eC's
/// `move5x6` (ri5x6.ec:1723-1931), whose own loop has no bound.
#[allow(
    clippy::too_many_lines,
    reason = "one eC function, kept whole so that it reads against its source line by line"
)]
fn move5x6_within(
    o: (f64, f64),
    mut dx: f64,
    mut dy: f64,
    n_rotations: i32,
    mut adjust: Option<&mut (f64, f64)>,
    final_cross: bool,
    max_turns: u32,
) -> Option<(f64, f64)> {
    const EPS: f64 = 1e-11;
    // `Sgn` as the eC's `int` (see [`sgn`]): the library compares with zero and takes one step
    // forward, one back, or none.
    let step = |d: f64| i32::from(d > 0.0) - i32::from(d < 0.0);

    let (mut x, mut y) = o;
    if x < 0.0 && y < ONE_PLUS_1E11 {
        x += 5.0;
        y += 5.0;
    }
    if dx.abs() < EPS && dy.abs() < EPS {
        return Some((x, y));
    }

    // The eC's `while(true)` (ri5x6.ec:1736), with the bound it lacks.
    for _ in 0..max_turns {
        // The point by which the rhombus is chosen: the walk's own point, nudged off the lines
        // of the layout. The eC's `north = !(cdy - cdx > 1)` (ri5x6.ec:1743) is read before
        // either coordinate is nudged (`0x39041` to `0x39054`).
        let (mut cdx, mut cdy) = (x, y);
        let south = y - x > 1.0;
        let iy = floor_i32(y + EPS);
        let ix = floor_i32(x + EPS);
        let mut do_nudge = true;

        // On a horizontal line of the layout, where it bounds a rhombus (ri5x6.ec:1749-1755).
        if (y - f64::from(iy)).abs() < EPS
            && (0..=5).contains(&iy)
            && x > f64::from(iy - 1) - EPS
            && x < f64::from(iy) + EPS
        {
            // `cdy += 2*Sgn(dy) * EPS` (ri5x6.ec:1752): the folded `2e-11` added, subtracted,
            // or nothing at all for an offset of zero (`0x38eda` to `0x38eea`, `0x39795`).
            if dy > 0.0 {
                cdy += TWICE_1E11;
            } else if dy < 0.0 {
                cdy -= TWICE_1E11;
            }
            if x > f64::from(iy - 1) + EPS && x < f64::from(iy) - EPS {
                do_nudge = false;
            }
        }
        // On a vertical line of the layout, likewise (ri5x6.ec:1756-1762).
        if (x - f64::from(ix)).abs() < EPS
            && (0..=5).contains(&ix)
            && y > f64::from(ix) - EPS
            && y < f64::from(ix + 1) + EPS
        {
            // `cdx += 2*Sgn(dx) * EPS` (ri5x6.ec:1759), likewise, in the two copies that the
            // library holds of it (`0x38fed` to `0x39008` with `0x39770`, and `0x38f62` to
            // `0x38f7d` with `0x39858`).
            if dx > 0.0 {
                cdx += TWICE_1E11;
            } else if dx < 0.0 {
                cdx -= TWICE_1E11;
            }
            if y > f64::from(ix) + EPS && y < f64::from(ix + 1) - EPS {
                do_nudge = false;
            }
        }

        if do_nudge {
            if south {
                cdx += EPS;
                cdy -= EPS;
            } else {
                cdx -= EPS;
                cdy += EPS;
            }
        }

        // The nudged point is folded into the layout, and the walk's point with it
        // (ri5x6.ec:1778-1792).
        if cdx < 0.0 && cdy < ONE_PLUS_1E11 {
            cdx += 5.0;
            cdy += 5.0;
            x += 5.0;
            y += 5.0;
        }
        if cdx > 5.0 && cdy > FIVE_MINUS_1E11 {
            cdx -= 5.0;
            cdy -= 5.0;
            x -= 5.0;
            y -= 5.0;
        }
        if cdy < 0.0 && cdx < EPS {
            cdx += 5.0;
            cdy += 5.0;
            x += 5.0;
            y += 5.0;
        }

        // faithful: site MV-8. The four bounds `if(cdx < 0) cdx = 0` and its fellows
        // (ri5x6.ec:1794-1797) are a `maxsd` with zero and a `minsd` with 5 or with 6
        // (`0x390cd` to `0x390e0`), and `(int)floor` of a number so bounded
        // (ri5x6.ec:1799-1800) is a bare `cvttsd2si` (`0x390e5`, `0x390ea`). The site is inert
        // on every finite input, where the two forms give the same rhombus; they part at a NaN
        // alone, which the library takes for zero.
        let cdx = if cdx > 0.0 { cdx } else { 0.0 };
        let cdy = if cdy > 0.0 { cdy } else { 0.0 };
        let cx = cvtt_i32(if cdx < 5.0 { cdx } else { 5.0 });
        let cy = cvtt_i32(if cdy < 6.0 { cdy } else { 6.0 });

        // The step of this turn: the offset, or as much of it as stays in the rhombus
        // (ri5x6.ec:1802-1803); `Max` and `Min` are the library's `maxsd` and `minsd`, which
        // keep the offset where the two are equal.
        let mut px = if dx < 0.0 {
            let to_edge = f64::from(cx) - x;
            if to_edge > dx { to_edge } else { dx }
        } else {
            let to_edge = f64::from(cx + 1) - x;
            if to_edge < dx { to_edge } else { dx }
        };
        let mut py = if dy < 0.0 {
            let to_edge = f64::from(cy) - y;
            if to_edge > dy { to_edge } else { dy }
        } else {
            let to_edge = f64::from(cy + 1) - y;
            if to_edge < dy { to_edge } else { dy }
        };
        // Cut short in one coordinate, the step is cut in the same proportion in the other
        // (ri5x6.ec:1805-1812), as written (`0x39135` to `0x3915d`, `0x39560` to `0x39570`).
        if dx != 0.0 && dy != 0.0 {
            let pkx = px / dx;
            let pky = py / dy;
            if pkx < pky {
                py = pkx * dy;
            } else if pky < pkx {
                px = pky * dx;
            }
        }

        x += px;
        y += py;
        // What is left of the offset, which the eC forms twice (ri5x6.ec:1819, :1877-1878) and
        // the library once (`0x3916b`, `0x39175`).
        let rest = (dx - px, dy - py);

        if !final_cross && rest.0.abs() < EPS && rest.1.abs() < EPS {
            return Some((x, y));
        }

        // The rhombus entered (ri5x6.ec:1823-1837), by the offset as it stood before this step.
        let at_vertex = (x - f64::from(cvtt_i32(x + 0.5))).abs() < EPS
            && (y - f64::from(cvtt_i32(y + 0.5))).abs() < EPS;
        // `floor(c.x + EPS * Sgn(dx))` (ri5x6.ec:1826-1827): the epsilon added, subtracted, or
        // absent (`0x39201` to `0x39216`, `0x39470`).
        let mut nx = floor_i32(if dx > 0.0 {
            x + EPS
        } else if dx < 0.0 {
            x - EPS
        } else {
            x
        });
        let mut ny = floor_i32(if dy > 0.0 {
            y + EPS
        } else if dy < 0.0 {
            y - EPS
        } else {
            y
        });
        if at_vertex && dx.abs() > EPS && dy.abs() <= EPS {
            ny = cy;
            nx = cx + step(dx);
        } else if at_vertex && dy.abs() > EPS && dx.abs() <= EPS {
            nx = cx;
            ny = cy + step(dy);
        } else if nx > cx && ny < cy {
            nx = cx;
        } else if nx < cx && ny > cy {
            ny = cy;
        }

        // The crossing of an interruption, and the turn it gives: -1 clockwise, 1 anticlockwise
        // (ri5x6.ec:1839-1875).
        let mut rotation = 0;
        if nx != cx || ny != cy {
            let in_north = (cx + cy) & 1 == 0;
            if in_north {
                if ny == cy && nx == cx + 1 {
                    // faithful: sites MV-4 and MV-5, crossing to the right in the north.
                    // `iy = (int)(c.x - 1 + 1E-11)` (ri5x6.ec:1848) is
                    // `(int)(c.x - fl(1 - 1e-11))`, the constant at `.rodata` `0x4a018`
                    // (`0x396a7` to `0x396bd`), and `c = { iy + 2 - (c.y - iy), c.x }`
                    // (ri5x6.ec:1849) is `{ (iy + 2) + (iy - c.y), c.x }` (`0x396c9` to
                    // `0x396dc`). The sum `iy + 2` wraps, as the `lea` at `0x396c9` does.
                    let iy = cvtt_i32(x - ONE_MINUS_1E11);
                    (x, y) = (f64::from(iy.wrapping_add(2)) + (f64::from(iy) - y), x);
                    rotation = -1;
                } else if nx == cx && ny == cy - 1 {
                    // faithful: site MV-1, crossing to the left in the north.
                    // `c = { c.y, ix - (c.x - ix) }` (ri5x6.ec:1855) is
                    // `{ c.y, (ix + ix) - c.x }` (`0x395db` to `0x395e3`); `ix` is as written
                    // (ri5x6.ec:1854, `0x395cf`, `0x395d3`).
                    let ix = f64::from(cvtt_i32(y + EPS));
                    (x, y) = (y, (ix + ix) - x);
                    rotation = 1;
                }
            } else if nx == cx && ny == cy + 1 {
                // faithful: sites MV-2 and MV-3, crossing to the right in the south.
                // `ix = (int)(c.y - 2 + 1E-11)` (ri5x6.ec:1864) is `(int)(c.y - fl(2 - 1e-11))`,
                // the constant at `.rodata` `0x4a258` (`0x3964e` to `0x39668`), and
                // `c = { c.y - 1, ix + 3 - (c.x - ix) }` (ri5x6.ec:1865) is
                // `{ c.y - 1, (ix + 3) + (ix - c.x) }` (`0x3967c` to `0x3968b`). The sum
                // `ix + 3` wraps, as the `lea` at `0x3967c` does.
                let ix = cvtt_i32(y - TWO_MINUS_1E11);
                (x, y) = (y - 1.0, f64::from(ix.wrapping_add(3)) + (f64::from(ix) - x));
                rotation = 1;
            } else if ny == cy && nx == cx - 1 {
                // faithful: sites MV-6 and MV-7, crossing to the left in the south.
                // `iy = (int)(c.x + 1 + 1E-11)` (ri5x6.ec:1870) is `(int)(fl(1 + 1e-11) + c.x)`,
                // the constant at `.rodata` `0x4a010` (`0x39720` to `0x39736`), and
                // `c = { iy - 1 - (c.y - iy), c.x + 1 }` (ri5x6.ec:1871) is
                // `{ (iy - 1) + (iy - c.y), 1 + c.x }` (`0x39742` to `0x3975d`). The sum
                // `iy - 1` wraps, as the `lea` at `0x39742` does: from the least `int`, which
                // the conversion gives where its operand is beyond an `int`, to the greatest.
                let iy = cvtt_i32(ONE_PLUS_1E11 + x);
                (x, y) = (f64::from(iy.wrapping_sub(1)) + (f64::from(iy) - y), 1.0 + x);
                rotation = -1;
            }
        }

        (dx, dy) = rest;

        // The point reached is folded into the layout (ri5x6.ec:1880-1885): forward by a period
        // from before it, to the left or below, which the eC tests in two arms, and back by
        // one from beyond it.
        if (x < -EPS && y < ONE_PLUS_1E11) || (y < 0.0 && x < EPS) {
            x += 5.0;
            y += 5.0;
        } else if x > FIVE_MINUS_1E11 && y > SIX_PLUS_1E11 {
            x -= 5.0;
            y -= 5.0;
        }

        // What is left of the offset is turned, and the caller's direction with it
        // (ri5x6.ec:1888-1922).
        if rotation != 0 {
            for _ in 0..n_rotations {
                if rotation == -1 {
                    // Sixty degrees clockwise.
                    (dx, dy) = (dx - dy, dx);
                    if let Some(a) = adjust.as_deref_mut() {
                        *a = (a.0 - a.1, a.0);
                    }
                } else {
                    // Sixty degrees anticlockwise.
                    (dx, dy) = (dy, dy - dx);
                    if let Some(a) = adjust.as_deref_mut() {
                        *a = (a.1, a.1 - a.0);
                    }
                }
            }
        }

        if dx.abs() < EPS {
            dx = 0.0;
        }
        if dy.abs() < EPS {
            dy = 0.0;
        }
        if dx == 0.0 && dy == 0.0 {
            return Some((x, y));
        }
    }
    None
}

#[cfg(test)]
#[allow(clippy::type_complexity, reason = "the tables of test vectors")]
mod tests {
    //! Every expected value below is the shipped engine's own answer, read by calling the
    //! function at its address in `libdggal.so`, BuildID
    //! `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, and written out as the shortest decimal that
    //! reads back to the same double; the comparisons are bit for bit. The tests of
    //! `add_intermediate_points` are the exception, and say where their values come from; those
    //! of `move5x6` hold their arguments and the engine's answers as bits.
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
        // `refine5x6`: the bounds of the layout, which are also the terms of sites R56-9 and
        // R56-10, and the epsilon either side of zero.
        assert_eq!(FIVE_PLUS_1E11.to_bits(), 0x4014_0000_0000_2bfb); // 0x4a188
        assert_eq!(SIX_PLUS_1E11.to_bits(), 0x4018_0000_0000_2bfb); // 0x4a2d0
        assert_eq!(FIVE_MINUS_1E11.to_bits(), 0x4013_ffff_ffff_d405); // 0x49f48
        assert_eq!(1e-11_f64.to_bits(), 0x3da5_fd7f_e179_6495); // 0x49ec0
        assert_eq!((-1e-11_f64).to_bits(), 0xbda5_fd7f_e179_6495); // 0x49f08
        // `move5x6`: twice the epsilon, beside the five it shares with the functions above
        // (`0x4a010`, `0x4a018`, `0x4a258`, `0x49f48` and `0x4a2d0`).
        assert_eq!(TWICE_1E11.to_bits(), 0x3db5_fd7f_e179_6495); // 0x49fa8
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

    /// One call of `add_intermediate_points`, as bits: the two ends of the edge, the number of
    /// divisions, the two points of the interruption where the edge crosses one, and the points
    /// the call adds. `crs84` is true in every one of them.
    ///
    /// The calls below were captured from a prototype of the refined boundary, each within the
    /// refined boundary of the zone and at the refinement its comment names, where the
    /// prototype's ring is the engine's `getZoneRefinedWGS84Vertices`, point for point and bit
    /// for bit, and where the site under test, put back in the eC's order alone, moves that
    /// ring away from the engine's. The engine's `addIntermediatePoints` takes an eC array and
    /// cannot be called at its address as the functions above were; the refined rings of the
    /// grids compare this function with the engine itself.
    ///
    /// Beyond these seven, the function was run once over every call the prototype made on a
    /// sample of the six grids on which each of its rings was the engine's, at refinements -7,
    /// -1, 0, 1, 2, 3, 7 and 50 and, on aperture 7, in the 5x6 CRS at 0, 1, 2, 12 and 20:
    /// 12,567,626 calls, 697,026 of them over an interruption, 1,908,026 with `crs84` false and
    /// 38,837 under the pole rule (3,956 of those on a plain edge), adding 114,983,836 points,
    /// and it added the same points at every one.
    struct Call {
        p: [u64; 2],
        n: [u64; 2],
        n_divisions: i32,
        interruption: Option<[[u64; 2]; 2]>,
        points: &'static [[u64; 2]],
    }

    fn point(b: [u64; 2]) -> (f64, f64) {
        (f64::from_bits(b[0]), f64::from_bits(b[1]))
    }

    fn assert_adds_the_captured_points(what: &str, call: &Call) {
        let mut points = Vec::new();
        add_intermediate_points(
            &mut points,
            point(call.p),
            point(call.n),
            call.n_divisions,
            call.interruption.map(|[i1, i2]| (point(i1), point(i2))),
            true,
        );
        assert_eq!(points.len(), call.points.len(), "{what}: how many points");
        for (k, (got, want)) in points.iter().zip(call.points).enumerate() {
            assert!(
                bits(*got) == (want[0], want[1]),
                "{what}: point {k} is ({:#018x}, {:#018x}), and the captured one is ({:#018x}, {:#018x})",
                got.0.to_bits(),
                got.1.to_bits(),
                want[0],
                want[1]
            );
        }
    }

    /// IGEO7 `06` (`0x6fffffffffffffff`), a pentagon of level 0, at refinement 7: a plain edge.
    const PLAIN_A7: Call = Call {
        p: [0x3fd5_5555_5555_5556, 0x3ffa_aaaa_aaaa_aaab],
        n: [0x3fe5_5555_5555_5556, 0x3ff5_5555_5555_5556],
        n_divisions: 7,
        interruption: None,
        points: &[
            [0x3fd5_5555_5555_5556, 0x3ffa_aaaa_aaaa_aaab],
            [0x3fd8_6186_1861_8619, 0x3ff9_e79e_79e7_9e7a],
            [0x3fdb_6db6_db6d_b6dc, 0x3ff9_2492_4924_924a],
            [0x3fde_79e7_9e79_e7a0, 0x3ff8_6186_1861_8619],
            [0x3fe0_c30c_30c3_0c31, 0x3ff7_9e79_e79e_79e8],
            [0x3fe2_4924_9249_2493, 0x3ff6_db6d_b6db_6db8],
            [0x3fe3_cf3c_f3cf_3cf4, 0x3ff6_1861_8618_6187],
        ],
    };

    /// ISEA3H `A0-0-D` (`0x0000000000000003`), a hexagon of level 1, at refinement 7: a plain edge.
    const PLAIN_A3: Call = Call {
        p: [0x3fe5_5555_5555_5555, 0x3ff0_0000_0000_0000],
        n: [0x3fe5_5555_5555_5555, 0x3fe5_5555_5555_5555],
        n_divisions: 7,
        interruption: None,
        points: &[
            [0x3fe5_5555_5555_5555, 0x3ff0_0000_0000_0000],
            [0x3fe5_5555_5555_5555, 0x3fee_79e7_9e79_e79e],
            [0x3fe5_5555_5555_5555, 0x3fec_f3cf_3cf3_cf3d],
            [0x3fe5_5555_5555_5555, 0x3feb_6db6_db6d_b6db],
            [0x3fe5_5555_5555_5555, 0x3fe9_e79e_79e7_9e7a],
            [0x3fe5_5555_5555_5555, 0x3fe8_6186_1861_8618],
            [0x3fe5_5555_5555_5555, 0x3fe6_db6d_b6db_6db6],
        ],
    };

    /// IGEO7 `0200` (`0x203fffffffffffff`), a pentagon of level 2, at refinement 0: the edge over the
    /// interruption.
    const BEFORE_A7: Call = Call {
        p: [0x3ff1_8618_6186_1861, 0x3ff0_c30c_30c3_0c30],
        n: [0x3fee_79e7_9e79_e79d, 0x3fec_f3cf_3cf3_cf3c],
        n_divisions: 20,
        interruption: Some([
            [0x3ff1_2492_4924_9249, 0x3ff0_0000_0000_0000],
            [0x3ff0_0000_0000_0000, 0x3fed_b6db_6db6_db6e],
        ]),
        points: &[
            [0x3ff1_8618_6186_1861, 0x3ff0_c30c_30c3_0c30],
            [0x3ff1_7c57_c57c_57c5, 0x3ff0_af8a_f8af_8af8],
            [0x3ff1_7297_2972_9729, 0x3ff0_9c09_c09c_09c0],
            [0x3ff1_68d6_8d68_d68d, 0x3ff0_8888_8888_8888],
            [0x3ff1_5f15_f15f_15f1, 0x3ff0_7507_5075_0750],
            [0x3ff1_5555_5555_5555, 0x3ff0_6186_1861_8618],
            [0x3ff1_4b94_b94b_94b9, 0x3ff0_4e04_e04e_04df],
            [0x3ff1_41d4_1d41_d41d, 0x3ff0_3a83_a83a_83a7],
            [0x3ff1_3813_8138_1381, 0x3ff0_2702_7027_026f],
            [0x3ff1_2e52_e52e_52e5, 0x3ff0_1381_3813_8137],
            [0x3ff1_2492_4924_9249, 0x3ff0_0000_0000_0000],
            [0x3ff0_0000_0000_0000, 0x3fed_b6db_6db6_db6e],
            [0x3fef_ffff_ffff_fffe, 0x3fed_b6db_6db6_db6d],
            [0x3fef_d8fd_8fd8_fd8e, 0x3fed_a35a_35a3_5a35],
            [0x3fef_b1fb_1fb1_fb1e, 0x3fed_8fd8_fd8f_d8fd],
            [0x3fef_8af8_af8a_f8ae, 0x3fed_7c57_c57c_57c5],
            [0x3fef_63f6_3f63_f63e, 0x3fed_68d6_8d68_d68d],
            [0x3fef_3cf3_cf3c_f3ce, 0x3fed_5555_5555_5555],
            [0x3fef_15f1_5f15_f15e, 0x3fed_41d4_1d41_d41c],
            [0x3fee_eeee_eeee_eeed, 0x3fed_2e52_e52e_52e4],
            [0x3fee_c7ec_7ec7_ec7d, 0x3fed_1ad1_ad1a_d1ac],
            [0x3fee_a0ea_0ea0_ea0d, 0x3fed_0750_7507_5074],
        ],
    };

    /// RTEA3H `F1-16C8-A` (`0x0a20000000005b20`), a hexagon of level 10 across an interruption, at
    /// refinement 0.
    const BEFORE_A3: Call = Call {
        p: [0x3f56_7980_e0bf_08c8, 0x3ff1_8eec_af95_3edc],
        n: [0x4013_9c44_d41a_b049, 0x4017_fe98_67f1_f410],
        n_divisions: 5,
        interruption: Some([
            [0x0000_0000_0000_0000, 0x3ff1_8c1d_7f79_26fb],
            [0x4013_9cf8_a021_b641, 0x4018_0000_0000_0000],
        ]),
        points: &[
            [0x3f56_7980_e0bf_08c8, 0x3ff1_8eec_af95_3edc],
            [0x3f4a_f834_40e5_412f, 0x3ff1_8dcd_02bd_354f],
            [0x3f31_facd_8098_e19c, 0x3ff1_8cad_55e5_2bc2],
            [0x0000_0000_0000_0000, 0x3ff1_8c1d_7f79_26fb],
            [0x4013_9cf8_a021_b641, 0x4018_0000_0000_0000],
            [0x4013_9cd4_aa86_b510, 0x4017_ffb8_14c9_fd9d],
            [0x4013_9c8c_bf50_b2ac, 0x4017_ff28_3e5d_f8d7],
        ],
    };

    /// IGEO7 `02323` (`0x269fffffffffffff`), a hexagon of level 3 across an interruption, at
    /// refinement 0.
    const AFTER_A7: Call = Call {
        p: [0x3ff7_f211_6a3b_35fd, 0x3ff0_37ba_5713_280e],
        n: [0x3fef_908b_51d9_afe4, 0x3fde_e95c_4ca0_37b9],
        n_divisions: 15,
        interruption: Some([
            [0x3ff8_0000_0000_0001, 0x3ff0_0000_0000_0000],
            [0x3ff0_0000_0000_0000, 0x3fdf_ffff_ffff_fffc],
        ]),
        points: &[
            [0x3ff7_f211_6a3b_35fd, 0x3ff0_37ba_5713_280e],
            [0x3ff7_f470_71c0_c937, 0x3ff0_2e3e_38fc_db28],
            [0x3ff7_f6cf_7946_5c70, 0x3ff0_24c2_1ae6_8e41],
            [0x3ff7_f92e_80cb_efaa, 0x3ff0_1b45_fcd0_415b],
            [0x3ff7_fb8d_8851_82e4, 0x3ff0_11c9_deb9_f474],
            [0x3ff7_fdec_8fd7_161e, 0x3ff0_084d_c0a3_a78e],
            [0x3ff8_0000_0000_0001, 0x3ff0_0000_0000_0000],
            [0x3ff0_0000_0000_0000, 0x3fdf_ffff_ffff_fffc],
            [0x3fef_fe7a_9a06_b64d, 0x3fdf_fc32_8110_c7bc],
            [0x3fef_f243_9201_b597, 0x3fdf_dda8_ed04_45f5],
            [0x3fef_e60c_89fc_b4e0, 0x3fdf_bf1f_58f7_c42d],
            [0x3fef_d9d5_81f7_b42a, 0x3fdf_a095_c4eb_4266],
            [0x3fef_cd9e_79f2_b374, 0x3fdf_820c_30de_c09e],
            [0x3fef_c167_71ed_b2bd, 0x3fdf_6382_9cd2_3ed7],
            [0x3fef_b530_69e8_b207, 0x3fdf_44f9_08c5_bd10],
            [0x3fef_a8f9_61e3_b151, 0x3fdf_266f_74b9_3b48],
            [0x3fef_9cc2_59de_b09a, 0x3fdf_07e5_e0ac_b980],
        ],
    };

    /// ISEA3H `B3-6-A` (`0x0260000000000018`), a hexagon of level 2 across an interruption, at
    /// refinement 7.
    const AFTER_A3: Call = Call {
        p: [0x3ff1_c71c_71c7_1c72, 0x4004_71c7_1c71_c71c],
        n: [0x3fdc_71c7_1c71_c71f, 0x3ffe_38e3_8e38_e38e],
        n_divisions: 7,
        interruption: Some([
            [0x3ff0_0000_0000_0000, 0x4004_0000_0000_0000],
            [0x3fe0_0000_0000_0000, 0x4000_0000_0000_0000],
        ]),
        points: &[
            [0x3ff1_c71c_71c7_1c72, 0x4004_71c7_1c71_c71c],
            [0x3ff1_4514_5145_1451, 0x4004_5145_1451_4514],
            [0x3ff0_c30c_30c3_0c31, 0x4004_30c3_0c30_c30c],
            [0x3ff0_4104_1041_0410, 0x4004_1041_0410_4104],
            [0x3ff0_0000_0000_0000, 0x4004_0000_0000_0000],
            [0x3fe0_0000_0000_0000, 0x4000_0000_0000_0000],
            [0x3fdf_7df7_df7d_f7e0, 0x3fff_befb_efbe_fbf0],
            [0x3fde_79e7_9e79_e7a0, 0x3fff_3cf3_cf3c_f3cf],
            [0x3fdd_75d7_5d75_d75f, 0x3ffe_baeb_aeba_ebae],
        ],
    };

    /// IGEO7 `00` (`0x0fffffffffffffff`), the pentagon of level 0 at the north pole, at refinement
    /// 3: the edge that crosses the interruption at the pole itself, from (5, 4.5) to (0.5, 0).
    const POLAR_A7: Call = Call {
        p: [0x4012_aaaa_aaaa_aaab, 0x4011_5555_5555_5555],
        n: [0x3fe5_5555_5555_5556, 0x3fd5_5555_5555_5555],
        n_divisions: 3,
        interruption: Some([
            [0x4014_0000_0000_0000, 0x4012_0000_0000_0000],
            [0x3fe0_0000_0000_0000, 0x0000_0000_0000_0000],
        ]),
        points: &[
            [0x4012_aaaa_aaaa_aaab, 0x4011_5555_5555_5555],
            [0x4012_b60b_60b6_0b61, 0x4011_5b05_b05b_05b0],
            [0x4013_9999_9999_999a, 0x4011_cccc_cccc_cccd],
            [0x4013_a4fa_4fa4_fa50, 0x4011_d27d_27d2_7d28],
            [0x4013_b05b_05b0_5b06, 0x4011_d82d_82d8_2d83],
            [0x4013_bbbb_bbbb_bbbc, 0x4011_dddd_dddd_ddde],
            [0x4013_c71c_71c7_1c72, 0x4011_e38e_38e3_8e39],
            [0x4013_d27d_27d2_7d28, 0x4011_e93e_93e9_3e94],
            [0x4013_dddd_dddd_ddde, 0x4011_eeee_eeee_eeef],
            [0x4013_e93e_93e9_3e94, 0x4011_f49f_49f4_9f4a],
            [0x4013_f49f_49f4_9f4a, 0x4011_fa4f_a4fa_4fa5],
            [0x4014_0000_0000_0000, 0x4012_0000_0000_0000],
            [0x3fe0_0000_0000_0000, 0x0000_0000_0000_0000],
            [0x3fe0_0000_0000_0000, 0x3c96_c16c_16c1_6c15],
            [0x3fe0_2d82_d82d_82d9, 0x3f86_c16c_16c1_6c43],
            [0x3fe0_5b05_b05b_05b1, 0x3f96_c16c_16c1_6c2c],
            [0x3fe0_8888_8888_8889, 0x3fa1_1111_1111_111b],
            [0x3fe0_b60b_60b6_0b61, 0x3fa6_c16c_16c1_6c20],
            [0x3fe0_e38e_38e3_8e39, 0x3fac_71c7_1c71_c726],
            [0x3fe1_1111_1111_1112, 0x3fb1_1111_1111_1115],
            [0x3fe1_3e93_e93e_93ea, 0x3fb3_e93e_93e9_3e98],
            [0x3fe1_6c16_c16c_16c2, 0x3fb6_c16c_16c1_6c1b],
            [0x3fe1_9999_9999_999a, 0x3fb9_9999_9999_999d],
            [0x3fe1_c71c_71c7_1c72, 0x3fbc_71c7_1c71_c720],
            [0x3fe1_f49f_49f4_9f4a, 0x3fbf_49f4_9f49_f4a3],
            [0x3fe2_2222_2222_2223, 0x3fc1_1111_1111_1113],
            [0x3fe2_4fa4_fa4f_a4fb, 0x3fc2_7d27_d27d_27d4],
            [0x3fe2_7d27_d27d_27d3, 0x3fc3_e93e_93e9_3e95],
            [0x3fe2_aaaa_aaaa_aaab, 0x3fc5_5555_5555_5557],
            [0x3fe2_d82d_82d8_2d83, 0x3fc6_c16c_16c1_6c18],
            [0x3fe3_05b0_5b05_b05c, 0x3fc8_2d82_d82d_82d9],
            [0x3fe3_3333_3333_3334, 0x3fc9_9999_9999_999a],
            [0x3fe3_60b6_0b60_b60c, 0x3fcb_05b0_5b05_b05c],
            [0x3fe3_8e38_e38e_38e4, 0x3fcc_71c7_1c71_c71d],
        ],
    };

    /// Site AIP-1 (`ri5x6.ec:702-703`, `0x370f6`): the points before an interruption are
    /// `p + j * (pi1 * (1 / t))`. In the eC's order, `p + j * pi1 / t`, a point of each of these
    /// edges falls a unit in the last place from the engine's.
    #[test]
    fn the_points_before_an_interruption_follow_the_compiled_order() {
        assert_adds_the_captured_points("IGEO7 0200", &BEFORE_A7);
        assert_adds_the_captured_points("RTEA3H F1-16C8-A", &BEFORE_A3);
    }

    /// Site AIP-2 (`ri5x6.ec:714-715`, `0x37105`): the points after an interruption are
    /// `i2 + (j - t) * (i2n * (1 / (nDivisions - t)))`, where the eC writes
    /// `i2 + (j - t) * i2n / (nDivisions - t)`.
    #[test]
    fn the_points_after_an_interruption_follow_the_compiled_order() {
        assert_adds_the_captured_points("IGEO7 02323", &AFTER_A7);
        assert_adds_the_captured_points("ISEA3H B3-6-A", &AFTER_A3);
    }

    /// Site AIP-3 (`ri5x6.ec:724`, `0x373f3`): the points of a plain edge are
    /// `p + j * (d * (1 / nDivisions))`, where the eC writes `p + j * d / nDivisions`.
    #[test]
    fn the_points_of_a_plain_edge_follow_the_compiled_order() {
        assert_adds_the_captured_points("IGEO7 06", &PLAIN_A7);
        assert_adds_the_captured_points("ISEA3H A0-0-D", &PLAIN_A3);
    }

    /// The pole rule (`ri5x6.ec:655-676`, `:698`): an interruption that ends on a pole's
    /// image makes the edge twenty times finer, 60 divisions for the 3 asked here, and the walk
    /// then takes them twenty at a time, save within twenty of the interruption, where it
    /// takes every one. So 34 points, and neither the 5 of the plain rule nor the 62 of every
    /// division.
    #[test]
    fn an_interruption_at_a_pole_is_refined_twenty_times_over() {
        assert_adds_the_captured_points("IGEO7 00", &POLAR_A7);
    }

    /// What the eC's `int` arithmetic gives where Rust's would give something else, and the
    /// branches none of the seven calls above enters. The expected values follow the compiled
    /// function as read at `0x36e70`, and are not captured answers. The run over every call of
    /// the prototype, described above, did enter these branches, save two that no sampled zone
    /// presents: a NaN `t` and the wrapped product.
    #[test]
    fn add_intermediate_points_keeps_the_compiled_integer_semantics() {
        let added = |p, n, n_divisions, interruption, crs84| {
            let mut points = Vec::new();
            add_intermediate_points(&mut points, p, n, n_divisions, interruption, crs84);
            points
        };
        let (p, n) = ((1.25, 1.5), (1.5, 1.75));
        // No division asked is one division (`ri5x6.ec:674`, `0x36fc0`), and one division of a
        // plain edge is its first end alone (`ri5x6.ec:719-720`).
        assert_eq!(added(p, n, 0, None, true), [p]);
        assert_eq!(added(p, n, 1, None, true), [p]);
        assert_eq!(added(p, n, 2, None, true), [p, (1.375, 1.625)]);
        // A negative number of divisions adds nothing to a plain edge (the `jle` at `0x373d3`),
        // and the first end alone to an interrupted one.
        assert_eq!(added(p, n, -7, None, true), []);
        let across = Some(((1.375, 1.5), (1.5, 1.625)));
        assert_eq!(added(p, n, -7, across, true), [p]);
        // One division of an interrupted edge: the interruption's two points are added in the
        // 5x6 CRS, and left out towards WGS84 (`ri5x6.ec:694-696`).
        assert_eq!(
            added(p, n, 1, across, false),
            [p, (1.375, 1.5), (1.5, 1.625)]
        );
        assert_eq!(added(p, n, 1, across, true), [p]);
        // `(int)t` of a NaN is `i32::MIN` (`cvttsd2si` at `0x370ee`), which no `j` equals, so
        // an edge of no length on either side of its interruption adds neither point of it;
        // Rust's `as` would give 0, and with it both points at `j == 1`.
        assert_eq!(added(p, n, 2, Some((p, n)), true), [p]);
        // Each of the four tests on the two ends alone (`ri5x6.ec:664-671`) refines twenty
        // times over, with or without an interruption; an edge that meets none does not.
        for (p, n) in [
            ((0.75, 0.005), (0.25, 0.005)),
            ((4.995, 4.75), (4.995, 4.25)),
            ((1.25, 2.995), (1.75, 2.995)),
            ((2.005, 3.25), (2.005, 3.75)),
        ] {
            assert_eq!(added(p, n, 1, None, true).len(), 20, "{p:?} to {n:?}");
            assert_eq!(added(p, n, 1, None, false).len(), 20, "{p:?} to {n:?}");
            assert_eq!(added(n, p, 1, None, true).len(), 1, "{n:?} to {p:?}");
        }
        // The product by twenty is a 32-bit one that wraps (`lea` and `shl` at `0x37398`):
        // 107,374,183 divisions become -2,147,483,636, and nothing is added.
        let (p, n) = ((0.75, 0.005), (0.25, 0.005));
        assert_eq!(added(p, n, 107_374_183, None, true), []);
    }

    /// One call of `move5x6` and the engine's answer to it, as bits: the origin, the offset, the
    /// turns a crossing, the caller's direction as it is given and as it is left (where one is
    /// passed), whether the last crossing is asked, and the point returned.
    ///
    /// The engine's function does not return from some arguments that its callers never give it
    /// (a point 1e-12 off an edge with a long offset, a corner of the layout, two turns a
    /// crossing from an arbitrary point), so none of these calls is invented. Each is a call
    /// that the aperture-7 sub-zone generator makes for the IGEO7 zone its comment names, at the
    /// depth it names (in `getFirstSubZoneCentroid`, in `iterateI7HSubZones` or, where the
    /// comment says so, in `zoneHasSubZone`), and the answer is the engine's own function's,
    /// read at its address in the library (`0x38d70`).
    ///
    /// Beyond these vectors, the function was run against the engine's over the calls made for
    /// every zone of levels 0 to 4 of IGEO7 and its pentagons of levels 5 and 6, at depths 1 to
    /// 3, and at depth 4 for the zones of levels 0 to 2 and the pentagons of levels 3 and 4
    /// (each order whole, and the test for a sub-zone at some seven of its entries): 16,528,523
    /// calls, and the engine's point and direction at every one.
    ///
    /// The vectors were chosen so that between them they tell the function from the same
    /// function with one arm spoilt (an arm removed, a nudge or an epsilon left out, a turn
    /// reversed), for every such variant that any call of that run tells apart. What no call of
    /// the run tells apart is pinned by no vector: the third fold of the nudged point
    /// (`ri5x6.ec:1788-1792`), three of the four bounds (`:1794`, `:1796`, `:1797`), the corner
    /// arm along `y` (`:1831-1832`) and a vertical line met with no offset along `x` (`:1759`,
    /// where `Sgn(dx)` is zero), which no call enters; the fold of the origin on entry
    /// (`:1727-1728`), which the first fold of the loop repeats wherever there is an offset; the
    /// bound of `:1795`, since a number between -1 and 0 truncates to 0 as it is; whether the
    /// nudge on a line is `2e-11` or `1e-11` moving up or moving left, and whether the diagonal
    /// nudge is withheld on a vertical line; the seven sites MV-1 to MV-7, which are inert for
    /// origins inside the layout, where every call of that run begins (outside it a site may
    /// move an answer, as MV-1 does from the origin given at [`move5x6`]); and site MV-8, which
    /// only a NaN could tell from the eC's form.
    struct Move {
        o: [u64; 2],
        d: [u64; 2],
        n_rotations: i32,
        direction: Option<[[u64; 2]; 2]>,
        final_cross: bool,
        v: [u64; 2],
    }

    /// The greatest number of turns that any walk took in the measurement described at
    /// [`MOVE5X6_MAX_TURNS`].
    const MOVE5X6_MEASURED_MAX_TURNS: u32 = 4;

    // The bound is at least sixteen times the greatest number of turns measured.
    const _: () = assert!(MOVE5X6_MAX_TURNS >= 16 * MOVE5X6_MEASURED_MAX_TURNS);

    /// Each vector is answered as the engine answers it, and within the greatest number of
    /// turns measured, of which the bound is at least sixteen times, so that a change which
    /// made the generator's walks longer would be seen here before it met the bound.
    fn assert_moves_as_the_engine(what: &str, moves: &[Move]) {
        for (k, m) in moves.iter().enumerate() {
            let walk = |max_turns: u32| {
                let mut direction = m.direction.map(|[given, _]| point(given));
                let v = move5x6_within(
                    point(m.o),
                    f64::from_bits(m.d[0]),
                    f64::from_bits(m.d[1]),
                    m.n_rotations,
                    direction.as_mut(),
                    m.final_cross,
                    max_turns,
                );
                (v.map(bits), direction.map(bits))
            };
            let mut direction = m.direction.map(|[given, _]| point(given));
            let got = move5x6(
                point(m.o),
                f64::from_bits(m.d[0]),
                f64::from_bits(m.d[1]),
                m.n_rotations,
                direction.as_mut(),
                m.final_cross,
            );
            let Some(got) = got else {
                panic!("{what}, vector {k}: the walk is given up")
            };
            assert!(
                bits(got) == (m.v[0], m.v[1]),
                "{what}, vector {k}: ({:#018x}, {:#018x}), the engine's ({:#018x}, {:#018x})",
                got.0.to_bits(),
                got.1.to_bits(),
                m.v[0],
                m.v[1]
            );
            assert_eq!(
                direction.map(bits),
                m.direction.map(|[_, left]| (left[0], left[1])),
                "{what}, vector {k}: the direction as it is left"
            );
            assert_eq!(
                walk(MOVE5X6_MEASURED_MAX_TURNS),
                (Some(bits(got)), direction.map(bits)),
                "{what}, vector {k}: the walk takes more turns than the greatest number measured"
            );
        }
    }

    const WITHIN: &[Move] = &[
        // `061`, depth 1: a step within one rhombus.
        Move {
            o: [0x4011_b6db_6db6_db6d, 0x4016_4924_9249_2492],
            d: [0x3fc2_4924_9249_2492, 0x3fc2_4924_9249_2492],
            n_rotations: 1,
            direction: None,
            final_cross: true,
            v: [0x4012_4924_9249_2492, 0x4016_db6d_b6db_6db7],
        },
        // `0004`, depth 2: a step along x alone.
        Move {
            o: [0x3fee_0a72_f053_9783, 0x3faf_58d0_fac6_87ee],
            d: [0x3f94_e5e0_a72f_0539, 0x0000_0000_0000_0000],
            n_rotations: 1,
            direction: None,
            final_cross: true,
            v: [0x3fee_b1a1_f58d_0fad, 0x3faf_58d0_fac6_87ee],
        },
        // `001`, depth 1: no offset, and the origin is returned.
        Move {
            o: [0x4010_0000_0000_0000, 0x4009_2492_4924_924a],
            d: [0x0000_0000_0000_0000, 0x0000_0000_0000_0000],
            n_rotations: 1,
            direction: None,
            final_cross: true,
            v: [0x4010_0000_0000_0000, 0x4009_2492_4924_924a],
        },
        // `001`, depth 1: no offset and a direction, both returned as given.
        Move {
            o: [0x4010_0000_0000_0000, 0x4009_2492_4924_924a],
            d: [0x8000_0000_0000_0000, 0x0000_0000_0000_0000],
            n_rotations: 1,
            direction: Some([
                [0x3fc2_4924_9249_2492, 0x0000_0000_0000_0000],
                [0x3fc2_4924_9249_2492, 0x0000_0000_0000_0000],
            ]),
            final_cross: true,
            v: [0x4010_0000_0000_0000, 0x4009_2492_4924_924a],
        },
        // `011`, depth 1: from one rhombus into the next, where no interruption lies between them.
        Move {
            o: [0x4011_2492_4924_9248, 0x4012_4924_9249_2492],
            d: [0x3fdb_6db6_db6d_b6db, 0x3fdb_6db6_db6d_b6db],
            n_rotations: 1,
            direction: None,
            final_cross: true,
            v: [0x4012_db6d_b6db_6db6, 0x4014_0000_0000_0000],
        },
        // `064141040441410`, depth 1: a step within one rhombus, deep in the hierarchy.
        Move {
            o: [0x4012_ae15_2d56_b4f6, 0x4017_ffff_ae83_1312],
            d: [0xbec4_5f3b_3bb0_829a, 0xbeb4_5f3b_3bb0_829a],
            n_rotations: 1,
            direction: Some([
                [0x3eb4_5f3b_3bb0_829a, 0x3eb4_5f3b_3bb0_829a],
                [0x3eb4_5f3b_3bb0_829a, 0x3eb4_5f3b_3bb0_829a],
            ]),
            final_cross: true,
            v: [0x4012_ae14_8a5c_db18, 0x4017_ffff_5d06_2623],
        },
        // `006`, depth 1: a step cut short along x at the edge of the rhombus, and along y in the
        // same proportion.
        Move {
            o: [0x3fff_3cf3_cf3c_f3d0, 0x3ff3_0c30_c30c_30c4],
            d: [0x3fa8_6186_1861_8618, 0xbfa8_6186_1861_8618],
            n_rotations: 1,
            direction: None,
            final_cross: false,
            v: [0x4000_0000_0000_0000, 0x3ff2_4924_9249_2494],
        },
        // `055`, depth 1: a step cut short along y at the edge of the rhombus, and along x in the
        // same proportion.
        Move {
            o: [0x400e_1861_8618_6186, 0x400f_9e79_e79e_79e8],
            d: [0x3fb8_6186_1861_8618, 0x3fa8_6186_1861_8618],
            n_rotations: 1,
            direction: None,
            final_cross: false,
            v: [0x400e_db6d_b6db_6db6, 0x4010_0000_0000_0000],
        },
        // `021`, depth 1: what is left of dx, below the epsilon, is taken for nothing.
        Move {
            o: [0x3fd2_4924_9249_2493, 0x3fe2_4924_9249_2494],
            d: [0x3fdb_6db6_db6d_b6db, 0x3fdb_6db6_db6d_b6db],
            n_rotations: 1,
            direction: None,
            final_cross: true,
            v: [0x3fe6_db6d_b6db_6db6, 0x3ff0_0000_0000_0000],
        },
        // `044`, depth 1: what is left of dy, below the epsilon, is taken for nothing.
        Move {
            o: [0x4006_db6d_b6db_6db8, 0x400a_4924_9249_2490],
            d: [0x3fc2_4924_9249_2492, 0x3fc2_4924_9249_2492],
            n_rotations: 1,
            direction: None,
            final_cross: true,
            v: [0x4008_0000_0000_0000, 0x400b_6db6_db6d_b6d8],
        },
    ];

    const CROSSINGS: &[Move] = &[
        // `00`, depth 3: north, to the left (MV-1).
        Move {
            o: [0x3ff8_7d63_43eb_1a1d, 0x3ff0_5397_829c_bc14],
            d: [0x3fb4_e5e0_a72f_0539, 0xbfa4_e5e0_a72f_053a],
            n_rotations: 1,
            direction: None,
            final_cross: true,
            v: [0x3fef_58d0_fac6_87d4, 0x3fd7_829c_bc14_e5e6],
        },
        // `074`, depth 1: south, to the right (MV-2 and MV-3).
        Move {
            o: [0x3feb_6db6_db6d_b6e4, 0x4000_0000_0000_0000],
            d: [0x0000_0000_0000_0000, 0x3fc2_4924_9249_2492],
            n_rotations: 1,
            direction: None,
            final_cross: true,
            v: [0x3ff2_4924_9249_2492, 0x4002_4924_9249_2490],
        },
        // `001`, depth 1: north, to the right (MV-4 and MV-5).
        Move {
            o: [0x400e_db6d_b6db_6db7, 0x400a_4924_9249_2493],
            d: [0x3fd2_4924_9249_2492, 0x0000_0000_0000_0000],
            n_rotations: 1,
            direction: None,
            final_cross: true,
            v: [0x4013_6db6_db6d_b6db, 0x4010_9249_2492_4925],
        },
        // `11`, depth 1: south, to the left (MV-6 and MV-7).
        Move {
            o: [0x4009_2492_4924_9249, 0x4011_b6db_6db6_db6e],
            d: [0xbfd2_4924_9249_2492, 0x3fc2_4924_9249_2492],
            n_rotations: 1,
            direction: None,
            final_cross: true,
            v: [0x4002_4924_9249_2492, 0x400e_db6d_b6db_6db7],
        },
        // `001313131313131`, depth 1, in the test for a sub-zone: north, to the left, deep in the
        // hierarchy (MV-1).
        Move {
            o: [0x4012_0000_cbb8_5055, 0x4010_0000_0000_0000],
            d: [0xbe9a_e41a_fcde_c0e1, 0xbeaa_e41a_fcde_c0e1],
            n_rotations: 1,
            direction: None,
            final_cross: false,
            v: [0x400f_ffff_946f_940d, 0x400b_fffe_32c7_295c],
        },
        // `111306640113066`, depth 1: south, to the right, deep in the hierarchy (MV-2 and MV-3).
        Move {
            o: [0x4003_3333_129a_d46c, 0x400f_ffff_5d06_2622],
            d: [0x0000_0000_0000_0000, 0x3eb4_5f3b_3bb0_829a],
            n_rotations: 1,
            direction: None,
            final_cross: true,
            v: [0x4008_0000_0000_0000, 0x4012_6666_76b2_95ca],
        },
        // `001313131313131`, depth 1: north, to the right, deep in the hierarchy (MV-4 and MV-5).
        Move {
            o: [0x400f_ffff_5d06_2622, 0x400b_ffff_0b89_3936],
            d: [0x3eb4_5f3b_3bb0_829a, 0x0000_0000_0000_0000],
            n_rotations: 1,
            direction: None,
            final_cross: true,
            v: [0x4012_0000_7a3b_6365, 0x4010_0000_0000_0000],
        },
        // `074104360341043`, depth 1, in the test for a sub-zone: south, to the left, deep in the
        // hierarchy (MV-6 and MV-7).
        Move {
            o: [0x3ff0_0000_0000_0000, 0x4002_6665_30be_e20e],
            d: [0xbeaa_e41a_fcde_c0e1, 0xbe9a_e41a_fcde_c0e1],
            n_rotations: 1,
            direction: None,
            final_cross: false,
            v: [0x3fe6_666a_65e3_9fe1, 0x3fff_ffff_28df_2819],
        },
    ];

    const TURNS: &[Move] = &[
        // `026`, depth 1: north, to the left.
        Move {
            o: [0x3ff9_2492_4924_924a, 0x3ff2_4924_9249_2492],
            d: [0xbfd2_4924_9249_2492, 0xbfc2_4924_9249_2492],
            n_rotations: 1,
            direction: Some([
                [0x3fc2_4924_9249_2492, 0x3fc2_4924_9249_2492],
                [0x3fc2_4924_9249_2492, 0x0000_0000_0000_0000],
            ]),
            final_cross: true,
            v: [0x3ff0_0000_0000_0000, 0x3fe6_db6d_b6db_6db4],
        },
        // `111`, depth 1: south, to the right.
        Move {
            o: [0x4003_6db6_db6d_b6db, 0x400d_b6db_6db6_db6c],
            d: [0xbfd2_4924_9249_2492, 0x3fd2_4924_9249_2492],
            n_rotations: 1,
            direction: Some([
                [0x3fc2_4924_9249_2492, 0x3fc2_4924_9249_2492],
                [0x3fc2_4924_9249_2492, 0x0000_0000_0000_0000],
            ]),
            final_cross: true,
            v: [0x4007_ffff_ffff_fffe, 0x4013_6db6_db6d_b6dc],
        },
        // `001`, depth 1: north, to the right.
        Move {
            o: [0x4010_0000_0000_0000, 0x4009_2492_4924_924a],
            d: [0x8000_0000_0000_0000, 0x3fd2_4924_9249_2492],
            n_rotations: 1,
            direction: Some([
                [0x3fc2_4924_9249_2492, 0x0000_0000_0000_0000],
                [0x3fc2_4924_9249_2492, 0x3fc2_4924_9249_2492],
            ]),
            final_cross: true,
            v: [0x4012_4924_9249_2492, 0x4010_0000_0000_0000],
        },
        // `084`, depth 1: south, to the left.
        Move {
            o: [0x4002_4924_9249_2493, 0x400a_4924_9249_2490],
            d: [0xbfd2_4924_9249_2492, 0x0000_0000_0000_0000],
            n_rotations: 1,
            direction: Some([
                [0x3fc2_4924_9249_2492, 0x3fc2_4924_9249_2492],
                [0x0000_0000_0000_0000, 0x3fc2_4924_9249_2492],
            ]),
            final_cross: true,
            v: [0x3ffb_6db6_db6d_b6e0, 0x4008_0000_0000_0001],
        },
    ];

    const FINAL: &[Move] = &[
        // `00`, depth 3: the offset ends on an interruption and the last crossing is made (north,
        // to the left).
        Move {
            o: [0x3ff8_d0fa_c687_d635, 0x3ff1_4e5e_0a72_f056],
            d: [0x3fc4_e5e0_a72f_0539, 0xbfb4_e5e0_a72f_053a],
            n_rotations: 1,
            direction: None,
            final_cross: true,
            v: [0x3ff0_0000_0000_0002, 0x3fd2_4924_9249_2490],
        },
        // `074`, depth 1: the offset ends on an interruption and the last crossing is made (south,
        // to the right).
        Move {
            o: [0x3fe2_4924_9249_249a, 0x3ffd_b6db_6db6_db6d],
            d: [0x0000_0000_0000_0000, 0x3fc2_4924_9249_2492],
            n_rotations: 1,
            direction: None,
            final_cross: true,
            v: [0x3fef_ffff_ffff_fffe, 0x4003_6db6_db6d_b6da],
        },
        // `005`, depth 1: the offset ends on an interruption and the last crossing is made (north,
        // to the right).
        Move {
            o: [0x4013_6db6_db6d_b6db, 0x4011_2492_4924_924a],
            d: [0x3fc2_4924_9249_2492, 0x0000_0000_0000_0000],
            n_rotations: 1,
            direction: None,
            final_cross: true,
            v: [0x4016_db6d_b6db_6db6, 0x4014_0000_0000_0000],
        },
        // `11`, depth 2: the offset ends on an interruption and the last crossing is made (south,
        // to the left).
        Move {
            o: [0x400b_6db6_db6d_b6dc, 0x4013_6db6_db6d_b6dc],
            d: [0xbfdb_6db6_db6d_b6db, 0x0000_0000_0000_0000],
            n_rotations: 1,
            direction: None,
            final_cross: true,
            v: [0x4001_2492_4924_9248, 0x4010_0000_0000_0000],
        },
        // `011`, depth 1: the first centroid, with no last crossing asked, stops within the
        // rhombus.
        Move {
            o: [0x4012_79e7_9e79_e79e, 0x4012_1861_8618_6186],
            d: [0x3fb8_6186_1861_8618, 0x3fa8_6186_1861_8618],
            n_rotations: 1,
            direction: None,
            final_cross: false,
            v: [0x4012_db6d_b6db_6db6, 0x4012_4924_9249_2492],
        },
        // `001`, depth 1: the first centroid, with no last crossing asked, stops on the
        // interruption it reaches.
        Move {
            o: [0x400f_9e79_e79e_79e8, 0x4009_8618_6186_1862],
            d: [0x3fa8_6186_1861_8618, 0xbfa8_6186_1861_8618],
            n_rotations: 1,
            direction: None,
            final_cross: false,
            v: [0x4010_0000_0000_0000, 0x4009_2492_4924_924a],
        },
        // `001`, depth 1, in the test for a sub-zone: no last crossing asked, and an interruption
        // crossed on the way (north, to the left).
        Move {
            o: [0x4012_db6d_b6db_6db7, 0x4010_0000_0000_0000],
            d: [0xbfa8_231b_cb56_4efe, 0xbfb8_231b_cb56_4efe],
            n_rotations: 1,
            direction: None,
            final_cross: false,
            v: [0x400f_3ee7_21a5_4d88, 0x4009_e898_231b_cb56],
        },
        // `0704`, depth 1, in the test for a sub-zone: no last crossing asked, and an interruption
        // crossed on the way (south, to the right).
        Move {
            o: [0x3fe8_29cb_c14e_5e0b, 0x3fff_ac68_7d63_43eb],
            d: [0x3f9b_95d6_9f3e_1122, 0x3fa1_3da6_2386_cab5],
            n_rotations: 1,
            direction: None,
            final_cross: false,
            v: [0x3ff0_3655_ae7f_7a41, 0x4001_d98c_39ad_a97c],
        },
        // `0003`, depth 1, in the test for a sub-zone: no last crossing asked, and an interruption
        // crossed on the way (north, to the right).
        Move {
            o: [0x4007_d634_3eb1_a1f6, 0x4001_cbc1_4e5e_0a73],
            d: [0x3fa1_3da6_2386_cab5, 0x3f7b_95d6_9f3e_1122],
            n_rotations: 1,
            direction: None,
            final_cross: false,
            v: [0x400e_419e_9d92_13a5, 0x4008_1b2a_d73f_bd21],
        },
        // `074`, depth 1, in the test for a sub-zone: no last crossing asked, and an interruption
        // crossed on the way (south, to the left).
        Move {
            o: [0x3ff0_0000_0000_0000, 0x4004_9249_2492_4924],
            d: [0xbfb8_231b_cb56_4efe, 0xbfa8_231b_cb56_4efe],
            n_rotations: 1,
            direction: None,
            final_cross: false,
            v: [0x3fd8_6953_6202_ed00, 0x3ffe_7dce_434a_9b10],
        },
    ];

    const LINES: &[Move] = &[
        // `07`, depth 1: on a horizontal line, moving up, where the nudge of 2e-11 settles the
        // rhombus.
        Move {
            o: [0x3fe2_4924_9249_2494, 0x3ffb_6db6_db6d_b6dd],
            d: [0x3feb_6db6_db6d_b6db, 0x3fe2_4924_9249_2492],
            n_rotations: 1,
            direction: None,
            final_cross: true,
            v: [0x3ff6_db6d_b6db_6db8, 0x4002_4924_9249_2493],
        },
        // `015`, depth 1: on a horizontal line, moving down, where the nudge of 2e-11 settles the
        // rhombus.
        Move {
            o: [0x4013_6db6_db6d_b6db, 0x4014_0000_0000_0000],
            d: [0xbfd2_4924_9249_2492, 0xbfc2_4924_9249_2492],
            n_rotations: 1,
            direction: Some([
                [0x3fc2_4924_9249_2492, 0x3fc2_4924_9249_2492],
                [0x3fc2_4924_9249_2492, 0x3fc2_4924_9249_2492],
            ]),
            final_cross: true,
            v: [0x4012_4924_9249_2492, 0x4013_6db6_db6d_b6db],
        },
        // `01000000000000003`, depth 1: on a horizontal line, moving down, where a nudge of 1e-11
        // would not do.
        Move {
            o: [0x3e8f_0b04_ec9e_dd69, 0x3e7f_0b04_ee9e_dd69],
            d: [0xbe97_4843_b1ee_4c1e, 0xbe87_4843_b1ee_4c1e],
            n_rotations: 1,
            direction: Some([
                [0x3e87_4843_b1ee_4c1e, 0x3e87_4843_b1ee_4c1e],
                [0x3e87_4843_b1ee_4c1e, 0x3e87_4843_b1ee_4c1e],
            ]),
            final_cross: true,
            v: [0x4013_ffff_f83d_3ec4, 0x4013_ffff_fc1e_9f63],
        },
        // `014`, depth 1: on a vertical line, moving right, where the nudge of 2e-11 settles the
        // rhombus.
        Move {
            o: [0x4013_6db6_db6d_b6db, 0x4015_2492_4924_9249],
            d: [0x3fd2_4924_9249_2492, 0x3fd2_4924_9249_2492],
            n_rotations: 1,
            direction: None,
            final_cross: true,
            v: [0x3fc2_4924_9249_2484, 0x3fe2_4924_9249_2491],
        },
        // `01000000000000003`, depth 1, in the test for a sub-zone: on a vertical line, moving
        // right, where a nudge of 1e-11 would not do.
        Move {
            o: [0x0000_0000_0000_0000, 0x0000_0000_0000_0000],
            d: [0x3e6e_bb8c_8eb5_6eb7, 0xbe6e_bb8c_8eb5_6eb7],
            n_rotations: 1,
            direction: None,
            final_cross: false,
            v: [0x4013_ffff_fc28_8e6e, 0x4013_ffff_f851_1cdc],
        },
        // `073`, depth 1: on a vertical line, moving left, where the nudge of 2e-11 settles the
        // rhombus.
        Move {
            o: [0x3ff0_0000_0000_0000, 0x3ff6_db6d_b6db_6db6],
            d: [0xbfd2_4924_9249_2492, 0x0000_0000_0000_0000],
            n_rotations: 1,
            direction: Some([
                [0x3fc2_4924_9249_2492, 0x3fc2_4924_9249_2492],
                [0x3fc2_4924_9249_2492, 0x3fc2_4924_9249_2492],
            ]),
            final_cross: true,
            v: [0x3fe6_db6d_b6db_6db7, 0x3ff6_db6d_b6db_6db6],
        },
        // `01`, depth 2: on a horizontal line and away from its ends, where the diagonal nudge is
        // withheld.
        Move {
            o: [0x4012_4924_9249_2493, 0x4013_ffff_ffff_ffff],
            d: [0x3fdb_6db6_db6d_b6db, 0x0000_0000_0000_0000],
            n_rotations: 1,
            direction: None,
            final_cross: true,
            v: [0x4014_0000_0000_0001, 0x4014_0000_0000_0000],
        },
        // `023`, depth 1: off the lines, in a northern rhombus, where the diagonal nudge settles
        // the rhombus.
        Move {
            o: [0x3ff0_0000_0000_0000, 0x3fdb_6db6_db6d_b6da],
            d: [0xbfd2_4924_9249_2492, 0x0000_0000_0000_0000],
            n_rotations: 1,
            direction: Some([
                [0x3fc2_4924_9249_2492, 0x3fc2_4924_9249_2492],
                [0x3fc2_4924_9249_2492, 0x3fc2_4924_9249_2492],
            ]),
            final_cross: true,
            v: [0x3fe6_db6d_b6db_6db7, 0x3fdb_6db6_db6d_b6da],
        },
        // `074`, depth 1: off the lines, in a southern rhombus, where the diagonal nudge settles
        // the rhombus.
        Move {
            o: [0x3feb_6db6_db6d_b6e4, 0x4000_0000_0000_0000],
            d: [0x0000_0000_0000_0000, 0x3fd2_4924_9249_2492],
            n_rotations: 1,
            direction: None,
            final_cross: true,
            v: [0x3ff4_9249_2492_4924, 0x4003_6db6_db6d_b6d9],
        },
    ];

    const WRAPS: &[Move] = &[
        // `01`, depth 1: the origin lies before the layout.
        Move {
            o: [0xbfd5_5555_5555_5555, 0xbfe5_5555_5555_5555],
            d: [0xbfce_79e7_9e79_e79e, 0xbfa8_6186_1861_8618],
            n_rotations: 1,
            direction: None,
            final_cross: false,
            v: [0x4011_b6db_6db6_db6e, 0x4011_2492_4924_9249],
        },
        // `014`, depth 1: the nudged point lies before the layout, and the walk is brought forward
        // by a period.
        Move {
            o: [0x3fd2_4924_9249_2492, 0x3fd2_4924_9249_2494],
            d: [0xbfdb_6db6_db6d_b6db, 0x0000_0000_0000_0000],
            n_rotations: 1,
            direction: Some([
                [0x3fc2_4924_9249_2492, 0x3fc2_4924_9249_2492],
                [0x3fc2_4924_9249_2492, 0x3fc2_4924_9249_2492],
            ]),
            final_cross: true,
            v: [0x4013_6db6_db6d_b6db, 0x4015_2492_4924_9249],
        },
        // `005`, depth 1: the nudged point lies beyond the layout, and the walk is brought back by
        // a period.
        Move {
            o: [0x4016_4924_9249_2492, 0x4014_0000_0000_0000],
            d: [0x3fc2_4924_9249_2492, 0x3fc2_4924_9249_2492],
            n_rotations: 1,
            direction: None,
            final_cross: true,
            v: [0x3fe6_db6d_b6db_6db4, 0x3fc2_4924_9249_2492],
        },
        // `064`, depth 1: the point reached lies before the layout, to the left, and is brought
        // forward by a period.
        Move {
            o: [0x3fd2_4924_9249_2492, 0x3ff4_9249_2492_4925],
            d: [0xbfd2_4924_9249_2492, 0x0000_0000_0000_0000],
            n_rotations: 1,
            direction: Some([
                [0x3fc2_4924_9249_2492, 0x3fc2_4924_9249_2492],
                [0x0000_0000_0000_0000, 0x3fc2_4924_9249_2492],
            ]),
            final_cross: true,
            v: [0x4012_db6d_b6db_6db7, 0x4018_0000_0000_0000],
        },
        // `016`, depth 1: the point reached lies below the layout and is brought forward by a
        // period.
        Move {
            o: [0x3fe2_4924_9249_2493, 0x3fc2_4924_9249_2490],
            d: [0xbfd2_4924_9249_2492, 0xbfc2_4924_9249_2492],
            n_rotations: 1,
            direction: Some([
                [0x3fc2_4924_9249_2492, 0x3fc2_4924_9249_2492],
                [0x3fc2_4924_9249_2492, 0x0000_0000_0000_0000],
            ]),
            final_cross: true,
            v: [0x4014_0000_0000_0000, 0x4012_db6d_b6db_6db7],
        },
        // `064`, depth 1: the point reached lies beyond the layout and is brought back by a period.
        Move {
            o: [0x4012_4924_9249_2492, 0x4017_6db6_db6d_b6db],
            d: [0x0000_0000_0000_0000, 0x3fc2_4924_9249_2492],
            n_rotations: 1,
            direction: None,
            final_cross: true,
            v: [0x0000_0000_0000_0000, 0x3ff6_db6d_b6db_6db8],
        },
    ];

    const CORNERS: &[Move] = &[
        // `00000`, depth 3: at a corner, moving along x alone, crossed north, to the right.
        Move {
            o: [0x3fff_8895_4569_3c76, 0x3fef_ffff_ffff_fffa],
            d: [0x3f9d_daae_a5b0_e2e4, 0x0000_0000_0000_0000],
            n_rotations: 1,
            direction: None,
            final_cross: true,
            v: [0x4008_0000_0000_0002, 0x4000_0000_0000_0000],
        },
        // `110`, depth 1: at a corner, moving along x alone, crossed south, to the left.
        Move {
            o: [0x4009_2492_4924_9248, 0x4014_0000_0000_0000],
            d: [0xbfc2_4924_9249_2492, 0x0000_0000_0000_0000],
            n_rotations: 1,
            direction: None,
            final_cross: true,
            v: [0x4000_0000_0000_0000, 0x4010_0000_0000_0000],
        },
        // `04`, depth 2: at a corner, moving along x alone, where no interruption lies.
        Move {
            o: [0x4004_9249_2492_4924, 0x4008_0000_0000_0000],
            d: [0x3fdb_6db6_db6d_b6db, 0x0000_0000_0000_0000],
            n_rotations: 1,
            direction: None,
            final_cross: true,
            v: [0x4007_ffff_ffff_ffff, 0x4008_0000_0000_0000],
        },
        // `00`, depth 1: forward in x and back in y at once, crossed north, to the left.
        Move {
            o: [0x3ffb_6db6_db6d_b6dd, 0x3ff2_4924_9249_2496],
            d: [0x3fd2_4924_9249_2492, 0xbfc2_4924_9249_2492],
            n_rotations: 1,
            direction: None,
            final_cross: true,
            v: [0x3ff0_0000_0000_0004, 0x0000_0000_0000_0000],
        },
        // `11`, depth 1: back in x and forward in y at once, crossed south, to the left.
        Move {
            o: [0x400a_4924_9249_2492, 0x4013_6db6_db6d_b6dc],
            d: [0xbfd2_4924_9249_2492, 0x3fc2_4924_9249_2492],
            n_rotations: 1,
            direction: None,
            final_cross: true,
            v: [0x4000_0000_0000_0000, 0x4010_0000_0000_0001],
        },
        // `013`, depth 1: the rhombus entered is read an epsilon ahead of the point, in x, moving
        // right.
        Move {
            o: [0x4012_db6d_b6db_6db6, 0x4011_2492_4924_9248],
            d: [0x3fd2_4924_9249_2492, 0x3fd2_4924_9249_2492],
            n_rotations: 1,
            direction: None,
            final_cross: true,
            v: [0x4015_b6db_6db6_db6f, 0x4013_ffff_ffff_ffff],
        },
        // `094`, depth 1: the rhombus entered is read an epsilon ahead of the point, in x, moving
        // left.
        Move {
            o: [0x400a_4924_9249_2493, 0x4011_2492_4924_9249],
            d: [0xbfd2_4924_9249_2492, 0x0000_0000_0000_0000],
            n_rotations: 1,
            direction: Some([
                [0x3fc2_4924_9249_2492, 0x3fc2_4924_9249_2492],
                [0x0000_0000_0000_0000, 0x3fc2_4924_9249_2492],
            ]),
            final_cross: true,
            v: [0x4005_b6db_6db6_db6e, 0x4010_0000_0000_0000],
        },
        // `081`, depth 1: the rhombus entered is read an epsilon ahead of the point, in y, moving
        // up.
        Move {
            o: [0x3ff4_9249_2492_4925, 0x4004_9249_2492_4924],
            d: [0x3fdb_6db6_db6d_b6db, 0x3fdb_6db6_db6d_b6db],
            n_rotations: 1,
            direction: None,
            final_cross: true,
            v: [0x3fff_ffff_ffff_fffe, 0x400a_4924_9249_2492],
        },
        // `036`, depth 1: the rhombus entered is read an epsilon ahead of the point, in y, moving
        // down.
        Move {
            o: [0x4004_9249_2492_4925, 0x4001_2492_4924_9248],
            d: [0xbfd2_4924_9249_2492, 0xbfc2_4924_9249_2492],
            n_rotations: 1,
            direction: Some([
                [0x3fc2_4924_9249_2492, 0x3fc2_4924_9249_2492],
                [0x3fc2_4924_9249_2492, 0x0000_0000_0000_0000],
            ]),
            final_cross: true,
            v: [0x4000_0000_0000_0000, 0x3ffb_6db6_db6d_b6d6],
        },
    ];

    const SEVERAL: &[Move] = &[
        // `001`, depth 1: two crossings.
        Move {
            o: [0x400e_db6d_b6db_6db7, 0x400a_4924_9249_2493],
            d: [0x3fdb_6db6_db6d_b6db, 0x0000_0000_0000_0000],
            n_rotations: 1,
            direction: None,
            final_cross: true,
            v: [0x4016_db6d_b6db_6db7, 0x4013_ffff_ffff_ffff],
        },
        // `00`, depth 3: three crossings.
        Move {
            o: [0x3ff8_d0fa_c687_d635, 0x3ff1_4e5e_0a72_f056],
            d: [0x3fdc_bc14_e5e0_a72e, 0xbfcc_bc14_e5e0_a730],
            n_rotations: 1,
            direction: None,
            final_cross: true,
            v: [0x4010_0000_0000_0000, 0x4009_2492_4924_9248],
        },
        // `00`, depth 3: three crossings in four turns of the loop.
        Move {
            o: [0x3ff8_7d63_43eb_1a1d, 0x3ff0_5397_829c_bc14],
            d: [0x3fdf_58d0_fac6_87d6, 0xbfcf_58d0_fac6_87d7],
            n_rotations: 1,
            direction: None,
            final_cross: true,
            v: [0x400f_d634_3eb1_a1f7, 0x4009_cbc1_4e5e_0a75],
        },
    ];

    const TWO_TURNS: &[Move] = &[
        // `11`, depth 2: crossed south, to the left, and what is left of the offset turned twice.
        Move {
            o: [0x400b_6db6_db6d_b6dc, 0x4014_0000_0000_0000],
            d: [0xbfe2_4924_9249_2492, 0x0000_0000_0000_0000],
            n_rotations: 2,
            direction: None,
            final_cross: true,
            v: [0x4000_0000_0000_0000, 0x400e_db6d_b6db_6db8],
        },
        // `1100`, depth 2: crossed south, to the left, and what is left of the offset turned twice.
        Move {
            o: [0x4008_7d63_43eb_1a20, 0x4014_0000_0000_0000],
            d: [0xbfb4_e5e0_a72f_0539, 0x0000_0000_0000_0000],
            n_rotations: 2,
            direction: None,
            final_cross: true,
            v: [0x4000_0000_0000_0000, 0x400f_d634_3eb1_a1f6],
        },
        // `000000`, depth 2: two turns asked at a last crossing, where nothing is left to turn
        // (north, to the left).
        Move {
            o: [0x3fff_dc2c_c805_f88a, 0x3fef_ffff_ffff_fffc],
            d: [0x3f67_e225_515a_4f1d, 0x0000_0000_0000_0000],
            n_rotations: 2,
            direction: None,
            final_cross: true,
            v: [0x3fef_ffff_ffff_fffc, 0x3f77_e225_515a_4e00],
        },
        // `00`, depth 2: two turns asked at a last crossing, where nothing is left to turn (north,
        // to the right).
        Move {
            o: [0x3ff9_2492_4924_924c, 0x3ff2_4924_9249_2497],
            d: [0x3fdb_6db6_db6d_b6db, 0x0000_0000_0000_0000],
            n_rotations: 2,
            direction: None,
            final_cross: true,
            v: [0x4006_db6d_b6db_6db4, 0x4000_0000_0000_0000],
        },
    ];

    const POLAR: &[Move] = &[
        // `00`, depth 1: the northern polar pentagon of level 0.
        Move {
            o: [0x4013_6db6_db6d_b6dc, 0x4012_4924_9249_2492],
            d: [0x3fe2_4924_9249_2492, 0xbfd2_4924_9249_2492],
            n_rotations: 1,
            direction: Some([
                [0xbfdb_6db6_db6d_b6db, 0xbfd2_4924_9249_2492],
                [0x3fd2_4924_9249_2492, 0xbfc2_4924_9249_2492],
            ]),
            final_cross: true,
            v: [0x3ffb_6db6_db6d_b6dd, 0x3ff2_4924_9249_2496],
        },
        // `11`, depth 1: the southern polar pentagon of level 0.
        Move {
            o: [0x3fc2_4924_9249_2492, 0x3ff6_db6d_b6db_6db7],
            d: [0xbfe2_4924_9249_2492, 0x3fd2_4924_9249_2492],
            n_rotations: 1,
            direction: Some([
                [0x3fdb_6db6_db6d_b6db, 0x3fd2_4924_9249_2492],
                [0xbfd2_4924_9249_2492, 0x3fc2_4924_9249_2492],
            ]),
            final_cross: true,
            v: [0x400a_4924_9249_2492, 0x4013_6db6_db6d_b6dc],
        },
        // `000`, depth 2: the northern polar pentagon of level 1.
        Move {
            o: [0x3fea_1f58_d0fa_c687, 0x3f94_e5e0_a72f_0568],
            d: [0x3fd8_d0fa_c687_d634, 0x3fc4_e5e0_a72f_0539],
            n_rotations: 1,
            direction: Some([
                [0xbf94_e5e0_a72f_053a, 0xbfaf_58d0_fac6_87d6],
                [0x3faf_58d0_fac6_87d6, 0x3fa4_e5e0_a72f_0539],
            ]),
            final_cross: true,
            v: [0x4006_87d6_343e_b1a1, 0x4000_29cb_c14e_5e0a],
        },
        // `110`, depth 2: the southern polar pentagon of level 1.
        Move {
            o: [0x4010_bc14_e5e0_a72f, 0x4017_eb1a_1f58_d0fb],
            d: [0xbfca_1f58_d0fa_c687, 0xbfc7_829c_bc14_e5e0],
            n_rotations: 1,
            direction: Some([
                [0x3f94_e5e0_a72f_053a, 0x3faf_58d0_fac6_87d6],
                [0xbfa4_e5e0_a72f_0539, 0x3f94_e5e0_a72f_053a],
            ]),
            final_cross: true,
            v: [0x4009_7829_cbc1_4e5d, 0x4013_eb1a_1f58_d0fb],
        },
    ];

    /// No interruption is met: a step, no step at all, a step from one rhombus into the next,
    /// and a step cut short at the edge of its rhombus.
    #[test]
    fn move5x6_within_the_layout_matches_the_engine() {
        assert_moves_as_the_engine("within the layout", WITHIN);
    }

    /// The four crossings, each performed as the library compiles it: site MV-1 north, to the
    /// left; MV-2 and MV-3 south, to the right; MV-4 and MV-5 north, to the right; MV-6 and MV-7
    /// south, to the left. All seven sites are inert for origins inside the layout, which are
    /// the arguments the sub-zone generator passes: with any one of them put back in the eC's
    /// order the function answers every vector of this module as before, as it did every call of
    /// the run described at [`Move`]. Outside the layout a site may move an answer, as MV-1 does
    /// from the origin given at [`move5x6`]. They are pinned by their constants, in
    /// `the_folded_constants_match_the_engines_rodata`, and performed as compiled all the same.
    #[test]
    fn move5x6_crosses_each_interruption_as_the_engine_does() {
        assert_moves_as_the_engine("a crossing", CROSSINGS);
    }

    /// The caller's direction is turned with what is left of the offset: sixty degrees clockwise
    /// at a crossing to the right in the north or to the left in the south, anticlockwise at
    /// the other two.
    #[test]
    fn move5x6_turns_the_callers_direction_at_a_crossing() {
        assert!(
            TURNS
                .iter()
                .all(|m| m.direction.is_some_and(|[given, left]| given != left))
        );
        assert_moves_as_the_engine("a crossing with a direction", TURNS);
    }

    /// With `final_cross` an offset that ends on an interruption is carried across it; without,
    /// the walk stops where the offset is spent, and still crosses what lies on its way.
    #[test]
    fn move5x6_makes_the_last_crossing_only_where_it_is_asked() {
        assert!(FINAL.iter().any(|m| m.final_cross) && FINAL.iter().any(|m| !m.final_cross));
        assert_moves_as_the_engine("the last crossing", FINAL);
    }

    /// The point by which the rhombus is chosen is nudged off the lines of the layout: by
    /// `2e-11` in the sense of the offset where the walk stands on a line, and diagonally by
    /// `1e-11` elsewhere (`ri5x6.ec:1749-1776`).
    #[test]
    fn move5x6_chooses_the_rhombus_by_the_nudged_point() {
        assert_moves_as_the_engine("a line of the layout", LINES);
    }

    /// The folds by a period of five at the two ends of the layout: of the origin, of the nudged
    /// point, and of the point reached.
    #[test]
    fn move5x6_wraps_at_the_ends_of_the_layout() {
        assert_moves_as_the_engine("an end of the layout", WRAPS);
    }

    /// Which rhombus is entered (`ri5x6.ec:1823-1837`): by rule at a corner, and elsewhere read
    /// an epsilon ahead of the point reached.
    #[test]
    fn move5x6_settles_the_rhombus_entered_as_the_engine_does() {
        assert_moves_as_the_engine("the rhombus entered", CORNERS);
    }

    #[test]
    fn move5x6_crosses_several_interruptions_in_one_call() {
        assert_moves_as_the_engine("several crossings", SEVERAL);
    }

    /// Two turns a crossing, which the generator asks at a pentagon of an even level and at an
    /// even depth (`RI7H.ec:3692`). Among the calls sampled only those of the southern polar
    /// pentagon tell two turns from one; at the others the crossing is the last, and nothing is
    /// left to turn.
    #[test]
    fn move5x6_turns_twice_a_crossing_where_it_is_asked() {
        assert!(TWO_TURNS.iter().all(|m| m.n_rotations == 2));
        assert_moves_as_the_engine("two turns a crossing", TWO_TURNS);
    }

    #[test]
    fn move5x6_at_the_polar_pentagons_matches_the_engine() {
        assert_moves_as_the_engine("a polar pentagon", POLAR);
    }

    /// The call that the eC's own test sets aside with the remark that it loops without end
    /// (`ri5x6.ec:1940-1941`). The engine's function does not return from it, so it is given to
    /// this crate alone, which gives the walk up at its bound; and the bound is what stops it,
    /// since the walk is no nearer its end with a thousand times as many turns.
    #[test]
    fn move5x6_gives_up_the_walk_that_the_ecs_own_test_sets_aside_as_endless() {
        let o = (4.7755102040816331, 5.8979591836734695);
        let (dx, dy) = (0.1836734693877551, 0.1224489795918367);
        assert_eq!(move5x6(o, dx, dy, 2, None, true), None);
        let long = MOVE5X6_MAX_TURNS * 1000;
        assert_eq!(move5x6_within(o, dx, dy, 2, None, true, long), None);
        // The same walk with one turn a crossing ends, at once, and at the engine's point.
        assert_eq!(
            move5x6(o, dx, dy, 1, None, true).map(bits),
            Some((0x3f94_e5e0_a72f_0550, 0x3ff0_fac6_87d6_343d))
        );
    }

    /// The integer sums of the crossings wrap, as the library's `lea` does. The origin is far
    /// beyond the layout, where no caller of the function stands, and the engine answers it all
    /// the same: the walk comes to the crossing to the left in the south at `x` = 2^31 - 1,
    /// where `(int)(c.x + 1 + 1E-11)` (`ri5x6.ec:1870`) is beyond an `int` and the library's
    /// conversion gives the least `int`; `iy - 1` (`ri5x6.ec:1871`) is then the greatest
    /// (`0x39742`). The answer is the engine's own function's, read at its address in the
    /// library.
    #[test]
    fn move5x6_wraps_its_integer_sums_as_the_library_does() {
        assert_eq!(
            move5x6((2_147_483_653.0, 5.0), -1.0, 0.0, 1, None, true).map(bits),
            Some((0xbff0_0000_0000_0000, 0x41e0_0000_0000_0000))
        );
    }
}
