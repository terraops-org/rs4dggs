//! Aperture-7 hexagonal topology on the 5x6 ISEA planar grid.
//!
//! This is the topology half of IGEO7/ISEA7H. It turns a projection's planar point
//! `(x, y)` in DGGAL's oblique 5x6 grid into a Z7 address (base cell plus direction
//! digits), rebuilds planar cell geometry (a centroid and five or six vertices)
//! from that address, and lists a zone's neighbours as the engine finds them, in the
//! same plane; the sibling module `hex_a7_subzones` orders its sub-zones. The geographic
//! lift belongs to the `Grid`, so this module stays independent of any particular
//! projection above the shared 5x6 layout.
//!
//! Provenance, cited inline against these sources:
//!   - `dggrs/RI7H.ec`: the I7H zone, the `fromCentroid` quantisation, the
//!     `centroid`/`getVertices` geometry, the neighbour search `getNeighbors`, the
//!     5x6 wrapping and interruption helpers, `isEdgeHex`, `calcCandidateParent`,
//!     `getPrimaryChildren`.
//!   - `dggrs/RI7H_Z7.ec`: the I7H to Z7 conversion, `from7H`/`to7H`,
//!     `getLevelRotationOffset`, `getChildPosition`, the pentagon child-position
//!     adjustments and the Z7 wrapper of `getZoneNeighbors`.
//!
//! The eC is the ground truth this module follows. py4dggs's `topologies/hex_a7.py`
//! serves as its readable guide: most items here correspond to a Python item, in the
//! same order, so that the two can be read side by side, and wherever the guide agrees
//! with the eC the arithmetic is held to it bit for bit, with the same expressions,
//! the same groupings, the same epsilons, the same evaluation order and the same
//! half-up `jsround`. Re-associating a single floating-point operation or merging a
//! branch would drift the last bit and fail the differential comparison against DGGAL.
//!
//! Where the guide contradicts the eC or the live engine, this module departs from it
//! and follows the engine, and a `// faithful:` comment at the site says so. The
//! departures lie mostly in the quantisation. `contains_point` carries the eC's
//! bounding boxes and compares against zero, where py4dggs runs a convex half-plane
//! test against an epsilon that is too coarse from resolution 13 upwards; the
//! odd-level candidate search answers the engine's null zone where it finds no cell,
//! rather than py4dggs's fallback one level finer; and `from_7h` refuses a parent
//! chain shorter than the zone's level, which py4dggs indexes past. Two boundary
//! tracers, `polar_pentagon_outline` and `add_non_polar_vertices_refined`, have no
//! counterpart in the guide at all: they reproduce the polygons the eC's own
//! containment test is given. The quantisation that results is checked against a
//! live DGGAL at every resolution.
//!
//! One class of departure goes further, and against the eC's own text. DGGAL v0.0.6
//! is compiled with `-O2 -ffast-math` (`Makefile.dggal:133`), and the eC marks only
//! some of its functions with the `-fno-unsafe-math-optimizations` attribute; in the
//! ones it does not mark, gcc was free to re-associate, and it did. At the sites
//! below this module performs the operation order read out of the shipped
//! `libdggal.so`, BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, rather than the
//! source's, and every one carries a `// faithful:` comment quoting its own address
//! and the eC line it replaces.
//!
//! The criterion is not uniform across them, and the difference matters when one is
//! read. The side test and the boundary tracing were ported because the rewriting
//! demonstrably changes answers: without them the quantiser disagreed with the engine
//! at thousands of boundary points. The five wrap comparisons were ported on the
//! principle that in a function the eC does not guard the shipped engine decides, and
//! not the eC's text: they are there so that the module performs what the engine
//! performs, not because an answer turns on them. Three of the five were afterwards measured to be
//! inert in place, on censuses that reach them: the bounding box on the `-5` path,
//! taken 2,031, 2,706 and 1,450 times in the three censuses of the boundary work, and
//! the corner test and the ordinate in `from_centroid`, each reached and taken 6,800
//! times in the census built for the corner, with no answer moved in any of them. The
//! remaining two, the same pair in `from_even_level_primary_child`, no census reaches
//! at all; they rest on the argument given at the identical sites in `from_centroid`,
//! that neither rewriting can change the integer it decides. The three rewritten
//! differences of the neighbour search were ported on the same principle, and measured
//! inert in the same way: over some 275,000 zones of each aperture-7 grid, where they are
//! taken some 900, 20 and 5,200 times, a variant in the source's order moved no
//! neighbour list, although it moves the point stepped to in its last bits. The sites
//! are:
//!
//! - the side test, `pointLineSide` (RI7H.ec:1360-1365), in `point_line_side` and
//!   `tie_break_side_is_negative`: it is asked for the sign of a quantity that is
//!   legitimately zero on a cell boundary, and the eC's grouping returns a tiny
//!   negative value there where the engine returns an exact zero;
//! - five wrap comparisons across three functions: the `-5` bounding box in
//!   `contains_point`, and the bottom-right corner test together with the truncation
//!   of the ordinate it wraps, the latter two appearing in both `from_centroid` and
//!   `from_even_level_primary_child`;
//! - the boundary tracing, in `add_non_polar_base_vertices`,
//!   `add_non_polar_vertices_refined`, `polar_pentagon_outline` and the polar fan of
//!   `get_vertices`, where the compiler cancelled the centroid out of the first
//!   direction, moved the whole number of each fan step inside the bracket and
//!   removed a pair of opposite rotations;
//! - the refined boundary of a polar pentagon, in `polar_pentagon_sides`, where the compiler
//!   moved the whole number of each turn about the pole inside the bracket, as in the outline
//!   above but not at every turn, and folded what it could of the two points of each
//!   interruption. Seven of its sites move a ring and are pinned by rings read from the engine;
//!   the others give the source's own doubles at every level, and a test records both. The
//!   walker of every other zone, `add_non_polar_vertices_refined`, serves the refined boundary
//!   too: the eC has it twice, and the library performs both copies in the same order;
//! - the centroid, in `zone_centroid`, where the divisions by seven became
//!   multiplications by the rounded reciprocal;
//! - the neighbour search, in `get_neighbors`, where the compiler formed each of the three
//!   offsets moved by five from the difference already taken, and moved the epsilon of the
//!   hemisphere test across the comparison.
//!
//! Everywhere else in this module the eC's source text still decides.

use super::hex_a7_subzones;
use crate::fivebysix::{
    Crossing, add_intermediate_points, cvtt_i32, move5x6_vertex, move5x6_vertex2, sgn,
};
use crate::interfaces::Topology;
use crate::interfaces::sealed::{Private, TopologyPlumbing};
use crate::math;
use crate::types::{Address, PlanarPoint};

/// The aperture-7 hexagonal topology on DGGAL's 5x6 ISEA planar grid.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct HexA7;

impl crate::interfaces::sealed::Sealed for HexA7 {}

// --------------------------------------------------------------------------- //
// I7H <-> Z7 mapping tables (RI7H_Z7.ec)
// --------------------------------------------------------------------------- //

/// I7H child position to Z7 digit.
const C_MAP: [i64; 7] = [0, 3, 1, 5, 4, 6, 2];
/// Z7 digit to I7H child position.
const INV_C_MAP: [i64; 7] = [0, 2, 6, 1, 4, 3, 5];
/// I7H root rhombus to Z7 base cell.
const ROOT_MAP: [u64; 12] = [1, 6, 2, 7, 3, 8, 4, 9, 5, 10, 0, 11];
/// Z7 base cell to I7H root rhombus.
const INV_ROOT_MAP: [i64; 12] = [10, 0, 2, 4, 6, 8, 1, 3, 5, 7, 9, 11];

/// Powers of seven, the count of cells along a rhombus edge at a given `l49r`.
/// The table covers every resolution the grid reaches; [`pow7`] falls back to a
/// rounded computation past it.
const POW7: [i64; 11] = [
    1, 7, 49, 343, 2401, 16807, 117649, 823543, 5764801, 40353607, 282475249,
];

/// The first Z7 level with no geometry. The 64-bit packing holds twenty direction
/// digits, so level 20 is representable, but `to7H` (RI7H_Z7.ec:348-352) "does not
/// support level 20 zones" and returns nullZone, which is why DGGAL's highest real
/// level is 19. The constant is kept here rather than taken from the indexing: it
/// is a property of the Z7 to 7H *geometry* conversion, which is this module's
/// concern. See [`HexA7::is_null_geometry`].
pub(super) const NULL_GEOMETRY_LEVEL: usize = 20;

// --------------------------------------------------------------------------- //
// Rounding and power helpers
// --------------------------------------------------------------------------- //

/// Rounds half towards positive infinity, that is `floor(x + 0.5)`, matching
/// JavaScript's `Math.round` and the eC. Rust's `f64::round` rounds half away from
/// zero and Python's built-in `round` is banker's rounding; either would pick a
/// different cell on an exact half-boundary, so the whole pipeline uses this. The eC writes it
/// inline as `(int64)(x * p + 0.5)` under a `Max(0, ...)`, in `fromCentroid`
/// (`RI7H.ec:1474-1475`) and in `fromEvenLevelPrimaryChild` (`RI7H.ec:1729-1730`).
fn jsround(x: f64) -> i64 {
    (x + 0.5).floor() as i64
}

/// `7**n`: a table lookup inside the table's range, a rounded computation beyond
/// it. Every caller passes a non-negative `l49r`, so the negative case cannot
/// arise; were it to, the cast below would send it to the fallback arm rather than
/// index out of bounds.
///
/// The fallback works in `f64` where Python's `round(7 ** n)` is exact at any size,
/// so it loses the last digits from `n == 19` and saturates the `i64` cast past
/// `n == 22`. The table covers `n` up to 10, which is `l49r` at Z7 level 20, one
/// past the deepest real level, so no resolution this crate admits reaches the
/// fallback at all. Were an absurd resolution ever to reach it, the saturated value
/// could make `calc_candidate_parent`'s `col + add_col` overflow and panic in a
/// debug build; that is a bound on the caller, not a case to handle here.
///
/// The eC's `POW7` (`RI7H.ec:34`) reads the same powers from `powersOf7` (`RI7H.ec:837-841`),
/// a longer table, and past its end adds `POW_EPSILON` to `pow(7, x)` rather than rounding.
pub(super) fn pow7(n: i64) -> i64 {
    if (n as usize) < POW7.len() {
        POW7[n as usize]
    } else {
        math::powi(7.0, n as i32).round_ties_even() as i64
    }
}

/// Python's two-argument `max` on floats: the second value wins only when it is
/// strictly greater, which keeps the first argument on a tie. `f64::max` is free to
/// return either operand when both are zero, so it is not used here.
///
/// It stands for the eC's `Max` in `crosses5x6InterruptionV2` (`ri5x6.ec:1645-1646`), which the
/// eC compiler expands as `(a > b) ? a : b` (`ectp/src/pass15.ec:9786-9796`), keeping the second
/// argument on a tie.
///
/// The library does not perform that text. At `0x38b00` (for `x`) and `0x38a98` (for `y`) in
/// `libdggal.so`, BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, gcc compiled it as
/// `(a <= b) ? b : a`, which keeps the first argument on NaN, as this function does. The two
/// differ only on a tie between zeros of opposite sign, and every caller passes a negative `b`,
/// so at every call site this function is the library's operation for every input.
fn pymax(a: f64, b: f64) -> f64 {
    if b > a { b } else { a }
}

/// The library's `Min` in `crosses5x6InterruptionV2` (`ri5x6.ec:1645-1646`): the first argument
/// when it is strictly less, and the second otherwise, so the second wins on a tie and whenever
/// either is NaN. gcc emitted each of the two as a single `minsd`, at `0x387d7` (for `x`) and
/// `0x387fa` (for `y`) in `libdggal.so`, BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, with
/// `a` in the destination and `b` in the source; `minsd` answers `dst < src ? dst : src`. That is
/// also the eC compiler's expansion of `Min(a, b)`, `(a < b) ? a : b`
/// (`ectp/src/pass15.ec:9786-9796`), so here the source text and the library agree.
///
/// The name is py4dggs's, kept so that the two files can be read side by side, but the function
/// is no longer Python's `min`, which keeps the first argument: the two part company on zeros of
/// opposite sign, where `pymin(0.0, -0.0)` is `-0.0`, and on NaN. Neither case reaches the call
/// sites, whose `b` is a component of a hexagon edge, never zero and never NaN.
fn pymin(a: f64, b: f64) -> f64 {
    if a < b { a } else { b }
}

// --------------------------------------------------------------------------- //
// The I7H zone and its derived predicates (RI7H.ec I7HZone)
// --------------------------------------------------------------------------- //

/// An I7H zone.
///
/// `l49r` is the aperture-49 resolution, which advances one step per two Z7
/// levels; `root` is the rhombus (0 to 9) or a polar sentinel (0xA north, 0xB
/// south); `row` and `col` index within the rhombus; `sub_hex` is the odd-level
/// child selector, zero at an even level. The Z7 level is
/// `2 * l49r + (sub_hex > 0)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) struct Z {
    pub(super) l49r: i64,
    pub(super) root: i64,
    row: i64,
    col: i64,
    pub(super) sub_hex: i64,
}

impl Z {
    /// No eC counterpart: the eC builds its zone with a brace initialiser at each site.
    fn new(l49r: i64, root: i64, row: i64, col: i64, sub_hex: i64) -> Z {
        Z {
            l49r,
            root,
            row,
            col,
            sub_hex,
        }
    }
}

/// The zone's Z7 level: two levels per aperture-49 step, plus one for a non-zero
/// sub-hex. Ports the `level` property (`RI7H.ec:889-892`).
pub(super) fn zone_level(z: &Z) -> i64 {
    2 * z.l49r + i64::from(z.sub_hex > 0)
}

/// The number of cell vertices: five for the rhombus-corner pentagon, six
/// otherwise. Ports the `nPoints` property (`RI7H.ec:894-904`).
pub(super) fn zone_npoints(z: &Z) -> i64 {
    if z.sub_hex > 1 {
        return 6;
    }
    if z.row == 0 && z.col == 0 { 5 } else { 6 }
}

/// Zone equality treating the absent zone as a distinct empty value, as py4dggs's
/// `zone_eq` does for `None`: two absent zones are equal, and an absent zone never
/// equals a present one. The derived `PartialEq` on `Option<Z>` already has exactly
/// that behaviour, so the body is a single comparison; the function is kept, and
/// called at every site the Python calls it, so that the correspondence between the
/// two files stays visible.
///
/// No eC counterpart: the eC compares zones with `==`; this is py4dggs's `zone_eq`.
fn zone_eq(a: Option<Z>, b: Option<Z>) -> bool {
    a == b
}

/// True for a hexagon lying on a rhombus seam, that is the north edge of a north
/// rhombus or the south edge of a south rhombus: the `isEdgeHex` property
/// (`RI7H.ec:906-930`). Edge hexagons need the rotation corrections applied during
/// child enumeration.
///
/// Every caller passes an even zone. The eC's parity test, `!(level & 0)` (`RI7H.ec:915`),
/// is always true, so the seam test applies at every level, and an odd zone is answered
/// by the row and column it shares with its even parent; this function does the same.
fn zone_is_edge_hex(z: &Z) -> bool {
    if zone_npoints(z) != 6 {
        return false;
    }
    if z.root & 1 != 0 {
        z.col == 0 // South
    } else {
        z.row == 0 // North
    }
}

// --------------------------------------------------------------------------- //
// 5x6 coordinate wrapping and interruption
// (the eC 5x6 helpers live in `projections/ri5x6.ec`, not in `dggrs/RI7H.ec`)
// --------------------------------------------------------------------------- //

/// Folds a 5x6 planar point into the canonical fundamental domain, handling the
/// polar diagonals and the rhombus interruptions. The epsilons disambiguate a point
/// sitting exactly on a seam.
///
/// This ports the *inline* domain fold that the eC performs inside `fromCentroid`
/// (RI7H.ec:1439-1458), and not the eC function that happens to carry the same name
/// (`canonicalize5x6`, ri5x6.ec:1562-1606), which is a different and more elaborate
/// routine: same name, different job. The engine applies that one to every planar
/// vertex in `getZoneWGS84Vertices` (RI7H.ec:461-462), before the inverse projection;
/// this crate ports it in `fivebysix.rs` for the aperture-3 grids but does not apply it
/// here, and `Grid::vertices` records why the omission changes no answer on the
/// aperture-7 grids.
///
/// The fold is inline except on its last arm, where the eC calls
/// `move5x6Vertex2(c, { 5, 5 }, c.x, c.y, false)` (RI7H.ec:1458, and again at :1719
/// inside the even-level parent search) rather than adding five to each coordinate.
/// The two agree on every input that reaches the arm, so the open-coded addition is
/// kept; ri5x6.ec:1470-1560 is the routine to read if that equivalence is ever
/// revisited.
fn canonicalize5x6(cx: f64, cy: f64) -> (f64, f64) {
    let (mut x, mut y) = (cx, cy);
    if (x - y - 1.0).abs() < 1e-10 {
        return (x, y); // north pole diagonal
    }
    if (y - x - 2.0).abs() < 1e-10 {
        return (x, y); // south pole diagonal
    }
    if y < -1e-11 && x > -1e-11 {
        x -= y;
        y = 0.0;
    } else if (x + 1e-11).floor() as i64 > (y + 1e-11).floor() as i64 {
        let iy = ((y + 1e-11).floor() as i64).min(5);
        x += iy as f64 + 1.0 - y;
        y = iy as f64 + 1.0;
    } else if ((y + 1e-11).floor() as i64) - ((x + 1e-11).floor() as i64) > 1 {
        let ix = ((x + 1e-11).floor() as i64).min(4);
        y += ix as f64 + 1.0 - x;
        x = ix as f64 + 1.0;
    } else if x < -1e-11 || y < -1e-11 {
        x += 5.0;
        y += 5.0;
    }
    (x, y)
}

/// The signed side of the point `(px, py)` relative to the directed line A to B;
/// positive means left. It is the half-plane test inside [`contains_point`], which
/// is the odd-level path of the quantiser.
///
/// faithful: this is not the eC's source order, and the reason is the compiler, not
/// the source. `pointLineSide` (RI7H.ec:1360-1365) reads `A*x + B*y + C` with
/// `A = dy`, `B = -dx` and `C = a.y*dx - dy*a.x`, and this crate evaluated it in
/// exactly that order until the order gcc actually emitted was read out of the
/// shipped library. DGGAL v0.0.6 is built with `-O2 -ffast-math`
/// (`Makefile.dggal:133`), and `containsPoint` (RI7H.ec:951) carries no
/// `__attribute__ ((optimize("-fno-unsafe-math-optimizations")))`, unlike the
/// eighteen functions in the tree that do, so the compiler was free to re-associate
/// the sum and it did. In `libdggal.so`, BuildID
/// `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, the inlining of the call at
/// RI7H.ec:1061 occupies `0x23241` to `0x23276` and computes the four products in
/// the source's own order but adds them in two halves:
///
/// ```text
///   23241  movapd xmm4,xmm2      ; b.x
///   23245  subsd  xmm3,xmm5      ; dy   = b.y - a.y
///   23249  subsd  xmm4,xmm1      ; dx   = b.x - a.x
///   2324d  mulsd  xmm4,xmm5      ; dx * a.y
///   23251  movapd xmm5,xmm7      ; v.x
///   23255  mulsd  xmm5,xmm3      ; v.x * dy
///   23259  mulsd  xmm3,xmm1      ; dy * a.x
///   2325d  addsd  xmm4,xmm5      ; s1 = dx*a.y + v.x*dy
///   23261  movapd xmm5,xmm1      ; a.x
///   23265  subsd  xmm5,xmm2      ; a.x - b.x
///   23269  mulsd  xmm5,xmm8      ; (a.x - b.x) * v.y
///   2326e  subsd  xmm5,xmm3      ; s2 = (a.x - b.x)*v.y - dy*a.x
///   23272  addsd  xmm4,xmm5      ; s  = s1 + s2
///   23276  comisd xmm15,xmm4     ; xmm15 is zero, set once by pxor at 0x231c2
///   2327b  ja     234d3          ; 0 > s, that is s < 0, rejects the zone
/// ```
///
/// The comparison itself is unchanged: still against an exact zero, so the caller
/// keeps `< 0.0`. No fused multiply-add appears, as expected on a baseline x86-64
/// target. `a.x - b.x` is written here rather than `-dx` because that is what the
/// instruction at `0x23265` does; the two are the same double, since IEEE
/// subtraction is exactly antisymmetric, but the site is quoted as it stands.
///
/// What this buys: on a point lying exactly on the edge two candidate cells share,
/// the source order returns a tiny negative value where the emitted order returns
/// an exact zero, so the eC's text rejects both cells and answers the null zone
/// where the engine names a cell. The order above reproduces the engine's answer
/// without any tolerance or epsilon.
fn point_line_side(px: f64, py: f64, ax: f64, ay: f64, bx: f64, by: f64) -> f64 {
    let dy = by - ay;
    let dx = bx - ax;
    let s1 = dx * ay + px * dy;
    let s2 = (ax - bx) * py - dy * ax;
    s1 + s2
}

/// The even-level tie-break, `pointLineSide(dx, dy, a, b) < 0` at RI7H.ec:1610,
/// 1616, 1626 and 1632, as the compiler emitted it.
///
/// faithful: at those four sites both endpoints are compile-time constants, so gcc
/// folded `A`, `B` and `C` and then re-associated the test itself, comparing the
/// running sum against `-C` instead of adding `C` and comparing against zero. The
/// four inlinings in `libdggal.so`, BuildID
/// `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, share one shape:
///
/// ```text
///   0x237a8  (RI7H.ec:1626)   0x239ea  (RI7H.ec:1610)
///   0x23d50  (RI7H.ec:1632)   0x23d80  (RI7H.ec:1616)
///
///   mulsd  xmm0,[A]           ; A * dx
///   movsd  xmm2,[-C]
///   mulsd  xmm5,[B]           ; B * dy
///   addsd  xmm0,xmm5          ; A*dx + B*dy
///   comisd xmm2,xmm0          ; -C against the sum
///   seta / ja / jbe           ; taken when -C > A*dx + B*dy
/// ```
///
/// The folded constants are in `.rodata` at `0x4a190` (`-1/3`), `0x4a198` (`1/6`),
/// `0x4a1a0` (`-1/4` to the last bit), `0x4a1a8` (`-1/6`), `0x4a1b0` (`1/3`),
/// `0x4a1b8` and `0x4a1c0` (`-1/12`), and the test below asserts that the values
/// this crate computes from the same endpoints match them bit for bit.
///
/// This rewriting cannot change an answer, and it is ported for the record rather
/// than for an effect: `fl(s + C) < 0` and `s < -C` can only disagree when `s + C`
/// is a non-zero value that underflows to zero, and where `s` is near `-C` the sum
/// is exact, so no such case exists at these magnitudes. The census confirms it:
/// not one even-level answer moved.
fn tie_break_side_is_negative(px: f64, py: f64, ax: f64, ay: f64, bx: f64, by: f64) -> bool {
    let dx = bx - ax;
    let dy = by - ay;
    let a = dy;
    let b = -dx;
    let c = ay * dx - dy * ax;
    a * px + b * py < -c
}

/// Rotates a 5x6 offset vector by sixty degrees. The oblique basis makes a sixty
/// degree turn a pure integer combination of the two components. Ports
/// `rotate5x6Offset` (`RI7H.ec:1933-1947`).
pub(super) fn rotate5x6_offset(dx: f64, dy: f64, clockwise: bool) -> (f64, f64) {
    if clockwise {
        (dx - dy, dx) // 60 degrees clockwise
    } else {
        (dy, dy - dx) // 60 degrees anticlockwise
    }
}

/// Walks an offset `(dx, dy)` from a 5x6 point and, when it crosses a rhombus
/// interruption (one of the icosahedron's cut edges), returns the crossing's in and
/// out points together with its hemisphere. Returns `None` when the walk stays
/// within a single rhombus. Ports `crosses5x6InterruptionV2` (ri5x6.ec:1611); it
/// drives the vertex tracing in [`add_non_polar_base_vertices`].
fn crosses5x6_interruption_v2(c_in_x: f64, c_in_y: f64, dx: f64, dy: f64) -> Option<Crossing> {
    let (mut cx0, mut cy0) = (c_in_x, c_in_y);
    if cx0 < 0.0 && cy0 < 1.0 + 1e-11 {
        cx0 += 5.0;
        cy0 += 5.0;
    }

    let (mut cdx, mut cdy) = (cx0, cy0);
    let north = cdx - cdy - 1e-11 > 0.0;
    if north {
        cdx -= 1e-11;
        cdy += 1e-11;
    } else {
        cdx += 1e-11;
        cdy -= 1e-11;
    }

    if cdx < 0.0 && cdy < 1.0 + 1e-11 {
        cdx += 5.0;
        cdy += 5.0;
        cx0 += 5.0;
        cy0 += 5.0;
    }
    if cdx > 5.0 && cdy > 5.0 - 1e-11 {
        cdx -= 5.0;
        cdy -= 5.0;
        cx0 -= 5.0;
        cy0 -= 5.0;
    }

    let icx = cdx.floor() as i64;
    let icy = cdy.floor() as i64;

    let mut px = if dx < 0.0 {
        pymax(icx as f64 - cx0, dx)
    } else {
        pymin(icx as f64 + 1.0 - cx0, dx)
    };
    let mut py = if dy < 0.0 {
        pymax(icy as f64 - cy0, dy)
    } else {
        pymin(icy as f64 + 1.0 - cy0, dy)
    };

    if dx != 0.0 && dy != 0.0 {
        let pkx = px / dx;
        let pky = py / dy;
        if pkx < pky {
            py = pkx * dy;
        } else if pky < pkx {
            px = pky * dx;
        }
    }

    cx0 += px;
    cy0 += py;

    if (dx - px).abs() < 1e-11 && (dy - py).abs() < 1e-11 {
        return None;
    }

    // faithful: `Sgn(0)` is 0, as the eC compiler expands it (see `sgn`), in
    // `nx = (int)floor(c.x + 1E-11 * Sgn(dx))` and its sibling (ri5x6.ec:1666-1667), and the
    // library tests for it here: in DGGAL's `-O2 -ffast-math` build, `libdggal.so`, BuildID
    // `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, `0x38946` compares `dx` with zero and
    // `0x3894a` takes no epsilon when they are equal, and `0x389be` and `0x389c2` do the same
    // for `dy`. This crate took it as 1 before, as py4dggs does; the answers are the same
    // either way here, since the walk steps along the edges of a hexagon and neither component
    // of any of its six edges is zero.
    let nx = (cx0 + 1e-11 * sgn(dx)).floor() as i64;
    let ny = (cy0 + 1e-11 * sgn(dy)).floor() as i64;

    if (nx != icx || ny != icy) && (nx > icx || (dx - px).abs() > 1e-11 || (dy - py).abs() > 1e-11)
    {
        let root = icx + icy;
        let mut i_src: Option<(f64, f64)> = None;
        let mut i_dst: (f64, f64) = (0.0, 0.0);
        let mut in_north = false;

        // faithful: the four `iy`/`ix` values below truncate in the eC,
        // `int iy = (int)(c.x - 1 + 1E-11)` and its three siblings at
        // ri5x6.ec:1678, 1687, 1700 and 1708, where py4dggs floors. py4dggs keeps
        // the difference on purpose, and it is inert on this crate's inputs: the
        // values are non-negative there, so flooring and truncating agree. That holds
        // of the points and offsets the two boundary tracers pass, not of the function
        // in general, and the crate's own inputs never reach the difference: over the
        // oracle suites and a sweep of every level, no input from the tracers reached a
        // value at which flooring and truncating differ. A synthetic sweep of 200,000
        // inputs, measured outside this test suite, found 3,105 answered differently
        // from the library, and all of them agree once these four values are truncated.
        // A new caller must not rely on it. The `icx`/`icy` and
        // `nx`/`ny` pairs above are not affected, since the eC floors there
        // (ri5x6.ec:1642-1643 and :1666-1667).
        if root & 1 == 0 {
            // North rhombus
            if ny == icy && nx == icx + 1 {
                let iy = (cx0 - 1.0 + 1e-11).floor();
                i_src = Some((cx0, cy0));
                i_dst = (iy + 2.0 - (cy0 - iy), cx0);
                in_north = true;
            } else if nx == icx && ny == icy - 1 {
                let ix = (cy0 + 1e-11).floor();
                i_src = Some((cx0, cy0));
                i_dst = (cy0, ix - (cx0 - ix));
                in_north = true;
            }
        } else {
            // South rhombus
            if nx == icx && ny == icy + 1 {
                let ix = (cy0 - 2.0 + 1e-11).floor();
                i_src = Some((cx0, cy0));
                i_dst = (cy0 - 1.0, ix + 3.0 - (cx0 - ix));
                in_north = false;
            } else if ny == icy && nx == icx - 1 {
                let iy = (cx0 + 1.0 + 1e-11).floor();
                i_src = Some((cx0, cy0));
                i_dst = (iy - 1.0 - (cy0 - iy), cx0 + 1.0);
                in_north = false;
            }
        }

        if let Some(i_src) = i_src {
            if i_dst.0 < 0.0 && i_dst.1 < 1.0 + 1e-11 {
                i_dst.0 += 5.0;
                i_dst.1 += 5.0;
            }
            return Some(Crossing {
                i_src,
                i_dst,
                in_north,
            });
        }
    }
    None
}

// --------------------------------------------------------------------------- //
// Hexagon vertex offset tables (RI7H.ec getVertices), in units of 1/(7*p)
// --------------------------------------------------------------------------- //

/// The rounded reciprocals the compiler substituted for the centroid property's
/// divisions by seven. They are the `.rodata` values at `0x49b48` and `0x49b40` in
/// `libdggal.so`, BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, and the tests
/// below pin them by bit pattern.
const ONE_SEVENTH: f64 = 1.0 / 7.0;
const TWO_SEVENTHS: f64 = 2.0 / 7.0;

pub(super) const EVEN_HEX_A: f64 = 7.0 / 3.0;
pub(super) const EVEN_HEX_B: f64 = 14.0 / 3.0;
pub(super) const ODD_HEX_A: f64 = 4.0 / 3.0;
pub(super) const ODD_HEX_B: f64 = 5.0 / 3.0;
pub(super) const ODD_HEX_C: f64 = 1.0 / 3.0;

const EVEN_HEX_VERTS: [(f64, f64); 6] = [
    (-EVEN_HEX_A, -EVEN_HEX_B),
    (-EVEN_HEX_B, -EVEN_HEX_A),
    (-EVEN_HEX_A, EVEN_HEX_A),
    (EVEN_HEX_A, EVEN_HEX_B),
    (EVEN_HEX_B, EVEN_HEX_A),
    (EVEN_HEX_A, -EVEN_HEX_A),
];
const ODD_HEX_VERTS: [(f64, f64); 6] = [
    (-ODD_HEX_A, -ODD_HEX_B),
    (-ODD_HEX_B, -ODD_HEX_C),
    (-ODD_HEX_C, ODD_HEX_A),
    (ODD_HEX_A, ODD_HEX_B),
    (ODD_HEX_B, ODD_HEX_C),
    (ODD_HEX_C, -ODD_HEX_A),
];

// --------------------------------------------------------------------------- //
// Zone centroid and containment (RI7H.ec centroid / containsPoint)
// --------------------------------------------------------------------------- //

/// The zone's centroid in 5x6 space: the `centroid` property of `I7HZone`
/// (`RI7H.ec:2895-2976`). An even level sits on
/// the rhombus lattice; an odd level is displaced by one of six sub-hexagon
/// offsets, re-expressed across an interruption through [`move5x6_vertex`]. An edge
/// hexagon and the south-pole corner rotate the sub-hexagon index first, matching
/// the children enumeration.
pub(super) fn zone_centroid(z: &Z) -> (f64, f64) {
    let l49r = z.l49r;
    let p = pow7(l49r);
    let oop = 1.0 / p as f64;
    let root = z.root;

    let (mut vx, mut vy);
    if root == 0xA {
        vx = 1.0;
        vy = 0.0;
        if z.sub_hex > 1 {
            // faithful: the two polar arms as the compiler emitted them. See the
            // citation on `oop_div7` below for the reciprocal; on top of it, the
            // compiler moved the bracket, so that the whole-number part is added to
            // the corner first and the fraction is applied to that sum, where the
            // eC forms `sh - 2 - 2*oop/7` and adds it to the corner in one go
            // (RI7H.ec:2918-2919). At `0x20429` in `libdggal.so`, BuildID
            // `e75d6ab18f8b460713ebcd21df6f07922fb909c3`:
            //
            // ```text
            //   20430  movsd  xmm3,[rip+0x29708]  ; 49b40 = 2/7
            //   20438  cvtsi2sd xmm1,ebp          ; (double)(sh - 2)
            //   2043c  mulsd  xmm3,xmm0           ; (2/7) * oop
            //   20440  mulsd  xmm0,[rip+0x29700]  ; 49b48 = 1/7, oop * (1/7)
            //   20448  addsd  xmm2,xmm1           ; 1 + (sh - 2)
            //   2044c  addsd  xmm0,xmm1           ; v.y = (1/7)*oop + (sh - 2)
            //   20450  subsd  xmm2,xmm3           ; v.x = (1 + (sh-2)) - (2/7)*oop
            // ```
            let k = (z.sub_hex - 2) as f64;
            vx = (vx + k) - TWO_SEVENTHS * oop;
            // The eC's `v.y += ...` starts from zero, which the compiler dropped;
            // adding a positive zero is exact, so the term stands on its own here.
            vy = ONE_SEVENTH * oop + k;
        }
    } else if root == 0xB {
        vx = 4.0;
        vy = 6.0;
        if z.sub_hex > 1 {
            // faithful: the south arm, the same rewriting, emitted as a pair of
            // packed operations at `0x20339`:
            //
            // ```text
            //   20340  movapd xmm0,[rip+0x29738]  ; 49a80 = (4.0, 6.0), the corner
            //   2034c  mulpd  xmm2,[rip+0x297ec]  ; 49b40 = (2/7, 1/7), times oop
            //   20354  unpcklpd xmm1,xmm1         ; ((sh-2), (sh-2))
            //   20358  subpd  xmm0,xmm1           ; corner - (sh - 2)
            //   20360  addpd  xmm1,xmm0           ; + the fractions, for v.x
            //   20364  subpd  xmm0,xmm2           ; - the fractions, for v.y
            //   20368  movsd  xmm0,xmm1           ; v.x from the first, v.y the second
            // ```
            let k = (z.sub_hex - 2) as f64;
            vx = (vx - k) + TWO_SEVENTHS * oop;
            vy = (vy - k) - ONE_SEVENTH * oop;
        }
    } else {
        let cx = root >> 1;
        let cy = cx + (root & 1);
        vx = cx as f64 + z.col as f64 * oop;
        vy = cy as f64 + z.row as f64 * oop;
    }

    if z.sub_hex != 0 && root < 10 {
        let mut sh = z.sub_hex;
        let south = root & 1 != 0;
        if z.row == 0 && z.col == 0 && south && sh >= 4 {
            sh += 1;
        } else if sh >= 2 {
            let parent_z = Z::new(z.l49r, z.root, z.row, z.col, 0);
            if zone_is_edge_hex(&parent_z) {
                // `sh >= 2` on this arm, so the left operand of the remainder is at
                // least one and a plain `%` matches Python's own result.
                sh = (sh + if south { -1 } else { 3 }) % 6 + 2;
            }
        }

        // faithful: the eC writes `oop /= 7` (RI7H.ec:2947) and this crate divided.
        // The compiler turned the division by the compile-time seven into a
        // multiplication by its rounded reciprocal, which `-freciprocal-math`, part
        // of `-ffast-math`, permits and which is not exact: `1/7` has no finite
        // binary expansion, and `oop * fl(1/7)` differs from `oop / 7` in the last
        // place at `levelI49R` 2, 3, 4 and 7, that is at resolutions 4 to 9 and 14
        // and 15. The centroid property (RI7H.ec:2895) carries no
        // `__attribute__ ((optimize("-fno-unsafe-math-optimizations")))`, and DGGAL
        // v0.0.6 is built with `-O2 -ffast-math` (`Makefile.dggal:133`). At
        // `0x20310` in `libdggal.so`, BuildID
        // `e75d6ab18f8b460713ebcd21df6f07922fb909c3`:
        //
        // ```text
        //   20310  mulsd  xmm0,[rip+0x29830]  ; 49b48 = 0.14285714285714285
        // ```
        //
        // The six sub-hexagon offsets below keep the eC's own `k * oop`, which the
        // compiler folded to a constant times the reciprocal-scaled `oop` at
        // `0x20560`, `0x20588`, `0x205b8`, `0x205e0`, `0x20610` and `0x20640`.
        let oop_div7 = oop * ONE_SEVENTH;
        let mut ddx = 0.0;
        let mut ddy = 0.0;
        // The Python writes `-1 * oop_div7` and `1 * oop_div7`; the unit factors are
        // dropped here because negation and multiplication by one are exact in
        // IEEE 754, so the values are bit-identical.
        match sh {
            2 => {
                ddx = -oop_div7;
                ddy = -3.0 * oop_div7;
            }
            3 => {
                ddx = -3.0 * oop_div7;
                ddy = -2.0 * oop_div7;
            }
            4 => {
                ddx = -2.0 * oop_div7;
                ddy = oop_div7;
            }
            5 => {
                ddx = oop_div7;
                ddy = 3.0 * oop_div7;
            }
            6 => {
                ddx = 3.0 * oop_div7;
                ddy = 2.0 * oop_div7;
            }
            7 => {
                ddx = 2.0 * oop_div7;
                ddy = -oop_div7;
            }
            _ => {}
        }
        if sh != 1 {
            let moved = move5x6_vertex(vx, vy, ddx, ddy);
            vx = moved.0;
            vy = moved.1;
        }
    }

    if vy.abs() < 1e-6 {
        vy = 0.0;
    } else if (vx - 5.0).abs() < 1e-6 {
        vx = 5.0;
    }
    if vx > 5.0 - 1e-6 || vy > 6.0 + 1e-6 {
        vx -= 5.0;
        vy -= 5.0;
    }
    if vx < -1e-6 {
        vx += 5.0;
        vy += 5.0;
    }

    (vx, vy)
}

/// The 5x6 outline of a polar pentagon, as the eC builds it for `containsPoint`
/// (`getBaseRefinedVerticesNoAlloc`, the `root > 9 && subHex == 1` arm at
/// RI7H.ec:2428-2480, reached through RI7H.ec:956).
///
/// A polar pentagon is not a small polygon around a centroid: in 5x6 space it is a
/// wedge that spans all five rhombi meeting at the pole, so neither the vertex
/// offset table nor the non-polar boundary tracing describes it. The eC walks its
/// five sides with `addIntermediatePointsNoAlloc`, which at `nDivisions == 1` and
/// `crs84 == false` emits exactly the side's start point followed by the two
/// interruption points it was handed (ri5x6.ec:758-765), and then closes the ring
/// with the far side, the pole itself and the first interruption point. That comes
/// to seventeen vertices.
///
/// The eC also computes `l1`, `l2` and `t` with two square roots before that loop
/// (ri5x6.ec:753-756). At `nDivisions == 1` the loop they feed never runs, so they
/// cannot affect the result and are not reproduced here; that also keeps this file
/// free of a square root, which the crate confines to `math.rs`.
fn polar_pentagon_outline(root: i64, oonp: f64) -> Vec<(f64, f64)> {
    let a_off = ODD_HEX_A;
    let b_off = ODD_HEX_B;
    let c_off = ODD_HEX_C;
    let r = 1.0 / 5.0;
    let mut out: Vec<(f64, f64)> = Vec::with_capacity(17);

    if root == 0xA {
        // North pole (RI7H.ec:2431-2454), in the order the compiler emitted it at
        // `0x21260` in `libdggal.so`, BuildID
        // `e75d6ab18f8b460713ebcd21df6f07922fb909c3`.
        //
        // faithful: three rewritings, of the kinds already met above.
        //
        // First, the whole number of each fan step moved inside the bracket, so
        // that `b.x + k` is `(1 + k) - oonp*C` and `ab.x + k` is
        // `(1 + k) + fl(-1.4)*oonp`, where the eC shifts the first corner. The
        // packed corners `(5,4)`, `(4,3)`, `(3,2)` and `(2,1)` are loaded at
        // `0x2126d`, `0x212b8`, `0x212e1` and `0x21302`, and the subtraction of the
        // offset follows each of them:
        //
        // ```text
        //   21260  mulpd  xmm1,[rip+0x28938]  ; 49ba0 = (1/3, 4/3) times oonp
        //   2129e  subpd  xmm4,xmm1           ; (5,4) - (C', A')
        //   212a6  addpd  xmm3,xmm7           ; (C', A') + (5,4)
        //   212b4  movsd  xmm3,xmm4           ; x from the first, y from the second
        // ```
        //
        // Second, `ab.x = a.x + (b.x - a.x)*r` was folded whole: every `oonp` term
        // collapses into one coefficient, `-B + (B - C)*r`, which rounds to exactly
        // `-1.4` and sits at `0x49bb0`. `ab.y` keeps the eC's own value, since
        // `b.y - a.y` is `A' - (-C')`, which the compiler wrote as `C' + A'`, an
        // exact identity.
        //
        // Third, `c` is never formed: `c.x - b.x` became `C' + A'` at `0x2134d` and
        // `(c.y - b.y)*r` became `fl((B - A)*r) * oonp`, the constant `0x4a160`
        // loaded at `0x21351`. `d.x` and `d.y` are then assembled at `0x213a7`, and
        // `d + k` is a plain packed addition of `(k, k)` at `0x213b5`, since `d`
        // carries no integer corner to move.
        let c_p = oonp * c_off;
        let a_p = oonp * a_off;
        let b = (1.0 - c_p, a_p);
        let ab_off_x = (-b_off + (b_off - c_off) * r) * oonp;
        let ab_off_y = (c_p + a_p) * r - c_p;
        let d = (b.0 + (c_p + a_p) * r, (b_off - a_off) * r * oonp + a_p);
        for k in 0..4 {
            let k = k as f64;
            out.push(((1.0 + k) - c_p, k + a_p));
            out.push((k + d.0, k + d.1));
            let k1 = k + 1.0;
            out.push(((1.0 + k1) + ab_off_x, k1 + ab_off_y));
        }
        out.push((5.0 - c_p, 4.0 + a_p));
        // This extends to the right border of the last triangle.
        out.push((4.0 + d.0, 4.0 + d.1));
        // These two are the north pole itself.
        out.push((5.0, 4.0));
        out.push((1.0, 0.0));
        out.push((1.0 + ab_off_x, ab_off_y));
    } else {
        // South pole (RI7H.ec:2456-2479), at `0x21608` in the same library.
        //
        // faithful: the compiler treated this arm differently from the north one,
        // and the difference is reproduced rather than tidied. Here `ab` is
        // computed as the eC writes it, `a + (b - a)*r`, with `b - a` a packed
        // subtraction at `0x216db`, and both `ab` and `d` are then shifted by a
        // plain packed addition of `(-k, -k)` at `0x21727` and `0x21761`. Only the
        // `b - k` fan moved, and it moved the same way as the north one:
        //
        // ```text
        //   216bb  addsd  xmm12,xmm3   ; 3 + C', that is b.x - 1
        //   2174e  subsd  xmm11,xmm4   ; 5 - A', that is b.y - 1
        //   2180e  addsd  xmm10,xmm3   ; 2 + C', that is b.x - 2
        //   218b2  subsd  xmm7,xmm4    ; 3 - A', that is b.y - 3
        //   21972  subsd  xmm2,xmm4    ; 2 - A', that is b.y - 4
        // ```
        //
        // `c` is not formed here either: `c.x - b.x` is `(4 - A') - b.x` at
        // `0x216f7`, which is the eC's own expression with the corner folded, and
        // `c.y - b.y` became `A' - B'` at `0x2170a`, a re-association. `d` is then
        // `(c - b)*r + b` at `0x2174a`, the eC's sum with its terms exchanged.
        let b_p = oonp * b_off;
        let c_p = oonp * c_off;
        let a_p = oonp * a_off;
        let a = (4.0 + b_p, 6.0 + c_p);
        let b = (4.0 + c_p, 6.0 - a_p);
        let ab = (a.0 + (b.0 - a.0) * r, a.1 + (b.1 - a.1) * r);
        let d = (((4.0 - a_p) - b.0) * r + b.0, (a_p - b_p) * r + b.1);
        for k in 0..4 {
            let k = k as f64;
            out.push(((4.0 - k) + c_p, (6.0 - k) - a_p));
            out.push((d.0 - k, d.1 - k));
            let k1 = k + 1.0;
            out.push((ab.0 - k1, ab.1 - k1));
        }
        out.push((0.0 + c_p, 2.0 - a_p));
        // This extends to the left wrapping point.
        out.push((d.0 - 4.0, d.1 - 4.0));
        // These two are the south pole itself.
        out.push((0.0, 2.0));
        out.push((4.0, 6.0));
        out.push(ab);
    }
    out
}

/// The polygon the eC hands to `containsPoint`: `getBaseRefinedVerticesNoAlloc`
/// with `crs84` false and `nDivisions` one (called at RI7H.ec:956, defined at
/// :2333).
///
/// Only an odd-level zone reaches this, because [`contains_point`]'s single caller
/// is `from_centroid`'s odd-level branch, which tests the primary children of an
/// even-level candidate parent and those are all odd. The eC's even-level arms are
/// therefore not ported, and the assertion below records that expectation rather
/// than leaving it to be rediscovered.
fn base_refined_vertices(z: &Z, oonp: f64) -> Vec<(f64, f64)> {
    debug_assert!(
        z.sub_hex != 0,
        "base_refined_vertices: only an odd-level zone reaches containsPoint"
    );
    if z.root > 9 && z.sub_hex == 1 {
        // Polar pentagons (RI7H.ec:2428).
        return polar_pentagon_outline(z.root, oonp);
    }
    // Odd level, every other zone (RI7H.ec:2482-2494).
    let v = scale_offsets(&ODD_HEX_VERTS, oonp);
    add_non_polar_vertices_refined(z, centroid5x6(z), &v, Emission::NoAlloc)
}

/// What the walker of a zone's boundary, [`add_non_polar_vertices_refined`], adds for each side
/// it walks. The eC has the walker twice, and the two copies differ in this alone:
/// `addNonPolarVerticesRefinedNoAlloc` (`RI7H.ec:2094-2165`) hands each side to
/// `addIntermediatePointsNoAlloc` (`ri5x6.ec:728-791`), and `addNonPolarVerticesRefined`
/// (`RI7H.ec:2021-2092`) to `addIntermediatePoints` (`ri5x6.ec:651-726`).
#[derive(Clone, Copy)]
enum Emission {
    /// The copy `containsPoint` runs, at `nDivisions` one and `crs84` false, where
    /// `addIntermediatePointsNoAlloc` writes the side's first point and, where the side crosses
    /// an interruption, the two points of the crossing (`ri5x6.ec:758-765`). It has no pole
    /// rule on this path: its flag is raised under `crs84` alone (`ri5x6.ec:732`).
    NoAlloc,
    /// The copy the refined boundary runs: each side divided by
    /// [`add_intermediate_points`], which has the pole rule whatever `crs84` is.
    Refined { crs84: bool, n_divisions: i32 },
}

/// Traces a non-polar cell's boundary side by side, for `containsPoint` and for the refined
/// boundary: it ports both `addNonPolarVerticesRefinedNoAlloc` (RI7H.ec:2094), at
/// `nDivisions` one, and `addNonPolarVerticesRefined` (`RI7H.ec:2021-2092`), which the eC
/// writes with the same text and the library performs in the same order, as the comments below
/// show for each copy. What each adds for a side is the `emission`.
///
/// This is deliberately a second tracer beside [`add_non_polar_base_vertices`],
/// because the eC keeps two: `getVertices` emits one point per side, whereas this
/// one also emits the two interruption points `i1` and `i2` whenever a side crosses
/// a rhombus seam, since `addIntermediatePointsNoAlloc` at `nDivisions` one and
/// `crs84` false writes `p`, `i1` and `i2` in that order (ri5x6.ec:758-765). A cell
/// lying on a seam therefore carries more vertices here than in its published
/// geometry, and the containment test needs them: without the pair the polygon
/// short-circuits across the interruption.
fn add_non_polar_vertices_refined(
    z: &Z,
    c: (f64, f64),
    v: &[(f64, f64); 6],
    emission: Emission,
) -> Vec<(f64, f64)> {
    let n_points = zone_npoints(z);
    let (cx, cy) = c;
    // faithful: as in [`add_non_polar_base_vertices`], and for the same reason.
    // `addNonPolarVerticesRefinedNoAlloc` (RI7H.ec:2093) is unguarded too, and gcc
    // compiled it separately, at `0x1f520` to `0x1fa22` in `libdggal.so`, BuildID
    // `e75d6ab18f8b460713ebcd21df6f07922fb909c3`. Its head is rewritten exactly as
    // the other walker's is:
    //
    // ```text
    //   1f566  addsd  xmm7,xmm4        ; xmm7 = 1e-11 + c.x, hoisted
    //   1f595  addsd  xmm3,xmm2        ; t.y  = c.y + v[i].y
    //   1f599  addsd  xmm1,xmm7        ; v[i].x + (c.x + 1e-11)
    //   1f5f7  addsd  xmm6,xmm2        ; point.y = c.y + v[start].y
    //   1f5ff  addsd  xmm4,xmm5        ; point.x = c.x + v[start].x
    //   1f60b  subsd  xmm5,[rdi]       ; dir.x = v[start].x - v[prev].x
    //   1f621  subsd  xmm2,[rdi+0x8]   ; dir.y = v[start].y - v[prev].y
    // ```
    //
    // This is the walker whose result `containsPoint` actually tests, so it is the
    // one that decides the quantiser's answers on a boundary.
    //
    // faithful: `addNonPolarVerticesRefined` (RI7H.ec:2021-2092), the allocating copy, is
    // unguarded as well and compiled on its own, at `0x28cc0` to `0x2920f` in the same library,
    // with the same head (RI7H.ec:2030-2041):
    //
    // ```text
    //   28d0a  addsd  xmm7,xmm4        ; xmm7 = 1e-11 + c.x, hoisted
    //   28d35  addsd  xmm3,xmm2        ; t.y  = c.y + v[i].y
    //   28d39  addsd  xmm1,xmm7        ; v[i].x + (c.x + 1e-11)
    //   28d97  addsd  xmm6,xmm2        ; point.y = c.y + v[start].y
    //   28d9f  addsd  xmm4,xmm5        ; point.x = c.x + v[start].x
    //   28dab  subsd  xmm5,[rsi]       ; dir.x = v[start].x - v[prev].x
    //   28dc1  subsd  xmm2,[rsi+0x8]   ; dir.y = v[start].y - v[prev].y
    // ```
    let cx_eps = cx + 1e-11;
    let mut start: usize = 0;
    for (i, vi) in v.iter().enumerate() {
        let ty = cy + vi.1;
        let itx = (vi.0 + cx_eps).floor();
        if !(ty - itx > 2.0 || ty < itx) {
            start = i;
            break;
        }
    }

    let mut point = (cx + v[start].0, cy + v[start].1);
    // `start` is in 0..6, so the left operand of the remainder is positive.
    let prev = (start + 5) % 6;
    let mut direction = (v[start].0 - v[prev].0, v[start].1 - v[prev].1);

    // Sized once. For `containsPoint` the eC's buffer is fixed (`Pointd v5x6[24]`,
    // RI7H.ec:955): three points a side at most, the start and the two interruption
    // points. For the refined boundary, see [`refined_capacity`].
    let mut out: Vec<(f64, f64)> = Vec::with_capacity(match emission {
        Emission::NoAlloc => 3 * n_points as usize,
        Emission::Refined { n_divisions, .. } => refined_capacity(n_divisions),
    });
    for _ in 0..n_points {
        let ndir = rotate5x6_offset(direction.0, direction.1, false); // 60 anticlockwise
        let mut n = (point.0 + ndir.0, point.1 + ndir.1);

        let mut p = point;
        if p.0 > 5.0 && p.1 > 5.0 {
            p.0 -= 5.0;
            p.1 -= 5.0;
        }
        if p.0 < 0.0 || p.1 < 0.0 {
            p.0 += 5.0;
            p.1 += 5.0;
        }

        if let Some(cross) = crosses5x6_interruption_v2(p.0, p.1, ndir.0, ndir.1) {
            let mut i1 = cross.i_src;
            let mut i2 = cross.i_dst;
            let north = cross.in_north;
            if point.0 - p.0 > 4.0 {
                i1.0 += 5.0;
                i1.1 += 5.0;
                i2.0 += 5.0;
                i2.1 += 5.0;
            }
            if p.0 - point.0 > 4.0 {
                i1.0 -= 5.0;
                i1.1 -= 5.0;
                i2.0 -= 5.0;
                i2.1 -= 5.0;
            }
            if i2.1 - i1.1 > 4.0 {
                i2.0 -= 5.0;
                i2.1 -= 5.0;
            }
            if i1.1 - i2.1 > 4.0 {
                i2.0 += 5.0;
                i2.1 += 5.0;
            }

            let crossing_left = if north { i2.0 < i1.0 } else { i2.0 > i1.0 };
            // faithful: the remainder of the step and the second rotation, as this
            // walker was compiled. The rewriting is of the same kind as the other
            // walker's and is **not the same rewriting**: gcc formed `ddx` and
            // `ddy` here by adding the rotated direction to the offset of `point`
            // from `i1`, where in `addNonPolarBaseVertices` it subtracted `i1` from
            // the uncrossed next point. The two brackets differ in the last bits,
            // so each is written as its own function was compiled. At `0x1f520`:
            //
            // ```text
            //   1f77b  subsd  xmm5,xmm6    ; point.x - i1.x
            //   1f78d  subsd  xmm3,xmm4    ; point.y - i1.y
            //   1f79e  addpd  xmm4,xmm7    ; ddx = (point.x - i1.x) + ndir.x
            //   1f7a2  addsd  xmm9,xmm1    ; ddy = (point.y - i1.y) + ndir.y
            //   -- crossingLeft false, the clockwise arm --
            //   1f7bb  addsd  xmm4,xmm0    ; n.y = ddx + i2.y
            //   1f7bf  subsd  xmm3,xmm9    ; d.x = ddx - ddy
            //   1f7c4  addsd  xmm3,xmm8    ; n.x = d.x + i2.x
            //   (dir is written nowhere on this arm)
            //   -- crossingLeft true, the anticlockwise arm --
            //   1f951  addpd  xmm2,xmm3    ; d.x = ndir.y + (point.y - i1.y)
            //   1f955  subpd  xmm3,xmm4    ; d.y = ddy - ddx
            //   1f963  movsd  [rsp+0x8],xmm1 ; dir.x := ndir.y
            //   1f969  xorpd  xmm2,...     ; dir.y := -dir.x
            //   1f971  addpd  xmm3,xmm0    ; n = i2 + d
            // ```
            //
            // Unlike the other walker, the two `n = i2 + d` sums keep the eC's own
            // bracketing here, so only `ddx`, `ddy` and the second rotation move.
            //
            // faithful: the allocating copy (RI7H.ec:2077-2079) is compiled to the same forms,
            // at `0x28cc0`:
            //
            // ```text
            //   28f13  subsd  xmm5,xmm6    ; point.x - i1.x
            //   28f24  subsd  xmm3,xmm4    ; point.y - i1.y
            //   28f35  addpd  xmm4,xmm7    ; ddx = (point.x - i1.x) + ndir.x, in the high half
            //   28f39  addsd  xmm9,xmm1    ; ddy = (point.y - i1.y) + ndir.y
            //   -- crossingLeft false, the clockwise arm --
            //   28f52  addsd  xmm0,xmm4    ; n.y = i2.y + ddx
            //   28f56  subsd  xmm3,xmm9    ; d.x = ddx - ddy
            //   28f5b  addsd  xmm3,xmm8    ; n.x = d.x + i2.x
            //   (dir is written nowhere on this arm: `[rsp]` and `[rsp+0x8]` keep the old pair)
            //   -- crossingLeft true, the anticlockwise arm --
            //   29131  addpd  xmm2,xmm3    ; d.x = ndir.y + (point.y - i1.y), in the low half
            //   29135  subpd  xmm3,xmm4    ; d.y = ddy - ddx, in the high half
            //   29142  movsd  [rsp],xmm1   ; dir.x := ndir.y
            //   29147  xorpd  xmm2,...     ; dir.y := -dir.x
            //   2914f  addpd  xmm3,xmm0    ; n = i2 + d
            // ```
            let ddx = (point.0 - i1.0) + ndir.0;
            let ddy = (point.1 - i1.1) + ndir.1;
            if crossing_left {
                n = (i2.0 + ddy, i2.1 + (ddy - ddx));
                direction = (ndir.1, -direction.0);
            } else {
                n = ((ddx - ddy) + i2.0, ddx + i2.1);
                // `direction` is deliberately left as it was: see above.
            }
            match emission {
                // `addIntermediatePointsNoAlloc(.., point, n, 1, i1, i2, false)`.
                Emission::NoAlloc => {
                    out.push(point);
                    out.push(i1);
                    out.push(i2);
                }
                // `addIntermediatePoints(vertices, point, n, nDivisions, i1, i2, crs84)`
                // (RI7H.ec:2081), called at `0x28f8e`: `nDivisions` goes as it is, a zero
                // too, which the callee counts as one.
                Emission::Refined { crs84, n_divisions } => {
                    add_intermediate_points(&mut out, point, n, n_divisions, Some((i1, i2)), crs84);
                }
            }
        } else {
            match emission {
                // `addIntermediatePointsNoAlloc(.., point, n, 1, null, null, false)`,
                // whose null arm writes the single point `p` at `nDivisions` one.
                Emission::NoAlloc => out.push(point),
                // `if(!nDivisions) vertices.Add(point)` (RI7H.ec:2085-2086): the test of
                // `r13d` at `0x29082`, and the `Add` called at `0x290bd`.
                Emission::Refined { n_divisions: 0, .. } => out.push(point),
                // `addIntermediatePoints(vertices, point, n, nDivisions, null, null, crs84)`
                // (RI7H.ec:2088), called at `0x2919a`.
                Emission::Refined { crs84, n_divisions } => {
                    add_intermediate_points(&mut out, point, n, n_divisions, None, crs84);
                }
            }
            direction = ndir;
        }
        point = n;
    }
    out
}

/// The points from which the five sides of a polar pentagon's refined boundary are built, at
/// each of the five turns about the pole: side `k` runs from `b[k]` to `b[k + 1]`, leaves its
/// rhombus at `d[k]` and enters the next at `ab[k + 1]`, and the fifth side comes back to
/// `b[0]` through `d[4]` and `ab[0]`. They are the eC's `b`, `d` and `ab` shifted by `k` along
/// the diagonal of the 5x6 layout (`RI7H.ec:2197-2202`), towards the south-east about the north
/// pole and towards the north-west about the south pole.
struct PolarSides {
    b: [(f64, f64); 5],
    d: [(f64, f64); 5],
    ab: [(f64, f64); 5],
    /// The two images of the pole between which the 5x6 CRS closes the ring.
    pole: [(f64, f64); 2],
}

/// The sides of the polar pentagon about `root` (`0xA` the north pole, `0xB` the south), at an
/// even level or at an odd one: the four polar arms of the allocating
/// `I7HZone::getBaseRefinedVertices` (`RI7H.ec:2187-2212`, `:2213-2239`, `:2264-2288` and
/// `:2289-2313`), in the order the library performs them. `None` for a root that is neither
/// pole, where the eC's odd arm adds nothing (`RI7H.ec:2261-2314`).
///
/// The function is not guarded, and gcc re-associated each arm in its own way. Every form below
/// was read from DGGAL v0.0.6 as built with `-O2 -ffast-math` (`libdggal.so`, BuildID
/// `e75d6ab18f8b460713ebcd21df6f07922fb909c3`), where `getBaseRefinedVertices` lies at
/// `0x29220` to `0x2a857`. The odd arms are compiled as those of the `NoAlloc` copy, which
/// [`polar_pentagon_outline`] ports for `containsPoint`, save at one site (ON-4); the even arms
/// have no counterpart there. The two copies are kept apart, each beside the addresses it was
/// read from.
///
/// Each arm shifts its points with `k` from zero. At zero the sums are exact (no coordinate
/// here is a negative zero), so the first of each five is the unshifted point, as in the
/// library, which drops the eC's `+ 0` and `- 0`.
fn polar_pentagon_sides(root: i64, odd: bool, oonp: f64) -> Option<PolarSides> {
    let turns = |shifted: &dyn Fn(f64) -> (f64, f64)| {
        [
            shifted(0.0),
            shifted(1.0),
            shifted(2.0),
            shifted(3.0),
            shifted(4.0),
        ]
    };
    let r = 1.0 / 5.0;
    match (root, odd) {
        // The north pole at an even level (RI7H.ec:2187-2212), compiled at `0x29990`.
        (0xA, false) => {
            // `oonp * A` and `oonp * B`, taken before the roots are told apart: `0x2930f`
            // (`.rodata` `0x4a110`) and `0x2930b` (`0x4a118`).
            let a_p = oonp * EVEN_HEX_A;
            let b_p = oonp * EVEN_HEX_B;
            // faithful: site EN-2. `ab.x = (a.x + b.x) / 2` (RI7H.ec:2191), with
            // `a.x = 1 - oonp * B` and `b.x = 1 - oonp * A`, is `(((1 - A') - B') + 1) * 0.5`:
            // `subsd` at `0x299aa` and `0x29a5b`, `addsd` at `0x29a7d`, `mulsd` by the half of
            // `.rodata` `0x4a028` at `0x29a89`. It is the source's double at each of the ten
            // even levels; it is ported as compiled all the same.
            let ab_x = (((1.0 - a_p) - b_p) + 1.0) * 0.5;
            // faithful: site EN-5. `d.y`, the second member of the rotated offset
            // `(b.y - ab.y) - (b.x - ab.x)` added to `b.y` (RI7H.ec:2194-2195), is
            // `(3 A' - 1) + ab.x`: `mulsd` by three at `0x29a75`, `subsd` at `0x29a9e`, `addsd`
            // at `0x29aaf`. The eC's order gives another double at every even level.
            let d_y = (3.0 * a_p - 1.0) + ab_x;
            Some(PolarSides {
                // faithful: site EN-1. The eC shifts the corner `b`, as `b.x + k` with
                // `b.x = 1 - A'` (RI7H.ec:2197-2200); the library adds `A'` to, and subtracts
                // it from, the packed whole numbers `(5, 4)`, `(4, 3)`, `(3, 2)` and `(2, 1)`
                // of `.rodata` `0x49b50` to `0x49b80`, so that `b.x + k` is `(1 + k) - A'`:
                // `subpd` at `0x299f3`, `0x29a10`, `0x29a31` and `0x29a4e`, with `A' + k` from
                // the `addpd` at `0x299d9`, `0x29a0c`, `0x29a28` and `0x29a52`. The eC's order
                // gives another double at levels 0, 4, 10 and 16.
                b: turns(&|k| ((1.0 + k) - a_p, a_p + k)),
                // faithful: site EN-4. `d.x`, which is `(b.y - ab.y) + b.x` (RI7H.ec:2194-2195), is folded to
                // one,
                // the word of `.rodata` `0x49a70` stored at `0x299ea`, and `d.x + k` to the
                // whole numbers stored at `0x29b5c`, `0x29bca`, `0x29c57` and `0x29d25` (or
                // `0x2a2f3`); `d.y + k` is the `addsd` at `0x29af1`, `0x29bea`, `0x29c77` and
                // `0x29d20` (or `0x2a2ee`).
                d: turns(&|k| (1.0 + k, d_y + k)),
                // faithful: site EN-3. `ab.y = (a.y + b.y) / 2` (RI7H.ec:2191) is folded to zero, stored at
                // `0x299af`, so that `ab.y + k` is `k`; `ab.x + k` is the `addsd` at `0x29ab3`
                // and the `addpd` at `0x29b4a`, `0x29bc1` and `0x29c4e`.
                ab: turns(&|k| (ab_x + k, k)),
                pole: [(5.0, 4.0), (1.0, 0.0)],
            })
        }
        // The south pole at an even level (RI7H.ec:2213-2239), compiled at `0x29db0`.
        (0xB, false) => {
            let a_p = oonp * EVEN_HEX_A;
            let b_p = oonp * EVEN_HEX_B;
            // `b.x = 4 + oonp * A` (`0x29de4`). `ab.x` is the eC's own `(a.x + b.x) / 2`
            // (RI7H.ec:2217), `a.x` being `B' + 4` (`0x29dc7`): the sum at `0x29e10` and the
            // product by the half at `0x29e14`, which is the quotient by two exactly.
            let b_x = 4.0 + a_p;
            let ab_x = ((b_p + 4.0) + b_x) * 0.5;
            // faithful: site ES-5. `d.y`, the rotated offset `(b.y - ab.y) - (b.x - ab.x)`
            // added to `b.y` (RI7H.ec:2220-2221), is `((6 - A') - A') + (ab.x - b.x)`: `subsd`
            // at `0x29e0c`, `0x29e24` and `0x29e46`, `addsd` at `0x29e5b`. It is the source's
            // double at each of the ten even levels; it is ported as compiled all the same.
            let d_y = ((6.0 - a_p) - a_p) + (ab_x - b_x);
            Some(PolarSides {
                // faithful: site ES-1. The eC shifts `b = (4 + A', 6 - A')` by `k`
                // (RI7H.ec:2223-2226); the library forms `(4 - k) + A'` and `(6 - k) - A'`
                // from the whole numbers: `addsd` at `0x29e5f`, `0x29f42` and `0x29ff7`,
                // `subsd` at `0x29e81`, `0x29f03`, `0x29fb2` and `0x2a098`, and at the fourth
                // turn `A'` itself (`0x2a09c`). The eC's order gives another double at every
                // even level.
                b: turns(&|k| ((4.0 - k) + a_p, (6.0 - k) - a_p)),
                // faithful: site ES-4. `d.x`, which is `(b.y - ab.y) + b.x` (RI7H.ec:2220-2221), is folded to
                // four,
                // the word of `.rodata` `0x49a80` stored at `0x29de8`, and `d.x - k` to the
                // whole numbers stored at `0x29f18`, `0x29fc5`, `0x2a054` and `0x2a144` (or
                // `0x2a251`); `d.y - k` is the `subsd` at `0x29f2d`, `0x29fe2`, `0x2a07f` and
                // `0x2a13c` (or `0x2a26a`).
                d: turns(&|k| (4.0 - k, d_y - k)),
                // faithful: site ES-3. `ab.y = (a.y + b.y) / 2` (RI7H.ec:2217) is folded to six, the word of
                // `.rodata` `0x49a88` stored at `0x29dfc`, and `ab.y - 1` to the five stored at
                // `0x29e78`; `ab - k` is the `subsd` at `0x29e42` and the `addpd` of `(-k, -k)`
                // at `0x29f0f`, `0x29fb6` and `0x2a061`.
                ab: turns(&|k| (ab_x - k, 6.0 - k)),
                pole: [(0.0, 2.0), (4.0, 6.0)],
            })
        }
        // The north pole at an odd level (RI7H.ec:2264-2288), compiled at `0x2a390`.
        (0xA, true) => {
            // `(C', A')`, the packed product of `oonp` and `.rodata` `0x49ba0` at `0x2a390`.
            let c_p = oonp * ODD_HEX_C;
            let a_p = oonp * ODD_HEX_A;
            // faithful: site ON-2. `ab = a + (b - a) * r` (RI7H.ec:2268) is folded whole: in
            // `ab.x` every term in `oonp` collapses into the one coefficient
            // `-B + (B - C) * r`, which rounds to `-1.4` (`.rodata` `0x49bb0`, the `mulpd` at
            // `0x2a48b`), and `ab.y` is `(C' + A') * r - C'` (`addpd` at `0x2a476`, the same
            // `mulpd`, `subpd` at `0x2a4a0`). The eC's order gives another double at levels
            // 1, 9 and 15. The offset is kept apart from the corner, for the first turn takes
            // it so (ON-4).
            let ab_off = (
                (-ODD_HEX_B + (ODD_HEX_B - ODD_HEX_C) * r) * oonp,
                (c_p + a_p) * r - c_p,
            );
            // `ab` itself, `(1 + ab_off.x, ab_off.y)`: the `addpd` at `0x2a4a5`, stored at
            // `0x2a4b5`.
            let ab = (1.0 + ab_off.0, ab_off.1);
            // faithful: site ON-3. `c` is never formed: in `d = b + (c - b) * r`
            // (RI7H.ec:2269-2270), `c.x - b.x` is `C' + A'`, and `(c.y - b.y) * r` is the
            // folded `fl((B - A) * r)` of `.rodata` `0x4a160` times `oonp` (`0x2a4d0`); the
            // `mulpd` at `0x2a493` and the `addpd` at `0x2a4f2` assemble both members. It is
            // the source's double at each of the ten odd levels; it is ported as compiled all
            // the same.
            let d = (
                (1.0 - c_p) + (c_p + a_p) * r,
                (ODD_HEX_B - ODD_HEX_A) * r * oonp + a_p,
            );
            Some(PolarSides {
                // faithful: site ON-1. As at the even level (EN-1), `b.x + k` (RI7H.ec:2272-2275) is
                // `(1 + k) - C'` and `b.y + k` is `A' + k`, from the same packed whole
                // numbers: `subpd` at `0x2a3cf`, `0x2a3f9`, `0x2a41a` and `0x2a42b`, `addpd`
                // at `0x2a3d7`, `0x2a3f5`, `0x2a411` and `0x2a3db`; `b` itself is the `subsd`
                // at `0x2a451`. The eC's order gives another double at levels 3, 9, 15 and 19.
                b: turns(&|k| ((1.0 + k) - c_p, a_p + k)),
                // `d + k`, a packed sum with `(k, k)`: `0x2a599`, `0x2a5f6`, `0x2a699` and
                // `0x2a4ff`.
                d: turns(&|k| (d.0 + k, d.1 + k)),
                // faithful: site ON-4. At the first turn alone the whole number goes inside
                // the bracket, `ab + 1` being `(1 + 1) + ab_off.x` and `1 + ab_off.y`, the
                // packed sum of the offset with `(2, 1)` at `0x2a4cb`; at the second, third
                // and fourth it is the eC's own `ab + k` (RI7H.ec:2273-2275), the stored `ab`
                // and `(k, k)` at `0x2a579`, `0x2a60c` and `0x2a658`. The eC's `ab + 1` gives
                // another double at levels 11, 13 and 19. The `NoAlloc` copy takes the
                // bracketed form at all four turns, and there the two copies part, at level 5
                // alone: at the fourth turn of `0000000`, by one unit in the last place.
                ab: turns(&|k| {
                    if k == 1.0 {
                        (2.0 + ab_off.0, 1.0 + ab_off.1)
                    } else {
                        (ab.0 + k, ab.1 + k)
                    }
                }),
                pole: [(5.0, 4.0), (1.0, 0.0)],
            })
        }
        // The south pole at an odd level (RI7H.ec:2289-2313), compiled at `0x29584`.
        (0xB, true) => {
            // `(B', C')`, the packed product of `oonp` and `.rodata` `0x49c00` at `0x29584`,
            // and `A'` at `0x295f7`.
            let a_p = oonp * ODD_HEX_A;
            let b_p = oonp * ODD_HEX_B;
            let c_p = oonp * ODD_HEX_C;
            let a = (4.0 + b_p, 6.0 + c_p);
            let b = (4.0 + c_p, 6.0 - a_p);
            // `ab = a + (b - a) * r` as the eC writes it (RI7H.ec:2293), the sum's terms
            // exchanged: `subpd` at `0x29653`, `mulpd` at `0x29669`, `addpd` at `0x29683`.
            let ab = (a.0 + (b.0 - a.0) * r, a.1 + (b.1 - a.1) * r);
            // faithful: site OS-2. `c` is not formed here either: in `d = b + (c - b) * r`
            // (RI7H.ec:2294-2295), `c.x - b.x` is the eC's own `(4 - A') - b.x` (`0x2962f`),
            // and `c.y - b.y`, which is `(6 - B') - (6 - A')`, is `A' - B'` (`0x29644`); the
            // `mulpd` at `0x29664` and the `addpd` at `0x2967e` follow. It is the source's
            // double at each of the ten odd levels; it is ported as compiled all the same.
            let d = (((4.0 - a_p) - b.0) * r + b.0, (a_p - b_p) * r + b.1);
            Some(PolarSides {
                // faithful: site OS-1. As at the even level (ES-1), `b - k` (RI7H.ec:2297-2300) is
                // `((4 - k) + C', (6 - k) - A')`: `addsd` at `0x296d1`, `0x2977a` and
                // `0x29821`, `subsd` at `0x2966e`, `0x2960f`, `0x2982d` and `0x298da`, and at
                // the fourth turn `C'` itself (`0x2989d`). The eC's order gives another double
                // at every odd level.
                b: turns(&|k| ((4.0 - k) + c_p, (6.0 - k) - a_p)),
                // `d - k` and `ab - k`, packed sums with `(-k, -k)`: `0x29710`, `0x297eb`,
                // `0x29899` and `0x29698` for `d`; `0x296bd`, `0x2976d`, `0x29814` and
                // `0x298b1` for `ab`.
                d: turns(&|k| (d.0 - k, d.1 - k)),
                ab: turns(&|k| (ab.0 - k, ab.1 - k)),
                pole: [(0.0, 2.0), (4.0, 6.0)],
            })
        }
        _ => None,
    }
}

/// The zone's refined boundary in the 5x6 plane: the allocating
/// `I7HZone::getBaseRefinedVertices(crs84, nDivisions)` (`RI7H.ec:2167-2331`), compiled at
/// `0x29220` to `0x2a857` in `libdggal.so`, BuildID
/// `e75d6ab18f8b460713ebcd21df6f07922fb909c3`. Every edge is divided by
/// [`add_intermediate_points`]: the edges of a polar pentagon from [`polar_pentagon_sides`], and
/// those of every other zone as [`add_non_polar_vertices_refined`] walks them.
///
/// It is the ring from which the engine's refined boundary in WGS84 is projected (`crs84`
/// true), and, with `crs84` false, the ring its planar accessor `getZoneRefinedCRSVertices`
/// answers in the 5x6 CRS (`RI7H.ec:3887`). The two differ where
/// [`add_intermediate_points`] says, and in the fifth side of a polar pentagon: towards WGS84
/// it is divided as the other four are, and in the 5x6 CRS it is closed through the two images
/// of the pole.
///
/// [`base_refined_vertices`] ports the other copy, `getBaseRefinedVerticesNoAlloc`, at one
/// division and odd levels, for `containsPoint`; the two share the walker and nothing else.
///
/// The zone must have geometry. A Z7 identifier of level 20 is the null zone to the eC's `to7H`
/// (`RI7H_Z7.ec:348-352`), which no `Z` of this crate stands for: the caller answers for such
/// an identifier before it asks for this ring.
fn get_base_refined_vertices(z: &Z, crs84: bool, n_divisions: i32) -> Vec<(f64, f64)> {
    // `1.0 / (7 * p)` (RI7H.ec:2174): `divsd` at `0x292dd`, the divisor converted as the
    // `uint64` it is at `0x292b6`. `POW7` is a table read of `levelI49R` (`0x292a7`).
    let oonp = 1.0 / (7 * pow7(z.l49r)) as f64;

    let polar = if z.sub_hex == 0 {
        // Even level (RI7H.ec:2181-2253): the test of `subHex` at `0x292fa`.
        polar_pentagon_sides(z.root, false, oonp)
    } else if z.root > 9 && z.sub_hex == 1 {
        // Polar pentagons at an odd level (RI7H.ec:2261-2314): the tests at `0x293f0` and
        // `0x293f5`. A root beyond the two poles adds nothing, as in the eC (`0x2957e`); this
        // crate builds no such zone.
        let Some(sides) = polar_pentagon_sides(z.root, true, oonp) else {
            return Vec::new();
        };
        Some(sides)
    } else {
        None
    };

    let Some(sides) = polar else {
        // Every other zone (RI7H.ec:2240-2252 and :2315-2328): the offsets of the six vertices
        // from the centroid, scaled by `oonp` (`0x2932b` to `0x293be` at an even level,
        // `0x293ff` to `0x29483` at an odd one, where the library scales by `-oonp` and takes
        // the opposite, the same doubles), and the walk about the centroid (`0x293cf`). The
        // centroid is folded as `getVertices` folds its own (RI7H.ec:2176-2179, the same text):
        // the comparisons at `0x292d1` and `0x294c0` with `fl(6 + 1e-9)` and `fl(5 + 1e-9)` of
        // `.rodata` `0x4a128` and `0x4a130`, and at `0x294d2` with zero.
        let table = if z.sub_hex == 0 {
            &EVEN_HEX_VERTS
        } else {
            &ODD_HEX_VERTS
        };
        return add_non_polar_vertices_refined(
            z,
            centroid5x6(z),
            &scale_offsets(table, oonp),
            Emission::Refined { crs84, n_divisions },
        );
    };

    let mut out: Vec<(f64, f64)> = Vec::with_capacity(refined_capacity(n_divisions));
    // The four sides that stay within the layout (RI7H.ec:2197-2200 and its three fellows),
    // each from `b + k` to `b + (k + 1)` across the interruption between `d + k` and
    // `ab + (k + 1)`.
    for k in 0..4 {
        add_intermediate_points(
            &mut out,
            sides.b[k],
            sides.b[k + 1],
            n_divisions,
            Some((sides.d[k], sides.ab[k + 1])),
            crs84,
        );
    }
    // The fifth side, which wraps to the first corner (RI7H.ec:2201-2211): the test of `crs84`
    // at `0x29c8b`, `0x2a0b2`, `0x2a6aa` and `0x298f4` in the four arms.
    if crs84 {
        add_intermediate_points(
            &mut out,
            sides.b[4],
            sides.b[0],
            n_divisions,
            Some((sides.d[4], sides.ab[0])),
            crs84,
        );
    } else {
        out.extend([
            sides.b[4],
            sides.d[4],
            sides.pole[0],
            sides.pole[1],
            sides.ab[0],
        ]);
    }
    out
}

/// The room the eC asks for in the array of a refined boundary before it fills it,
/// `Max(1, nDivisions) * 6` points (`RI7H.ec:2169`). It is a hint and no more: the array grows
/// where the pole rule or an interruption adds points. The hint stops at the 100,000 divisions
/// beyond which the public surface refuses a refinement, so that no argument asks for more
/// memory than its ring will fill.
fn refined_capacity(n_divisions: i32) -> usize {
    6 * n_divisions.clamp(1, crate::grid::MAX_EDGE_REFINEMENT as i32) as usize
}

/// The inside-or-boundary test for the zone's polygon (`containsPoint`,
/// RI7H.ec:951-1080).
///
/// The polygon comes from [`base_refined_vertices`], not from the centroid plus the
/// vertex offset table. Two bounding boxes then do the seam work: the zone's own
/// box, shifted by five when the point lies a whole band away from it, rejects the
/// point outright; and the point's integer rhombus cell is used to skip any side
/// both of whose endpoints fall outside it. Only the surviving sides are tested,
/// and `result` starts false, so a zone all of whose sides were skipped does not
/// contain the point.
///
/// faithful: py4dggs re-expresses this routine as a plain convex half-plane test
/// over the centroid plus the offset table, comparing the side against `-1e-11`
/// rather than against zero and carrying neither bounding box. That is the defect
/// this port no longer reproduces; see the note on [`from_centroid`].
fn contains_point(z: &Z, vx: f64, vy: f64) -> bool {
    let level = zone_level(z);
    // `level` is non-negative, so the integer division matches Python's
    // `math.floor(level / 2)`.
    let p = pow7(level / 2);
    let oonp = 1.0 / (7 * p) as f64;

    let poly = base_refined_vertices(z, oonp);
    if poly.is_empty() {
        return false;
    }

    // The point's own integer rhombus cell. The eC truncates with
    // `(int)(v.x + 1E-11)`; the coordinate is non-negative here, so flooring agrees.
    let pbb_tl = ((vx + 1e-11).floor(), (vy + 1e-11).floor());
    let pbb_br = (pbb_tl.0 + 1.0, pbb_tl.1 + 1.0);

    let mut bbox_tl = (100.0f64, 100.0f64);
    let mut bbox_br = (-100.0f64, -100.0f64);
    for &(x, y) in &poly {
        if x > bbox_br.0 {
            bbox_br.0 = x;
        }
        if y > bbox_br.1 {
            bbox_br.1 = y;
        }
        if x < bbox_tl.0 {
            bbox_tl.0 = x;
        }
        if y < bbox_tl.1 {
            bbox_tl.1 = y;
        }
    }

    if vx - bbox_br.0 > 3.0 && vy - bbox_br.1 > 3.0 && vx - bbox_tl.0 > 3.0 && vy - bbox_tl.1 > 3.0
    {
        bbox_tl.0 += 5.0;
        bbox_tl.1 += 5.0;
        bbox_br.0 += 5.0;
        bbox_br.1 += 5.0;
    }
    // The threshold the engine compares `v.x` against. On every path but the one
    // below it is the plain `bbox.tl.x - 1e-11` of RI7H.ec:1003.
    let mut tl_x_limit = None;
    if vx - bbox_br.0 < -3.0
        && vy - bbox_br.1 < -3.0
        && vx - bbox_tl.0 < -3.0
        && vy - bbox_tl.1 < -3.0
    {
        // faithful: on this one path, and for `bbox.tl.x` alone, the compiler fused
        // the wrap of RI7H.ec:997 with the threshold of RI7H.ec:1003, so that
        // `(tl.x - 5) - 1e-11` is evaluated as `tl.x - fl(5 + 1e-11)`. That is a
        // re-association, which needs `-fassociative-math`, and DGGAL v0.0.6 is
        // built with `-O2 -ffast-math` (`Makefile.dggal:133`) while `containsPoint`
        // (RI7H.ec:951) carries no
        // `__attribute__ ((optimize("-fno-unsafe-math-optimizations")))`. In
        // `libdggal.so`, BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`:
        //
        // ```text
        //   2346d  movsd  xmm0,[rip+0x266db]  ; 49b50 = 5.0
        //   23475  subsd  xmm2,[rip+0x26d0b]  ; 4a188 = 5.00000000001, tl.x
        //   2347d  subsd  xmm5,xmm0           ; tl.y -= 5
        //   23481  subsd  xmm4,xmm0           ; br.x -= 5
        //   23485  subsd  xmm3,xmm0           ; br.y -= 5
        //   23489  jmp    23148               ; past the `subsd xmm2,1e-11` at 23140
        // ```
        //
        // The asymmetry is the engine's and is reproduced as it stands: only
        // `tl.x`, and only on this path. The `+5` wrap at `0x23418` adds five to
        // all four coordinates with a plain `addsd` and then falls into the
        // ordinary threshold test, and `tl.y`, `br.x` and `br.y` keep the plain
        // subtraction and addition of `1e-11` at `0x23154` to `0x23179`. The moved
        // threshold is worth up to one unit in the last place of five, about
        // 8.9e-16, and the constant is pinned by bit pattern in the tests below.
        tl_x_limit = Some(bbox_tl.0 - (5.0 + 1e-11));
        bbox_tl.0 -= 5.0;
        bbox_tl.1 -= 5.0;
        bbox_br.0 -= 5.0;
        bbox_br.1 -= 5.0;
    }

    if vx < tl_x_limit.unwrap_or(bbox_tl.0 - 1e-11)
        || vy < bbox_tl.1 - 1e-11
        || vx > bbox_br.0 + 1e-11
        || vy > bbox_br.1 + 1e-11
    {
        return false;
    }

    let mut result = false;
    for i in 0..poly.len() {
        let j = if i < poly.len() - 1 { i + 1 } else { 0 };
        let (mut ax, mut ay) = poly[i];
        let (mut bx, mut by) = poly[j];
        if (ax - vx).abs() > 3.0 && (ay - vy).abs() > 3.0 {
            if ax > 3.0 && ay > 3.0 {
                ax -= 5.0;
                ay -= 5.0;
            } else {
                ax += 5.0;
                ay += 5.0;
            }
        }
        if (bx - vx).abs() > 3.0 && (by - vy).abs() > 3.0 {
            if bx > 3.0 && by > 3.0 {
                bx -= 5.0;
                by -= 5.0;
            } else {
                bx += 5.0;
                by += 5.0;
            }
        }

        // Skip a side both of whose endpoints lie outside the point's own cell.
        if (ax - 1e-11 < pbb_tl.0
            || ay - 1e-11 < pbb_tl.1
            || ax + 1e-11 > pbb_br.0
            || ay + 1e-11 > pbb_br.1)
            && (bx - 1e-11 < pbb_tl.0
                || by - 1e-11 < pbb_tl.1
                || bx + 1e-11 > pbb_br.0
                || by + 1e-11 > pbb_br.1)
        {
            continue;
        }

        if point_line_side(vx, vy, ax, ay, bx, by) < 0.0 {
            result = false;
            break;
        }
        // At least one side has to have been tested; the bounding box alone yields
        // false positives.
        result = true;
    }
    result
}

/// Traces the cell boundary of a non-polar hexagon or pentagon, rotating the offset
/// around the centroid and re-expressing it across every rhombus interruption it
/// crosses (RI7H.ec `getVertices`, the non-polar path). `v` is the scaled vertex
/// offset table; the boundary points come back in 5x6 space.
///
/// The crossing branch ports the eC's `addNonPolarBaseVertices` faithfully,
/// including the four i1/i2 reconciliation blocks (RI7H.ec:1991-2004) that the
/// igeo7-py oracle from which py4dggs descends had dropped. Without them the
/// far-side vertices of an interruption-spanning cell land outside the 5x6 band and
/// degenerate on the inverse projection. This is py4dggs's one deliberate
/// divergence from that oracle, made in order to match DGGAL itself.
fn add_non_polar_base_vertices(
    cx: f64,
    cy: f64,
    n_points: i64,
    v: &[(f64, f64); 6],
) -> Vec<(f64, f64)> {
    // faithful: the eC forms `t.x = c.x + v[i].x` and then floors `t.x + 1E-11`
    // (RI7H.ec:1959-1960). The compiler hoisted `c.x + 1e-11` out of the loop and
    // added the offset to it, so the sum is bracketed the other way round;
    // `t.x` itself is never formed, since nothing else uses it. See the citation
    // on the direction below for why the engine, and not the eC's text, decides
    // this function.
    let cx_eps = cx + 1e-11;
    let mut start: usize = 0;
    for (i, vi) in v.iter().enumerate() {
        let ty = cy + vi.1;
        let itx = (vi.0 + cx_eps).floor();
        if !(ty - itx > 2.0 || ty < itx) {
            start = i;
            break;
        }
    }

    let mut point = (cx + v[start].0, cy + v[start].1);
    // `start` is in 0..6, so the left operand of the remainder is positive.
    let prev = (start + 5) % 6;
    // faithful: the eC writes `dir = point - (c + v[prev])` (RI7H.ec:1970), with
    // `point` already `c + v[start]`, so the centroid appears on both sides; the
    // compiler cancelled it and subtracted the two offsets directly. That is a
    // re-association, which needs `-fassociative-math`, and DGGAL v0.0.6 is built
    // with `-O2 -ffast-math` (`Makefile.dggal:133`) while `addNonPolarBaseVertices`
    // (RI7H.ec:1949) carries no
    // `__attribute__ ((optimize("-fno-unsafe-math-optimizations")))`. In
    // `libdggal.so`, BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, the
    // function occupies `0x1ebe0` to `0x1efe4`, and its head reads:
    //
    // ```text
    //   1ec10  addsd  xmm7,xmm2        ; xmm7 = 1e-11 + c.x, hoisted
    //   1ec35  addsd  xmm3,xmm4        ; t.y  = c.y + v[i].y
    //   1ec39  addsd  xmm1,xmm7        ; v[i].x + (c.x + 1e-11)
    //   1ec91  addsd  xmm6,xmm4        ; point.y = c.y + v[start].y
    //   1ec9c  addsd  xmm2,xmm5        ; point.x = c.x + v[start].x
    //   1eca0  subsd  xmm5,[rsi]       ; dir.x = v[start].x - v[prev].x
    //   1eca4  subsd  xmm4,[rsi+0x8]   ; dir.y = v[start].y - v[prev].y
    // ```
    //
    // The difference is one unit in the last place at a time, and it enters the
    // walk before the first turn, so every later vertex carries it. This is the
    // arithmetic whose absence left the boundary class open after the side test
    // was ported.
    let mut direction = (v[start].0 - v[prev].0, v[start].1 - v[prev].1);

    // Sized once, as the eC's buffer is fixed (`Pointd v5x6[6]` in
    // `getZoneWGS84Vertices`, RI7H.ec:456): one point a side.
    let mut vertices: Vec<(f64, f64)> = Vec::with_capacity(n_points as usize);
    for i in start..start + n_points as usize {
        let ndir = rotate5x6_offset(direction.0, direction.1, false); // 60 anticlockwise
        let mut n = (point.0 + ndir.0, point.1 + ndir.1);

        let mut p = point;
        if p.0 > 5.0 && p.1 > 5.0 {
            p.0 -= 5.0;
            p.1 -= 5.0;
        }
        if p.0 < 0.0 || p.1 < 0.0 {
            p.0 += 5.0;
            p.1 += 5.0;
        }

        let cross = crosses5x6_interruption_v2(p.0, p.1, ndir.0, ndir.1);
        if let Some(cross) = cross {
            // Reconcile the crossing in and out points (i1 is i_src, i2 is i_dst)
            // into a consistent frame before using them: first shift both into the
            // unwrapped `point` frame, then bring i2 into i1's frame. Dropping these
            // four corrections (eC RI7H.ec:1991-2004) is the igeo7-py port defect
            // that gave an interruption-spanning cell a wrong or degenerate
            // far-side vertex, because the remainder direction and the new point
            // `i2 + d` were then computed in a frame off by five.
            let mut i1 = cross.i_src;
            let mut i2 = cross.i_dst;
            let north = cross.in_north;
            if point.0 - p.0 > 4.0 {
                i1.0 += 5.0;
                i1.1 += 5.0;
                i2.0 += 5.0;
                i2.1 += 5.0;
            }
            if p.0 - point.0 > 4.0 {
                i1.0 -= 5.0;
                i1.1 -= 5.0;
                i2.0 -= 5.0;
                i2.1 -= 5.0;
            }
            if i2.1 - i1.1 > 4.0 {
                i2.0 -= 5.0;
                i2.1 -= 5.0;
            }
            if i1.1 - i2.1 > 4.0 {
                i2.0 += 5.0;
                i2.1 += 5.0;
            }

            let crossing_left = if north { i2.0 < i1.0 } else { i2.0 > i1.0 };
            // faithful: the eC rotates the remainder of the step,
            // `rotate5x6Offset(d, dir.x - (i1.x - point.x), dir.y - (i1.y - point.y),
            // !crossingLeft)` (RI7H.ec:2008), then sets `n = i2 + d` and rotates
            // `dir` a second time (RI7H.ec:2010). The compiler re-associated the
            // remainder, forming it from the uncrossed next point rather than from
            // the direction, and then removed the second rotation altogether, since
            // under `-fassociative-math` a turn one way followed by a turn the
            // other is the identity. At `0x1ebe0` in `libdggal.so`, BuildID
            // `e75d6ab18f8b460713ebcd21df6f07922fb909c3`:
            //
            // ```text
            //   1edc1  movapd xmm9,xmm1    ; ndir.y
            //   1edc6  addsd  xmm9,xmm3    ; n.y = ndir.y + point.y (from 1ee33 for x)
            //   -- crossingLeft false, the clockwise arm --
            //   1edd3  movapd xmm1,xmm8    ; n.x
            //   1edd8  subsd  xmm5,xmm9    ; i1.y - n.y, that is -ddy
            //   1eddd  subsd  xmm1,xmm7    ; ddx = n.x - i1.x
            //   1ede1  addsd  xmm5,xmm1    ; d.x = (-ddy) + ddx
            //   1ede5  addsd  xmm0,xmm1    ; n.y = i2.y + ddx
            //   1ede9  addsd  xmm5,xmm6    ; n.x = d.x + i2.x
            //   (dir is written nowhere on this arm)
            //   -- crossingLeft true, the anticlockwise arm --
            //   1ef40  subsd  xmm9,xmm5    ; ddy = n.y - i1.y
            //   1ef45  subsd  xmm7,xmm8    ; i1.x - n.x, that is -ddx
            //   1ef4a  movsd  xmm4,[rsp]   ; dir.y := old dir.x, then negated at 1ef54
            //   1ef4f  movsd  [rsp],xmm1   ; dir.x := ndir.y
            //   1ef5c  addsd  xmm0,xmm9    ; i2.y + ddy
            //   1ef66  addsd  xmm8,xmm6    ; n.x = ddy + i2.x
            //   1ef6b  addsd  xmm0,xmm7    ; n.y = (i2.y + ddy) - ddx
            // ```
            //
            // So `ddx` and `ddy` are the offsets of the uncrossed next point from
            // `i1`, not the direction less the offset of `i1` from `point`; the
            // clockwise arm leaves `dir` untouched where the source text would
            // recompute it as `dir.y - (dir.y - dir.x)`; and the anticlockwise arm
            // writes `(ndir.y, -dir.x)` where the source text would compute
            // `(dir.y - dir.x) - dir.y`. Each of those three is exact in real
            // arithmetic and not in IEEE.
            let ddx = n.0 - i1.0;
            if crossing_left {
                let ddy = n.1 - i1.1;
                n = (ddy + i2.0, (i2.1 + ddy) - ddx);
                direction = (ndir.1, -direction.0);
            } else {
                // `-ddy` as the instruction at `0x1edd8` forms it. IEEE subtraction
                // is exactly antisymmetric, so this is the exact negation of the
                // `ddy` of the other arm, and `(-ddy) + ddx` is exactly `ddx - ddy`.
                let minus_ddy = i1.1 - n.1;
                n = ((minus_ddy + ddx) + i2.0, i2.1 + ddx);
                // `direction` is deliberately left as it was: see above.
            }
            vertices.push(point);
            point = n;
        } else {
            // Inert by construction, since the loop header already bounds `i`, and
            // faithful: the eC carries the identical dead guard,
            // `else if(i < start + nPoints)` at RI7H.ec:2014 against the loop
            // `for(i = start + 0; i < start + nPoints; i++)` at :1973. Kept by the
            // byte-faithfulness convention; it must not be simplified away.
            if i < start + n_points as usize {
                vertices.push(point);
            }
            direction = ndir;
            point = n;
        }
    }
    vertices
}

// --------------------------------------------------------------------------- //
// Hierarchy: candidate parents, children, the parent walk (RI7H.ec)
// --------------------------------------------------------------------------- //

/// Maps a perturbed `(row, col)` back to the even-level zone that contains it,
/// wrapping across rhombus seams and poles (`calcCandidateParent`, `RI7H.ec:1367-1429`).
/// Returns `None` when the perturbed cell falls outside every rhombus.
fn calc_candidate_parent(
    l49r: i64,
    root: i64,
    row: i64,
    col: i64,
    add_col: i64,
    add_row: i64,
) -> Option<Z> {
    let p = pow7(l49r);
    let mut r = root;
    let mut c = col + add_col;
    let mut rw = row + add_row;
    let south = r & 1 != 0;

    if c == p && rw < p && !south {
        c = p - rw;
        rw = 0;
        r += 2;
    } else if rw == p && c < p && south {
        rw = p - c;
        c = 0;
        r += 2;
    } else if rw < 0 && c < 0 {
        rw += p;
        c += p;
        r -= 2;
    } else if rw < 0 {
        rw += p;
        r -= 1;
    } else if c < 0 {
        c += p;
        r -= 1;
    } else if c >= p && rw >= p {
        rw -= p;
        c -= p;
        r += 2;
    } else if rw >= p {
        rw -= p;
        r += 1;
    } else if c >= p {
        c -= p;
        r += 1;
    }
    if r < 0 {
        r += 10;
    } else if r > 9 {
        r -= 10;
    }

    let south = r & 1 != 0;
    if !south && rw == 0 && c == p {
        return Some(Z::new(l49r, 0xA, 0, 0, 0));
    }
    if south && rw == p && c == 0 {
        return Some(Z::new(l49r, 0xB, 0, 0, 0));
    }

    if rw < 0 || rw >= p || c < 0 || c >= p {
        return None;
    }
    Some(Z::new(l49r, r, rw, c, 0))
}

/// The root, row and column of the even-level centroid child of an odd-level zone,
/// that is the aperture-49 child sitting on the parent's own centroid
/// (`getOddLevelCentroidChildRootRowCol`, RI7H.ec). It handles the poles, the south
/// corner and the edge-hexagon sub-hexagon rotations, and then the seam wrapping at
/// the finer `cp` scale.
#[derive(Clone, Copy, Debug)]
struct CentroidChild {
    root: i64,
    row: i64,
    col: i64,
    cp: i64,
}

/// Ports `getOddLevelCentroidChildRootRowCol` (`RI7H.ec:2501-2608`); see `CentroidChild`.
fn get_odd_level_centroid_child_root_row_col(z: &Z) -> Option<CentroidChild> {
    let p = pow7(z.l49r);
    let cp = p * 7;
    let root = z.root;

    if root == 0xA {
        if z.sub_hex == 1 {
            return Some(CentroidChild {
                root: 0xA,
                row: 0,
                col: 0,
                cp,
            });
        }
        return Some(CentroidChild {
            root: 2 * (z.sub_hex - 2),
            row: 1,
            col: cp - 2,
            cp,
        });
    }
    if root == 0xB {
        if z.sub_hex == 1 {
            return Some(CentroidChild {
                root: 0xB,
                row: 0,
                col: 0,
                cp,
            });
        }
        return Some(CentroidChild {
            root: 9 - 2 * (z.sub_hex - 2),
            row: cp - 1,
            col: 2,
            cp,
        });
    }

    let mut sh = z.sub_hex;
    let mut c_rhombus = root;
    let south = c_rhombus & 1 != 0;
    let mut row = 7 * z.row;
    let mut col = 7 * z.col;

    if z.row == 0 && z.col == 0 && south && sh >= 4 {
        sh += 1;
    } else if sh >= 2 {
        let parent_z = Z::new(z.l49r, z.root, z.row, z.col, 0);
        if zone_is_edge_hex(&parent_z) {
            // `sh >= 2` on this arm, so the left operand of the remainder is at
            // least one and a plain `%` matches Python's own result.
            sh = (sh + if south { -1 } else { 3 }) % 6 + 2;
        }
    }

    match sh {
        2 => {
            row -= 3;
            col -= 1;
        }
        3 => {
            row -= 2;
            col -= 3;
        }
        4 => {
            row += 1;
            col -= 2;
        }
        5 => {
            row += 3;
            col += 1;
        }
        6 => {
            row += 2;
            col += 3;
        }
        7 => {
            row -= 1;
            col += 2;
        }
        _ => {}
    }

    if col == cp && row < cp && !south {
        col = cp - row;
        row = 0;
        c_rhombus += 2;
    } else if row == cp && col < cp && south {
        row = cp - col;
        col = 0;
        c_rhombus += 2;
    } else if col > 0 && col < cp && row < 0 && !south {
        let ncol = cp + row;
        let nrow = cp - col;
        col = ncol;
        row += nrow;
        c_rhombus -= 2;
    } else if row > 0 && row < cp && col < 0 && south {
        let nrow = cp + col;
        let ncol = cp - row;
        row = nrow;
        col += ncol;
        c_rhombus -= 2;
    }

    // Python runs this as a second, separate `if` chain rather than as the `else`
    // arm of the one above; the other three wrapping sites in this file do use an
    // `else`, so the difference is preserved deliberately.
    if row < 0 && col < 0 {
        row += cp;
        col += cp;
        c_rhombus -= 2;
    } else if row < 0 {
        row += cp;
        c_rhombus -= 1;
    } else if col < 0 {
        col += cp;
        c_rhombus -= 1;
    } else if col >= cp && row >= cp {
        row -= cp;
        col -= cp;
        c_rhombus += 2;
    } else if row >= cp {
        row -= cp;
        c_rhombus += 1;
    } else if col >= cp {
        col -= cp;
        c_rhombus += 1;
    }

    if c_rhombus < 0 {
        c_rhombus += 10;
    } else if c_rhombus > 9 {
        c_rhombus -= 10;
    }

    if (0..cp).contains(&row) && (0..cp).contains(&col) {
        return Some(CentroidChild {
            root: c_rhombus,
            row,
            col,
            cp,
        });
    }
    None
}

/// A zone's primary children, seven at most, held inline as the eC holds them: every
/// caller of `getPrimaryChildren(I7HZone children[7])` (`RI7H.ec:2721`) passes an
/// array of seven on its own stack (`RI7H.ec:1535`, `RI7H.ec:1747`,
/// `RI7H_Z7.ec:140`, `RI7H_Z7.ec:235`, `RI7H_Z7.ec:356`). The geometric parent search
/// asks for some eight such lists on every call, and allocating each on the heap was a
/// large part of its cost. It reads as a slice of the zones pushed, in the order in
/// which they were pushed, and iterates by value, as the vector it replaces did.
#[derive(Debug)]
struct Children {
    zones: [Z; 7],
    len: usize,
}

impl Children {
    const EMPTY: Children = Children {
        zones: [Z {
            l49r: 0,
            root: 0,
            row: 0,
            col: 0,
            sub_hex: 0,
        }; 7],
        len: 0,
    };

    /// No eC counterpart: the eC appends with `children[count++] = ...` at each site.
    fn push(&mut self, z: Z) {
        debug_assert!(
            self.len < 7,
            "a zone has at most seven primary children, and an eighth was pushed"
        );
        self.zones[self.len] = z;
        self.len += 1;
    }
}

impl std::ops::Deref for Children {
    type Target = [Z];

    /// No eC counterpart: the eC reads the array together with the count returned beside it.
    fn deref(&self) -> &[Z] {
        &self.zones[..self.len]
    }
}

impl IntoIterator for Children {
    type Item = Z;
    type IntoIter = std::iter::Take<std::array::IntoIter<Z, 7>>;

    /// No eC counterpart: the eC loops over the array up to the count returned beside it.
    fn into_iter(self) -> Self::IntoIter {
        self.zones.into_iter().take(self.len)
    }
}

/// The six or seven primary children of a zone (`getPrimaryChildren`, `RI7H.ec:2721-2870`).
/// Even to odd simply enumerates the sub-hexagons one to seven, six of them for the central
/// pentagon. Odd to even places the centroid child plus the ring of six offset
/// children, wrapping each across the rhombus seams, with the polar fan-out and the
/// edge-hexagon rotation.
fn get_primary_children(z: &Z) -> Children {
    let (l49r, root, z_row, z_col, sub_hex) = (z.l49r, z.root, z.row, z.col, z.sub_hex);
    // Into an array of seven, as the eC's callers pass one, and never on the heap.
    let mut children = Children::EMPTY;
    if l49r > 9 || (l49r == 9 && sub_hex != 0) {
        return children;
    }

    if sub_hex == 0 {
        for sh in 1..7 {
            children.push(Z::new(l49r, root, z_row, z_col, sh));
        }
        if z_row != 0 || z_col != 0 {
            children.push(Z::new(l49r, root, z_row, z_col, 7));
        }
        return children;
    }

    let Some(cc) = get_odd_level_centroid_child_root_row_col(z) else {
        return children;
    };
    let (c_root, crow, ccol, cp) = (cc.root, cc.row, cc.col, cc.cp);
    let c_child = Z::new(l49r + 1, c_root, crow, ccol, 0);
    children.push(c_child);

    if c_root == 0xA {
        for i in 0..5 {
            children.push(Z::new(l49r + 1, 2 * i, 0, cp - 1, 0));
        }
    } else if c_root == 0xB {
        for i in 0..5 {
            children.push(Z::new(l49r + 1, 9 - 2 * i, cp - 1, 0, 0));
        }
    } else {
        let n_points = zone_npoints(z);
        let south = c_root & 1 != 0;
        let c_offsets: [(i64, i64); 6] = [(-1, 0), (-1, -1), (0, -1), (1, 0), (1, 1), (0, 1)];

        let mut edge_hex_fix = false;
        if n_points == 6 && sub_hex == 1 {
            let parent_z = Z::new(l49r, root, z_row, z_col, 0);
            if zone_is_edge_hex(&parent_z) {
                edge_hex_fix = true;
            }
        }

        for i in 0..6usize {
            // `i` is in 0..6 and the addend is positive, so the left operand of the
            // remainder is positive.
            let ii = if edge_hex_fix {
                (i + if south { 1 } else { 5 }) % 6
            } else {
                i
            };

            if n_points == 5 {
                if south && i == 2 {
                    continue;
                }
                if !south && i == 0 {
                    continue;
                }
            }

            let mut row = crow + c_offsets[ii].0;
            let mut col = ccol + c_offsets[ii].1;
            let mut c_rhombus = c_root;

            if col == cp && row < cp && !south {
                col = cp - row;
                row = 0;
                c_rhombus += 2;
            } else if row == cp && col < cp && south {
                row = cp - col;
                col = 0;
                c_rhombus += 2;
            } else if col > 0 && col < cp && row < 0 && !south {
                let ncol = cp + row;
                let nrow = cp - col;
                col = ncol;
                row += nrow;
                c_rhombus -= 2;
            } else if row > 0 && row < cp && col < 0 && south {
                let nrow = cp + col;
                let ncol = cp - row;
                row = nrow;
                col += ncol;
                c_rhombus -= 2;
            } else if row < 0 && col < 0 {
                row += cp;
                col += cp;
                c_rhombus -= 2;
            } else if row < 0 {
                row += cp;
                c_rhombus -= 1;
            } else if col < 0 {
                col += cp;
                c_rhombus -= 1;
            } else if col >= cp && row >= cp {
                row -= cp;
                col -= cp;
                c_rhombus += 2;
            } else if row >= cp {
                row -= cp;
                c_rhombus += 1;
            } else if col >= cp {
                col -= cp;
                c_rhombus += 1;
            }

            if c_rhombus < 0 {
                c_rhombus += 10;
            } else if c_rhombus > 9 {
                c_rhombus -= 10;
            }

            if (0..cp).contains(&row) && (0..cp).contains(&col) {
                children.push(Z::new(l49r + 1, c_rhombus, row, col, 0));
            }
        }
    }
    children
}

/// The single child sharing the parent's centroid: sub-hexagon one for an even
/// zone, the even-level centroid child for an odd zone. Used by the rotation
/// offset. Ports the `centroidChild` property (`RI7H.ec:2610-2632`), with its refusal
/// (`RI7H.ec:2614`): an odd zone at `l49r` 9 is at level 19, the deepest, and has no
/// child. No caller asks at that level.
fn get_centroid_child(z: &Z) -> Option<Z> {
    if z.l49r == 9 && z.sub_hex != 0 {
        return None;
    }
    if z.sub_hex == 0 {
        return Some(Z::new(z.l49r, z.root, z.row, z.col, 1));
    }
    let result = get_odd_level_centroid_child_root_row_col(z)?;
    Some(Z::new(z.l49r + 1, result.root, result.row, result.col, 0))
}

/// Finds the odd-level parent of an even-level zone by inverting the centroid
/// geometry: quantise the child's centroid back to a rhombus cell, enumerate the
/// candidate parents, and return the one whose grandchildren include `child` (the
/// even-level parent search in RI7H.ec). Returns `None` when no parent matches.
fn from_even_level_primary_child(child: &Z) -> Option<Z> {
    let l49r = child.l49r - 1;
    if child.sub_hex != 0 || l49r < 0 {
        return None;
    }

    let (cx, cy) = zone_centroid(child);
    let (mut x, mut y) = (cx, cy);
    let p = pow7(l49r);
    let oop = 1.0 / p as f64;

    if (x - y - 1.0).abs() < 1e-10 {
        // north pole diagonal: leave the point alone
    } else if (y - x - 2.0).abs() < 1e-10 {
        // south pole diagonal: leave the point alone
    } else if y < -1e-11 && x > -1e-11 {
        x -= y;
        y = 0.0;
    } else if (x + 1e-11).floor() as i64 > (y + 1e-11).floor() as i64 {
        let iy = ((y + 1e-11).floor() as i64).min(5);
        x += iy as f64 + 1.0 - y;
        y = iy as f64 + 1.0;
    } else if ((y + 1e-11).floor() as i64) - ((x + 1e-11).floor() as i64) > 1 {
        let ix = ((x + 1e-11).floor() as i64).min(4);
        y += ix as f64 + 1.0 - x;
        x = ix as f64 + 1.0;
    } else if x < -1e-11 || y < -1e-11 {
        x += 5.0;
        y += 5.0;
    }

    // faithful: the bottom-right corner test, as the engine evaluates it. The eC
    // writes `c.x + c.y > 5.0 + 5.0 - oop - 1E-11` (RI7H.ec:1722), which is
    // `(10 - oop) - 1e-11`; the compiler folded the two constants together and
    // subtracted the reciprocal last, giving `fl(10 - 1e-11) - oop`. The
    // reasoning and the citation are those of the identical site in
    // [`from_centroid`]; this inlining sits at `0x1fc54` in `libdggal.so`, BuildID
    // `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, and loads the same
    // `0x49f58 = 9.99999999999`.
    // faithful: the ordinate the truncation below is given, which on this path the
    // compiler computes in one subtraction. The reasoning and the citation are those
    // of the identical site in [`from_centroid`], including the demonstration that
    // it cannot change an answer, the guard excluding the one double at which the
    // two forms part company. This inlining sits at `0x1fc82`, and the abscissa
    // keeps the plain `1e-11 + c.x` at `0x1fdc0` on both paths, as it does there.
    // This copy is never reached by any census taken so far: the corner wrap fires
    // only in `from_centroid`, and `from_even_level_primary_child` is entered from
    // the even-level parent walk, which no sample has carried to the corner.
    let mut wrapped_y = None;
    if x > 5.0 - 1e-11 && y > 5.0 - 1e-11 && x + y > (10.0 - 1e-11) - oop {
        wrapped_y = Some(y - (5.0 - 1e-11));
        x -= 5.0;
        y -= 5.0;
    }

    // faithful: the eC truncates here, `Min(4, (int)(c.x + 1E-11))` at
    // RI7H.ec:1726, where py4dggs floors. py4dggs keeps the difference on purpose,
    // and it is inert: the values are non-negative at this site, so flooring and
    // truncating agree, and on the wrapped ordinate above they agree too, since the
    // guard `y > 5 - 1e-11` makes `y - fl(5 - 1e-11)` strictly positive. The inlined
    // domain fold above is not affected, since the eC floors there (RI7H.ec:1445,
    // 1448, 1451, 1454).
    let ix = ((x + 1e-11).floor() as i64).min(4);
    let iy = (wrapped_y.unwrap_or(y + 1e-11).floor() as i64).min(5);
    let root = ix + iy;
    let fx = x - ix as f64;
    let fy = y - iy as f64;
    let col = jsround(fx * p as f64).max(0);
    let row = jsround(fy * p as f64).max(0);
    let south = y - x - 1e-11 > 1.0;
    let north = x - y - 1e-11 > 0.0;
    let north_pole = north && (x - y - 1.0).abs() < 1e-11;
    let south_pole = south && (y - x - 2.0).abs() < 1e-11;

    if north_pole {
        return Some(Z::new(l49r, 0xA, 0, 0, 1));
    }
    if south_pole {
        return Some(Z::new(l49r, 0xB, 0, 0, 1));
    }

    let mut candidate_parents: [Option<Z>; 7] = [None; 7];
    if north && row == 0 && col == p {
        candidate_parents[0] = Some(Z::new(l49r, 0xA, 0, 0, 0));
    } else if south && row == p && col == 0 {
        candidate_parents[0] = Some(Z::new(l49r, 0xB, 0, 0, 0));
    } else {
        candidate_parents[0] = calc_candidate_parent(l49r, root, row, col, 0, 0);
    }
    candidate_parents[1] = calc_candidate_parent(l49r, root, row, col, 0, -1);
    candidate_parents[2] = calc_candidate_parent(l49r, root, row, col, 0, 1);
    candidate_parents[3] = calc_candidate_parent(l49r, root, row, col, 1, 0);
    candidate_parents[4] = calc_candidate_parent(l49r, root, row, col, -1, 0);
    candidate_parents[5] = calc_candidate_parent(l49r, root, row, col, -1, -1);
    candidate_parents[6] = calc_candidate_parent(l49r, root, row, col, 1, 1);

    for cand in candidate_parents.into_iter().flatten() {
        for ch in get_primary_children(&cand) {
            for gch in get_primary_children(&ch) {
                if zone_eq(Some(gch), Some(*child)) {
                    return Some(ch);
                }
            }
        }
    }
    None
}

/// The immediate parent of a zone: drop the sub-hexagon for an odd zone, otherwise
/// find the odd parent of an even zone through
/// [`from_even_level_primary_child`]. Returns `None` at the base level. Ports the
/// `parent0` property (`RI7H.ec:1099-1119`).
fn get_parent0(z: &Z) -> Option<Z> {
    if z.sub_hex > 0 {
        return Some(Z::new(z.l49r, z.root, z.row, z.col, 0));
    }
    let level = zone_level(z);
    if level == 0 {
        return None;
    }
    from_even_level_primary_child(z)
}

/// The full parent chain of `zone`, from its parent to the base zone, as the eC's
/// `computeParents` builds it (`RI7H_Z7.ec:286-296`).
fn compute_parents(zone: &Z) -> Vec<Z> {
    let mut parents: Vec<Z> = Vec::new();
    let mut z = *zone;
    while zone_level(&z) > 0 {
        let Some(parent) = get_parent0(&z) else { break };
        parents.push(parent);
        z = parent;
    }
    parents
}

// --------------------------------------------------------------------------- //
// fromCentroid: quantise a 5x6 point to an I7H zone (RI7H.ec)
// --------------------------------------------------------------------------- //

/// Quantises a 5x6 planar point to the I7H zone whose cell contains it, at `level`
/// (`I7HZone::fromCentroid`).
///
/// An even level snaps to the rhombus lattice, then nudges the `(row, col)` by one
/// cell using the sub-cell offset `(dx, dy)` and the four hexagon-edge half-plane
/// tests. An odd level enumerates the candidate even parents' children by brute
/// force and picks the one whose polygon contains the point; when no candidate
/// contains it, the absent value comes back, which is the eC's nullZone. There is
/// no fallback beyond that, as the eC has none; the comment on that arm explains
/// the recovery py4dggs adds there and why it is not reproduced. The even level
/// gives the absent value too, where the eC does, should its wrapped row or column
/// still fall outside the rhombus.
pub(super) fn from_centroid(level: i64, cx: f64, cy: f64) -> Option<Z> {
    // `level` is non-negative, so the integer division matches Python's
    // `math.floor(level / 2)`.
    let l49r = level / 2;
    let p = pow7(l49r);
    let oop = 1.0 / p as f64;

    let (mut x, mut y) = canonicalize5x6(cx, cy);

    // faithful: the eC writes `c.x + c.y > 5.0 + 5.0 - oop - 1E-11` (RI7H.ec:1461),
    // which is `(10 - oop) - 1e-11`; the compiler re-associated it into
    // `fl(10 - 1e-11) - oop`, folding the two constants together and subtracting
    // the reciprocal last. `fromCentroid` (RI7H.ec:1431) carries no
    // `__attribute__ ((optimize("-fno-unsafe-math-optimizations")))`, and DGGAL
    // v0.0.6 is built with `-O2 -ffast-math` (`Makefile.dggal:133`). In
    // `libdggal.so`, BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`:
    //
    // ```text
    //   235cb  movapd xmm8,xmm4           ; 1.0
    //   235d0  movsd  xmm5,[rip+0x26980]  ; 49f58 = 9.99999999999
    //   235d8  divsd  xmm8,xmm7           ; oop = 1.0 / p
    //   235dd  subsd  xmm5,xmm8           ; (10 - 1e-11) - oop
    //   235e2  movapd xmm8,xmm0
    //   235e7  addsd  xmm8,xmm1           ; c.x + c.y
    //   235ec  comisd xmm8,xmm5
    //   235f1  ja     23a58
    // ```
    //
    // The two `5 - 1E-11` tests above it load the folded `0x49f48`, which is
    // ordinary constant folding and therefore exact. The constant of this site is
    // pinned by bit pattern in the tests below.
    //
    // It cannot change an answer, and the port is for the record. The first two
    // conjuncts force `c.x + c.y > 10 - 2e-11`, while `oop` is `7^-l49r` with `l49r`
    // at most 9, so `oop >= 2.47e-8` and neither form of the threshold ever rises
    // above `10 - 2.47e-8`. The third conjunct therefore cannot fail once the other
    // two hold, at any level this grid admits, and the two groupings of it decide
    // alike. The site is exercised rather than merely present: in the census built at
    // the corner it is reached 6,800 times and taken 6,800 times, and no answer of
    // that census moved when it was ported.
    // faithful: the ordinate that the truncation below is given, on the wrap path
    // alone. The eC takes `Min(5, (int)(c.y + 1E-11))` (RI7H.ec:1472) of the already
    // wrapped `c.y`, that is of `(c.y - 5) + 1e-11`; the compiler fused the two into
    // `c.y - fl(5 - 1e-11)`, a re-association that moves the truncation's argument
    // by up to one unit in the last place of five. At `0x23a58` in `libdggal.so`,
    // BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`:
    //
    // ```text
    //   23a58  movapd xmm11,xmm1   ; keep c.y
    //   23a5d  subsd  xmm1,xmm2    ; xmm2 = 0x49f48 = 4.99999999999
    //   23a7b  subsd  xmm11,xmm5   ; c.y -= 5, with xmm5 = 5.0
    //   23a80  subsd  xmm0,xmm5    ; c.x -= 5
    //   23a84  cvttsd2si edx,xmm1  ; the truncation, on the fused value
    //   23a9e  cmp    edx,eax      ; Min(5, ...) against eax = 5
    // ```
    //
    // The asymmetry is the engine's and is reproduced as it stands: the abscissa
    // keeps the plain `c.x + 1e-11` at `0x23710`, on the wrap path as on every
    // other, and the ordinate is fused on the wrap path alone, the plain form being
    // at `0x23617`. The same fusion sits at `0x1fc82` in `fromEvenLevelPrimaryChild`.
    //
    // It cannot change an answer, and the port is for the record. The two forms
    // were compared at every integer crossing this domain admits, and they flip at
    // the same double at all of them but one: the crossing to zero, where the plain
    // form flips at `0x1.3ffffffffd406p+2` and the fused one at
    // `0x1.3ffffffffd405p+2`, which is `fl(5 - 1e-11)` itself. The guard above,
    // `y > 5 - 1e-11`, excludes exactly that double by a strict comparison against
    // it. The values do part company: over the corner census the fused form differs
    // in its bits from the plain one in 4,780 of 6,800 evaluations. Only the integer
    // they feed never does.
    let mut wrapped_y = None;
    if x > 5.0 - 1e-11 && y > 5.0 - 1e-11 && x + y > (10.0 - 1e-11) - oop {
        wrapped_y = Some(y - (5.0 - 1e-11));
        x -= 5.0;
        y -= 5.0;
    }

    // faithful: the eC truncates here, `Min(4, (int)(c.x + 1E-11))` at
    // RI7H.ec:1471, where py4dggs floors. py4dggs keeps the difference on purpose,
    // and it is inert: the values are non-negative at this site, so flooring and
    // truncating agree, and on the wrapped ordinate above they agree too, since the
    // guard `y > 5 - 1e-11` makes `y - fl(5 - 1e-11)` strictly positive. The
    // `canonicalize5x6` fold above is not affected, since the eC floors there
    // (RI7H.ec:1445, 1448, 1451, 1454).
    let ix = ((x + 1e-11).floor() as i64).min(4);
    let iy = (wrapped_y.unwrap_or(y + 1e-11).floor() as i64).min(5);
    let mut root = ix + iy;
    let fx = x - ix as f64;
    let fy = y - iy as f64;
    let mut col = jsround(fx * p as f64).max(0);
    let mut row = jsround(fy * p as f64).max(0);
    let dx = fx * p as f64 + 0.5 - col as f64;
    let dy = fy * p as f64 + 0.5 - row as f64;
    let south_rhombus = root & 1 != 0;
    let south = y - x - 1e-11 > 1.0;
    let north = x - y - 1e-11 > 0.0;
    let north_pole = north && (x - y - 1.0).abs() < 1e-11;
    let south_pole = south && (y - x - 2.0).abs() < 1e-11;

    if level & 1 != 0 {
        if north_pole {
            return Some(Z::new(l49r, 0xA, 0, 0, 1));
        }
        if south_pole {
            return Some(Z::new(l49r, 0xB, 0, 0, 1));
        }

        let mut candidate_parents: [Option<Z>; 7] = [None; 7];
        if north && row == 0 && col == p {
            candidate_parents[0] = Some(Z::new(l49r, 0xA, 0, 0, 0));
        } else if south && row == p && col == 0 {
            candidate_parents[0] = Some(Z::new(l49r, 0xB, 0, 0, 0));
        } else {
            candidate_parents[0] = calc_candidate_parent(l49r, root, row, col, 0, 0);
        }
        candidate_parents[1] = calc_candidate_parent(l49r, root, row, col, 0, -1);
        candidate_parents[2] = calc_candidate_parent(l49r, root, row, col, 0, 1);
        candidate_parents[3] = calc_candidate_parent(l49r, root, row, col, 1, 0);
        candidate_parents[4] = calc_candidate_parent(l49r, root, row, col, -1, 0);
        candidate_parents[5] = calc_candidate_parent(l49r, root, row, col, -1, -1);
        candidate_parents[6] = calc_candidate_parent(l49r, root, row, col, 1, 1);

        for cand in candidate_parents.into_iter().flatten() {
            for ch in get_primary_children(&cand) {
                if contains_point(&ch, x, y) {
                    return Some(ch);
                }
            }
        }

        // faithful: the eC has no fallback here. `zone` is initialised to nullZone
        // (RI7H.ec:1489) and, when the candidate search finds nothing, is returned
        // unchanged (RI7H.ec:1594). py4dggs instead adds a two-step recovery, first
        // quantising at `level + 1` and then walking back up through
        // `from_even_level_primary_child`, and falling back again to the rhombus
        // corner. That recovery is dropped here, for two reasons.
        //
        // It is wrong: DGGAL really does return nullZone for these points, which
        // lie within roughly a quarter of a degree of a pole at odd resolutions of
        // 15 and above, and py4dggs answers them with a fabricated zone. It is also
        // unsafe: on such a point the recovery reaches `from_7h` with a zone whose
        // parent chain is empty and panics on the index.
        //
        // Both were demonstrated, not assumed. With this arm replaced by
        // `unreachable!`, a sweep of 480,120 quantisations never entered it, but a
        // larger one of 1,500,120 did, and the oracle answers `(null)` for exactly
        // those inputs.
        return None;
    }

    if north_pole {
        return Some(Z::new(l49r, 0xA, 0, 0, 0));
    }
    if south_pole {
        return Some(Z::new(l49r, 0xB, 0, 0, 0));
    }

    // The four tie-breaks below are RI7H.ec:1610, 1616, 1626 and 1632, each of them
    // a `pointLineSide(...) < 0` whose endpoints are constants. See
    // [`tie_break_side_is_negative`] for the order the compiler emitted and the
    // address of each inlining.
    if dx > 1.0 - dy {
        if dx > dy {
            // faithful: 0x239ea in libdggal.so e75d6ab1, replacing RI7H.ec:1610.
            if tie_break_side_is_negative(dx, dy, 1.0, 0.5, 5.0 / 6.0, 1.0 / 6.0) {
                col += 1;
            }
            // faithful: 0x23d80 in libdggal.so e75d6ab1, replacing RI7H.ec:1616.
        } else if tie_break_side_is_negative(dx, dy, 1.0 / 6.0, 5.0 / 6.0, 0.5, 1.0) {
            row += 1;
        }
    } else if dx > dy {
        // faithful: 0x237a8 in libdggal.so e75d6ab1, replacing RI7H.ec:1626.
        if tie_break_side_is_negative(dx, dy, 5.0 / 6.0, 1.0 / 6.0, 0.5, 0.0) {
            row -= 1;
        }
        // faithful: 0x23d50 in libdggal.so e75d6ab1, replacing RI7H.ec:1632.
    } else if tie_break_side_is_negative(dx, dy, 0.0, 0.5, 1.0 / 6.0, 5.0 / 6.0) {
        col -= 1;
    }

    if north && col == p && row == 0 {
        return Some(Z::new(l49r, 0xA, 0, 0, 0));
    }
    if south && col == 0 && row == p {
        return Some(Z::new(l49r, 0xB, 0, 0, 0));
    }

    if col == p && row < p && !south_rhombus {
        col = p - row;
        row = 0;
        root += 2;
    } else if row == p && col < p && south_rhombus {
        row = p - col;
        col = 0;
        root += 2;
    } else if row < 0 && col < 0 {
        row += p;
        col += p;
        root -= 2;
    } else if row < 0 {
        row += p;
        root -= 1;
    } else if col < 0 {
        col += p;
        root -= 1;
    } else if col >= p && row >= p {
        row -= p;
        col -= p;
        root += 2;
    } else if row >= p {
        row -= p;
        root += 1;
    } else if col >= p {
        col -= p;
        root += 1;
    }
    if root < 0 {
        root += 10;
    } else if root > 9 {
        root -= 10;
    }

    if row < 0 || row >= p || col < 0 || col >= p {
        // faithful: the eC's `fromCentroid` returns nullZone here (RI7H.ec:1680),
        // and so does this, as the absent value. py4dggs returns a corner cell of
        // the rhombus instead, a difference it keeps on purpose as inert. The arm is
        // unreachable in every census taken so far: no hits over two hundred
        // thousand random points at every resolution in py4dggs, nor over the
        // 1,500,120 quantisations this port has been measured across. Should it ever
        // be reached, the engine's answer is the one given.
        return None;
    }

    Some(Z::new(l49r, root, row, col, 0))
}

// --------------------------------------------------------------------------- //
// I7H <-> Z7 conversion (RI7H_Z7.ec)
// --------------------------------------------------------------------------- //

/// The index of `zone` among `parent`'s primary children, counting from zero. An
/// odd level reads it straight off the sub-hexagon; an even level scans the
/// children list.
///
/// faithful: on a miss the eC returns the child count, since its loop index
/// survives the loop and `return i;` follows it (RI7H_Z7.ec:143-146), and so does
/// this. py4dggs returns zero instead, and does not record the difference. The miss
/// is real: deep under the polar pentagon of base cell 0, the geometric parent search
/// that [`from_7h`] makes through [`compute_parents`] can name a parent whose primary
/// children do not include the zone. There the child count, carried through the
/// rotation offset and the digit map, gives the engine's identifier, where zero gave
/// another zone's: the quantisation of the engine's own centroid of twelve zones per
/// grid, at resolutions 17 and 19, depended on it.
fn get_child_position(parent: &Z, zone: &Z) -> i64 {
    let level = zone_level(zone);
    if level & 1 != 0 {
        return zone.sub_hex - 1;
    }
    let children = get_primary_children(parent);
    for (i, ch) in children.iter().enumerate() {
        if zone_eq(Some(*ch), Some(*zone)) {
            return i as i64;
        }
    }
    children.len() as i64
}

/// Re-indexes a child position around a pentagon base cell so that the Z7 digit
/// ordering stays consistent (`adjustZ7PentagonChildPosition`, `RI7H_Z7.ec:150-167`).
fn adjust_z7_pentagon_child_position(i: i64, level: i64, p_root: i64) -> i64 {
    if i == 0 {
        return 0;
    }
    let mut i = i;
    let south_p_rhombus = p_root & 1 != 0;
    let odd_level = level & 1 != 0;
    // `i` is a non-zero child position, hence positive, so every left operand below
    // is non-negative and a plain `%` matches Python's own result.
    if p_root == 10 {
        i = (i + 1) % 5 + 1;
    } else if p_root == 11 {
        i = (i + if odd_level { 3 } else { 4 }) % 5 + 1;
    } else if !odd_level && !south_p_rhombus {
        i = (i + 5) % 5 + 1;
    }
    if south_p_rhombus && i >= 3 {
        i += 1;
    }
    i
}

/// The inverse of [`adjust_z7_pentagon_child_position`]
/// (`deadjustZ7PentagonChildPosition`, `RI7H_Z7.ec:169-187`), used when reconstructing the
/// I7H child position from a Z7 digit.
fn deadjust_z7_pentagon_child_position(i: i64, level: i64, p_root: i64) -> i64 {
    if i == 0 {
        return 0;
    }
    let mut i = i;
    let south_p_rhombus = p_root & 1 != 0;
    let odd_level = level & 1 != 0;
    if south_p_rhombus && i >= 4 {
        i -= 1;
    }
    // `i` is at least one on every arm, so `i - 1` plus a positive addend is
    // non-negative and a plain `%` matches Python's own result.
    if p_root == 10 {
        i = (i - 1 + 3) % 5 + 1;
    } else if p_root == 11 {
        i = (i - 1 + if odd_level { 1 } else { 5 }) % 5 + 1;
    } else if !odd_level && !south_p_rhombus {
        i = (i - 1 + 4) % 5 + 1;
    }
    i
}

/// The per-level rotation correction accumulated while walking from parent to
/// child, so that the Z7 digit ring stays aligned through pentagons and edge
/// hexagons (`getLevelRotationOffset`, RI7H_Z7.ec). The cascade of cases mirrors
/// the eC exactly, including its child-position sentinel: the eC opens with
/// `if(i == -1) i = getChildPosition(parent, zone);` (RI7H_Z7.ec:200-201), so
/// `prev_i` is an input. A caller that already knows the position passes it, and
/// only `-1` asks this function to look it up. The supplied value is a cache of
/// exactly the same lookup, so the branch preserves the structure rather than
/// changing the result.
fn get_level_rotation_offset(
    l: i64,
    prev_i: i64,
    zone: &Z,
    parent: Option<&Z>,
    grand_parent: Option<&Z>,
) -> i64 {
    let Some(parent) = parent else { return 0 };
    let mut i = prev_i;
    if i == -1 {
        i = get_child_position(parent, zone);
    }
    if i == 0 {
        return 0;
    }

    let mut offset = 0;
    let p_root = parent.root;
    let pn_points = zone_npoints(parent);
    let odd_level = l & 1 != 0;
    let south_p_rhombus = p_root & 1 != 0;
    let is_edge_hex = !odd_level && zone_is_edge_hex(zone);
    let p_edge_hex = odd_level && zone_is_edge_hex(parent);
    let gp_edge_hex = !odd_level && grand_parent.is_some_and(zone_is_edge_hex);

    if pn_points == 5 {
        i = adjust_z7_pentagon_child_position(i, l, p_root);
    }

    if p_root >= 10 {
        if pn_points == 5 {
            offset += i + if odd_level {
                if south_p_rhombus { 0 } else { 3 }
            } else if south_p_rhombus {
                5
            } else {
                2
            };
        } else if is_edge_hex
            && (!south_p_rhombus || !zone_eq(Some(*zone), get_centroid_child(parent)))
        {
            offset += 5;
        }
    }

    if south_p_rhombus && is_edge_hex {
        offset += 1;
    }

    if p_edge_hex {
        if !south_p_rhombus && i >= 4 {
            offset += 1;
        } else if south_p_rhombus && (i == 0 || (3..=5).contains(&i)) {
            offset += 5;
        }
    } else if gp_edge_hex {
        let grand_parent = grand_parent.expect("gp_edge_hex implies a grandparent");
        let pc = get_primary_children(grand_parent);
        let c = get_primary_children(parent);
        if south_p_rhombus {
            if pc.len() > 1
                && zone_eq(Some(pc[1]), Some(*parent))
                && c.len() > 2
                && c.len() > 5
                && c[2].root != c[5].root
                && (i == 4 || i == 5)
            {
                offset += 5;
            }
        } else if pc.len() > 4 && zone_eq(Some(pc[4]), Some(*parent)) && (i == 1 || i == 2) {
            offset += 5;
        }

        if zone_eq(Some(*parent), get_centroid_child(grand_parent)) {
            if south_p_rhombus {
                if i > 2 {
                    offset += 5;
                }
            } else if i == 5 || i == 6 {
                offset += 1;
            }
        }

        if south_p_rhombus && is_edge_hex && (i == 4 || i == 5) {
            offset += 5;
        }
    }
    offset
}

/// The Z7 address of an I7H zone: the base pentagon plus a twenty-level packed
/// ancestry of three-bit digits. Replaces py4dggs's dictionary return; the `level`
/// key it also carries is not read by any caller here, so it is not reproduced.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Z7Parts {
    root_pentagon: u64,
    ancestry: u64,
}

/// Converts an I7H zone to its Z7 address (`Z7Zone::from7H`, RI7H_Z7.ec). It walks
/// the parent chain from the top down, applying the rotation offset and the
/// pentagon re-indexing at each level, and then maps every I7H child position
/// through `C_MAP`. Returns `None` above the representable level.
///
/// faithful: the eC rejects `level > 19` and returns nullZone (RI7H_Z7.ec:384),
/// whereas this rejects only `level > 20`, so a level-20 zone converts here where
/// DGGAL declines to. py4dggs does the same on purpose, on the ground that Z7 level
/// 20 packs but has no geometry: the packed identifier stays bit-identical to
/// DGGAL's, and [`HexA7::is_null_geometry`] is the contract that keeps a level-20
/// address from being given real geometry. See [`zone_from_steps`] for the two
/// matching guards on the inverse conversion.
fn from_7h(zone: &Z) -> Option<Z7Parts> {
    let level = zone_level(zone);
    if level > 20 {
        return None;
    }

    let parents = compute_parents(zone);
    if parents.len() < level as usize {
        // faithful: with a stated limit. The eC's `computeParents` has no early
        // exit: its `while(l > 0)` loop (RI7H_Z7.ec:289-294) always writes exactly
        // `level` entries, storing nullZone wherever `parent0` finds no parent and
        // then continuing to take `parent0` OF that nullZone. `from7H` goes on to
        // read those null entries, and DGGAL duly answers with a deterministic but
        // meaningless identifier: all five sample points that reach this state at
        // resolution 18 give the same `01222222222222222220`.
        //
        // py4dggs instead breaks out of the walk (hex_a7.py:821-822), leaving a
        // short list that its own `from_7h` then indexes out of range. This port
        // reproduces neither: matching the eC would mean emulating the bit layout
        // of the `I7HZone` packing and of every predicate's behaviour on the null
        // pattern, which is a piece of work in its own right and yields garbage in
        // any case, while following py4dggs means a panic. The absent value is
        // returned instead, and `quantize` turns it into its documented nullZone
        // placeholder.
        //
        // Reachable, but barely: five of 1,500,000 sample quantisations, all at
        // resolution 18 and all within a quarter of a degree of the north pole.
        // A later census at the bottom-right corner of the 5x6 layout reaches it
        // away from any pole as well, at `(5.0, 5.000001)` and the eighteen doubles
        // just below that abscissa, at resolutions 18 and 19: the corner wrap
        // carries the point to `(0, 1.2e-6)`, `from_centroid` gives the zone
        // `l49r 9, root 0, row 40, col 0`, and its chain is 0 and 1 long against a
        // level of 18 and 19. The engine answers `01222222222222222220` there, the
        // same meaningless identifier, whose centroid lies some twenty million cell
        // widths from the point it was asked about.
        return None;
    }
    let mut ancestry: u64 = 0;
    let mut offset: i64 = 0;
    let mut prev_i: i64 = 0;

    for l in 1..=level {
        let p_index = (level - l) as usize;
        let z = if l == level {
            *zone
        } else {
            parents[p_index - 1]
        };
        // The guard above has established that the chain holds `level` entries, so
        // `p_index`, and `p_index + 1` and `+ 2` under their own conditions, are all
        // in range.
        let parent = parents[p_index];
        let grand_parent = if l > 1 {
            Some(parents[p_index + 1])
        } else {
            None
        };
        let great_grand_parent = if l > 2 {
            Some(parents[p_index + 2])
        } else {
            None
        };

        let mut i = get_child_position(&parent, &z);
        // The offset and the rotation contribution are both non-negative, so a
        // plain `%` matches Python's own result.
        offset = (offset
            + get_level_rotation_offset(
                l - 1,
                prev_i,
                &parent,
                grand_parent.as_ref(),
                great_grand_parent.as_ref(),
            ))
            % 6;
        prev_i = i;

        if i != 0 {
            if zone_npoints(&parent) == 5 {
                i = adjust_z7_pentagon_child_position(i, l, parent.root);
            }
            // `i >= 1` and `offset >= 0`, so the left operand is non-negative.
            i = (i - 1 + offset) % 6 + 1;
        }

        let shift = (19 - (l - 1)) * 3;
        ancestry |= (C_MAP[i as usize] as u64) << shift;
    }

    for l in level + 1..21 {
        let shift = (19 - (l - 1)) * 3;
        ancestry |= 7u64 << shift;
    }

    let top_parent = if level == 0 {
        *zone
    } else {
        *parents.last().expect("a non-zero level has a parent chain")
    };
    // Bounded, where the eC indexes `rootMap[... .rootRhombus]` unchecked
    // (RI7H_Z7.ec:412). For a root of 0 to 11 this is the same lookup, and for any other
    // the unbounded one panicked, so no answer given before changes. A point far outside
    // the 5x6 layout, which `HexA7::quantize` can be handed directly, drives the root to
    // any value at all; such a point has no cell, and the engine offers no answer to
    // follow for it, since it keeps the root in a four-bit field (RI7H.ec:869) and so
    // reads roots of 12 to 15 past the end of the twelve-entry `rootMap`
    // (RI7H_Z7.ec:10), which is undefined. The null zone is answered instead of an index
    // out of bounds.
    let root_pentagon = *ROOT_MAP.get(usize::try_from(top_parent.root).ok()?)?;

    Some(Z7Parts {
        root_pentagon,
        ancestry,
    })
}

/// Reconstructs the I7H zone from a Z7 base cell and its direction digits
/// (`Z7Zone::to7H`, RI7H_Z7.ec), the inverse of [`from_7h`]. It walks the digits
/// from the top down, undoing the rotation offset and the pentagon re-indexing, and
/// steps even to odd (a sub-hexagon) and odd to even (a primary child) alternately.
///
/// faithful: the eC guards this conversion twice against a level-20 zone and this
/// port carries neither guard. It opens with `if((this & 7) != 7) return nullZone;`
/// ("I7HZone are only valid up to level 19", RI7H_Z7.ec:310-311), which rejects any
/// identifier whose twentieth digit slot is not the terminator, and it breaks with
/// `zone = nullZone` on `else if(level == 19)` ("7H does not support level 20
/// zones", RI7H_Z7.ec:348-352) when a twentieth digit would step odd to even.
/// py4dggs drops both, deliberately, since a level-20 zone is kept from any
/// geometry elsewhere, and this port follows it. The consequence is that a
/// twenty-digit address reconstructs to a real I7H zone and
/// `HexA7::quantize(p, 20)` yields a real address where DGGAL yields nullZone.
/// [`HexA7::is_null_geometry`] is the contract that rescues it: `Grid` consults it
/// before asking for any geometry, and it refuses every address of twenty digits or
/// more.
///
/// The parent and the grandparent handed to the rotation offset are the walk's own two
/// previous zones, as `to7H` reads them back from the array it fills on the way down,
/// and are not found again through [`get_parent0`]. Each zone of the walk is a primary
/// child of the one before it, and the geometric search that [`get_parent0`] makes at
/// an even level almost always names that same zone. Deep under the pentagon of base
/// cell 0, the north pole of the 5x6 layout, it need not: there the primary children of
/// neighbouring zones overlap, and the search, which takes the first odd zone that lists
/// the child, can name the polar zone's own centre, or no zone at all, where the walk
/// came from a neighbouring sub-hexagon. At such zones, `0000000000000000131` among
/// them, only the walk's own zone gives the engine's centroid and vertices.
pub(super) fn zone_from_steps(base_cell: u64, directions: &[u8]) -> Z {
    let z7_root = INV_ROOT_MAP[base_cell as usize];
    let mut zone = Z::new(0, z7_root, 0, 0, 0);
    let mut offset: i64 = 0;
    let mut prev_cix: i64 = 0;
    // faithful: the rotation offset is given the walk's own two previous zones, which
    // the eC stores at every level, `parents[pStart] = zone` (RI7H_Z7.ec:321), and
    // reads back as `parents[pStart + 1]` and `parents[pStart + 2]`
    // (RI7H_Z7.ec:332-333). py4dggs recomputes both through `get_parent0`, whose
    // geometric search can name another zone deep under the north polar pentagon.
    let mut parent: Option<Z> = None;
    let mut grand_parent: Option<Z> = None;

    for l in 0..directions.len() as i64 {
        let b = directions[l as usize];
        if b == 7 {
            break;
        }
        let mut cix = INV_C_MAP[b as usize];
        let n_points = zone_npoints(&zone);

        if cix != 0 || l < 19 {
            // Both terms are non-negative, so a plain `%` matches Python's result.
            offset = (offset
                + get_level_rotation_offset(
                    l,
                    prev_cix,
                    &zone,
                    parent.as_ref(),
                    grand_parent.as_ref(),
                ))
                % 6;
        }
        if cix != 0 {
            cix = cix - 1 - offset;
            if cix < 0 {
                cix += 6;
            }
            cix += 1;
            if n_points == 5 {
                cix = deadjust_z7_pentagon_child_position(cix, l + 1, zone.root);
            }
        }
        prev_cix = cix;

        let child = if l & 1 == 0 {
            Z::new(zone.l49r, zone.root, zone.row, zone.col, 1 + cix)
        } else {
            let children = get_primary_children(&zone);
            if (cix as usize) < children.len() {
                children[cix as usize]
            } else {
                break;
            }
        };
        grand_parent = parent;
        parent = Some(zone);
        zone = child;
    }
    zone
}

// --------------------------------------------------------------------------- //
// Public planar geometry: centroid and vertices in 5x6 space (RI7H.ec)
// --------------------------------------------------------------------------- //

/// The zone centroid in 5x6 space, folded into the canonical domain: the eC's `centroid`
/// property (`RI7H.ec:2895-2976`), then the fold at the head of `getVertices`
/// (`RI7H.ec:1806-1809`).
fn centroid5x6(zone: &Z) -> (f64, f64) {
    let (cx, cy) = zone_centroid(zone);
    fold_centroid(cx, cy)
}

/// The fold `getVertices` applies to the centroid it is handed (`RI7H.ec:1806-1809`).
pub(super) fn fold_centroid(mut cx: f64, mut cy: f64) -> (f64, f64) {
    if cy > 6.0 + 1e-9 || cx > 5.0 + 1e-9 {
        cx -= 5.0;
        cy -= 5.0;
    } else if cx < 0.0 {
        cx += 5.0;
        cy += 5.0;
    }
    (cx, cy)
}

/// The zone's five or six vertices in 5x6 space around the centroid `(cx, cy)`
/// (`I7HZone::getVertices`, `RI7H.ec:1799-1895`, ported as [`get_vertices`]). A polar
/// base cell (root 0xA or 0xB) fans five vertices along the pole diagonal, or traces the
/// odd-level hexagon; every other cell traces its boundary through
/// [`add_non_polar_base_vertices`]. A final seam wrap brings each vertex into the
/// canonical domain.
fn verts5x6(zone: &Z, cx: f64, cy: f64, level: i64) -> Vec<(f64, f64)> {
    let mut out = get_vertices(zone.root, zone.sub_hex, zone_npoints(zone), cx, cy, level);
    for (vx, vy) in out.iter_mut() {
        if *vx > 5.0 && *vy > 5.0 {
            *vx -= 5.0;
            *vy -= 5.0;
        } else if *vx < 0.0 && *vy < 1.0 {
            *vx += 5.0;
            *vy += 5.0;
        }
    }
    out
}

/// The vertices of a zone of level `level`, rhombus `root` and sub-hex `sub_hex` about the
/// folded centroid `(cx, cy)`, as `I7HZone::getVertices` (`RI7H.ec:1799-1895`, compiled at
/// `0x1eff0`) returns them: without the final wrap of [`verts5x6`], and with `n_pts` the
/// caller's, as `getNeighbors` passes its own zone's count for its centroid child
/// (`RI7H.ec:1146`).
pub(super) fn get_vertices(
    root: i64,
    sub_hex: i64,
    n_pts: i64,
    cx: f64,
    cy: f64,
    level: i64,
) -> Vec<(f64, f64)> {
    // `level` is non-negative, so the integer division matches Python's
    // `math.floor(level / 2)`.
    let p = pow7(level / 2);
    let oonp = 1.0 / (7 * p) as f64;
    let is_odd = level & 1 != 0;
    let mut verts: Option<Vec<(f64, f64)>> = None;

    if root == 0xA || root == 0xB {
        let is_pole = sub_hex <= 1;
        // The two polar arms differ only in the sign of the fan and in the corner
        // they start from; they are kept apart as the Python keeps them. The two
        // scaled offsets stay optional, as Python's `bx = by = None` does, so that
        // a fan built without its corner having been set fails loudly rather than
        // quietly emitting five points at the origin.
        //
        // faithful: the eC forms the corner once, `b = { 1 - oonp*A, 0 + oonp*A }`
        // for the north and `b = { 4 + oonp*A, 6 - oonp*A }` for the south, and
        // then writes the fan as `b.x + k` or `b.x - k` (RI7H.ec:1816-1833 and
        // :1858-1879). The compiler moved the whole number inside the bracket, so
        // that each vertex is the integer corner of its own rhombus plus or minus
        // the offset, rather than the first corner shifted. In `libdggal.so`,
        // BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, the even north arm
        // at `0x1f2b0` reads:
        //
        // ```text
        //   1f2b4  movsd  xmm3,[rip+0x2a8bc]  ; 49b78 = 2.0
        //   1f2bc  addsd  xmm2,xmm1           ; 1 + oonp*A, the second y
        //   1f2c0  subsd  xmm0,xmm1           ; 1 - oonp*A, the first x
        //   1f2cc  movapd xmm0,xmm3           ; reload the corner, here 2.0
        //   1f2d0  addsd  xmm3,xmm1           ; 2 + oonp*A
        //   1f2d4  subsd  xmm0,xmm1           ; 2 - oonp*A, the second x
        //   ... and so on through 3.0, 4.0 and 5.0 at 49b70, 49a80 and 49b50
        // ```
        //
        // so the k-th vertex is `(1 + k) - oonp*A`, not `(1 - oonp*A) + k`. The
        // odd north arm at `0x1f480`, the even south arm at `0x1f0bf` and the odd
        // south arm follow the same shape, subtracting `k` from the corner. The
        // fraction is added last in every case, which is why the offsets are kept
        // here rather than the corner.
        let (mut fan_x, mut fan_y): (Option<f64>, Option<f64>) = (None, None);
        if root == 0xA {
            if !is_odd || is_pole {
                let a = if is_odd { ODD_HEX_C } else { EVEN_HEX_A };
                fan_x = Some(oonp * a);
                fan_y = Some(oonp * (if is_odd { ODD_HEX_A } else { EVEN_HEX_A }));
            } else {
                let v = scale_offsets(&ODD_HEX_VERTS, oonp);
                verts = Some(add_non_polar_base_vertices(cx, cy, n_pts, &v));
            }
            if verts.is_none() {
                let (fx, fy) = (fan_x.expect(POLAR_FAN), fan_y.expect(POLAR_FAN));
                verts = Some(
                    (0..5)
                        .map(|i| {
                            let k = i as f64;
                            ((1.0 + k) - fx, k + fy)
                        })
                        .collect(),
                );
            }
        } else {
            if !is_odd || is_pole {
                let a = if is_odd { ODD_HEX_C } else { EVEN_HEX_A };
                fan_x = Some(oonp * a);
                fan_y = Some(oonp * (if is_odd { ODD_HEX_A } else { EVEN_HEX_A }));
            } else {
                let v = scale_offsets(&ODD_HEX_VERTS, oonp);
                verts = Some(add_non_polar_base_vertices(cx, cy, n_pts, &v));
            }
            if verts.is_none() {
                let (fx, fy) = (fan_x.expect(POLAR_FAN), fan_y.expect(POLAR_FAN));
                verts = Some(
                    (0..5)
                        .map(|i| {
                            let k = i as f64;
                            ((4.0 - k) + fx, (6.0 - k) - fy)
                        })
                        .collect(),
                );
            }
        }
    } else {
        let hex_verts = if is_odd {
            &ODD_HEX_VERTS
        } else {
            &EVEN_HEX_VERTS
        };
        let v = scale_offsets(hex_verts, oonp);
        verts = Some(add_non_polar_base_vertices(cx, cy, n_pts, &v));
    }

    verts.expect("every branch above assigns the vertex list")
}

/// The message behind the two `expect` calls on the polar fan in [`get_vertices`]: the
/// corner is set on exactly the branch that does not trace a boundary instead, so
/// reaching a fan without one would be a logic error in that pairing.
const POLAR_FAN: &str =
    "get_vertices: the polar fan needs its corner, which only the non-tracing branch sets";

/// Scales a vertex offset table by `oonp`, the Python list comprehension that
/// precedes each call to [`add_non_polar_base_vertices`]. The eC writes the scaled
/// offsets out in `getVertices` (`RI7H.ec:1839-1844`, `RI7H.ec:1884-1889`).
fn scale_offsets(table: &[(f64, f64); 6], oonp: f64) -> [(f64, f64); 6] {
    let mut v = [(0.0, 0.0); 6];
    for i in 0..6 {
        v[i] = (table[i].0 * oonp, table[i].1 * oonp);
    }
    v
}

// --------------------------------------------------------------------------- //
// Neighbours, as the engine finds them (RI7H.ec getNeighbors)
// --------------------------------------------------------------------------- //

/// The zone's neighbours as `I7HZone::getNeighbors` lists them (`RI7H.ec:1121-1218`),
/// in its order, one per vertex of the zone's centroid child, `None` where
/// `fromCentroid` answers nullZone. A neighbour is found by stepping from the zone's
/// centroid three times the offset of each of its centroid child's vertices, and
/// quantising the point reached; the level-20 centroid child of a level-19 zone, which
/// `I7HZone` cannot represent, is never built, only its vertices.
///
/// `getNeighbors` carries no `__attribute__ ((optimize("-fno-unsafe-math-optimizations")))`,
/// and it is compiled at `0x241d0` to `0x245ff` in `libdggal.so`, BuildID
/// `e75d6ab18f8b460713ebcd21df6f07922fb909c3`. Every function it calls is called out of
/// line there, at the address this crate's port of it follows: the `centroid` property
/// (`0x20200`), `getOddLevelCentroidChildRootRowCol` (`0x221c0`), `getVertices`
/// (`0x1eff0`), `cross5x6Interruption` (`0x360d0`), `move5x6Vertex2` (`0x37be0`),
/// `canonicalize5x6` (`0x3a080`) and `fromCentroid` (`0x234e0`).
fn get_neighbors(z: &Z) -> Vec<Option<Z>> {
    let n_level = zone_level(z);
    if n_level > 19 {
        return Vec::new();
    }
    // The centroid property, unfolded, as `c = centroid` (RI7H.ec:1132); `getVertices`
    // folds its own copy (RI7H.ec:1806-1809).
    let c = zone_centroid(z);
    let (root, c_l49r, c_sh) = if n_level & 1 != 0 {
        // The eC passes the -1 of an out-of-range centroid child on as the root
        // (RI7H.ec:1136), and `getVertices`, given an even child, compares the root only
        // with 0xA and 0xB to choose a polar fan, so it traces a boundary all the same.
        let root = get_odd_level_centroid_child_root_row_col(z).map_or(-1, |cc| cc.root);
        (root, z.l49r + 1, 0)
    } else {
        (z.root, z.l49r, 1)
    };
    let (fx, fy) = fold_centroid(c.0, c.1);
    let c_verts = get_vertices(
        root,
        c_sh,
        zone_npoints(z),
        fx,
        fy,
        2 * c_l49r + i64::from(c_sh > 0),
    );

    // `Min(4, (int)(c.x + 1E-11))` and `Min(5, (int)(c.y + 1E-11))` (RI7H.ec:1152), each a
    // 32-bit `cvttsd2si` at `0x242e0` and `0x242f7`.
    let cx = f64::from(cvtt_i32(c.0 + 1e-11).min(4));
    let cy = f64::from(cvtt_i32(c.1 + 1e-11).min(5));
    // faithful: `c.x - c.y - 1E-11 > 0` (RI7H.ec:1153) is compiled as `c.x - c.y > 1e-11`,
    // at `0x24416` and again at `0x24468` in `libdggal.so`, BuildID
    // `e75d6ab18f8b460713ebcd21df6f07922fb909c3`. The two are the same test in IEEE
    // arithmetic, since a difference of two doubles is zero only when they are equal, so
    // either form would do; the compiled one is written.
    let north = c.0 - c.1 > 1e-11;

    let mut out = Vec::with_capacity(c_verts.len());
    for &(vx, vy) in &c_verts {
        let mut cc = c;
        let mut dx = vx - cc.0;
        let mut dy = vy - cc.1;
        // faithful: `dx = cVerts[i].x - 5 - cc.x` and its ordinate (RI7H.ec:1164-1165)
        // are compiled as `(cVerts[i].x - cc.x) - 5`, at `0x243cc` and `0x243d4`
        // (`subsd` of 5.0 at `0x49b50` from the difference already taken), and
        // `cVerts[i].x + 5 - cc.x` and its ordinate (RI7H.ec:1169-1170) as
        // `(cVerts[i].x - cc.x) + 5`, at `0x24574` and `0x2457c`, in `libdggal.so`, BuildID
        // `e75d6ab18f8b460713ebcd21df6f07922fb909c3`. The rounding falls at the magnitude
        // of five rather than of the offset.
        if dx > 3.0 && dy > 3.0 {
            dx -= 5.0;
            dy -= 5.0;
        } else if dx < -3.0 && dy < -3.0 {
            dx += 5.0;
            dy += 5.0;
        }

        if dx.abs() < 1.0
            && dy.abs() < 1.0
            && ((north && (c.1 - cy).abs() < 1e-11) || (!north && (c.0 - cx).abs() < 1e-11))
        {
            let ci = crate::fivebysix::cross5x6_interruption(c.0, c.1, !north, true);
            let mut x = vx - ci.0;
            let mut y = vy - ci.1;
            // faithful: the eC's first arm, `if(x > 3 && dy > 3)` (RI7H.ec:1186), tests
            // `dy`, which is below one in magnitude on this path, so it can never be taken;
            // the compiler removed it, the first comparison at `0x244aa` being the second
            // arm's against -3.0, and it is omitted here too. The second arm,
            // `cVerts[i].x + 5 - ci.x` and its ordinate (RI7H.ec:1193-1194), is compiled as
            // `(cVerts[i].x - ci.x) + 5`, at `0x244ed` and `0x244f1` in `libdggal.so`,
            // BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, as the pair above.
            if x < -3.0 && y < -3.0 {
                x += 5.0;
                y += 5.0;
            }
            if x.abs() < dx.abs() && y.abs() < dy.abs() {
                cc = ci;
                dx = x;
                dy = y;
            }
        }

        let v = move5x6_vertex2(cc.0, cc.1, dx * 3.0, dy * 3.0, false);
        let v = crate::fivebysix::canonicalize5x6(v.0, v.1);
        out.push(from_centroid(n_level, v.0, v.1));
    }
    out
}

// --------------------------------------------------------------------------- //
// Parents and children, as the engine finds them (RI7H.ec getParents, getChildren)
// --------------------------------------------------------------------------- //

/// The factor by which `getParents` shortens the step from a centroid towards a vertex, the
/// `.99` of `RI7H.ec:1338`: the double that the library holds at `0x4a120`.
const PARENT_STEP_FACTOR: f64 = 0.99;

/// The factor by which `getChildren` lengthens the step from the centroid of the centroid
/// child towards one of its vertices, the `3` of `RI7H.ec:2706`: the double that the library
/// holds at `0x49b70`.
const CHILD_STEP_FACTOR: f64 = 3.0;

/// The zone's parents as `I7HZone::getParents` lists them (`RI7H.ec:1261-1358`), in its
/// order and before anything is left out: none where `parent0` answers nullZone, which is at
/// level 0 and wherever the search for the parent of an even zone finds none; the primary
/// parent alone for a zone that is its centroid child; and otherwise the primary parent and
/// the first other zone of its level found by stepping from the zone's centroid towards
/// each of the zone's own vertices in turn, a hundredth short of it, and quantising the
/// point reached. Where no vertex leads to another zone the primary parent stands alone, as
/// in the eC, whose only word on it is a message of its debug build.
///
/// `getParents` carries no `__attribute__ ((optimize("-fno-unsafe-math-optimizations")))`,
/// and it is compiled at `0x24610` to `0x24a10` in `libdggal.so`, BuildID
/// `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, in the shape of `getNeighbors`. Every function
/// it calls is called out of line there: the `parent0` property (`0x201b0`), the
/// `centroidChild` property (`0x227b0`), the `centroid` property (`0x20200`), `getVertices`
/// (`0x1eff0`), `cross5x6Interruption` (`0x360d0`), `move5x6Vertex2` (`0x37be0`) and
/// `fromCentroid` (`0x234e0`). The five sites at which the compiled order departs from the
/// source, GP7-1 to GP7-5, are followed below. They are inert: GP7-4 is the same test either
/// way, and with the others put back to the source's order no list changed at any of the
/// 62,876 zones compared.
fn get_parents(z: &Z) -> Vec<Z> {
    let Some(parent0) = get_parent0(z) else {
        return Vec::new();
    };
    let mut out = vec![parent0];
    if get_centroid_child(&parent0) == Some(*z) {
        return out;
    }
    // The centroid property, unfolded, as `c = centroid` (RI7H.ec:1275); `getVertices`
    // folds its own copy (RI7H.ec:1806-1809).
    let c = zone_centroid(z);
    let (fx, fy) = fold_centroid(c.0, c.1);
    let vertices = get_vertices(z.root, z.sub_hex, zone_npoints(z), fx, fy, zone_level(z));
    let p_level = zone_level(&parent0);

    for &(vx, vy) in &vertices {
        let mut acc = c;
        let mut dx = vx - acc.0;
        let mut dy = vy - acc.1;
        // faithful (GP7-1, GP7-2): `dx = vertices[i].x - 5 - acc.x` and its ordinate
        // (RI7H.ec:1292-1293) are compiled as `(vertices[i].x - acc.x) - 5`, at `0x24950`
        // and `0x24958` (`subsd` of 5.0 at `0x49b50` from the differences taken at
        // `0x2477e` and `0x2478a`), and `vertices[i].x + 5 - acc.x` and its ordinate
        // (RI7H.ec:1297-1298) as `(vertices[i].x - acc.x) + 5`, at `0x24968` and `0x24970`,
        // in `libdggal.so`, BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`. The rounding
        // falls at the magnitude of five rather than of the offset. Either coordinate
        // suffices here (RI7H.ec:1290, :1295), where `getNeighbors` and `getChildren` ask
        // for both.
        if dx > 3.0 || dy > 3.0 {
            dx -= 5.0;
            dy -= 5.0;
        } else if dx < -3.0 || dy < -3.0 {
            dx += 5.0;
            dy += 5.0;
        }

        if dx.abs() < 1.0 && dy.abs() < 1.0 {
            // faithful (GP7-4): `acc.x - acc.y - 1E-11 > 0` (RI7H.ec:1303) is compiled as
            // `acc.x - acc.y > 1e-11`, the difference at `0x24801` compared at `0x24806` and
            // again at `0x2485a`. The two are the same test in IEEE arithmetic, as at the
            // like site of `get_neighbors`; the compiled one is written.
            let north = acc.0 - acc.1 > 1e-11;
            // `(int)(acc.y + 1E-11)` and `(int)(acc.x + 1E-11)` (RI7H.ec:1304-1305), each a
            // 32-bit `cvttsd2si`, at `0x24821` and `0x2498c`; the library takes only the one
            // its arm compares, and neither is bounded as `getNeighbors` bounds its own.
            let on_interruption = if north {
                (acc.1 - f64::from(cvtt_i32(acc.1 + 1e-11))).abs() < 1e-11
            } else {
                (acc.0 - f64::from(cvtt_i32(acc.0 + 1e-11))).abs() < 1e-11
            };
            if on_interruption {
                let ci = crate::fivebysix::cross5x6_interruption(c.0, c.1, !north, true);
                let mut x = vx - ci.0;
                let mut y = vy - ci.1;
                // faithful (GP7-3, GP7-5): the eC's first arm, `if(x > 3 && dy > 3)`
                // (RI7H.ec:1318), tests `dy`, which is below one in magnitude on this path,
                // so it can never be taken; the compiler removed it, the first comparison
                // after the call, at `0x2489b`, being the second arm's against -3.0, and it
                // is omitted here too. The second arm, `vertices[i].x + 5 - ci.x` and its
                // ordinate (RI7H.ec:1325-1326), is compiled as `(vertices[i].x - ci.x) + 5`,
                // at `0x248db` and `0x248df`, as the pair above.
                if x < -3.0 && y < -3.0 {
                    x += 5.0;
                    y += 5.0;
                }
                if x.abs() < dx.abs() && y.abs() < dy.abs() {
                    acc = ci;
                    dx = x;
                    dy = y;
                }
            }
        }

        // `.99 * dx` and `.99 * dy` (RI7H.ec:1338), each one product with the 0.99 at
        // `0x4a120`, at `0x24710` and `0x24720`.
        let v = move5x6_vertex2(
            acc.0,
            acc.1,
            dx * PARENT_STEP_FACTOR,
            dy * PARENT_STEP_FACTOR,
            false,
        );
        if let Some(q) = from_centroid(p_level, v.0, v.1) {
            if q != parent0 {
                out.push(q);
                return out;
            }
        }
    }
    out
}

/// The zone's children as `I7HZone::getChildren` lists them (`RI7H.ec:2634-2719`), in its
/// order and before anything is left out: the primary children, the centroid child first,
/// and then one zone for each vertex of the centroid child, found by stepping from the
/// centroid child's centroid three times the offset of the vertex and quantising the point
/// reached, `None` where `fromCentroid` answers nullZone. A zone without primary children
/// has none at all: a zone of level 19, and an odd zone for which
/// [`get_primary_children`] finds no centroid child within the rhombi and answers none,
/// where the eC goes on from its null zone (`RI7H.ec:2745-2750`).
///
/// `getChildren` carries no `__attribute__ ((optimize("-fno-unsafe-math-optimizations")))`
/// either, and it is compiled at `0x24a20` to `0x24de0` in `libdggal.so`, BuildID
/// `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, in the shape of `getNeighbors` once more,
/// calling out of line `getPrimaryChildren` (`0x22890`), the `centroid` property
/// (`0x20200`), `getVertices` (`0x1eff0`), `cross5x6Interruption` (`0x360d0`),
/// `move5x6Vertex2` (`0x37be0`), `canonicalize5x6` (`0x3a080`) and `fromCentroid`
/// (`0x234e0`). The five sites at which the compiled order departs from the source, GC7-1
/// to GC7-5, are followed below. They are inert: GC7-4 is the same test either way, and with
/// the others put back to the source's order no list changed at any of the 62,876 zones
/// compared.
fn get_children(z: &Z) -> Vec<Option<Z>> {
    let primary = get_primary_children(z);
    let mut out: Vec<Option<Z>> = primary.iter().copied().map(Some).collect();
    let Some(&c0) = primary.first() else {
        return out;
    };
    // The centroid property of the centroid child, unfolded, as `c = children[0].centroid`
    // (RI7H.ec:2639); `getVertices` folds its own copy (RI7H.ec:1806-1809).
    let c = zone_centroid(&c0);
    let (fx, fy) = fold_centroid(c.0, c.1);
    let c_level = zone_level(&c0);
    let c_verts = get_vertices(c0.root, c0.sub_hex, zone_npoints(&c0), fx, fy, c_level);

    // faithful (GC7-4): `c.x - c.y - 1E-11 > 0` (RI7H.ec:2644) is compiled as
    // `c.x - c.y > 1e-11`, the difference taken once at `0x24b23` and compared at `0x24c67`
    // and again at `0x24cbb`. The two are the same test in IEEE arithmetic, as at the like
    // site of `get_neighbors`; the compiled one is written.
    let north = c.0 - c.1 > 1e-11;
    // `(int)(c.y + 1E-11)` and `(int)(c.x + 1E-11)` (RI7H.ec:2645-2646), each a 32-bit
    // `cvttsd2si`, at `0x24b43` and `0x24b38`, and neither bounded as `getNeighbors` bounds
    // its own.
    let cy = f64::from(cvtt_i32(c.1 + 1e-11));
    let cx = f64::from(cvtt_i32(c.0 + 1e-11));

    for &(vx, vy) in &c_verts {
        let mut cc = c;
        let mut dx = vx - cc.0;
        let mut dy = vy - cc.1;
        // faithful (GC7-1, GC7-2): `dx = cVerts[i].x - 5 - cc.x` and its ordinate
        // (RI7H.ec:2664-2665) are compiled as `(cVerts[i].x - cc.x) - 5`, at `0x24c1e` and
        // `0x24c26` (`subsd` of 5.0 at `0x49b50` from the differences taken at `0x24bf2` and
        // `0x24c04`), and `cVerts[i].x + 5 - cc.x` and its ordinate (RI7H.ec:2669-2670) as
        // `(cVerts[i].x - cc.x) + 5`, at `0x24dac` and `0x24db4`, in `libdggal.so`, BuildID
        // `e75d6ab18f8b460713ebcd21df6f07922fb909c3`. The rounding falls at the magnitude of
        // five rather than of the offset.
        if dx > 3.0 && dy > 3.0 {
            dx -= 5.0;
            dy -= 5.0;
        } else if dx < -3.0 && dy < -3.0 {
            dx += 5.0;
            dy += 5.0;
        }

        if dx.abs() < 1.0
            && dy.abs() < 1.0
            && ((north && (c.1 - cy).abs() < 1e-11) || (!north && (c.0 - cx).abs() < 1e-11))
        {
            let ci = crate::fivebysix::cross5x6_interruption(c.0, c.1, !north, true);
            let mut x = vx - ci.0;
            let mut y = vy - ci.1;
            // faithful (GC7-3, GC7-5): the eC's first arm, `if(x > 3 && dy > 3)`
            // (RI7H.ec:2686), tests `dy`, which is below one in magnitude on this path, so
            // it can never be taken; the compiler removed it, the first comparison after
            // the call, at `0x24cfd`, being the second arm's against -3.0, and it is omitted
            // here too. The second arm, `cVerts[i].x + 5 - ci.x` and its ordinate
            // (RI7H.ec:2693-2694), is compiled as `(cVerts[i].x - ci.x) + 5`, at `0x24d40`
            // and `0x24d44`, as the pair above.
            if x < -3.0 && y < -3.0 {
                x += 5.0;
                y += 5.0;
            }
            if x.abs() < dx.abs() && y.abs() < dy.abs() {
                cc = ci;
                dx = x;
                dy = y;
            }
        }

        // `dx * 3` and `dy * 3` (RI7H.ec:2706), each one product with the 3.0 at `0x49b70`,
        // at `0x24b90` and `0x24ba0`.
        let v = move5x6_vertex2(
            cc.0,
            cc.1,
            dx * CHILD_STEP_FACTOR,
            dy * CHILD_STEP_FACTOR,
            false,
        );
        let v = crate::fivebysix::canonicalize5x6(v.0, v.1);
        out.push(from_centroid(c_level, v.0, v.1));
    }
    out
}

/// The Z7 address of `z7` at resolution `res`: the packed ancestry unpacked into its
/// first `res` digits, stopping at the level-7 terminator, as the `level` property reads
/// the packing (`RI7H_Z7.ec:46-50`). Shared by `quantize`, `neighbors`, `parents` and
/// `children`.
fn z7_address(z7: Z7Parts, res: i64) -> Address {
    let mut address = Address::new(z7.root_pentagon, &[]);
    for l in 0..res {
        let shift = (19 - l) * 3;
        let d = ((z7.ancestry >> shift) & 7) as u8;
        if d == 7 {
            break;
        }
        address.push(d);
    }
    address
}

// --------------------------------------------------------------------------- //
// The zones of a level as a lattice, in the order of listZones (RI7H.ec)
// --------------------------------------------------------------------------- //

/// The value of the zone as an `I7HZone`, as the eC packs its bit class (`RI7H.ec:866-871`):
/// `levelI49R` in four bits from bit 58, `rootRhombus` in four from bit 54, `rhombusIX`, which
/// is `row * 7^levelI49R + col`, in fifty-one from bit 3, and `subHex` in the lowest three.
///
/// It is the key by which the engine orders the zones of one level. `listZones` sorts its
/// answer (`RI7H.ec:829`) by `I7HZone::OnCompare` (`RI7H.ec:873-884`), which compares the
/// levels and then the values, and the Z7 class then converts each zone to its own identifier
/// where it stands (`RI7H_Z7.ec:674-695`): the order of a list is that of these values and
/// not that of its Z7 identifiers. An odd level shares `levelI49R` with the even level below
/// it, so the value orders the zones of one level alone.
///
/// For a zone of the lattice, its row and its column within `7^levelI49R`, the greatest
/// index, `7^18 - 1`, fits the fifty-one bits. A zone that is not of the lattice is packed as
/// the eC's bit class packs any value, each member cut to its own bits, so that the function
/// answers every zone.
fn i7h_value(z: &Z) -> u64 {
    let p = pow7(z.l49r) as u64;
    let index = (z.row as u64).wrapping_mul(p).wrapping_add(z.col as u64);
    ((z.l49r as u64 & 0xf) << 58)
        | ((z.root as u64 & 0xf) << 54)
        | ((index & ((1 << 51) - 1)) << 3)
        | (z.sub_hex as u64 & 7)
}

/// Whether `a` names the child a pentagon does not have: its first digit that is not 0 is 2
/// under the base cells 0 to 5, or 5 under 6 to 11, which the eC's `fromTextID` refuses
/// (`RI7H_Z7.ec:430-438`).
fn names_deleted_child(a: &Address) -> bool {
    let deleted = if a.base <= 5 { 2 } else { 5 };
    a.digits().iter().find(|&&d| d != 0) == Some(&deleted)
}

/// The Z7 address of the I7H zone `z` of `level`, with its key, where `z` is a zone of the
/// lattice; `None` where it is not.
///
/// A cell holds seven values of `subHex` at an odd level (`RI7H.ec:866-871`), and the cell of
/// a pentagon has six zones: the seventh value converts to the address of the child a pentagon
/// does not have, which is refused by its digit, or to the address of one of the six. So the
/// address must also be of the level and read back to `z` through [`zone_from_steps`]: no cell
/// hosts a zone under an address that names another.
///
/// From level 16 that leaves cells that host nothing, in bands beside the two edges of the
/// icosahedron along which the engine's own answers are not consistent, from base cell 0 to 1
/// and from 1 to 6. There the conversion gives a zone no address of its own: either none,
/// where the engine prints an identifier that it cannot read back, or the address of another
/// zone. From level 16 the lattice therefore holds fewer zones than the level counts.
fn lattice_zone(level: u8, z: Z) -> Option<(Address, u64)> {
    let p = pow7(z.l49r);
    if !(0..12).contains(&z.root) || !(0..p).contains(&z.row) || !(0..p).contains(&z.col) {
        return None;
    }
    let a = z7_address(from_7h(&z)?, i64::from(level));
    if names_deleted_child(&a) {
        return None;
    }
    (a.len() == usize::from(level) && zone_from_steps(a.base, a.digits()) == z)
        .then(|| (a, i7h_value(&z)))
}

/// The zones of one cell of the root `root`, a rhombus or a polar root, in the ascending order
/// of their keys: the cell's own zone at an even level, its sub-hexagons B to H at an odd one
/// (`RI7H.ec:866-871`), those that are zones.
fn lattice_cell(level: u8, root: u8, row: u64, col: u64) -> Vec<(Address, u64)> {
    let sub_hexes: &[i64] = if level % 2 == 0 {
        &[0]
    } else {
        &[1, 2, 3, 4, 5, 6, 7]
    };
    sub_hexes
        .iter()
        .filter_map(|&sub_hex| {
            let z = Z::new(
                i64::from(level / 2),
                i64::from(root),
                row as i64,
                col as i64,
                sub_hex,
            );
            lattice_zone(level, z)
        })
        .collect()
}

// --------------------------------------------------------------------------- //
// The compaction of a set of zones (RI7H.ec)
// --------------------------------------------------------------------------- //

/// The compaction of a set of zones: every parent whose children are all in the set stands
/// for those of them that have no other parent, or whose other parent was taken too, one
/// level a pass from the finest level of the set upwards, in the ascending order of the
/// zones' values as `I7HZone`s. Ports `compactI7HZones` (`RI7H.ec:3754-3862`) with the part
/// of `RhombicIcosahedral7H::compactZones` that reads the array into a tree
/// (`RI7H.ec:172-197`), statement for statement; the eC's local names are kept.
///
/// The eC's tree is ordered by the value of the zone, so that a zone given twice counts
/// once, and it leaves the null zone out; a `BTreeMap` from that value, [`i7h_value`], to the
/// zone is that tree, and the caller hands no null zone. A parent is taken when every child
/// that `getChildren` lists for it is in the set, the null entries apart: thirteen for a
/// hexagon and eleven for a pentagon, the primary children and the zones across the vertices
/// of the centroid child. A zone leaves the set when every parent that `getParents` lists for
/// it was taken, which is one for a centroid child and two for any other: so the children of
/// a lone parent give way at its centroid child alone, and the twelve others stay beside the
/// parent and overlap it. The congruent hierarchy of the Z7 digits plays no part.
///
/// Five things in the eC look like oversights and are the engine's answer, so all are kept.
/// `output` is not emptied between passes. `next` is not emptied after a pass that ends the
/// ascent, and the loop goes on after it (the eC's `break` is commented out). A zone for which
/// `getParents` lists no parent leaves the set as soon as a pass meets it, the test that all
/// its parents were taken being true of none: a zone of level 0 beside a finer zone, and, in
/// the two broken seams, a zone for which the search of a parent finds none. And the short
/// cut that follows the loop, meant to put the twelve zones of level 0 in the place of the
/// seventy-two of level 1, builds each of its twelve zones as `{ r, 0 }`, which sets
/// `levelI49R` to `r` and `rootRhombus` to nought: the zone at the origin of the first
/// rhombus at twelve even levels, the last two of which are beyond the grid. No set of one
/// level reaches the short cut, since a pass that holds the seventy-two takes the twelve and
/// drops them all; a set of two levels can.
///
/// The function is meant for zones of one level, which `Grid` checks. The cost is that of
/// the set: as many passes as the finest level of the set counts, at most nineteen, each a
/// bounded amount of work for every zone it holds (its parents twice, the children of each,
/// and a search of a tree for each child); and the zones a pass holds are zones of the set
/// or parents of zones of the pass before, at most two for each.
fn compact(input: &[Z]) -> Vec<Z> {
    use std::collections::BTreeMap;
    let mut max_level = 0;
    let mut zones: BTreeMap<u64, Z> = BTreeMap::new();
    for zone in input {
        max_level = max_level.max(zone_level(zone));
        zones.insert(i7h_value(zone), *zone);
    }

    let mut output: BTreeMap<u64, Z> = BTreeMap::new();
    let mut next: BTreeMap<u64, Z> = BTreeMap::new();
    let mut l = max_level - 1;
    while l >= 0 {
        for zone in zones.values() {
            for c_parent in get_parents(zone) {
                let key = i7h_value(&c_parent);
                if next.contains_key(&key) {
                    continue;
                }
                let parent_all_in = get_children(&c_parent)
                    .iter()
                    .flatten()
                    .all(|c| zones.contains_key(&i7h_value(c)));
                if parent_all_in {
                    next.insert(key, c_parent);
                }
            }
        }

        for (&key, zone) in &zones {
            let all_in = get_parents(zone)
                .iter()
                .all(|c_parent| next.contains_key(&i7h_value(c_parent)));
            if !all_in {
                output.insert(key, *zone);
            }
        }

        if l >= 1 && !next.is_empty() {
            // Not done: the next level becomes the zones to compact.
            zones = std::mem::take(&mut next);
        } else {
            // Done: next is combined with output into the final zones.
            zones = output.clone();
            zones.extend(next.iter().map(|(&key, &zone)| (key, zone)));
        }
        l -= 1;
    }

    if zones.len() >= 72
        && zones
            .values()
            .next()
            .is_some_and(|zone| zone_level(zone) == 1)
    {
        let mut n_l1 = 0;
        for zone in zones.values() {
            match zone_level(zone) {
                1 => n_l1 += 1,
                0 => {}
                _ => break,
            }
        }
        if n_l1 == 72 {
            // Simplifying the full globe to the zones of level 0 (the eC's own note); the
            // zones built are those of the eC's initialiser, as the note above says.
            zones = (0..12)
                .map(|r| Z::new(r, 0, 0, 0, 0))
                .map(|zone| (i7h_value(&zone), zone))
                .collect();
        }
    }
    zones.into_values().collect()
}

// --------------------------------------------------------------------------- //
// The Topology facade
// --------------------------------------------------------------------------- //

impl Topology for HexA7 {
    const APERTURE: u8 = 7;

    /// The twelve icosahedral base cells, and digits within the aperture.
    ///
    /// Both bounds are the ones the tables below are built to: `ROOT_MAP` has
    /// one entry per base cell, and a direction digit indexes the seven-entry
    /// child maps. A `ZoneId` that was never produced by this crate decodes to
    /// a base cell of 12 to 15, which is representable in Z7's four base bits
    /// and has no icosahedral cell behind it, so `Grid` filters such an
    /// address out here before any of the lookups below would run off the end
    /// of a table. They are also the bounds on a base cell and a digit that the
    /// eC's `fromTextID` checks (`RI7H_Z7.ec:428-438`).
    fn is_valid_address(a: &Address) -> bool {
        a.base < ROOT_MAP.len() as u64 && a.digits().iter().all(|&d| d < Self::APERTURE)
    }

    /// Quantises the planar point `p` to a Z7 cell at resolution `res`. As in
    /// DGGAL, `p.face` is not read: the cell follows from `(p.x, p.y)` alone
    /// through `fromCentroid` and then the I7H to Z7 conversion, and the packed
    /// ancestry is unpacked into the first `res` direction digits, stopping at the
    /// level-7 terminator.
    ///
    /// A resolution above the packing's twenty levels gives `None`, as the eC's
    /// `from7H` gives nullZone for any level above 19 (RI7H_Z7.ec:384). Level 20
    /// itself is still quantised, for the reason given at `from_7h`.
    ///
    /// A point with a coordinate that is not finite, or that is a billion or more in
    /// magnitude, is answered `None` before any arithmetic, and so is not a question this
    /// function answers. `p` should lie in or near the layout, as every point the
    /// projection's `forward` returns does: [`crate::Grid`] refuses a coordinate that is
    /// not finite before projecting it, and the three projections' `forward` placed every
    /// finite coordinate tried within the layout, 790,314 pairs up to `f64::MAX` in
    /// magnitude, all within 0 to 5 and 0 to 6. The bound is there for a direct call
    /// through this trait, since the integer row and column arithmetic is sized for the
    /// layout: without it, arithmetic with overflow checks overflowed once a coordinate was
    /// large enough, the threshold falling as the resolution rises, from a magnitude of
    /// about 5.6e18 at resolution 0 to 1.8e11 at resolution 19 and 1.8e10 at the packing's
    /// level 20, over magnitudes a quarter of a decade apart from 1 to `f64::MAX`, and never
    /// below them. The bound lies more than a decade below the lowest.
    ///
    /// Between the layout and the bound a point has no cell either, and the eC's own casts
    /// and table lookups are undefined there; what this function does was sampled, not
    /// proved. It answers `None` at most such points and a meaningless cell at some, and
    /// it does not panic. The one lookup that did panic there, that of the root rhombus in
    /// `from_7h`, is bounded and answers `None`.
    fn quantize(p: PlanarPoint, res: u8) -> Option<Address> {
        // The packing holds twenty levels, the last of them the level with no
        // geometry, so nothing finer can be represented. This also keeps a large
        // `res` away from `from_centroid`, whose arithmetic is sized for the
        // packing's levels.
        if usize::from(res) > NULL_GEOMETRY_LEVEL {
            return None;
        }
        // A point far outside the layout, or not a point at all: see above. The negated
        // form refuses NaN as well.
        if !(p.x.abs() < 1e9 && p.y.abs() < 1e9) {
            return None;
        }

        // Both arms below are DGGAL's nullZone, and the absent value carries it
        // faithfully to the caller.
        //
        // `from_centroid` yields it where the eC returns nullZone (RI7H.ec:1489
        // and :1594), which is reachable: in a narrow strip along the icosahedron
        // edges from base cell 0 to 1 and from 1 to 6, at odd resolutions of 15
        // and above, the candidate search finds no containing cell and DGGAL
        // answers `(null)` too (see `ZoneId::NULL`). It also
        // yields it on the even-level path where the eC does (RI7H.ec:1680), an
        // arm that no census has yet reached. It used to yield it at a few
        // further points as well, where its own candidate search accepted
        // neither of the two cells related by a half-turn about an exact
        // boundary tie while DGGAL's accepted one; those are gone, now that the
        // side test and the geometry both follow the compiled engine.
        //
        // `from_7h` yields it in two cases. Above the representable level
        // (level > 20) it is out of reach from here, since the guard above and
        // `Grid::zone_from_geo`, which rejects a resolution above 19, both come
        // first. On a parent chain shorter than the zone's level it is reachable,
        // and there the eC answers with a meaningless but non-null identifier
        // instead; the comment on `from_7h` explains why that meaningless
        // identifier is not reproduced. It was once thought to be reachable only at
        // a pole, and it is not: a census at the bottom-right corner reaches it at
        // the rhombus-0 corner as well, at `(5.0, 5.000001)`, which the corner wrap
        // carries to `(0, 1.2e-6)`. There `from_centroid` gives the zone
        // `l49r 9, root 0, row 40, col 0` at resolutions 18 and 19, whose parent
        // chain is 0 and 1 long against a level of 18 and 19, and the engine answers
        // `01222222222222222220`, a zone whose centroid lies some twenty million
        // cell widths from the point it was asked about.
        let i7h = from_centroid(i64::from(res), p.x, p.y)?;
        let z7 = from_7h(&i7h)?;

        Some(z7_address(z7, i64::from(res)))
    }

    /// The zone's neighbours as the engine's `getZoneNeighbors` lists them
    /// (`RI7H_Z7.ec:532-538`): the zone taken to its I7H form, as `to7H` takes it; its
    /// neighbours found by `I7HZone::getNeighbors` (`RI7H.ec:1121-1218`), one for each
    /// vertex of its centroid child, in that order, by quantising a point stepped to
    /// across each; and each converted back to Z7, as `from7H` converts it. There are six
    /// for a hexagon and five for a pentagon, fewer where an entry is dropped, in the
    /// engine's order, which is the order in which [`crate::Disk`] walks them.
    ///
    /// A departure from the engine: four kinds of entry that the engine lists are dropped,
    /// as this crate never hands them out. They are its null zone, where the point stepped
    /// to lies in no cell; the zone itself; a repeat of an entry already listed; and an
    /// identifier that the engine cannot read back, which it gives a cell whose parent
    /// chain is too short and for which `from_7h` answers nothing here. All four occur only
    /// in the broken seams, along the edges from base cell 0 to 1 and from 1 to 6 at
    /// resolutions 15 to 19, where the engine's lists are not consistent with themselves;
    /// the list is otherwise the engine's, entry for entry. The oracle suites compare every
    /// list they sample with the engine's less those four kinds, and re-check each entry
    /// dropped as one of them. An entry kept there is the engine's own even where it lies
    /// far from the zone: on IGEO7 the list of `00515151515151510530`, at resolution 18,
    /// holds `00454545454545405216`, `...212` and `...213`, each 35.88 degrees away.
    ///
    /// An address that fails [`HexA7::is_valid_address`], and one of the twentieth level,
    /// which has no geometry (see [`HexA7::is_null_geometry`]), have no neighbours, as
    /// `Grid` answers them and as the engine answers a zone of the twentieth level; so the
    /// method answers every address, on a direct call as well.
    fn neighbors(a: &Address) -> Option<Vec<Address>> {
        if !Self::is_valid_address(a) || Self::is_null_geometry(a) {
            return Some(Vec::new());
        }
        let zone = zone_from_steps(a.base, a.digits());
        let res = a.len() as i64;
        let mut out: Vec<Address> = Vec::with_capacity(6);
        for n in get_neighbors(&zone) {
            let Some(z7) = n.and_then(|z| from_7h(&z)) else {
                continue;
            };
            let nb = z7_address(z7, res);
            if nb != *a && !out.contains(&nb) {
                out.push(nb);
            }
        }
        Some(out)
    }

    /// The zone's parents as the engine's `getZoneParents` lists them
    /// (`RI7H_Z7.ec:550-556`): the zone taken to its I7H form, as `to7H` takes it; its
    /// parents found by `I7HZone::getParents` (`RI7H.ec:1261-1358`); and each converted back
    /// to Z7, as `from7H` converts it. The first is the engine's primary parent, `parent0`,
    /// which is the zone the identifier names with its last digit dropped. A zone that is
    /// the centroid child of that parent, which is one whose last digit is 0, has it alone;
    /// any other zone lies across the boundary of two zones of the level above and has
    /// both, the second found by quantising a point a hundredth short of one of the zone's
    /// own vertices. A zone of level 0 has none.
    ///
    /// A departure from the engine: two kinds of entry that the engine lists are dropped,
    /// as this crate never hands them out. They are an identifier that the engine cannot
    /// read back, for which `from_7h` answers nothing here, and a repeat of an entry already
    /// listed, where two I7H zones are given one identifier. Both occur only in the broken
    /// seams, along the edges from base cell 0 to 1 and from 1 to 6 at resolutions 15 to
    /// 19, where the engine's lists are not consistent with themselves; the list is
    /// otherwise the engine's, entry for entry. There the account above does not hold of
    /// every zone: the engine may find no second parent, its primary parent need not be the
    /// zone the identifier names with its last digit dropped, and a zone may have no parent
    /// at all. On the three grids the parents of `0000000000000001644` are
    /// `000000000000000005` and `000000000000000000`, and `00000000000000001311` has none.
    /// Neither of those two identifiers is the zone found at its own centroid, and the same
    /// held of every zone with no parent, or with another primary parent, that was compared
    /// with the engine: as measured, they are identifiers that a text can name and that no
    /// position was seen to quantise to, whereas a zone with one parent that is no centroid
    /// child is met at zones obtained from positions (see [`crate::Grid::parents`] for what
    /// was measured).
    ///
    /// An address that fails [`HexA7::is_valid_address`], and one of the twentieth level,
    /// which has no geometry (see [`HexA7::is_null_geometry`]), have no parents, as `Grid`
    /// answers them and as the engine answers a zone of the twentieth level; so the method
    /// answers every address, on a direct call as well, and never `None`.
    fn parents(a: &Address) -> Option<Vec<Address>> {
        if !Self::is_valid_address(a) || Self::is_null_geometry(a) || a.is_empty() {
            return Some(Vec::new());
        }
        let zone = zone_from_steps(a.base, a.digits());
        let res = a.len() as i64 - 1;
        let mut out: Vec<Address> = Vec::with_capacity(2);
        for p in get_parents(&zone) {
            let Some(z7) = from_7h(&p) else {
                continue;
            };
            let parent = z7_address(z7, res);
            if !out.contains(&parent) {
                out.push(parent);
            }
        }
        Some(out)
    }

    /// The zone's children as the engine's `getZoneChildren` lists them
    /// (`RI7H_Z7.ec:558-564`): the zone taken to its I7H form; its children found by
    /// `I7HZone::getChildren` (`RI7H.ec:2634-2719`); and each converted back to Z7. They are
    /// thirteen for a hexagon and eleven for a pentagon, in the engine's order: first the
    /// primary children, seven or six, the centroid child at their head, which are the
    /// zones the identifier names with one digit more, though not in the order of the
    /// digits; and then one zone across each vertex of the centroid child, each of which
    /// lies across this zone's boundary and is a child of a neighbour too. A zone of level
    /// 19, the finest, has none.
    ///
    /// A departure from the engine: three kinds of entry that the engine lists are dropped,
    /// as this crate never hands them out. They are its null zone, where the point stepped
    /// to lies in no cell; an identifier that the engine cannot read back, for which
    /// `from_7h` answers nothing here; and a repeat of an entry already listed. All three
    /// occur only in the broken seams, along the edges from base cell 0 to 1 and from 1 to
    /// 6, in the lists of zones of resolutions 14 to 18, where the engine's lists are not
    /// consistent with themselves; the list is otherwise the engine's, entry for entry.
    /// There a list may be shorter, a zone need not be among the children of its parents,
    /// nor a zone that the identifier names with one digit more among the children: on the
    /// three grids the engine's children of `0000000000000005` hold its null zone four
    /// times, and those of `00052626050026015` hold `012222222222222220`, which it cannot
    /// read back, four times, twice among the primary children.
    ///
    /// An address that fails [`HexA7::is_valid_address`], and one of the twentieth level,
    /// have no children, as for [`HexA7::parents`]; the method answers every address and
    /// never `None`.
    fn children(a: &Address) -> Option<Vec<Address>> {
        if !Self::is_valid_address(a) || Self::is_null_geometry(a) {
            return Some(Vec::new());
        }
        let zone = zone_from_steps(a.base, a.digits());
        let res = a.len() as i64 + 1;
        let mut out: Vec<Address> = Vec::with_capacity(13);
        for c in get_children(&zone) {
            let Some(z7) = c.and_then(|z| from_7h(&z)) else {
                continue;
            };
            let child = z7_address(z7, res);
            if !out.contains(&child) {
                out.push(child);
            }
        }
        Some(out)
    }

    /// The first of the zone's parents, as [`HexA7::parents`] lists them, that is itself a
    /// centroid child, and `Some(None)` where none is: what
    /// `RhombicIcosahedral7H::getZoneCentroidParent` defines (`RI7H.ec:129-138`).
    ///
    /// The rule applied to each parent is [`HexA7::is_centroid_child`], the last digit of its
    /// Z7 identifier, which is what the engine's `isZoneCentroidChild` answers on these grids
    /// (`RI7H_Z7.ec:81-103`). The engine's `getZoneCentroidParent` tests instead the geometric
    /// property of the I7H zone, that it is the centroid child of its own primary parent
    /// (`RI7H.ec:2978-2992`). The two rules coincide away from the broken seams and may part
    /// within them, where the last digit is followed all the same: of the parents of
    /// `0000000000000001644`, `000000000000000005` and `000000000000000000`, the second is
    /// answered. An identifier that the engine cannot read back is not among the parents
    /// here, and so is never the centroid parent, though it may end in 0 and stand first in
    /// the engine's list.
    ///
    /// A departure from the engine: its `getZoneCentroidParent` answered the null zone at
    /// every zone of these grids at which it was asked. That function asks for the zone's
    /// parents through the class of the grid (`RI7H.ec:132`), whose own `getZoneParents`
    /// then reads the I7H identifier it is handed as a Z7 one (`RI7H_Z7.ec:540-543`,
    /// `:550-556`). The definition is followed here, not that outcome.
    fn centroid_parent(a: &Address) -> Option<Option<Address>> {
        let parents = Self::parents(a)?;
        Some(
            parents
                .into_iter()
                .find(|p| Self::is_centroid_child(p) == Some(true)),
        )
    }

    /// Whether the zone is a centroid child, as `Z7Zone::isCentroidChild` reads the
    /// identifier (`RI7H_Z7.ec:81-103`): its last digit is 0. A zone of level 0 has no
    /// digit and is none, and an address that fails [`HexA7::is_valid_address`] is none
    /// either. An address of the twentieth level is read by its last digit as any other,
    /// as the engine reads it, although it has no parent.
    fn is_centroid_child(a: &Address) -> Option<bool> {
        Some(Self::is_valid_address(a) && a.digits().last() == Some(&0))
    }

    /// The aperture-7 grids order their sub-zones as `I7HZone::iterateI7HSubZones` generates
    /// them (`RI7H.ec:3260-3705`): in scanlines across the zone, from the vertex of it that
    /// stands uppermost on the icosahedral net.
    const HAS_SUB_ZONE_ORDER: bool = true;

    /// The number of sub-zones `depth` levels below the zone, as the engine's `countSubZones`
    /// answers (`RI7H_Z7.ec:492-500`, the closed form of `getSubZonesCount` at
    /// `RI7H.ec:3063-3116`): for a hexagon `7^d + 5 * 7^((d - 1) / 2) + 1` at an odd depth `d`
    /// and `7^d + 7^(d / 2) - 1` at an even one, which is 13 at depth 1 and 55 at depth 2; for
    /// a pentagon five sixths of that, raised to the next integer, 11 and 46; and 1 at depth 0.
    ///
    /// `None` where the sub-zones would lie beyond level 19, the finest that has geometry,
    /// which `Grid` refuses first, and for an address that fails
    /// [`HexA7::is_valid_address`]: so this method and the three below answer every address at
    /// every depth, on a direct call as well.
    fn count_sub_zones(a: &Address, depth: u8) -> Option<u64> {
        hex_a7_subzones::count(a, depth)
    }

    /// The first sub-zone `depth` levels below the zone, in the engine's order, as its
    /// `getFirstSubZone` answers (`RI7H_Z7.ec:576-579`, `RI7H.ec:3118-3124`): the centroid
    /// beside the zone's uppermost vertex, quantised at the sub-zone level. `None` as for
    /// [`HexA7::count_sub_zones`].
    ///
    /// Where no zone holds that centroid, which happens only in the broken seams, the answer
    /// is `Some` of the address that the Z7 packing writes as its null identifier (base cell
    /// 15 and no digit), never `None`: the order has an entry there, and the entry is no zone.
    ///
    /// Two departures from the engine. At depth 0 this answers the zone itself, as `Grid` does
    /// before it asks, where the engine names a neighbour at most zones. And at the pentagon
    /// of the south pole (base cell 11, every digit 0) at an odd level, at depth 2, where the
    /// engine's call ends the process, turning a direction that it was not given
    /// (`RI7H.ec:3246-3255`), this answers the first entry of the order.
    fn first_sub_zone(a: &Address, depth: u8) -> Option<Address> {
        hex_a7_subzones::first_sub_zone(a, depth)
    }

    /// Every sub-zone `depth` levels below the zone, in the engine's order, as its
    /// `getSubZones` answers (`dggrs.ec:195-222`, through `getSubZoneCRSCentroids` at
    /// `RI7H.ec:615-644`): the sub-zone centroids, generated as the engine generates them,
    /// each quantised at the sub-zone level. The order is the engine's at every zone, entry
    /// for entry; nothing in it is corrected.
    ///
    /// An entry is the address that the Z7 packing writes as its null identifier where the
    /// engine's entry is its null zone or an identifier that it cannot read back, which
    /// happens only in the broken seams, along the edges from base cell 0 to 1 and from 1 to
    /// 6, among sub-zones of level 15 and finer. The length of the order and the place of
    /// every other entry are the engine's, and whatever else its order holds there is kept:
    /// a zone named twice, a zone that is no descendant of this one.
    ///
    /// `None` as for [`HexA7::count_sub_zones`], and above
    /// [`crate::grid::max_materialised_sub_zones`], which `Grid` refuses before it asks, so
    /// that a direct caller of this trait method cannot materialise gigabytes of addresses
    /// either. The engine's own list ends the process from `2^28` centroids, and answers
    /// nothing from `2^32`.
    fn sub_zones(a: &Address, depth: u8) -> Option<Vec<Address>> {
        hex_a7_subzones::sub_zones(a, depth)
    }

    /// The sub-zone at `index`, `depth` levels below the zone, found through the generator's
    /// own search, without building [`HexA7::sub_zones`]'s whole sequence: it reckons the
    /// width of every scanline before the one that holds the index, some
    /// `(10 / 3) * 7^(depth / 2)` of them at most, and generates one centroid. Ports
    /// `getSubZoneAtIndex` (`RI7H.ec:241-256`), and is not limited by the length of the order.
    ///
    /// One departure from the engine: at a pentagon, at an odd depth, the engine's function
    /// answers another zone than its `getSubZones` lists at that index, from the scanline
    /// after the pentagon's own (it reads a width that it recorded only while generating that
    /// scanline, `RI7H.ec:3401-3402`). Here the answer is the entry of the order at every
    /// index. A second, at depth 0, where this answers the zone itself: see
    /// [`HexA7::first_sub_zone`].
    ///
    /// `None` as for [`HexA7::count_sub_zones`], and at an `index` at or beyond the count; a
    /// direct caller of this trait method must tell those apart from "no override" itself, as
    /// [`crate::Topology::sub_zone_at_index`] documents. An entry that is no zone is answered
    /// as in [`HexA7::sub_zones`].
    fn sub_zone_at_index(a: &Address, depth: u8, index: u64) -> Option<Address> {
        hex_a7_subzones::sub_zone_at_index(a, depth, index)
    }

    /// The index of `sub` in the order of the sub-zones of `parent`, by the engine's own walk,
    /// `getSubZoneIndex` (`RI7H_Z7.ec:599-602`, `RI7H.ec:224-239`): the zone must pass the
    /// engine's test of a sub-zone, `zoneHasSubZone` (`RI7H.ec:258-308`), and the index is that
    /// of the first centroid of the order, generated from its head, that lies within 1e-11 of
    /// the zone's own in each coordinate. That index is answered where the entry of the order
    /// there is `sub`. Where the test accepts the zone and no centroid lies so near, or the
    /// order holds another zone at the walk's index, which happens in the broken seams, `sub` is
    /// sought in the order itself, generated entry by entry, and its first place is the index,
    /// or `Some(None)` where the order does not hold it. Where the test refuses the zone the
    /// answer is `Some(None)`, as the engine's -1.
    ///
    /// Three departures from the engine: the search in the order, where the engine answers -1
    /// or the index of another zone for a zone that its test accepts; and two of the crate's
    /// own rule at `Grid`: of two zones of one level a zone is its own sub-zone at index 0 and
    /// no other is, where the engine answers 0 for any two; and the walk is not made over an
    /// order longer than [`crate::grid::max_materialised_sub_zones`], since it costs one
    /// centroid for each entry before the one sought, and the search a quantisation for each.
    /// `None` there, for an address that fails [`HexA7::is_valid_address`], and for a `sub`
    /// beyond level 19.
    fn sub_zone_index(parent: &Address, sub: &Address) -> Option<Option<u64>> {
        hex_a7_subzones::sub_zone_index(parent, sub)
    }

    /// True when this address has no geometry at all, DGGAL's `nullZone` outcome.
    ///
    /// Z7 level 20 is representable, since the 64-bit packing has twenty
    /// direction-digit slots, but it is not a valid DGGAL zone: the Z7 to 7H
    /// conversion (`to7H`, RI7H_Z7.ec:348-352, "does not support level 20 zones")
    /// returns nullZone, so the WGS84 centroid and vertices degenerate to the zero
    /// point. Verified against DGGAL for all three aperture-7 grids: a level-20
    /// zone yields a centroid of (0, 0) and six (0, 0) vertices, on both the
    /// hexagon and the pentagon path, while level 19 behaves normally and a
    /// pentagon-path level-19 zone correctly yields five vertices.
    ///
    /// `Grid` consults this before calling `planar_centroid` or
    /// `planar_vertices`; without it, this topology would compute a
    /// plausible coordinate for a zone that DGGAL declares does not exist.
    fn is_null_geometry(a: &Address) -> bool {
        a.len() >= NULL_GEOMETRY_LEVEL
    }

    /// The cell's centroid as a planar point in 5x6 space. `face` carries the `-1`
    /// sentinel, "derive it from x and y"; the `Grid` runs the projection's inverse
    /// to obtain a latitude and longitude.
    ///
    /// **This is not DGGAL's `getZoneCRSCentroid`, and must not be read as it.** The
    /// coordinates returned here are canonically folded, through the internal
    /// `centroid5x6`, which is the fold the eC performs at the head of `getVertices`
    /// (RI7H.ec:1806-1809) rather than in the `centroid` property itself. DGGAL's
    /// accessor returns the property unfolded, so on a cell whose centroid wraps the
    /// two differ by exactly five in both coordinates: the engine reports a small
    /// negative abscissa where this reports the same point five units along. It was
    /// measured over the 16,163 zones of the boundary census used to prove the
    /// vertex port, and it happens at 80 of them, all at resolutions 15, 17 and 19;
    /// at every other zone the two agree bit for bit, and in all 80 the difference
    /// is exactly five and not a rounding. That was measured on x86-64 Linux with
    /// glibc, where the engine and this crate share one `libm`; it says nothing of
    /// wasm32 or of any target whose floating-point functions come from another
    /// library.
    ///
    /// The fold is deliberate and is not a divergence to be closed. This function is
    /// a primitive that feeds the inverse projection, and the geographic centroid it
    /// leads to is the engine's, which is what the oracle suite enforces; the two
    /// planar forms name the same point on the icosahedral net. A caller who wants
    /// the engine's own CRS coordinates should not take them from here.
    ///
    /// # Panics
    ///
    /// `a` must satisfy [`HexA7::is_valid_address`]. A base cell above 11 or a
    /// digit above 7 indexes past the end of a lookup table and panics; a digit
    /// of 7, the terminator, ends the walk there without a panic, and the answer
    /// is then meaningless. `Grid` checks the address on every path before
    /// calling this.
    fn planar_centroid(a: &Address) -> PlanarPoint {
        let (cx, cy) = centroid5x6(&zone_from_steps(a.base, a.digits()));
        PlanarPoint {
            face: -1,
            x: cx,
            y: cy,
        }
    }

    /// The cell's boundary vertices, six for a hexagon and five for a pentagon, as
    /// planar points in 5x6 space, each carrying the `-1` face sentinel. The
    /// resolution, that is the Z7 level, is the number of digits.
    ///
    /// They are the vertices of `I7HZone::getVertices` (`RI7H.ec:1799-1895`), taken about
    /// the zone's centroid as the engine's `getZoneWGS84Vertices` takes them
    /// (`RI7H.ec:454-467`), without the `canonicalize5x6` it then applies to each;
    /// `Grid::vertices` records why that omission changes no answer on these grids.
    ///
    /// # Panics
    ///
    /// `a` must satisfy [`HexA7::is_valid_address`], on the same terms as for
    /// [`HexA7::planar_centroid`]. `Grid` checks the address on every path before
    /// calling this.
    fn planar_vertices(a: &Address) -> Vec<PlanarPoint> {
        let zone = zone_from_steps(a.base, a.digits());
        let (cx, cy) = centroid5x6(&zone);
        let level = a.len() as i64;
        verts5x6(&zone, cx, cy, level)
            .into_iter()
            .map(|(vx, vy)| PlanarPoint {
                face: -1,
                x: vx,
                y: vy,
            })
            .collect()
    }
}

impl TopologyPlumbing for HexA7 {
    /// The ring of [`HexA7::planar_refined_ring`] as `getRefinedVertices` asks it of the zone
    /// for WGS84, with `crs84` true (`RI7H.ec:544`).
    fn planar_refined_vertices(
        _: Private,
        a: &Address,
        n_divisions: i32,
    ) -> Option<Vec<(f64, f64)>> {
        Some(Self::planar_refined_ring(a, true, n_divisions))
    }

    /// The engine's `compactZones` on these grids (`RI7H_Z7.ec:581-597`): each zone taken to
    /// its I7H form, the set compacted there by the function `compact` of this module, and
    /// each zone of the answer converted back to Z7, in the order of the I7H values, which is
    /// that of the lattice within a level and not that of the Z7 identifiers.
    ///
    /// The two conversions are the engine's, so that the answer names each zone as the engine
    /// names it: an address that the engine reads as another zone, the child a pentagon does
    /// not have or an address of a broken seam, comes back as the address of that zone. An
    /// address of twenty digits is the null zone to the eC's `to7H` (`RI7H_Z7.ec:348-352`),
    /// which the set leaves out, and it is left out here.
    ///
    /// A departure from the engine: a zone of the answer for which `from_7h` answers nothing,
    /// where the engine prints an identifier that it cannot read back, is left out, as this
    /// crate never hands such an identifier out. It was met in the broken seams alone, at
    /// addresses of levels 17 and 19 built from digits, which are not the zone found at their
    /// own centroid, and at no set made of zones found from positions. Two zones are besides
    /// beyond the grid, the last two of the twelve that the short cut of that function builds,
    /// which the engine converts to its null zone and no set of one level reaches.
    fn compact(_: Private, zones: &[Address]) -> Option<Vec<Address>> {
        let zones: Vec<Z> = zones
            .iter()
            .filter(|a| !Self::is_null_geometry(a))
            .map(|a| zone_from_steps(a.base, a.digits()))
            .collect();
        Some(
            compact(&zones)
                .iter()
                .filter_map(|z| Some(z7_address(from_7h(z)?, zone_level(z))))
                .collect(),
        )
    }

    /// `7^(level / 2)` cells along a rhombus: the grid of the even level at or below `level`,
    /// whose index within a rhombus is `row * 7^levelI49R + col` (`RI7H.ec:866-871`). Nought
    /// beyond level 19.
    fn lattice_edge(_: Private, level: u8) -> u64 {
        if usize::from(level) >= NULL_GEOMETRY_LEVEL {
            return 0;
        }
        pow7(i64::from(level / 2)) as u64
    }

    /// The zones of one cell of a rhombus: one at an even level; at an odd level its seven
    /// sub-hexagons, or six in the cell at row 0 and column 0, which is a pentagon's
    /// (`RI7H.ec:866-871`). Each comes with the value of its `I7HZone`, the key of the engine's
    /// order, which is not that of the Z7 identifiers: see `i7h_value`.
    fn cell_zones(_: Private, level: u8, root: u8, row: u64, col: u64) -> Vec<(Address, u64)> {
        let p = Self::lattice_edge(Private, level);
        if root > 9 || row >= p || col >= p {
            return Vec::new();
        }
        lattice_cell(level, root, row, col)
    }

    /// The zones of a polar root: the polar pentagon at an even level, its six sub-hexagons at
    /// an odd one, the roots 10 and 11 of `rootRhombus` (`RI7H.ec:866-871`).
    fn polar_zones(_: Private, level: u8, root: u8) -> Vec<(Address, u64)> {
        if usize::from(level) >= NULL_GEOMETRY_LEVEL || !(10..=11).contains(&root) {
            return Vec::new();
        }
        lattice_cell(level, root, 0, 0)
    }

    /// The root, the row and the column of the I7H zone that the address reads to, as the eC's
    /// `to7H` reads it (`RI7H_Z7.ec:298-370`), and the value of that zone.
    ///
    /// `None` for an address that [`HexA7::is_valid_address`] refuses, for one of twenty digits,
    /// for one that names the child a pentagon does not have, and for any other that the cell it
    /// reads to does not host under that address: an answer is given exactly where
    /// `cell_zones` or `polar_zones` at the place answered holds `a`.
    fn locate(_: Private, a: &Address) -> Option<(u8, u64, u64, u64)> {
        if !Self::is_valid_address(a) || Self::is_null_geometry(a) {
            return None;
        }
        let z = zone_from_steps(a.base, a.digits());
        let (hosted, key) = lattice_zone(a.len() as u8, z)?;
        (hosted == *a).then_some((z.root as u8, z.row as u64, z.col as u64, key))
    }
}

impl HexA7 {
    /// The zone's refined boundary in the 5x6 plane, before any projection: the ring of
    /// `I7HZone::getBaseRefinedVertices(crs84, nDivisions)` (`RI7H.ec:2167-2331`), ported as
    /// [`get_base_refined_vertices`], for the zone the address names. `crs84` is true where the
    /// ring is bound for WGS84 and false where it stays in the 5x6 CRS, and `n_divisions` is the
    /// eC's `int`, the number of parts each edge is divided into.
    ///
    /// The address must name a zone with geometry. A Z7 identifier of level 20 is the null zone
    /// to the eC's `to7H` (`RI7H_Z7.ec:348-352`), which this function cannot be handed: the
    /// caller answers for it first, as it does before [`HexA7::planar_vertices`].
    ///
    /// # Panics
    ///
    /// `a` must satisfy [`HexA7::is_valid_address`], on the same terms as for
    /// [`HexA7::planar_centroid`].
    pub(crate) fn planar_refined_ring(
        a: &Address,
        crs84: bool,
        n_divisions: i32,
    ) -> Vec<(f64, f64)> {
        get_base_refined_vertices(&zone_from_steps(a.base, a.digits()), crs84, n_divisions)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Address, Topology};

    #[test]
    fn aperture_is_7() {
        assert_eq!(HexA7::APERTURE, 7);
    }
    #[test]
    fn lisbon_res5_planar_centroid_matches_py4dggs() {
        let c = HexA7::planar_centroid(&Address::new(0, &[6, 4, 1, 5, 6]));
        assert_eq!(c.face, -1);
        assert!(
            (c.x - 1.9620991253644315).abs() < 1e-13 && (c.y - 1.3556851311953353).abs() < 1e-13,
            "{c:?}"
        );
    }
    #[test]
    fn vertex_counts_hexagon_and_pentagon() {
        assert_eq!(
            HexA7::planar_vertices(&Address::new(0, &[6, 4, 1, 5, 6])).len(),
            6
        );
        assert_eq!(
            HexA7::planar_vertices(&Address::new(0, &[0, 0, 0, 0, 0])).len(),
            5
        );
    }
    /// The bounds the tables are built to: twelve base cells and seven
    /// directions. A base cell of 12 is representable in Z7's four base bits
    /// and must be refused here, or `ROOT_MAP[12]` would end the process.
    #[test]
    fn valid_address_bounds_the_base_cell_and_the_digits() {
        assert!(HexA7::is_valid_address(&Address::new(11, &[6, 0, 6])));
        assert!(!HexA7::is_valid_address(&Address::new(12, &[])));
        assert!(!HexA7::is_valid_address(&Address::new(15, &[])));
        assert!(!HexA7::is_valid_address(&Address::new(0, &[1, 7])));
    }
    #[test]
    fn level_20_has_null_geometry_level_19_does_not() {
        assert!(HexA7::is_null_geometry(&Address::new(1, &[1; 20])));
        assert!(!HexA7::is_null_geometry(&Address::new(1, &[1; 19])));
    }
    #[test]
    fn quantize_anchor_lisbon() {
        use crate::projections::Isea;
        use crate::{GridConfig, Projection};
        let g = Isea::build_geometry(&GridConfig::default());
        let p = Isea::forward(&g, 38.7223, -9.1393);
        assert_eq!(
            HexA7::quantize(p, 5),
            Some(Address::new(0, &[6, 4, 1, 5, 6]))
        );
        assert_eq!(HexA7::quantize(p, 1), Some(Address::new(0, &[6])));
        assert_eq!(HexA7::quantize(p, 0), Some(Address::new(0, &[])));
    }
    /// Points far outside the 5x6 layout at which the root rhombus that `from_7h`
    /// looks up left the twelve-entry table, and the lookup panicked, in a release
    /// build and a debug one alike. Each has no cell, and answers the null zone.
    #[test]
    fn a_root_outside_the_table_is_the_null_zone() {
        for (x, y) in [
            (-11.0, -11.0),
            (-20.0, -7.0),
            (-1e3, -1e2),
            (-1e6, -1e6),
            (-1e12, f64::NAN),
            (-1e18, -1e18),
        ] {
            let p = PlanarPoint { face: -1, x, y };
            assert_eq!(HexA7::quantize(p, 0), None, "({x}, {y}) at resolution 0");
        }
    }
    #[test]
    fn quantize_of_centroid_is_the_cell() {
        let a = Address::new(0, &[6, 4, 1, 5, 6]);
        assert_eq!(HexA7::quantize(HexA7::planar_centroid(&a), 5), Some(a));
    }

    /// Odd-resolution anchors, the case py4dggs quantises wrongly from resolution
    /// 13 upwards. Every expected value here was taken from live DGGAL through
    /// `dggal-oracle`, never from py4dggs, which disagrees with all twenty of them.
    ///
    /// The five points are chosen to cover the branches that the Lisbon anchor
    /// alone does not: the two poles, where `containsPoint` uses the eC's polar
    /// pentagon outline; the antipode of the orientation axis, on a south rhombus;
    /// and an ordinary southern-hemisphere point.
    #[test]
    fn odd_resolution_anchors_match_dggal() {
        use crate::projections::Isea;
        use crate::{GridConfig, Projection};
        // (latitude, longitude, resolution, expected digits below base cell)
        let cases: &[(f64, f64, u8, u64, &str)] = &[
            (38.7223, -9.1393, 13, 0, "6415654630655"),
            (38.7223, -9.1393, 15, 0, "641565463065523"),
            (38.7223, -9.1393, 17, 0, "64156546306552355"),
            (38.7223, -9.1393, 19, 0, "6415654630655235512"),
            (89.9, 45.0, 13, 0, "5151513414556"),
            (89.9, 45.0, 15, 0, "515151341455626"),
            (89.9, 45.0, 17, 0, "51515134145562643"),
            (89.9, 45.0, 19, 0, "5151513414556264316"),
            (
                -89.996_017_488_839,
                161.037_209_077_084_53,
                13,
                11,
                "4545454556604",
            ),
            (
                -89.996_017_488_839,
                161.037_209_077_084_53,
                15,
                11,
                "454545455660423",
            ),
            (
                -89.996_017_488_839,
                161.037_209_077_084_53,
                17,
                11,
                "45454545566042304",
            ),
            (
                -89.996_017_488_839,
                161.037_209_077_084_53,
                19,
                11,
                "4545454556604230464",
            ),
            (-31.717_474_411_461_1, 168.8, 13, 11, "3230136021525"),
            (-31.717_474_411_461_1, 168.8, 15, 11, "323013602152524"),
            (-31.717_474_411_461_1, 168.8, 17, 11, "32301360215252465"),
            (-31.717_474_411_461_1, 168.8, 19, 11, "3230136021525246523"),
            (-33.9, 151.2, 13, 10, "1550516413015"),
            (-33.9, 151.2, 15, 10, "155051641301555"),
            (-33.9, 151.2, 17, 10, "15505164130155520"),
            (-33.9, 151.2, 19, 10, "1550516413015552042"),
        ];
        let g = Isea::build_geometry(&GridConfig::default());
        for &(lat, lon, res, base, digits) in cases {
            let expected =
                Address::new(base, &digits.bytes().map(|b| b - b'0').collect::<Vec<u8>>());
            let got = HexA7::quantize(Isea::forward(&g, lat, lon), res);
            assert_eq!(got, Some(expected), "({lat}, {lon}) at resolution {res}");
        }
    }

    /// The null zone is a real answer, not a failure.
    ///
    /// Within roughly a quarter of a degree of the north pole, at odd resolutions
    /// of 15 and above, the odd-level candidate search finds no containing cell,
    /// and the eC returns nullZone there (RI7H.ec:1489 and :1594). DGGAL prints
    /// `(null)` for this very point at resolution 17, so the absent value is the
    /// faithful answer rather than a gap.
    ///
    /// The same point at an even resolution is an ordinary zone, which is what
    /// keeps this from being a blanket "the poles do not work".
    #[test]
    fn near_pole_odd_resolution_is_the_null_zone() {
        use crate::projections::Isea;
        use crate::{GridConfig, Projection};
        let g = Isea::build_geometry(&GridConfig::default());
        let p = Isea::forward(&g, 89.999_129_147_255_16, 9.333_761_199_597_859);
        assert_eq!(HexA7::quantize(p, 17), None);
        assert_eq!(
            HexA7::quantize(p, 16),
            Some(Address::new(
                0,
                &[5, 1, 5, 1, 5, 1, 5, 1, 5, 1, 5, 3, 0, 2, 5, 4]
            ))
        );
    }

    /// Containment inside and around a polar pentagon, at odd resolutions, which
    /// is the path that uses [`polar_pentagon_outline`] rather than a traced
    /// boundary.
    ///
    /// The outline test below pins that polygon's shape; this pins what the shape
    /// is for, namely that `from_centroid`'s candidate search picks the right cell
    /// near a pole. The five points around each pentagon's centre fan out in
    /// latitude and longitude so that several of its half-planes decide the
    /// answer, not just the one nearest the centre. Every expected value was taken
    /// from live DGGAL, so this covers the path with the engine absent.
    #[test]
    fn polar_pentagon_containment_anchors_match_dggal() {
        use crate::projections::Isea;
        use crate::{GridConfig, Projection};
        // (latitude, longitude, resolution, base cell, digits below it)
        let cases: &[(f64, f64, u8, u64, &str)] = &[
            // The north polar pentagon's own centre, then four points around it.
            (58.397_145_907_431_117, 11.2, 3, 0, "000"),
            (58.397_145_907_431_117, 11.2, 5, 0, "00000"),
            (58.397_145_907_431_117, 11.2, 7, 0, "0000000"),
            (59.5, 11.2, 5, 0, "00050"),
            (59.5, 11.2, 7, 0, "0005025"),
            (57.3, 11.2, 5, 0, "00063"),
            (57.3, 11.2, 7, 0, "0006336"),
            (58.4, 13.8, 5, 0, "00012"),
            (58.4, 13.8, 7, 0, "0001251"),
            (58.4, 8.6, 5, 0, "00046"),
            (58.4, 8.6, 7, 0, "0004614"),
            // The south polar pentagon, the other arm of the outline.
            (-58.397_145_907_431_117, -168.8, 3, 11, "000"),
            (-58.397_145_907_431_117, -168.8, 5, 11, "00000"),
            (-58.397_145_907_431_117, -168.8, 7, 11, "0000000"),
            (-59.5, -168.8, 5, 11, "00040"),
            (-59.5, -168.8, 7, 11, "0004034"),
            (-57.3, -168.8, 5, 11, "00021"),
            (-57.3, -168.8, 7, 11, "0002112"),
            (-58.4, -166.2, 5, 11, "00062"),
            (-58.4, -166.2, 7, 11, "0006256"),
            (-58.4, -171.4, 5, 11, "00012"),
            (-58.4, -171.4, 7, 11, "0001251"),
        ];
        let g = Isea::build_geometry(&GridConfig::default());
        for &(lat, lon, res, base, digits) in cases {
            let expected =
                Address::new(base, &digits.bytes().map(|b| b - b'0').collect::<Vec<u8>>());
            let got = HexA7::quantize(Isea::forward(&g, lat, lon), res);
            assert_eq!(got, Some(expected), "({lat}, {lon}) at resolution {res}");
        }
    }

    /// The polar pentagon outline is the only polygon in `containsPoint` that is
    /// not traced from the offset table, so its shape is pinned here: seventeen
    /// vertices, both poles present, and the ring closing on `ab`.
    #[test]
    fn polar_pentagon_outline_has_the_ecs_seventeen_vertices() {
        for root in [0xA, 0xB] {
            let out = polar_pentagon_outline(root, 1.0 / 7.0);
            assert_eq!(out.len(), 17, "root {root:#x}");
            // The pole pair the eC appends at RI7H.ec:2451-2452 and :2476-2477.
            let poles: [(f64, f64); 2] = if root == 0xA {
                [(5.0, 4.0), (1.0, 0.0)]
            } else {
                [(0.0, 2.0), (4.0, 6.0)]
            };
            assert_eq!(&out[14..16], &poles[..], "root {root:#x}");
        }
    }

    /// The side test in the order the compiled engine performs it, on the case that
    /// forced the change.
    ///
    /// The planar point (4.5, 5.5), where base cells 6 and 10 meet, at resolution 9
    /// on any of the three aperture-7 grids. The point lies exactly on the edge two
    /// candidate cells share. The endpoints below are the last side of the refined
    /// polygon of zone `06131313131`, taken bit for bit from this crate and checked
    /// against the engine's own `getZoneCRSVertices` in the 5x6 CRS, where they are
    /// identical, so nothing but the arithmetic separates the two sides here.
    ///
    /// In the eC's written order the side value is -5.421010862427522e-20, so both
    /// candidates are rejected and the quantiser answers the null zone; in the order
    /// `libdggal.so` performs it the value is an exact zero and the cell is
    /// accepted, which is the engine's own answer, `06131313131`. Restoring the
    /// source order fails the first assertion.
    #[test]
    fn the_side_test_is_zero_at_the_boundary_in_the_compiled_order() {
        let (px, py) = (
            f64::from_bits(0x4012_0000_0000_0000),
            f64::from_bits(0x4016_0000_0000_0000),
        );
        let (ax, ay) = (
            f64::from_bits(0x4012_000c_ff6d_2016),
            f64::from_bits(0x4016_0002_997c_399e),
        );
        let (bx, by) = (
            f64::from_bits(0x4011_fff3_0092_dfeb),
            f64::from_bits(0x4015_fffd_6683_c662),
        );
        assert_eq!((px, py), (4.5, 5.5));

        assert_eq!(
            point_line_side(px, py, ax, ay, bx, by).to_bits(),
            0.0f64.to_bits(),
            "the compiled order must return an exact positive zero here"
        );

        // The eC's written order, `A*x + B*y + C`, kept here so that the value it
        // returns is on the record and the difference is not a matter of belief.
        fn source_order(px: f64, py: f64, ax: f64, ay: f64, bx: f64, by: f64) -> f64 {
            let dx = bx - ax;
            let dy = by - ay;
            let a = dy;
            let b = -dx;
            let c = ay * dx - dy * ax;
            a * px + b * py + c
        }
        assert_eq!(
            source_order(px, py, ax, ay, bx, by),
            -5.421010862427522e-20,
            "the eC's written order still returns its tiny negative value"
        );
    }

    /// The constants gcc folded for the four even-level tie-breaks, against the
    /// values this crate computes from the same endpoints.
    ///
    /// Each triple is `(A, B, -C)` as `libdggal.so` holds it in `.rodata`, read at
    /// the addresses named beside it. If Rust's folding of the endpoint arithmetic
    /// ever drifted from the compiler's, this test would say so rather than let a
    /// tie-break quietly answer on different coefficients.
    #[test]
    fn the_tie_break_coefficients_match_the_engines_folded_constants() {
        fn coefficients(ax: f64, ay: f64, bx: f64, by: f64) -> (u64, u64, u64) {
            let dx = bx - ax;
            let dy = by - ay;
            let c = ay * dx - dy * ax;
            (dy.to_bits(), (-dx).to_bits(), (-c).to_bits())
        }
        // RI7H.ec:1610, inlined at 0x239ea: 0x4a190, 0x4a198, 0x4a1a0.
        assert_eq!(
            coefficients(1.0, 0.5, 5.0 / 6.0, 1.0 / 6.0),
            (
                0xbfd5_5555_5555_5556,
                0x3fc5_5555_5555_5554,
                0xbfd0_0000_0000_0001
            )
        );
        // RI7H.ec:1616, inlined at 0x23d80: 0x4a198, 0x4a190, 0x4a1a0.
        assert_eq!(
            coefficients(1.0 / 6.0, 5.0 / 6.0, 0.5, 1.0),
            (
                0x3fc5_5555_5555_5554,
                0xbfd5_5555_5555_5556,
                0xbfd0_0000_0000_0001
            )
        );
        // RI7H.ec:1626, inlined at 0x237a8: 0x4a1a8, 0x4a1b0, 0x4a1b8.
        assert_eq!(
            coefficients(5.0 / 6.0, 1.0 / 6.0, 0.5, 0.0),
            (
                0xbfc5_5555_5555_5555,
                0x3fd5_5555_5555_5556,
                0xbfb5_5555_5555_5556
            )
        );
        // RI7H.ec:1632, inlined at 0x23d50: 0x4a1b0, 0x4a1a8, 0x4a1c0.
        assert_eq!(
            coefficients(0.0, 0.5, 1.0 / 6.0, 5.0 / 6.0),
            (
                0x3fd5_5555_5555_5556,
                0xbfc5_5555_5555_5555,
                0xbfb5_5555_5555_5555
            )
        );
    }

    /// The two constants gcc folded for the wrap comparisons, against the values
    /// this crate computes.
    ///
    /// `contains_point`'s `-5` wrap path subtracts `fl(5 + 1e-11)` from
    /// `bbox.tl.x` in one instruction, and `from_centroid`'s bottom-right corner
    /// test compares against `fl(10 - 1e-11) - oop`. Both folded values are in
    /// `.rodata`, read at the addresses named beside them; the assertions hold
    /// Rust's own folding to the same bits.
    #[test]
    fn the_wrap_constants_match_the_engines_folded_constants() {
        // 0x4a188, loaded by the `subsd` at 0x23475.
        assert_eq!((5.0f64 + 1e-11).to_bits(), 0x4014_0000_0000_2bfb);
        // 0x49f58, loaded at 0x235d0 and again at 0x1fc54.
        assert_eq!((10.0f64 - 1e-11).to_bits(), 0x4023_ffff_ffff_ea03);
        // 0x49f48, the ordinary folding of `5 - 1E-11` that guards both corner
        // tests and that this crate already writes as the eC does.
        assert_eq!((5.0f64 - 1e-11).to_bits(), 0x4013_ffff_ffff_d405);
    }

    /// The constants the compiler folded for the geometry, against the values this
    /// crate computes from the same eC terms.
    ///
    /// The vertices this crate computes are bit-identical to the engine's only so
    /// long as Rust's folding of these expressions agrees with gcc's. If it ever
    /// drifted, the boundaries would quietly stop matching and only a census would
    /// say so; this test says so at once. Each value is the `.rodata` word at the
    /// address named beside it in `libdggal.so`, BuildID
    /// `e75d6ab18f8b460713ebcd21df6f07922fb909c3`.
    #[test]
    fn the_geometry_constants_match_the_engines_folded_constants() {
        // 0x49b48 and 0x49b40, the rounded reciprocals that replaced the centroid
        // property's divisions by seven.
        assert_eq!(ONE_SEVENTH.to_bits(), 0x3fc2_4924_9249_2492);
        assert_eq!(TWO_SEVENTHS.to_bits(), 0x3fd2_4924_9249_2492);
        // And that the reciprocal really is not the division, which is the whole
        // reason the site is ported: they part company at `levelI49R` 2, 3, 4 and 7.
        let parted: Vec<i64> = (0..10)
            .filter(|&k| {
                let oop = 1.0 / pow7(k) as f64;
                oop * ONE_SEVENTH != oop / 7.0
            })
            .collect();
        assert_eq!(parted, vec![2, 3, 4, 7]);

        // 0x49bc0, the eC's `r = 1 / 5.0` in the polar refined arms.
        assert_eq!((1.0f64 / 5.0).to_bits(), 0x3fc9_9999_9999_999a);
        // 0x49bb0, the whole `oonp` coefficient of the north arm's `ab.x`, which the
        // compiler folded out of `a.x + (b.x - a.x) * r`.
        assert_eq!(
            (-ODD_HEX_B + (ODD_HEX_B - ODD_HEX_C) * (1.0 / 5.0)).to_bits(),
            0xbff6_6666_6666_6666
        );
        // 0x4a160, the same arm's `(c.y - b.y) * r`, folded whole. Note the last
        // digit: it is not `fl(1/15)`, which is `0x3fb1111111111111`, but the
        // rounding of this particular expression.
        assert_eq!(
            ((ODD_HEX_B - ODD_HEX_A) * (1.0 / 5.0)).to_bits(),
            0x3fb1_1111_1111_1113
        );
        assert_ne!(
            ((ODD_HEX_B - ODD_HEX_A) * (1.0 / 5.0)).to_bits(),
            (1.0f64 / 15.0).to_bits()
        );

        // 0x49ba0, 0x49ba8 and 0x49c00, the odd offsets the polar arms scale by.
        assert_eq!(ODD_HEX_C.to_bits(), 0x3fd5_5555_5555_5555);
        assert_eq!(ODD_HEX_A.to_bits(), 0x3ff5_5555_5555_5555);
        assert_eq!(ODD_HEX_B.to_bits(), 0x3ffa_aaaa_aaaa_aaab);

        // 0x4a110 and 0x4a118, the even offsets `getBaseRefinedVertices` scales by, at
        // 0x2930f and 0x2930b.
        assert_eq!(EVEN_HEX_A.to_bits(), 0x4002_aaaa_aaaa_aaab);
        assert_eq!(EVEN_HEX_B.to_bits(), 0x4012_aaaa_aaaa_aaab);
        // 0x4a128 and 0x4a130, the two bounds of the fold it applies to the centroid, compared
        // at 0x292d1 and 0x294c0.
        assert_eq!((6.0f64 + 1e-9).to_bits(), 0x4018_0000_0011_2e0c);
        assert_eq!((5.0f64 + 1e-9).to_bits(), 0x4014_0000_0011_2e0c);
    }

    /// Every site of the four polar arms, by the levels at which the member as the eC writes
    /// it is another double than the one the library computes, and [`polar_pentagon_sides`]
    /// with it. `oonp` takes ten values, one to each level of an arm, so the levels below are
    /// the whole of the matter: a site with none cannot change an answer, and is ported as
    /// compiled all the same; a site with some is pinned by the engine's rings at one of them,
    /// in `POLAR_RINGS`. The joint switch for ON-2 and ON-3 is told apart here: `ab`
    /// alone parts from the eC's order, and `d` nowhere.
    ///
    /// The members the compiler folded to whole numbers are compared with the words of the
    /// library as well: one at `.rodata` `0x49a70`, four at `0x49a80`, six at `0x49a88`, and
    /// the zero stored at `0x299af`.
    #[test]
    fn the_polar_arms_part_from_the_ecs_order_at_the_recorded_levels_alone() {
        use std::collections::BTreeMap;
        type Points = [(f64, f64)];
        let r = 1.0 / 5.0;
        let turns = |shifted: &dyn Fn(f64) -> (f64, f64)| [0.0, 1.0, 2.0, 3.0, 4.0].map(shifted);
        let mut parted: BTreeMap<&str, Vec<i64>> = BTreeMap::new();
        let mut compare = |site: &'static str, level: i64, compiled: &Points, source: &Points| {
            let levels = parted.entry(site).or_default();
            assert_eq!(compiled.len(), source.len());
            if compiled
                .iter()
                .zip(source)
                .any(|(a, b)| a.0.to_bits() != b.0.to_bits() || a.1.to_bits() != b.1.to_bits())
            {
                levels.push(level);
            }
        };
        // A member alone, the other held at zero.
        let xs = |points: &Points| points.iter().map(|p| (p.0, 0.0)).collect::<Vec<_>>();
        let ys = |points: &Points| points.iter().map(|p| (0.0, p.1)).collect::<Vec<_>>();

        for l49r in 0..=9 {
            let oonp = 1.0 / (7 * pow7(l49r)) as f64;
            let (even, odd) = (2 * l49r, 2 * l49r + 1);

            // Even north, as the eC writes it (RI7H.ec:2189-2202).
            let (a_p, b_p) = (oonp * EVEN_HEX_A, oonp * EVEN_HEX_B);
            let a = (1.0 - b_p, 0.0 - a_p);
            let b = (1.0 - a_p, 0.0 + a_p);
            let ab = ((a.0 + b.0) / 2.0, (a.1 + b.1) / 2.0);
            let offset = rotate5x6_offset(b.0 - ab.0, b.1 - ab.1, false);
            let d = (offset.0 + b.0, offset.1 + b.1);
            let north = polar_pentagon_sides(0xA, false, oonp).expect("the north pole");
            compare("EN-1", even, &north.b, &turns(&|k| (b.0 + k, b.1 + k)));
            compare(
                "EN-2",
                even,
                &xs(&north.ab),
                &xs(&turns(&|k| (ab.0 + k, 0.0))),
            );
            compare(
                "EN-3",
                even,
                &ys(&north.ab),
                &ys(&turns(&|k| (0.0, ab.1 + k))),
            );
            compare(
                "EN-4",
                even,
                &xs(&north.d),
                &xs(&turns(&|k| (d.0 + k, 0.0))),
            );
            compare(
                "EN-5",
                even,
                &ys(&north.d),
                &ys(&turns(&|k| (0.0, d.1 + k))),
            );
            assert_eq!(north.ab[0].1.to_bits(), 0, "EN-3 at level {even}");
            assert_eq!(north.d[0].0.to_bits(), 0x3ff0_0000_0000_0000, "EN-4");

            // Even south (RI7H.ec:2215-2229).
            let a = (4.0 + b_p, 6.0 + a_p);
            let b = (4.0 + a_p, 6.0 - a_p);
            let ab = ((a.0 + b.0) / 2.0, (a.1 + b.1) / 2.0);
            let offset = rotate5x6_offset(b.0 - ab.0, b.1 - ab.1, false);
            let d = (offset.0 + b.0, offset.1 + b.1);
            let south = polar_pentagon_sides(0xB, false, oonp).expect("the south pole");
            compare("ES-1", even, &south.b, &turns(&|k| (b.0 - k, b.1 - k)));
            compare(
                "ES-2",
                even,
                &xs(&south.ab),
                &xs(&turns(&|k| (ab.0 - k, 0.0))),
            );
            compare(
                "ES-3",
                even,
                &ys(&south.ab),
                &ys(&turns(&|k| (0.0, ab.1 - k))),
            );
            compare(
                "ES-4",
                even,
                &xs(&south.d),
                &xs(&turns(&|k| (d.0 - k, 0.0))),
            );
            compare(
                "ES-5",
                even,
                &ys(&south.d),
                &ys(&turns(&|k| (0.0, d.1 - k))),
            );
            assert_eq!(south.ab[0].1.to_bits(), 0x4018_0000_0000_0000, "ES-3");
            assert_eq!(south.d[0].0.to_bits(), 0x4010_0000_0000_0000, "ES-4");

            // Odd north (RI7H.ec:2266-2277).
            let (a_p, b_p, c_p) = (oonp * ODD_HEX_A, oonp * ODD_HEX_B, oonp * ODD_HEX_C);
            let a = (1.0 - b_p, 0.0 - c_p);
            let b = (1.0 - c_p, 0.0 + a_p);
            let ab = (a.0 + (b.0 - a.0) * r, a.1 + (b.1 - a.1) * r);
            let c = (1.0 + a_p, 0.0 + b_p);
            let d = (b.0 + (c.0 - b.0) * r, b.1 + (c.1 - b.1) * r);
            let north = polar_pentagon_sides(0xA, true, oonp).expect("the north pole");
            compare("ON-1", odd, &north.b, &turns(&|k| (b.0 + k, b.1 + k)));
            compare("ON-2", odd, &north.ab[..1], &[ab]);
            compare("ON-3", odd, &north.d, &turns(&|k| (d.0 + k, d.1 + k)));
            // ON-4 is of the shift alone: the eC's `ab + k`, on the library's own `ab`.
            let shifted = turns(&|k| (north.ab[0].0 + k, north.ab[0].1 + k));
            compare("ON-4", odd, &north.ab[1..2], &shifted[1..2]);
            compare(
                "ab + k beyond the first turn",
                odd,
                &north.ab[2..],
                &shifted[2..],
            );

            // Odd south (RI7H.ec:2291-2302).
            let a = (4.0 + b_p, 6.0 + c_p);
            let b = (4.0 + c_p, 6.0 - a_p);
            let ab = (a.0 + (b.0 - a.0) * r, a.1 + (b.1 - a.1) * r);
            let c = (4.0 - a_p, 6.0 - b_p);
            let d = (b.0 + (c.0 - b.0) * r, b.1 + (c.1 - b.1) * r);
            let south = polar_pentagon_sides(0xB, true, oonp).expect("the south pole");
            compare("OS-1", odd, &south.b, &turns(&|k| (b.0 - k, b.1 - k)));
            compare("OS-2", odd, &south.d, &turns(&|k| (d.0 - k, d.1 - k)));
            compare("ab - k", odd, &south.ab, &turns(&|k| (ab.0 - k, ab.1 - k)));
        }

        let every_even: Vec<i64> = (0..=18).step_by(2).collect();
        let every_odd: Vec<i64> = (1..=19).step_by(2).collect();
        let recorded: BTreeMap<&str, Vec<i64>> = [
            ("EN-1", vec![0, 4, 10, 16]),
            ("EN-2", vec![]),
            ("EN-3", vec![]),
            ("EN-4", vec![]),
            ("EN-5", every_even.clone()),
            ("ES-1", every_even),
            ("ES-2", vec![]),
            ("ES-3", vec![]),
            ("ES-4", vec![]),
            ("ES-5", vec![]),
            ("ON-1", vec![3, 9, 15, 19]),
            ("ON-2", vec![1, 9, 15]),
            ("ON-3", vec![]),
            ("ON-4", vec![11, 13, 19]),
            ("ab + k beyond the first turn", vec![]),
            ("OS-1", every_odd),
            ("OS-2", vec![]),
            ("ab - k", vec![]),
        ]
        .into_iter()
        .collect();
        assert_eq!(parted, recorded);
    }

    /// The emission towards WGS84, which the engine's planar accessor never asks for (it
    /// passes `crs84` false): the fifth side of a polar pentagon is then divided as the other
    /// four are, from the fifth corner back to the first across the interruption between
    /// `d + 4` and `ab` (RI7H.ec:2201-2202 and its three fellows). The ring is rebuilt here
    /// from the points of the ring in the 5x6 CRS, which the engine's own rings pin, one
    /// zone to an arm; the comparison with the engine is that of the projected ring.
    #[test]
    fn towards_wgs84_the_fifth_side_of_a_polar_pentagon_is_divided_as_the_others() {
        use crate::Indexing as _;
        for text in ["0000", "00000", "1100", "11000"] {
            let id = crate::indexings::Z7::from_text(text).expect("a canonical identifier");
            let address = crate::indexings::Z7::decode(id);
            // `b`, `d` and the next `ab` of each of four sides, then `b + 4`, `d + 4`, the two
            // images of the pole, and `ab`.
            let planar = HexA7::planar_refined_ring(&address, false, 1);
            assert_eq!(planar.len(), 17, "{text}");
            for n in [1, 2, 12] {
                let mut expected: Vec<(f64, f64)> = Vec::new();
                for k in 0..4 {
                    add_intermediate_points(
                        &mut expected,
                        planar[3 * k],
                        planar[3 * k + 3],
                        n,
                        Some((planar[3 * k + 1], planar[3 * k + 2])),
                        true,
                    );
                }
                add_intermediate_points(
                    &mut expected,
                    planar[12],
                    planar[0],
                    n,
                    Some((planar[13], planar[16])),
                    true,
                );
                let ring = HexA7::planar_refined_ring(&address, true, n);
                assert_eq!(ring.len(), expected.len(), "{text} at {n} divisions");
                for (i, (a, b)) in ring.iter().zip(&expected).enumerate() {
                    assert!(
                        a.0.to_bits() == b.0.to_bits() && a.1.to_bits() == b.1.to_bits(),
                        "{text} at {n} divisions: point {i} is {a:?} for {b:?}"
                    );
                }
            }
            // At one division the ring towards WGS84 is the five corners and nothing else.
            assert_eq!(
                HexA7::planar_refined_ring(&address, true, 1),
                [planar[0], planar[3], planar[6], planar[9], planar[12]],
                "{text}"
            );
        }
    }

    /// `pymin` is the library's `minsd`, which answers its second operand on a tie
    /// and on NaN, and `pymax` the library's `(a <= b) ? b : a`, which answers its
    /// first on NaN. The zeros are compared by their bits, since `+0.0 == -0.0`.
    #[test]
    fn pymin_and_pymax_are_the_librarys_operations() {
        assert_eq!(pymin(1.0, 2.0), 1.0);
        assert_eq!(pymin(2.0, 1.0), 1.0);
        assert_eq!(pymin(0.0, -0.0).to_bits(), (-0.0f64).to_bits());
        assert_eq!(pymin(-0.0, 0.0).to_bits(), 0.0f64.to_bits());
        assert_eq!(pymin(f64::NAN, 1.0), 1.0);
        assert!(pymin(1.0, f64::NAN).is_nan());

        assert!(pymax(f64::NAN, -1.0).is_nan());
        assert_eq!(pymax(-2.0, f64::NAN), -2.0);
    }

    /// The level-19 refusal of `centroidChild` (RI7H.ec:2614): an odd zone at
    /// `l49r` 9 has no centroid child, while the even zone above it has one.
    #[test]
    fn a_level_19_zone_has_no_centroid_child() {
        assert_eq!(get_centroid_child(&Z::new(9, 0, 1, 1, 2)), None);
        assert_eq!(get_centroid_child(&Z::new(9, 0, 1, 1, 1)), None);
        assert_eq!(
            get_centroid_child(&Z::new(9, 0, 1, 1, 0)),
            Some(Z::new(9, 0, 1, 1, 1))
        );
        assert!(get_centroid_child(&Z::new(8, 0, 1, 1, 2)).is_some());
    }

    /// `isEdgeHex` applies its seam test at an odd zone too (`RI7H.ec:915`), which the
    /// function used to refuse before it matched the eC: a zone with `sub_hex` 2 (six
    /// points, since `sub_hex > 1`), a north root and row 0 has `rhombusIX = 0 * 7 + 3 =
    /// 3`, less than `POW7(1) = 7`, so the eC answers true; the removed branch, which
    /// returned false whenever `sub_hex != 0`, answered false here.
    #[test]
    fn zone_is_edge_hex_answers_an_odd_zone_by_its_even_seam_test() {
        assert!(zone_is_edge_hex(&Z::new(1, 0, 0, 3, 2)));
    }

    /// The refined polygon this crate builds, against the nearest thing DGGAL exposes
    /// to the one its own containment test is given.
    ///
    /// `containsPoint` is handed `getBaseRefinedVerticesNoAlloc(false, 1, ...)`
    /// (RI7H.ec:956), for which the binding exposes no accessor at all, so until now
    /// [`base_refined_vertices`] and the two tracers beneath it rested on the
    /// identical-input census alone: an instrument built by hand for one task rather
    /// than a check that stands. `getZoneRefinedCRSVertices` in the 5x6 CRS routes
    /// through `getIcoNetRefinedVertices` to `getBaseRefinedVertices`
    /// (RI7H.ec:2167), the allocating twin of that routine, which the eC writes with
    /// the same text. It is therefore a close proxy and not the thing itself, and
    /// this test is written to say exactly where the two part company rather than to
    /// paper over it.
    ///
    /// Two places, and only two, and both are asserted rather than tolerated.
    ///
    /// **The zones whose interruption lies near a pole.** The two intermediate-point
    /// routines are not the same function either. `addIntermediatePointsNoAlloc`
    /// raises `interruptionNearPole` only under `crs84` (ri5x6.ec:732), which is
    /// false on this path, so a crossing side emits exactly `p`, `i1`, `i2` and
    /// nothing else. The allocating `addIntermediatePoints` carries four further
    /// clauses that do not test `crs84` at all (ri5x6.ec:664-671) and multiply
    /// `nDivisions` by twenty, so a side crossing near one of the four points they
    /// name comes back as a twenty-fold subdivision: 29 points where this crate's
    /// polygon has 10. Those zones are a different polygon, not a differently
    /// spelled one, and a coordinate comparison of them would be meaningless. They
    /// are therefore counted and the count asserted exactly, so that a zone may
    /// never drop out of this test merely by the shape of its answer; if the count
    /// moves, the tracing has changed and the test says so. Every other
    /// interruption-crossing zone is compared in full, and 127 of the 360 below
    /// are of that kind.
    ///
    /// **One vertex, and its cause.** At the north polar pentagon `0000000`, vertex
    /// 11, the abscissa differs by one unit in the last place, and the reason is that
    /// gcc compiled the two twins' north polar arm differently. The eC writes the
    /// fifth interruption point as `ab.x + 4` (RI7H.ec:2443). The allocating twin
    /// keeps that, at `0x2a390` onwards in `libdggal.so`, BuildID
    /// `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, giving `4.995918367346938`. The
    /// no-alloc twin moves the whole number inside the bracket, adding the packed
    /// corner pair `(5, 4)` from `0x49b50` to the offsets form of `ab` at `0x2152a`,
    /// giving `(1 + 4) + fl(-1.4)*oonp` and so `4.995918367346939`. This crate ports
    /// the no-alloc twin, because that is the one `containsPoint` uses, so **this
    /// crate's value is the right one for the purpose and the accessor's is not**.
    /// The exception is pinned by zone, index and both bit patterns, so that it
    /// cannot silently become two.
    ///
    /// That difference is the site the port of the allocating copy names ON-4, in
    /// [`polar_pentagon_sides`]: there the library puts the whole number inside the bracket
    /// at the first turn about the pole alone (`0x2a4cb`), and keeps the eC's `ab + k` at the
    /// second, third and fourth (`0x2a579`, `0x2a60c`, `0x2a658`). At this zone the two forms
    /// give the same double at the second and third turns and part at the fourth, which is
    /// vertex 11. Since the refined boundary was ported this crate ports both copies: [`get_base_refined_vertices`]
    /// answers the accessor's own value here, and the test that follows this one compares
    /// that copy with the accessor exactly, at every level and division. This test is kept
    /// beside it as the guard of the `NoAlloc` copy, on which
    /// `containsPoint` rests, and as the record of the one site at which the two copies part.
    ///
    /// What the test actually guards was checked by reverting each ported order in
    /// turn and seeing whether this test failed. Four of the seven in this path are
    /// caught: the cancellation of the centroid out of the first direction, the
    /// north polar arm's folded `ab`, its `b` fan, and the centroid's reciprocal.
    /// The other three are not, and the reason is not a gap in the sample: the
    /// hoisted `c.x + 1e-11` in the start search, the dropped second rotation and
    /// the bracketing of `ddx` and `ddy` were each computed both ways over 8,257
    /// odd-level zones and **never differed in a single bit**, so no comparison
    /// against this accessor, or any other, could distinguish them. They are ported
    /// because the engine performs them that way, not because they are observable.
    ///
    /// Gated behind the `oracle` feature, meaningful only within this repository: it
    /// calls the live engine through the test-only `dggal-oracle` crate, a
    /// dev-dependency that names only a path and no registry version, so `cargo
    /// package` strips it, and this test is never even compiled in the published
    /// crate.
    #[test]
    #[cfg(feature = "oracle")]
    fn the_refined_polygon_matches_the_engines_refined_accessor() {
        // Six base cells, including both polar pentagons (Z7 base 0 is I7H root 0xA,
        // base 11 is root 0xB), and six digit patterns whose first non-zero digit
        // avoids the omitted pentagon child under either rule, so that every
        // truncation below is a canonical identifier.
        const BASES: [&str; 6] = ["00", "01", "04", "06", "09", "11"];
        const PATTERNS: [&str; 6] = [
            "0000000000000000000",
            "1313131313131313131",
            "3232323232323232323",
            "4141414141414141414",
            "6060606060606060606",
            "1456013456013456013",
        ];
        // Only an odd level reaches the refined polygon: an even one is a whole
        // sub-hexagon and `containsPoint` is never asked about it.
        let zones: Vec<String> = BASES
            .iter()
            .flat_map(|b| {
                PATTERNS
                    .iter()
                    .flat_map(move |p| (1..=19).step_by(2).map(move |r| format!("{b}{}", &p[..r])))
            })
            .collect();
        assert_eq!(zones.len(), 360, "the zone set is fixed by construction");

        for oracle in [
            dggal_oracle::IGEO7,
            dggal_oracle::IVEA7H,
            dggal_oracle::RTEA7H,
        ] {
            let mut coords = 0usize;
            let mut seam_crossing = 0usize;
            let mut exceptions: Vec<(String, usize, u64, u64)> = Vec::new();
            for text in &zones {
                use crate::Indexing as _;
                let id = crate::indexings::Z7::from_text(text)
                    .unwrap_or_else(|e| panic!("{text} is not a canonical Z7 identifier: {e}"));
                let address = crate::indexings::Z7::decode(id);
                let zone = zone_from_steps(address.base, address.digits());
                let level = zone_level(&zone);
                let oonp = 1.0 / (7 * pow7(level / 2)) as f64;
                let ours = base_refined_vertices(&zone, oonp);

                let theirs_id = dggal_oracle::zone_from_text(oracle, text);
                assert_ne!(
                    theirs_id,
                    dggal_oracle::NULL_ZONE,
                    "{oracle} does not know {text}"
                );
                let theirs = dggal_oracle::refined_crs_vertices_5x6(oracle, theirs_id, 1);
                assert!(!theirs.is_empty(), "{oracle} {text}: no refined vertices");

                if ours.len() != theirs.len() {
                    seam_crossing += 1;
                    continue;
                }
                for (i, (a, b)) in ours.iter().zip(theirs.iter()).enumerate() {
                    coords += 1;
                    if a.0.to_bits() != b.0.to_bits() || a.1.to_bits() != b.1.to_bits() {
                        exceptions.push((text.clone(), i, a.0.to_bits(), b.0.to_bits()));
                    }
                }
            }

            // Measured, and asserted exactly so that the test cannot quietly compare
            // less than it claims: 17 of the 360 zones cross an interruption, and the
            // remaining 343 contribute 2,666 vertices each of which is compared.
            assert_eq!(
                seam_crossing, 17,
                "{oracle}: the count of interruption-crossing zones moved"
            );
            assert_eq!(
                coords, 2666,
                "{oracle}: the number of vertices compared moved"
            );

            // The single documented exception, pinned by zone, index and both bit
            // patterns. Its cause is in this test's documentation.
            assert_eq!(
                exceptions,
                vec![(
                    "0000000".to_string(),
                    11usize,
                    0x4013_fbd2_0644_f699_u64,
                    0x4013_fbd2_0644_f698_u64
                )],
                "{oracle}: the refined polygon differs from the engine's accessor \
                 somewhere other than the one documented vertex"
            );
        }
    }

    /// The edge divisions at which the planar ring is compared with the engine's: none, one,
    /// two, and the three that `getRefinedVertices` chooses by level when it is asked for the
    /// automatic refinement (`RI7H.ec:542-543`: twenty below level 3, fifteen below level 5,
    /// twelve from there on).
    #[cfg(feature = "oracle")]
    const PLANAR_REFINEMENTS: [i32; 6] = [0, 1, 2, 12, 15, 20];

    /// Where this crate's planar ring of the zone `id`, at `n` divisions and in the 5x6 CRS,
    /// leaves the engine's: `None` where the two are the same sequence of the same doubles.
    /// The engine's is `getZoneRefinedCRSVertices` in the 5x6 CRS, which is
    /// `getBaseRefinedVertices(false, n)` itself (`RI7H.ec:3887`).
    #[cfg(feature = "oracle")]
    fn planar_ring_difference(oracle: &str, id: u64, n: i32) -> Option<String> {
        use crate::Indexing as _;
        let zone = crate::ZoneId(id);
        let ours = HexA7::planar_refined_ring(&crate::indexings::Z7::decode(zone), false, n);
        let theirs = dggal_oracle::refined_crs_vertices_5x6(oracle, id, n);
        let text = crate::indexings::Z7::to_text(zone);
        if ours.len() != theirs.len() {
            return Some(format!(
                "{oracle} {text} ({id:#018x}) at {n} divisions: {} points, the engine has {}",
                ours.len(),
                theirs.len()
            ));
        }
        ours.iter()
            .zip(&theirs)
            .position(|(a, b)| a.0.to_bits() != b.0.to_bits() || a.1.to_bits() != b.1.to_bits())
            .map(|i| {
                format!(
                    "{oracle} {text} ({id:#018x}) at {n} divisions: point {i} of {} is {:?}, \
                     the engine has {:?}",
                    ours.len(),
                    ours[i],
                    theirs[i]
                )
            })
    }

    /// The planar ring of the refined boundary against the engine's own, before any
    /// projection: a stage oracle.
    ///
    /// `getZoneRefinedCRSVertices` in the 5x6 CRS returns `getBaseRefinedVertices(false, n)`
    /// as it stands (`RI7H.ec:3887`), so the walker, the even-level tables, the four polar arms
    /// and `add_intermediate_points` beneath them are compared here bit for bit, as sequences,
    /// with nothing of the inverse projection between this crate and the engine. The ring is
    /// the same on the three grids, which share the 5x6 plane; all three are asked, since each
    /// quantises the sampled points to zones of its own.
    ///
    /// The zones: every zone to level 3, taken from the engine's own lists of children; the
    /// zones that pin each compiled-order site and one zone of each kind, the degenerate ones
    /// among them (`PINNED_ZONES`); the six digit patterns of the test above under six base
    /// cells at every level to 19, which holds both polar pentagons at all twenty levels; and,
    /// at every level, the zones of both poles and of points beside them, of the antimeridian,
    /// of the meridian of the icosahedron's first vertex and of a fixed pseudo-random set.
    ///
    /// What the flag `crs84` changes (the fifth side of a polar pentagon, the interruption
    /// points at one division, the pole rule's first clause) this accessor never asks, since
    /// it passes false; that emission is compared through the projected ring, in the suites.
    #[test]
    #[cfg(feature = "oracle")]
    fn the_planar_refined_ring_is_the_engines_at_every_level_and_refinement() {
        use std::collections::BTreeSet;

        const BASES: [&str; 6] = ["00", "01", "04", "06", "09", "11"];
        const PATTERNS: [&str; 6] = [
            "0000000000000000000",
            "1313131313131313131",
            "3232323232323232323",
            "4141414141414141414",
            "6060606060606060606",
            "1456013456013456013",
        ];
        // Both poles and points beside them, the antimeridian, and the meridian of the
        // icosahedron's first vertex with its antimeridian, which together are the edge that
        // runs over the north pole: after the suites' `adversarial_points`.
        const SPECIAL: [(f64, f64); 24] = [
            (90.0, 0.0),
            (-90.0, 0.0),
            (89.9999999, 0.0),
            (-89.9999999, 0.0),
            (89.9999999, 11.2),
            (-89.9999999, -168.8),
            (89.99, 101.2),
            (-89.99, -78.8),
            (0.0, 180.0),
            (0.0, -180.0),
            (45.0, 180.0),
            (-45.0, -180.0),
            (-80.0, 11.2),
            (-45.0, 11.2),
            (-10.0, 11.2),
            (10.0, 11.2),
            (45.0, 11.2),
            (80.0, 11.2),
            (-80.0, -168.8),
            (-45.0, -168.8),
            (-10.0, -168.8),
            (10.0, -168.8),
            (45.0, -168.8),
            (80.0, -168.8),
        ];
        const RANDOM_POINTS_A_LEVEL: usize = 150;

        // (grid, zones compared, coordinates compared)
        let mut compared: Vec<(&str, usize, usize)> = Vec::new();
        for oracle in [
            dggal_oracle::IGEO7,
            dggal_oracle::IVEA7H,
            dggal_oracle::RTEA7H,
        ] {
            let mut zones: BTreeSet<u64> = BTreeSet::new();

            // Every zone to level 3: the twelve base cells and three generations of the
            // engine's own children.
            let mut generation: BTreeSet<u64> = (0..12)
                .map(|base| dggal_oracle::zone_from_text(oracle, &format!("{base:02}")))
                .collect();
            for _ in 0..3 {
                zones.extend(&generation);
                generation = generation
                    .iter()
                    .flat_map(|&zone| dggal_oracle::children(oracle, zone))
                    .collect();
            }
            zones.extend(&generation);
            assert!(!zones.contains(&dggal_oracle::NULL_ZONE));
            assert_eq!(
                zones.len(),
                12 + 72 + 492 + 3432,
                "{oracle}: the zones to level 3"
            );

            zones.extend(PINNED_ZONES.iter().flat_map(|(_, pinned)| pinned.iter()));

            for base in BASES {
                for pattern in PATTERNS {
                    for level in 0..=19 {
                        let text = format!("{base}{}", &pattern[..level]);
                        let id = dggal_oracle::zone_from_text(oracle, &text);
                        assert_ne!(id, dggal_oracle::NULL_ZONE, "{oracle} cannot read {text}");
                        zones.insert(id);
                    }
                }
            }

            // A fixed pseudo-random set: a linear congruential generator, whose high bits
            // give a latitude and a longitude in degrees.
            let mut state: u64 = 0x2026_1001;
            let mut uniform = || {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                (state >> 11) as f64 / (1u64 << 53) as f64
            };
            for level in 0..=19 {
                let random: Vec<(f64, f64)> = (0..RANDOM_POINTS_A_LEVEL)
                    .map(|_| (180.0 * uniform() - 90.0, 360.0 * uniform() - 180.0))
                    .collect();
                for &(lat, lon) in SPECIAL.iter().chain(&random) {
                    let id = dggal_oracle::zone_from_geo(oracle, lat, lon, level);
                    assert_ne!(
                        id,
                        dggal_oracle::NULL_ZONE,
                        "{oracle}: no zone at ({lat}, {lon}), level {level}"
                    );
                    zones.insert(id);
                }
            }

            let mut coordinates = 0usize;
            let mut levels = [0usize; 20];
            // The kinds, counted on the ring at one division: an edge across an interruption
            // adds its two points, and an edge under the pole rule is divided in twenty.
            let (mut across_an_interruption, mut under_the_pole_rule) = (0usize, 0usize);
            for &id in &zones {
                use crate::Indexing as _;
                let address = crate::indexings::Z7::decode(crate::ZoneId(id));
                levels[address.len()] += 1;
                for n in PLANAR_REFINEMENTS {
                    if let Some(difference) = planar_ring_difference(oracle, id, n) {
                        panic!("{difference}");
                    }
                    let points = HexA7::planar_refined_ring(&address, false, n).len();
                    coordinates += 2 * points;
                    if n == 1 {
                        across_an_interruption += usize::from(points > 6);
                        under_the_pole_rule += usize::from(points >= 20);
                    }
                }
            }

            println!(
                "{oracle}: {} zones, {coordinates} coordinates, by level {levels:?}, \
                 {across_an_interruption} across an interruption, {under_the_pole_rule} under \
                 the pole rule",
                zones.len()
            );
            assert!(
                levels[4..].iter().all(|&count| count >= 200),
                "{oracle}: zones by level {levels:?}"
            );
            assert!(
                across_an_interruption >= 700 && under_the_pole_rule >= 130,
                "{oracle}: {across_an_interruption} rings across an interruption, \
                 {under_the_pole_rule} under the pole rule"
            );
            compared.push((oracle, zones.len(), coordinates));
        }

        // Asserted exactly, so that the test cannot quietly compare less than it claims.
        assert_eq!(
            compared,
            [
                (dggal_oracle::IGEO7, 7298, 4_535_378),
                (dggal_oracle::IVEA7H, 7298, 4_535_576),
                (dggal_oracle::RTEA7H, 7297, 4_534_602),
            ]
        );
    }

    /// The pinning zones on aperture 7, by the site each pins: reverted
    /// alone to the eC's order, the site moves the planar ring of each of its zones at some
    /// division of `PLANAR_REFINEMENTS`. The last entry holds one zone of each kind, the
    /// degenerate ones among them, and the identifier the engine writes for a neighbour and
    /// cannot read back.
    #[cfg(feature = "oracle")]
    const PINNED_ZONES: [(&str, &[u64]); 11] = [
        // `addIntermediatePoints`, before an interruption: `0200`, `0704`, `0741`.
        (
            "AIP-1",
            &[
                0x203f_ffff_ffff_ffff,
                0x713f_ffff_ffff_ffff,
                0x787f_ffff_ffff_ffff,
            ],
        ),
        // After an interruption: `10`, `02323`, `06040`, `00453`, `01320`.
        (
            "AIP-2",
            &[
                0xafff_ffff_ffff_ffff,
                0x269f_ffff_ffff_ffff,
                0x6107_ffff_ffff_ffff,
                0x095f_ffff_ffff_ffff,
                0x1687_ffff_ffff_ffff,
            ],
        ),
        // A plain edge: `06`, `010`, `014`, `016`.
        (
            "AIP-3",
            &[
                0x6fff_ffff_ffff_ffff,
                0x11ff_ffff_ffff_ffff,
                0x19ff_ffff_ffff_ffff,
                0x1dff_ffff_ffff_ffff,
            ],
        ),
        // The even north polar arm: levels 0, 4 and 10; then 6, 12 and 18.
        (
            "EN-1",
            &[
                0x0fff_ffff_ffff_ffff,
                0x0000_ffff_ffff_ffff,
                0x0000_0000_3fff_ffff,
            ],
        ),
        (
            "EN-5",
            &[
                0x0000_03ff_ffff_ffff,
                0x0000_0000_00ff_ffff,
                0x0000_0000_0000_003f,
            ],
        ),
        // The even south polar arm: levels 2, 4, 10 and 12.
        (
            "ES-1",
            &[
                0xb03f_ffff_ffff_ffff,
                0xb000_ffff_ffff_ffff,
                0xb000_0000_3fff_ffff,
                0xb000_0000_00ff_ffff,
            ],
        ),
        // The odd north polar arm: levels 3, 9 and 15; 1, 9 and 15; 11, 13 and 19. ON-3
        // moves no ring, ON-2 alone moves those of levels 1, 9 and 15, and the zone of level
        // 11 is ON-4's.
        (
            "ON-1",
            &[
                0x0007_ffff_ffff_ffff,
                0x0000_0001_ffff_ffff,
                0x0000_0000_0000_7fff,
            ],
        ),
        (
            "ON-2",
            &[
                0x01ff_ffff_ffff_ffff,
                0x0000_0001_ffff_ffff,
                0x0000_0000_0000_7fff,
            ],
        ),
        (
            "ON-4",
            &[
                0x0000_0000_07ff_ffff,
                0x0000_0000_001f_ffff,
                0x0000_0000_0000_0007,
            ],
        ),
        // The odd south polar arm: levels 3, 5 and 9.
        (
            "OS-1",
            &[
                0xb007_ffff_ffff_ffff,
                0xb000_1fff_ffff_ffff,
                0xb000_0001_ffff_ffff,
            ],
        ),
        // `0313522`, `0100000`, `0000000`, `110000`, `0100033`, `0132336`, `0132323`,
        // `013232323232323232323`, `005151515151515151515`, `0000000000000000131`,
        // `01323232323232322` and `01222222222222222220`.
        (
            "one of each kind",
            &[
                0x32ea_5fff_ffff_ffff,
                0x1000_1fff_ffff_ffff,
                0x0000_1fff_ffff_ffff,
                0xb000_ffff_ffff_ffff,
                0x1003_7fff_ffff_ffff,
                0x169b_dfff_ffff_ffff,
                0x169a_7fff_ffff_ffff,
                0x169a_69a6_9a69_a69f,
                0x0a69_a69a_69a6_9a6f,
                0x0000_0000_0000_b3ff,
                0x169a_69a6_9a69_7fff,
                0x1492_4924_9249_243f,
            ],
        ),
    ];

    /// Each compiled-order site of the planar ring, at the zones that pin it, against the
    /// engine. Every difference is gathered before the test fails, so that a site put back in
    /// the eC's order names its zones in the failure.
    #[test]
    #[cfg(feature = "oracle")]
    fn the_compiled_order_sites_hold_at_the_zones_that_pin_them() {
        let mut differences: Vec<String> = Vec::new();
        for oracle in [
            dggal_oracle::IGEO7,
            dggal_oracle::IVEA7H,
            dggal_oracle::RTEA7H,
        ] {
            for (site, zones) in PINNED_ZONES {
                for &id in zones {
                    for n in PLANAR_REFINEMENTS {
                        if let Some(difference) = planar_ring_difference(oracle, id, n) {
                            differences.push(format!("{site}: {difference}"));
                        }
                    }
                }
            }
        }
        assert!(
            differences.is_empty(),
            "{} rings leave the engine's:\n{}",
            differences.len(),
            differences.join("\n")
        );
    }

    /// The engine's planar rings of six polar pentagons at one division, in the 5x6 CRS
    /// (`getZoneRefinedCRSVertices`, that is `getBaseRefinedVertices(false, 1)`), read from
    /// DGGAL v0.0.6, BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, and written out as the
    /// shortest decimals that read back to the same doubles; they are the same on the three
    /// grids. At one division such a ring is the arm's own points and nothing else: for each of
    /// the first four sides `b`, `d` and the next side's `ab`, then the fifth `b` and `d`, the
    /// two images of the pole, and the first `ab`.
    const POLAR_RINGS: [(&str, [(f64, f64); 17]); 6] = [
        // level 4, even, north: EN-1 and EN-5.
        (
            "000000",
            [
                (0.9931972789115646, 0.006802721088435375),
                (1.0, 0.010204081632653073),
                (1.989795918367347, 1.0),
                (1.9931972789115646, 1.0068027210884354),
                (2.0, 1.010204081632653),
                (2.989795918367347, 2.0),
                (2.993197278911565, 2.006802721088435),
                (3.0, 2.010204081632653),
                (3.989795918367347, 3.0),
                (3.993197278911565, 3.006802721088435),
                (4.0, 3.010204081632653),
                (4.989795918367347, 4.0),
                (4.993197278911564, 4.006802721088436),
                (5.0, 4.010204081632653),
                (5.0, 4.0),
                (1.0, 0.0),
                (0.9897959183673469, 0.0),
            ],
        ),
        // level 6, even, north: EN-5.
        (
            "00000000",
            [
                (0.9990281827016521, 0.0009718172983479106),
                (1.0, 0.001457725947521804),
                (1.998542274052478, 1.0),
                (1.999028182701652, 1.000971817298348),
                (2.0, 1.001457725947522),
                (2.998542274052478, 2.0),
                (2.999028182701652, 2.000971817298348),
                (3.0, 2.001457725947522),
                (3.998542274052478, 3.0),
                (3.999028182701652, 3.000971817298348),
                (4.0, 3.001457725947522),
                (4.998542274052478, 4.0),
                (4.999028182701652, 4.000971817298348),
                (5.0, 4.001457725947522),
                (5.0, 4.0),
                (1.0, 0.0),
                (0.9985422740524781, 0.0),
            ],
        ),
        // level 2, even, south: ES-1.
        (
            "1100",
            [
                (4.0476190476190474, 5.9523809523809526),
                (4.0, 5.928571428571429),
                (3.071428571428571, 5.0),
                (3.0476190476190474, 4.9523809523809526),
                (3.0, 4.928571428571429),
                (2.071428571428571, 4.0),
                (2.0476190476190474, 3.9523809523809526),
                (2.0, 3.928571428571429),
                (1.0714285714285712, 3.0),
                (1.0476190476190477, 2.9523809523809526),
                (1.0, 2.928571428571429),
                (0.07142857142857117, 2.0),
                (0.047619047619047616, 1.9523809523809523),
                (0.0, 1.9285714285714288),
                (0.0, 2.0),
                (4.0, 6.0),
                (4.071428571428571, 6.0),
            ],
        ),
        // level 9, odd, north: ON-1 and ON-2.
        (
            "00000000000",
            [
                (0.9999801669939112, 7.933202435493147e-5),
                (1.0, 8.329862557267805e-5),
                (1.9999167013744272, 1.0),
                (1.9999801669939112, 1.0000793320243548),
                (2.0, 1.0000832986255728),
                (2.9999167013744272, 2.0),
                (2.9999801669939115, 2.000079332024355),
                (3.0, 2.0000832986255728),
                (3.9999167013744272, 3.0),
                (3.9999801669939115, 3.000079332024355),
                (4.0, 3.0000832986255728),
                (4.999916701374428, 4.0),
                (4.999980166993911, 4.000079332024355),
                (5.0, 4.000083298625572),
                (5.0, 4.0),
                (1.0, 0.0),
                (0.9999167013744273, 0.0),
            ],
        ),
        // level 11, odd, north: ON-4.
        (
            "0000000000000",
            [
                (0.9999971667134159, 1.1333146336418782e-5),
                (1.0, 1.1899803653239721e-5),
                (1.9999881001963467, 1.0),
                (1.999997166713416, 1.0000113331463365),
                (2.0, 1.0000118998036533),
                (2.999988100196347, 2.0),
                (2.999997166713416, 2.0000113331463365),
                (3.0, 2.000011899803653),
                (3.999988100196347, 3.0),
                (3.999997166713416, 3.0000113331463365),
                (4.0, 3.000011899803653),
                (4.999988100196346, 4.0),
                (4.999997166713416, 4.000011333146336),
                (5.0, 4.000011899803654),
                (5.0, 4.0),
                (1.0, 0.0),
                (0.9999881001963468, 0.0),
            ],
        ),
        // level 3, odd, south: OS-1.
        (
            "11000",
            [
                (4.006802721088436, 5.9727891156462585),
                (4.0, 5.9714285714285715),
                (3.0285714285714285, 5.0),
                (3.006802721088435, 4.9727891156462585),
                (3.0, 4.9714285714285715),
                (2.0285714285714285, 4.0),
                (2.006802721088435, 3.9727891156462585),
                (2.0, 3.9714285714285715),
                (1.0285714285714285, 3.0),
                (1.0068027210884354, 2.9727891156462585),
                (1.0, 2.9714285714285715),
                (0.02857142857142847, 2.0),
                (0.006802721088435373, 1.9727891156462585),
                (0.0, 1.9714285714285715),
                (0.0, 2.0),
                (4.0, 6.0),
                (4.0285714285714285, 6.0),
            ],
        ),
    ];

    /// The four polar arms of [`get_base_refined_vertices`] in the order the library performs
    /// them, on the engine's own rings. Each of the seven sites of the arms that move a ring
    /// (EN-1, EN-5, ES-1, ON-1, ON-2, ON-4 and OS-1), put back in the eC's order, fails this
    /// test at the ring named for it in `POLAR_RINGS`; each was so put back, and seen to fail,
    /// when the arms were ported. The test of the pinned zones says the same against the live
    /// engine, at three zones or more a site.
    #[test]
    fn the_polar_arms_follow_the_compiled_order() {
        use crate::Indexing as _;
        let mut differences: Vec<String> = Vec::new();
        for (text, engine) in &POLAR_RINGS {
            let id = crate::indexings::Z7::from_text(text).expect("a canonical identifier");
            let ours = HexA7::planar_refined_ring(&crate::indexings::Z7::decode(id), false, 1);
            if ours.len() != engine.len() {
                differences.push(format!("{text}: {} points for 17", ours.len()));
                continue;
            }
            for (i, (a, b)) in ours.iter().zip(engine).enumerate() {
                if a.0.to_bits() != b.0.to_bits() || a.1.to_bits() != b.1.to_bits() {
                    differences.push(format!("{text}: point {i} is {a:?}, the engine has {b:?}"));
                }
            }
        }
        assert!(
            differences.is_empty(),
            "{} differences from the engine's rings:\n{}",
            differences.len(),
            differences.join("\n")
        );
    }

    /// Six zones deep under the polar pentagon of base cell 0, two at each of
    /// resolutions 17, 18 and 19, pinned to the engine's own doubles on each
    /// aperture-7 grid: the centroid and then the vertices, in the engine's order, as
    /// `getZoneWGS84Centroid` and `getZoneWGS84Vertices` answer them in degrees. Near
    /// the pole the geometric parent search can name a zone other than the one the
    /// walk came from, and the walk reaches these answers only by handing the rotation
    /// offset its own two previous zones, as the eC's does.
    ///
    /// Each zone is read from its text, and from its text alone: these six lie in the broken
    /// seam over the north pole, where the engine's hierarchy does not follow the digits.
    /// None of them is among the engine's children of the zone that its text less the last
    /// digit names, and so none is among this crate's; and the engine's `getZoneParents`
    /// answers for them, the same on the three grids, the lists of `PARENTS`, which are
    /// this crate's: for `0000000000000000136` the one zone `000000000000000005`, which the
    /// engine lists twice, and for `00000000000000001311` no parent at all.
    #[test]
    fn deep_polar_zones_answer_as_the_engine() {
        const PARENTS: [(&str, &[&str]); 6] = [
            (
                "0000000000000001644",
                &["000000000000000005", "000000000000000000"],
            ),
            ("0000000000000000136", &["000000000000000005"]),
            (
                "00000000000000016522",
                &["0000000000000000053", "0000000000000001601"],
            ),
            ("00000000000000001311", &[]),
            (
                "000000000000000013033",
                &["00000000000000000501", "00000000000000000512"],
            ),
            (
                "000000000000000000133",
                &["00000000000000000005", "00000000000000000166"],
            ),
        ];
        type Pinned = [(&'static str, [(f64, f64); 7]); 6];
        const IGEO7: Pinned = [
            (
                "0000000000000001644",
                [
                    (58.397148892705964, 11.200000000000001),
                    (58.39714714422544, 11.200003779577845),
                    (58.397149785150866, 11.20000536966907),
                    (58.397151160887354, 11.200001073933842),
                    (58.3971505022444, 11.199995704264643),
                    (58.39714791976935, 11.199995146488341),
                    (58.39714668297508, 11.199998926066272),
                ],
            ),
            (
                "0000000000000000136",
                [
                    (58.39715285998272, 11.200051208943597),
                    (58.39715029130645, 11.200050042651654),
                    (58.397150456323665, 11.200054485957805),
                    (58.39715302499988, 11.200055652250157),
                    (58.39715544672103, 11.200052481186624),
                    (58.397155299765906, 11.200048143830406),
                    (58.39715269496541, 11.200046765637165),
                ],
            ),
            (
                "00000000000000016522",
                [
                    (58.397155267628996, 11.20009660999695),
                    (58.39715451014794, 11.200095642011275),
                    (58.3971542139153, 11.200097379671922),
                    (58.39715497139635, 11.200098347657606),
                    (58.397156025110064, 11.200097577982689),
                    (58.39715632134269, 11.200095840321945),
                    (58.39715556386164, 11.200094872336239),
                ],
            ),
            (
                "00000000000000001311",
                [
                    (58.397147400068526, 11.200000000000003),
                    (58.39714668297509, 11.200001073933725),
                    (58.39714745851905, 11.200002147867496),
                    (58.39714817561247, 11.200001073933755),
                    (58.39714817561249, 11.199998926066224),
                    (58.397147458519065, 11.199997852132508),
                    (58.39714668297508, 11.199998926066293),
                ],
            ),
            (
                "000000000000000013033",
                [
                    (58.3971539239244, 11.200060428783281),
                    (58.39715355697068, 11.200060262169982),
                    (58.39715358054452, 11.200060896928102),
                    (58.39715394749825, 11.200061063541405),
                    (58.397154290878134, 11.20006059539658),
                    (58.39715426730431, 11.200059960638466),
                    (58.39715390035056, 11.200059794025162),
                ],
            ),
            (
                "000000000000000000133",
                [
                    (58.39714697137622, 11.200009219835607),
                    (58.397146604422424, 11.200009053223047),
                    (58.39714662799649, 11.200009687980943),
                    (58.397146994950305, 11.200009854593521),
                    (58.39714733833005, 11.20000938644816),
                    (58.397147314755955, 11.20000875169022),
                    (58.397146947802135, 11.20000858507765),
                ],
            ),
        ];
        const IVEA7H: Pinned = [
            (
                "0000000000000001644",
                [
                    (58.39714900515731, 11.200000000000001),
                    (58.397147136639106, 11.200003867709864),
                    (58.39714958069271, 11.200005086545648),
                    (58.39715132266484, 11.200001034580486),
                    (58.397150446849864, 11.199995891967394),
                    (58.39714788060774, 11.199995094849791),
                    (58.39714664208346, 11.19999898269095),
                ],
            ),
            (
                "0000000000000000136",
                [
                    (58.39715320981185, 11.200052881381522),
                    (58.397150653319386, 11.200051266734318),
                    (58.39715085193676, 11.200055776627805),
                    (58.39715340953324, 11.200057400917675),
                    (58.39715575671461, 11.200054425777099),
                    (58.397155554093324, 11.20004988836872),
                    (58.39715300976769, 11.200048359516526),
                ],
            ),
            (
                "00000000000000016522",
                [
                    (58.39715596029178, 11.200099137385932),
                    (58.397155200525084, 11.200098023083372),
                    (58.39715491983607, 11.200099719696787),
                    (58.397155679810155, 11.200100836162735),
                    (58.39715671974425, 11.200100248588035),
                    (58.397156999839076, 11.20009854611689),
                    (58.39715624061926, 11.200097437089935),
                ],
            ),
            (
                "00000000000000001311",
                [
                    (58.39714745629419, 11.2),
                    (58.39714664208346, 11.200001017309047),
                    (58.397147376735774, 11.200002034618132),
                    (58.397148217245665, 11.200001032948023),
                    (58.39714821724565, 11.199998967051975),
                    (58.397147376735745, 11.199997965381856),
                    (58.39714664208343, 11.199998982690955),
                ],
            ),
            (
                "000000000000000013033",
                [
                    (58.39715433863966, 11.200062373089201),
                    (58.39715397390715, 11.200062146472977),
                    (58.39715400241319, 11.200062791927873),
                    (58.39715436717328, 11.200063018752987),
                    (58.39715470320677, 11.200062598481582),
                    (58.39715467463741, 11.20006195255857),
                    (58.39715431010095, 11.200061727387567),
                ],
            ),
            (
                "000000000000000000133",
                [
                    (58.39714703612952, 11.200009490684032),
                    (58.397146670519504, 11.200009256409519),
                    (58.39714669887709, 11.200009900487856),
                    (58.39714706457941, 11.20001013568022),
                    (58.39714740073354, 11.200009716426145),
                    (58.39714737205291, 11.200009069761224),
                    (58.39714700765885, 11.20000884551337),
                ],
            ),
        ];
        const RTEA7H: Pinned = [
            (
                "0000000000000001644",
                [
                    (58.39714912772701, 11.2),
                    (58.39714712393803, 11.200003914207423),
                    (58.39714950215034, 11.200004977784204),
                    (58.39715145681878, 11.200000995556925),
                    (58.39715039335442, 11.199996017772525),
                    (58.397147842881836, 11.19999509023568),
                    (58.39714662637497, 11.19999900444326),
                ],
            ),
            (
                "0000000000000000136",
                [
                    (58.397153542155095, 11.200054449604016),
                    (58.3971509510183, 11.20005213702691),
                    (58.39715117396726, 11.200056666522938),
                    (58.39715376510397, 11.200058979100351),
                    (58.397156080060704, 11.200056449932259),
                    (58.397155803880516, 11.200051608186381),
                    (58.39715331920608, 11.200049920107856),
                ],
            ),
            (
                "00000000000000016522",
                [
                    (58.39715654381035, 11.200101052377338),
                    (58.39715577163617, 11.200099744569151),
                    (58.39715549702272, 11.200101355413391),
                    (58.39715626919686, 11.20010266322162),
                    (58.397157315984494, 11.200102360185602),
                    (58.39715759059797, 11.200100749341313),
                    (58.39715681842381, 11.200099441533059),
                ],
            ),
            (
                "00000000000000001311",
                [
                    (58.397147517579064, 11.200000000000003),
                    (58.39714662637498, 11.200000995556783),
                    (58.39714734531882, 11.200001991113574),
                    (58.397148236522916, 11.200000995556834),
                    (58.39714823652288, 11.199999004443205),
                    (58.3971473453188, 11.199998008886446),
                    (58.39714662637495, 11.199999004443248),
                ],
            ),
            (
                "000000000000000013033",
                [
                    (58.39715472837744, 11.200064169333222),
                    (58.39715435821508, 11.200063838964939),
                    (58.39715439006487, 11.200064486035844),
                    (58.39715476022724, 11.200064816404158),
                    (58.397155098539834, 11.200064499701524),
                    (58.397155066690026, 11.200063852630596),
                    (58.39715469652766, 11.200063522262296),
                ],
            ),
            (
                "000000000000000000133",
                [
                    (58.39714709365774, 11.200009719725982),
                    (58.39714672349523, 11.200009389358248),
                    (58.39714675534532, 11.200010036429022),
                    (58.39714712550783, 11.200010366796768),
                    (58.39714746382026, 11.200010050093697),
                    (58.397147431970176, 11.200009403022923),
                    (58.397147061807665, 11.200009072655195),
                ],
            ),
        ];
        let bits = |p: &(f64, f64)| (p.0.to_bits(), p.1.to_bits());
        let mut departures = Vec::new();
        for (name, pinned) in [("IGEO7", IGEO7), ("IVEA7H", IVEA7H), ("RTEA7H", RTEA7H)] {
            let grid = crate::registry::get_grid(name).unwrap();
            for ((text, engine), (parents_of, parents)) in pinned.into_iter().zip(PARENTS) {
                assert_eq!(text, parents_of, "the two tables name the same zones");
                let id = grid.zone_from_text(text).unwrap();
                let by_digits = grid.zone_from_text(&text[..text.len() - 1]).unwrap();
                assert!(
                    !grid.children(by_digits).contains(&id),
                    "{name} {text} is among the children of the zone one digit shorter, where \
                     the engine's list does not hold it"
                );
                let ours: Vec<String> = grid
                    .parents(id)
                    .into_iter()
                    .map(|p| grid.text_id(p))
                    .collect();
                assert_eq!(ours, parents, "{name}: the parents of {text}");
                let c = grid.centroid(id);
                let ours: Vec<(f64, f64)> = std::iter::once((c.lat, c.lon))
                    .chain(grid.vertices(id).iter().map(|v| (v.lat, v.lon)))
                    .collect();
                if !ours.iter().map(bits).eq(engine.iter().map(bits)) {
                    departures.push(format!(
                        "{name} {text}: ours {ours:?}, the engine's {engine:?}"
                    ));
                }
            }
        }
        assert!(
            departures.is_empty(),
            "{} of 18 zones depart from the engine:\n{}",
            departures.len(),
            departures.join("\n")
        );
    }

    /// The engine's own centroid of twelve zones deep under the polar pentagon of base cell 0,
    /// six at each of resolutions 17 and 19, quantised at the zone's resolution, pinned on each
    /// aperture-7 grid to the zone the engine names for it: the zone, its centroid as the engine
    /// answers it in degrees, and the engine's quantisation of that centroid. The engine does
    /// not name the zone itself there but another (`...053` for the centroid of `...122`), and
    /// so must this crate. Converting the cell the point falls in to its Z7 identifier, the
    /// geometric parent search names a parent whose primary children do not include that
    /// cell, once for each of these points, and the child position the conversion then uses
    /// is the one the eC answers for a child it does not find, the child count
    /// (RI7H_Z7.ec:143-146). Every literal was read from the engine through `dggal-oracle`.
    #[test]
    fn the_engines_centroid_quantises_as_the_engine_where_the_child_scan_misses() {
        type Pinned = [(&'static str, (f64, f64), &'static str); 12];
        const IGEO7: Pinned = [
            (
                "0000000000000000122",
                (58.39714739259525, 11.200039989748737),
                "0000000000000000053",
            ),
            (
                "0000000000000000123",
                (58.397150126289134, 11.200045599345536),
                "0000000000000000051",
            ),
            (
                "0000000000000000126",
                (58.39714963123623, 11.200032269427513),
                "0000000000000000052",
            ),
            (
                "0000000000000000132",
                (58.397150621340714, 11.200058929263998),
                "0000000000000000053",
            ),
            (
                "0000000000000000133",
                (58.3971533550338, 11.200064538863442),
                "0000000000000000051",
            ),
            (
                "0000000000000000136",
                (58.39715285998272, 11.200051208943597),
                "0000000000000000052",
            ),
            (
                "000000000000000000122",
                (58.397146119598126, 11.200005712821035),
                "000000000000000000053",
            ),
            (
                "000000000000000000123",
                (58.39714651012607, 11.200006514191472),
                "000000000000000000051",
            ),
            (
                "000000000000000000126",
                (58.39714643940375, 11.200004609917743),
                "000000000000000000052",
            ),
            (
                "000000000000000000132",
                (58.39714658084831, 11.20000841846512),
                "000000000000000000053",
            ),
            (
                "000000000000000000133",
                (58.39714697137622, 11.200009219835607),
                "000000000000000000051",
            ),
            (
                "000000000000000000136",
                (58.39714690065396, 11.20000731556186),
                "000000000000000000052",
            ),
        ];
        const IVEA7H: Pinned = [
            (
                "0000000000000000122",
                (58.39714769244801, 11.200040554962078),
                "0000000000000000053",
            ),
            (
                "0000000000000000123",
                (58.39715045463654, 11.200046756136695),
                "0000000000000000051",
            ),
            (
                "0000000000000000126",
                (58.39714985787237, 11.200033217396685),
                "0000000000000000052",
            ),
            (
                "0000000000000000132",
                (58.39715105050506, 11.200060285970732),
                "0000000000000000053",
            ),
            (
                "0000000000000000133",
                (58.397153808305056, 11.20006643480083),
                "0000000000000000051",
            ),
            (
                "0000000000000000136",
                (58.39715320981185, 11.200052881381522),
                "0000000000000000052",
            ),
            (
                "000000000000000000122",
                (58.39714616243432, 11.200005793565763),
                "000000000000000000053",
            ),
            (
                "000000000000000000123",
                (58.397146557032954, 11.200006679447409),
                "000000000000000000051",
            ),
            (
                "000000000000000000126",
                (58.397146471780424, 11.200004745341975),
                "000000000000000000052",
            ),
            (
                "000000000000000000132",
                (58.397146642157686, 11.200008612280435),
                "000000000000000000053",
            ),
            (
                "000000000000000000133",
                (58.39714703612952, 11.200009490684032),
                "000000000000000000051",
            ),
            (
                "000000000000000000136",
                (58.397146950629725, 11.200007554481752),
                "000000000000000000052",
            ),
        ];
        const RTEA7H: Pinned = [
            (
                "0000000000000000122",
                (58.397147913982955, 11.200040765458855),
                "0000000000000000053",
            ),
            (
                "0000000000000000123",
                (58.39715072806923, 11.200047607530921),
                "0000000000000000051",
            ),
            (
                "0000000000000000126",
                (58.397150059220856, 11.200034019043217),
                "0000000000000000052",
            ),
            (
                "0000000000000000132",
                (58.39715139691604, 11.200061196019005),
                "0000000000000000053",
            ),
            (
                "0000000000000000133",
                (58.397154211001215, 11.200068038093134),
                "0000000000000000051",
            ),
            (
                "0000000000000000136",
                (58.397153542155095, 11.200054449604016),
                "0000000000000000052",
            ),
            (
                "000000000000000000122",
                (58.397146194082225, 11.200005823636767),
                "000000000000000000053",
            ),
            (
                "000000000000000000123",
                (58.397146596094856, 11.200006801075252),
                "000000000000000000051",
            ),
            (
                "000000000000000000126",
                (58.397146500544544, 11.200004859862952),
                "000000000000000000052",
            ),
            (
                "000000000000000000132",
                (58.39714669164512, 11.200008742287455),
                "000000000000000000053",
            ),
            (
                "000000000000000000133",
                (58.39714709365774, 11.200009719725982),
                "000000000000000000051",
            ),
            (
                "000000000000000000136",
                (58.397146998107466, 11.200007778513667),
                "000000000000000000052",
            ),
        ];
        let mut departures = Vec::new();
        for (name, pinned) in [("IGEO7", IGEO7), ("IVEA7H", IVEA7H), ("RTEA7H", RTEA7H)] {
            let grid = crate::registry::get_grid(name).unwrap();
            for (text, (lat, lon), engine) in pinned {
                let res = (text.len() - 2) as u8;
                let ours = grid
                    .zone_from_geo(lat, lon, res)
                    .map(|id| grid.text_id(id))
                    .unwrap_or_else(|e| format!("{e:?}"));
                if ours != engine {
                    departures.push(format!(
                        "{name}: the centroid of {text} quantises to {ours}, the engine's to {engine}"
                    ));
                }
            }
        }
        assert!(
            departures.is_empty(),
            "{} of 36 quantisations depart from the engine:\n{}",
            departures.len(),
            departures.join("\n")
        );
    }

    /// The neighbours of nine zones, pinned on each aperture-7 grid to the engine's own list,
    /// in its order, less the four kinds of entry this crate never hands out: a hexagon; a
    /// pentagon of each kind, the polar one of base cell 0, whose vertices are a fan, and that
    /// of base cell 2, whose boundary is traced; the zones about both poles at resolution 13;
    /// a zone whose neighbours the engine finds by carrying its centroid across a rhombus
    /// interruption; and three zones in the broken seam deep under the pentagon of base cell
    /// 0. The engine's list for the last three holds entries dropped here: for
    /// `0000000000000000043`, `[...005, ...004, ...000, ...000, (null), ...005]`, its null zone
    /// and two repeats; for `000000000000000014`, `[...016, ...014, ...051, ...041, ...041,
    /// ...041]`, the zone itself and two repeats; and for `000000000000000104`, `[...100,
    /// ...105, ...152, 012222222222222220, 012222222222222220, ...106]`, twice an identifier the
    /// engine cannot read back. The lists are the same on the three grids, which share the
    /// topology and differ only in the projection, and each was read from each grid's engine
    /// through `dggal-oracle`.
    #[test]
    fn the_neighbours_are_the_engines_lists_in_order() {
        const PINNED: [(&str, &[&str]); 9] = [
            (
                "0064156",
                &[
                    "0064150", "0064154", "0064141", "0064143", "0064105", "0064152",
                ],
            ),
            (
                "0000000",
                &["0000005", "0000003", "0000004", "0000001", "0000006"],
            ),
            (
                "0200000",
                &["0200003", "0200001", "0200005", "0200004", "0200006"],
            ),
            (
                "013232323232323",
                &[
                    "005151515151511",
                    "013232323232336",
                    "013232323232321",
                    "013232323232320",
                    "013232323232322",
                    "005151515151515",
                ],
            ),
            (
                "084545454545454",
                &[
                    "084545454545450",
                    "084545454545455",
                    "114545454545454",
                    "114545454545455",
                    "084545454545441",
                    "084545454545456",
                ],
            ),
            (
                "02320400360",
                &[
                    "02320400364",
                    "02320400366",
                    "02320400362",
                    "02320400363",
                    "02320400361",
                    "02320400365",
                ],
            ),
            (
                "0000000000000000043",
                &[
                    "0000000000000000005",
                    "0000000000000000004",
                    "0000000000000000000",
                ],
            ),
            (
                "000000000000000014",
                &[
                    "000000000000000016",
                    "000000000000000051",
                    "000000000000000041",
                ],
            ),
            (
                "000000000000000104",
                &[
                    "000000000000000100",
                    "000000000000000105",
                    "000000000000000152",
                    "000000000000000106",
                ],
            ),
        ];
        let mut departures = Vec::new();
        for name in ["IGEO7", "IVEA7H", "RTEA7H"] {
            let grid = crate::registry::get_grid(name).unwrap();
            for (text, engine) in PINNED {
                let id = grid.zone_from_text(text).unwrap();
                let ours: Vec<String> = grid
                    .neighbors(id)
                    .into_iter()
                    .map(|n| grid.text_id(n))
                    .collect();
                if ours != engine {
                    departures.push(format!(
                        "{name} {text}: ours {ours:?}, the engine's {engine:?}"
                    ));
                }
            }
        }
        assert!(
            departures.is_empty(),
            "{} of 27 lists depart from the engine's:\n{}",
            departures.len(),
            departures.join("\n")
        );
    }

    /// `quantize` answers the null zone, and does not panic, for a planar point with a
    /// coordinate that is not finite or lies a billion units or more from the origin, at
    /// every resolution the packing holds; the thirteen values below, taken in every pair,
    /// include such coordinates and others inside the layout. `Grid` never passes such a
    /// point, since it refuses a coordinate that is not finite and every projection places a
    /// finite one within the layout, but this is a public method of the trait, and a debug
    /// build overflowed an integer on those points.
    #[test]
    fn a_point_far_outside_the_layout_is_the_null_zone() {
        let values = [
            f64::NAN,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::MAX,
            f64::MIN,
            -0.0,
            0.0,
            1e-300,
            5.0,
            6.0,
            1e9,
            -1e9,
            1e18,
        ];
        for x in values {
            for y in values {
                let far = !(x.abs() < 1e9 && y.abs() < 1e9);
                for res in 0..=20 {
                    let p = PlanarPoint { face: -1, x, y };
                    let got = HexA7::quantize(p, res);
                    assert!(!far || got.is_none(), "({x}, {y}) at {res}: {got:?}");
                }
            }
        }
        // The bound lies no lower than it must. Inside the layout a point has a cell at every
        // resolution; and just under the bound, the point `(-9.99e8, 0.5)` is still answered
        // at resolution 0, with the meaningless base cell 11, rather than refused, so that a
        // bound lowered anywhere below it fails here.
        for res in 0..=19 {
            let p = PlanarPoint {
                face: -1,
                x: 1.3,
                y: 2.2,
            };
            assert!(HexA7::quantize(p, res).is_some(), "(1.3, 2.2) at {res}");
        }
        let p = PlanarPoint {
            face: -1,
            x: -9.99e8,
            y: 0.5,
        };
        assert_eq!(HexA7::quantize(p, 0), Some(Address::new(11, &[])));
    }

    /// `neighbors` is total on a direct call through the trait, as `HexA3`'s is: an address
    /// this topology has no cell for (a base cell of 12 or 15, a digit of 9) and one of the
    /// twentieth level, which has no geometry, answer no neighbours, as `Grid` answers them
    /// and as the engine answers a zone of the twentieth level.
    #[test]
    fn neighbours_are_total_on_a_direct_call() {
        for a in [
            Address::new(3, &[1; 20]),
            Address::new(12, &[]),
            Address::new(15, &[1, 2, 3]),
            Address::new(4, &[1, 9, 2]),
        ] {
            let got = std::panic::catch_unwind(|| HexA7::neighbors(&a));
            assert!(
                matches!(&got, Ok(Some(ns)) if ns.is_empty()),
                "{a:?}: {:?}",
                got.map_err(|_| "a panic")
            );
        }
    }

    /// The parents and the children of thirteen zones, pinned on each aperture-7 grid to the
    /// engine's own lists, `getZoneParents` and `getZoneChildren`, in its order, less the three
    /// kinds of entry this crate never hands out: its null zone, an identifier it cannot read
    /// back, and a repeat. The lists are the same on the three grids, which share the topology
    /// and differ only in the projection, and each was read from each grid's engine.
    ///
    /// They are a hexagon, with two parents and thirteen children, the first seven of which
    /// are the zones its text names with one digit more, in the engine's order and not in the
    /// digits'; the polar pentagons of base cells 0 and 11 at level 0, with no parent and
    /// eleven children, and the second at level 1; the pentagon of base cell 3 at level 2; a
    /// centroid child, which has one parent, and one of its children; a zone of level 19,
    /// which has no children; and five zones in the broken seam over the north pole, from
    /// whose lists entries are left out. There the engine's lists are:
    ///
    /// - the children of `00052626050026015`: `[...150, ...152, ...153, ...151, ...155,
    ///   012222222222222220, 012222222222222220, 000526260500260105, 000526260500260114,
    ///   000526260500261626, 012222222222222220, 000526260500260533, 012222222222222220]`,
    ///   four times an identifier the engine cannot read back, two of them among its primary
    ///   children;
    /// - the children of `0000000000000005`: `[...050, ...052, ...053, ...051, ...055, ...054,
    ///   ...056, (null), (null), 00000000000000005, (null), (null), 00000000000000532]`, four
    ///   times its null zone; and the child `00000000000000005` that is kept does not name
    ///   the zone among its own parents, which are the one zone `0000000000000000`;
    /// - the children of `00000000000000005`: `[...050, ...053, ...051, ...055, ...054,
    ///   ...056, ...052, ...055, ...056, 000000000000000004, 000000000000000001,
    ///   000000000000000010, 000000000000000522]`, two repeats;
    /// - the children of the pentagon `0000000000000000`: `[...000, ...005, ...004, ...006,
    ///   ...003, ...001, (null), 00000000000000016, 00000000000000061, 00000000000000052,
    ///   00000000000000034]`, its null zone once;
    /// - the parents of `0000000000000000136`: `[000000000000000005, 000000000000000005]`, a
    ///   repeat, and its children: `[...520, ...522, ...523, ...521, ...525,
    ///   00000000000000000005, 00000000000000000005, 00000000000000000515,
    ///   00000000000000003555, 00000000000000000536, 00000000000000000005,
    ///   01222222222222222220, 00000000000000000005]`, three repeats and an identifier the
    ///   engine cannot read back.
    #[test]
    fn the_parents_and_children_are_the_engines_lists_in_order() {
        const PINNED: [(&str, &[&str], &[&str]); 13] = [
            (
                "006415654636",
                &["00641565463", "00641565462"],
                &[
                    "0064156546360",
                    "0064156546361",
                    "0064156546365",
                    "0064156546364",
                    "0064156546366",
                    "0064156546362",
                    "0064156546363",
                    "0064156546342",
                    "0064156546033",
                    "0064156546251",
                    "0064156546215",
                    "0064156546324",
                    "0064156546306",
                ],
            ),
            (
                "00",
                &[],
                &[
                    "000", "005", "004", "006", "003", "001", "023", "053", "033", "013", "043",
                ],
            ),
            (
                "11",
                &[],
                &[
                    "110", "112", "113", "111", "114", "116", "104", "074", "094", "064", "084",
                ],
            ),
            (
                "110",
                &["11"],
                &[
                    "1100", "1103", "1101", "1104", "1106", "1102", "1125", "1143", "1134", "1161",
                    "1116",
                ],
            ),
            (
                "0300",
                &["030"],
                &[
                    "03000", "03003", "03001", "03005", "03004", "03006", "03016", "03052",
                    "03043", "03061", "03034",
                ],
            ),
            (
                "0064156546360",
                &["006415654636"],
                &[
                    "00641565463600",
                    "00641565463601",
                    "00641565463605",
                    "00641565463604",
                    "00641565463606",
                    "00641565463602",
                    "00641565463603",
                    "00641565463616",
                    "00641565463652",
                    "00641565463643",
                    "00641565463661",
                    "00641565463625",
                    "00641565463634",
                ],
            ),
            (
                "00641565463601",
                &["0064156546360", "0064156546361"],
                &[
                    "006415654636010",
                    "006415654636011",
                    "006415654636015",
                    "006415654636014",
                    "006415654636016",
                    "006415654636012",
                    "006415654636013",
                    "006415654636162",
                    "006415654636053",
                    "006415654636001",
                    "006415654636035",
                    "006415654636344",
                    "006415654636126",
                ],
            ),
            (
                "006415654636011111111",
                &["00641565463601111111", "00641565463612424242"],
                &[],
            ),
            (
                "00052626050026015",
                &["0005262605002601", "0005262605002616"],
                &[
                    "000526260500260150",
                    "000526260500260152",
                    "000526260500260153",
                    "000526260500260151",
                    "000526260500260155",
                    "000526260500260105",
                    "000526260500260114",
                    "000526260500261626",
                    "000526260500260533",
                ],
            ),
            (
                "0000000000000005",
                &["000000000000000", "000000000000005"],
                &[
                    "00000000000000050",
                    "00000000000000052",
                    "00000000000000053",
                    "00000000000000051",
                    "00000000000000055",
                    "00000000000000054",
                    "00000000000000056",
                    "00000000000000005",
                    "00000000000000532",
                ],
            ),
            (
                "00000000000000005",
                &["0000000000000000"],
                &[
                    "000000000000000050",
                    "000000000000000053",
                    "000000000000000051",
                    "000000000000000055",
                    "000000000000000054",
                    "000000000000000056",
                    "000000000000000052",
                    "000000000000000004",
                    "000000000000000001",
                    "000000000000000010",
                    "000000000000000522",
                ],
            ),
            (
                "0000000000000000",
                &["000000000000000"],
                &[
                    "00000000000000000",
                    "00000000000000005",
                    "00000000000000004",
                    "00000000000000006",
                    "00000000000000003",
                    "00000000000000001",
                    "00000000000000016",
                    "00000000000000061",
                    "00000000000000052",
                    "00000000000000034",
                ],
            ),
            (
                "0000000000000000136",
                &["000000000000000005"],
                &[
                    "00000000000000000520",
                    "00000000000000000522",
                    "00000000000000000523",
                    "00000000000000000521",
                    "00000000000000000525",
                    "00000000000000000005",
                    "00000000000000000515",
                    "00000000000000003555",
                    "00000000000000000536",
                ],
            ),
        ];
        let mut departures = Vec::new();
        for name in ["IGEO7", "IVEA7H", "RTEA7H"] {
            let grid = crate::registry::get_grid(name).unwrap();
            let texts = |ids: Vec<crate::ZoneId>| -> Vec<String> {
                ids.into_iter().map(|z| grid.text_id(z)).collect()
            };
            for (text, parents, children) in PINNED {
                let id = grid.zone_from_text(text).unwrap();
                let ours = texts(grid.parents(id));
                if ours != parents {
                    departures.push(format!(
                        "{name}, the parents of {text}: ours {ours:?}, the engine's {parents:?}"
                    ));
                }
                let ours = texts(grid.children(id));
                if ours != children {
                    departures.push(format!(
                        "{name}, the children of {text}: ours {ours:?}, the engine's {children:?}"
                    ));
                }
            }
        }
        assert!(
            departures.is_empty(),
            "{} of 78 lists depart from the engine's:\n{}",
            departures.len(),
            departures.join("\n")
        );
    }

    /// `get_parents` and `get_children` answer the engine's lists before anything is left
    /// out of them, entry for entry: `I7HZone::getParents` and `getChildren` as the engine's
    /// `getZoneParents` and `getZoneChildren` then convert them (`RI7H_Z7.ec:550-564`). Each
    /// entry is shown as the engine prints it, save an identifier it cannot read back, for
    /// which `from_7h` answers nothing here: the engine prints `012222222222222220` for the
    /// four of the second list, which was read, as the others, from the engine.
    #[test]
    fn the_lists_before_anything_is_left_out_are_the_engines() {
        use crate::Indexing;
        use crate::indexings::Z7;

        const UNREADABLE: &str = "an identifier the engine cannot read back";
        fn zone(text: &str) -> Z {
            let digits: Vec<u8> = text.bytes().skip(2).map(|b| b - b'0').collect();
            zone_from_steps(text[..2].parse().unwrap(), &digits)
        }
        fn shown(z: Option<Z>, res: usize) -> String {
            match z.map(|z| from_7h(&z)) {
                None => crate::NULL_TEXT.to_string(),
                Some(None) => UNREADABLE.to_string(),
                Some(Some(z7)) => Z7::to_text(Z7::encode(&z7_address(z7, res as i64))),
            }
        }
        let children = |text: &str| -> Vec<String> {
            get_children(&zone(text))
                .into_iter()
                .map(|c| shown(c, text.len() - 1))
                .collect()
        };

        assert_eq!(
            children("0000000000000005"),
            [
                "00000000000000050",
                "00000000000000052",
                "00000000000000053",
                "00000000000000051",
                "00000000000000055",
                "00000000000000054",
                "00000000000000056",
                "(null)",
                "(null)",
                "00000000000000005",
                "(null)",
                "(null)",
                "00000000000000532",
            ]
        );
        assert_eq!(
            children("00052626050026015"),
            [
                "000526260500260150",
                "000526260500260152",
                "000526260500260153",
                "000526260500260151",
                "000526260500260155",
                UNREADABLE,
                UNREADABLE,
                "000526260500260105",
                "000526260500260114",
                "000526260500261626",
                UNREADABLE,
                "000526260500260533",
                UNREADABLE,
            ]
        );
        assert_eq!(
            children("00000000000000005"),
            [
                "000000000000000050",
                "000000000000000053",
                "000000000000000051",
                "000000000000000055",
                "000000000000000054",
                "000000000000000056",
                "000000000000000052",
                "000000000000000055",
                "000000000000000056",
                "000000000000000004",
                "000000000000000001",
                "000000000000000010",
                "000000000000000522",
            ]
        );
        // A zone of level 19 has no primary children, and so no children at all.
        assert!(children("006415654636011111111").is_empty());

        // The two parents the engine lists for this zone are two zones to `getParents`, which
        // the conversion gives one identifier.
        let text = "0000000000000000136";
        let parents = get_parents(&zone(text));
        assert_eq!(parents.len(), 2);
        assert_ne!(parents[0], parents[1]);
        let shown_parents: Vec<String> = parents
            .into_iter()
            .map(|p| shown(Some(p), text.len() - 3))
            .collect();
        assert_eq!(shown_parents, ["000000000000000005", "000000000000000005"]);
        // A centroid child has its one parent, and a zone of level 0 none.
        assert_eq!(get_parents(&zone("0064156546360")).len(), 1);
        assert!(get_parents(&zone("00")).is_empty());
    }

    /// The factor by which `getParents` shortens the step towards a vertex, as the library
    /// holds it at `0x4a120`, where the two products at `0x24710` and `0x24720` read it; the
    /// three by which `getChildren` lengthens it is the 3.0 at `0x49b70`. They are the
    /// constants that `get_parents` and `get_children` multiply by.
    #[test]
    fn the_hierarchy_factors_match_the_engines_folded_constants() {
        assert_eq!(PARENT_STEP_FACTOR.to_bits(), 0x3fef_ae14_7ae1_47ae); // 0x4a120
        assert_eq!(CHILD_STEP_FACTOR.to_bits(), 0x4008_0000_0000_0000); // 0x49b70
    }

    /// A zone is a centroid child where the last digit of its identifier is 0, as
    /// `Z7Zone::isCentroidChild` reads it (`RI7H_Z7.ec:81-103`), and a zone of level 0, which
    /// has no digit, is none. The engine's `isZoneCentroidChild` answers so for each
    /// identifier below on the three grids, the two of twenty digits among them, which it
    /// reads by their twentieth digit although it gives them no parent.
    ///
    /// The centroid parent is the first of the zone's parents that is itself a centroid
    /// child. The engine's own `getZoneCentroidParent` on these grids answers its null zone
    /// for every zone, which is not followed; each answer below is the first of the engine's
    /// parents for which its `isZoneCentroidChild` holds: the first parent of
    /// `00641565463601`, the second of `00641565463616`, whose parents are `0064156546361`
    /// and `0064156546360`, and none for `00641565463611`, whose parents are `0064156546361`
    /// and `0064156546304`, for the centroid child `0064156546360`, whose one parent
    /// `006415654636` is none, and at level 0.
    #[test]
    fn the_centroid_child_and_the_centroid_parent_are_read_from_the_last_digit() {
        use crate::ZoneId;
        for name in ["IGEO7", "IVEA7H", "RTEA7H"] {
            let grid = crate::registry::get_grid(name).unwrap();
            let id = |text: &str| grid.zone_from_text(text).unwrap();
            for (text, is) in [
                ("0064156546360", true),
                ("006415654636", false),
                ("0064100", true),
                ("0064156", false),
                ("110", true),
                ("11", false),
                ("00", false),
            ] {
                assert_eq!(grid.is_centroid_child(id(text)), is, "{name} {text}");
            }
            // `0064156546360111111110` and `0064156546360111111111`.
            assert!(
                grid.is_centroid_child(ZoneId(0x0d0d_d667_8124_9248)),
                "{name}"
            );
            assert!(
                !grid.is_centroid_child(ZoneId(0x0d0d_d667_8124_9249)),
                "{name}"
            );

            for (text, parent) in [
                ("00641565463601", Some("0064156546360")),
                ("00641565463616", Some("0064156546360")),
                ("0064100", Some("006410")),
                ("00641565463611", None),
                ("0064156546360", None),
                ("0064156", None),
                ("00", None),
            ] {
                assert_eq!(
                    grid.centroid_parent(id(text)).map(|p| grid.text_id(p)),
                    parent.map(str::to_string),
                    "{name} {text}"
                );
            }
        }
    }

    /// The four methods of the hierarchy are total on a direct call through the trait, and
    /// each answers for itself: never `None`, on which `Grid` would fall back to the digits
    /// of the identifier. An address this topology has no cell for (a base cell of 12 or 15,
    /// a digit of 9) and one of the twentieth level, which has no geometry, have no parents,
    /// no children and no centroid parent, as the engine answers a zone of the twentieth
    /// level; and a zone of level 0 has no parent and is no centroid child.
    #[test]
    fn the_hierarchy_is_total_on_a_direct_call() {
        for a in [
            Address::new(3, &[1; 20]),
            Address::new(12, &[]),
            Address::new(15, &[1, 2, 3]),
            Address::new(4, &[1, 9, 2]),
        ] {
            let got = std::panic::catch_unwind(|| {
                (
                    HexA7::parents(&a),
                    HexA7::children(&a),
                    HexA7::centroid_parent(&a),
                    HexA7::is_centroid_child(&a),
                )
            });
            assert!(
                matches!(
                    &got,
                    Ok((Some(ps), Some(cs), Some(None), Some(_))) if ps.is_empty() && cs.is_empty()
                ),
                "{a:?}: {:?}",
                got.map_err(|_| "a panic")
            );
        }
        for a in [Address::new(12, &[0]), Address::new(4, &[1, 9, 0])] {
            assert_eq!(HexA7::is_centroid_child(&a), Some(false), "{a:?}");
        }
        let base = Address::new(0, &[]);
        assert_eq!(HexA7::parents(&base), Some(Vec::new()));
        assert_eq!(HexA7::centroid_parent(&base), Some(None));
        assert_eq!(HexA7::is_centroid_child(&base), Some(false));
        assert_eq!(HexA7::children(&base).map(|cs| cs.len()), Some(11));
    }

    /// The compaction of sets of several levels, which `Grid` refuses and the engine answers,
    /// each read from the engine's `compactZones` on ISEA7H_Z7. A zone of level 0 has no
    /// parent, and leaves the set as soon as a pass meets it, which the first three sets show.
    ///
    /// The fourth reaches the short cut that follows the loop: of the seventy-two zones of
    /// level 1, those whose digit is 1, 3 or 4 as they are, and the children of the others.
    /// No zone of level 0 has all its children in either part, so the seventy-two all stay,
    /// with zones of level 2 beside them; and the short cut then puts in their place not the
    /// twelve zones of level 0 but the zone at the origin of the first rhombus at twelve even
    /// levels, ten of which the grid has. The engine's answer for that set is those ten
    /// identifiers and its null zone twice.
    ///
    /// The fifth shows that the loop goes on after a pass that takes nothing: the children of
    /// `0064` but the one at its centre, `00640`, with the thirteen children of that one. The
    /// first pass takes `00640`, the second takes nothing and ends the ascent, and the third,
    /// which then finds every child of `0064` among the zones kept, takes `0064`. Were the
    /// loop left at the second pass, the answer would lack `0064`.
    #[test]
    fn the_compaction_of_several_levels_is_the_engines() {
        use crate::Indexing;
        use crate::indexings::Z7;

        fn address(text: &str) -> Address {
            let digits: Vec<u8> = text.bytes().skip(2).map(|b| b - b'0').collect();
            Address::new(text[..2].parse().unwrap(), &digits)
        }
        fn compacted(set: &[Address]) -> String {
            let answer = HexA7::compact(Private, set).unwrap();
            let texts: Vec<String> = answer.iter().map(|a| Z7::to_text(Z7::encode(a))).collect();
            texts.join(" ")
        }
        let of = |texts: &str| -> Vec<Address> { texts.split(' ').map(address).collect() };

        for (set, answer) in [
            ("00 0064", "0064"),
            ("05 011", "011"),
            ("0064 00641 00", "0064 00641"),
            (
                "00641 00645 00644 00646 00642 00643 00665 00604 00656 00422 03332 03323 \
                 006400 006401 006405 006404 006406 006402 006403 006461 006425 006434 006416 \
                 006452 006443",
                "00422 03332 03323 0064 00640 00641 00645 00644 00646 00642 00643 00656 00665 \
                 00604 006416 006452 006405 006404 006443 006406 006401 006461 006402 006403 \
                 006434 006425",
            ),
        ] {
            assert_eq!(compacted(&of(set)), answer, "{set}");
        }

        let mut set: Vec<Z> = Vec::new();
        for base in 0..12u64 {
            for digit in 0..7u8 {
                let level_1 = Address::new(base, &[digit]);
                if names_deleted_child(&level_1) {
                    continue;
                }
                let zone = zone_from_steps(base, &[digit]);
                if matches!(digit, 1 | 3 | 4) {
                    set.push(zone);
                } else {
                    set.extend(get_children(&zone).into_iter().flatten());
                }
            }
        }
        let level_1 = set.iter().filter(|z| zone_level(z) == 1).count();
        assert_eq!((level_1, set.len()), (36, 36 + 36 * 13 - 2 * 12));
        let answer = compact(&set);
        let built: Vec<Z> = (0..12).map(|r| Z::new(r, 0, 0, 0, 0)).collect();
        assert_eq!(answer, built);
        let set: Vec<Address> = set
            .iter()
            .map(|z| z7_address(from_7h(z).unwrap(), zone_level(z)))
            .collect();
        assert_eq!(
            compacted(&set),
            "01 0100 010000 01000000 0100000000 010000000000 01000000000000 0100000000000000 \
             010000000000000000 01000000000000000000"
        );

        // A zone for which the engine finds no parent leaves the set, with nothing in its
        // place. This one, of level 18 in the broken seam over the pole of base cell 0, is
        // asked of the compaction itself: its identifier is one the engine cannot read back,
        // so that the conversion of the answer would leave it out whatever the compaction did.
        let orphan = zone_from_steps(0, address("00000000000000001311").digits());
        assert_eq!(zone_level(&orphan), 18);
        assert!(get_parents(&orphan).is_empty());
        assert!(from_7h(&orphan).is_none());
        assert_eq!(compact(&[orphan]), []);

        // No zone gives no zone; a zone given twice counts once; an address of twenty digits,
        // which the engine reads as its null zone, is left out of the set.
        assert_eq!(HexA7::compact(Private, &[]), Some(Vec::new()));
        assert_eq!(compacted(&of("0064 0064 00")), "0064");
        assert_eq!(
            compacted(&[address("0064"), Address::new(0, &[0; 20])]),
            "0064"
        );
        assert_eq!(compacted(&[Address::new(3, &[1; 20])]), "");
    }
}
