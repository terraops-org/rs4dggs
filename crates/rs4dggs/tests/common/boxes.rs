//! Boxes for the suites that put a bounding box to the engine's `listZones`, and the test of which
//! of them the engine can be asked.
//!
//! A box is made in degrees, `[ll.lat, ll.lon, ur.lat, ur.lon]`, as a caller gives it, and is
//! handed to the engine in radians by the product with `Pi / 180` that this crate makes everywhere
//! ([`radians`]): degrees to radians and back is not the same double, and a box made in radians
//! could part the crate from the engine at a tie. Nine kinds of box are made, each seeded and
//! so repeatable; each is a way in which a box meets the structure of the grid.
//!
//! The engine does not survive every box, and for some it costs more than a test may spend. The
//! suites ask [`engine_survives`] before they call it, and count the boxes it refuses.
use std::collections::BTreeMap;
use std::f64::consts::PI;
use std::sync::Mutex;

use dggal_oracle::{
    IGEO7, ISEA3H, IVEA3H, IVEA7H, NULL_ZONE, RTEA3H, RTEA7H, centroid, list_zones, max_level,
    vertices, zone_from_geo,
};

/// The nine kinds of box.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    /// Uniform centre; sides from a tenth of a zone to sixteen zones.
    Rand,
    /// Sides from a millionth to a tenth of a zone.
    Tiny,
    /// Sides up to most of the sphere.
    Big,
    /// About the antimeridian, so that the box is cut there.
    Anti,
    /// One side on a pole.
    Pole,
    /// About one of the twelve vertices of the icosahedron, slightly displaced.
    Vertex,
    /// No width and no height.
    Point,
    /// No width or no height.
    Line,
    /// One corner exactly on a vertex or the centroid of a zone.
    Corner,
}

/// Every kind, in the order the suites report them.
pub const KINDS: [Kind; 9] = [
    Kind::Rand,
    Kind::Tiny,
    Kind::Big,
    Kind::Anti,
    Kind::Pole,
    Kind::Vertex,
    Kind::Point,
    Kind::Line,
    Kind::Corner,
];

impl Kind {
    pub fn name(self) -> &'static str {
        match self {
            Kind::Rand => "rand",
            Kind::Tiny => "tiny",
            Kind::Big => "big",
            Kind::Anti => "anti",
            Kind::Pole => "pole",
            Kind::Vertex => "vertex",
            Kind::Point => "point",
            Kind::Line => "line",
            Kind::Corner => "corner",
        }
    }
}

/// A repeatable stream of numbers (splitmix64).
#[derive(Clone, Debug)]
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Rng {
        Rng(seed)
    }
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }
    /// Uniform in `[0, 1)`.
    pub fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    pub fn range(&mut self, a: f64, b: f64) -> f64 {
        a + (b - a) * self.unit()
    }
    /// A point uniform on the sphere, `(lat, lon)` in degrees.
    pub fn point(&mut self) -> (f64, f64) {
        let lat = (2.0 * self.unit() - 1.0).asin().to_degrees();
        (lat, self.range(-180.0, 180.0))
    }
    fn pick(&mut self, from: &[(f64, f64)]) -> (f64, f64) {
        from[(self.next_u64() % from.len() as u64) as usize]
    }
}

/// The places some kinds are made about, in degrees `(lat, lon)`: the centres of the twelve
/// base cells, and the vertices and centroids of a sample of the zones of the level.
#[derive(Clone, Debug)]
pub struct Sites {
    pub roots: Vec<(f64, f64)>,
    pub corners: Vec<(f64, f64)>,
}

/// The centres of the twelve base cells of a grid, in degrees `(lat, lon)`: the pentagons of
/// level 0, which lie at the vertices of the icosahedron. Asked of the engine once for each grid.
pub fn roots(grid: &str) -> Vec<(f64, f64)> {
    static CACHE: Mutex<BTreeMap<String, Vec<(f64, f64)>>> = Mutex::new(BTreeMap::new());
    let mut cache = CACHE.lock().unwrap_or_else(|e| e.into_inner());
    cache
        .entry(grid.to_string())
        .or_insert_with(|| {
            list_zones(grid, 0, None)
                .into_iter()
                .map(|z| centroid(grid, z))
                .collect()
        })
        .clone()
}

/// The sites for a grid and a level, the sample drawn from `rng`.
pub fn sites(grid: &str, level: i32, rng: &mut Rng) -> Sites {
    let roots = roots(grid);
    let mut corners = Vec::new();
    for _ in 0..200 {
        let (lat, lon) = rng.point();
        let z = zone_from_geo(grid, lat, lon, level);
        if z != NULL_ZONE {
            corners.push(centroid(grid, z));
            // The engine's longitudes lie within 180 degrees. Were one ever beyond, a clamp here
            // would make of it, in silence, a box of no width on the antimeridian.
            let ring = vertices(grid, z);
            assert!(
                ring.iter().all(|(_, lo)| (-180.0..=180.0).contains(lo)),
                "{grid}: a vertex of the engine beyond 180 degrees of longitude: {ring:?}"
            );
            corners.extend(ring);
        }
    }
    Sites { roots, corners }
}

/// Whether a grid is one of aperture 7.
pub fn is_aperture7(grid: &str) -> bool {
    [IGEO7, IVEA7H, RTEA7H].contains(&grid)
}

/// The number of zones at `level`, as a float: `10 * a^level + 2`.
pub fn zone_count(grid: &str, level: i32) -> f64 {
    let a: f64 = if is_aperture7(grid) { 7.0 } else { 3.0 };
    10.0 * a.powi(level) + 2.0
}

/// The width of a zone at `level` in degrees: the side of a square of the zone's area.
pub fn zone_width_deg(grid: &str, level: i32) -> f64 {
    (4.0 * PI / zone_count(grid, level)).sqrt().to_degrees()
}

/// A box in degrees as the engine receives it, in radians: each coordinate times `Pi / 180`.
pub fn radians(b: &[f64; 4]) -> [f64; 4] {
    b.map(|x| x * (PI / 180.0))
}

/// One box of a kind, in degrees. `size` is the width of a zone at the level, in degrees
/// ([`zone_width_deg`]). Longitudes stay within +-180 degrees, as the engine requires.
pub fn make_box(kind: Kind, rng: &mut Rng, size: f64, sites: &Sites) -> [f64; 4] {
    let clamp = |x: f64| x.clamp(-90.0, 90.0);
    let lon_in = |x: f64| x.clamp(-180.0, 180.0);
    // A box about a centre with half-height `h` and half-width `w`, cut at the antimeridian.
    let about = |lat: f64, lon: f64, h: f64, w: f64| -> [f64; 4] {
        let wrap = |x: f64| lon_in(x - 360.0 * ((x + 180.0) / 360.0).floor());
        let (lo, hi) = if w >= 180.0 {
            (-180.0, 180.0)
        } else {
            (wrap(lon - w), wrap(lon + w))
        };
        [clamp(lat - h), lo, clamp(lat + h), hi]
    };
    // A half-width in degrees of longitude that spans `w` on the ground at latitude `lat`.
    let widen = |lat: f64, w: f64| (w / lat.to_radians().cos().max(0.05)).min(180.0 * 0.999);
    match kind {
        Kind::Rand => {
            let (lat, lon) = rng.point();
            let h = size * 10f64.powf(rng.range(-1.0, 1.2));
            let w = size * 10f64.powf(rng.range(-1.0, 1.2));
            about(lat, lon, h.min(85.0), widen(lat, w))
        }
        Kind::Tiny => {
            let (lat, lon) = rng.point();
            let h = size * 10f64.powf(rng.range(-6.0, -1.0));
            let w = size * 10f64.powf(rng.range(-6.0, -1.0));
            about(lat, lon, h, widen(lat, w))
        }
        Kind::Big => {
            let (lat, lon) = rng.point();
            about(lat, lon, rng.range(3.0, 80.0), rng.range(3.0, 179.0))
        }
        Kind::Anti => {
            let lat = (2.0 * rng.unit() - 1.0).asin().to_degrees();
            let lon = if rng.unit() < 0.5 {
                180.0 - rng.range(0.0, size)
            } else {
                -180.0 + rng.range(0.0, size)
            };
            let h = size * 10f64.powf(rng.range(-1.0, 1.0));
            let w = size * 10f64.powf(rng.range(-0.5, 1.0));
            about(lat, lon, h.min(85.0), widen(lat, w))
        }
        Kind::Pole => {
            let north = rng.unit() < 0.5;
            let h = (size * 10f64.powf(rng.range(-2.0, 1.0))).min(85.0);
            let (lo, hi) = if rng.unit() < 0.4 {
                (-180.0, 180.0)
            } else {
                (rng.range(-180.0, 180.0), rng.range(-180.0, 180.0))
            };
            if north {
                [90.0 - h, lo, 90.0, hi]
            } else {
                [-90.0, lo, -90.0 + h, hi]
            }
        }
        Kind::Vertex => {
            let (lat, lon) = rng.pick(&sites.roots);
            let h = size * 10f64.powf(rng.range(-1.5, 1.0));
            let w = size * 10f64.powf(rng.range(-1.5, 1.0));
            let dlat = size * rng.range(-0.5, 0.5);
            let dlon = size * rng.range(-0.5, 0.5);
            about(clamp(lat + dlat), lon + dlon, h.min(85.0), widen(lat, w))
        }
        Kind::Point => {
            let (lat, lon) = rng.point();
            [lat, lon, lat, lon]
        }
        Kind::Line => {
            let (lat, lon) = rng.point();
            let l = size * 10f64.powf(rng.range(-1.0, 1.0));
            if rng.unit() < 0.5 {
                about(lat, lon, 0.0, widen(lat, l))
            } else {
                about(lat, lon, l.min(85.0), 0.0)
            }
        }
        Kind::Corner => {
            let (lat, lon) = rng.pick(&sites.corners);
            let h = size * 10f64.powf(rng.range(-1.0, 0.7));
            let w = widen(lat, size * 10f64.powf(rng.range(-1.0, 0.7)));
            let lon = lon_in(lon);
            match rng.next_u64() % 4 {
                0 => [lat, lon, clamp(lat + h), lon_in(lon + w)],
                1 => [clamp(lat - h), lon_in(lon - w), lat, lon],
                2 => [lat, lon_in(lon - w), clamp(lat + h), lon],
                _ => [clamp(lat - h), lon, lat, lon_in(lon + w)],
            }
        }
    }
}

/// The greatest ratio of a box's area to the area of a zone of its level at which the engine is
/// asked: beyond it, on aperture 3 the engine visits every cell of the box's rectangle in the
/// plane, and on aperture 7 its descent costs by the answer, which no test may spend.
pub const MAX_ZONES_IN_BOX: f64 = 20_000.0;

/// The finest level at which the engine is asked for the zones of a box on aperture 7: the level
/// at which this crate's own box query ends on those grids, so that nothing is compared beyond
/// it. At levels 18 and 19, besides, the engine answers even a box of six zone widths upon an
/// edge of the icosahedron with a null array, which no test survives.
pub const MAX_BOX_LEVEL_7: i32 = 14;

/// The margin, in radians, by which the engine demands that the extent of a zone overlap the
/// box, in each direction, for the zone to be in its answer.
pub const ENGINE_MARGIN: f64 = 2e-15;

/// The pieces of an extent `[ll.lat, ll.lon, ur.lat, ur.lon]`, in radians: itself, or the two
/// it is cut into at the antimeridian where its west lies above its east.
fn pieces(e: &[f64; 4]) -> Vec<[f64; 4]> {
    if e[1] > e[3] {
        vec![[e[0], e[1], e[2], PI], [e[0], -PI, e[2], e[3]]]
    } else {
        vec![*e]
    }
}

/// The engine's test of two extents, in radians, written here anew and apart from the crate's:
/// each piece of the one against each piece of the other, strictly, by the engine's margin.
/// It is the rule by which a zone belongs to the answer of a box, put to the engine's own
/// extent of the zone.
pub fn meets(a: &[f64; 4], b: &[f64; 4]) -> bool {
    pieces(a).iter().any(|a| {
        pieces(b).iter().any(|b| {
            b[2] - ENGINE_MARGIN > a[0]
                && a[2] - ENGINE_MARGIN > b[0]
                && b[3] - ENGINE_MARGIN > a[1]
                && a[3] - ENGINE_MARGIN > b[1]
        })
    })
}

/// The area of a box in degrees, in zones of the level.
pub fn area_in_zones(grid: &str, level: i32, b: &[f64; 4]) -> f64 {
    let [s, w, n, e] = radians(b);
    let dlon = if e < w { e - w + 2.0 * PI } else { e - w };
    dlon * (n.sin() - s.sin()) / (4.0 * PI) * zone_count(grid, level)
}

/// On aperture 3 the engine visits every cell of the rectangle that the box spans in the plane,
/// and a long thin box, which holds few zones, spans a rectangle as wide as it is long: the
/// longest side of a box, in zone widths, is kept to this.
pub const MAX_SIDE_IN_ZONE_WIDTHS: f64 = 100.0;

/// On aperture 3, the level beyond which a box upon the seam of the engine's layout is kept
/// from it.
///
/// The engine lays the icosahedron out in a plane, ten squares along a staircase, and visits
/// every cell of the rectangle that the images of the box span there. Eleven of the thirty edges
/// of the icosahedron are cut in that layout: the two sides of such an edge lie apart in the
/// plane, and a box across it spans the rectangle between them, however small the box.
///
/// Three of the eleven are the seam ([`seam`]), where the staircase ends and begins anew: a box
/// across it spans the whole layout, wherever on the seam it lies. Measured, on the three grids
/// alike, with a box of four zone widths: 0.07 s and 12 MB at level 8, 0.2 s at level 9, 0.5 s at
/// level 10, 1.6 s and 160 MB at level 11, 4.7 s and 360 MB at level 12, 14 s and 1.2 GB at level
/// 13, 43 s and 4 GB at level 14, three times as much for each level more. A box that lies a
/// hundredth of a zone width from the seam and does not touch it is answered in under 0.1 s at
/// every level to the finest.
pub const SEAM_LEVEL_3: i32 = 8;

/// On aperture 3, how far along one of the eight cut edges that are not of the seam a box may
/// lie, in zone widths from the end of the cut at which the layout is whole.
///
/// Those eight run, four from each, from the two vertices at which the cuts meet (the two ends
/// of the seam) to the vertices about them. A box across one of them spans a square whose side
/// is its distance from the other end of the cut, where the two sides come together again.
/// Measured, on the three grids alike, with a box of four zone widths: 0.04 s at 3,000 zone
/// widths from that end, 0.17 s at 10,000, 0.9 s at 30,000, 8 s at 100,000 and no answer within
/// 20 s at 300,000, at every level at which the edge is that long (it is 6,480 zone widths long
/// at level 16 and 58,300 at level 20). The south pole lies in the middle of one of the eight: a
/// box upon it costs 0.03 s at level 16, 0.12 s at level 18, 0.75 s at level 20, 6 s at level 22,
/// and is not answered within 20 s at level 26. The cost grows with the box as well: the largest
/// box allowed (a hundred zone widths a side, some ten thousand zones), across such an edge,
/// costs up to 0.4 s at 1,000 zone widths from the whole end, up to 1 s at 3,400 and up to 1.3 s
/// at 4,900.
///
/// Across the other nineteen edges the layout is continuous, and so it is at the eight vertices
/// that are not of the seam: there a box of up to a hundred zone widths is answered in under
/// 0.2 s at every level to the finest.
pub const CUT_REACH_IN_ZONE_WIDTHS: f64 = 3_500.0;

/// How near a box must come to a cut edge, in zone widths, to be taken as lying across it. The
/// engine itself is fast a hundredth of a zone width from a cut; the margin is there for the way
/// the nearness is found, on a grid of 33 by 33 points over the box: the longest side allowed
/// being [`MAX_SIDE_IN_ZONE_WIDTHS`], a point of a box that touches a cut lies within 2.3 zone
/// widths of a point of that grid.
pub const CUT_MARGIN_IN_ZONE_WIDTHS: f64 = 4.0;

/// The four vertices of the seam of the engine's layout on a grid of aperture 3, in degrees
/// `(lat, lon)`, in the order in which the seam runs through them: its three edges join each to
/// the next, the first of them over the north pole. They are the engine's own: it lists the
/// twelve base cells by their place in the layout, the ten of the staircase first, then the two
/// vertices at which the cuts meet.
pub fn seam(grid: &str) -> [(f64, f64); 4] {
    let r = roots(grid);
    [r[10], r[0], r[1], r[11]]
}

type Vector = [f64; 3];

/// The unit vector of a point `(lat, lon)` in degrees on the sphere of the ellipsoid's area
/// (WGS 84), upon which the edges of the icosahedron are arcs of great circles: the latitude is
/// taken to that sphere by the closed formula.
fn on_sphere((lat, lon): (f64, f64)) -> Vector {
    let f: f64 = 1.0 / 298.257_223_563;
    let e2 = f * (2.0 - f);
    let e = e2.sqrt();
    let q = |phi: f64| -> f64 {
        let s = phi.sin();
        (1.0 - e2) * (s / (1.0 - e2 * s * s) - ((1.0 - e * s) / (1.0 + e * s)).ln() / (2.0 * e))
    };
    let beta = (q(lat.to_radians()) / q(PI / 2.0)).clamp(-1.0, 1.0).asin();
    let l = lon.to_radians();
    [beta.cos() * l.cos(), beta.cos() * l.sin(), beta.sin()]
}

fn dot(a: &Vector, b: &Vector) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: &Vector, b: &Vector) -> Vector {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// The angle in degrees between two unit vectors, by their chord: the arc cosine of the scalar
/// product loses a small angle, and at the finest levels a zone is a hundred-millionth of a
/// radian wide.
fn angle_deg(a: &Vector, b: &Vector) -> f64 {
    let chord = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    (2.0 * (dot(&chord, &chord).sqrt() / 2.0).min(1.0).asin()).to_degrees()
}

/// The angular distance in degrees from the unit vector `p` to the arc of great circle from `a`
/// to `b`: to the great circle where the foot of `p` falls between the two, to the nearer end
/// elsewhere.
fn from_arc_deg(p: &Vector, a: &Vector, b: &Vector) -> f64 {
    let n = cross(a, b);
    if dot(&cross(a, p), &n) >= 0.0 && dot(&cross(p, b), &n) >= 0.0 {
        (dot(p, &n) / dot(&n, &n).sqrt())
            .abs()
            .min(1.0)
            .asin()
            .to_degrees()
    } else {
        angle_deg(p, a).min(angle_deg(p, b))
    }
}

/// Whether the box (degrees), on a grid of aperture 3, lies across a cut of the engine's layout
/// at a place where the rectangle it then spans is more than a test may spend: across the seam
/// beyond [`SEAM_LEVEL_3`], or across one of the eight other cut edges more than
/// [`CUT_REACH_IN_ZONE_WIDTHS`] from the end at which the layout is whole. A box lies across a cut
/// when a point of the grid of 33 by 33 points over it is within [`CUT_MARGIN_IN_ZONE_WIDTHS`] of
/// the cut.
fn across_a_costly_cut(grid: &str, level: i32, b: &[f64; 4]) -> bool {
    let r: Vec<Vector> = roots(grid).into_iter().map(on_sphere).collect();
    let [top, v0, v1, bottom] = seam(grid).map(on_sphere);
    let seam = [[top, v0], [v0, v1], [v1, bottom]];
    let width = zone_width_deg(grid, level);
    let margin = CUT_MARGIN_IN_ZONE_WIDTHS * width;
    let reach = CUT_REACH_IN_ZONE_WIDTHS * width;
    let span = if b[3] < b[1] {
        b[3] - b[1] + 360.0
    } else {
        b[3] - b[1]
    };
    (0..=32).any(|i| {
        let lat = b[0] + (b[2] - b[0]) * f64::from(i) / 32.0;
        (0..=32).any(|j| {
            let p = on_sphere((lat, b[1] + span * f64::from(j) / 32.0));
            let on_the_seam =
                level > SEAM_LEVEL_3 && seam.iter().any(|[a, c]| from_arc_deg(&p, a, c) < margin);
            // The other cuts: from the corners of the staircase but the two of the seam, the
            // third to the tenth base cell, alternately to the one meeting vertex and the other.
            on_the_seam
                || (2..10).any(|k| {
                    let meeting = if k % 2 == 0 { &top } else { &bottom };
                    angle_deg(&p, &r[k]) > reach && from_arc_deg(&p, &r[k], meeting) < margin
                })
        })
    })
}

/// Whether the engine may be asked for the zones at `level` of the box `b`, in degrees.
/// Conservative: a box wrongly let through hangs the gate or kills it, a box wrongly kept back
/// only hides a comparison, which the suites count. Each rule is stated where it is applied,
/// with the measurement it rests on.
///
/// **What stays open.** On aperture 7 the engine answers a box in which it finds no zone with a
/// null array, at any level, and the asking process aborts. The boxes that are so answered by
/// rule are kept back; a box of positive area in which the engine finds no zone at a level up
/// to [`MAX_BOX_LEVEL_7`] would still abort the process, and none has been met.
///
/// The whole world, which the engine takes in a short way, is asked as `None` of
/// [`list_zones`] and needs no test.
pub fn engine_survives(grid: &str, level: i32, b: &[f64; 4]) -> bool {
    let seven = is_aperture7(grid);
    // A level the grid has not; on aperture 7, a level beyond the one at which the box query
    // ends.
    let finest = if seven {
        MAX_BOX_LEVEL_7
    } else {
        max_level(grid)
    };
    if !(0..=finest).contains(&level) {
        return false;
    }
    // A coordinate that is not finite, a latitude beyond a pole, a south above its north: the
    // engine has no answer for them. A longitude beyond 180 degrees: in a box over the
    // antimeridian the engine cuts the box there and calls itself on the pieces without end.
    if !b.iter().all(|x| x.is_finite()) {
        return false;
    }
    let r = radians(b);
    let [s, w, n, e] = r;
    let half = PI / 2.0;
    if !(-half <= s && s <= n && n <= half) || w.abs() > PI || e.abs() > PI {
        return false;
    }
    // A sliver, a side under 1e-11 radian that is not nought: the engine steps over the box by
    // a hundredth of each side, and its walk makes no progress there.
    let dlat = n - s;
    let dlon = if e < w { e - w + 2.0 * PI } else { e - w };
    if (dlat != 0.0 && dlat < 1e-11) || (dlon != 0.0 && dlon < 1e-11) {
        return false;
    }
    // A point on a pole, which the engine answers with a null array.
    if dlat == 0.0 && dlon == 0.0 && s.abs() == half {
        return false;
    }
    // On aperture 7, a box that no extent can overlap by the engine's margin, for which it
    // answers a null array at every level: a box whose south lies within the margin of the north
    // pole or whose north within the margin of the south pole, and a box whose west lies within
    // the margin of +180 degrees or whose east within the margin of -180 degrees (over the
    // antimeridian, where the engine cuts the box in two, both at once). With no sliver left,
    // these are the boxes of no height on a pole and of no width on the antimeridian. Measured
    // on the three grids, each call in a process of its own: such a box aborts the process at
    // levels 0, 2, 3, 4, 5 and 14, on the pole or the antimeridian itself and one or two units
    // in the last place of a degree from it; a box of no height 1e-12 degree from a pole, a box
    // of no width 1e-12 degree from the antimeridian or at any other longitude, and a point
    // anywhere else are answered.
    if seven {
        let (west_out, east_out) = (PI - ENGINE_MARGIN <= w, e - ENGINE_MARGIN <= -PI);
        let lon_out = if w > e {
            west_out && east_out
        } else {
            west_out || east_out
        };
        if half - ENGINE_MARGIN <= s || n - ENGINE_MARGIN <= -half || lon_out {
            return false;
        }
    }
    // A box of more zones' area than a test may spend, at any level, on either aperture.
    let area = dlon * (n.sin() - s.sin());
    if area / (4.0 * PI) * zone_count(grid, level) > MAX_ZONES_IN_BOX {
        return false;
    }
    if seven {
        return true;
    }
    // On aperture 3, a box whose longest side on the ground, the east-west one taken at the
    // latitude nearest the equator, is beyond what the engine's rectangle allows: a band of 100
    // degrees that holds 15,000 zones costs it 10 s at level 14.
    let nearest = if s <= 0.0 && 0.0 <= n {
        0.0
    } else {
        s.abs().min(n.abs())
    };
    let side = dlat.to_degrees().max(dlon.to_degrees() * nearest.cos());
    if side > MAX_SIDE_IN_ZONE_WIDTHS * zone_width_deg(grid, level) {
        return false;
    }
    // On aperture 3, a box across a cut of the engine's layout where that costs it the layout.
    !across_a_costly_cut(grid, level, b)
}

/// The six grids' names, aperture 7 first.
pub const GRIDS: [&str; 6] = [IGEO7, IVEA7H, RTEA7H, ISEA3H, IVEA3H, RTEA3H];
