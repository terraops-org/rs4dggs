//! The census of the cleaning, compiled for the tests alone: a sample of zones about the poles,
//! the antimeridian and a vertex of the icosahedron, and a judge of the validity of a shape that
//! needs no geometry library.
//!
//! The sample is, on each of the six grids and at every level, the zone at each of fifteen
//! points and its neighbours, each zone once. The points are both poles; four points beside
//! them (89.99 N 11.2 E, 89.9 N 120 W, 89.99 S 168.8 W, 89.9 S 60 E); seven on the antimeridian
//! (on the equator, at 45 N, at 45 S, at 60 N a thousandth of a degree to its west, at 20 S a
//! thousandth to its east, at 75 N and at 75 S); the vertex of the icosahedron at
//! 58.397145907431 N 11.2 E; and Lisbon.

use std::collections::{BTreeSet, HashSet};

use rs4dggs::{AnyGrid, GeoPoint, ZoneId};

use super::{Shape, clean, prepared};

/// The fifteen points of the sample, as (latitude, longitude).
const POINTS: [(f64, f64); 15] = [
    (90.0, 0.0),
    (-90.0, 0.0),
    (89.99, 11.2),
    (89.9, -120.0),
    (-89.99, -168.8),
    (-89.9, 60.0),
    (0.0, 180.0),
    (45.0, 180.0),
    (-45.0, -180.0),
    (60.0, 179.999),
    (-20.0, -179.999),
    (75.0, 180.0),
    (-75.0, 180.0),
    (58.397145907431, 11.2),
    (38.7223, -9.1393),
];

/// The zones within `radius` steps of neighbours of the zone at each of the fifteen points, at
/// every level of `grid`, each once, in the order first met.
fn around(grid: &AnyGrid, radius: usize) -> Vec<ZoneId> {
    let mut seen = HashSet::new();
    let mut out = Vec::new();
    for level in 0..=grid.max_resolution() {
        for (lat, lon) in POINTS {
            // The disk about this point, walked by itself, so that a zone met before from
            // another point still leads on to its neighbours.
            let centre = grid.zone_from_geo(lat, lon, level).unwrap();
            let mut disk = HashSet::from([centre]);
            let mut ring = vec![centre];
            for step in 0..=radius {
                let mut next = Vec::new();
                for z in ring {
                    if z == ZoneId::NULL {
                        continue;
                    }
                    if seen.insert(z) {
                        out.push(z);
                    }
                    if step < radius {
                        next.extend(grid.neighbors(z).into_iter().filter(|n| disk.insert(*n)));
                    }
                }
                ring = next;
            }
        }
    }
    out
}

/// The zones of the census on `grid`: at every level, the zone at each of the fifteen points and
/// its neighbours, each zone once.
pub(crate) fn sample(grid: &AnyGrid) -> Vec<ZoneId> {
    around(grid, 1)
}

/// The rings of a shape.
fn rings(shape: &Shape) -> Vec<&Vec<[f64; 2]>> {
    match shape {
        Shape::Polygon(r) => vec![r],
        Shape::MultiPolygon(v) => v.iter().collect(),
    }
}

/// Twice the signed area of the triangle `a`, `b`, `c`: positive if it turns anticlockwise.
fn orient(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
}

/// Whether `p`, known to lie on the line through `a` and `b`, lies on the segment between them.
fn within(a: [f64; 2], b: [f64; 2], p: [f64; 2]) -> bool {
    a[0].min(b[0]) <= p[0]
        && p[0] <= a[0].max(b[0])
        && a[1].min(b[1]) <= p[1]
        && p[1] <= a[1].max(b[1])
}

/// Whether the closed segments `a`-`b` and `c`-`d` share a point.
fn meet(a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2]) -> bool {
    if a[0].max(b[0]) < c[0].min(d[0])
        || c[0].max(d[0]) < a[0].min(b[0])
        || a[1].max(b[1]) < c[1].min(d[1])
        || c[1].max(d[1]) < a[1].min(b[1])
    {
        return false;
    }
    let (d1, d2) = (orient(c, d, a), orient(c, d, b));
    let (d3, d4) = (orient(a, b, c), orient(a, b, d));
    if ((d1 > 0.0 && d2 < 0.0) || (d1 < 0.0 && d2 > 0.0))
        && ((d3 > 0.0 && d4 < 0.0) || (d3 < 0.0 && d4 > 0.0))
    {
        return true;
    }
    (d1 == 0.0 && within(c, d, a))
        || (d2 == 0.0 && within(c, d, b))
        || (d3 == 0.0 && within(a, b, c))
        || (d4 == 0.0 && within(a, b, d))
}

/// Whether the closed ring `ring` (its last position its first) is simple: no two edges that
/// are not neighbours share a point, and no two neighbours overlap, that is, the second does not
/// turn back along the first. A position repeated at once is one position, as GeoJSON readers
/// take it.
fn simple(ring: &[[f64; 2]]) -> Result<(), String> {
    let mut r: Vec<[f64; 2]> = Vec::with_capacity(ring.len());
    for &p in ring {
        if r.last() != Some(&p) {
            r.push(p);
        }
    }
    let n = r.len() - 1; // edges: r[i] to r[i + 1]
    for i in 0..n {
        let (a, b, c) = (r[i], r[i + 1], r[(i + 2) % n]);
        if orient(a, b, c) == 0.0
            && (c[0] - b[0]) * (a[0] - b[0]) + (c[1] - b[1]) * (a[1] - b[1]) > 0.0
        {
            return Err(format!("the edges at {b:?} turn back along each other"));
        }
        // Every later edge but the two that share a position with this one.
        for j in i + 2..n {
            if i == 0 && j == n - 1 {
                continue;
            }
            if meet(a, b, r[j], r[j + 1]) {
                return Err(format!(
                    "the edge {a:?} to {b:?} meets the edge {:?} to {:?}",
                    r[j],
                    r[j + 1]
                ));
            }
        }
    }
    Ok(())
}

/// Whether `shape` is a valid GeoJSON geometry: every ring closed, of four positions or more
/// once a position repeated at once is counted one, anticlockwise (its signed area positive),
/// every coordinate finite, longitudes within -180 to 180 and latitudes within -90 to 90, and
/// simple.
///
/// Each ring is judged by itself: that the polygons of a `MultiPolygon` do not overlap is not
/// judged here. It holds by construction for a ring that is simple before the cut, whose pieces
/// lie on the two sides of the meridian or between distinct pairs of its crossings.
pub(crate) fn valid(shape: &Shape) -> std::result::Result<(), String> {
    let rings = rings(shape);
    if rings.is_empty() {
        return Err("no ring".into());
    }
    for r in rings {
        if r.first() != r.last() {
            return Err("a ring is not closed".into());
        }
        let mut distinct = r.clone();
        distinct.dedup();
        if distinct.len() < 4 {
            return Err(format!("a ring of {} positions", distinct.len()));
        }
        for p in r {
            if !p[0].is_finite() || !p[1].is_finite() {
                return Err(format!("{p:?} is not finite"));
            }
            if !(-180.0..=180.0).contains(&p[0]) || !(-90.0..=90.0).contains(&p[1]) {
                return Err(format!("{p:?} is out of range"));
            }
        }
        if super::area(r) <= 0.0 {
            return Err(format!("a ring of area {}", super::area(r)));
        }
        simple(r)?;
    }
    Ok(())
}

/// `shape` as it is written with `decimals` places: every coordinate rounded.
pub(crate) fn as_written(shape: &Shape, decimals: u8) -> Shape {
    let ring = |r: &Vec<[f64; 2]>| -> Vec<[f64; 2]> {
        r.iter()
            .map(|p| {
                [
                    super::rounded(p[0], decimals),
                    super::rounded(p[1], decimals),
                ]
            })
            .collect()
    };
    match shape {
        Shape::Polygon(r) => Shape::Polygon(ring(r)),
        Shape::MultiPolygon(v) => Shape::MultiPolygon(v.iter().map(ring).collect()),
    }
}

/// The ring `prepared` gives, with the turn it makes, closed as a ring of the plane: as it is
/// for a ring that makes no turn, and for a ring about a pole carried one turn on and closed
/// along the pole, as the cleaning closes it.
fn closure(u: &[[f64; 2]], turn: f64) -> Vec<[f64; 2]> {
    let mut r = u.to_vec();
    if turn.abs() >= 180.0 {
        let pole = if turn > 0.0 { 90.0 } else { -90.0 };
        let end = [u[0][0] + turn.signum() * 360.0, u[0][1]];
        r.push(end);
        r.push([end[0], pole]);
        r.push([u[0][0], pole]);
    }
    r.push(u[0]);
    r
}

/// Whether `p` lies within the extent of the zone, by at most `slack` degrees; a longitude is
/// read eastward from the western end, through 180 degrees where the western end is the greater.
fn bracketed(e: &rs4dggs::Extent, p: [f64; 2], slack: f64) -> bool {
    let lat = e.ll.lat - slack <= p[1] && p[1] <= e.ur.lat + slack;
    let lon = if e.ll.lon <= e.ur.lon {
        e.ll.lon - slack <= p[0] && p[0] <= e.ur.lon + slack
    } else {
        p[0] >= e.ll.lon - slack || p[0] <= e.ur.lon + slack
    };
    lat && lon
}

/// What the census finds on one grid.
#[derive(Debug, Default)]
struct Census {
    zones: usize,
    /// Zones whose ring has no area to draw (every point on a pole, or fewer than three distinct).
    no_area: usize,
    /// The zones whose shape is refused, by name.
    refused: BTreeSet<String>,
    /// Of those, the ones whose ring as the library gives it holds a point on a pole.
    refused_with_pole: usize,
    /// The zones whose ring, as it stands before the cut, is not simple.
    not_simple_before: BTreeSet<String>,
    /// The zones cut into more than two rings.
    more_than_two: BTreeSet<String>,
    /// The positions of valid shapes outside the zone's extent by more than 1e-9 degrees.
    outside: Vec<String>,
}

fn census(
    grid: &AnyGrid,
    ring_of: &dyn Fn(ZoneId) -> Vec<GeoPoint>,
    clean: &dyn Fn(&[GeoPoint]) -> Shape,
) -> Census {
    let mut c = Census::default();
    for z in sample(grid) {
        c.zones += 1;
        let ring = ring_of(z);
        let Some((u, turn)) = prepared(&ring) else {
            c.no_area += 1;
            continue;
        };
        let name = grid.text_id(z);
        if simple(&closure(&u, turn)).is_err() {
            c.not_simple_before.insert(name.clone());
        }
        let shape = clean(&ring);
        if rings(&shape).len() > 2 {
            c.more_than_two.insert(name.clone());
        }
        if valid(&shape).is_err() {
            if ring.iter().any(|p| p.lat.abs() == 90.0) {
                c.refused_with_pole += 1;
            }
            c.refused.insert(name);
            continue;
        }
        let e = grid.extent(z).unwrap();
        for p in rings(&shape).into_iter().flatten() {
            if !bracketed(&e, *p, 1e-9) {
                c.outside.push(format!("{name} {p:?} {e:?}"));
            }
        }
    }
    c
}

/// What the census must find on one grid.
struct Expected {
    grid: &'static str,
    zones: usize,
    /// Zones whose ring has no area to draw, plain or refined alike: the three zones of the
    /// finest level of a Z7 grid whose every point lies on a pole.
    no_area: usize,
    /// The refined rings whose shape is refused, every one a ring that is not simple before the
    /// cut, as the engine gives it with its points on a pole replaced.
    refused: &'static [&'static str],
    /// How many of those hold a point on a pole.
    with_pole: usize,
    /// The refined rings cut into more than two polygons.
    more_than_two: &'static [&'static str],
}

/// The finest zones on the aperture-3 grids whose ring holds 45 of its 47 points, or 246 of
/// its 248, on a pole.
const A3: [&str; 4] = [
    "Q0-1486BA0-C",
    "Q3-6954FE0D5D2E0-D",
    "Q5-34AA7EFC6B3A0-D",
    "Q8-34AA7F2578AE0-C",
];

/// Rings of levels 15 to 18 on the Z7 grids that cross themselves: two beside the north pole on
/// the antimeridian, two beside the vertex of the icosahedron at 58.4 N 11.2 E.
const Z7: [&str; 4] = [
    "00000000000000000001",
    "000000000000000001",
    "01323232323232322",
    "0132323232323232322",
];

const CUT_IN_THREE: [&str; 2] = ["01323232323232322", "0132323232323232322"];

const EXPECTED: [Expected; 6] = [
    Expected {
        grid: "ISEA3H",
        zones: 2991,
        no_area: 0,
        refused: &A3,
        with_pole: 4,
        more_than_two: &[],
    },
    Expected {
        grid: "IVEA3H",
        zones: 2991,
        no_area: 0,
        refused: &A3,
        with_pole: 4,
        more_than_two: &[],
    },
    Expected {
        grid: "RTEA3H",
        zones: 2985,
        no_area: 0,
        refused: &[
            "N5-127E95CC238-A",
            "N5-127E975160B-A",
            "Q0-1486BA0-C",
            "Q3-6954FE0D5D2E0-D",
            "Q5-34AA7ED35DC5F-A",
            "Q5-34AA7EFC6B3A0-D",
            "Q8-34AA7F2578AE0-C",
        ],
        with_pole: 5,
        more_than_two: &[],
    },
    Expected {
        grid: "IGEO7",
        zones: 1757,
        no_area: 3,
        refused: &Z7,
        with_pole: 0,
        more_than_two: &["01012", "01323232323232322", "0132323232323232322"],
    },
    Expected {
        grid: "IVEA7H",
        zones: 1757,
        no_area: 3,
        refused: &Z7,
        with_pole: 0,
        more_than_two: &CUT_IN_THREE,
    },
    Expected {
        grid: "RTEA7H",
        zones: 1757,
        no_area: 3,
        refused: &[
            "000000000000000000005",
            "00000000000000000001",
            "000000000000000001",
            "01323232323232322",
            "0132323232323232322",
            "084545454545454541",
            "084545454545454545",
            "114545454545454545",
            "11454545454545454545",
        ],
        with_pole: 3,
        more_than_two: &CUT_IN_THREE,
    },
];

fn names(list: &[&str]) -> BTreeSet<String> {
    list.iter().map(|s| s.to_string()).collect()
}

fn refined(grid: &AnyGrid) -> impl Fn(ZoneId) -> Vec<GeoPoint> + '_ {
    |z| grid.refined_vertices(z, 0).unwrap()
}

fn plain(grid: &AnyGrid) -> impl Fn(ZoneId) -> Vec<GeoPoint> + '_ {
    |z| grid.vertices(z)
}

#[test]
fn the_sample_holds_the_zones_it_held() {
    for e in &EXPECTED {
        let grid = rs4dggs::get_grid(e.grid).unwrap();
        assert_eq!(sample(&grid).len(), e.zones, "{}", e.grid);
    }
}

#[test]
fn every_refined_ring_simple_before_the_cut_gives_a_valid_shape_within_its_extent() {
    for e in &EXPECTED {
        let grid = rs4dggs::get_grid(e.grid).unwrap();
        let c = census(&grid, &refined(&grid), &clean);
        let what = e.grid;
        assert_eq!(c.zones, e.zones, "{what}");
        assert_eq!(c.no_area, e.no_area, "{what}");
        assert_eq!(c.refused, names(e.refused), "{what}: the shapes refused");
        assert_eq!(c.refused_with_pole, e.with_pole, "{what}");
        // Every shape refused is that of a ring not simple before the cut, and every ring not
        // simple before the cut gives a shape refused.
        assert_eq!(
            c.not_simple_before, c.refused,
            "{what}: the rings not simple before the cut"
        );
        assert_eq!(
            c.more_than_two,
            names(e.more_than_two),
            "{what}: cut in more than two"
        );
        assert!(
            c.outside.is_empty(),
            "{what}: outside the extent: {:?}",
            c.outside
        );
    }
}

/// What the census must find with eight decimals, the writers' default, on each grid: the
/// shapes refused, and the zones written without geometry.
///
/// It differs from the census in the shortest form in two ways. Rounding brings the crossings of
/// eight of the engine's rings on the RTEA grids within half a step, where the cleaning removes
/// them as spikes, so that those shapes are valid. And the four rings of each aperture-3 grid
/// that the engine draws out and back along the pole enclose nothing once rounded, and are
/// written without geometry. No other shape is refused.
const EIGHT_DECIMALS: [(&str, &[&str], &[&str]); 6] = [
    ("ISEA3H", &[], &A3),
    ("IVEA3H", &[], &A3),
    ("RTEA3H", &[], &A3),
    ("IGEO7", &Z7, &[]),
    ("IVEA7H", &Z7, &[]),
    (
        "RTEA7H",
        &[
            "00000000000000000001",
            "000000000000000001",
            "01323232323232322",
            "0132323232323232322",
        ],
        &[],
    ),
];

#[test]
fn with_eight_decimals_only_rings_the_engine_draws_crossing_are_refused() {
    for (name, refused, undrawn) in EIGHT_DECIMALS {
        let grid = rs4dggs::get_grid(name).unwrap();
        let (mut bad, mut none) = (BTreeSet::new(), BTreeSet::new());
        for z in sample(&grid) {
            let ring = grid.refined_vertices(z, 0).unwrap();
            if prepared(&ring).is_none() {
                continue;
            }
            match super::drawn(&ring, Some(8)) {
                None => none.insert(grid.text_id(z)),
                Some(s) if valid(&as_written(&s, 8)).is_err() => bad.insert(grid.text_id(z)),
                Some(_) => false,
            };
        }
        assert_eq!(bad, names(refused), "{name}: the shapes refused");
        assert_eq!(none, names(undrawn), "{name}: the zones without geometry");
    }
}

#[test]
fn every_plain_ring_gives_a_valid_shape() {
    for e in &EXPECTED {
        let grid = rs4dggs::get_grid(e.grid).unwrap();
        let c = census(&grid, &plain(&grid), &clean);
        assert_eq!((c.zones, c.no_area), (e.zones, e.no_area), "{}", e.grid);
        assert!(c.refused.is_empty(), "{}: {:?}", e.grid, c.refused);
        assert!(c.not_simple_before.is_empty(), "{}", e.grid);
    }
}

#[test]
fn the_judge_accepts_simple_rings_and_refuses_the_others() {
    let square = vec![
        [0.0, 0.0],
        [10.0, 0.0],
        [10.0, 10.0],
        [0.0, 10.0],
        [0.0, 0.0],
    ];
    assert_eq!(valid(&Shape::Polygon(square.clone())), Ok(()));
    let far = square.iter().map(|p| [p[0] + 20.0, p[1]]).collect();
    assert_eq!(
        valid(&Shape::MultiPolygon(vec![square.clone(), far])),
        Ok(())
    );
    // Closed about the north pole along it, as the cleaning closes a ring about a pole.
    let cap = vec![
        [-180.0, 80.0],
        [-90.0, 80.0],
        [0.0, 80.0],
        [90.0, 80.0],
        [180.0, 80.0],
        [180.0, 90.0],
        [-180.0, 90.0],
        [-180.0, 80.0],
    ];
    assert_eq!(valid(&Shape::Polygon(cap)), Ok(()));
    // A position repeated at once is one position.
    let repeated = vec![
        [0.0, 0.0],
        [10.0, 0.0],
        [10.0, 0.0],
        [10.0, 10.0],
        [0.0, 0.0],
    ];
    assert_eq!(valid(&Shape::Polygon(repeated)), Ok(()));

    let refused = |r: Vec<[f64; 2]>| valid(&Shape::Polygon(r)).is_err();
    // Two edges that cross: a bow tie.
    assert!(refused(vec![
        [0.0, 0.0],
        [10.0, 10.0],
        [10.0, 0.0],
        [0.0, 10.0],
        [0.0, 0.0]
    ]));
    // A vertex that touches an edge that is not its neighbour.
    assert!(refused(vec![
        [0.0, 0.0],
        [10.0, 0.0],
        [10.0, 10.0],
        [5.0, 0.0],
        [0.0, 10.0],
        [0.0, 0.0]
    ]));
    // A vertex met twice: two loops that touch at one point.
    assert!(refused(vec![
        [0.0, 0.0],
        [5.0, 5.0],
        [10.0, 0.0],
        [10.0, 10.0],
        [5.0, 5.0],
        [0.0, 10.0],
        [0.0, 0.0]
    ]));
    // A spike: an edge that turns back along the one before it.
    assert!(refused(vec![
        [0.0, 0.0],
        [10.0, 0.0],
        [15.0, 0.0],
        [10.0, 0.0],
        [10.0, 10.0],
        [0.0, 0.0]
    ]));
    // A triangle whose third point turns back along its first edge.
    assert!(simple(&[[0.0, 0.0], [10.0, 0.0], [5.0, 0.0], [0.0, 0.0]]).is_err());
    // Two edges that lie along each other, not neighbours.
    assert!(refused(vec![
        [0.0, 0.0],
        [10.0, 0.0],
        [10.0, 5.0],
        [4.0, 0.0],
        [2.0, 0.0],
        [0.0, 5.0],
        [0.0, 0.0]
    ]));
    // Clockwise, open, too short, out of range, not finite.
    assert!(refused(square.iter().rev().copied().collect()));
    assert!(refused(square[..4].to_vec()));
    assert!(refused(vec![[0.0, 0.0], [10.0, 0.0], [0.0, 0.0]]));
    assert!(refused(vec![
        [0.0, 0.0],
        [10.0, 0.0],
        [10.0, 0.0],
        [0.0, 0.0]
    ]));
    assert!(refused(vec![
        [170.0, 0.0],
        [190.0, 0.0],
        [190.0, 10.0],
        [170.0, 0.0]
    ]));
    assert!(refused(vec![
        [0.0, 80.0],
        [10.0, 80.0],
        [10.0, 91.0],
        [0.0, 80.0]
    ]));
    assert!(refused(vec![
        [0.0, 0.0],
        [f64::NAN, 0.0],
        [10.0, 10.0],
        [0.0, 0.0]
    ]));
    assert!(valid(&Shape::MultiPolygon(Vec::new())).is_err());
}
