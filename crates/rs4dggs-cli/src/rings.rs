//! A zone's vertices made into GeoJSON polygons that QGIS, GDAL and web maps draw correctly:
//! split at the antimeridian, closed about the poles (RFC 7946, section 3.1.9). A pure function.
//!
//! The library gives its rings anticlockwise as seen from outside the sphere, on every grid (its
//! `Grid::vertices` names the few engine rings that are not). The direction is judged again here
//! all the same, as a guard, on the ring unwrapped in longitude, so that a ring crossing the
//! antimeridian or touching a pole is judged by its true shape.

use rs4dggs::GeoPoint;

/// One polygon or two, each a closed exterior ring of `[lon, lat]` points, anticlockwise.
#[derive(Debug, Clone, PartialEq)]
pub enum Shape {
    Polygon(Vec<[f64; 2]>),
    MultiPolygon(Vec<Vec<[f64; 2]>>),
}

/// How far, in degrees of longitude, an edge may fall short of half a turn and still run over a
/// pole. The engine's own edges over a pole miss it by up to 7e-6; every other edge of every grid
/// lies at least 21 degrees from it.
const HALF_TURN_TOL: f64 = 1e-3;

fn on_pole(p: &[f64; 2]) -> bool {
    p[1].abs() == 90.0
}

/// `d` brought within [-180, 180] by whole turns.
fn wrap(d: f64) -> f64 {
    d - 360.0 * (d / 360.0).round()
}

fn closed(mut ring: Vec<[f64; 2]>) -> Vec<[f64; 2]> {
    if let (Some(first), Some(last)) = (ring.first().copied(), ring.last().copied()) {
        if first != last {
            ring.push(first);
        }
    }
    ring
}

/// `ring` without consecutive repeats, the last compared with the first.
fn dedup(ring: Vec<[f64; 2]>) -> Vec<[f64; 2]> {
    let mut out: Vec<[f64; 2]> = Vec::with_capacity(ring.len());
    for p in ring {
        if out.last() != Some(&p) {
            out.push(p);
        }
    }
    while out.len() > 1 && out.first() == out.last() {
        out.pop();
    }
    out
}

/// The point where the straight edge from `a` to `b` meets longitude `m`.
fn crossing(a: [f64; 2], b: [f64; 2], m: f64) -> [f64; 2] {
    let t = (m - a[0]) / (b[0] - a[0]);
    [m, a[1] + t * (b[1] - a[1])]
}

/// The part of `ring` with longitude at most `m` (`below`) or at least `m` (Sutherland and
/// Hodgman against one half-plane, which keeps the ring's direction).
fn clip(ring: &[[f64; 2]], m: f64, below: bool) -> Vec<[f64; 2]> {
    let inside = |p: &[f64; 2]| if below { p[0] <= m } else { p[0] >= m };
    let n = ring.len();
    let mut out = Vec::with_capacity(n + 2);
    for i in 0..n {
        let (prev, cur) = (ring[(i + n - 1) % n], ring[i]);
        match (inside(&prev), inside(&cur)) {
            (true, true) => out.push(cur),
            (true, false) => out.push(crossing(prev, cur, m)),
            (false, true) => {
                out.push(crossing(prev, cur, m));
                out.push(cur);
            }
            (false, false) => {}
        }
    }
    dedup(out)
}

/// `pts` with each longitude within 180 of the one before, and the turn the ring makes in
/// longitude on the way back to its start (0 for a ring that does not surround a pole). An edge
/// from one pole point to another, half a turn apart, is ambiguous: it is taken the way that
/// leaves the ring with no turn, as the ring then only touches the pole. Each longitude is the
/// engine's own, plus whole turns.
fn unwrapped(pts: &[[f64; 2]]) -> (Vec<[f64; 2]>, f64) {
    let n = pts.len();
    let mut deltas: Vec<f64> = (0..n)
        .map(|i| wrap(pts[(i + 1) % n][0] - pts[i][0]))
        .collect();
    let half = |i: usize| {
        on_pole(&pts[i])
            && on_pole(&pts[(i + 1) % n])
            && (deltas[i].abs() - 180.0).abs() < HALF_TURN_TOL
    };
    if let Some(k) = (0..n).find(|&i| half(i)) {
        let rest: f64 = deltas.iter().sum::<f64>() - deltas[k];
        if (rest + deltas[k]).abs() > 180.0 {
            deltas[k] -= 360.0 * deltas[k].signum();
        }
    }
    let mut lon = pts[0][0];
    let mut u = Vec::with_capacity(n);
    for (i, p) in pts.iter().enumerate() {
        u.push([p[0] + 360.0 * ((lon - p[0]) / 360.0).round(), p[1]]);
        lon += deltas[i];
    }
    (u, deltas.iter().sum())
}

/// The signed area of `ring` by the shoelace formula, about its first point.
fn area(ring: &[[f64; 2]]) -> f64 {
    let o = ring[0];
    let n = ring.len();
    (0..n)
        .map(|i| {
            let (a, b) = (ring[i], ring[(i + 1) % n]);
            (a[0] - o[0]) * (b[1] - o[1]) - (b[0] - o[0]) * (a[1] - o[1])
        })
        .sum::<f64>()
        / 2.0
}

pub fn clean(vertices: &[GeoPoint]) -> Shape {
    let raw: Vec<[f64; 2]> = vertices.iter().map(|v| [v.lon, v.lat]).collect();
    if raw.iter().flatten().any(|c| !c.is_finite()) {
        return Shape::Polygon(closed(raw));
    }
    let pts = dedup(raw.clone());
    // No area: every point on a pole (written as it comes), or fewer than three distinct points.
    if raw.iter().all(on_pole) {
        return Shape::Polygon(closed(raw));
    }
    if pts.len() < 3 {
        return Shape::Polygon(closed(pts));
    }
    // A great-circle edge half a turn of longitude long runs over a pole: a pole point is put in
    // its way, and then treated as any other. (Were both ends on the equator, the edge would be
    // antipodal and go to the south; no edge of the library is.)
    let n = pts.len();
    let mut over: Vec<[f64; 2]> = Vec::with_capacity(n + 2);
    for i in 0..n {
        let (a, b) = (pts[i], pts[(i + 1) % n]);
        over.push(a);
        if !on_pole(&a) && !on_pole(&b) && ((b[0] - a[0]).abs() - 180.0).abs() < HALF_TURN_TOL {
            over.push([a[0], if a[1] + b[1] > 0.0 { 90.0 } else { -90.0 }]);
        }
    }
    let pts = over;
    let n = pts.len();
    // A vertex on a pole becomes two points on it, at the longitudes of the nearest vertices
    // before and after it that are not on a pole; a run of them becomes one edge along it.
    let nearest = |i: usize, step: usize| -> f64 {
        let mut j = i;
        loop {
            j = (j + step) % n;
            if !on_pole(&pts[j]) {
                return pts[j][0];
            }
        }
    };
    let mut expanded = Vec::with_capacity(n + 4);
    for (i, p) in pts.iter().enumerate() {
        if on_pole(p) {
            if !on_pole(&pts[(i + n - 1) % n]) {
                expanded.push([nearest(i, n - 1), p[1]]);
                expanded.push([nearest(i, 1), p[1]]);
            }
        } else {
            expanded.push(*p);
        }
    }
    let mut pts = dedup(expanded);
    let (mut u, mut turn) = unwrapped(&pts);
    // A guard: the library's rings are anticlockwise, and one that is not is turned about, so that
    // every exterior ring is. A ring about a pole is clockwise if it turns against that pole.
    let mean_lat = pts.iter().map(|p| p[1]).sum::<f64>() / pts.len() as f64;
    let clockwise = if turn.abs() < 180.0 {
        area(&u) < 0.0
    } else {
        (turn > 0.0) != (mean_lat > 0.0)
    };
    if clockwise {
        pts.reverse();
        (u, turn) = unwrapped(&pts);
    }
    if turn.abs() < 180.0 {
        let min = u.iter().map(|p| p[0]).fold(f64::INFINITY, f64::min);
        let max = u.iter().map(|p| p[0]).fold(f64::NEG_INFINITY, f64::max);
        if min >= -180.0 && max <= 180.0 {
            return Shape::Polygon(closed(u));
        }
        let m = if max > 180.0 { 180.0 } else { -180.0 };
        let shift = if m > 0.0 { -360.0 } else { 360.0 };
        let near = clip(&u, m, m > 0.0);
        let far: Vec<[f64; 2]> = clip(&u, m, m < 0.0)
            .into_iter()
            .map(|p| [p[0] + shift, p[1]])
            .collect();
        return match (near.len() >= 3, far.len() >= 3) {
            (true, true) => Shape::MultiPolygon(vec![closed(near), closed(far)]),
            (true, false) => Shape::Polygon(closed(near)),
            _ => Shape::Polygon(closed(far)),
        };
    }
    // About a pole: north if the ring turns eastward (+360), south if westward (-360).
    let north = turn > 0.0;
    let mut path = u.clone();
    path.push([u[0][0] + if north { 360.0 } else { -360.0 }, u[0][1]]);
    // The antimeridian the path crosses: 180 + 360k strictly inside its span, or at its start.
    let (lo, hi) = if north {
        (path[0][0], path[path.len() - 1][0])
    } else {
        (path[path.len() - 1][0], path[0][0])
    };
    let m = 180.0 + 360.0 * ((lo - 180.0) / 360.0).ceil();
    let m = if m > hi { m - 360.0 } else { m };
    let i = (0..path.len() - 1)
        .find(|&i| {
            let (a, b) = (path[i][0], path[i + 1][0]);
            a != b && a.min(b) <= m && m <= a.max(b)
        })
        .expect("a ring about a pole crosses every meridian");
    let x = crossing(path[i], path[i + 1], m);
    let turn_back = if north { 360.0 } else { -360.0 };
    // From the crossing round the ring and back to the crossing, one turn on.
    let mut edge = vec![x];
    edge.extend_from_slice(&path[i + 1..]);
    edge.extend(path[1..=i].iter().map(|p| [p[0] + turn_back, p[1]]));
    edge.push([x[0] + turn_back, x[1]]);
    // Moved so that it runs from -180 to 180 (north) or from 180 to -180 (south).
    let offset = if north { -180.0 - m } else { 180.0 - m };
    let mut ring: Vec<[f64; 2]> = edge.into_iter().map(|p| [p[0] + offset, p[1]]).collect();
    let pole = if north { 90.0 } else { -90.0 };
    let (end, start) = if north {
        (180.0, -180.0)
    } else {
        (-180.0, 180.0)
    };
    ring.push([end, pole]);
    ring.push([start, pole]);
    Shape::Polygon(closed(dedup(ring)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::{HashSet, VecDeque};

    fn gp(lon: f64, lat: f64) -> GeoPoint {
        GeoPoint { lat, lon }
    }

    fn rings(s: &Shape) -> Vec<&Vec<[f64; 2]>> {
        match s {
            Shape::Polygon(r) => vec![r],
            Shape::MultiPolygon(v) => v.iter().collect(),
        }
    }

    /// The census's checks on one result; `raw` are the vertices it was made from.
    fn assert_clean(s: &Shape, raw: &[GeoPoint], what: &str) {
        if raw.iter().all(|v| v.lat.abs() == 90.0) {
            return;
        }
        for r in rings(s) {
            assert_eq!(r.first(), r.last(), "{what}: ring not closed");
            assert!(
                area(r) > 0.0,
                "{what}: area {} not positive: {r:?} from {raw:?}",
                area(r)
            );
            for p in r {
                assert!(
                    (-180.0..=180.0).contains(&p[0]) && (-90.0..=90.0).contains(&p[1]),
                    "{what}: {p:?}"
                );
            }
            let body = &r[..r.len() - 1];
            for (i, p) in body.iter().enumerate() {
                assert!(
                    !body[i + 1..].contains(p),
                    "{what}: {p:?} repeated in {r:?}"
                );
            }
            let flat: Vec<(f64, f64)> = r
                .windows(2)
                .filter(|w| on_pole(&w[0]) && w[0][1] == w[1][1])
                .map(|w| (w[0][0].min(w[1][0]), w[0][0].max(w[1][0])))
                .collect();
            for (i, a) in flat.iter().enumerate() {
                for b in &flat[i + 1..] {
                    assert!(
                        a.0.max(b.0) >= a.1.min(b.1),
                        "{what}: edges along the pole overlap in {r:?}"
                    );
                }
            }
            for p in r {
                // Only the cut at 180 is new; every other longitude is the engine's, up to the
                // rounding of a shift by 360.
                assert!(
                    p[0].abs() == 180.0 || raw.iter().any(|v| (v.lon - p[0]).abs() < 1e-12),
                    "{what}: {p:?} is not a longitude of {raw:?}"
                );
            }
            for w in r.windows(2) {
                assert!(
                    (w[0][0] - w[1][0]).abs() <= 180.0 || (on_pole(&w[0]) && on_pole(&w[1])),
                    "{what}: long edge {:?} {:?}",
                    w[0],
                    w[1]
                );
            }
        }
        // The comparison stands only where the ring is not altered at a pole: a vertex on it, or
        // an edge over it (half a turn of longitude long).
        let n = raw.len();
        let altered = (0..n).any(|i| {
            let d = (raw[(i + 1) % n].lon - raw[i].lon).abs();
            raw[i].lat.abs() == 90.0 || (d - 180.0).abs() < HALF_TURN_TOL
        });
        if let (Shape::MultiPolygon(v), false) = (s, altered) {
            let mut u: Vec<[f64; 2]> = Vec::new();
            for v in raw {
                let lon = match u.last() {
                    None => v.lon,
                    Some(q) => q[0] + wrap(v.lon - q[0]),
                };
                u.push([lon, v.lat]);
            }
            let total: f64 = v.iter().map(|r| area(r)).sum();
            let whole = area(&closed(dedup(u))).abs();
            assert!(
                // The floor is the rounding of coordinates near 180 on cells of a millionth of a
                // degree.
                (total - whole).abs() <= 1e-9 * whole.abs() + 1e-12,
                "{what}: {total} vs {whole}"
            );
        }
    }

    #[test]
    fn a_plain_square_is_only_closed() {
        let s = clean(&[gp(0.0, 0.0), gp(10.0, 0.0), gp(10.0, 10.0), gp(0.0, 10.0)]);
        assert_eq!(
            s,
            Shape::Polygon(vec![
                [0.0, 0.0],
                [10.0, 0.0],
                [10.0, 10.0],
                [0.0, 10.0],
                [0.0, 0.0]
            ])
        );
    }

    #[test]
    fn squares_across_the_antimeridian_are_split() {
        // Eastward across 180, anticlockwise: lat 0 from 170 to -170 then back at lat 10.
        let s = clean(&[
            gp(170.0, 0.0),
            gp(-170.0, 0.0),
            gp(-170.0, 10.0),
            gp(170.0, 10.0),
        ]);
        let Shape::MultiPolygon(v) = &s else {
            panic!("{s:?}")
        };
        assert_eq!(v.len(), 2);
        let near = v.iter().find(|r| r.iter().any(|p| p[0] == 170.0)).unwrap();
        let far = v.iter().find(|r| r.iter().any(|p| p[0] == -170.0)).unwrap();
        assert!(near.iter().all(|p| p[0] >= 170.0 && p[0] <= 180.0));
        assert!(far.iter().all(|p| p[0] <= -170.0 && p[0] >= -180.0));
        assert!(near.contains(&[180.0, 0.0]) && near.contains(&[180.0, 10.0]));
        assert!(far.contains(&[-180.0, 0.0]) && far.contains(&[-180.0, 10.0]));
        assert_clean(
            &s,
            &[
                gp(170.0, 0.0),
                gp(-170.0, 0.0),
                gp(-170.0, 10.0),
                gp(170.0, 10.0),
            ],
            "east",
        );
        // The other way: anticlockwise westward across -180.
        let raw = [
            gp(-170.0, 0.0),
            gp(170.0, 0.0),
            gp(170.0, -10.0),
            gp(-170.0, -10.0),
        ];
        let s = clean(&raw);
        assert!(
            matches!(&s, Shape::MultiPolygon(v) if v.len() == 2),
            "{s:?}"
        );
        assert_clean(&s, &raw, "west");
        // Reversed order is clockwise: it is turned about.
        let raw = [
            gp(170.0, 0.0),
            gp(170.0, 10.0),
            gp(-170.0, 10.0),
            gp(-170.0, 0.0),
        ];
        let s = clean(&raw);
        assert!(
            matches!(&s, Shape::MultiPolygon(v) if v.len() == 2),
            "{s:?}"
        );
        assert_clean(&s, &raw, "clockwise");
    }

    #[test]
    fn a_slanted_edge_is_cut_at_its_interpolated_latitude() {
        let raw = [
            gp(170.0, 0.0),
            gp(-170.0, 10.0),
            gp(-170.0, 20.0),
            gp(170.0, 20.0),
        ];
        let Shape::MultiPolygon(v) = clean(&raw) else {
            panic!()
        };
        assert!(v.iter().flatten().any(|p| *p == [180.0, 5.0]));
        assert!(v.iter().flatten().any(|p| *p == [-180.0, 5.0]));
    }

    fn hexagon(lat: f64) -> Vec<GeoPoint> {
        (0..6)
            .map(|k| gp(-180.0 + 60.0 * k as f64 + 30.0, lat))
            .collect()
    }

    #[test]
    fn a_hexagon_about_the_north_pole_is_closed_about_it() {
        let raw = hexagon(80.0);
        let s = clean(&raw);
        let Shape::Polygon(r) = &s else {
            panic!("{s:?}")
        };
        assert!(r.contains(&[180.0, 90.0]) && r.contains(&[-180.0, 90.0]));
        assert_clean(&s, &raw, "north");
    }

    #[test]
    fn a_hexagon_about_the_south_pole_is_closed_about_it() {
        let mut raw = hexagon(-80.0);
        raw.reverse();
        let s = clean(&raw);
        let Shape::Polygon(r) = &s else {
            panic!("{s:?}")
        };
        assert!(r.contains(&[180.0, -90.0]) && r.contains(&[-180.0, -90.0]));
        assert_clean(&s, &raw, "south");
    }

    #[test]
    fn a_vertex_on_the_pole_becomes_two_points() {
        let raw = [gp(0.0, 80.0), gp(40.0, 90.0), gp(100.0, 80.0)];
        let s = clean(&raw);
        let Shape::Polygon(r) = &s else {
            panic!("{s:?}")
        };
        assert!(
            r.contains(&[0.0, 90.0]) && r.contains(&[100.0, 90.0]),
            "{r:?}"
        );
        assert_eq!(r.first(), r.last());
    }

    #[test]
    fn an_edge_over_the_pole_takes_in_the_cap() {
        let raw = [gp(90.0, 60.0), gp(-90.0, 60.0), gp(0.0, 30.0)];
        let s = clean(&raw);
        let Shape::Polygon(r) = &s else {
            panic!("{s:?}")
        };
        assert!(
            r.contains(&[90.0, 90.0]) && r.contains(&[-90.0, 90.0]),
            "{r:?}"
        );
        assert!((area(r) - 8100.0).abs() < 1e-9, "{}", area(r));
    }

    #[test]
    fn a_run_of_pole_vertices_is_one_edge_along_the_pole() {
        let raw = [
            gp(0.0, 80.0),
            gp(20.0, 90.0),
            gp(60.0, 90.0),
            gp(100.0, 80.0),
        ];
        let s = clean(&raw);
        assert_clean(&s, &raw, "run");
        let Shape::Polygon(r) = &s else {
            panic!("{s:?}")
        };
        assert_eq!(r.iter().filter(|p| p[1] == 90.0).count(), 2, "{r:?}");
    }

    #[test]
    fn a_coordinate_that_is_not_finite_does_not_panic() {
        for bad in [f64::NAN, f64::INFINITY] {
            let raw = [gp(bad, 0.0), gp(10.0, 0.0), gp(10.0, 10.0)];
            assert!(matches!(clean(&raw), Shape::Polygon(r) if r.len() == 4));
        }
    }

    #[test]
    fn repeated_vertices_are_dropped() {
        let raw = [
            gp(0.0, 0.0),
            gp(0.0, 0.0),
            gp(10.0, 0.0),
            gp(10.0, 10.0),
            gp(10.0, 10.0),
        ];
        let Shape::Polygon(r) = clean(&raw) else {
            panic!()
        };
        assert_eq!(r, vec![[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 0.0]]);
    }

    #[test]
    fn a_degenerate_pole_ring_is_written_as_it_comes() {
        let raw = [gp(0.0, 90.0), gp(0.0, 90.0), gp(0.0, 90.0)];
        assert_eq!(clean(&raw), Shape::Polygon(vec![[0.0, 90.0]; 3]));
    }

    fn aperture(name: &str) -> u64 {
        if matches!(name, "IGEO7" | "IVEA7H" | "RTEA7H") {
            7
        } else {
            3
        }
    }

    /// Every zone of `grid` at `res`, by a breadth-first walk over the neighbours.
    fn walk(grid: &rs4dggs::AnyGrid, res: u8) -> Vec<rs4dggs::ZoneId> {
        let first = grid.zone_from_geo(0.0, 0.0, res).unwrap();
        let mut seen = HashSet::from([first]);
        let mut queue = VecDeque::from([first]);
        let mut all = Vec::new();
        while let Some(z) = queue.pop_front() {
            all.push(z);
            for n in grid.neighbors(z) {
                if seen.insert(n) {
                    queue.push_back(n);
                }
            }
        }
        all
    }

    #[test]
    fn every_zone_to_resolution_3_is_clean() {
        for name in rs4dggs::GRID_NAMES {
            let grid = rs4dggs::get_grid(name).unwrap();
            let mut any_multi = false;
            for res in 0..=3 {
                let zones = walk(&grid, res);
                assert_eq!(
                    zones.len() as u64,
                    10 * aperture(name).pow(u32::from(res)) + 2,
                    "{name} {res}"
                );
                let (mut poly, mut multi) = (0, 0);
                let mut tiled = 0.0;
                for &id in &zones {
                    let raw = grid.vertices(id);
                    let s = clean(&raw);
                    assert_clean(&s, &raw, &format!("{name} {res} {}", grid.text_id(id)));
                    tiled += rings(&s).into_iter().map(|r| area(r)).sum::<f64>();
                    match s {
                        Shape::Polygon(_) => poly += 1,
                        Shape::MultiPolygon(_) => multi += 1,
                    }
                }
                assert!(
                    (tiled - 64800.0).abs() < 1e-6,
                    "{name} {res}: the pieces sum to {tiled}, not 64800"
                );
                let mut poles = String::new();
                for lat in [90.0, -90.0] {
                    let id = grid.zone_from_geo(lat, 0.0, res).unwrap();
                    let shape = clean(&grid.vertices(id));
                    let on: Vec<[f64; 2]> = rings(&shape)
                        .into_iter()
                        .flatten()
                        .copied()
                        .filter(|p| p[1] == lat)
                        .collect();
                    assert!(!on.is_empty(), "{name} {res}: no point on latitude {lat}");
                    let raw = grid.vertices(id);
                    let n = raw.len();
                    let kind = if raw.iter().any(|p| p.lat == lat) {
                        "vertex"
                    } else if (0..n).any(|i| {
                        ((raw[(i + 1) % n].lon - raw[i].lon).abs() - 180.0).abs() < HALF_TURN_TOL
                    }) {
                        "edge over the pole"
                    } else {
                        "surrounds"
                    };
                    poles.push_str(&format!(
                        " [{kind}, {}]",
                        if matches!(shape, Shape::Polygon(_)) {
                            "Polygon"
                        } else {
                            "MultiPolygon"
                        }
                    ));
                }
                println!("{name} res {res}: {poly} Polygon, {multi} MultiPolygon; poles:{poles}");
                any_multi |= multi > 0;
            }
            assert!(any_multi, "{name}: no MultiPolygon at any resolution");
        }
    }

    #[test]
    fn the_polar_cells_at_the_finest_resolutions() {
        let grid = rs4dggs::get_grid("IGEO7").unwrap();
        for res in 15..=19 {
            let north = grid.zone_from_geo(90.0, 0.0, res).unwrap();
            let mut ids = vec![north, grid.zone_from_geo(89.999999, 0.0, res).unwrap()];
            ids.extend(grid.neighbors(north));
            for id in ids {
                let raw = grid.vertices(id);
                let s = clean(&raw);
                assert_clean(&s, &raw, &format!("IGEO7 {res} {}", grid.text_id(id)));
            }
        }
    }

    /// The zones of `id` and of its neighbours and theirs.
    fn disk2(grid: &rs4dggs::AnyGrid, id: rs4dggs::ZoneId) -> Vec<rs4dggs::ZoneId> {
        let mut all = vec![id];
        for depth in 0..2 {
            let from = all.clone();
            for z in from {
                for n in grid.neighbors(z) {
                    if !all.contains(&n) {
                        all.push(n);
                    }
                }
            }
            let _ = depth;
        }
        all
    }

    #[test]
    fn the_zone_at_each_pole_reaches_it_at_every_resolution() {
        for name in rs4dggs::GRID_NAMES {
            let grid = rs4dggs::get_grid(name).unwrap();
            for res in 0..=grid.max_resolution() {
                for lat in [90.0, -90.0] {
                    let id = grid.zone_from_geo(lat, 0.0, res).unwrap();
                    let shape = clean(&grid.vertices(id));
                    assert!(
                        rings(&shape).into_iter().flatten().any(|p| p[1] == lat),
                        "{name} {res}: the zone {} does not reach latitude {lat}",
                        grid.text_id(id)
                    );
                    for z in disk2(&grid, id) {
                        let raw = grid.vertices(z);
                        assert_clean(
                            &clean(&raw),
                            &raw,
                            &format!("{name} {res} {}", grid.text_id(z)),
                        );
                    }
                }
            }
        }
    }

    /// A zone by name: simple, anticlockwise, in range, and reaching a pole or split in two.
    fn pinned(grid: &str, text: &str, pole: bool) -> Shape {
        let grid = rs4dggs::get_grid(grid).unwrap();
        let raw = grid.vertices(grid.zone_from_text(text).unwrap());
        let s = clean(&raw);
        assert_clean(&s, &raw, text);
        if pole {
            assert!(
                rings(&s).into_iter().flatten().any(|p| p[1].abs() == 90.0),
                "{text}: no point on a pole"
            );
        }
        s
    }

    #[test]
    fn the_pole_zones_at_resolution_1_by_name() {
        for (grid, text) in [
            ("IGEO7", "013"),
            ("IGEO7", "084"),
            ("ISEA3H", "A3-0-D"),
            ("ISEA3H", "A8-0-C"),
        ] {
            pinned(grid, text, true);
        }
    }

    #[test]
    fn a_zone_on_the_antimeridian_by_name() {
        let Shape::MultiPolygon(v) = pinned("IGEO7", "114", false) else {
            panic!("not a MultiPolygon")
        };
        assert_eq!(v.len(), 2);
        let meets = |x: f64| v.iter().filter(|r| r.iter().any(|p| p[0] == x)).count();
        assert_eq!((meets(180.0), meets(-180.0)), (1, 1));
    }
}
