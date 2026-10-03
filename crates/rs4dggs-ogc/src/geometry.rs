//! A zone's ring made into a GeoJSON geometry that maps draw correctly: cut at the antimeridian
//! and closed about the poles (RFC 7946, section 3.1.9). A pure function, [`clean`].
//!
//! It takes the ring as the library gives it, the plain vertices of `vertices` or the refined
//! ring of `refined_vertices`, whose edges may cross the antimeridian several times. The library
//! gives its rings anticlockwise as seen from outside the sphere, on every grid (its
//! `Grid::vertices` names the few engine rings that are not). The direction is judged again here
//! all the same, as a guard, on the ring unwrapped in longitude, so that a ring crossing the
//! antimeridian or touching a pole is judged by its true shape.
//!
//! A ring that the engine itself draws crossing or touching itself, as it does for some zones of
//! the finest levels beside the poles and beside the vertices of the icosahedron, is cut and
//! closed as any other and is not mended: the shape then holds the same fault, and a strict
//! reader of GeoJSON may refuse it. Of some fourteen thousand zones about the poles and the
//! antimeridian, at every level of the six grids, every ring that does not cross or touch itself
//! once its points on a pole are replaced, as described under [`clean`], gave a valid shape.

use rs4dggs::GeoPoint;

#[cfg(test)]
pub(crate) mod census;

/// One polygon or several, each a closed exterior ring of `[lon, lat]` positions, anticlockwise.
///
/// A ring that crosses the antimeridian is cut there into as many polygons as it has pieces on
/// either side of it, each within -180 to 180 degrees of longitude.
///
/// The type may gain variants in a later release, so that a `match` on it outside this crate
/// needs an arm for the others.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum Shape {
    /// One closed ring.
    Polygon(Vec<[f64; 2]>),
    /// Two closed rings or more, one for each polygon.
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

/// The point where the straight edge from `a` to `b` meets longitude `m`: an end of the edge
/// itself where it lies on `m`, so that two crossings at one vertex fall at one latitude.
fn crossing(a: [f64; 2], b: [f64; 2], m: f64) -> [f64; 2] {
    if b[0] == m {
        return b;
    }
    let t = (m - a[0]) / (b[0] - a[0]);
    [m, a[1] + t * (b[1] - a[1])]
}

/// A ring of `[lon, lat]` positions.
type Ring = Vec<[f64; 2]>;

/// The pieces of `ring` on either side of the meridian `m`, which is 180 or -180: first those
/// on the near side (the side of zero, a point on the meridian among them), then those beyond
/// it. Each piece keeps the ring's direction, is closed along the meridian, and begins at the
/// first of its points on the ring, the crossings put in their places.
///
/// A ring may cross the meridian any even number of times (each crossing takes it to the other
/// side, and it ends where it began), and a refined ring may do so more than twice. Along the
/// meridian the inside of a simple ring lies between the first crossing and the second, by
/// latitude, between the third and the fourth, and so on; a piece that leaves its side by one
/// crossing therefore continues from the crossing paired with it.
///
/// A vertex on the meridian both of whose neighbours lie on the near side is where the ring
/// touches the meridian from the near side and turns back: it is counted beyond, as if the ring
/// crossed there and returned at once, so that the near side is cut there into two pieces that
/// meet at that point, and not left as one piece that touches itself.
///
/// An edge may lie along the meridian, and so may a run of edges, whether the ring then crosses
/// it, touches it from either side, or crosses it elsewhere too. Such a run bounds the inside of
/// the ring on one side only, and is given to that side: on an anticlockwise ring the inside lies
/// to the left, to the west of a run that goes northward and to the east of one that goes
/// southward. The run then continues that side's closing edge along the meridian, and does not
/// turn back along it.
fn split(ring: &[[f64; 2]], m: f64) -> (Vec<Ring>, Vec<Ring>) {
    let beyond = |p: &[f64; 2]| if m > 0.0 { p[0] > m } else { p[0] < m };
    let on = |p: &[f64; 2]| p[0] == m;
    let strictly_near = |p: &[f64; 2]| !beyond(p) && !on(p);
    let n = ring.len();
    // Whether the run of points on the meridian through `i` goes northward, or `None` where it
    // goes neither way (its ends at one latitude, or the whole ring on the meridian).
    let northward = |i: usize| -> Option<bool> {
        let (mut first, mut last) = (i, i);
        for _ in 0..n {
            if !on(&ring[(first + n - 1) % n]) {
                break;
            }
            first = (first + n - 1) % n;
        }
        for _ in 0..n {
            if !on(&ring[(last + 1) % n]) {
                break;
            }
            last = (last + 1) % n;
        }
        let rise = ring[last][1] - ring[first][1];
        if on(&ring[(first + n - 1) % n]) || rise == 0.0 {
            None
        } else {
            Some(rise > 0.0)
        }
    };
    let near: Vec<bool> = (0..n)
        .map(|i| {
            let (prev, p, next) = (&ring[(i + n - 1) % n], &ring[i], &ring[(i + 1) % n]);
            if !on(p) {
                return !beyond(p);
            }
            if on(prev) || on(next) {
                if let Some(north) = northward(i) {
                    // The near side is to the west of 180, and to the east of -180.
                    return north == (m > 0.0);
                }
            }
            !(strictly_near(prev) && strictly_near(next))
        })
        .collect();
    // The ring with its crossings put in, each before the point that ends its edge, marked, and
    // each point with its side (a crossing with the side of the point after it).
    let mut aug: Vec<([f64; 2], bool, bool)> = Vec::with_capacity(n + 8);
    for i in 0..n {
        let prev = (i + n - 1) % n;
        if near[prev] != near[i] {
            aug.push((crossing(ring[prev], ring[i], m), true, near[i]));
        }
        aug.push((ring[i], false, near[i]));
    }
    let k = aug.len();
    let mut cross: Vec<usize> = (0..k).filter(|&i| aug[i].1).collect();
    if cross.is_empty() {
        let whole = vec![dedup(ring.to_vec())];
        return if near[0] {
            (whole, Vec::new())
        } else {
            (Vec::new(), whole)
        };
    }
    // Two crossings at one latitude are a vertex on the meridian at which the ring touches it
    // and turns back, from beyond or (as counted above) from the near side. Each must be paired
    // with the crossing on its own side of that vertex: on an anticlockwise ring, the crossing by
    // which the ring leaves the near side is taken as the lower of the two at -180 and as the
    // higher at 180.
    cross.sort_by(|&a, &b| {
        aug[a].0[1].total_cmp(&aug[b].0[1]).then_with(|| {
            let (ea, eb) = (aug[a].2, aug[b].2);
            if m > 0.0 { eb.cmp(&ea) } else { ea.cmp(&eb) }
        })
    });
    let mut partner = vec![usize::MAX; k];
    for pair in cross.chunks_exact(2) {
        partner[pair[0]] = pair[1];
        partner[pair[1]] = pair[0];
    }
    let mut sides = (Vec::new(), Vec::new());
    for side in [true, false] {
        let mut seen = vec![false; k];
        let mut pieces: Vec<Vec<usize>> = Vec::new();
        for &start in &cross {
            // A crossing by which the ring enters this side: the point after it lies there.
            let next = &aug[(start + 1) % k];
            if seen[start] || next.1 || next.2 != side {
                continue;
            }
            let (mut piece, mut c) = (Vec::new(), start);
            loop {
                seen[c] = true;
                piece.push(c);
                let mut i = (c + 1) % k;
                while !aug[i].1 {
                    piece.push(i);
                    i = (i + 1) % k;
                }
                // The crossing by which it leaves, and along the meridian to its pair.
                piece.push(i);
                c = partner[i];
                if c == usize::MAX || c == start || seen[c] {
                    break;
                }
            }
            let first = (0..piece.len()).min_by_key(|&j| piece[j]).unwrap_or(0);
            piece.rotate_left(first);
            pieces.push(piece);
        }
        pieces.sort_by_key(|p| p[0]);
        let out = if side { &mut sides.0 } else { &mut sides.1 };
        out.extend(
            pieces
                .into_iter()
                .map(|p| dedup(p.into_iter().map(|i| aug[i].0).collect())),
        );
    }
    sides
}

/// `pts` with each longitude within 180 of the one before, and the turn the ring makes in
/// longitude on the way back to its start (0 for a ring that does not surround a pole). An edge
/// from one pole point to another runs along the pole, where longitude has no meaning, and may
/// be taken either way: it is taken the way that leaves the ring with no turn, as the ring then
/// passes through the pole and does not go about it. (The two points of a vertex on the pole lie
/// at the longitudes of its neighbours on the ring, which on a refined ring may be more than half
/// a turn apart.) Each longitude is the engine's own, plus whole turns.
fn unwrapped(pts: &[[f64; 2]]) -> (Vec<[f64; 2]>, f64) {
    let n = pts.len();
    let mut deltas: Vec<f64> = (0..n)
        .map(|i| wrap(pts[(i + 1) % n][0] - pts[i][0]))
        .collect();
    let half = |i: usize| on_pole(&pts[i]) && on_pole(&pts[(i + 1) % n]);
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

/// The ring's `[lon, lat]` positions, and whether it has no area to draw: a coordinate that is
/// not finite, or every point on a pole. Such a ring is written as it comes.
fn positions(ring: &[GeoPoint]) -> (Vec<[f64; 2]>, bool) {
    let raw: Vec<[f64; 2]> = ring.iter().map(|v| [v.lon, v.lat]).collect();
    let as_it_comes = raw.iter().flatten().any(|c| !c.is_finite()) || raw.iter().all(on_pole);
    (raw, as_it_comes)
}

/// The ring as it stands before any cut, and the turn it makes in longitude on the way back to
/// its start (0, or a whole turn for a ring about a pole): its repeated points dropped, a point
/// on a pole put in an edge that runs over one, each run of points on a pole replaced by two
/// points on it, its longitudes continuous and its direction anticlockwise. `None` for a ring
/// with no area to draw, which [`clean`] writes as it comes.
pub(crate) fn prepared(ring: &[GeoPoint]) -> Option<(Vec<[f64; 2]>, f64)> {
    let (raw, as_it_comes) = positions(ring);
    if as_it_comes {
        return None;
    }
    let pts = dedup(raw);
    if pts.len() < 3 {
        return None;
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
    Some((u, turn))
}

/// A zone's ring, as the library gives it, made into one polygon or several that GeoJSON holds
/// and maps draw: closed, anticlockwise, every longitude within -180 to 180 degrees.
///
/// A ring that crosses the antimeridian is cut there, and each piece on either side becomes a
/// polygon of its own. A point on a pole becomes two points on it, at the longitudes of its
/// neighbours, and a ring about a pole is closed along it from -180 to 180 degrees. A ring with
/// no area to draw (a coordinate that is not finite, every point on a pole, fewer than three
/// distinct points) is returned as one ring, closed and otherwise as it comes.
pub fn clean(ring: &[GeoPoint]) -> Shape {
    drawn(ring, None).unwrap_or_else(|| {
        let (raw, as_it_comes) = positions(ring);
        Shape::Polygon(closed(if as_it_comes { raw } else { dedup(raw) }))
    })
}

/// `x` rounded to `decimals` places, as the writers write it.
pub(crate) fn rounded(x: f64, decimals: u8) -> f64 {
    let p = usize::from(decimals);
    format!("{x:.p$}").parse().unwrap_or(x)
}

/// The shape [`clean`] makes of `ring`, or `None` for a ring with no area to draw (a coordinate
/// that is not finite, every point on a pole, fewer than three distinct points). With
/// `decimals`, the ring is rounded to that many places before it is made a shape, so that points
/// the rounding makes equal are one point, and a vertex at which the rounded ring doubles back on
/// itself to within half a rounding step is removed; a ring then left with fewer than three
/// distinct points has no area to draw.
pub(crate) fn drawn(ring: &[GeoPoint], decimals: Option<u8>) -> Option<Shape> {
    let Some(d) = decimals else {
        return cut(ring, None);
    };
    let ring: Vec<GeoPoint> = ring
        .iter()
        .map(|p| GeoPoint {
            lat: rounded(p.lat, d),
            lon: rounded(p.lon, d),
        })
        .collect();
    cut(&ring, decimals)
}

/// Whether the ring, coming from `a` to `b` and going on to `c`, doubles back on itself at `b`
/// to within `tol`: `a` and `c` are one point, or the turn is a reversal (the ring goes back
/// towards where it came from) and the nearer of `a` and `c` lies within `tol` of the line of the
/// longer edge, so that the two edges cannot be told apart at that precision.
fn spike(a: [f64; 2], b: [f64; 2], c: [f64; 2], tol: f64) -> bool {
    let (ab, bc) = ([b[0] - a[0], b[1] - a[1]], [c[0] - b[0], c[1] - b[1]]);
    let cross = (ab[0] * bc[1] - ab[1] * bc[0]).abs();
    let longer = ab[0].hypot(ab[1]).max(bc[0].hypot(bc[1]));
    a == c || (ab[0] * bc[0] + ab[1] * bc[1] < 0.0 && cross <= tol * longer)
}

/// `u`, the ring as [`prepared`] gives it with the turn it makes, without its spikes narrower
/// than `tol`: each vertex at which the ring doubles back on itself to within `tol` is removed,
/// and again, until none is left. With `tol` half the rounding step, such a spike cannot be
/// written at that precision, and the sliver its removal takes is narrower than the step.
fn without_spikes(mut u: Vec<[f64; 2]>, turn: f64, tol: f64) -> Vec<[f64; 2]> {
    // The ring closes one whole turn on (or none), whatever the rounding of the sum.
    let turn = 360.0 * (turn / 360.0).round();
    loop {
        let n = u.len();
        if n < 3 {
            return u;
        }
        let at = |i: usize| -> [f64; 2] {
            let p = u[i % n];
            if i < n { p } else { [p[0] + turn, p[1]] }
        };
        let before = |i: usize| -> [f64; 2] {
            if i == 0 {
                [u[n - 1][0] - turn, u[n - 1][1]]
            } else {
                u[i - 1]
            }
        };
        let Some(i) = (0..n).find(|&i| spike(before(i), at(i), at(i + 1), tol)) else {
            return u;
        };
        u.remove(i);
        // The two neighbours of a vertex removed may now be one point.
        u = dedup(u);
    }
}

/// The shape [`clean`] makes of `ring`, or `None` for a ring with no area to draw. With
/// `decimals`, every point made by the cut (a crossing of the meridian, a point moved by a whole
/// turn) is brought to that many places before repeated points are dropped and any piece left
/// with fewer than three points is set aside, so that the shape is the one written.
/// A crossing so brought to the grid is not proved unable to make a new turn back within a piece
/// (the turns back are removed before the cut); the census at eight decimals, and a sample of
/// some twelve thousand further zones, found none.
fn cut(ring: &[GeoPoint], decimals: Option<u8>) -> Option<Shape> {
    let snapped = |r: Vec<[f64; 2]>| -> Vec<[f64; 2]> {
        match decimals {
            Some(d) => dedup(
                r.into_iter()
                    .map(|p| [rounded(p[0], d), rounded(p[1], d)])
                    .collect(),
            ),
            None => r,
        }
    };
    let (u, turn) = prepared(ring)?;
    let u = if let Some(d) = decimals {
        // Half the rounding step.
        without_spikes(u, turn, 0.5 * 10f64.powi(-i32::from(d)))
    } else {
        u
    };
    if u.len() < 3 {
        return None;
    }
    if turn.abs() < 180.0 {
        let min = u.iter().map(|p| p[0]).fold(f64::INFINITY, f64::min);
        let max = u.iter().map(|p| p[0]).fold(f64::NEG_INFINITY, f64::max);
        if min >= -180.0 && max <= 180.0 {
            return Some(Shape::Polygon(closed(u)));
        }
        let m = if max > 180.0 { 180.0 } else { -180.0 };
        let shift = if m > 0.0 { -360.0 } else { 360.0 };
        let (near, far) = split(&u, m);
        let far = far.into_iter().map(|r| {
            r.into_iter()
                .map(|p| [p[0] + shift, p[1]])
                .collect::<Vec<_>>()
        });
        // A piece of fewer than three points is where the ring only touches the meridian. (A
        // ring of three distinct points or more always leaves a piece of three on one side, so
        // that the case of no piece at all does not arise; were it to, there would be nothing
        // to draw.)
        let mut pieces: Vec<Vec<[f64; 2]>> = near
            .into_iter()
            .chain(far)
            .map(snapped)
            .filter(|r| r.len() >= 3)
            .map(closed)
            .collect();
        return match pieces.len() {
            0 => None,
            1 => Some(Shape::Polygon(pieces.remove(0))),
            _ => Some(Shape::MultiPolygon(pieces)),
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
    // A ring about a pole crosses every meridian; one that seems to cross none has longitudes so
    // great that a turn is lost in their rounding, and is returned as it stands.
    let Some(i) = (0..path.len() - 1).find(|&i| {
        let (a, b) = (path[i][0], path[i + 1][0]);
        a != b && a.min(b) <= m && m <= a.max(b)
    }) else {
        return Some(Shape::Polygon(closed(u)));
    };
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
    Some(Shape::Polygon(closed(dedup(snapped(ring)))))
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

    // The refined rings: the engine's own refinement of each edge.

    use super::census::{as_written, valid};

    /// The refined ring of a zone by name, and the zone's extent.
    fn refined(grid: &str, text: &str) -> (Vec<GeoPoint>, rs4dggs::Extent) {
        let grid = rs4dggs::get_grid(grid).unwrap();
        let id = grid.zone_from_text(text).unwrap();
        (
            grid.refined_vertices(id, 0).unwrap(),
            grid.extent(id).unwrap(),
        )
    }

    #[test]
    fn a_zone_with_the_north_pole_as_a_vertex_is_one_polygon_that_touches_it() {
        let (raw, e) = refined("IGEO7", "005");
        assert_eq!(raw.len(), 145);
        let on: Vec<usize> = (0..raw.len()).filter(|&i| raw[i].lat == 90.0).collect();
        assert_eq!(on, [117]);
        // Its neighbours on the ring lie a little more than half a turn apart.
        assert_eq!(
            format!("{:.6} {:.6}", raw[116].lon, raw[118].lon),
            "63.630949 -116.370269"
        );
        assert_eq!(
            format!("{:.6} {:.6}", e.ll.lon, e.ur.lon),
            "-116.633082 63.630949"
        );
        let s = clean(&raw);
        let Shape::Polygon(r) = &s else {
            panic!("{s:?}")
        };
        valid(&s).unwrap();
        for p in r {
            assert!(e.ll.lon <= p[0] && p[0] <= e.ur.lon, "{p:?} outside {e:?}");
        }
        assert_eq!(r.iter().filter(|p| p[1] == 90.0).count(), 2, "{r:?}");
    }

    #[test]
    fn the_other_zones_of_level_1_with_a_pole_as_a_vertex_are_valid() {
        let (raw, _) = refined("IGEO7", "114");
        let s = clean(&raw);
        assert!(matches!(s, Shape::MultiPolygon(_)), "{s:?}");
        valid(&s).unwrap();
        for text in ["013", "084"] {
            let (raw, _) = refined("IVEA7H", text);
            let s = clean(&raw);
            valid(&s).unwrap_or_else(|e| panic!("{text}: {e}"));
        }
    }

    #[test]
    fn a_zone_with_the_pole_on_an_edge_keeps_its_shape() {
        let (raw, _) = refined("ISEA3H", "A8-0-C");
        assert_eq!(raw.len(), 500);
        let n = raw.len();
        let i = (0..n).find(|&i| raw[i].lat.abs() == 90.0).unwrap();
        // Its neighbours on the ring lie half a turn apart, to some 1e-11 degrees.
        let (a, b) = (raw[(i + n - 1) % n], raw[(i + 1) % n]);
        assert!(((a.lon - b.lon).abs() - 180.0).abs() < 1e-9, "{a:?} {b:?}");
        let s = clean(&raw);
        assert!(matches!(s, Shape::MultiPolygon(_)), "{s:?}");
        valid(&s).unwrap();
    }

    #[test]
    fn a_refined_ring_that_crosses_the_antimeridian_four_times_is_cut_in_three() {
        let (raw, _) = refined("IGEO7", "01012");
        assert_eq!(raw.len(), 90);
        let lons = raw
            .iter()
            .map(|p| if p.lon < 0.0 { p.lon + 360.0 } else { p.lon });
        let (lo, hi) = lons.fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), x| {
            (lo.min(x), hi.max(x))
        });
        assert_eq!(format!("{lo:.4} {hi:.4}"), "172.0025 180.3305");
        let s = clean(&raw);
        let Shape::MultiPolygon(v) = &s else {
            panic!("{s:?}")
        };
        assert_eq!(v.len(), 3, "{v:?}");
        valid(&s).unwrap();
        // One piece on this side of the meridian, and the two beyond it, each between two of
        // the four crossings.
        assert!(v[0].iter().all(|p| p[0] > 0.0));
        let cut = |r: &Vec<[f64; 2]>| -> Vec<f64> {
            let mut l: Vec<f64> = r[..r.len() - 1]
                .iter()
                .filter(|p| p[0] == -180.0)
                .map(|p| p[1])
                .collect();
            l.sort_by(f64::total_cmp);
            l
        };
        assert_eq!(cut(&v[1]), [63.24971908462399, 64.45383071289578]);
        assert_eq!(cut(&v[2]), [64.84881572393388, 65.58641528280835]);
        let mut near: Vec<f64> = v[0]
            .iter()
            .filter(|p| p[0] == 180.0)
            .map(|p| p[1])
            .collect();
        near.sort_by(f64::total_cmp);
        near.dedup();
        assert_eq!(
            near,
            [
                63.24971908462399,
                64.45383071289578,
                64.84881572393388,
                65.58641528280835
            ]
        );
        assert!(v[1..].iter().flatten().all(|p| p[0] < 0.0));
    }

    /// A ring of `[lon, lat]` pairs as the library gives one.
    fn ring(points: &[[f64; 2]]) -> Vec<GeoPoint> {
        points.iter().map(|p| gp(p[0], p[1])).collect()
    }

    #[test]
    fn a_vertex_exactly_on_the_antimeridian_is_cut_in_two() {
        // Eastward across 180 with a vertex on it, and westward across -180 with one.
        for raw in [
            ring(&[
                [170.0, 0.0],
                [180.0, 0.0],
                [-170.0, 0.0],
                [-170.0, 10.0],
                [170.0, 10.0],
            ]),
            ring(&[
                [-170.0, 0.0],
                [-170.0, 10.0],
                [-180.0, 10.0],
                [170.0, 10.0],
                [170.0, 0.0],
            ]),
            ring(&[
                [170.0, 0.0],
                [-170.0, 0.0],
                [-170.0, 10.0],
                [180.0, 5.0],
                [170.0, 10.0],
            ]),
        ] {
            let s = clean(&raw);
            assert!(
                matches!(&s, Shape::MultiPolygon(v) if v.len() == 2),
                "{s:?}"
            );
            valid(&s).unwrap();
        }
    }

    #[test]
    fn a_ring_that_touches_the_antimeridian_without_crossing_it_is_one_polygon() {
        // Beyond 180 but for its first vertex, which lies on it.
        let raw = ring(&[[180.0, 0.0], [-170.0, 0.0], [-170.0, 10.0]]);
        let s = clean(&raw);
        assert_eq!(
            s,
            Shape::Polygon(vec![
                [-180.0, 0.0],
                [-170.0, 0.0],
                [-170.0, 10.0],
                [-180.0, 0.0]
            ])
        );
        // Within 180 and touching it: nothing to cut.
        let raw = ring(&[[170.0, 0.0], [180.0, 5.0], [170.0, 10.0]]);
        assert!(matches!(clean(&raw), Shape::Polygon(r) if r.len() == 4));
    }

    #[test]
    fn a_ring_that_touches_the_antimeridian_from_the_near_side_is_cut_there() {
        // Across 180 between 0 and 10 N, and back to touch it at 5 N from the west: the western
        // side is two pieces that meet at that point, the eastern one piece.
        let points = [
            [170.0, 0.0],
            [-170.0, 0.0],
            [-170.0, 10.0],
            [170.0, 10.0],
            [180.0, 5.0],
        ];
        for start in 0..points.len() {
            let mut p = points.to_vec();
            p.rotate_left(start);
            let s = clean(&ring(&p));
            valid(&s).unwrap_or_else(|e| panic!("begun at {start}: {s:?}: {e}"));
            let Shape::MultiPolygon(v) = &s else {
                panic!("begun at {start}: {s:?}")
            };
            assert_eq!(v.len(), 3, "begun at {start}: {s:?}");
            for r in v {
                let west = r.iter().all(|q| q[0] >= 170.0);
                let east = r.iter().all(|q| q[0] <= -170.0);
                assert!(west || east, "begun at {start}: {r:?}");
            }
        }
    }

    #[test]
    fn a_ring_that_touches_the_antimeridian_from_beyond_is_cut_within_range() {
        // A ring between 170 W and 170 E across the antimeridian, notched back to touch it at
        // 5 N from the far side, the touching vertex written as 180 and as -180; and its mirror
        // across the prime meridian. Two crossings meet at 5 N; each ring is begun at every one
        // of its points. At 1.01 N the edge from 6 N ends where its interpolation, 6 + (1.01 - 6),
        // does not fall: the crossing there must be the vertex itself.
        for (t, lat, after) in [(180.0, 5.0, 4.0), (-180.0, 5.0, 4.0), (180.0, 1.01, 0.5)] {
            let west = [
                [-170.0, 0.0],
                [-170.0, 10.0],
                [170.0, 10.0],
                [170.0, 6.0],
                [t, lat],
                [170.0, after],
                [170.0, 0.0],
            ];
            // Mirrored, and taken backwards so that it stays anticlockwise.
            let east: Vec<[f64; 2]> = west.iter().rev().map(|p| [-p[0], p[1]]).collect();
            for (which, points) in [("west", west.to_vec()), ("east", east)] {
                for start in 0..points.len() {
                    let mut p = points.clone();
                    p.rotate_left(start);
                    let s = clean(&ring(&p));
                    let what = format!("{which}, {t} at {lat}, begun at {start}: {s:?}");
                    valid(&s).unwrap_or_else(|e| panic!("{what}: {e}"));
                    assert!(
                        rings(&s)
                            .iter()
                            .all(|r| r.iter().all(|q| q[0].abs() <= 180.0)),
                        "{what}"
                    );
                }
            }
        }
    }

    #[test]
    fn a_comb_whose_three_teeth_cross_the_antimeridian_gives_four_polygons() {
        // A back along 170 E from 0 to 50 N, and three teeth eastward across 180, to 175 W.
        let raw = ring(&[
            [170.0, 0.0],
            [-175.0, 0.0],
            [-175.0, 10.0],
            [175.0, 10.0],
            [175.0, 20.0],
            [-175.0, 20.0],
            [-175.0, 30.0],
            [175.0, 30.0],
            [175.0, 40.0],
            [-175.0, 40.0],
            [-175.0, 50.0],
            [170.0, 50.0],
        ]);
        let s = clean(&raw);
        let Shape::MultiPolygon(v) = &s else {
            panic!("{s:?}")
        };
        assert_eq!(v.len(), 4, "{v:?}");
        valid(&s).unwrap();
        // The near side whole, then the three teeth from the south.
        assert!(v[0].iter().all(|p| p[0] >= 170.0));
        for (k, tooth) in v[1..].iter().enumerate() {
            let lat = 20.0 * k as f64;
            let mut box_: Vec<[f64; 2]> = tooth[..tooth.len() - 1].to_vec();
            box_.sort_by(|a, b| a[0].total_cmp(&b[0]).then(a[1].total_cmp(&b[1])));
            assert_eq!(
                box_,
                [
                    [-180.0, lat],
                    [-180.0, lat + 10.0],
                    [-175.0, lat],
                    [-175.0, lat + 10.0]
                ]
            );
        }
        let area: f64 = v.iter().map(|r| area(r)).sum();
        assert_eq!(area, 10.0 * 50.0 - 2.0 * 5.0 * 10.0 + 3.0 * 5.0 * 10.0);
    }

    #[test]
    fn a_ring_with_an_edge_along_the_antimeridian_is_cut_within_range() {
        // Each ring is anticlockwise, written here with longitudes beyond 180 where it lies past
        // the antimeridian, and holds a run of edges along 180 between two points on either
        // side of it; each is cut into the number of polygons given.
        let cases: [(&str, usize, &[[f64; 2]]); 8] = [
            (
                "across, southward along it",
                2,
                &[
                    [170.0, 0.0],
                    [180.0, 0.0],
                    [180.0, -2.0],
                    [190.0, 0.0],
                    [190.0, 10.0],
                    [170.0, 10.0],
                ],
            ),
            (
                "across, northward along it",
                2,
                &[
                    [170.0, 0.0],
                    [180.0, 0.0],
                    [180.0, 2.0],
                    [190.0, 2.0],
                    [190.0, 10.0],
                    [170.0, 10.0],
                ],
            ),
            (
                "back, northward along it",
                2,
                &[
                    [170.0, 0.0],
                    [190.0, 0.0],
                    [190.0, 10.0],
                    [180.0, 10.0],
                    [180.0, 12.0],
                    [170.0, 12.0],
                ],
            ),
            (
                "back, southward along it",
                2,
                &[
                    [170.0, 0.0],
                    [190.0, 0.0],
                    [190.0, 10.0],
                    [180.0, 10.0],
                    [180.0, 8.0],
                    [170.0, 8.0],
                ],
            ),
            (
                "touching from this side, the inside beyond",
                3,
                &[
                    [170.0, 0.0],
                    [190.0, 0.0],
                    [190.0, 10.0],
                    [170.0, 10.0],
                    [170.0, 6.0],
                    [180.0, 6.0],
                    [180.0, 4.0],
                    [170.0, 4.0],
                ],
            ),
            (
                "touching from this side, the inside here",
                2,
                &[
                    [170.0, 0.0],
                    [190.0, 0.0],
                    [190.0, 10.0],
                    [175.0, 10.0],
                    [175.0, 20.0],
                    [180.0, 20.0],
                    [180.0, 25.0],
                    [180.0, 30.0],
                    [170.0, 30.0],
                ],
            ),
            (
                "touching from beyond, the inside here",
                3,
                &[
                    [170.0, 0.0],
                    [190.0, 0.0],
                    [190.0, 4.0],
                    [180.0, 4.0],
                    [180.0, 6.0],
                    [190.0, 6.0],
                    [190.0, 10.0],
                    [170.0, 10.0],
                ],
            ),
            (
                "touching from beyond, the inside beyond",
                2,
                &[
                    [170.0, 0.0],
                    [190.0, 0.0],
                    [190.0, 30.0],
                    [180.0, 30.0],
                    [180.0, 25.0],
                    [180.0, 20.0],
                    [185.0, 20.0],
                    [185.0, 10.0],
                    [170.0, 10.0],
                ],
            ),
        ];
        for (what, pieces, points) in cases {
            let wrapped: Vec<[f64; 2]> = points.iter().map(|p| [wrap(p[0]), p[1]]).collect();
            // The same about -180: mirrored across the prime meridian and taken backwards, so
            // that it stays anticlockwise; and with its points on the meridian written as -180.
            let mirrored: Vec<[f64; 2]> = wrapped.iter().rev().map(|p| [-p[0], p[1]]).collect();
            let other_sign: Vec<[f64; 2]> = wrapped
                .iter()
                .map(|p| [if p[0] == 180.0 { -180.0 } else { p[0] }, p[1]])
                .collect();
            for (how, points) in [
                ("", wrapped),
                (", mirrored", mirrored),
                (", as -180", other_sign),
            ] {
                for start in 0..points.len() {
                    let mut p = points.clone();
                    p.rotate_left(start);
                    let raw = ring(&p);
                    let s = clean(&raw);
                    let what = format!("{what}{how}, begun at {start}: {s:?}");
                    valid(&s).unwrap_or_else(|e| panic!("{what}: {e}"));
                    assert_clean(&s, &raw, &what);
                    assert_eq!(rings(&s).len(), pieces, "{what}");
                }
            }
        }
    }

    #[test]
    fn a_spike_made_by_the_rounding_is_removed() {
        // A square with a sliver three billionths of a degree wide on its eastern side: with
        // eight decimals the sliver is a spike that goes out to 15 E and back the same way.
        let points = [
            [0.0, 0.0],
            [10.0, 0.0],
            [10.0, 5.0],
            [15.0, 5.000000004],
            [10.0, 5.000000003],
            [10.0, 10.0],
            [0.0, 10.0],
        ];
        for start in 0..points.len() {
            let mut p = points.to_vec();
            p.rotate_left(start);
            let raw = ring(&p);
            valid(&clean(&raw)).unwrap_or_else(|e| panic!("begun at {start}, unrounded: {e}"));
            let s = as_written(&drawn(&raw, Some(8)).unwrap(), 8);
            valid(&s).unwrap_or_else(|e| panic!("begun at {start}: {e}: {s:?}"));
            let Shape::Polygon(r) = &s else {
                panic!("begun at {start}: {s:?}")
            };
            assert_eq!(area(r), 100.0, "begun at {start}: {r:?}");
            assert!(!r.iter().any(|q| q[0] == 15.0), "begun at {start}: {r:?}");
        }
    }

    #[test]
    fn zones_beside_a_pole_whose_rounded_ring_doubles_back_give_valid_shapes() {
        // On the RTEA grids' seam at 11.2 E and 168.8 W, beside the poles, the engine's ring
        // holds two points a few billionths of a degree apart, which eight decimals bring onto
        // one line: the ring then goes out and back the same way, and the spike is removed.
        let zones: [(&str, &[&str]); 2] = [
            (
                "RTEA3H",
                &[
                    "P5-5DA0E09CD0BC-A",
                    "P5-5DA0E177C327-A",
                    "P5-5DA0E252B592-A",
                    "Q5-34AA7EFC6B3A0-A",
                    "Q5-34AA7F2578AE1-A",
                ],
            ),
            (
                "RTEA7H",
                &[
                    "0132323232323232320",
                    "013232323232323236",
                    "08454545454545450",
                    "08454545454545454",
                    "0845454545454545450",
                    "0845454545454545454",
                    "08454545454545454541",
                    "08454545454545454545",
                    "084545454545454545450",
                    "1145454545454545454",
                ],
            ),
        ];
        for (grid, texts) in zones {
            for &text in texts {
                let (raw, _) = refined(grid, text);
                valid(&clean(&raw)).unwrap_or_else(|e| panic!("{text}, unrounded: {e}"));
                let s = as_written(&drawn(&raw, Some(8)).unwrap(), 8);
                valid(&s).unwrap_or_else(|e| panic!("{text}: {e}"));
                let rounded_ring: Vec<GeoPoint> = raw
                    .iter()
                    .map(|p| gp(rounded(p.lon, 8), rounded(p.lat, 8)))
                    .collect();
                let (u, turn) = prepared(&rounded_ring).unwrap();
                let (before, after) = (area(&u), area(&without_spikes(u.clone(), turn, 5e-9)));
                // The sliver taken is narrower than half a step: here below 1e-14 square degrees.
                assert!(
                    (before - after).abs() < 1e-14,
                    "{text}: {before} then {after}"
                );
            }
        }
    }

    #[test]
    fn zones_beside_a_pole_whose_rounded_ring_crosses_by_less_than_a_step_give_valid_shapes() {
        // On the same seam, the two points come one step apart in longitude once rounded, and
        // the edges that leave them cross each other some 1e-15 degrees from the second: the
        // ring turns back within half a step, and the vertex at the turn is removed.
        let zones: [(&str, &[&str]); 2] = [
            (
                "RTEA3H",
                &[
                    "M5-20E1023287-A",
                    "M5-20E10A4E78-A",
                    "M5-20E1126A69-A",
                    "N5-127E98D69DE-A",
                    "O5-A67351DC663-A",
                    "O5-A673566C1DC-A",
                    "O5-A6735AFBD55-A",
                ],
            ),
            (
                "RTEA7H",
                &["0845454545454541", "0845454545454545", "1145454545454545"],
            ),
        ];
        for (grid, texts) in zones {
            for &text in texts {
                let (raw, _) = refined(grid, text);
                valid(&clean(&raw)).unwrap_or_else(|e| panic!("{text}, unrounded: {e}"));
                let s = as_written(&drawn(&raw, Some(8)).unwrap(), 8);
                valid(&s).unwrap_or_else(|e| panic!("{text}: {e}"));
                let rounded_ring: Vec<GeoPoint> = raw
                    .iter()
                    .map(|p| gp(rounded(p.lon, 8), rounded(p.lat, 8)))
                    .collect();
                let (u, turn) = prepared(&rounded_ring).unwrap();
                let (before, after) = (area(&u), area(&without_spikes(u.clone(), turn, 5e-9)));
                // The sliver taken is narrower than half a step: here below 1e-14 square degrees.
                assert!(
                    (before - after).abs() < 1e-14,
                    "{text}: {before} then {after}"
                );
            }
        }
    }

    #[test]
    fn two_zones_rounded_coarsely_give_valid_shapes() {
        // Rounded to whole degrees, the ring of the first has an edge along the antimeridian
        // and crosses it elsewhere; rounded to five places, a crossing of the second falls
        // between two places, beside a point of the ring.
        for (grid, text, decimals) in [("ISEA3H", "C9-7-A", 0), ("IVEA3H", "M8-2E8C638937-A", 5)] {
            let (raw, _) = refined(grid, text);
            let s = as_written(&drawn(&raw, Some(decimals)).unwrap(), decimals);
            valid(&s).unwrap_or_else(|e| panic!("{text}: {e}: {s:?}"));
        }
    }

    #[test]
    fn coordinates_finite_but_extreme_do_not_panic() {
        for raw in [
            ring(&[[f64::MAX, 0.0], [-f64::MAX, 0.0], [0.0, 10.0]]),
            ring(&[[1e300, 89.0], [-1e300, 89.0], [0.0, 89.5], [5e-324, 90.0]]),
            ring(&[[0.0, -f64::MAX], [10.0, f64::MAX], [5.0, 0.0]]),
        ] {
            let _ = clean(&raw);
        }
    }
}
