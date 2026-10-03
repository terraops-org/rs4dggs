//! A zone as GeoJSON (RFC 7946): its geometry alone, the zone as a feature, and a list of zones
//! as a feature collection. Every feature carries a number as its `id` and the zone's text
//! identifier as the property `zoneID`, and its members are written in the order `type`, `id`,
//! `geometry`, `properties`. No white space is written, nor a final line feed.

use std::fmt;
use std::io::{self, Write};

use rs4dggs::{AnyGrid, GeoPoint, ZoneId};

use crate::geometry::{self, Shape};
use crate::{Error, Result, json};

/// How a zone's geometry is written: the standard's `geometry` parameter.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ZoneGeometry {
    /// The zone's region: its boundary as the library's refined ring, made a `Polygon`, or a
    /// `MultiPolygon` where the antimeridian cuts it, by [`geometry::clean`].
    Region {
        /// The parts into which each edge of the zone is divided, 0 for the engine's own choice;
        /// one beyond the library's limit is refused with [`Error::Grid`].
        edge_refinement: u32,
    },
    /// The zone's centroid, as a `Point`.
    Centroid,
    /// No geometry: a feature's `geometry` is `null`.
    None,
}

/// How the coordinates are written: by default with eight decimals, about a millimetre on the
/// ground.
///
/// Eight is also the fewest: a caller who asks fewer decimals gets eight. Rounded more coarsely,
/// the rings of fine zones were written invalid in great numbers (with six decimals, 17 to 27 per
/// cent of the aperture-7 zones of levels 17 to 19, and 16 to 23 per cent of the aperture-3 zones
/// of levels 31 to 33, in a sample at middle latitudes; with seven, 19 per cent of the aperture-7
/// zones of level 19), the rounding making their rings cross or touch themselves.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Coordinates {
    /// The number of decimals of every coordinate, never fewer than eight (a smaller number is
    /// taken as eight), or `None` for the shortest form that reads back to the same double.
    ///
    /// With decimals, the ring is rounded before it is made a polygon, so that points the
    /// rounding makes equal are one point, and a vertex at which the rounded ring turns back on
    /// itself to within half a rounding step is removed, taking a sliver narrower than the
    /// precision asked; a zone whose ring is then left with fewer than three distinct points has
    /// no region to draw. On a sample of some fourteen thousand zones about the poles, the
    /// antimeridian and a vertex of the icosahedron, the only shapes a strict reader refuses at
    /// eight decimals are those of rings the engine itself draws crossing, and twelve zones of the
    /// finest aperture-3 levels at the poles are written without geometry, their rounded ring
    /// enclosing nothing.
    pub decimals: Option<u8>,
}

/// The fewest decimals with which a coordinate is written.
const MIN_DECIMALS: u8 = 8;

impl Default for Coordinates {
    /// Eight decimals.
    fn default() -> Self {
        Self {
            decimals: Some(MIN_DECIMALS),
        }
    }
}

impl Coordinates {
    /// The decimals with which the coordinates are written: those asked, but never fewer than
    /// eight; `None` for the shortest form.
    fn applied(self) -> Option<u8> {
        self.decimals.map(|d| d.max(MIN_DECIMALS))
    }
}

/// What a list of zones as GeoJSON holds, as [`ZoneListGeoJson::finish`] counts it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Written {
    /// The features written, one for each zone given.
    pub features: u64,
    /// The features written with `"geometry":null`: every feature under [`ZoneGeometry::None`],
    /// and otherwise those of the zones that cannot be drawn.
    pub without_geometry: u64,
    /// The null zones passed over. A list of zones refuses the null zone and leaves this at
    /// nought; a writer that walks a sub-zone order itself passes them over and counts them here.
    pub null_zones_skipped: u64,
}

/// A geometry made, before any byte of it is written.
enum Drawn {
    Point([f64; 2]),
    Shape(Shape),
}

/// The geometry of `zone`, or `None` where there is none to write: under [`ZoneGeometry::None`],
/// or for a zone that cannot be drawn (a ring that is empty, wholly on a pole, of fewer than three
/// distinct points, or with a coordinate that is not finite; a centroid that is not finite).
fn drawn(
    grid: AnyGrid,
    zone: ZoneId,
    geometry: ZoneGeometry,
    coords: Coordinates,
) -> Result<Option<Drawn>> {
    if zone == ZoneId::NULL {
        return Err(Error::NullZone);
    }
    Ok(match geometry {
        ZoneGeometry::None => None,
        ZoneGeometry::Centroid => {
            let c = grid.centroid(zone);
            let p = [c.lon, c.lat];
            p.iter().all(|x| x.is_finite()).then_some(Drawn::Point(p))
        }
        ZoneGeometry::Region { edge_refinement } => {
            return region(grid.refined_vertices(zone, edge_refinement)?, coords);
        }
    })
}

/// The region of a zone whose ring is `ring`, or `None` where it cannot be drawn; refused with
/// [`Error::NotFinite`] where the cleaning, from finite coordinates, makes one that is not, which
/// no ring of the library's grids does.
fn region(ring: Vec<GeoPoint>, coords: Coordinates) -> Result<Option<Drawn>> {
    let Some(shape) = geometry::drawn(&ring, coords.applied()) else {
        return Ok(None);
    };
    let finite = |r: &Vec<[f64; 2]>| r.iter().flatten().all(|x| x.is_finite());
    let all_finite = match &shape {
        Shape::Polygon(r) => finite(r),
        Shape::MultiPolygon(rs) => rs.iter().all(finite),
    };
    if !all_finite {
        return Err(Error::NotFinite);
    }
    Ok(Some(Drawn::Shape(shape)))
}

fn number<W: Write>(out: &mut W, x: f64, coords: Coordinates) -> io::Result<()> {
    match coords.applied() {
        Some(d) => json::fixed(out, x, d),
        None => json::f64(out, x),
    }
}

fn position<W: Write>(out: &mut W, p: [f64; 2], coords: Coordinates) -> io::Result<()> {
    out.write_all(b"[")?;
    number(out, p[0], coords)?;
    out.write_all(b",")?;
    number(out, p[1], coords)?;
    out.write_all(b"]")
}

fn ring<W: Write>(out: &mut W, ring: &[[f64; 2]], coords: Coordinates) -> io::Result<()> {
    out.write_all(b"[")?;
    for (i, &p) in ring.iter().enumerate() {
        if i > 0 {
            out.write_all(b",")?;
        }
        position(out, p, coords)?;
    }
    out.write_all(b"]")
}

fn write_drawn<W: Write>(out: &mut W, drawn: &Drawn, coords: Coordinates) -> io::Result<()> {
    match drawn {
        Drawn::Point(p) => {
            out.write_all(br#"{"type":"Point","coordinates":"#)?;
            position(out, *p, coords)?;
        }
        Drawn::Shape(Shape::Polygon(r)) => {
            out.write_all(br#"{"type":"Polygon","coordinates":["#)?;
            ring(out, r, coords)?;
            out.write_all(b"]")?;
        }
        Drawn::Shape(Shape::MultiPolygon(rs)) => {
            out.write_all(br#"{"type":"MultiPolygon","coordinates":["#)?;
            for (i, r) in rs.iter().enumerate() {
                out.write_all(if i > 0 { b",[" } else { b"[" })?;
                ring(out, r, coords)?;
                out.write_all(b"]")?;
            }
            out.write_all(b"]")?;
        }
    }
    out.write_all(b"}")
}

/// Writes `zone` as a feature; its properties begin with what `properties` writes (members, each
/// followed by a comma) and end with `zoneID`.
fn write_feature<W: Write>(
    out: &mut W,
    grid: AnyGrid,
    zone: ZoneId,
    id: u64,
    drawn: Option<&Drawn>,
    coords: Coordinates,
    properties: impl FnOnce(&mut W) -> io::Result<()>,
) -> io::Result<()> {
    write!(out, r#"{{"type":"Feature","id":{id},"geometry":"#)?;
    match drawn {
        Some(d) => write_drawn(out, d, coords)?,
        None => out.write_all(b"null")?,
    }
    out.write_all(br#","properties":{"#)?;
    properties(out)?;
    out.write_all(br#""zoneID":"#)?;
    json::string(out, &grid.text_id(zone))?;
    out.write_all(b"}}")
}

/// Writes the geometry of `zone` as a GeoJSON geometry object and answers `true`; where there is
/// none to write, under [`ZoneGeometry::None`] or for a zone that cannot be drawn, it writes
/// nothing at all and answers `false`.
///
/// It writes a member of the caller's own document, such as the `geometry` of a zone's
/// information, which the caller leaves out when the answer is `false`. A zone cannot be drawn
/// when its ring is empty, lies wholly on a pole, has fewer than three distinct points (with
/// `decimals`, once rounded) or a coordinate that is not finite: on ISEA7H_Z7 at its finest
/// level, a few zones beside a pole have a ring that is a single point.
///
/// The coordinates are written as `coords` says: with eight decimals by default and never fewer,
/// or in the shortest form that reads back to the same double (see [`Coordinates`]).
///
/// The null zone is refused with [`Error::NullZone`], a refinement beyond the library's limit
/// with [`Error::Grid`], and a geometry with a coordinate that is not finite with
/// [`Error::NotFinite`], before anything is written. The sink is not flushed.
///
/// # Examples
///
/// ```
/// use rs4dggs_ogc::{Coordinates, ZoneGeometry, write_zone_geometry};
///
/// let grid = rs4dggs::get_grid("ISEA3H")?;
/// let zone = grid.zone_from_text("E6-317-A")?;
/// let mut out = Vec::new();
/// let coords = Coordinates::default();
/// assert!(write_zone_geometry(&mut out, grid, zone, ZoneGeometry::Centroid, coords)?);
/// assert_eq!(out, br#"{"type":"Point","coordinates":[34.78016915,45.42937742]}"#);
///
/// out.clear();
/// assert!(!write_zone_geometry(&mut out, grid, zone, ZoneGeometry::None, coords)?);
/// assert!(out.is_empty());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn write_zone_geometry<W: Write>(
    mut out: W,
    grid: AnyGrid,
    zone: ZoneId,
    geometry: ZoneGeometry,
    coords: Coordinates,
) -> Result<bool> {
    let Some(d) = drawn(grid, zone, geometry, coords)? else {
        return Ok(false);
    };
    write_drawn(&mut out, &d, coords)?;
    Ok(true)
}

/// Writes `zone` as one GeoJSON feature, with `id` as its `id` and the zone's text identifier as
/// its property `zoneID`, and answers whether it has a geometry: where there is none to write
/// (see [`write_zone_geometry`]), its `geometry` is `null` and the answer is `false`.
///
/// The coordinates are written as `coords` says: with eight decimals by default and never fewer,
/// or in the shortest form that reads back to the same double (see [`Coordinates`]).
///
/// The null zone is refused with [`Error::NullZone`], a refinement beyond the library's limit
/// with [`Error::Grid`], and a geometry with a coordinate that is not finite with
/// [`Error::NotFinite`], before anything is written. The sink is flushed.
pub fn write_zone_feature<W: Write>(
    mut out: W,
    grid: AnyGrid,
    zone: ZoneId,
    id: u64,
    geometry: ZoneGeometry,
    coords: Coordinates,
) -> Result<bool> {
    let d = drawn(grid, zone, geometry, coords)?;
    write_feature(&mut out, grid, zone, id, d.as_ref(), coords, |_| Ok(()))?;
    out.flush()?;
    Ok(d.is_some())
}

/// The writer of a list of zones as GeoJSON, the document a zone query answers with
/// `application/geo+json`: `{"type":"FeatureCollection","features":[...]}`, one feature for each
/// zone, as [`write_zone_feature`] writes it.
///
/// The zones are given one at a time, in the order they are to appear, and each is written as
/// soon as it is given, the collection being opened with the first; the features are numbered
/// from `first_id` on, so that a service that pages a list continues the count of the page
/// before. A zone that cannot be drawn is written with `"geometry":null` and counted, so that the
/// list holds a feature for every zone given.
///
/// A failure of the sink, or a zone the library refuses, in the middle of the list leaves the
/// document truncated where it arose.
///
/// # Examples
///
/// ```
/// use rs4dggs_ogc::{Coordinates, ZoneGeometry, ZoneListGeoJson};
///
/// let grid = rs4dggs::get_grid("ISEA3H")?;
/// let region = ZoneGeometry::Region { edge_refinement: 0 };
/// let mut list = ZoneListGeoJson::new(Vec::new(), grid, region, Coordinates::default(), 1);
/// list.zone(grid.zone_from_text("E6-317-A")?)?;
/// let written = list.finish()?;
/// assert_eq!((written.features, written.without_geometry), (1, 0));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub struct ZoneListGeoJson<W: Write> {
    out: W,
    grid: AnyGrid,
    geometry: ZoneGeometry,
    coords: Coordinates,
    next_id: u64,
    written: Written,
}

impl<W: Write> ZoneListGeoJson<W> {
    /// A list of zones of `grid`, to be written to `out`, the first feature numbered `first_id`.
    /// Nothing is written yet. The coordinates are written as `coords` says: with eight decimals
    /// by default and never fewer, or in the shortest form that reads back to the same double
    /// (see [`Coordinates`]).
    pub fn new(
        out: W,
        grid: AnyGrid,
        geometry: ZoneGeometry,
        coords: Coordinates,
        first_id: u64,
    ) -> Self {
        Self {
            out,
            grid,
            geometry,
            coords,
            next_id: first_id,
            written: Written::default(),
        }
    }

    /// Writes one zone as a feature.
    ///
    /// The null zone is refused with [`Error::NullZone`] before anything is written: it names no
    /// zone, and a list that may hold it is filtered by the caller. A refinement beyond the
    /// library's limit is refused with [`Error::Grid`], and a geometry with a coordinate that is
    /// not finite with [`Error::NotFinite`], also before anything is written.
    pub fn zone(&mut self, zone: ZoneId) -> Result<()> {
        self.feature(zone, |_| Ok(()))
    }

    /// Writes one zone as a feature, as [`ZoneListGeoJson::zone`] does, its properties beginning
    /// with what `properties` writes: members, each followed by a comma.
    pub(crate) fn feature(
        &mut self,
        zone: ZoneId,
        properties: impl FnOnce(&mut W) -> io::Result<()>,
    ) -> Result<()> {
        let d = drawn(self.grid, zone, self.geometry, self.coords)?;
        self.out.write_all(if self.written.features == 0 {
            br#"{"type":"FeatureCollection","features":["#
        } else {
            b","
        })?;
        write_feature(
            &mut self.out,
            self.grid,
            zone,
            self.next_id,
            d.as_ref(),
            self.coords,
            properties,
        )?;
        // No list reaches the end of the numbers; a caller who begins near it sees them wrap.
        self.next_id = self.next_id.wrapping_add(1);
        self.written.features += 1;
        if d.is_none() {
            self.written.without_geometry += 1;
        }
        Ok(())
    }

    /// Ends the document and returns what it holds. The sink is flushed.
    pub fn finish(mut self) -> Result<Written> {
        self.out.write_all(if self.written.features == 0 {
            br#"{"type":"FeatureCollection","features":[]}"#
        } else {
            b"]}"
        })?;
        self.out.flush()?;
        Ok(self.written)
    }
}

impl<W: Write> fmt::Debug for ZoneListGeoJson<W> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("ZoneListGeoJson")
            .field("grid", &self.grid.name())
            .field("geometry", &self.geometry)
            .field("coords", &self.coords)
            .field("next_id", &self.next_id)
            .field("written", &self.written)
            .finish_non_exhaustive()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{census, prepared};

    const REGION: ZoneGeometry = ZoneGeometry::Region { edge_refinement: 0 };
    const SHORTEST: Coordinates = Coordinates { decimals: None };

    fn grid(name: &str) -> AnyGrid {
        rs4dggs::get_grid(name).unwrap()
    }

    fn zone(grid: AnyGrid, text: &str) -> ZoneId {
        grid.zone_from_text(text).unwrap()
    }

    /// What `write` wrote, as text, and what it answered.
    fn written<T>(write: impl FnOnce(&mut Vec<u8>) -> Result<T>) -> (String, T) {
        let mut out = Vec::new();
        let answer = write(&mut out).unwrap();
        (String::from_utf8(out).unwrap(), answer)
    }

    /// The list of `zones`, as text, and what it holds.
    fn list(
        grid: AnyGrid,
        zones: &[ZoneId],
        geometry: ZoneGeometry,
        coords: Coordinates,
        first_id: u64,
    ) -> (String, Written) {
        written(|out| {
            let mut list = ZoneListGeoJson::new(out, grid, geometry, coords, first_id);
            for &z in zones {
                list.zone(z)?;
            }
            list.finish()
        })
    }

    /// The zones of the census sample whose refined ring has no area to draw.
    fn undrawable(grid: AnyGrid) -> Vec<ZoneId> {
        census::sample(&grid)
            .into_iter()
            .filter(|&z| prepared(&grid.refined_vertices(z, 0).unwrap()).is_none())
            .collect()
    }

    #[test]
    fn the_list_over_the_census_counts_every_zone_and_those_that_cannot_be_drawn() {
        for (name, cannot) in [
            ("ISEA3H", 0),
            ("IVEA3H", 0),
            ("RTEA3H", 0),
            ("IGEO7", 3),
            ("IVEA7H", 3),
            ("RTEA7H", 3),
        ] {
            let g = grid(name);
            let sample = census::sample(&g);
            let mut list = ZoneListGeoJson::new(io::sink(), g, REGION, SHORTEST, 1);
            for &z in &sample {
                list.zone(z).unwrap();
            }
            let w = list.finish().unwrap();
            assert_eq!(w.features, sample.len() as u64, "{name}");
            assert_eq!(w.without_geometry, cannot, "{name}");
            assert_eq!(w.without_geometry, undrawable(g).len() as u64, "{name}");
            assert_eq!(w.null_zones_skipped, 0, "{name}");
        }
    }

    #[test]
    fn the_single_point_rings_of_the_finest_level_have_no_geometry() {
        for name in ["IGEO7", "IVEA7H", "RTEA7H"] {
            let g = grid(name);
            let zones = undrawable(g);
            assert_eq!(zones.len(), 3, "{name}");
            for &z in &zones {
                let text = format!("{name} {}", g.text_id(z));
                assert_eq!(g.resolution(z), 19, "{text}");
                let ring = g.refined_vertices(z, 0).unwrap();
                assert!(!ring.is_empty(), "{text}");
                assert!(ring.iter().all(|p| *p == ring[0]), "{text}: one point");
                // What the cleaning returns for such a ring is no valid geometry, and is never
                // written.
                assert!(census::valid(&geometry::clean(&ring)).is_err(), "{text}");
                let (out, drawn) = written(|o| write_zone_geometry(o, g, z, REGION, SHORTEST));
                assert!(!drawn, "{text}");
                assert!(out.is_empty(), "{text}: {out}");
                let (out, drawn) = written(|o| write_zone_feature(o, g, z, 7, REGION, SHORTEST));
                assert!(!drawn, "{text}");
                assert_eq!(
                    out,
                    format!(
                        r#"{{"type":"Feature","id":7,"geometry":null,"properties":{{"zoneID":"{}"}}}}"#,
                        g.text_id(z)
                    )
                );
            }
            let (out, w) = list(g, &zones, REGION, SHORTEST, 1);
            assert_eq!((w.features, w.without_geometry), (3, 3), "{name}");
            assert_eq!(out.matches(r#""geometry":null"#).count(), 3, "{name}");
        }
    }

    #[test]
    fn no_geometry_is_a_null_member_of_a_feature_and_nothing_alone() {
        let g = grid("ISEA3H");
        let z = zone(g, "E6-317-A");
        let none = ZoneGeometry::None;
        let (out, drawn) = written(|o| write_zone_feature(o, g, z, 1, none, SHORTEST));
        assert!(!drawn);
        assert_eq!(
            out,
            r#"{"type":"Feature","id":1,"geometry":null,"properties":{"zoneID":"E6-317-A"}}"#
        );
        let (out, drawn) = written(|o| write_zone_geometry(o, g, z, none, SHORTEST));
        assert!(!drawn);
        assert!(out.is_empty());
        let (out, w) = list(g, &[z, z], none, SHORTEST, 1);
        assert_eq!((w.features, w.without_geometry), (2, 2));
        assert_eq!(out.matches(r#""geometry":null"#).count(), 2);
    }

    #[test]
    fn the_geometry_alone_is_the_geometry_member_of_the_feature() {
        let g = grid("ISEA3H");
        let z = zone(g, "E6-317-A");
        for geometry in [REGION, ZoneGeometry::Centroid] {
            let (geom, drawn) = written(|o| write_zone_geometry(o, g, z, geometry, SHORTEST));
            assert!(drawn);
            let (feature, drawn) = written(|o| write_zone_feature(o, g, z, 1, geometry, SHORTEST));
            assert!(drawn);
            assert_eq!(
                feature,
                format!(
                    r#"{{"type":"Feature","id":1,"geometry":{geom},"properties":{{"zoneID":"E6-317-A"}}}}"#
                )
            );
            let (listed, _) = list(g, &[z], geometry, SHORTEST, 1);
            assert_eq!(
                listed,
                format!(r#"{{"type":"FeatureCollection","features":[{feature}]}}"#)
            );
        }
    }

    #[test]
    fn the_features_are_numbered_from_the_first_id() {
        let g = grid("ISEA3H");
        let zones = [zone(g, "E6-317-A"), zone(g, "A8-0-C")];
        let (out, w) = list(g, &zones, ZoneGeometry::None, SHORTEST, 41);
        assert_eq!(w.features, 2);
        assert!(
            out.starts_with(
                r#"{"type":"FeatureCollection","features":[{"type":"Feature","id":41,"#
            ),
            "{out}"
        );
        assert!(out.contains(r#"},{"type":"Feature","id":42,"#), "{out}");
    }

    /// The numbers of the coordinates of every geometry in `doc`.
    fn coordinates(doc: &str) -> Vec<&str> {
        doc.split(r#""coordinates":"#)
            .skip(1)
            .flat_map(|rest| {
                let end = rest.find('}').unwrap();
                rest[..end].split(['[', ']', ',']).filter(|t| !t.is_empty())
            })
            .collect()
    }

    /// Whether every number of the coordinates in `doc` has `n` decimals.
    fn every_coordinate_has(doc: &str, n: usize) -> bool {
        let numbers = coordinates(doc);
        !numbers.is_empty()
            && numbers.iter().all(|t| {
                let Some((whole, decimals)) = t.split_once('.') else {
                    return false;
                };
                let whole = whole.strip_prefix('-').unwrap_or(whole);
                !whole.is_empty()
                    && whole.bytes().all(|b| b.is_ascii_digit())
                    && decimals.len() == n
                    && decimals.bytes().all(|b| b.is_ascii_digit())
            })
    }

    #[test]
    fn by_default_every_coordinate_has_eight_decimals_and_never_fewer() {
        assert_eq!(Coordinates::default(), Coordinates { decimals: Some(8) });
        for (name, texts) in [
            ("ISEA3H", ["E6-317-A", "A8-0-C"]),
            ("IGEO7", ["01012", "005"]),
        ] {
            let g = grid(name);
            let zones: Vec<ZoneId> = texts.iter().map(|t| zone(g, t)).collect();
            for geometry in [REGION, ZoneGeometry::Centroid] {
                let (eight, w) = list(g, &zones, geometry, Coordinates::default(), 1);
                assert_eq!(w.without_geometry, 0, "{name}");
                assert!(every_coordinate_has(&eight, 8), "{name}: {eight}");
                // Fewer decimals asked are eight.
                for fewer in [0, 3, 6, 7] {
                    let (out, _) = list(
                        g,
                        &zones,
                        geometry,
                        Coordinates {
                            decimals: Some(fewer),
                        },
                        1,
                    );
                    assert_eq!(out, eight, "{name}, {fewer} decimals asked");
                }
                let (ten, _) = list(g, &zones, geometry, Coordinates { decimals: Some(10) }, 1);
                assert!(every_coordinate_has(&ten, 10), "{name}: {ten}");
            }
        }
        // A cut across the antimeridian is still a MultiPolygon once rounded.
        let g = grid("IGEO7");
        let (out, drawn) = written(|o| {
            write_zone_geometry(o, g, zone(g, "01012"), REGION, Coordinates::default())
        });
        assert!(drawn);
        assert!(out.starts_with(r#"{"type":"MultiPolygon","#), "{out}");
    }

    #[test]
    fn a_zone_of_the_finest_level_asked_with_no_decimals_keeps_its_geometry() {
        // Rounded to whole degrees, its ring would be one point; eight decimals keep it.
        let g = grid("IGEO7");
        let lisbon = g.zone_from_geo(38.7223, -9.1393, 19).unwrap();
        let ring = g.refined_vertices(lisbon, 0).unwrap();
        assert!(geometry::drawn(&ring, Some(0)).is_none());
        let (out, w) = list(g, &[lisbon], REGION, Coordinates { decimals: Some(0) }, 1);
        assert_eq!((w.features, w.without_geometry), (1, 0));
        assert!(every_coordinate_has(&out, 8), "{out}");
        census::valid(&geometry::drawn(&ring, Some(8)).unwrap()).unwrap();
    }

    #[test]
    fn a_refinement_beyond_the_limit_is_refused_before_anything_is_written() {
        let g = grid("ISEA3H");
        let z = zone(g, "E6-317-A");
        let beyond = ZoneGeometry::Region {
            edge_refinement: 100_001,
        };
        let mut out = Vec::new();
        let e = write_zone_geometry(&mut out, g, z, beyond, SHORTEST).unwrap_err();
        assert!(matches!(e, Error::Grid(_)), "{e}");
        let e = write_zone_feature(&mut out, g, z, 1, beyond, SHORTEST).unwrap_err();
        assert!(matches!(e, Error::Grid(_)), "{e}");
        let mut list = ZoneListGeoJson::new(&mut out, g, beyond, SHORTEST, 1);
        assert!(matches!(list.zone(z), Err(Error::Grid(_))));
        assert_eq!(list.finish().unwrap(), Written::default());
        assert_eq!(out, br#"{"type":"FeatureCollection","features":[]}"#);
    }

    #[test]
    fn the_null_zone_is_refused_before_anything_is_written() {
        let g = grid("IGEO7");
        let null = ZoneId::NULL;
        let mut out = Vec::new();
        for geometry in [REGION, ZoneGeometry::Centroid, ZoneGeometry::None] {
            let e = write_zone_geometry(&mut out, g, null, geometry, SHORTEST);
            assert!(matches!(e, Err(Error::NullZone)));
            let e = write_zone_feature(&mut out, g, null, 1, geometry, SHORTEST);
            assert!(matches!(e, Err(Error::NullZone)));
            let mut list = ZoneListGeoJson::new(&mut out, g, geometry, SHORTEST, 1);
            assert!(matches!(list.zone(null), Err(Error::NullZone)));
        }
        assert!(out.is_empty());
    }

    #[test]
    fn a_coordinate_made_not_finite_by_the_cleaning_is_refused() {
        // Finite, but so great that the cleaning makes a longitude that is not a number: GeoJSON
        // has no form for it, and it is neither written nor taken for a zone without geometry.
        let extreme: Vec<GeoPoint> = [[f64::MAX, 0.0], [-f64::MAX, 0.0], [0.0, 10.0]]
            .iter()
            .map(|p| GeoPoint {
                lon: p[0],
                lat: p[1],
            })
            .collect();
        for coords in [SHORTEST, Coordinates { decimals: Some(3) }] {
            let e = region(extreme.clone(), coords);
            assert!(matches!(e, Err(Error::NotFinite)), "{coords:?}");
        }
        // A ring with a coordinate that is not finite as it comes is a zone that cannot be drawn.
        let nan = vec![
            GeoPoint {
                lon: f64::NAN,
                lat: 0.0
            };
            3
        ];
        assert!(matches!(region(nan, SHORTEST), Ok(None)));
    }

    #[test]
    fn an_empty_list_is_an_empty_collection() {
        let (out, w) = list(grid("ISEA3H"), &[], REGION, SHORTEST, 1);
        assert_eq!(out, r#"{"type":"FeatureCollection","features":[]}"#);
        assert_eq!(w, Written::default());
    }

    /// A sink that takes `left` bytes and then refuses.
    struct Failing {
        left: usize,
    }

    impl Write for Failing {
        fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
            if self.left == 0 {
                return Err(io::Error::other("full"));
            }
            let n = buf.len().min(self.left);
            self.left -= n;
            Ok(n)
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn a_sink_that_fails_in_the_middle_is_an_error() {
        let g = grid("IGEO7");
        let zones = [zone(g, "01012"), zone(g, "005")];
        for left in [0, 1, 40, 700, 2000] {
            let e = write_zone_geometry(Failing { left }, g, zones[0], REGION, SHORTEST);
            assert!(matches!(e, Err(Error::Io(_))), "{left}");
            let e = write_zone_feature(Failing { left }, g, zones[0], 1, REGION, SHORTEST);
            assert!(matches!(e, Err(Error::Io(_))), "{left}");
            let mut list = ZoneListGeoJson::new(Failing { left }, g, REGION, SHORTEST, 1);
            let e = zones.iter().try_for_each(|&z| list.zone(z));
            assert!(matches!(e, Err(Error::Io(_))), "{left}");
        }
    }
}
