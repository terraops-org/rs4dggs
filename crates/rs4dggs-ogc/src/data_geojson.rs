//! The data of a zone as GeoJSON: each sub-zone a feature, with its values among its properties.

use std::io::{self, Write};

use rs4dggs::ZoneId;

use crate::data::{RESERVED, Values, ZoneData, check};
use crate::{Coordinates, Error, Result, Written, ZoneGeometry, ZoneListGeoJson, json};

/// Writes the data of a zone as GeoJSON, the document a request for zone data answers with
/// `application/geo+json`: one `FeatureCollection` holding, for each relative depth in the order
/// of [`ZoneData::depths`], one feature for each sub-zone at that depth, in the order of the
/// library's `sub_zones`. Each feature is written as [`ZoneListGeoJson`] writes a zone, its
/// `properties` holding one member for each property, in the order given, and then `zoneID`:
/// `{"type":"Feature","id":1,"geometry":null,"properties":{"Elevation":-14.2,"zoneID":"E0-8C1-A"}}`
/// for a band `Elevation`, without geometry. The features are numbered from 1 through all the
/// depths, and a sub-zone's level is read from its identifier.
///
/// Values are written as in DGGS-JSON: integers as integers, a value of `f32` in the shortest
/// form that reads back to the same single, one of `f64` in the shortest form that reads back to
/// the same double, and a missing value, or one that is not finite, as `null`. The coordinates
/// of a geometry are written as `coords` says: with eight decimals by default and never fewer, or
/// in the shortest form that reads back to the same double (see [`Coordinates`]). A zone that
/// cannot be drawn is written with `"geometry":null`, as is every zone under
/// [`ZoneGeometry::None`], and
/// counted in [`Written::without_geometry`]. An entry of the sub-zone order that is the null zone
/// (on an aperture-7 grid, at a broken seam) names no zone: it is passed over, its value with it,
/// without taking a number, and counted in [`Written::null_zones_skipped`]. The sink is flushed.
///
/// The data are held first to what the document requires and refused before anything is written,
/// as [`write_dggs_json`](crate::write_dggs_json) refuses them (a property named `zoneID`, which
/// would stand beside the member of that name, among the faults). The identifiers of the
/// sub-zones at every depth are then obtained from the library, again before anything is written:
/// a list beyond the library's limit on the sub-zones it materialises is refused with
/// [`Error::Grid`](crate::Error::Grid), as is a refinement beyond its limit, which the first zone
/// to be drawn meets. Each value is written to the sink as it comes, in many small writes: a
/// service passes a buffered sink. A failure of the sink, or a geometry with a coordinate that is
/// not finite ([`Error::NotFinite`](crate::Error::NotFinite), which no zone of the library's grids
/// gives), leaves the document truncated where it arose.
///
/// # Examples
///
/// ```
/// use rs4dggs_ogc::{
///     Coordinates, Dggrs, Property, Values, ZoneData, ZoneGeometry, write_zone_data_geojson,
/// };
///
/// let dggrs = Dggrs::from_id("ISEA3H").unwrap();
/// let zone = dggrs.grid().zone_from_text("C2-23-C")?;
/// let band = [Some(1), Some(2), None, Some(4), Some(5), Some(6), Some(7)];
/// let data = ZoneData {
///     dggrs,
///     zone,
///     depths: &[1],
///     properties: &[Property { name: "class", depths: &[Values::U8(&band)] }],
/// };
/// let mut out = Vec::new();
/// let coords = Coordinates::default();
/// let written = write_zone_data_geojson(&mut out, &data, ZoneGeometry::None, coords)?;
/// assert_eq!((written.features, written.without_geometry), (7, 7));
/// let out = String::from_utf8(out)?;
/// assert!(out.starts_with(r#"{"type":"FeatureCollection","features":[{"type":"Feature","#));
/// assert!(out.contains(r#""id":3,"geometry":null,"properties":{"class":null,"zoneID":"#));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn write_zone_data_geojson<W: Write>(
    out: W,
    data: &ZoneData,
    geometry: ZoneGeometry,
    coords: Coordinates,
) -> Result<Written> {
    check(data, RESERVED)?;
    let grid = data.dggrs.grid();
    let orders = data
        .depths
        .iter()
        .map(|&depth| grid.sub_zones(data.zone, depth))
        .collect::<rs4dggs::Result<Vec<_>>>()?;
    // The library makes the length of each order the count that `check` held the values to, so
    // that every position of an order has its value; were the two to differ, the document is
    // refused here, before its first byte, and not left unfinished at the first missing value.
    let mut at_depth = Vec::with_capacity(orders.len());
    for (i, (order, &depth)) in orders.iter().zip(data.depths).enumerate() {
        let mut values = Vec::with_capacity(data.properties.len());
        for property in data.properties {
            let v = property.depths.get(i);
            let found = v.map_or(0, Values::len);
            match v {
                Some(&v) if found == order.len() => values.push(v),
                _ => {
                    return Err(Error::Shape {
                        property: property.name.into(),
                        depth,
                        expected: order.len() as u64,
                        found: found as u64,
                    });
                }
            }
        }
        at_depth.push(values);
    }
    let mut list = ZoneListGeoJson::new(out, grid, geometry, coords, 1);
    let mut skipped = 0;
    for (order, values) in orders.iter().zip(&at_depth) {
        for (k, &zone) in order.iter().enumerate() {
            if zone == ZoneId::NULL {
                skipped += 1;
                continue;
            }
            list.feature(zone, |out| {
                for (property, &v) in data.properties.iter().zip(values) {
                    json::string(out, property.name)?;
                    out.write_all(b":")?;
                    value(out, v, k)?;
                    out.write_all(b",")?;
                }
                Ok(())
            })?;
        }
    }
    let mut written = list.finish()?;
    written.null_zones_skipped = skipped;
    Ok(written)
}

/// The value at position `k` of `values`, or `null` where it is missing or not finite. (The
/// values were held to the length of the order before the first byte, so that `k` is always
/// within them.)
fn value<W: Write>(out: &mut W, values: Values, k: usize) -> io::Result<()> {
    fn one<W: Write, T: Copy>(
        out: &mut W,
        value: Option<T>,
        write: impl Fn(&mut W, T) -> io::Result<()>,
    ) -> io::Result<()> {
        match value {
            Some(x) => write(out, x),
            None => out.write_all(b"null"),
        }
    }
    match values {
        Values::U8(s) => one(out, s.get(k).copied().flatten(), |o, x| write!(o, "{x}")),
        Values::I8(s) => one(out, s.get(k).copied().flatten(), |o, x| write!(o, "{x}")),
        Values::U16(s) => one(out, s.get(k).copied().flatten(), |o, x| write!(o, "{x}")),
        Values::I16(s) => one(out, s.get(k).copied().flatten(), |o, x| write!(o, "{x}")),
        Values::U32(s) => one(out, s.get(k).copied().flatten(), |o, x| write!(o, "{x}")),
        Values::I32(s) => one(out, s.get(k).copied().flatten(), |o, x| write!(o, "{x}")),
        Values::F32(s) => one(out, s.get(k).copied().flatten(), json::f32),
        Values::F64(s) => one(out, s.get(k).copied().flatten(), json::f64),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Dggrs, Error, Property};

    fn isea3h() -> &'static Dggrs {
        Dggrs::from_id("ISEA3H").unwrap()
    }

    /// The data of `C2-23-C` at `depths`, with the properties given.
    fn data<'a>(depths: &'a [u8], properties: &'a [Property<'a>]) -> ZoneData<'a> {
        ZoneData {
            dggrs: isea3h(),
            zone: isea3h().grid().zone_from_text("C2-23-C").unwrap(),
            depths,
            properties,
        }
    }

    /// The document of `data` with `geometry`, and what it holds.
    fn written(data: &ZoneData, geometry: ZoneGeometry) -> (String, Written) {
        let mut out = Vec::new();
        let w = write_zone_data_geojson(&mut out, data, geometry, Coordinates::default()).unwrap();
        (String::from_utf8(out).unwrap(), w)
    }

    /// The error of writing `data`, having asserted that nothing reached the sink.
    fn refused(data: &ZoneData, geometry: ZoneGeometry) -> Error {
        let mut out = Vec::new();
        let e =
            write_zone_data_geojson(&mut out, data, geometry, Coordinates::default()).unwrap_err();
        assert!(out.is_empty(), "{} bytes written before {e}", out.len());
        e
    }

    /// The members `"id":N` of the features, in order.
    fn ids(doc: &str) -> Vec<u64> {
        doc.split(r#"{"type":"Feature","id":"#)
            .skip(1)
            .map(|f| f[..f.find(',').unwrap()].parse().unwrap())
            .collect()
    }

    /// The properties of each feature, as written.
    fn properties(doc: &str) -> Vec<&str> {
        doc.split(r#""properties":{"#)
            .skip(1)
            .map(|f| &f[..f.find('}').unwrap()])
            .collect()
    }

    const SEVEN: [Option<f64>; 7] = [
        Some(-14.209265189621),
        Some(57.2747655762701),
        None,
        Some(f64::NAN),
        Some(f64::INFINITY),
        Some(-0.0),
        Some(1e17),
    ];

    #[test]
    fn each_sub_zone_is_a_feature_with_its_value_and_then_its_identifier() {
        let band = [Property {
            name: "Elevation",
            depths: &[Values::F64(&SEVEN)],
        }];
        let (doc, w) = written(&data(&[1], &band), ZoneGeometry::None);
        let grid = isea3h().grid();
        let zone = grid.zone_from_text("C2-23-C").unwrap();
        let texts: Vec<String> = grid
            .sub_zones(zone, 1)
            .unwrap()
            .into_iter()
            .map(|z| grid.text_id(z))
            .collect();
        let values = [
            "-14.209265189621",
            "57.2747655762701",
            "null",
            "null",
            "null",
            "-0",
            "1e17",
        ];
        let expected: Vec<String> = values
            .iter()
            .zip(&texts)
            .map(|(v, t)| format!(r#""Elevation":{v},"zoneID":"{t}""#))
            .collect();
        assert_eq!(properties(&doc), expected);
        assert_eq!(ids(&doc), (1..=7).collect::<Vec<_>>());
        let first = r#"{"type":"Feature","id":1,"geometry":null,"properties":{"#;
        assert!(doc.starts_with(&format!(
            r#"{{"type":"FeatureCollection","features":[{first}"#
        )));
        assert!(doc.ends_with("}}]}"), "{doc}");
        assert_eq!(
            w,
            Written {
                features: 7,
                without_geometry: 7,
                null_zones_skipped: 0
            }
        );
    }

    #[test]
    fn the_bands_stand_in_the_order_given_each_in_its_own_form() {
        let class = [
            Some(0u8),
            Some(1),
            Some(2),
            Some(3),
            Some(4),
            Some(5),
            Some(255),
        ];
        let t = [
            Some(0.1f32),
            None,
            Some(f32::NAN),
            Some(1.5),
            Some(2.0),
            Some(3.0),
            Some(4.0),
        ];
        let i = [
            Some(-1i16),
            Some(i16::MIN),
            None,
            Some(0),
            Some(1),
            Some(2),
            Some(i16::MAX),
        ];
        let bands = [
            Property {
                name: "t",
                depths: &[Values::F32(&t)],
            },
            Property {
                name: "class",
                depths: &[Values::U8(&class)],
            },
            Property {
                name: "a \"quoted\" name",
                depths: &[Values::I16(&i)],
            },
        ];
        let (doc, _) = written(&data(&[1], &bands), ZoneGeometry::None);
        let p = properties(&doc);
        assert!(
            p[0].starts_with(r#""t":0.1,"class":0,"a \"quoted\" name":-1,"zoneID":"#),
            "{}",
            p[0]
        );
        assert!(
            p[1].starts_with(r#""t":null,"class":1,"a \"quoted\" name":-32768,"zoneID":"#),
            "{}",
            p[1]
        );
        assert!(
            p[2].starts_with(r#""t":null,"class":2,"a \"quoted\" name":null,"zoneID":"#),
            "{}",
            p[2]
        );
        assert!(
            p[6].starts_with(r#""t":4,"class":255,"a \"quoted\" name":32767,"zoneID":"#),
            "{}",
            p[6]
        );
    }

    #[test]
    fn every_depth_follows_the_one_before_and_the_numbers_run_through_them_all() {
        let seven = [Some(7u32); 7];
        let thirteen = [Some(13u32); 13];
        let band = [Property {
            name: "n",
            depths: &[Values::U32(&seven), Values::U32(&thirteen)],
        }];
        let (doc, w) = written(&data(&[1, 2], &band), ZoneGeometry::None);
        assert_eq!(ids(&doc), (1..=20).collect::<Vec<_>>());
        let p = properties(&doc);
        assert!(p[..7].iter().all(|p| p.starts_with(r#""n":7,"#)));
        assert!(p[7..].iter().all(|p| p.starts_with(r#""n":13,"#)));
        let grid = isea3h().grid();
        let zone = grid.zone_from_text("C2-23-C").unwrap();
        let order: Vec<String> = [1, 2]
            .iter()
            .flat_map(|&d| grid.sub_zones(zone, d).unwrap())
            .map(|z| format!(r#""zoneID":"{}""#, grid.text_id(z)))
            .collect();
        assert!(p.iter().zip(&order).all(|(p, z)| p.ends_with(z.as_str())));
        assert_eq!((w.features, w.null_zones_skipped), (20, 0));
    }

    #[test]
    fn a_null_entry_of_the_order_is_passed_over_without_a_number() {
        let z7 = Dggrs::from_id("ISEA7H_Z7").unwrap();
        let seam = z7.grid().zone_from_text("00055353260226021").unwrap();
        let positions: Vec<Option<i32>> = (0..13).map(Some).collect();
        let band = [Property {
            name: "i",
            depths: &[Values::I32(&positions)],
        }];
        let data = ZoneData {
            dggrs: z7,
            zone: seam,
            depths: &[1],
            properties: &band,
        };
        let (doc, w) = written(&data, ZoneGeometry::None);
        assert_eq!(
            w,
            Written {
                features: 11,
                without_geometry: 11,
                null_zones_skipped: 2
            }
        );
        assert_eq!(ids(&doc), (1..=11).collect::<Vec<_>>());
        // The values of the entries 8 and 12, the null zone, are not written.
        let values: Vec<&str> = properties(&doc)
            .iter()
            .map(|p| &p[r#""i":"#.len()..p.find(',').unwrap()])
            .collect();
        assert_eq!(
            values,
            ["0", "1", "2", "3", "4", "5", "6", "7", "9", "10", "11"]
        );
    }

    #[test]
    fn a_zone_drawn_carries_its_geometry_and_counts_none_without() {
        let positions: Vec<Option<u8>> = (0..7).map(Some).collect();
        let band = [Property {
            name: "i",
            depths: &[Values::U8(&positions)],
        }];
        let (doc, w) = written(&data(&[1], &band), ZoneGeometry::Centroid);
        assert_eq!(w.without_geometry, 0);
        assert_eq!(
            doc.matches(r#""geometry":{"type":"Point","coordinates":["#)
                .count(),
            7
        );
        let (doc, w) = written(
            &data(&[1], &band),
            ZoneGeometry::Region { edge_refinement: 0 },
        );
        assert_eq!(w.without_geometry, 0);
        assert_eq!(
            doc.matches(r#""geometry":{"type":"Polygon","coordinates":[[["#)
                .count(),
            7
        );
    }

    #[test]
    fn a_band_named_zone_id_is_refused_before_anything_is_written() {
        let positions: Vec<Option<u8>> = (0..7).map(Some).collect();
        let band = [Property {
            name: "zoneID",
            depths: &[Values::U8(&positions)],
        }];
        assert!(matches!(
            refused(&data(&[1], &band), ZoneGeometry::None),
            Error::ReservedName(n) if n == "zoneID"
        ));
    }

    #[test]
    fn values_one_short_are_refused_before_anything_is_written() {
        let six: Vec<Option<u8>> = (0..6).map(Some).collect();
        let band = [Property {
            name: "i",
            depths: &[Values::U8(&six)],
        }];
        assert!(matches!(
            refused(&data(&[1], &band), ZoneGeometry::None),
            Error::Shape {
                expected: 7,
                found: 6,
                ..
            }
        ));
    }

    #[test]
    fn a_refinement_beyond_the_limit_is_refused_before_anything_is_written() {
        let positions: Vec<Option<u8>> = (0..7).map(Some).collect();
        let band = [Property {
            name: "i",
            depths: &[Values::U8(&positions)],
        }];
        let region = ZoneGeometry::Region {
            edge_refinement: 100_001,
        };
        assert!(matches!(
            refused(&data(&[1], &band), region),
            Error::Grid(_)
        ));
    }

    #[test]
    fn no_depth_is_an_empty_collection() {
        let (doc, w) = written(&data(&[], &[]), ZoneGeometry::None);
        assert_eq!(doc, r#"{"type":"FeatureCollection","features":[]}"#);
        assert_eq!(w, Written::default());
    }

    #[test]
    fn a_sink_that_fails_in_the_middle_is_an_error() {
        struct Failing(usize);
        impl Write for Failing {
            fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
                if self.0 < buf.len() {
                    return Err(io::Error::other("full"));
                }
                self.0 -= buf.len();
                Ok(buf.len())
            }
            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }
        let positions: Vec<Option<f64>> = (0..7).map(|i| Some(f64::from(i))).collect();
        let band = [Property {
            name: "i",
            depths: &[Values::F64(&positions)],
        }];
        for at in [0, 1, 40, 300, 600] {
            let e = write_zone_data_geojson(
                Failing(at),
                &data(&[1], &band),
                ZoneGeometry::Centroid,
                Coordinates::default(),
            )
            .unwrap_err();
            assert!(matches!(e, Error::Io(_)), "{at}: {e}");
        }
    }
}
