//! The data of a zone, the values of its sub-zones, as DGGS-JSON.

use std::io::{self, Write};

use rs4dggs::ZoneId;

use crate::{Dggrs, Error, Result, json};

/// The values of one property at one relative depth, one for each sub-zone, in the order of the
/// zone's sub-zones at that depth (the library's `sub_zones`). `None` is a missing value.
///
/// A value of `f32` or `f64` that is not finite is written as a missing value, since neither
/// encoding of the standard has a form for it.
///
/// The type may gain variants in a later release, so that a `match` on it outside this crate
/// needs an arm for the others.
#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub enum Values<'a> {
    /// Unsigned integers of 8 bits.
    U8(&'a [Option<u8>]),
    /// Signed integers of 8 bits.
    I8(&'a [Option<i8>]),
    /// Unsigned integers of 16 bits.
    U16(&'a [Option<u16>]),
    /// Signed integers of 16 bits.
    I16(&'a [Option<i16>]),
    /// Unsigned integers of 32 bits.
    U32(&'a [Option<u32>]),
    /// Signed integers of 32 bits.
    I32(&'a [Option<i32>]),
    /// Floating-point numbers of single precision.
    F32(&'a [Option<f32>]),
    /// Floating-point numbers of double precision.
    F64(&'a [Option<f64>]),
}

impl Values<'_> {
    /// The number of values, missing ones included.
    pub fn len(&self) -> usize {
        match self {
            Values::U8(s) => s.len(),
            Values::I8(s) => s.len(),
            Values::U16(s) => s.len(),
            Values::I16(s) => s.len(),
            Values::U32(s) => s.len(),
            Values::I32(s) => s.len(),
            Values::F32(s) => s.len(),
            Values::F64(s) => s.len(),
        }
    }

    /// Whether there are no values at all.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// One property of the data (a band of a raster, for instance), with its values at each
/// relative depth asked.
///
/// Built by [`Property::new`]; the type may gain fields in a later release, so that it is not
/// built from its fields outside this crate.
///
/// Built from its fields outside this crate, it does not compile:
///
/// ```compile_fail,E0639
/// let property = rs4dggs_ogc::Property { name: "t", depths: &[] };
/// ```
#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub struct Property<'a> {
    /// The name of the property, as the document carries it.
    pub name: &'a str,
    /// The values at each relative depth, in the order of [`ZoneData::depths`]: one entry for
    /// each depth.
    pub depths: &'a [Values<'a>],
}

/// The data of one zone: the values of its sub-zones at one or more relative depths, for one or
/// more properties.
///
/// Built by [`ZoneData::new`]; the type may gain fields in a later release, so that it is not
/// built from its fields outside this crate. The arguments borrow their slices: a slice built in
/// the call itself lives only to the end of its statement, so that data used in later statements
/// take their properties and values from slices bound beforehand, as the example of
/// [`write_dggs_json`] does.
///
/// Built from its fields outside this crate, it does not compile:
///
/// ```compile_fail,E0639
/// let dggrs = rs4dggs_ogc::Dggrs::from_id("ISEA3H").unwrap();
/// let zone = dggrs.grid().zone_from_text("C2-23-C").unwrap();
/// let data = rs4dggs_ogc::ZoneData { dggrs, zone, depths: &[1], properties: &[] };
/// ```
#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub struct ZoneData<'a> {
    /// The grid, whose URI the document carries.
    pub dggrs: &'static Dggrs,
    /// The zone whose sub-zones the values describe.
    pub zone: ZoneId,
    /// The relative depths, in the order the document lists them. At depth 0 the one sub-zone
    /// is the zone itself.
    pub depths: &'a [u8],
    /// The properties, in the order the document lists them.
    pub properties: &'a [Property<'a>],
}

impl<'a> Property<'a> {
    /// The property `name`, with its values at each relative depth: `depths` holds one entry
    /// for each depth of the data, in their order.
    pub const fn new(name: &'a str, depths: &'a [Values<'a>]) -> Self {
        Self { name, depths }
    }
}

impl<'a> ZoneData<'a> {
    /// The data of `zone` of `dggrs`: at each relative depth of `depths`, in that order, the
    /// values of every property of `properties`.
    pub const fn new(
        dggrs: &'static Dggrs,
        zone: ZoneId,
        depths: &'a [u8],
        properties: &'a [Property<'a>],
    ) -> Self {
        Self {
            dggrs,
            zone,
            depths,
            properties,
        }
    }
}

/// The names the data documents reserve: a property of this name would be confused, in zone
/// data as GeoJSON, with the member that holds each zone's identifier, and is refused in every
/// encoding so that one set of data is accepted by all of them or by none.
pub(crate) const RESERVED: &[&str] = &["zoneID"];

/// Holds `data` to what a document of it requires, and answers the number of sub-zones at each
/// of its depths, in their order. Nothing is written: every fault found here refuses the
/// document before its first byte.
///
/// Refused: the null zone ([`Error::NullZone`]); a zone or a depth the library refuses
/// ([`Error::Grid`]); a property whose values are not given once for each depth
/// ([`Error::Depths`]); values whose number is not that of the sub-zones at their depth
/// ([`Error::Shape`]); a name among `reserved` ([`Error::ReservedName`]); a name given twice
/// ([`Error::RepeatedName`]).
pub(crate) fn check(data: &ZoneData, reserved: &[&str]) -> Result<Vec<u64>> {
    if data.zone == ZoneId::NULL {
        return Err(Error::NullZone);
    }
    let grid = data.dggrs.grid();
    // The zone is validated even when no depth is asked, which the counts alone would not do.
    grid.count_sub_zones(data.zone, 0)?;
    let counts = data
        .depths
        .iter()
        .map(|&depth| grid.count_sub_zones(data.zone, depth))
        .collect::<rs4dggs::Result<Vec<u64>>>()?;
    for (i, property) in data.properties.iter().enumerate() {
        let name = property.name;
        if reserved.contains(&name) {
            return Err(Error::ReservedName(name.into()));
        }
        if data.properties[..i].iter().any(|p| p.name == name) {
            return Err(Error::RepeatedName(name.into()));
        }
        if property.depths.len() != counts.len() {
            return Err(Error::Depths {
                property: name.into(),
                expected: counts.len(),
                found: property.depths.len(),
            });
        }
        for ((&depth, &expected), values) in data.depths.iter().zip(&counts).zip(property.depths) {
            let found = values.len() as u64;
            if found != expected {
                return Err(Error::Shape {
                    property: name.into(),
                    depth,
                    expected,
                    found,
                });
            }
        }
    }
    Ok(counts)
}

/// Writes the data of a zone as DGGS-JSON, the document a request for zone data answers with
/// `application/json`:
/// `{"dggrs":URI,"zoneId":TEXT,"depths":[...],"values":{NAME:[{"depth":D,"shape":{"count":N,
/// "subZones":N},"data":[...]},...],...}}`, with one entry for each depth under each property.
///
/// No white space is written, nor a final line feed. Integers are written as integers; a value
/// of `f32` in the shortest form that reads back to the same single (`0.1`, not
/// `0.10000000149011612`), and one of `f64` in the shortest form that reads back to the same
/// double; a missing value, and one that is not finite, as `null`. The document carries neither
/// a schema of the properties nor further dimensions. The sink is flushed.
///
/// The data are held first to what the document requires, and refused before anything is
/// written: the null zone ([`Error::NullZone`]); a zone or a depth the library refuses
/// ([`Error::Grid`]); a property whose values are not given once for each depth
/// ([`Error::Depths`]); values whose number is not that of the sub-zones at their depth
/// ([`Error::Shape`]), which a client would otherwise read in silence with each value on another
/// zone; a property named `zoneID` ([`Error::ReservedName`]), which zone data as GeoJSON cannot
/// carry; two properties of one name ([`Error::RepeatedName`]).
///
/// Each value is written to the sink as it comes, in many small writes: a service passes a
/// buffered sink. A failure of the sink leaves the document truncated where it arose.
///
/// # Examples
///
/// ```
/// use rs4dggs_ogc::{Dggrs, Property, Values, ZoneData, write_dggs_json};
///
/// let dggrs = Dggrs::from_id("ISEA3H").unwrap();
/// let zone = dggrs.grid().zone_from_text("C2-23-C")?;
/// let band = [Some(0.5), Some(1.5), None, Some(0.1), Some(22.8), Some(f32::NAN), Some(6.0)];
/// let values = [Values::F32(&band)];
/// let properties = [Property::new("t", &values)];
/// let data = ZoneData::new(dggrs, zone, &[1], &properties);
/// let mut out = Vec::new();
/// write_dggs_json(&mut out, &data)?;
/// assert!(out.ends_with(br#""data":[0.5,1.5,null,0.1,22.8,null,6]}]}}"#));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn write_dggs_json<W: Write>(mut out: W, data: &ZoneData) -> Result<()> {
    let counts = check(data, RESERVED)?;
    let out = &mut out;
    out.write_all(b"{\"dggrs\":")?;
    json::string(out, data.dggrs.uri())?;
    out.write_all(b",\"zoneId\":")?;
    json::string(out, &data.dggrs.grid().text_id(data.zone))?;
    out.write_all(b",\"depths\":[")?;
    for (i, depth) in data.depths.iter().enumerate() {
        write!(out, "{}{depth}", if i == 0 { "" } else { "," })?;
    }
    out.write_all(b"],\"values\":{")?;
    for (i, property) in data.properties.iter().enumerate() {
        if i > 0 {
            out.write_all(b",")?;
        }
        json::string(out, property.name)?;
        out.write_all(b":[")?;
        for (j, ((depth, count), values)) in data
            .depths
            .iter()
            .zip(&counts)
            .zip(property.depths)
            .enumerate()
        {
            write!(
                out,
                "{}{{\"depth\":{depth},\"shape\":{{\"count\":{count},\"subZones\":{count}}},\
                 \"data\":[",
                if j == 0 { "" } else { "," }
            )?;
            match values {
                Values::U8(s) => list(out, s, |o, x| write!(o, "{x}"))?,
                Values::I8(s) => list(out, s, |o, x| write!(o, "{x}"))?,
                Values::U16(s) => list(out, s, |o, x| write!(o, "{x}"))?,
                Values::I16(s) => list(out, s, |o, x| write!(o, "{x}"))?,
                Values::U32(s) => list(out, s, |o, x| write!(o, "{x}"))?,
                Values::I32(s) => list(out, s, |o, x| write!(o, "{x}"))?,
                Values::F32(s) => list(out, s, json::f32)?,
                Values::F64(s) => list(out, s, json::f64)?,
            }
            out.write_all(b"]}")?;
        }
        out.write_all(b"]")?;
    }
    out.write_all(b"}}")?;
    out.flush()?;
    Ok(())
}

/// The values, separated by commas, each by `one` or as `null` where it is missing.
fn list<W: Write, T: Copy>(
    out: &mut W,
    values: &[Option<T>],
    one: impl Fn(&mut W, T) -> io::Result<()>,
) -> io::Result<()> {
    for (i, value) in values.iter().enumerate() {
        if i > 0 {
            out.write_all(b",")?;
        }
        match value {
            Some(x) => one(out, *x)?,
            None => out.write_all(b"null")?,
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn isea3h() -> &'static Dggrs {
        Dggrs::from_id("ISEA3H").unwrap()
    }

    fn c2_23_c() -> ZoneId {
        isea3h().grid().zone_from_text("C2-23-C").unwrap()
    }

    /// The data of `C2-23-C` at `depths`, with the properties given.
    fn data<'a>(depths: &'a [u8], properties: &'a [Property<'a>]) -> ZoneData<'a> {
        ZoneData {
            dggrs: isea3h(),
            zone: c2_23_c(),
            depths,
            properties,
        }
    }

    /// The error of writing `data`, having asserted that nothing reached the sink.
    fn refused(data: &ZoneData) -> Error {
        let mut out = Vec::new();
        let e = write_dggs_json(&mut out, data).unwrap_err();
        assert!(out.is_empty(), "{} bytes written before {e}", out.len());
        e
    }

    const SEVEN: [Option<i32>; 7] = [
        Some(0),
        Some(1),
        Some(2),
        Some(3),
        Some(4),
        Some(5),
        Some(6),
    ];

    #[test]
    fn values_tell_their_length() {
        assert_eq!(Values::I32(&SEVEN).len(), 7);
        assert!(!Values::I32(&SEVEN).is_empty());
        assert!(Values::F64(&[]).is_empty());
        assert_eq!(Values::U8(&[None, Some(1)]).len(), 2);
    }

    #[test]
    fn six_values_where_seven_are_due_are_refused() {
        let six = &SEVEN[..6];
        let e = refused(&data(
            &[1],
            &[Property {
                name: "i",
                depths: &[Values::I32(six)],
            }],
        ));
        assert!(
            matches!(&e, Error::Shape { property, depth: 1, expected: 7, found: 6 } if property == "i"),
            "{e:?}"
        );
    }

    #[test]
    fn a_property_with_one_entry_for_two_depths_is_refused() {
        let e = refused(&data(
            &[1, 2],
            &[Property {
                name: "i",
                depths: &[Values::I32(&SEVEN)],
            }],
        ));
        assert!(
            matches!(&e, Error::Depths { property, expected: 2, found: 1 } if property == "i"),
            "{e:?}"
        );
    }

    #[test]
    fn a_depth_beyond_the_finest_level_is_refused() {
        let e = refused(&data(
            &[40],
            &[Property {
                name: "i",
                depths: &[Values::I32(&SEVEN)],
            }],
        ));
        assert!(matches!(e, Error::Grid(_)), "{e:?}");
    }

    #[test]
    fn the_null_zone_is_refused() {
        let mut d = data(
            &[1],
            &[Property {
                name: "i",
                depths: &[Values::I32(&SEVEN)],
            }],
        );
        d.zone = ZoneId::NULL;
        assert!(matches!(refused(&d), Error::NullZone));
    }

    #[test]
    fn a_zone_the_grid_does_not_hold_is_refused() {
        let mut d = data(&[1], &[]);
        d.zone = ZoneId(0x0123_4567_89ab_cdef);
        assert!(matches!(refused(&d), Error::Grid(_)));
    }

    #[test]
    fn a_zone_the_grid_does_not_hold_is_refused_with_no_depth_asked() {
        let mut d = data(&[], &[]);
        d.zone = ZoneId(0x0123_4567_89ab_cdef);
        assert!(matches!(refused(&d), Error::Grid(_)));
    }

    #[test]
    fn two_properties_of_one_name_are_refused() {
        let e = refused(&data(
            &[1],
            &[
                Property {
                    name: "t",
                    depths: &[Values::I32(&SEVEN)],
                },
                Property {
                    name: "t",
                    depths: &[Values::F64(&[Some(0.5); 7])],
                },
            ],
        ));
        assert!(matches!(&e, Error::RepeatedName(n) if n == "t"), "{e:?}");
    }

    #[test]
    fn a_property_named_as_the_identifier_of_a_zone_is_refused() {
        let e = refused(&data(
            &[1],
            &[Property {
                name: "zoneID",
                depths: &[Values::I32(&SEVEN)],
            }],
        ));
        assert!(
            matches!(&e, Error::ReservedName(n) if n == "zoneID"),
            "{e:?}"
        );
    }

    #[test]
    fn a_name_that_needs_escaping_is_accepted_and_escaped() {
        let mut out = Vec::new();
        write_dggs_json(
            &mut out,
            &data(
                &[0],
                &[Property {
                    name: "a\"b",
                    depths: &[Values::U8(&[Some(9)])],
                }],
            ),
        )
        .unwrap();
        let text = String::from_utf8(out).unwrap();
        assert!(
            text.ends_with(
                r#""values":{"a\"b":[{"depth":0,"shape":{"count":1,"subZones":1},"data":[9]}]}}"#
            ),
            "{text}"
        );
    }

    #[test]
    fn the_counts_of_each_depth_are_answered_in_order() {
        let d = data(&[2, 0, 1], &[]);
        assert_eq!(check(&d, &[]).unwrap(), [13, 1, 7]);
        let d = data(
            &[1],
            &[Property {
                name: "band",
                depths: &[Values::I32(&SEVEN)],
            }],
        );
        assert!(matches!(check(&d, &["band"]), Err(Error::ReservedName(n)) if n == "band"));
        assert!(check(&d, &[]).is_ok());
    }

    #[test]
    fn data_with_no_property_or_no_depth_is_a_document_still() {
        let mut out = Vec::new();
        write_dggs_json(&mut out, &data(&[], &[])).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            r#"{"dggrs":"https://www.opengis.net/def/dggrs/OGC/1.0/ISEA3H","zoneId":"C2-23-C","depths":[],"values":{}}"#
        );
    }

    #[test]
    fn every_integer_band_writes_its_extremes_as_integers() {
        let mut out = Vec::new();
        let p = [
            Property {
                name: "u8",
                depths: &[Values::U8(&[Some(u8::MAX)])],
            },
            Property {
                name: "i8",
                depths: &[Values::I8(&[Some(i8::MIN)])],
            },
            Property {
                name: "u16",
                depths: &[Values::U16(&[Some(u16::MAX)])],
            },
            Property {
                name: "i16",
                depths: &[Values::I16(&[Some(i16::MIN)])],
            },
            Property {
                name: "u32",
                depths: &[Values::U32(&[Some(u32::MAX)])],
            },
            Property {
                name: "i32",
                depths: &[Values::I32(&[Some(i32::MIN)])],
            },
            Property {
                name: "f64",
                depths: &[Values::F64(&[Some(f64::NEG_INFINITY)])],
            },
        ];
        write_dggs_json(&mut out, &data(&[0], &p)).unwrap();
        let text = String::from_utf8(out).unwrap();
        for (name, value) in [
            ("u8", "255"),
            ("i8", "-128"),
            ("u16", "65535"),
            ("i16", "-32768"),
            ("u32", "4294967295"),
            ("i32", "-2147483648"),
            ("f64", "null"),
        ] {
            let entry = format!(
                r#""{name}":[{{"depth":0,"shape":{{"count":1,"subZones":1}},"data":[{value}]}}]"#
            );
            assert!(text.contains(&entry), "{entry} not in {text}");
        }
    }

    /// A sink that refuses every write.
    struct Refusing;

    impl Write for Refusing {
        fn write(&mut self, _: &[u8]) -> io::Result<usize> {
            Err(io::Error::other("refused"))
        }
        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn a_sink_that_fails_is_an_error_not_a_panic() {
        let e = write_dggs_json(Refusing, &data(&[1], &[])).unwrap_err();
        assert!(matches!(e, Error::Io(_)), "{e:?}");
    }
}
