//! The data of a zone, the values of its sub-zones, as DGGS-UBJSON: the document of DGGS-JSON in
//! the binary form of UBJSON.

use std::io::{self, Write};

use crate::data::{RESERVED, Values, ZoneData, check};
use crate::{Result, ubjson};

/// Writes the data of a zone as DGGS-UBJSON, the document a request for zone data answers with
/// `application/ubjson`: the document of [`write_dggs_json`](crate::write_dggs_json), member
/// for member and in the same order, in UBJSON (draft 12), whose numbers are big endian.
///
/// `dggrs` and `zoneId` are strings; each depth is an unsigned integer of 8 bits; `count` and
/// `subZones` take the smallest integer type that holds them. The values of one property at one
/// depth are written in one container, of one type:
///
/// - a band of `f32` as float32 (`d`), and one of `f64` as float64 (`D`);
/// - a band of integers as the smallest signed type that holds every value present in it: int8
///   (`i`), int16 (`I`), int32 (`l`) or int64 (`L`). The type of one band may therefore differ
///   from one zone to another; the values do not. The unsigned type of 8 bits is never used for
///   values, since some readers take a typed container of it for a string of bytes.
///
/// When every value is present and finite, the container is typed (`[$`, the type, `#`, the
/// count, then the values without their markers); otherwise it is counted (`[#`, the count, then
/// each value with its marker, and `Z`, the null of UBJSON, for a missing value or one that is
/// not finite). A band with no value present is a counted container of `Z` alone.
///
/// The data are held to what the document requires, and refused before anything is written, as
/// by [`write_dggs_json`](crate::write_dggs_json). Each value is written to the sink as it
/// comes, in many small writes: a service passes a buffered sink. The sink is flushed. A failure
/// of the sink leaves the document truncated where it arose.
///
/// # Examples
///
/// ```
/// use rs4dggs_ogc::{Dggrs, Property, Values, ZoneData, write_dggs_ubjson};
///
/// let dggrs = Dggrs::from_id("ISEA3H").unwrap();
/// let zone = dggrs.grid().zone_from_text("C2-23-C")?;
/// let band = [Some(10u16), Some(200), Some(300), Some(400), Some(500), Some(600), Some(700)];
/// let values = [Values::U16(&band)];
/// let properties = [Property::new("h", &values)];
/// let data = ZoneData::new(dggrs, zone, &[1], &properties);
/// let mut out = Vec::new();
/// write_dggs_ubjson(&mut out, &data)?;
/// // The data: a typed container of seven int16 (`[$I#i\x07`), and the closing marks.
/// assert!(out.ends_with(b"data[$I#i\x07\x00\x0a\x00\xc8\x01\x2c\x01\x90\x01\xf4\x02\x58\x02\xbc}]}}"));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn write_dggs_ubjson<W: Write>(mut out: W, data: &ZoneData) -> Result<()> {
    let counts = check(data, RESERVED)?;
    let out = &mut out;
    out.write_all(b"{")?;
    ubjson::key(out, "dggrs")?;
    ubjson::string(out, data.dggrs.uri())?;
    ubjson::key(out, "zoneId")?;
    ubjson::string(out, &data.dggrs.grid().text_id(data.zone))?;
    ubjson::key(out, "depths")?;
    out.write_all(b"[#")?;
    ubjson::count(out, data.depths.len() as u64)?;
    for &depth in data.depths {
        out.write_all(&[b'U', depth])?;
    }
    ubjson::key(out, "values")?;
    out.write_all(b"{")?;
    for property in data.properties {
        ubjson::key(out, property.name)?;
        out.write_all(b"[")?;
        for ((&depth, &count), values) in data.depths.iter().zip(&counts).zip(property.depths) {
            out.write_all(b"{")?;
            ubjson::key(out, "depth")?;
            out.write_all(&[b'U', depth])?;
            ubjson::key(out, "shape")?;
            out.write_all(b"{")?;
            ubjson::key(out, "count")?;
            ubjson::count(out, count)?;
            ubjson::key(out, "subZones")?;
            ubjson::count(out, count)?;
            out.write_all(b"}")?;
            ubjson::key(out, "data")?;
            band(out, values)?;
            out.write_all(b"}")?;
        }
        out.write_all(b"]")?;
    }
    out.write_all(b"}}")?;
    out.flush()?;
    Ok(())
}

/// One band's values at one depth, as the member `data`.
fn band<W: Write>(out: &mut W, values: &Values) -> io::Result<()> {
    match values {
        Values::U8(s) => integers(out, s),
        Values::I8(s) => integers(out, s),
        Values::U16(s) => integers(out, s),
        Values::I16(s) => integers(out, s),
        Values::U32(s) => integers(out, s),
        Values::I32(s) => integers(out, s),
        Values::F32(s) => container(out, s, b'd', f32::is_finite, |o, x| {
            o.write_all(&x.to_be_bytes())
        }),
        Values::F64(s) => container(out, s, b'D', f64::is_finite, |o, x| {
            o.write_all(&x.to_be_bytes())
        }),
    }
}

/// Integers, in the smallest signed type that holds every value present: one pass finds the
/// least and the greatest.
fn integers<W: Write, T: Copy + Into<i64>>(out: &mut W, values: &[Option<T>]) -> io::Result<()> {
    let (least, greatest) = values
        .iter()
        .flatten()
        .map(|&v| v.into())
        .fold((0, 0), |(lo, hi), v| (v.min(lo), v.max(hi)));
    let within = |lo: i64, hi: i64| least >= lo && greatest <= hi;
    let (marker, width) = if within(i8::MIN.into(), i8::MAX.into()) {
        (b'i', 1)
    } else if within(i16::MIN.into(), i16::MAX.into()) {
        (b'I', 2)
    } else if within(i32::MIN.into(), i32::MAX.into()) {
        (b'l', 4)
    } else {
        (b'L', 8)
    };
    container(
        out,
        values,
        marker,
        |_| true,
        |o, v| o.write_all(&v.into().to_be_bytes()[8 - width..]),
    )
}

/// The values in a container of the type `marker`: typed when every value is present and
/// `finite`, counted otherwise, with `Z` for each value missing or not finite. `payload` writes a
/// value without its marker.
fn container<W: Write, T: Copy>(
    out: &mut W,
    values: &[Option<T>],
    marker: u8,
    finite: impl Fn(T) -> bool,
    payload: impl Fn(&mut W, T) -> io::Result<()>,
) -> io::Result<()> {
    let present = |v: &Option<T>| v.filter(|&x| finite(x));
    if values.iter().all(|v| present(v).is_some()) {
        out.write_all(&[b'[', b'$', marker, b'#'])?;
        ubjson::count(out, values.len() as u64)?;
        for &v in values.iter().flatten() {
            payload(out, v)?;
        }
    } else {
        out.write_all(b"[#")?;
        ubjson::count(out, values.len() as u64)?;
        for v in values {
            match present(v) {
                Some(x) => {
                    out.write_all(&[marker])?;
                    payload(out, x)?;
                }
                None => out.write_all(b"Z")?,
            }
        }
    }
    Ok(())
}

#[cfg(test)]
#[path = "../tests/common/data_cases.rs"]
mod data_cases;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ubjson::read::to_json;
    use crate::{Dggrs, Error, Property};

    /// The member `data` of one band.
    fn bytes(values: Values) -> Vec<u8> {
        let mut out = Vec::new();
        band(&mut out, &values).unwrap();
        out
    }

    /// The bytes of a typed container of `count` values (below 128) of the type `marker`.
    fn typed(marker: u8, count: u8, payload: &[u8]) -> Vec<u8> {
        [&[b'[', b'$', marker, b'#', b'i', count][..], payload].concat()
    }

    #[test]
    fn every_document_reads_back_as_the_golden_of_dggs_json() {
        let goldens = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/goldens");
        let mut seen = 0;
        data_cases::each(|name, data| {
            let golden = std::fs::read_to_string(goldens.join(name)).unwrap();
            let mut out = Vec::new();
            write_dggs_ubjson(&mut out, data).unwrap();
            assert_eq!(to_json(&out), golden, "{name}");
            seen += 1;
        });
        assert_eq!(seen, 8);
    }

    #[test]
    fn a_u8_band_to_127_is_a_typed_container_of_int8() {
        let band = [Some(0), Some(5), Some(127)];
        assert_eq!(bytes(Values::U8(&band)), typed(b'i', 3, &[0, 5, 127]));
    }

    #[test]
    fn a_u8_band_holding_200_is_a_typed_container_of_int16_never_of_uint8() {
        let band = [Some(0), Some(200), Some(7)];
        let out = bytes(Values::U8(&band));
        assert_eq!(out, typed(b'I', 3, &[0, 0, 0, 200, 0, 7]));
        assert!(!out.starts_with(b"[$U"));
    }

    #[test]
    fn a_u8_band_with_a_value_missing_is_counted_with_one_null() {
        assert_eq!(
            bytes(Values::U8(&[Some(1), None, Some(2)])),
            b"[#i\x03i\x01Zi\x02"
        );
        assert_eq!(bytes(Values::U8(&[Some(200), None])), b"[#i\x02I\x00\xc8Z");
    }

    #[test]
    fn a_u32_band_takes_the_type_its_values_need() {
        assert_eq!(
            bytes(Values::U32(&[Some(1000), Some(32_767)])),
            typed(b'I', 2, &[0x03, 0xe8, 0x7f, 0xff])
        );
        assert_eq!(
            bytes(Values::U32(&[Some(1000), Some(70_000)])),
            typed(b'l', 2, &[0, 0, 0x03, 0xe8, 0, 0x01, 0x11, 0x70])
        );
        assert_eq!(
            bytes(Values::U32(&[Some(1000), Some(3_000_000_000)])),
            typed(
                b'L',
                2,
                &[
                    0, 0, 0, 0, 0, 0, 0x03, 0xe8, 0, 0, 0, 0, 0xb2, 0xd0, 0x5e, 0x00
                ]
            )
        );
    }

    #[test]
    fn a_negative_value_widens_the_type_too() {
        assert_eq!(
            bytes(Values::I16(&[Some(-129), Some(5)])),
            typed(b'I', 2, &[0xff, 0x7f, 0x00, 0x05])
        );
        assert_eq!(
            bytes(Values::I8(&[Some(-128), Some(127)])),
            typed(b'i', 2, &[0x80, 0x7f])
        );
        assert_eq!(
            bytes(Values::I32(&[Some(i32::MIN)])),
            typed(b'l', 1, &[0x80, 0, 0, 0])
        );
    }

    #[test]
    fn a_band_with_every_value_missing_is_a_counted_container_of_nulls() {
        assert_eq!(bytes(Values::I32(&[None, None, None])), b"[#i\x03ZZZ");
        assert_eq!(bytes(Values::F64(&[None, Some(f64::NAN)])), b"[#i\x02ZZ");
    }

    #[test]
    fn a_float_that_is_not_finite_is_a_null_in_a_counted_container() {
        assert_eq!(
            bytes(Values::F32(&[Some(1.0), Some(f32::NAN)])),
            b"[#i\x02d\x3f\x80\x00\x00Z"
        );
        assert_eq!(
            bytes(Values::F64(&[
                Some(f64::INFINITY),
                Some(-0.0),
                Some(f64::NEG_INFINITY)
            ])),
            b"[#i\x03ZD\x80\x00\x00\x00\x00\x00\x00\x00Z"
        );
    }

    #[test]
    fn floats_keep_their_own_width_and_their_extremes() {
        assert_eq!(
            bytes(Values::F32(&[Some(0.1)])),
            typed(b'd', 1, &[0x3d, 0xcc, 0xcc, 0xcd])
        );
        assert_eq!(
            bytes(Values::F64(&[Some(0.1)])),
            typed(b'D', 1, &0.1f64.to_be_bytes())
        );
        let extremes = [Some(f32::MAX), Some(f32::MIN_POSITIVE), Some(-0.0)];
        assert_eq!(
            to_json(&bytes(Values::F32(&extremes))),
            "[3.4028235e38,1.1754944e-38,-0]"
        );
        let extremes = [Some(f64::MAX), Some(f64::MIN_POSITIVE), Some(5e-324)];
        assert_eq!(
            to_json(&bytes(Values::F64(&extremes))),
            "[1.7976931348623157e308,2.2250738585072014e-308,5e-324]"
        );
    }

    #[test]
    fn an_f32_band_of_5536_values_is_four_bytes_a_value_typed_and_five_counted() {
        let mut band = vec![Some(22.8f32); 5536];
        let out = bytes(Values::F32(&band));
        assert_eq!(&out[..7], b"[$d#I\x15\xa0");
        assert_eq!(out.len() - 7, 22_144);
        for i in [0, 100, 5535] {
            band[i] = None;
        }
        let out = bytes(Values::F32(&band));
        assert_eq!(&out[..5], b"[#I\x15\xa0");
        assert_eq!(out.len() - 5, 27_668);
        assert_eq!(out.iter().filter(|&&b| b == b'Z').count(), 3);
    }

    #[test]
    fn the_data_are_refused_before_a_byte_as_in_dggs_json() {
        let dggrs = Dggrs::from_id("ISEA3H").unwrap();
        let zone = dggrs.grid().zone_from_text("C2-23-C").unwrap();
        let six = [Some(0.5); 6];
        let mut out = Vec::new();
        let e = write_dggs_ubjson(
            &mut out,
            &ZoneData {
                dggrs,
                zone,
                depths: &[1],
                properties: &[Property {
                    name: "t",
                    depths: &[Values::F64(&six)],
                }],
            },
        )
        .unwrap_err();
        assert!(
            matches!(
                e,
                Error::Shape {
                    expected: 7,
                    found: 6,
                    ..
                }
            ),
            "{e:?}"
        );
        assert!(out.is_empty());
    }
}
