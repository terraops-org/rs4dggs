//! The primitives of JSON text: strings, and numbers in their shortest exact form or to a fixed
//! number of decimals.

use std::fmt::{Display, LowerExp};
use std::io::{self, Write};

/// `s` as a JSON string, quoted. The quotation mark and the reverse solidus are escaped with a
/// reverse solidus, and every character below U+0020 as `\u00` and two lower-case hexadecimal
/// digits; every other character is written as it is, in UTF-8.
pub(crate) fn string<W: Write>(out: &mut W, s: &str) -> io::Result<()> {
    out.write_all(b"\"")?;
    let bytes = s.as_bytes();
    let mut start = 0;
    for (i, &b) in bytes.iter().enumerate() {
        if b == b'"' || b == b'\\' || b < 0x20 {
            out.write_all(&bytes[start..i])?;
            match b {
                b'"' => out.write_all(b"\\\"")?,
                b'\\' => out.write_all(b"\\\\")?,
                _ => write!(out, "\\u{b:04x}")?,
            }
            start = i + 1;
        }
    }
    out.write_all(&bytes[start..])?;
    out.write_all(b"\"")
}

/// A number in the shortest form that reads back to the same value: plain decimals where
/// `plain`, the exponent form otherwise, `null` where it is not finite (JSON has no form for an
/// infinity or for a value that is not a number).
fn number<W: Write, T: Display + LowerExp>(
    out: &mut W,
    x: T,
    finite: bool,
    plain: bool,
) -> io::Result<()> {
    if !finite {
        out.write_all(b"null")
    } else if plain {
        write!(out, "{x}")
    } else {
        write!(out, "{x:e}")
    }
}

/// `x` in the shortest form that reads back to the same double: plain decimals where its
/// magnitude is nought or lies from 1e-5 to below 1e17 (`0.00001`, `-0`, `31170676883138.758`),
/// the exponent form beyond (`9.999e-6`, `1e17`), and `null` where it is not finite.
pub(crate) fn f64<W: Write>(out: &mut W, x: f64) -> io::Result<()> {
    let m = x.abs();
    number(out, x, x.is_finite(), m == 0.0 || (1e-5..1e17).contains(&m))
}

/// `x` in the shortest form that reads back to the same single (`0.1`, not the double nearest
/// to it), by the rule of [`f64`].
pub(crate) fn f32<W: Write>(out: &mut W, x: f32) -> io::Result<()> {
    let m = x.abs();
    number(out, x, x.is_finite(), m == 0.0 || (1e-5..1e17).contains(&m))
}

/// `x` rounded to `decimals` places, never in the exponent form, and `null` where it is not
/// finite. A value that rounds to nought is written without a sign: `-0.0000001` to three places
/// is `0.000`, not `-0.000`.
pub(crate) fn fixed<W: Write>(out: &mut W, x: f64, decimals: u8) -> io::Result<()> {
    if !x.is_finite() {
        return out.write_all(b"null");
    }
    let p = usize::from(decimals);
    let s = format!("{x:.p$}");
    let s = match s.strip_prefix('-') {
        Some(t) if t.bytes().all(|b| b == b'0' || b == b'.') => t,
        _ => &s,
    };
    out.write_all(s.as_bytes())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn text(write: impl FnOnce(&mut Vec<u8>) -> io::Result<()>) -> String {
        let mut out = Vec::new();
        write(&mut out).unwrap();
        String::from_utf8(out).unwrap()
    }

    fn s(x: &str) -> String {
        text(|o| string(o, x))
    }

    fn d(x: f64) -> String {
        text(|o| f64(o, x))
    }

    fn f(x: f32) -> String {
        text(|o| f32(o, x))
    }

    #[test]
    fn quote_and_backslash_are_escaped() {
        assert_eq!(s(r#"a"b\c"#), r#""a\"b\\c""#);
        assert_eq!(s(""), r#""""#);
    }

    #[test]
    fn every_control_character_is_written_as_a_unicode_escape() {
        assert_eq!(s("\n"), r#""\u000a""#);
        for b in 0u8..0x20 {
            let c = char::from(b);
            assert_eq!(
                s(&format!("x{c}y")),
                format!("\"x\\u00{b:02x}y\""),
                "byte {b:#04x}"
            );
        }
    }

    #[test]
    fn other_characters_are_written_as_they_are() {
        assert_eq!(s("é/"), "\"é/\"");
        assert_eq!(s("\u{1F30D}"), "\"\u{1F30D}\"");
        assert_eq!(s("\u{7f}"), "\"\u{7f}\"");
    }

    #[test]
    fn doubles_take_the_plain_form_from_1e_minus_5_to_below_1e17() {
        for (x, want) in [
            (0.0, "0"),
            (-0.0, "-0"),
            (3.0, "3"),
            (0.1, "0.1"),
            (1e-5, "0.00001"),
            (9.999e-6, "9.999e-6"),
            (1e16, "10000000000000000"),
            (1e17, "1e17"),
            (1e300, "1e300"),
            (5e-324, "5e-324"),
            (-1.5e-7, "-1.5e-7"),
            (f64::MAX, "1.7976931348623157e308"),
            (31170676883138.758, "31170676883138.758"),
            (34.78016915103061, "34.78016915103061"),
        ] {
            assert_eq!(d(x), want, "{x:e}");
        }
    }

    #[test]
    fn a_number_that_is_not_finite_is_null() {
        for x in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert_eq!(d(x), "null");
        }
        for x in [f32::NAN, f32::INFINITY, f32::NEG_INFINITY] {
            assert_eq!(f(x), "null");
        }
        assert_eq!(text(|o| fixed(o, f64::NAN, 3)), "null");
        assert_eq!(text(|o| fixed(o, f64::INFINITY, 3)), "null");
    }

    #[test]
    fn singles_are_written_as_themselves() {
        assert_eq!(f(0.1), "0.1");
        assert_eq!(f(22.8), "22.8");
        assert_eq!(f(f32::MAX), "3.4028235e38");
        assert_eq!(f(-0.0), "-0");
        assert_eq!(f(1e-5), "0.00001");
    }

    #[test]
    fn fixed_decimals_never_write_a_negative_zero() {
        assert_eq!(text(|o| fixed(o, -0.0000001, 3)), "0.000");
        assert_eq!(text(|o| fixed(o, -0.0, 2)), "0.00");
        assert_eq!(text(|o| fixed(o, -0.0006, 3)), "-0.001");
        assert_eq!(text(|o| fixed(o, 12.34567, 2)), "12.35");
        assert_eq!(text(|o| fixed(o, 7.0, 0)), "7");
    }

    /// SplitMix64, so that the sample is the same everywhere.
    fn patterns(seed: u64) -> impl Iterator<Item = u64> {
        let mut state = seed;
        std::iter::repeat_with(move || {
            state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
            let mut z = state;
            z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
            z ^ (z >> 31)
        })
    }

    fn double_reads_back(x: f64) {
        let t = d(x);
        let y: f64 = t.parse().unwrap_or_else(|e| panic!("{t}: {e}"));
        assert_eq!(y.to_bits(), x.to_bits(), "{t}");
    }

    fn single_reads_back(x: f32) {
        let t = f(x);
        let y: f32 = t.parse().unwrap_or_else(|e| panic!("{t}: {e}"));
        assert_eq!(y.to_bits(), x.to_bits(), "{t}");
    }

    #[test]
    fn every_double_written_reads_back_to_the_same_bits() {
        let below = |x: f64| f64::from_bits(x.to_bits() - 1);
        let extremes = [
            0.0,
            f64::MIN_POSITIVE,
            below(f64::MIN_POSITIVE),
            f64::from_bits(1),
            f64::MAX,
            f64::EPSILON,
            1e-5,
            below(1e-5),
            1e17,
            below(1e17),
            1.0,
            below(1.0),
            0.1,
            1e300,
        ];
        for x in extremes {
            double_reads_back(x);
            double_reads_back(-x);
        }
        let mut n = 0;
        for bits in patterns(0x5eed_0001) {
            let x = f64::from_bits(bits);
            if x.is_finite() {
                double_reads_back(x);
                n += 1;
                if n == 100_000 {
                    break;
                }
            }
        }
    }

    #[test]
    fn every_single_written_reads_back_to_the_same_bits() {
        let below = |x: f32| f32::from_bits(x.to_bits() - 1);
        let extremes = [
            0.0,
            f32::MIN_POSITIVE,
            below(f32::MIN_POSITIVE),
            f32::from_bits(1),
            f32::MAX,
            f32::EPSILON,
            1e-5,
            below(1e-5),
            1e17,
            below(1e17),
            0.1,
            22.8,
        ];
        for x in extremes {
            single_reads_back(x);
            single_reads_back(-x);
        }
        let mut n = 0;
        for bits in patterns(0x5eed_0002) {
            let x = f32::from_bits(bits as u32);
            if x.is_finite() {
                single_reads_back(x);
                n += 1;
                if n == 100_000 {
                    break;
                }
            }
        }
    }
}
