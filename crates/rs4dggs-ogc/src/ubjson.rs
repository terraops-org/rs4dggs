//! The primitives of UBJSON (Universal Binary JSON, draft 12, <https://ubjson.org>): lengths and
//! counts, the keys of an object, and strings. Every number of UBJSON is big endian.

use std::io::{self, Write};

/// `n`, a length or a count, as a marked integer of the smallest type that holds it: `i` (int8)
/// to 127, `U` (uint8) to 255, `I` (int16), `l` (int32), `L` (int64). UBJSON has no unsigned
/// 64-bit type, so that a number beyond the range of int64 is refused, with nothing written.
pub(crate) fn count<W: Write>(out: &mut W, n: u64) -> io::Result<()> {
    if let Ok(v) = i8::try_from(n) {
        out.write_all(&[b'i', v.to_be_bytes()[0]])
    } else if let Ok(v) = u8::try_from(n) {
        out.write_all(&[b'U', v])
    } else if let Ok(v) = i16::try_from(n) {
        out.write_all(b"I")?;
        out.write_all(&v.to_be_bytes())
    } else if let Ok(v) = i32::try_from(n) {
        out.write_all(b"l")?;
        out.write_all(&v.to_be_bytes())
    } else if let Ok(v) = i64::try_from(n) {
        out.write_all(b"L")?;
        out.write_all(&v.to_be_bytes())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("{n} exceeds the largest integer UBJSON can write"),
        ))
    }
}

/// `s` as the key of an object: its length in bytes, then its bytes, with no `S` marker.
pub(crate) fn key<W: Write>(out: &mut W, s: &str) -> io::Result<()> {
    count(out, s.len() as u64)?;
    out.write_all(s.as_bytes())
}

/// `s` as a string value: the marker `S`, then as a key.
pub(crate) fn string<W: Write>(out: &mut W, s: &str) -> io::Result<()> {
    out.write_all(b"S")?;
    key(out, s)
}

/// A reader of UBJSON for the tests, which turns a document back into JSON text: integers as
/// integers, float32 and float64 in the shortest form of their own width, `Z` as `null`, strings
/// escaped as the writer of JSON escapes them. It panics on what it cannot read: an unknown
/// marker, a document cut short, bytes after the document.
#[cfg(test)]
pub(crate) mod read {
    use crate::json;

    /// `bytes`, one UBJSON document, as JSON text.
    pub(crate) fn to_json(bytes: &[u8]) -> String {
        let mut reader = Reader { bytes, at: 0 };
        let mut out = Vec::new();
        let marker = reader.byte();
        reader.value(marker, &mut out);
        assert_eq!(reader.at, bytes.len(), "bytes after the document");
        String::from_utf8(out).expect("UTF-8")
    }

    struct Reader<'a> {
        bytes: &'a [u8],
        at: usize,
    }

    impl Reader<'_> {
        fn take<const N: usize>(&mut self) -> [u8; N] {
            let taken = self.bytes[self.at..self.at + N].try_into().unwrap();
            self.at += N;
            taken
        }

        fn byte(&mut self) -> u8 {
            self.take::<1>()[0]
        }

        fn peek(&self) -> u8 {
            self.bytes[self.at]
        }

        fn integer(&mut self, marker: u8) -> Option<i64> {
            Some(match marker {
                b'i' => i8::from_be_bytes(self.take()).into(),
                b'U' => self.byte().into(),
                b'I' => i16::from_be_bytes(self.take()).into(),
                b'l' => i32::from_be_bytes(self.take()).into(),
                b'L' => i64::from_be_bytes(self.take()),
                _ => return None,
            })
        }

        fn length(&mut self) -> usize {
            let marker = self.byte();
            let n = self.integer(marker).expect("a length is an integer");
            usize::try_from(n).expect("a length is not negative")
        }

        fn text(&mut self) -> &str {
            let n = self.length();
            let text = std::str::from_utf8(&self.bytes[self.at..self.at + n]).expect("UTF-8");
            self.at += n;
            text
        }

        fn value(&mut self, marker: u8, out: &mut Vec<u8>) {
            if let Some(n) = self.integer(marker) {
                out.extend(n.to_string().bytes());
                return;
            }
            match marker {
                b'Z' => out.extend(b"null"),
                b'd' => json::f32(out, f32::from_be_bytes(self.take())).unwrap(),
                b'D' => json::f64(out, f64::from_be_bytes(self.take())).unwrap(),
                b'S' => json::string(out, self.text()).unwrap(),
                b'{' => {
                    out.push(b'{');
                    let mut first = true;
                    while self.peek() != b'}' {
                        if !first {
                            out.push(b',');
                        }
                        first = false;
                        json::string(out, self.text()).unwrap();
                        out.push(b':');
                        let marker = self.byte();
                        self.value(marker, out);
                    }
                    self.at += 1;
                    out.push(b'}');
                }
                b'[' => {
                    let typed = (self.peek() == b'$').then(|| {
                        self.at += 1;
                        self.byte()
                    });
                    let count = (self.peek() == b'#').then(|| {
                        self.at += 1;
                        self.length()
                    });
                    assert!(typed.is_none() || count.is_some(), "a type without a count");
                    out.push(b'[');
                    let mut i = 0;
                    while count.map_or_else(|| self.peek() != b']', |n| i < n) {
                        if i > 0 {
                            out.push(b',');
                        }
                        let marker = typed.unwrap_or_else(|| self.byte());
                        self.value(marker, out);
                        i += 1;
                    }
                    if count.is_none() {
                        self.at += 1;
                    }
                    out.push(b']');
                }
                _ => panic!("an unknown marker {marker:#04x} before byte {}", self.at),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_reader_gives_json_text() {
        // An object of a string, a plain array of an int16 and a null, a counted array of a
        // float64, and a typed array of float32.
        let doc = b"{i\x01aSi\x02\"\\i\x01b[I\x01\x00Z]i\x01c[#i\x01D\x3f\xf8\x00\x00\x00\x00\x00\x00i\x01d[$d#i\x02\x3d\xcc\xcc\xcd\xc0\x00\x00\x00}";
        assert_eq!(
            read::to_json(doc),
            r#"{"a":"\"\\","b":[256,null],"c":[1.5],"d":[0.1,-2]}"#
        );
    }

    #[test]
    #[should_panic(expected = "an unknown marker 0x54")]
    fn the_reader_refuses_a_marker_it_does_not_know() {
        read::to_json(b"[T]");
    }

    fn bytes(write: impl FnOnce(&mut Vec<u8>) -> io::Result<()>) -> Vec<u8> {
        let mut out = Vec::new();
        write(&mut out).unwrap();
        out
    }

    #[test]
    fn a_count_takes_the_smallest_integer_type_that_holds_it() {
        for (n, want) in [
            (0u64, &[0x69, 0x00][..]),
            (5, &[0x69, 0x05]),
            (127, &[0x69, 0x7f]),
            (128, &[0x55, 0x80]),
            (255, &[0x55, 0xff]),
            (256, &[0x49, 0x01, 0x00]),
            (32_767, &[0x49, 0x7f, 0xff]),
            (32_768, &[0x6c, 0x00, 0x00, 0x80, 0x00]),
            (2_147_483_647, &[0x6c, 0x7f, 0xff, 0xff, 0xff]),
            (
                2_147_483_648,
                &[0x4c, 0x00, 0x00, 0x00, 0x00, 0x80, 0x00, 0x00, 0x00],
            ),
            (
                i64::MAX as u64,
                &[0x4c, 0x7f, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff],
            ),
        ] {
            assert_eq!(bytes(|o| count(o, n)), want, "{n}");
        }
    }

    #[test]
    fn a_count_beyond_the_signed_64_bit_range_is_refused() {
        let mut out = Vec::new();
        let e = count(&mut out, i64::MAX as u64 + 1).unwrap_err();
        assert_eq!(e.kind(), io::ErrorKind::InvalidInput);
        assert!(out.is_empty());
    }

    #[test]
    fn a_key_is_its_length_and_its_bytes() {
        assert_eq!(
            bytes(|o| key(o, "dggrs")),
            [0x69, 0x05, 0x64, 0x67, 0x67, 0x72, 0x73]
        );
    }

    #[test]
    fn a_string_is_marked_then_written_as_a_key() {
        assert_eq!(
            bytes(|o| string(o, "C2-23-C")),
            [0x53, 0x69, 0x07, 0x43, 0x32, 0x2d, 0x32, 0x33, 0x2d, 0x43]
        );
        // The length counts bytes, not characters.
        assert_eq!(bytes(|o| string(o, "é")), [0x53, 0x69, 0x02, 0xc3, 0xa9]);
    }
}
