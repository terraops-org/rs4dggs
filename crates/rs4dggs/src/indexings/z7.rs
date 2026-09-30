//! Z7 indexing: a 64-bit id = 4-bit base cell (0-11) + twenty 3-bit direction digits (0-6;
//! 7 = unused). Resolution is implicit: the first 7 digit. Port of py4dggs `indexings/z7.py`;
//! bit layout verified against `dggal-v0.06/src/dggrs/RI7H_Z7.ec`.
use crate::error::quoted;
use crate::{Address, Error, Indexing, NULL_TEXT, Result, ZoneId};

/// Base cell 0 and all twenty digit slots = 7.
const DEFAULT: u64 = 0x0FFF_FFFF_FFFF_FFFF;
/// Direction-digit slots in the packing. Level 20 is representable but not a DGGAL zone.
pub const PACKING_LEVELS: u8 = 20;

/// DGGAL's Z7 indexing of the aperture-7 grids: a base cell and up to twenty direction
/// digits packed into one 64-bit identifier, with the resolution implicit in the first
/// unused digit slot, and its text form, the base cell in two decimal digits followed by
/// one digit per resolution level (`RI7H_Z7.ec`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct Z7;

impl crate::interfaces::sealed::Sealed for Z7 {}

/// The shift of the digit of resolution `res`, as the eC's loops count it down from `19 * 3`,
/// three bits a level (`RI7H_Z7.ec:42-46`).
#[inline]
fn bit_offset(res: u8) -> u32 {
    (u32::from(PACKING_LEVELS) - u32::from(res)) * 3
}
/// The digit of resolution `res`, read as the eC reads it,
/// `(int)((ancestry & (7LL << shift)) >> shift)` (`RI7H_Z7.ec:48`).
#[inline]
fn get_dir(bits: u64, res: u8) -> u8 {
    ((bits >> bit_offset(res)) & 0x7) as u8
}
/// `bits` with the digit of resolution `res` replaced by `d`. The eC writes a digit into a
/// cleared slot, `ancestry |= (uint64)(c - '0') << shift` (`RI7H_Z7.ec:444`), and clears one
/// with `ancestry & ~(7LL << shift)` (`RI7H_Z7.ec:76`).
#[inline]
fn set_dir(bits: u64, res: u8, d: u8) -> u64 {
    let off = bit_offset(res);
    (bits & !(0x7u64 << off)) | (u64::from(d) << off)
}
/// The direction digit deleted under a pentagon rooted at `base`: 2 north (0-5), 5 south (6-11).
/// eC `fromTextID`: `c == (south ? '5' : '2')` with `south = root >= 6` (`RI7H_Z7.ec:430-438`).
#[inline]
fn pentagon_omit(base: u64) -> u8 {
    if base <= 5 { 2 } else { 5 }
}

impl Indexing for Z7 {
    // faithful: DGGAL getMaxDGGRSZoneLevel() == 19; level 20 is nullZone (RI7H_Z7.ec:348-352).
    const MAX_RESOLUTION: u8 = 19;

    /// The base cell, the eC's field `rootPentagon` (`RI7H_Z7.ec:16`).
    fn base_cell(id: ZoneId) -> u64 {
        (id.0 >> 60) & 0xF
    }

    /// faithful: with one deliberate exception. The eC's `Z7Zone::level`
    /// special-cases the sentinel and returns -1 (RI7H_Z7.ec:37-38), whereas this
    /// answers 0 for it, because every digit slot of the all-ones identifier holds
    /// the terminator.
    ///
    /// The signature stays `u8` rather than widening to carry -1. A negative
    /// resolution is not a value a caller can act on, and carrying it would push a
    /// signed type through `Grid::resolution` and `Zone::resolution` for no gain,
    /// while the null zone is already identifiable by its identifier, which is
    /// [`ZoneId::NULL`], and by its text, which is [`NULL_TEXT`]. A caller that
    /// needs to tell a resolution-0 zone from the sentinel compares the identifier.
    fn resolution(id: ZoneId) -> u8 {
        for r in 1..=PACKING_LEVELS {
            if get_dir(id.0, r) == 7 {
                return r - 1;
            }
        }
        PACKING_LEVELS
    }

    fn is_pentagon(id: ZoneId) -> bool {
        // The sentinel is not a pentagon. Without this it would reach the `res == 0`
        // arm below, since its resolution reads as 0, and claim to be one. The eC's
        // `Z7Zone::nPoints` special-cases it and returns 0 (RI7H_Z7.ec:109-110),
        // which the live engine confirms as `countZoneEdges` of 0: the sentinel has
        // no edges at all, so it is neither a pentagon nor a hexagon.
        if id == ZoneId::NULL {
            return false;
        }
        let res = Self::resolution(id);
        res == 0 || (1..=res).all(|r| get_dir(id.0, r) == 0)
    }

    /// The digits a child may append: 0 to 6, less the one a pentagon omits, as the eC's
    /// `fromTextID` accepts them (`RI7H_Z7.ec:435-448`).
    fn child_digits(id: ZoneId) -> Vec<u8> {
        if Self::is_pentagon(id) {
            let omit = pentagon_omit(Self::base_cell(id));
            (0..7).filter(|&d| d != omit).collect()
        } else {
            (0..7).collect()
        }
    }

    /// Congruent digit-path parent (drop the last digit).
    ///
    /// No eC counterpart: deliberately so, since DGGAL's Z7 hierarchy is the geometric one.
    fn parent(id: ZoneId) -> Option<ZoneId> {
        let res = Self::resolution(id);
        if res == 0 {
            None
        } else {
            Some(ZoneId(set_dir(id.0, res, 7)))
        }
    }

    /// Packs the address into the Z7 bit layout, as the eC's `fromTextID` packs a zone
    /// (`RI7H_Z7.ec:434-455`): the base cell in the top four bits, one digit a slot from the
    /// top, and the terminator 7 in every slot after the last. Nothing is validated here:
    /// `Grid` only ever encodes addresses that the topology produced or accepted.
    ///
    /// It does not panic, in a debug build or a release one. An [`Address`] holds
    /// at most [`crate::types::ADDRESS_CAPACITY`] digits, twenty, which is the
    /// number of digit slots in the packing, so every digit has a slot; a
    /// twenty-first is refused by [`Address::new`] and [`Address::push`], which
    /// panic in every build.
    ///
    /// Nor are the fields checked, and an address that does not fit them yields a
    /// well-formed identifier of another cell rather than an error. A base cell
    /// above 15 keeps only its lowest four bits, so that 16 encodes as base cell 0
    /// and 17 as base cell 1. A digit above 7 keeps its lowest three bits in its own
    /// slot and carries the rest into the next coarser slot, or into the base cell
    /// for the first digit. A digit of 7 is the terminator, so the path reads as
    /// ending before it. This follows from the arithmetic, and was confirmed on a
    /// sample of such addresses in both builds.
    fn encode(a: &Address) -> ZoneId {
        let mut bits = (DEFAULT & !(0xFu64 << 60)) | (a.base << 60);
        for (i, &d) in a.digits().iter().enumerate() {
            bits = set_dir(bits, i as u8 + 1, d);
        }
        ZoneId(bits)
    }

    /// The base cell and the digits up to the first terminator, as the eC's `getTextID`
    /// reads them (`RI7H_Z7.ec:471-479`).
    fn decode(id: ZoneId) -> Address {
        let res = Self::resolution(id);
        let mut a = Address::new(Self::base_cell(id), &[]);
        for r in 1..=res {
            a.push(get_dir(id.0, r));
        }
        a
    }

    /// The text id, as the eC's `getTextID` prints it (`RI7H_Z7.ec:461-482`).
    fn to_text(id: ZoneId) -> String {
        // The null sentinel prints as "(null)", exactly as the eC's
        // `getZoneTextID` does: it special-cases nullZone before it decodes any
        // bit field. Without this branch the all-ones sentinel stringifies into
        // the plausible-looking "15", the base field's own all-ones value, which
        // hides the one signal DGGAL gives that a zone does not exist. Verified
        // against the live engine, which prints "(null)" for this identifier on
        // the Z7 grid; it is what py4dggs already does for the aperture-3 grids
        // (`indexings/i3h.py`), and it becomes load-bearing here now that
        // `HexA7::quantize` can produce the sentinel near the poles.
        if id == ZoneId::NULL {
            return NULL_TEXT.to_string();
        }
        let a = Self::decode(id);
        // Sized once, and the base cell written as the eC's `sprintf(zoneID, "%02d",
        // rootPentagon)` writes it (`RI7H_Z7.ec:471`): the base field has four bits,
        // so it is at most 15 and always two decimal digits.
        let mut s = String::with_capacity(2 + a.len());
        s.push(char::from(b'0' + (a.base / 10) as u8));
        s.push(char::from(b'0' + (a.base % 10) as u8));
        for &d in a.digits() {
            s.push(char::from(b'0' + d));
        }
        s
    }

    /// Reads a text id, as the eC's `fromTextID` does (`RI7H_Z7.ec:418-459`), with the
    /// departures noted below.
    fn from_text(text: &str) -> Result<ZoneId> {
        // "(null)" round-trips back to the sentinel, so `from_text(to_text(v))`
        // holds for it too. Only that exact spelling does: DGGAL is laxer and
        // answers nullZone for any unparseable string at all, where this port
        // reports an error instead. That is the same deviation py4dggs records
        // for the aperture-3 grids, where the engine signals failure with a
        // sentinel return and this library signals it with a typed error.
        if text == NULL_TEXT {
            return Ok(ZoneId::NULL);
        }
        let bytes = text.as_bytes();
        // A deliberate departure from the engine. DGGAL's parser accepts twenty
        // digits, the whole packing, and so reads a level-20 identifier, which it
        // then gives no geometry: the zero point as its centroid, six zero vertices
        // and no neighbours, answers a caller cannot tell from a genuine cell at
        // (0, 0). This refuses any text longer than an identifier at
        // `MAX_RESOLUTION`, the two characters of the base cell and nineteen
        // digits, so that no text names a zone the engine gives no geometry. The
        // length is bounded here, before the text is scanned or any digit is
        // collected, so that the cost of refusing a long input does not grow with it.
        let longest = 2 + usize::from(Self::MAX_RESOLUTION);
        if bytes.len() > longest {
            return Err(Error::InvalidZone(format!(
                "bad Z7 text id {}: longer than the {longest} characters of an \
                 identifier at the finest resolution, {}",
                quoted(text),
                Self::MAX_RESOLUTION
            )));
        }
        let bad = || Error::InvalidZone(format!("bad Z7 text id {}", quoted(text)));
        if bytes.len() < 2 || !bytes.iter().all(u8::is_ascii_digit) {
            return Err(bad());
        }
        let base = u64::from((bytes[0] - b'0') * 10 + (bytes[1] - b'0'));
        let digits: Vec<u8> = bytes[2..].iter().map(|c| c - b'0').collect();
        if base > 11 || digits.iter().any(|&d| d > 6) {
            return Err(bad());
        }
        let omit = pentagon_omit(base);
        for &d in &digits {
            if d == 0 {
                continue;
            }
            if d == omit {
                return Err(Error::InvalidZone(format!(
                    "non-canonical pentagon child: digit {d} invalid under base {base}"
                )));
            }
            break;
        }
        Ok(Self::encode(&Address::new(base, &digits)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Address, Error, Indexing, ZoneId};

    #[test]
    fn max_resolution_matches_dggal() {
        assert_eq!(Z7::MAX_RESOLUTION, 19);
        assert_eq!(PACKING_LEVELS, 20);
    }
    #[test]
    fn encode_decode_roundtrip() {
        let a = Address::new(3, &[6, 5, 3, 1]);
        assert_eq!(Z7::decode(Z7::encode(&a)), a);
    }
    #[test]
    fn lisbon_res5_value_and_text() {
        let id = Z7::encode(&Address::new(0, &[6, 4, 1, 5, 6]));
        assert_eq!(id, ZoneId(0x0d0d_dfff_ffff_ffff));
        assert_eq!(Z7::to_text(id), "0064156");
        assert_eq!(Z7::from_text("0064156").unwrap(), id);
    }
    #[test]
    fn resolution_and_pentagon() {
        let v = Z7::encode(&Address::new(0, &[0, 0]));
        assert_eq!(Z7::resolution(v), 2);
        assert!(Z7::is_pentagon(v));
        let h = Z7::encode(&Address::new(0, &[6, 4, 1, 5, 6]));
        assert!(!Z7::is_pentagon(h));
        assert_eq!(Z7::child_digits(h).len(), 7);
        assert_eq!(Z7::child_digits(v).len(), 6);
        assert!(Z7::is_pentagon(Z7::encode(&Address::new(7, &[]))));
    }
    #[test]
    fn pentagon_child_digits_omit_hemisphere_digit() {
        assert_eq!(
            Z7::child_digits(Z7::encode(&Address::new(0, &[0]))),
            vec![0, 1, 3, 4, 5, 6]
        );
        assert_eq!(
            Z7::child_digits(Z7::encode(&Address::new(6, &[]))),
            vec![0, 1, 2, 3, 4, 6]
        );
        assert_eq!(
            Z7::child_digits(Z7::encode(&Address::new(6, &[1]))),
            vec![0, 1, 2, 3, 4, 5, 6]
        );
    }
    /// The null sentinel prints as "(null)" and parses back, as DGGAL does.
    /// Without the branch it would print "15", the base field's all-ones value,
    /// which reads like an ordinary identifier.
    #[test]
    fn null_zone_prints_and_parses_as_null_text() {
        assert_eq!(Z7::to_text(ZoneId::NULL), NULL_TEXT);
        assert_eq!(Z7::from_text(NULL_TEXT).unwrap(), ZoneId::NULL);
        // Only the exact canonical spelling is the sentinel. DGGAL is laxer and
        // answers nullZone for any garbage; this port reports an error, the same
        // deviation py4dggs records for the aperture-3 grids.
        for t in ["null", "(NULL)", "(null", "", "  (null)  "] {
            assert!(
                Z7::from_text(t).is_err(),
                "{t:?} should not be the sentinel"
            );
        }
    }

    #[test]
    fn parent_drops_last_digit() {
        let id = Z7::from_text("0064156").unwrap();
        assert_eq!(Z7::to_text(Z7::parent(id).unwrap()), "006415");
        assert_eq!(Z7::parent(Z7::from_text("05").unwrap()), None);
    }
    #[test]
    fn bit63_base_cells_round_trip() {
        for base in 8..=11u64 {
            let id = Z7::encode(&Address::new(base, &[1, 2, 3]));
            assert_eq!(id.0 >> 63, 1);
            assert_eq!(Z7::base_cell(id), base);
            assert_eq!(Z7::from_text(&Z7::to_text(id)).unwrap(), id);
        }
    }
    #[test]
    fn from_text_rejects_malformed_and_non_canonical() {
        for bad in [
            "",
            "0",
            "12",
            "0a",
            "007",
            "00123456789012345678901",
            "-01",
            "００",
        ] {
            assert!(
                matches!(Z7::from_text(bad), Err(Error::InvalidZone(_))),
                "{bad:?}"
            );
        }
        assert!(Z7::from_text("0002").is_err()); // north pentagon, deleted digit 2
        assert!(Z7::from_text("06005").is_err()); // south pentagon, deleted digit 5
        assert!(Z7::from_text("0012").is_ok()); // first non-zero digit 1: hexagon, 2 fine after
        assert!(Z7::from_text("06002").is_ok()); // 2 is valid under a south base
        assert_eq!(Z7::resolution(Z7::from_text("00").unwrap()), 0);
    }
    /// The packing represents a twentieth level, and prints it, but the parser
    /// refuses to read it back: the engine gives such a zone no geometry, and this
    /// crate names no such zone from text. So `from_text(to_text(id))` fails for a
    /// level-20 identifier, deliberately, and the refusal names the bound.
    #[test]
    fn twenty_digits_representable() {
        let a = Address::new(1, &[1; 20]);
        let id = Z7::encode(&a);
        assert_eq!(Z7::resolution(id), 20);
        let text = Z7::to_text(id);
        assert_eq!(text, format!("01{}", "1".repeat(20)));
        match Z7::from_text(&text) {
            Err(Error::InvalidZone(msg)) => assert!(
                msg.contains("longer than the 21 characters") && msg.contains(&text),
                "{msg}"
            ),
            other => panic!("level-20 text {text} gave {other:?}"),
        }
        // Nineteen digits, the finest resolution, still read back.
        let finest = Z7::encode(&Address::new(1, &[1; 19]));
        assert_eq!(Z7::from_text(&Z7::to_text(finest)).unwrap(), finest);
    }

    /// A rejected text is quoted in the error at bounded length, however long it is,
    /// and the length check comes before the text is scanned.
    #[test]
    fn a_rejected_text_is_quoted_at_bounded_length() {
        let long = "0".repeat(200_002);
        let Err(Error::InvalidZone(msg)) = Z7::from_text(&long) else {
            panic!("200,002 digits accepted");
        };
        assert!(msg.len() < 200, "{} bytes: {msg}", msg.len());
        assert!(msg.contains("(200002 bytes in all)"), "{msg}");
        // Short enough to pass the length check, and still quoted whole.
        let Err(Error::InvalidZone(msg)) = Z7::from_text("0a") else {
            panic!("0a accepted");
        };
        assert_eq!(msg, "bad Z7 text id \"0a\"");
    }
}
