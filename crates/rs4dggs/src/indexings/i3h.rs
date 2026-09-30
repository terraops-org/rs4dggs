//! I3H indexing: DGGAL's `I3HZone` bit class (`RI3H.ec:855-861`), one 64-bit identifier per
//! cell, `subHex` in bits 0 to 1, `rhombusIX` in bits 2 to 52, `rootRhombus` in bits 53 to 56
//! and `levelI9R` in bits 57 to 61; bits 62 and 63 belong to no field. The whole cell travels in
//! [`Address::base`] with no digits, as py4dggs's `indexings/i3h.py` carries it, since the
//! aperture-3 hierarchy is not a digit path and lives on the topology.
use crate::error::quoted;
use crate::{Address, Error, Indexing, NULL_TEXT, Result, ZoneId};

/// The finest I9R level, the level letter `Q`: `3^16 * 3^16` indices still fit the 51 bits of
/// `rhombusIX` (`RI3H.ec:858`).
const MAX_LEVEL_I9R: u64 = 16;

/// The longest text of a zone the indexing can read: a level letter, one root digit, thirteen
/// hexadecimal digits of index at level letter `Q`, a sub-hex letter and two hyphens.
const LONGEST_TEXT: usize = 18;

/// DGGAL's I3H indexing of the aperture-3 grids: a level of the rhombic I9R grid, a root
/// rhombus, an index within it and a sub-hexagon, packed into one 64-bit identifier, and its
/// text form, for instance `B6-5-A` (`RI3H.ec:855-949`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct I3h;

impl crate::interfaces::sealed::Sealed for I3h {}

/// The four fields of an identifier, `(levelI9R, rootRhombus, rhombusIX, subHex)`, read masked,
/// as the eC reads a bit class: those of `I3HZone` (`RI3H.ec:858-861`).
pub(crate) fn fields(id: u64) -> (u64, u64, u64, u64) {
    (
        (id >> 57) & 0x1F,
        (id >> 53) & 0xF,
        (id >> 2) & ((1 << 51) - 1),
        id & 0x3,
    )
}

/// An identifier built from its four fields as the eC builds an `I3HZone`.
///
/// faithful: construction is not masked. The eC compiles a bit-class constructor as unmasked
/// shifts ORed together (`ectp/src/pass16.ec:317-331`), so a value too wide for its field spills
/// into the fields above it rather than being truncated to its width, as py4dggs's `pack_i3h`
/// truncates it; the engine's children of its own null zone are built this way (`RI3H.ec:2072`).
/// With every field in range the two agree. A value shifted past bit 63 is lost, which is what
/// the eC's 64-bit shift does too.
pub(crate) fn pack(level_i9r: u64, root: u64, ix: u64, sub_hex: u64) -> u64 {
    (level_i9r << 57) | (root << 53) | (ix << 2) | sub_hex
}

/// Whether `id` names a cell: true exactly for the identifiers [`I3h::from_text`] reads back
/// from their own text.
///
/// Those are the zones the eC's `fromI9R` can build (`RI3H.ec:951-963`): an I9R level of 16 at
/// most, and either a rhombus of roots 0 to 9 with an index below `3^level * 3^level` and any
/// sub-hex, or one of the two polar pentagons, roots 10 and 11, at index 0. Bits 62 and 63 must
/// be clear, since no text can set them. Three classes that the engine answers are excluded,
/// each a documented departure: the null zone, whose text `"(null)"` reads back to it but which
/// names no cell; the sub-hexes C and D of a polar pentagon, which the engine reads from text
/// although no relation of any zone reaches them; and the identifiers of roots 10 and 11 with a
/// non-zero index, which the engine's quantiser answers along two seams of the layout and its
/// own parser refuses: root 10 along the icosahedron edge that runs over the geographic north
/// pole, root 11 along an edge in the southern hemisphere that ends at the southern polar
/// vertex, at every resolution from 2 to 21 (see `HexA3::quantize`).
pub(crate) fn is_readable(id: u64) -> bool {
    if id >> 62 != 0 {
        return false;
    }
    let (level_i9r, root, ix, sub_hex) = fields(id);
    if level_i9r > MAX_LEVEL_I9R {
        return false;
    }
    match root {
        0..=9 => ix < 3u64.pow(2 * level_i9r as u32),
        10 | 11 => ix == 0 && sub_hex <= 1,
        _ => false,
    }
}

impl Indexing for I3h {
    // faithful: DGGAL's getMaxDGGRSZoneLevel() is 33 (`RI3H.ec:42`): I9R level 16 with a sub-hex
    // of B, C or D.
    const MAX_RESOLUTION: u8 = 33;

    /// The identity: the whole cell is the address's base, with no digits.
    ///
    /// No eC counterpart: the eC keeps the cell packed and has no address to carry it in.
    fn encode(a: &Address) -> ZoneId {
        ZoneId(a.base)
    }

    /// The whole cell as the address's base, with no digits.
    ///
    /// No eC counterpart: the eC keeps the cell packed and has no address to carry it in.
    fn decode(id: ZoneId) -> Address {
        Address::new(id.0, &[])
    }

    /// faithful: `levelI9R * 2 + (subHex > 0)` (`RI3H.ec:878-881`), with one deliberate
    /// exception. The engine reads the null zone's fields as they stand and reports 63; this
    /// answers 0 for it, as for the null zone of the aperture-7 grids, so that the null zone
    /// reads alike on every grid. A caller that needs to tell a resolution-0 zone from it
    /// compares the identifier with [`ZoneId::NULL`].
    fn resolution(id: ZoneId) -> u8 {
        if id == ZoneId::NULL {
            return 0;
        }
        let (level_i9r, _, _, sub_hex) = fields(id.0);
        (level_i9r * 2 + u64::from(sub_hex > 0)) as u8
    }

    /// The root rhombus: 0 to 9, or 10 and 11 for the polar pentagons, the eC's field
    /// `rootRhombus` (`RI3H.ec:859`).
    fn base_cell(id: ZoneId) -> u64 {
        fields(id.0).1
    }

    /// faithful: `nPoints` is 5 when `rhombusIX == 0 && subHex <= 1` (`RI3H.ec:883-892`), the
    /// twelve pentagons of every level, the polar two among them.
    fn is_pentagon(id: ZoneId) -> bool {
        let (_, _, ix, sub_hex) = fields(id.0);
        ix == 0 && sub_hex <= 1
    }

    /// faithful: `getZoneID` (`RI3H.ec:935-949`), `{'A' + levelI9R}{root:X}-{rhombusIX:X}-{'A' +
    /// subHex}`, and `"(null)"` for the null zone. The level letter is the I9R level, not the
    /// resolution: `B6-5-A` is at resolution 2. Every identifier prints, one with no cell behind
    /// it included, as the eC prints its fields whatever they hold.
    fn to_text(id: ZoneId) -> String {
        if id == ZoneId::NULL {
            return NULL_TEXT.to_string();
        }
        let (level_i9r, root, ix, sub_hex) = fields(id.0);
        format!(
            "{}{root:X}-{ix:X}-{}",
            char::from(b'A' + level_i9r as u8),
            char::from(b'A' + sub_hex as u8)
        )
    }

    /// Reads a canonical text id: exactly the texts [`I3h::to_text`] prints for the identifiers
    /// that name a cell, and `"(null)"` for the null zone.
    ///
    /// The eC's `fromZoneID` (`RI3H.ec:908-931`) scans the text, rebuilds the zone and answers
    /// its null zone unless the zone prints back as the input, so it too accepts canonical text
    /// only. Two departures, each deliberate. Text that names no zone is an error here, where the
    /// engine answers its null zone, as on the aperture-7 grids and in py4dggs. And the sub-hexes
    /// C and D of a polar pentagon, such as `AA-0-C`, are refused, where the engine reads them:
    /// it leaves the pentagon's row and column unset and validates nothing for roots 10 and 11
    /// (`RI3H.ec:915-924`, `:971-972`), and no neighbour, parent, child or sub-zone of any zone
    /// is such a cell.
    fn from_text(text: &str) -> Result<ZoneId> {
        if text == NULL_TEXT {
            return Ok(ZoneId::NULL);
        }
        // The length is bounded before the text is scanned, so that the cost of refusing a long
        // input does not grow with it.
        if text.len() > LONGEST_TEXT {
            return Err(Error::InvalidZone(format!(
                "bad I3H text id {}: longer than the {LONGEST_TEXT} characters of the longest \
                 identifier",
                quoted(text)
            )));
        }
        let bad = || Error::InvalidZone(format!("bad I3H text id {}", quoted(text)));
        let mut parts = text.split('-');
        let (Some(head), Some(ix), Some(sub_hex), None) =
            (parts.next(), parts.next(), parts.next(), parts.next())
        else {
            return Err(bad());
        };
        let (Some(&letter), Some(root), &[sub_letter]) =
            (head.as_bytes().first(), head.get(1..), sub_hex.as_bytes())
        else {
            return Err(bad());
        };
        let (Ok(root), Ok(ix)) = (u64::from_str_radix(root, 16), u64::from_str_radix(ix, 16))
        else {
            return Err(bad());
        };
        let level_i9r = u64::from(letter.wrapping_sub(b'A'));
        let sub_hex = u64::from(sub_letter.wrapping_sub(b'A'));
        if level_i9r > MAX_LEVEL_I9R || sub_hex > 3 {
            return Err(bad());
        }
        if (root == 10 || root == 11) && ix == 0 && sub_hex > 1 {
            return Err(Error::InvalidZone(format!(
                "bad I3H text id {}: a polar pentagon has no sub-hex {}",
                quoted(text),
                char::from(sub_letter)
            )));
        }
        // A field too wide for its bits spills in `pack`, and the zone it then names does not
        // print back as the input, so the comparison refuses it along with every non-canonical
        // spelling (a leading zero, a lower-case digit, a sign).
        let id = pack(level_i9r, root, ix, sub_hex);
        if is_readable(id) && Self::to_text(ZoneId(id)) == text {
            Ok(ZoneId(id))
        } else {
            Err(bad())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A small, fixed pseudo-random sequence (splitmix64), so that the sweep is the same on
    /// every run.
    fn splitmix(state: &mut u64) -> u64 {
        *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = *state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Whether the text of `id` reads back to `id`.
    fn round_trips(id: u64) -> bool {
        I3h::from_text(&I3h::to_text(ZoneId(id))) == Ok(ZoneId(id))
    }

    /// Every cell of I9R levels 0 and 1, that is of resolutions 0 to 3.
    fn cells_to_resolution_3() -> Vec<u64> {
        let mut ids = Vec::new();
        for level_i9r in 0..=1 {
            for root in 0..=9 {
                for ix in 0..3u64.pow(2 * level_i9r as u32) {
                    for sub_hex in 0..=3 {
                        ids.push(pack(level_i9r, root, ix, sub_hex));
                    }
                }
            }
            for root in [10, 11] {
                for sub_hex in 0..=1 {
                    ids.push(pack(level_i9r, root, 0, sub_hex));
                }
            }
        }
        ids
    }

    /// Identifiers the engine reads from its own text, with the values it gives them.
    #[test]
    fn layout_matches_the_engine() {
        for (text, id, res) in [
            ("A0-0-A", 0x0000_0000_0000_0000, 0),
            ("A0-0-C", 0x0000_0000_0000_0002, 1),
            ("B6-5-A", 0x02C0_0000_0000_0014, 2),
            ("Q9-6949F8B1A8B78-A", 0x213A_527E_2C6A_2DE0, 32),
            ("QA-0-B", 0x2140_0000_0000_0001, 33),
        ] {
            assert_eq!(I3h::from_text(text), Ok(ZoneId(id)), "{text}");
            assert_eq!(I3h::to_text(ZoneId(id)), text);
            assert_eq!(I3h::resolution(ZoneId(id)), res, "{text}");
        }
        assert_eq!(I3h::base_cell(ZoneId(0x02C0_0000_0000_0014)), 6);
    }

    #[test]
    fn the_null_zone() {
        assert_eq!(I3h::from_text("(null)"), Ok(ZoneId::NULL));
        assert_eq!(I3h::to_text(ZoneId::NULL), "(null)");
        assert_eq!(I3h::resolution(ZoneId::NULL), 0);
        assert!(!I3h::is_pentagon(ZoneId::NULL));
        // The one identifier on which reading back and naming a cell differ.
        assert!(round_trips(ZoneId::NULL.0));
        assert!(!is_readable(ZoneId::NULL.0));
    }

    /// Resolution `r` has `10 * 3^r + 2` zones, and each prints and reads back.
    #[test]
    fn every_cell_to_resolution_3_round_trips() {
        let mut per_resolution = [0u64; 4];
        for id in cells_to_resolution_3() {
            assert!(round_trips(id), "{}", I3h::to_text(ZoneId(id)));
            assert!(is_readable(id));
            per_resolution[usize::from(I3h::resolution(ZoneId(id)))] += 1;
        }
        assert_eq!(per_resolution, [12, 32, 92, 272]);
    }

    #[test]
    fn pentagons() {
        for level_i9r in 0..=MAX_LEVEL_I9R {
            for root in 0..=11 {
                for sub_hex in 0..=1 {
                    let id = ZoneId(pack(level_i9r, root, 0, sub_hex));
                    assert!(I3h::is_pentagon(id));
                    assert!(round_trips(id.0) && is_readable(id.0));
                }
                assert!(!I3h::is_pentagon(ZoneId(pack(level_i9r, root, 0, 2))));
            }
        }
        assert!(!I3h::is_pentagon(ZoneId(pack(1, 6, 5, 0))));
    }

    #[test]
    fn malformed_and_non_canonical_text_is_refused() {
        for text in [
            "null", "", "junk", "a0-0-a", "A0-00-A", " A0-0-A", "A0-0-AX", "R0-0-A", "QC-0-A",
            "AC-0-A", "A0-1-A", "A00-0-A", "A+0-0-A", "A0-+0-A", "A0-0-E", "B6-5", "B6-5-A-A",
            "B6-5-a", "b6-5-A", "B6--5-A", "é0-0-A", "Bé-0-A",
        ] {
            assert!(
                matches!(I3h::from_text(text), Err(Error::InvalidZone(_))),
                "{text:?} accepted"
            );
        }
    }

    /// The engine reads these; this crate refuses them.
    #[test]
    fn phantom_polar_sub_hexes_are_refused() {
        for level in ['A', 'B', 'H', 'Q'] {
            for tail in ["A-0-C", "A-0-D", "B-0-C", "B-0-D"] {
                let text = format!("{level}{tail}");
                let Err(Error::InvalidZone(msg)) = I3h::from_text(&text) else {
                    panic!("{text} accepted");
                };
                assert!(msg.contains("polar pentagon"), "{msg}");
            }
        }
    }

    #[test]
    fn a_rejected_text_is_quoted_at_bounded_length() {
        let long = "A".repeat(200_000);
        let Err(Error::InvalidZone(msg)) = I3h::from_text(&long) else {
            panic!("200,000 characters accepted");
        };
        assert!(msg.len() < 200, "{} bytes: {msg}", msg.len());
        assert!(msg.contains("(200000 bytes in all)"), "{msg}");
        // The longest real text is within the bound.
        let longest = "Q9-6949F8B1A8B78-B";
        assert_eq!(longest.len(), LONGEST_TEXT);
        assert!(I3h::from_text(longest).is_ok());
    }

    /// faithful: the engine's children of its null zone are built with unmasked fields
    /// (`RI3H.ec:2072`, measured values), an index above 51 bits spilling into the root and
    /// level fields.
    #[test]
    fn pack_spills_as_the_engine_does() {
        assert_eq!(pack(0, 0, 1 << 51, 0), 1 << 53);
        let p = 3u64.pow(31);
        let ix = (3 * p).wrapping_mul(3 * p - 1);
        for (root, id) in [
            (1, 0xE4AE_48F2_5D81_FA00),
            (3, 0xE4EE_48F2_5D81_FA00),
            (5, 0xE4AE_48F2_5D81_FA00),
            (7, 0xE4EE_48F2_5D81_FA00),
            (9, 0xE5AE_48F2_5D81_FA00),
        ] {
            assert_eq!(pack(32, root, ix, 0), id, "root {root}");
        }
    }

    /// `is_readable` is the text round trip, on every class, the null zone apart.
    #[test]
    fn readable_is_the_round_trip() {
        let mut ids = cells_to_resolution_3();
        for level_i9r in 0..=MAX_LEVEL_I9R {
            for root in 0..=15 {
                for sub_hex in 0..=3 {
                    ids.push(pack(level_i9r, root, 0, sub_hex));
                }
            }
            // The unreadable north-pole class, `BA-2-A` to `KA-7355-B`, and its like.
            ids.push(pack(level_i9r, 10, 2, 0));
            ids.push(pack(level_i9r, 11, 5, 1));
            // Past the last index of the level.
            ids.push(pack(level_i9r, 3, 3u64.pow(2 * level_i9r as u32), 0));
        }
        // Level 34 and beyond, and bits 62 and 63.
        for level_i9r in 17..=31 {
            ids.push(pack(level_i9r, 4, 0, 1));
        }
        ids.push(pack(1, 6, 5, 0) | 1 << 62);
        ids.push(pack(1, 6, 5, 0) | 1 << 63);
        // Half the random values as they come, half with bits 62 and 63 clear and a level
        // letter of `A` to `Q`, so that the root and the index decide.
        let mut state = 0x5EED;
        for _ in 0..50_000 {
            let r = splitmix(&mut state);
            ids.push(r);
            ids.push((r & ((1 << 57) - 1)) | ((r >> 57) % (MAX_LEVEL_I9R + 1)) << 57);
        }
        for id in ids.into_iter().filter(|&id| id != ZoneId::NULL.0) {
            assert_eq!(
                is_readable(id),
                round_trips(id),
                "{id:#018x} {}",
                I3h::to_text(ZoneId(id))
            );
        }
        assert!(!is_readable(pack(1, 10, 2, 0)));
        assert_eq!(I3h::to_text(ZoneId(0x0340_0000_0000_0008)), "BA-2-A");
    }
}
