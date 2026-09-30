//! Errors this crate returns. Ported from py4dggs's exceptions, which are a
//! handful of stdlib exception types (`ValueError`, `KeyError`, `IndexError`,
//! `NotImplementedError`) each carrying its own formatted message; here they
//! become one `Error` enum so a caller can match on the failure kind instead
//! of parsing a string.
//!
//! Source: this crate's own API, in place of py4dggs's exception types.

use std::fmt;

/// Every way an operation in this crate can fail.
///
/// Marked `#[non_exhaustive]`: a caller's `match` needs a wildcard arm, so that a
/// variant added in a later release, for another kind of input the engine answers
/// only meaninglessly, is not a breaking change.
///
/// No error carries caller input at unbounded length. Where an error reproduces
/// text a caller supplied, as [`Error::InvalidZone`] does for a rejected text
/// identifier and [`Error::UnknownGrid`] for a grid name, it keeps at most the
/// first thirty-two characters, so that an input of any size yields an error of
/// bounded size, which a service can log, by its message or its `Debug` form,
/// without the input choosing how much is written.
#[derive(Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum Error {
    /// A zone id, text id or digit path that could not be decoded: malformed,
    /// non-canonical, or otherwise not a real zone. Mirrors py4dggs's
    /// `InvalidZoneError`, whose call sites already build their own detailed
    /// message (for instance, which text id was rejected and why), so the
    /// message travels as a plain `String` rather than a fixed template.
    InvalidZone(String),
    /// A requested resolution outside `0..=max`, where `max` is the relevant
    /// `Indexing::MAX_RESOLUTION`. `res` is the same `u8` a caller passed in,
    /// carried as it stands rather than widened, since the two fields name
    /// the one quantity, a resolution level, and read most plainly of one type.
    ResolutionOutOfRange {
        /// The resolution requested.
        res: u8,
        /// The finest resolution the grid offers.
        max: u8,
    },
    /// A grid name with no matching entry in the registry, matched against
    /// [`crate::registry::GRID_NAMES`] without regard to case. Carries the name as
    /// the caller spelt it, which the message quotes beside the names the registry
    /// does offer; a name longer than thirty-two characters is carried as its first
    /// thirty-two followed by `...`, so that neither the message nor the error's
    /// `Debug` form grows with the name.
    UnknownGrid(String),
    /// The topology has no ordered sub-zone enumeration, so a sub-zone
    /// request cannot be honoured at all. Unlike `neighbors` or the congruent
    /// hierarchy, it has no grid-agnostic substitute, since no grid-agnostic
    /// method can enumerate an arbitrary refinement in a canonical order.
    NoSubZoneOrder,
    /// Materialising every sub-zone at the requested depth would produce
    /// `count` zones, past the `limit` this crate is willing to allocate in
    /// one call.
    TooManySubZones {
        /// How many sub-zones the request would produce.
        count: u64,
        /// The most this crate materialises in one call,
        /// [`crate::grid::MAX_MATERIALISED_SUB_ZONES`].
        limit: u64,
    },
    /// `index` is not a valid position among `count` sub-zones.
    IndexOutOfRange {
        /// The index requested.
        index: u64,
        /// How many sub-zones there are at the depth requested.
        count: u64,
    },
    /// A number the caller supplied is NaN or infinite, where only a finite one
    /// has a meaning: a coordinate given to `zone_from_geo`, or a field of the
    /// [`crate::GridConfig`] given to [`crate::Grid::new`]. `name` names the
    /// parameter or field at fault, `"latitude"` or `"orientation_lat_deg"` for
    /// instance; the value itself is not carried, since a non-finite value is
    /// one of only three, and the parameter is what the caller must correct.
    ///
    /// DGGAL does not refuse a non-finite coordinate: it answers one fixed cell
    /// whatever the input, the pentagon of base cell 1, in the northern Pacific,
    /// and this crate refuses instead, so that no caller mistakes that cell for
    /// the place a coordinate names. The engine takes no configuration from a
    /// caller, so the refusal of a non-finite field has no counterpart in it.
    NonFinite {
        /// The parameter or field whose value is not finite.
        name: &'static str,
    },
}

/// The most characters of caller input an error message quotes; see [`Error`].
/// Kept private, and stated in words there, so that the bound is documented without
/// becoming a constant of the public surface.
///
/// Thirty-two exceeds the longest text a grid of this crate accepts, the
/// twenty-one characters of a Z7 identifier at resolution 19, so every input that
/// comes close to being valid is quoted whole, and only an input that no grid
/// could accept is cut short.
pub(crate) const QUOTED_INPUT_LIMIT: usize = 32;

/// `text` quoted for an error message: whole if it has at most
/// [`QUOTED_INPUT_LIMIT`] characters, and otherwise its first
/// [`QUOTED_INPUT_LIMIT`] characters followed by the length of the whole.
///
/// The cut is made on characters, never on bytes, so it cannot fall inside a
/// multi-byte character. The quoted part is escaped as Rust's `Debug` escapes a
/// string, which may lengthen each character to at most ten, so the message stays
/// bounded whatever the input.
pub(crate) fn quoted(text: &str) -> String {
    let mut chars = text.chars();
    let head: String = chars.by_ref().take(QUOTED_INPUT_LIMIT).collect();
    if chars.next().is_none() {
        format!("{head:?}")
    } else {
        format!("{head:?}... ({} bytes in all)", text.len())
    }
}

/// `text` itself if it has at most [`QUOTED_INPUT_LIMIT`] characters, and otherwise
/// its first [`QUOTED_INPUT_LIMIT`] characters followed by `...`: the bounded copy
/// of caller input that an error carries as data rather than as a finished message.
pub(crate) fn truncated(text: &str) -> String {
    let mut chars = text.chars();
    let mut head: String = chars.by_ref().take(QUOTED_INPUT_LIMIT).collect();
    if chars.next().is_some() {
        head.push_str("...");
    }
    head
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::InvalidZone(msg) => write!(f, "{msg}"),
            Error::ResolutionOutOfRange { res, max } => {
                write!(f, "resolution {res} out of range 0 to {max}")
            }
            Error::UnknownGrid(name) => write!(
                f,
                "unknown grid {name:?}; this crate offers {}",
                crate::registry::GRID_NAMES.join(", ")
            ),
            Error::NoSubZoneOrder => write!(f, "this grid has no sub-zone order"),
            Error::TooManySubZones { count, limit } => write!(
                f,
                "{count} sub-zones exceeds the {limit} materialisation limit; \
                 use count_sub_zones() for the count, or refine in smaller steps"
            ),
            Error::IndexOutOfRange { index, count } => {
                write!(
                    f,
                    "index {index} out of range: {count} sub-zone(s) available"
                )
            }
            Error::NonFinite { name } => write!(f, "{name} is not a finite number"),
        }
    }
}

impl std::error::Error for Error {}

/// This crate's own `Result`, fixed to [`Error`].
pub type Result<T> = std::result::Result<T, Error>;

#[cfg(test)]
mod tests {
    use super::*;

    /// Input up to the limit is quoted whole; beyond it, only the first
    /// characters are quoted, and the full length is stated. A cut that fell
    /// inside a multi-byte character would panic, which is why it is made on
    /// characters and tried here on one of three bytes.
    #[test]
    fn quoted_input_is_bounded() {
        let at_limit = "0".repeat(QUOTED_INPUT_LIMIT);
        assert_eq!(quoted(&at_limit), format!("{at_limit:?}"));
        let long = "0".repeat(200_002);
        let q = quoted(&long);
        assert_eq!(
            q,
            format!(
                "{:?}... (200002 bytes in all)",
                "0".repeat(QUOTED_INPUT_LIMIT)
            )
        );
        let wide = "\u{ff10}".repeat(QUOTED_INPUT_LIMIT + 1);
        assert_eq!(
            quoted(&wide),
            format!(
                "{:?}... ({} bytes in all)",
                "\u{ff10}".repeat(QUOTED_INPUT_LIMIT),
                3 * (QUOTED_INPUT_LIMIT + 1)
            )
        );
    }

    /// An unknown grid's name is kept within the same bound, and a short one whole.
    #[test]
    fn unknown_grid_name_is_bounded() {
        assert_eq!(truncated("NOPE"), "NOPE");
        assert_eq!(
            truncated(&"x".repeat(QUOTED_INPUT_LIMIT)).len(),
            QUOTED_INPUT_LIMIT
        );
        let long = truncated(&"\u{ff10}".repeat(200_000));
        assert_eq!(
            long,
            format!("{}...", "\u{ff10}".repeat(QUOTED_INPUT_LIMIT))
        );
    }

    #[test]
    fn non_finite_names_the_parameter() {
        let e = Error::NonFinite { name: "latitude" };
        assert_eq!(e.to_string(), "latitude is not a finite number");
    }
}
