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
    ///
    /// It also answers a zone that an enumeration of the zones of a level is asked to be
    /// entered after ([`crate::Zones::after`], [`crate::ZonesInBox::after`]) and that the
    /// level does not hold: a zone of another level among them, which is a real zone and
    /// still not one of that sequence.
    InvalidZone(String),
    /// A requested resolution outside `0..=max`, where `max` is the relevant
    /// `Indexing::MAX_RESOLUTION`. `res` is the same `u8` a caller passed in,
    /// carried as it stands rather than widened, since the two fields name
    /// the one quantity, a resolution level, and read most plainly of one type.
    ///
    /// From an enumeration of zones ([`crate::Grid::zones`], [`crate::Grid::zones_in_box`])
    /// `max` is the limit of the enumeration, [`crate::Grid::max_box_level`], and not the
    /// grid's last level: on the aperture-7 grids it is 14, where the grids themselves reach
    /// 19, so that a level from 15 to 19 is refused there though the grid has it.
    ResolutionOutOfRange {
        /// The resolution requested.
        res: u8,
        /// The finest resolution the grid offers, or, from an enumeration of zones, the
        /// finest level enumerated.
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
    ///
    /// No grid that this crate ships returns it: the aperture-3 and the aperture-7
    /// topologies both order their sub-zones, as the engine does. The variant stays for
    /// a topology that does not.
    NoSubZoneOrder,
    /// Materialising every sub-zone at the requested depth would produce
    /// `count` zones, past the `limit` this crate is willing to allocate in
    /// one call; or the index of a sub-zone was asked in an order so long. The
    /// limit is [`crate::grid::max_materialised_sub_zones`]:
    /// [`crate::grid::MAX_MATERIALISED_SUB_ZONES`], four million, unless the
    /// environment variable `RS4DGGS_MAX_MATERIALISED_SUB_ZONES` lowers it.
    TooManySubZones {
        /// How many sub-zones the request would produce.
        count: u64,
        /// The most this crate materialises in one call: the limit in force,
        /// [`crate::grid::max_materialised_sub_zones`].
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
    /// has a meaning: a coordinate given to `zone_from_geo`, a coordinate of the
    /// bounding box given to [`crate::Grid::zones_in_box`] or
    /// [`crate::Grid::estimate_zones_in_box`], or a field of the
    /// [`crate::GridConfig`] given to [`crate::Grid::new`]. `name` names the
    /// parameter or field at fault, `"latitude"`, `"south"` or `"orientation_lat_deg"`
    /// for instance; the value itself is not carried, since a non-finite value is
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
    /// An edge refinement above `max`, the limit in force, asked of `refined_vertices`. The
    /// limit is [`crate::grid::max_edge_refinement`]: [`crate::grid::MAX_EDGE_REFINEMENT`],
    /// 100,000, unless the environment variable `RS4DGGS_MAX_EDGE_REFINEMENT` lowers it. A
    /// refinement of 0 is never refused. The two fields are of the type the caller passed, as
    /// in [`Error::ResolutionOutOfRange`].
    ///
    /// DGGAL sets no such bound: its refined boundary takes any `int`, a hexagon at a
    /// refinement of 100,000 having 600,000 points, and beyond some 107 million the engine's
    /// own count of divisions wraps. This crate refuses instead, above the greatest
    /// refinement at which its boundary was compared with the engine's, and an operator may
    /// lower that bound further.
    EdgeRefinementOutOfRange {
        /// The edge refinement requested.
        refinement: u32,
        /// The greatest edge refinement accepted.
        max: u32,
    },
    /// A latitude of a bounding box beyond a pole, outside -90 to 90 degrees. `name` is the
    /// coordinate at fault, `"south"` or `"north"`, as in [`Error::NonFinite`]; the caller holds
    /// the box, and so its value.
    ///
    /// DGGAL answers a box of this kind with something, or with nothing; this crate refuses it
    /// instead, so that no caller takes a box beyond the poles for a box at them.
    LatitudeOutOfRange {
        /// The coordinate of the box at fault: `"south"` or `"north"`.
        name: &'static str,
    },
    /// A longitude of a bounding box outside -180 to 180 degrees. `name` is the coordinate at
    /// fault, `"west"` or `"east"`, as in [`Error::LatitudeOutOfRange`].
    ///
    /// DGGAL answers such a box if its west lies below its east, and overflows its stack if
    /// not; this crate refuses it in both cases, and a service that holds longitudes beyond
    /// that range brings them within it before it asks.
    LongitudeOutOfRange {
        /// The coordinate of the box at fault: `"west"` or `"east"`.
        name: &'static str,
    },
    /// A bounding box whose south lies above its north. A west above the east is not an error:
    /// the box then runs eastwards over the antimeridian.
    InvertedBox,
    /// Zones of more than one level, given to [`crate::Grid::compact_zones`], which compacts
    /// the zones of one level. `level` is the level of the first zone of the slice and
    /// `other` that of the first zone found of another.
    ///
    /// DGGAL's `compactZones` takes such a set and loses zones of it: a zone that has no
    /// coarser zone to give way to is dropped as soon as the set holds a finer one, so that
    /// its own answer, which is of several levels, is not a set that it compacts to itself.
    /// This crate refuses the set instead, so that no caller loses a zone in silence.
    MixedLevels {
        /// The level of the first zone.
        level: u8,
        /// The level of the first zone of another level.
        other: u8,
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
            Error::EdgeRefinementOutOfRange { refinement, max } => {
                write!(f, "edge refinement {refinement} out of range 0 to {max}")
            }
            Error::LatitudeOutOfRange { name } => {
                write!(f, "the {name} latitude is outside -90 to 90 degrees")
            }
            Error::LongitudeOutOfRange { name } => {
                write!(f, "the {name} longitude is outside -180 to 180 degrees")
            }
            Error::InvertedBox => {
                write!(f, "the box's south latitude lies above its north latitude")
            }
            Error::MixedLevels { level, other } => write!(
                f,
                "the zones are of more than one level: {level} and {other}"
            ),
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

    #[test]
    fn an_edge_refinement_out_of_range_names_both_numbers() {
        let e = Error::EdgeRefinementOutOfRange {
            refinement: 100_001,
            max: 100_000,
        };
        assert_eq!(
            e.to_string(),
            "edge refinement 100001 out of range 0 to 100000"
        );
    }

    #[test]
    fn the_box_refusals_name_the_coordinate_at_fault() {
        let e = Error::LatitudeOutOfRange { name: "north" };
        assert_eq!(
            e.to_string(),
            "the north latitude is outside -90 to 90 degrees"
        );
        let e = Error::LongitudeOutOfRange { name: "west" };
        assert_eq!(
            e.to_string(),
            "the west longitude is outside -180 to 180 degrees"
        );
        assert_eq!(
            Error::InvertedBox.to_string(),
            "the box's south latitude lies above its north latitude"
        );
    }
}
