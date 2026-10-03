//! What can prevent a document from being written, or from being written whole.

use std::fmt;

/// Why a document was refused, or left unfinished.
///
/// A refusal of the values or of the names ([`Error::Shape`], [`Error::Depths`],
/// [`Error::ReservedName`], [`Error::RepeatedName`]) comes before the first byte is written. So
/// does the null zone ([`Error::NullZone`]) in the writers of zone data, of a zone's geometry or
/// feature, and of the list of 64-bit integers, and a geometry with a coordinate that is not
/// finite ([`Error::NotFinite`]) in the writers of a zone's geometry or feature. In a list given
/// one zone at a time either comes before any byte of the zone concerned, and so does
/// [`Error::NotFinite`] in zone data as GeoJSON, which draws one feature at a time: the zones
/// before it are written and the document is left unfinished. An error of the sink
/// ([`Error::Io`]) or of the library ([`Error::Grid`]) may come in the middle of a document,
/// which is then left truncated. The message of these two says what failed, and
/// [`std::error::Error::source`] gives the cause, so that a chain of errors names it once.
#[derive(Debug)]
#[non_exhaustive]
pub enum Error {
    /// The sink refused a write.
    Io(std::io::Error),
    /// The grid library refused a zone or a request.
    Grid(rs4dggs::Error),
    /// The null zone was given where a zone is required.
    NullZone,
    /// The values of a property at one relative depth do not match the number of sub-zones.
    Shape {
        /// The name of the property.
        property: String,
        /// The relative depth at which the count differs.
        depth: u8,
        /// The number of sub-zones at that depth.
        expected: u64,
        /// The number of values given.
        found: u64,
    },
    /// A property does not carry one set of values for each relative depth asked.
    Depths {
        /// The name of the property.
        property: String,
        /// The number of depths asked.
        expected: usize,
        /// The number of sets of values given.
        found: usize,
    },
    /// A property bears a name the encoding reserves for itself.
    ReservedName(String),
    /// Two properties bear the same name.
    RepeatedName(String),
    /// A coordinate of a geometry made from a zone's ring is not finite, and GeoJSON has no form
    /// for it.
    NotFinite,
}

/// This crate's own `Result`, fixed to [`Error`].
pub type Result<T> = std::result::Result<T, Error>;

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Io(_) => write!(f, "the document could not be written"),
            Error::Grid(_) => write!(f, "the grid refused the request"),
            Error::NullZone => write!(f, "the null zone stands where a zone is required"),
            Error::Shape {
                property,
                depth,
                expected,
                found,
            } => write!(
                f,
                "property {property:?} at relative depth {depth}: {found} values for \
                 {expected} sub-zones"
            ),
            Error::Depths {
                property,
                expected,
                found,
            } => write!(
                f,
                "property {property:?}: {found} sets of values for {expected} relative depths"
            ),
            Error::ReservedName(name) => {
                write!(f, "the property name {name:?} is reserved by the encoding")
            }
            Error::RepeatedName(name) => {
                write!(f, "the property name {name:?} is given more than once")
            }
            Error::NotFinite => write!(f, "a coordinate of the geometry is not finite"),
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Error::Io(e) => Some(e),
            Error::Grid(e) => Some(e),
            _ => None,
        }
    }
}

impl From<std::io::Error> for Error {
    fn from(e: std::io::Error) -> Self {
        Error::Io(e)
    }
}

impl From<rs4dggs::Error> for Error {
    fn from(e: rs4dggs::Error) -> Self {
        Error::Grid(e)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::error::Error as _;

    #[test]
    fn the_sink_and_the_library_are_the_sources() {
        let io: Error = std::io::Error::other("full").into();
        assert!(matches!(io, Error::Io(_)));
        assert!(io.source().is_some());
        let grid: Error = rs4dggs::Error::NoSubZoneOrder.into();
        assert!(matches!(grid, Error::Grid(_)));
        assert!(grid.source().is_some());
        assert!(Error::NullZone.source().is_none());
        // The cause is the source's to tell, not the message's.
        assert_eq!(io.to_string(), "the document could not be written");
        assert!(
            !grid
                .to_string()
                .contains(&rs4dggs::Error::NoSubZoneOrder.to_string())
        );
    }

    #[test]
    fn a_name_is_shown_escaped() {
        let e = Error::RepeatedName("a\"b\n".into());
        assert_eq!(
            e.to_string(),
            "the property name \"a\\\"b\\n\" is given more than once"
        );
    }
}
