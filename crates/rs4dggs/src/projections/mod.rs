//! Projections: geographic lat/lon <-> planar (face, x, y) in DGGAL's oblique 5x6 layout.
mod icovertex;
pub mod isea;
pub mod ivea;
pub mod rtea;
pub use icovertex::IcoGeometry;
pub use isea::Isea;
pub use ivea::Ivea;
pub use rtea::Rtea;

// The eC's `wrapLonAt` belongs to no projection: it lives beside `wrapLon`, whose constants it
// shares, and the rest of the crate reaches it by this name. `Grid` calls it where it builds the
// refined boundary of a zone.
pub(crate) use icovertex::wrap_lon_at;
