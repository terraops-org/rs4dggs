//! rs4dggs-ogc: the encodings of OGC API - DGGS over the grids of `rs4dggs`.
//!
//! Each encoding is a pure function that writes one document to a sink (anything that implements
//! [`std::io::Write`]): it reads no environment, opens no connection and builds no document in
//! memory. A service passes the sink of its response; a test passes a `Vec<u8>`.
//!
//! Every public function returns [`Result`]. A fault that can be known before the first byte is
//! written (the shape of the values, the names of the properties) refuses the document before
//! anything reaches the sink; a fault that arises while writing (the sink itself, a zone the
//! library refuses) leaves the document truncated where it arose.
#![forbid(unsafe_code)]
#![deny(missing_docs)]

mod error;
pub mod geometry;
mod zone_geojson;
pub use zone_geojson::{
    Coordinates, Written, ZoneGeometry, ZoneListGeoJson, write_zone_feature, write_zone_geometry,
};
// The writers of primitives, which the writers of documents call.
mod json;
mod ubjson;

pub use error::{Error, Result};

mod link;
mod zones;
pub use link::{Link, LinkTemplate, rel};
pub use zones::{ZoneListJson, write_zone_list_uint64};

mod dggrs;
pub use dggrs::Dggrs;

mod data;
pub use data::{Property, Values, ZoneData, write_dggs_json};

mod data_ubjson;
pub use data_ubjson::write_dggs_ubjson;

mod data_geojson;
pub use data_geojson::write_zone_data_geojson;

// The README's examples, compiled and run with the doctests.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
struct ReadmeExamples;

/// The media types of the encodings, as a service names them in `Content-Type`.
pub mod media_type {
    /// JSON: a list of zones, DGGS-JSON, a definition.
    pub const JSON: &str = "application/json";
    /// GeoJSON: a list of zones with their geometry, or zone data as features.
    pub const GEO_JSON: &str = "application/geo+json";
    /// UBJSON: DGGS-UBJSON, the binary form of DGGS-JSON.
    pub const UBJSON: &str = "application/ubjson";
    /// A list of zones as 64-bit integers.
    pub const UINT64: &str = "application/x-binary";
}
