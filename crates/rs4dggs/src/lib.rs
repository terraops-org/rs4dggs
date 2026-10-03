//! rs4dggs: DGGAL's discrete global grid systems in pure Rust, ported from py4dggs.
#![doc = include_str!("../README.md")]
#![forbid(unsafe_code)]
#![deny(missing_docs)]

mod disk;
pub mod error;
mod facts;
mod fivebysix;
pub mod grid;
pub mod indexings;
pub mod interfaces;
pub(crate) mod math;
pub mod projections;
pub mod registry;
pub mod topologies;
pub mod types;
pub mod zone;
mod zones;

pub use disk::{Disk, Rings};
pub use error::{Error, Result};
pub use grid::Grid;
pub use interfaces::{Indexing, Projection, Topology};
pub use registry::{
    AnyGrid, GRID_NAMES, Igeo7, Isea3h, Ivea3h, Ivea7h, Rtea3h, Rtea7h, get_grid, igeo7, isea3h,
    ivea3h, ivea7h, rtea3h, rtea7h,
};
pub use types::{Address, Extent, GeoPoint, GridConfig, NULL_TEXT, PlanarPoint, ZoneId};
pub use zone::Zone;
pub use zones::{Zones, ZonesInBox};

// The suites that ask the crate's own trait methods, which no code outside the crate can reach,
// run as modules of the library's own tests. Their files stay under `tests/`, and each names the
// crate `rs4dggs`, as code outside it would.
#[cfg(test)]
extern crate self as rs4dggs;
#[cfg(test)]
#[path = "../tests/compact_zones.rs"]
mod compact_zones_suite;
#[cfg(test)]
#[path = "../tests/lattice.rs"]
mod lattice_suite;
#[cfg(test)]
#[path = "../tests/zones_in_box.rs"]
mod zones_in_box_suite;
#[cfg(test)]
#[path = "../tests/zones.rs"]
mod zones_suite;

// The oracle suites, one for each grid, and what they share.
#[cfg(all(test, feature = "oracle"))]
#[path = "../tests/common/mod.rs"]
mod common;
#[cfg(all(test, feature = "oracle"))]
#[path = "../tests/oracle_igeo7.rs"]
mod oracle_igeo7;
#[cfg(all(test, feature = "oracle"))]
#[path = "../tests/oracle_isea3h.rs"]
mod oracle_isea3h;
#[cfg(all(test, feature = "oracle"))]
#[path = "../tests/oracle_ivea3h.rs"]
mod oracle_ivea3h;
#[cfg(all(test, feature = "oracle"))]
#[path = "../tests/oracle_ivea7h.rs"]
mod oracle_ivea7h;
#[cfg(all(test, feature = "oracle"))]
#[path = "../tests/oracle_rtea3h.rs"]
mod oracle_rtea3h;
#[cfg(all(test, feature = "oracle"))]
#[path = "../tests/oracle_rtea7h.rs"]
mod oracle_rtea7h;
