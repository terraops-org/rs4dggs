//! rs4dggs: DGGAL's discrete global grid systems in pure Rust, ported from py4dggs.
#![doc = include_str!("../README.md")]
#![forbid(unsafe_code)]
#![deny(missing_docs)]

mod disk;
pub mod error;
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

pub use disk::{Disk, Rings};
pub use error::{Error, Result};
pub use grid::Grid;
pub use interfaces::{Indexing, Projection, Topology};
pub use registry::{
    AnyGrid, GRID_NAMES, Igeo7, Isea3h, Ivea3h, Ivea7h, Rtea3h, Rtea7h, get_grid, igeo7, isea3h,
    ivea3h, ivea7h, rtea3h, rtea7h,
};
pub use types::{Address, GeoPoint, GridConfig, NULL_TEXT, PlanarPoint, ZoneId};
pub use zone::Zone;
