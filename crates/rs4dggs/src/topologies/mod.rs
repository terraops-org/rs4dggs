//! Topologies: planar point <-> (base, digits), plus planar cell geometry.
mod hex_a3;
mod hex_a3_subzones;
mod hex_a7;
mod hex_a7_subzones;
pub use hex_a3::HexA3;
pub(crate) use hex_a3::pow3;
pub use hex_a7::HexA7;
