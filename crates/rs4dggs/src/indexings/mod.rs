//! Indexings: (base, digits) <-> packed u64 <-> text id.
mod i3h;
mod z7;
pub use i3h::I3h;
pub(crate) use i3h::{fields, is_readable, pack};
pub use z7::Z7;
