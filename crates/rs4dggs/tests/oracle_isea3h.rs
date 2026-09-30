//! The ISEA3H grid compared operation by operation against the live DGGAL engine.
mod common;
common::oracle_suite!(common::subject::Isea3hSubject, aperture 3);
