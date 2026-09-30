//! The IVEA7H grid compared operation by operation against the live DGGAL engine.
mod common;
common::oracle_suite!(common::subject::Ivea7hSubject, aperture 7);
