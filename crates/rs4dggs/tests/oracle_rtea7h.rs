//! The RTEA7H grid compared operation by operation against the live DGGAL engine.
mod common;
common::oracle_suite!(common::subject::Rtea7hSubject, aperture 7);
