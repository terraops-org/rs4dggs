//! The RTEA3H grid compared operation by operation against the live DGGAL engine.
mod common;
common::oracle_suite!(common::subject::Rtea3hSubject, aperture 3);
