//! What the oracle suite needs to know about a grid under test: the grid itself, its
//! name in this crate, and the class name the engine gives it. Its aperture is its
//! topology's, `Topology::APERTURE`, and the suite keys everything that differs between
//! the apertures on it.
use rs4dggs::indexings::{I3h, Z7};
use rs4dggs::topologies::{HexA3, HexA7};
use rs4dggs::{Grid, Indexing, Projection, Topology};

pub trait Subject {
    type P: Projection + 'static;
    type T: Topology + 'static;
    type I: Indexing + 'static;
    /// The grid's name in this crate, as `get_grid` accepts it.
    const NAME: &'static str;
    /// The engine's class name for the same grid.
    const ORACLE: &'static str;
    fn grid() -> &'static Grid<Self::P, Self::T, Self::I>;
}

#[derive(Debug)]
pub struct Igeo7Subject;
impl Subject for Igeo7Subject {
    type P = rs4dggs::projections::Isea;
    type T = HexA7;
    type I = Z7;
    const NAME: &'static str = "IGEO7";
    const ORACLE: &'static str = dggal_oracle::IGEO7;
    fn grid() -> &'static Grid<Self::P, HexA7, Z7> {
        rs4dggs::igeo7()
    }
}

#[derive(Debug)]
pub struct Ivea7hSubject;
impl Subject for Ivea7hSubject {
    type P = rs4dggs::projections::Ivea;
    type T = HexA7;
    type I = Z7;
    const NAME: &'static str = "IVEA7H";
    const ORACLE: &'static str = dggal_oracle::IVEA7H;
    fn grid() -> &'static Grid<Self::P, HexA7, Z7> {
        rs4dggs::ivea7h()
    }
}

#[derive(Debug)]
pub struct Rtea7hSubject;
impl Subject for Rtea7hSubject {
    type P = rs4dggs::projections::Rtea;
    type T = HexA7;
    type I = Z7;
    const NAME: &'static str = "RTEA7H";
    const ORACLE: &'static str = dggal_oracle::RTEA7H;
    fn grid() -> &'static Grid<Self::P, HexA7, Z7> {
        rs4dggs::rtea7h()
    }
}

#[derive(Debug)]
pub struct Isea3hSubject;
impl Subject for Isea3hSubject {
    type P = rs4dggs::projections::Isea;
    type T = HexA3;
    type I = I3h;
    const NAME: &'static str = "ISEA3H";
    const ORACLE: &'static str = dggal_oracle::ISEA3H;
    fn grid() -> &'static Grid<Self::P, HexA3, I3h> {
        rs4dggs::isea3h()
    }
}

#[derive(Debug)]
pub struct Ivea3hSubject;
impl Subject for Ivea3hSubject {
    type P = rs4dggs::projections::Ivea;
    type T = HexA3;
    type I = I3h;
    const NAME: &'static str = "IVEA3H";
    const ORACLE: &'static str = dggal_oracle::IVEA3H;
    fn grid() -> &'static Grid<Self::P, HexA3, I3h> {
        rs4dggs::ivea3h()
    }
}

#[derive(Debug)]
pub struct Rtea3hSubject;
impl Subject for Rtea3hSubject {
    type P = rs4dggs::projections::Rtea;
    type T = HexA3;
    type I = I3h;
    const NAME: &'static str = "RTEA3H";
    const ORACLE: &'static str = dggal_oracle::RTEA3H;
    fn grid() -> &'static Grid<Self::P, HexA3, I3h> {
        rs4dggs::rtea3h()
    }
}

/// The names of all the subjects, against which the grid named by each listed divergence
/// is checked, so that an entry naming no subject cannot escape every re-check.
pub const NAMES: [&str; 6] = [
    Igeo7Subject::NAME,
    Ivea7hSubject::NAME,
    Rtea7hSubject::NAME,
    Isea3hSubject::NAME,
    Ivea3hSubject::NAME,
    Rtea3hSubject::NAME,
];
