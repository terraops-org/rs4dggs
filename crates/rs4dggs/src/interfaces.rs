//! The three pluggable parts a `Grid` composes: one `Projection`, one
//! `Topology` and one `Indexing` make a concrete DGGRS (IGEO7, IVEA7H, and so
//! on). A new grid adds a new combination of parts; none of the existing
//! parts change to accommodate it. Ported from py4dggs's `interfaces.py`,
//! whose three `Protocol` classes play the same role.
//!
//! All three traits are called through a generic parameter (`Grid<P, T, I>`),
//! never through a trait object, which is why every method here is an
//! associated function rather than one taking `&self`: there is no instance,
//! only a type that stands for one grid's projection, topology or indexing.
//!
//! `Topology` and `Indexing` each expose a few optional capabilities as
//! trait methods with a default body: `None` for the methods returning an
//! `Option`, `false` for `is_null_geometry` and for the `HAS_SUB_ZONE_ORDER`
//! constant, and an empty `Vec` for `child_digits` (which cannot return
//! `None`, so an empty list is the only way to spell "no children" once
//! `Grid` already knows there is no override to ask instead). This mirrors
//! py4dggs's own dispatch, which reads an optional method off a topology
//! with `getattr(topology, name, None)`: there, a missing attribute means
//! "this concrete topology declares no override for this capability", never
//! "the answer happens to be empty", and the same reading applies to every
//! default below. `Grid` is the one place that interprets these
//! defaults: for some capabilities it falls back to a grid-agnostic
//! implementation (the digit-path hierarchy built from
//! `Indexing::parent`/`child_digits` for `parents`/`children`/
//! `centroid_parent`/`is_centroid_child`), which answers for a topology
//! that declares no hierarchy of its own, and so for none of those this
//! crate ships, each of which ports the engine's; for `neighbors` it
//! answers none, every topology here listing its own; for others no
//! grid-agnostic answer exists, and `Grid` reports `Error::NoSubZoneOrder`
//! instead of guessing, which happens for none of those this crate ships,
//! each of which orders its sub-zones as the engine does.
//!
//! Source: py4dggs's `interfaces.py`.

use crate::error::Result;
use crate::types::{Address, GeoPoint, GridConfig, PlanarPoint, ZoneId};

/// Bars implementing [`Projection`], [`Topology`] or [`Indexing`] outside this crate.
///
/// The usual private-supertrait pattern: `sealed::Sealed` is nameable only from within
/// this crate, so only this crate's own types can satisfy it, and so only they can
/// implement the three traits above, which each carry it as a supertrait. The traits
/// themselves stay public and nameable as bounds, as [`crate::Grid`]'s own parameters
/// and the test suite's `Subject` do; a caller composes a grid from the parts this
/// crate ships, rather than supplying a fourth kind of any of the three.
pub(crate) mod sealed {
    pub trait Sealed {}
}

/// Converts between geographic coordinates and a projection's planar face
/// coordinates.
///
/// `Geometry` is the opaque, precomputed state a projection derives once
/// from a `GridConfig`: icosahedron vertex coordinates, authalic
/// coefficients, and whatever else that projection needs but does not want
/// to recompute on every call. `Grid` builds it once and passes it back into
/// `forward`/`inverse`, which keeps both pure functions of their explicit
/// arguments rather than methods on a stateful object.
///
/// Sealed by a private supertrait: only this crate's own projections implement it.
pub trait Projection: sealed::Sealed {
    /// The state this projection precomputes once for a configuration, and reads on
    /// every call thereafter.
    type Geometry: Send + Sync + 'static;

    /// Precomputes this projection's geometry for the given orientation and
    /// datum handling.
    fn build_geometry(config: &GridConfig) -> Self::Geometry;

    /// Geographic coordinates to this projection's planar point.
    fn forward(geom: &Self::Geometry, lat: f64, lon: f64) -> PlanarPoint;

    /// This projection's planar point back to geographic coordinates.
    ///
    /// `odd_grid` is the eC's own `oddGrid` argument to `inverse` (`ri5x6.ec:511`):
    /// whether the point belongs to a grid of odd resolution. The icosahedral
    /// projections read it only where they move a point that lies within a hair of a
    /// pole onto the pole itself (`fixPoles`), to choose the longitude it is given
    /// there. The caller decides it as the engine does for its grid; `Grid` documents
    /// the rule for the aperture-7 grids at the two call sites.
    fn inverse(geom: &Self::Geometry, p: PlanarPoint, odd_grid: bool) -> GeoPoint;

    /// [`Projection::inverse`] as the eC's `inverse` itself answers (`ri5x6.ec:511-556`):
    /// the latitude and the longitude in radians, in that order, or `None` where the eC's
    /// function returns `false`, which is where [`Projection::inverse`] answers the zero
    /// point. `Grid` builds the refined boundary of a zone on it, as the engine does: there
    /// a point whose inverse fails is left out, and neighbouring points are compared in
    /// radians.
    ///
    /// This is the crate's own plumbing and no part of its public interface, which answers
    /// in degrees throughout; it may change or go in any release. It is a method of this
    /// trait only because `Grid` reaches its projection through the trait alone.
    #[doc(hidden)]
    fn inverse_radians(geom: &Self::Geometry, p: PlanarPoint, odd_grid: bool)
    -> Option<(f64, f64)>;

    /// The latitude, on the sphere this projection works on, of the geodetic latitude `lat`,
    /// both in radians: the authalic latitude, by the eC's `latGeodeticToAuthalic`
    /// (`authalic.ec:9`), where the configuration asks for that conversion, exactly as
    /// [`Projection::forward`] makes it (`ri5x6.ec:252-255`, `:397`), and `lat` itself where
    /// it does not. The zones of a grid are equal in area on that sphere, so that `Grid` takes the
    /// area of a bounding box there when it counts in zones.
    ///
    /// This is the crate's own plumbing, as [`Projection::inverse_radians`] is, and no part of
    /// its public interface.
    #[doc(hidden)]
    fn sphere_latitude(geom: &Self::Geometry, lat: f64) -> f64;
}

/// Converts between a projection's planar coordinates and a topology's own
/// (base cell, direction digits) address, and supplies that address's
/// planar geometry. Carries the tessellation's aperture and cell shape; a
/// `Grid` lifts the planar geometry this trait returns into geographic
/// coordinates through its `Projection`.
///
/// Every method takes an `Address` and never a geometry parameter: unlike
/// py4dggs's `Topology` protocol, which accepts a `geom` argument on every
/// method for parity with `Projection`'s own signature, no topology here
/// actually reads it, so it is left out rather than threaded through
/// unused.
///
/// Sealed by a private supertrait: only this crate's own topologies implement it.
pub trait Topology: sealed::Sealed {
    /// The tessellation's aperture: how many children one cell splits into
    /// per resolution level (7 for the hexagonal aperture-7 grids, 3 for the
    /// hexagonal aperture-3 grids).
    const APERTURE: u8;

    /// Whether `a` lies within the bounds this topology's tables are built for:
    /// a base cell the tables have an entry for, and digits within the aperture.
    ///
    /// This is a bounds check, not a test that a cell exists behind `a`. On the
    /// aperture-7 topology, base cell 0 with the single digit 2 passes, although
    /// that digit names the deleted child of a pentagon: a `Grid` then gives the
    /// address real geometry and six neighbours, and prints it as "002", which
    /// [`Indexing::from_text`] rejects. What a true answer guarantees is only that
    /// the other methods of this trait can be called on `a` without indexing past
    /// the end of a table.
    ///
    /// Every other method here trusts the address it is given and indexes its
    /// own fixed-size tables with it, so a `Grid` asks this first on every
    /// path that reaches this trait. It has to: a [`ZoneId`] is a bare 64-bit
    /// integer, and a caller may hand over one this crate never produced,
    /// whether fabricated, corrupted or read from a stale store. The Z7
    /// packing, for one, takes the base cell from four bits and so decodes
    /// such an identifier to a base of 12 to 15, which the twelve-entry tables
    /// of an icosahedral topology have no entry for; this crate's contract,
    /// that caller input returns an error rather than panicking, covers
    /// [`crate::Grid`], [`crate::Zone`] and [`crate::AnyGrid`]. This trait and
    /// the topologies this crate implements it for are its own internals,
    /// made visible so that a caller can compose a new grid from them; the
    /// contract does not itself bind a `Topology` implementation. `Grid`'s
    /// asking this method first is what keeps the ones this crate ships
    /// inside it.
    ///
    /// Deliberately without a default body, unlike the optional capabilities
    /// below: a default of "every address is valid" would silently reopen that
    /// hole for the next topology added, and the bounds are the topology's own
    /// knowledge in any case. They cannot be hoisted into `Grid`, because the
    /// aperture-3 indexing packs an entire cell into the base field, where a
    /// bound of eleven would be meaningless.
    ///
    /// This is a different question from `is_null_geometry` below, which is
    /// about an address this topology does recognise and DGGAL still gives no
    /// geometry for.
    ///
    /// On the aperture-3 topology the answer is true exactly for the identifiers the
    /// I3H indexing reads back from its own text. The I3H null zone is not among them,
    /// so `Grid` answers it as it answers the null zone of the aperture-7 grids, with no
    /// geometry and no relations. That departs on purpose from the engine, which gives
    /// the I3H null zone a meaningless centroid and vertices and a parent chain on which
    /// it divides by zero, and from py4dggs, which mirrors the engine there.
    fn is_valid_address(a: &Address) -> bool;

    /// The address of the cell at resolution `res` that contains `p`, or `None`
    /// when there is no such cell.
    ///
    /// `None` is DGGAL's `nullZone`, and a topology returns it where the eC
    /// does and, in this crate's aperture-7 topology, in the two further cases
    /// described below. It is not an error and not a way of
    /// reporting bad input: it is the engine's own answer that the point falls
    /// in no cell at that resolution, and a caller that reproduces DGGAL has to
    /// be able to give the same answer. `Grid` turns it into [`ZoneId::NULL`],
    /// which the indexing prints as [`crate::NULL_TEXT`].
    ///
    /// It is genuinely reachable rather than defensive. An aperture-7 topology
    /// returns `None` wherever the eC's own candidate search finds no containing
    /// cell (`fromCentroid`, RI7H.ec:1489 and :1594), which over the sphere is a
    /// narrow strip along two of the icosahedron's edges, from base cell 0 to
    /// base cell 1 over the north pole and from base cell 1 to base cell 6, at
    /// the odd resolutions from 15 upwards, and also at isolated stations on the
    /// grid's mirror planes, where the input lies exactly on a cell boundary and
    /// the engine finds no cell for it. It returns `None` at resolutions 16 and
    /// 18 in the same strip, and at a few points of the strip's outer margin at
    /// 19, where the engine answers an identifier
    /// its own text-identifier parser cannot read back; the null zone is what
    /// the engine itself answers at the neighbouring odd resolutions (see "Known
    /// limitations" in the crate README). The aperture-3 topology returns it where
    /// the eC does, which a geographic point reaches near the poles at resolutions 0
    /// and 1 and beside some icosahedron edges at the coarsest resolutions, and along
    /// two seams of the layout where the engine answers an identifier it cannot read
    /// back, as [`crate::topologies::HexA3`] describes at its `quantize`.
    ///
    /// The class in which this topology answered `None` at an exact boundary tie
    /// while the engine named one of the two cells is closed. Its cause was not
    /// the side test but the polygon that test was given: the vertices this
    /// crate computed differed from the engine's by a unit or two in the last
    /// place. The vertex generation now follows the order the shipped library
    /// performs, the two sides are handed the same polygon, and at each of the
    /// three points that class held, the midpoint of the icosahedron edge from
    /// base cell 6 to base cell 10 at resolutions 9 and 19 and the south pole at
    /// resolution 15, both sides now name a cell on each of the three aperture-7
    /// grids. The polygons were found to be the same bit for bit on x86-64 Linux with glibc,
    /// where the engine and this crate share one `libm`; that says nothing of
    /// wasm32 or of any target whose floating-point functions come from another
    /// library.
    ///
    /// One case of that shape remains, and it is the engine's. Fed straight to
    /// this planar quantiser, at the bottom-right corner of the planar layout at
    /// resolutions 18 and 19, an aperture-7 topology returns `None` where the
    /// engine names a cell whose centroid lies some twenty million cell widths
    /// from the point it was asked about. No point of any oracle sample reaches
    /// it.
    ///
    /// An `Option` rather than a sentinel address, because a topology's address
    /// space has no spare value to spend: for Z7, base cell 0 with no digits is
    /// a real zone, namely base cell 0 at resolution 0, so returning that would
    /// quietly name the wrong cell instead of naming none.
    fn quantize(p: PlanarPoint, res: u8) -> Option<Address>;

    /// The planar centroid of the cell at `a`.
    ///
    /// `a` must satisfy [`Topology::is_valid_address`]; an implementation may
    /// panic otherwise, and `Grid` never calls it on any other address.
    fn planar_centroid(a: &Address) -> PlanarPoint;

    /// The planar vertices of the cell at `a`, in the topology's own winding
    /// order. Both topologies this crate ships give them in the order that runs
    /// anticlockwise as seen from outside the sphere once projected (see
    /// [`crate::Grid::vertices`]).
    ///
    /// `a` must satisfy [`Topology::is_valid_address`], as for
    /// [`Topology::planar_centroid`].
    fn planar_vertices(a: &Address) -> Vec<PlanarPoint>;

    /// The refined boundary of the cell at `a` in the projection's plane, before any
    /// projection: the ring of points `(x, y)` that the engine's `getRefinedVertices` builds
    /// first (`RI7H.ec:544`, `RI3H.ec:507`) and then takes through the inverse projection one
    /// by one, with each edge divided into `n_divisions` parts, in the engine's own sequence,
    /// as the engine builds it for a ring bound for WGS84 (the eC's `crs84` is true).
    /// `n_divisions` is the eC's `int` as it reaches that ring: `Grid` has already put, in
    /// place of an automatic refinement, the number the engine chooses for the level of the
    /// zone.
    ///
    /// `None` means "no override": this topology has no refined boundary, and `Grid` then
    /// answers the empty ring and no extent, for there is no grid-agnostic way to refine an
    /// edge as the engine does.
    ///
    /// `a` must satisfy [`Topology::is_valid_address`] and must not be one for which
    /// [`Topology::is_null_geometry`] holds; an implementation may panic otherwise, and `Grid`
    /// asks both before it calls this.
    ///
    /// This is the crate's own plumbing and no part of its public interface, which answers
    /// in degrees throughout; it may change or go in any release. It is a method of this
    /// trait only because `Grid` reaches its topology through the trait alone.
    #[doc(hidden)]
    fn planar_refined_vertices(_a: &Address, _n_divisions: i32) -> Option<Vec<(f64, f64)>> {
        None
    }

    /// Whether `a`'s geometry is degenerate: the address is representable in
    /// the packing, but there is no real geometry behind it (DGGAL's own
    /// `nullZone` outcome). Defaults to `false`, i.e. "this topology has no
    /// degenerate level, every representable address has real geometry",
    /// which is the right default for a topology that never produces one.
    /// When `true`, `Grid` returns DGGAL's null geometry (latitude/longitude
    /// 0, 0, and six zero vertices) instead of calling `planar_centroid`/
    /// `planar_vertices`. A hexagonal aperture-7 topology overrides this for
    /// Z7's twentieth packing level, which the packing can represent but
    /// which the DGGAL indexing it mirrors refuses to turn into a real
    /// zone. A hexagonal aperture-3 topology does not need to: every I3H
    /// identifier with no cell behind it, its null zone included, already
    /// fails [`Topology::is_valid_address`], and `Grid` gives it the same null
    /// geometry by that route.
    fn is_null_geometry(_a: &Address) -> bool {
        false
    }

    /// The neighbours of `a`, in the topology's own order. `None` means "no
    /// override", and `Grid` then answers no neighbours: there is no grid-agnostic
    /// neighbour search to fall back to. Both topologies this crate ships override it
    /// with a port of the engine's own neighbour search, which lists the neighbours in
    /// the engine's order (see [`crate::Grid::neighbors`]).
    fn neighbors(_a: &Address) -> Option<Vec<Address>> {
        None
    }

    /// Exact geometric parents of `a`, the primary parent first, when the
    /// hierarchy is not a simple digit-path truncation. `None` means "no
    /// override": `Grid` then falls back to the digit-path parent built from
    /// `Indexing::parent`. Both topologies this crate ships override it with
    /// a port of the engine's own search, for neither hierarchy is a tree: on
    /// the hexagonal aperture-3 topology one cell can have three parents, and
    /// on the hexagonal aperture-7 topology two, the first of which is the
    /// digit-path parent, save at some identifiers of the broken seams (see
    /// [`crate::Grid::parents`]).
    fn parents(_a: &Address) -> Option<Vec<Address>> {
        None
    }

    /// Exact geometric children of `a`, in the topology's own order. `None`
    /// means "no override": `Grid` then answers one child for each digit of
    /// `Indexing::child_digits`. See `parents` for which topologies override
    /// it, and [`crate::Grid::children`].
    fn children(_a: &Address) -> Option<Vec<Address>> {
        None
    }

    /// DGGAL's centroid parent of `a`, when the hierarchy is geometric: the
    /// first of its parents that is itself a centroid child, `Some(None)`
    /// where none is. That is not the parent of which `a` is the centroid
    /// child; see [`crate::Grid::centroid_parent`] and `parents`. `None`
    /// means "no override": `Grid` then answers `Indexing::parent`.
    fn centroid_parent(_a: &Address) -> Option<Option<Address>> {
        None
    }

    /// Whether `a` is a centroid child, when the hierarchy is geometric. See
    /// `parents`. `None` means "no override": `Grid` then answers that it is,
    /// as every zone of a digit-path hierarchy is its one parent's.
    fn is_centroid_child(_a: &Address) -> Option<bool> {
        None
    }

    /// Whether this topology overrides `count_sub_zones`, `first_sub_zone` and `sub_zones`
    /// below, the three methods a sub-zone order must supply. `Grid` reads this before validating a
    /// sub-zone request's depth, so a depth-dependent generator is never asked about depth 0 (a
    /// case it is not designed for) or about a depth so deep it would exceed the indexing's own
    /// maximum resolution; both are rejected at the `Grid` level instead. A topology without an
    /// ordered sub-zone enumeration leaves this `false`. `sub_zone_at_index`,
    /// the fourth sub-zone method below, is a further optimisation a topology may leave
    /// unoverridden even where this is `true`; it is not read here.
    const HAS_SUB_ZONE_ORDER: bool = false;

    /// The number of sub-zones `depth` resolution levels below `a`. `None`
    /// when `HAS_SUB_ZONE_ORDER` is `false`, which `Grid` reports as
    /// `Error::NoSubZoneOrder` rather than approximating: there is no
    /// grid-agnostic way to enumerate an arbitrary refinement.
    fn count_sub_zones(_a: &Address, _depth: u8) -> Option<u64> {
        None
    }

    /// The first (index 0), in this topology's own ordering, sub-zone
    /// `depth` levels below `a`. See `count_sub_zones`.
    ///
    /// `None` is "no order", never "no zone". A topology whose order has an
    /// entry that is no zone, as the aperture-7 one has in the broken seams,
    /// answers `Some` of an address that its indexing encodes as the null
    /// identifier, here and in `sub_zones` and `sub_zone_at_index`.
    fn first_sub_zone(_a: &Address, _depth: u8) -> Option<Address> {
        None
    }

    /// Every sub-zone `depth` levels below `a`, in this topology's own
    /// order. See `count_sub_zones`.
    ///
    /// A topology that bounds the list on a direct call reads the limit in force,
    /// [`crate::grid::max_materialised_sub_zones`], and no other: `Grid` reads the same
    /// one before it asks, and takes a `None` here for "no order".
    fn sub_zones(_a: &Address, _depth: u8) -> Option<Vec<Address>> {
        None
    }

    /// The sub-zone at `index`, `depth` levels below `a`, in this topology's own order, found
    /// without materialising [`Topology::sub_zones`]'s whole sequence. Optional even where
    /// `HAS_SUB_ZONE_ORDER` is `true`: `None` means "no override", and `Grid` then falls back to
    /// indexing into `sub_zones(a, depth)`, which its own cap,
    /// [`crate::grid::max_materialised_sub_zones`], still guards. A topology overrides this when it can jump to `index` directly, as the
    /// aperture-3 and the aperture-7 hexagons do through the engine's own generators, at a cost
    /// measured in the number of scanlines rather than in `index` itself.
    ///
    /// Where a topology does override it, `None` also answers an `index` at or beyond the count,
    /// or a `depth` reaching past the maximum resolution, exactly as `first_sub_zone` and
    /// `sub_zones` do; `Grid` checks both before it ever asks, so its own answers never meet that
    /// case, but a direct caller of this trait method cannot tell "no override" apart from those
    /// boundary answers, and must not read either meaning into the other.
    fn sub_zone_at_index(_a: &Address, _depth: u8, _index: u64) -> Option<Address> {
        None
    }

    /// Where `sub` stands in the order of the sub-zones of `parent`, found without
    /// materialising [`Topology::sub_zones`]'s whole sequence. Optional even where
    /// `HAS_SUB_ZONE_ORDER` is `true`: `None` means "no override", and `Grid` then builds the
    /// order and searches it, as the aperture-3 hexagons have it do. The aperture-7 hexagons
    /// override it with the engine's own walk, which generates the centroids of the order from
    /// its head until it meets that of `sub`, at a cost that grows with the index.
    ///
    /// `Some(None)` is "no sub-zone of it, or one that the walk does not find". `Some(Some(i))`
    /// is the walk's index and no more: `Grid` answers it only once it has confirmed that the
    /// entry of the order at `i` is `sub`, and a direct caller of this trait method must do
    /// the same before it relies on the index.
    ///
    /// Where a topology does override it, `None` also answers an address that names no zone, a
    /// `sub` finer than the maximum resolution, and an order longer than
    /// [`crate::grid::max_materialised_sub_zones`], so that the walk is bounded as the list
    /// is; `Grid` refuses each of them before it ever asks, and a direct caller cannot tell
    /// them apart from "no override", as with `sub_zone_at_index`.
    fn sub_zone_index(_parent: &Address, _sub: &Address) -> Option<Option<u64>> {
        None
    }

    /// The number of cells along the side of a root rhombus at `level`. The zones of a level
    /// form a lattice: ten root rhombi, 0 to 9, each of `lattice_edge(level)` rows by as many
    /// columns of cells, and two polar roots, 10 and 11, of one cell each. Taken rhombus by
    /// rhombus, row by row and cell by cell, and then the two polar roots, the zones those cells
    /// host are every zone of the level, each once, in the order of DGGAL's `listZones`.
    ///
    /// A topology that does not order its zones so leaves this and the three methods below as
    /// they are: a side of nought, no zone in any cell and no place for any zone, so that a walk
    /// of its lattice visits nothing. Nought also answers a level beyond the topology's finest.
    #[doc(hidden)]
    fn lattice_edge(_level: u8) -> u64 {
        0
    }

    /// The zones that the cell at `row` and `col` of the root rhombus `root` (0 to 9) hosts at
    /// `level`, in the order of DGGAL's `listZones`, each with its key in that order: the keys
    /// ascend along the walk that [`Topology::lattice_edge`] describes. A key orders the zones
    /// of one level alone and is not to be compared across levels. Nothing for a cell the
    /// lattice does not have: a root above 9, a row or a column at or beyond the side, a level
    /// beyond the finest.
    #[doc(hidden)]
    fn cell_zones(_level: u8, _root: u8, _row: u64, _col: u64) -> Vec<(Address, u64)> {
        Vec::new()
    }

    /// The zones of the polar root `root`, 10 or 11, at `level`, as [`Topology::cell_zones`]
    /// gives those of a cell. Nothing for any other root.
    #[doc(hidden)]
    fn polar_zones(_level: u8, _root: u8) -> Vec<(Address, u64)> {
        Vec::new()
    }

    /// Where the lattice holds the zone `a`: its root, the row and the column of the cell that
    /// hosts it (both nought under a polar root), and its key, as [`Topology::cell_zones`] and
    /// [`Topology::polar_zones`] give it. `None` for an address the lattice does not hold.
    #[doc(hidden)]
    fn locate(_a: &Address) -> Option<(u8, u64, u64, u64)> {
        None
    }

    /// DGGAL's `compactZones` of a set of zones: the zones that stand for the set, a coarser
    /// zone in place of those of its sub-zones that it alone holds wherever the set has all
    /// of them, in the engine's order of answer. A zone given twice counts once. `None` for
    /// a topology that has no compaction.
    ///
    /// Every address must satisfy [`Topology::is_valid_address`], and the addresses are meant
    /// to be of one level; `Grid` checks both before it asks. The cost is that of the slice.
    #[doc(hidden)]
    fn compact(_zones: &[Address]) -> Option<Vec<Address>> {
        None
    }
}

/// Converts between a topology's (base cell, digits) address, a packed
/// [`ZoneId`], and that id's canonical text form.
///
/// Sealed by a private supertrait: only this crate's own indexings implement it.
pub trait Indexing: sealed::Sealed {
    /// The finest resolution this indexing can pack into a `ZoneId`. A
    /// `Grid` rejects any request for a coarser-than-0 or finer-than-this
    /// resolution before it ever reaches a `Projection` or `Topology`.
    const MAX_RESOLUTION: u8;

    /// Packs an address into its zone id.
    fn encode(a: &Address) -> ZoneId;

    /// Unpacks a zone id back into its address.
    fn decode(id: ZoneId) -> Address;

    /// The resolution level encoded in `id`.
    fn resolution(id: ZoneId) -> u8;

    /// The base cell encoded in `id`.
    fn base_cell(id: ZoneId) -> u64;

    /// Whether `id` names a pentagon (an icosahedron vertex cell, which has
    /// five neighbours and children rather than the hexagonal six/seven).
    fn is_pentagon(id: ZoneId) -> bool;

    /// The canonical text form of `id`.
    fn to_text(id: ZoneId) -> String;

    /// Parses a canonical text id back into a `ZoneId`, rejecting anything
    /// malformed or non-canonical.
    ///
    /// It also rejects text that names a resolution finer than
    /// [`Indexing::MAX_RESOLUTION`], even where the packing could hold it, so that
    /// no text names a zone the engine gives no geometry; for Z7 that is the
    /// twentieth level, which [`Indexing::to_text`] still prints and this does not
    /// read back. It quotes the rejected text at bounded length, never whole, since
    /// the text is the caller's and may be of any size.
    fn from_text(text: &str) -> Result<ZoneId>;

    /// The digit-path parent of `id`, when the identifier is a digit path
    /// (dropping the last digit is exactly one level up). `None` means "this
    /// indexing has no digit path to offer". `Grid` reads this only for a
    /// topology that declares no hierarchy of its own: a geometric hierarchy
    /// lives on the `Topology`, never here, precisely because it does not
    /// fit this digit-truncation shape. That is so of an aperture-3 grid's,
    /// whose identifiers are no digit paths, and of an aperture-7 grid's,
    /// where the Z7 digit path gives the first of a zone's parents and not
    /// the second, save at some identifiers of the broken seams (see
    /// [`crate::Grid::parents`]). `None` here and `id` already being at resolution 0 both
    /// mean "no parent"; telling those two apart, if it matters, is `Grid`'s
    /// job, not this trait's.
    fn parent(_id: ZoneId) -> Option<ZoneId> {
        None
    }

    /// The child digits `id` can be extended with, one resolution level
    /// down, when the identifier is a digit path. An empty `Vec` (rather
    /// than `None`, which this method cannot return) is reached only once
    /// `Grid` already knows that the topology declares no children of its
    /// own, so there is no "no override" case left to signal here; it omits
    /// any digit whose child would be the collapsed pentagon child, matching
    /// py4dggs's own `child_digits`, which never returns a digit with no
    /// zone behind it.
    fn child_digits(_id: ZoneId) -> Vec<u8> {
        Vec::new()
    }
}
