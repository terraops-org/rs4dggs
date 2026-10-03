//! A `Grid` is one DGGRS: a `Projection`, a `Topology` and an `Indexing`
//! composed into a working grid. Every operation that needs all three at once
//! lives here.
//!
//! Ported from py4dggs's `grid.py`, keeping its fallback structure: a topology
//! may expose an exact hierarchy or an ordered sub-zone enumeration, and where
//! it does not, the `Grid` either falls back to something grid-agnostic (the
//! digit path of the identifier) or refuses the request outright
//! (`Error::NoSubZoneOrder`), never approximating. Neighbours are the
//! topology's own: every topology here ports the engine's neighbour search.
//! So is the hierarchy: every topology here ports the engine's parents and
//! children, and the digit path answers for none of the grids this crate ships.
//! So is the sub-zone order: every topology here ports the engine's generator
//! of it, and the refusal is returned by none of the grids this crate ships.
//!
//! Source: py4dggs's `grid.py`.

use std::marker::PhantomData;
use std::sync::OnceLock;

use crate::disk::Disk;
use crate::facts;
use crate::math;
use crate::projections::wrap_lon_at;
use crate::zone::Zone;
use crate::{
    Address, Error, Extent, GeoPoint, GridConfig, Indexing, PlanarPoint, Projection, Result,
    Topology, ZoneId,
};

/// Ceiling on how many sub-zones [`Grid::sub_zones`] will materialise into a
/// `Vec`, carried over from py4dggs's `MAX_MATERIALISED_SUB_ZONES`: the greatest
/// number this crate will ever list, whatever the environment. The limit in force,
/// which the environment variable `RS4DGGS_MAX_MATERIALISED_SUB_ZONES` can lower
/// and never raise, is [`max_materialised_sub_zones`].
///
/// `count_sub_zones` is a closed form and stays cheap at any depth, but
/// `sub_zones` builds the whole sequence, and the depth was validated only
/// against the maximum resolution, never against the cardinality of the
/// result. An aperture-3 cell at resolution 0 counts some 4.6e15 sub-zones at
/// depth 33 and passes every other guard, which turned a caller-supplied depth
/// into a one-line denial of service. Four million keeps every realistic
/// tiling request working while refusing the unservable ones; a caller who
/// wants the number without the cost still calls `count_sub_zones`.
pub const MAX_MATERIALISED_SUB_ZONES: u64 = 4_000_000;

/// The name of the environment variable that lowers the limit on a list of sub-zones.
const MAX_MATERIALISED_SUB_ZONES_VARIABLE: &str = "RS4DGGS_MAX_MATERIALISED_SUB_ZONES";

/// The limit in force on a list of sub-zones: [`MAX_MATERIALISED_SUB_ZONES`], unless the
/// environment variable `RS4DGGS_MAX_MATERIALISED_SUB_ZONES` lowers it.
///
/// The variable is read once, at the first use, and the answer is kept for the life of the
/// process. Its value is read as Rust reads a `u64` (`u64::from_str`, so a leading `+` and leading
/// zeros are taken) once surrounding whitespace is trimmed; a number above
/// [`MAX_MATERIALISED_SUB_ZONES`] gives that ceiling, so the variable can only lower the limit.
/// A value that is empty, not a number, negative or too large for a `u64` is ignored in
/// silence, as a library prints nothing, and the limit is the ceiling. On a target without an
/// environment (`wasm32-unknown-unknown`, for one) the variable is absent and the limit is the
/// ceiling.
///
/// The limit bounds the list that [`Grid::sub_zones`] builds, and with it the search of
/// [`Grid::sub_zone_index`], which builds that list or walks the order to the zone sought: a
/// request whose order is longer is refused with [`Error::TooManySubZones`], which carries
/// the limit in force. [`Grid::count_sub_zones`], [`Grid::first_sub_zone`] and
/// [`Grid::sub_zone_at_index`] are not bounded by it, and a depth of 0, whose list is the zone
/// itself, is answered whatever the limit, a limit of 0 included.
pub fn max_materialised_sub_zones() -> u64 {
    static LIMIT: OnceLock<u64> = OnceLock::new();
    *LIMIT.get_or_init(|| {
        materialised_sub_zones_limit(
            std::env::var(MAX_MATERIALISED_SUB_ZONES_VARIABLE)
                .ok()
                .as_deref(),
        )
    })
}

/// The reading of the variable, apart from the environment: see
/// [`max_materialised_sub_zones`].
fn materialised_sub_zones_limit(value: Option<&str>) -> u64 {
    value
        .and_then(|v| v.trim().parse::<u64>().ok())
        .map_or(MAX_MATERIALISED_SUB_ZONES, |n| {
            n.min(MAX_MATERIALISED_SUB_ZONES)
        })
}

/// The geometry of a zone that has none: the zero point as its centroid, and six
/// zero vertices. It is DGGAL's own on the aperture-7 grids, where
/// `getZoneWGS84Centroid` and `getZoneWGS84Vertices` on a `nullZone` give the zero
/// point and six zero vertices, six regardless of whether the address sits on a
/// pentagon path. On the aperture-3 grids this crate gives the same, a deliberate
/// departure: there the engine gives its null zone a centroid that belongs to no
/// cell and six vertices about it.
const NULL_GEOPOINT: GeoPoint = GeoPoint { lat: 0.0, lon: 0.0 };

/// How many vertices a zone with no geometry reports: six on every grid. On the
/// three aperture-7 grids that is DGGAL's own count, on the pentagon path as well
/// as the hexagon one, verified against DGGAL for each. On the three aperture-3
/// grids the engine gives its null zone six vertices too, as the oracle suites
/// re-check on each, though not at the zero point: there the count is DGGAL's and
/// the zero vertices are this crate's departure.
const NULL_VERTEX_COUNT: usize = 6;

/// The ceiling of the edge refinement [`Grid::refined_vertices`] accepts: the greatest one this
/// crate will ever take, whatever the environment. A greater one is refused with
/// [`Error::EdgeRefinementOutOfRange`]. The limit in force, which the environment variable
/// `RS4DGGS_MAX_EDGE_REFINEMENT` can lower and never raise, is [`max_edge_refinement`].
///
/// DGGAL sets no bound. A refinement divides every edge of a zone into that many parts, so a
/// hexagon at 100,000 has 600,000 points, and one with an edge beside a pole, which the engine
/// divides twenty times as finely, 2,500,000, each of them an inverse projection. That is the
/// greatest refinement at which this crate's boundary was compared with the engine's; beyond
/// some 107 million the engine's own count of divisions, which it multiplies by twenty beside
/// a pole, wraps.
pub const MAX_EDGE_REFINEMENT: u32 = 100_000;

/// The name of the environment variable that lowers the edge-refinement limit.
const MAX_EDGE_REFINEMENT_VARIABLE: &str = "RS4DGGS_MAX_EDGE_REFINEMENT";

/// The edge-refinement limit in force: [`MAX_EDGE_REFINEMENT`], unless the environment variable
/// `RS4DGGS_MAX_EDGE_REFINEMENT` lowers it.
///
/// The variable is read once, at the first use, and the answer is kept for the life of the
/// process. Its value is read as Rust reads a `u32` (`u32::from_str`, so a leading `+` and leading
/// zeros are taken) once surrounding whitespace is trimmed; a number above
/// [`MAX_EDGE_REFINEMENT`] gives that ceiling, so the variable can only lower the limit. A
/// value that is empty, not a number, negative or too large for a `u32` is ignored in silence,
/// as a library prints nothing, and the limit is the ceiling. On a target without an
/// environment (`wasm32-unknown-unknown`, for one) the variable is absent and the limit is the ceiling.
///
/// The limit bounds the refinement a caller asks of [`Grid::refined_vertices`]. A refinement of
/// 0, the engine's own choice, is always accepted, whatever the limit: that choice is
/// at most twenty points to an edge on an ordinary edge, though an edge beside a pole's image
/// takes twenty times as many, and [`Grid::extent`] uses it. A limit below 20 therefore does
/// not bound the automatic ring.
pub fn max_edge_refinement() -> u32 {
    static LIMIT: OnceLock<u32> = OnceLock::new();
    *LIMIT.get_or_init(|| {
        edge_refinement_limit(std::env::var(MAX_EDGE_REFINEMENT_VARIABLE).ok().as_deref())
    })
}

/// The reading of the variable, apart from the environment: see [`max_edge_refinement`].
fn edge_refinement_limit(value: Option<&str>) -> u32 {
    value
        .and_then(|v| v.trim().parse::<u32>().ok())
        .map_or(MAX_EDGE_REFINEMENT, |n| n.min(MAX_EDGE_REFINEMENT))
}

/// The eC's `Degrees { 0.05 }` of `getRefinedVertices` (`RI7H.ec:573`, `RI3H.ec:517`), in
/// radians: the shift of the longitude about which the engine unwraps the longitudes of a
/// refined boundary. It is pinned as bits because neither of its source forms gives it: the
/// eC compiler folded the angle through a literal of sixteen digits, `0.0008726646259972`,
/// which is the double the library holds (`.rodata` at `0x49fe0`), where `0.05 * Pi / 180`
/// is `0x3f4c987103b761f5`.
const REFINED_LON_SHIFT: f64 = f64::from_bits(0x3f4c_9871_03b7_633a);

/// `fl(Pi + 1e-9)`, beyond which `getRefinedVertices` brings the longitude of the centroid
/// back by a turn (`RI7H.ec:559-563`, `RI3H.ec:481-485`): the library holds it folded at
/// `0x49fd0`, and its negation, the folded `-Pi - 1E-9`, at `0x49fc8`.
const CENTROID_LON_LIMIT: f64 = math::PI + 1e-9;

/// What `GeoExtent::clear` writes (`GeoExtent.ec:102-106`, compiled at `0x42840`): the eC's
/// `MAXDOUBLE` taken as degrees and folded into radians (`.rodata` at `0x4a4c8`, and its
/// negation at `0x4a4d0`), which is not `f64::MAX`. The lower left corner of a cleared extent
/// holds it in both members and the upper right its negation, so that the first point of a
/// ring replaces all four.
const CLEARED_EXTENT: f64 = f64::from_bits(0x7f91_df46_a252_9d38);

/// The latitude, in degrees, beyond which `RI3H::getRefinedVertices` takes a point of a refined
/// boundary to lie on a pole, north or south (`RI3H.ec:529`): the library holds it at `0x49ff8`.
const POLE_LATITUDE: f64 = 89.999999;

/// The distance in longitude, in degrees, from the centroid of the zone beyond which
/// `RI3H::getRefinedVertices` turns a pole of a refined boundary by half a turn
/// (`RI3H.ec:601-604`, `:624-627`): the library holds it at `0x4a008`, and its negation at
/// `0x4a000`. An ordinary point is turned beyond 120 degrees (`RI3H.ec:522-525`). The turn
/// is taken at zones of every sample compared with the engine, but at none of some 87,000
/// zones a grid does a pole lie between 95 and 120 degrees from the centroid: the bound is
/// held by this constant against the library's, and by no ring.
const POLE_HALF_TURN_BEYOND: f64 = 95.0;

/// The number of parts into which the engine divides each edge of a zone of `level` where it
/// is asked for the automatic refinement, an `edgeRefinement` of 0: the eC's `nDivisions`
/// (`RI7H.ec:542-543`, `RI3H.ec:505-506`), by the aperture, which is all that separates the
/// two classes here.
fn automatic_edge_divisions(aperture: u8, level: u8) -> i32 {
    // faithful: the library takes the chain of conditionals as written, on the level of the
    // zone: on aperture 7 at `0x2a9c0` to `0x2a9de` and `0x2ac90` to `0x2ac98`, where the eC's
    // last two arms, both 12, are one; on aperture 3 at `0xf6e2` to `0xf6fa` and `0x102f0` to
    // `0x10309`, where its last two, both 5, are one.
    match (aperture, level) {
        (_, 0..=2) => 20,
        (_, 3..=4) => 15,
        (7, _) => 12,
        (_, 5..=7) => 10,
        (_, 8..=9) => 8,
        _ => 5,
    }
}

/// The `oddGrid` that `getRefinedVertices` passes to the inverse projection of every point of
/// a refined boundary, and under which alone it applies its correction of half a turn. On
/// aperture 7 it is true whatever the level (`RI7H.ec:552`, where the eC's comment says that
/// it mends a south pole otherwise flipped; the library loads the constant at `0x2ac26`); on
/// aperture 3 it is `zone.subHex > 0`, that is an odd level (`RI3H.ec:473`, read from the
/// stack at `0xfc85`). [`Grid::vertices`] and [`Grid::centroid`] pass the parity of the level
/// on both.
///
/// On aperture 3 it follows the eC and the library, and no ring of the engine shows it: with
/// the flag always true, or never, the ring of each of some 87,000 zones a grid is the same.
/// The inverse projection reads the flag only for a point on a pole's image, which the pole
/// branch then replaces, and no other point of those rings lies beyond 120 degrees of its
/// centroid.
fn refined_odd_grid(aperture: u8, level: u8) -> bool {
    aperture == 7 || level % 2 == 1
}

/// One discrete global grid reference system: a projection, a topology and an
/// indexing, plus the orientation they share.
///
/// The three parts are type parameters rather than fields, so a grid costs
/// nothing at run time beyond its precomputed projection geometry, and a new
/// grid is a new combination of existing parts rather than a change to any of
/// them.
pub struct Grid<P: Projection, T: Topology, I: Indexing> {
    geom: P::Geometry,
    config: GridConfig,
    name: &'static str,
    _m: PhantomData<fn() -> (T, I)>,
}

impl<P: Projection, T: Topology, I: Indexing> std::fmt::Debug for Grid<P, T, I> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Grid").field("name", &self.name).finish()
    }
}

/// `Ok` if `value` is finite, and otherwise [`Error::NonFinite`] naming it.
fn finite(name: &'static str, value: f64) -> Result<()> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(Error::NonFinite { name })
    }
}

impl<P: Projection, T: Topology, I: Indexing> Grid<P, T, I> {
    /// Builds the grid, precomputing the projection's geometry once.
    ///
    /// Refuses a configuration whose `orientation_lat_deg`, `orientation_lon_deg` or
    /// `azimuth_deg` is NaN or infinite, with [`Error::NonFinite`] naming the first
    /// such field in that order. A grid built on one would place every point at
    /// NaN, and so answer each question with the same fixed cell or with NaN
    /// coordinates, rather than fail. This check has no counterpart in the engine,
    /// which fixes each grid's orientation in its own source and takes none from
    /// a caller; it is this crate's own, as the configuration is.
    ///
    /// Only finiteness is checked. A finite orientation outside the usual ranges
    /// of latitude and longitude is taken as it stands, as the angle it denotes.
    pub fn new(config: GridConfig, name: &'static str) -> Result<Self> {
        finite("orientation_lat_deg", config.orientation_lat_deg)?;
        finite("orientation_lon_deg", config.orientation_lon_deg)?;
        finite("azimuth_deg", config.azimuth_deg)?;
        Ok(Grid {
            geom: P::build_geometry(&config),
            config,
            name,
            _m: PhantomData,
        })
    }

    /// The grid's DGGRS name, as the registry knows it.
    pub fn name(&self) -> &'static str {
        self.name
    }

    /// The orientation and datum handling this grid was built with.
    pub fn config(&self) -> &GridConfig {
        &self.config
    }

    /// The finest resolution this grid's indexing can pack.
    pub fn max_resolution(&self) -> u8 {
        I::MAX_RESOLUTION
    }

    /// Decodes `id` and refuses an address this grid has no cell for.
    ///
    /// A [`ZoneId`] is a bare 64-bit integer, and [`Grid::zone`] wraps any of
    /// them without asking questions, so an identifier this crate never
    /// produced, whether fabricated, corrupted in transit or read from a
    /// stale store, reaches every operation below. Nothing in the packing
    /// stops it: `Z7::decode` reads the base cell out of the top four bits, so
    /// such an identifier decodes to a base cell of 12 to 15, and the topology
    /// then indexes its own twelve-entry tables with it. This crate's contract,
    /// that caller input to [`Grid`], [`crate::Zone`] and [`crate::AnyGrid`]
    /// returns an error rather than panicking, is why every operation that
    /// reaches the topology decodes through here first. [`Projection`],
    /// [`Topology`] and [`Indexing`], and the implementations of them this
    /// crate ships, are its own internals made visible so that a caller can
    /// compose a new grid from them; the contract binds `Grid`, `Zone` and
    /// `AnyGrid`, the three types a caller actually calls, not a trait or an
    /// implementation of one taken on its own.
    ///
    /// Which addresses have a cell behind them is the topology's own
    /// knowledge, not this type's: the aperture-3 indexing packs an entire
    /// cell into the base field, where a bound of eleven would be meaningless,
    /// so the question goes to [`Topology::is_valid_address`].
    ///
    /// Note that this is a different question from
    /// [`Topology::is_null_geometry`]. There, the address is one the topology
    /// recognises and DGGAL still refuses to give geometry for, Z7's twentieth
    /// packing level being the case in point. Here, there is no such address
    /// at all. The two give a caller the same answers, for different reasons.
    fn address(&self, id: ZoneId) -> Result<Address> {
        let a = I::decode(id);
        if T::is_valid_address(&a) {
            Ok(a)
        } else {
            Err(Error::InvalidZone(format!(
                "zone id {:#018x} decodes to an address {} has no cell for",
                id.0, self.name
            )))
        }
    }

    // --- indexing entry points ---

    /// The zone containing the given geographic point at resolution `res`, or
    /// the null zone where the engine says the point falls in no cell. The point
    /// is a latitude and a longitude in decimal degrees.
    ///
    /// A null answer is a successful call, not an error: [`ZoneId::NULL`] comes
    /// back inside `Ok`, exactly as DGGAL returns its `nullZone` from
    /// `getZoneFromWGS84Centroid` rather than failing. It is reachable. On the
    /// three aperture-7 grids it arises in a narrow strip along two of the
    /// icosahedron's edges, one of which passes over the north pole, though the
    /// strip itself is not polar: at the odd resolutions from 15 the engine
    /// answers the null zone there as well, and at resolutions 16 and 18, and at a
    /// few points of the strip's margin at 19, the engine answers an identifier
    /// that its own parser cannot read back, where this crate answers the null
    /// zone by design. On the three aperture-3 grids the engine answers it near
    /// both poles and beside some icosahedron edges at the coarsest resolutions,
    /// and along two seams of the layout, at resolutions 2 to 21, this crate
    /// answers it where the engine answers an identifier it cannot read back.
    /// [`ZoneId::NULL`] describes where, and [`Topology::quantize`] and
    /// [`crate::topologies::HexA3`] why. No sample the oracle suites take
    /// contains a point at which this crate answers the null zone while the
    /// engine names a cell it can read back. The `Err` arm is reserved for what the caller got
    /// wrong: a resolution outside the indexing's range, which is
    /// [`Error::ResolutionOutOfRange`], or a coordinate that is not finite, which
    /// is [`Error::NonFinite`] naming `"latitude"` or `"longitude"`. The
    /// resolution is checked first, then the latitude, then the longitude, and the
    /// first at fault is the one reported.
    ///
    /// The bound on `res` is still load-bearing rather than defensive: it is
    /// what keeps `quantize` from being asked about a resolution the packing
    /// cannot represent, so it must not be weakened.
    ///
    /// The refusal of a non-finite coordinate is a deliberate departure from
    /// DGGAL. Given NaN or an infinity for either coordinate, the engine answers
    /// one fixed cell whatever the input, a pentagon at the requested resolution
    /// centred at about (58.4, -168.8) in the northern Pacific, on all six grids:
    /// that of base cell 1 on the aperture-7 grids, and that of root rhombus 0 on
    /// the aperture-3 grids. A cell is an answer about a place, and a service that
    /// passed such a coordinate through would report that pentagon as the place it
    /// named, so the input is refused instead.
    ///
    /// A finite coordinate is not otherwise validated. A latitude outside -90 to
    /// 90 or a longitude outside -180 to 180 is projected as it stands, as the
    /// engine projects it, and yields some zone, or the null zone, rather than an
    /// error; a caller who needs such input refused checks it first.
    ///
    /// # Examples
    ///
    /// ```
    /// let grid = rs4dggs::igeo7();
    /// let zone = grid.zone_from_geo(38.7223, -9.1393, 5)?; // Lisbon
    /// assert_eq!(zone.text_id(), "0064156");
    /// assert_eq!(zone.resolution(), 5);
    ///
    /// // A resolution above the grid's maximum, or a coordinate that is not finite, is refused.
    /// assert!(grid.zone_from_geo(38.7223, -9.1393, 20).is_err());
    /// assert!(grid.zone_from_geo(f64::NAN, -9.1393, 5).is_err());
    /// # Ok::<(), rs4dggs::Error>(())
    /// ```
    pub fn zone_from_geo(&self, lat: f64, lon: f64, res: u8) -> Result<Zone<'_, Self>> {
        if res > I::MAX_RESOLUTION {
            return Err(Error::ResolutionOutOfRange {
                res,
                max: I::MAX_RESOLUTION,
            });
        }
        finite("latitude", lat)?;
        finite("longitude", lon)?;
        let p = P::forward(&self.geom, lat, lon);
        let id = match T::quantize(p, res) {
            Some(a) => I::encode(&a),
            None => ZoneId::NULL,
        };
        Ok(self.zone(id))
    }

    /// The zone named by its canonical text id, which the indexing validates.
    ///
    /// Text is refused with [`Error::InvalidZone`] when it is malformed, when it
    /// is not canonical, as a digit path through a pentagon's deleted child is not,
    /// and when it names a resolution finer than [`Grid::max_resolution`]. The
    /// error quotes the text at bounded length: at most its first thirty-two
    /// characters, and then its length.
    ///
    /// Two of these are deliberate departures from DGGAL, which answers its null
    /// zone for any text it cannot parse where this crate names the failure, and
    /// which reads the twenty-digit text of a Z7 level-20 identifier, a level its
    /// packing can hold but it gives no geometry: the zero point as its centroid,
    /// six zero vertices and no neighbours. A caller could not tell that answer
    /// from a genuine cell at (0, 0), so this crate refuses the text instead.
    /// [`Grid::text_id`] still prints such an identifier, if one is bound with
    /// [`Grid::zone`], and that text does not read back. On the aperture-3 grids
    /// the parser departs from DGGAL a third time, as [`crate::indexings::I3h`]
    /// describes: the text of a polar pentagon's sub-hexagon C or D, such as
    /// `AA-0-C`, is refused, where the engine reads it although no relation of any
    /// zone reaches such a cell.
    ///
    /// The text, like the identifier it spells, carries no grid. The three
    /// aperture-7 grids share one text form, and so do the three aperture-3 grids,
    /// so the text of a zone of one grid is read on the others of its family without
    /// complaint, and there it names another cell, commonly several kilometres away
    /// on the aperture-7 grids; see [`Grid::zone`].
    ///
    /// # Examples
    ///
    /// ```
    /// let grid = rs4dggs::igeo7();
    /// let zone = grid.zone_from_text("0064156")?;
    /// assert_eq!(zone, grid.zone_from_geo(38.7223, -9.1393, 5)?);
    /// assert_eq!(zone.text_id(), "0064156");
    ///
    /// // A digit outside 0 to 6 names no zone.
    /// assert!(grid.zone_from_text("0064158").is_err());
    /// # Ok::<(), rs4dggs::Error>(())
    /// ```
    pub fn zone_from_text(&self, text: &str) -> Result<Zone<'_, Self>> {
        Ok(self.zone(I::from_text(text)?))
    }

    /// Binds an identifier to this grid without validating it.
    ///
    /// This is the cheap way back from a stored identifier to a usable handle,
    /// so it stays free of checks. An identifier this crate never produced is
    /// not rejected here: it answers as DGGAL's null zone on the geometry and
    /// the hierarchy, namely the zero point, six zero vertices and no
    /// neighbours, parents or children, and it is reported as
    /// [`Error::InvalidZone`] by the sub-zone methods, which can say so.
    ///
    /// Nor is the text of such an identifier guaranteed to lead back to it.
    /// [`Grid::text_id`] prints what the bits spell, and that can be text which
    /// [`Grid::zone_from_text`] rejects, such as "002", whose digit 2 is the
    /// deleted child of a pentagon, or text that parses to a different
    /// identifier, such as "03" for an identifier whose unused digit slots do
    /// not all hold the terminator, or text longer than the finest resolution
    /// allows, such as the twenty-two characters of a Z7 level-20 identifier.
    ///
    /// A level-20 identifier bound here is one the topology recognises, and it
    /// keeps its resolution and its text; it answers the engine's null geometry,
    /// the zero point as its centroid, six zero vertices and no neighbours, parents
    /// or children, exactly as DGGAL answers it. The refusal of such a
    /// zone's text by [`Grid::zone_from_text`] therefore closes the text route to
    /// it only, and a caller that binds raw identifiers from outside should compare
    /// their resolution with [`Grid::max_resolution`].
    ///
    /// Nor are `resolution` and `is_pentagon` guarded: for an unvalidated
    /// identifier they simply read the bit fields as they stand, so an
    /// identifier with base cell 12, which names no real base cell, still
    /// reports a pentagon, which is also what the engine answers for it.
    ///
    /// An identifier carries no grid. The three aperture-7 grids share the Z7
    /// packing, so an identifier taken from one of them is bound here without
    /// complaint, and names another cell. Measured over 60,000 identifiers drawn
    /// from all three grids, the centroid of the cell it names on a second grid
    /// lies on average some 7 to 14 km from its centroid on the first, according
    /// to the pair of grids (IVEA7H and RTEA7H being the closest, IGEO7 and RTEA7H
    /// the farthest), and up to 56 km. That holds at every resolution from 1, the
    /// average barely changing from resolution 3 onwards, so at the finer
    /// resolutions the cell named lies many thousands of cells from the one
    /// intended. The pentagons, the base cells of resolution 0 among them, are
    /// centred on the same twelve points on every grid, and they do coincide.
    /// The three aperture-3 grids likewise share the I3H packing: the zone that
    /// contains Lisbon at resolution 5 is `C2-23-C` on each of them, but at
    /// resolution 15 each names a different zone for it. Nothing in the identifier
    /// reveals the mistake, and a caller that keeps identifiers of more than one
    /// grid must keep each grid's name beside them.
    pub fn zone(&self, id: ZoneId) -> Zone<'_, Self> {
        Zone::new(self, id)
    }

    /// The identifier's canonical text form.
    pub fn text_id(&self, id: ZoneId) -> String {
        I::to_text(id)
    }

    /// The resolution encoded in the identifier.
    ///
    /// One deliberate divergence, at the null zone, whose resolution reads as 0 on
    /// every grid. On the aperture-7 grids DGGAL reports a level of -1 for its
    /// sentinel, special-cased before any field is decoded (`Z7Zone::level`,
    /// RI7H_Z7.ec:37-38), whereas this reports 0, because the sentinel's digit slots
    /// all hold the terminator. On the aperture-3 grids DGGAL reads the sentinel's
    /// fields as they stand and reports 63 (`getZoneLevel`, RI3H.ec:53-56, and the
    /// `level` property, :878-881), whereas this reports 0 there too, so that the
    /// null zone reads alike on every grid. The return type stays unsigned rather
    /// than widening to carry a value no caller can act on: a caller that has to
    /// tell a resolution-0 zone from the null zone compares the identifier against
    /// [`ZoneId::NULL`], or its text against [`crate::NULL_TEXT`]. [`ZoneId::NULL`]
    /// carries the same note.
    pub fn resolution(&self, id: ZoneId) -> u8 {
        I::resolution(id)
    }

    /// Whether the identifier names a pentagon, that is an icosahedron vertex
    /// cell.
    pub fn is_pentagon(&self, id: ZoneId) -> bool {
        I::is_pentagon(id)
    }

    // --- geometry: planar (topology) -> geographic (projection) ---

    /// Whether the zone belongs to an odd grid, the `oddGrid` argument the engine
    /// passes to the inverse projection, which uses it only to choose the longitude of
    /// a point it snaps onto a pole.
    ///
    /// The aperture-7 engine asks two different questions at the two call sites this
    /// port mirrors: the WGS84 centroid passes `zone.subHex > 0` (`RI7H.ec:349`) and
    /// the WGS84 vertices pass `zone.level & 1` (`RI7H.ec:458`). They are the same
    /// question, because an `I7HZone`'s level is `2 * levelI49R + (subHex > 0)`
    /// (`RI7H.ec:891`), so `subHex > 0` exactly when the level is odd. A Z7 zone
    /// reaches those call sites through `to7H`, which walks one 7H child per Z7 digit
    /// (`RI7H_Z7.ec:298-345`), so the 7H level is the Z7 resolution, and both reduce to
    /// an odd resolution. A third call site, the centroid in a WGS84 CRS, passes
    /// `false` whatever the zone (`RI7H.ec:338`); this port has no counterpart to it and
    /// follows the WGS84 centroid.
    ///
    /// The same rule serves the aperture-3 grids. They pass `zone.subHex > 0`
    /// at both of theirs (`RI3H.ec:268` and `:377`), and an `I3HZone`'s level is
    /// likewise `levelI9R * 2 + (subHex > 0)` (`RI3H.ec:880`), with `subHex` zero
    /// exactly at even levels (`RI3H.ec:861`).
    fn odd_grid(&self, id: ZoneId) -> bool {
        I::resolution(id) % 2 == 1
    }

    /// The zone's centroid in WGS84 degrees.
    pub fn centroid(&self, id: ZoneId) -> GeoPoint {
        let Ok(a) = self.address(id) else {
            return NULL_GEOPOINT;
        };
        if T::is_null_geometry(&a) {
            return NULL_GEOPOINT;
        }
        P::inverse(&self.geom, T::planar_centroid(&a), self.odd_grid(id))
    }

    /// The zone's boundary vertices in WGS84 degrees, anticlockwise as seen from
    /// outside the sphere: six for a hexagon, five for a pentagon. No closing vertex is
    /// appended, though the engine's own ring can end on the point it began at, as
    /// set out below. A zone with no geometry, the null zone or a Z7 level-20 zone,
    /// gives six zero points, as the engine does; so does an identifier with no
    /// cell behind it.
    ///
    /// Every ring that encloses any area runs anticlockwise, on every grid, save the
    /// few aperture-7 rings that hold the point (0, 0), set out below. On the
    /// aperture-7 grids that is the engine's own order. It was measured on some
    /// 35,700 zones of each of the three, every zone to resolution 3 and, at every
    /// resolution, the pentagons, the cells at both poles, cells on the antimeridian
    /// and a uniform sample, and the oracle suites find every ring they compare in the
    /// engine's order and anticlockwise, save those set out below. On the aperture-3
    /// grids the engine's ring runs clockwise, save at the two polar pentagons, and
    /// this crate departs from its order on purpose: the ring starts at the engine's
    /// first vertex and takes the others in reverse, `[v0, v(n-1), ..., v1]` for the
    /// engine's `[v0, v1, ..., v(n-1)]`, while the polar pentagons, whose rings already
    /// run anticlockwise, keep the engine's order. Either way every coordinate is the
    /// engine's own, and the oracle suites compare each with the engine's, bit for
    /// bit, in that order.
    ///
    /// Each vertex's longitude is given on its own, between -180 and 180, as the
    /// inverse projection yields it, and the ring is not unwrapped: where a cell
    /// straddles the antimeridian, consecutive vertices differ in longitude by
    /// almost 360 degrees, and a ring about a pole runs through every longitude.
    /// A caller that draws or intersects rings on a plane must split or shift such
    /// a ring itself.
    ///
    /// A ring can also name one point more than once. Near a pole at resolutions
    /// 18 and 19 of the aperture-7 grids, the engine's pole snap, which this crate
    /// reproduces, moves some vertices onto the pole, so that consecutive vertices
    /// coincide, either exactly or on the pole at different longitudes, and at
    /// resolution 19 some cells have every vertex on the pole and so enclose no
    /// area at all. At resolutions 16 and 18, beside the pentagons of base cells 0
    /// and 1, a ring of six can end on the point it began at, its sixth vertex
    /// repeating its first exactly, so that a caller testing the first and last
    /// vertices for equality would take it for a closed pentagon. Under the pentagon
    /// of base cell 0, at resolutions 15, 17 and 19, the engine puts the point (0, 0)
    /// among the vertices of some cells, one or three of them: a planar vertex there
    /// lies on the diagonal of the 5x6 layout, which the inverse projection answers
    /// with the zero point. The rest of such a ring runs anticlockwise. Each of these
    /// is the engine's own answer.
    pub fn vertices(&self, id: ZoneId) -> Vec<GeoPoint> {
        let Ok(a) = self.address(id) else {
            return vec![NULL_GEOPOINT; NULL_VERTEX_COUNT];
        };
        if T::is_null_geometry(&a) {
            return vec![NULL_GEOPOINT; NULL_VERTEX_COUNT];
        }
        let odd = self.odd_grid(id);
        // The engine's `getZoneWGS84Vertices` passes every planar vertex through
        // `canonicalize5x6` (`ri5x6.ec:1562-1606`) before the inverse projection
        // (`RI7H.ec:461-462`), folding a vertex on a seam of the 5x6 layout, or beyond
        // one of its edges, back into the canonical domain. This crate omits the call, as
        // py4dggs does for the aperture-7 grids, and the omission changes no answer on
        // them: the engine's own `canonicalize5x6`, applied to each of the 10,917,150
        // planar vertices of 1,819,645 zones on the three grids (every zone to resolution
        // 4, and at every resolution from 1 to 19 the zones that hold 13,386 points on and
        // beside the thirty edges of the icosahedron and near both poles, with their
        // neighbours), moved none of them, and every centroid and vertex of those zones
        // was the engine's to the last bit. On the aperture-3 grids it does move vertices,
        // on cells whose boundary crosses an interruption, and `HexA3::planar_vertices`
        // performs it (`RI3H.ec:380-381`), so the vertices reach this line in the engine's form.
        T::planar_vertices(&a)
            .into_iter()
            .map(|p| P::inverse(&self.geom, p, odd))
            .collect()
    }

    /// The zone's area in square metres, as DGGAL's `getZoneArea` computes it
    /// (`RI7H_Z7.ec:514-525`; `RI3H.ec:66-88`), or `None` for a zone without geometry.
    ///
    /// It is a nominal area and not one measured from the boundary: the engine takes every
    /// hexagon of a resolution as the same share of the area of the whole ellipsoid, that area
    /// divided by the count of hexagons, and a pentagon as five sixths of one. It is so on all
    /// six grids, whatever the projection, and the doubles are the engine's to the last bit.
    ///
    /// A zone without geometry has no area. That is the null zone, a Z7 identifier of twenty
    /// digits, and an identifier that no cell of this grid has behind it, the very ones for which
    /// [`Grid::vertices`] gives six zero points. The engine answers `+inf` for its null zone, and
    /// for a twenty-digit identifier the area of that level's formula, a finite number for a zone
    /// it cannot draw; this crate answers neither. An aperture-7 identifier that has an address
    /// but does not read back as text keeps its ordinary area, as it keeps its ordinary vertices.
    ///
    /// # Examples
    ///
    /// ```
    /// let grid = rs4dggs::igeo7();
    /// let hexagon = grid.zone_from_text("0064156")?;
    /// assert!((grid.area(hexagon.id()).unwrap() - 3.03484e9).abs() < 1e4);
    ///
    /// // A pentagon is five sixths of a hexagon's share.
    /// let pentagon = grid.zone_from_text("00")?;
    /// assert!((pentagon.area().unwrap() - 4.250547e13).abs() < 1e8);
    ///
    /// // The null zone has none.
    /// assert_eq!(grid.area(rs4dggs::ZoneId::NULL), None);
    /// # Ok::<(), rs4dggs::Error>(())
    /// ```
    pub fn area(&self, id: ZoneId) -> Option<f64> {
        self.geometry_address(id)?;
        Some(facts::zone_area(
            T::APERTURE,
            I::resolution(id),
            I::is_pentagon(id),
        ))
    }

    /// The address of a zone that has geometry, or `None` for one that has none: an identifier
    /// with no cell of this grid behind it, the null zone among them, or an address whose
    /// geometry the topology declares null, which is a Z7 identifier of twenty digits. They are
    /// the very zones for which [`Grid::vertices`] gives its six zero points, by the same two
    /// tests.
    ///
    /// No eC counterpart: the engine answers each such zone in its own way, an empty ring, a
    /// cleared extent and an area all the same, and this crate answers them alike.
    pub(crate) fn geometry_address(&self, id: ZoneId) -> Option<Address> {
        let a = self.address(id).ok()?;
        (!T::is_null_geometry(&a)).then_some(a)
    }

    /// The zone's centroid as the engine's `getZoneWGS84Centroid` answers it (`RI7H.ec:347-350`,
    /// `RI3H.ec:266-269`): the latitude and the longitude in radians, and the zero point where
    /// the inverse projection fails, since the eC's `inverse` clears its result before it
    /// returns false (`ri5x6.ec:554`). [`Grid::centroid`] is this point in degrees.
    fn centroid_radians(&self, id: ZoneId, a: &Address) -> (f64, f64) {
        P::inverse_radians(&self.geom, T::planar_centroid(a), self.odd_grid(id))
            .unwrap_or((0.0, 0.0))
    }

    /// The refined boundary of a zone with geometry as the engine's `getRefinedVertices`
    /// builds it for WGS84 (`RI7H.ec:536-612`, compiled at `0x2a960`; `RI3H.ec:451-673`, at
    /// `0xf580`): each point a latitude and a longitude in radians, in that order.
    /// `edge_refinement` is the eC's `edgeRefinement`, 0 for the engine's own choice. The ring
    /// is empty where the topology has no refined boundary.
    pub(crate) fn refined_ring_radians(
        &self,
        id: ZoneId,
        a: &Address,
        edge_refinement: i32,
    ) -> Vec<(f64, f64)> {
        // faithful: DGGAL v0.0.6 as built with `-O2 -ffast-math` (`Makefile.dggal:133`),
        // BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`. `getRefinedVertices` carries no
        // attribute, so the library decides its order, and this function follows the library
        // from `0x2a960` to `0x2ae4c` and, for aperture 3, from `0xf580` to `0x1031f`. One sum
        // is re-associated in either, site GR-1 below; the rest is performed as the eC writes
        // it.
        let level = I::resolution(id);
        let n_divisions = if edge_refinement != 0 {
            edge_refinement
        } else {
            automatic_edge_divisions(T::APERTURE, level)
        };
        let Some(ring) = T::planar_refined_vertices(a, n_divisions) else {
            return Vec::new();
        };
        let odd_grid = refined_odd_grid(T::APERTURE, level);

        // faithful: `RI7H.ec:557-563`, at `0x2aa59` to `0x2aa79` and `0x2ae30` to `0x2ae46`:
        // the centroid's longitude, brought back by a turn where it lies beyond
        // `Pi + 1e-9` of zero. The eC's two tests are compiled as one chain, which is the
        // same, since a longitude that the first has raised cannot pass the second. Neither
        // is ever entered on the grids of this crate: the inverse projection wraps its
        // longitude to within `Pi + 10 * DBL_EPSILON` of zero, which is nearer.
        let mut c_lon = self.centroid_radians(id, a).1;
        if c_lon < -CENTROID_LON_LIMIT {
            c_lon += math::TWO_PI;
        } else if c_lon > CENTROID_LON_LIMIT {
            c_lon -= math::TWO_PI;
        }

        // The eC's two tests of nearness, `RI7H.ec:587-589` and `:591-593`: both members
        // within 1e-11, in radians, on the pair as stored (`0x2ab75` to `0x2abb7`, `0x2ad78`
        // to `0x2adb7`).
        let near =
            |p: (f64, f64), q: (f64, f64)| (p.0 - q.0).abs() < 1e-11 && (p.1 - q.1).abs() < 1e-11;
        // The first of the two, which drops a point that repeats the one before it. A ring of
        // one point skips it (`0x2ab61`).
        let drop_a_repeat = |points: &mut Vec<(f64, f64)>| {
            let count = points.len();
            if count >= 2 && near(points[count - 1], points[count - 2]) {
                points.pop();
            }
        };
        // faithful: `RI3H.ec:474`, at `0xf647` to `0xf677`, and `:553` with `:592-593`, hoisted
        // to `0xf733` to `0xf75c`: the eC's `poleOffset`, a thousandth of two to the power of
        // the I9R level, and the distance by which the pole branch below moves a planar point
        // to either side, `(val * poleOffset) * poleOffset` with `val` 0.00001. The library
        // forms `poleOffset * -1e-5` once, turns its sign, and multiplies either by
        // `poleOffset` where an arm asks it (`0xf85b`, `0xfeca`, `0xfece`, `0x100be`,
        // `0x101c8`, `0x101cc`): the same double as this, to the sign. Aperture 7 has no such
        // branch.
        let pole_nudge = if T::APERTURE == 3 {
            let pole_offset = f64::from(1u32 << (level / 2)) * 0.001;
            pole_offset * 1e-5 * pole_offset
        } else {
            0.0
        };
        let mut points: Vec<(f64, f64)> = Vec::with_capacity(ring.len());
        for (i, &(x, y)) in ring.iter().enumerate() {
            // A point whose inverse fails is left out (`RI7H.ec:569`, `0x2ac33` to `0x2ac37`),
            // and with it both tests below: so the ring is closed against its first point only
            // where the inverse of the last planar point succeeds.
            let Some((lat, mut lon)) =
                P::inverse_radians(&self.geom, PlanarPoint { face: -1, x, y }, odd_grid)
            else {
                continue;
            };

            // faithful: site GR-1, `RI7H.ec:573` at `0x2aa98` to `0x2aacf` (`RI3H.ec:517` at
            // `0xf79d` to `0xf7ab`). The eC adds the centroid's longitude to the offset that
            // `wrapLonAt` answers and then takes the shift from the sum,
            // `(w + centroid.lon) - D`; the library takes the shift from the centroid's
            // longitude first, a second time (`0x2aac1`, `0x2aac7`), and adds the offset to
            // that (`0x2aacf`): `(centroid.lon - D) + w`.
            let about = c_lon - REFINED_LON_SHIFT;
            lon = about + wrap_lon_at(lon, about);

            // faithful: `RI7H.ec:576-582`, at `0x2aad3` to `0x2ab15` and `0x2acb8` to
            // `0x2acf4`: the correction of half a turn, taken under `oddGrid` alone. Both
            // longitudes go to degrees by the product with `fl(180 / Pi)` (the two calls
            // through `0x56d30`) before they are subtracted, and the eC's `180` degrees is
            // `Pi` at `0x49ac8`. The library converts the pair again for the second test, to
            // the same two doubles.
            //
            // ACTIVE on aperture 7, INERT on aperture 3, as measured: with the test
            // `odd_grid` taken out, the zones `0132323` and `00515151515151515151` of the
            // aperture-7 grids no longer have the engine's ring and extent, while on
            // the three aperture-3 grids no ring of the suites' samples (odd levels included)
            // and no ring of the censuses of the library's own tests moves. The branch is kept
            // on aperture 3 all the same, as the eC writes it for every grid, since an inert
            // site is ported in the compiled order and its absence from the rings sampled proves
            // nothing of the rings not sampled.
            if odd_grid {
                let apart = lon * math::RAD2DEG - c_lon * math::RAD2DEG;
                if apart < -120.0 {
                    lon += math::PI;
                } else if apart > 120.0 {
                    lon -= math::PI;
                }
            }

            // faithful: the pole branch of aperture 3, `RI3H.ec:529-638`, at `0xf7ff` to
            // `0xfbad`, which the library performs as the eC writes it. Where the latitude,
            // taken to degrees by the product with `fl(180 / Pi)` (`0xf805`), lies beyond
            // 89.999999 either way, the engine gives in place of the point the pole as seen
            // from either side: the points of two planar positions, moved from `(x, y)` by
            // `pole_nudge` one way and the other, along `x` or along `y` by the pole and by
            // where the point lies in the layout (`RI3H.ec:557-593`; the four arms at `0xf81e`
            // to `0xf88e`, `0xfeb0` to `0xfee0`, `0x100b8` to `0x100d4` and `0x101b8` to
            // `0x101da`). Of the terms that the eC multiplies by zero the library forms none
            // for the first position, and for the second adds a zero to the coordinate that
            // does not move (`addpd` at `0xf88a`).
            if T::APERTURE == 3 && (lat * math::RAD2DEG).abs() > POLE_LATITUDE {
                let (first, second) = if lat < 0.0 {
                    if y > 3.0 {
                        ((x, y - pole_nudge), (x + 0.0, y + pole_nudge))
                    } else {
                        ((x - pole_nudge, y), (x + pole_nudge, y + 0.0))
                    }
                } else if x < 1.0 {
                    ((x + pole_nudge, y), (x - pole_nudge, y + 0.0))
                } else {
                    ((x, y + pole_nudge), (x + 0.0, y - pole_nudge))
                };
                for (x, y) in [first, second] {
                    // Both inverses pass `oddGrid` true, whatever the level (`RI3H.ec:594`,
                    // `:617`; `ecx` set to 1 at `0xf8b8` and `0xfa49`), and a failure leaves
                    // that side out.
                    let Some((lat, mut lon)) =
                        P::inverse_radians(&self.geom, PlanarPoint { face: -1, x, y }, true)
                    else {
                        continue;
                    };
                    // `Sgn(out.lat) * 90`, in degrees, taken to radians by the product with
                    // `fl(Pi / 180)` (`RI3H.ec:596`, `:619`; `0xf8c7` to `0xf8f8` and `0xfa58`
                    // to `0xfa89`, through `0x56e08`): the pole itself. The library answers
                    // zero where the latitude is not above zero nor below it.
                    let sign = if lat > 0.0 {
                        1.0
                    } else if lat < 0.0 {
                        -1.0
                    } else {
                        0.0
                    };
                    let pole = sign * 90.0 * math::DEG2RAD;
                    // `RI3H.ec:601-604` and `:624-627`, at `0xf90c` to `0xf950` and `0x10010`
                    // to `0x1004c`, and at `0xfa9d` to `0xfae1` and `0xffa8` to `0xffe4`: the
                    // longitude as the inverse gives it, with neither `wrapLonAt` nor the
                    // sum of site GR-1, which the eC has commented out here, and a correction
                    // of half a turn beyond 95 degrees of the centroid's, whatever the level.
                    let apart = lon * math::RAD2DEG - c_lon * math::RAD2DEG;
                    if apart < -POLE_HALF_TURN_BEYOND {
                        lon += math::PI;
                    } else if apart > POLE_HALF_TURN_BEYOND {
                        lon -= math::PI;
                    }
                    // `RI3H.ec:607-613` and `:630-636`, at `0xf956` to `0xfa19` and `0xfae7`
                    // to `0xfbad`.
                    points.push((pole, lon));
                    drop_a_repeat(&mut points);
                }
            } else {
                points.push((lat, lon));
            }

            // `RI7H.ec:587-594`, `RI3H.ec:642-649`. On aperture 3 the first test runs here
            // again after the pole branch, whether or not that added a point (`0xfbaf` to
            // `0xfc39`). The count is read again after a point is dropped (`0x2abd8`), and a
            // ring of one point skips both tests (`0x2ab61`, `0x2abdd`).
            drop_a_repeat(&mut points);
            let count = points.len();
            if count >= 2 && i == ring.len() - 1 && near(points[0], points[count - 1]) {
                points.pop();
            }
        }
        points
    }

    /// The extent of a zone with geometry as the engine's `getZoneWGS84Extent` computes it
    /// (`RI7H.ec:383-416`, compiled at `0x2ae60`; `RI3H.ec:302-335`, at `0x10320`, instruction
    /// for instruction the same): `[ll.lat, ll.lon, ur.lat, ur.lon]` in radians, or `None`
    /// where the ring is empty, for which the engine leaves the extent as `clear()` wrote it.
    fn extent_radians(&self, id: ZoneId, a: &Address) -> Option<[f64; 4]> {
        // The ring at the automatic refinement (`RI7H.ec:388`, `0x2ae62` to `0x2ae87`).
        let ring = self.refined_ring_radians(id, a, 0);
        if ring.is_empty() {
            return None;
        }
        // faithful: the library performs this function in the eC's own order, re-associating
        // nothing (`0x2ae60` to `0x2aff5`). The centroid is the one `getZoneWGS84Centroid`
        // answers (`RI7H.ec:391`), which the ring's own loop has not adjusted.
        let c_lon = self.centroid_radians(id, a).1;
        let (mut ll_lat, mut ll_lon) = (CLEARED_EXTENT, CLEARED_EXTENT);
        let (mut ur_lat, mut ur_lon) = (-CLEARED_EXTENT, -CLEARED_EXTENT);
        let (mut least, mut greatest) = (99999.0_f64, -99999.0_f64);
        for &(lat, lon) in &ring {
            // `RI7H.ec:397-400`, at `0x2af61` to `0x2af7a` and `0x2af10`: the eC's two tests,
            // compiled as one chain, which is the same, since an offset that the first has
            // lowered cannot pass the second.
            let mut apart = lon - c_lon;
            if apart > math::PI {
                apart -= math::TWO_PI;
            } else if apart < -math::PI {
                apart += math::TWO_PI;
            }
            // `RI7H.ec:402-408`, at `0x2af14` to `0x2af44`. The library's four tests are the
            // eC's on every latitude and offset that is a number, and no ring of the engine
            // that the tests of this module have met holds one that is not.
            if lat > ur_lat {
                ur_lat = lat;
            }
            if lat < ll_lat {
                ll_lat = lat;
            }
            if apart > greatest {
                greatest = apart;
                ur_lon = lon;
            }
            if apart < least {
                least = apart;
                ll_lon = lon;
            }
        }
        // `RI7H.ec:410-413`, at `0x2afb0` to `0x2afdc`.
        if ll_lon < -math::PI {
            ll_lon += math::TWO_PI;
        }
        if ur_lon > math::PI {
            ur_lon -= math::TWO_PI;
        }
        Some([ll_lat, ll_lon, ur_lat, ur_lon])
    }

    /// The zone's boundary in WGS84 degrees with its edges refined, as DGGAL's
    /// `getZoneRefinedWGS84Vertices` gives it (`RI7H.ec:477-480`, which the Z7 grids reach
    /// through `RI7H_Z7.ec:654-657`; `RI3H.ec:396-399`), or the empty ring for a zone without
    /// geometry.
    ///
    /// An edge of a zone is a straight line in the plane of the projection and a curve on the
    /// ellipsoid, which the vertices of [`Grid::vertices`] alone do not follow. Here every
    /// edge is divided into `edge_refinement` parts in that plane and the points of division
    /// are projected with the vertices, so that the ring follows the zone as the engine draws
    /// it. An `edge_refinement` of 0 asks for the engine's own choice, which depends on the
    /// resolution: on the aperture-7 grids twenty parts below resolution 3, fifteen below
    /// resolution 5 and twelve from there on, a hexagon then having 120, 90 and 72 points; on
    /// the aperture-3 grids twenty below resolution 3, fifteen below resolution 5, ten below
    /// resolution 8, eight below resolution 10 and five from there on, a hexagon then having
    /// 120, 90, 60, 48 and 30.
    /// A refinement above [`max_edge_refinement`] is refused with
    /// [`Error::EdgeRefinementOutOfRange`], where the engine sets no bound; the argument is
    /// examined before the zone, so the refusal does not depend on the identifier. That limit
    /// is [`MAX_EDGE_REFINEMENT`], 100,000, unless the environment variable
    /// `RS4DGGS_MAX_EDGE_REFINEMENT` lowers it (read once, never raises it, a malformed value
    /// ignored in silence, absent on a target without an environment, `wasm32-unknown-unknown`
    /// for one). A refinement of 0 is always accepted,
    /// whatever the limit.
    ///
    /// The ring is the engine's own sequence, every coordinate the engine's to the last bit:
    /// it runs anticlockwise as seen from outside the sphere and it is not closed, its last
    /// point not repeating its first. On the aperture-7 grids it begins at the first vertex
    /// of [`Grid::vertices`] and passes the others in their order, save at the few zones on a
    /// pole whose rings degenerate, set out below. On the aperture-3 grids it begins at the
    /// second vertex of [`Grid::vertices`] and passes the others in their order, the first
    /// coming last; the two polar pentagons alone begin elsewhere, the northern at its fifth
    /// vertex and the southern at its first. The engine computes the ring in radians, and it
    /// is given in degrees by the same product with `180 / Pi` that gives [`Grid::vertices`]
    /// theirs.
    ///
    /// **Away from the poles the longitudes are continuous**, and so may lie beyond 180
    /// degrees either way. The
    /// engine unwraps the longitude of every point about that of the zone's centroid, so
    /// that the ring of a zone over the antimeridian does not jump by a turn from one point
    /// to the next, as the vertices of [`Grid::vertices`] do, each of which lies between -180
    /// and 180. A caller that needs each longitude within that interval wraps the points
    /// itself; the continuous ring cannot be had back from a wrapped one.
    ///
    /// The count of points is not always the count of edges times the refinement. An edge
    /// that crosses an interruption of the planar layout is divided on either side of it and
    /// gains a point there, from a refinement of 2. An edge that passes close by the image of
    /// a pole is divided twenty times as finely, since it is there that a straight edge of
    /// the plane lies furthest from its image on the ellipsoid. Two consecutive points
    /// within 1e-11 radians of one another in both coordinates are given once. And a point
    /// whose inverse projection fails is left out: the vertex that [`Grid::vertices`] reports
    /// as (0, 0) is absent here, so that the zone `0000000000000000131`, a hexagon, has five
    /// points at a refinement of 1. A refinement of 1 therefore need not give the vertices
    /// of [`Grid::vertices`], nor as many points; and where it gives the same points, a
    /// longitude may differ from theirs in its last bits, or by a turn. On the aperture-3
    /// grids, at an odd resolution, two corners of a zone that is not centred on a cell of
    /// the coarser grid are drawn 2e-11 of a rhombus's side from where [`Grid::vertices`] has
    /// them, some 1e-8 degrees.
    ///
    /// A vertex on a pole, its latitude exactly 90 or -90 degrees, carries no meaningful
    /// longitude, and the ring may jump by half a turn across it. On the aperture-3 grids a
    /// point of the ring that falls on a pole is not given as it stands: the engine gives in
    /// its place one or several points of exactly that latitude, whose longitude is the
    /// engine's raw one and is often that of neither neighbouring side. Of the zones
    /// examined, on ISEA3H at the engine's own refinement, `A8-0-C` runs from
    /// (89.946, 191.2) to (90, 11.2) and on to (89.946, 11.2); `F0-79-A` has five points on
    /// the pole, at -168.8, between sides at -78.8 and -258.8; and on IGEO7 the zone `013`
    /// has (90, 191.2) between longitudes of 243.6 and 63.6. A caller that reads the ring as
    /// planar coordinates replaces each run of pole vertices by two points at that pole, with
    /// the longitudes of the neighbouring vertices before and after the run. On the rings
    /// examined near the poles the pole was always a vertex of the ring, never enclosed by it.
    ///
    /// Beside the poles the engine's rings degenerate, and these are the engine's rings. At
    /// resolution 19 of the aperture-7 grids a zone centred on a pole has a ring of a single
    /// point, the pole itself, every one of its planar points being moved onto it, and a ring
    /// of fewer than three points is not a polygon. Of the main sample and the fixed list of
    /// the test suites, on the sphere with the pole taken as one point, no ring of the
    /// aperture-3 grids crosses itself, and three on each aperture-7 grid do, at resolutions
    /// 15 to 19. Read as planar longitudes and latitudes, as a drawing or a GeoJSON reader
    /// reads them, 17, 12 and 19 rings on ISEA3H, IVEA3H and RTEA3H cross themselves, from
    /// resolution 10, and 3, 3 and 6 on IGEO7, IVEA7H and RTEA7H. Four rings on each
    /// aperture-3 grid, at resolution 33 and centred on or beside a pole, have no area.
    ///
    /// A zone without geometry gives the empty ring: the null zone, a Z7 identifier of
    /// twenty digits, and an identifier that no cell of this grid has behind it, the very
    /// ones for which [`Grid::vertices`] gives six zero points and [`Grid::area`] nothing.
    /// On the aperture-3 grids the last include the identifiers of a polar root with a
    /// non-zero index, such as `BA-2-A`, which the engine's own quantiser answers beside two
    /// edges of the icosahedron and for which the engine has no ring either; and the
    /// sub-hexes C and D of a polar pentagon, such as `CA-0-C`, which the engine reads from
    /// their text and draws, though no relation of any zone reaches them.
    /// An aperture-7 identifier that has an address but does not read back as text keeps its
    /// ordinary ring.
    ///
    /// # Examples
    ///
    /// ```
    /// let grid = rs4dggs::igeo7();
    /// let zone = grid.zone_from_text("0064156")?;
    ///
    /// // The engine's own choice at resolution 5: twelve points to an edge.
    /// let ring = grid.refined_vertices(zone.id(), 0)?;
    /// assert_eq!(ring.len(), 72);
    /// assert_eq!(ring[0], zone.vertices()[0]);
    ///
    /// // A zone without geometry has no ring, and too fine a refinement is refused.
    /// assert!(grid.refined_vertices(rs4dggs::ZoneId::NULL, 0)?.is_empty());
    /// assert!(grid.refined_vertices(zone.id(), 100_001).is_err());
    ///
    /// // On an aperture-3 grid the engine's choice at resolution 5 is ten points to an edge.
    /// let zone = rs4dggs::isea3h().zone_from_text("C4-49-B")?;
    /// assert_eq!(zone.refined_vertices(0)?.len(), 60);
    /// # Ok::<(), rs4dggs::Error>(())
    /// ```
    pub fn refined_vertices(&self, id: ZoneId, edge_refinement: u32) -> Result<Vec<GeoPoint>> {
        let max = max_edge_refinement();
        if edge_refinement > max {
            return Err(Error::EdgeRefinementOutOfRange {
                refinement: edge_refinement,
                max,
            });
        }
        Ok(self.refined_vertices_within_limit(id, edge_refinement))
    }

    /// The ring of [`Grid::refined_vertices`] once its refinement has passed the limit, which
    /// is at most [`MAX_EDGE_REFINEMENT`] and so within the eC's `int`.
    pub(crate) fn refined_vertices_within_limit(
        &self,
        id: ZoneId,
        edge_refinement: u32,
    ) -> Vec<GeoPoint> {
        debug_assert!(edge_refinement <= MAX_EDGE_REFINEMENT);
        let Some(a) = self.geometry_address(id) else {
            return Vec::new();
        };
        self.refined_ring_radians(id, &a, edge_refinement as i32)
            .into_iter()
            .map(|(lat, lon)| GeoPoint {
                lat: lat * math::RAD2DEG,
                lon: lon * math::RAD2DEG,
            })
            .collect()
    }

    /// The zone's geographic extent in WGS84 degrees, as DGGAL's `getZoneWGS84Extent` gives
    /// it (`RI7H.ec:383-416`, which the Z7 grids reach through `RI7H_Z7.ec:634-637`;
    /// `RI3H.ec:302-335`), or `None` for a zone without geometry.
    ///
    /// It is the extent of the ring of [`Grid::refined_vertices`] at the engine's own choice
    /// of refinement, and not of the vertices alone: the least and the greatest latitude of
    /// that ring, and the longitudes of the two points of it that lie furthest to the west
    /// and to the east of the zone's centroid. Every one of the four is a coordinate of a
    /// point of the ring, the engine's to the last bit, save that a longitude the ring
    /// carries beyond 180 degrees is brought back by a turn: both ends lie between -180 and
    /// 180.
    ///
    /// **Over the antimeridian the western end is the greater**, `ll.lon > ur.lon`, as
    /// [`Extent`] describes. A zone that touches a pole has `ur.lat` of exactly 90, or
    /// `ll.lat` of exactly -90, and its two ends lie some 180 degrees apart: the engine does
    /// not widen such an extent to the whole circle, and the longitudes are those of the
    /// ring, where the pole is one point among the others, or on the aperture-3 grids two.
    /// At resolution 19 of the aperture-7 grids the extent of a zone centred on a pole is
    /// that single point.
    ///
    /// A zone without geometry has no extent, where the engine answers a sentinel, its
    /// cleared extent of `1.7976931348623155e308` and its negation: the null zone, a Z7
    /// identifier of twenty digits, and an identifier that no cell of this grid has behind
    /// it, the very ones for which [`Grid::refined_vertices`] gives the empty ring.
    ///
    /// # Examples
    ///
    /// ```
    /// let grid = rs4dggs::igeo7();
    /// let extent = grid.zone_from_text("0064156")?.extent().unwrap();
    /// assert!(extent.ll.lat < 38.7223 && 38.7223 < extent.ur.lat);
    /// assert!(extent.ll.lon < -9.1393 && -9.1393 < extent.ur.lon);
    ///
    /// // A zone that touches the north pole, and whose extent passes the antimeridian.
    /// let polar = grid.extent(grid.zone_from_text("013")?.id()).unwrap();
    /// assert_eq!(polar.ur.lat, 90.0);
    /// assert!(polar.ll.lon > polar.ur.lon);
    ///
    /// // The null zone has none.
    /// assert_eq!(grid.extent(rs4dggs::ZoneId::NULL), None);
    /// # Ok::<(), rs4dggs::Error>(())
    /// ```
    pub fn extent(&self, id: ZoneId) -> Option<Extent> {
        let a = self.geometry_address(id)?;
        let [ll_lat, ll_lon, ur_lat, ur_lon] = self.extent_radians(id, &a)?;
        Some(Extent {
            ll: GeoPoint {
                lat: ll_lat * math::RAD2DEG,
                lon: ll_lon * math::RAD2DEG,
            },
            ur: GeoPoint {
                lat: ur_lat * math::RAD2DEG,
                lon: ur_lon * math::RAD2DEG,
            },
        })
    }

    // --- neighbours: the topology's own ---

    /// The zone's edge neighbours, as the topology lists them: on every grid here the
    /// engine's own list, in the engine's order, which is the order in which [`Disk`]
    /// walks them. Each topology ports the engine's neighbour search: on an aperture-3
    /// grid `I3HZone::getNeighbors`, and on an aperture-7 grid `I7HZone::getNeighbors`
    /// through the Z7 conversion, which drops four kinds of entry the engine lists and
    /// this crate never hands out (its null zone, the zone itself, a repeat and an
    /// identifier it cannot read back; see [`crate::topologies::HexA7`]'s `neighbors`).
    ///
    /// Along the two icosahedron edges where the aperture-7 engine's odd resolutions
    /// leave a strip with no cells, from base cell 0 to 1 over the north pole and from 1
    /// to 6, from resolution 15, the engine's lists are not consistent with themselves,
    /// and the lists here are the engine's there too, less those four kinds. There the
    /// relation is not symmetric: a zone may list a neighbour that does not list it back,
    /// as in the engine. A kept entry there is the engine's own even where it lies far
    /// from the zone: on IGEO7 the list of `00515151515151510530`, at resolution 18,
    /// holds `00454545454545405216`, `...212` and `...213`, each 35.88 degrees away.
    /// [`Disk`] remembers every zone it has yielded for this reason.
    ///
    /// A zone with no geometry has no neighbours, which is what DGGAL answers for a Z7
    /// level-20 zone and for its null zone; so has an identifier with no zone behind it.
    /// A topology that listed no neighbours of its own would give none here; every
    /// topology this crate ships lists them.
    pub fn neighbors(&self, id: ZoneId) -> Vec<ZoneId> {
        let Ok(a) = self.address(id) else {
            return Vec::new();
        };
        if T::is_null_geometry(&a) {
            return Vec::new();
        }
        T::neighbors(&a)
            .map(|ns| ns.iter().map(I::encode).collect())
            .unwrap_or_default()
    }

    /// Every zone within `k` neighbour steps of the zone, the zone included, as a
    /// lazy iterator of identifiers: the zone first, then ring after ring outward,
    /// each ring in the order [`Disk`] defines, which is the same on every call.
    /// Nothing is computed until the disk is consumed, so the caller controls the
    /// cost by how much it takes, and `k` may be as large as `u32::MAX`; see
    /// [`Disk`] for the order, the cost and the memory, and [`Disk::rings`] for the
    /// same walk a ring at a time. [`crate::Zone::disk`] is the same walk yielding
    /// zones.
    pub fn disk(&self, id: ZoneId, k: u32) -> Disk<&Self> {
        Disk::new(self, id, k)
    }

    // --- hierarchy: exact topology override if provided, else the digits ---

    /// The zone's parents as DGGAL lists them (`getZoneParents`), its primary parent
    /// first: none at resolution 0, and above it one or two on an aperture-7 grid and
    /// one or three on an aperture-3 grid.
    ///
    /// On an aperture-7 grid the primary parent is the zone the identifier names with
    /// its last digit dropped. A centroid child, whose last digit is 0, has that parent
    /// alone; any other zone lies across the boundary between two zones of the
    /// resolution above and has both: `006415654636` has `00641565463` and
    /// `00641565462`. On an aperture-3 grid a centroid child has one parent and every
    /// other zone, a vertex child, three (see [`crate::topologies::HexA3`]'s `parents`).
    ///
    /// Along the two icosahedron edges where the aperture-7 engine's lists are not
    /// consistent with themselves (see [`Grid::neighbors`]), from resolution 15, the list
    /// is the engine's less an identifier it cannot read back and a repeat, and the
    /// account above does not hold of every zone: there a zone that is no centroid child
    /// may have one parent, the first parent may be another zone than the one the
    /// identifier names with its last digit dropped, and a zone may have none (see
    /// [`crate::topologies::HexA7`]'s `parents`).
    ///
    /// The first of those three is met at zones obtained from positions. The other two
    /// were met only at identifiers built by text. On each of the three aperture-7 grids,
    /// of 2,276 identifiers of a broken seam built by text at resolutions 14 to 19 and
    /// compared with the engine, 107 have no parent and 98 a first parent that the
    /// identifier does not name, and not one of those is the zone found at its own
    /// centroid, by this crate or by the engine; and among some 5,600 zones of each grid
    /// obtained from positions, at fifteen resolutions between 0 and 19, there is none of
    /// either kind. That is what was measured of those samples: it is no proof that no position
    /// quantises to such an identifier.
    ///
    /// An identifier with no zone behind it has no parents, nor has a Z7 identifier of
    /// level 20. A topology that exposed no hierarchy of its own would be answered from
    /// the digits of the identifier, the parent being the address with its last digit
    /// dropped; every topology this crate ships exposes the engine's.
    pub fn parents(&self, id: ZoneId) -> Vec<ZoneId> {
        let Ok(a) = self.address(id) else {
            return Vec::new();
        };
        if let Some(ps) = T::parents(&a) {
            return ps.iter().map(I::encode).collect();
        }
        I::parent(id).into_iter().collect()
    }

    /// The zone's children as DGGAL lists them (`getZoneChildren`), in its order, the
    /// centroid child first; none at the finest resolution.
    ///
    /// On an aperture-7 grid a hexagon has thirteen and a pentagon eleven: first the
    /// seven, or six, that the identifier names with one digit more, in the engine's
    /// order and not in the digits', and then six, or five, that lie across the zone's
    /// boundary, each of them a child of a neighbour as well. On an aperture-3 grid a
    /// hexagon has seven and a pentagon six, the centroid child and one on each vertex.
    ///
    /// Along the two icosahedron edges where the aperture-7 engine's lists are not
    /// consistent with themselves (see [`Grid::neighbors`]), in the lists of zones of
    /// resolutions 14 to 18, the list is the engine's less its null zone, an identifier
    /// it cannot read back and a repeat. There a list may be shorter, and a zone need
    /// not be among the children of its parents, nor a zone that the identifier names
    /// with one digit more among the children (see [`crate::topologies::HexA7`]'s
    /// `children`).
    ///
    /// An identifier with no zone behind it has no children, nor has a Z7 identifier of
    /// level 20. A topology that exposed no hierarchy of its own would be answered with
    /// one child for each digit its indexing allows.
    pub fn children(&self, id: ZoneId) -> Vec<ZoneId> {
        let Ok(a) = self.address(id) else {
            return Vec::new();
        };
        if let Some(cs) = T::children(&a) {
            return cs.iter().map(I::encode).collect();
        }
        // A zone at the maximum resolution has no children: a further digit
        // does not fit the packing. This mirrors the parent of a resolution-0
        // zone being absent, and it is what keeps `Indexing::encode` within
        // its packing, since encode bounds nothing itself and a digit past the
        // last slot would shift off the end of the identifier.
        if I::resolution(id) >= I::MAX_RESOLUTION {
            return Vec::new();
        }
        I::child_digits(id)
            .into_iter()
            .map(|d| {
                let mut c: Address = a;
                c.push(d);
                I::encode(&c)
            })
            .collect()
    }

    /// The zone's centroid parent as DGGAL defines it (`getZoneCentroidParent`): the first of
    /// its parents that is itself a centroid child, or `None` where none is, at resolution 0
    /// among others.
    ///
    /// On an aperture-7 grid a parent is a centroid child where its last digit is 0, so a zone
    /// has a centroid parent where one of its one or two parents ends so: `00641565463601`,
    /// whose parents are `0064156546360` and `0064156546361`, has the first, and
    /// `006415654636`, whose parents are `00641565463` and `00641565462`, has none. This is
    /// DGGAL's definition (`RI7H.ec:129-138`) and not its answer on these grids, which was the
    /// null zone at every zone asked (see [`crate::topologies::HexA7`]'s `centroid_parent`).
    ///
    /// On an aperture-3 grid it is `parent0`, the first parent, if that is a centroid child,
    /// and otherwise the first of the two others that is one (`RI3H.ec:1203-1223`). This is not
    /// the parent of which the zone is the centroid child, and the two coincide only for a
    /// centroid child whose one parent is itself a centroid child: a centroid child whose
    /// parent is not has no centroid parent, and a vertex child, the centroid child of no zone,
    /// has one whenever one of its three parents is a centroid child. Over every zone of ISEA3H
    /// at resolutions 0 to 6, 10,944 in all, the engine answers a centroid parent for 8,512,
    /// which is the parent of which the zone is the centroid child for 1,232; and 2,432 of the
    /// 3,664 centroid children have none.
    pub fn centroid_parent(&self, id: ZoneId) -> Option<ZoneId> {
        let a = self.address(id).ok()?;
        match T::centroid_parent(&a) {
            Some(r) => r.as_ref().map(I::encode),
            None => I::parent(id),
        }
    }

    /// Whether the zone is a centroid child, as DGGAL answers (`isZoneCentroidChild`): above
    /// resolution 0, the child that shares the centre of its primary parent.
    ///
    /// On an aperture-7 grid it is a zone whose identifier ends in the digit 0, and so no
    /// zone of resolution 0; a Z7 identifier of level 20 is read by its last digit as any
    /// other, as the engine reads it. On an aperture-3 grid it is as
    /// [`crate::topologies::HexA3`]'s `is_centroid_child` describes, every zone of
    /// resolution 0 among them. An identifier with no zone behind it is none.
    pub fn is_centroid_child(&self, id: ZoneId) -> bool {
        let Ok(a) = self.address(id) else {
            // Not a zone at all, so not any parent's centroid child.
            return false;
        };
        T::is_centroid_child(&a).unwrap_or(true)
    }

    // --- sub-zones: required topology override, no congruent default ---

    /// The guard the sub-zone methods share, in py4dggs's own order: refuse a
    /// grid with no sub-zone order, which none of those this crate ships is,
    /// then refuse a depth that would reach past the maximum resolution. Both
    /// cases would otherwise reach the
    /// depth-dependent geometric generators, which are not written for them
    /// and misbehave quietly: only `zone_from_geo` bounds a resolution, so an
    /// over-deep request would produce zones past the maximum with no error.
    ///
    /// The identifier is validated first, ahead of both, because an identifier
    /// with no cell behind it is not a question about depth at all, and
    /// because a `Result` is one of the few places where this crate can name
    /// that failure rather than answer as the null zone.
    fn check_sub_zone_depth(&self, id: ZoneId, depth: u8) -> Result<()> {
        self.address(id)?;
        if !T::HAS_SUB_ZONE_ORDER {
            return Err(Error::NoSubZoneOrder);
        }
        let target = u32::from(I::resolution(id)) + u32::from(depth);
        if target > u32::from(I::MAX_RESOLUTION) {
            return Err(Error::InvalidZone(format!(
                "relative_depth {depth} would reach resolution {target}, past max_resolution {}",
                I::MAX_RESOLUTION
            )));
        }
        Ok(())
    }

    /// How many sub-zones lie `depth` resolution levels below this zone.
    ///
    /// On an aperture-7 grid it is the engine's closed form, the length of the order that
    /// [`Grid::sub_zones`] describes: 13, 55, 379 and 2,449 at depths 1 to 4 for a hexagon,
    /// and 11, 46, 316 and 2,041 for a pentagon. That is more than the `7^depth` zones that
    /// the identifier's digits name below the zone, since the order holds every zone that
    /// lies across the zone's boundary as well. An entry of the order that is
    /// [`ZoneId::NULL`] is counted.
    ///
    /// A depth that would reach past `max_resolution` is refused with [`Error::InvalidZone`],
    /// where the engine's own `getSubZonesCount` answers its ordinary formula regardless, for a
    /// level past its own maximum (see Departures from the engine in the crate's README).
    pub fn count_sub_zones(&self, id: ZoneId, depth: u8) -> Result<u64> {
        self.check_sub_zone_depth(id, depth)?;
        if depth == 0 {
            // Depth 0 is the zone itself, which the aperture-3 topology's own
            // override says too (eC `I3HSubZones.ec:1787-1790`). It is settled
            // here rather than there because it holds for any DGGS, and
            // because it keeps depth 0 away from the depth-dependent
            // generators, which are not written for it.
            return Ok(1);
        }
        T::count_sub_zones(&I::decode(id), depth).ok_or(Error::NoSubZoneOrder)
    }

    /// The first sub-zone, in the topology's own order, `depth` levels below
    /// this zone.
    ///
    /// One deliberate departure, at depth 0, where this answers the zone itself. The engine's own
    /// `getFirstSubZone` has no depth-0 case: it quantises a vertex of the zone at the zone's own
    /// level there and so names a neighbour, for about two thirds of zones on aperture 3 (1,970,
    /// 1,986 and 1,983 of 3,000 measured live), while its `getSubZones` answers the zone itself
    /// at that depth.
    ///
    /// At the zones on the edges of the rhombi where [`Grid::sub_zones`] corrects the engine's
    /// order on aperture 3, this is the first of the corrected order, which need not be the
    /// engine's own `getFirstSubZone`: the correction can move the whole order, its start
    /// included.
    ///
    /// On an aperture-7 grid it is the first entry of the order that [`Grid::sub_zones`]
    /// describes: the sub-zone beside the vertex of the zone that stands uppermost on the
    /// icosahedral net, and never the centroid descendant, which stands in the middle. It is
    /// the engine's `getFirstSubZone` wherever that answers. At the pentagon of the south
    /// pole (base cell 11, every digit 0) at an odd level, at depth 2, the engine's call ends
    /// the process, and this answers the first of the list that the engine's `getSubZones`
    /// gives there. The answer is [`ZoneId::NULL`], and no error, where the first entry of
    /// the order is no zone, as [`Grid::sub_zones`] explains.
    ///
    /// At depth 0 this answers the zone itself on these grids too, and the engine's
    /// `getFirstSubZone` does not. Of the 5,607, 5,579 and 5,574 zones compared on IGEO7,
    /// IVEA7H and RTEA7H it names a neighbour of the zone at 5,594, 5,570 and 5,566, its
    /// null zone at 6, 5 and 4, the zone itself at 3, 1 and 1 and, in the broken seams,
    /// another identifier that is no neighbour at 4, 3 and 3. In those seams the engine's
    /// own `getSubZones` at depth 0 holds another identifier than the zone at 10, 5 and 5
    /// of the zones compared, where this crate's list is the zone.
    ///
    /// A depth that would reach past `max_resolution` is refused with [`Error::InvalidZone`],
    /// where the engine's own `getFirstSubZone` quantises a vertex at a level past its own
    /// maximum and names one of the level-34 identifiers its own parser cannot read back (see
    /// Departures from the engine in the crate's README).
    pub fn first_sub_zone(&self, id: ZoneId, depth: u8) -> Result<ZoneId> {
        self.check_sub_zone_depth(id, depth)?;
        if depth == 0 {
            return Ok(id); // See the note in `count_sub_zones`.
        }
        T::first_sub_zone(&I::decode(id), depth)
            .map(|a| I::encode(&a))
            .ok_or(Error::NoSubZoneOrder)
    }

    /// Every sub-zone `depth` levels below this zone, in the topology's own
    /// order, up to [`max_materialised_sub_zones`] of them: [`MAX_MATERIALISED_SUB_ZONES`],
    /// four million, unless the environment variable `RS4DGGS_MAX_MATERIALISED_SUB_ZONES`
    /// lowers it (read once, never raises it, a malformed value ignored). A longer list is
    /// refused with [`Error::TooManySubZones`], which carries the limit in force.
    ///
    /// On an aperture-3 grid the order is the engine's own, save at some zones on the edges of
    /// the rhombi, where the engine's order is at fault: a rounding of one unit in the last
    /// place decides whether a point of the generator's scan is carried across an interruption
    /// of the layout that it should only reach, and where it is carried wrongly the engine's
    /// order names a zone more than once, names a zone that is not a descendant of this one, or
    /// omits a descendant. There this crate answers each descendant once, as many as the
    /// engine's own count, equal to the children at depth 1, in the generator's own scanline
    /// order; [`Grid::first_sub_zone`], [`Grid::sub_zone_at_index`] and [`Grid::sub_zone_index`]
    /// follow the same corrected order. Every other order is the engine's, down to the bits of
    /// the generator's centroids.
    ///
    /// Measured over the edges of every rhombus, the fault falls on them alone: at resolutions 13,
    /// 19 and 21, on some sub-hexes along the top edges of the northern rhombi 2, 4 and 8 and the
    /// left edges of the southern rhombi 3, 5 and 9, at every depth measured (1 to 8); on the
    /// left-edge hexagons of rhombus 5 at resolutions 16, 20, 22, 24 and 26, and of rhombi 1 and 9
    /// at 24, 26 and 28, at some even depths; and on the southern polar pentagon at resolutions 11
    /// and 25, at depths 3, 5 and 7. The correction follows the mechanism, not this list.
    ///
    /// On an aperture-7 grid the order is the engine's own at every zone, entry for entry, and
    /// nothing in it is corrected. The engine generates the sub-zones as centroids, in
    /// scanlines across the zone that begin beside the vertex of it that stands uppermost on
    /// the icosahedral net, and quantises each at the sub-zone level: at depth 1 a hexagon has
    /// scanlines of 1, 4, 3, 4 and 1 zones, and at depth 2 of 3, 6, 7, 8, 7, 8, 7, 6 and 3.
    /// The centroid descendant, which is the zone that the identifier names with `depth`
    /// zeros more, stands in the middle of the order and never at its head: at index 6 of
    /// 13, 27 of 55, 189 of 379, and at the same indices of a pentagon's 11, 46 and 316. The
    /// order has no relation to the digits of the identifiers. The first three sub-zones of
    /// `006415654636` at depth 1 are `0064156546306`, `0064156546342` and `0064156546361`,
    /// which the digits make children of three different zones. At depth 1 the sub-zones are
    /// the zone's [`Grid::children`], in another order, save at some zones of the broken
    /// seams, among sub-zones of resolutions 16 and 18 as measured, where the engine's
    /// `getSubZones` and its `getZoneChildren` part: the order of `0000000000000001644`
    /// holds thirteen zones, and its children are four. From depth 3 some of the zones that
    /// the digits name below the zone are no sub-zones of it at all, and the order holds as
    /// many zones that the digits name below its neighbours.
    ///
    /// On an aperture-7 grid an entry may be [`ZoneId::NULL`]. In the broken seams, along the
    /// edges of the icosahedron from base cell 0 to 1 and from 1 to 6, among sub-zones of
    /// resolution 15 and finer, the engine's quantiser gives some centroids of an order its
    /// null zone, or an identifier that the engine itself cannot read back: the first was met
    /// at the resolutions 15, 17 and 19, the second at 16 and 18 and, in an order, at 19 too.
    /// Each such entry is the null zone here, at the engine's place, so that the length of the
    /// order is [`Grid::count_sub_zones`] and the place of every other entry is the engine's: of
    /// the 17,053 sub-zones of `010004000400` at depth 5, 99 are the null zone. Whatever else the
    /// engine's order holds there is kept as it is: a zone named twice, a zone that is no
    /// descendant of this one, a descendant left out. A caller who walks an order skips the
    /// null zone, which has no geometry and is refused by every method that validates an
    /// identifier.
    ///
    /// A depth that would reach past `max_resolution` is refused the same way as on
    /// [`Grid::count_sub_zones`] and [`Grid::first_sub_zone`], with [`Error::InvalidZone`],
    /// where the engine's own `getSubZones` answers nothing there.
    pub fn sub_zones(&self, id: ZoneId, depth: u8) -> Result<Vec<ZoneId>> {
        let n = self.count_sub_zones(id, depth)?;
        if depth == 0 {
            return Ok(vec![id]); // See the note in `count_sub_zones`.
        }
        let limit = max_materialised_sub_zones();
        if n > limit {
            return Err(Error::TooManySubZones { count: n, limit });
        }
        let subs = T::sub_zones(&I::decode(id), depth).ok_or(Error::NoSubZoneOrder)?;
        Ok(subs.iter().map(I::encode).collect())
    }

    /// Where `sub` sits in this zone's sub-zone order, or `None` if it is no
    /// sub-zone of it, or one that has no index. On an aperture-3 grid it is found as the
    /// eC's own generic method finds it (`dggrs.ec:115-131`): the ordered sequence is built
    /// and the position sought in it. On an aperture-7 grid it is found by the engine's own
    /// walk of the order, described below. On every grid `Some(i)` is answered only where
    /// [`Grid::sub_zone_at_index`] at `i` is `sub`.
    ///
    /// Both identifiers are validated first, as the other sub-zone methods validate
    /// theirs, and one with no cell behind it, the null zone among them, is refused
    /// with [`Error::InvalidZone`]. That departs from DGGAL, which on all three
    /// aperture-7 grids answers 0 for the null zone, or any fabricated identifier,
    /// within itself, and -1 for one against a real zone: it compares the two
    /// levels before anything else, and finds them equal in the first case. On the
    /// aperture-3 grids, asked for the index of the null zone within a real zone,
    /// the engine never returns.
    ///
    /// It departs from DGGAL a second time, also on purpose, for two different zones
    /// at one level. The engine answers 0 whenever the two share a level, whatever
    /// they are (`getSubZoneIndex`, `RI7H.ec:229-230`, which the Z7 grids reach
    /// through `RI7H_Z7.ec:599-601`), and so reports a zone's neighbour as its
    /// sub-zone at index 0. This crate answers `None`: at depth 0 a zone's only
    /// sub-zone is the zone itself, and index 0 names that zone, not the one asked
    /// about.
    ///
    /// A known inconsistency, kept on purpose. When `sub` is at this zone's own
    /// resolution or coarser, the answer is settled here, above the topology, before
    /// the sub-zone order is consulted: a zone is its own sub-zone at index 0, and a
    /// coarser zone is none. So on a grid with no sub-zone order, which none of those
    /// this crate ships is, `sub_zone_index(z, z)` would answer `Ok(Some(0))` while
    /// `count_sub_zones(z, 0)` and the other sub-zone methods return
    /// [`Error::NoSubZoneOrder`]. Only a finer `sub` reaches the order and the
    /// error. py4dggs's `grid.py` has exactly the same inconsistency, answering
    /// depth 0 above the topology here and checking for an order first in the
    /// others. What remains of it on the shipped grids is a Z7 identifier of twenty
    /// digits, which is its own sub-zone at index 0 here and has no sub-zones by the
    /// other methods, which refuse it as lying past the finest resolution.
    ///
    /// On an aperture-7 grid the index is found as the engine's `getSubZoneIndex` finds it
    /// (`RI7H.ec:224-239`), without building the order. The engine first asks whether `sub`
    /// is a sub-zone of the zone at all, by its `zoneHasSubZone`: it is one where a point 99
    /// hundredths of the way from its centroid to one of its six vertices lies in the zone.
    /// It then generates the centroids of the order from its head, and the index is that of
    /// the first which lies within 1e-11 of the centroid of `sub` in each coordinate. This
    /// crate answers that index only once it has confirmed it: `Some(i)` where the entry of
    /// the order at `i` is `sub`, and `None` otherwise. So `None` is "no sub-zone of it, or
    /// one that has no index", and every index answered is one that
    /// [`Grid::sub_zone_at_index`] gives back as `sub`.
    ///
    /// Away from the broken seams every entry of an order has its index, which is its place,
    /// at the pentagons too, where the engine's `getSubZoneAtIndex` is at fault and its walk
    /// is not. In the broken seams, along the edges of the icosahedron from base cell 0 to 1
    /// and from 1 to 6, among sub-zones of resolution 15 and finer, an order may name a zone
    /// that has none: the engine's test of a sub-zone refuses it, or no centroid of the
    /// order lies within 1e-11 of its own, and the engine answers -1 there, as this crate
    /// answers `None`. Of the order of `0000000000000005` at depth 1 the entry at 4,
    /// `00000000000000005`, and the entry at 9, `00000000000000055`, have no index. The null
    /// zone, which such an order may hold, is refused as `sub` like any identifier with no
    /// cell behind it, and so has none either.
    ///
    /// The walk costs one centroid for each entry before the one sought, so its cost grows
    /// with the index: in a release build the middle of the 825,259 sub-zones of a hexagon
    /// at depth 7 is some twelve milliseconds away, and the last twice that, while a zone
    /// that is no sub-zone is refused in microseconds. The walk is therefore bound as a list
    /// is: where the order is longer than [`max_materialised_sub_zones`], which under the
    /// ceiling is any depth from 8 below a hexagon or a pentagon, the request is refused with
    /// [`Error::TooManySubZones`] before any walk is begun, where the engine walks an order
    /// of any length. A caller who pages through a deeper order with
    /// [`Grid::sub_zone_at_index`] holds the index already.
    ///
    /// A third departure, within the cap that bounds this search: a sub-zone at resolution 33 has
    /// its index here, where DGGAL's own `getSubZoneIndex` answers -1 for some 83 to 84 per cent
    /// of them. The oracle suites, which compare every sub-zone at resolution 33 of the zones they
    /// sample at resolutions 31 and 32, find it answering -1 for 7,300 of 8,725 on ISEA3H, 7,277
    /// of 8,725 on IVEA3H and 7,270 of 8,711 on RTEA3H. Its `zoneHasSubZone` test fails there on
    /// the needle's level-34 centroid child, whose index overflows 51 bits, so the walk that would
    /// find it is never entered; this crate's search does not go through `zoneHasSubZone` and
    /// answers the index it finds.
    ///
    /// A fourth departure, a consequence of the corrected order: at the zones where
    /// [`Grid::sub_zones`] corrects the engine's order, the index answered here is the sub-zone's
    /// place in the corrected order, which is unique. DGGAL's own `getSubZoneIndex` searches its
    /// own faulty order. For a descendant that order omits, its `zoneHasSubZone` accepts the zone,
    /// but its search of the order's points, to within 1e-11, meets none of them, and it answers
    /// -1. For any other zone it answers the zone's place in its own order, which need not be this
    /// one.
    pub fn sub_zone_index(&self, id: ZoneId, sub: ZoneId) -> Result<Option<u64>> {
        self.address(id)?;
        self.address(sub)?;
        let (sub_level, level) = (I::resolution(sub), I::resolution(id));
        if sub_level == level {
            // Depth 0 is the zone itself, so a zone is its own sub-zone at
            // index 0. Treating this as the "not a descendant" case would
            // report a real sub-zone as absent.
            return Ok(if sub == id { Some(0) } else { None });
        }
        if sub_level < level {
            return Ok(None);
        }
        let depth = sub_level - level;
        // The bound of a list is the bound of the search too, whichever way the topology
        // makes it: a walk of the order costs as its length, as building it does.
        let count = self.count_sub_zones(id, depth)?;
        let limit = max_materialised_sub_zones();
        if count > limit {
            return Err(Error::TooManySubZones { count, limit });
        }
        if let Some(found) = T::sub_zone_index(&I::decode(id), &I::decode(sub)) {
            // The topology's index is its walk's and no more: it is answered where the
            // entry of the order there is the zone asked about.
            return Ok(found.filter(|&i| self.sub_zone_at_index(id, depth, i) == Ok(sub)));
        }
        let subs = self.sub_zones(id, depth)?;
        Ok(subs.iter().position(|&s| s == sub).map(|i| i as u64))
    }

    /// The sub-zone at `index`, `depth` levels below this zone. Generic, as
    /// the eC's own `dggrs.ec:133-149`.
    ///
    /// At depth 0 (so index 0, its only one) this answers the zone itself, where the engine's
    /// `getSubZoneAtIndex` inherits `getFirstSubZone`'s neighbour there; see
    /// [`Grid::first_sub_zone`] for the departure and its figures. For any other index, a topology
    /// that orders its sub-zones without materialising them (both of those this crate ships,
    /// through the engine's own generators and their index skip) is asked directly, which reaches
    /// any depth to the maximum resolution without building the order; only a topology with no
    /// such shortcut falls back to indexing into [`Grid::sub_zones`], which
    /// [`max_materialised_sub_zones`] still guards.
    ///
    /// At the zones on the edges of the rhombi where [`Grid::sub_zones`] corrects the engine's
    /// order on aperture 3, this follows the corrected order, through the same generators;
    /// elsewhere on aperture 3 it is the engine's own index skip, bit for bit.
    ///
    /// On an aperture-7 grid it is the entry at `index` of the order that [`Grid::sub_zones`]
    /// describes, at every index, [`ZoneId::NULL`] where that entry is no zone, with no error.
    /// It is not limited by the length of the order: the search reckons the width of each
    /// scanline before the one that holds the index and generates one sub-zone, so that its
    /// cost grows with the number of those scanlines, some `(10 / 3) * 7^(depth / 2)` at most
    /// (134 million below a zone of resolution 0 at depth 19, a third of a second in a release
    /// build), and never with the number of sub-zones.
    ///
    /// The answer departs from the engine at the pentagons. At an odd depth the engine's
    /// `getSubZoneAtIndex` answers another zone than its `getSubZones` lists at that index,
    /// from the scanline after the pentagon's own: for `00` at depth 1 it answers `030` and
    /// `034` at the indices 9 and 10, where its list holds `006` and `043`, and many of its
    /// answers are no sub-zones of the pentagon at all. Here the two agree at every zone.
    pub fn sub_zone_at_index(&self, id: ZoneId, depth: u8, index: u64) -> Result<ZoneId> {
        let count = self.count_sub_zones(id, depth)?;
        if index >= count {
            return Err(Error::IndexOutOfRange { index, count });
        }
        if index == 0 {
            return self.first_sub_zone(id, depth);
        }
        if let Some(a) = T::sub_zone_at_index(&I::decode(id), depth, index) {
            return Ok(I::encode(&a));
        }
        // `index < count` is established, but `count` comes from the topology's
        // closed form and the sequence from its enumeration; should the two ever
        // disagree, that is reported rather than allowed to panic.
        let subs = self.sub_zones(id, depth)?;
        usize::try_from(index)
            .ok()
            .and_then(|i| subs.get(i))
            .copied()
            .ok_or(Error::IndexOutOfRange { index, count })
    }

    // --- the grid's facts, as DGGAL's `DGGRS` gives them ---

    /// Refuses a level beyond the finest resolution, as [`Error::ResolutionOutOfRange`].
    fn check_level(&self, level: u8) -> Result<()> {
        if level > I::MAX_RESOLUTION {
            return Err(Error::ResolutionOutOfRange {
                res: level,
                max: I::MAX_RESOLUTION,
            });
        }
        Ok(())
    }

    /// The factor by which a level has more zones than the one above it, DGGAL's
    /// `getRefinementRatio` (`RI7H.ec:48`; `RI3H.ec:43`): 7 on the aperture-7 grids, 3 on the
    /// aperture-3 grids.
    pub fn refinement_ratio(&self) -> u8 {
        T::APERTURE
    }

    /// The most parents a zone has, DGGAL's `getMaxParents` (`RI7H.ec:49`; `RI3H.ec:44`): 2 on the
    /// aperture-7 grids, 3 on the aperture-3 grids.
    pub fn max_parents(&self) -> u8 {
        facts::max_parents(T::APERTURE)
    }

    /// The most children a zone has, DGGAL's `getMaxChildren` (`RI7H.ec:51`; `RI3H.ec:46`): 13 on
    /// the aperture-7 grids, 7 on the aperture-3 grids.
    pub fn max_children(&self) -> u8 {
        facts::max_children(T::APERTURE)
    }

    /// The most edge neighbours a zone has, DGGAL's `getMaxNeighbors` (`RI7H.ec:50`; `RI3H.ec:45`):
    /// 6 on every grid.
    pub fn max_neighbors(&self) -> u8 {
        facts::max_neighbors(T::APERTURE)
    }

    /// The depth at which a zone has about 65,536 sub-zones, DGGAL's `get64KDepth`
    /// (`dggrs.ec:224-227`): 6 on the aperture-7 grids, 10 on the aperture-3 grids.
    pub fn depth_64k(&self) -> u8 {
        facts::depth_64k(T::APERTURE)
    }

    /// The deepest relative depth at which DGGAL lists a zone's sub-zones, its `getMaxDepth`
    /// (`dggrs.ec:230-241`): the deepest at which a hexagon has at most `2^28` of them, 9 on
    /// the aperture-7 grids and 17 on the aperture-3 grids.
    pub fn max_depth(&self) -> u8 {
        facts::max_depth(T::APERTURE)
    }

    /// The number of zones at `level`, DGGAL's `countZones` (`RI7H.ec:42-45`; `RI3H.ec:37-40`):
    /// ten times the refinement ratio to the power of the level, plus two, that is twelve
    /// pentagons and the hexagons. It is [`Error::ResolutionOutOfRange`] beyond
    /// [`Grid::max_resolution`], where the engine goes on answering from tables that run out,
    /// with counts that wrap and then settle at 2.
    ///
    /// # Examples
    ///
    /// ```
    /// let grid = rs4dggs::igeo7();
    /// assert_eq!(grid.count_zones(0)?, 12);
    /// assert_eq!(grid.count_zones(1)?, 72);
    /// assert!(grid.count_zones(20).is_err());
    /// # Ok::<(), rs4dggs::Error>(())
    /// ```
    pub fn count_zones(&self, level: u8) -> Result<u64> {
        self.check_level(level)?;
        Ok(facts::count_zones(T::APERTURE, level))
    }

    /// The reference area of a zone at `level` in square metres, DGGAL's `getRefZoneArea`
    /// (`dggrs.ec:266-284`): the
    /// area of the whole ellipsoid divided by the count of zones, which is neither a hexagon's
    /// area nor a pentagon's, and is the figure the engine converts levels from and to sizes
    /// with. It is [`Error::ResolutionOutOfRange`] beyond [`Grid::max_resolution`].
    ///
    /// # Examples
    ///
    /// ```
    /// let grid = rs4dggs::igeo7();
    /// // Twelve pentagons share the ellipsoid at level 0.
    /// assert!((grid.ref_zone_area(0)? - 4.250547e13).abs() < 1e8);
    /// # Ok::<(), rs4dggs::Error>(())
    /// ```
    pub fn ref_zone_area(&self, level: u8) -> Result<f64> {
        self.check_level(level)?;
        Ok(facts::ref_zone_area(T::APERTURE, level))
    }

    /// The size in metres of the sub-zones `depth` levels below a zone of `level`, DGGAL's
    /// `getMetersPerSubZoneFromLevel` (`dggrs.ec:305-308`): the square root of the reference
    /// area of the level `level + depth`. It is [`Error::ResolutionOutOfRange`] where that level
    /// is beyond [`Grid::max_resolution`], carrying the sum, saturated at 255.
    pub fn meters_per_sub_zone(&self, level: u8, depth: u8) -> Result<f64> {
        let target = level.saturating_add(depth);
        self.check_level(target)?;
        Ok(facts::meters_per_sub_zone(T::APERTURE, target))
    }

    /// The level whose reference area an area of `square_metres` reaches, DGGAL's
    /// `getLevelFromRefZoneArea` (`dggrs.ec:250-265`): the first level at which the count of
    /// zones is at least the area of the ellipsoid over it.
    ///
    /// It mirrors the engine for every input, and answers up to `max_resolution + depth_64k`,
    /// 25 on the aperture-7 grids and 43 on the aperture-3 grids, which is beyond the finest
    /// resolution: an area of 0 or NaN answers that, and a negative or infinite area answers 0.
    /// Nor is it the inverse of [`Grid::ref_zone_area`]: on ISEA3H the reference area of level
    /// 9 answers 10, since the quotient rounds above the count. A caller that needs a resolution
    /// the grid has takes the lesser of the answer and [`Grid::max_resolution`].
    pub fn level_from_ref_zone_area(&self, square_metres: f64) -> u8 {
        facts::level_from_ref_zone_area(T::APERTURE, square_metres)
    }

    /// The level whose sub-zones `depth` levels below it measure `metres`, DGGAL's
    /// `getLevelFromMetersPerSubZone` (`dggrs.ec:300-303`): [`Grid::level_from_ref_zone_area`]
    /// of the square of `metres`, less `depth`, and not below 0. It mirrors the engine for every input, and
    /// answers beyond the finest resolution where that does, up to `max_resolution + depth_64k`.
    ///
    /// # Examples
    ///
    /// ```
    /// let grid = rs4dggs::igeo7();
    /// // The level to ask for, so that the zones are no coarser than a 100-metre pixel.
    /// let level = grid.level_from_meters_per_sub_zone(100.0, 0);
    /// assert_eq!(level, 12);
    /// assert!(grid.meters_per_sub_zone(level, 0)? <= 100.0);
    /// # Ok::<(), rs4dggs::Error>(())
    /// ```
    pub fn level_from_meters_per_sub_zone(&self, metres: f64, depth: u8) -> u8 {
        let level = i32::from(facts::level_from_ref_zone_area(
            T::APERTURE,
            metres * metres,
        ));
        // The eC's `Max(0, level - relativeDepth)`.
        (level - i32::from(depth)).max(0) as u8
    }

    // --- the zones of a level, enumerated ---

    /// The finest level at which this crate enumerates zones: every zone of a level
    /// ([`Grid::zones`]), and the zones of a bounding box ([`Grid::zones_in_box`]). It is 33
    /// on the aperture-3 grids, their last level, and 14 on the aperture-7 grids, whose last
    /// level is 19; either enumeration refuses a finer level with
    /// [`Error::ResolutionOutOfRange`], which carries this limit. A service that offers the
    /// zones of a box advertises this level, and not [`Grid::max_resolution`], as the finest
    /// it answers at.
    ///
    /// The crate enumerates by the lattice in which the engine itself orders the zones of a
    /// level. On the aperture-7 grids that lattice is not whole from level 16: along two edges
    /// of the icosahedron, from base cell 0 to base cell 1 over the north pole and from base
    /// cell 1 to base cell 6, there are cells that host no zone under an identifier of its
    /// own, and the engine's own answers there are not consistent with one another. The limit
    /// is set at 14, below the first level at which anything of the kind was met.
    ///
    /// The zones of a box have a second reason for it. They are found by a search of the
    /// lattice that takes a zone to lie within one side of a cell of the corner of the cell
    /// that hosts it: to level 14 no zone was found beyond 0.44 of a side, and from level 15
    /// there are zones along those two edges that lie farther than a side, which the search
    /// would not find. Within the limit the answer is whole wherever it was compared.
    pub fn max_box_level(&self) -> u8 {
        if T::APERTURE == 7 {
            crate::zones::APERTURE_7_ENUMERATION_LIMIT.min(I::MAX_RESOLUTION)
        } else {
            I::MAX_RESOLUTION
        }
    }

    /// The walk of the lattice of `level` from its first zone, or
    /// [`Error::ResolutionOutOfRange`] for a level beyond [`Grid::max_box_level`], carrying
    /// that limit as its `max`.
    pub(crate) fn lattice_walk(&self, level: u8) -> Result<crate::zones::Walk> {
        let max = self.max_box_level();
        if level > max {
            return Err(Error::ResolutionOutOfRange { res: level, max });
        }
        Ok(crate::zones::Walk::new(level, T::lattice_edge(level)))
    }

    /// The zones that one cell of the lattice of `level` hosts, each with its key in the
    /// order of the level: the cell at `row` and `col` of a rhombus, `root` 0 to 9, or the one
    /// cell of a polar root, 10 or 11. Nothing for a cell the lattice does not have.
    pub(crate) fn lattice_cell(
        &self,
        level: u8,
        root: u8,
        row: u64,
        col: u64,
    ) -> Vec<(ZoneId, u64)> {
        let hosted = if root < 10 {
            T::cell_zones(level, root, row, col)
        } else {
            T::polar_zones(level, root)
        };
        hosted
            .into_iter()
            .map(|(a, key)| (I::encode(&a), key))
            .collect()
    }

    /// Where the lattice of `level` holds `zone`: its root, its row, its column and its key,
    /// as [`Topology::locate`] gives them. [`Error::InvalidZone`] unless `zone` is an
    /// identifier that the enumeration of `level` itself yields: one of that level, which the
    /// lattice places, and which is the identifier of the address it decodes to (a Z7
    /// identifier with anything but padding after its last digit decodes to a zone and is not
    /// that zone's identifier).
    pub(crate) fn lattice_place(&self, level: u8, zone: ZoneId) -> Result<(u8, u64, u64, u64)> {
        // The level is compared first, so that a zone of another level is refused before the
        // lattice is asked to place it.
        let place = (I::resolution(zone) == level)
            .then(|| I::decode(zone))
            .filter(|a| I::encode(a) == zone)
            .and_then(|a| T::locate(&a));
        place.ok_or_else(|| {
            Error::InvalidZone(format!(
                "zone id {:#018x} is no zone of level {level} of {}, so the zones of that \
                 level cannot be entered after it",
                zone.0, self.name
            ))
        })
    }

    /// Every zone of `level`, as a lazy iterator of identifiers in the order of DGGAL's
    /// `listZones` for the whole world: see [`crate::Zones`] for that order, which on the
    /// aperture-7 grids is not the order of the Z7 identifiers, and for the cost.
    /// [`crate::Zones::after`] enters the same sequence after any zone of the level, at a
    /// cost that does not grow with the zone's place in it, which is how the sequence is
    /// taken in pages.
    ///
    /// Nothing is computed until the iterator is consumed, and no list is ever held: a level
    /// of 5.6e16 zones is as safe to ask for as level 0, and costs what the caller takes.
    ///
    /// It is [`Error::ResolutionOutOfRange`] beyond [`Grid::max_box_level`], the error
    /// carrying that limit: on the aperture-7 grids levels 15 to 19 are refused, though the
    /// grid has them.
    ///
    /// # Examples
    ///
    /// ```
    /// let grid = rs4dggs::igeo7();
    /// // The twelve zones of level 0, in the engine's order, which is not that of the
    /// // identifiers.
    /// let level_0: Vec<String> = grid.zones(0)?.map(|z| grid.text_id(z)).collect();
    /// assert_eq!(
    ///     level_0,
    ///     ["01", "06", "02", "07", "03", "08", "04", "09", "05", "10", "00", "11"]
    /// );
    /// // Every zone of a level, each once.
    /// assert_eq!(grid.zones(3)?.count() as u64, grid.count_zones(3)?);
    /// // A deep level is asked for as a coarse one is: only what is taken is computed.
    /// assert_eq!(grid.zones(14)?.take(5).count(), 5);
    /// assert!(grid.zones(15).is_err());
    /// # Ok::<(), rs4dggs::Error>(())
    /// ```
    pub fn zones(&self, level: u8) -> Result<crate::Zones<&Self>> {
        Ok(crate::Zones::new(self, self.lattice_walk(level)?))
    }

    // --- the zones of a level whose extent meets a bounding box ---

    /// The point of the sphere, a latitude and a longitude in radians, that the inverse
    /// projection gives for the point `(x, y)` of the plane of the ten rhombi, or `None` where
    /// it gives none. It is the point the search for the cells near a box measures a block of
    /// the lattice from.
    pub(crate) fn lattice_point(&self, x: f64, y: f64) -> Option<(f64, f64)> {
        P::inverse_radians(&self.geom, PlanarPoint { face: -1, x, y }, false)
    }

    /// Whether the extent of `zone` meets the box `bbox`, `[ll.lat, ll.lon, ur.lat, ur.lon]`
    /// in radians, by [`crate::zones::intersects`]: the rule of [`Grid::zones_in_box`] for one
    /// zone. The extent, which is the dearer part, is computed only for a zone whose centroid
    /// lies within `reach` of the box, `reach` being a bound on how far a point of a zone's
    /// boundary lies from its centroid: beyond it the extent cannot meet the box. A zone
    /// without geometry meets nothing.
    pub(crate) fn zone_meets_box(&self, zone: ZoneId, bbox: &[f64; 4], reach: f64) -> bool {
        let Some(a) = self.geometry_address(zone) else {
            return false;
        };
        // A centroid that the inverse projection cannot place prunes nothing: the zone is
        // then decided by its extent alone.
        let near = P::inverse_radians(&self.geom, T::planar_centroid(&a), self.odd_grid(zone))
            .is_none_or(|(lat, lon)| crate::zones::near_box(lat, lon, reach, bbox));
        near && self
            .extent_radians(zone, &a)
            .is_some_and(|extent| crate::zones::intersects(&extent, bbox))
    }

    /// The walk of the lattice of `level` and the search of it for the cells near `bbox`,
    /// after the two refusals of a query of a box, made before any work and in this order:
    /// [`Error::ResolutionOutOfRange`] for a level beyond [`Grid::max_box_level`], then the
    /// refusals of [`crate::zones::box_radians`] for a box that is none.
    pub(crate) fn box_search(
        &self,
        level: u8,
        bbox: &Extent,
    ) -> Result<(crate::zones::Walk, crate::zones::Search)> {
        let walk = self.lattice_walk(level)?;
        let bbox = crate::zones::box_radians(bbox)?;
        let search =
            crate::zones::Search::new(bbox, T::lattice_edge(level), self.count_zones(level)?);
        Ok((walk, search))
    }

    /// The zones of `level` whose extent meets the bounding box `bbox`, as a lazy iterator of
    /// identifiers in the order of DGGAL's `listZones`: a zone is yielded if and only if its
    /// extent ([`Grid::extent`], a bounding box of the zone) meets the box by the engine's own
    /// test of two extents. See [`crate::ZonesInBox`] for the rule in full, for how the answer
    /// differs from the engine's own, for the order and for the cost.
    ///
    /// The box is an [`Extent`], in degrees, as [`Grid::extent`] gives one: `ll` its south and
    /// its west, `ur` its north and its east. A west above the east is a box that runs
    /// eastwards over the antimeridian. A box may be a segment or a point, and may be the
    /// whole world, for which [`Grid::zones`] is the cheaper call.
    ///
    /// [`crate::ZonesInBox::after`] enters the answer after any zone of the level, which is
    /// how it is taken in pages; nothing is computed until the iterator is consumed, and no
    /// list is ever held.
    ///
    /// # What a caller must know
    ///
    /// - The extent of a zone is a bounding box of it, so the answer covers a little more than
    ///   the box: a zone whose extent meets the box at a corner is yielded although the zone
    ///   itself does not touch the box.
    /// - A point exactly on a pole yields no zone, since the test is strict and no extent
    ///   reaches beyond a pole; and a point on the equator, or on one of the four meridians
    ///   along which the rhombi are joined, may yield none, where the extents of the zones on
    ///   either side end at the line or short of it. [`Grid::zone_from_geo`] gives the zone
    ///   at a point.
    /// - On the aperture-7 grids the answer is given to level 14, [`Grid::max_box_level`], and
    ///   not to the grids' last level, 19.
    /// - A step of the iterator never costs by its place in the answer nor by the area of the
    ///   box: it costs by the length, in cells, of the stretch of the rim of the box that
    ///   runs along a row of a rhombus, which is a few cells for most boxes and the whole of
    ///   a side that lies on one of those four meridians (11.2 degrees east, 168.8 west, 78.8
    ///   west and 101.2 east on the grids of the registry). A side of `n` cells puts at least
    ///   `n` zones in the answer, so a caller that refuses a box whose
    ///   [`Grid::estimate_zones_in_box`] exceeds a limit of its own bounds a step by that
    ///   many cells. [`crate::ZonesInBox`] gives the measured cost of a cell.
    ///
    /// # Errors
    ///
    /// Refused before any work, in this order: a level beyond [`Grid::max_box_level`]
    /// ([`Error::ResolutionOutOfRange`], carrying that limit); a coordinate that is not finite
    /// ([`Error::NonFinite`]); a latitude outside -90 to 90 degrees
    /// ([`Error::LatitudeOutOfRange`]); a longitude outside -180 to 180 degrees
    /// ([`Error::LongitudeOutOfRange`]); a south above the north ([`Error::InvertedBox`]).
    ///
    /// # Examples
    ///
    /// ```
    /// use rs4dggs::{Extent, GeoPoint};
    /// let grid = rs4dggs::isea3h();
    /// // The box about Lisbon: from 38 to 39 degrees north, from 10 to 9 degrees west.
    /// let lisbon = Extent {
    ///     ll: GeoPoint { lat: 38.0, lon: -10.0 },
    ///     ur: GeoPoint { lat: 39.0, lon: -9.0 },
    /// };
    /// // At level 2 one zone holds the whole box; at level 3 two zones share it.
    /// let texts = |level| -> rs4dggs::Result<Vec<String>> {
    ///     Ok(grid.zones_in_box(level, &lisbon)?.map(|z| grid.text_id(z)).collect())
    /// };
    /// assert_eq!(texts(2)?, ["B4-2-A"]);
    /// assert_eq!(texts(3)?, ["B2-5-C", "B4-2-B"]);
    /// // A deep level is asked for as a coarse one is: only what is taken is computed.
    /// assert_eq!(grid.zones_in_box(33, &lisbon)?.take(5).count(), 5);
    /// // A box beyond a pole is refused, and so is a level the grid has not.
    /// let beyond = Extent { ur: GeoPoint { lat: 91.0, lon: -9.0 }, ..lisbon };
    /// assert!(grid.zones_in_box(3, &beyond).is_err());
    /// assert!(grid.zones_in_box(34, &lisbon).is_err());
    /// # Ok::<(), rs4dggs::Error>(())
    /// ```
    pub fn zones_in_box(&self, level: u8, bbox: &Extent) -> Result<crate::ZonesInBox<&Self>> {
        let (walk, search) = self.box_search(level, bbox)?;
        Ok(crate::ZonesInBox::new(self, walk, search))
    }

    // --- an upper bound on the zones of a box, by arithmetic ---

    /// An upper bound on how many zones [`Grid::zones_in_box`] yields for `level` and `bbox`,
    /// by arithmetic alone: no zone is visited, and the cost is the same for every box and
    /// every level. It is what a service refuses a request by before it does any work, or
    /// sizes a page with.
    ///
    /// # What it is
    ///
    /// With `n` the count of zones of the level ([`Grid::count_zones`]) and `w` the width of
    /// a zone, the square root of the sphere's area over `n`, it is the area of the box grown
    /// by `1.75 w` on every side, over the area of a zone, rounded up, plus 3; and never more
    /// than `n + 3`. The grown box is cut at the poles and is nowhere wider than the whole
    /// circle. The area is taken on the sphere on which the zones are equal in area: the two
    /// latitudes of the box are first converted from geodetic to authalic, as the grid itself
    /// converts a latitude ([`crate::GridConfig`]'s `authalic`). A cap above 70 degrees of
    /// latitude is some 0.81 per cent larger on that sphere than its geodetic latitudes say,
    /// which on a box of some two hundred thousand zones is more than the margin about its
    /// rim.
    ///
    /// # What it rests on
    ///
    /// It is a bound by measurement, and it is not proved. The measurement: it was not under
    /// the count on any box on which the two were compared, at every level at which
    /// [`Grid::zones_in_box`] answers, beside the poles, over the antimeridian and about the
    /// vertices of the icosahedron, on boxes from a single point to caps and bands of half a
    /// million and of a million and a half zones.
    ///
    /// Away from the poles an argument supports the margin, where the zones of the level are
    /// disjoint and of one area: a zone is yielded when its extent meets the box; a corner
    /// of an extent lies at most `sqrt(2)` radii of the zone from its centroid, and the zone
    /// within one radius of it; no zone was found with a radius above 0.7224 `w`; so a zone
    /// of the answer lies wholly within `(sqrt(2) + 1) * 0.7224 w`, which is `1.744 w`, of
    /// the box.
    ///
    /// **Upon a pole the bound is measured and not argued.** The box is grown in longitude
    /// by `1.75 w` over the cosine of the latitude, which beside a pole is less than the
    /// zones reach by arc; and at the finest level of the aperture-3 grids the extents of
    /// ten zones reach the north pole, from centroids up to `1.7 w` away. The thinnest margin
    /// was measured there: a wedge upon the north pole at level 33, a millionth of a zone
    /// width high and one degree wide, holds 8 zones and is estimated at 9, on the three
    /// aperture-3 grids alike. That margin is one zone, and the estimate holds it by the 3
    /// that it adds; at levels 31 and 32 the thinnest margin measured is 6 zones, and on the
    /// aperture-7 grids 7.
    ///
    /// A caller that must not be surprised keeps a limit of its own on what it takes from
    /// the iterator.
    ///
    /// # How far above the count it lies
    ///
    /// The margin is a rim of some one zone about the box, so that it weighs by the size of
    /// the answer, and most on a box that is long and thin. As measured: an answer of under
    /// ten zones may be estimated at seventeen to twenty-one times its size (a point away
    /// from the poles, which one to three zones hold, is estimated at 16 zones, and at 12 to
    /// 16 at the two coarsest levels); one of ten to a hundred zones at
    /// up to some five and a half times; one of a hundred to a thousand at up to some
    /// twice; one of a thousand to ten thousand at no more than a third above it; one of ten
    /// thousand or more at no more than an eighth; and on a large box the excess shrinks
    /// towards one row of zones about its rim, eight parts in a thousand on a band about the
    /// world of half a million zones.
    ///
    /// # Errors
    ///
    /// Those of [`Grid::zones_in_box`], in the same order and before any arithmetic: a level
    /// beyond [`Grid::max_box_level`] ([`Error::ResolutionOutOfRange`], carrying that limit),
    /// then a box that is none ([`Error::NonFinite`], [`Error::LatitudeOutOfRange`],
    /// [`Error::LongitudeOutOfRange`], [`Error::InvertedBox`]).
    ///
    /// # Examples
    ///
    /// ```
    /// use rs4dggs::{Extent, GeoPoint};
    /// let grid = rs4dggs::isea3h();
    /// let iberia = Extent {
    ///     ll: GeoPoint { lat: 36.0, lon: -10.0 },
    ///     ur: GeoPoint { lat: 44.0, lon: 3.0 },
    /// };
    /// // Sixty zones of level 7 meet the box; the estimate is above that, and costs nothing.
    /// let estimate = grid.estimate_zones_in_box(7, &iberia)?;
    /// assert_eq!(grid.zones_in_box(7, &iberia)?.count(), 60);
    /// assert!((60..=120).contains(&estimate));
    /// // At the finest level the same box holds some two hundred million million zones,
    /// // which the estimate says without visiting one.
    /// assert!(grid.estimate_zones_in_box(33, &iberia)? > 100_000_000_000_000);
    /// // The whole world is never estimated above the zones of the level and three.
    /// let world = Extent {
    ///     ll: GeoPoint { lat: -90.0, lon: -180.0 },
    ///     ur: GeoPoint { lat: 90.0, lon: 180.0 },
    /// };
    /// assert_eq!(grid.estimate_zones_in_box(7, &world)?, grid.count_zones(7)? + 3);
    /// # Ok::<(), rs4dggs::Error>(())
    /// ```
    pub fn estimate_zones_in_box(&self, level: u8, bbox: &Extent) -> Result<u64> {
        self.lattice_walk(level)?;
        let [south, west, north, east] = crate::zones::box_radians(bbox)?;
        let on_the_sphere = [
            P::sphere_latitude(&self.geom, south),
            west,
            P::sphere_latitude(&self.geom, north),
            east,
        ];
        Ok(crate::zones::estimate(
            &on_the_sphere,
            self.count_zones(level)?,
        ))
    }

    // --- the compaction of a set of zones ---

    /// DGGAL's `compactZones` of the zones of one level (`RI3H.ec:161-186` with `:2250-2386`;
    /// `RI7H_Z7.ec:581-597` and `RI7H.ec:172-197` with `:3754-3862`): a shorter list of zones
    /// that stands for the same set, in which a coarser zone takes the place of those of its
    /// sub-zones that it alone holds, wherever the set has all of its sub-zones of the level.
    /// It is the engine's algorithm and the engine's answer, zone for zone and in its order,
    /// on the six grids. The two apertures compact by different rules.
    ///
    /// # What the answer is on the aperture-3 grids
    ///
    /// - The compaction goes two levels at a time, from the level of the set upwards: a zone
    ///   two levels up takes the place of the seven sub-zones that lie wholly within it,
    ///   when the set has all thirteen of its sub-zones of that depth or, for each of the six
    ///   on its vertices that is missing, has already kept the zone at that one's centre. So
    ///   the zones of the answer are of the level of the set or an even number of levels
    ///   above it; and every zone of level 1, all thirty-two, become the twelve of level 0.
    /// - **The zones of the answer overlap.** A zone on a vertex of a coarser zone is shared
    ///   by three such zones, and stays in the answer unless all three are in it: the
    ///   thirteen sub-zones of one zone, two levels down, compact to that zone and the six
    ///   on its vertices, whose areas sum to fifteen thirteenths of the set's. The area that
    ///   a compacted list covers is that of the set, not the sum over the list.
    /// - The sub-zones, at the level of the set, of the zones of the answer
    ///   ([`Grid::sub_zones`]) hold every zone of the set, and no other where the set has no
    ///   gap: the zones of a box, the sub-zones of a zone, a whole level.
    /// - **A set with a gap may be covered beyond itself, by two routes.** The engine takes
    ///   a coarser zone when the six neighbours of the sub-zone at its centre are in the
    ///   set, and does not look for that sub-zone itself: the thirteen sub-zones above
    ///   without the one at the centre compact to the same seven zones, which cover it. And
    ///   from the second pass on, a sub-zone on a vertex of the coarser zone that the pass
    ///   lacks counts as present where the answer already holds the one zone, two levels
    ///   finer, at its centre: the coarser zone is then taken, and covers the part of itself
    ///   about that vertex, whatever of it the set held. A page of a longer answer is a set
    ///   with a gap where it ends before a zone whose neighbours it holds; the zone on a
    ///   pole, which comes last in the order of [`Grid::zones_in_box`], is one.
    /// - The answer is in ascending order of the identifiers, which puts coarser zones
    ///   before finer ones, and does not depend on the order of `zones`. A zone given twice
    ///   counts once, and no zone gives no zone.
    ///
    /// # What the answer is on the aperture-7 grids
    ///
    /// - The compaction goes one level at a time, from the level of the set upwards: a zone
    ///   one level up is taken when the set has every one of its children (thirteen for a
    ///   hexagon, eleven for a pentagon), and a zone gives way when every one of its parents
    ///   was taken. The lists read are the engine's own (those of [`Grid::children`] and
    ///   [`Grid::parents`], before anything this crate leaves out of them). A zone at the
    ///   centre of its parent has that parent alone and gives way to it; any other has two,
    ///   and stays unless both were taken. So the zones of the answer are of the level of the
    ///   set or of any level above it, and every zone of a level becomes the twelve of
    ///   level 0. The zones that an identifier names with one digit more are seven of the
    ///   thirteen children, and do not compact by themselves.
    /// - **The zones of the answer overlap, more than on aperture 3.** The thirteen children
    ///   of `0064` compact to thirteen zones: `0064` in the place of `00640`, the child at its
    ///   centre, and the twelve others, which it overlaps. Their areas sum to nineteen
    ///   thirteenths of the set's.
    /// - Away from the two broken seams, the sub-zones, at the level of the set, of the zones
    ///   of the answer are the zones of the set, and **no set is covered beyond itself**: a
    ///   coarser zone is taken only where the set has every one of its children, the one at
    ///   its centre among them, so the thirteen above without any one of them come back as
    ///   they are.
    /// - **The answer is in the order of the engine's own identifiers**, which within a level
    ///   is the order of [`Grid::zones`] and is not that of the Z7 identifiers. A level and
    ///   the odd level above it share their cells there, so that the zones of the two
    ///   alternate in the list and coarser zones do not come first: in the answer above
    ///   `0064`, of level 2, stands between `03323` and `00641`, of level 3. A caller that
    ///   wants the coarser zones first sorts the answer by level, stably. The answer does not
    ///   depend on the order of `zones`; a zone given twice counts once, and no zone gives
    ///   no zone.
    /// - **An identifier that the engine reads as another zone comes back as that zone's.**
    ///   The engine compacts the zones that the identifiers are read to and names each zone
    ///   of the answer anew. An identifier given by its bits that names the child a pentagon
    ///   does not have comes back as the engine reads it, a child that the pentagon has:
    ///   `0x05FF_FFFF_FFFF_FFFF`, whose bits spell `002` (a text that
    ///   [`Grid::zone_from_text`] refuses), comes back as `003`. So do some identifiers of
    ///   the two broken seams.
    /// - **In the two broken seams**, along the edges of the icosahedron from base cell 0 to
    ///   1 and from 1 to 6 from resolution 15, the engine's lists of parents and children are
    ///   not consistent with themselves, and the answer is the engine's all the same: a zone
    ///   for which the engine finds no parent (`00000000000000001311`) leaves the answer and
    ///   no zone stands for it; a zone is taken although the engine lists its null zone among
    ///   the children; and the sub-zones of the answer need not be the set. One departure is
    ///   made there: where the engine's answer holds an identifier that it cannot itself
    ///   read back, this crate's answer lacks it. That was met only for identifiers of
    ///   resolutions 17 and 19 that are not the zone found at their own centroid, and for no
    ///   set made of zones found from positions.
    ///
    /// The cost is that of the slice, in time and in memory: a bounded amount of work for
    /// each zone of a pass, at most seventeen passes on aperture 3 and nineteen on aperture
    /// 7, and three sets of the zones of the slice and of the coarser zones taken for them.
    ///
    /// # Errors
    ///
    /// Checked before any work, in this order:
    ///
    /// - [`Error::InvalidZone`] for an identifier that is no zone of the grid, wherever it
    ///   stands in the slice. The null zone is one, and so is a Z7 identifier of twenty
    ///   digits, which the engine reads as its null zone: the engine leaves both out of the
    ///   set, and this crate refuses them, as it refuses every identifier it reads no zone
    ///   for.
    /// - [`Error::MixedLevels`] for zones of more than one level, where the engine loses
    ///   zones: an answer of this method is therefore not a slice to hand it again.
    ///
    /// # Examples
    ///
    /// ```
    /// let grid = rs4dggs::isea3h();
    /// // The thirteen sub-zones, two levels down, of the zone of level 2 that holds Lisbon.
    /// let zone = grid.zone_from_text("B4-2-A")?.id();
    /// let thirteen = grid.sub_zones(zone, 2)?;
    /// // The zone takes the place of seven of them; the six on its vertices stay.
    /// let compact: Vec<String> = grid
    ///     .compact_zones(&thirteen)?
    ///     .into_iter()
    ///     .map(|z| grid.text_id(z))
    ///     .collect();
    /// assert_eq!(
    ///     compact,
    ///     ["B4-2-A", "C2-11-A", "C2-19-A", "C2-2C-A", "C4-E-A", "C4-11-A", "C4-19-A"]
    /// );
    /// // Every zone of a level compacts to the twelve zones of level 0.
    /// let level_4: Vec<_> = grid.zones(4)?.collect();
    /// assert_eq!(grid.compact_zones(&level_4)?.len(), 12);
    ///
    /// // On an aperture-7 grid the thirteen children of a zone give way at its centre alone.
    /// let grid = rs4dggs::igeo7();
    /// let zone = grid.zone_from_text("0064")?.id();
    /// let compact: Vec<String> = grid
    ///     .compact_zones(&grid.children(zone))?
    ///     .into_iter()
    ///     .map(|z| grid.text_id(z))
    ///     .collect();
    /// assert_eq!(compact.len(), 13);
    /// assert_eq!(compact[2..5], ["03323", "0064", "00641"]);
    /// assert!(!compact.contains(&"00640".to_string()));
    /// # Ok::<(), rs4dggs::Error>(())
    /// ```
    pub fn compact_zones(&self, zones: &[ZoneId]) -> Result<Vec<ZoneId>> {
        let mut addresses = Vec::with_capacity(zones.len());
        for &zone in zones {
            let a = self.address(zone)?;
            if T::is_null_geometry(&a) {
                return Err(Error::InvalidZone(format!(
                    "zone id {:#018x} is of a level at which {} has no zone, so it cannot be \
                     compacted",
                    zone.0, self.name
                )));
            }
            addresses.push(a);
        }
        if let Some(&first) = zones.first() {
            let level = I::resolution(first);
            if let Some(other) = zones
                .iter()
                .map(|&z| I::resolution(z))
                .find(|&l| l != level)
            {
                return Err(Error::MixedLevels { level, other });
            }
        }
        let compacted = T::compact(&addresses).ok_or(Error::NoSubZoneOrder)?;
        Ok(compacted.iter().map(I::encode).collect())
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::FRAC_PI_2;

    use crate::indexings::I3h;
    use crate::indexings::Z7;
    use crate::math;
    use crate::projections::Isea;
    use crate::registry::{get_grid, igeo7, ivea7h, rtea7h};
    use crate::topologies::HexA3;
    #[cfg(feature = "oracle")]
    use crate::topologies::HexA7;
    use crate::{
        Error, Extent, GeoPoint, Grid, GridConfig, Indexing, Projection, Topology, ZoneId,
    };

    use super::{
        CENTROID_LON_LIMIT, CLEARED_EXTENT, MAX_EDGE_REFINEMENT, MAX_MATERIALISED_SUB_ZONES,
        POLE_HALF_TURN_BEYOND, POLE_LATITUDE, REFINED_LON_SHIFT, automatic_edge_divisions,
        edge_refinement_limit, materialised_sub_zones_limit, max_edge_refinement,
        max_materialised_sub_zones, refined_odd_grid,
    };

    /// Builds ISEA3H directly, as `hex_a3.rs`'s own tests do.
    fn isea3h() -> Grid<Isea, HexA3, I3h> {
        Grid::new(GridConfig::default(), "ISEA3H").unwrap()
    }

    fn a3_id(grid: &Grid<Isea, HexA3, I3h>, text: &str) -> ZoneId {
        grid.zone_from_text(text).unwrap().id()
    }

    /// A fabricated identifier: base cell 12, which the four-bit base field
    /// can hold and the icosahedron has no cell for.
    const BASE_12: ZoneId = ZoneId(0xCFFF_FFFF_FFFF_FFFF);

    #[test]
    fn zone_from_geo_lisbon() {
        let z = igeo7().zone_from_geo(38.7223, -9.1393, 5).unwrap();
        assert_eq!(z.text_id(), "0064156");
        assert_eq!(z.id(), ZoneId(0x0d0d_dfff_ffff_ffff));
    }
    #[test]
    fn resolution_out_of_range_is_an_error() {
        assert!(matches!(
            igeo7().zone_from_geo(0.0, 0.0, 20),
            Err(Error::ResolutionOutOfRange { res: 20, max: 19 })
        ));
    }
    /// The centroid is DGGAL v0.0.6's own, `getZoneWGS84Centroid` on `ISEA7H_Z7` in
    /// degrees, which py4dggs, whence it was first taken, also gives.
    #[test]
    fn centroid_and_vertices() {
        let z = igeo7().zone_from_text("0064156").unwrap();
        let c = z.centroid();
        assert_eq!((c.lat, c.lon), (38.61687542349096, -8.874508932649084));
        assert_eq!(z.vertices().len(), 6);
    }
    /// The neighbours of a hexagon and of a pentagon, in the engine's order: one for each
    /// vertex of the zone's centroid child, in the order the engine traces them.
    #[test]
    fn neighbours_are_in_the_engines_order() {
        let n: Vec<String> = igeo7()
            .zone_from_text("0064156")
            .unwrap()
            .neighbors()
            .iter()
            .map(|z| z.text_id())
            .collect();
        assert_eq!(
            n,
            [
                "0064150", "0064154", "0064141", "0064143", "0064105", "0064152"
            ]
        );
        let p: Vec<String> = igeo7()
            .zone_from_text("0000000")
            .unwrap()
            .neighbors()
            .iter()
            .map(|z| z.text_id())
            .collect();
        assert_eq!(p, ["0000005", "0000003", "0000004", "0000001", "0000006"]);
    }
    /// A zone of the finest resolution has no children, as the engine answers for it; one
    /// level above, a hexagon has the thirteen the engine lists, in its order: the seven its
    /// text names with one digit more, and six others.
    #[test]
    fn children_at_max_resolution_are_empty() {
        let g = igeo7();
        let z = g.zone_from_geo(38.7223, -9.1393, 19).unwrap();
        assert_eq!(z.text_id(), "006415654630655235512");
        assert!(z.children().is_empty());
        let z = g.zone_from_geo(38.7223, -9.1393, 18).unwrap();
        assert_eq!(z.text_id(), "00641565463065523551");
        let children: Vec<String> = z.children().iter().map(|c| c.text_id()).collect();
        assert_eq!(
            children,
            [
                "006415654630655235510",
                "006415654630655235511",
                "006415654630655235515",
                "006415654630655235514",
                "006415654630655235516",
                "006415654630655235512",
                "006415654630655235513",
                "006415654630655213262",
                "006415654630655235553",
                "006415654630655235501",
                "006415654630655235535",
                "006415654630655235144",
                "006415654630655213226",
            ]
        );
    }

    /// The hierarchy about `0064156` on IGEO7, each answer the engine's own: its two parents,
    /// `006415` and `006414`, neither a centroid child, so that it has no centroid parent and,
    /// its last digit being 6, is no centroid child itself; its ancestors two levels up, among
    /// them `00422`, the second parent of `006415`, and `00645`, the second of `006414`, but
    /// not `00642`; and the centroid child `00641560`, whose one parent it is.
    #[test]
    fn geometric_hierarchy() {
        let g = igeo7();
        let zone = |text: &str| g.zone_from_text(text).unwrap();
        let z = zone("0064156");
        assert_eq!(z.parent().unwrap().text_id(), "006415");
        let parents: Vec<String> = z.parents().iter().map(|p| p.text_id()).collect();
        assert_eq!(parents, ["006415", "006414"]);
        assert!(!z.is_centroid_child());
        assert_eq!(z.centroid_parent(), None);
        for p in z.parents() {
            assert!(z.is_immediate_child_of(&p) && !p.is_immediate_child_of(&z));
            assert!(p.is_ancestor_of(&z) && !z.is_ancestor_of(&p));
        }
        assert!(z.parent().unwrap().parent().unwrap().is_ancestor_of(&z));
        for (text, is) in [
            ("00641", true),
            ("00422", true),
            ("00645", true),
            ("00642", false),
        ] {
            assert_eq!(zone(text).is_ancestor_of(&z), is, "{text}");
        }
        let c = zone("00641560");
        assert_eq!(c.parents(), vec![z]);
        assert!(c.is_centroid_child() && c.is_immediate_child_of(&z));
        assert_eq!(z.children()[0], c);
        // A parent that is itself a centroid child is the centroid parent.
        assert_eq!(
            zone("00641565463601").centroid_parent(),
            Some(zone("0064156546360"))
        );
        assert!(zone("00").parent().is_none());
        assert!(!zone("00").is_centroid_child());
    }

    /// The digit path lies within the engine's hierarchy, at every zone of resolutions 0 to
    /// 4 of IGEO7, 28,020 in all: the first parent is the zone the text names with its last
    /// digit dropped, and the zones the text names with one digit more are the first seven
    /// children, or the first six under a pentagon, as a set, for the engine lists them in
    /// another order than the digits'. A hexagon has thirteen children and a pentagon eleven.
    #[test]
    fn the_digit_path_lies_within_the_hierarchy() {
        use std::collections::BTreeSet;
        let g = igeo7();
        let mut level: Vec<ZoneId> = (0..12)
            .map(|base| Z7::encode(&crate::Address::new(base, &[])))
            .collect();
        let mut zones = 0usize;
        for res in 0..=4u8 {
            let mut next = Vec::new();
            for &id in &level {
                let z = g.zone(id);
                let text = z.text_id();
                assert_eq!(z.resolution(), res, "{text}");
                assert_eq!(
                    z.parent().map(|p| p.text_id()),
                    (res > 0).then(|| text[..text.len() - 1].to_string()),
                    "the first parent of {text}"
                );
                let by_digits: BTreeSet<ZoneId> = Z7::child_digits(id)
                    .into_iter()
                    .map(|d| g.zone_from_text(&format!("{text}{d}")).unwrap().id())
                    .collect();
                let children = g.children(id);
                assert_eq!(
                    (by_digits.len(), children.len()),
                    if z.is_pentagon() { (6, 11) } else { (7, 13) },
                    "the children of {text}"
                );
                let first: BTreeSet<ZoneId> = children[..by_digits.len()].iter().copied().collect();
                assert_eq!(first, by_digits, "the first children of {text}");
                next.extend(by_digits);
                zones += 1;
            }
            level = next;
        }
        assert_eq!(zones, 12 + 72 + 492 + 3_432 + 24_012);
    }
    #[test]
    fn null_geometry_at_level_20() {
        let g = igeo7();
        let id =
            <crate::indexings::Z7 as crate::Indexing>::encode(&crate::Address::new(1, &[1; 20]));
        assert_eq!(g.centroid(id), crate::GeoPoint { lat: 0.0, lon: 0.0 });
        assert_eq!(
            g.vertices(id),
            vec![crate::GeoPoint { lat: 0.0, lon: 0.0 }; 6]
        );
        assert!(g.neighbors(id).is_empty());
    }
    /// The texts of an order, the null zone as its own text.
    fn texts_of(g: &crate::AnyGrid, order: &[ZoneId]) -> Vec<String> {
        order.iter().map(|&s| g.text_id(s)).collect()
    }

    /// The three aperture-7 grids, which share every answer below: the order is generated in
    /// the 5x6 layout, before any projection.
    const APERTURE_7: [&str; 3] = ["IGEO7", "IVEA7H", "RTEA7H"];

    /// The sub-zone order of an aperture-7 zone is the engine's: its scanlines across the zone,
    /// 1, 4, 3, 4 and 1 zones at depth 1, the centroid descendant in the middle and never at
    /// the head, in no relation to the digits. Every list is DGGAL's own `getSubZones`.
    #[test]
    fn aperture7_sub_zones_are_in_the_engines_order() {
        for name in APERTURE_7 {
            let g = get_grid(name).unwrap();
            let z = g.zone_from_text("006415654636").unwrap();
            assert_eq!(g.count_sub_zones(z, 1), Ok(13), "{name}");
            let order = g.sub_zones(z, 1).unwrap();
            assert_eq!(
                texts_of(&g, &order),
                [
                    "0064156546306",
                    "0064156546342",
                    "0064156546361",
                    "0064156546363",
                    "0064156546324",
                    "0064156546365",
                    "0064156546360",
                    "0064156546362",
                    "0064156546033",
                    "0064156546364",
                    "0064156546366",
                    "0064156546215",
                    "0064156546251",
                ],
                "{name}"
            );
            assert_eq!(g.first_sub_zone(z, 1), Ok(order[0]), "{name}");
            for (index, &sub) in order.iter().enumerate() {
                assert_eq!(g.sub_zone_at_index(z, 1, index as u64), Ok(sub), "{name}");
            }
            assert_eq!(
                g.sub_zone_at_index(z, 1, 13),
                Err(Error::IndexOutOfRange {
                    index: 13,
                    count: 13
                }),
                "{name}"
            );
            // Depth 0 is the zone itself, as on aperture 3.
            assert_eq!(g.count_sub_zones(z, 0), Ok(1), "{name}");
            assert_eq!(g.sub_zones(z, 0), Ok(vec![z]), "{name}");
            assert_eq!(g.first_sub_zone(z, 0), Ok(z), "{name}");
            assert_eq!(g.sub_zone_at_index(z, 0, 0), Ok(z), "{name}");

            // The centroid descendant stands in the middle: index 27 of 55, 189 of 379.
            let z = g.zone_from_text("0064156").unwrap();
            for (depth, count, middle, descendant) in [
                (1, 13, 6, "00641560"),
                (2, 55, 27, "006415600"),
                (3, 379, 189, "0064156000"),
                (4, 2_449, 1_224, "00641560000"),
            ] {
                assert_eq!(g.count_sub_zones(z, depth), Ok(count), "{name}");
                let order = g.sub_zones(z, depth).unwrap();
                assert_eq!(order.len() as u64, count, "{name}");
                assert_eq!(g.text_id(order[middle]), descendant, "{name}");
                assert_eq!(g.first_sub_zone(z, depth), Ok(order[0]), "{name}");
            }
            assert_eq!(
                g.text_id(g.first_sub_zone(z, 4).unwrap()),
                "00641561311",
                "{name}"
            );
            assert_eq!(
                g.text_id(g.sub_zone_at_index(z, 4, 2_448).unwrap()),
                "00641566466",
                "{name}"
            );
        }
    }

    /// At a pentagon, at an odd depth, the engine's `getSubZoneAtIndex` answers another zone
    /// than its `getSubZones` lists at that index, from the scanline after the pentagon's own:
    /// for `00` at depth 1 it answers `030` and `034` at the indices 9 and 10, where the list
    /// holds `006` and `043`. Here the zone at an index is the entry of the list.
    #[test]
    fn aperture7_sub_zone_at_index_is_the_entry_of_the_order_at_a_pentagon() {
        for name in APERTURE_7 {
            let g = get_grid(name).unwrap();
            let z = g.zone_from_text("00").unwrap();
            assert_eq!(g.count_sub_zones(z, 1), Ok(11), "{name}");
            let order = g.sub_zones(z, 1).unwrap();
            assert_eq!(
                texts_of(&g, &order),
                [
                    "013", "023", "005", "001", "053", "004", "000", "003", "033", "006", "043"
                ],
                "{name}"
            );
            for (index, &sub) in order.iter().enumerate() {
                assert_eq!(g.sub_zone_at_index(z, 1, index as u64), Ok(sub), "{name}");
            }
            assert_eq!(g.text_id(g.sub_zone_at_index(z, 1, 9).unwrap()), "006");
            assert_eq!(g.text_id(g.sub_zone_at_index(z, 1, 10).unwrap()), "043");
            // At depths 3 and 5 as well, at every index.
            for depth in [3, 5] {
                let order = g.sub_zones(z, depth).unwrap();
                for (index, &sub) in order.iter().enumerate() {
                    assert_eq!(
                        g.sub_zone_at_index(z, depth, index as u64),
                        Ok(sub),
                        "{name}, depth {depth}, index {index}"
                    );
                }
            }
        }
    }

    /// At the pentagon of the south pole at an odd level, at depth 2, the engine's
    /// `getFirstSubZone` ends the process; its `getSubZones` answers, and the first sub-zone
    /// here is the first of that list.
    #[test]
    fn aperture7_first_sub_zone_where_the_engines_call_ends_the_process() {
        for name in APERTURE_7 {
            let g = get_grid(name).unwrap();
            for (zone, first, last) in [("110", "11022", "11166"), ("11000", "1100022", "1100166")]
            {
                let z = g.zone_from_text(zone).unwrap();
                assert_eq!(g.text_id(g.first_sub_zone(z, 2).unwrap()), first, "{name}");
                let order = g.sub_zones(z, 2).unwrap();
                assert_eq!(order.len(), 46, "{name}");
                assert_eq!(g.first_sub_zone(z, 2), Ok(order[0]), "{name}");
                assert_eq!(g.sub_zone_at_index(z, 2, 0), Ok(order[0]), "{name}");
                assert_eq!(g.text_id(order[45]), last, "{name}");
            }
        }
    }

    /// In the broken seams the engine's quantiser gives some centroids of an order no zone:
    /// the order of `010004000400` at depth 5 holds the engine's null zone at 99 of its 17,053
    /// places. Each is the null zone here, at the engine's place, and no error: the length of
    /// the order and the place of every other entry are the engine's.
    #[test]
    fn aperture7_seam_order_holds_the_null_zone_at_the_engines_places() {
        #[rustfmt::skip]
        const NULL_PLACES: [usize; 99] = [
            204, 316, 380, 523, 603, 774, 862, 1_040, 1_131, 1_316, 1_410, 1_601, 1_699, 1_895, 1_996, 2_200, 2_303, 2_513, 2_620, 2_836, 2_946, 3_168, 3_282, 3_510, 3_626, 3_862, 3_981, 4_222, 4_345, 4_592, 4_718, 4_972, 5_101, 5_361, 5_493, 5_760, 5_895, 6_167, 6_306, 6_583, 6_722, 6_999, 7_138, 7_415, 7_554, 7_831, 7_970, 8_247, 8_386, 8_663, 8_802, 9_079, 9_218, 9_495, 9_634, 9_911, 10_050, 10_327, 10_466, 10_743, 10_882, 11_154, 11_289, 11_556, 11_688, 11_948, 12_077, 12_331, 12_457, 12_704, 12_827, 13_068, 13_187, 13_423, 13_539, 13_767, 13_881, 14_103, 14_213, 14_429, 14_536, 14_746, 14_849, 15_053, 15_154, 15_350, 15_448, 15_639, 15_733, 15_918, 16_009, 16_187, 16_275, 16_446, 16_526, 16_669, 16_733, 16_845, 16_893
        ];
        for name in APERTURE_7 {
            let g = get_grid(name).unwrap();
            let z = g.zone_from_text("010004000400").unwrap();
            assert_eq!(g.count_sub_zones(z, 5), Ok(17_053), "{name}");
            let order = g.sub_zones(z, 5).unwrap();
            assert_eq!(order.len(), 17_053, "{name}");
            let nulls: Vec<usize> = (0..order.len())
                .filter(|&i| order[i] == ZoneId::NULL)
                .collect();
            assert_eq!(nulls, NULL_PLACES, "{name}");
            // The whole of the engine's list, the null zone at its places, folded into one
            // number (FNV-1a over the identifiers).
            let folded = order.iter().fold(0xcbf2_9ce4_8422_2325_u64, |h, s| {
                (h ^ s.0).wrapping_mul(0x0100_0000_01b3)
            });
            assert_eq!(folded, 0x5f3f_61ea_e636_a552, "{name}");
            assert_eq!(g.text_id(order[0]), "01000400040346464", "{name}");
            assert_eq!(g.text_id(order[203]), "01000400040032330", "{name}");
            assert_eq!(g.text_id(order[8_526]), "01000400040000000", "{name}");
            assert_eq!(g.text_id(order[17_052]), "01000400040431313", "{name}");
            // The entry at such a place is the null zone, and no error.
            assert_eq!(g.sub_zone_at_index(z, 5, 204), Ok(ZoneId::NULL), "{name}");
            assert_eq!(g.text_id(order[204]), crate::NULL_TEXT, "{name}");
        }
        let g = igeo7();
        let zone = g.zone_from_text("010004000400").unwrap();
        let order = zone.sub_zones(5).unwrap();
        assert_eq!(order.len(), 17_053);
        assert_eq!(order[204].id(), ZoneId::NULL);
        assert_eq!(order[204].text_id(), crate::NULL_TEXT);
        assert_eq!(zone.sub_zone_at_index(5, 204).unwrap().id(), ZoneId::NULL);
    }

    /// The zone at an index is found without building the order, at any index of any depth to
    /// the finest level: the last of the 678,223,896,391 sub-zones of a hexagon at depth 14,
    /// which is DGGAL's own answer, and the last of the 9,499,079,489,284,316 of a pentagon of
    /// level 0 at depth 19, the greatest index that there is.
    #[test]
    fn aperture7_sub_zone_at_index_reaches_the_last_index_without_materialising() {
        for name in APERTURE_7 {
            let g = get_grid(name).unwrap();
            let z = g.zone_from_text("0064156").unwrap();
            let count = g.count_sub_zones(z, 14).unwrap();
            assert_eq!(count, 678_223_896_391, "{name}");
            for (index, sub) in [
                (1, "006415613131313131313"),
                (339_111_948_195, "006415600000000000000"),
                (678_223_896_390, "006415664646464646466"),
            ] {
                assert_eq!(
                    g.text_id(g.sub_zone_at_index(z, 14, index).unwrap()),
                    sub,
                    "{name}, index {index}"
                );
            }
            assert_eq!(
                g.sub_zone_at_index(z, 14, count),
                Err(Error::IndexOutOfRange {
                    index: count,
                    count
                }),
                "{name}"
            );
            // The list itself is refused, above the limit on a materialised list.
            assert_eq!(
                g.sub_zones(z, 14),
                Err(Error::TooManySubZones {
                    count,
                    limit: max_materialised_sub_zones()
                }),
                "{name}"
            );

            // The first of the greatest order, at depth 19 below a pentagon of level 0.
            let z = g.zone_from_text("00").unwrap();
            assert_eq!(
                g.count_sub_zones(z, 19),
                Ok(9_499_079_489_284_316),
                "{name}"
            );
            assert_eq!(
                g.text_id(g.first_sub_zone(z, 19).unwrap()),
                "013131313131313131313",
                "{name}"
            );
        }

        // The last of it, on one grid: the three share the topology, and the search reckons
        // the 134 million scanlines before the last. No engine answers this one: its own search
        // by index is at fault at a pentagon at an odd depth (it answers
        // `034503320445033204454` here, a zone of another base cell), and its list is answered
        // to depth 10 alone. The answer is this crate's own, by the search that agrees with
        // the engine's list at every index of every order compared; at depths 1, 3, 5 and 7
        // the last of the engine's list is `043`, `03622`, `0066201` and `006620131`, as here.
        // The engine's `zoneHasSubZone` holds the zone answered here for a sub-zone of `00`,
        // and not the zone that its own search answers.
        let g = get_grid("IGEO7").unwrap();
        let z = g.zone_from_text("00").unwrap();
        for (depth, last) in [
            (1, "043"),
            (3, "03622"),
            (5, "0066201"),
            (7, "006620131"),
            (19, "006620115066201156311"),
        ] {
            let count = g.count_sub_zones(z, depth).unwrap();
            let sub = g.sub_zone_at_index(z, depth, count - 1).unwrap();
            assert_eq!(g.text_id(sub), last, "depth {depth}");
            assert_eq!(g.resolution(sub), depth, "depth {depth}");
        }
        assert_eq!(
            g.sub_zone_at_index(z, 19, 9_499_079_489_284_316),
            Err(Error::IndexOutOfRange {
                index: 9_499_079_489_284_316,
                count: 9_499_079_489_284_316
            })
        );
    }

    /// The limit on a materialised list holds on aperture 7 as on aperture 3: a hexagon's
    /// 5,767,201 sub-zones at depth 8 and a pentagon's 4,806,001 are refused, while the count,
    /// the first sub-zone and the zone at an index are answered at any depth.
    #[test]
    fn aperture7_sub_zones_above_the_limit_are_refused() {
        let g = igeo7();
        for (zone, count) in [("0064156", 5_767_201), ("00", 4_806_001)] {
            let z = g.zone_from_text(zone).unwrap().id();
            assert_eq!(g.count_sub_zones(z, 8), Ok(count));
            assert!(count > MAX_MATERIALISED_SUB_ZONES);
            assert_eq!(
                g.sub_zones(z, 8),
                Err(Error::TooManySubZones {
                    count,
                    limit: max_materialised_sub_zones()
                })
            );
            assert!(g.first_sub_zone(z, 8).is_ok());
            assert!(g.sub_zone_at_index(z, 8, count - 1).is_ok());
        }
        // A depth that would pass the finest level is refused before anything is counted.
        let z = g.zone_from_text("0064156").unwrap().id();
        for depth in [15, 16, 255] {
            assert!(matches!(
                g.count_sub_zones(z, depth),
                Err(Error::InvalidZone(_))
            ));
            assert!(matches!(g.sub_zones(z, depth), Err(Error::InvalidZone(_))));
            assert!(matches!(
                g.first_sub_zone(z, depth),
                Err(Error::InvalidZone(_))
            ));
            assert!(matches!(
                g.sub_zone_at_index(z, depth, 1),
                Err(Error::InvalidZone(_))
            ));
        }
    }

    /// The index of a sub-zone on aperture 7 is the engine's own, `getSubZoneIndex`, read from
    /// the live engine at every row; where the engine answers -1 the answer is `None`.
    #[test]
    fn aperture7_sub_zone_index_is_the_engines_walk() {
        for name in APERTURE_7 {
            let g = get_grid(name).unwrap();
            let id = |text: &str| g.zone_from_text(text).unwrap();
            for (zone, sub, index) in [
                ("0064156", "00641524", Some(0)),
                ("0064156", "00641565", Some(5)),
                ("0064156", "00641560", Some(6)),
                // A sub-zone of the neighbour `0064154` alone.
                ("0064156", "00641545", None),
                // A pentagon at an odd depth, where the engine's `getSubZoneAtIndex` answers
                // `030` and `034` at the indices 9 and 10, and its walk answers these.
                ("00", "006", Some(9)),
                ("00", "043", Some(10)),
                ("00", "030", None),
                // The pentagon of the south pole.
                ("11", "110", Some(6)),
                ("110", "11022", Some(0)),
                ("110", "11000", Some(27)),
                ("11000", "1100022", Some(0)),
                // The middle of an order in a broken seam.
                ("010004000400", "01000400040000000", Some(8_526)),
            ] {
                assert_eq!(
                    g.sub_zone_index(id(zone), id(sub)),
                    Ok(index),
                    "{name}: {sub} in {zone}"
                );
                if let Some(i) = index {
                    let depth = (sub.len() - zone.len()) as u8;
                    assert_eq!(
                        g.sub_zone_at_index(id(zone), depth, i),
                        Ok(id(sub)),
                        "{name}: {zone} at index {i}"
                    );
                }
            }
        }
    }

    /// In the broken seams an order may name a zone that has no index. The engine's walk
    /// first asks whether the zone is a sub-zone at all, by a test that the seams defeat, and
    /// then seeks the zone's own centroid among those of the order, which need not hold it.
    /// Each row is the engine's `getSubZones` and its `getSubZoneIndex`, the same on the
    /// three grids.
    #[test]
    fn aperture7_seam_entries_may_have_no_index() {
        for name in APERTURE_7 {
            let g = get_grid(name).unwrap();
            let id = |text: &str| g.zone_from_text(text).unwrap();
            for (zone, entries) in [
                (
                    "0000000000000005",
                    &[
                        (2_usize, "00000000000000053", Some(2_u64)),
                        (3, "00000000000000052", Some(3)),
                        // Refused by the engine's test of a sub-zone.
                        (4, "00000000000000005", None),
                        (5, "00000000000000051", Some(5)),
                        (6, "00000000000000050", Some(6)),
                        (8, "00000000000000532", Some(8)),
                        // Accepted by that test, and met by no centroid of the order.
                        (9, "00000000000000055", None),
                    ][..],
                ),
                (
                    "0000000000000000136",
                    // The order names the first of these at the places 0, 3, 7 and 11.
                    &[
                        (0, "00000000000000000005", None),
                        (1, "00000000000000000515", None),
                        (8, "00000000000000003555", Some(8)),
                    ][..],
                ),
                (
                    "00052626050026015",
                    &[
                        (1, "000526260500260105", Some(1)),
                        (4, "000526260500260533", None),
                        (12, "000526260500261626", Some(12)),
                    ][..],
                ),
            ] {
                let order = g.sub_zones(id(zone), 1).unwrap();
                for &(place, sub, index) in entries {
                    assert_eq!(order[place], id(sub), "{name}: entry {place} of {zone}");
                    assert_eq!(
                        g.sub_zone_index(id(zone), id(sub)),
                        Ok(index),
                        "{name}: {sub} in {zone}"
                    );
                }
            }
        }
    }

    /// An index is answered only where the entry of the order there is the zone asked about.
    /// In the broken seams the engine's walk, which compares centroids, may answer an index at
    /// which its own order, which holds what its quantiser makes of them, names another zone:
    /// for `000000000000000001` below `00000000000000000` the engine's `getSubZoneIndex`
    /// answers 6, where its `getSubZones` holds `000000000000000000`, and the zone asked about
    /// stands at another place of that order. Such a zone has no index here.
    #[test]
    fn aperture7_sub_zone_index_is_confirmed_by_the_order() {
        for name in APERTURE_7 {
            let g = get_grid(name).unwrap();
            let zone = g.zone_from_text("00000000000000000").unwrap();
            let sub = g.zone_from_text("000000000000000001").unwrap();
            let order = g.sub_zones(zone, 1).unwrap();
            assert!(order.contains(&sub), "{name}");
            assert_eq!(g.text_id(order[6]), "000000000000000000", "{name}");
            assert_eq!(g.sub_zone_index(zone, sub), Ok(None), "{name}");
            // Its centroid child, which the order holds at 6, has that index.
            assert_eq!(g.sub_zone_index(zone, order[6]), Ok(Some(6)), "{name}");
        }
    }

    /// The index is found by walking the order from its head, without building it, and the
    /// walk is bounded as a list is: the centroid descendant in the middle of the 825,259
    /// sub-zones of a hexagon at depth 7 has its index, 412,629, which is the engine's; one
    /// level deeper the order is longer than a list may be, and the walk is refused before
    /// it begins. The limit is the one in force, which the environment of the run may have
    /// lowered: below 825,259 the index at depth 7 is refused as well.
    #[test]
    fn aperture7_sub_zone_index_is_bounded_by_the_limit_of_a_list() {
        let g = igeo7();
        let limit = max_materialised_sub_zones();
        let z = g.zone_from_text("0064156").unwrap().id();
        let middle = g.zone_from_text("00641560000000").unwrap().id();
        assert_eq!(g.count_sub_zones(z, 7), Ok(825_259));
        assert_eq!(
            g.sub_zone_index(z, middle),
            if limit >= 825_259 {
                Ok(Some(412_629))
            } else {
                Err(Error::TooManySubZones {
                    count: 825_259,
                    limit,
                })
            }
        );
        assert_eq!(g.sub_zone_at_index(z, 7, 412_629), Ok(middle));
        for (zone, sub, count) in [
            ("0064156", "006415600000000", 5_767_201),
            ("00", "0000000000", 4_806_001),
            ("00", "000000000000000000000", 9_499_079_489_284_316),
        ] {
            let (z, sub) = (
                g.zone_from_text(zone).unwrap().id(),
                g.zone_from_text(sub).unwrap().id(),
            );
            assert!(count > MAX_MATERIALISED_SUB_ZONES);
            assert_eq!(
                g.sub_zone_index(z, sub),
                Err(Error::TooManySubZones { count, limit }),
                "{zone}"
            );
        }
    }

    /// The disk of radius 1 is the zone and then its neighbours, in the order
    /// `neighbors` gives them; of radius 0, the zone alone. The grid's own walk
    /// yields the identifiers of the zone's walk.
    #[test]
    fn disk_k1_is_self_plus_neighbours() {
        let g = igeo7();
        let z = g.zone_from_text("0064156").unwrap();
        let d: Vec<ZoneId> = z.disk(1).map(|z| z.id()).collect();
        let mut want = vec![z.id()];
        want.extend(g.neighbors(z.id()));
        assert_eq!(d, want);
        assert_eq!(g.disk(z.id(), 1).collect::<Vec<_>>(), want);
        assert_eq!(z.disk(0).map(|z| z.id()).collect::<Vec<_>>(), [z.id()]);
    }
    #[test]
    fn get_grid_by_name() {
        assert_eq!(get_grid("IGEO7").unwrap().name(), "IGEO7");
        assert!(matches!(get_grid("NOPE"), Err(Error::UnknownGrid(_))));
    }

    /// A fabricated identifier, and DGGAL's own null sentinel, must answer as
    /// the null zone everywhere rather than panicking inside the topology's
    /// fixed-size tables, which is what an unvalidated base cell of 12 to 15
    /// would do.
    #[test]
    fn fabricated_identifiers_answer_as_the_null_zone() {
        let g = igeo7();
        for id in [BASE_12, ZoneId::NULL] {
            assert_eq!(g.centroid(id), GeoPoint { lat: 0.0, lon: 0.0 });
            assert_eq!(g.vertices(id), vec![GeoPoint { lat: 0.0, lon: 0.0 }; 6]);
            assert!(g.neighbors(id).is_empty());
            assert!(g.parents(id).is_empty());
            assert!(g.children(id).is_empty());
            assert_eq!(g.centroid_parent(id), None);
            assert!(!g.is_centroid_child(id));
            assert_eq!(g.zone(id).disk(1).count(), 1);
            assert_eq!(g.disk(id, u32::MAX).collect::<Vec<_>>(), [id]);
        }
    }

    /// The sub-zone methods return a `Result`, so they name the failure
    /// instead of answering as the null zone, and they do so before the
    /// grid's own lack of a sub-zone order, which is a question about a real
    /// zone.
    #[test]
    fn fabricated_identifiers_are_reported_by_the_sub_zone_methods() {
        let g = igeo7();
        assert!(matches!(
            g.count_sub_zones(BASE_12, 1),
            Err(Error::InvalidZone(_))
        ));
        assert!(matches!(
            g.sub_zones(BASE_12, 0),
            Err(Error::InvalidZone(_))
        ));
        assert!(matches!(
            g.first_sub_zone(ZoneId::NULL, 1),
            Err(Error::InvalidZone(_))
        ));
        assert!(matches!(
            g.sub_zone_at_index(BASE_12, 1, 0),
            Err(Error::InvalidZone(_))
        ));
    }

    /// A zone with no geometry is not the same thing as an identifier with no
    /// zone behind it: the level-20 address is one this topology recognises,
    /// and, bound as a raw identifier, it keeps its resolution and its text id
    /// even though DGGAL gives it no geometry. It has no parents and no children,
    /// as the engine answers for it.
    ///
    /// Its text no longer leads back to it. DGGAL reads the twenty-two characters
    /// as that level-20 zone, and so names a zone without geometry from text; this
    /// crate refuses the text instead, as longer than the finest resolution allows.
    /// The round trip from identifier to text and back is thereby broken for this
    /// one level, on purpose, and nothing else in the crate relies on it.
    #[test]
    fn null_geometry_is_not_an_invalid_identifier() {
        let g = igeo7();
        let id =
            <crate::indexings::Z7 as crate::Indexing>::encode(&crate::Address::new(1, &[1; 20]));
        assert_eq!(g.resolution(id), 20);
        assert!(g.parents(id).is_empty());
        assert!(g.children(id).is_empty());
        assert_eq!(g.centroid_parent(id), None);
        assert!(g.zone(id).parent().is_none());
        let text = g.text_id(id);
        assert_eq!(text, format!("01{}", "1".repeat(20)));
        assert!(matches!(
            g.zone_from_text(&text),
            Err(Error::InvalidZone(_))
        ));
        // It has no sub-zones: any depth below it lies past the finest level, the depth of
        // the zone itself among them.
        for depth in [0, 1] {
            assert!(matches!(
                g.count_sub_zones(id, depth),
                Err(Error::InvalidZone(_))
            ));
            assert!(matches!(g.sub_zones(id, depth), Err(Error::InvalidZone(_))));
        }
        // Validated, it is a zone, and it is its own sub-zone at index 0, as the
        // engine also answers.
        assert_eq!(g.sub_zone_index(id, id), Ok(Some(0)));
    }

    /// A NaN or infinite coordinate is refused, naming the coordinate at fault,
    /// where the engine answers the pentagon of base cell 1 whatever the input.
    /// The resolution is checked before either coordinate, and the latitude before
    /// the longitude.
    #[test]
    fn a_non_finite_coordinate_is_refused() {
        let g = igeo7();
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert_eq!(
                g.zone_from_geo(bad, 0.0, 5).err(),
                Some(Error::NonFinite { name: "latitude" })
            );
            assert_eq!(
                g.zone_from_geo(0.0, bad, 5).err(),
                Some(Error::NonFinite { name: "longitude" })
            );
            assert_eq!(
                g.zone_from_geo(bad, bad, 5).err(),
                Some(Error::NonFinite { name: "latitude" })
            );
            assert_eq!(
                g.zone_from_geo(bad, 0.0, 20).err(),
                Some(Error::ResolutionOutOfRange { res: 20, max: 19 })
            );
        }
        // The extreme finite coordinates are still answered, each with a real cell
        // and not the pentagon a non-finite one would have been given.
        for (lat, lon) in [(f64::MAX, 0.0), (0.0, f64::MIN), (f64::MIN_POSITIVE, 0.0)] {
            let z = g.zone_from_geo(lat, lon, 5).unwrap();
            assert_ne!(z.id(), ZoneId::NULL, "({lat:e}, {lon:e})");
            assert_ne!(z.text_id(), "0100000", "({lat:e}, {lon:e})");
        }
    }

    /// A configuration with a NaN or infinite angle is refused, naming the field,
    /// rather than built into a grid that answers NaN.
    #[test]
    fn a_non_finite_configuration_is_refused() {
        use crate::GridConfig;
        use crate::registry::Igeo7;
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let c = GridConfig {
                orientation_lat_deg: bad,
                ..GridConfig::default()
            };
            assert_eq!(
                Igeo7::new(c, "IGEO7").err(),
                Some(Error::NonFinite {
                    name: "orientation_lat_deg"
                })
            );
            let c = GridConfig {
                orientation_lon_deg: bad,
                ..GridConfig::default()
            };
            assert_eq!(
                Igeo7::new(c, "IGEO7").err(),
                Some(Error::NonFinite {
                    name: "orientation_lon_deg"
                })
            );
            let c = GridConfig {
                azimuth_deg: bad,
                ..GridConfig::default()
            };
            assert_eq!(
                Igeo7::new(c, "IGEO7").err(),
                Some(Error::NonFinite {
                    name: "azimuth_deg"
                })
            );
        }
        // A finite configuration, the canonical one included, builds the same grid
        // as the registry's.
        let g = Igeo7::new(GridConfig::default(), "IGEO7").unwrap();
        let ours = g.zone_from_geo(38.7223, -9.1393, 5).unwrap().id();
        assert_eq!(
            ours,
            igeo7().zone_from_geo(38.7223, -9.1393, 5).unwrap().id()
        );
    }

    /// `sub_zone_index` validates both identifiers before it answers, as the other
    /// sub-zone methods do; the engine answers 0 for each of these within itself.
    #[test]
    fn sub_zone_index_refuses_an_identifier_with_no_cell() {
        let g = igeo7();
        let real = g.zone_from_text("0064156").unwrap().id();
        for bad in [BASE_12, ZoneId::NULL] {
            assert!(matches!(
                g.sub_zone_index(bad, bad),
                Err(Error::InvalidZone(_))
            ));
            assert!(matches!(
                g.sub_zone_index(real, bad),
                Err(Error::InvalidZone(_))
            ));
            assert!(matches!(
                g.sub_zone_index(bad, real),
                Err(Error::InvalidZone(_))
            ));
        }
        assert_eq!(g.sub_zone_index(real, real), Ok(Some(0)));
    }

    /// A point the engine places in no cell comes back as the null zone, and
    /// the call succeeds: DGGAL returns its `nullZone` from
    /// `getZoneFromWGS84Centroid` rather than failing, so this does too.
    ///
    /// The coordinates are a point 0.0009 degrees from the north pole, verified
    /// against live DGGAL, which answers `(null)` at resolution 17 and an
    /// ordinary zone at 16.
    #[test]
    fn a_point_in_no_cell_is_the_null_zone_and_not_an_error() {
        let g = igeo7();
        let z = g
            .zone_from_geo(89.999_129_147_255_16, 9.333_761_199_597_859, 17)
            .expect("the null zone is a successful answer, not an error");
        assert_eq!(z.id(), ZoneId::NULL);
        assert_eq!(z.text_id(), crate::NULL_TEXT);
        // The same point one resolution coarser is an ordinary zone.
        let ok = g
            .zone_from_geo(89.999_129_147_255_16, 9.333_761_199_597_859, 16)
            .unwrap();
        assert_ne!(ok.id(), ZoneId::NULL);
    }

    /// The null zone's geometry, hierarchy, neighbours, pentagon flag and text
    /// match DGGAL; its resolution does not, and that one is named here.
    ///
    /// Five of these answers, the centroid, the vertices, the neighbours, the
    /// parents and the children, run through `address`, which rejects the
    /// sentinel's decoded base cell of fifteen, so they fall out of that guard
    /// rather than from any special case. The other two do not: `Grid::is_pentagon`
    /// and `Grid::text_id` never call `address`, and their answers, false and
    /// `(null)`, are deliberate special cases for the sentinel in `Z7::is_pentagon`
    /// and `Z7::to_text`. All are asserted rather than assumed. Each expected value
    /// was read from the live engine for its own null zone: the zero centroid, six
    /// zero vertices, and no neighbours, parents or children.
    ///
    /// The one divergence this does not claim to cover is `resolution`. DGGAL
    /// reports -1 for the sentinel and this port reports 0, deliberately, for the
    /// reasons given at [`Grid::resolution`]. The assertion at the end pins the
    /// value this port gives, so the difference is recorded rather than hidden.
    #[test]
    fn the_null_zone_matches_dggal_but_for_its_resolution() {
        let g = igeo7();
        let z = g.zone(ZoneId::NULL);
        assert_eq!(z.centroid(), GeoPoint { lat: 0.0, lon: 0.0 });
        assert_eq!(z.vertices().len(), 6);
        assert!(
            z.vertices()
                .iter()
                .all(|v| *v == GeoPoint { lat: 0.0, lon: 0.0 })
        );
        assert!(z.neighbors().is_empty());
        assert!(z.parents().is_empty());
        assert!(z.children().is_empty());
        assert_eq!(z.text_id(), crate::NULL_TEXT);
        // And it survives a trip through its own text form.
        assert_eq!(
            g.zone_from_text(crate::NULL_TEXT).unwrap().id(),
            ZoneId::NULL
        );
        // DGGAL's `nPoints` is 0 for the sentinel, so it is neither a pentagon nor
        // a hexagon; what matters to a caller is that it is not claimed to be one.
        assert!(!z.is_pentagon());
        // The known divergence, pinned: DGGAL says -1 here.
        assert_eq!(z.resolution(), 0);
    }

    /// Every neighbour of these zones has the zone back as a neighbour, which
    /// the neighbour search does not enforce and would break if a step across an
    /// edge landed one cell too far.
    ///
    /// It does not hold everywhere. At resolutions 15 to 19, along the edges from
    /// base cell 0 to 1 and from 1 to 6, where the engine's odd resolutions leave a
    /// strip with no cells, the engine's own lists, and so these, may name a
    /// neighbour that does not name the zone back. That is why [`crate::Disk`]
    /// remembers every zone it has yielded rather than the last two rings alone.
    #[test]
    fn neighbours_are_reciprocal() {
        let g = igeo7();
        for text in ["0064156", "0000000", "0611111", "05"] {
            let z = g.zone_from_text(text).unwrap();
            for n in z.neighbors() {
                assert!(
                    n.neighbors().contains(&z),
                    "{} is not a neighbour of its neighbour {}",
                    text,
                    n.text_id()
                );
            }
        }
    }

    /// `sub_zone_at_index` agrees with the materialised order at every index, for a sample of
    /// zones and depths small enough to materialise safely.
    #[test]
    fn sub_zone_at_index_matches_the_materialised_order() {
        let g = isea3h();
        for text in ["B6-5-A", "A2-0-A", "A6-0-C", "B6-1-B", "B5-3-B"] {
            let z = a3_id(&g, text);
            for depth in 1..=3 {
                let subs = g.sub_zones(z, depth).unwrap();
                for (i, &sub) in subs.iter().enumerate() {
                    assert_eq!(
                        g.sub_zone_at_index(z, depth, i as u64),
                        Ok(sub),
                        "{text} at depth {depth}, index {i}"
                    );
                }
            }
        }
    }

    /// `sub_zone_at_index` reaches resolution 33 from resolution 0, through the generators' own
    /// index skip rather than by materialising 4.6e15 sub-zones, for the first, the last and a
    /// middle index. The slowest takes about 0.3 s unoptimised on a laptop and 1.1 s on a
    /// shared CI runner running the other tests beside it; the bound of ten seconds still
    /// separates the skip from the materialisation, which would never finish.
    #[test]
    fn sub_zone_at_index_reaches_depth_33_from_resolution_0_promptly() {
        let g = isea3h();
        let z = a3_id(&g, "A0-0-A");
        let count = g.count_sub_zones(z, 33).unwrap();
        assert_eq!(count, 4_632_550_579_746_406);
        for index in [0, count / 2, count - 1] {
            let start = std::time::Instant::now();
            let sub = g.sub_zone_at_index(z, 33, index).unwrap();
            let elapsed = start.elapsed();
            assert!(
                elapsed < std::time::Duration::from_secs(10),
                "index {index} took {elapsed:?}"
            );
            assert_eq!(g.resolution(sub), 33);
        }
    }

    /// An index at or beyond the count is `Error::IndexOutOfRange`, never a panic and never the
    /// materialised fallback's own out-of-bounds.
    #[test]
    fn sub_zone_at_index_beyond_the_count_is_index_out_of_range() {
        let g = isea3h();
        let z = a3_id(&g, "A0-0-A");
        let count = g.count_sub_zones(z, 1).unwrap();
        assert_eq!(count, 6);
        for index in [count, count + 1000] {
            assert_eq!(
                g.sub_zone_at_index(z, 1, index),
                Err(Error::IndexOutOfRange { index, count })
            );
        }
    }

    /// `sub_zone_index` from resolution 0 of a resolution-33 sub-zone refuses the search rather
    /// than run it: the sub-zone's own index is found first, promptly, through
    /// `sub_zone_at_index`'s index skip, but asking `sub_zone_index` to find it back searches by
    /// materialising, which the cap of [`max_materialised_sub_zones`] refuses at that depth, within a
    /// second, as `Grid::sub_zone_index` keeps its materialising search and its cap.
    #[test]
    fn sub_zone_index_of_a_resolution_33_sub_zone_is_too_many_sub_zones() {
        let g = isea3h();
        let z = a3_id(&g, "A0-0-A");
        let count = g.count_sub_zones(z, 33).unwrap();
        let sub = g.sub_zone_at_index(z, 33, 0).unwrap();
        assert_eq!(g.resolution(sub), 33);
        let start = std::time::Instant::now();
        assert_eq!(
            g.sub_zone_index(z, sub),
            Err(Error::TooManySubZones {
                count,
                limit: max_materialised_sub_zones(),
            })
        );
        assert!(start.elapsed() < std::time::Duration::from_secs(1));
    }

    /// The round trip `sub_zone_index(z, sub_zone_at_index(z, d, i)) == Some(i)` over a sample of
    /// zones at depths 1 to 6, reaching resolution 33 three ways: deep (depth 6 from resolution
    /// 27), shallow (depth 1 from resolution 32) and in between (depth 2 from resolution 31), the
    /// index this crate finds each time where the engine's own `zoneHasSubZone` would refuse.
    #[test]
    fn sub_zone_index_and_at_index_round_trip() {
        let g = isea3h();
        let a0 = a3_id(&g, "A0-0-A");
        let cases = [
            (a3_id(&g, "B6-5-A"), 6),
            (a3_id(&g, "A2-0-A"), 6),
            (a3_id(&g, "A6-0-C"), 6),
            (g.first_sub_zone(a0, 27).unwrap(), 6),
            (g.first_sub_zone(a0, 31).unwrap(), 2),
            (g.first_sub_zone(a0, 32).unwrap(), 1),
        ];
        for (z, max_depth) in cases {
            assert!(g.resolution(z) + max_depth <= 33);
            for depth in 1..=max_depth {
                let count = g.count_sub_zones(z, depth).unwrap();
                for index in 0..count {
                    let sub = g.sub_zone_at_index(z, depth, index).unwrap();
                    assert_eq!(
                        g.sub_zone_index(z, sub),
                        Ok(Some(index)),
                        "{z:?} at depth {depth}, index {index}"
                    );
                }
            }
        }
    }

    // --- the refined boundary and the extent ---

    /// A Z7 identifier of twenty digits on the pentagon path, `0000000000000000000000`, and one
    /// on a hexagon path, `0312345601234560123456`: the packing holds them and the engine gives
    /// them no geometry.
    const LEVEL_20: [ZoneId; 2] = [ZoneId(0), ZoneId(0x329c_b814_e5c0_a72e)];

    /// `01222222222222222220`, which the engine writes for a neighbour and cannot read back,
    /// and to which it gives an ordinary hexagon's ring and extent.
    const UNREADABLE: ZoneId = ZoneId(0x1492_4924_9249_243f);

    /// What DGGAL v0.0.6, BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, answers of one
    /// zone: its `getZoneRefinedWGS84Vertices` and `getZoneWGS84Extent`, in radians as the
    /// engine gives them, written as the shortest decimals that read back to the same doubles.
    struct EngineZone {
        text: &'static str,
        id: u64,
        /// The count of points of the ring at the automatic refinement.
        points_at_0: usize,
        /// `[ll.lat, ll.lon, ur.lat, ur.lon]`.
        extent: [f64; 4],
        /// The ring at a refinement of 1: latitude, then longitude.
        ring_at_1: &'static [(f64, f64)],
    }

    /// Where a ring of this crate leaves the engine's, as a sequence of doubles compared by
    /// their bits: `None` where the two are the same.
    fn ring_difference(ours: &[(f64, f64)], theirs: &[(f64, f64)]) -> Option<String> {
        if ours.len() != theirs.len() {
            return Some(format!(
                "{} points, the engine has {}",
                ours.len(),
                theirs.len()
            ));
        }
        ours.iter()
            .zip(theirs)
            .position(|(a, b)| a.0.to_bits() != b.0.to_bits() || a.1.to_bits() != b.1.to_bits())
            .map(|i| {
                format!(
                    "point {i} of {} is {:?}, the engine has {:?}",
                    ours.len(),
                    ours[i],
                    theirs[i]
                )
            })
    }

    /// A point of the engine, in radians, as the public methods give it: in degrees, by the
    /// product with `180 / Pi`.
    fn degrees((lat, lon): (f64, f64)) -> GeoPoint {
        GeoPoint {
            lat: lat * math::RAD2DEG,
            lon: lon * math::RAD2DEG,
        }
    }

    /// Every way in which this crate's ring at a refinement of 1, its count of points at the
    /// automatic refinement and its extent leave the engine's at the zones of `engine`, in
    /// radians and then through the public methods, in degrees.
    fn differences_from_the_engine<P: Projection, T: Topology, I: Indexing>(
        grid: &Grid<P, T, I>,
        engine: &[EngineZone],
    ) -> Vec<String> {
        let mut differences = Vec::new();
        for zone in engine {
            let id = ZoneId(zone.id);
            assert_eq!(grid.text_id(id), zone.text);
            let name = format!("{} {}", grid.name(), zone.text);
            let a = grid.geometry_address(id).expect("a zone with geometry");

            let ring = grid.refined_ring_radians(id, &a, 1);
            if let Some(difference) = ring_difference(&ring, zone.ring_at_1) {
                differences.push(format!("{name} at a refinement of 1: {difference}"));
            }
            let points = grid.refined_ring_radians(id, &a, 0).len();
            if points != zone.points_at_0 {
                differences.push(format!(
                    "{name} at the automatic refinement: {points} points, the engine has {}",
                    zone.points_at_0
                ));
            }
            let extent = grid.extent_radians(id, &a);
            if extent.map(|e| e.map(f64::to_bits)) != Some(zone.extent.map(f64::to_bits)) {
                differences.push(format!(
                    "{name}: the extent is {extent:?}, the engine has {:?}",
                    zone.extent
                ));
            }

            let public: Vec<GeoPoint> = zone.ring_at_1.iter().copied().map(degrees).collect();
            if grid.refined_vertices_within_limit(id, 1) != public {
                differences.push(format!("{name}: the public ring is not the engine's"));
            }
            let [ll_lat, ll_lon, ur_lat, ur_lon] = zone.extent;
            let public = Extent {
                ll: degrees((ll_lat, ll_lon)),
                ur: degrees((ur_lat, ur_lon)),
            };
            if grid.extent(id) != Some(public) || grid.zone(id).extent() != Some(public) {
                differences.push(format!("{name}: the public extent is not the engine's"));
            }
        }
        differences
    }

    /// The differences of [`differences_from_the_engine`] on the three aperture-7 grids, at the
    /// zones of each grid's table that `zones` selects.
    fn differences_on_the_three_grids(
        zones: impl Fn(&'static [EngineZone; 8]) -> &'static [EngineZone],
    ) -> Vec<String> {
        let mut differences = differences_from_the_engine(igeo7(), zones(&IGEO7_ENGINE));
        differences.extend(differences_from_the_engine(ivea7h(), zones(&IVEA7H_ENGINE)));
        differences.extend(differences_from_the_engine(rtea7h(), zones(&RTEA7H_ENGINE)));
        differences
    }

    /// Site GR-1, the sum that follows `wrapLonAt` in `getRefinedVertices`, at two of the zones
    /// that pin it on each of the three grids: the hexagon `0010` and `0003`, a hexagon across
    /// an interruption, both of level 2. Put back in the eC's order, `(w + centroid.lon) - D`,
    /// the site moves a longitude of each of these six rings in its last bits, and with it an
    /// end of each extent; it was so put back, and this test seen to fail at every one.
    #[test]
    fn the_sum_after_wrap_lon_at_follows_the_compiled_order() {
        let differences = differences_on_the_three_grids(|engine| &engine[..2]);
        assert!(
            differences.is_empty(),
            "{} answers leave the engine's:\n{}",
            differences.len(),
            differences.join("\n")
        );
    }

    /// One zone of each kind, on the three grids, against the engine's ring at a refinement of
    /// 1, its count of points at the automatic refinement and its extent: a hexagon, a pentagon,
    /// a pentagon at a pole of the icosahedron, a hexagon across an interruption, a hexagon
    /// that touches the north pole, and the zone centred on that pole at resolution 19, whose
    /// ring is the pole alone.
    #[test]
    fn a_zone_of_each_kind_has_the_engines_ring_and_extent() {
        let differences = differences_on_the_three_grids(|engine| &engine[2..]);
        assert!(
            differences.is_empty(),
            "{} answers leave the engine's:\n{}",
            differences.len(),
            differences.join("\n")
        );

        // What the crate promises of these zones, on the public answers.
        for grid in [get_grid("IGEO7"), get_grid("IVEA7H"), get_grid("RTEA7H")] {
            let grid = grid.unwrap();
            // The zone that touches the pole reaches a latitude of exactly 90, and its extent
            // passes the antimeridian: the western end is the greater.
            let polar = grid.zone_from_text("0132323").unwrap();
            let extent = grid.extent(polar).unwrap();
            assert_eq!(extent.ur.lat, 90.0);
            assert!(extent.ll.lon > extent.ur.lon, "{extent:?}");
            let ring = grid.refined_vertices(polar, 0).unwrap();
            assert!(ring.iter().any(|p| p.lon > 180.0), "a continuous ring");
            // The zone centred on the pole: a ring of the pole alone, at any refinement, and an
            // extent that is that point.
            let pole = grid.zone_from_text("013232323232323232323").unwrap();
            for refinement in [0, 1, 2, 3] {
                let ring = grid.refined_vertices_within_limit(pole, refinement);
                assert_eq!(ring.len(), 1);
                assert_eq!(ring[0].lat, 90.0);
            }
            let extent = grid.extent(pole).unwrap();
            assert_eq!(extent.ll, extent.ur);
            assert_eq!(extent.ur.lat, 90.0);
            // A ring is never closed, and a hexagon's at the automatic refinement of its
            // resolution has twelve points to an edge.
            let hexagon = grid.zone_from_text("0313522").unwrap();
            let ring = grid.refined_vertices(hexagon, 0).unwrap();
            assert_eq!(ring.len(), 72);
            assert_ne!(ring.first(), ring.last());
        }
    }

    /// The differences of [`differences_from_the_engine`] on the three aperture-3 grids, at the
    /// zones of each grid's table that `zones` selects.
    fn differences_on_the_three_grids_of_aperture_3(
        zones: impl Fn(&'static [EngineZone; 15]) -> &'static [EngineZone],
    ) -> Vec<String> {
        let mut differences =
            differences_from_the_engine(crate::registry::isea3h(), zones(&ISEA3H_ENGINE));
        differences.extend(differences_from_the_engine(
            crate::registry::ivea3h(),
            zones(&IVEA3H_ENGINE),
        ));
        differences.extend(differences_from_the_engine(
            crate::registry::rtea3h(),
            zones(&RTEA3H_ENGINE),
        ));
        differences
    }

    /// Site GR-1 on aperture 3, where the loop of `getRefinedVertices` is `RI3H`'s own
    /// (`RI3H.ec:517`, the sum at `0xf79d` to `0xf7ab`), at two of the zones that pin it on
    /// each of the three grids: the hexagon `B4-4-A` and the pentagon `B4-0-A`, both of level
    /// 2. Put back in the eC's order the site moves a longitude of each of these six rings;
    /// it was so put back, and this test seen to fail at every one.
    #[test]
    fn the_sum_after_wrap_lon_at_follows_the_compiled_order_on_aperture_3() {
        let differences = differences_on_the_three_grids_of_aperture_3(|engine| &engine[..2]);
        assert!(
            differences.is_empty(),
            "{} answers leave the engine's:\n{}",
            differences.len(),
            differences.join("\n")
        );
    }

    /// One zone of each kind on the three aperture-3 grids, against the engine's ring at a
    /// refinement of 1, its count of points at the automatic refinement and its extent: a
    /// hexagon, a pentagon, a pentagon at either pole of the icosahedron, a hexagon across an
    /// interruption, a hexagon with an edge beside a pole, which is divided twenty times as
    /// finely, and seven zones that touch a pole and so pass through the pole branch. Five are
    /// on the north pole: `C0-4-C` of level 5; `F0-79-A` of level 10, whose ring names the pole
    /// twice; `M0-40DF8-A` of level 24, whose ring crosses itself; and of level 33 the zone
    /// centred on the pole, `Q0-1486BA0-C`, and `Q8-34AA7F2578AE0-C` beside it. The sixth is
    /// `A0-0-A` of level 0, an even level: there the two positions of the pole branch lie so
    /// near the pole's image that their inverses answer the pole itself, at a longitude that
    /// the `oddGrid` passed to them decides, and the engine passes true whatever the level.
    /// The last, `Q5-34AA7F2578AE1-A` of level 32, is on the south pole, and holds the two
    /// southern arms of the branch as the others hold the two northern: each of the four,
    /// with its two positions exchanged, was seen to fail this test.
    #[test]
    fn a_zone_of_each_kind_has_the_engines_ring_and_extent_on_aperture_3() {
        let differences = differences_on_the_three_grids_of_aperture_3(|engine| &engine[2..]);
        assert!(
            differences.is_empty(),
            "{} answers leave the engine's:\n{}",
            differences.len(),
            differences.join("\n")
        );

        // What the crate promises of these zones, on the public answers.
        for name in ["ISEA3H", "IVEA3H", "RTEA3H"] {
            let grid = get_grid(name).unwrap();
            // A zone that touches the pole reaches a latitude of exactly 90, and its ring
            // names the pole twice, as it is seen from either side; the two ends of its extent
            // lie half a turn apart.
            let polar = grid.zone_from_text("C0-4-C").unwrap();
            let extent = grid.extent(polar).unwrap();
            assert_eq!(extent.ur.lat, 90.0, "{name}");
            let span = (extent.ur.lon - extent.ll.lon).rem_euclid(360.0);
            assert!((span - 180.0).abs() < 1.0, "{name}: {extent:?}");
            let ring = grid.refined_vertices(polar, 0).unwrap();
            assert_eq!(ring.len(), 251, "{name}");
            assert_eq!(ring.iter().filter(|p| p.lat == 90.0).count(), 2, "{name}");
            // The zone centred on the pole at the finest level: its extent reaches the pole.
            let pole = grid.zone_from_text("Q0-1486BA0-C").unwrap();
            assert_eq!(grid.extent(pole).unwrap().ur.lat, 90.0, "{name}");
            // And a zone on the south pole reaches a latitude of exactly -90.
            let southern = grid.zone_from_text("Q5-34AA7F2578AE1-A").unwrap();
            assert_eq!(grid.extent(southern).unwrap().ll.lat, -90.0, "{name}");
            // A ring is never closed, and a hexagon's at the automatic refinement of level 5
            // has ten points to an edge.
            let hexagon = grid.zone_from_text("C4-49-B").unwrap();
            let ring = grid.refined_vertices(hexagon, 0).unwrap();
            assert_eq!(ring.len(), 60, "{name}");
            assert_ne!(ring.first(), ring.last(), "{name}");
        }
    }

    /// On aperture 3 the ring is the engine's sequence, which begins at the point
    /// `vertices()[1]` and passes the others in the order of [`Grid::vertices`]; at the two
    /// polar pentagons, whose plain ring is the engine's own, it begins at the engine's `v4`
    /// under root 10 and at its `v0` under root 11. The points are the same, and not the same
    /// doubles: a longitude of the refined ring may differ from the plain one in its last
    /// bits or by a turn, and on a sub-hex C or D two corners of the refined ring are drawn
    /// 2e-11 of a rhombus off the edge on which the plain ring has them (`RI3H.ec:1995-1996`,
    /// `:1999`, `:2004`).
    #[test]
    fn the_refined_ring_of_aperture_3_begins_at_the_second_plain_vertex() {
        // Whether `a` and `b` are one point of the sphere, to 1e-9 radians between them.
        let same = |a: GeoPoint, b: GeoPoint| {
            let unit = |p: GeoPoint| {
                let (lat, lon) = (p.lat * math::DEG2RAD, p.lon * math::DEG2RAD);
                [
                    math::cos(lat) * math::cos(lon),
                    math::cos(lat) * math::sin(lon),
                    math::sin(lat),
                ]
            };
            let (a, b) = (unit(a), unit(b));
            (0..3).all(|i| (a[i] - b[i]).abs() < 1e-9)
        };
        for name in ["ISEA3H", "IVEA3H", "RTEA3H"] {
            let grid = get_grid(name).unwrap();
            // A hexagon and a pentagon of either parity, sub-hexes C and D, a hexagon across
            // an interruption, and zones of the coarsest and of the finest levels.
            for text in [
                "C4-49-B",
                "C8-2C-A",
                "C0-0-A",
                "C0-0-B",
                "C0-3-C",
                "C5-1B-D",
                "C0-6-A",
                "A1-0-A",
                "A0-0-D",
                "Q4-32B90539F2412-B",
            ] {
                let id = grid.zone_from_text(text).unwrap();
                let plain = grid.vertices(id);
                let ring = grid.refined_vertices_within_limit(id, 1);
                assert_eq!(ring.len(), plain.len(), "{name} {text}");
                // The same points in the same order, from the second.
                for (i, &point) in ring.iter().enumerate() {
                    let vertex = plain[(i + 1) % plain.len()];
                    assert!(same(point, vertex), "{name} {text}: {point:?}, {vertex:?}");
                }
            }
            // The polar pentagons, at an even level and at an odd one.
            for (text, first) in [("CA-0-A", 4), ("CA-0-B", 4), ("CB-0-A", 0), ("CB-0-B", 0)] {
                let id = grid.zone_from_text(text).unwrap();
                let plain = grid.vertices(id);
                let ring = grid.refined_vertices_within_limit(id, 1);
                assert_eq!((ring.len(), plain.len()), (5, 5), "{name} {text}");
                for (i, &point) in ring.iter().enumerate() {
                    let vertex = plain[(first + i) % 5];
                    assert!(same(point, vertex), "{name} {text}: {point:?}, {vertex:?}");
                }
            }
        }
    }

    /// The thirty identifiers of a polar root with a non-zero index that the engine's own
    /// quantisation handed out for the points of a sample, at levels 2 to 15, on
    /// each of the three grids: `BA-2-A`, `BB-3-A`, `BB-6-A` and the rest. The
    /// engine gives a null array as the ring of each, and its parser refuses their text.
    const NULL_RING_CLASS: [u64; 30] = [
        0x0340_0000_0000_0008,
        0x0360_0000_0000_000c,
        0x0360_0000_0000_0018,
        0x0340_0000_0000_0009,
        0x0360_0000_0000_000d,
        0x0360_0000_0000_0019,
        0x0540_0000_0000_0014,
        0x0560_0000_0000_0048,
        0x0560_0000_0000_0120,
        0x0560_0000_0000_0024,
        0x0540_0000_0000_0015,
        0x0560_0000_0000_0049,
        0x0560_0000_0000_0121,
        0x0560_0000_0000_0025,
        0x0740_0000_0000_0038,
        0x0760_0000_0000_006c,
        0x0760_0000_0000_0af8,
        0x0760_0000_0000_00d8,
        0x0740_0000_0000_0039,
        0x0760_0000_0000_006d,
        0x0760_0000_0000_0af9,
        0x0760_0000_0000_00d9,
        0x0940_0000_0000_00a4,
        0x0940_0000_0000_00a5,
        0x0b40_0000_0000_01e8,
        0x0b40_0000_0000_01e9,
        0x0d40_0000_0000_05b4,
        0x0d40_0000_0000_05b5,
        0x0f40_0000_0000_1118,
        0x0f40_0000_0000_1119,
    ];

    /// The sub-hexes C and D of the two polar roots at level 5, `CA-0-C` and its three
    /// fellows. The engine reads the text of each and draws it a ring; this crate has no cell
    /// for any (`HexA3::is_valid_address`).
    const POLAR_SUB_HEXES: [(&str, u64); 4] = [
        ("CA-0-C", 0x0540_0000_0000_0002),
        ("CA-0-D", 0x0540_0000_0000_0003),
        ("CB-0-C", 0x0560_0000_0000_0002),
        ("CB-0-D", 0x0560_0000_0000_0003),
    ];

    /// On aperture 3, an identifier that has no cell answers the empty ring, no
    /// extent and no area, and is the very identifier for which `vertices` gives its six zero
    /// points. Such are the null zone; the identifiers of a polar root with a non-zero index,
    /// for which the engine's ring is a null array; the sub-hexes C and D of a polar root,
    /// for which the engine has a ring; a level beyond the thirty-third; and an index beyond
    /// its rhombus.
    #[test]
    fn an_identifier_with_no_cell_has_no_refined_boundary_and_no_extent_on_aperture_3() {
        let others = [
            ZoneId::NULL,
            ZoneId(0x2200_0000_0000_0000), // I9R level 17
            ZoneId(0x0400_0000_0000_0144), // `C0-51-A`: the rhombus holds 81 cells
        ];
        for name in ["ISEA3H", "IVEA3H", "RTEA3H"] {
            let grid = get_grid(name).unwrap();
            let polar = POLAR_SUB_HEXES.map(|(_, id)| ZoneId(id));
            for id in NULL_RING_CLASS
                .map(ZoneId)
                .into_iter()
                .chain(others)
                .chain(polar)
            {
                assert_eq!(
                    grid.refined_vertices(id, 0),
                    Ok(Vec::new()),
                    "{name} {id:?}"
                );
                for refinement in [1, 2, 3, MAX_EDGE_REFINEMENT] {
                    assert_eq!(
                        grid.refined_vertices_within_limit(id, refinement),
                        Vec::new(),
                        "{name} {id:?}"
                    );
                }
                assert_eq!(grid.extent(id), None, "{name} {id:?}");
                assert_eq!(grid.area(id), None, "{name} {id:?}");
                assert_eq!(grid.vertices(id), vec![GeoPoint { lat: 0.0, lon: 0.0 }; 6]);
            }
            // The refusal of a refinement does not depend on the identifier.
            assert_eq!(
                grid.refined_vertices(ZoneId(NULL_RING_CLASS[0]), max_edge_refinement() + 1),
                Err(Error::EdgeRefinementOutOfRange {
                    refinement: max_edge_refinement() + 1,
                    max: max_edge_refinement()
                })
            );
        }
        let zone = crate::registry::isea3h().zone(ZoneId(NULL_RING_CLASS[0]));
        assert_eq!(zone.refined_vertices(0), Ok(Vec::new()));
        assert_eq!(zone.extent(), None);
    }

    /// No identifier makes the refined boundary or the extent panic on aperture 3: 3,000
    /// arbitrary identifiers on each grid, of every level to the seventeenth of I9R, of
    /// thirteen roots and of every sub-hex, a third of them with an index within its rhombus;
    /// and a ring is empty exactly where there is no extent.
    #[test]
    fn no_identifier_makes_the_refined_boundary_of_aperture_3_panic() {
        for name in ["ISEA3H", "IVEA3H", "RTEA3H"] {
            let grid = get_grid(name).unwrap();
            let mut state: u64 = 0x2026_1001;
            let mut draw = |below: u64| {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                (state >> 11) % below
            };
            let (mut with_a_ring, mut without) = (0, 0);
            for turn in 0..3_000 {
                let level = draw(18);
                let cells = 9u64.pow(level.min(16) as u32);
                let index = if turn % 3 == 0 {
                    draw(cells)
                } else {
                    draw(1 << 51)
                };
                let id = ZoneId((level << 57) | (draw(13) << 53) | (index << 2) | draw(4));
                let extent = grid.extent(id);
                for refinement in [0, 1, 3] {
                    let ring = grid.refined_vertices_within_limit(id, refinement);
                    assert_eq!(ring.is_empty(), extent.is_none(), "{name} {id:?}");
                }
                if extent.is_some() {
                    with_a_ring += 1;
                } else {
                    without += 1;
                }
            }
            assert!(
                with_a_ring >= 600 && without >= 1_500,
                "{name}: {with_a_ring}, {without}"
            );
        }
    }

    /// The constants of the refined boundary and of the extent against the doubles the library
    /// holds, and the two rules that the aperture decides.
    #[test]
    fn the_constants_of_the_refined_boundary_match_the_engines_rodata() {
        // `Degrees { 0.05 }`, at `0x49fe0`: the double of the eC compiler's own literal, and
        // not that of either expression of the angle.
        assert_eq!(REFINED_LON_SHIFT.to_bits(), 0x3f4c_9871_03b7_633a);
        assert_eq!(REFINED_LON_SHIFT, 0.000_872_664_625_997_2);
        assert_eq!((0.05 * math::PI / 180.0).to_bits(), 0x3f4c_9871_03b7_61f5);
        assert_eq!((0.05 * math::DEG2RAD).to_bits(), 0x3f4c_9871_03b7_61f5);
        assert_eq!(CENTROID_LON_LIMIT.to_bits(), 0x4009_21fb_5466_8930); // 0x49fd0
        assert_eq!((-CENTROID_LON_LIMIT).to_bits(), 0xc009_21fb_5466_8930); // 0x49fc8
        assert_eq!(math::PI.to_bits(), 0x4009_21fb_5444_2d18); // 0x49ac8
        assert_eq!((-math::PI).to_bits(), 0xc009_21fb_5444_2d18); // 0x49ab8
        assert_eq!(math::TWO_PI.to_bits(), 0x4019_21fb_5444_2d18); // 0x49f90
        // The cleared extent, at `0x4a4c8`: in degrees, one unit below the greatest double.
        assert_eq!(CLEARED_EXTENT.to_bits(), 0x7f91_df46_a252_9d38);
        assert_eq!(
            (CLEARED_EXTENT * math::RAD2DEG).to_bits(),
            f64::MAX.to_bits() - 1
        );
        // A pole, in degrees.
        assert_eq!(FRAC_PI_2 * math::RAD2DEG, 90.0);

        // The pole branch of aperture 3 (`RI3H.ec:529-638`): the latitude beyond which a
        // point is on a pole, the two bounds of its correction of half a turn, and the two
        // factors of `poleOffset` and of the distance it steps aside.
        assert_eq!(POLE_LATITUDE.to_bits(), 0x4056_7fff_fbce_4218); // 0x49ff8
        assert_eq!((-POLE_HALF_TURN_BEYOND).to_bits(), 0xc057_c000_0000_0000); // 0x4a000
        assert_eq!(POLE_HALF_TURN_BEYOND.to_bits(), 0x4057_c000_0000_0000); // 0x4a008
        assert_eq!(0.001_f64.to_bits(), 0x3f50_624d_d2f1_a9fc); // 0x49fc0
        assert_eq!((-1e-5_f64).to_bits(), 0xbee4_f8b5_88e3_68f1); // 0x49fd8
        // `Sgn(lat) * 90` taken to radians (0x49fb0, 0x49fb8, through 0x56e08) is the pole.
        assert_eq!((90.0 * math::DEG2RAD).to_bits(), FRAC_PI_2.to_bits());
        assert_eq!((-90.0 * math::DEG2RAD).to_bits(), (-FRAC_PI_2).to_bits());
        // The distance as this crate forms it is the library's, which multiplies
        // `poleOffset * -1e-5`, or its negation, by `poleOffset`, at every level.
        for level in 0..=33_u32 {
            let pole_offset = f64::from(1u32 << (level / 2)) * 0.001;
            let hoisted = pole_offset * -1e-5;
            let nudge = pole_offset * 1e-5 * pole_offset;
            assert_eq!((pole_offset * -hoisted).to_bits(), nudge.to_bits());
            assert_eq!((hoisted * pole_offset).to_bits(), (-nudge).to_bits());
        }

        let divisions = |aperture| -> Vec<i32> {
            (0..=12)
                .map(|level| automatic_edge_divisions(aperture, level))
                .collect()
        };
        assert_eq!(
            divisions(7),
            [20, 20, 20, 15, 15, 12, 12, 12, 12, 12, 12, 12, 12]
        );
        assert_eq!(
            divisions(3),
            [20, 20, 20, 15, 15, 10, 10, 10, 8, 8, 5, 5, 5]
        );
        assert_eq!(automatic_edge_divisions(7, 19), 12);
        assert_eq!(automatic_edge_divisions(3, 33), 5);
        assert!((0..=19).all(|level| refined_odd_grid(7, level)));
        assert!((0..=33).all(|level| refined_odd_grid(3, level) == (level % 2 == 1)));
    }

    /// A zone without geometry answers the empty ring and no extent, and it is the
    /// very zone for which `vertices` gives its six zero points and `area` nothing. The
    /// identifier that the engine cannot read back keeps its ordinary geometry.
    #[test]
    fn a_zone_without_geometry_has_no_refined_boundary_and_no_extent() {
        for grid in [get_grid("IGEO7"), get_grid("IVEA7H"), get_grid("RTEA7H")] {
            let grid = grid.unwrap();
            for id in [ZoneId::NULL, BASE_12, LEVEL_20[0], LEVEL_20[1]] {
                assert_eq!(grid.refined_vertices(id, 0), Ok(Vec::new()));
                for refinement in [1, 2, 3, MAX_EDGE_REFINEMENT] {
                    assert_eq!(
                        grid.refined_vertices_within_limit(id, refinement),
                        Vec::new()
                    );
                }
                assert_eq!(grid.extent(id), None);
                assert_eq!(grid.area(id), None);
                assert_eq!(grid.vertices(id), vec![GeoPoint { lat: 0.0, lon: 0.0 }; 6]);
            }

            assert_eq!(grid.text_id(UNREADABLE), "01222222222222222220");
            let points: Vec<usize> = [0, 1, 2, 3, 7]
                .iter()
                .map(|&n| grid.refined_vertices_within_limit(UNREADABLE, n).len())
                .collect();
            assert_eq!(points, [72, 6, 12, 18, 42], "{}", grid.name());
            let extent = grid.extent(UNREADABLE).unwrap();
            assert!(extent.ll.lat > 85.4 && extent.ur.lat < 85.5, "{extent:?}");
            assert!(extent.ll.lon > -88.8 && extent.ur.lon < -88.5, "{extent:?}");
            assert!(extent.ll.lat < extent.ur.lat && extent.ll.lon < extent.ur.lon);
            assert!(grid.area(UNREADABLE).is_some());
        }
        let zone = igeo7().zone(ZoneId::NULL);
        assert_eq!(zone.refined_vertices(0), Ok(Vec::new()));
        assert_eq!(zone.extent(), None);
    }

    /// A refinement above the limit in force is refused, whatever the
    /// zone, and the ceiling itself is taken by the ring, a hexagon then having 600,000 points,
    /// as the engine's has. The tests derive the limit from `max_edge_refinement()`, since the
    /// variable that lowers it may be exported and is read once per process.
    #[test]
    fn an_edge_refinement_above_the_maximum_is_refused() {
        assert_eq!(MAX_EDGE_REFINEMENT, 100_000);
        let max = max_edge_refinement();
        assert!(max <= MAX_EDGE_REFINEMENT);
        let g = igeo7();
        let hexagon = g.zone_from_text("0313522").unwrap();
        let beyond = max + 1;
        for id in [hexagon.id(), ZoneId::NULL, BASE_12, LEVEL_20[1]] {
            assert!(g.refined_vertices(id, max).is_ok());
            for refinement in [beyond, 107_374_183, u32::MAX] {
                assert_eq!(
                    g.refined_vertices(id, refinement),
                    Err(Error::EdgeRefinementOutOfRange { refinement, max })
                );
            }
            assert!(g.refined_vertices(id, 0).is_ok());
        }
        assert_eq!(
            get_grid("RTEA7H")
                .unwrap()
                .refined_vertices(hexagon.id(), beyond),
            Err(Error::EdgeRefinementOutOfRange {
                refinement: beyond,
                max
            })
        );
        assert_eq!(
            g.refined_vertices_within_limit(hexagon.id(), MAX_EDGE_REFINEMENT)
                .len(),
            600_000
        );
    }

    /// The reading of `RS4DGGS_MAX_MATERIALISED_SUB_ZONES`: it only lowers the ceiling, and a
    /// value that is not a `u64` is ignored.
    #[test]
    fn the_limit_on_a_list_of_sub_zones_is_read_from_the_variable() {
        for (value, limit) in [
            (None, 4_000_000),
            (Some("4000000"), 4_000_000),
            (Some("4000001"), 4_000_000),
            (Some("18446744073709551615"), 4_000_000),
            (Some("18446744073709551616"), 4_000_000),
            (Some("12"), 12),
            (Some(" 12 "), 12),
            (Some("+12"), 12),
            (Some("0012"), 12),
            (Some("0"), 0),
            (Some(""), 4_000_000),
            (Some("abc"), 4_000_000),
            (Some("-1"), 4_000_000),
            (Some("1e3"), 4_000_000),
            (Some("12 000"), 4_000_000),
        ] {
            assert_eq!(materialised_sub_zones_limit(value), limit, "{value:?}");
        }
    }

    /// A list of sub-zones longer than the limit in force is refused on every grid, with its
    /// length and that limit, and so is the index of a zone in an order so long; the count,
    /// the first sub-zone and the zone at an index are not limited, and a depth of 0, whose
    /// list is the zone itself, is answered whatever the limit. The limit is derived from
    /// `max_materialised_sub_zones()`, since the variable that lowers it may be exported and
    /// is read once per process.
    #[test]
    fn a_list_of_sub_zones_above_the_limit_in_force_is_refused() {
        assert_eq!(MAX_MATERIALISED_SUB_ZONES, 4_000_000);
        let limit = max_materialised_sub_zones();
        assert!(limit <= MAX_MATERIALISED_SUB_ZONES);
        for (name, zone) in [
            ("IGEO7", "0064156"),
            ("IVEA7H", "0064156"),
            ("RTEA7H", "0064156"),
            ("ISEA3H", "A4-0-A"),
            ("IVEA3H", "A4-0-A"),
            ("RTEA3H", "A4-0-A"),
        ] {
            let g = get_grid(name).unwrap();
            let z = g.zone_from_text(zone).unwrap();
            // The first depth at which the list is longer than the limit: under the ceiling
            // it is 8 below a hexagon of aperture 7 and 14 below one of aperture 3.
            let deepest = g.max_resolution() - g.resolution(z);
            let depth = (1..=deepest)
                .find(|&depth| g.count_sub_zones(z, depth).unwrap() > limit)
                .unwrap();
            let count = g.count_sub_zones(z, depth).unwrap();
            let too_many = Error::TooManySubZones { count, limit };
            assert_eq!(g.sub_zones(z, depth), Err(too_many.clone()), "{name}");
            let sub = g.sub_zone_at_index(z, depth, 1).unwrap();
            assert_eq!(g.sub_zone_index(z, sub), Err(too_many), "{name}");
            assert!(g.first_sub_zone(z, depth).is_ok(), "{name}");
            assert_eq!(g.sub_zones(z, 0), Ok(vec![z]), "{name}");
            // One level above, the list is within the limit: it is built, where it is
            // short, and never refused for want of an order.
            if depth > 1 {
                let within = g.count_sub_zones(z, depth - 1).unwrap();
                assert!(within <= limit, "{name}");
                if within <= 60_000 {
                    assert_eq!(
                        g.sub_zones(z, depth - 1).map(|order| order.len() as u64),
                        Ok(within),
                        "{name}"
                    );
                }
            }
        }
        if limit == MAX_MATERIALISED_SUB_ZONES {
            let g = igeo7();
            let z = g.zone_from_text("0064156").unwrap().id();
            assert_eq!(
                g.sub_zones(z, 8),
                Err(Error::TooManySubZones {
                    count: 5_767_201,
                    limit: 4_000_000
                })
            );
        }
    }

    /// The reading of `RS4DGGS_MAX_EDGE_REFINEMENT`: it only lowers the ceiling, and a value
    /// that is not a `u32` is ignored.
    #[test]
    fn the_limit_is_read_from_the_variable() {
        for (value, limit) in [
            (None, 100_000),
            (Some("100000"), 100_000),
            (Some("100001"), 100_000),
            (Some("4294967295"), 100_000),
            (Some("4294967296"), 100_000),
            (Some("12"), 12),
            (Some(" 12 "), 12),
            (Some("0"), 0),
            (Some(""), 100_000),
            (Some("abc"), 100_000),
            (Some("-1"), 100_000),
            (Some("1e3"), 100_000),
        ] {
            assert_eq!(edge_refinement_limit(value), limit, "{value:?}");
        }
    }

    /// No identifier makes the refined boundary or the extent panic: 2,000 arbitrary
    /// identifiers on each grid, most of which this crate never hands out, and among them
    /// addresses of every level to the twentieth and base cells that the icosahedron has not.
    #[test]
    fn no_identifier_makes_the_refined_boundary_panic() {
        for grid in [get_grid("IGEO7"), get_grid("IVEA7H"), get_grid("RTEA7H")] {
            let grid = grid.unwrap();
            let mut state: u64 = 0x2026_1001;
            let (mut with_a_ring, mut without) = (0, 0);
            for _ in 0..2_000 {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                // The high bits of the generator, and every digit slot from a level drawn
                // first, so that the deep levels are as frequent as the shallow.
                let level = (state >> 59) % 21;
                let unused = u64::MAX.checked_shr(4 + 3 * level as u32).unwrap_or(0);
                let id = ZoneId((state << 3) | unused);
                let extent = grid.extent(id);
                for refinement in [0, 1, 3] {
                    let ring = grid.refined_vertices_within_limit(id, refinement);
                    assert_eq!(ring.is_empty(), extent.is_none(), "{id:?}");
                }
                if extent.is_some() {
                    with_a_ring += 1;
                } else {
                    without += 1;
                }
            }
            assert!(
                with_a_ring >= 1_000 && without >= 300,
                "{with_a_ring}, {without}"
            );
        }
    }

    /// The pinning zones on aperture 7 that concern the ring and the
    /// extent, each with the count of points the engine gives it at the automatic refinement
    /// and at a refinement of 1 on IGEO7, IVEA7H and RTEA7H, in that order.
    ///
    /// The first four pin site GR-1 (`050` on IGEO7, at the automatic refinement alone; `0003`
    /// and `0010` on the three grids, and `0013` on IVEA7H and RTEA7H, at every refinement
    /// tried): put back in the eC's order alone, the site moves the ring of each. The rest is
    /// one zone of each kind: a hexagon, a pentagon, a pentagon at either pole of the
    /// icosahedron, a hexagon across an interruption, a hexagon with an edge beside a pole,
    /// which is divided twenty times as finely, a hexagon that touches a pole, the zone centred
    /// on a pole at resolution 19 and one beside it, a zone that loses a vertex whose inverse
    /// fails, a ring that crosses itself, and the identifier the engine writes for a neighbour
    /// and cannot read back.
    #[cfg(feature = "oracle")]
    const PINNED_ZONES: [(&str, u64, [usize; 3], [usize; 3]); 16] = [
        ("050", 0x51ff_ffff_ffff_ffff, [101; 3], [5; 3]),
        ("0003", 0x00ff_ffff_ffff_ffff, [120; 3], [6; 3]),
        ("0010", 0x023f_ffff_ffff_ffff, [120; 3], [6; 3]),
        ("0013", 0x02ff_ffff_ffff_ffff, [120; 3], [6; 3]),
        ("0313522", 0x32ea_5fff_ffff_ffff, [72; 3], [6; 3]),
        ("0100000", 0x1000_1fff_ffff_ffff, [61; 3], [5; 3]),
        ("0000000", 0x0000_1fff_ffff_ffff, [65; 3], [5; 3]),
        ("110000", 0xb000_ffff_ffff_ffff, [80; 3], [5; 3]),
        ("0100033", 0x1003_7fff_ffff_ffff, [74; 3], [6; 3]),
        ("0132336", 0x169b_dfff_ffff_ffff, [300; 3], [25; 3]),
        ("0132323", 0x169a_7fff_ffff_ffff, [108; 3], [26; 3]),
        (
            "013232323232323232323",
            0x169a_69a6_9a69_a69f,
            [1; 3],
            [1; 3],
        ),
        (
            "005151515151515151515",
            0x0a69_a69a_69a6_9a6f,
            [67, 65, 68],
            [24, 24, 25],
        ),
        (
            "0000000000000000131",
            0x0000_0000_0000_b3ff,
            [58; 3],
            [5; 3],
        ),
        (
            "01323232323232322",
            0x169a_69a6_9a69_7fff,
            [149, 149, 150],
            [46, 46, 47],
        ),
        ("01222222222222222220", UNREADABLE.0, [72; 3], [6; 3]),
    ];

    /// What a comparison with the engine has met.
    #[cfg(feature = "oracle")]
    #[derive(Debug, Default, PartialEq)]
    struct Census {
        zones: usize,
        /// The doubles compared, of rings and of extents.
        coordinates: usize,
        /// The zones with neither ring nor extent, on either side.
        without_geometry: usize,
        /// The rings, at the automatic refinement, that hold a longitude beyond half a turn.
        beyond_half_a_turn: usize,
        /// The extents whose western end is the greater.
        western_end_the_greater: usize,
        /// The rings, at the automatic refinement, that hold a pole.
        at_a_pole: usize,
        /// The rings, at the automatic refinement, of a single point.
        of_one_point: usize,
        /// The zones of level 15 and beyond that the engine cannot read back from its own text.
        unreadable: usize,
        /// The identifiers whose ring is not asked of the engine: on aperture 3, those the
        /// engine cannot read.
        ring_not_asked: usize,
    }

    /// Every way in which this crate's ring of the zone `id`, at each of `refinements`, and its
    /// extent leave the engine's: both in radians, as the engine gives them and as this crate
    /// computes them, compared as sequences of doubles by their bits. A zone without geometry
    /// must have the engine's empty ring at every refinement and the engine's cleared extent.
    #[cfg(feature = "oracle")]
    fn differences_from_the_live_engine<P: Projection, T: Topology, I: Indexing>(
        oracle: &str,
        grid: &Grid<P, T, I>,
        id: u64,
        refinements: &[i32],
        census: &mut Census,
    ) -> Vec<String> {
        let zone = ZoneId(id);
        let name = format!("{oracle} {} ({id:#018x})", grid.text_id(zone));
        let address = grid.geometry_address(zone);
        let mut differences = Vec::new();
        census.zones += 1;

        // On aperture 3 the engine answers a null array as the ring of an identifier of a
        // polar root with a non-zero index, on which the vendored binding crashes, and the
        // oracle refuses to ask the ring of anything the engine cannot read: that class, the
        // null zone and the levels beyond the maximum, whose rings the engine gives empty.
        // This crate has no address for any of them, and so answers the empty ring and no
        // extent: that is asserted here, the ring is not asked, and the zone is counted. Its
        // extent is asked of the engine below all the same.
        //
        // The converse does not hold: the engine reads the sub-hexes C and D of a polar root,
        // for which this crate has no address either. The ring of one is asked, and this
        // crate's empty ring is then reported as a difference, which it is.
        let asked = T::APERTURE != 3 || dggal_oracle::engine_can_read(oracle, id);
        if !asked {
            assert!(
                address.is_none(),
                "{name}: an address the engine cannot read"
            );
            census.ring_not_asked += 1;
        }

        for &refinement in refinements.iter().filter(|_| asked) {
            let ours = address
                .map(|a| grid.refined_ring_radians(zone, &a, refinement))
                .unwrap_or_default();
            let theirs = dggal_oracle::refined_vertices(oracle, id, refinement);
            if let Some(difference) = ring_difference(&ours, &theirs) {
                differences.push(format!(
                    "{name} at a refinement of {refinement}: {difference}"
                ));
            }
            census.coordinates += 2 * theirs.len();
            assert!(
                theirs.iter().all(|p| !p.0.is_nan() && !p.1.is_nan()),
                "{name}: the engine's ring holds a NaN"
            );
            if refinement == 0 {
                census.beyond_half_a_turn +=
                    usize::from(theirs.iter().any(|p| p.1.abs() > math::PI));
                census.at_a_pole += usize::from(theirs.iter().any(|p| p.0.abs() == FRAC_PI_2));
                census.of_one_point += usize::from(theirs.len() == 1);
            }
        }

        let ours = address.and_then(|a| grid.extent_radians(zone, &a));
        let theirs = dggal_oracle::extent(oracle, id);
        let cleared = [
            CLEARED_EXTENT,
            CLEARED_EXTENT,
            -CLEARED_EXTENT,
            -CLEARED_EXTENT,
        ];
        if ours.unwrap_or(cleared).map(f64::to_bits) != theirs.map(f64::to_bits) {
            differences.push(format!(
                "{name}: the extent is {ours:?}, the engine has {theirs:?}"
            ));
        }
        match ours {
            Some([_, ll_lon, _, ur_lon]) => {
                census.coordinates += 4;
                census.western_end_the_greater += usize::from(ll_lon > ur_lon);
                // What `Extent` promises of the two ends.
                assert!(
                    ll_lon.abs() <= math::PI && ur_lon.abs() <= math::PI,
                    "{name}: an end of the extent beyond half a turn, {ours:?}"
                );
            }
            None => census.without_geometry += 1,
        }
        differences
    }

    /// The zones of `PINNED_ZONES`, and the zones without geometry, against the live engine on
    /// the three grids: the ring at the automatic refinement and at 1, 2, 3 and 7, and the
    /// extent. Every difference is gathered before the test fails, so that site GR-1 put back
    /// in the eC's order names its zones in the failure. The counts of points are the
    /// engine's, asserted, so that each zone is seen to be of its kind.
    #[test]
    #[cfg(feature = "oracle")]
    fn the_pinned_zones_have_the_engines_ring_and_extent() {
        fn on<P: Projection>(
            oracle: &str,
            grid: &Grid<P, HexA7, Z7>,
            column: usize,
        ) -> Vec<String> {
            let mut census = Census::default();
            let mut differences = Vec::new();
            for (text, id, at_0, at_1) in PINNED_ZONES {
                assert_eq!(dggal_oracle::text_id(oracle, id), text);
                assert_eq!(
                    dggal_oracle::refined_vertices(oracle, id, 0).len(),
                    at_0[column],
                    "{oracle} {text} at the automatic refinement"
                );
                assert_eq!(
                    dggal_oracle::refined_vertices(oracle, id, 1).len(),
                    at_1[column],
                    "{oracle} {text} at a refinement of 1"
                );
                differences.extend(differences_from_the_live_engine(
                    oracle,
                    grid,
                    id,
                    &[0, 1, 2, 3, 7],
                    &mut census,
                ));
                // The public answers are the engine's, in degrees.
                let ring: Vec<GeoPoint> = dggal_oracle::refined_vertices(oracle, id, 0)
                    .into_iter()
                    .map(degrees)
                    .collect();
                let [ll_lat, ll_lon, ur_lat, ur_lon] = dggal_oracle::extent(oracle, id);
                let extent = Extent {
                    ll: degrees((ll_lat, ll_lon)),
                    ur: degrees((ur_lat, ur_lon)),
                };
                if grid.refined_vertices(ZoneId(id), 0) != Ok(ring) {
                    differences.push(format!(
                        "{oracle} {text}: the public ring is not the engine's"
                    ));
                }
                if grid.extent(ZoneId(id)) != Some(extent) {
                    differences.push(format!(
                        "{oracle} {text}: the public extent is not the engine's"
                    ));
                }
            }
            assert_eq!(census.without_geometry, 0);

            // The engine reads only the last of these back, as itself.
            assert!(!dggal_oracle::engine_can_read(oracle, UNREADABLE.0));
            for id in [LEVEL_20[0].0, LEVEL_20[1].0, dggal_oracle::NULL_ZONE] {
                differences.extend(differences_from_the_live_engine(
                    oracle,
                    grid,
                    id,
                    &[0, 1, 2, 3, 7],
                    &mut census,
                ));
            }
            assert_eq!(census.without_geometry, 3);
            differences
        }

        let mut differences = on(dggal_oracle::IGEO7, igeo7(), 0);
        differences.extend(on(dggal_oracle::IVEA7H, ivea7h(), 1));
        differences.extend(on(dggal_oracle::RTEA7H, rtea7h(), 2));
        assert!(
            differences.is_empty(),
            "{} answers leave the engine's:\n{}",
            differences.len(),
            differences.join("\n")
        );
    }

    /// The greatest refinement accepted, against the engine: a hexagon, a hexagon that touches
    /// a pole and one with an edge beside a pole, which is divided twenty times as finely, on
    /// IGEO7, each edge divided into 100,000 parts.
    #[test]
    #[cfg(feature = "oracle")]
    fn the_greatest_refinement_is_the_engines_too() {
        let mut census = Census::default();
        for (id, points) in [
            (0x32ea_5fff_ffff_ffff, 600_000),
            (0x169a_7fff_ffff_ffff, 600_029),
            (0x169b_dfff_ffff_ffff, 2_500_000),
        ] {
            let differences = differences_from_the_live_engine(
                dggal_oracle::IGEO7,
                igeo7(),
                id,
                &[MAX_EDGE_REFINEMENT as i32],
                &mut census,
            );
            assert!(differences.is_empty(), "{differences:?}");
            assert_eq!(
                igeo7()
                    .refined_vertices_within_limit(ZoneId(id), MAX_EDGE_REFINEMENT)
                    .len(),
                points
            );
        }
    }

    /// The zones at which a census below compares this crate with the engine, on one grid.
    ///
    /// Every zone to the level `generations`, taken from the engine's own lists of children
    /// from the twelve zones `bases` of level 0, and counted; the zones `pinned`; and, at
    /// every level of the grid, the zones in which the engine places both poles at three
    /// longitudes, points beside the poles from 1e-7 degrees away, points on and beside the
    /// antimeridian, `points_along_an_edge` points along each of the thirty edges of the
    /// icosahedron, and a fixed pseudo-random set of `random_points_a_level`, uniform over the
    /// sphere. Whatever the engine answers is kept: its null zone, which on aperture 7 it
    /// answers on two of those edges at the odd levels from 15, and on aperture 3 the
    /// identifiers it hands out and cannot read.
    #[cfg(feature = "oracle")]
    fn census_zones(
        oracle: &str,
        bases: &[u64],
        generations: usize,
        zones_to_that_level: usize,
        pinned: &[u64],
        points_along_an_edge: usize,
        random_points_a_level: usize,
    ) -> std::collections::BTreeSet<u64> {
        assert_eq!(bases.len(), 12);
        let mut zones: std::collections::BTreeSet<u64> = std::collections::BTreeSet::new();
        let mut generation: std::collections::BTreeSet<u64> = bases.iter().copied().collect();
        for _ in 0..generations {
            zones.extend(&generation);
            generation = generation
                .iter()
                .flat_map(|&zone| dggal_oracle::children(oracle, zone))
                .collect();
        }
        zones.extend(&generation);
        assert_eq!(
            zones.len(),
            zones_to_that_level,
            "{oracle}: to level {generations}"
        );
        zones.extend(pinned);

        let mut points: Vec<(f64, f64)> = Vec::new();
        for pole in [90.0, -90.0] {
            for lon in [0.0, 57.3, -143.2] {
                points.push((pole, lon));
            }
            let mut away = 3e-7;
            while away < 30.0 {
                for lon in [17.2, 114.6, -68.8, -166.2] {
                    points.push((pole - away * f64::signum(pole), lon));
                }
                away *= 3.0;
            }
        }
        for step in -14..=14 {
            let lat = 6.0 * f64::from(step);
            points.extend([(lat, 180.0), (lat, -180.0), (lat, 180.0 - 5e-6)]);
        }
        // The vertices of the icosahedron are the centroids of the twelve zones of level 0, and
        // its edges join those some 63.4 degrees apart.
        let unit = |(lat, lon): (f64, f64)| {
            let (lat, lon) = (lat * math::DEG2RAD, lon * math::DEG2RAD);
            [
                math::cos(lat) * math::cos(lon),
                math::cos(lat) * math::sin(lon),
                math::sin(lat),
            ]
        };
        let vertices: Vec<[f64; 3]> = bases
            .iter()
            .map(|&base| unit(dggal_oracle::centroid(oracle, base)))
            .collect();
        let mut edges = 0;
        for (i, a) in vertices.iter().enumerate() {
            for b in &vertices[i + 1..] {
                let angle = math::acos((a[0] * b[0] + a[1] * b[1] + a[2] * b[2]).clamp(-1.0, 1.0));
                if angle >= 1.2 {
                    continue;
                }
                edges += 1;
                for part in 1..=points_along_an_edge {
                    let t = part as f64 / (points_along_an_edge + 1) as f64;
                    let (wa, wb) = (
                        math::sin((1.0 - t) * angle) / math::sin(angle),
                        math::sin(t * angle) / math::sin(angle),
                    );
                    let p = [
                        wa * a[0] + wb * b[0],
                        wa * a[1] + wb * b[1],
                        wa * a[2] + wb * b[2],
                    ];
                    points.push((
                        math::asin(p[2].clamp(-1.0, 1.0)) * math::RAD2DEG,
                        math::atan2(p[1], p[0]) * math::RAD2DEG,
                    ));
                }
            }
        }
        assert_eq!(edges, 30, "{oracle}: the edges of the icosahedron");

        // A linear congruential generator, whose high bits give a point uniform over the
        // sphere.
        let mut state: u64 = 0x2026_1001;
        let mut uniform = || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        for level in 0..=dggal_oracle::max_level(oracle) {
            let random: Vec<(f64, f64)> = (0..random_points_a_level)
                .map(|_| {
                    (
                        math::asin(2.0 * uniform() - 1.0) * math::RAD2DEG,
                        360.0 * uniform() - 180.0,
                    )
                })
                .collect();
            for &(lat, lon) in points.iter().chain(&random) {
                zones.insert(dggal_oracle::zone_from_geo(oracle, lat, lon, level));
            }
        }
        zones
    }

    /// The census: the ring at the automatic refinement and at 1, 2 and 3, and the extent, of
    /// every zone of `census_zones`, against the live engine on the three aperture-7 grids,
    /// bit for bit and as sequences. No answer differs; the counts are asserted exactly, so
    /// that the test cannot quietly compare less than it claims.
    #[test]
    #[cfg(feature = "oracle")]
    fn the_refined_ring_and_the_extent_are_the_engines_at_every_level() {
        fn on<P: Projection>(oracle: &str, grid: &Grid<P, HexA7, Z7>) -> Census {
            let bases: Vec<u64> = (0..12)
                .map(|base| dggal_oracle::zone_from_text(oracle, &format!("{base:02}")))
                .collect();
            let pinned: Vec<u64> = PINNED_ZONES.iter().map(|pinned| pinned.1).collect();
            let zones = census_zones(oracle, &bases, 3, 12 + 72 + 492 + 3432, &pinned, 7, 60);
            let mut census = Census::default();
            let mut levels = [0usize; 20];
            for &id in &zones {
                if id != dggal_oracle::NULL_ZONE {
                    let level = usize::from(grid.resolution(ZoneId(id)));
                    levels[level] += 1;
                    // The engine writes an identifier that it cannot read back only in the
                    // broken seams, from level 15.
                    if level >= 15 && !dggal_oracle::engine_can_read(oracle, id) {
                        census.unreadable += 1;
                    }
                }
                let differences =
                    differences_from_the_live_engine(oracle, grid, id, &[0, 1, 2, 3], &mut census);
                assert!(differences.is_empty(), "{}", differences.join("\n"));
            }
            println!("{oracle}: {census:?}, by level {levels:?}");
            assert!(
                levels[4..].iter().all(|&count| count >= 300),
                "{oracle}: zones by level {levels:?}"
            );
            // The one zone without geometry is the engine's null zone.
            assert!(zones.contains(&dggal_oracle::NULL_ZONE));
            census
        }

        let compared = [
            on(dggal_oracle::IGEO7, igeo7()),
            on(dggal_oracle::IVEA7H, ivea7h()),
            on(dggal_oracle::RTEA7H, rtea7h()),
        ];
        // On each grid the one zone without geometry is the null zone and the one identifier
        // that the engine cannot read back is `01222222222222222220`; ninety-four rings hold a
        // pole, and three, at level 19, are a pole alone; and a ring holds a longitude beyond
        // half a turn as often as the western end of an extent is the greater.
        let expected = |zones, coordinates, beyond_half_a_turn| Census {
            zones,
            coordinates,
            without_geometry: 1,
            beyond_half_a_turn,
            western_end_the_greater: beyond_half_a_turn,
            at_a_pole: 94,
            of_one_point: 3,
            unreadable: 1,
            ring_not_asked: 0,
        };
        assert_eq!(
            compared,
            [
                expected(10_107, 2_533_434, 653),
                expected(10_090, 2_529_576, 641),
                expected(10_103, 2_533_310, 647),
            ]
        );
    }

    /// The pinning zones on aperture 3 that concern what the two
    /// topologies share and the ring itself, each with the site it pins and with the count of
    /// points the engine gives it at the automatic refinement and at a refinement of 1 on
    /// ISEA3H, IVEA3H and RTEA3H, in that order.
    ///
    /// The sites are those of `addIntermediatePoints`, AIP-1 to AIP-3, which meet the engine
    /// on aperture 3 here, and GR-1: put back in the eC's order alone, each moves the ring of
    /// its zones, on one of the three grids at least. The rest is one zone of each kind, as
    /// `a_zone_of_each_kind_has_the_engines_ring_and_extent_on_aperture_3` sets them out.
    #[cfg(feature = "oracle")]
    #[allow(clippy::type_complexity, reason = "the table of pinned zones")]
    const PINNED_ZONES_OF_APERTURE_3: [(&str, &str, u64, [usize; 3], [usize; 3]); 32] = [
        (
            "AIP-1",
            "BB-0-A",
            0x0360_0000_0000_0000,
            [100, 100, 100],
            [5, 5, 5],
        ),
        (
            "AIP-1",
            "E0-4-A",
            0x0800_0000_0000_0010,
            [48, 48, 48],
            [6, 6, 6],
        ),
        (
            "AIP-1",
            "Q5-34AA7F2578AE1-A",
            0x20ad_2a9f_c95e_2b84,
            [130, 129, 131],
            [64, 64, 65],
        ),
        (
            "AIP-1",
            "Q5-34AA7F4E86222-A",
            0x20ad_2a9f_d3a1_8888,
            [90, 89, 92],
            [44, 44, 46],
        ),
        (
            "AIP-1",
            "Q5-34AA7FA0A10A4-A",
            0x20ad_2a9f_e828_4290,
            [90, 89, 92],
            [44, 44, 46],
        ),
        (
            "AIP-1",
            "F1-16C8-A",
            0x0a20_0000_0000_5b20,
            [32, 32, 32],
            [6, 6, 6],
        ),
        (
            "AIP-2",
            "B3-6-A",
            0x0260_0000_0000_0018,
            [120, 120, 120],
            [6, 6, 6],
        ),
        (
            "AIP-2",
            "B6-1-A",
            0x02c0_0000_0000_0004,
            [120, 120, 120],
            [6, 6, 6],
        ),
        (
            "AIP-2",
            "C3-3F-A",
            0x0460_0000_0000_00fc,
            [92, 92, 92],
            [6, 6, 6],
        ),
        (
            "AIP-2",
            "E0-23-A",
            0x0800_0000_0000_008c,
            [48, 48, 48],
            [6, 6, 6],
        ),
        (
            "AIP-2",
            "C0-3-A",
            0x0400_0000_0000_000c,
            [92, 92, 92],
            [6, 6, 6],
        ),
        (
            "AIP-2",
            "C3-2D-A",
            0x0460_0000_0000_00b4,
            [92, 92, 92],
            [6, 6, 6],
        ),
        (
            "AIP-3",
            "A1-0-A",
            0x0020_0000_0000_0000,
            [100, 100, 100],
            [5, 5, 5],
        ),
        (
            "AIP-3",
            "A0-0-B",
            0x0000_0000_0000_0001,
            [100, 100, 100],
            [5, 5, 5],
        ),
        (
            "AIP-3",
            "A0-0-D",
            0x0000_0000_0000_0003,
            [120, 120, 120],
            [6, 6, 6],
        ),
        (
            "AIP-3",
            "A1-0-B",
            0x0020_0000_0000_0001,
            [100, 100, 100],
            [5, 5, 5],
        ),
        (
            "GR-1",
            "B4-0-A",
            0x0280_0000_0000_0000,
            [100, 100, 100],
            [5, 5, 5],
        ),
        (
            "GR-1",
            "B4-4-A",
            0x0280_0000_0000_0010,
            [120, 120, 120],
            [6, 6, 6],
        ),
        (
            "GR-1",
            "B4-5-A",
            0x0280_0000_0000_0014,
            [120, 120, 120],
            [6, 6, 6],
        ),
        (
            "a kind",
            "C4-49-B",
            0x0480_0000_0000_0125,
            [60, 60, 60],
            [6, 6, 6],
        ),
        (
            "a kind",
            "C0-0-B",
            0x0400_0000_0000_0001,
            [50, 50, 50],
            [5, 5, 5],
        ),
        (
            "a kind",
            "CA-0-B",
            0x0540_0000_0000_0001,
            [50, 50, 50],
            [5, 5, 5],
        ),
        (
            "a kind",
            "CB-0-A",
            0x0560_0000_0000_0000,
            [80, 80, 80],
            [5, 5, 5],
        ),
        (
            "a kind",
            "C0-6-A",
            0x0400_0000_0000_0018,
            [92, 92, 92],
            [6, 6, 6],
        ),
        (
            "a kind",
            "E0-7A-A",
            0x0800_0000_0000_01e8,
            [200, 200, 200],
            [25, 25, 25],
        ),
        (
            "a kind",
            "C0-4-C",
            0x0400_0000_0000_0012,
            [251, 251, 251],
            [26, 26, 26],
        ),
        (
            "a kind",
            "F0-79-A",
            0x0a00_0000_0000_01e4,
            [64, 64, 65],
            [29, 29, 30],
        ),
        (
            "a kind",
            "M0-40DF8-A",
            0x1800_0000_0010_37e0,
            [62, 62, 62],
            [25, 25, 25],
        ),
        (
            "a kind",
            "Q0-1486BA0-C",
            0x2000_0000_0521_ae82,
            [248, 47, 248],
            [48, 7, 48],
        ),
        (
            "a kind",
            "Q8-34AA7F2578AE0-C",
            0x210d_2a9f_c95e_2b82,
            [47, 47, 47],
            [7, 7, 7],
        ),
        (
            "a kind",
            "A0-0-A",
            0x0000_0000_0000_0000,
            [139, 139, 139],
            [24, 24, 24],
        ),
        // A hexagon of level 2 whose boundary reaches the south pole (its extent has `ll.lat` of
        // -90), kept as a case of the pole branch on the southern side.
        (
            "a kind",
            "B5-6-A",
            0x02a0_0000_0000_0018,
            [158, 158, 158],
            [25, 25, 25],
        ),
    ];

    /// The zones of `PINNED_ZONES_OF_APERTURE_3` against the live engine on the three grids:
    /// the ring at the automatic refinement and at 1, 2, 3 and 7, and the extent. Every
    /// difference is gathered before the test fails, so that a site put back in the eC's
    /// order names its zones in the failure. Then the identifiers the engine cannot read: the
    /// thirty of `NULL_RING_CLASS` and the null zone, whose ring is not asked of the engine
    /// and whose extent is the engine's cleared one.
    #[test]
    #[cfg(feature = "oracle")]
    fn the_pinned_zones_of_aperture_3_have_the_engines_ring_and_extent() {
        fn on<P: Projection>(
            oracle: &str,
            grid: &Grid<P, HexA3, I3h>,
            column: usize,
        ) -> Vec<String> {
            let mut census = Census::default();
            let mut differences = Vec::new();
            for (site, text, id, at_0, at_1) in PINNED_ZONES_OF_APERTURE_3 {
                assert_eq!(dggal_oracle::text_id(oracle, id), text);
                assert_eq!(
                    dggal_oracle::refined_vertices(oracle, id, 0).len(),
                    at_0[column],
                    "{oracle} {text} at the automatic refinement"
                );
                assert_eq!(
                    dggal_oracle::refined_vertices(oracle, id, 1).len(),
                    at_1[column],
                    "{oracle} {text} at a refinement of 1"
                );
                let found = differences_from_the_live_engine(
                    oracle,
                    grid,
                    id,
                    &[0, 1, 2, 3, 7],
                    &mut census,
                );
                differences.extend(found.iter().map(|found| format!("{site}: {found}")));
                // The public answers are the engine's, in degrees.
                let ring: Vec<GeoPoint> = dggal_oracle::refined_vertices(oracle, id, 0)
                    .into_iter()
                    .map(degrees)
                    .collect();
                let [ll_lat, ll_lon, ur_lat, ur_lon] = dggal_oracle::extent(oracle, id);
                let extent = Extent {
                    ll: degrees((ll_lat, ll_lon)),
                    ur: degrees((ur_lat, ur_lon)),
                };
                if grid.refined_vertices(ZoneId(id), 0) != Ok(ring) {
                    differences.push(format!(
                        "{site}: {oracle} {text}: the public ring is not the engine's"
                    ));
                }
                if grid.extent(ZoneId(id)) != Some(extent) {
                    differences.push(format!(
                        "{site}: {oracle} {text}: the public extent is not the engine's"
                    ));
                }
            }
            assert_eq!((census.without_geometry, census.ring_not_asked), (0, 0));

            for id in NULL_RING_CLASS.into_iter().chain([dggal_oracle::NULL_ZONE]) {
                assert!(!dggal_oracle::engine_can_read(oracle, id));
                differences.extend(differences_from_the_live_engine(
                    oracle,
                    grid,
                    id,
                    &[0, 1, 2, 3, 7],
                    &mut census,
                ));
            }
            assert_eq!((census.without_geometry, census.ring_not_asked), (31, 31));
            // To each identifier of the class the engine gives the area of a zone of its
            // level, as it gives the polar pentagon's geometry to its plain vertices; this
            // crate gives it none.
            for id in NULL_RING_CLASS {
                let area = dggal_oracle::area(oracle, id);
                assert!(
                    area.is_finite() && area > 0.0,
                    "{oracle} {id:#018x}: {area}"
                );
                assert_eq!(grid.area(ZoneId(id)), None);
            }

            // The identifiers that the engine can read and this crate has no cell for: the
            // sub-hexes C and D of a polar root, at every odd level. The engine reads their
            // text back and draws each a hexagon's ring, an extent and an area, though no
            // relation of any zone reaches them; this crate, which takes them for no cell
            // (`HexA3::is_valid_address`), answers the empty ring, no extent and no area. A departure from the engine, characterised here and not mirrored.
            for (text, id) in POLAR_SUB_HEXES {
                assert_eq!(dggal_oracle::text_id(oracle, id), text);
                assert!(dggal_oracle::engine_can_read(oracle, id), "{oracle} {text}");
                assert_eq!(dggal_oracle::refined_vertices(oracle, id, 0).len(), 60);
                assert_eq!(dggal_oracle::refined_vertices(oracle, id, 1).len(), 6);
                let extent = dggal_oracle::extent(oracle, id);
                assert!(
                    extent.iter().all(|member| member.abs() < 4.0),
                    "{oracle} {text}: {extent:?}"
                );
                assert!(dggal_oracle::area(oracle, id) > 0.0, "{oracle} {text}");
                assert!(
                    grid.geometry_address(ZoneId(id)).is_none(),
                    "{oracle} {text}"
                );
                assert_eq!(grid.refined_vertices(ZoneId(id), 0), Ok(Vec::new()));
                assert_eq!(grid.extent(ZoneId(id)), None);
                assert_eq!(grid.area(ZoneId(id)), None);
            }
            differences
        }

        let mut differences = on(dggal_oracle::ISEA3H, crate::registry::isea3h(), 0);
        differences.extend(on(dggal_oracle::IVEA3H, crate::registry::ivea3h(), 1));
        differences.extend(on(dggal_oracle::RTEA3H, crate::registry::rtea3h(), 2));
        assert!(
            differences.is_empty(),
            "{} answers leave the engine's:\n{}",
            differences.len(),
            differences.join("\n")
        );
    }

    /// The greatest refinement accepted, against the engine on aperture 3: a hexagon, a
    /// hexagon with an edge beside a pole and one that touches a pole, on ISEA3H, each edge
    /// divided into 100,000 parts.
    #[test]
    #[cfg(feature = "oracle")]
    fn the_greatest_refinement_is_the_engines_on_aperture_3() {
        let mut census = Census::default();
        let grid = crate::registry::isea3h();
        for (id, points) in [
            (0x0480_0000_0000_0125, 600_000),
            (0x0800_0000_0000_01e8, 2_500_000),
            (0x0400_0000_0000_0012, 2_500_001),
        ] {
            let differences = differences_from_the_live_engine(
                dggal_oracle::ISEA3H,
                grid,
                id,
                &[MAX_EDGE_REFINEMENT as i32],
                &mut census,
            );
            assert!(differences.is_empty(), "{differences:?}");
            assert_eq!(
                grid.refined_vertices_within_limit(ZoneId(id), MAX_EDGE_REFINEMENT)
                    .len(),
                points
            );
        }
    }

    /// The census on aperture 3: the ring at the automatic refinement and at 1, 2 and 3, and
    /// the extent, of every zone of `census_zones`, against the live engine on the three
    /// grids, bit for bit and as sequences. No answer differs; the counts are asserted
    /// exactly, so that the test cannot quietly compare less than it claims.
    ///
    /// The identifiers the engine hands out and cannot read are kept and counted: this crate
    /// answers each the empty ring and no extent, the engine its cleared extent, and its ring
    /// is not asked (see `differences_from_the_live_engine`).
    #[test]
    #[cfg(feature = "oracle")]
    fn the_refined_ring_and_the_extent_of_aperture_3_are_the_engines_at_every_level() {
        fn on<P: Projection>(oracle: &str, grid: &Grid<P, HexA3, I3h>) -> Census {
            let bases: Vec<u64> = (0..12)
                .map(|root| dggal_oracle::zone_from_text(oracle, &format!("A{root:X}-0-A")))
                .collect();
            let pinned: Vec<u64> = PINNED_ZONES_OF_APERTURE_3
                .iter()
                .map(|pinned| pinned.2)
                .chain(NULL_RING_CLASS)
                .collect();
            let zones = census_zones(
                oracle,
                &bases,
                CENSUS_GENERATIONS,
                CENSUS_ZONES_TO_THAT_LEVEL,
                &pinned,
                CENSUS_POINTS_ALONG_AN_EDGE,
                CENSUS_RANDOM_POINTS_A_LEVEL,
            );
            let mut census = Census::default();
            let mut levels = [0usize; 34];
            for &id in &zones {
                if grid.geometry_address(ZoneId(id)).is_some() {
                    levels[usize::from(grid.resolution(ZoneId(id)))] += 1;
                }
                let differences =
                    differences_from_the_live_engine(oracle, grid, id, &[0, 1, 2, 3], &mut census);
                assert!(differences.is_empty(), "{}", differences.join("\n"));
            }
            println!("{oracle}: {census:?}, by level {levels:?}");
            assert!(
                levels[5..].iter().all(|&count| count >= 150),
                "{oracle}: zones by level {levels:?}"
            );
            // The engine's null zone is in the sample of itself: the engine answers it for
            // points of the sample. With the thirty of `NULL_RING_CLASS` it makes the
            // thirty-one identifiers without geometry that the counts below assert, so that
            // the sample meets no other identifier that the engine cannot read.
            assert!(zones.contains(&dggal_oracle::NULL_ZONE));
            census
        }
        const CENSUS_GENERATIONS: usize = 4;
        const CENSUS_ZONES_TO_THAT_LEVEL: usize = 12 + 32 + 92 + 272 + 812;
        const CENSUS_POINTS_ALONG_AN_EDGE: usize = 3;
        const CENSUS_RANDOM_POINTS_A_LEVEL: usize = 24;

        let compared = [
            on(dggal_oracle::ISEA3H, crate::registry::isea3h()),
            on(dggal_oracle::IVEA3H, crate::registry::ivea3h()),
            on(dggal_oracle::RTEA3H, crate::registry::rtea3h()),
        ];
        // On each grid thirty-one identifiers have neither ring nor extent, and theirs are the
        // rings not asked of the engine; 148 rings hold a pole, and none is a single point;
        // and a ring holds a longitude beyond half a turn as often as the western end of an
        // extent is the greater.
        let expected = |zones, coordinates, beyond_half_a_turn| Census {
            zones,
            coordinates,
            without_geometry: 31,
            beyond_half_a_turn,
            western_end_the_greater: beyond_half_a_turn,
            at_a_pole: 148,
            of_one_point: 0,
            unreadable: 0,
            ring_not_asked: 31,
        };
        assert_eq!(
            compared,
            [
                expected(7_613, 1_340_292, 1_049),
                expected(7_616, 1_337_242, 1_065),
                expected(7_626, 1_344_734, 1_046),
            ]
        );
    }

    // The engine's answers of two zones that pin site GR-1 and of one zone of each kind, on each
    // of the three grids.
    const IGEO7_ENGINE: [EngineZone; 8] = [
        EngineZone {
            text: "0010",
            id: 0x023f_ffff_ffff_ffff,
            points_at_0: 120,
            extent: [
                1.0646009269447458,
                0.8512256007476848,
                1.2579496219837465,
                1.314096656025106,
            ],
            ring_at_1: &[
                (1.0646009269447458, 0.9853152791746168),
                (1.0767933034139172, 1.1829183308202427),
                (1.1647322081531561, 1.314096656025106),
                (1.2579065616078107, 1.2300208113847044),
                (1.2466011469937759, 0.9381757640163388),
                (1.140948820615925, 0.8512256007476848),
            ],
        },
        EngineZone {
            text: "0003",
            id: 0x00ff_ffff_ffff_ffff,
            points_at_0: 120,
            extent: [
                0.79543500221503,
                0.1954768762233649,
                0.9833030093411923,
                0.5035376783907152,
            ],
            ring_at_1: &[
                (0.8444268725659635, 0.46684954748521873),
                (0.9369652874902633, 0.5035376783907152),
                (0.9833030093411923, 0.3600170411994272),
                (0.9234723105095722, 0.19547687622336493),
                (0.8289399532286812, 0.19547687622336496),
                (0.79543500221503, 0.3251779503316056),
            ],
        },
        EngineZone {
            text: "0313522",
            id: 0x32ea_5fff_ffff_ffff,
            points_at_0: 72,
            extent: [
                0.2354736313080461,
                -0.7736012078236408,
                0.24641263036766073,
                -0.7631588486829763,
            ],
            ring_at_1: &[
                (0.2445533532421138, -0.772101813357245),
                (0.23908424958821012, -0.7736012078236408),
                (0.2354736313080461, -0.769884780944429),
                (0.2373298154682325, -0.764670978008093),
                (0.2427996620206654, -0.7631588486829763),
                (0.24641263036766073, -0.7668732236159013),
            ],
        },
        EngineZone {
            text: "0100000",
            id: 0x1000_1fff_ffff_ffff,
            points_at_0: 61,
            extent: [
                1.0146312230643104,
                -2.9544187618658073,
                1.0238794028591711,
                -2.9369865992760915,
            ],
            ring_at_1: &[
                (1.0217476386100264, -2.9539105827109085),
                (1.0161224498252872, -2.9530814821709406),
                (1.0147718995768933, -2.94270418748007),
                (1.0195440465762722, -2.9369865992760915),
                (1.0238794028591711, -2.943896035153389),
            ],
        },
        EngineZone {
            text: "0000000",
            id: 0x0000_1fff_ffff_ffff,
            points_at_0: 65,
            extent: [
                1.0146371165807027,
                0.18722861598748797,
                1.0238794028591682,
                0.20460605431370169,
            ],
            ring_at_1: &[
                (1.0217476386100233, 0.18768207087888522),
                (1.0161224498252845, 0.18851117141885287),
                (1.0147718995768904, 0.19888846610972272),
                (1.019544046576269, 0.20460605431370169),
                (1.0238794028591682, 0.19769661843640396),
            ],
        },
        EngineZone {
            text: "0100033",
            id: 0x1003_7fff_ffff_ffff,
            points_at_0: 74,
            extent: [
                1.0439968318023491,
                -2.9645021444703623,
                1.0531686500358184,
                -2.9415066541079073,
            ],
            ring_at_1: &[
                (1.0518013564667712, -2.962308587113436),
                (1.0472157794867838, -2.9645021444703623),
                (1.0439968318023491, -2.955256737733471),
                (1.0454096839005271, -2.9438275269661442),
                (1.0500809106551199, -2.9415066541079073),
                (1.0531686500358184, -2.95074610164365),
            ],
        },
        EngineZone {
            text: "0132323",
            id: 0x169a_7fff_ffff_ffff,
            points_at_0: 108,
            extent: [
                1.5609065075830497,
                1.1104356612552402,
                FRAC_PI_2,
                -2.0310138127084336,
            ],
            ring_at_1: &[
                (1.5633026229608513, 1.7662732030180608),
                (1.5609065075830497, 2.290060193787843),
                (1.5615089610942188, 2.857778830566318),
                (1.5641191403972365, 3.498141692078455),
                (1.5680941371315655, 4.25202831484609),
                (1.5684390736367453, 4.252046516998501),
                (1.568784008938155, 4.252064741301931),
                (1.5691289430647055, 4.252082987696969),
                (1.569473876045275, 4.252101256125398),
                (1.5698188079087034, 4.252119546528103),
                (1.5701637386837908, 4.252137858847174),
                (1.5705086683993035, 4.252156193025904),
                (FRAC_PI_2, 3.337069529813158),
                (1.5707594492431998, 1.1105768829555696),
                (1.5705373423867857, 1.1105650645307603),
                (1.570315235103937, 1.1105532551024804),
                (1.5700931273869867, 1.110541454725338),
                (1.5698710192282612, 1.1105296634160235),
                (1.56964891062008, 1.1105178811898386),
                (1.5694268015547566, 1.11050610806312),
                (1.569204692024597, 1.1104943440503179),
                (1.5689825820219025, 1.1104825891677614),
                (1.5687604715389671, 1.1104708434311859),
                (1.5685383605680796, 1.1104591068567051),
                (1.56831624910152, 1.1104473794593046),
                (1.5680941371315633, 1.1104356612552402),
            ],
        },
        EngineZone {
            text: "013232323232323232323",
            id: 0x169a_69a6_9a69_a69f,
            points_at_0: 1,
            extent: [FRAC_PI_2, -2.946115777366428, FRAC_PI_2, -2.946115777366428],
            ring_at_1: &[(FRAC_PI_2, -2.946115777366428)],
        },
    ];
    const IVEA7H_ENGINE: [EngineZone; 8] = [
        EngineZone {
            text: "0010",
            id: 0x023f_ffff_ffff_ffff,
            points_at_0: 120,
            extent: [
                1.0646612280779808,
                0.8542021938068779,
                1.257945516221589,
                1.3167176020555231,
            ],
            ring_at_1: &[
                (1.0646612280779808, 0.9887676584135472),
                (1.07738239745752, 1.1852844734590953),
                (1.1659347120254, 1.3167176020555231),
                (1.257945516221589, 1.2068656477985042),
                (1.244302353530911, 0.9124953566343073),
                (1.1410502336699215, 0.8542021938068779),
            ],
        },
        EngineZone {
            text: "0003",
            id: 0x00ff_ffff_ffff_ffff,
            points_at_0: 120,
            extent: [
                0.7934694379530692,
                0.1954768762233649,
                0.9851790870047552,
                0.49334260778055367,
            ],
            ring_at_1: &[
                (0.8423781916893086, 0.46824326645460007),
                (0.9407979280545962, 0.49334260778055367),
                (0.9851790870047552, 0.35300531591383577),
                (0.9277898733945906, 0.19547687622336496),
                (0.8360928407434417, 0.195476876223365),
                (0.7934694379530692, 0.32673811853074847),
            ],
        },
        EngineZone {
            text: "0313522",
            id: 0x32ea_5fff_ffff_ffff,
            points_at_0: 72,
            extent: [
                0.2354577339632168,
                -0.7742693017922222,
                0.24645330731635065,
                -0.7638118007208565,
            ],
            ring_at_1: &[
                (0.24455692673999357, -0.7727030730187372),
                (0.23905908042860402, -0.7742693017922222),
                (0.2354577339632168, -0.7706152766614258),
                (0.23735474509221613, -0.7653947356192589),
                (0.24285241493439252, -0.7638118007208565),
                (0.24645330731635065, -0.7674661505000296),
            ],
        },
        EngineZone {
            text: "0100000",
            id: 0x1000_1fff_ffff_ffff,
            points_at_0: 61,
            extent: [
                1.0147283313059376,
                -2.954222344162003,
                1.0239642583616924,
                -2.9368571756866926,
            ],
            ring_at_1: &[
                (1.021732316479856, -2.954090122311022),
                (1.0160285004172136, -2.953108961737678),
                (1.0147283313059376, -2.9425443183754383),
                (1.0196106063165031, -2.9368571756866926),
                (1.0239642583616924, -2.9439783086785276),
            ],
        },
        EngineZone {
            text: "0000000",
            id: 0x0000_1fff_ffff_ffff,
            points_at_0: 65,
            extent: [
                1.0147283313059352,
                0.18737093716457637,
                1.0239642583616897,
                0.20473547790309995,
            ],
            ring_at_1: &[
                (1.021732316479853, 0.18750253127876926),
                (1.016028500417211, 0.1884836918521156),
                (1.0147283313059352, 0.1990483352143539),
                (1.0196106063164998, 0.20473547790309995),
                (1.0239642583616897, 0.1976143449112653),
            ],
        },
        EngineZone {
            text: "0100033",
            id: 0x1003_7fff_ffff_ffff,
            points_at_0: 74,
            extent: [
                1.0445636558956015,
                -2.9638784001281806,
                1.0541841485527195,
                -2.9416335891534815,
            ],
            ring_at_1: &[
                (1.0523675719817585, -2.9618135134840053),
                (1.0474814201094913, -2.9638784001281806),
                (1.0445636558956015, -2.954973190741774),
                (1.0462476521583448, -2.943892757590509),
                (1.051002675320016, -2.9416335891534815),
                (1.0541841485527195, -2.950623006269021),
            ],
        },
        EngineZone {
            text: "0132323",
            id: 0x169a_7fff_ffff_ffff,
            points_at_0: 108,
            extent: [
                1.5609057502167343,
                1.1096959569693379,
                FRAC_PI_2,
                -2.0317283258487153,
            ],
            ring_at_1: &[
                (1.5633069233206947, 1.7662732030180608),
                (1.5609057502167343, 2.2898993999368837),
                (1.5615085093325194, 2.857579887885085),
                (1.5641168477605598, 3.4980955335327595),
                (1.5680945199668195, 4.251456981330871),
                (1.5684394352424427, 4.251435464252733),
                (1.5687843412421252, 4.25141391097788),
                (1.569129237986206, 4.2513923215511396),
                (1.569474125495033, 4.25137069601844),
                (1.5698190037889568, 4.251349034424355),
                (1.5701638728883296, 4.251327336814687),
                (1.5705087328135119, 4.251305603235791),
                (FRAC_PI_2, 3.337069529813158),
                (1.5707594578169286, 1.1096971256239985),
                (1.5705374006337343, 1.1097111413054692),
                (1.5703153396557519, 1.1097251420597032),
                (1.570093274877545, 1.1097391279069215),
                (1.5698712062936768, 1.1097530988349025),
                (1.569649133898712, 1.109767054832119),
                (1.5694270576872165, 1.1097809958869678),
                (1.5692049776537544, 1.1097949219864534),
                (1.5689828937928934, 1.1098088331190068),
                (1.5687608060992007, 1.1098227292727254),
                (1.568538714567245, 1.1098366104357758),
                (1.5683166191915938, 1.1098504765954893),
                (1.5680945199668177, 1.109864327740151),
            ],
        },
        EngineZone {
            text: "013232323232323232323",
            id: 0x169a_69a6_9a69_a69f,
            points_at_0: 1,
            extent: [FRAC_PI_2, -2.946115777366428, FRAC_PI_2, -2.946115777366428],
            ring_at_1: &[(FRAC_PI_2, -2.946115777366428)],
        },
    ];
    const RTEA7H_ENGINE: [EngineZone; 8] = [
        EngineZone {
            text: "0010",
            id: 0x023f_ffff_ffff_ffff,
            points_at_0: 120,
            extent: [
                1.0652938011712971,
                0.8508670502271447,
                1.2579503088972461,
                1.316474845330918,
            ],
            ring_at_1: &[
                (1.0652938011712971, 0.9887985174165007),
                (1.0783242763269543, 1.1842900485166463),
                (1.1659570522534408, 1.316474845330918),
                (1.2579503088972461, 1.213835936091201),
                (1.244673895044674, 0.9164762675126693),
                (1.1414815201131843, 0.8508670502271447),
            ],
        },
        EngineZone {
            text: "0003",
            id: 0x00ff_ffff_ffff_ffff,
            points_at_0: 120,
            extent: [
                0.7938389966008388,
                0.19547687622336438,
                0.9856604141025156,
                0.4920547240979091,
            ],
            ring_at_1: &[
                (0.8424607668643684, 0.4676457025401466),
                (0.9412757285794003, 0.4920547240979091),
                (0.9856604141025156, 0.3511861351427969),
                (0.9289060795711012, 0.19547687622336443),
                (0.8369911077926054, 0.19547687622336454),
                (0.7938389966008388, 0.32698144538942225),
            ],
        },
        EngineZone {
            text: "0313522",
            id: 0x32ea_5fff_ffff_ffff,
            points_at_0: 72,
            extent: [
                0.23469886401119266,
                -0.7744510107752828,
                0.24596605938655264,
                -0.7642600673316566,
            ],
            ring_at_1: &[
                (0.24392299011335727, -0.7729924601168706),
                (0.23830546190987906, -0.7744510107752828),
                (0.23469886401119266, -0.7708041674063861),
                (0.2367165522974506, -0.7657067082025604),
                (0.2423653440530303, -0.7642600673316566),
                (0.24596605938655264, -0.7678990862681049),
            ],
        },
        EngineZone {
            text: "0100000",
            id: 0x1000_1fff_ffff_ffff,
            points_at_0: 61,
            extent: [
                1.01472032548627,
                -2.95446601807496,
                1.0240027528277493,
                -2.9368100447100822,
            ],
            ring_at_1: &[
                (1.0217062408209738, -2.9541831096477087),
                (1.0159742200618251, -2.9530905014579005),
                (1.01472032548627, -2.9424410799161973),
                (1.0196599044409143, -2.9368100447100822),
                (1.0240027528277493, -2.944054151057197),
            ],
        },
        EngineZone {
            text: "0000000",
            id: 0x0000_1fff_ffff_ffff,
            points_at_0: 65,
            extent: [
                1.0147203254862676,
                0.18707001229289083,
                1.0240027528277467,
                0.2047826088797105,
            ],
            ring_at_1: &[
                (1.0217062408209714, 0.18740954394208334),
                (1.0159742200618225, 0.18850215213189356),
                (1.0147203254862676, 0.1991515736735955),
                (1.0196599044409114, 0.2047826088797105),
                (1.0240027528277467, 0.19753850253259667),
            ],
        },
        EngineZone {
            text: "0100033",
            id: 0x1003_7fff_ffff_ffff,
            points_at_0: 74,
            extent: [
                1.0449079245240926,
                -2.9633294184475467,
                1.0551401397620093,
                -2.9417770283039313,
            ],
            ring_at_1: &[
                (1.0525924409439182, -2.9613228446197244),
                (1.0474718061779824, -2.9633294184475467),
                (1.0449079245240926, -2.9546815343409505),
                (1.0470905141952849, -2.943965644109549),
                (1.0518503888610276, -2.9417770283039313),
                (1.0551401397620093, -2.9504815409896934),
            ],
        },
        EngineZone {
            text: "0132323",
            id: 0x169a_7fff_ffff_ffff,
            points_at_0: 108,
            extent: [
                1.5608035282624262,
                1.112264393763418,
                FRAC_PI_2,
                -2.0293282532154304,
            ],
            ring_at_1: &[
                (1.5634244201932699, 1.7662732030180575),
                (1.5608035282624262, 2.2909308628855314),
                (1.5614576749773215, 2.846074452726508),
                (1.5642321327679056, 3.506983652267898),
                (1.5680592920758163, 4.253857047390701),
                (1.5684086727986106, 4.25385704739085),
                (1.56875805344265, 4.253857047390955),
                (1.569107434019456, 4.253857047390753),
                (1.5694568145405547, 4.253857047391051),
                (1.5698061950174715, 4.253857047391335),
                (1.5701555754617278, 4.253857047392199),
                (1.5705049558848492, 4.253857047395218),
                (FRAC_PI_2, 3.337069529813158),
                (1.5707589732622798, 1.1122643937721473),
                (1.5705339999062164, 1.112264393798666),
                (1.570309026546565, 1.1122643937990364),
                (1.5700840531802487, 1.1122643937992416),
                (1.5698590798041911, 1.1122643937993222),
                (1.569634106415315, 1.1122643937993504),
                (1.5694091330105446, 1.1122643937999532),
                (1.5691841595868012, 1.112264393799931),
                (1.5689591861410093, 1.1122643937998746),
                (1.5687342126700914, 1.1122643937998302),
                (1.5685092391709712, 1.1122643938001597),
                (1.5682842656405709, 1.112264393800108),
                (1.5680592920758143, 1.1122643938000634),
            ],
        },
        EngineZone {
            text: "013232323232323232323",
            id: 0x169a_69a6_9a69_a69f,
            points_at_0: 1,
            extent: [FRAC_PI_2, -2.946115777366428, FRAC_PI_2, -2.946115777366428],
            ring_at_1: &[(FRAC_PI_2, -2.946115777366428)],
        },
    ];

    // The engine's answers on the three aperture-3 grids: of two zones that pin site GR-1, a
    // hexagon and a pentagon of level 2, and then of one zone of each kind.
    const ISEA3H_ENGINE: [EngineZone; 15] = [
        EngineZone {
            text: "B4-4-A",
            id: 0x0280_0000_0000_0010,
            points_at_0: 120,
            extent: [
                -0.25322127440706776,
                -0.1787194949679124,
                0.25322127440706477,
                0.19547687622336501,
            ],
            ring_at_1: &[
                (0.25322127440706477, 0.004370231930614905),
                (0.13024648689955778, -0.1776358979876636),
                (-0.13024648689956064, -0.17763589798766377),
                (-0.25322127440706776, 0.004370231930614743),
                (-0.12291855091939194, 0.19547687622336496),
                (0.1229185509193889, 0.19547687622336496),
            ],
        },
        EngineZone {
            text: "B4-0-A",
            id: 0x0280_0000_0000_0000,
            points_at_0: 100,
            extent: [
                -0.21170005999557479,
                -0.5799031971730955,
                0.21170005999557218,
                -0.1776358979876629,
            ],
            ring_at_1: &[
                (0.21170005999557218, -0.42767333555827247),
                (-1.1404657533354003e-15, -0.5799031971730955),
                (-0.21170005999557479, -0.4276733355582727),
                (-0.13024648689956111, -0.17763589798766302),
                (0.13024648689955826, -0.1776358979876629),
            ],
        },
        EngineZone {
            text: "C4-49-B",
            id: 0x0480_0000_0000_0125,
            points_at_0: 60,
            extent: [
                -0.839337269837361,
                0.12740730179486015,
                -0.761278867168083,
                0.26354645065186977,
            ],
            ring_at_1: &[
                (-0.7663554681913034, 0.16262437968791152),
                (-0.8066541662138593, 0.12740730179486015),
                (-0.839337269837361, 0.1603347563757197),
                (-0.8393372698373609, 0.23061899607101033),
                (-0.8066541662138595, 0.26354645065186977),
                (-0.7663554681913036, 0.22832937275881846),
            ],
        },
        EngineZone {
            text: "C0-0-B",
            id: 0x0400_0000_0000_0001,
            points_at_0: 50,
            extent: [
                0.9818240465680356,
                -3.018394191713056,
                1.058356474923492,
                -2.873837363019799,
            ],
            ring_at_1: &[
                (1.058356474923492, -2.946115777366428),
                (1.0301695008222396, -3.018394191713056),
                (0.9871479190776499, -2.987818812248788),
                (0.9871479190776499, -2.904412742484068),
                (1.0301695008222396, -2.873837363019799),
            ],
        },
        EngineZone {
            text: "CA-0-B",
            id: 0x0540_0000_0000_0001,
            points_at_0: 50,
            extent: [
                0.9818240465680332,
                0.12319846187673542,
                1.0583564749234884,
                0.26775529056999503,
            ],
            ring_at_1: &[
                (1.0583564749234884, 0.1954768762233649),
                (1.0301695008222367, 0.12319846187673542),
                (0.9871479190776469, 0.15377384134100527),
                (0.9871479190776469, 0.23717991110572456),
                (1.0301695008222367, 0.26775529056999503),
            ],
        },
        EngineZone {
            text: "CB-0-A",
            id: 0x0560_0000_0000_0000,
            points_at_0: 80,
            extent: [
                -1.0783132479880575,
                -3.0759086867160654,
                -0.9446353185123617,
                -2.8163228680167913,
            ],
            ring_at_1: &[
                (-0.9922793780831645, -3.0759086867160654),
                (-1.0777861725813531, -3.03863115083069),
                (-1.0777861725813531, -2.8536004039021665),
                (-0.9922793780831644, -2.8163228680167913),
                (-0.9446353185123617, -2.9461157773664284),
            ],
        },
        EngineZone {
            text: "C0-6-A",
            id: 0x0400_0000_0000_0018,
            points_at_0: 92,
            extent: [
                1.3169102797676526,
                -0.2180434975973176,
                1.4453724243299435,
                0.608997250044049,
            ],
            ring_at_1: &[
                (1.3169102797676526, 0.3613875016988786),
                (1.3646144623791008, 0.608997250044049),
                (1.4384300521252502, 0.5106899169083345),
                (1.4384300521252507, -0.11973616446160565),
                (1.3646144623791003, -0.2180434975973176),
                (1.3169102797676526, 0.02956625074785149),
            ],
        },
        EngineZone {
            text: "E0-7A-A",
            id: 0x0800_0000_0000_01e8,
            points_at_0: 200,
            extent: [
                1.548145618941151,
                -2.031610648835633,
                1.5662624843353457,
                -0.7190282523072226,
            ],
            ring_at_1: &[
                (1.5593531313536695, -0.7190282523072226),
                (1.5662624843353457, -1.3753194505711999),
                (1.5660230669341808, -1.4484969041765752),
                (1.565760596093337, -1.514372918146923),
                (1.5654784818073826, -1.5734972837029835),
                (1.5651796809553176, -1.6265032088613978),
                (1.564866713378447, -1.674038622122016),
                (1.56454170324212, -1.716724682019766),
                (1.564206430114049, -1.755133818940524),
                (1.5638623806592642, -1.789780558086172),
                (1.5635107963110786, -1.8211200658909834),
                (1.5631527150170292, -1.8495509760547424),
                (1.5627890066856862, -1.8754203045800228),
                (1.5624204027228623, -1.8990291388740415),
                (1.5620475203646031, -1.9206383583856912),
                (1.5616708825942462, -1.940474000543992),
                (1.5612909343933683, -1.9587320983258814),
                (1.560908055989569, -1.9755829370634268),
                (1.5605225736634691, -1.991174743414816),
                (1.5601347685803775, -2.005636851479176),
                (1.559744884026123, -2.0190824039774604),
                (1.5593531313536713, -2.031610648835633),
                (1.5513711293993275, -1.7433280211026314),
                (1.548145618941151, -1.3753194505714652),
                (1.5513711293993264, -1.007310880040285),
            ],
        },
        EngineZone {
            text: "C0-4-C",
            id: 0x0400_0000_0000_0012,
            points_at_0: 251,
            extent: [
                1.4867442066092442,
                -2.9461157773664284,
                FRAC_PI_2,
                0.19547687622336474,
            ],
            ring_at_1: &[
                (1.512252049057587, -2.1756995352715767),
                (1.4867442066092451, -1.628975367331603),
                (1.4867442066092442, -1.1216635338114256),
                (1.5122520490575848, -0.5749393658714503),
                (1.5498672963521878, 0.19547687517012768),
                (1.5519600767686321, 0.1954768750531528),
                (1.5540528959079065, 0.19547687490691956),
                (1.5561457494686433, 0.1954768747189075),
                (1.5582386331490616, 0.1954768744682145),
                (1.56033154264702, 0.1954768741172399),
                (1.5624244736600676, 0.19547687359075328),
                (1.564517421885496, 0.19547687271326808),
                (1.5666103830203908, 0.19547687095828592),
                (1.568703352761682, 0.1954768656932413),
                (FRAC_PI_2, 0.19547687622336474),
                (FRAC_PI_2, -2.9461157773664284),
                (1.5687033527390806, -2.9461157668363045),
                (1.5666103829977898, -2.946115772101349),
                (1.5645174218628952, -2.9461157738563313),
                (1.5624244736374668, -2.9461157747338165),
                (1.5603315426244193, -2.9461157752602993),
                (1.558238633126461, -2.9461157756112777),
                (1.5561457494460431, -2.9461157758619705),
                (1.5540528958853068, -2.946115776049985),
                (1.5519600767460326, -2.946115776196216),
                (1.549867296329589, -2.946115776313191),
            ],
        },
        EngineZone {
            text: "F0-79-A",
            id: 0x0a00_0000_0000_01e4,
            points_at_0: 64,
            extent: [
                1.5659052175735526,
                1.7662732029685255,
                FRAC_PI_2,
                -1.3753194505217956,
            ],
            ring_at_1: &[
                (1.5692848356260396, -4.516912104162319),
                (1.5669821037428857, -3.861016032843341),
                (1.5659052175735526, -3.260289857076885),
                (1.5659052175735526, -2.631941697656023),
                (1.5669821037428853, -2.031215521889581),
                (1.5692848356260398, -1.3753194505705368),
                (1.5694359752450657, -1.3753194505704263),
                (1.5695871169720033, -1.375319450570288),
                (1.5697382608078738, -1.3753194505701105),
                (1.5698894067536984, -1.3753194505698736),
                (1.5700405548104988, -1.375319450569542),
                (1.5701917049792966, -1.3753194505690447),
                (1.5703428572611133, -1.3753194505682158),
                (1.5704940116569708, -1.3753194505665578),
                (1.5706451681678912, -1.3753194505615842),
                (FRAC_PI_2, -2.9461157859242113),
                (FRAC_PI_2, -2.9461157688086423),
                (FRAC_PI_2, -2.946115785924211),
                (FRAC_PI_2, -2.9461157688086423),
                (FRAC_PI_2, -2.946115785924211),
                (1.570645168167891, -4.516912104171272),
                (1.5704940116569701, -4.516912104166298),
                (1.5703428572611133, -4.516912104164641),
                (1.5701917049792964, -4.516912104163811),
                (1.5700405548104983, -4.516912104163314),
                (1.5698894067536984, -4.516912104162984),
                (1.5697382608078734, -4.516912104162746),
                (1.5695871169720026, -4.516912104162568),
                (1.5694359752450657, -4.516912104162429),
            ],
        },
        EngineZone {
            text: "M0-40DF8-A",
            id: 0x1800_0000_0010_37e0,
            points_at_0: 62,
            extent: [
                1.570794090614931,
                1.766273148632095,
                FRAC_PI_2,
                -1.3753193961853651,
            ],
            ring_at_1: &[
                (1.5707956356212318, -4.516912106336772),
                (1.5707945827825367, -3.861217654954196),
                (1.5707940906149311, -3.2603484204031825),
                (1.570794090614931, -2.631883134207023),
                (1.5707945827825363, -2.0310138995708935),
                (1.5707956356212314, -1.3753194483960847),
                (1.570795704738596, -1.3753194481543685),
                (1.5707957738559608, -1.375319447852223),
                (1.5707958429733262, -1.3753194474637505),
                (1.5707959120906922, -1.375319446945787),
                (1.5707959812080585, -1.375319446220638),
                (1.5707960503254252, -1.3753194451329147),
                (1.5707961194427922, -1.3753194433200426),
                (1.57079618856016, -1.375319439694298),
                (1.570796257677528, -1.3753194288170647),
                (FRAC_PI_2, -2.94611577736695),
                (1.5707962576775285, -4.516912125915791),
                (1.5707961885601598, -4.516912115038558),
                (1.5707961194427924, -4.516912111412813),
                (1.5707960503254257, -4.516912109599941),
                (1.5707959812080583, -4.5169121085122175),
                (1.5707959120906922, -4.5169121077870695),
                (1.5707958429733266, -4.516912107269105),
                (1.5707957738559617, -4.516912106880634),
                (1.570795704738596, -4.516912106578488),
            ],
        },
        EngineZone {
            text: "Q0-1486BA0-C",
            id: 0x2000_0000_0521_ae82,
            points_at_0: 248,
            extent: [
                1.570796309176799,
                -1.6263214253769345,
                FRAC_PI_2,
                0.1954772282041839,
            ],
            ring_at_1: &[
                (FRAC_PI_2, 0.19547670023293187),
                (FRAC_PI_2, 0.19547705221373413),
                (1.5707963091767996, -1.6263214253769345),
                (1.570796309176799, -1.1243173134488027),
                (FRAC_PI_2, 0.1954767002329961),
                (FRAC_PI_2, 0.19547705221379807),
                (FRAC_PI_2, 0.1954768757688144),
                (FRAC_PI_2, 0.1954768766779158),
                (FRAC_PI_2, 0.19547687576881356),
                (FRAC_PI_2, 0.1954768766779158),
                (FRAC_PI_2, 0.1954768757688136),
                (FRAC_PI_2, 0.1954768766779158),
                (FRAC_PI_2, 0.1954768757688144),
                (FRAC_PI_2, 0.1954768766779158),
                (FRAC_PI_2, 0.19547687576881437),
                (FRAC_PI_2, 0.1954768766779167),
                (FRAC_PI_2, 0.19547687576881437),
                (FRAC_PI_2, 0.1954768766779158),
                (FRAC_PI_2, 0.19547687576881356),
                (FRAC_PI_2, 0.1954768766779167),
                (FRAC_PI_2, 0.19547687576881434),
                (FRAC_PI_2, 0.1954768766779158),
                (FRAC_PI_2, 0.19547687576881353),
                (FRAC_PI_2, 0.1954768766779167),
                (FRAC_PI_2, 0.19547687576881437),
                (FRAC_PI_2, 0.1954768766779158),
                (FRAC_PI_2, 0.19547687576881434),
                (FRAC_PI_2, 0.1954768766779158),
                (FRAC_PI_2, 0.19547687576881434),
                (FRAC_PI_2, 0.1954768766779158),
                (FRAC_PI_2, 0.19547687576881434),
                (FRAC_PI_2, 0.1954768766779158),
                (FRAC_PI_2, 0.1954768757688143),
                (FRAC_PI_2, 0.1954768766779158),
                (FRAC_PI_2, 0.19547687576881434),
                (FRAC_PI_2, 0.1954768766779158),
                (FRAC_PI_2, 0.1954768757688143),
                (FRAC_PI_2, 0.1954768766779158),
                (FRAC_PI_2, 0.1954768757688143),
                (FRAC_PI_2, 0.1954768766779158),
                (FRAC_PI_2, 0.1954768757688143),
                (FRAC_PI_2, 0.1954768766779158),
                (FRAC_PI_2, 0.1954768757688143),
                (FRAC_PI_2, 0.19547687667791624),
                (FRAC_PI_2, 0.1954768757688143),
                (FRAC_PI_2, 0.1954768766779158),
                (FRAC_PI_2, 0.19547687576881434),
                (FRAC_PI_2, 0.1954768766779158),
            ],
        },
        EngineZone {
            text: "Q8-34AA7F2578AE0-C",
            id: 0x210d_2a9f_c95e_2b82,
            points_at_0: 47,
            extent: [
                1.5707963091767996,
                0.1954765242425447,
                FRAC_PI_2,
                2.017275161543644,
            ],
            ring_at_1: &[
                (1.5707963091767996, 2.017275161543644),
                (FRAC_PI_2, 0.19547670023300734),
                (FRAC_PI_2, 0.19547705221378608),
                (FRAC_PI_2, 0.19547687622336696),
                (FRAC_PI_2, 0.19547670023294383),
                (FRAC_PI_2, 0.1954770522137228),
                (1.5707963091823995, 1.516514247647949),
            ],
        },
        EngineZone {
            text: "A0-0-A",
            id: 0x0000_0000_0000_0000,
            points_at_0: 139,
            extent: [
                0.36635950360263975,
                1.7662732030174428,
                FRAC_PI_2,
                -1.3753194505707127,
            ],
            ring_at_1: &[
                (1.2074224434760321, -4.516912104161329),
                (0.6175926646873967, -3.731513940763877),
                (0.36635950360263975, -2.9461157773664284),
                (0.6175926646873967, -2.160717613968979),
                (1.207422443476032, -1.3753194505715276),
                (1.2435838967460207, -1.375319450571527),
                (1.2797501142146526, -1.3753194505715265),
                (1.3159337228545505, -1.3753194505715258),
                (1.3521474480612032, -1.3753194505715247),
                (1.3884041486830376, -1.3753194505715234),
                (1.4247168529542358, -1.3753194505715214),
                (1.4610987956473578, -1.3753194505715178),
                (1.4975634567828358, -1.375319450571511),
                (1.5341246022584634, -1.3753194505714905),
                (FRAC_PI_2, -2.9461157773664284),
                (1.5341246022584638, -4.516912104161366),
                (1.4975634567828362, -4.516912104161345),
                (1.4610987956473585, -4.516912104161339),
                (1.4247168529542367, -4.516912104161335),
                (1.3884041486830385, -4.5169121041613325),
                (1.3521474480612035, -4.516912104161332),
                (1.3159337228545507, -4.516912104161331),
                (1.279750114214653, -4.51691210416133),
                (1.2435838967460209, -4.516912104161329),
            ],
        },
        EngineZone {
            text: "Q5-34AA7F2578AE1-A",
            id: 0x20ad_2a9f_c95e_2b84,
            points_at_0: 130,
            extent: [
                -FRAC_PI_2,
                2.366015781344457,
                -1.5707962991877369,
                -1.9750620395448706,
            ],
            ring_at_1: &[
                (-FRAC_PI_2, -2.946115601376036),
                (-FRAC_PI_2, -2.9461159533568204),
                (-FRAC_PI_2, -2.946115618975079),
                (-FRAC_PI_2, -2.9461159357577773),
                (-FRAC_PI_2, -2.9461156365741115),
                (-FRAC_PI_2, -2.9461159181587453),
                (-FRAC_PI_2, -2.946115654173154),
                (-FRAC_PI_2, -2.946115900559702),
                (-FRAC_PI_2, -2.9461156717721977),
                (-FRAC_PI_2, -2.9461158829606586),
                (-FRAC_PI_2, -2.946115689371231),
                (-FRAC_PI_2, -2.9461158653616257),
                (-FRAC_PI_2, -2.9461157069702737),
                (-FRAC_PI_2, -2.9461158477625826),
                (-FRAC_PI_2, -2.9461157245693075),
                (-FRAC_PI_2, -2.946115830163549),
                (-FRAC_PI_2, -2.9461157421683506),
                (-FRAC_PI_2, -2.9461158125645057),
                (-FRAC_PI_2, -2.946115759767394),
                (-FRAC_PI_2, -2.946115794965462),
                (-FRAC_PI_2, -2.9461157773664275),
                (-FRAC_PI_2, -2.946115759767386),
                (-FRAC_PI_2, -2.94611579496547),
                (-FRAC_PI_2, -2.9461157421683417),
                (-FRAC_PI_2, -2.946115812564514),
                (-FRAC_PI_2, -2.9461157245693093),
                (-FRAC_PI_2, -2.946115830163547),
                (-FRAC_PI_2, -2.9461157069702657),
                (-FRAC_PI_2, -2.946115847762591),
                (-FRAC_PI_2, -2.9461156893712235),
                (-FRAC_PI_2, -2.946115865361633),
                (-FRAC_PI_2, -2.9461156717721897),
                (-FRAC_PI_2, -2.9461158829606666),
                (-FRAC_PI_2, -2.9461156541731475),
                (-FRAC_PI_2, -2.9461159005597097),
                (-FRAC_PI_2, -2.946115636574103),
                (-FRAC_PI_2, -2.9461159181587533),
                (-FRAC_PI_2, -2.946115618975071),
                (-FRAC_PI_2, -2.9461159357577857),
                (-FRAC_PI_2, -2.946115601376028),
                (-FRAC_PI_2, -2.9461159533568284),
                (-1.5707963052638791, -2.0310138809662517),
                (-1.5707962991877369, -2.6318831186680898),
                (-1.5707962994394413, -2.661553145408573),
                (-1.5707962996666296, -2.6917470802550088),
                (-1.5707962998686824, -2.7224225126548514),
                (-1.5707963000450287, -2.7535320171275894),
                (-1.570796300195158, -2.785023454245568),
                (-1.5707963003186245, -2.8168404398881504),
                (-1.5707963004150536, -2.8489227320656894),
                (-1.5707963004841476, -2.881206910109472),
                (-1.570796300525691, -2.9136271235956253),
                (-1.5707963005395535, -2.946115778711427),
                (-1.5707963005256913, -2.9786044150800928),
                (-1.5707963004841479, -3.011024631087399),
                (-1.5707963004150531, -3.043308807406452),
                (-1.5707963003186245, -3.075391101051478),
                (-1.5707963001951584, -3.107208086234565),
                (-1.5707963000450285, -3.1386995205882897),
                (-1.5707962998686824, -3.169809025714121),
                (-1.57079629966663, -3.200484462987536),
                (-1.5707962994394409, -3.2306783898810294),
                (-1.570796299187737, -3.260348424932014),
                (-1.5707963052638791, -3.8612176870697366),
            ],
        },
    ];
    const IVEA3H_ENGINE: [EngineZone; 15] = [
        EngineZone {
            text: "B4-4-A",
            id: 0x0280_0000_0000_0010,
            points_at_0: 120,
            extent: [
                -0.248923076827407,
                -0.18403351866790815,
                0.248923076827404,
                0.195476876223365,
            ],
            ring_at_1: &[
                (0.248923076827404, -0.0023344543308559994),
                (0.12572552476898874, -0.184033518667908),
                (-0.1257255247689916, -0.18403351866790815),
                (-0.248923076827407, -0.002334454330856166),
                (-0.12328000274069219, 0.19547687622336476),
                (0.12328000274068945, 0.19547687622336476),
            ],
        },
        EngineZone {
            text: "B4-0-A",
            id: 0x0280_0000_0000_0000,
            points_at_0: 100,
            extent: [
                -0.20428760985147687,
                -0.5721230160797547,
                0.20428760985147415,
                -0.1756669246730134,
            ],
            ring_at_1: &[
                (0.20428760985147415, -0.42516298502123356),
                (-1.1484140412133812e-15, -0.5721230160797547),
                (-0.20428760985147687, -0.4251629850212339),
                (-0.12572552476899207, -0.18403351866790757),
                (0.12572552476898918, -0.1840335186679074),
            ],
        },
        EngineZone {
            text: "C4-49-B",
            id: 0x0480_0000_0000_0125,
            points_at_0: 60,
            extent: [
                -0.8420906198265579,
                0.1261609349726837,
                -0.7695902523587982,
                0.2647928174740465,
            ],
            ring_at_1: &[
                (-0.7703615172362234, 0.16209198888816745),
                (-0.8072438085312589, 0.1261609349726837),
                (-0.8420906198265579, 0.1593872742661208),
                (-0.8420906198265579, 0.23156647818060944),
                (-0.8072438085312589, 0.2647928174740465),
                (-0.7703615172362233, 0.22886176355856272),
            ],
        },
        EngineZone {
            text: "C0-0-B",
            id: 0x0400_0000_0000_0001,
            points_at_0: 50,
            extent: [
                0.9836841516899438,
                -3.02077287852043,
                1.0596213568111716,
                -2.871458676212425,
            ],
            ring_at_1: &[
                (1.0596213568111716, -2.946115777366428),
                (1.0304842741073945, -3.02077287852043),
                (0.9860981484276611, -2.9890988782570203),
                (0.9860981484276617, -2.9031326764758365),
                (1.0304842741073947, -2.871458676212425),
            ],
        },
        EngineZone {
            text: "CA-0-B",
            id: 0x0540_0000_0000_0001,
            points_at_0: 50,
            extent: [
                0.9836841516899425,
                0.12081977506936153,
                1.059621356811168,
                0.270133977377369,
            ],
            ring_at_1: &[
                (1.059621356811168, 0.1954768762233649),
                (1.0304842741073914, 0.12081977506936153),
                (0.9860981484276589, 0.15249377533277364),
                (0.9860981484276587, 0.23845997711395622),
                (1.0304842741073919, 0.270133977377369),
            ],
        },
        EngineZone {
            text: "CB-0-A",
            id: 0x0560_0000_0000_0000,
            points_at_0: 80,
            extent: [
                -1.0798212532740166,
                -3.0700896145240244,
                -0.9481241136953423,
                -2.822141940208832,
            ],
            ring_at_1: &[
                (-0.9937054235302549, -3.0700896145240244),
                (-1.0751349226589186, -3.0338680713221136),
                (-1.075134922658919, -2.8583634834107423),
                (-0.993705423530255, -2.822141940208832),
                (-0.9481241136953423, -2.9461157773664284),
            ],
        },
        EngineZone {
            text: "C0-6-A",
            id: 0x0400_0000_0000_0018,
            points_at_0: 92,
            extent: [
                1.3192222619887986,
                -0.22537123461814296,
                1.4461240889731666,
                0.6163249870648754,
            ],
            ring_at_1: &[
                (1.3192222619887992, 0.3644919528857882),
                (1.3658248375902748, 0.6163249870648754),
                (1.4392916616077007, 0.5167111785303711),
                (1.4392916616077005, -0.1257574260836416),
                (1.3658248375902742, -0.22537123461814296),
                (1.3192222619887986, 0.026461799560942106),
            ],
        },
        EngineZone {
            text: "E0-7A-A",
            id: 0x0800_0000_0000_01e8,
            points_at_0: 200,
            extent: [
                1.548143962118078,
                -2.0311953128393427,
                1.5662656880295158,
                -0.7194435883035099,
            ],
            ring_at_1: &[
                (1.5593513853200704, -0.7194435883035099),
                (1.5662656880295158, -1.3753194505711868),
                (1.5660261253829983, -1.4485774129497728),
                (1.5657634543320538, -1.5145100551108872),
                (1.5654811001375497, -1.5736698834615834),
                (1.5651820321020344, -1.6266932758418813),
                (1.5648687798616023, -1.6742312837814395),
                (1.5645434751404252, -1.7169078848909067),
                (1.5642079032569263, -1.7552979261969195),
                (1.5638635552104672, -1.7899179387004152),
                (1.5635116756812664, -1.8212247195557452),
                (1.5631533050392592, -1.8496182107837564),
                (1.5627893149948586, -1.8754464693302038),
                (1.5624204382905071, -1.899011406653687),
                (1.5620472931506981, -1.9205745527319122),
                (1.5616704032866575, -1.9403624578676184),
                (1.5612902142127068, -1.9585715592444377),
                (1.5609071065428488, -1.975372460830291),
                (1.5605214068341948, -1.9909136405912546),
                (1.560133396445875, -2.0053246309440516),
                (1.5597433187952907, -2.0187187311772035),
                (1.5593513853200724, -2.0311953128393427),
                (1.551359408801181, -1.7426582593092472),
                (1.548143962118078, -1.3753194505714652),
                (1.55135940880118, -1.0079806418336719),
            ],
        },
        EngineZone {
            text: "C0-4-C",
            id: 0x0400_0000_0000_0012,
            points_at_0: 251,
            extent: [
                1.4863964494896615,
                -2.9461157773664284,
                FRAC_PI_2,
                0.19547687622336474,
            ],
            ring_at_1: &[
                (1.5122293439935903, -2.170265298066484),
                (1.4863964494896624, -1.6264848978945448),
                (1.4863964494896615, -1.1241540032484838),
                (1.5122293439935883, -0.5803736030765426),
                (1.5498782342213222, 0.1954768762233603),
                (1.5519673118865007, 0.19547687622335985),
                (1.5540569885804607, 0.19547687622335919),
                (1.5561472672877845, 0.1954768762233583),
                (1.558238151000091, 0.1954768762233572),
                (1.5603296427160753, 0.19547687622335563),
                (1.5624217454415525, 0.1954768762233532),
                (1.564514462189499, 0.19547687622334942),
                (1.566607795980093, 0.19547687622334142),
                (1.5687017498407585, 0.19547687622331789),
                (FRAC_PI_2, 0.19547687622336474),
                (FRAC_PI_2, -2.9461157773664284),
                (1.5687017498181435, -2.946115777366381),
                (1.566607795957485, -2.9461157773664044),
                (1.5645144621668976, -2.9461157773664124),
                (1.562421745418958, -2.9461157773664164),
                (1.560329642693487, -2.9461157773664186),
                (1.5582381509775092, -2.9461157773664204),
                (1.5561472672652092, -2.9461157773664213),
                (1.554056988557892, -2.946115777366422),
                (1.5519673118639383, -2.946115777366423),
                (1.5498782341987662, -2.9461157773664235),
            ],
        },
        EngineZone {
            text: "F0-79-A",
            id: 0x0a00_0000_0000_01e4,
            points_at_0: 64,
            extent: [
                1.5659033323123248,
                1.7662732029511652,
                FRAC_PI_2,
                -1.3753194505180908,
            ],
            ring_at_1: &[
                (1.5692861118256802, -4.516912104161975),
                (1.5669825001618098, -3.860572349812576),
                (1.5659033323123248, -3.260103295754979),
                (1.565903332312325, -2.632128258977972),
                (1.5669825001618098, -2.031659204920374),
                (1.569286111825681, -1.375319450570459),
                (1.5694371332788388, -1.3753194505703825),
                (1.5695881547458168, -1.375319450570191),
                (1.5697391762250792, -1.3753194505699993),
                (1.5698901977150908, -1.3753194505697501),
                (1.5700412192143156, -1.3753194505692474),
                (1.5701922407212185, -1.3753194505688595),
                (1.5703432622342635, -1.375319450567969),
                (1.5704942837519151, -1.3753194505658215),
                (1.5706453052726381, -1.3753194505608433),
                (FRAC_PI_2, -2.9461157859163896),
                (FRAC_PI_2, -2.9461157688164645),
                (FRAC_PI_2, -2.9461157859163896),
                (FRAC_PI_2, -2.9461157688164645),
                (FRAC_PI_2, -2.9461157859163896),
                (1.570645305272638, -4.516912104171281),
                (1.5704942837519147, -4.516912104166304),
                (1.5703432622342632, -4.5169121041636675),
                (1.5701922407212179, -4.516912104163081),
                (1.5700412192143147, -4.516912104162731),
                (1.569890197715091, -4.516912104162984),
                (1.5697391762250792, -4.516912104162747),
                (1.569588154745816, -4.516912104161655),
                (1.5694371332788386, -4.516912104162046),
            ],
        },
        EngineZone {
            text: "M0-40DF8-A",
            id: 0x1800_0000_0010_37e0,
            points_at_0: 62,
            extent: [
                1.5707940889598995,
                1.7662731485823393,
                FRAC_PI_2,
                -1.3753193921341738,
            ],
            ring_at_1: &[
                (1.570795636252956, -4.516912107139049),
                (1.5707945831900982, -3.860333808925805),
                (1.5707940889598995, -3.2598113787361727),
                (1.5707940889598995, -2.632420175970583),
                (1.5707945831900982, -2.0318977456518246),
                (1.5707956362529554, -1.3753194483940947),
                (1.5707957053071495, -1.3753194479743156),
                (1.5707957743613434, -1.3753194478497355),
                (1.5707958434155376, -1.375319447689561),
                (1.5707959124697317, -1.37531944694247),
                (1.570795981523926, -1.3753194465367726),
                (1.5707960505781202, -1.3753194451279394),
                (1.5707961196323144, -1.3753194433134086),
                (1.5707961886865083, -1.37531943888406),
                (1.5707962577407024, -1.3753194303977372),
                (FRAC_PI_2, -2.9461157773669497),
                (1.5707962577407029, -4.5169121259356935),
                (1.5707961886865078, -4.516912111847359),
                (1.5707961196323146, -4.516912113553547),
                (1.5707960505781207, -4.516912111205491),
                (1.5707959815239259, -4.516912108516198),
                (1.570795912469732, -4.516912107790386),
                (1.570795843415538, -4.516912107271949),
                (1.5707957743613443, -4.516912108883838),
                (1.5707957053071497, -4.516912107292065),
            ],
        },
        EngineZone {
            text: "Q0-1486BA0-C",
            id: 0x2000_0000_0521_ae82,
            points_at_0: 47,
            extent: [
                1.5707963091909127,
                -1.62676179010905,
                FRAC_PI_2,
                0.19547723079897805,
            ],
            ring_at_1: &[
                (FRAC_PI_2, 0.19547669893553396),
                (FRAC_PI_2, 0.19547705351113365),
                (1.5707963091909136, -1.62676179010905),
                (1.5707963091909127, -1.1238769532924773),
                (FRAC_PI_2, 0.19547669893559722),
                (FRAC_PI_2, 0.19547705351119582),
                (FRAC_PI_2, 0.19547687622336304),
            ],
        },
        EngineZone {
            text: "Q8-34AA7F2578AE0-C",
            id: 0x210d_2a9f_c95e_2b82,
            points_at_0: 47,
            extent: [
                1.5707963091909132,
                0.19547652164774032,
                FRAC_PI_2,
                2.0177155240754203,
            ],
            ring_at_1: &[
                (1.5707963091909132, 2.0177155240754203),
                (FRAC_PI_2, 0.19547669893560293),
                (FRAC_PI_2, 0.19547705351118336),
                (FRAC_PI_2, 0.19547687622336696),
                (FRAC_PI_2, 0.19547669893553987),
                (FRAC_PI_2, 0.19547705351112055),
                (1.5707963091965282, 1.5160759064131395),
            ],
        },
        EngineZone {
            text: "A0-0-A",
            id: 0x0000_0000_0000_0000,
            points_at_0: 139,
            extent: [
                0.36635950360263975,
                1.766273203017442,
                FRAC_PI_2,
                -1.375319450570712,
            ],
            ring_at_1: &[
                (1.2074224434760321, -4.516912104161329),
                (0.6175926646873967, -3.731513940763877),
                (0.36635950360263975, -2.9461157773664284),
                (0.6175926646873967, -2.160717613968979),
                (1.207422443476032, -1.3753194505715276),
                (1.243151499054315, -1.3753194505715267),
                (1.279069621142875, -1.3753194505715265),
                (1.3151574693548627, -1.3753194505715256),
                (1.3513951321611892, -1.3753194505715247),
                (1.387762187701894, -1.3753194505715232),
                (1.4242377686673735, -1.3753194505715214),
                (1.4608006307823067, -1.3753194505715178),
                (1.4974292243744907, -1.375319450571511),
                (1.5341017684668892, -1.3753194505714905),
                (FRAC_PI_2, -2.9461157773664284),
                (1.5341017684668896, -4.516912104161366),
                (1.497429224374491, -4.516912104161345),
                (1.4608006307823074, -4.516912104161337),
                (1.4242377686673742, -4.516912104161334),
                (1.3877621877018949, -4.516912104161338),
                (1.3513951321611892, -4.5169121041613325),
                (1.3151574693548629, -4.516912104161331),
                (1.2790696211428751, -4.516912104161331),
                (1.2431514990543158, -4.51691210416133),
            ],
        },
        EngineZone {
            text: "Q5-34AA7F2578AE1-A",
            id: 0x20ad_2a9f_c95e_2b84,
            points_at_0: 129,
            extent: [
                -FRAC_PI_2,
                2.366868196440341,
                -1.5707962991673,
                -1.9759144468930274,
            ],
            ring_at_1: &[
                (-FRAC_PI_2, -2.946115600078641),
                (-FRAC_PI_2, -2.946115954654217),
                (-FRAC_PI_2, -2.9461156178074233),
                (-FRAC_PI_2, -2.9461159369254344),
                (-FRAC_PI_2, -2.946115635536197),
                (-FRAC_PI_2, -2.9461159191966617),
                (-FRAC_PI_2, -2.94611565326498),
                (-FRAC_PI_2, -2.9461159014678793),
                (-FRAC_PI_2, -2.946115670993763),
                (-FRAC_PI_2, -2.9461158837390964),
                (-FRAC_PI_2, -2.9461156887225375),
                (-FRAC_PI_2, -2.9461158660103224),
                (-FRAC_PI_2, -2.94611570645132),
                (-FRAC_PI_2, -2.946115848281539),
                (-FRAC_PI_2, -2.9461157241800917),
                (-FRAC_PI_2, -2.9461158305527664),
                (-FRAC_PI_2, -2.946115741908875),
                (-FRAC_PI_2, -2.946115812823984),
                (-FRAC_PI_2, -2.946115759637658),
                (-FRAC_PI_2, -2.946115795095199),
                (-FRAC_PI_2, -2.9461157773664275),
                (-FRAC_PI_2, -2.9461157596376495),
                (-FRAC_PI_2, -2.9461157950952086),
                (-FRAC_PI_2, -2.946115741908866),
                (-FRAC_PI_2, -2.946115812823994),
                (-FRAC_PI_2, -2.9461157241800935),
                (-FRAC_PI_2, -2.9461158305527646),
                (-FRAC_PI_2, -2.94611570645131),
                (-FRAC_PI_2, -2.946115848281549),
                (-FRAC_PI_2, -2.9461156887225273),
                (-FRAC_PI_2, -2.9461158660103326),
                (-FRAC_PI_2, -2.9461156709937546),
                (-FRAC_PI_2, -2.946115883739104),
                (-FRAC_PI_2, -2.9461156532649713),
                (-FRAC_PI_2, -2.9461159014678895),
                (-FRAC_PI_2, -2.9461156355361884),
                (-FRAC_PI_2, -2.946115919196672),
                (-FRAC_PI_2, -2.9461156178074153),
                (-FRAC_PI_2, -2.9461159369254437),
                (-FRAC_PI_2, -2.9461156000786324),
                (-FRAC_PI_2, -2.946115954654228),
                (-1.5707963052689122, -2.0318979303555977),
                (-1.5707962991673001, -2.632420323686162),
                (-1.5707962994183549, -2.6620455575788724),
                (-1.5707962996449503, -2.6921921288418775),
                (-1.5707962998464697, -2.7228177885641918),
                (-1.570796300022347, -2.7538752979905397),
                (-1.5707963001720737, -2.785312754975098),
                (-1.5707963002952066, -2.817074040741738),
                (-1.5707963003913734, -2.8490992162667617),
                (-1.5707963004602792, -2.881325185317659),
                (-1.5707963005017094, -2.9136864519925),
                (-1.5707963005155336, -2.946115778710198),
                (-1.5707963005017092, -2.978545087735218),
                (-1.5707963004602794, -3.01090635841417),
                (-1.570796300391373, -3.043132318426572),
                (-1.5707963002952066, -3.075157496607383),
                (-1.5707963001720737, -3.106918779007823),
                (-1.5707963000223466, -3.1383562395881537),
                (-1.5707962998464697, -3.1694137468337487),
                (-1.5707962996449505, -3.200039410128618),
                (-1.5707962994183544, -3.230185973294495),
                (-1.5707962991673, -3.2598112140534012),
                (-1.570796305268912, -3.860333629893638),
            ],
        },
    ];
    const RTEA3H_ENGINE: [EngineZone; 15] = [
        EngineZone {
            text: "B4-4-A",
            id: 0x0280_0000_0000_0010,
            points_at_0: 120,
            extent: [
                -0.2500094092468502,
                -0.1845444184936117,
                0.250009409246847,
                0.19547687622336504,
            ],
            ring_at_1: &[
                (0.250009409246847, -0.0006428856160125025),
                (0.12536404336291887, -0.1845444184936117),
                (-0.12536404336292198, -0.18454441849361156),
                (-0.2500094092468502, -0.000642885616012303),
                (-0.12155967180226708, 0.195476876223365),
                (0.12155967180226401, 0.195476876223365),
            ],
        },
        EngineZone {
            text: "B4-0-A",
            id: 0x0280_0000_0000_0000,
            points_at_0: 100,
            extent: [
                -0.2036952386447587,
                -0.5715013235358124,
                0.2036952386447558,
                -0.17092461819013663,
            ],
            ring_at_1: &[
                (0.2036952386447558, -0.42496272333822005),
                (-1.0801225347490178e-15, -0.5715013235358124),
                (-0.2036952386447587, -0.4249627233382203),
                (-0.12536404336292203, -0.18454441849361092),
                (0.12536404336291893, -0.1845444184936112),
            ],
        },
        EngineZone {
            text: "C4-49-B",
            id: 0x0480_0000_0000_0125,
            points_at_0: 60,
            extent: [
                -0.8436747564743802,
                0.12640177966343075,
                -0.7698053966201996,
                0.26455197278329934,
            ],
            ring_at_1: &[
                (-0.7715350927245634, 0.162424311166478),
                (-0.8087806104949777, 0.12640177966343075),
                (-0.8436747564743802, 0.15948654241898874),
                (-0.8436747564743801, 0.23146721002774123),
                (-0.8087806104949776, 0.26455197278329934),
                (-0.7715350927245634, 0.2285294412802519),
            ],
        },
        EngineZone {
            text: "C0-0-B",
            id: 0x0400_0000_0000_0001,
            points_at_0: 50,
            extent: [
                0.984309630561419,
                -3.023547354800163,
                1.0610950228559666,
                -2.868684199932694,
            ],
            ring_at_1: &[
                (1.0610950228559666, -2.946115777366428),
                (1.0308478558390992, -3.023547354800163),
                (0.9848741281568848, -2.990585150872912),
                (0.9848741281568851, -2.901646403859944),
                (1.0308478558390992, -2.868684199932694),
            ],
        },
        EngineZone {
            text: "CA-0-B",
            id: 0x0540_0000_0000_0001,
            points_at_0: 50,
            extent: [
                0.9843096305614176,
                0.11804529878963071,
                1.0610950228559632,
                0.27290845365709915,
            ],
            ring_at_1: &[
                (1.0610950228559632, 0.1954768762233649),
                (1.030847855839096, 0.11804529878963071),
                (0.9848741281568818, 0.1510075027168809),
                (0.9848741281568819, 0.239946249729849),
                (1.0308478558390963, 0.27290845365709915),
            ],
        },
        EngineZone {
            text: "CB-0-A",
            id: 0x0560_0000_0000_0000,
            points_at_0: 80,
            extent: [
                -1.0819391811619674,
                -3.0684140768167585,
                -0.9491258447912547,
                -2.823817477916098,
            ],
            ring_at_1: &[
                (-0.994112017553241, -3.0684140768167585),
                (-1.0743719065374713, -3.0325090446002894),
                (-1.0743719065374715, -2.8597225101325665),
                (-0.9941120175532407, -2.823817477916098),
                (-0.9491258447912547, -2.9461157773664284),
            ],
        },
        EngineZone {
            text: "C0-6-A",
            id: 0x0400_0000_0000_0018,
            points_at_0: 92,
            extent: [
                1.3199184530455632,
                -0.22547598046439873,
                1.4498073274394094,
                0.616429732911133,
            ],
            ring_at_1: &[
                (1.3199184530455634, 0.36539099867805225),
                (1.3646571882851195, 0.616429732911133),
                (1.439150168364611, 0.5216243364831021),
                (1.4391501683646113, -0.13067058403636997),
                (1.3646571882851188, -0.22547598046439873),
                (1.3199184530455632, 0.025562753768679503),
            ],
        },
        EngineZone {
            text: "E0-7A-A",
            id: 0x0800_0000_0000_01e8,
            points_at_0: 200,
            extent: [
                1.54849813987436,
                -2.0293282597892306,
                1.5663367852773795,
                -0.72131064135362,
            ],
            ring_at_1: &[
                (1.559206101219953, -0.72131064135362),
                (1.5663367852773795, -1.3753194505711943),
                (1.5660771656275025, -1.4505469725708977),
                (1.5657963442326301, -1.5176325904887031),
                (1.5654974471184144, -1.577370827525319),
                (1.5651832216722836, -1.6305780850461387),
                (1.5648560181864193, -1.678035996017196),
                (1.5645178166656344, -1.7204613777667466),
                (1.564170271538959, -1.7584940275890815),
                (1.5638147601722723, -1.7926952899515738),
                (1.5634524285234233, -1.8235524714620897),
                (1.5630842312859063, -1.851485979506548),
                (1.562710965901824, -1.8768573432241638),
                (1.5623333007739924, -1.8999771104387166),
                (1.5619517983833977, -1.9211121246781822),
                (1.5615669341082197, -1.9404919825451612),
                (1.5611791114988487, -1.9583146341532331),
                (1.5607886746698219, -1.9747511719641153),
                (1.5603959183638914, -1.989949890108837),
                (1.5600010961435091, -2.004039707763288),
                (1.559604427077776, -2.0171330485097774),
                (1.5592061012199547, -2.0293282597892306),
                (1.551268251967315, -1.7462683257258176),
                (1.54849813987436, -1.375319450571464),
                (1.551268251967314, -1.0043705754171026),
            ],
        },
        EngineZone {
            text: "C0-4-C",
            id: 0x0400_0000_0000_0012,
            points_at_0: 251,
            extent: [
                1.486572969103528,
                -2.9461157773664284,
                FRAC_PI_2,
                0.19547687622336474,
            ],
            ring_at_1: &[
                (1.5115453676909434, -2.167969108190714),
                (1.486572969103529, -1.6305780850464127),
                (1.486572969103528, -1.1200608160966166),
                (1.5115453676909412, -0.5826697929523135),
                (1.5506444122892926, 0.19547687508745848),
                (1.5526596667662491, 0.19547687496124655),
                (1.5546749013396235, 0.19547687480348164),
                (1.556690118221338, 0.19547687460064123),
                (1.5587053196231684, 0.19547687433018734),
                (1.5607205077567634, 0.19547687395155178),
                (1.5627356848336624, 0.19547687338359832),
                (1.5647508530653138, 0.19547687243700929),
                (1.5667660146630935, 0.19547687054383145),
                (1.5687811718383222, 0.19547686486429794),
                (FRAC_PI_2, 0.19547687622336474),
                (FRAC_PI_2, -2.9461157773664284),
                (1.5687811718165614, -2.946115766007361),
                (1.566766014641333, -2.946115771686905),
                (1.5647508530435532, -2.9461157735800723),
                (1.5627356848119016, -2.9461157745266613),
                (1.5607205077350024, -2.946115775094619),
                (1.5587053196014073, -2.946115775473254),
                (1.5566901181995767, -2.9461157757437073),
                (1.554674901317862, -2.9461157759465446),
                (1.5526596667444874, -2.9461157761043095),
                (1.550644412267531, -2.9461157762305215),
            ],
        },
        EngineZone {
            text: "F0-79-A",
            id: 0x0a00_0000_0000_01e4,
            points_at_0: 65,
            extent: [
                1.565924341331054,
                1.7662732029505204,
                FRAC_PI_2,
                -1.3753194505123727,
            ],
            ring_at_1: &[
                (1.5693098141392021, -4.516912104162337),
                (1.566932939046034, -3.862903294943743),
                (1.565924341331054, -3.272263237626309),
                (1.5659243413310544, -2.6199683171065646),
                (1.566932939046034, -2.0293282597889934),
                (1.5693098141392028, -1.37531945057052),
                (1.5694584654300692, -1.3753194505704076),
                (1.5696071167129468, -1.375319450570267),
                (1.5697557679887235, -1.3753194505700865),
                (1.5699044192582872, -1.3753194505698456),
                (1.570053070522525, -1.3753194505695086),
                (1.5702017217823252, -1.3753194505690027),
                (1.5703503730385746, -1.3753194505681599),
                (1.570499024292161, -1.375319450566474),
                (1.5706476755439727, -1.3753194505614164),
                (FRAC_PI_2, -2.946115786254691),
                (FRAC_PI_2, -2.946115776080306),
                (FRAC_PI_2, -2.9461157684781627),
                (FRAC_PI_2, -2.946115752045054),
                (FRAC_PI_2, -2.9461157684781627),
                (FRAC_PI_2, -2.946115752045054),
                (1.5706476755439722, -4.516912104171439),
                (1.5704990242921604, -4.516912104166382),
                (1.5703503730385746, -4.516912104162264),
                (1.570201721782325, -4.51691210416203),
                (1.570053070522524, -4.516912104163348),
                (1.569904419258287, -4.51691210416301),
                (1.569755767988723, -4.5169121041627704),
                (1.5696071167129453, -4.516912104162589),
                (1.5694584654300685, -4.516912104162449),
            ],
        },
        EngineZone {
            text: "M0-40DF8-A",
            id: 0x1800_0000_0010_37e0,
            points_at_0: 62,
            extent: [
                1.570794099094735,
                1.7662731477143492,
                FRAC_PI_2,
                -1.3753193929210914,
            ],
            ring_at_1: &[
                (1.5707956470909605, -4.516912107124369),
                (1.5707945602721942, -3.8629032951402675),
                (1.570794099094735, -3.2722632376604968),
                (1.5707940990947353, -2.6199683168745342),
                (1.5707945602721942, -2.0293282590797017),
                (1.5707956470909605, -1.3753194482655138),
                (1.5707957150613543, -1.3753194479049997),
                (1.5707957830317478, -1.375319448040988),
                (1.5707958510021414, -1.3753194472772206),
                (1.570795918972535, -1.3753194467281689),
                (1.5707959869429287, -1.3753194461472182),
                (1.5707960549133222, -1.3753194445718349),
                (1.5707961228837157, -1.3753194431976765),
                (1.5707961908541095, -1.375319438572138),
                (1.5707962588245028, -1.3753194284499661),
                (FRAC_PI_2, -2.9461157773669706),
                (1.570796258824503, -4.5169121262828895),
                (1.5707961908541093, -4.516912107713219),
                (1.5707961228837155, -4.516912114038142),
                (1.5707960549133222, -4.516912111568937),
                (1.5707959869429282, -4.516912108585638),
                (1.570795918972535, -4.516912107848253),
                (1.5707958510021414, -4.516912107321549),
                (1.5707957830317483, -4.5169121088037425),
                (1.5707957150613538, -4.516912107453596),
            ],
        },
        EngineZone {
            text: "Q0-1486BA0-C",
            id: 0x2000_0000_0521_ae82,
            points_at_0: 248,
            extent: [
                1.5707963091915085,
                -1.6305780044326654,
                FRAC_PI_2,
                0.19547725547355377,
            ],
            ring_at_1: &[
                (FRAC_PI_2, 0.19547668659824016),
                (FRAC_PI_2, 0.19547706584842128),
                (1.5707963091915094, -1.6305780044326654),
                (1.5707963091915085, -1.1200607383241052),
                (FRAC_PI_2, 0.19547668659830866),
                (FRAC_PI_2, 0.19547706584849056),
                (FRAC_PI_2, 0.19547687573359873),
                (FRAC_PI_2, 0.19547687671313163),
                (FRAC_PI_2, 0.19547687573359876),
                (FRAC_PI_2, 0.19547687671313074),
                (FRAC_PI_2, 0.19547687573359876),
                (FRAC_PI_2, 0.19547687671313074),
                (FRAC_PI_2, 0.19547687573359873),
                (FRAC_PI_2, 0.19547687671313074),
                (FRAC_PI_2, 0.19547687573359876),
                (FRAC_PI_2, 0.19547687671313163),
                (FRAC_PI_2, 0.19547687573359873),
                (FRAC_PI_2, 0.19547687671313163),
                (FRAC_PI_2, 0.19547687573359876),
                (FRAC_PI_2, 0.19547687671313074),
                (FRAC_PI_2, 0.19547687573359873),
                (FRAC_PI_2, 0.1954768767131303),
                (FRAC_PI_2, 0.19547687573359876),
                (FRAC_PI_2, 0.1954768767131303),
                (FRAC_PI_2, 0.19547687573359873),
                (FRAC_PI_2, 0.19547687671313163),
                (FRAC_PI_2, 0.1954768757335988),
                (FRAC_PI_2, 0.19547687671313163),
                (FRAC_PI_2, 0.1954768757335988),
                (FRAC_PI_2, 0.1954768767131303),
                (FRAC_PI_2, 0.1954768757335988),
                (FRAC_PI_2, 0.1954768767131303),
                (FRAC_PI_2, 0.1954768757335988),
                (FRAC_PI_2, 0.1954768767131303),
                (FRAC_PI_2, 0.1954768757335988),
                (FRAC_PI_2, 0.19547687671313163),
                (FRAC_PI_2, 0.19547687573359876),
                (FRAC_PI_2, 0.1954768767131303),
                (FRAC_PI_2, 0.19547687573359873),
                (FRAC_PI_2, 0.19547687671313163),
                (FRAC_PI_2, 0.1954768757335988),
                (FRAC_PI_2, 0.1954768767131303),
                (FRAC_PI_2, 0.1954768757335987),
                (FRAC_PI_2, 0.1954768767131303),
                (FRAC_PI_2, 0.19547687573359873),
                (FRAC_PI_2, 0.1954768767131312),
                (FRAC_PI_2, 0.1954768757335988),
                (FRAC_PI_2, 0.1954768767131303),
            ],
        },
        EngineZone {
            text: "Q8-34AA7F2578AE0-C",
            id: 0x210d_2a9f_c95e_2b82,
            points_at_0: 47,
            extent: [
                1.5707963091915091,
                0.19547649697317482,
                FRAC_PI_2,
                2.0215317700863347,
            ],
            ring_at_1: &[
                (1.5707963091915091, 2.0215317700863347),
                (FRAC_PI_2, 0.19547668659832107),
                (FRAC_PI_2, 0.19547706584848576),
                (FRAC_PI_2, 0.19547687622336696),
                (FRAC_PI_2, 0.19547668659825268),
                (FRAC_PI_2, 0.19547706584841748),
                (1.5707963091982757, 1.5122598623613945),
            ],
        },
        EngineZone {
            text: "A0-0-A",
            id: 0x0000_0000_0000_0000,
            points_at_0: 139,
            extent: [
                0.36635950360263975,
                1.7662732030174286,
                FRAC_PI_2,
                -1.375319450570699,
            ],
            ring_at_1: &[
                (1.2074224434760321, -4.516912104161329),
                (0.6175926646873967, -3.731513940763877),
                (0.36635950360263975, -2.9461157773664284),
                (0.6175926646873967, -2.160717613968979),
                (1.207422443476032, -1.3753194505715276),
                (1.2441317451217546, -1.3753194505715272),
                (1.2807221560503335, -1.3753194505715265),
                (1.3172076991104953, -1.3753194505715256),
                (1.353602101704864, -1.3753194505715247),
                (1.389918838186159, -1.3753194505715232),
                (1.426171169973816, -1.3753194505715212),
                (1.4623721837848103, -1.3753194505715176),
                (1.4985348283394815, -1.3753194505715107),
                (1.5346719498769865, -1.3753194505714899),
                (FRAC_PI_2, -2.9461157773664284),
                (1.5346719498769867, -4.516912104161366),
                (1.4985348283394817, -4.516912104161346),
                (1.4623721837848107, -4.516912104161339),
                (1.4261711699738167, -4.516912104161335),
                (1.3899188381861591, -4.5169121041613325),
                (1.3536021017048636, -4.516912104161332),
                (1.3172076991104948, -4.516912104161331),
                (1.2807221560503335, -4.51691210416133),
                (1.2441317451217546, -4.516912104161329),
            ],
        },
        EngineZone {
            text: "Q5-34AA7F2578AE1-A",
            id: 0x20ad_2a9f_c95e_2b84,
            points_at_0: 131,
            extent: [
                -FRAC_PI_2,
                2.3657049958715146,
                -1.5707962992924243,
                -1.9747512286321363,
            ],
            ring_at_1: &[
                (-FRAC_PI_2, -2.946115587741347),
                (-FRAC_PI_2, -2.946115966991506),
                (-FRAC_PI_2, -2.946115606703859),
                (-FRAC_PI_2, -2.946115948028994),
                (-FRAC_PI_2, -2.9461156256663603),
                (-FRAC_PI_2, -2.946115929066493),
                (-FRAC_PI_2, -2.9461156446288723),
                (-FRAC_PI_2, -2.946115910103981),
                (-FRAC_PI_2, -2.9461156635913848),
                (-FRAC_PI_2, -2.9461158911414684),
                (-FRAC_PI_2, -2.946115682553886),
                (-FRAC_PI_2, -2.946115872178968),
                (-FRAC_PI_2, -2.946115701516398),
                (-FRAC_PI_2, -2.9461158532164546),
                (-FRAC_PI_2, -2.9461157204789),
                (-FRAC_PI_2, -2.946115834253953),
                (-FRAC_PI_2, -2.946115739441413),
                (-FRAC_PI_2, -2.9461158152914404),
                (-FRAC_PI_2, -2.9461157584039253),
                (-FRAC_PI_2, -2.9461157963289275),
                (-FRAC_PI_2, -2.9461157773664275),
                (-FRAC_PI_2, -2.9461157584039164),
                (-FRAC_PI_2, -2.946115796328937),
                (-FRAC_PI_2, -2.946115739441404),
                (-FRAC_PI_2, -2.9461158152914493),
                (-FRAC_PI_2, -2.9461157204789017),
                (-FRAC_PI_2, -2.9461158342539493),
                (-FRAC_PI_2, -2.9461157015163892),
                (-FRAC_PI_2, -2.9461158532164617),
                (-FRAC_PI_2, -2.946115682553877),
                (-FRAC_PI_2, -2.9461158721789764),
                (-FRAC_PI_2, -2.9461156635913754),
                (-FRAC_PI_2, -2.9461158911414778),
                (-FRAC_PI_2, -2.946115644628863),
                (-FRAC_PI_2, -2.9461159101039898),
                (-FRAC_PI_2, -2.946115625666351),
                (-FRAC_PI_2, -2.946115929066502),
                (-FRAC_PI_2, -2.9461156067038496),
                (-FRAC_PI_2, -2.946115948029004),
                (-FRAC_PI_2, -2.9461155877413376),
                (-FRAC_PI_2, -2.946115966991516),
                (-1.5707963049859737, -2.0293283204540282),
                (-1.5707962992924245, -2.6199683246622207),
                (-1.5707962996027671, -2.6499294230290773),
                (-1.5707962998930696, -2.6805595369950668),
                (-1.5707963001631735, -2.7118342708125796),
                (-1.570796300413066, -2.743725467807729),
                (-1.570796300642913, -2.7762016503912137),
                (-1.5707963008531, -2.809228652184506),
                (-1.5707963010442807, -2.8427702941042834),
                (-1.570796301217439, -2.87678947107924),
                (-1.5707963013739683, -2.9112495261067246),
                (-1.5707963015157695, -2.946115778763368),
                (-1.5707963015157693, -2.946115768710331),
                (-1.5707963013739679, -2.9809820042301376),
                (-1.5707963012174389, -3.0154420595672393),
                (-1.57079630104428, -3.0494612346148),
                (-1.5707963008530998, -3.0830028756528014),
                (-1.570796300642913, -3.116029880389401),
                (-1.5707963004130656, -3.1485060593955683),
                (-1.570796300163173, -3.180397256609059),
                (-1.5707962998930693, -3.211671994232417),
                (-1.5707962996027667, -3.2423021046477767),
                (-1.5707962992924243, -3.2722632041921798),
                (-1.5707963049859737, -3.86290322274219),
            ],
        },
    ];
}
