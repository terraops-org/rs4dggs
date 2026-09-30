//! A `Grid` is one DGGRS: a `Projection`, a `Topology` and an `Indexing`
//! composed into a working grid. Every operation that needs all three at once
//! lives here.
//!
//! Ported from py4dggs's `grid.py`, keeping its fallback structure: a topology
//! may expose an exact hierarchy or an ordered sub-zone enumeration, and where
//! it does not, the `Grid` either falls back to something grid-agnostic (the
//! congruent digit path) or refuses the request outright
//! (`Error::NoSubZoneOrder`), never approximating. Neighbours are the
//! topology's own: every topology here ports the engine's neighbour search.
//!
//! Source: py4dggs's `grid.py`.

use std::marker::PhantomData;

use crate::disk::Disk;
use crate::zone::Zone;
use crate::{Address, Error, GeoPoint, GridConfig, Indexing, Projection, Result, Topology, ZoneId};

/// Ceiling on how many sub-zones [`Grid::sub_zones`] will materialise into a
/// `Vec`, carried over from py4dggs's `MAX_MATERIALISED_SUB_ZONES`.
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
    /// keeps its resolution, its text and its congruent parent; it answers the
    /// engine's null geometry, the zero point as its centroid, six zero vertices
    /// and no neighbours, exactly as DGGAL answers it. The refusal of such a
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

    /// The zone's parents: one for a congruent grid, one or three for a
    /// non-congruent aperture-3 one.
    ///
    /// Where the topology exposes a geometric hierarchy it is used; otherwise
    /// the hierarchy is the congruent digit path, where the parent is the
    /// address with its last digit dropped.
    pub fn parents(&self, id: ZoneId) -> Vec<ZoneId> {
        let Ok(a) = self.address(id) else {
            return Vec::new();
        };
        if let Some(ps) = T::parents(&a) {
            return ps.iter().map(I::encode).collect();
        }
        I::parent(id).into_iter().collect()
    }

    /// The zone's children, in the topology's own order for a geometric
    /// hierarchy, or one per available digit for a congruent one.
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
    /// among others. On a congruent grid that is the single parent.
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

    /// Whether this zone is its parent's centroid child, which on a congruent
    /// grid every zone is, having exactly one parent.
    pub fn is_centroid_child(&self, id: ZoneId) -> bool {
        let Ok(a) = self.address(id) else {
            // Not a zone at all, so not any parent's centroid child.
            return false;
        };
        T::is_centroid_child(&a).unwrap_or(true)
    }

    // --- sub-zones: required topology override, no congruent default ---

    /// The guard the three sub-zone methods share, in py4dggs's own order:
    /// refuse a grid with no sub-zone order, then refuse a depth that would
    /// reach past the maximum resolution. Both cases would otherwise reach the
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
    /// order, this is the first of the corrected order, which need not be the engine's own
    /// `getFirstSubZone`: the correction can move the whole order, its start included.
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
    /// order, up to [`MAX_MATERIALISED_SUB_ZONES`] of them.
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
    /// A depth that would reach past `max_resolution` is refused the same way as on
    /// [`Grid::count_sub_zones`] and [`Grid::first_sub_zone`], with [`Error::InvalidZone`],
    /// where the engine's own `getSubZones` answers nothing there.
    pub fn sub_zones(&self, id: ZoneId, depth: u8) -> Result<Vec<ZoneId>> {
        let n = self.count_sub_zones(id, depth)?;
        if depth == 0 {
            return Ok(vec![id]); // See the note in `count_sub_zones`.
        }
        if n > MAX_MATERIALISED_SUB_ZONES {
            return Err(Error::TooManySubZones {
                count: n,
                limit: MAX_MATERIALISED_SUB_ZONES,
            });
        }
        let subs = T::sub_zones(&I::decode(id), depth).ok_or(Error::NoSubZoneOrder)?;
        Ok(subs.iter().map(I::encode).collect())
    }

    /// Where `sub` sits in this zone's sub-zone order, or `None` if it is no
    /// sub-zone of it. Generic, as the eC's own `dggrs.ec:115-131`: build the
    /// ordered sequence and find the position.
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
    /// coarser zone is none. So on a grid with no sub-zone order, an aperture-7 grid
    /// among them, `sub_zone_index(z, z)` answers `Ok(Some(0))` while
    /// `count_sub_zones(z, 0)` and the other sub-zone methods return
    /// [`Error::NoSubZoneOrder`]. Only a finer `sub` reaches `sub_zones` and the
    /// error. py4dggs's `grid.py` has exactly the same inconsistency, answering
    /// depth 0 above the topology here and checking for an order first in the
    /// others.
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
        let subs = self.sub_zones(id, sub_level - level)?;
        Ok(subs.iter().position(|&s| s == sub).map(|i| i as u64))
    }

    /// The sub-zone at `index`, `depth` levels below this zone. Generic, as
    /// the eC's own `dggrs.ec:133-149`.
    ///
    /// At depth 0 (so index 0, its only one) this answers the zone itself, where the engine's
    /// `getSubZoneAtIndex` inherits `getFirstSubZone`'s neighbour there; see
    /// [`Grid::first_sub_zone`] for the departure and its figures. For any other index, a topology
    /// that orders its sub-zones without materialising them (an aperture-3 hexagon, through the
    /// engine's own generators and their index skip) is asked directly, which reaches any depth to
    /// the maximum resolution cheaply; only a topology with no such shortcut falls back to
    /// indexing into [`Grid::sub_zones`], which [`MAX_MATERIALISED_SUB_ZONES`] still guards.
    ///
    /// At the zones on the edges of the rhombi where [`Grid::sub_zones`] corrects the engine's
    /// order, this follows the corrected order, through the same generators; elsewhere it is
    /// the engine's own index skip, bit for bit.
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
}

#[cfg(test)]
mod tests {
    use crate::indexings::I3h;
    use crate::projections::Isea;
    use crate::registry::{get_grid, igeo7};
    use crate::topologies::HexA3;
    use crate::{Error, GeoPoint, Grid, GridConfig, ZoneId};

    use super::MAX_MATERIALISED_SUB_ZONES;

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
    #[test]
    fn children_at_max_resolution_are_empty() {
        let g = igeo7();
        let z = g.zone_from_geo(38.7223, -9.1393, 19).unwrap();
        assert!(z.children().is_empty());
        assert_eq!(
            g.zone_from_geo(38.7223, -9.1393, 18)
                .unwrap()
                .children()
                .len(),
            7
        );
    }
    #[test]
    fn congruent_hierarchy() {
        let z = igeo7().zone_from_text("0064156").unwrap();
        assert_eq!(z.parent().unwrap().text_id(), "006415");
        assert_eq!(z.parents().len(), 1);
        assert!(z.is_centroid_child());
        assert_eq!(z.centroid_parent(), z.parent());
        assert!(z.is_immediate_child_of(&z.parent().unwrap()));
        assert!(z.parent().unwrap().parent().unwrap().is_ancestor_of(&z));
        assert!(igeo7().zone_from_text("00").unwrap().parent().is_none());
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
    #[test]
    fn aperture7_has_no_sub_zone_order() {
        let z = igeo7().zone_from_text("0064156").unwrap();
        assert!(matches!(z.sub_zones(1), Err(Error::NoSubZoneOrder)));
        assert!(matches!(z.count_sub_zones(1), Err(Error::NoSubZoneOrder)));
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
    /// and, bound as a raw identifier, it keeps its resolution, its text id and
    /// its congruent parent even though DGGAL gives it no geometry.
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
        assert_eq!(g.parents(id).len(), 1);
        let text = g.text_id(id);
        assert_eq!(text, format!("01{}", "1".repeat(20)));
        assert!(matches!(
            g.zone_from_text(&text),
            Err(Error::InvalidZone(_))
        ));
        assert!(matches!(
            g.count_sub_zones(id, 1),
            Err(Error::NoSubZoneOrder)
        ));
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
    /// index skip rather than by materialising 4.6e15 sub-zones, well within a second for the
    /// first, the last and a middle index.
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
                elapsed < std::time::Duration::from_secs(1),
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
    /// materialising, which the [`MAX_MATERIALISED_SUB_ZONES`] cap refuses at that depth, within a
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
                limit: MAX_MATERIALISED_SUB_ZONES,
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
}
