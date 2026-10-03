//! The zones of a level, enumerated: every zone of a level in the order of DGGAL's `listZones`,
//! and the zones of a level whose extent meets a bounding box, in that same order, each as an
//! iterator that can be entered after any zone of the level; with them the test of a zone's
//! extent against the box, the validation of the box and the search of the lattice for the
//! cells near it.
//!
//! Source: `GeoExtent::intersects` (`GeoExtent.ec:130-149`), as DGGAL v0.0.6 compiles it
//! (`0x42860`); the order of `listZones` (`RI3H.ec:824`, `RI7H.ec:829`), which the topologies
//! give as a lattice; the iterators, the entry after a zone, the rule for what a box may be and
//! the search of the lattice by its distance from the box are this crate's own API.

use crate::error::{Error, Result};
use crate::grid::Grid;
use crate::math;
use crate::registry::AnyGrid;
use crate::types::Extent;
use crate::{Indexing, Projection, Topology, ZoneId};
use std::f64::consts::{FRAC_PI_2, PI};
use std::iter::FusedIterator;

/// The finest level at which the zones of the aperture-7 grids are enumerated; see
/// [`Grid::max_box_level`], which says why it lies below the grids' last level.
pub(crate) const APERTURE_7_ENUMERATION_LIMIT: u8 = 14;

/// Every zone of one level of a grid, as a lazy iterator of identifiers, in the order in which
/// DGGAL's `listZones` gives the whole world.
///
/// It is taken from a [`Grid`] with [`Grid::zones`] or from an [`AnyGrid`] with
/// [`AnyGrid::zones`]; the type parameter is the handle it was taken from, and it yields
/// identifiers from either.
///
/// # Order
///
/// The zones of a level stand in a lattice: ten root rhombi, numbered 0 to 9, each a square of
/// rows and columns of cells, and two polar roots, 10 and 11, which are the two vertices of the
/// icosahedron that lie in no rhombus. A cell hosts one zone at an even level; at an odd level
/// it hosts three on the aperture-3 grids, and on the aperture-7 grids seven, or six where the
/// cell is a pentagon's. The order is that of the engine's own list: rhombus after rhombus;
/// within a rhombus row after row, from the top; within a row cell after cell, from the left;
/// within a cell the zones it hosts, in the order of the engine's sub-hexagons; and the zones of
/// the two polar roots last. It is the same on every call, in every run and on every platform.
///
/// On the aperture-3 grids that is the ascending order of the identifiers, so that sorting a
/// list of identifiers of one level gives it.
///
/// On the aperture-7 grids it is **not** the order of the Z7 identifiers, and no sorting of them
/// gives it. The engine sorts by the value of its internal identifier (the rhombus, the row, the
/// column, the sub-hexagon) and converts each zone to its Z7 identifier where it stands: level 0
/// comes as `01 06 02 07 03 08 04 09 05 10 00 11`.
///
/// # Entering after a zone
///
/// [`Zones::after`] gives the zones that follow a given zone of the level, which is how a
/// sequence is taken in pages: the last zone of one page is where the next begins. A page is
/// reached without a visit to anything that lies before it, so that its cost does not grow
/// with its distance from the first zone.
///
/// # Cost and memory
///
/// Nothing is computed until the iterator is consumed. One step asks the lattice for the zones
/// of one cell, wherever every cell hosts a zone, whatever the level and the position:
/// arithmetic on the aperture-3 grids, and up to seven conversions of an identifier on the
/// aperture-7 grids. A cell that hosted no zone would be passed over, the step going on to the
/// next, and a step would then ask for more than one. On the aperture-3 grids every cell hosts
/// one zone or three, by arithmetic. On the aperture-7 grids no cell that hosts nothing was
/// found at the levels enumerated: none in the engine's own list of the whole world to level
/// 5; none in a census of the lattice to level 9, where the zones number exactly the count of
/// the level (1,176,492, 8,235,432, 57,648,012 and 403,536,072 at levels 6 to 9) and no cell
/// hosts other than one zone, seven or, at a pentagon, six; and none in samples along the
/// borders of the rhombi at levels 11 to 14. The iterator holds its place and the zones of the
/// cell in hand, seven at most, and no list. A level has up to 5.6e16 zones on the aperture-3
/// grids and 6.8e12 at the finest level enumerated on the aperture-7 grids, so a caller takes
/// what it needs and no more.
///
/// Once the last zone has been yielded the iterator stays ended: `Zones` is a
/// [`FusedIterator`].
#[must_use = "the zones of a level are walked lazily, and nothing happens unless they are consumed"]
#[derive(Debug, Clone)]
pub struct Zones<S> {
    source: S,
    walk: Walk,
}

impl<S> Zones<S> {
    pub(crate) fn new(source: S, walk: Walk) -> Self {
        Zones { source, walk }
    }
}

impl<P: Projection, T: Topology, I: Indexing> Zones<&Grid<P, T, I>> {
    /// The zones of the level that come after `zone` in the order, `zone` itself excluded,
    /// wherever this iterator stood before: nothing already taken from it is counted.
    ///
    /// `zone` must be a zone of the level this iterator enumerates, and the sequence is then
    /// what remains of [`Grid::zones`] once `zone` has been yielded; after the last zone of the
    /// level it is empty. Entering costs no more the farther `zone` lies from the first zone:
    /// its place in the lattice is computed from the identifier, and nothing before it is
    /// visited.
    ///
    /// Since a zone to enter after commonly comes back from a client, it is refused, with
    /// [`crate::Error::InvalidZone`], unless it is an identifier that the enumeration of this
    /// level itself yields: a zone of another level, the null zone and an identifier that
    /// names no zone of the grid are all refused, where a walk begun at any of them would give
    /// a meaningless page in silence.
    ///
    /// # Examples
    ///
    /// ```
    /// let grid = rs4dggs::isea3h();
    /// let whole: Vec<_> = grid.zones(2)?.collect();
    /// // Two pages of forty zones, the second entered after the last zone of the first.
    /// let first: Vec<_> = grid.zones(2)?.take(40).collect();
    /// let second: Vec<_> = grid.zones(2)?.after(first[39])?.take(40).collect();
    /// assert_eq!([first, second].concat(), whole[..80]);
    /// // A zone of another level is refused.
    /// assert!(grid.zones(3)?.after(whole[0]).is_err());
    /// # Ok::<(), rs4dggs::Error>(())
    /// ```
    pub fn after(mut self, zone: ZoneId) -> Result<Self> {
        let place = self.source.lattice_place(self.walk.level, zone)?;
        self.walk.enter_after(place);
        Ok(self)
    }
}

impl Zones<AnyGrid> {
    /// The zones of the level that come after `zone` in the order, `zone` itself excluded; see
    /// [`Zones::after`] on a [`Grid`], which this forwards to, for what is refused and why.
    pub fn after(mut self, zone: ZoneId) -> Result<Self> {
        let place = self.source.lattice_place(self.walk.level, zone)?;
        self.walk.enter_after(place);
        Ok(self)
    }
}

impl<P: Projection, T: Topology, I: Indexing> Iterator for Zones<&Grid<P, T, I>> {
    type Item = ZoneId;

    fn next(&mut self) -> Option<ZoneId> {
        let (g, level) = (self.source, self.walk.level);
        self.walk
            .next_zone(|root, row, col| g.lattice_cell(level, root, row, col))
    }
}

impl Iterator for Zones<AnyGrid> {
    type Item = ZoneId;

    fn next(&mut self) -> Option<ZoneId> {
        let (g, level) = (self.source, self.walk.level);
        self.walk
            .next_zone(|root, row, col| g.lattice_cell(level, root, row, col))
    }
}

// Fused, because the walk's end is a sink: `next_zone` returns `None` only once the place has
// passed the last root, and from there it returns `None` at once and asks for nothing.
impl<P: Projection, T: Topology, I: Indexing> FusedIterator for Zones<&Grid<P, T, I>> {}

impl FusedIterator for Zones<AnyGrid> {}

/// The roots of the lattice: the ten rhombi, 0 to 9, and the two polar roots, 10 and 11.
const ROOTS: u8 = 12;

/// The walk of the lattice of one level, cell by cell, with the zones of a cell supplied at
/// each step, so that the two handles share one implementation and differ only in how they ask
/// for a cell.
///
/// It holds a place, the cell to be asked for next, and the zones of the cell asked for last
/// that are still to be yielded. A cell may host nothing, and is then passed over.
#[derive(Debug, Clone)]
pub(crate) struct Walk {
    level: u8,
    /// The number of cells along the side of a rhombus at `level`.
    edge: u64,
    /// The cell to be asked for next: its root, and under a rhombus its row and its column.
    /// A root of [`ROOTS`] is the end of the walk.
    root: u8,
    row: u64,
    col: u64,
    /// The zones of the cell asked for last, each with its key in the order, and the position
    /// of the next of them to be considered.
    pending: Vec<(ZoneId, u64)>,
    next_pending: usize,
    /// The key of the zone the walk was entered after: only zones of a greater key are yielded.
    floor: Option<u64>,
}

impl Walk {
    /// The walk of the level from its first zone, over a lattice of `edge` cells a side.
    pub(crate) fn new(level: u8, edge: u64) -> Self {
        Walk {
            level,
            edge,
            root: 0,
            row: 0,
            col: 0,
            pending: Vec::new(),
            next_pending: 0,
            floor: None,
        }
    }

    /// Moves the walk to the cell of `place`, a root, a row, a column and a key as
    /// [`Topology::locate`] gives them, so that the next zone is the first whose key exceeds
    /// that key: the zones of that cell that follow, then the cells after it.
    fn enter_after(&mut self, (root, row, col, key): (u8, u64, u64, u64)) {
        self.root = root;
        self.row = row;
        self.col = col;
        self.pending = Vec::new();
        self.next_pending = 0;
        self.floor = Some(key);
    }

    /// The next zone of the walk, `cell` giving the zones that the cell of a root, a row and a
    /// column hosts, each with its key, in ascending order of key.
    fn next_zone(
        &mut self,
        mut cell: impl FnMut(u8, u64, u64) -> Vec<(ZoneId, u64)>,
    ) -> Option<ZoneId> {
        loop {
            while let Some(&(zone, key)) = self.pending.get(self.next_pending) {
                self.next_pending += 1;
                if self.floor.is_none_or(|floor| key > floor) {
                    return Some(zone);
                }
            }
            if self.root >= ROOTS {
                self.pending = Vec::new();
                return None;
            }
            self.pending = cell(self.root, self.row, self.col);
            self.next_pending = 0;
            self.step();
        }
    }

    /// Moves the place to the cell after it: the next column; past the last column, the first
    /// column of the next row; past the last row, the first cell of the next root. A polar root
    /// is one cell. The comparisons are with `>=`, so that a side of nought, which no level
    /// enumerated has, still ends the walk after one cell a root.
    fn step(&mut self) {
        if self.root < 10 {
            self.col += 1;
            if self.col < self.edge {
                return;
            }
            self.col = 0;
            self.row += 1;
            if self.row < self.edge {
                return;
            }
            self.row = 0;
        }
        self.root += 1;
    }
}

/// The margin by which two extents must overlap to intersect, in radians. The eC writes
/// `10 * DBL_EPSILON` (`GeoExtent.ec:14`), which is 2.22e-15; the library holds 2e-15
/// (`0x3ce203af9ee75616`, `.rodata` at `0x4a4d8`, used at `0x428a5`), and a zone whose extent
/// overlaps a box by between the two is in the engine's answer, so the library's value is the
/// one that gives the engine's answers.
const RAD_EPSILON: f64 = f64::from_bits(0x3ce2_03af_9ee7_5616);

/// The value the cleared extent carries for its west end and, negated, for its east end and its
/// latitudes: `MAXDOUBLE` degrees, in radians (`0x7f91df46a2529d38`, at `0x4a4c8`). An extent
/// whose west end is not below it has no geometry, and is never cut at the antimeridian.
const CLEARED: f64 = f64::from_bits(0x7f91_df46_a252_9d38);

/// Whether the extent `a` intersects the extent `b`, both `[ll.lat, ll.lon, ur.lat, ur.lon]` in
/// radians, as `GeoExtent::intersects` answers it, in its order and strictly: two extents that
/// touch do not intersect, and each must overlap the other by more than [`RAD_EPSILON`] in both
/// directions.
///
/// An extent whose `ll.lon` exceeds its `ur.lon` runs eastwards over the antimeridian, and is
/// cut there at `+-Pi` into two, either of which meeting the other extent is a meeting; `a` is
/// cut before `b`. The engine makes the cut by calling itself on the pieces, which does not end
/// for a longitude beyond `Pi` in magnitude; this function compares the pieces, of `a` and of
/// `b` if both are cut, one with another, which ends, and for any longitude within `+-Pi` gives
/// the engine's answer.
///
/// faithful: the four comparisons are those of the eC (`GeoExtent.ec:145-148`), each after a
/// subtraction of the margin from the far end, in the order the library performs them
/// (`0x428ad` to `0x428ed`); the guard against the cleared extent is `ll.lon < MAXDOUBLE`
/// (`0x4286d`), the cut `0x42960` for `a` and `0x428f8` for `b`.
pub(crate) fn intersects(a: &[f64; 4], b: &[f64; 4]) -> bool {
    match (cut(a), cut(b)) {
        (None, None) => plain(a, b),
        (Some([west, east]), None) => plain(&west, b) || plain(&east, b),
        (None, Some([west, east])) => plain(a, &west) || plain(a, &east),
        (Some([a_west, a_east]), Some([b_west, b_east])) => {
            plain(&a_west, &b_west)
                || plain(&a_west, &b_east)
                || plain(&a_east, &b_west)
                || plain(&a_east, &b_east)
        }
    }
}

/// The two pieces, `[west, Pi]` and `[-Pi, east]`, of an extent over the antimeridian
/// (`GeoExtent.ec:132-135` and `:138-141`); `None` for any other, and for the cleared extent.
fn cut(e: &[f64; 4]) -> Option<[[f64; 4]; 2]> {
    (e[1] < CLEARED && e[1] > e[3]).then(|| [[e[0], e[1], e[2], PI], [e[0], -PI, e[2], e[3]]])
}

/// The last branch of `GeoExtent::intersects` (`GeoExtent.ec:145-148`): neither extent crosses
/// the antimeridian.
fn plain(a: &[f64; 4], b: &[f64; 4]) -> bool {
    b[2] - RAD_EPSILON > a[0]
        && a[2] - RAD_EPSILON > b[0]
        && b[3] - RAD_EPSILON > a[1]
        && a[3] - RAD_EPSILON > b[1]
}

/// The box `bbox`, in degrees as [`crate::Grid::extent`] gives an extent, as the radians that
/// [`intersects`] compares: each coordinate by the product with `Pi / 180`, made once, since
/// compared in degrees a tie can fall the other way.
///
/// What a box may be: any coordinate that is finite; latitudes between -90 and 90 with the
/// south not above the north; longitudes between -180 and 180, the west above the east
/// meaning a box that runs eastwards over the antimeridian. A box of no width or no height, or
/// of a single point, is a box; so is the whole world. Anything else is refused with the error
/// that names the coordinate at fault.
pub(crate) fn box_radians(bbox: &Extent) -> Result<[f64; 4]> {
    let (south, west, north, east) = (bbox.ll.lat, bbox.ll.lon, bbox.ur.lat, bbox.ur.lon);
    for (name, value) in [
        ("south", south),
        ("west", west),
        ("north", north),
        ("east", east),
    ] {
        if !value.is_finite() {
            return Err(Error::NonFinite { name });
        }
    }
    for (name, lat) in [("south", south), ("north", north)] {
        if !(-90.0..=90.0).contains(&lat) {
            return Err(Error::LatitudeOutOfRange { name });
        }
    }
    for (name, lon) in [("west", west), ("east", east)] {
        if !(-180.0..=180.0).contains(&lon) {
            return Err(Error::LongitudeOutOfRange { name });
        }
    }
    if south > north {
        return Err(Error::InvertedBox);
    }
    Ok([south, west, north, east].map(math::radians))
}

/// The zones of one level of a grid whose extent meets a bounding box, as a lazy iterator of
/// identifiers, in the order in which DGGAL's `listZones` gives the zones of a level.
///
/// It is taken from a [`Grid`] with [`Grid::zones_in_box`] or from an [`AnyGrid`] with
/// [`AnyGrid::zones_in_box`]; the type parameter is the handle it was taken from, and it yields
/// identifiers from either.
///
/// # Which zones
///
/// A zone of the level is yielded if and only if its extent meets the box, by the engine's own
/// test of two extents. The extent of a zone is what [`Grid::extent`] gives: the least and the
/// greatest latitude of its boundary and its westernmost and easternmost longitudes, a bounding
/// box of the zone. The test is strict: two extents meet where they overlap by more than 2e-15
/// radian both in latitude and in longitude, and two that only touch do not meet.
///
/// That is every zone that DGGAL's `listZones` returns for the box, and those that its own
/// search does not reach: on the aperture-3 grids the engine looks for zones about a sample of
/// points of the box, and leaves out a zone that belongs to another rhombus than the one its
/// samples fall in, the zone that holds the whole of a small box among them, and every zone of
/// a small cap about a pole; on the aperture-7 grids it descends from the twelve base cells
/// and, rarely, leaves out a zone whose extent reaches beyond the extents of its parents.
///
/// A caller must know three things of the rule.
///
/// - **The answer covers a little more than the box.** The extent is a bounding box, not the
///   hexagon: a zone whose extent meets the box at a corner is yielded although the zone itself
///   does not touch the box.
/// - **A box may run over the antimeridian.** A west above the east is no error: the box runs
///   eastwards from its west, through 180 degrees, to its east, so that a west of 10 and an
///   east of 5 make a box 355 degrees wide.
/// - **A box may have no area.** A segment yields the zones whose extents it strictly crosses,
///   and a point those whose extents strictly hold it: from one to three, or none where the
///   point lies on a line at which the extents of the zones on either side end. **A box of
///   no height exactly on a pole, a point on a pole among them, and a box of no width
///   exactly on the antimeridian (a west and an east both of 180 degrees, or both of -180)
///   yield no zone at all**, on every grid and at every level, since no extent reaches
///   beyond a pole and the pieces of an extent cut at the antimeridian end there: a segment
///   along the antimeridian crosses no extent strictly, whereas one a thousandth of a degree
///   from it yields zones. Points of the lines about which the grid is symmetric, on which
///   vertices and edges of zones lie, may yield none as well. On the equator, from a tenth to a quarter of the points tried yielded no zone (on
///   the aperture-7 grids at the even levels alone). On the four meridians along which the
///   rhombi are joined (see the cost, below), a third of them yielded none at the odd levels
///   of the aperture-3 grids, where the engine's extent of the zone that holds the point ends
///   up to 2.6e-10 radian before the meridian. [`Grid::zone_from_geo`] is the call for the
///   zone at a point.
///
/// # Order, and entering after a zone
///
/// The order is that of [`Zones`], of which this sequence is a part: rhombus after rhombus, row
/// after row, cell after cell, and the two polar zones last. On the aperture-3 grids that is
/// the ascending order of the identifiers; on the aperture-7 grids it is not the order of the
/// Z7 identifiers, and no sorting of them gives it.
///
/// [`ZonesInBox::after`] gives the zones of the answer that follow a given zone of the level,
/// which is how an answer is taken in pages: the last zone of one page is where the next
/// begins, and the pages together are the whole answer, each zone once.
///
/// # Cost and memory
///
/// Nothing is computed until the iterator is consumed, and no list is ever held: the iterator
/// keeps the box, its place in the lattice and the zones of the cell in hand, seven at most.
///
/// A step looks for the next cell of the lattice that may host a zone of the answer, by halving
/// blocks of cells: a block is left out when the point of the sphere at its middle lies farther
/// from the box than the block and a zone together can reach. The zones of a cell so found are
/// tested by their centroids, and those that pass by their extents, which is the dearer part.
/// Within a box that is one cell and one extent a zone. A step never costs by the place in the
/// answer, so that a page far into an answer costs what the first page costs, nor by the area
/// of the box; a box smaller than a zone is answered in some six hundred to three thousand
/// tests of a block, by the level, and an ordinary step takes a fraction of a millisecond.
///
/// **What a step costs by.** By the length, in cells, of the stretch of the rim of the box
/// that runs along a row of a rhombus. The search cannot leave out the cells along such a
/// stretch, which are near the box, and few of them host a zone of the answer: the step that
/// meets the stretch searches the whole of it. Most of that cost is the search itself, some
/// forty tests of a block for each cell of the stretch, and the lesser part the cells passed
/// over that host no zone of the answer: nine tenths against one to level 26 of the
/// aperture-3 grids, the cells and their extents weighing more at their finest level and on
/// the aperture-7 grids. Measured on an optimised build, a cell of the stretch costs from 7
/// to 47 microseconds: 7 to 9 on the aperture-3 grids to level 26, some 18 at level 32 and 35
/// at level 33; 11 to 21 on the aperture-7 grids at the even levels and 40 to 47 at the odd
/// ones.
///
/// - Where the rim crosses the rows of the lattice, the stretch is a few cells.
/// - Where the rim turns back, at the first or the last row that the box reaches in a rhombus,
///   the rows just beyond it follow the rim for a stretch that grows as the square root of the
///   side of the lattice. That is some hundreds of cells at the levels a map uses, and tens of
///   thousands at the three finest levels of the aperture-3 grids, where one such step was
///   measured at up to seven tenths of a second at the finest (0.10, 0.23 and 0.68 seconds
///   at levels 31, 32 and 33).
/// - **The worst case is a side of the box that lies on one of the meridians along which the
///   rhombi are joined**, to within some two cells: the stretch is then the whole of that
///   side. On the grids of the registry, whose orientation is DGGAL's own, there are four:
///   11.2 degrees east and 168.8 degrees west, beyond 58.4 degrees of latitude in either
///   hemisphere, and 78.8 degrees west and 101.2 degrees east, between 31.7 degrees south and
///   north. They are values of one decimal, which a box typed by hand, or a tiling by fifths
///   of a degree, lands on exactly. Measured, for a segment of one degree of the meridian of
///   11.2 degrees east: 7 milliseconds for the dearest step at level 20 of ISEA3H, 228 at
///   level 26, and 286 at level 14 of IGEO7; a third of a degree away from the meridian the
///   same step costs 0.3 millisecond. A grid built in another orientation has its rhombi
///   elsewhere, and the same holds of their sides wherever one lies along a meridian or a
///   parallel.
///
/// **The guard.** A side of `n` cells puts at least `n` zones in the answer, so that the
/// number of zones of an answer bounds every stretch of its rim. A caller that refuses a box
/// whose estimate ([`Grid::estimate_zones_in_box`]) exceeds `N` zones thereby bounds one step
/// by `N` cells, at the cost of a cell given above. A caller that sets no such limit can be
/// made to spend seconds on one step by a box with a side on one of those meridians at a deep
/// level.
///
/// Once the last zone has been yielded the iterator stays ended: `ZonesInBox` is a
/// [`FusedIterator`].
#[must_use = "the zones of a box are walked lazily, and nothing happens unless they are consumed"]
#[derive(Debug, Clone)]
pub struct ZonesInBox<S> {
    source: S,
    walk: Walk,
    search: Search,
}

impl<S> ZonesInBox<S> {
    pub(crate) fn new(source: S, walk: Walk, search: Search) -> Self {
        ZonesInBox {
            source,
            walk,
            search,
        }
    }
}

impl<P: Projection, T: Topology, I: Indexing> ZonesInBox<&Grid<P, T, I>> {
    /// The zones of the answer that come after `zone` in the order, `zone` itself excluded,
    /// wherever this iterator stood before: nothing already taken from it is counted.
    ///
    /// `zone` must be a zone of the level, and need not be a zone of the answer: the sequence
    /// is what [`Grid::zones_in_box`] yields beyond the place of `zone` in the order of the
    /// level, and is empty after the last zone of the answer. Entering costs no more the
    /// farther `zone` lies from the first zone: its place in the lattice is computed from the
    /// identifier, and nothing before it is visited.
    ///
    /// Since a zone to enter after commonly comes back from a client, it is refused, with
    /// [`crate::Error::InvalidZone`], unless it is an identifier that the enumeration of the
    /// level itself yields, exactly as [`Zones::after`] refuses it: a zone of another level,
    /// the null zone and an identifier that names no zone of the grid are all refused.
    ///
    /// # Examples
    ///
    /// ```
    /// use rs4dggs::{Extent, GeoPoint};
    /// let grid = rs4dggs::igeo7();
    /// let iberia = Extent {
    ///     ll: GeoPoint { lat: 36.0, lon: -10.0 },
    ///     ur: GeoPoint { lat: 44.0, lon: 3.0 },
    /// };
    /// let whole: Vec<_> = grid.zones_in_box(5, &iberia)?.collect();
    /// // Pages of ten zones, each entered after the last zone of the one before.
    /// let mut paged: Vec<_> = grid.zones_in_box(5, &iberia)?.take(10).collect();
    /// while let Some(&last) = paged.last() {
    ///     let page: Vec<_> = grid.zones_in_box(5, &iberia)?.after(last)?.take(10).collect();
    ///     if page.is_empty() {
    ///         break;
    ///     }
    ///     paged.extend(page);
    /// }
    /// assert_eq!(paged, whole);
    /// # Ok::<(), rs4dggs::Error>(())
    /// ```
    pub fn after(mut self, zone: ZoneId) -> Result<Self> {
        let place = self.source.lattice_place(self.walk.level, zone)?;
        self.walk.enter_after(place);
        Ok(self)
    }
}

impl ZonesInBox<AnyGrid> {
    /// The zones of the answer that come after `zone` in the order, `zone` itself excluded;
    /// see [`ZonesInBox::after`] on a [`Grid`], which this forwards to, for what is refused
    /// and why.
    pub fn after(mut self, zone: ZoneId) -> Result<Self> {
        let place = self.source.lattice_place(self.walk.level, zone)?;
        self.walk.enter_after(place);
        Ok(self)
    }
}

impl<P: Projection, T: Topology, I: Indexing> Iterator for ZonesInBox<&Grid<P, T, I>> {
    type Item = ZoneId;

    fn next(&mut self) -> Option<ZoneId> {
        let (g, level, search) = (self.source, self.walk.level, &self.search);
        self.walk.next_zone_in_box(
            search,
            |x, y| g.lattice_point(x, y),
            |root, row, col| g.lattice_cell(level, root, row, col),
            |zone| g.zone_meets_box(zone, &search.bbox, search.reach),
        )
    }
}

impl Iterator for ZonesInBox<AnyGrid> {
    type Item = ZoneId;

    fn next(&mut self) -> Option<ZoneId> {
        let (g, level, search) = (self.source, self.walk.level, &self.search);
        self.walk.next_zone_in_box(
            search,
            |x, y| g.lattice_point(x, y),
            |root, row, col| g.lattice_cell(level, root, row, col),
            |zone| g.zone_meets_box(zone, &search.bbox, search.reach),
        )
    }
}

// Fused for the reason [`Zones`] is: the end of the walk is a sink.
impl<P: Projection, T: Topology, I: Indexing> FusedIterator for ZonesInBox<&Grid<P, T, I>> {}

impl FusedIterator for ZonesInBox<AnyGrid> {}

impl Walk {
    /// The next zone of the walk that `meets` the box, the walk passing over every cell that
    /// `search` finds too far from the box to host such a zone: `point` gives the point of the
    /// sphere at a point of the plane of the rhombi, `cell` the zones a cell hosts, as in
    /// [`Walk::next_zone`]. The two polar roots are one cell each, which is always asked for.
    fn next_zone_in_box(
        &mut self,
        search: &Search,
        point: impl Fn(f64, f64) -> Option<(f64, f64)>,
        mut cell: impl FnMut(u8, u64, u64) -> Vec<(ZoneId, u64)>,
        mut meets: impl FnMut(ZoneId) -> bool,
    ) -> Option<ZoneId> {
        loop {
            while let Some(&(zone, key)) = self.pending.get(self.next_pending) {
                self.next_pending += 1;
                if self.floor.is_none_or(|floor| key > floor) && meets(zone) {
                    return Some(zone);
                }
            }
            if self.root >= ROOTS {
                self.pending = Vec::new();
                return None;
            }
            if self.root < 10 && !self.seek(search, &point) {
                continue;
            }
            self.pending = cell(self.root, self.row, self.col);
            self.next_pending = 0;
            self.step();
        }
    }

    /// Moves the place, a cell of a rhombus, to the first cell at or after it, in the order of
    /// the walk, that may host a zone of the answer: in the row in hand, or else in the first
    /// row after it that has such a cell. Where the rhombus has none left, the place becomes
    /// the first cell of the next root, and the answer is `false`.
    fn seek(&mut self, search: &Search, point: &impl Fn(f64, f64) -> Option<(f64, f64)>) -> bool {
        loop {
            if let Some(col) = search.next_cell(point, self.root, self.row, self.col) {
                self.col = col;
                return true;
            }
            self.col = 0;
            match search.next_row(point, self.root, self.row + 1) {
                Some(row) => self.row = row,
                None => {
                    self.row = 0;
                    self.root += 1;
                    return false;
                }
            }
        }
    }
}

/// Half the side of a triangle of the icosahedral net, on the unit sphere, where the net is
/// drawn so that a triangle has the area of a face, a twentieth of the sphere's: the side of an
/// equilateral triangle of area `Pi / 5` is `sqrt(4 Pi / (5 sqrt(3)))`.
const HALF_SIDE: f64 = 0.602_295_502_927_627_3;

/// The bound on how much the inverse projection stretches a length of the net: the ratio of a
/// length on the sphere to the length of its image in the net. The greatest measured is 1.1655
/// on the ISEA projection, 1.107 on IVEA and 1.131 on RTEA; this is a measured greatest with a
/// margin, and not a proved bound.
const STRETCH: f64 = 1.25;

/// The bound on how far a point of a zone's boundary lies from the zone's centroid, in zone
/// widths, the width of a zone being the square root of the sphere's area over the number of
/// zones of the level. The greatest measured is 0.72233, to level 30 on the aperture-3 grids
/// and to level 14 on the aperture-7 grids; a measured greatest with a margin, not a proved
/// bound.
const RADIUS_IN_ZONE_WIDTHS: f64 = 0.75;

/// What is added to that radius, in radians, for the three finest levels of the aperture-3
/// grids, where a zone is some 1.5e-8 radian wide and the rounding of its boundary's own
/// coordinates is no longer small beside it: the greatest excess measured there, over the
/// 0.7224 widths that hold the radius at every coarser level, is 1.5e-8 radian.
const RADIUS_FLOOR: f64 = 3e-8;

/// The bound on how far a zone's centroid lies from the top left corner of the cell that hosts
/// it, in sides of a cell. The greatest measured is 0.672 on the aperture-3 grids to level 32
/// and 0.915 at level 33, and 0.44 on the aperture-7 grids to level 14; a measured greatest
/// with a margin, not a proved bound. (On the aperture-7 grids from level 15 there are zones
/// beyond it, which is one reason for which the enumeration ends at level 14 there.)
const HOSTING_IN_SIDES: f64 = 1.0;

/// The length in the net, on the unit sphere, of a displacement `(dx, dy)` of the plane of the
/// ten rhombi, in which a rhombus is a unit square whose sides and whose diagonal from the top
/// left to the bottom right are sides of triangles (`ri5x6.ec:639-647`, `toIcosahedronNet`).
fn net_length(dx: f64, dy: f64) -> f64 {
    HALF_SIDE * math::sqrt((dx + dy) * (dx + dy) + 3.0 * (dx - dy) * (dx - dy))
}

/// Whether the point of latitude `lat` and longitude `lon` may lie within the angle `rho` of
/// the box, all in radians: true wherever it does, and at some points where it does not.
///
/// A point within `rho` of this one differs from it in latitude by `rho` at most, and, where
/// `phi`, the greater of the two latitudes in magnitude, lies short of the pole by more than
/// `rho`, in longitude by `asin(sin(rho) / cos(phi))` at most; nearer the pole than that it
/// may have any longitude.
pub(crate) fn near_box(lat: f64, lon: f64, rho: f64, bbox: &[f64; 4]) -> bool {
    let [south, west, north, east] = *bbox;
    if lat < south - rho || lat > north + rho {
        return false;
    }
    let phi = lat.abs() + rho;
    if rho + phi >= PI / 2.0 {
        return true;
    }
    let spread = math::asin((math::sin(rho) / math::cos(phi)).min(1.0));
    // The box's width, eastwards from its west, and how far east of that west the point lies,
    // each within one turn.
    let turn = |x: f64| x - math::TWO_PI * (x / math::TWO_PI).floor();
    let width = if east >= west {
        east - west
    } else {
        east - west + math::TWO_PI
    };
    let ahead = turn(lon - west);
    ahead <= width + spread || ahead >= math::TWO_PI - spread
}

/// The search of the lattice of one level for the cells that may host a zone whose extent
/// meets a box.
///
/// A cell hosts its zones about its top left corner, none farther from it than
/// [`HOSTING_IN_SIDES`]; a zone lies within `reach` of its centroid, and its extent with it;
/// and a length of the net is at most [`STRETCH`] times as long on the sphere. So a block of
/// cells can host a zone of the answer only if the point of the sphere at its middle lies
/// within the sum of those three of the box: its half-diagonal and a side, stretched, and the
/// reach.
#[derive(Debug, Clone)]
pub(crate) struct Search {
    /// The box, `[ll.lat, ll.lon, ur.lat, ur.lon]` in radians.
    bbox: [f64; 4],
    /// The number of cells along the side of a rhombus at the level.
    edge: u64,
    /// How far from a zone's centroid a point of its boundary may lie, in radians.
    reach: f64,
}

impl Search {
    /// The search for the box `bbox`, in radians, at a level of `zones` zones whose lattice has
    /// `edge` cells a side.
    pub(crate) fn new(bbox: [f64; 4], edge: u64, zones: u64) -> Self {
        let width = math::sqrt(4.0 * PI / zones as f64);
        Search {
            bbox,
            edge,
            reach: RADIUS_IN_ZONE_WIDTHS * width + RADIUS_FLOOR,
        }
    }

    /// Whether a cell of the block of rows `r0..r1` and columns `c0..c1` of the rhombus `root`
    /// may host a zone of the answer. A block is measured from its middle, a single cell from
    /// its top left corner, where its zones sit. Where the inverse projection gives no point
    /// the block is kept.
    fn may_hold(
        &self,
        point: &impl Fn(f64, f64) -> Option<(f64, f64)>,
        root: u8,
        (r0, r1): (u64, u64),
        (c0, c1): (u64, u64),
    ) -> bool {
        let p = self.edge as f64;
        let (hx, hy) = if r1 - r0 == 1 && c1 - c0 == 1 {
            (0.0, 0.0)
        } else {
            ((c1 - c0) as f64 / 2.0 / p, (r1 - r0) as f64 / 2.0 / p)
        };
        // The rhombus `root` is the unit square whose top left corner is at `(root / 2,
        // (root + 1) / 2)`, the divisions whole.
        let x = f64::from(root / 2) + c0 as f64 / p + hx;
        let y = f64::from(root.div_ceil(2)) + r0 as f64 / p + hy;
        let Some((lat, lon)) = point(x, y) else {
            return true;
        };
        let half_diagonal = net_length(hx, hy).max(net_length(hx, -hy));
        let rho =
            STRETCH * (half_diagonal + HOSTING_IN_SIDES * net_length(1.0 / p, 0.0)) + self.reach;
        near_box(lat, lon, rho, &self.bbox)
    }

    /// The first column from `from` on, in the row `row`, whose cell may host a zone of the
    /// answer.
    fn next_cell(
        &self,
        point: &impl Fn(f64, f64) -> Option<(f64, f64)>,
        root: u8,
        row: u64,
        from: u64,
    ) -> Option<u64> {
        if from >= self.edge {
            return None;
        }
        if self.may_hold(point, root, (row, row + 1), (from, from + 1)) {
            return Some(from);
        }
        self.first_col(point, root, row, (0, self.edge), from)
    }

    /// The first column of `lo..hi` that is not before `from`, in the row `row`, whose cell
    /// may host a zone of the answer, by halving the columns.
    fn first_col(
        &self,
        point: &impl Fn(f64, f64) -> Option<(f64, f64)>,
        root: u8,
        row: u64,
        (lo, hi): (u64, u64),
        from: u64,
    ) -> Option<u64> {
        if hi <= from || !self.may_hold(point, root, (row, row + 1), (lo, hi)) {
            return None;
        }
        if hi - lo == 1 {
            return Some(lo);
        }
        let mid = lo + (hi - lo) / 2;
        self.first_col(point, root, row, (lo, mid), from)
            .or_else(|| self.first_col(point, root, row, (mid, hi), from))
    }

    /// The first row from `from` on that has a cell that may host a zone of the answer.
    ///
    /// The row `from` itself is tried first, since within a box it is the one. Failing it,
    /// the rows after it are halved, the upper half being kept wherever it has such a cell:
    /// the first row is then found in as many halvings as the side has binary digits at most,
    /// each of them a search for any cell of a band of rows, which ends at the first it finds.
    /// (A single search of the whole rhombus for its first row, upper and left halves first,
    /// costs by the length of the box's rim wherever the rim rises to the right: it finds a row
    /// at the left, and then a better one in every block along the rim.)
    fn next_row(
        &self,
        point: &impl Fn(f64, f64) -> Option<(f64, f64)>,
        root: u8,
        from: u64,
    ) -> Option<u64> {
        if from >= self.edge {
            return None;
        }
        let columns = (0, self.edge);
        if self.first_col(point, root, from, columns, 0).is_some() {
            return Some(from);
        }
        // The first row lies in `lo..=hi`, and the row `hi` has such a cell.
        let mut lo = from + 1;
        let mut hi = self.any_cell(point, root, (lo, self.edge), columns)?;
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            match self.any_cell(point, root, (lo, mid + 1), columns) {
                Some(row) => hi = row,
                None => lo = mid + 1,
            }
        }
        Some(lo)
    }

    /// The row of a cell of the block that may host a zone of the answer, if the block has
    /// one: the first found by halving the block along its longer side, the upper or the left
    /// half first. It is a row of the block that has such a cell, and not always the first.
    fn any_cell(
        &self,
        point: &impl Fn(f64, f64) -> Option<(f64, f64)>,
        root: u8,
        (r0, r1): (u64, u64),
        (c0, c1): (u64, u64),
    ) -> Option<u64> {
        if r1 <= r0 || !self.may_hold(point, root, (r0, r1), (c0, c1)) {
            return None;
        }
        if r1 - r0 == 1 && c1 - c0 == 1 {
            return Some(r0);
        }
        if r1 - r0 >= c1 - c0 {
            let mid = r0 + (r1 - r0) / 2;
            self.any_cell(point, root, (r0, mid), (c0, c1))
                .or_else(|| self.any_cell(point, root, (mid, r1), (c0, c1)))
        } else {
            let mid = c0 + (c1 - c0) / 2;
            self.any_cell(point, root, (r0, r1), (c0, mid))
                .or_else(|| self.any_cell(point, root, (r0, r1), (mid, c1)))
        }
    }
}

// --- an upper bound on the zones of a box, by arithmetic ---

/// By how much a box is grown on every side before its area is counted in zones, in widths of
/// a zone, the width being the square root of the sphere's area over the count of zones.
///
/// A zone is counted when its extent meets the box. A corner of an extent lies at most
/// `sqrt(2)` radii of the zone from its centroid, and the zone within one radius of it, so
/// that a counted zone lies wholly within `sqrt(2) + 1` radii of the box; the radius of a
/// zone is at most 0.7224 widths (see [`RADIUS_IN_ZONE_WIDTHS`]), which makes 1.744 widths.
/// Where the zones are disjoint and of one area, the area of the box grown by that much, over
/// the area of a zone, is therefore not under their number.
const ESTIMATE_MARGIN_IN_ZONE_WIDTHS: f64 = 1.75;

/// What the estimate adds to the quotient of the areas, for three reasons. The twelve
/// pentagons have five sixths of a hexagon's area. At the coarsest levels the zones are too
/// few for the quotient to stand alone. And upon a pole the quotient alone is under the count:
/// the box is grown in longitude by the margin over the cosine of the latitude, which beside
/// a pole is less than the zones reach by arc, and at the finest level of the aperture-3 grids
/// extents reach the pole from centroids up to 1.7 widths away. A narrow wedge upon the north
/// pole there holds 8 zones where the quotient, rounded up, is 6: the estimate is 9, and
/// holds by this allowance alone, with one zone to spare, which is the thinnest margin
/// measured. With it the estimate is never nought.
const ESTIMATE_ALLOWANCE: u64 = 3;

/// The room left for the rounding of the arithmetic, as a factor of the quotient of the
/// areas: where a box is the whole sphere there is no rim to absorb it, and at the finest
/// levels of aperture 3 the count of zones is itself beyond what a `f64` holds exactly.
const ESTIMATE_ROUNDING: f64 = 1.0 + 16.0 * f64::EPSILON;

/// `sin(b) - sin(a)`, taken as `2 cos((a + b) / 2) sin((b - a) / 2)`. Beside a pole at a fine
/// level the two sines are the same number, or differ by a unit in the last place, where the
/// band between the two latitudes still holds zones; the product loses nothing there.
fn sine_difference(a: f64, b: f64) -> f64 {
    2.0 * math::cos((a + b) / 2.0) * math::sin((b - a) / 2.0)
}

/// The area, on the unit sphere, of the box `bbox` grown by `margin` radians of arc on every
/// side: `[south, west, north, east]` in radians, the west above the east where the box runs
/// over the antimeridian.
///
/// The grown box reaches from `south - margin` to `north + margin`, cut at the poles, and at
/// the latitude `lat` it is `width + 2 margin / cos(lat)` wide, or the whole circle where that
/// exceeds it, which is from the latitude `turn` towards either pole. Its area is the integral
/// of that width times `cos(lat)`, in closed form over the three bands the two latitudes
/// `-turn` and `turn` make between the poles, which is where the cut at the poles is made.
fn grown_area(bbox: &[f64; 4], margin: f64) -> f64 {
    let (south, north) = (bbox[0] - margin, bbox[2] + margin);
    let mut width = bbox[3] - bbox[1];
    if width < 0.0 {
        width += 2.0 * PI;
    }
    let rest = 2.0 * PI - width;
    let turn = if 2.0 * margin < rest {
        math::acos(2.0 * margin / rest)
    } else {
        0.0
    };
    let mut area = 0.0;
    for (from, to, whole) in [
        (-FRAC_PI_2, -turn, true),
        (-turn, turn, false),
        (turn, FRAC_PI_2, true),
    ] {
        let (a, b) = (south.max(from), north.min(to));
        if b > a {
            area += if whole {
                2.0 * PI * sine_difference(a, b)
            } else {
                width * sine_difference(a, b) + 2.0 * margin * (b - a)
            };
        }
    }
    area
}

/// An upper bound on the number of zones, of a level of `zones` zones, whose extent meets the
/// box `bbox`: the area of the box grown by [`ESTIMATE_MARGIN_IN_ZONE_WIDTHS`] widths of a
/// zone on every side, over the area of a zone, rounded up, plus [`ESTIMATE_ALLOWANCE`], and
/// never above `zones` plus that allowance.
///
/// `bbox` is `[south, west, north, east]` in radians, as [`box_radians`] gives it, with its
/// two latitudes already taken on the sphere on which the zones are equal in area.
pub(crate) fn estimate(bbox: &[f64; 4], zones: u64) -> u64 {
    let width = math::sqrt(4.0 * PI / zones as f64);
    let share = grown_area(bbox, ESTIMATE_MARGIN_IN_ZONE_WIDTHS * width) / (4.0 * PI);
    share_in_zones(share, zones)
        .saturating_add(ESTIMATE_ALLOWANCE)
        .min(zones + ESTIMATE_ALLOWANCE)
}

/// The share `share` of the sphere, in zones of a level of `zones` zones, rounded up, with
/// the room of [`ESTIMATE_ROUNDING`]. The cast saturates, and a share is not negative.
fn share_in_zones(share: f64, zones: u64) -> u64 {
    (zones as f64 * share * ESTIMATE_ROUNDING).ceil() as u64
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::interfaces::sealed::Private;
    use crate::types::GeoPoint;

    fn bits(b: [u64; 4]) -> [f64; 4] {
        b.map(f64::from_bits)
    }

    fn degrees(s: f64, w: f64, n: f64, e: f64) -> Extent {
        Extent {
            ll: GeoPoint { lat: s, lon: w },
            ur: GeoPoint { lat: n, lon: e },
        }
    }

    /// The margin is the library's 2e-15, which is what `.rodata` holds at `0x4a4d8` (the
    /// bytes `16 56 e7 9e af 03 e2 3c`), and not the 2.22e-15 of `10 * DBL_EPSILON`.
    #[test]
    fn the_margin_is_the_librarys() {
        assert_eq!(RAD_EPSILON.to_bits(), 0x3ce2_03af_9ee7_5616);
        assert_eq!(
            u64::from_le_bytes([0x16, 0x56, 0xe7, 0x9e, 0xaf, 0x03, 0xe2, 0x3c]),
            RAD_EPSILON.to_bits()
        );
        assert_eq!(RAD_EPSILON, 2e-15);
        assert_eq!(CLEARED.to_bits(), 0x7f91_df46_a252_9d38);
    }

    /// A zone of level 3 on ISEA7H_Z7 (`10625`) and a box whose south edge lies 37 units in the
    /// last place below the zone's north edge: 2.054e-15 of overlap, which is enough for the
    /// engine's 2e-15 and not for `10 * DBL_EPSILON`.
    const ZONE_10625: [u64; 4] = [
        0x3fcb_05fa_009b_54ae,
        0xc008_9e26_2d94_0998,
        0x3fd2_8226_8af3_5a68,
        0xc008_0d2a_918b_2e42,
    ];
    const BOX_AT_10625: [u64; 4] = [
        0x3fd2_8226_8af3_5a43,
        0xc008_6e21_7aba_75f3,
        0x3fd9_4cba_74ff_0784,
        0xc008_2b5a_1e3d_22ac,
    ];

    #[test]
    fn an_overlap_of_37_units_in_the_last_place_is_enough() {
        let (zone, bx) = (bits(ZONE_10625), bits(BOX_AT_10625));
        assert_eq!(ZONE_10625[2] - BOX_AT_10625[0], 37);
        assert!(intersects(&zone, &bx));
        assert!(intersects(&bx, &zone));
        // With 10 * DBL_EPSILON as the margin the zone would be out.
        let lower = 10.0 * f64::EPSILON;
        assert!(zone[2] - lower <= bx[0]);
        assert!(zone[2] - RAD_EPSILON > bx[0]);
    }

    #[test]
    fn extents_that_touch_do_not_intersect() {
        let a = [0.0, 0.0, 1.0, 1.0];
        // Sharing an edge, in each of the four directions.
        assert!(!intersects(&a, &[1.0, 0.0, 2.0, 1.0]));
        assert!(!intersects(&a, &[-1.0, 0.0, 0.0, 1.0]));
        assert!(!intersects(&a, &[0.0, 1.0, 1.0, 2.0]));
        assert!(!intersects(&a, &[0.0, -1.0, 1.0, 0.0]));
        // Overlapping by less than the margin, and by more.
        assert!(!intersects(&a, &[1.0 - 1e-15, 0.0, 2.0, 1.0]));
        assert!(intersects(&a, &[1.0 - 1e-14, 0.0, 2.0, 1.0]));
        // A point strictly inside, and the extent against itself.
        assert!(intersects(&a, &[0.5, 0.5, 0.5, 0.5]));
        assert!(intersects(&a, &a));
        // A point does not intersect what it lies on the edge of, nor a point another.
        assert!(!intersects(&a, &[1.0, 0.5, 1.0, 0.5]));
        assert!(!intersects(&[0.5, 0.5, 0.5, 0.5], &[0.5, 0.5, 0.5, 0.5]));
    }

    /// Each of the four comparisons is strict: an overlap of the margin exactly, in one
    /// direction alone, is no intersection, and an overlap of one and a half times the margin
    /// is one. The four are, in the order of the source: the second's north against the first's
    /// south, the first's north against the second's south, the second's east against the
    /// first's west, the first's east against the second's west.
    #[test]
    fn an_overlap_of_the_margin_exactly_is_no_intersection() {
        let a = [0.0, 0.0, 1.0, 1.0];
        let more = 1.5 * RAD_EPSILON;
        assert!(!intersects(&a, &[-1.0, 0.0, RAD_EPSILON, 1.0]));
        assert!(intersects(&a, &[-1.0, 0.0, more, 1.0]));
        assert!(!intersects(&[-1.0, 0.0, RAD_EPSILON, 1.0], &a));
        assert!(intersects(&[-1.0, 0.0, more, 1.0], &a));
        assert!(!intersects(&a, &[0.0, -1.0, 1.0, RAD_EPSILON]));
        assert!(intersects(&a, &[0.0, -1.0, 1.0, more]));
        assert!(!intersects(&[0.0, -1.0, 1.0, RAD_EPSILON], &a));
        assert!(intersects(&[0.0, -1.0, 1.0, more], &a));
    }

    #[test]
    fn the_cleared_extent_intersects_nothing() {
        let cleared = [CLEARED, CLEARED, -CLEARED, -CLEARED];
        let world = [-FRAC_PI_2, -PI, FRAC_PI_2, PI];
        assert!(!intersects(&cleared, &world));
        assert!(!intersects(&world, &cleared));
        assert!(!intersects(&cleared, &cleared));
    }

    /// An extent whose west end is the cleared one's is never cut at the antimeridian, although
    /// that end lies above its east end: cut, its piece from `-Pi` to its east end would meet
    /// the other extent.
    #[test]
    fn a_cleared_west_end_is_not_cut_at_the_antimeridian() {
        let a = [0.0, 0.0, 1.0, 1.0];
        let west_cleared = [0.0, CLEARED, 1.0, 0.5];
        assert!(!intersects(&west_cleared, &a));
        assert!(!intersects(&a, &west_cleared));
        // The same extent with a west end just above its east end is cut, and meets it.
        assert!(intersects(&[0.0, 3.0, 1.0, 0.5], &a));
        assert!(intersects(&a, &[0.0, 3.0, 1.0, 0.5]));
    }

    /// Each of the four cases of the cut at the antimeridian: neither extent over it, the first,
    /// the second, both. The pieces are `[west, Pi]` and `[-Pi, east]`.
    #[test]
    fn the_cut_at_the_antimeridian() {
        let d = PI / 180.0;
        // Over the antimeridian from 170 to -170 degrees, between 0 and 10 degrees of latitude.
        let over = [0.0, 170.0 * d, 10.0 * d, -170.0 * d];
        let east_side = [2.0 * d, 175.0 * d, 4.0 * d, 179.0 * d];
        let west_side = [2.0 * d, -179.0 * d, 4.0 * d, -175.0 * d];
        let outside = [2.0 * d, 0.0, 4.0 * d, 10.0 * d];
        let both_pieces = [2.0 * d, 179.0 * d, 4.0 * d, -179.0 * d];
        let another_over = [5.0 * d, 160.0 * d, 15.0 * d, -160.0 * d];
        let another_apart = [20.0 * d, 160.0 * d, 30.0 * d, -160.0 * d];

        // Neither over the antimeridian.
        assert!(intersects(
            &east_side,
            &[3.0 * d, 176.0 * d, 5.0 * d, 178.0 * d]
        ));
        // The extent over it, the box not.
        assert!(intersects(&over, &east_side));
        assert!(intersects(&over, &west_side));
        assert!(!intersects(&over, &outside));
        // The box over it, the extent not.
        assert!(intersects(&east_side, &over));
        assert!(intersects(&west_side, &over));
        assert!(!intersects(&outside, &over));
        // Both over it.
        assert!(intersects(&over, &another_over));
        assert!(!intersects(&over, &another_apart));
        assert!(intersects(&over, &both_pieces));
        assert!(intersects(&both_pieces, &over));
        // Touching at the antimeridian from the two sides: the pieces share only the line.
        let touching_west = [2.0 * d, 100.0 * d, 4.0 * d, 170.0 * d];
        let over_east = [0.0, 170.0 * d, 10.0 * d, -100.0 * d];
        assert!(!intersects(&touching_west, &over_east));
    }

    #[test]
    fn a_box_is_made_into_radians_by_the_product() {
        let r = box_radians(&degrees(38.0, -10.0, 39.0, -9.0)).unwrap();
        let k = PI / 180.0;
        assert_eq!(
            r.map(f64::to_bits),
            [38.0 * k, -10.0 * k, 39.0 * k, -9.0 * k].map(f64::to_bits)
        );
    }

    #[test]
    fn what_a_box_may_be() {
        // The whole world is the box of the poles and the antimeridian, to the bit.
        let world = box_radians(&degrees(-90.0, -180.0, 90.0, 180.0)).unwrap();
        assert_eq!(
            world.map(f64::to_bits),
            [-FRAC_PI_2, -PI, FRAC_PI_2, PI].map(f64::to_bits)
        );
        // West above east is a box over the antimeridian, 355 degrees wide.
        assert!(box_radians(&degrees(30.0, 10.0, 40.0, 5.0)).is_ok());
        // No width, no height, one point, one point on a pole.
        assert!(box_radians(&degrees(30.0, 10.0, 40.0, 10.0)).is_ok());
        assert!(box_radians(&degrees(30.0, 10.0, 30.0, 20.0)).is_ok());
        assert!(box_radians(&degrees(38.7223, -9.1393, 38.7223, -9.1393)).is_ok());
        assert!(box_radians(&degrees(90.0, 0.0, 90.0, 0.0)).is_ok());
    }

    #[test]
    fn a_box_with_a_coordinate_that_is_not_finite_is_refused() {
        for (i, name) in ["south", "west", "north", "east"].into_iter().enumerate() {
            for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
                let mut c = [10.0, 20.0, 30.0, 40.0];
                c[i] = bad;
                assert_eq!(
                    box_radians(&degrees(c[0], c[1], c[2], c[3])),
                    Err(Error::NonFinite { name })
                );
            }
        }
    }

    #[test]
    fn a_latitude_beyond_a_pole_is_refused() {
        assert_eq!(
            box_radians(&degrees(-90.000001, 0.0, 10.0, 10.0)),
            Err(Error::LatitudeOutOfRange { name: "south" })
        );
        assert_eq!(
            box_radians(&degrees(80.0, 0.0, 100.0, 10.0)),
            Err(Error::LatitudeOutOfRange { name: "north" })
        );
    }

    #[test]
    fn a_south_above_the_north_is_refused() {
        assert_eq!(
            box_radians(&degrees(40.0, 10.0, 30.0, 20.0)),
            Err(Error::InvertedBox)
        );
    }

    #[test]
    fn a_longitude_beyond_the_antimeridian_is_refused() {
        // With the west below the east, and with it above: the engine answers the first and
        // overflows its stack on the second.
        for (w, e, bad) in [
            (170.0, 181.0, "east"),
            (-190.0, -170.0, "west"),
            (170.0, -181.0, "east"),
            (181.0, 170.0, "west"),
        ] {
            assert_eq!(
                box_radians(&degrees(10.0, w, 20.0, e)),
                Err(Error::LongitudeOutOfRange { name: bad })
            );
        }
    }

    /// The literals above are the engine's own extent of the zone, and its answer for the box.
    #[cfg(feature = "oracle")]
    #[test]
    fn the_literals_are_the_engines() {
        let zone = dggal_oracle::zone_from_text(dggal_oracle::IGEO7, "10625");
        assert_eq!(
            dggal_oracle::extent(dggal_oracle::IGEO7, zone).map(f64::to_bits),
            ZONE_10625
        );
        let listed = dggal_oracle::list_zones(dggal_oracle::IGEO7, 3, Some(bits(BOX_AT_10625)));
        assert!(listed.contains(&zone));
    }

    /// A lattice of the test's own, of two cells a side: a cell hosts two zones, or none where
    /// its row and its column are of unlike parity, and a polar root one. The key of a zone is
    /// its identifier, which ascends along the walk. Every cell asked for is recorded.
    fn toy(asked: &mut Vec<(u8, u64, u64)>, root: u8, row: u64, col: u64) -> Vec<(ZoneId, u64)> {
        asked.push((root, row, col));
        let cell = u64::from(root) * 100 + row * 20 + col * 10;
        let hosted = if root >= 10 {
            vec![cell]
        } else if (row + col) % 2 == 0 {
            vec![cell, cell + 1]
        } else {
            Vec::new()
        };
        hosted.into_iter().map(|id| (ZoneId(id), id)).collect()
    }

    /// The walk asks for every cell once, in the order of the roots, the rows and the columns,
    /// passes over a cell that hosts nothing, and once ended asks for nothing more.
    #[test]
    fn the_walk_visits_every_cell_once_in_order_and_stays_ended() {
        let mut asked = Vec::new();
        let mut walk = Walk::new(2, 2);
        let mut yielded = Vec::new();
        while let Some(z) = walk.next_zone(|root, row, col| toy(&mut asked, root, row, col)) {
            yielded.push(z.0);
        }
        let mut cells = Vec::new();
        for root in 0..10 {
            cells.extend([(root, 0, 0), (root, 0, 1), (root, 1, 0), (root, 1, 1)]);
        }
        cells.extend([(10, 0, 0), (11, 0, 0)]);
        assert_eq!(asked, cells);
        let mut zones = Vec::new();
        for root in 0..10 {
            zones.extend([0, 1, 30, 31].map(|z| root * 100 + z));
        }
        zones.extend([1000, 1100]);
        assert_eq!(yielded, zones);
        for _ in 0..3 {
            assert_eq!(
                walk.next_zone(|root, row, col| toy(&mut asked, root, row, col)),
                None
            );
        }
        assert_eq!(asked.len(), 42);
    }

    /// Entered after a zone, the walk begins at that zone's cell, yields what the cell hosts
    /// beyond the zone and goes on; no cell before it is asked for.
    #[test]
    fn the_walk_entered_after_a_zone_begins_at_its_cell() {
        for (place, first_asked, expected) in [
            ((3, 1, 1, 330), (3, 1, 1), vec![331, 400, 401]),
            ((3, 1, 1, 331), (3, 1, 1), vec![400, 401, 430]),
            ((9, 1, 1, 931), (9, 1, 1), vec![1000, 1100]),
            ((10, 0, 0, 1000), (10, 0, 0), vec![1100]),
            ((11, 0, 0, 1100), (11, 0, 0), vec![]),
        ] {
            let mut asked = Vec::new();
            let mut walk = Walk::new(2, 2);
            // What the walk gave before it is entered is not counted.
            walk.next_zone(|root, row, col| toy(&mut asked, root, row, col));
            asked.clear();
            walk.enter_after(place);
            let mut yielded = Vec::new();
            while let Some(z) = walk.next_zone(|root, row, col| toy(&mut asked, root, row, col)) {
                yielded.push(z.0);
                if yielded.len() == 3 {
                    break;
                }
            }
            assert_eq!(yielded, expected, "{place:?}");
            assert_eq!(asked[0], first_asked, "{place:?}");
        }
    }

    /// A lattice of no side, which no level enumerated has, is one cell a root, and the walk
    /// ends.
    #[test]
    fn the_walk_of_a_lattice_of_no_side_ends() {
        let mut asked = 0;
        let mut walk = Walk::new(0, 0);
        let next = walk.next_zone(|_, _, _| {
            asked += 1;
            Vec::new()
        });
        assert_eq!(next, None);
        assert_eq!(asked, 12);
    }

    /// A repeatable stream of numbers (splitmix64), uniform in `[0, 1)`.
    struct Rng(u64);

    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
            z ^ (z >> 31)
        }
        fn unit(&mut self) -> f64 {
            (self.next() >> 11) as f64 / (1u64 << 53) as f64
        }
        fn range(&mut self, a: f64, b: f64) -> f64 {
            a + (b - a) * self.unit()
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n
        }
        /// A size spread evenly over the decades from `10^lo` to `10^hi`.
        fn decades(&mut self, lo: i32, hi: i32) -> f64 {
            let decade = lo + self.below((hi - lo) as u64) as i32;
            self.range(1.0, 10.0) * math::powi(10.0, decade)
        }
        /// A point uniform on the sphere, a latitude and a longitude in degrees.
        fn point(&mut self) -> (f64, f64) {
            let lat = math::asin(2.0 * self.unit() - 1.0) * math::RAD2DEG;
            (lat, self.range(-180.0, 180.0))
        }
    }

    /// The angle between two points of the sphere, each a latitude and a longitude in radians,
    /// by the haversine, which keeps a small angle.
    fn arc(a: (f64, f64), b: (f64, f64)) -> f64 {
        let s1 = math::sin((b.0 - a.0) / 2.0);
        let s2 = math::sin((b.1 - a.1) / 2.0);
        2.0 * math::asin(math::sqrt(s1 * s1 + math::cos(a.0) * math::cos(b.0) * s2 * s2).min(1.0))
    }

    /// The point at the angle `d` from the point `a`, along the bearing `bearing`, all in
    /// radians.
    fn displaced(a: (f64, f64), bearing: f64, d: f64) -> (f64, f64) {
        let lat = math::asin(
            (math::sin(a.0) * math::cos(d) + math::cos(a.0) * math::sin(d) * math::cos(bearing))
                .clamp(-1.0, 1.0),
        );
        let lon = a.1
            + math::atan2(
                math::sin(bearing) * math::sin(d) * math::cos(a.0),
                math::cos(d) - math::sin(a.0) * math::sin(lat),
            );
        (lat, lon)
    }

    #[test]
    fn the_net_is_measured_in_sides_of_a_triangle_of_a_twentieth_of_the_sphere() {
        // The side of an equilateral triangle of area Pi / 5.
        let side = math::sqrt(4.0 * PI / (5.0 * math::sqrt(3.0)));
        assert!((2.0 * HALF_SIDE - side).abs() < 1e-15);
        // The sides of a rhombus and its diagonal from the top left to the bottom right are
        // sides of triangles; the other diagonal is twice the height of one.
        for (dx, dy) in [(1.0, 0.0), (0.0, 1.0), (1.0, 1.0)] {
            assert!((net_length(dx, dy) - side).abs() < 1e-15);
        }
        assert!((net_length(1.0, -1.0) - math::sqrt(3.0) * side).abs() < 1e-15);
        assert_eq!(net_length(0.25, 0.0), net_length(1.0, 0.0) / 4.0);
    }

    /// The test of a point against a box grown by an angle never refuses a point that lies
    /// within that angle of the box: points of boxes of every sort, displaced by less than the
    /// angle in every direction, all pass.
    #[test]
    fn no_point_within_the_angle_of_a_box_is_found_far_from_it() {
        let mut rng = Rng(0x00b0_c5e5);
        let mut tried = 0;
        for i in 0..4000 {
            let (lat, lon) = rng.point();
            let half_height = rng.decades(-9, 2).min(50.0);
            let half_width = rng.decades(-9, 3).min(179.0);
            let wrap = |x: f64| x - 360.0 * ((x + 180.0) / 360.0).floor();
            let mut b = [
                (lat - half_height).max(-90.0),
                wrap(lon - half_width),
                (lat + half_height).min(90.0),
                wrap(lon + half_width),
            ];
            // Every fourth box has no height, every fifth no width, and some both; every
            // seventh has a side on a pole.
            if i % 4 == 0 {
                b[2] = b[0];
            }
            if i % 5 == 0 {
                b[3] = b[1];
            }
            if i % 7 == 0 {
                b[2] = 90.0;
            }
            let r = b.map(math::radians);
            let width = if r[3] >= r[1] {
                r[3] - r[1]
            } else {
                r[3] - r[1] + math::TWO_PI
            };
            for _ in 0..12 {
                let inside = (r[0] + (r[2] - r[0]) * rng.unit(), r[1] + width * rng.unit());
                let rho = rng.decades(-9, 0);
                let p = displaced(inside, rng.range(0.0, math::TWO_PI), rho * rng.unit());
                assert!(
                    near_box(p.0, p.1, rho, &r),
                    "the point {p:?}, within {rho} of {inside:?} of the box {r:?}"
                );
                tried += 1;
            }
        }
        assert_eq!(tried, 48_000);
    }

    /// And it refuses points far from the box: in latitude, in longitude on either side, and
    /// on the far side of a box that runs over the antimeridian.
    #[test]
    fn a_point_far_from_a_box_is_found_far_from_it() {
        let d = PI / 180.0;
        let b = [10.0 * d, 20.0 * d, 20.0 * d, 30.0 * d];
        let rho = d;
        assert!(near_box(15.0 * d, 25.0 * d, rho, &b));
        assert!(near_box(20.5 * d, 25.0 * d, rho, &b));
        assert!(!near_box(21.5 * d, 25.0 * d, rho, &b));
        assert!(!near_box(8.5 * d, 25.0 * d, rho, &b));
        assert!(near_box(15.0 * d, 30.9 * d, rho, &b));
        assert!(near_box(15.0 * d, 19.1 * d, rho, &b));
        assert!(!near_box(15.0 * d, 32.0 * d, rho, &b));
        assert!(!near_box(15.0 * d, 18.0 * d, rho, &b));
        assert!(!near_box(15.0 * d, -155.0 * d, rho, &b));
        // The longitude of the point is taken within one turn of the box's west.
        assert!(near_box(15.0 * d, (25.0 - 360.0) * d, rho, &b));
        assert!(!near_box(15.0 * d, (32.0 + 360.0) * d, rho, &b));
        // Over the antimeridian, from 170 to -170 degrees.
        let over = [10.0 * d, 170.0 * d, 20.0 * d, -170.0 * d];
        assert!(near_box(15.0 * d, 175.0 * d, rho, &over));
        assert!(near_box(15.0 * d, -175.0 * d, rho, &over));
        assert!(near_box(15.0 * d, 180.0 * d, rho, &over));
        assert!(near_box(15.0 * d, -169.1 * d, rho, &over));
        assert!(!near_box(15.0 * d, -168.0 * d, rho, &over));
        assert!(!near_box(15.0 * d, 168.0 * d, rho, &over));
        assert!(!near_box(15.0 * d, 0.0, rho, &over));
        // Beside a pole every longitude is near.
        let cap = [89.0 * d, 0.0, 90.0 * d, 10.0 * d];
        assert!(near_box(88.5 * d, 180.0 * d, rho, &cap));
        assert!(!near_box(87.5 * d, 180.0 * d, rho, &cap));
    }

    /// The greatest and the least ratio of a length on the sphere to the length of its image
    /// in the net, over `n * n` points of each rhombus and six directions at each.
    fn stretch<P: Projection, T: Topology, I: Indexing>(g: &Grid<P, T, I>, n: u32) -> (f64, f64) {
        let (mut least, mut greatest) = (f64::MAX, 0.0_f64);
        let h = 1e-5;
        for root in 0..10u32 {
            let (x0, y0) = (f64::from(root / 2), f64::from(root.div_ceil(2)));
            for i in 0..n {
                for j in 0..n {
                    let x = x0 + (f64::from(i) + 0.5) / f64::from(n);
                    let y = y0 + (f64::from(j) + 0.5) / f64::from(n);
                    let a = g.lattice_point(x, y).unwrap();
                    for (dx, dy) in [
                        (h, 0.0),
                        (0.0, h),
                        (h, h),
                        (h, -h),
                        (2.0 * h, h),
                        (h, 2.0 * h),
                    ] {
                        if x + dx >= x0 + 1.0 || y + dy >= y0 + 1.0 || y + dy <= y0 {
                            continue;
                        }
                        let q = g.lattice_point(x + dx, y + dy).unwrap();
                        let ratio = arc(a, q) / net_length(dx, dy);
                        least = least.min(ratio);
                        greatest = greatest.max(ratio);
                    }
                }
            }
        }
        (least, greatest)
    }

    /// The bound on the stretch of the inverse projection stands on these measured greatest
    /// values, one for each projection (1.165455, 1.106966 and 1.131024): a sample that exceeds
    /// one fails here, before the search that rests on the bound can lose a zone.
    #[test]
    fn the_stretch_of_each_projection_stays_under_its_measured_greatest() {
        let n = 60;
        for (name, (least, greatest), ceiling) in [
            ("ISEA", stretch(crate::isea3h(), n), 1.1655),
            ("IVEA", stretch(crate::ivea3h(), n), 1.1070),
            ("RTEA", stretch(crate::rtea3h(), n), 1.1311),
        ] {
            assert!(greatest <= ceiling, "{name}: a stretch of {greatest}");
            // The sample reaches the greatest to a thousandth.
            assert!(
                greatest > ceiling - 0.001,
                "{name}: a stretch of {greatest}"
            );
            assert!(least > 0.85, "{name}: a stretch of {least}");
            assert!(greatest < STRETCH, "{name}: a stretch of {greatest}");
        }
        // The projection is one on the grids of both apertures.
        assert_eq!(stretch(crate::igeo7(), 12), stretch(crate::isea3h(), 12));
    }

    /// The cells of a lattice of `edge` cells a side that the measurements below sample at a
    /// deep level: the four corners of every rhombus and the cells beside them, cells along
    /// its four sides, in the two rows and columns nearest each, and cells within it.
    fn sampled_cells(edge: u64, rng: &mut Rng) -> Vec<(u8, u64, u64)> {
        let mut cells = Vec::new();
        let last = edge - 1;
        let near = |k: u64| k.min(last);
        for root in 0..10u8 {
            for (row, col) in [(0, 0), (0, last), (last, 0), (last, last)] {
                for (dr, dc) in [(0, 0), (1, 0), (0, 1), (1, 1)] {
                    let r = if row == 0 { near(dr) } else { last - near(dr) };
                    let c = if col == 0 { near(dc) } else { last - near(dc) };
                    cells.push((root, r, c));
                }
            }
            for _ in 0..3 {
                for side in 0..2 {
                    let (k, along) = (near(side), rng.below(edge));
                    cells.extend([
                        (root, k, along),
                        (root, last - k, along),
                        (root, along, k),
                        (root, along, last - k),
                    ]);
                }
            }
            for _ in 0..4 {
                cells.push((root, rng.below(edge), rng.below(edge)));
            }
        }
        cells.sort_unstable();
        cells.dedup();
        cells
    }

    /// The zones the measurements below are taken over at `level`: every zone of the level
    /// where `whole`; else the zones of the sampled cells and of the two polar roots, the zones
    /// at forty points within four zone widths of each geographic pole, and the zones at four
    /// points where the zones of the ISEA grids were found widest.
    fn measured_zones<P: Projection, T: Topology, I: Indexing>(
        g: &Grid<P, T, I>,
        level: u8,
        whole: bool,
        rng: &mut Rng,
    ) -> Vec<ZoneId> {
        if whole {
            return g.zones(level).unwrap().collect();
        }
        let mut zones = Vec::new();
        for (root, row, col) in sampled_cells(T::lattice_edge(Private, level), rng) {
            zones.extend(g.lattice_cell(level, root, row, col));
        }
        for root in [10, 11] {
            zones.extend(g.lattice_cell(level, root, 0, 0));
        }
        let mut zones: Vec<ZoneId> = zones.into_iter().map(|(zone, _)| zone).collect();
        let width = math::sqrt(4.0 * PI / g.count_zones(level).unwrap() as f64) * math::RAD2DEG;
        let mut at = |lat: f64, lon: f64| {
            let zone = g.zone_from_geo(lat, lon, level).unwrap().id();
            assert_ne!(
                zone,
                ZoneId::NULL,
                "{} level {level}: ({lat}, {lon})",
                g.name()
            );
            zones.push(zone);
        };
        for pole in [90.0, -90.0_f64] {
            for _ in 0..40 {
                let away = width * rng.range(0.0, 4.0);
                at(pole - pole.signum() * away, rng.range(-180.0, 180.0));
            }
        }
        for (lat, lon) in [
            (15.4257, -178.2479),
            (-16.4049, 19.0267),
            (-12.2145, 176.7370),
            (-15.5317, 20.5334),
        ] {
            at(lat, lon);
        }
        zones.sort_unstable();
        zones.dedup();
        zones
    }

    /// The greatest radius of a zone and the greatest distance of a zone from the corner of
    /// its cell, over the measured zones of every level of `levels`, every zone of a level to
    /// `whole_to` and a sample beyond. The radius is the greatest angle from the zone's
    /// centroid to a point of the boundary its extent is taken from: it is given in zone
    /// widths, and as an excess in radians over `radius` zone widths. The hosting distance is
    /// given in sides of a cell. Last, the number of zones measured.
    fn measure<P: Projection, T: Topology, I: Indexing>(
        g: &Grid<P, T, I>,
        levels: std::ops::RangeInclusive<u8>,
        whole_to: u8,
        radius: f64,
    ) -> (f64, f64, f64, usize) {
        let mut rng = Rng(0x5eed_0000 + u64::from(*levels.start()));
        let (mut widths, mut excess, mut sides, mut count) = (0.0_f64, f64::MIN, 0.0_f64, 0);
        for level in levels {
            let width = math::sqrt(4.0 * PI / g.count_zones(level).unwrap() as f64);
            let edge = T::lattice_edge(Private, level) as f64;
            let side = net_length(1.0 / edge, 0.0);
            for zone in measured_zones(g, level, level <= whole_to, &mut rng) {
                let centre = g.centroid(zone);
                let centre = (math::radians(centre.lat), math::radians(centre.lon));
                let ring = g.refined_vertices(zone, 0).unwrap();
                assert!(!ring.is_empty(), "{}", g.text_id(zone));
                let r = ring
                    .iter()
                    .map(|p| arc(centre, (math::radians(p.lat), math::radians(p.lon))))
                    .fold(0.0, f64::max);
                widths = widths.max(r / width);
                excess = excess.max(r - radius * width);
                let (root, row, col, _) = T::locate(Private, &I::decode(zone)).unwrap();
                if root < 10 {
                    let x = f64::from(root / 2) + col as f64 / edge;
                    let y = f64::from(root.div_ceil(2)) + row as f64 / edge;
                    let corner = g.lattice_point(x, y).unwrap();
                    sides = sides.max(arc(corner, centre) / side);
                }
                count += 1;
            }
        }
        (widths, excess, sides, count)
    }

    /// The greatest radius measured, in zone widths: 0.72233, on the ISEA grids of both
    /// apertures, where the projection stretches most.
    const GREATEST_RADIUS: f64 = 0.7224;

    /// The greatest excess measured over that radius at the three finest levels of the
    /// aperture-3 grids, in radians: 1.474e-8, at the zones upon the geographic poles.
    const GREATEST_EXCESS: f64 = 1.5e-8;

    /// What was measured lies under the bounds the search prunes by: a radius in zone widths,
    /// its excess in radians where the radius is taken at [`GREATEST_RADIUS`], and a distance
    /// from the corner of the cell in sides.
    fn under_the_bounds_of_the_search(name: &str, widths: f64, excess: f64, sides: f64) {
        assert!(
            widths < RADIUS_IN_ZONE_WIDTHS,
            "{name}: a radius of {widths}"
        );
        assert!(excess < RADIUS_FLOOR, "{name}: an excess of {excess}");
        assert!(sides < HOSTING_IN_SIDES, "{name}: {sides} sides");
    }

    /// To level 30 a zone of the aperture-3 grids lies within [`GREATEST_RADIUS`] widths of
    /// its centroid, and to level 32 within 0.672 sides of the corner of its cell. At the
    /// three finest levels a zone is 4.5e-8 to 1.5e-8 radian wide, and the zones upon the two
    /// geographic poles, whose coordinates the projection rounds by about as much, lie
    /// farther: by [`GREATEST_EXCESS`] at most, which the floor of the radius carries; at
    /// level 33 one of them lies 0.915 sides from the corner of its cell. `floor` is the
    /// radius that the lists of this grid must reach, so that the measurement is one.
    fn aperture_3_stays_under_what_was_measured<P: Projection, T: Topology, I: Indexing>(
        g: &Grid<P, T, I>,
        floor: f64,
    ) {
        let name = g.name();
        let (widths, _, sides, count) = measure(g, 0..=30, 5, GREATEST_RADIUS);
        under_the_bounds_of_the_search(name, widths, 0.0, sides);
        assert!(widths <= GREATEST_RADIUS, "{name}: a radius of {widths}");
        assert!(widths > floor, "{name}: a radius of {widths}");
        assert!(sides <= 0.672, "{name}: {sides} sides from the corner");
        assert!(sides > 0.63, "{name}: {sides} sides from the corner");
        assert!(count >= 26_000, "{name}: {count} zones");
        let (_, excess, sides, count) = measure(g, 31..=32, 0, GREATEST_RADIUS);
        under_the_bounds_of_the_search(name, GREATEST_RADIUS, excess, sides);
        assert!(excess <= GREATEST_EXCESS, "{name}: an excess of {excess}");
        assert!(excess > 6e-9, "{name}: an excess of {excess}");
        assert!(sides <= 0.672, "{name}: {sides} sides from the corner");
        assert!(count >= 500, "{name}: {count} zones");
        let (_, excess, sides, count) = measure(g, 33..=33, 0, GREATEST_RADIUS);
        under_the_bounds_of_the_search(name, GREATEST_RADIUS, excess, sides);
        assert!(excess <= GREATEST_EXCESS, "{name}: an excess of {excess}");
        assert!(excess > 1.4e-8, "{name}: an excess of {excess}");
        assert!(sides <= 0.915, "{name}: {sides} sides from the corner");
        assert!(sides > 0.89, "{name}: {sides} sides from the corner");
        assert!(count >= 500, "{name}: {count} zones");
    }

    /// The bound on the radius of a zone and the bound on its distance from the corner of its
    /// cell stand on these measured greatest values: a zone of the lists that exceeds one fails
    /// here, before the search that rests on the bounds can lose a zone.
    #[test]
    fn the_zones_of_aperture_3_stay_under_the_measured_radius_and_hosting_distance() {
        aperture_3_stays_under_what_was_measured(crate::isea3h(), 0.722);
        aperture_3_stays_under_what_was_measured(crate::ivea3h(), 0.685);
        aperture_3_stays_under_what_was_measured(crate::rtea3h(), 0.70);
    }

    /// On the aperture-7 grids, to the level at which their zones are enumerated: within
    /// [`GREATEST_RADIUS`] widths of the centroid and within 0.44 sides of the corner of the
    /// cell.
    fn aperture_7_stays_under_what_was_measured<P: Projection, T: Topology, I: Indexing>(
        g: &Grid<P, T, I>,
        floor: f64,
    ) {
        let name = g.name();
        assert_eq!(g.max_box_level(), 14);
        let (widths, _, sides, count) = measure(g, 0..=14, 3, GREATEST_RADIUS);
        under_the_bounds_of_the_search(name, widths, 0.0, sides);
        assert!(widths <= GREATEST_RADIUS, "{name}: a radius of {widths}");
        assert!(widths > floor, "{name}: a radius of {widths}");
        assert!(sides <= 0.44, "{name}: {sides} sides from the corner");
        assert!(sides > 0.41, "{name}: {sides} sides from the corner");
        assert!(count >= 22_000, "{name}: {count} zones");
    }

    #[test]
    fn the_zones_of_aperture_7_stay_under_the_measured_radius_and_hosting_distance() {
        aperture_7_stays_under_what_was_measured(crate::igeo7(), 0.722);
        aperture_7_stays_under_what_was_measured(crate::ivea7h(), 0.685);
        aperture_7_stays_under_what_was_measured(crate::rtea7h(), 0.70);
    }

    /// The first zones of a large box at a deep level cost a search of the lattice, and not
    /// its rim: the first three hundred zones of boxes whose rim rises to the right across a
    /// rhombus, for hundreds of thousands of cells, are had in some ten to twenty thousand tests
    /// of a block. (A search of a rhombus for its first row that takes the upper and the left
    /// halves first, and keeps the best row found, costs here six million tests for the first
    /// zone alone, and as many again at each row.)
    #[test]
    fn the_first_zones_of_a_large_box_cost_a_search_and_not_its_rim() {
        fn tests_of_a_block<P: Projection, T: Topology, I: Indexing>(
            g: &Grid<P, T, I>,
            level: u8,
            b: [f64; 4],
        ) -> usize {
            let (mut walk, search) = g
                .box_search(level, &degrees(b[0], b[1], b[2], b[3]))
                .unwrap();
            let tests = std::cell::Cell::new(0);
            let mut zones = 0;
            while zones < 300
                && walk
                    .next_zone_in_box(
                        &search,
                        |x, y| {
                            tests.set(tests.get() + 1);
                            g.lattice_point(x, y)
                        },
                        |root, row, col| g.lattice_cell(level, root, row, col),
                        |zone| g.zone_meets_box(zone, &search.bbox, search.reach),
                    )
                    .is_some()
            {
                zones += 1;
            }
            assert_eq!(zones, 300, "{} level {level}, the box {b:?}", g.name());
            tests.get()
        }
        let south = [
            -90.0,
            56.478892231624016,
            -42.34771473210071,
            149.54716528164403,
        ];
        let east = [
            -48.2217911335482,
            131.17959902814255,
            8.496604404680816,
            168.33992068981183,
        ];
        let north = [
            8.770311054338151,
            -150.56808348646433,
            90.0,
            -134.3430224673964,
        ];
        for (name, tests) in [
            ("RTEA3H 26", tests_of_a_block(crate::rtea3h(), 26, south)),
            ("RTEA3H 33", tests_of_a_block(crate::rtea3h(), 33, south)),
            ("ISEA3H 32", tests_of_a_block(crate::isea3h(), 32, east)),
            ("IVEA7H 14", tests_of_a_block(crate::ivea7h(), 14, north)),
        ] {
            assert!(tests < 60_000, "{name}: {tests} tests of a block");
            assert!(tests > 300, "{name}: {tests} tests of a block");
        }
    }

    /// A box of the test's own, in degrees: `[south, west, north, east]`, of `kind` 0 to 7,
    /// about a point of the sphere, its sides measured in `width`, the width of a zone in
    /// degrees.
    fn a_box(kind: u64, width: f64, rng: &mut Rng) -> [f64; 4] {
        let (lat, lon) = rng.point();
        let wrap = |x: f64| (x - 360.0 * ((x + 180.0) / 360.0).floor()).clamp(-180.0, 180.0);
        let about = |lat: f64, lon: f64, h: f64, w: f64| {
            let w = (w / math::cos(math::radians(lat)).max(0.05)).min(179.0);
            [
                (lat - h).max(-90.0),
                wrap(lon - w),
                (lat + h).min(90.0),
                wrap(lon + w),
            ]
        };
        let size = |rng: &mut Rng, lo: i32, hi: i32| (width * rng.decades(lo, hi)).min(85.0);
        match kind {
            // Sides from a tenth of a zone to ten zones.
            0 => {
                let (h, w) = (size(rng, -1, 1), size(rng, -1, 1));
                about(lat, lon, h, w)
            }
            // Far smaller than a zone.
            1 => {
                let (h, w) = (size(rng, -6, -1), size(rng, -6, -1));
                about(lat, lon, h, w)
            }
            // A large part of the sphere.
            2 => about(lat, lon, rng.range(3.0, 80.0), rng.range(3.0, 179.0)),
            // Over the antimeridian.
            3 => {
                let (h, w) = (size(rng, -1, 1), size(rng, -1, 1));
                about(lat, 180.0 - rng.range(-1.0, 1.0) * width, h, w)
            }
            // A side on a pole, about the whole circle or a part of it.
            4 => {
                let h = size(rng, -2, 1);
                let (west, east) = if rng.unit() < 0.4 {
                    (-180.0, 180.0)
                } else {
                    (rng.range(-180.0, 180.0), rng.range(-180.0, 180.0))
                };
                if rng.unit() < 0.5 {
                    [90.0 - h, west, 90.0, east]
                } else {
                    [-90.0, west, -90.0 + h, east]
                }
            }
            // A point.
            5 => [lat, lon, lat, lon],
            // A segment of a parallel, or of a meridian.
            6 => {
                let l = size(rng, -1, 1);
                if rng.unit() < 0.5 {
                    about(lat, lon, 0.0, l)
                } else {
                    about(lat, lon, l, 0.0)
                }
            }
            // A side upon one of the meridians and the parallel along which sides of rhombi
            // run: the meridians of 11.2 and -168.8 degrees, and those of -78.8 and 101.2.
            _ => {
                let meridian = [11.2, -168.8, -78.8, 101.2][rng.below(4) as usize];
                let (h, w) = (size(rng, -1, 1), size(rng, -1, 1));
                let mut b = about(lat, meridian, h, w);
                if rng.unit() < 0.5 {
                    b[1] = meridian;
                } else {
                    b[3] = meridian;
                }
                b
            }
        }
    }

    /// The search leaves out no zone of the answer: for boxes of every sort, at the levels at
    /// which every zone of a level can be tested, the walk of the box gives exactly the zones
    /// of the level whose extent meets it, in the order of the level, where every zone is put
    /// to the test of its extent and none is left out by its distance.
    fn the_search_leaves_out_nothing<P: Projection, T: Topology, I: Indexing>(
        g: &Grid<P, T, I>,
        levels: std::ops::RangeInclusive<u8>,
        boxes_a_kind: u64,
    ) -> usize {
        let mut rng = Rng(0x0b0c_5000 + u64::from(T::APERTURE));
        let mut zones_met = 0;
        for level in levels {
            // Every zone of the level with its extent as the public surface gives it, grown
            // by far more than the conversion to degrees and back can move it: a zone whose
            // grown extent misses the box cannot meet it, and is spared the exact test.
            let whole: Vec<(ZoneId, [f64; 4])> = g
                .zones(level)
                .unwrap()
                .map(|z| {
                    let e = g.extent(z).unwrap();
                    let [s, w, n, e] = [e.ll.lat, e.ll.lon, e.ur.lat, e.ur.lon].map(math::radians);
                    (z, [s - 1e-12, w - 1e-12, n + 1e-12, e + 1e-12])
                })
                .collect();
            let width = math::sqrt(4.0 * PI / whole.len() as f64) * math::RAD2DEG;
            for kind in 0..8 {
                for _ in 0..boxes_a_kind {
                    let b = a_box(kind, width, &mut rng);
                    let bbox = degrees(b[0], b[1], b[2], b[3]);
                    let r = box_radians(&bbox).unwrap();
                    let expected: Vec<ZoneId> = whole
                        .iter()
                        .filter(|(z, grown)| {
                            intersects(grown, &r) && g.zone_meets_box(*z, &r, f64::INFINITY)
                        })
                        .map(|(z, _)| *z)
                        .collect();
                    let got: Vec<ZoneId> = g.zones_in_box(level, &bbox).unwrap().collect();
                    assert_eq!(got, expected, "{} level {level}, the box {b:?}", g.name());
                    zones_met += got.len();
                }
            }
        }
        zones_met
    }

    #[test]
    fn the_search_leaves_out_no_zone_of_a_box_on_aperture_3() {
        let met = the_search_leaves_out_nothing(crate::isea3h(), 0..=5, 12)
            + the_search_leaves_out_nothing(crate::ivea3h(), 0..=5, 12)
            + the_search_leaves_out_nothing(crate::rtea3h(), 0..=5, 12);
        assert!(met >= 56_000, "{met}");
    }

    #[test]
    fn the_search_leaves_out_no_zone_of_a_box_on_aperture_7() {
        let met = the_search_leaves_out_nothing(crate::igeo7(), 0..=3, 12)
            + the_search_leaves_out_nothing(crate::ivea7h(), 0..=3, 12)
            + the_search_leaves_out_nothing(crate::rtea7h(), 0..=3, 12);
        assert!(met >= 63_000, "{met}");
    }

    /// The zones of the finest level of aperture 3, the level at which the arithmetic of the
    /// estimate has least room.
    const ZONES_OF_LEVEL_33: u64 = 55_590_605_665_555_232;

    /// Beside a pole the two sines are one number where the band between the two latitudes
    /// is not empty: the difference taken as a product keeps it, to the last digits, and
    /// away from the poles it is the difference of the sines.
    #[test]
    fn the_difference_of_two_sines_does_not_cancel_beside_a_pole() {
        // A band of half a zone of the finest level of aperture 3, on the pole.
        assert_eq!(
            math::sin(FRAC_PI_2) - math::sin(FRAC_PI_2 - 7.5e-9),
            0.0,
            "the two sines are one number"
        );
        for x in [7.5e-9, 1e-8, 3e-9, 2.5e-8, 1e-6] {
            // The latitude `x` from the pole, as a `f64` holds it, and its distance from the
            // pole, which is exact; `1 - cos(d)` is `d * d / 2` at this size. The mean of the
            // two latitudes is held to a unit in its last place, 2.2e-16, which is that part
            // of half the distance.
            let near = FRAC_PI_2 - x;
            let d = FRAC_PI_2 - near;
            let expected = d * d / 2.0;
            let got = sine_difference(near, FRAC_PI_2);
            assert!((got / expected - 1.0).abs() < 1e-15 / d, "{x}: {got}");
            let got = sine_difference(-FRAC_PI_2, -near);
            assert!((got / expected - 1.0).abs() < 1e-15 / d, "{x}: {got}");
        }
        for (a, b) in [(-1.2, 0.3), (0.0, FRAC_PI_2), (0.5, 0.5), (-0.7, -0.2)] {
            let (got, expected) = (sine_difference(a, b), math::sin(b) - math::sin(a));
            assert!((got - expected).abs() < 1e-15, "{a} {b}: {got}");
        }
    }

    /// The area of the grown box in closed form is the integral it stands for: the width of
    /// the grown box at each latitude, `width + 2 margin / cos(lat)` and at most the whole
    /// circle, times `cos(lat)`, summed here over thin bands. Boxes of every sort: plain,
    /// over the antimeridian, on a pole, of no width, of no height, the whole circle, and
    /// with a margin beside which the box is small or large.
    #[test]
    fn the_area_of_the_grown_box_is_the_integral_of_its_width() {
        let boxes = [
            [0.66, -0.17, 0.68, -0.15],
            [-0.3, 3.0, 0.4, -3.0],
            [1.2, -PI, FRAC_PI_2, PI],
            [-FRAC_PI_2, 0.2, -1.5, 2.9],
            [0.1, 1.0, 0.9, 1.0],
            [-0.5, -2.0, -0.5, 2.0],
            [1.57, 0.3, 1.57, 0.3],
            [-1.0, -3.1, 1.0, 3.1],
            [-FRAC_PI_2, -PI, FRAC_PI_2, PI],
            [-0.2, 0.5, 0.2, 0.4],
        ];
        for b in &boxes {
            for margin in [0.0, 1e-4, 0.01, 0.3, 2.0] {
                let south = (b[0] - margin).max(-FRAC_PI_2);
                let north = (b[2] + margin).min(FRAC_PI_2);
                let width = if b[3] < b[1] {
                    b[3] - b[1] + 2.0 * PI
                } else {
                    b[3] - b[1]
                };
                let steps = 200_000;
                let step = (north - south) / f64::from(steps);
                let summed: f64 = (0..steps)
                    .map(|i| {
                        let cos = math::cos(south + (f64::from(i) + 0.5) * step);
                        (width * cos + 2.0 * margin).min(2.0 * PI * cos) * step
                    })
                    .sum();
                let got = grown_area(b, margin);
                assert!(
                    (got - summed).abs() <= 1e-6 * (1.0 + summed),
                    "{b:?} grown by {margin}: {got}, summed {summed}"
                );
                assert!(got <= 4.0 * PI * (1.0 + 4.0 * f64::EPSILON));
            }
        }
    }

    /// A cap about a pole that is smaller than a zone, at the finest level: the estimate is
    /// the area of the cap grown by the margin, a disc of `radius + 1.75` widths of a zone,
    /// rounded up, and the allowance. The sines of the two latitudes of such a cap differ by
    /// three units in the last place or fewer, each worth three zones.
    #[test]
    fn a_cap_about_a_pole_smaller_than_a_zone_is_estimated_by_its_area() {
        let width = math::sqrt(4.0 * PI / ZONES_OF_LEVEL_33 as f64);
        for tenths in 0..=30 {
            let radius = f64::from(tenths) / 10.0;
            let disc = PI * (radius + ESTIMATE_MARGIN_IN_ZONE_WIDTHS) * (radius + 1.75);
            for cap in [
                [FRAC_PI_2 - radius * width, -PI, FRAC_PI_2, PI],
                [-FRAC_PI_2, -PI, radius * width - FRAC_PI_2, PI],
            ] {
                let over = estimate(&cap, ZONES_OF_LEVEL_33) as f64 - disc - 3.0;
                assert!(
                    (-1e-6..1.0 + 1e-6).contains(&over),
                    "a cap of {radius} widths: {over} above its disc"
                );
            }
        }
    }

    /// The whole sphere is every zone of the level and the allowance, at every level of
    /// either aperture; and a share of the sphere that the rounding of the arithmetic left
    /// a unit in the last place under the whole still counts every zone, at the levels whose
    /// count a `f64` does not hold to the unit.
    #[test]
    fn the_whole_sphere_is_every_zone_of_the_level() {
        let sphere = [-FRAC_PI_2, -PI, FRAC_PI_2, PI];
        for (aperture, levels) in [(3, 0..=33), (7, 0..=19)] {
            for level in levels {
                let zones = crate::facts::count_zones(aperture, level);
                assert_eq!(estimate(&sphere, zones), zones + ESTIMATE_ALLOWANCE);
                assert!(share_in_zones(1.0 - f64::EPSILON / 2.0, zones) >= zones);
                assert!(share_in_zones(1.0 - 2.0 * f64::EPSILON, zones) >= zones);
            }
        }
        assert_eq!(crate::facts::count_zones(3, 33), ZONES_OF_LEVEL_33);
        // Without that room the share one unit under the whole loses eight zones there.
        assert_eq!(
            (ZONES_OF_LEVEL_33 as f64 * (1.0 - f64::EPSILON / 2.0)) as u64,
            ZONES_OF_LEVEL_33 - 8
        );
    }

    /// A box of no width exactly on the antimeridian, of either sign, and a box of no height
    /// exactly on either pole yield no zone: the pieces of an extent cut at the antimeridian
    /// end there, and no extent reaches beyond a pole, so that no overlap exceeds the margin.
    /// A segment a thousandth of a degree west of the antimeridian yields zones.
    #[test]
    fn a_box_of_no_area_on_the_antimeridian_or_on_a_pole_yields_no_zone() {
        fn count<P: Projection, T: Topology, I: Indexing>(
            g: &Grid<P, T, I>,
            s: f64,
            w: f64,
            n: f64,
            e: f64,
        ) -> usize {
            g.zones_in_box(6, &degrees(s, w, n, e)).unwrap().count()
        }
        fn none_on_the_lines<P: Projection, T: Topology, I: Indexing>(g: &Grid<P, T, I>) {
            for lon in [180.0, -180.0] {
                assert_eq!(count(g, -90.0, lon, 90.0, lon), 0, "{} at {lon}", g.name());
            }
            for lat in [90.0, -90.0] {
                assert_eq!(
                    count(g, lat, -180.0, lat, 180.0),
                    0,
                    "{} at {lat}",
                    g.name()
                );
            }
            assert!(count(g, -90.0, 179.999, 90.0, 179.999) > 0, "{}", g.name());
        }
        none_on_the_lines(crate::isea3h());
        none_on_the_lines(crate::igeo7());
    }
}
