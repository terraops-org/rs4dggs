//! The grids this crate ships, built once on first use and looked up by name.
//! Ported from py4dggs's `registry.py`, which builds all six eagerly at import
//! time; here each is a `OnceLock`, so a caller pays for the projection
//! geometry of the grid it actually asks for.
//!
//! Source: py4dggs's `registry.py`; `AnyGrid`, the handle on a grid chosen by name at run
//! time, is this crate's own.

use std::sync::OnceLock;

use crate::disk::Disk;
use crate::grid::Grid;
use crate::indexings::{I3h, Z7};
use crate::projections::{Isea, Ivea, Rtea};
use crate::topologies::{HexA3, HexA7};
use crate::{Error, Extent, GeoPoint, GridConfig, Result, ZoneId};

/// A registry grid on the canonical orientation.
///
/// `Grid::new` refuses only a configuration with a NaN or infinite angle, and
/// `GridConfig::default` holds three finite literals, which
/// `types::tests::grid_config_defaults_are_igeo7` pins; so this cannot fail, and
/// the `expect` states an invariant of this crate rather than handling input.
fn canonical<P: crate::Projection, T: crate::Topology, I: crate::Indexing>(
    name: &'static str,
) -> Grid<P, T, I> {
    Grid::new(GridConfig::default(), name).expect("the canonical orientation is finite")
}

/// IGEO7: the ISEA projection, the aperture-7 hexagonal topology and the Z7
/// indexing, on DGGAL's canonical orientation.
pub type Igeo7 = Grid<Isea, HexA7, Z7>;

/// The IGEO7 grid, built on first use and shared thereafter.
pub fn igeo7() -> &'static Igeo7 {
    static G: OnceLock<Igeo7> = OnceLock::new();
    G.get_or_init(|| canonical("IGEO7"))
}

/// IVEA7H: the IVEA projection, the aperture-7 hexagonal topology and the Z7
/// indexing, on DGGAL's canonical orientation.
pub type Ivea7h = Grid<Ivea, HexA7, Z7>;

/// The IVEA7H grid, built on first use and shared thereafter.
pub fn ivea7h() -> &'static Ivea7h {
    static G: OnceLock<Ivea7h> = OnceLock::new();
    G.get_or_init(|| canonical("IVEA7H"))
}

/// RTEA7H: the RTEA projection, the aperture-7 hexagonal topology and the Z7
/// indexing, on DGGAL's canonical orientation.
pub type Rtea7h = Grid<Rtea, HexA7, Z7>;

/// The RTEA7H grid, built on first use and shared thereafter.
pub fn rtea7h() -> &'static Rtea7h {
    static G: OnceLock<Rtea7h> = OnceLock::new();
    G.get_or_init(|| canonical("RTEA7H"))
}

/// ISEA3H: the ISEA projection, the aperture-3 hexagonal topology and the I3H
/// indexing, on DGGAL's canonical orientation.
pub type Isea3h = Grid<Isea, HexA3, I3h>;

/// The ISEA3H grid, built on first use and shared thereafter.
pub fn isea3h() -> &'static Isea3h {
    static G: OnceLock<Isea3h> = OnceLock::new();
    G.get_or_init(|| canonical("ISEA3H"))
}

/// IVEA3H: the IVEA projection, the aperture-3 hexagonal topology and the I3H
/// indexing, on DGGAL's canonical orientation.
pub type Ivea3h = Grid<Ivea, HexA3, I3h>;

/// The IVEA3H grid, built on first use and shared thereafter.
pub fn ivea3h() -> &'static Ivea3h {
    static G: OnceLock<Ivea3h> = OnceLock::new();
    G.get_or_init(|| canonical("IVEA3H"))
}

/// RTEA3H: the RTEA projection, the aperture-3 hexagonal topology and the I3H
/// indexing, on DGGAL's canonical orientation.
pub type Rtea3h = Grid<Rtea, HexA3, I3h>;

/// The RTEA3H grid, built on first use and shared thereafter.
pub fn rtea3h() -> &'static Rtea3h {
    static G: OnceLock<Rtea3h> = OnceLock::new();
    G.get_or_init(|| canonical("RTEA3H"))
}

/// A grid chosen at run time, by name rather than by type.
///
/// A `Grid` composes its three parts as type parameters, so two grids are two
/// types and a caller that only learns which grid it wants at run time cannot
/// name one. This enum is that caller's handle, and it offers the grid's own
/// operations, taking and returning identifiers rather than zones. It names the
/// six grids of this crate: the three of aperture 7, IGEO7, IVEA7H and RTEA7H,
/// and the three of aperture 3, ISEA3H, IVEA3H and RTEA3H. The enum is
/// `non_exhaustive`, so that a grid added later does not break a caller's `match`.
#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub enum AnyGrid {
    /// IGEO7, as [`igeo7`] gives it.
    Igeo7(&'static Igeo7),
    /// IVEA7H, as [`ivea7h`] gives it.
    Ivea7h(&'static Ivea7h),
    /// RTEA7H, as [`rtea7h`] gives it.
    Rtea7h(&'static Rtea7h),
    /// ISEA3H, as [`isea3h`] gives it.
    Isea3h(&'static Isea3h),
    /// IVEA3H, as [`ivea3h`] gives it.
    Ivea3h(&'static Ivea3h),
    /// RTEA3H, as [`rtea3h`] gives it.
    Rtea3h(&'static Rtea3h),
}

/// The canonical names of the grids this crate offers, in the order [`get_grid`]
/// tries them: what a grid's own `name()` returns, and what [`Error::UnknownGrid`]'s
/// message lists. [`get_grid`] itself accepts any of these regardless of case; this
/// constant carries only the canonical spelling.
pub const GRID_NAMES: &[&str] = &["IGEO7", "IVEA7H", "RTEA7H", "ISEA3H", "IVEA3H", "RTEA3H"];

/// The grid of this name, matched without regard to case, so `"igeo7"`, `"IGEO7"`
/// and `"Igeo7"` all name the same grid, or [`Error::UnknownGrid`] naming the
/// grids this crate does offer, [`GRID_NAMES`]. A grid's own `name()` always
/// answers in its canonical form regardless of how the caller spelt it here.
///
/// # Examples
///
/// ```
/// let grid = rs4dggs::get_grid("igeo7")?; // any casing
/// assert_eq!(grid.name(), "IGEO7");
/// assert_eq!(
///     rs4dggs::GRID_NAMES,
///     ["IGEO7", "IVEA7H", "RTEA7H", "ISEA3H", "IVEA3H", "RTEA3H"]
/// );
///
/// // An unknown name is refused, and the message lists the names there are.
/// let err = rs4dggs::get_grid("H3").unwrap_err();
/// assert!(err.to_string().contains("IGEO7"));
/// # Ok::<(), rs4dggs::Error>(())
/// ```
pub fn get_grid(name: &str) -> Result<AnyGrid> {
    if name.eq_ignore_ascii_case("IGEO7") {
        Ok(AnyGrid::Igeo7(igeo7()))
    } else if name.eq_ignore_ascii_case("IVEA7H") {
        Ok(AnyGrid::Ivea7h(ivea7h()))
    } else if name.eq_ignore_ascii_case("RTEA7H") {
        Ok(AnyGrid::Rtea7h(rtea7h()))
    } else if name.eq_ignore_ascii_case("ISEA3H") {
        Ok(AnyGrid::Isea3h(isea3h()))
    } else if name.eq_ignore_ascii_case("IVEA3H") {
        Ok(AnyGrid::Ivea3h(ivea3h()))
    } else if name.eq_ignore_ascii_case("RTEA3H") {
        Ok(AnyGrid::Rtea3h(rtea3h()))
    } else {
        Err(Error::UnknownGrid(crate::error::truncated(name)))
    }
}

/// Forwards one operation to whichever grid the variant holds. Written once
/// here so that adding a grid is one more arm in one place rather than one
/// more arm in every method below.
macro_rules! dispatch {
    ($self:ident, $g:ident => $e:expr) => {
        match $self {
            AnyGrid::Igeo7($g) => $e,
            AnyGrid::Ivea7h($g) => $e,
            AnyGrid::Rtea7h($g) => $e,
            AnyGrid::Isea3h($g) => $e,
            AnyGrid::Ivea3h($g) => $e,
            AnyGrid::Rtea3h($g) => $e,
        }
    };
}

impl AnyGrid {
    /// The grid's DGGRS name.
    pub fn name(&self) -> &'static str {
        dispatch!(self, g => g.name())
    }

    /// The finest resolution the grid's indexing can pack.
    pub fn max_resolution(&self) -> u8 {
        dispatch!(self, g => g.max_resolution())
    }

    /// The identifier of the zone containing the given point at resolution
    /// `res`. Identifiers, rather than zones, because a `Zone` borrows the
    /// grid it was made from and this handle hides which grid that is.
    ///
    /// Refuses a NaN or infinite coordinate, where the engine answers one fixed
    /// pentagon; see [`Grid::zone_from_geo`].
    pub fn zone_from_geo(&self, lat: f64, lon: f64, res: u8) -> Result<ZoneId> {
        dispatch!(self, g => g.zone_from_geo(lat, lon, res).map(|z| z.id()))
    }

    /// The identifier named by a canonical text id.
    ///
    /// Refuses text past the finest resolution, which the engine reads as a zone
    /// without geometry, and quotes rejected text at bounded length; see
    /// [`Grid::zone_from_text`]. The text carries no grid, and read on another of
    /// the grids names another cell; see [`Grid::zone`].
    pub fn zone_from_text(&self, text: &str) -> Result<ZoneId> {
        dispatch!(self, g => g.zone_from_text(text).map(|z| z.id()))
    }

    /// The identifier's canonical text form.
    pub fn text_id(&self, id: ZoneId) -> String {
        dispatch!(self, g => g.text_id(id))
    }

    /// The resolution encoded in the identifier; see [`Grid::resolution`] for
    /// the null zone.
    pub fn resolution(&self, id: ZoneId) -> u8 {
        dispatch!(self, g => g.resolution(id))
    }

    /// Whether the identifier names a pentagon.
    pub fn is_pentagon(&self, id: ZoneId) -> bool {
        dispatch!(self, g => g.is_pentagon(id))
    }

    /// The zone's centroid in WGS84 degrees.
    pub fn centroid(&self, id: ZoneId) -> GeoPoint {
        dispatch!(self, g => g.centroid(id))
    }

    /// The zone's boundary vertices in WGS84 degrees, anticlockwise as seen from outside
    /// the sphere on every grid; see [`Grid::vertices`] for the order, the antimeridian
    /// and the vertices that can coincide.
    pub fn vertices(&self, id: ZoneId) -> Vec<GeoPoint> {
        dispatch!(self, g => g.vertices(id))
    }

    /// The zone's area in square metres, or `None` for a zone without geometry; see
    /// [`Grid::area`].
    pub fn area(&self, id: ZoneId) -> Option<f64> {
        dispatch!(self, g => g.area(id))
    }

    /// The zone's boundary in WGS84 degrees with every edge divided into `edge_refinement`
    /// parts, 0 for the engine's own choice, or the empty ring for a zone without geometry;
    /// see [`Grid::refined_vertices`] for the sequence, the longitudes, which are continuous
    /// beyond 180 degrees, and the refinement that is refused: one above
    /// [`crate::grid::max_edge_refinement`], which `RS4DGGS_MAX_EDGE_REFINEMENT` may lower from
    /// 100,000 and which never refuses 0.
    pub fn refined_vertices(&self, id: ZoneId, edge_refinement: u32) -> Result<Vec<GeoPoint>> {
        dispatch!(self, g => g.refined_vertices(id, edge_refinement))
    }

    /// For the tests of the ring itself, which hold whatever the edge-refinement limit is.
    #[cfg(test)]
    pub(crate) fn refined_vertices_within_limit(
        &self,
        id: ZoneId,
        edge_refinement: u32,
    ) -> Vec<GeoPoint> {
        dispatch!(self, g => g.refined_vertices_within_limit(id, edge_refinement))
    }

    /// The zone's geographic extent in WGS84 degrees, or `None` for a zone without geometry;
    /// see [`Grid::extent`], and [`Extent`] for an extent over the antimeridian.
    pub fn extent(&self, id: ZoneId) -> Option<Extent> {
        dispatch!(self, g => g.extent(id))
    }

    /// The zone's edge neighbours; see [`Grid::neighbors`], which also says where the
    /// relation is not symmetric.
    pub fn neighbors(&self, id: ZoneId) -> Vec<ZoneId> {
        dispatch!(self, g => g.neighbors(id))
    }

    /// The zone's parents as DGGAL lists them, the primary parent first; see
    /// [`Grid::parents`].
    pub fn parents(&self, id: ZoneId) -> Vec<ZoneId> {
        dispatch!(self, g => g.parents(id))
    }

    /// The zone's children as DGGAL lists them, in its order; see [`Grid::children`].
    pub fn children(&self, id: ZoneId) -> Vec<ZoneId> {
        dispatch!(self, g => g.children(id))
    }

    /// The zone's centroid parent as DGGAL defines it, the first of its parents that is itself
    /// a centroid child; see [`Grid::centroid_parent`].
    pub fn centroid_parent(&self, id: ZoneId) -> Option<ZoneId> {
        dispatch!(self, g => g.centroid_parent(id))
    }

    /// Whether the zone is a centroid child; see [`Grid::is_centroid_child`].
    pub fn is_centroid_child(&self, id: ZoneId) -> bool {
        dispatch!(self, g => g.is_centroid_child(id))
    }

    /// How many sub-zones lie `depth` levels below the zone.
    pub fn count_sub_zones(&self, id: ZoneId, depth: u8) -> Result<u64> {
        dispatch!(self, g => g.count_sub_zones(id, depth))
    }

    /// The first sub-zone, in the grid's own order, `depth` levels below the
    /// zone.
    pub fn first_sub_zone(&self, id: ZoneId, depth: u8) -> Result<ZoneId> {
        dispatch!(self, g => g.first_sub_zone(id, depth))
    }

    /// Every sub-zone `depth` levels below the zone, in the grid's own order, which on
    /// every grid is the engine's, in scanlines across the zone; see [`Grid::sub_zones`],
    /// which also says where an entry on an aperture-7 grid is [`ZoneId::NULL`].
    pub fn sub_zones(&self, id: ZoneId, depth: u8) -> Result<Vec<ZoneId>> {
        dispatch!(self, g => g.sub_zones(id, depth))
    }

    /// Where `sub` sits in the zone's sub-zone order, or `None` if it is no
    /// sub-zone of it, or one that has no index. Both identifiers are read on this one grid, and each is
    /// validated first; see [`Grid::sub_zone_index`], which also gives the ways
    /// in which the answer departs from the engine's.
    pub fn sub_zone_index(&self, id: ZoneId, sub: ZoneId) -> Result<Option<u64>> {
        dispatch!(self, g => g.sub_zone_index(id, sub))
    }

    /// The sub-zone at `index`, `depth` levels below the zone, found without building
    /// the order; see [`Grid::sub_zone_at_index`].
    pub fn sub_zone_at_index(&self, id: ZoneId, depth: u8, index: u64) -> Result<ZoneId> {
        dispatch!(self, g => g.sub_zone_at_index(id, depth, index))
    }

    /// The identifier of the zone's primary parent, as [`crate::Zone::parent`]; `None` at
    /// the root.
    pub fn parent(&self, id: ZoneId) -> Option<ZoneId> {
        dispatch!(self, g => g.zone(id).parent().map(|z| z.id()))
    }

    /// Whether `a` is an ancestor of `b`, as [`crate::Zone::is_ancestor_of`].
    pub fn is_ancestor_of(&self, a: ZoneId, b: ZoneId) -> bool {
        dispatch!(self, g => g.zone(a).is_ancestor_of(&g.zone(b)))
    }

    /// Whether `a` and `b` are siblings, as [`crate::Zone::is_sibling_of`].
    pub fn is_sibling_of(&self, a: ZoneId, b: ZoneId) -> bool {
        dispatch!(self, g => g.zone(a).is_sibling_of(&g.zone(b)))
    }

    /// Whether `parent` is one of `child`'s parents.
    pub fn is_immediate_child_of(&self, child: ZoneId, parent: ZoneId) -> bool {
        dispatch!(self, g => g.zone(child).is_immediate_child_of(&g.zone(parent)))
    }

    /// The factor by which a level has more zones than the one above it; see
    /// [`Grid::refinement_ratio`].
    pub fn refinement_ratio(&self) -> u8 {
        dispatch!(self, g => g.refinement_ratio())
    }

    /// The most parents a zone has; see [`Grid::max_parents`].
    pub fn max_parents(&self) -> u8 {
        dispatch!(self, g => g.max_parents())
    }

    /// The most children a zone has; see [`Grid::max_children`].
    pub fn max_children(&self) -> u8 {
        dispatch!(self, g => g.max_children())
    }

    /// The most edge neighbours a zone has; see [`Grid::max_neighbors`].
    pub fn max_neighbors(&self) -> u8 {
        dispatch!(self, g => g.max_neighbors())
    }

    /// The depth at which a zone has about 65,536 sub-zones; see [`Grid::depth_64k`].
    pub fn depth_64k(&self) -> u8 {
        dispatch!(self, g => g.depth_64k())
    }

    /// The deepest relative depth at which DGGAL lists a zone's sub-zones; see
    /// [`Grid::max_depth`].
    pub fn max_depth(&self) -> u8 {
        dispatch!(self, g => g.max_depth())
    }

    /// The number of zones at `level`; see [`Grid::count_zones`], which also says why a level
    /// beyond the finest resolution is refused.
    pub fn count_zones(&self, level: u8) -> Result<u64> {
        dispatch!(self, g => g.count_zones(level))
    }

    /// The reference area of a zone at `level` in square metres; see [`Grid::ref_zone_area`].
    pub fn ref_zone_area(&self, level: u8) -> Result<f64> {
        dispatch!(self, g => g.ref_zone_area(level))
    }

    /// The size in metres of the sub-zones `depth` levels below a zone of `level`; see
    /// [`Grid::meters_per_sub_zone`].
    pub fn meters_per_sub_zone(&self, level: u8, depth: u8) -> Result<f64> {
        dispatch!(self, g => g.meters_per_sub_zone(level, depth))
    }

    /// The level whose reference area an area of `square_metres` reaches; see
    /// [`Grid::level_from_ref_zone_area`], which also says how far the answer reaches.
    pub fn level_from_ref_zone_area(&self, square_metres: f64) -> u8 {
        dispatch!(self, g => g.level_from_ref_zone_area(square_metres))
    }

    /// The level whose sub-zones `depth` levels below it measure `metres`; see
    /// [`Grid::level_from_meters_per_sub_zone`].
    pub fn level_from_meters_per_sub_zone(&self, metres: f64, depth: u8) -> u8 {
        dispatch!(self, g => g.level_from_meters_per_sub_zone(metres, depth))
    }

    /// Every zone within `k` neighbour steps of the zone, the zone included, as a
    /// lazy iterator of identifiers in the order [`Disk`] defines, the same walk as
    /// [`Grid::disk`] on the grid this handle names. Nothing is computed until the
    /// disk is consumed, so the caller controls the cost by how much it takes.
    ///
    /// The iterator is a concrete type rather than a boxed trait object: it holds
    /// this handle, which is `Copy` and `'static`, and asks the grid behind it for
    /// each zone's neighbours as the walk needs them. So it allocates nothing beyond
    /// the walk itself, it is `Clone`, `Debug`, `Send` and `Sync`, it offers
    /// [`Disk::rings`], and it keeps its type as further grids join the enum.
    pub fn disk(&self, id: ZoneId, k: u32) -> Disk<AnyGrid> {
        Disk::new(*self, id, k)
    }

    /// The finest level at which this crate enumerates zones, 33 on the aperture-3 grids and
    /// 14 on the aperture-7 grids; see [`Grid::max_box_level`], which says why.
    pub fn max_box_level(&self) -> u8 {
        dispatch!(self, g => g.max_box_level())
    }

    /// The zones that one cell of the lattice of `level` hosts; see `Grid::lattice_cell`.
    pub(crate) fn lattice_cell(
        &self,
        level: u8,
        root: u8,
        row: u64,
        col: u64,
    ) -> Vec<(ZoneId, u64)> {
        dispatch!(self, g => g.lattice_cell(level, root, row, col))
    }

    /// Where the lattice of `level` holds `zone`; see `Grid::lattice_place`.
    pub(crate) fn lattice_place(&self, level: u8, zone: ZoneId) -> Result<(u8, u64, u64, u64)> {
        dispatch!(self, g => g.lattice_place(level, zone))
    }

    /// Every zone of `level`, as a lazy iterator of identifiers in the order of DGGAL's
    /// `listZones` for the whole world, the same sequence as [`Grid::zones`] on the grid this
    /// handle names; see [`crate::Zones`] for the order and the cost, and
    /// [`crate::Zones::after`] to enter the sequence after one of its zones. It is
    /// [`Error::ResolutionOutOfRange`] beyond [`AnyGrid::max_box_level`].
    ///
    /// As with [`AnyGrid::disk`], the iterator is a concrete type that holds this handle, and
    /// keeps its type as further grids join the enum.
    pub fn zones(&self, level: u8) -> Result<crate::Zones<AnyGrid>> {
        let walk = dispatch!(self, g => g.lattice_walk(level))?;
        Ok(crate::Zones::new(*self, walk))
    }

    /// The point of the sphere at a point of the plane of the rhombi; see
    /// `Grid::lattice_point`.
    pub(crate) fn lattice_point(&self, x: f64, y: f64) -> Option<(f64, f64)> {
        dispatch!(self, g => g.lattice_point(x, y))
    }

    /// Whether the extent of `zone` meets the box; see `Grid::zone_meets_box`.
    pub(crate) fn zone_meets_box(&self, zone: ZoneId, bbox: &[f64; 4], reach: f64) -> bool {
        dispatch!(self, g => g.zone_meets_box(zone, bbox, reach))
    }

    /// The zones of `level` whose extent meets the bounding box `bbox`, as a lazy iterator of
    /// identifiers in the order of DGGAL's `listZones`, the same sequence as
    /// [`Grid::zones_in_box`] on the grid this handle names, with the same refusals; see that
    /// method for what a box may be and what a caller must know, [`crate::ZonesInBox`] for the
    /// rule, the order and the cost, and [`crate::ZonesInBox::after`] to enter the answer
    /// after a zone of the level.
    ///
    /// As with [`AnyGrid::disk`], the iterator is a concrete type that holds this handle, and
    /// keeps its type as further grids join the enum.
    pub fn zones_in_box(&self, level: u8, bbox: &Extent) -> Result<crate::ZonesInBox<AnyGrid>> {
        let (walk, search) = dispatch!(self, g => g.box_search(level, bbox))?;
        Ok(crate::ZonesInBox::new(*self, walk, search))
    }

    /// An upper bound, by arithmetic alone, on how many zones [`AnyGrid::zones_in_box`] yields
    /// for `level` and `bbox`, with the same refusals; see [`Grid::estimate_zones_in_box`] for
    /// what it is, what it rests on and how far above the count it lies.
    pub fn estimate_zones_in_box(&self, level: u8, bbox: &Extent) -> Result<u64> {
        dispatch!(self, g => g.estimate_zones_in_box(level, bbox))
    }

    /// DGGAL's `compactZones` of the zones of one level, on the grid this handle names: a
    /// shorter list that stands for the same set, in the engine's order; see
    /// [`Grid::compact_zones`] for what the answer is, that its zones overlap, and the
    /// refusals.
    pub fn compact_zones(&self, zones: &[ZoneId]) -> Result<Vec<ZoneId>> {
        dispatch!(self, g => g.compact_zones(zones))
    }
}

#[cfg(test)]
mod tests {
    use super::{GRID_NAMES, get_grid, igeo7, isea3h, ivea3h, ivea7h, rtea3h, rtea7h};
    use crate::{GridConfig, ZoneId};

    #[test]
    fn the_singleton_is_built_once() {
        assert!(std::ptr::eq(igeo7(), igeo7()));
        assert_eq!(igeo7().name(), "IGEO7");
        assert_eq!(*igeo7().config(), GridConfig::default());
        assert_eq!(igeo7().max_resolution(), 19);
    }

    /// `GRID_NAMES` and `get_grid` agree: every name the constant lists is one
    /// `get_grid` accepts, and answers with a grid whose own `name()` is that same
    /// canonical spelling.
    #[test]
    fn grid_names_matches_the_registry() {
        assert_eq!(
            GRID_NAMES,
            ["IGEO7", "IVEA7H", "RTEA7H", "ISEA3H", "IVEA3H", "RTEA3H"]
        );
        for name in GRID_NAMES {
            assert_eq!(get_grid(name).unwrap().name(), *name);
        }
    }

    /// `get_grid` matches a name without regard to case; a grid's own `name()`
    /// keeps the canonical form regardless of how the caller spelt it.
    #[test]
    fn get_grid_is_case_insensitive() {
        for name in ["igeo7", "IGEO7", "Igeo7", "iGeO7"] {
            assert_eq!(get_grid(name).unwrap().name(), "IGEO7");
        }
        assert_eq!(get_grid("ivea7h").unwrap().name(), "IVEA7H");
        assert_eq!(get_grid("RtEa7H").unwrap().name(), "RTEA7H");
    }

    /// The error names the grids this crate does offer, from `GRID_NAMES`, beside
    /// the name the caller actually gave.
    #[test]
    fn unknown_grid_message_lists_the_registered_names() {
        let err = get_grid("NOPE").unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("\"NOPE\""), "{msg}");
        for name in GRID_NAMES {
            assert!(msg.contains(name), "{msg}: missing {name}");
        }
        // A name of any length yields an error of bounded size, by its message and
        // by its `Debug` form alike.
        let err = get_grid(&"x".repeat(200_000)).unwrap_err();
        assert!(err.to_string().len() < 200, "{err}");
        assert!(format!("{err:?}").len() < 200, "{err:?}");
    }

    /// Each forwarding method added to give the handle the grid's full set of
    /// operations, called once, with `disk` also compared with the typed
    /// grid's own answer.
    #[test]
    fn any_grid_forwards_the_rest_of_the_grid() {
        let g = get_grid("IGEO7").unwrap();
        let typed = igeo7();
        let id = g.zone_from_text("0064156").unwrap();
        let pentagon = g.zone_from_text("0000").unwrap();
        let parent = g.zone_from_text("006415").unwrap();
        let sibling = g.zone_from_text("0064154").unwrap();
        let root = g.zone_from_text("00").unwrap();
        assert_eq!(g.max_resolution(), 19);
        assert_eq!(g.resolution(id), 5);
        assert!(!g.is_pentagon(id));
        assert!(g.is_pentagon(pentagon));
        // `0064156` is no centroid child, and neither of its parents, `006415` and `006414`,
        // is one, so that it has no centroid parent; its first child, `00641560`, is one, and
        // is therefore the centroid parent of its own child `006415601`.
        assert_eq!(g.centroid_parent(id), None);
        assert!(!g.is_centroid_child(id));
        let first_child = g.children(id)[0];
        assert_eq!(g.text_id(first_child), "00641560");
        assert!(g.is_centroid_child(first_child));
        assert_eq!(
            g.centroid_parent(g.zone_from_text("006415601").unwrap()),
            Some(first_child)
        );
        // The sub-zones of `0064156` at depth 1, in the engine's order, which is the typed
        // grid's.
        assert_eq!(g.count_sub_zones(id, 1), Ok(13));
        let order = g.sub_zones(id, 1).unwrap();
        assert_eq!(Ok(&order), typed.sub_zones(id, 1).as_ref());
        let texts: Vec<String> = order.iter().map(|&s| g.text_id(s)).collect();
        assert_eq!(
            texts,
            [
                "00641524", "00641506", "00641561", "00641563", "00641055", "00641565", "00641560",
                "00641562", "00641542", "00641564", "00641566", "00641431", "00641413",
            ]
        );
        assert_eq!(g.first_sub_zone(id, 1), Ok(order[0]));
        assert_eq!(g.sub_zone_at_index(id, 1, 0), Ok(order[0]));
        assert_eq!(g.sub_zone_at_index(id, 1, 12), Ok(order[12]));
        assert_eq!(g.sub_zone_index(id, id), Ok(Some(0)));
        // The first child is the centroid child, in the middle of the order.
        assert_eq!(g.sub_zone_index(id, g.children(id)[0]), Ok(Some(6)));
        assert_eq!(g.parent(id), Some(parent));
        assert_eq!(g.parent(id), typed.zone(id).parent().map(|z| z.id()));
        assert_eq!(g.parent(root), None);
        assert!(g.is_ancestor_of(parent, id) && !g.is_ancestor_of(id, parent));
        assert!(g.is_sibling_of(id, sibling) && g.is_sibling_of(sibling, id));
        assert!(!g.is_sibling_of(id, root));
        // A zone is not its own sibling, as DGGAL has it (`dggrs.ec:360`), pinned here
        // rather than inherited from the typed grid's own tests.
        assert!(!g.is_sibling_of(id, id));
        assert!(!g.is_sibling_of(root, root));
        assert!(g.is_immediate_child_of(id, parent) && !g.is_immediate_child_of(parent, id));
        let through_any: Vec<ZoneId> = g.disk(id, 2).collect();
        let typed_ids: Vec<ZoneId> = typed.zone(id).disk(2).map(|z| z.id()).collect();
        assert_eq!(through_any, typed_ids);
        assert_eq!(g.disk(id, 1).count(), 7);
        let rings: Vec<usize> = g.disk(id, 2).rings().map(|r| r.len()).collect();
        assert_eq!(rings, [1, 6, 12]);
    }

    /// What `AnyGrid::disk` promises of its iterator's type: nameable, cloneable,
    /// printable, and free to cross threads and to outlive the call that made it.
    #[test]
    fn any_grid_disk_is_a_concrete_static_iterator() {
        fn promised<T: Clone + std::fmt::Debug + Send + Sync + 'static>(_: &T) {}
        let g = get_grid("RTEA7H").unwrap();
        let disk: crate::Disk<super::AnyGrid> = g.disk(g.zone_from_text("0064156").unwrap(), 3);
        promised(&disk);
        promised(&disk.clone().rings());
        assert_eq!(disk.clone().count(), 37);
        assert_eq!(disk.rings().flatten().count(), 37);
    }

    /// The centroids are DGGAL v0.0.6's own, `getZoneWGS84Centroid` on `IVEA7H_Z7` and
    /// `RTEA7H_Z7` converted to degrees with `180 / PI`. They were taken from py4dggs
    /// until the projection was put into the order the compiled engine performs it, and
    /// under a tolerance of 1e-9 degrees; against the engine both longitudes had moved, by
    /// 3.6e-15 degrees each, two units in the last place, and the latitudes had not.
    #[test]
    fn ivea7h_and_rtea7h_anchors() {
        let z = ivea7h().zone_from_geo(38.7223, -9.1393, 5).unwrap();
        assert_eq!(z.text_id(), "0064156");
        let c = z.centroid();
        assert_eq!((c.lat, c.lon), (38.559613393079225, -8.943504906609359));
        assert_eq!(
            ivea7h()
                .zone_from_geo(38.7223, -9.1393, 15)
                .unwrap()
                .text_id(),
            "00641565005045433"
        );

        let z = rtea7h().zone_from_geo(38.7223, -9.1393, 5).unwrap();
        assert_eq!(z.text_id(), "0064156");
        let c = z.centroid();
        assert_eq!((c.lat, c.lon), (38.56406766826461, -8.989341241083546));
        assert_eq!(
            rtea7h()
                .zone_from_geo(38.7223, -9.1393, 15)
                .unwrap()
                .text_id(),
            "00641565034206566"
        );

        assert_eq!(get_grid("IVEA7H").unwrap().name(), "IVEA7H");
        assert_eq!(get_grid("RTEA7H").unwrap().name(), "RTEA7H");
    }

    #[test]
    fn any_grid_forwards_to_the_named_grid() {
        let g = get_grid("IGEO7").unwrap();
        let id = g.zone_from_geo(38.7223, -9.1393, 5).unwrap();
        assert_eq!(g.text_id(id), "0064156");
        assert_eq!(g.zone_from_text("0064156").unwrap(), id);
        assert_eq!(g.neighbors(id).len(), 6);
        assert_eq!(g.children(id).len(), 13);
        assert_eq!(
            g.parents(id),
            [
                g.zone_from_text("006415").unwrap(),
                g.zone_from_text("006414").unwrap()
            ]
        );
        assert_eq!(g.vertices(id).len(), 6);
        // DGGAL v0.0.6's own centroid latitude, as in `ivea7h_and_rtea7h_anchors`.
        assert_eq!(g.centroid(id).lat, 38.61687542349096);

        // Every aperture-7 grid, through the same handle, once: the relational operations
        // agree on the same texts across all three, since the hierarchy is found in the
        // plane the three share, and only the projection differs between them.
        for (name, centroid_lat) in [
            ("IGEO7", 38.61687542349096),
            ("IVEA7H", 38.559613393079225),
            ("RTEA7H", 38.56406766826461),
        ] {
            let g = get_grid(name).unwrap();
            let id = g.zone_from_geo(38.7223, -9.1393, 5).unwrap();
            let parent = g.zone_from_text("006415").unwrap();
            let sibling = g.zone_from_text("0064154").unwrap();
            assert_eq!(g.text_id(id), "0064156");
            assert_eq!(g.centroid(id).lat, centroid_lat, "{name}");
            assert_eq!(g.parent(id), Some(parent));
            assert_eq!(g.parent(g.zone_from_text("00").unwrap()), None);
            assert!(g.is_ancestor_of(parent, id) && !g.is_ancestor_of(id, parent));
            assert!(g.is_sibling_of(id, sibling) && g.is_sibling_of(sibling, id));
            assert!(
                !g.is_sibling_of(id, id),
                "{name}: a zone is its own sibling"
            );
            assert!(g.is_immediate_child_of(id, parent) && !g.is_immediate_child_of(parent, id));
        }
    }

    /// ISEA3H through the registry: its singleton, its name in any casing, its
    /// finest resolution, and the zone containing Lisbon at resolution 5 with its
    /// relatives, as DGGAL v0.0.6 answers them (`C2-23-C`, parents `C2-23-A`,
    /// `C4-5-A` and `C4-6-A`, six neighbours and seven children).
    #[test]
    fn isea3h_through_the_registry() {
        assert!(std::ptr::eq(isea3h(), isea3h()));
        assert_eq!(isea3h().name(), "ISEA3H");
        assert_eq!(isea3h().max_resolution(), 33);
        let g = get_grid("isea3h").unwrap();
        assert_eq!(g.name(), "ISEA3H");
        assert_eq!(g.max_resolution(), 33);
        let id = g.zone_from_geo(38.7223, -9.1393, 5).unwrap();
        assert_eq!(g.text_id(id), "C2-23-C");
        assert_eq!(id, ZoneId(0x0440_0000_0000_008E));
        let parents: Vec<String> = g.parents(id).into_iter().map(|p| g.text_id(p)).collect();
        assert_eq!(parents, ["C2-23-A", "C4-5-A", "C4-6-A"]);
        assert_eq!(g.neighbors(id).len(), 6);
        assert_eq!(g.children(id).len(), 7);
        assert_eq!(g.vertices(id).len(), 6);
    }

    /// IVEA3H and RTEA3H through the registry: their singletons, their names in any
    /// casing, their finest resolution, and the zone containing Lisbon, as DGGAL v0.0.6
    /// answers it on each. At resolution 5 it is ISEA3H's `C2-23-C`, with the same parents,
    /// the three grids sharing their topology and indexing; at 15 the three projections
    /// place the point in three different zones: `H2-19D3E7-B` on ISEA3H, `H2-19C2D2-B` on
    /// IVEA3H and `H2-19C2D3-D` on RTEA3H.
    #[test]
    fn ivea3h_and_rtea3h_through_the_registry() {
        assert!(std::ptr::eq(ivea3h(), ivea3h()));
        assert!(std::ptr::eq(rtea3h(), rtea3h()));
        assert_eq!(ivea3h().name(), "IVEA3H");
        assert_eq!(rtea3h().name(), "RTEA3H");
        assert_eq!(ivea3h().max_resolution(), 33);
        assert_eq!(rtea3h().max_resolution(), 33);
        for (spelt, name, at_15) in [
            ("isea3h", "ISEA3H", "H2-19D3E7-B"),
            ("Ivea3h", "IVEA3H", "H2-19C2D2-B"),
            ("rTeA3H", "RTEA3H", "H2-19C2D3-D"),
        ] {
            let g = get_grid(spelt).unwrap();
            assert_eq!(g.name(), name);
            assert_eq!(g.max_resolution(), 33, "{name}");
            let id = g.zone_from_geo(38.7223, -9.1393, 5).unwrap();
            assert_eq!(g.text_id(id), "C2-23-C", "{name}");
            assert_eq!(id, ZoneId(0x0440_0000_0000_008E), "{name}");
            let parents: Vec<String> = g.parents(id).into_iter().map(|p| g.text_id(p)).collect();
            assert_eq!(parents, ["C2-23-A", "C4-5-A", "C4-6-A"], "{name}");
            assert_eq!(g.neighbors(id).len(), 6, "{name}");
            assert_eq!(g.children(id).len(), 7, "{name}");
            let id = g.zone_from_geo(38.7223, -9.1393, 15).unwrap();
            assert_eq!(g.text_id(id), at_15, "{name}");
        }
    }
}
