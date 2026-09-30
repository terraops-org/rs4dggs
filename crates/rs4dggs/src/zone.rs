//! A `Zone` is a handle: a grid reference and an identifier, with no logic of
//! its own. Every question it answers it forwards to its [`Grid`], so that
//! there is one implementation of each operation and a zone is free to copy
//! and cheap to keep. Ported from py4dggs's `zone.py`.
//!
//! Source: py4dggs's `zone.py`.

use std::collections::HashSet;

use crate::disk::Disk;
use crate::grid::Grid;
use crate::{GeoPoint, Indexing, Projection, Result, Topology, ZoneId};

/// One cell of one grid.
///
/// A zone carries its grid, so it can answer any question about itself, and
/// two zones are the same zone only when they belong to the same grid
/// instance. Zones of different grid types cannot be compared at all, since
/// they are different types; but nothing stops a caller from building a second
/// grid of the same type, with [`Grid::new`], perhaps on another orientation,
/// and there the same identifier names a different cell, so identifiers alone
/// cannot tell a cell of one instance from a cell of the other.
pub struct Zone<'g, G> {
    grid: &'g G,
    id: ZoneId,
}

impl<G> Clone for Zone<'_, G> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<G> Copy for Zone<'_, G> {}
impl<G> PartialEq for Zone<'_, G> {
    fn eq(&self, o: &Self) -> bool {
        std::ptr::eq(self.grid, o.grid) && self.id == o.id
    }
}
impl<G> Eq for Zone<'_, G> {}
impl<G> std::hash::Hash for Zone<'_, G> {
    fn hash<H: std::hash::Hasher>(&self, h: &mut H) {
        self.id.hash(h)
    }
}

impl<P: Projection, T: Topology, I: Indexing> std::fmt::Debug for Zone<'_, Grid<P, T, I>> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Zone({:?})", self.text_id())
    }
}

impl<'g, P: Projection, T: Topology, I: Indexing> Zone<'g, Grid<P, T, I>> {
    pub(crate) fn new(grid: &'g Grid<P, T, I>, id: ZoneId) -> Self {
        Zone { grid, id }
    }

    fn wrap(&self, ids: Vec<ZoneId>) -> Vec<Self> {
        ids.into_iter().map(|id| self.grid.zone(id)).collect()
    }

    /// The zone's packed identifier.
    pub fn id(&self) -> ZoneId {
        self.id
    }

    /// The grid this zone belongs to.
    pub fn grid(&self) -> &'g Grid<P, T, I> {
        self.grid
    }

    /// The zone's canonical text id.
    pub fn text_id(&self) -> String {
        self.grid.text_id(self.id)
    }

    /// The zone's resolution.
    pub fn resolution(&self) -> u8 {
        self.grid.resolution(self.id)
    }

    /// Whether the zone is a pentagon, that is an icosahedron vertex cell.
    pub fn is_pentagon(&self) -> bool {
        self.grid.is_pentagon(self.id)
    }

    /// The zone's centroid in WGS84 degrees.
    pub fn centroid(&self) -> GeoPoint {
        self.grid.centroid(self.id)
    }

    /// The zone's boundary vertices in WGS84 degrees, anticlockwise as seen from outside
    /// the sphere on every grid; see [`Grid::vertices`] for the order, the antimeridian
    /// and the vertices that can coincide.
    ///
    /// # Examples
    ///
    /// ```
    /// let zone = rs4dggs::igeo7().zone_from_text("0064156")?;
    /// let ring = zone.vertices();
    /// assert_eq!(ring.len(), 6); // five for a pentagon
    ///
    /// // The ring is not closed: formats such as GeoJSON need the first vertex repeated.
    /// let mut closed = ring.clone();
    /// closed.push(ring[0]);
    /// assert_eq!(closed.len(), 7);
    /// # Ok::<(), rs4dggs::Error>(())
    /// ```
    pub fn vertices(&self) -> Vec<GeoPoint> {
        self.grid.vertices(self.id)
    }

    /// The zone's edge neighbours; see [`Grid::neighbors`], which also says where the
    /// relation is not symmetric.
    ///
    /// # Examples
    ///
    /// ```
    /// let grid = rs4dggs::igeo7();
    /// assert_eq!(grid.zone_from_text("0064156")?.neighbors().len(), 6);
    ///
    /// // Every base cell is a pentagon at resolution 0, with five neighbours.
    /// let base = grid.zone_from_text("00")?;
    /// assert!(base.is_pentagon());
    /// assert_eq!(base.neighbors().len(), 5);
    /// # Ok::<(), rs4dggs::Error>(())
    /// ```
    pub fn neighbors(&self) -> Vec<Self> {
        self.wrap(self.grid.neighbors(self.id))
    }

    /// The primary parent, that is the first of `parents`, or `None` at the
    /// root. On a congruent grid this is the single parent; on an aperture-3
    /// grid it is DGGAL's own `parent0`.
    ///
    /// # Examples
    ///
    /// ```
    /// let zone = rs4dggs::igeo7().zone_from_text("0064156")?;
    /// let parent = zone.parent().expect("a zone above resolution 0 has a parent");
    /// assert_eq!(parent.text_id(), "006415");
    /// assert!(parent.is_ancestor_of(&zone));
    ///
    /// // A base cell has none.
    /// assert!(rs4dggs::igeo7().zone_from_text("00")?.parent().is_none());
    /// # Ok::<(), rs4dggs::Error>(())
    /// ```
    pub fn parent(&self) -> Option<Self> {
        self.grid
            .parents(self.id)
            .first()
            .map(|&p| self.grid.zone(p))
    }

    /// Every parent: one on a congruent grid, one or three on a non-congruent
    /// aperture-3 one.
    pub fn parents(&self) -> Vec<Self> {
        self.wrap(self.grid.parents(self.id))
    }

    /// The zone's children, empty at the maximum resolution.
    ///
    /// # Examples
    ///
    /// ```
    /// let grid = rs4dggs::igeo7();
    /// let zone = grid.zone_from_text("0064156")?;
    /// let children = zone.children();
    /// assert_eq!(children.len(), 7);
    /// assert!(children.iter().all(|child| child.is_immediate_child_of(&zone)));
    ///
    /// // A pentagon has six.
    /// assert_eq!(grid.zone_from_text("00")?.children().len(), 6);
    /// # Ok::<(), rs4dggs::Error>(())
    /// ```
    pub fn children(&self) -> Vec<Self> {
        self.wrap(self.grid.children(self.id))
    }

    /// The zone's centroid parent as DGGAL defines it: the first of its parents
    /// that is itself a centroid child, or `None` where none is; the single
    /// parent on a congruent grid. On an aperture-3 grid this is not the parent
    /// of which the zone is the centroid child; see [`Grid::centroid_parent`].
    pub fn centroid_parent(&self) -> Option<Self> {
        self.grid
            .centroid_parent(self.id)
            .map(|p| self.grid.zone(p))
    }

    /// Whether the zone is its parent's centroid child, which every zone of a
    /// congruent grid is.
    pub fn is_centroid_child(&self) -> bool {
        self.grid.is_centroid_child(self.id)
    }

    /// Whether `other` is one of this zone's parents.
    pub fn is_immediate_child_of(&self, other: &Self) -> bool {
        self.parents().contains(other)
    }

    /// Every sub-zone `depth` levels below this one, in the grid's own order.
    pub fn sub_zones(&self, depth: u8) -> Result<Vec<Self>> {
        Ok(self.wrap(self.grid.sub_zones(self.id, depth)?))
    }

    /// How many sub-zones lie `depth` levels below this one.
    pub fn count_sub_zones(&self, depth: u8) -> Result<u64> {
        self.grid.count_sub_zones(self.id, depth)
    }

    /// The first sub-zone `depth` levels below this one.
    pub fn first_sub_zone(&self, depth: u8) -> Result<Self> {
        Ok(self.grid.zone(self.grid.first_sub_zone(self.id, depth)?))
    }

    /// Where `sub` sits in this zone's sub-zone order, or `None` if it is no
    /// sub-zone of it; see [`Grid::sub_zone_index`] for the ways in which the
    /// answer departs from the engine's.
    ///
    /// The grid check is not redundant: two instances of the same grid type
    /// share one packing, so comparing identifiers alone cannot tell a genuine
    /// sub-zone from the same number on the other instance, and a foreign zone
    /// would be handed a plausible index. Equality already compares grids, and
    /// this keeps the two consistent.
    pub fn sub_zone_index(&self, sub: &Self) -> Result<Option<u64>> {
        if !std::ptr::eq(self.grid, sub.grid) {
            return Ok(None);
        }
        self.grid.sub_zone_index(self.id, sub.id)
    }

    /// The sub-zone at `index`, `depth` levels below this one.
    pub fn sub_zone_at_index(&self, depth: u8, index: u64) -> Result<Self> {
        Ok(self
            .grid
            .zone(self.grid.sub_zone_at_index(self.id, depth, index)?))
    }

    /// Whether this zone is an ancestor of `other`.
    ///
    /// The walk goes up the parent graph rather than down the children,
    /// because a non-congruent hierarchy is a directed acyclic graph and not a
    /// tree: a zone may have three parents, whose own ancestors overlap, so
    /// the frontier is kept as a set and every zone already visited is
    /// dropped.
    ///
    /// A zone of another grid instance is never a descendant, for the reason
    /// given at [`Zone::sub_zone_index`], which the same check guards.
    pub fn is_ancestor_of(&self, other: &Self) -> bool {
        if !std::ptr::eq(self.grid, other.grid) {
            return false;
        }
        let mut seen: HashSet<ZoneId> = HashSet::new();
        let mut frontier: HashSet<ZoneId> = HashSet::from([other.id]);
        while !frontier.is_empty() {
            let mut next = HashSet::new();
            for z in frontier {
                for p in self.grid.parents(z) {
                    if p == self.id {
                        return true;
                    }
                    if seen.insert(p) {
                        next.insert(p);
                    }
                }
            }
            frontier = next;
        }
        false
    }

    /// Whether this zone and `other` share a parent, and with it a
    /// resolution. A zone of another grid instance is never a sibling, for the
    /// reason given at [`Zone::sub_zone_index`].
    ///
    /// A zone is not its own sibling. DGGAL's `areZonesSiblings` excludes it
    /// (`dggrs.ec:363`); py4dggs does not, since it only tests for a shared
    /// parent, which a zone trivially shares with itself. This crate follows
    /// DGGAL.
    pub fn is_sibling_of(&self, other: &Self) -> bool {
        if !std::ptr::eq(self.grid, other.grid) {
            return false;
        }
        if self.id == other.id {
            return false;
        }
        if self.resolution() != other.resolution() {
            return false;
        }
        let mine: HashSet<ZoneId> = self.grid.parents(self.id).into_iter().collect();
        self.grid.parents(other.id).iter().any(|p| mine.contains(p))
    }

    /// Every zone within `k` neighbour steps, this one included, as a lazy
    /// iterator: this zone first, then ring after ring outward, in the order
    /// [`Disk`] defines, which is the same on every call.
    ///
    /// The disk is composed from `neighbors` step by step rather than from any
    /// closed form, so it inherits whatever the grid's own neighbour list is,
    /// the engine's own on every grid here, and it stops growing where the grid
    /// does: once a ring comes back empty, the walk ends, however large `k` is. Nothing is
    /// computed until the disk is consumed, so the caller controls the cost by
    /// how much it takes; [`Disk::rings`] gives the same walk a ring at a time.
    ///
    /// ```
    /// let zone = rs4dggs::igeo7().zone_from_text("0064156")?;
    /// let within_two: Vec<_> = zone.disk(2).collect();
    /// assert_eq!(within_two.len(), 19);
    /// assert_eq!(within_two[0], zone);
    /// let first_ten: Vec<_> = zone.disk(u32::MAX).take(10).collect();
    /// assert_eq!(first_ten[..], within_two[..10]);
    /// # Ok::<(), rs4dggs::Error>(())
    /// ```
    pub fn disk(&self, k: u32) -> Disk<Self> {
        Disk::new(*self, self.id, k)
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use crate::registry::igeo7;

    #[test]
    fn equality_and_debug_read_from_the_identifier() {
        let g = igeo7();
        let a = g.zone_from_text("0064156").unwrap();
        let b = g.zone_from_text("0064156").unwrap();
        assert_eq!(a, b);
        assert_ne!(a, g.zone_from_text("0064154").unwrap());
        assert_eq!(format!("{a:?}"), "Zone(\"0064156\")");
    }

    /// A second instance of the same grid type may name different cells by the
    /// same identifiers, since it may be built on another orientation, so its
    /// zones are neither equal to, nor ancestors, siblings or sub-zones of, the
    /// first instance's. That holds even when, as here, the two instances share
    /// an orientation.
    #[test]
    fn zones_of_another_instance_are_not_related() {
        let other = crate::registry::Igeo7::new(crate::GridConfig::default(), "IGEO7").unwrap();
        let g = igeo7();
        let parent = g.zone_from_text("006415").unwrap();
        let child = g.zone_from_text("0064156").unwrap();
        let sibling = g.zone_from_text("0064154").unwrap();
        assert!(parent.is_ancestor_of(&child));
        assert!(child.is_sibling_of(&sibling));
        let foreign_parent = other.zone(parent.id());
        let foreign_child = other.zone(child.id());
        let foreign_sibling = other.zone(sibling.id());
        assert_ne!(parent, foreign_parent);
        assert!(!parent.is_ancestor_of(&foreign_child));
        assert!(!foreign_parent.is_ancestor_of(&child));
        assert!(!child.is_sibling_of(&foreign_sibling));
        assert!(!foreign_child.is_sibling_of(&sibling));
        assert_eq!(child.sub_zone_index(&foreign_child), Ok(None));
    }

    /// The walk ends once a ring comes back empty, so a vast `k` costs no more
    /// than the grid it covers. Resolution 0 has 12 zones, all within three
    /// steps of any one, and each is yielded once.
    #[test]
    fn disk_stops_when_the_grid_is_covered() {
        let z = igeo7().zone_from_text("00").unwrap();
        let all: Vec<_> = z.disk(u32::MAX).collect();
        assert_eq!(all.len(), 12);
        let distinct: HashSet<_> = all.iter().map(|z| z.id()).collect();
        assert_eq!(distinct.len(), 12);
    }

    #[test]
    fn siblings_share_a_parent() {
        let g = igeo7();
        let a = g.zone_from_text("0064156").unwrap();
        let b = g.zone_from_text("0064154").unwrap();
        assert!(a.is_sibling_of(&b));
        assert!(!a.is_sibling_of(&g.zone_from_text("0064141").unwrap()));
        assert!(!a.is_sibling_of(&a.parent().unwrap()));
    }

    /// A zone is not its own sibling, following DGGAL's `areZonesSiblings`
    /// (`dggrs.ec:363`), which excludes it; py4dggs answers true here, since it
    /// only tests for a shared parent.
    #[test]
    fn a_zone_is_not_its_own_sibling() {
        let g = igeo7();
        let a = g.zone_from_text("0064156").unwrap();
        assert!(!a.is_sibling_of(&a));
        let b = g.zone_from_text("0064156").unwrap();
        assert!(!a.is_sibling_of(&b));
    }
}
