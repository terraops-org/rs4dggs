//! The disk of radius `k` about a zone, walked lazily and outward, ring by ring.
//!
//! DGGAL has no disk of its own: it is composed here from the grid's neighbour
//! relation, one step at a time, so that it inherits whatever that relation is and
//! stops growing where the grid does. The walk is an iterator, so that a caller pays
//! for the zones it takes and for no others: a radius arriving from outside, however
//! large, costs nothing until the disk is consumed, and costs only as much as is
//! consumed. See [`Disk`] for the order it yields and for what it holds in memory.
//!
//! Source: this crate's own API.

use std::collections::HashSet;
use std::iter::FusedIterator;

use crate::grid::Grid;
use crate::registry::AnyGrid;
use crate::zone::Zone;
use crate::{Indexing, Projection, Topology, ZoneId};

/// Every zone within `k` neighbour steps of a centre, the centre included, as a lazy
/// iterator.
///
/// A disk is taken from a [`Zone`] with [`Zone::disk`], and then yields zones; or
/// from a [`Grid`] with [`Grid::disk`], or from an [`AnyGrid`] with
/// [`AnyGrid::disk`], and then yields identifiers. The type parameter is the handle
/// the disk was taken from, and the three `Iterator` implementations below state
/// what each yields. [`Disk::rings`] turns the same walk into one of whole rings.
///
/// # Order
///
/// The order is defined by this crate, since the engine has no disk and so no order
/// for one, and it is the same on every call, in every run and on every platform: no
/// hashed collection decides it, sets serving only to recognise a zone already
/// yielded. The centre comes first, as ring 0. Ring `n + 1` follows ring `n`, and it
/// is built from ring `n`'s zones taken in the order in which they were themselves
/// yielded: for each of them in turn, its neighbours in the order in which
/// [`Grid::neighbors`] returns them, each yielded on first arrival and passed over on
/// any later one. So a zone belongs to the ring of the fewest neighbour steps that
/// reach it, and within a ring the zones stand in the order in which the walk first
/// came upon them.
///
/// # Cost and termination
///
/// Nothing is computed until the disk is consumed. Ring `n + 1` costs one neighbour
/// query for each zone of ring `n`, and each query is made only when the walk needs
/// the zones it may contribute, so a caller controls the cost by how much it takes:
/// `disk(u32::MAX)` is safe to construct, and taking the centre and its neighbours
/// from it asks for the neighbours of the centre alone. The walk stops after ring
/// `k`, without asking for the neighbours of ring `k`, and earlier when a ring comes
/// back empty, which happens once the grid is exhausted: at resolution 0, twelve
/// zones in all, within three steps of any one.
///
/// Once the walk has ended, it stays ended: every later call of `next` returns `None`
/// and asks for nothing, so `Disk` and [`Rings`] are each a [`FusedIterator`].
///
/// The centre is yielded whatever it is. An identifier with no zone behind it, the
/// null zone among them, has no neighbours (see [`Grid::neighbors`]), and its disk is
/// that identifier alone.
///
/// # Memory
///
/// The walk remembers every zone it has yielded, together with the ring from which
/// the ring in progress is being built, so what it holds grows with what the caller
/// takes, and a caller that consumes a disk of radius `k` whole holds of the order of
/// `3k²` identifiers by the end.
///
/// It cannot keep only the last two rings, which would bound it by the order of `k`.
/// That economy is sound only for a symmetric neighbour relation, in which a zone's
/// neighbours lie in the ring before its own, its own ring or the ring after it, and
/// this crate's relation is not symmetric everywhere. On the aperture-7 grids, at
/// resolutions 15 to 19, within some 1.3e-5 degrees of the two icosahedron edges where
/// the engine's odd resolutions leave a strip with no cells (from base cell 0 to 1 over
/// the north pole, and from 1 to 6), a zone may list a neighbour that does not list it
/// back, and there a walk that remembers two rings yields some zones twice. Over some
/// 950,000 zones of each of those three grids, every zone of resolutions 0 to 4, the
/// pentagons and both poles, the two broken edges and every other icosahedron edge, and
/// the two-step neighbourhoods of those zones at every resolution among them, no such
/// pair was found further out. Among the 2,860 digit paths deep under the two polar
/// pentagons, zeros followed by at most three other digits, the reach is greater, every
/// such pair lying under the pentagon of base cell 0: up to 3.57e-5 degrees, over 1,673
/// one-way pairs on each grid.
#[must_use = "a disk is walked lazily, and does nothing unless it is consumed"]
#[derive(Debug, Clone)]
pub struct Disk<S> {
    source: S,
    walk: Walk,
}

/// The same walk as a [`Disk`], yielding one ring at a time rather than one zone:
/// ring 0, the centre alone, then ring 1, and so on, each a `Vec` in the order the
/// disk would have yielded its zones. The ring number is the position in the
/// sequence, so `.enumerate()` pairs each ring with it.
///
/// A ring is yielded once it is complete, which costs the neighbour queries of every
/// zone of the ring inside it. No ring is ever empty: the walk ends instead.
#[must_use = "a disk is walked lazily, and does nothing unless it is consumed"]
#[derive(Debug, Clone)]
pub struct Rings<S> {
    source: S,
    walk: Walk,
}

impl<S> Disk<S> {
    pub(crate) fn new(source: S, centre: ZoneId, k: u32) -> Self {
        Disk {
            source,
            walk: Walk::new(centre, k),
        }
    }

    /// The same walk, one ring at a time. On a disk not yet consumed the first item
    /// is ring 0; on one partly consumed it is whatever remains of the ring in
    /// progress, and the rings that follow are whole.
    ///
    /// # Examples
    ///
    /// ```
    /// let zone = rs4dggs::igeo7().zone_from_text("0064156")?;
    /// let sizes: Vec<usize> = zone.disk(2).rings().map(|ring| ring.len()).collect();
    /// assert_eq!(sizes, [1, 6, 12]);
    /// # Ok::<(), rs4dggs::Error>(())
    /// ```
    pub fn rings(self) -> Rings<S> {
        Rings {
            source: self.source,
            walk: self.walk,
        }
    }
}

impl<'g, P: Projection, T: Topology, I: Indexing> Iterator for Disk<Zone<'g, Grid<P, T, I>>> {
    type Item = Zone<'g, Grid<P, T, I>>;

    fn next(&mut self) -> Option<Self::Item> {
        let g = self.source.grid();
        self.walk.next_zone(|z| g.neighbors(z)).map(|z| g.zone(z))
    }
}

impl<P: Projection, T: Topology, I: Indexing> Iterator for Disk<&Grid<P, T, I>> {
    type Item = ZoneId;

    fn next(&mut self) -> Option<ZoneId> {
        let g = self.source;
        self.walk.next_zone(|z| g.neighbors(z))
    }
}

impl Iterator for Disk<AnyGrid> {
    type Item = ZoneId;

    fn next(&mut self) -> Option<ZoneId> {
        let g = self.source;
        self.walk.next_zone(|z| g.neighbors(z))
    }
}

// Fused, because the walk's ended state is a sink: `next_zone` and `next_ring` return
// `None` only once the walk has ended, and on an ended walk they return `None` at once.
impl<'g, P: Projection, T: Topology, I: Indexing> FusedIterator for Disk<Zone<'g, Grid<P, T, I>>> {}

impl<P: Projection, T: Topology, I: Indexing> FusedIterator for Disk<&Grid<P, T, I>> {}

impl FusedIterator for Disk<AnyGrid> {}

impl<'g, P: Projection, T: Topology, I: Indexing> Iterator for Rings<Zone<'g, Grid<P, T, I>>> {
    type Item = Vec<Zone<'g, Grid<P, T, I>>>;

    fn next(&mut self) -> Option<Self::Item> {
        let g = self.source.grid();
        self.walk
            .next_ring(|z| g.neighbors(z))
            .map(|ring| ring.into_iter().map(|z| g.zone(z)).collect())
    }
}

impl<P: Projection, T: Topology, I: Indexing> Iterator for Rings<&Grid<P, T, I>> {
    type Item = Vec<ZoneId>;

    fn next(&mut self) -> Option<Vec<ZoneId>> {
        let g = self.source;
        self.walk.next_ring(|z| g.neighbors(z))
    }
}

impl Iterator for Rings<AnyGrid> {
    type Item = Vec<ZoneId>;

    fn next(&mut self) -> Option<Vec<ZoneId>> {
        let g = self.source;
        self.walk.next_ring(|z| g.neighbors(z))
    }
}

impl<'g, P: Projection, T: Topology, I: Indexing> FusedIterator for Rings<Zone<'g, Grid<P, T, I>>> {}

impl<P: Projection, T: Topology, I: Indexing> FusedIterator for Rings<&Grid<P, T, I>> {}

impl FusedIterator for Rings<AnyGrid> {}

/// Where a walk stands.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum State {
    /// Nothing yielded yet: the centre waits in `outer`.
    Fresh,
    /// The centre yielded, and ring `ring` in progress.
    Walking,
    /// The walk ended, at ring `k` or at an empty ring, and its memory released.
    Ended,
}

/// The walk itself, over identifiers, with the neighbour relation supplied at each
/// step, so that the three handles share one implementation and differ only in how
/// they ask for neighbours and in what they yield.
#[derive(Debug, Clone)]
struct Walk {
    k: u32,
    state: State,
    /// The ring in progress, whose zones are being yielded.
    ring: u32,
    /// Every zone yielded so far. It serves to recognise a zone already yielded and
    /// never decides the order.
    seen: HashSet<ZoneId>,
    /// Ring `ring - 1`, whose zones' neighbours make up the ring in progress, and the
    /// position of the next of them whose neighbours are to be asked for.
    inner: Vec<ZoneId>,
    next_inner: usize,
    /// The neighbours of the zone of `inner` asked last, and the position of the next
    /// of them to be considered.
    pending: Vec<ZoneId>,
    next_pending: usize,
    /// The ring in progress, as far as it has been yielded.
    outer: Vec<ZoneId>,
}

impl Walk {
    fn new(centre: ZoneId, k: u32) -> Self {
        Walk {
            k,
            state: State::Fresh,
            ring: 0,
            seen: HashSet::new(),
            inner: Vec::new(),
            next_inner: 0,
            pending: Vec::new(),
            next_pending: 0,
            outer: vec![centre],
        }
    }

    /// Yields the centre, which is ring 0 whole.
    fn start(&mut self) -> ZoneId {
        let centre = self.outer[0];
        self.seen.insert(centre);
        self.state = State::Walking;
        centre
    }

    /// The next zone of the ring in progress, or `None` once that ring is complete.
    fn next_in_ring(
        &mut self,
        neighbors: &mut impl FnMut(ZoneId) -> Vec<ZoneId>,
    ) -> Option<ZoneId> {
        loop {
            while let Some(&w) = self.pending.get(self.next_pending) {
                self.next_pending += 1;
                if self.seen.insert(w) {
                    self.outer.push(w);
                    return Some(w);
                }
            }
            let &z = self.inner.get(self.next_inner)?;
            self.next_inner += 1;
            self.pending = neighbors(z);
            self.next_pending = 0;
        }
    }

    /// Moves on from a complete ring to the next, or ends the walk, when the ring just
    /// completed was ring `k` or came back empty. Returns whether the walk goes on.
    /// The comparison with `k` comes before the increment, so that `k = u32::MAX`
    /// cannot overflow the ring number.
    fn advance(&mut self) -> bool {
        if self.outer.is_empty() || self.ring == self.k {
            self.state = State::Ended;
            self.seen = HashSet::new();
            self.inner = Vec::new();
            self.pending = Vec::new();
            self.outer = Vec::new();
            return false;
        }
        self.inner = std::mem::take(&mut self.outer);
        self.next_inner = 0;
        self.ring += 1;
        true
    }

    /// The next zone of the walk.
    fn next_zone(&mut self, mut neighbors: impl FnMut(ZoneId) -> Vec<ZoneId>) -> Option<ZoneId> {
        match self.state {
            State::Fresh => return Some(self.start()),
            State::Ended => return None,
            State::Walking => {}
        }
        loop {
            if let Some(z) = self.next_in_ring(&mut neighbors) {
                return Some(z);
            }
            if !self.advance() {
                return None;
            }
        }
    }

    /// What remains of the ring in progress, or the next ring whole if nothing does.
    fn next_ring(
        &mut self,
        mut neighbors: impl FnMut(ZoneId) -> Vec<ZoneId>,
    ) -> Option<Vec<ZoneId>> {
        match self.state {
            State::Fresh => return Some(vec![self.start()]),
            State::Ended => return None,
            State::Walking => {}
        }
        loop {
            let from = self.outer.len();
            while self.next_in_ring(&mut neighbors).is_some() {}
            if self.outer.len() > from {
                return Some(self.outer[from..].to_vec());
            }
            if !self.advance() {
                return None;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use std::cell::Cell;

    use super::Walk;
    use crate::ZoneId;
    use crate::registry::{get_grid, igeo7};

    /// The hexagonal lattice in axial coordinates, packed into an identifier: a plane
    /// without end, in which ring `n` has exactly `6n` cells, and whose neighbour
    /// order is fixed by the list below.
    fn hex(q: i32, r: i32) -> ZoneId {
        ZoneId((u64::from(q as u32) << 32) | u64::from(r as u32))
    }
    fn hex_neighbours(z: ZoneId) -> Vec<ZoneId> {
        let (q, r) = ((z.0 >> 32) as u32 as i32, z.0 as u32 as i32);
        [(1, 0), (1, -1), (0, -1), (-1, 0), (-1, 1), (0, 1)]
            .iter()
            .map(|&(dq, dr)| hex(q + dq, r + dr))
            .collect()
    }

    fn zones(walk: &mut Walk, f: impl Fn(ZoneId) -> Vec<ZoneId>) -> Vec<ZoneId> {
        std::iter::from_fn(|| walk.next_zone(&f)).collect()
    }

    /// Once a walk has ended, by reaching ring `k` or by coming back empty, every later
    /// call answers `None` and asks for no neighbours, whether zones or rings are
    /// taken; and every public shape of the walk declares it.
    #[test]
    fn an_ended_walk_stays_ended() {
        use std::iter::FusedIterator;

        let calls = Cell::new(0usize);
        let counted = |z| {
            calls.set(calls.get() + 1);
            hex_neighbours(z)
        };
        let mut walk = Walk::new(hex(0, 0), 2);
        assert_eq!(zones(&mut walk, counted).len(), 19);
        let spent = calls.get();
        for _ in 0..5 {
            assert_eq!(walk.next_zone(counted), None);
            assert_eq!(walk.next_ring(counted), None);
        }
        assert_eq!(calls.get(), spent, "an ended walk asked for neighbours");

        let mut walk = Walk::new(hex(0, 0), 2);
        assert_eq!(std::iter::from_fn(|| walk.next_ring(counted)).count(), 3);
        for _ in 0..5 {
            assert_eq!(walk.next_ring(counted), None);
            assert_eq!(walk.next_zone(counted), None);
        }

        // Ended by an empty ring rather than by `k`: resolution 0 has twelve zones.
        let g = igeo7();
        let base = g.zone_from_text("00").unwrap();
        let mut disk = base.disk(u32::MAX);
        assert_eq!(disk.by_ref().count(), 12);
        assert!((0..5).all(|_| disk.next().is_none()));
        let mut rings = g.disk(base.id(), u32::MAX).rings();
        assert_eq!(rings.by_ref().count(), 4);
        assert!((0..5).all(|_| rings.next().is_none()));
        let any = get_grid("IGEO7").unwrap();
        let mut disk = any.disk(base.id(), 1);
        assert_eq!(disk.by_ref().count(), 6);
        assert!((0..5).all(|_| disk.next().is_none()));

        fn fused<I: FusedIterator>(_: &I) {}
        fused(&base.disk(1));
        fused(&base.disk(1).rings());
        fused(&g.disk(base.id(), 1));
        fused(&g.disk(base.id(), 1).rings());
        fused(&any.disk(base.id(), 1));
        fused(&any.disk(base.id(), 1).rings());
    }

    /// A radius of `u32::MAX` on a plane without end: taking a few zones asks for the
    /// neighbours of only the zones those few are reached from, however large `k`.
    #[test]
    fn taking_from_an_unbounded_radius_costs_what_is_taken() {
        let calls = Cell::new(0usize);
        let counted = |z| {
            calls.set(calls.get() + 1);
            hex_neighbours(z)
        };
        let mut walk = Walk::new(hex(0, 0), u32::MAX);
        let first: Vec<ZoneId> = std::iter::from_fn(|| walk.next_zone(counted))
            .take(7)
            .collect();
        assert_eq!(first.len(), 7);
        assert_eq!(
            calls.get(),
            1,
            "the centre and its six neighbours need one query"
        );

        // Rings 0 to 2 whole, 1 + 6 + 12 zones, from the centre's query and the six of
        // ring 1, and no more.
        let mut walk = Walk::new(hex(0, 0), u32::MAX);
        calls.set(0);
        let n = std::iter::from_fn(|| walk.next_zone(counted))
            .take(19)
            .count();
        assert_eq!(n, 19);
        assert!(calls.get() <= 7, "{} queries for three rings", calls.get());
        assert_eq!(walk.ring, 2);

        // A walk consumed to its end at ring `k` does not ask for the neighbours of
        // ring `k`: the disk of radius 1, consumed whole, needs the centre's alone.
        let mut walk = Walk::new(hex(0, 0), 1);
        calls.set(0);
        assert_eq!(std::iter::from_fn(|| walk.next_zone(counted)).count(), 7);
        assert_eq!(calls.get(), 1);
    }

    /// On the lattice ring `n` has `6n` cells, which the walk recovers, ring by ring,
    /// and the flat walk and the rings agree zone for zone.
    #[test]
    fn rings_and_the_flat_walk_agree() {
        let mut walk = Walk::new(hex(0, 0), 12);
        let rings: Vec<Vec<ZoneId>> =
            std::iter::from_fn(|| walk.next_ring(hex_neighbours)).collect();
        let sizes: Vec<usize> = rings.iter().map(Vec::len).collect();
        assert_eq!(sizes[0], 1);
        for (n, &s) in sizes.iter().enumerate().skip(1) {
            assert_eq!(s, 6 * n);
        }
        assert_eq!(sizes.len(), 13);
        let flat = zones(&mut Walk::new(hex(0, 0), 12), hex_neighbours);
        assert_eq!(rings.concat(), flat);
        let distinct: std::collections::HashSet<ZoneId> = flat.iter().copied().collect();
        assert_eq!(distinct.len(), flat.len());
    }

    /// A walk partly consumed zone by zone and then taken ring by ring gives the rest
    /// of the ring in progress first, and whole rings after it.
    #[test]
    fn rings_resume_a_partly_consumed_walk() {
        let flat = zones(&mut Walk::new(hex(0, 0), 3), hex_neighbours);
        let mut walk = Walk::new(hex(0, 0), 3);
        let head: Vec<ZoneId> = (0..4)
            .map(|_| walk.next_zone(hex_neighbours).unwrap())
            .collect();
        let rest: Vec<Vec<ZoneId>> =
            std::iter::from_fn(|| walk.next_ring(hex_neighbours)).collect();
        assert_eq!(
            rest[0].len(),
            3,
            "ring 1 had six zones, three of them taken"
        );
        assert_eq!(rest[1].len(), 12);
        assert_eq!(rest[2].len(), 18);
        assert_eq!([head, rest.concat()].concat(), flat);
    }

    /// The order on a graph small enough to spell out, with a one-way edge: zone 5 of
    /// ring 2 lists zone 1 of ring 1, which does not list it back. A walk that
    /// remembered only the last two rings would yield zone 1 again in ring 3; this one
    /// remembers every zone it has yielded, and does not.
    #[test]
    fn the_order_is_first_arrival_ring_by_ring() {
        let graph = |z: ZoneId| -> Vec<ZoneId> {
            let ns: &[u64] = match z.0 {
                0 => &[2, 1, 3],
                1 => &[0, 4],
                2 => &[5, 0, 4],
                3 => &[0],
                4 => &[1, 2, 6],
                5 => &[2, 1, 7],
                6 => &[4],
                7 => &[5, 1],
                _ => &[],
            };
            ns.iter().map(|&n| ZoneId(n)).collect()
        };
        let got: Vec<u64> = zones(&mut Walk::new(ZoneId(0), u32::MAX), graph)
            .iter()
            .map(|z| z.0)
            .collect();
        assert_eq!(got, [0, 2, 1, 3, 5, 4, 7, 6]);
        let mut walk = Walk::new(ZoneId(0), u32::MAX);
        let rings: Vec<Vec<u64>> = std::iter::from_fn(|| walk.next_ring(graph))
            .map(|r| r.iter().map(|z| z.0).collect())
            .collect();
        assert_eq!(rings, [vec![0], vec![2, 1, 3], vec![5, 4], vec![7, 6]]);
        let at_one: Vec<u64> = zones(&mut Walk::new(ZoneId(0), 1), graph)
            .iter()
            .map(|z| z.0)
            .collect();
        assert_eq!(at_one, [0, 2, 1, 3]);
        assert_eq!(zones(&mut Walk::new(ZoneId(0), 0), graph), [ZoneId(0)]);
    }

    /// The order on a real grid, pinned zone by zone: the centre, its six neighbours in
    /// the order `neighbors` gives them, then ring 2 in first-arrival order. Ring 2 was
    /// derived by hand from the six neighbour lists of ring 1, in ring 1's order, as the
    /// engine gives them: `0064150` contributes `0064151`, `0064155` and `0064153`;
    /// `0064154` contributes `0042223` and `0042222`, its `0064155` having arrived
    /// already; `0064141` gives `0064145` and `0064140`; `0064143` gives `0064142` and
    /// `0064104`; `0064105` gives `0064100` and `0064101`; and `0064152` gives `0064114`,
    /// the rest of its list having arrived already.
    #[test]
    fn the_order_on_igeo7_is_pinned() {
        let z = igeo7().zone_from_text("0064156").unwrap();
        let got: Vec<String> = z.disk(2).map(|z| z.text_id()).collect();
        assert_eq!(
            got,
            [
                "0064156", // ring 0
                "0064150", "0064154", "0064141", "0064143", "0064105", "0064152", // ring 1
                "0064151", "0064155", "0064153", "0042223", "0042222", "0064145", // ring 2
                "0064140", "0064142", "0064104", "0064100", "0064101", "0064114",
            ]
        );
        let pentagon = igeo7().zone_from_text("0000000").unwrap();
        let got: Vec<String> = pentagon.disk(1).map(|z| z.text_id()).collect();
        assert_eq!(
            got,
            [
                "0000000", "0000005", "0000003", "0000004", "0000001", "0000006"
            ]
        );
    }

    /// At resolution 0 the grid is exhausted after ring 3, twelve zones in all, and
    /// the walk ends there whatever the radius.
    #[test]
    fn the_walk_ends_when_the_grid_is_exhausted() {
        let z = igeo7().zone_from_text("00").unwrap();
        let sizes: Vec<usize> = z.disk(u32::MAX).rings().map(|r| r.len()).collect();
        assert_eq!(sizes, [1, 5, 5, 1]);
        assert_eq!(z.disk(u32::MAX).count(), 12);
    }

    /// At resolution 19 a radius of `u32::MAX` covers the sphere, some 1.4e16 zones,
    /// and taking a hundred of them is prompt all the same.
    #[test]
    fn taking_a_little_from_a_vast_disk_on_a_real_grid_is_prompt() {
        let g = igeo7();
        let z = g.zone_from_geo(38.7223, -9.1393, 19).unwrap();
        let first: Vec<ZoneId> = g.disk(z.id(), u32::MAX).take(100).collect();
        assert_eq!(first.len(), 100);
        let whole_two: Vec<ZoneId> = g.disk(z.id(), 2).collect();
        assert_eq!(&first[..whole_two.len()], &whole_two[..]);
        let any = get_grid("IGEO7").unwrap();
        let through_any: Vec<ZoneId> = any.disk(z.id(), u32::MAX).take(100).collect();
        assert_eq!(through_any, first);
    }
}
