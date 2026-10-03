//! The tests of an aperture-3 grid that the aperture-7 suites have no counterpart for, or
//! whose counterpart there asks the engine what it cannot be asked here.
//!
//! On aperture 3 the engine's hierarchy is the one this crate implements, so its parents,
//! children, centroid parent and predicates are the oracle for this crate's, and so are its
//! sub-zone counts, first sub-zones, sub-zone orders and sub-zone indices. And five
//! departures were decided for aperture 3, each tested here with the engine's side of it
//! re-checked live on every run: the null zone answers as on aperture 7; the identifiers
//! the engine answers along two seams and cannot read back are the null zone here; a zone
//! at the finest resolution has no children, and a polar pentagon's phantom sub-hexes are
//! refused; the sub-zone answers the engine gives in contradiction of its own are
//! answered consistently; and the sub-zone orders the engine gets wrong at some zones on
//! the edges of the rhombi are corrected. The refusal of a non-finite coordinate, which the
//! aperture-7 suites test together with refusals that have no aperture-3 counterpart, is
//! tested here against the aperture-3 engine's own answer.
//!
//! Every relational question to the engine passes the guard `engine_may_be_asked` first;
//! the generators' centroids, bit for bit, are checked by the crate's own oracle test in
//! `hex_a3_subzones.rs`, which this suite does not repeat.

use std::collections::{BTreeMap, BTreeSet};

use super::expected_divergences as div;
use super::subject::Subject;
use super::suite::{SEED, engine_may_be_asked, engine_neighbours, zones};
use super::*;
use dggal_oracle as o;
use rs4dggs::indexings::I3h;
use rs4dggs::topologies::HexA3;
use rs4dggs::{Error, GridConfig, PlanarPoint, Projection, ZoneId};

/// A subject of aperture 3: the bound of every test here.
pub trait Aperture3: Subject<T = HexA3, I = I3h> {}
impl<S: Subject<T = HexA3, I = I3h>> Aperture3 for S {}

/// What has been measured of one aperture-3 grid against its engine by the tests of this
/// module: the floors of the samples drawn from the grid's own zones, which each grid
/// quantises under its own projection, and the counts of the engine's departures, which are
/// held exactly, so that a change in the engine's answers or in the sample shows. The
/// samples that do not depend on the grid, the pentagons of resolution 0, the zones of the
/// great-depth sample, the seam and polar points of the unreadable answers, the phantom
/// texts and the non-finite coordinates, keep their floors in the tests.
#[derive(Debug)]
struct Measured {
    /// `hierarchy_and_predicates_match_the_engine`: the zones compared.
    hierarchy_zones: usize,
    /// The same: the predicates the engine answers false and true, each distinct pair
    /// asked once; exactly.
    predicate_answers: [usize; 2],
    /// `sub_zones_and_paging_match_the_engine`: the orders compared.
    sub_zone_orders: usize,
    /// The same: the sub-zone indices asked of the engine.
    sub_zone_indices: usize,
    /// The same: the orders in which this crate departs from the engine's, the engine's
    /// failing its own descendant test and this crate's passing it; exactly.
    corrected_orders: usize,
    /// `unreadable_answers_along_two_seams_are_the_null_zone`: the answers the rule admits,
    /// at most.
    band_max_admitted: usize,
    /// `children_at_resolution_33_are_none`: the hexagons and the pentagons at resolution
    /// 33; exactly.
    finest_zones: [usize; 2],
    /// `sub_zone_answers_depart_from_the_engine_as_recorded`: the zones compared at depth 0,
    /// and those whose first sub-zone at depth 0 the engine answers with a neighbour;
    /// exactly.
    depth_0: [usize; 2],
    /// The same: the sub-zones at resolution 33 compared, and those to which the engine
    /// gives no index; exactly.
    finest_sub_zones: [usize; 2],
    /// The same: the neighbours compared at one level.
    same_level: usize,
}

/// ISEA3H, as measured against its own engine. Every value is the observation itself.
static ISEA3H_MEASURED: Measured = Measured {
    hierarchy_zones: 5_877,
    // At 1,176 zones.
    predicate_answers: [18_737, 18_700],
    // Holding 27,158 sub-zones, at the zones of the sample and at those of `EDGE_ZONES`, of
    // which the main sample holds `K5-67E99A74-D`.
    sub_zone_orders: 1_503,
    sub_zone_indices: 4_294,
    // All at the zones of `EDGE_ZONES`: `K5-67E99A74-D` and `K5-89336A36-D` at depths 1, 2 and
    // 3, `K5-8855C3D0-A` at depth 2 and `FB-0-B` at depth 3; none elsewhere.
    corrected_orders: 8,
    // The rule is met on both seams at every resolution from 2 to 21; the widest reach from
    // the edge is 2.24e-5 degrees on the northern seam and 2.05e-5 on the southern one.
    band_max_admitted: 2_437,
    finest_zones: [425, 12],
    depth_0: [1_176, 763],
    // 83.7 per cent.
    finest_sub_zones: [8_725, 7_300],
    same_level: 7_015,
};

/// IVEA3H, as measured against its own engine, on the zones that its projection makes of the
/// main sample. Every value is the observation itself.
static IVEA3H_MEASURED: Measured = Measured {
    hierarchy_zones: 5_876,
    // At 1,176 zones.
    predicate_answers: [18_720, 18_617],
    // Holding 27,179 sub-zones, at the zones of the sample and at those of `EDGE_ZONES`, of
    // which the main sample holds `K5-8855C3D0-A`.
    sub_zone_orders: 1_505,
    sub_zone_indices: 4_294,
    // As on ISEA3H, the same eight orders, all at the zones of `EDGE_ZONES`; none elsewhere.
    corrected_orders: 8,
    // The rule is met on both seams at every resolution from 2 to 21; the widest reach from
    // the edge is 2.18e-5 degrees on the northern seam and 2.00e-5 on the southern one.
    band_max_admitted: 2_256,
    finest_zones: [425, 12],
    depth_0: [1_175, 760],
    // 83.4 per cent.
    finest_sub_zones: [8_725, 7_277],
    same_level: 7_012,
};

/// RTEA3H, as measured against its own engine, on the zones that its projection makes of the
/// main sample. Every value is the observation itself.
static RTEA3H_MEASURED: Measured = Measured {
    hierarchy_zones: 5_872,
    // At 1,175 zones.
    predicate_answers: [18_705, 18_439],
    // Holding 27,149 sub-zones, at the zones of the sample and at those of `EDGE_ZONES`, of
    // which the main sample holds `K5-89336A36-D` and `FB-0-B`.
    sub_zone_orders: 1_509,
    sub_zone_indices: 4_279,
    // As on ISEA3H, the same eight orders, all at the zones of `EDGE_ZONES`; none elsewhere.
    corrected_orders: 8,
    // The rule is met on both seams at every resolution from 2 to 21; the widest reach from
    // the edge is 2.19e-5 degrees on the northern seam and 1.74e-5 on the southern one.
    band_max_admitted: 2_287,
    finest_zones: [425, 12],
    depth_0: [1_175, 757],
    // 83.5 per cent.
    finest_sub_zones: [8_711, 7_270],
    same_level: 7_001,
};

/// What has been measured of the aperture-3 grid under test. A grid of which nothing has
/// been measured fails here rather than borrowing another grid's counts; each test prints
/// what it observed before it asserts, so that the table can be filled from the run.
fn measured<S: Aperture3>() -> &'static Measured {
    match S::NAME {
        "ISEA3H" => &ISEA3H_MEASURED,
        "IVEA3H" => &IVEA3H_MEASURED,
        "RTEA3H" => &RTEA3H_MEASURED,
        other => panic!(
            "nothing has been measured of {other} against the engine by the aperture-3 tests"
        ),
    }
}

/// The pentagons of resolution 0, indexed by the base cell of the aperture-7 grids whose
/// vertex each is centred on: the ten root rhombi's, `A0-0-A` on base cell 1, `A1-0-A` on 6,
/// `A2-0-A` on 2 and so on alternately, and the two polar ones, `AA-0-A` on base cell 0 and
/// `AB-0-A` on base cell 11.
pub const BASE_PENTAGONS: [&str; 12] = [
    "AA-0-A", "A0-0-A", "A2-0-A", "A4-0-A", "A6-0-A", "A8-0-A", "A1-0-A", "A3-0-A", "A5-0-A",
    "A7-0-A", "A9-0-A", "AB-0-A",
];

/// The seams and the icosahedron's vertices of the sample rest on the twelve points of
/// [`icosahedron_vertices`], which are IGEO7's base cells. That is sound for an aperture-3
/// grid only if its engine centres its twelve pentagons of resolution 0 on the same points,
/// as the shared orientation implies: checked here pentagon by pentagon, bit for bit, and
/// this crate's centroids with them.
pub fn the_engine_places_its_pentagons_on_the_icosahedron_vertices<S: Aperture3>() {
    let mut floor = Floor::new("pentagons on the icosahedron vertices", 12);
    for (base, (&text, &vertex)) in BASE_PENTAGONS
        .iter()
        .zip(icosahedron_vertices().iter())
        .enumerate()
    {
        let z = o::zone_from_text(S::ORACLE, text);
        assert_ne!(z, o::NULL_ZONE, "{} reads no {text}", S::ORACLE);
        let bits = |p: (f64, f64)| (p.0.to_bits(), p.1.to_bits());
        let theirs = o::centroid(S::ORACLE, z);
        assert_eq!(
            bits(theirs),
            bits(vertex),
            "{text} of {} is centred at {theirs:?}, where base cell {base:02} of {} is at \
             {vertex:?}",
            S::ORACLE,
            o::IGEO7
        );
        let ours = S::grid().zone_from_text(text).unwrap();
        assert_eq!(ours.id().0, z, "{text} read by this crate");
        let c = ours.centroid();
        assert_eq!(
            bits((c.lat, c.lon)),
            bits(vertex),
            "{text}, centred by this crate"
        );
        floor.hit();
    }
    eprintln!(
        "pentagons: {} centres its twelve pentagons of resolution 0 on the icosahedron's \
         vertices, bit for bit",
        S::ORACLE
    );
}

/// Asserts that the engine may be asked about `z`, which is one of its own answers about a
/// zone it can read, and fails naming it otherwise: an engine answer that it cannot itself
/// be asked about would be a finding to characterise, never a case to pass over.
fn askable<S: Aperture3>(z: u64, whose: &str) {
    assert!(
        engine_may_be_asked::<S>(z),
        "the engine lists {} ({z:#018x}) among the {whose}, and cannot read it back \
         (seed {SEED:#x})",
        o::text_id(S::ORACLE, z)
    );
}

/// A hierarchy predicate the engine answers, of an ordered pair `(a, b)`: whether `a` is an
/// ancestor of `b` (`isZoneAncestorOf`, with no bound on the depth), whether `a` is an
/// immediate child of `b` (`isZoneImmediateChildOf`), and whether the two are siblings
/// (`areZonesSiblings`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Predicate {
    Ancestor,
    ImmediateChild,
    Sibling,
}

/// The identifiers of a list of zones.
fn raw(zs: &[ZoneId]) -> Vec<u64> {
    zs.iter().map(|z| z.0).collect()
}

/// The hierarchy against the engine's, at every zone of the sample: the parents, as a list in
/// the engine's order; the children below the finest resolution, likewise; the centroid
/// parent; and whether the zone is a centroid child. And, at every fifth zone, the three
/// predicates against the engine's own `isZoneAncestorOf`, `isZoneImmediateChildOf` and
/// `areZonesSiblings`: each parent and each parent of a parent as an ancestor, and the zone
/// as none of theirs; each neighbour of the first parent, which may be another parent or
/// none; and the zone itself, each of its neighbours and each child of its first parent, as
/// siblings or not. Each pair is asked once, however many of these routes lead to it.
///
/// The engine's ancestry predicate recurses through its `getZoneParents` from the descendant
/// up to the level below the ancestor, so it is asked about ancestors at most two levels up,
/// where it passes through the zone and its parents alone, each found askable first. The null
/// zone, which the sample holds where the engine answers its own or an identifier it cannot
/// read back, has no relations here, and the engine is not asked about it. At the finest
/// resolution this crate answers no children, which `children_at_resolution_33_are_none`
/// sets against the engine's.
pub fn hierarchy_and_predicates_match_the_engine<S: Aperture3>() {
    let g = S::grid();
    let m = measured::<S>();
    let mut floor = Floor::new("hierarchy", m.hierarchy_zones);
    let mut predicates = Floor::new(
        "hierarchy predicates",
        m.predicate_answers[0] + m.predicate_answers[1],
    );
    let (mut null_zone, mut finest) = (0usize, 0usize);
    // The predicates the engine answered false and true, which must both be reached.
    let mut answers = [0usize; 2];
    for (k, id) in zones::<S>().into_iter().enumerate() {
        let z = g.zone(id);
        let text = z.text_id();
        let at = format!("{text} (seed {SEED:#x})");
        if !engine_may_be_asked::<S>(id.0) {
            assert!(
                id == ZoneId::NULL
                    && g.parents(id).is_empty()
                    && g.children(id).is_empty()
                    && g.centroid_parent(id).is_none()
                    && !g.is_centroid_child(id),
                "the engine cannot read {at} back, or it has relations"
            );
            null_zone += 1;
            floor.hit();
            continue;
        }
        let parents = o::parents(S::ORACLE, id.0);
        assert_eq!(raw(&g.parents(id)), parents, "parents of {at}");
        for &p in &parents {
            askable::<S>(p, &format!("parents of {text}"));
        }
        if z.resolution() < g.max_resolution() {
            assert_eq!(
                raw(&g.children(id)),
                o::children(S::ORACLE, id.0),
                "children of {at}"
            );
        } else {
            assert!(g.children(id).is_empty(), "{at} has children");
            finest += 1;
        }
        assert_eq!(
            g.centroid_parent(id).map_or(o::NULL_ZONE, |p| p.0),
            o::centroid_parent(S::ORACLE, id.0),
            "centroid parent of {at}"
        );
        assert_eq!(
            g.is_centroid_child(id),
            o::is_centroid_child(S::ORACLE, id.0),
            "whether {at} is a centroid child"
        );
        floor.hit();
        if k % 5 != 0 {
            continue;
        }

        // The pairs to ask about, each once, however many routes lead to it: the other
        // parents are also neighbours of the first, and the zone is among its first parent's
        // children, where it is asked about as its own sibling already.
        let mut pairs: Vec<(Predicate, u64, u64)> = Vec::new();
        let mut ask = |q: (Predicate, u64, u64)| {
            if !pairs.contains(&q) {
                pairs.push(q);
            }
        };
        let mut grandparents = BTreeSet::new();
        for &p in &parents {
            ask((Predicate::Ancestor, p, id.0));
            ask((Predicate::Ancestor, id.0, p));
            ask((Predicate::ImmediateChild, id.0, p));
            grandparents.extend(o::parents(S::ORACLE, p));
        }
        for gp in grandparents {
            ask((Predicate::Ancestor, gp, id.0));
            ask((Predicate::ImmediateChild, id.0, gp));
        }
        ask((Predicate::Sibling, id.0, id.0));
        if let Some(&first) = parents.first() {
            for n in engine_neighbours::<S>(first) {
                askable::<S>(
                    n,
                    &format!("neighbours of {}", o::text_id(S::ORACLE, first)),
                );
                ask((Predicate::Ancestor, n, id.0));
                ask((Predicate::ImmediateChild, id.0, n));
            }
            for c in o::children(S::ORACLE, first) {
                askable::<S>(c, &format!("children of {}", o::text_id(S::ORACLE, first)));
                ask((Predicate::Sibling, id.0, c));
            }
        }
        for n in engine_neighbours::<S>(id.0) {
            askable::<S>(n, &format!("neighbours of {text}"));
            ask((Predicate::Sibling, id.0, n));
        }
        for (predicate, a, b) in pairs {
            let (za, zb) = (g.zone(ZoneId(a)), g.zone(ZoneId(b)));
            let (ours, theirs, what) = match predicate {
                Predicate::Ancestor => (
                    za.is_ancestor_of(&zb),
                    o::is_ancestor_of(S::ORACLE, a, b),
                    "an ancestor of",
                ),
                Predicate::ImmediateChild => (
                    za.is_immediate_child_of(&zb),
                    o::is_immediate_child_of(S::ORACLE, a, b),
                    "an immediate child of",
                ),
                Predicate::Sibling => (
                    za.is_sibling_of(&zb),
                    o::are_siblings(S::ORACLE, a, b),
                    "a sibling of",
                ),
            };
            assert_eq!(
                ours,
                theirs,
                "whether {} is {what} {}, ours against the engine's (seed {SEED:#x})",
                o::text_id(S::ORACLE, a),
                o::text_id(S::ORACLE, b)
            );
            answers[usize::from(theirs)] += 1;
            predicates.hit();
        }
    }
    eprintln!(
        "hierarchy: the null zone {null_zone} time(s), not asked of the engine; {finest} zones \
         at the finest resolution; predicates the engine answers false and true: {answers:?}"
    );
    assert!(
        null_zone > 0,
        "the null-zone branch never ran (seed {SEED:#x})"
    );
    // Asserted exactly, since the predicates are asked once for each distinct pair and a
    // pair asked twice would swell a floor.
    assert_eq!(
        answers, m.predicate_answers,
        "the predicates the engine answers false and true (seed {SEED:#x})"
    );
}

/// Zones on the edges of the rhombi at which the engine's sub-zone order is not the zone's
/// descendants, on every aperture-3 grid, for each of the three ways in which a unit in the
/// last place carries a point across an interruption of the layout that it should only reach:
/// a vertex of the zone, from which the first centroid is read, at `K5-67E99A74-D` and
/// `K5-89336A36-D` (sub-hex D on the left edge of rhombus 5 at resolution 21, at depths 1 to
/// 3); a step of the generator of an even-resolution zone, then carried across a second time,
/// at `K5-8855C3D0-A` (on the same edge at resolution 20, at depth 2); and the same, of the
/// generator of an odd-resolution zone, at `FB-0-B` (the southern polar pentagon at resolution
/// 11, at depth 3). Found in the samples of ISEA3H (`K5-67E99A74-D`), IVEA3H (`K5-8855C3D0-A`)
/// and RTEA3H (`K5-89336A36-D` and `FB-0-B`), and added to the sub-zone sample of every grid,
/// so that every grid meets all three ways whatever its projection makes of the main sample.
pub const EDGE_ZONES: [&str; 4] = ["K5-67E99A74-D", "K5-89336A36-D", "K5-8855C3D0-A", "FB-0-B"];

/// The sub-zone order and its paging against the engine's.
///
/// At every tenth zone of the sample and at the zones of [`EDGE_ZONES`], to depth 3: the count
/// against the engine's `countSubZones`, and the whole order against its `getSubZones`. Where
/// the two orders are the same, the first sub-zone against the engine's `getFirstSubZone`.
/// Where they differ, the orders the engine gets wrong at some zones on the edges of the
/// rhombi, the departure is re-checked live by
/// `expected_divergences::engine_sub_zone_order_is_faulty`: the engine's order fails its own
/// descendant test, and this crate's passes it; the engine's first sub-zone is re-checked to be
/// the first of its own order, as this crate's is of its own; and the number of such orders is
/// held exactly. Either way, this crate's order names no zone twice, and its first sub-zone is
/// its first entry. The sub-zone at every index against that order, as the engine's
/// `getSubZoneAtIndex` cannot be reached through the vendored binding, whose call drops its
/// depth, and was measured equal to `getSubZones` at every pair tried; the index one past the
/// last refused; and the index of every sub-zone in this crate, its place in the order. In the
/// engine's `getSubZoneIndex`, where the orders are the same and the sub-zones lie at
/// resolution 32 or coarser (at 33 the engine departs: see
/// `sub_zone_answers_depart_from_the_engine_as_recorded`), the index of the first, second and
/// last sub-zone, and at the zones of [`EDGE_ZONES`] of every sub-zone. Where the orders
/// differ, the engine's index refers to its own order, and is not asked.
///
/// And at great depth, from each of the twelve pentagons of resolution 0 and a hexagon of
/// each root rhombus at resolution 1, at every depth to the finest resolution: the count and
/// the first sub-zone
/// against the engine's, and the sub-zones at indices 1 to 3, reached through the
/// generators' skip without the order being built, against the engine's own index of each,
/// whose search costs as much as the index and so stays cheap.
pub fn sub_zones_and_paging_match_the_engine<S: Aperture3>() {
    the_limit_on_a_list_of_sub_zones_is_its_ceiling();
    let g = S::grid();
    let m = measured::<S>();
    // The deep sample, the same on every grid: 716 depths from its twenty-two zones.
    let mut orders = Floor::new("sub-zone orders", m.sub_zone_orders);
    let mut indices = Floor::new("sub-zone indices asked of the engine", m.sub_zone_indices);
    let mut deep = Floor::new("sub-zones at great depth", 716);
    let mut sub_zones_compared = 0usize;
    // The orders in which this crate departs from the engine's, each re-checked live, with
    // what is wrong with the engine's.
    let mut corrected: Vec<String> = Vec::new();
    let mut sample: Vec<ZoneId> = zones::<S>()
        .into_iter()
        .step_by(10)
        .filter(|&id| id != ZoneId::NULL)
        .collect();
    let edge: Vec<ZoneId> = EDGE_ZONES
        .iter()
        .map(|t| g.zone_from_text(t).unwrap().id())
        .collect();
    for &id in &edge {
        if !sample.contains(&id) {
            sample.push(id);
        }
    }
    for id in sample {
        let text = g.text_id(id);
        assert!(engine_may_be_asked::<S>(id.0), "{text} (seed {SEED:#x})");
        let res = g.resolution(id);
        for depth in 1..=3u8 {
            if res + depth > g.max_resolution() {
                continue;
            }
            let at = format!("{text} at depth {depth} (seed {SEED:#x})");
            let count = g.count_sub_zones(id, depth).unwrap();
            assert_eq!(
                count,
                o::count_sub_zones(S::ORACLE, id.0, i32::from(depth)),
                "count of {at}"
            );
            let subs = g.sub_zones(id, depth).unwrap();
            assert_eq!(subs.len() as u64, count, "count of {at} against its order");
            let distinct: BTreeSet<ZoneId> = subs.iter().copied().collect();
            assert_eq!(distinct.len(), subs.len(), "a zone named twice in {at}");
            let first = g.first_sub_zone(id, depth).unwrap();
            assert_eq!(first, subs[0], "first sub-zone of {at} against its order");
            let theirs = o::sub_zones(S::ORACLE, id.0, i32::from(depth));
            let same = raw(&subs) == theirs;
            if same {
                assert_eq!(
                    first.0,
                    o::first_sub_zone(S::ORACLE, id.0, i32::from(depth)),
                    "first sub-zone of {at}"
                );
            } else {
                let fault = div::engine_sub_zone_order_is_faulty::<S>(id.0, depth, &raw(&subs))
                    .unwrap_or_else(|why| panic!("{why} (seed {SEED:#x})"));
                assert_eq!(
                    o::first_sub_zone(S::ORACLE, id.0, i32::from(depth)),
                    theirs[0],
                    "the engine's first sub-zone of {at} against its own order"
                );
                corrected.push(format!("{at}: {fault:?}"));
            }
            for (i, &s) in subs.iter().enumerate() {
                assert_eq!(
                    g.sub_zone_at_index(id, depth, i as u64),
                    Ok(s),
                    "sub-zone {i} of {at}"
                );
                assert_eq!(
                    g.sub_zone_index(id, s),
                    Ok(Some(i as u64)),
                    "index of {} in {at}",
                    g.text_id(s)
                );
            }
            assert_eq!(
                g.sub_zone_at_index(id, depth, count),
                Err(Error::IndexOutOfRange {
                    index: count,
                    count
                }),
                "the index past the last of {at}"
            );
            sub_zones_compared += subs.len();
            if same && res + depth < g.max_resolution() {
                // At an edge zone every sub-zone; elsewhere the first, second and last.
                let places: BTreeSet<usize> = if edge.contains(&id) {
                    (0..subs.len()).collect()
                } else {
                    BTreeSet::from([0, 1, subs.len() - 1])
                };
                for i in places {
                    let s = subs[i];
                    askable::<S>(s.0, &format!("sub-zones of {text}"));
                    assert_eq!(
                        o::sub_zone_index(S::ORACLE, id.0, s.0),
                        i as i64,
                        "the engine's index of {} in {at}",
                        g.text_id(s)
                    );
                    indices.hit();
                }
            }
            orders.hit();
        }
    }

    let mut from: Vec<ZoneId> = BASE_PENTAGONS
        .iter()
        .map(|t| g.zone_from_text(t).unwrap().id())
        .collect();
    from.extend((0..10).map(|root| g.zone_from_text(&format!("A{root:X}-0-C")).unwrap().id()));
    for id in from {
        let text = g.text_id(id);
        let span = g.max_resolution() - g.resolution(id);
        for depth in 1..=span {
            let at = format!("{text} at depth {depth}");
            let count = g.count_sub_zones(id, depth).unwrap();
            assert_eq!(
                count,
                o::count_sub_zones(S::ORACLE, id.0, i32::from(depth)),
                "count of {at}"
            );
            assert_eq!(
                g.first_sub_zone(id, depth).unwrap().0,
                o::first_sub_zone(S::ORACLE, id.0, i32::from(depth)),
                "first sub-zone of {at}"
            );
            if depth < span {
                let first: Vec<ZoneId> = (0..=3u64)
                    .map(|i| g.sub_zone_at_index(id, depth, i).unwrap())
                    .collect();
                for (i, &s) in first.iter().enumerate().skip(1) {
                    askable::<S>(s.0, &format!("sub-zones of {text}"));
                    assert_eq!(
                        o::sub_zone_index(S::ORACLE, id.0, s.0),
                        i as i64,
                        "the engine's index of {} in {at}",
                        g.text_id(s)
                    );
                }
            }
            deep.hit();
        }
    }
    eprintln!(
        "sub-zones: {} orders, {sub_zones_compared} sub-zones in them; {} indices asked of the \
         engine; {} counts and first sub-zones at great depth",
        orders.n, indices.n, deep.n
    );
    eprintln!(
        "sub-zones: {} orders corrected, the engine's failing its own descendant test and this \
         crate's passing it:\n  {}",
        corrected.len(),
        corrected.join("\n  ")
    );
    // Held exactly, so that the departure can neither grow nor vanish unnoticed.
    assert_eq!(
        corrected.len(),
        m.corrected_orders,
        "the orders corrected (seed {SEED:#x})"
    );
}

/// The null zone answers as it does on aperture 7, and as this crate decided: its text
/// `(null)`, read back to it; resolution 0; the point (0, 0) and six zero vertices; no
/// neighbours, parents, children or centroid parent, and no centroid child; a disk of itself
/// alone; no ancestor, descendant or sibling, a real zone tried both ways; and every sub-zone
/// method refused as asked about no zone.
///
/// The engine's side, re-checked live, is the departure: it prints and reads the same text,
/// but reads its null zone at level 63, gives it a centroid and six vertices computed from
/// its bits (see `expected_divergences::engine_null_zone_geometry`), and cannot read it at a
/// level of 33 or less, which is why no relational question is put to it: its parents are
/// identifiers on which the engine's own parent arithmetic divides by zero.
pub fn the_null_zone_departs_from_the_engine_as_recorded<S: Aperture3>() {
    let g = S::grid();
    let mut floor = Floor::new("null zone", 5);
    let null = ZoneId::NULL;
    let z = g.zone(null);

    assert_eq!(g.text_id(null), rs4dggs::NULL_TEXT);
    assert_eq!(o::text_id(S::ORACLE, null.0), rs4dggs::NULL_TEXT);
    assert_eq!(
        g.zone_from_text(rs4dggs::NULL_TEXT).map(|z| z.id()),
        Ok(null)
    );
    floor.hit();

    assert_eq!(g.resolution(null), 0);
    let c = g.centroid(null);
    assert_eq!((c.lat, c.lon), (0.0, 0.0), "the null zone's centroid");
    let vs: Vec<(f64, f64)> = g.vertices(null).iter().map(|p| (p.lat, p.lon)).collect();
    assert_eq!(vs, vec![(0.0, 0.0); 6], "the null zone's vertices");
    assert!(!g.is_pentagon(null));
    div::engine_null_zone_geometry::<S>().unwrap_or_else(|why| panic!("{why}"));
    floor.hit();

    assert!(
        g.neighbors(null).is_empty()
            && g.parents(null).is_empty()
            && g.children(null).is_empty()
            && g.centroid_parent(null).is_none()
            && !g.is_centroid_child(null),
        "the null zone has relations"
    );
    assert_eq!(z.disk(3).map(|z| z.id()).collect::<Vec<_>>(), vec![null]);
    assert!(
        !o::engine_can_read(S::ORACLE, null.0),
        "the engine reads its null zone"
    );
    floor.hit();

    let real = g.zone_from_text("A0-0-A").unwrap();
    let child = real.children()[1];
    for other in [real, child] {
        assert!(
            !z.is_ancestor_of(&other)
                && !other.is_ancestor_of(&z)
                && !z.is_sibling_of(&other)
                && !other.is_sibling_of(&z)
                && !z.is_immediate_child_of(&other)
                && !other.is_immediate_child_of(&z),
            "the null zone is related to {}",
            other.text_id()
        );
    }
    assert!(!z.is_sibling_of(&z) && !z.is_ancestor_of(&z));
    floor.hit();

    let invalid = |r: rs4dggs::Result<()>| matches!(r, Err(Error::InvalidZone(_)));
    for depth in [0u8, 1, 2, 33] {
        assert!(
            invalid(g.count_sub_zones(null, depth).map(drop)),
            "count at {depth}"
        );
        assert!(
            invalid(g.first_sub_zone(null, depth).map(drop)),
            "first at {depth}"
        );
        assert!(
            invalid(g.sub_zones(null, depth).map(drop)),
            "sub-zones at {depth}"
        );
        assert!(
            invalid(g.sub_zone_at_index(null, depth, 0).map(drop)),
            "sub-zone at index at {depth}"
        );
    }
    for (a, b) in [(null, null), (null, real.id()), (real.id(), null)] {
        assert!(
            invalid(g.sub_zone_index(a, b).map(drop)),
            "sub-zone index of {b:?} in {a:?}"
        );
    }
    floor.hit();
}

/// Where the engine answers an identifier it cannot read back, along the seam over the
/// north pole and along the seam that ends at the southern polar vertex, at resolutions 2 to
/// 21, this crate answers the null zone. Checked on a sample built to meet both seams at
/// every one of those resolutions: points `c / 3^(res / 2)` inside each seam of the planar
/// layout, `c` from 1e-9 to 1e-6, at nineteen stations along it, lifted to the sphere by this
/// crate's inverse projection; and 200 points a resolution with colatitudes from 1e-10 to
/// 1e-2 degrees about the north pole, where the seam over the pole narrows fastest. Every
/// resolution from 0 to 33 is sampled.
///
/// At every point the two answers are the same, or the class rule
/// `expected_divergences::engine_answers_an_unreadable_polar_root` holds, live: the engine's
/// answer names a polar root with a non-zero index, neither side reads it back, and the
/// engine draws it as that polar pentagon. The rule must be met on both seams at every
/// resolution from 2 to 21, and it admits nothing at any other.
pub fn unreadable_answers_along_two_seams_are_the_null_zone<S: Aperture3>() {
    let g = S::grid();
    let geom = S::P::build_geometry(&GridConfig::default());
    // Measured: 352 points at each of 34 resolutions.
    let mut floor = Floor::new("seams of unreadable answers", 11_968);
    let mut met: BTreeMap<(u64, u8), usize> = BTreeMap::new();
    let (mut agreed, mut nulls_agreed) = (0usize, 0usize);
    let mut widest = [0.0f64; 2];
    let v = icosahedron_vertices();
    for res in 0..=g.max_resolution() {
        let mut points = Vec::new();
        for c in [1e-6, 5e-7, 1e-7, 1e-9] {
            let delta = c / 3u64.pow(u32::from(res / 2)) as f64;
            for i in 1..20 {
                let t = f64::from(i) / 20.0;
                for (x, y) in [(5.0 - delta, 4.0 + t), (4.0 + t, 6.0 - delta)] {
                    let p = S::P::inverse(&geom, PlanarPoint { face: -1, x, y }, res % 2 == 1);
                    points.push((p.lat, p.lon));
                }
            }
        }
        let mut r = Rng(SEED ^ u64::from(res));
        for _ in 0..200 {
            let colatitude = 10f64.powf(-10.0 + 8.0 * r.unit());
            points.push((90.0 - colatitude, 360.0 * r.unit() - 180.0));
        }
        for (lat, lon) in points {
            let ours = g.zone_from_geo(lat, lon, res).unwrap().id().0;
            let theirs = o::zone_from_geo(S::ORACLE, lat, lon, i32::from(res));
            floor.hit();
            if ours == theirs {
                agreed += 1;
                nulls_agreed += usize::from(ours == o::NULL_ZONE);
                continue;
            }
            div::engine_answers_an_unreadable_polar_root::<S>(lat, lon, res, ours, theirs)
                .unwrap_or_else(|why| {
                    panic!("zone_from_geo({lat}, {lon}, {res}): {why} (seed {SEED:#x})")
                });
            let root = (theirs >> 53) & 0xF;
            *met.entry((root, res)).or_default() += 1;
            let (edge, i) = if root == 10 {
                ((0, 1), 0)
            } else {
                ((6, 11), 1)
            };
            widest[i] = widest[i].max(dist_to_arc_deg((lat, lon), v[edge.0], v[edge.1]));
        }
    }
    eprintln!(
        "seams of unreadable answers: {agreed} answers agreed ({nulls_agreed} null), and the \
         rule met {} times; the widest reach from the edge {widest:?} degrees",
        met.values().sum::<usize>()
    );
    eprintln!("  by root and resolution: {met:?}");
    // The admissions are held under the measured count, so that a loosened rule shows.
    let admitted: usize = met.values().sum();
    assert!(
        admitted <= measured::<S>().band_max_admitted,
        "the rule admits {admitted} answers (seed {SEED:#x})"
    );
    for root in [10, 11] {
        let missing: Vec<u8> = div::UNREADABLE_POLAR_RESOLUTIONS
            .filter(|&res| !met.contains_key(&(root, res)))
            .collect();
        assert!(
            missing.is_empty(),
            "the sample no longer meets the unreadable answers of root {root} at resolutions \
             {missing:?} (seed {SEED:#x})"
        );
    }
}

/// A zone at the finest resolution, 33, has no children here. The engine answers seven for
/// a hexagon and six for a pentagon, whose indices overflow their 51 bits into the root and
/// level fields: each is re-checked live to be at level 34 and unreadable by the engine,
/// and refused by this crate, and nothing else is asked about it, since the engine's parent
/// arithmetic divides by zero on it. Sampled at every zone of the main sample at resolution
/// 33, the twelve pentagons among them, since the sample holds the icosahedron's vertices.
pub fn children_at_resolution_33_are_none<S: Aperture3>() {
    let g = S::grid();
    let m = measured::<S>();
    let mut floor = Floor::new(
        "children at resolution 33",
        m.finest_zones[0] + m.finest_zones[1],
    );
    let ids: Vec<ZoneId> = zones::<S>()
        .into_iter()
        .filter(|&id| id != ZoneId::NULL && g.resolution(id) == 33)
        .collect();
    let mut theirs_counted = [0usize; 2];
    for id in ids {
        let z = g.zone(id);
        let text = z.text_id();
        assert!(z.children().is_empty(), "{text} has children");
        assert!(g.children(id).is_empty(), "{text} has children in the grid");
        assert!(
            engine_may_be_asked::<S>(id.0),
            "the engine cannot read {text} back"
        );
        let theirs = o::children(S::ORACLE, id.0);
        let want = if z.is_pentagon() { 6 } else { 7 };
        assert_eq!(theirs.len(), want, "the engine's children of {text}");
        for c in theirs {
            let ct = o::text_id(S::ORACLE, c);
            assert_eq!(
                o::level(S::ORACLE, c),
                34,
                "the engine's child {ct} of {text}"
            );
            assert_eq!(
                o::zone_from_text(S::ORACLE, &ct),
                o::NULL_ZONE,
                "the engine's own parser reads its child {ct} of {text} back"
            );
            assert!(
                g.zone_from_text(&ct).is_err(),
                "this crate reads {ct}, the engine's child of {text}"
            );
        }
        theirs_counted[usize::from(z.is_pentagon())] += 1;
        floor.hit();
    }
    eprintln!(
        "children at resolution 33: none here, where the engine answers level-34 identifiers \
         for {} hexagons and {} pentagons",
        theirs_counted[0], theirs_counted[1]
    );
    assert_eq!(
        theirs_counted, m.finest_zones,
        "the hexagons and pentagons at resolution 33 (seed {SEED:#x})"
    );
}

/// The sub-hexes C and D of a polar pentagon, `AA-0-C` and its like at every level, are
/// refused here. The engine reads them, since it validates nothing for the polar roots, and
/// prints them back as read, which is re-checked live for every one; and none of them is a
/// neighbour, parent or child, in the engine's own lists, of a polar pentagon at any level,
/// which is where such a cell would be met.
pub fn phantom_polar_texts_are_refused<S: Aperture3>() {
    let g = S::grid();
    let mut floor = Floor::new("phantom polar texts", 17 * 2 * 2);
    for level in 0..=16u8 {
        let letter = char::from(b'A' + level);
        for root in ['A', 'B'] {
            for sub_hex in ['C', 'D'] {
                let text = format!("{letter}{root}-0-{sub_hex}");
                assert!(
                    matches!(g.zone_from_text(&text), Err(Error::InvalidZone(_))),
                    "this crate reads {text}"
                );
                let e = o::zone_from_text(S::ORACLE, &text);
                assert_ne!(e, o::NULL_ZONE, "the engine no longer reads {text}");
                assert_eq!(
                    o::text_id(S::ORACLE, e),
                    text,
                    "the engine prints {text} back"
                );
                floor.hit();
            }
        }
    }
    let phantom = |z: u64| {
        let (root, index, sub_hex) = ((z >> 53) & 0xF, (z >> 2) & ((1 << 51) - 1), z & 3);
        z != o::NULL_ZONE && (root == 10 || root == 11) && index == 0 && sub_hex > 1
    };
    let mut relations = Floor::new("relations of the polar pentagons", 68);
    for level in 0..=16u8 {
        let letter = char::from(b'A' + level);
        for root in ['A', 'B'] {
            for sub_hex in ['A', 'B'] {
                let text = format!("{letter}{root}-0-{sub_hex}");
                let p = o::zone_from_text(S::ORACLE, &text);
                assert!(
                    engine_may_be_asked::<S>(p),
                    "the engine cannot read {text} back"
                );
                let mut related = engine_neighbours::<S>(p);
                related.extend(o::parents(S::ORACLE, p));
                // At resolution 33 the engine's children are identifiers of level 34 (see
                // `children_at_resolution_33_are_none`), and are not asked for here.
                if !(level == 16 && sub_hex == 'B') {
                    related.extend(o::children(S::ORACLE, p));
                }
                if let Some(&z) = related.iter().find(|&&z| phantom(z)) {
                    panic!("the engine relates {text} to {}", o::text_id(S::ORACLE, z));
                }
                relations.hit();
            }
        }
    }
}

/// The sub-zone answers in which the engine contradicts itself, answered consistently here,
/// each with the engine's side re-checked live.
///
/// - Depth 0, at every fifth zone of the sample: the zone is its own first sub-zone and its
///   sub-zone at index 0 here, as it is the whole of its sub-zones at that depth on both
///   sides, `getSubZones` answering the zone alone and `countSubZones` one. The engine's
///   `getFirstSubZone` has no case for depth 0 and quantises a vertex of the zone at the
///   zone's own level, so that for most zones it names another: re-checked to be one of the
///   zone's neighbours in the engine's own list, at a floor of zones.
/// - A sub-zone at resolution 33, at every sub-zone of the zones of the sample at
///   resolutions 31 and 32, to depth 2 and 1: it has its index here, its place in the
///   order, where the engine's `getSubZoneIndex` answers -1 for most of them (some 83 to 84
///   per cent on each grid), having first asked `zoneHasSubZone`, which fails on the
///   needle's level-34 centroid child. Re-checked: wherever the engine answers -1, its own
///   `zoneHasSubZone` is false; wherever it answers an index, it is this crate's. A zone at
///   33 itself is not taken, since at depth 0 it is its own only sub-zone and both sides
///   answer 0 at one level; the number of sub-zones compared is asserted exactly with the
///   departures, so that such comparisons cannot swell it.
/// - Two different zones at one level, a zone and each of its neighbours: `None` here, where
///   the engine answers 0, having compared the two levels only, as on aperture 7.
pub fn sub_zone_answers_depart_from_the_engine_as_recorded<S: Aperture3>() {
    the_limit_on_a_list_of_sub_zones_is_its_ceiling();
    let g = S::grid();
    let zones: Vec<ZoneId> = zones::<S>()
        .into_iter()
        .filter(|&id| id != ZoneId::NULL)
        .collect();

    // The engine's departures are counted too, and asserted exactly, with the zones and
    // sub-zones compared, so that a change in the engine's answers shows.
    let m = measured::<S>();
    let mut depth_0 = Floor::new("sub-zones at depth 0", m.depth_0[0]);
    let mut elsewhere = 0usize;
    for &id in zones.iter().step_by(5) {
        let text = g.text_id(id);
        assert_eq!(
            g.first_sub_zone(id, 0),
            Ok(id),
            "first sub-zone of {text} at depth 0"
        );
        assert_eq!(
            g.sub_zone_at_index(id, 0, 0),
            Ok(id),
            "sub-zone 0 of {text} at depth 0"
        );
        assert_eq!(
            g.sub_zones(id, 0),
            Ok(vec![id]),
            "sub-zones of {text} at depth 0"
        );
        assert_eq!(
            g.count_sub_zones(id, 0),
            Ok(1),
            "count of {text} at depth 0"
        );
        assert!(
            engine_may_be_asked::<S>(id.0),
            "the engine cannot read {text} back"
        );
        assert_eq!(
            o::count_sub_zones(S::ORACLE, id.0, 0),
            1,
            "the engine's count of {text}"
        );
        assert_eq!(
            o::sub_zones(S::ORACLE, id.0, 0),
            vec![id.0],
            "the engine's sub-zones of {text} at depth 0"
        );
        let first = o::first_sub_zone(S::ORACLE, id.0, 0);
        if first != id.0 {
            assert!(
                engine_neighbours::<S>(id.0).contains(&first),
                "the engine's first sub-zone of {text} at depth 0, {}, is no neighbour of it",
                o::text_id(S::ORACLE, first)
            );
            elsewhere += 1;
        }
        depth_0.hit();
    }

    let mut finest = Floor::new("sub-zones at resolution 33", m.finest_sub_zones[0]);
    let mut engine_has_none = 0usize;
    // The zones at resolutions 31 and 32, whose sub-zones at depths 2 and 1 lie at 33; a zone
    // at 33 is its own only sub-zone at depth 0, where both sides answer 0 at one level, and
    // is not taken.
    let parents = zones
        .iter()
        .filter(|&&id| (31..=32).contains(&g.resolution(id)));
    for &id in parents {
        let text = g.text_id(id);
        let depth = 33 - g.resolution(id);
        let subs = g.sub_zones(id, depth).unwrap();
        for (i, &s) in subs.iter().enumerate() {
            let st = g.text_id(s);
            assert_eq!(
                g.sub_zone_index(id, s),
                Ok(Some(i as u64)),
                "index of {st}, at {i}, in {text}"
            );
            assert!(
                engine_may_be_asked::<S>(s.0),
                "the engine cannot read {st} back"
            );
            match o::sub_zone_index(S::ORACLE, id.0, s.0) {
                -1 => {
                    assert!(
                        !o::zone_has_sub_zone(S::ORACLE, id.0, s.0),
                        "the engine finds {st} in {text} and gives it no index"
                    );
                    engine_has_none += 1;
                }
                e => assert_eq!(e, i as i64, "the engine's index of {st} in {text}"),
            }
            finest.hit();
        }
    }

    let mut same_level = Floor::new("sub-zone index at one level", m.same_level);
    for &id in zones.iter().step_by(5) {
        let text = g.text_id(id);
        for n in g.neighbors(id) {
            let nt = g.text_id(n);
            assert_eq!(g.sub_zone_index(id, n), Ok(None), "index of {nt} in {text}");
            assert!(
                engine_may_be_asked::<S>(n.0),
                "the engine cannot read {nt} back"
            );
            assert_eq!(
                o::sub_zone_index(S::ORACLE, id.0, n.0),
                0,
                "the engine's index of {nt}, a neighbour, in {text}"
            );
            same_level.hit();
        }
    }
    eprintln!(
        "sub-zone departures: at depth 0 the engine's first sub-zone is a neighbour for {elsewhere} \
         of {} zones; at resolution 33 it has no index for {engine_has_none} of {} sub-zones",
        depth_0.n, finest.n
    );
    assert_eq!(
        [depth_0.n, elsewhere, finest.n, engine_has_none],
        [
            m.depth_0[0],
            m.depth_0[1],
            m.finest_sub_zones[0],
            m.finest_sub_zones[1]
        ],
        "the engine's departures at depth 0 and at resolution 33, and the sub-zones at 33 \
         compared (seed {SEED:#x})"
    );
}

/// A non-finite latitude or longitude is refused here, naming the coordinate at fault, at
/// every resolution. The engine answers it with the pentagon of root rhombus 0 at the
/// resolution asked for, one fixed, real cell centred on the icosahedron's vertex of base
/// cell 1, whatever the input: re-checked live at every case.
pub fn non_finite_coordinates_are_refused<S: Aperture3>() {
    let g = S::grid();
    let mut floor = Floor::new("non-finite coordinates", 7 * 34);
    let vertex = icosahedron_vertices()[1];
    for (lat, lon, name) in div::NON_FINITE_COORDINATES {
        for res in 0..=g.max_resolution() {
            let at = format!("({lat}, {lon}) at {res}");
            assert_eq!(
                g.zone_from_geo(lat, lon, res).err(),
                Some(Error::NonFinite { name }),
                "our answer at {at}"
            );
            let e = o::zone_from_geo(S::ORACLE, lat, lon, i32::from(res));
            assert_eq!(
                o::text_id(S::ORACLE, e),
                div::aperture_3_engine_answer_to_a_non_finite_coordinate(res),
                "the engine no longer answers as recorded at {at}"
            );
            let c = o::centroid(S::ORACLE, e);
            assert!(
                arc_deg(c, vertex) < TOL_DEG,
                "the engine's answer at {at} is centred at {c:?}, not at the vertex {vertex:?}"
            );
            floor.hit();
        }
    }
}
