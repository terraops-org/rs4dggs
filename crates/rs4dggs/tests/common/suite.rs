//! A grid checked against a live DGGAL, operation by operation.
//!
//! Every test here is generic over the grid under test, a `Subject`, and an oracle test
//! file runs the suite for one grid through `oracle_suite!`, which lists the tests that
//! apply to the grid's aperture. A test that holds on both apertures takes any subject; a
//! test written for one aperture binds its subject's topology and indexing, so that it
//! cannot be run on the other. The samples and their seeds are the same for every grid of
//! an aperture. What differs is the grid itself, the engine's class for it, the divergences
//! listed for it, its unresolved regions, and what has been measured of it against the
//! engine (see `measured`); and, between the apertures, what the suite takes from each
//! (see the helpers keyed on `Topology::APERTURE`).
//!
//! Nothing here skips an input. A missing engine fails the build, every test asserts a
//! floor on the number of cases it compared, and every failure message carries the seed
//! of the sample. Identifiers, text identifiers, neighbour lists and the hierarchy compare
//! exactly, and so do coordinates: every latitude and longitude this crate gives must be
//! the engine's own double, bit for bit, a claim made for x86-64 Linux with glibc, where
//! the suites run (see `assert_bit_identical`). Each geometric test still prints the
//! largest gap it measured, which is zero while the claim holds. A neighbour list compares
//! as a sequence with the engine's, which on aperture 7 is taken less the four kinds of
//! entry this crate drops, every one of them re-checked against the engine as it is left
//! out (see `engine_list`); nothing else is set aside, anywhere.
//!
//! A disagreement is allowed only in two ways. Either it is a characterised divergence,
//! where our answer is shown to be the defensible one and the evidence is re-checked
//! against the engine at every instance on every run (`common::expected_divergences`:
//! the witnesses listed for the grid, six for IGEO7, two each for IVEA7H and RTEA7H, and
//! four for each aperture-3 grid, the class rules of the grid's aperture, the input this
//! crate refuses where the engine answers it meaninglessly, in
//! `refusals_depart_from_the_engine_as_recorded`, and the sub-zone index of two different
//! zones at one level, in `sub_zones_at_depth_1_as_the_engine_lists_them`); or it lies inside one
//! of the grid's measured regions in `common::unresolved`, for each aperture-7 grid the
//! broken seams, which are open questions rather than rulings, and whose cases are
//! counted, printed and held under a ceiling. Anything else fails, with the seed.
//!
//! The engine is asked a relational question, about neighbours, parents, children,
//! ancestry, siblings or a sub-zone index, only about an identifier it reads back, as
//! `engine_may_be_asked` allows: on aperture 3 some of those questions take the test
//! process down, and no Rust mechanism can catch that. On aperture 7 `engine_can_read`
//! holds for the null zone too, which the engine reads back as `(null)` at level -1, so
//! the ordinary neighbour comparison asks it for the null zone's neighbours as it would
//! any other zone's, and it answers none, as this crate does. The null-zone tests ask
//! the engine directly for the same, and for the neighbours of a zone at the twentieth
//! level, without going through `engine_may_be_asked`; the engine answers none there
//! too.
//!
//! On aperture 7 the parents and the children are DGGAL's own, its geometric hierarchy, and
//! are compared with its `getZoneParents` and `getZoneChildren` directly, as sequences: the
//! engine's lists less the three kinds of entry this crate drops, every one of them
//! re-checked against the engine as it is left out (see `engine_relatives`). With them go
//! the centroid child, the centroid parent as the engine defines it, and the three
//! predicates against the engine's own, on pairs drawn from the lists (see `Hierarchy`).
//! Two samples are compared so: the zones that positions give, in
//! `parents_and_children_as_the_engine_lists_them`, and identifiers of a broken seam that
//! are built by text, in `the_hierarchy_at_seam_identifiers_built_by_text`. What the broken
//! seams do to the hierarchy is no divergence, for this crate answers there as the engine
//! does: it is counted by class and by level and held, exactly, to what is recorded of each
//! grid, and every zone so counted must lie within the band about the seams. The anomalies
//! of the parents, a zone with none and a first parent that the identifier does not name,
//! are met at the identifiers built by text alone, and each of those is asserted not to be
//! the zone found at its own centroid.
//!
//! The sub-zone order of the aperture-7 grids is DGGAL's own too, and is compared with its
//! `getSubZones` entry by entry: at depth 1 on a part of this file's sample, in
//! `sub_zones_at_depth_1_as_the_engine_lists_them`, with the refusals and the index of a
//! sub-zone; and to depth 5, with the first sub-zone, the zone at every index and what the
//! broken seams do to an order, in `common::sub_zone_order`.

use std::collections::{BTreeMap, BTreeSet, HashSet};

use super::expected_divergences::{self as div, DIVERGENCES, Evidence};
use super::subject::Subject;
use super::unresolved;
use super::*;
use dggal_oracle as o;
use rs4dggs::indexings::Z7;
use rs4dggs::topologies::HexA7;
use rs4dggs::{Error, Indexing, Topology, ZoneId};

pub(super) const SEED: u64 = 0x16E0_7000_0000_0001;

/// What has been measured of one grid against the engine. These are findings about the
/// grid, not properties of the samples: the ceilings under which admitted disagreements
/// are held, so that a class rule or a region loosened by a later edit shows as a jump,
/// and the phenomena a sample must reach, so that the rules which admit them are
/// exercised rather than merely described. The floors on compared cases describe the
/// samples, and stay with the tests.
#[derive(Debug)]
struct Measured {
    /// `zone_from_geo_text_id_and_resolution`: inputs at which both sides answer the null
    /// zone, at least.
    main_min_nulls_agreed: usize,
    /// The same: quantisations inside a region that neither answer contains, at most.
    main_max_unresolved: usize,
    /// The same: ties admitted, at most.
    main_max_ties: usize,
    /// The same: unreadable answers admitted, at most.
    main_max_unreadable: usize,
    /// The same: ties in which our answer is the null zone, exactly.
    main_ties_ours_null: usize,
    /// `seams_between_faces_quantise_like_dggal`: inputs at which both sides answer the
    /// null zone, at least.
    seam_min_nulls_agreed: usize,
    /// The same: resolutions at which the unreadable-answer rule must admit something.
    seam_unreadable_at: &'static [u8],
    /// The same: quantisations inside a region that neither answer contains, exactly.
    seam_unresolved: usize,
    /// The same: ties admitted, at most.
    seam_max_ties: usize,
    /// The same: unreadable answers admitted, at most.
    seam_max_unreadable: usize,
    /// The same: ties in which our answer is the null zone, at most.
    seam_max_ties_ours_null: usize,
    /// `listed_divergences_still_hold_with_their_evidence`: entries listed for the grid
    /// in `DIVERGENCES`, exactly. It is a constant of the grid and not the length of the
    /// table, which that test would otherwise both take its floor from and check: with
    /// the length, deleting every entry for a grid would leave a floor of nothing and
    /// pass on no iterations.
    listed_divergences: usize,
    /// `quantize_matches_dggal_at_every_resolution`: null answers agreed, at least.
    topology_min_nulls_agreed: usize,
    /// The same: listed divergences met, exactly.
    topology_carve_outs: usize,
}

/// IGEO7, as measured when this suite was written.
static IGEO7_MEASURED: Measured = Measured {
    main_min_nulls_agreed: 10,
    // Measured: none. The one case this line carried, at (89.99997981407529,
    // 101.19999974473858) at resolution 16, 2e-5 degrees from the north pole, in which
    // neither side's answer contained the input, went when the projection was put into the
    // order the compiled engine performs it: this crate's planar point there is now the
    // engine's own, and the two sides name the same cell. The ceiling follows the
    // observation down to nothing.
    main_max_unresolved: 0,
    // Measured: no ties and 3 unreadable answers. The ties fell from 29 to none when the
    // projection was put into the order the compiled engine performs it: every one of them
    // was an input on a cell boundary, the quantiser was already the engine's own, and so
    // the answers could part only on the planar point, which is now the engine's own too;
    // both sides name the same cell at all 29. The ceiling follows the observation down
    // to nothing; a tie now fails here.
    main_max_ties: 0,
    main_max_unreadable: 4,
    // Measured: none. The one inferior tie this line carried, at (0, 11.2) at resolution
    // 1, is gone. Its cause was never the side test but the polygon the test was given:
    // this crate's vertices for `041` and `036` differed from `getZoneCRSVertices` by one
    // and two units in the last place. The vertex generation has since been put into the
    // order the compiled engine performs it, and every vertex this crate computes in the
    // 5x6 CRS is now bit-identical to the engine's over a census of 16,163 zones, so the
    // two sides answer over the same polygon and the tie is an ordinary one. The bound
    // follows the observation down to nothing.
    main_ties_ours_null: 0,
    // The sample must reach the strip and the pole, or the class rules go untested.
    seam_min_nulls_agreed: 50,
    seam_unreadable_at: &[16, 18],
    seam_unresolved: 0,
    // Measured: no ties and 46 unreadable answers. The ties fell from 91 to 87 when the
    // side test was put into the compiled engine's order, from 87 to 84 when the geometry
    // followed, and from 84 to none when the projection did: the 84 lay at nine stations
    // on the mirror planes, and at every one the two sides' planar points had differed in
    // their last bits, where they are now the same. The ceiling follows the observation
    // down to nothing.
    seam_max_ties: 0,
    seam_max_unreadable: 51,
    // The ties where our answer is the inferior one, null: measured none. The one that
    // remained, at the south pole at resolution 15, went when the vertex generation was
    // put into the compiled engine's order, as the two at the midpoint of the edge from
    // base cell 6 to 10 had gone with the side test before it. The ceiling follows the
    // observation down to nothing. See `boundary_tie_on_a_mirror_plane`.
    seam_max_ties_ours_null: 0,
    // The six entries of `DIVERGENCES` that name IGEO7. There were eight: the two tie
    // witnesses at resolution 14 no longer diverge, both sides naming the engine's cell
    // there once the projection was ported.
    listed_divergences: 6,
    // `NULL_PROBE` gives the null zone at resolutions 15, 17 and 19.
    topology_min_nulls_agreed: 3,
    // The resolution-18 divergence at `NULL_PROBE`, and nothing else.
    topology_carve_outs: 1,
};

/// IVEA7H, as measured against its own engine. Every disagreement below was
/// characterised on this grid before its count was entered: the unreadable answers lie
/// in the broken seams along the edges from base cell 0 to 1 and from 1 to 6, as on
/// IGEO7, and the ties lie on the mirror planes and at the axis points, where the grid's
/// symmetries were measured to hold as on IGEO7.
/// Ceilings allow a tenth more than was measured, as IGEO7's do.
static IVEA7H_MEASURED: Measured = Measured {
    // Measured: 10, at resolutions 15, 17 and 19 at the two points listed for unreadable
    // answers, and at 17 and 19 at the two adversarial points a hair from the north pole.
    main_min_nulls_agreed: 5,
    // Measured: none. The phenomenon exists in the strip, between 1.8e-5 and 2.5e-5
    // degrees from the north pole at resolution 16, but no point of the main sample meets
    // it: at IGEO7's own witness point the two sides of IVEA7H name one cell.
    main_max_unresolved: 0,
    // Measured: no ties and 3 unreadable answers. The 15 ties, at (0, 11.2) and at
    // (10, 11.2), went when the projection was put into the order the compiled engine
    // performs it, as IGEO7's did. The ceiling follows the observation down to nothing.
    main_max_ties: 0,
    main_max_unreadable: 4,
    // Measured: none. The two that stood here, (0, -168.8) at resolutions 9 and 19, are
    // gone: the face the eC takes carries that input to the planar point of the seam
    // station, (4.5, 5.5), exactly, and the side test in the order the compiled engine
    // performs it returns an exact zero there and accepts the engine's cell.
    main_ties_ours_null: 0,
    // Measured: 135 null answers agreed, two more than the 133 measured when this table
    // was first filled. The engine's answers have not moved, so each rise is on this
    // side: a seam station at which the engine answers its null zone, and this crate
    // named a cell, now answers the null zone as well. The first came with the side test;
    // the second, at the station on the north pole at resolution 1, (89.99999999999991,
    // 6.536633648797066), came with the projection.
    seam_min_nulls_agreed: 60,
    seam_unreadable_at: &[16, 18],
    seam_unresolved: 0,
    // Measured: no ties and 41 unreadable answers, 11 at 16 and 30 at 18. The ties fell
    // from 68 to 61 when the side test was put into the compiled engine's order, from 61
    // to 58 when the geometry followed, and from 58 to none when the projection did: the
    // 58 lay at seven stations, among them the north pole, carried to face 0 a few units
    // in the last place from the engine's planar point, and the station on the equator at
    // longitude 101.19999994270422, carried by its sub-triangle to a cell boundary in the
    // same way. The ceiling follows the observation down to nothing.
    seam_max_ties: 0,
    seam_max_unreadable: 46,
    // Measured: none. The one that stood here, the south pole at resolution 15, went when
    // the vertex generation was put into the compiled engine's order, as the two at the
    // midpoint of the edge from base cell 6 to 10 had gone with the side test before it.
    // The ceiling follows the observation down to nothing.
    seam_max_ties_ours_null: 0,
    // The two entries of `DIVERGENCES` that name IVEA7H. There were three: the tie
    // witness at (0, 11.2) at resolution 0 no longer diverges.
    listed_divergences: 2,
    // `NULL_PROBE` gives the null zone at resolutions 15, 17 and 19, as on IGEO7.
    topology_min_nulls_agreed: 3,
    // The resolution-18 divergence at `NULL_PROBE`, listed for IVEA7H, and nothing else.
    topology_carve_outs: 1,
};

/// RTEA7H, as measured against its own engine. Every disagreement below was
/// characterised on this grid before its count was entered: the unreadable answers lie
/// in the broken seams along the edges from base cell 0 to 1 and from 1 to 6, as on
/// IGEO7, and the ties lie on the mirror planes and at the axis points, where the grid's
/// symmetries were measured to hold as on IGEO7.
/// Ceilings allow a tenth more than was measured, as IGEO7's do.
static RTEA7H_MEASURED: Measured = Measured {
    // Measured: 10, at resolutions 15, 17 and 19 at the two points listed for unreadable
    // answers, and at 17 and 19 at the two adversarial points a hair from the north pole.
    main_min_nulls_agreed: 5,
    // Measured: none.
    main_max_unresolved: 0,
    // Measured: no ties and 3 unreadable answers. The 13 ties, at the two poles, at
    // (0, 11.2) and at three other points on or beside the mirror planes, went when the
    // projection was put into the order the compiled engine performs it, as IGEO7's did;
    // at the north pole at resolution 19 this crate's planar point had lain nine units in
    // the last place from the engine's. The ceiling follows the observation down to
    // nothing.
    main_max_ties: 0,
    main_max_unreadable: 4,
    // Measured: none. The two that stood here, the south pole at resolution 15 and the
    // north pole at 11, were decided not by the side test but by vertices a unit or two
    // in the last place from the engine's. The vertex generation has since been put into
    // the order the compiled engine performs it, and both are gone. See
    // `boundary_tie_on_a_mirror_plane`.
    main_ties_ours_null: 0,
    // Measured: 131 null answers agreed, one more than the 130 measured when this table
    // was first filled, and for the reason given for IVEA7H: the engine's answers have
    // not moved, so a seam station at which it answers its null zone, and this crate
    // named a cell, now answers the null zone as well.
    seam_min_nulls_agreed: 65,
    seam_unreadable_at: &[16, 18],
    seam_unresolved: 0,
    // Measured: no ties and 37 unreadable answers, 11 at 16 and 26 at 18. The ties did
    // not move when the side test was put into the compiled engine's order, fell from 94
    // to 92 when the geometry followed, and from 92 to none when the projection did: the
    // 92 lay at ten stations on the mirror planes and the half-turn axes, among them one
    // 8.8e-14 degrees from the midpoint of the edge from base cell 6 to 10, which tied at
    // all twenty resolutions. The ceiling follows the observation down to nothing.
    seam_max_ties: 0,
    seam_max_unreadable: 41,
    // Measured: none. The south pole at resolution 15 and the station on the north pole
    // at 11 both went when the vertex generation was put into the compiled engine's
    // order, as the two at the midpoint of the edge from base cell 6 to 10 had gone with
    // the side test before them. The ceiling follows the observation down to nothing.
    seam_max_ties_ours_null: 0,
    // The two entries of `DIVERGENCES` that name RTEA7H. There were three: the tie
    // witness at (0, 11.2) at resolution 3 no longer diverges.
    listed_divergences: 2,
    // `NULL_PROBE` gives the null zone at resolutions 15, 17 and 19, as on IGEO7.
    topology_min_nulls_agreed: 3,
    // The resolution-18 divergence at `NULL_PROBE`, listed for RTEA7H, and nothing else.
    topology_carve_outs: 1,
};

/// ISEA3H, as measured against its own engine. Its one class of disagreement is the
/// identifiers the engine answers along two seams and cannot read back, where this crate
/// answers the null zone (see `expected_divergences`); it has no unresolved region, and its
/// neighbour lists, disks and hierarchy agree with the engine's wherever the engine can be
/// asked. Every value is the observation itself, the counts being fixed by the samples.
static ISEA3H_MEASURED: Measured = Measured {
    // Measured: 12, at resolutions 0 and 1, at the point a hair from the south pole, at the
    // two a hair from the north pole on the meridian 101.2, and at the three listed witnesses
    // that are not adversarial inputs, where the engine answers its null zone on the seams.
    main_min_nulls_agreed: 12,
    main_max_unresolved: 0,
    // The tie rule is not offered on aperture 3.
    main_max_ties: 0,
    // Measured: 48, at resolutions 2 to 20, from the two points a hair from the north pole on
    // the meridian 101.2 and from the listed witnesses at the resolutions other than their
    // own.
    main_max_unreadable: 48,
    main_ties_ours_null: 0,
    // Measured: 624.
    seam_min_nulls_agreed: 624,
    // Measured: the rule admits answers at every resolution from 2 to 13, along the seam over
    // the north pole at all of them and along the southern one at 2 to 5, at the seam
    // stations and at 1e-9 radians either side of them.
    seam_unreadable_at: &[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
    seam_unresolved: 0,
    seam_max_ties: 0,
    // Measured: 182.
    seam_max_unreadable: 182,
    seam_max_ties_ours_null: 0,
    // The four entries of `DIVERGENCES` that name ISEA3H.
    listed_divergences: 4,
    // `APERTURE_3_NULL_PROBE` answers the null zone at resolutions 0 and 1 on both sides.
    topology_min_nulls_agreed: 2,
    topology_carve_outs: 0,
};

/// IVEA3H, as measured against its own engine, on its own samples: the main sample and the
/// seam sample quantised under its own projection, and its own four witnesses. Its one class
/// of disagreement is ISEA3H's, the identifiers the engine answers along two seams and
/// cannot read back, where this crate answers the null zone; it has no unresolved region,
/// and its neighbour lists, disks and hierarchy agree with the engine's wherever the engine
/// can be asked. Every value is the observation itself.
static IVEA3H_MEASURED: Measured = Measured {
    // Measured: 12, at resolutions 0 and 1, at the point a hair from the south pole, at the
    // two a hair from the north pole on the meridian 101.2, and at the three listed witnesses
    // that are not adversarial inputs, where the engine answers its null zone on the seams.
    main_min_nulls_agreed: 12,
    main_max_unresolved: 0,
    // The tie rule is not offered on aperture 3.
    main_max_ties: 0,
    // Measured: 44, at resolutions 2 to 20: from the two points a hair from the north pole on
    // the meridian 101.2, at 3 from the one that is the witness at 2 and at 2 to 5 from the
    // other; from the southern witness of resolution 2, at 3; and from the two witnesses of
    // resolution 21, at 2 to 20.
    main_max_unreadable: 44,
    main_ties_ours_null: 0,
    // Measured: 628, at resolutions 0 to 5.
    seam_min_nulls_agreed: 628,
    // Measured: as on ISEA3H, the rule admits answers at every resolution from 2 to 13, along
    // the seam over the north pole at all of them and along the southern one at 2 to 5.
    seam_unreadable_at: &[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
    seam_unresolved: 0,
    seam_max_ties: 0,
    // Measured: 182.
    seam_max_unreadable: 182,
    seam_max_ties_ours_null: 0,
    // The four entries of `DIVERGENCES` that name IVEA3H.
    listed_divergences: 4,
    // `APERTURE_3_NULL_PROBE` answers the null zone at resolutions 0 and 1 on both sides.
    topology_min_nulls_agreed: 2,
    topology_carve_outs: 0,
};

/// RTEA3H, as measured against its own engine, on its own samples: the main sample and the
/// seam sample quantised under its own projection, and its own four witnesses. Its one class
/// of disagreement is ISEA3H's, the identifiers the engine answers along two seams and
/// cannot read back, where this crate answers the null zone; it has no unresolved region,
/// and its neighbour lists, disks and hierarchy agree with the engine's wherever the engine
/// can be asked. Every value is the observation itself.
static RTEA3H_MEASURED: Measured = Measured {
    // Measured: 12, at resolutions 0 and 1, at the point a hair from the south pole, at the
    // two a hair from the north pole on the meridian 101.2, and at the three listed witnesses
    // that are not adversarial inputs, where the engine answers its null zone on the seams.
    main_min_nulls_agreed: 12,
    main_max_unresolved: 0,
    // The tie rule is not offered on aperture 3.
    main_max_ties: 0,
    // Measured: 44, at resolutions 2 to 20: from the two points a hair from the north pole on
    // the meridian 101.2, at 3 from the one that is the witness at 2 and at 2 to 5 from the
    // other; from the southern witness of resolution 2, at 3; and from the two witnesses of
    // resolution 21, at 2 to 20.
    main_max_unreadable: 44,
    main_ties_ours_null: 0,
    // Measured: 636, at resolutions 0 to 5.
    seam_min_nulls_agreed: 636,
    // Measured: as on ISEA3H, the rule admits answers at every resolution from 2 to 13, along
    // the seam over the north pole at all of them and along the southern one at 2 to 5.
    seam_unreadable_at: &[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
    seam_unresolved: 0,
    seam_max_ties: 0,
    // Measured: 182.
    seam_max_unreadable: 182,
    seam_max_ties_ours_null: 0,
    // The four entries of `DIVERGENCES` that name RTEA3H.
    listed_divergences: 4,
    // `APERTURE_3_NULL_PROBE` answers the null zone at resolutions 0 and 1 on both sides.
    topology_min_nulls_agreed: 2,
    topology_carve_outs: 0,
};

/// What has been measured of the grid under test. A grid of which nothing has been
/// measured fails here, once its sweep is done, rather than borrowing the findings of
/// another grid: a ceiling taken from IGEO7 would be a meaningless baseline for another
/// grid, against which a loosened rule would not show, and a required phenomenon could
/// demand one the grid does not have. Every test that consults this first prints what it
/// observed, with [`print_observed`], so that a failure here carries the values the
/// table is to be filled from.
fn measured<S: Subject>() -> &'static Measured {
    match S::NAME {
        "IGEO7" => &IGEO7_MEASURED,
        "IVEA7H" => &IVEA7H_MEASURED,
        "RTEA7H" => &RTEA7H_MEASURED,
        "ISEA3H" => &ISEA3H_MEASURED,
        "IVEA3H" => &IVEA3H_MEASURED,
        "RTEA3H" => &RTEA3H_MEASURED,
        other => panic!(
            "nothing has been measured of {other} against the engine, so its ceilings and \
             required phenomena are unknown. The line \"observed for {other}: ...\" printed \
             above gives the values this run met, field by field; the table is to be filled \
             from them only once each has been characterised, and never from another grid"
        ),
    }
}

/// Prints, under the names of the [`Measured`] fields a test consults, the values the
/// test observed, with the number of cases it compared. A test calls this before it
/// consults [`measured`], so that the values are on record whether or not the grid has
/// been measured, and whether or not they are within the table.
fn print_observed<S: Subject>(fields: &[(&str, &dyn std::fmt::Debug)], cases: usize) {
    let fields: Vec<String> = fields.iter().map(|(k, v)| format!("{k} = {v:?}")).collect();
    eprintln!(
        "observed for {}: {}; compared cases = {cases}",
        S::NAME,
        fields.join(", ")
    );
}

/// Drops repeated points, keeping the first, so that a point listed twice (an
/// adversarial input that is also a divergence witness) is compared once.
fn distinct(points: Vec<(f64, f64)>) -> Vec<(f64, f64)> {
    let mut seen = HashSet::new();
    points
        .into_iter()
        .filter(|p| seen.insert((p.0.to_bits(), p.1.to_bits())))
        .collect()
}

/// The main sample: the adversarial inputs, every divergence witness, and 400 seeded
/// uniform points.
fn inputs<S: Subject>() -> Vec<(f64, f64)> {
    let mut v = adversarial_points();
    v.extend(
        DIVERGENCES
            .iter()
            .filter(|d| d.grid == S::NAME)
            .map(|d| (d.lat, d.lon)),
    );
    v.extend(sphere_points(SEED, 400));
    distinct(v)
}

/// Resolutions at which the main sample's zones are taken for the geometric, neighbour
/// and hierarchy tests. On aperture 7: the coarse end, a spread through the middle, and
/// every resolution from 13, where the odd-level quantisation and the strip begin. On
/// aperture 3: the coarse end, where the engine answers its null zone near the poles;
/// both parities through the middle; 20, 21 and 22, the last two resolutions at which it
/// answers an identifier it cannot read back and the first at which it does not; and 31,
/// 32 and 33, the finest three, whose sub-zones lie at 33 at depths 2, 1 and 0, and whose
/// last level has no children here.
fn zone_resolutions<S: Subject>() -> &'static [u8] {
    match S::T::APERTURE {
        7 => &[0, 1, 2, 3, 5, 8, 10, 11, 13, 14, 15, 16, 17, 18, 19],
        3 => &[0, 1, 2, 3, 4, 5, 8, 11, 14, 17, 20, 21, 22, 27, 31, 32, 33],
        a => panic!("no zone resolutions are recorded for aperture {a}"),
    }
}

/// The texts of the zones added to the main sample's: every base cell, and a pentagon at
/// resolution 4 under each. On aperture 3 the twelve pentagons of a level are the ten root
/// rhombi's and the two polar ones, roots 10 and 11, at index 0.
fn fixed_zone_texts<S: Subject>() -> Vec<String> {
    match S::T::APERTURE {
        7 => (0..12)
            .flat_map(|base| [format!("{base:02}"), format!("{base:02}0000")])
            .collect(),
        3 => (0..12)
            .flat_map(|root| [format!("A{root:X}-0-A"), format!("C{root:X}-0-A")])
            .collect(),
        a => panic!("no fixed zones are recorded for aperture {a}"),
    }
}

/// Zones to check geometry, topology and hierarchy on: the main sample at the resolutions
/// of [`zone_resolutions`], and the zones of [`fixed_zone_texts`].
pub(super) fn zones<S: Subject>() -> Vec<ZoneId> {
    let mut ids = Vec::new();
    for (lat, lon) in inputs::<S>() {
        for &res in zone_resolutions::<S>() {
            ids.push(S::grid().zone_from_geo(lat, lon, res).unwrap().id());
        }
    }
    for text in fixed_zone_texts::<S>() {
        ids.push(S::grid().zone_from_text(&text).unwrap().id());
    }
    ids.sort();
    ids.dedup();
    ids
}

/// The class rules that can admit a disagreement, by name for the tallies.
const UNREADABLE: &str = "engine cannot read its own answer";
const TIE: &str = "boundary tie on a mirror plane";
/// A tie admitted by the same rule in which our answer is the null zone: the one kind
/// of tie in which our answer is the inferior one, counted apart so that a new instance
/// shows.
const TIE_OURS_NULL: &str = "boundary tie on a mirror plane, our answer null (inferior)";

/// The signature of a class rule in `common::expected_divergences`.
type ClassRule = fn(f64, f64, u8, u64, u64) -> Result<(), String>;

/// The class rules that may admit a quantisation disagreement on the grid under test, in
/// the order they are tried. Each is a finding about its aperture: the mirror planes and
/// half-turn axes of the tie rule, and the strip of the rule for unreadable answers, were
/// measured on the aperture-7 grids. On aperture 3 the one rule is that for the
/// unreadable identifiers along two seams; no tie has been met there, and the symmetries
/// the tie rule rests on have not been measured, so it is not offered.
fn class_rules<S: Subject>() -> Vec<(&'static str, ClassRule)> {
    match S::T::APERTURE {
        7 => vec![
            (UNREADABLE, div::engine_cannot_read_its_own_answer::<S>),
            (TIE, div::boundary_tie_on_a_mirror_plane::<S>),
        ],
        3 => vec![(
            UNREADABLE,
            div::engine_answers_an_unreadable_polar_root::<S>,
        )],
        a => panic!("no class rules are recorded for aperture {a}"),
    }
}

/// Whether this crate's geometry for `id` departs from the engine's by decision: only the
/// null zone on aperture 3, which this crate answers with the point (0, 0) and six zero
/// vertices, as it does on aperture 7, where the engine does too, and to which the
/// aperture-3 engine gives a centroid and vertices computed from the bits of its
/// identifier. The geometric tests re-check the engine's side of it at every instance, by
/// [`div::engine_null_zone_geometry`].
fn null_geometry_departs<S: Subject>(id: ZoneId) -> bool {
    id == ZoneId::NULL && S::T::APERTURE == 3
}

/// The level the engine reports for its null zone, which this crate reports as resolution
/// 0 on every grid (see `Grid::resolution`): -1 on aperture 7, where the eC's `getLevel`
/// answers so for it, and 63 on aperture 3, where the eC reads the level out of the null
/// zone's bits as out of any other identifier's (`RI3H.ec:878-881`).
fn engine_null_level<S: Subject>() -> i32 {
    match S::T::APERTURE {
        7 => -1,
        3 => 63,
        a => panic!("no null-zone level is recorded for aperture {a}"),
    }
}

/// Whether the engine may be asked a relational question about `z`: its neighbours, its
/// parents or children, an ancestry or sibling predicate, or a sub-zone index. Every such
/// call in this suite is made only when this holds.
///
/// It is `engine_can_read`, on both apertures. On aperture 3 the engine dies with `SIGFPE`
/// asking for the parents of the identifiers it cannot read back and of level-34
/// identifiers, writes seven neighbours into six-slot arrays for the unreadable ones, and
/// never returns from a sub-zone index with a sub-zone above resolution 33, and no Rust
/// mechanism can catch any of it. On aperture 7 the engine answers such questions, but an
/// identifier it cannot read back is no zone this crate names, and nothing compared here
/// needs its answers: the engine's walk of a disk is composed from its lists less the
/// entries this crate drops, which include every identifier it cannot read back.
pub(super) fn engine_may_be_asked<S: Subject>(z: u64) -> bool {
    o::engine_can_read(S::ORACLE, z)
}

/// The engine's neighbours of `z`, asked only where [`engine_may_be_asked`] allows it, and
/// otherwise a failure naming the identifier: an identifier the engine lists and cannot
/// itself be asked about is a finding to characterise, not a sample to skip.
pub(super) fn engine_neighbours<S: Subject>(z: u64) -> Vec<u64> {
    assert!(
        engine_may_be_asked::<S>(z),
        "the engine cannot read {} ({z:#018x}) back, and asking it for its neighbours would \
         bring the test process down (seed {SEED:#x})",
        o::text_id(S::ORACLE, z)
    );
    o::neighbors(S::ORACLE, z)
}

/// The four kinds of entry in the engine's neighbour lists that the aperture-7 grids
/// drop, in the order in which [`engine_list`] tries them, for the tallies.
const DROPPED: [&str; 4] = [
    "the engine's null zone",
    "the zone itself",
    "a repeat of an entry already listed",
    "an identifier the engine cannot read back",
];

/// The engine's neighbour list for `z` as this crate must give it, and how many entries of
/// each kind of [`DROPPED`] it left out, asked only where [`engine_may_be_asked`] allows.
///
/// On aperture 7 it is the engine's list, in its order, less its null zone, the zone
/// itself, a repeat of an entry already kept or left out, and an identifier the engine
/// cannot read back, each entry classified live by the first of those four kinds that fits
/// it, so that every entry left out is re-checked against the engine as one of the four on
/// every run. On aperture 3 it is the engine's list verbatim, which this crate gives
/// entry for entry.
fn engine_list<S: Subject>(z: u64) -> (Vec<u64>, [usize; 4]) {
    let list = engine_neighbours::<S>(z);
    if S::T::APERTURE != 7 {
        return (list, [0; 4]);
    }
    let mut dropped = [0usize; 4];
    let mut seen = HashSet::new();
    let mut kept = Vec::new();
    for t in list {
        let kind = if t == o::NULL_ZONE {
            0
        } else if t == z {
            1
        } else if !seen.insert(t) {
            2
        } else if !o::engine_can_read(S::ORACLE, t) {
            3
        } else {
            kept.push(t);
            continue;
        };
        dropped[kind] += 1;
    }
    (kept, dropped)
}

/// How the disagreements of one quantisation sweep were accounted for.
#[derive(Default)]
struct Accounted {
    /// Hits per listed divergence, by its index in `DIVERGENCES`.
    listed: BTreeMap<usize, usize>,
    /// Disagreements admitted by a class rule, by rule and resolution.
    by_class: BTreeMap<(&'static str, u8), usize>,
    /// Inputs at which both sides answered the null zone.
    nulls_agreed: usize,
    /// Disagreements inside the broken seams where neither answer contains the input.
    unresolved: usize,
}

impl Accounted {
    /// Compares one quantisation. A disagreement must be a listed divergence, whose
    /// answers are checked here, or satisfy one of the two class rules, whose evidence
    /// is checked live, or be the unresolved case in which neither answer contains the
    /// input; anything else fails with the seed.
    fn quantisation<S: Subject>(
        &mut self,
        sample: &str,
        lat: f64,
        lon: f64,
        res: u8,
        ours: ZoneId,
        theirs: u64,
    ) {
        if ours.0 == theirs {
            if theirs == o::NULL_ZONE {
                self.nulls_agreed += 1;
            }
            return;
        }
        if let Some(d) = div::lookup(S::NAME, lat, lon, res) {
            assert_eq!(
                (
                    S::grid().text_id(ours).as_str(),
                    o::text_id(S::ORACLE, theirs).as_str()
                ),
                (d.ours, d.dggal),
                "listed divergence at ({lat}, {lon}, {res}) now answers differently (seed {SEED:#x})"
            );
            let i = DIVERGENCES.iter().position(|e| std::ptr::eq(e, d)).unwrap();
            *self.listed.entry(i).or_default() += 1;
            return;
        }
        // The class rules, in turn; the first whose evidence holds admits it.
        let mut refusals = Vec::new();
        for (name, rule) in class_rules::<S>() {
            match rule(lat, lon, res, ours.0, theirs) {
                Ok(()) => {
                    let name = if name == TIE && ours == ZoneId::NULL {
                        eprintln!(
                            "zone_from_geo({lat}, {lon}, {res}): a tie where we answer null \
                             and DGGAL {}",
                            o::text_id(S::ORACLE, theirs)
                        );
                        TIE_OURS_NULL
                    } else {
                        name
                    };
                    *self.by_class.entry((name, res)).or_default() += 1;
                    return;
                }
                Err(why) => refusals.push(format!("{name}: {why}")),
            }
        }
        // Not a divergence, but possibly the one measured unresolved case.
        match unresolved::neither_answer_contains_the_input::<S>(lat, lon, res, ours.0, theirs) {
            Ok(()) => {
                self.unresolved += 1;
                eprintln!(
                    "zone_from_geo({lat}, {lon}, {res}): unresolved, neither answer contains the input \
                     (ours {}, DGGAL {})",
                    S::grid().text_id(ours),
                    o::text_id(S::ORACLE, theirs)
                );
                return;
            }
            Err(why) => refusals.push(format!("unresolved, neither contains it: {why}")),
        }
        panic!(
            "zone_from_geo({lat}, {lon}, {res}) in the {sample} sample: ours {}, DGGAL {}, \
             and no characterised divergence accounts for it:\n  {}\n(seed {SEED:#x})",
            S::grid().text_id(ours),
            o::text_id(S::ORACLE, theirs),
            refusals.join("\n  ")
        );
    }

    fn class_count(&self, rule: &str) -> usize {
        self.by_class
            .iter()
            .filter(|((r, _), _)| *r == rule)
            .map(|(_, n)| n)
            .sum()
    }
}

/// The text identifier, its reading back on both sides, and the resolution, for one
/// identifier. The null zone's resolution is the one answer this port deliberately gives
/// differently (0 against DGGAL's -1 on aperture 7 and 63 on aperture 3, documented at
/// `Grid::resolution`), so it is pinned to both values rather than compared.
fn check_identifier<S: Subject>(id: ZoneId, res: u8) {
    let text = S::grid().text_id(id);
    assert_eq!(
        text,
        o::text_id(S::ORACLE, id.0),
        "text of {:#018x} (seed {SEED:#x})",
        id.0
    );
    assert_eq!(
        S::grid().zone_from_text(&text).map(|z| z.id()),
        Ok(id),
        "reading back {text} (seed {SEED:#x})"
    );
    assert_eq!(
        o::zone_from_text(S::ORACLE, &text),
        id.0,
        "DGGAL reading back {text} (seed {SEED:#x})"
    );
    let level = o::level(S::ORACLE, id.0);
    if id == ZoneId::NULL {
        assert_eq!(
            (S::grid().resolution(id), level),
            (0, engine_null_level::<S>()),
            "the null zone's resolution"
        );
    } else {
        assert_eq!(
            i32::from(S::grid().resolution(id)),
            level,
            "level of {text} (seed {SEED:#x})"
        );
        assert_eq!(
            S::grid().resolution(id),
            res,
            "resolution of {text} (seed {SEED:#x})"
        );
    }
}

pub fn max_resolution_matches<S: Subject>() {
    assert_eq!(
        i32::from(S::grid().max_resolution()),
        o::max_level(S::ORACLE)
    );
}

/// The suite takes the icosahedron's vertices from IGEO7's engine for every aperture-7
/// grid, and with them the thirty seams it samples and the strips of the unresolved
/// regions. That is sound only if the three grids place their base cells on the same
/// twelve points, as their shared orientation implies; this checks it against the
/// engines themselves, base cell by base cell.
///
/// The comparison is bit for bit, not within [`TOL_DEG`]: the three engines were
/// measured to give the same doubles, so a tolerance here would admit a difference the
/// grids do not have. It runs in each grid's suite rather than in one of them, because
/// `seam_points` and every region's `near` rest on these twelve points on every grid,
/// and a suite should check what it rests on.
pub fn the_three_engines_agree_on_the_icosahedron_vertices<S: Subject<T = HexA7, I = Z7>>() {
    use super::subject::{Igeo7Subject, Ivea7hSubject, Rtea7hSubject};
    let mut floor = Floor::new("icosahedron vertices", 24);
    let v = icosahedron_vertices();
    for oracle in [Ivea7hSubject::ORACLE, Rtea7hSubject::ORACLE] {
        for (base, &vertex) in v.iter().enumerate() {
            let z = o::zone_from_text(oracle, &format!("{base:02}"));
            assert_ne!(z, o::NULL_ZONE, "{oracle} reads no base cell {base:02}");
            let theirs = o::centroid(oracle, z);
            assert_eq!(
                (theirs.0.to_bits(), theirs.1.to_bits()),
                (vertex.0.to_bits(), vertex.1.to_bits()),
                "base cell {base:02} of {oracle} is {theirs:?}, where {} gives {vertex:?}",
                Igeo7Subject::ORACLE
            );
            floor.hit();
        }
    }
    eprintln!(
        "icosahedron vertices: the three engines agree bit for bit, from the {} suite",
        S::NAME
    );
}

/// The floor of the main sample, which differs by grid: the part of the sample that is
/// the same for every grid, the 45 adversarial inputs and the 400 seeded points at all
/// twenty resolutions, 8,900 quantisations, and twenty more for each point the grid's
/// listed witnesses add to it, that is each distinct witness point that is not already
/// an adversarial input. For IGEO7 that is six points and 9,020 quantisations, and for
/// IVEA7H and RTEA7H two points each and 8,940. On aperture 3 every point is taken at
/// thirty-four resolutions, so that the shared part is 15,130 quantisations, and each of
/// ISEA3H, IVEA3H and RTEA3H lists four witnesses, of which one is an adversarial input:
/// three points more and 15,232. The whole sample, so that the loss of a single point
/// fails.
///
/// A constant for each grid, and not a count taken from `DIVERGENCES`: the floor guards
/// the witnesses that table supplies to the sample, and read from the table it would fall
/// in step with any witness deleted from it, and so never fail on that deletion.
fn main_sample_floor<S: Subject>() -> usize {
    match S::NAME {
        "IGEO7" => 9_020,
        "IVEA7H" | "RTEA7H" => 8_940,
        "ISEA3H" | "IVEA3H" | "RTEA3H" => 15_232,
        other => panic!(
            "no main-sample floor is recorded for {other}; count its sample, the 8,900 \
             quantisations every grid shares and twenty for each witness point its listed \
             divergences add, and record it here"
        ),
    }
}

pub fn zone_from_geo_text_id_and_resolution<S: Subject>() {
    // See `main_sample_floor`; each listed witness must also be met, which is asserted
    // below.
    let mut floor = Floor::new("zone_from_geo", main_sample_floor::<S>());
    let mut acc = Accounted::default();
    for (lat, lon) in inputs::<S>() {
        for res in 0..=S::grid().max_resolution() {
            let ours = S::grid().zone_from_geo(lat, lon, res).unwrap_or_else(|e| {
                panic!("zone_from_geo({lat}, {lon}, {res}): {e} (seed {SEED:#x})")
            });
            let theirs = o::zone_from_geo(S::ORACLE, lat, lon, i32::from(res));
            acc.quantisation::<S>("main", lat, lon, res, ours.id(), theirs);
            check_identifier::<S>(ours.id(), res);
            floor.hit();
        }
    }
    eprintln!(
        "zone_from_geo: {} nulls agreed, {} listed divergences met, class-rule divergences by rule and resolution {:?}",
        acc.nulls_agreed,
        acc.listed.values().sum::<usize>(),
        acc.by_class
    );
    print_observed::<S>(
        &[
            ("main_min_nulls_agreed", &acc.nulls_agreed),
            ("main_max_unresolved", &acc.unresolved),
            ("main_max_ties", &acc.class_count(TIE)),
            ("main_max_unreadable", &acc.class_count(UNREADABLE)),
            ("main_ties_ours_null", &acc.class_count(TIE_OURS_NULL)),
        ],
        floor.n,
    );
    assert_main_accounting::<S>(&acc);
}

/// What the main sample's accounting must show, whichever route its quantisations took:
/// every listed divergence met exactly once, the null path exercised, and the admissions
/// under their ceilings. Shared by `zone_from_geo_text_id_and_resolution` and by
/// `any_grid_answers_as_the_typed_grid_and_the_engine`, which quantise the same sample.
fn assert_main_accounting<S: Subject>(acc: &Accounted) {
    // Every listed divergence is in the sample and must be met exactly once, or the
    // table has gone stale.
    for (i, d) in DIVERGENCES
        .iter()
        .enumerate()
        .filter(|(_, d)| d.grid == S::NAME)
    {
        assert_eq!(
            acc.listed.get(&i).copied().unwrap_or(0),
            1,
            "listed divergence at ({}, {}, {}) was not met exactly once",
            d.lat,
            d.lon,
            d.res
        );
    }
    let m = measured::<S>();
    // The null path has to be exercised, not merely described.
    assert!(
        acc.nulls_agreed >= m.main_min_nulls_agreed,
        "only {} null answers agreed",
        acc.nulls_agreed
    );
    assert!(
        acc.unresolved <= m.main_max_unresolved,
        "{} unresolved quantisations",
        acc.unresolved
    );
    // Admissions are held under ceilings too, so that a class rule loosened by a later
    // edit shows as a jump here.
    let (ties, unreadable) = (acc.class_count(TIE), acc.class_count(UNREADABLE));
    assert!(
        ties <= m.main_max_ties && unreadable <= m.main_max_unreadable,
        "{ties} ties, {unreadable} unreadable"
    );
    assert_eq!(
        acc.class_count(TIE_OURS_NULL),
        m.main_ties_ours_null,
        "ties where we answer null in the main sample"
    );
}

/// Every one of the thirty seams between projection faces, at every resolution.
pub fn seams_between_faces_quantise_like_dggal<S: Subject>() {
    let mut floor = Floor::new("seam zone_from_geo", 30_000);
    let mut acc = Accounted::default();
    for (lat, lon) in seam_points() {
        for res in 0..=S::grid().max_resolution() {
            let ours = S::grid().zone_from_geo(lat, lon, res).unwrap();
            let theirs = o::zone_from_geo(S::ORACLE, lat, lon, i32::from(res));
            acc.quantisation::<S>("seam", lat, lon, res, ours.id(), theirs);
            check_identifier::<S>(ours.id(), res);
            floor.hit();
        }
    }
    eprintln!(
        "seams: {} nulls agreed, class-rule divergences by rule and resolution {:?}",
        acc.nulls_agreed, acc.by_class
    );
    let unreadable_at: Vec<u8> = acc
        .by_class
        .keys()
        .filter(|(rule, _)| *rule == UNREADABLE)
        .map(|&(_, res)| res)
        .collect();
    print_observed::<S>(
        &[
            ("seam_min_nulls_agreed", &acc.nulls_agreed),
            ("seam_unreadable_at", &unreadable_at),
            ("seam_unresolved", &acc.unresolved),
            ("seam_max_ties", &acc.class_count(TIE)),
            ("seam_max_unreadable", &acc.class_count(UNREADABLE)),
            ("seam_max_ties_ours_null", &acc.class_count(TIE_OURS_NULL)),
        ],
        floor.n,
    );
    let m = measured::<S>();
    // The sample must reach what was measured along the seams, or the class rules go
    // untested.
    assert!(
        acc.nulls_agreed >= m.seam_min_nulls_agreed,
        "the seam sample no longer reaches the strip"
    );
    assert!(
        m.seam_unreadable_at
            .iter()
            .all(|&r| acc.by_class.contains_key(&(UNREADABLE, r))),
        "the seam sample no longer meets the unreadable answers at resolutions {:?}",
        m.seam_unreadable_at
    );
    assert_eq!(
        acc.unresolved, m.seam_unresolved,
        "unresolved quantisations along the seams"
    );
    let (ties, unreadable) = (acc.class_count(TIE), acc.class_count(UNREADABLE));
    assert!(
        ties <= m.seam_max_ties && unreadable <= m.seam_max_unreadable,
        "{ties} ties, {unreadable} unreadable"
    );
    // The ties where our answer is the inferior one, null; see
    // `boundary_tie_on_a_mirror_plane`.
    let inferior = acc.class_count(TIE_OURS_NULL);
    assert!(
        inferior <= m.seam_max_ties_ours_null,
        "{inferior} ties where we answer null"
    );
}

pub fn centroids<S: Subject>() {
    let mut floor = Floor::new("centroid", 5000);
    let mut worst = 0.0f64;
    let mut null_zone = 0usize;
    for id in zones::<S>() {
        let c = S::grid().centroid(id);
        if null_geometry_departs::<S>(id) {
            assert_eq!((c.lat, c.lon), (0.0, 0.0), "the null zone's centroid");
            div::engine_null_zone_geometry::<S>().unwrap_or_else(|why| panic!("{why}"));
            null_zone += 1;
            floor.hit();
            continue;
        }
        assert_bit_identical(
            &format!("centroid of {} (seed {SEED:#x})", S::grid().text_id(id)),
            (c.lat, c.lon),
            o::centroid(S::ORACLE, id.0),
            &mut worst,
        );
        floor.hit();
    }
    eprintln!("centroid: worst gap {worst:e} degrees");
    if null_zone > 0 {
        eprintln!("centroid: the null zone {null_zone} time(s), at (0, 0) as decided");
    }
}

pub fn vertices_and_pentagons<S: Subject>() {
    let mut floor = Floor::new("vertices", 5000);
    let mut worst = 0.0f64;
    let mut pentagons = 0;
    let mut null_zone = 0usize;
    let mut windings = Windings::default();
    for id in zones::<S>() {
        let text = S::grid().text_id(id);
        let ours = S::grid().vertices(id);
        if null_geometry_departs::<S>(id) {
            let zeros: Vec<(f64, f64)> = ours.iter().map(|p| (p.lat, p.lon)).collect();
            assert_eq!(zeros, vec![(0.0, 0.0); 6], "the null zone's vertices");
            assert!(!S::grid().is_pentagon(id), "the null zone is a pentagon");
            div::engine_null_zone_geometry::<S>().unwrap_or_else(|why| panic!("{why}"));
            null_zone += 1;
            floor.hit();
            continue;
        }
        let theirs = windings.engine_ring::<S>(id);
        assert_eq!(
            ours.len(),
            theirs.len(),
            "vertex count of {text} (seed {SEED:#x})"
        );
        for (a, b) in ours.iter().zip(&theirs) {
            assert_bit_identical(
                &format!("vertex of {text} (seed {SEED:#x})"),
                (a.lat, a.lon),
                *b,
                &mut worst,
            );
        }
        // The oracle exposes no pentagon predicate, so the engine's vertex count serves
        // as one; the null zone's six zero vertices make it no pentagon on either side.
        assert_eq!(
            S::grid().is_pentagon(id),
            theirs.len() == 5,
            "pentagon flag of {text} (seed {SEED:#x})"
        );
        pentagons += usize::from(theirs.len() == 5);
        floor.hit();
    }
    assert!(pentagons >= 24, "only {pentagons} pentagons compared");
    eprintln!("vertices: worst gap {worst:e} degrees; {pentagons} pentagons");
    windings.finish::<S>("vertices", 34, 2, 0);
    if null_zone > 0 {
        eprintln!("vertices: the null zone {null_zone} time(s), six zero vertices as decided");
    }
}

/// The eC's `fixPoles` (ri5x6.ec:461-509), which moves any coordinate within 1e-8 of a
/// pole's 5x6 position onto the pole itself, is matched coordinate for coordinate. Every
/// centroid and vertex of the cells at and around both poles, at resolutions 17 to 19,
/// is the engine's to the last bit, so that wherever the engine puts a coordinate exactly
/// on a pole, ours is on the pole too, at the same longitude. The
/// longitude is the part that exercises the odd-grid flag: `fixPoles` gives -168.8 or
/// 11.2 on an odd grid and -78.8 or 101.2 on an even one, and `Grid` decides the flag
/// as the engine does.
pub fn the_poles_snap_as_the_engine_does<S: Subject>() {
    let (resolutions, min, snapped_at) = pole_snap_sample::<S>();
    let mut floor = Floor::new("pole snap", min);
    let mut on_pole = BTreeMap::<u8, usize>::new();
    let mut worst = 0.0f64;
    let mut windings = Windings::default();
    let mut ids = BTreeSet::new();
    for res in resolutions {
        for sign in [1.0, -1.0] {
            for j in 0..12 {
                // From 1e-5 degrees down to the pole itself, a decade at a time.
                let d = if j == 11 { 0.0 } else { 1e-5 * 10f64.powi(-j) };
                for i in 0..12 {
                    let lon = -180.0 + 30.0 * f64::from(i);
                    ids.insert(
                        S::grid()
                            .zone_from_geo(sign * (90.0 - d), lon, res)
                            .unwrap()
                            .id(),
                    );
                }
            }
        }
    }
    ids.remove(&ZoneId::NULL);
    for id in ids {
        let text = S::grid().text_id(id);
        let res = S::grid().resolution(id);
        let c = S::grid().centroid(id);
        let mut pairs = vec![((c.lat, c.lon), o::centroid(S::ORACLE, id.0))];
        let theirs = windings.engine_ring::<S>(id);
        let ours = S::grid().vertices(id);
        assert_eq!(
            ours.len(),
            theirs.len(),
            "vertex count of {text} (seed {SEED:#x})"
        );
        pairs.extend(ours.iter().map(|v| (v.lat, v.lon)).zip(theirs));
        for (a, b) in pairs {
            assert_bit_identical(
                &format!("{text} at the pole (seed {SEED:#x})"),
                a,
                b,
                &mut worst,
            );
            if b.0.abs() == 90.0 {
                *on_pole.entry(res).or_default() += 1;
            }
        }
        floor.hit();
    }
    eprintln!(
        "pole snap: worst gap {worst:e} degrees; snapped coordinates by resolution {on_pole:?}"
    );
    windings.finish::<S>("pole snap", 2, 4, 0);
    // Each parity must be exercised, or the longitude choice for it goes untested.
    assert!(
        snapped_at
            .iter()
            .all(|r| on_pole.get(r).is_some_and(|&n| n >= 4)),
        "too few snapped coordinates: {on_pole:?}"
    );
}

/// The pole-snap test's resolutions, its floor on the zones compared, and the resolutions at
/// which at least four coordinates must lie exactly on a pole. On aperture 7 the snap reaches
/// one vertex of each polar cell at resolution 18 and more at 19. On aperture 3 it reaches
/// the pole only at the two finest resolutions, 12 coordinates at 32 and 28 at 33, and the
/// sample takes every resolution, 286 zones in all, since the polar cells of the coarser
/// ones are few.
fn pole_snap_sample<S: Subject>() -> (std::ops::RangeInclusive<u8>, usize, &'static [u8]) {
    match S::T::APERTURE {
        7 => (17..=19, 80, &[18, 19]),
        3 => (0..=33, 286, &[32, 33]),
        a => panic!("no pole-snap sample is recorded for aperture {a}"),
    }
}

/// Deep digit paths under the two polar pentagons, base cells 0 and 11: at every
/// resolution from 15 to 19, the pentagon's own chain of zeros and every tail of up to
/// three digits after it, 3,430 texts, of which 570 pass through a deleted pentagon
/// child. Under the pentagon of base cell 0 the geometric parent search can name a zone
/// other than the one a walk down the digits came from, and the walk answers as the
/// engine only by keeping its own two previous zones, as the eC's does: a walk that
/// recomputed them gave twelve zones of this sample on each grid a centroid and
/// vertices other than the engine's.
///
/// Both sides refuse the same texts. Every zone read is one the engine reads back, and
/// its centroid and vertices are the engine's to the last bit. Its neighbour list is
/// compared as in `neighbours_as_the_engine_lists_them`. Under base cell 0 the paths lie
/// in the broken seam over the north pole, where the engine's lists hold each of the four
/// kinds of entry this crate drops, and the sample must meet every kind, so that each
/// arm of the rule is exercised; under base cell 11 the engine's lists must be taken
/// whole.
pub fn deep_polar_digit_paths_answer_as_the_engine<S: Subject<T = HexA7, I = Z7>>() {
    let g = S::grid();
    let mut floor = Floor::new("deep polar paths", 2_860);
    let mut refused = 0usize;
    let mut worst = 0.0f64;
    let mut windings = Windings::default();
    let (mut north, mut south) = (Neighbours::default(), Neighbours::default());
    for base in [0u8, 11] {
        for res in 15..=19usize {
            for tail in 0..343usize {
                let text = format!(
                    "{base:02}{}{}{}{}",
                    "0".repeat(res - 3),
                    tail / 49,
                    tail / 7 % 7,
                    tail % 7
                );
                let engine_id = o::zone_from_text(S::ORACLE, &text);
                let Ok(zone) = g.zone_from_text(&text) else {
                    assert_eq!(
                        engine_id,
                        o::NULL_ZONE,
                        "{text}: the engine reads a zone this crate refuses"
                    );
                    refused += 1;
                    continue;
                };
                let id = zone.id();
                assert_eq!(id.0, engine_id, "{text}: the engine reads another zone");
                assert!(
                    o::engine_can_read(S::ORACLE, id.0),
                    "the engine cannot read {text} back"
                );
                let c = g.centroid(id);
                assert_bit_identical(
                    &format!("centroid of {text}"),
                    (c.lat, c.lon),
                    o::centroid(S::ORACLE, id.0),
                    &mut worst,
                );
                let ours = g.vertices(id);
                let theirs = windings.engine_ring::<S>(id);
                assert_eq!(ours.len(), theirs.len(), "vertex count of {text}");
                for (a, b) in ours.iter().zip(&theirs) {
                    assert_bit_identical(
                        &format!("vertex of {text}"),
                        (a.lat, a.lon),
                        *b,
                        &mut worst,
                    );
                }
                let u = if base == 0 { &mut north } else { &mut south };
                u.compare::<S>(id);
                floor.hit();
            }
        }
    }
    assert_eq!(
        refused, 570,
        "texts refused, each passing through a deleted pentagon child, on both sides"
    );
    eprintln!("deep polar paths: worst gap {worst:e} degrees");
    windings.finish::<S>("deep polar paths", 0, 0, 16);
    north.print::<S>("deep polar paths under base cell 0");
    south.print::<S>("deep polar paths under base cell 11");
    assert_eq!(
        south.zones_dropping, 0,
        "the engine's lists under the south polar pentagon have entries left out"
    );
    for (k, n) in DROPPED.iter().zip(north.dropped) {
        assert!(
            n > 0,
            "no entry of the kind \"{k}\" was left out under base cell 0, so that arm of \
             the rule goes unexercised"
        );
    }
}

/// What our own neighbour list must satisfy wherever the zone is, including inside the
/// broken seams: real zones only, at the zone's resolution, each once, and never
/// the zone itself.
fn check_neighbour_list<S: Subject>(id: ZoneId, ns: &[ZoneId]) {
    let text = S::grid().text_id(id);
    assert!(
        ns.len() <= 6,
        "{text} has {} neighbours (seed {SEED:#x})",
        ns.len()
    );
    let distinct: HashSet<ZoneId> = ns.iter().copied().collect();
    assert_eq!(
        distinct.len(),
        ns.len(),
        "{text} lists a neighbour twice (seed {SEED:#x})"
    );
    for &n in ns {
        let nt = S::grid().text_id(n);
        assert!(
            n != ZoneId::NULL && n != id,
            "{text} lists {nt} (seed {SEED:#x})"
        );
        assert_eq!(
            S::grid().resolution(n),
            S::grid().resolution(id),
            "{text} lists {nt} (seed {SEED:#x})"
        );
        assert_eq!(
            S::grid().zone_from_text(&nt).map(|z| z.id()),
            Ok(n),
            "{text} lists {nt}, which is not a readable zone (seed {SEED:#x})"
        );
    }
}

/// Tallies of the neighbour lists compared with the engine's: the entries left out of the
/// engine's lists by kind, over every zone compared, and this crate's null zone.
#[derive(Default)]
struct Neighbours {
    /// Entries of the engine's lists left out, by kind, in the order of [`DROPPED`].
    dropped: [usize; 4],
    /// Zones whose engine list had an entry left out.
    zones_dropping: usize,
    /// This crate's null zone, met where the engine may not be asked about it (see
    /// [`engine_may_be_asked`]): no neighbours on this side, as on every grid. Only ever
    /// counted on aperture 3, the one aperture where the engine cannot read its own null
    /// zone back; on aperture 7 it can, so the null zone's neighbours are asked of it and
    /// compared like any other zone's.
    null_zone: usize,
}

impl Neighbours {
    /// Compares one zone's neighbour list with the engine's, as a sequence: it must be
    /// [`engine_list`], the engine's list in its order, less on aperture 7 the four kinds
    /// of entry this crate drops, each re-checked live as one of them.
    fn compare<S: Subject>(&mut self, id: ZoneId) {
        self.compare_list::<S>(id, &S::grid().neighbors(id));
    }

    /// As [`Neighbours::compare`], for a neighbour list obtained by another route than
    /// the typed grid, `AnyGrid`'s among them, so that it meets the engine under exactly
    /// the same rule.
    fn compare_list<S: Subject>(&mut self, id: ZoneId, ns: &[ZoneId]) {
        check_neighbour_list::<S>(id, ns);
        if !engine_may_be_asked::<S>(id.0) {
            // Only this crate's null zone may stand here: every other zone it names is one
            // the engine reads back. It has no neighbours, as `check_neighbour_list` and the
            // null-zone test require; the engine's are not asked for.
            assert_eq!(
                (id, ns.len()),
                (ZoneId::NULL, 0),
                "the engine cannot read {} back (seed {SEED:#x})",
                S::grid().text_id(id)
            );
            self.null_zone += 1;
            return;
        }
        let (want, dropped) = engine_list::<S>(id.0);
        let seq: Vec<u64> = ns.iter().map(|z| z.0).collect();
        assert_eq!(
            seq,
            want,
            "neighbours of {}: ours {:?}, the engine's less the four kinds {:?}, from the \
             engine's list {:?} (seed {SEED:#x})",
            S::grid().text_id(id),
            ns.iter().map(|&z| S::grid().text_id(z)).collect::<Vec<_>>(),
            want.iter()
                .map(|&z| o::text_id(S::ORACLE, z))
                .collect::<Vec<_>>(),
            o::neighbors(S::ORACLE, id.0)
                .into_iter()
                .map(|z| o::text_id(S::ORACLE, z))
                .collect::<Vec<_>>()
        );
        if dropped != [0; 4] {
            self.zones_dropping += 1;
        }
        for (n, d) in self.dropped.iter_mut().zip(dropped) {
            *n += d;
        }
    }

    /// Prints the tallies: the entries left out, on aperture 7, where the engine's lists
    /// are taken less the four kinds, and the null zone, where it was met.
    fn print<S: Subject>(&self, what: &str) {
        if S::T::APERTURE == 7 {
            eprintln!(
                "{what}: the engine's list had entries left out at {} zones",
                self.zones_dropping
            );
            for (k, n) in DROPPED.iter().zip(self.dropped) {
                eprintln!("  {k}: {n}");
            }
        }
        if self.null_zone > 0 {
            eprintln!(
                "{what}: the null zone {} time(s), where the engine cannot read its own",
                self.null_zone
            );
        }
    }
}

/// Every zone of [`zones`], its neighbour list compared with the engine's as a sequence:
/// on aperture 7 the engine's list less the four kinds of entry this crate drops, each
/// re-checked live, and on aperture 3 the engine's list itself (see [`Neighbours`]).
pub fn neighbours_as_the_engine_lists_them<S: Subject>() {
    let mut floor = Floor::new("neighbours", 5000);
    let mut bit63 = 0;
    let mut u = Neighbours::default();
    for id in zones::<S>() {
        u.compare::<S>(id);
        if id != ZoneId::NULL && id.0 >> 63 == 1 {
            bit63 += 1;
        }
        floor.hit();
    }
    u.print::<S>("neighbours");
    eprintln!("neighbours: {bit63} with bit 63 set");
    // On aperture 7 the base cells 8 to 11 set bit 63 of every identifier under them, which
    // must be reached; no zone of aperture 3 sets it.
    if S::T::APERTURE == 7 {
        assert!(bit63 >= 1000, "only {bit63} bit-63 zones compared");
    }
}

/// The measurement asked for: neighbour lists near both poles and along the seams between
/// the faces, at the resolutions of [`poles_and_seams_sample`] (13 to 19 on aperture 7),
/// each compared with the engine's as in [`neighbours_as_the_engine_lists_them`], with the
/// entries left out of the engine's lists tallied and printed. Near the south pole, which
/// lies on no broken seam, the engine's lists must be taken whole.
pub fn neighbours_near_the_poles_and_along_the_seams<S: Subject>() {
    let (resolutions, samples) = poles_and_seams_sample::<S>();
    let mut sample: BTreeMap<&'static str, BTreeSet<ZoneId>> = BTreeMap::new();
    let mut lons: Vec<f64> = (0..18).map(|i| -180.0 + 20.0 * f64::from(i)).collect();
    lons.extend([11.2, -168.8, 101.2, -78.8]);
    for &res in resolutions {
        for (pole, sign) in [("north pole", 1.0), ("south pole", -1.0)] {
            for j in 0..28 {
                // From 0.3 degrees down to about 1e-9, a third of a decade at a time.
                let d = 0.3 * 10f64.powf(-f64::from(j) / 3.0);
                for &lon in &lons {
                    let z = S::grid()
                        .zone_from_geo(sign * (90.0 - d), lon, res)
                        .unwrap()
                        .id();
                    sample.entry(pole).or_default().insert(z);
                }
            }
        }
        for (lat, lon) in seam_points() {
            let z = S::grid().zone_from_geo(lat, lon, res).unwrap().id();
            sample.entry("seams").or_default().insert(z);
        }
    }
    let mut south = usize::MAX;
    for (what, name, min) in samples {
        let mut floor = Floor::new(name, min);
        let ids = sample.remove(what).unwrap_or_default();
        let mut u = Neighbours::default();
        for &id in &ids {
            u.compare::<S>(id);
            floor.hit();
        }
        u.print::<S>(&format!("{what} ({} zones)", ids.len()));
        if what == "south pole" {
            south = u.zones_dropping;
        }
    }
    assert_eq!(
        south, 0,
        "the engine's lists near the south pole have entries left out (seed {SEED:#x})"
    );
}

/// The resolutions of `neighbours_near_the_poles_and_along_the_seams`, and its three samples,
/// each with the name of its floor and the floor: one for each sample, so that neither pole's
/// sample can vanish unnoticed behind the others. On aperture 7, measured: 1,788 zones near
/// the north pole, 2,020 near the south pole and 6,435 along the seams. On aperture 3 the
/// resolutions are those where the engine answers its null zone near the poles, 0 and 1;
/// the first two and the last two at which it answers identifiers it cannot read back along
/// the seams, 2 and 3 and 20 and 21, and one between; and the finest two, where the pole
/// snap reaches the vertices. The zones are those each grid's projection makes of the
/// points, and so differ by grid. Measured, and the floors: 1,074, 1,075 and 3,536 zones
/// on ISEA3H, 1,072, 1,073 and 3,547 on IVEA3H, and 1,077, 1,078 and 3,542 on RTEA3H.
type PolesAndSeams = (&'static [u8], [(&'static str, &'static str, usize); 3]);
fn poles_and_seams_sample<S: Subject>() -> PolesAndSeams {
    match S::T::APERTURE {
        7 => (
            &[13, 14, 15, 16, 17, 18, 19],
            [
                ("north pole", "north-pole neighbours", 1700),
                ("south pole", "south-pole neighbours", 1900),
                ("seams", "seam neighbours", 6000),
            ],
        ),
        3 => {
            let [north, south, seams] = match S::NAME {
                "ISEA3H" => [1_074, 1_075, 3_536],
                "IVEA3H" => [1_072, 1_073, 3_547],
                "RTEA3H" => [1_077, 1_078, 3_542],
                other => panic!("no pole and seam floors are recorded for {other}"),
            };
            (
                &[0, 1, 2, 3, 12, 20, 21, 32, 33],
                [
                    ("north pole", "north-pole neighbours", north),
                    ("south pole", "south-pole neighbours", south),
                    ("seams", "seam neighbours", seams),
                ],
            )
        }
        a => panic!("no pole and seam sample is recorded for aperture {a}"),
    }
}

/// The three kinds of entry in the engine's lists of parents and of children that the
/// aperture-7 grids drop, in the order in which [`engine_relatives`] tries them, for the
/// tallies.
const DROPPED_RELATIVES: [&str; 3] = [
    "the engine's null zone",
    "an identifier the engine cannot read back",
    "a repeat of an entry already listed",
];

/// One of the engine's lists of parents or of children as an aperture-7 grid must give it,
/// and how many entries of each kind of [`DROPPED_RELATIVES`] it left out: the list in its
/// order, less the engine's null zone, an identifier the engine cannot read back, and a
/// repeat of an entry already kept, each entry classified live by the first of the three
/// kinds that fits it, so that every entry left out is re-checked against the engine as one
/// of the three on every run. The kinds are tried in the order in which this crate drops
/// them: an identifier that does not read back is never kept, so a second instance of it is
/// counted as unreadable and not as a repeat.
fn engine_relatives<S: Subject>(list: Vec<u64>) -> (Vec<u64>, [usize; 3]) {
    let mut dropped = [0usize; 3];
    let mut kept: Vec<u64> = Vec::new();
    for t in list {
        let kind = if t == o::NULL_ZONE {
            0
        } else if !o::engine_can_read(S::ORACLE, t) {
            1
        } else if kept.contains(&t) {
            2
        } else {
            kept.push(t);
            continue;
        };
        dropped[kind] += 1;
    }
    (kept, dropped)
}

/// The engine's parents and children of `id`, each as [`engine_relatives`] leaves it, with
/// the entries left out of the two lists together, asked only where [`engine_may_be_asked`]
/// allows it and otherwise a failure naming the identifier.
fn engine_parents_and_children<S: Subject>(id: ZoneId) -> (Vec<u64>, Vec<u64>, [usize; 3]) {
    assert!(
        engine_may_be_asked::<S>(id.0),
        "the engine cannot read {} ({:#018x}) back, and is not asked for its parents and \
         children (seed {SEED:#x})",
        S::grid().text_id(id),
        id.0
    );
    let (parents, left_out) = engine_relatives::<S>(o::parents(S::ORACLE, id.0));
    let (children, more) = engine_relatives::<S>(o::children(S::ORACLE, id.0));
    let dropped = [0, 1, 2].map(|k| left_out[k] + more[k]);
    (parents, children, dropped)
}

/// The classes that the comparison of the hierarchy counts by level (see [`Hierarchy`]).
/// Every one of them was met in the broken seams alone.
///
/// The entries left out of the engine's lists of parents and of children, by the kind of
/// entry, in the order of [`DROPPED_RELATIVES`].
const LEFT_OUT_OF_PARENTS: [&str; 3] = [
    "parents: the engine's null zone left out",
    "parents: an identifier the engine cannot read back left out",
    "parents: a repeat left out",
];
const LEFT_OUT_OF_CHILDREN: [&str; 3] = [
    "children: the engine's null zone left out",
    "children: an identifier the engine cannot read back left out",
    "children: a repeat left out",
];
/// A zone above level 0 with no parent: the engine names none, or none that it reads back.
const NO_PARENT: &str = "no parent above level 0";
/// A zone whose first parent is not the zone its identifier names with the last digit
/// dropped.
const FIRST_PARENT_ELSEWHERE: &str = "a first parent that is not the digit parent";
/// A zone of whose parents the engine lists one identifier twice.
const PARENT_TWICE: &str = "a parent listed twice by the engine";
/// A zone with one parent that is no centroid child; away from the broken seams the zones
/// with one parent are the centroid children, and no others.
const ONE_PARENT_NO_CENTROID_CHILD: &str = "one parent and no centroid child";
/// A zone at which the first of the engine's parents that the engine calls a centroid child
/// is an identifier it cannot read back, which this crate does not hand out: the centroid
/// parent here is the first centroid child among the parents that are kept, or none.
const CENTROID_PARENT_UNREADABLE: &str =
    "the engine's first centroid child among the parents cannot be read back";
/// A zone below the finest level with no children at all.
const NO_CHILDREN: &str = "no children below the finest level";
/// A child, counted at the level of the zone that lists it, which does not name that zone
/// among its own parents. Every child is asked from level 14, where the seams begin to
/// shorten the lists; below it two or three of the thirteen are asked, in rotation, so that
/// the count there is of the children asked and not of every child.
const CHILD_NAMES_ANOTHER: &str = "a child that does not name the zone among its parents";
/// A zone whose identifier with one digit more does not name the first of its children,
/// taken as a set.
const DIGIT_CHILDREN_ELSEWHERE: &str = "the digit children not the first of the children";
/// A zone of whose children an entry is left out and which is all the same the zone found
/// at its own centroid: the shortened lists are met by one who quantises positions.
const CHILDREN_LEFT_OUT_AT_A_ZONE_OF_A_POSITION: &str =
    "children left out at a zone found at its own centroid";

/// The three predicates of the hierarchy, in the order of `Hierarchy::predicates`.
const PREDICATES: [&str; 3] = ["is_immediate_child_of", "is_sibling_of", "is_ancestor_of"];

/// The coarsest level at which the broken seams were measured to shorten a list of
/// children; the lists of parents depart from 15.
const SEAMS_FROM_LEVEL: u8 = 14;

/// What is recorded of one sample of one grid: every class the comparison counted, with its
/// count at each level at which it was met. The comparison is of the whole table, so a class
/// that appears, vanishes or moves by one zone at one level fails it.
type HierarchyCounts = &'static [(&'static str, &'static [(u8, usize)])];

/// How far a point lies from the nearest broken seam of the grid under test, in degrees.
pub(super) fn broken_seam_distance_deg<S: Subject>(lat: f64, lon: f64) -> f64 {
    let v = icosahedron_vertices();
    unresolved::regions::<S>()
        .iter()
        .flat_map(|r| r.edges)
        .map(|&(i, j)| dist_to_arc_deg((lat, lon), v[i], v[j]))
        .fold(f64::INFINITY, f64::min)
}

/// Whether `id` is the zone found at its own centroid, at its own level: the engine's zone
/// at the engine's centroid of it, and this crate's zone at this crate's, which must give
/// one answer to the question. Where it is not, no one who obtains zones from positions is
/// handed `id` at its own centre.
fn found_at_its_own_centroid<S: Subject>(id: ZoneId) -> bool {
    let g = S::grid();
    let (lat, lon) = o::centroid(S::ORACLE, id.0);
    let theirs = o::zone_from_geo(S::ORACLE, lat, lon, o::level(S::ORACLE, id.0));
    let c = g.centroid(id);
    let ours = g
        .zone_from_geo(c.lat, c.lon, g.resolution(id))
        .unwrap()
        .id();
    assert_eq!(
        ours == id,
        theirs == id.0,
        "{}: at its own centroid this crate finds {} and the engine {} (seed {SEED:#x})",
        g.text_id(id),
        g.text_id(ours),
        o::text_id(S::ORACLE, theirs)
    );
    theirs == id.0
}

/// Tallies of the aperture-7 hierarchy compared with the engine's, over one sample.
#[derive(Default)]
struct Hierarchy {
    /// Zones compared.
    zones: usize,
    /// Every class met, by level: see the constants above and [`HierarchyCounts`].
    classes: BTreeMap<&'static str, BTreeMap<u8, usize>>,
    /// Whether the zone under comparison was counted in a class.
    in_a_class: bool,
    /// Zones counted in a class, and the farthest that one lies from a broken seam, in
    /// degrees.
    zones_in_a_class: usize,
    reach_deg: f64,
    /// Zones with one parent and with two; hexagons with thirteen children and pentagons
    /// with eleven; zones of the finest level; centroid children; zones with a centroid
    /// parent.
    one_parent: usize,
    two_parents: usize,
    thirteen: usize,
    eleven: usize,
    finest: usize,
    centroid_children: usize,
    with_centroid_parent: usize,
    /// Answers compared with the engine's for each of [`PREDICATES`]: `false`, then `true`.
    predicates: [[usize; 2]; 3],
}

impl Hierarchy {
    fn count(&mut self, class: &'static str, level: u8, n: usize) {
        if n > 0 {
            *self
                .classes
                .entry(class)
                .or_default()
                .entry(level)
                .or_default() += n;
            self.in_a_class = true;
        }
    }

    /// One predicate of [`PREDICATES`] on the pair `(a, b)`, this crate's answer against the
    /// engine's, exactly; returns the answer.
    fn predicate<S: Subject<T = HexA7, I = Z7>>(&mut self, which: usize, a: u64, b: u64) -> bool {
        let g = S::grid();
        let (za, zb) = (g.zone(ZoneId(a)), g.zone(ZoneId(b)));
        let (ours, theirs) = match which {
            0 => (
                za.is_immediate_child_of(&zb),
                o::is_immediate_child_of(S::ORACLE, a, b),
            ),
            1 => (za.is_sibling_of(&zb), o::are_siblings(S::ORACLE, a, b)),
            _ => (za.is_ancestor_of(&zb), o::is_ancestor_of(S::ORACLE, a, b)),
        };
        assert_eq!(
            ours,
            theirs,
            "{} of {} and {}, against the engine's (seed {SEED:#x})",
            PREDICATES[which],
            za.text_id(),
            zb.text_id()
        );
        self.predicates[which][usize::from(ours)] += 1;
        ours
    }

    /// Compares one zone with the engine (see [`Hierarchy::lists_and_classes`] and
    /// [`Hierarchy::predicates_on_pairs`]), and holds every zone counted in a class to the
    /// band about the broken seams that the suites use, [`div::BROKEN_SEAM_BAND_DEG`].
    fn compare<S: Subject<T = HexA7, I = Z7>>(&mut self, id: ZoneId) {
        self.in_a_class = false;
        self.zones += 1;
        if let Some((parents, children)) = self.lists_and_classes::<S>(id) {
            self.predicates_on_pairs::<S>(id, &parents, &children);
        }
        if self.in_a_class {
            let g = S::grid();
            let c = g.centroid(id);
            assert!(
                div::near_broken_seam::<S>(c.lat, c.lon),
                "{} is counted in a class of the broken seams and lies {:e} degrees from \
                 them (seed {SEED:#x})",
                g.text_id(id),
                broken_seam_distance_deg::<S>(c.lat, c.lon)
            );
            self.zones_in_a_class += 1;
            self.reach_deg = self
                .reach_deg
                .max(broken_seam_distance_deg::<S>(c.lat, c.lon));
        }
    }

    /// The zone's parents and its children as sequences against the engine's lists less the
    /// three kinds, with every entry left out re-checked live (see [`engine_relatives`]) and
    /// counted by kind and by level; whether it is a centroid child; its centroid parent;
    /// its primary parent; and the classes that the broken seams produce, each counted by
    /// level. Returns the two lists, or nothing for the null zone, which has no relatives.
    fn lists_and_classes<S: Subject<T = HexA7, I = Z7>>(
        &mut self,
        id: ZoneId,
    ) -> Option<(Vec<u64>, Vec<u64>)> {
        let g = S::grid();
        let text = g.text_id(id);
        let texts = |list: &[u64]| -> Vec<String> {
            list.iter().map(|&z| o::text_id(S::ORACLE, z)).collect()
        };
        assert!(
            engine_may_be_asked::<S>(id.0),
            "the engine cannot read {text} ({:#018x}) back, and is not asked for its parents \
             and children (seed {SEED:#x})",
            id.0
        );
        let raw_parents = o::parents(S::ORACLE, id.0);
        let raw_children = o::children(S::ORACLE, id.0);
        let (parents, left_out_of_parents) = engine_relatives::<S>(raw_parents.clone());
        let (children, left_out_of_children) = engine_relatives::<S>(raw_children.clone());
        let ours: Vec<u64> = g.parents(id).iter().map(|z| z.0).collect();
        assert_eq!(
            ours,
            parents,
            "parents of {text}, against the engine's less the three kinds, from the engine's \
             list {:?} (seed {SEED:#x})",
            texts(&raw_parents)
        );
        let ours: Vec<u64> = g.children(id).iter().map(|z| z.0).collect();
        // Said first of the one case that an empty list could hide: where the engine names a
        // child that is kept, this crate names children.
        assert!(
            children.is_empty() || !ours.is_empty(),
            "{text} has no children here, and the engine's are {:?} (seed {SEED:#x})",
            texts(&children)
        );
        assert_eq!(
            ours,
            children,
            "children of {text}, against the engine's less the three kinds, from the engine's \
             list {:?} (seed {SEED:#x})",
            texts(&raw_children)
        );
        assert_eq!(
            g.zone(id).parent().map(|p| p.id().0),
            parents.first().copied(),
            "the primary parent of {text} is not the first of its parents (seed {SEED:#x})"
        );

        // The centroid child, and the centroid parent: the engine's own function answers its
        // null zone everywhere, and what it defines is taken from its own list.
        let centroid_child = g.is_centroid_child(id);
        assert_eq!(
            centroid_child,
            o::is_centroid_child(S::ORACLE, id.0),
            "whether {text} is a centroid child (seed {SEED:#x})"
        );
        assert_eq!(
            o::centroid_parent(S::ORACLE, id.0),
            o::NULL_ZONE,
            "the engine now names a centroid parent of {text} (seed {SEED:#x})"
        );
        let centroid_parent = g.centroid_parent(id).map(|z| z.0);
        let mut first_of_the_engines = None;
        for &p in &raw_parents {
            // Asked of the engine's raw parents, of which some are identifiers that it cannot
            // read back: on purpose, and safe, since its `isZoneCentroidChild` reads a digit
            // of the identifier and follows no relation of the zone.
            if o::is_centroid_child(S::ORACLE, p) {
                first_of_the_engines = Some(p);
                break;
            }
        }
        if id == ZoneId::NULL {
            // The null zone has no level and no relatives, on either side.
            assert!(
                raw_parents.is_empty() && raw_children.is_empty() && centroid_parent.is_none(),
                "the null zone has relatives (seed {SEED:#x})"
            );
            return None;
        }
        let res = g.resolution(id);
        if centroid_parent != first_of_the_engines {
            // The one way in which the two part: the engine's first is an identifier it
            // cannot read back, which is left out of the parents, and this crate answers
            // the first centroid child of those that remain.
            let unreadable = first_of_the_engines.expect("the engine names one where we do");
            let mut first_kept = None;
            for &p in &parents {
                if o::is_centroid_child(S::ORACLE, p) {
                    first_kept = Some(p);
                    break;
                }
            }
            assert!(
                !o::engine_can_read(S::ORACLE, unreadable) && centroid_parent == first_kept,
                "centroid parent of {text}: ours {:?}, the first centroid child of the \
                 engine's list {}, of {:?} (seed {SEED:#x})",
                centroid_parent.map(|z| g.text_id(ZoneId(z))),
                o::text_id(S::ORACLE, unreadable),
                texts(&raw_parents)
            );
            self.count(CENTROID_PARENT_UNREADABLE, res, 1);
        }
        self.centroid_children += usize::from(centroid_child);
        self.with_centroid_parent += usize::from(centroid_parent.is_some());

        // The entries left out, by kind and by level.
        for k in 0..3 {
            self.count(LEFT_OUT_OF_PARENTS[k], res, left_out_of_parents[k]);
            self.count(LEFT_OUT_OF_CHILDREN[k], res, left_out_of_children[k]);
        }

        // The parents.
        match parents.len() {
            1 => self.one_parent += 1,
            2 => self.two_parents += 1,
            _ => {}
        }
        if parents.len() == 1 && !centroid_child {
            self.count(ONE_PARENT_NO_CENTROID_CHILD, res, 1);
        }
        let twice = (1..raw_parents.len()).any(|i| raw_parents[..i].contains(&raw_parents[i]));
        self.count(PARENT_TWICE, res, usize::from(twice));
        if res == 0 {
            assert!(
                raw_parents.is_empty(),
                "{text} has a parent (seed {SEED:#x})"
            );
        } else {
            let shorter = &text[..text.len() - 1];
            let digit_parent = g.zone_from_text(shorter).map(|z| z.id().0).ok();
            assert_eq!(
                digit_parent.unwrap_or(o::NULL_ZONE),
                o::zone_from_text(S::ORACLE, shorter),
                "{shorter}, read on both sides (seed {SEED:#x})"
            );
            let class = match parents.first() {
                None => Some(NO_PARENT),
                Some(&p) if Some(p) != digit_parent => Some(FIRST_PARENT_ELSEWHERE),
                Some(_) => None,
            };
            if let Some(class) = class {
                self.count(class, res, 1);
                // What matters to a caller: such an identifier is not the zone of its own
                // centre, so that no one who quantises positions is handed it there.
                assert!(
                    !found_at_its_own_centroid::<S>(id),
                    "{text} ({class}) is the zone found at its own centroid: the anomalies \
                     of the parents are then met by positions too (seed {SEED:#x})"
                );
            }
        }

        // The children. The engine's own list is never short: what shortens a list is what
        // is left out of it.
        if res == g.max_resolution() {
            assert!(
                raw_children.is_empty(),
                "{text} has children (seed {SEED:#x})"
            );
            self.finest += 1;
        } else {
            let pentagon = g.is_pentagon(id);
            assert_eq!(
                raw_children.len(),
                if pentagon { 11 } else { 13 },
                "the entries of the engine's list of the children of {text} (seed {SEED:#x})"
            );
            match (children.len(), pentagon) {
                (13, false) => self.thirteen += 1,
                (11, true) => self.eleven += 1,
                _ => {}
            }
            self.count(NO_CHILDREN, res, usize::from(children.is_empty()));
            let digit_children: BTreeSet<u64> = (0..7)
                .filter_map(|d| g.zone_from_text(&format!("{text}{d}")).ok())
                .map(|z| z.id().0)
                .collect();
            let first: BTreeSet<u64> = children
                .iter()
                .take(digit_children.len())
                .copied()
                .collect();
            self.count(
                DIGIT_CHILDREN_ELSEWHERE,
                res,
                usize::from(first != digit_children),
            );
            if left_out_of_children != [0; 3] && found_at_its_own_centroid::<S>(id) {
                self.count(CHILDREN_LEFT_OUT_AT_A_ZONE_OF_A_POSITION, res, 1);
            }
        }
        Some((parents, children))
    }

    /// The three predicates against the engine's `isZoneImmediateChildOf`,
    /// `areZonesSiblings` and `isZoneAncestorOf`, exactly, on pairs drawn from the lists of
    /// the zone, which [`Hierarchy::lists_and_classes`] has compared, and from this crate's
    /// lists of those zones, each of which the engine must read back before it is asked.
    ///
    /// The pairs are drawn in rotation, so that the comparison costs little at each zone and
    /// reaches every position of a list over the sample. At every zone: the zone with each
    /// parent, and with its children, each of them from the level at which the broken seams
    /// begin to shorten the lists and two or three below it. At one zone in four, about one
    /// parent: the parent as an ancestor, and as no child and no sibling; a sibling through
    /// it, both ways; a zone that shares a parent with that sibling, and perhaps none with
    /// this zone; each grandparent; and a child as a descendant. At one zone in sixteen, the
    /// two questions that are answered only at the top of the hierarchy: whether the zone
    /// is an ancestor of its parent, and whether a child of a grandparent is one of the zone.
    fn predicates_on_pairs<S: Subject<T = HexA7, I = Z7>>(
        &mut self,
        id: ZoneId,
        parents: &[u64],
        children: &[u64],
    ) {
        let g = S::grid();
        let text = g.text_id(id);
        let res = g.resolution(id);
        let asked = |z: ZoneId| {
            assert!(
                engine_may_be_asked::<S>(z.0),
                "the engine cannot read {} back, a relative of a relative of {text} (seed \
                 {SEED:#x})",
                g.text_id(z)
            );
            z.0
        };
        let n = self.zones;
        for &p in parents {
            assert!(
                self.predicate::<S>(0, id.0, p),
                "{text} is no child of its parent (seed {SEED:#x})"
            );
        }
        for (j, &c) in children.iter().enumerate() {
            if res >= SEAMS_FROM_LEVEL || j % 6 == n % 6 {
                let named = self.predicate::<S>(0, c, id.0);
                self.count(CHILD_NAMES_ANOTHER, res, usize::from(!named));
            }
        }
        if n % 4 != 0 {
            return;
        }
        let turn = n / 4;
        if let Some(&c) = children.get(turn % children.len().max(1)) {
            let child = self.predicate::<S>(0, c, id.0);
            assert_eq!(
                self.predicate::<S>(2, id.0, c),
                child,
                "{text} and its child {} (seed {SEED:#x})",
                g.text_id(ZoneId(c))
            );
        }
        let Some(&p) = parents.get(turn % parents.len().max(1)) else {
            return;
        };
        assert!(
            self.predicate::<S>(2, p, id.0),
            "the parent of {text} is no ancestor of it (seed {SEED:#x})"
        );
        assert!(
            !self.predicate::<S>(0, p, id.0) && !self.predicate::<S>(1, id.0, p),
            "the parent of {text} is its child or its sibling (seed {SEED:#x})"
        );
        let siblings = g.children(ZoneId(p));
        if let Some(&s) = siblings.get(turn % siblings.len().max(1)) {
            let s = asked(s);
            assert_eq!(
                self.predicate::<S>(1, id.0, s),
                self.predicate::<S>(1, s, id.0),
                "{text} and {} are siblings one way only (seed {SEED:#x})",
                g.text_id(ZoneId(s))
            );
            if let Some(&q) = g.parents(ZoneId(s)).last() {
                let cousins = g.children(q);
                if let Some(&t) = cousins.get(turn / 13 % cousins.len().max(1)) {
                    self.predicate::<S>(1, id.0, asked(t));
                }
            }
        }
        let grandparents = g.parents(ZoneId(p));
        for &gp in &grandparents {
            let gp = asked(gp);
            assert!(
                self.predicate::<S>(2, gp, id.0) && !self.predicate::<S>(0, id.0, gp),
                "the grandparent {} of {text} (seed {SEED:#x})",
                g.text_id(ZoneId(gp))
            );
        }
        if turn % 4 != 0 {
            return;
        }
        assert!(
            !self.predicate::<S>(2, id.0, p),
            "{text} is an ancestor of its parent (seed {SEED:#x})"
        );
        if let Some(&gp) = grandparents.get(turn / 4 % grandparents.len().max(1)) {
            let uncles = g.children(gp);
            if let Some(&u) = uncles.get(turn / 4 % uncles.len().max(1)) {
                assert_eq!(
                    self.predicate::<S>(2, asked(u), id.0),
                    parents.contains(&u.0),
                    "{}, a child of a grandparent of {text}, as its ancestor (seed {SEED:#x})",
                    g.text_id(u)
                );
            }
        }
    }

    /// Prints the tallies of one sample; holds them to the floors given, in the order of
    /// the tallies as they are printed and then of the predicates, `false` before `true`;
    /// and holds the classes to what is recorded of them: every class, at every level,
    /// exactly.
    fn finish<S: Subject>(&self, what: &str, floors: [usize; 13], recorded: HierarchyCounts) {
        eprintln!(
            "{what}: {} zones; {} with one parent and {} with two; {} hexagons with thirteen \
             children and {} pentagons with eleven; {} of the finest level; {} centroid \
             children and {} zones with a centroid parent; {} zones counted in a class of \
             the broken seams, the farthest {:e} degrees from them",
            self.zones,
            self.one_parent,
            self.two_parents,
            self.thirteen,
            self.eleven,
            self.finest,
            self.centroid_children,
            self.with_centroid_parent,
            self.zones_in_a_class,
            self.reach_deg
        );
        for (name, [no, yes]) in PREDICATES.iter().zip(self.predicates) {
            eprintln!("  {name}: {yes} pairs true and {no} false, as the engine answers");
        }
        for (class, by_level) in &self.classes {
            eprintln!(
                "  {class}: {} {by_level:?}",
                by_level.values().sum::<usize>()
            );
        }
        let [[c0, c1], [s0, s1], [a0, a1]] = self.predicates;
        let tallies = [
            self.one_parent,
            self.two_parents,
            self.thirteen,
            self.eleven,
            self.finest,
            self.centroid_children,
            self.with_centroid_parent,
            c0,
            c1,
            s0,
            s1,
            a0,
            a1,
        ];
        assert!(
            tallies.iter().zip(floors).all(|(&n, min)| n >= min),
            "{what} on {}: the tallies {tallies:?} are below the floors {floors:?} (seed \
             {SEED:#x})",
            S::NAME
        );
        let recorded: BTreeMap<&str, BTreeMap<u8, usize>> = recorded
            .iter()
            .map(|(class, by_level)| (*class, by_level.iter().copied().collect()))
            .collect();
        assert_eq!(
            self.classes,
            recorded,
            "{what} on {}: the classes counted by level, against what is recorded of the \
             grid (seed {SEED:#x})",
            S::NAME
        );
    }
}

/// What is recorded of the hierarchy of the grid under test over the zones of [`zones`], the
/// sample that positions give. The zones are those each grid's projection makes of the
/// points, and so the counts differ by grid. No list of parents has an entry left out, no
/// zone lacks a parent or has another first parent than its identifier names, and the
/// children that the identifier names are the first of the list everywhere; what the sample
/// does meet, at 49 zones on IGEO7, 39 on IVEA7H and 38 on RTEA7H, is children left out,
/// children that name other parents, and a zone that is no centroid child with one parent.
fn main_hierarchy_counts<S: Subject>() -> HierarchyCounts {
    match S::NAME {
        "IGEO7" => &[
            (LEFT_OUT_OF_CHILDREN[0], &[(14, 42), (16, 41), (18, 21)]),
            (LEFT_OUT_OF_CHILDREN[2], &[(15, 14), (16, 2), (17, 8)]),
            (
                CHILDREN_LEFT_OUT_AT_A_ZONE_OF_A_POSITION,
                &[(14, 12), (15, 4), (16, 7), (17, 3), (18, 7)],
            ),
            (
                CHILD_NAMES_ANOTHER,
                &[(14, 10), (15, 27), (16, 25), (17, 20), (18, 5)],
            ),
            (
                ONE_PARENT_NO_CENTROID_CHILD,
                &[(15, 1), (16, 1), (17, 1), (19, 1)],
            ),
        ],
        "IVEA7H" => &[
            (LEFT_OUT_OF_CHILDREN[0], &[(14, 26), (16, 29), (18, 21)]),
            (LEFT_OUT_OF_CHILDREN[2], &[(15, 2), (17, 6)]),
            (
                CHILDREN_LEFT_OUT_AT_A_ZONE_OF_A_POSITION,
                &[(14, 8), (15, 1), (16, 7), (17, 3), (18, 7)],
            ),
            (
                CHILD_NAMES_ANOTHER,
                &[(14, 6), (15, 18), (16, 15), (17, 22), (18, 5)],
            ),
            (ONE_PARENT_NO_CENTROID_CHILD, &[(17, 1), (19, 2)]),
        ],
        "RTEA7H" => &[
            (LEFT_OUT_OF_CHILDREN[0], &[(14, 26), (16, 29), (18, 21)]),
            (LEFT_OUT_OF_CHILDREN[2], &[(15, 2), (17, 6)]),
            (
                CHILDREN_LEFT_OUT_AT_A_ZONE_OF_A_POSITION,
                &[(14, 8), (15, 1), (16, 7), (17, 3), (18, 7)],
            ),
            (
                CHILD_NAMES_ANOTHER,
                &[(14, 6), (15, 18), (16, 15), (17, 22), (18, 5)],
            ),
            (ONE_PARENT_NO_CENTROID_CHILD, &[(17, 1), (19, 1)]),
        ],
        other => panic!(
            "nothing is recorded of the hierarchy of {other}: the classes printed above are \
             what this run met, and the table is to be filled from them only once each has \
             been characterised, and never from another grid"
        ),
    }
}

/// What is recorded of the hierarchy over the identifiers of [`seam_texts`]: 2,276 zones,
/// the same texts on every grid, and, as measured on IGEO7, IVEA7H and RTEA7H, the same
/// count of every class at every level. A grid of which nothing has been measured is not
/// handed these counts. The rustdoc of `Grid::parents` quotes three of these numbers (the
/// 2,276 identifiers, the 107 with no parent and the 98 with another first parent), and the
/// crate's README the same: a change of this table is a change of both.
fn seam_hierarchy_counts<S: Subject>() -> HierarchyCounts {
    match S::NAME {
        "IGEO7" | "IVEA7H" | "RTEA7H" => &[
            (LEFT_OUT_OF_PARENTS[1], &[(17, 8), (19, 83)]),
            (LEFT_OUT_OF_PARENTS[2], &[(17, 1)]),
            (LEFT_OUT_OF_CHILDREN[0], &[(14, 11), (16, 84), (18, 364)]),
            (
                LEFT_OUT_OF_CHILDREN[1],
                &[(15, 18), (16, 42), (17, 325), (18, 617)],
            ),
            (
                LEFT_OUT_OF_CHILDREN[2],
                &[(15, 12), (16, 20), (17, 55), (18, 120)],
            ),
            (NO_PARENT, &[(16, 6), (17, 2), (18, 88), (19, 11)]),
            (FIRST_PARENT_ELSEWHERE, &[(17, 8), (18, 5), (19, 85)]),
            (PARENT_TWICE, &[(17, 1)]),
            (
                ONE_PARENT_NO_CENTROID_CHILD,
                &[(15, 3), (16, 8), (17, 8), (18, 9), (19, 72)],
            ),
            (CENTROID_PARENT_UNREADABLE, &[(17, 8), (19, 83)]),
            (
                CHILD_NAMES_ANOTHER,
                &[(14, 2), (15, 31), (16, 41), (17, 273), (18, 295)],
            ),
            (
                DIGIT_CHILDREN_ELSEWHERE,
                &[(15, 3), (16, 6), (17, 48), (18, 93)],
            ),
            (
                CHILDREN_LEFT_OUT_AT_A_ZONE_OF_A_POSITION,
                &[(14, 4), (15, 4), (16, 11), (17, 3), (18, 5)],
            ),
        ],
        other => panic!(
            "nothing is recorded of the hierarchy of {other} at the identifiers built by \
             text: the classes printed above are what this run met"
        ),
    }
}

/// Every zone of [`zones`] compared with the engine as [`Hierarchy::compare`] does: its
/// parents and its children as sequences with the engine's `getZoneParents` and
/// `getZoneChildren`, in its order, less the three kinds of entry this crate drops, each
/// re-checked live as it is left out (see [`engine_relatives`]); whether it is a centroid
/// child, its centroid parent and its primary parent; and the three predicates on pairs
/// drawn from the lists. The null zone, which the sample holds at odd resolutions of 15 and
/// above, has no parents and no children on either side, and neither has a zone of
/// resolution 0 a parent nor a zone of the finest resolution a child.
///
/// What the sample must reach has a floor, so that no arm of the comparison falls silent,
/// and what the broken seams do to the lists is counted by class and by level and held to
/// what is recorded of the grid, exactly. The sample is of zones that positions gave, and
/// at every one of them the primary parent is the zone that the identifier names with its
/// last digit dropped: the two anomalies of the parents are not met here.
pub fn parents_and_children_as_the_engine_lists_them<S: Subject<T = HexA7, I = Z7>>() {
    let mut floor = Floor::new("parents and children", 5000);
    let mut h = Hierarchy::default();
    for id in zones::<S>() {
        h.compare::<S>(id);
        floor.hit();
    }
    // Measured, on IGEO7, IVEA7H and RTEA7H: 982, 889 and 906 zones with one parent and
    // 4,613, 4,678 and 4,656 with two; 4,961, 4,945 and 4,940 hexagons with thirteen
    // children, and 171 pentagons with eleven and 434 zones of the finest level on each; 978,
    // 886 and 904 centroid children, and 1,573, 1,563 and 1,580 zones with a centroid
    // parent; and of the predicates, some 50,000 pairs for the immediate child, 5,600 for
    // the sibling and 5,900 for the ancestor, each answered both ways.
    h.finish::<S>(
        "parents and children",
        [
            850, 4_500, 4_900, 170, 430, 850, 1_500, 3_900, 45_000, 2_000, 3_300, 600, 5_000,
        ],
        main_hierarchy_counts::<S>(),
    );
    for class in [NO_PARENT, FIRST_PARENT_ELSEWHERE] {
        assert!(
            !h.classes.contains_key(class),
            "a zone that a position gave has {class} (seed {SEED:#x})"
        );
    }
}

/// The stride at which [`seam_texts`] takes the tails of five digits: prime to 7, so that
/// every digit takes every value, and wide enough for the sample to cost a few seconds.
const SEAM_TEXT_STRIDE: usize = 37;

/// Identifiers of the broken seam under base cell 0, built by text: at each level from 14
/// to 19, `00`, zeros and a tail of five digits, the tails taken at a stride of
/// [`SEAM_TEXT_STRIDE`] through the 16,807; and eight identifiers known for what the
/// engine answers of them, among them one of each anomaly of the parents. No position is
/// behind any of them, and the same texts are asked of every grid.
pub(super) fn seam_texts() -> Vec<String> {
    let mut texts: Vec<String> = [
        // No parent; a first parent that is not the digit parent; a parent listed twice.
        "00000000000000001311",
        "0000000000000001644",
        "0000000000000000136",
        // The zone found at the centroid of another identifier; and four whose children
        // have entries left out: null zones, repeats, the polar pentagon's own, and four
        // entries that the engine cannot read back.
        "000000000000000140",
        "0000000000000005",
        "00000000000000005",
        "0000000000000000",
        "00052626050026015",
    ]
    .map(String::from)
    .to_vec();
    for res in 14..=19usize {
        for tail in (res..16_807).step_by(SEAM_TEXT_STRIDE) {
            texts.push(format!(
                "00{}{}{}{}{}{}",
                "0".repeat(res - 5),
                tail / 2401,
                tail / 343 % 7,
                tail / 49 % 7,
                tail / 7 % 7,
                tail % 7
            ));
        }
    }
    let mut seen = HashSet::new();
    texts.retain(|t| seen.insert(t.clone()));
    texts
}

/// The hierarchy at identifiers of the broken seam that are built by text, [`seam_texts`],
/// each compared with the engine as [`Hierarchy::compare`] does, like the zones that
/// positions give. Here the anomalies of the parents are met: a zone above level 0 with no
/// parent, a first parent that is not the zone the identifier names with its last digit
/// dropped, a parent listed twice. Each is counted by level, exactly, and of the first two
/// it is asserted at every instance that the identifier is not the zone found at its own
/// centroid, on either side: they are identifiers that a text can name and that no
/// position of the sample quantises to.
///
/// A text that this crate refuses passes through a deleted child of the pentagon, and the
/// engine must answer its null zone for it.
pub fn the_hierarchy_at_seam_identifiers_built_by_text<S: Subject<T = HexA7, I = Z7>>() {
    let g = S::grid();
    let mut h = Hierarchy::default();
    let mut refused = 0usize;
    for text in seam_texts() {
        let engine_id = o::zone_from_text(S::ORACLE, &text);
        let Ok(zone) = g.zone_from_text(&text) else {
            assert_eq!(
                engine_id,
                o::NULL_ZONE,
                "{text}: the engine reads a zone this crate refuses"
            );
            refused += 1;
            continue;
        };
        assert_eq!(
            zone.id().0,
            engine_id,
            "{text}: the engine reads another zone"
        );
        h.compare::<S>(zone.id());
    }
    assert_eq!(
        (h.zones, refused),
        (2_276, 456),
        "identifiers compared, and texts refused on both sides"
    );
    // Measured, the same on the three grids: 398 zones with one parent and 1,771 with two;
    // 1,684 hexagons with thirteen children, and no pentagon with eleven, for the one
    // pentagon of the sample, of level 14, has an entry left out; 378 zones of the finest
    // level; 326
    // centroid children and 572 zones with a centroid parent; and of the predicates 28,968
    // pairs for the immediate child, 2,172 for the sibling and 2,298 for the ancestor.
    h.finish::<S>(
        "seam identifiers built by text",
        [
            390, 1_750, 1_650, 0, 370, 320, 560, 2_200, 26_000, 840, 1_300, 280, 2_000,
        ],
        seam_hierarchy_counts::<S>(),
    );
    // A floor of none cannot fail, so the count is held exactly beside it: a pentagon with
    // eleven children appearing in this sample must be noticed.
    assert_eq!(
        h.eleven, 0,
        "pentagons with eleven children among the seam identifiers built by text, where the \
         one pentagon of the sample has an entry left out"
    );
}

/// The four methods of the hierarchy and its three predicates, asked of hostile identifiers:
/// the null zone, the identifiers at the ends of the range, those with no cell behind them,
/// and 3,000 drawn at random from every level. None panics or fails to end; the typed grid,
/// `AnyGrid` and `Zone` give one answer; a list holds real zones alone, each once, one level
/// from the zone; an identifier with no geometry has no relatives; and the predicates agree
/// with the lists. The engine is not asked: these are not identifiers it may be handed.
pub fn the_hierarchy_does_not_panic_on_hostile_input<S: Subject<T = HexA7, I = Z7>>() {
    let g = S::grid();
    let any = rs4dggs::get_grid(S::NAME).unwrap();
    let ids = super::geometry::hostile_identifiers::<S>();
    let floor = ids.len() / 20;
    let (mut with_relatives, mut without) = (0usize, 0usize);
    let mut previous = ZoneId::NULL;
    for (i, id) in ids.into_iter().enumerate() {
        let id = ZoneId(id);
        let at = format!("{:#018x}", id.0);
        let zone = g.zone(id);
        let parents = g.parents(id);
        let children = g.children(id);
        let centroid_parent = g.centroid_parent(id);
        let centroid_child = g.is_centroid_child(id);
        let ids_of = |zones: Vec<rs4dggs::Zone<'_, _>>| -> Vec<ZoneId> {
            zones.iter().map(|z| z.id()).collect()
        };
        assert_eq!(any.parents(id), parents, "{at}");
        assert_eq!(ids_of(zone.parents()), parents, "{at}");
        assert_eq!(any.children(id), children, "{at}");
        assert_eq!(ids_of(zone.children()), children, "{at}");
        assert_eq!(any.centroid_parent(id), centroid_parent, "{at}");
        assert_eq!(
            zone.centroid_parent().map(|z| z.id()),
            centroid_parent,
            "{at}"
        );
        assert_eq!(any.is_centroid_child(id), centroid_child, "{at}");
        assert_eq!(zone.is_centroid_child(), centroid_child, "{at}");
        assert_eq!(any.parent(id), parents.first().copied(), "{at}");
        assert_eq!(
            zone.parent().map(|z| z.id()),
            parents.first().copied(),
            "{at}"
        );

        assert!(parents.len() <= 2 && children.len() <= 13, "{at}");
        let res = i32::from(g.resolution(id));
        for (list, level) in [(&parents, res - 1), (&children, res + 1)] {
            let distinct: HashSet<ZoneId> = list.iter().copied().collect();
            assert_eq!(distinct.len(), list.len(), "{at} lists a zone twice");
            for &r in list {
                let text = g.text_id(r);
                assert_eq!(
                    g.zone_from_text(&text).map(|z| z.id()),
                    Ok(r),
                    "{at} lists {text}, which is no zone"
                );
                assert_eq!(i32::from(g.resolution(r)), level, "{at} lists {text}");
            }
        }
        if let Some(p) = centroid_parent {
            assert!(
                parents.contains(&p) && g.is_centroid_child(p),
                "{at}: the centroid parent {}",
                g.text_id(p)
            );
        }
        if g.extent(id).is_none() {
            assert!(
                parents.is_empty() && children.is_empty() && centroid_parent.is_none(),
                "{at} has no geometry and has relatives"
            );
        }

        // The predicates: a zone with itself, with each parent, with two children in
        // rotation, and with the identifier drawn before it, which is of any level or of
        // none.
        assert!(
            !zone.is_immediate_child_of(&zone)
                && !zone.is_sibling_of(&zone)
                && !zone.is_ancestor_of(&zone),
            "{at} is its own relative"
        );
        for &p in &parents {
            let parent = g.zone(p);
            assert!(
                zone.is_immediate_child_of(&parent) && parent.is_ancestor_of(&zone),
                "{at} and its parent {}",
                g.text_id(p)
            );
            assert!(any.is_immediate_child_of(id, p) && any.is_ancestor_of(p, id));
        }
        for &c in children.iter().skip(i % 13).take(2) {
            assert_eq!(
                g.zone(c).is_immediate_child_of(&zone),
                g.parents(c).contains(&id),
                "{at} and its child {}",
                g.text_id(c)
            );
        }
        let other = g.zone(previous);
        assert_eq!(
            zone.is_immediate_child_of(&other),
            any.is_immediate_child_of(id, previous),
            "{at} and {:#018x}",
            previous.0
        );
        assert_eq!(
            zone.is_sibling_of(&other),
            any.is_sibling_of(id, previous),
            "{at} and {:#018x}",
            previous.0
        );
        // One in eight, for the walk to the top of the hierarchy is the dearest of them.
        if i % 8 == 0 {
            assert_eq!(
                zone.is_ancestor_of(&other),
                any.is_ancestor_of(id, previous),
                "{at} and {:#018x}",
                previous.0
            );
        }
        previous = id;

        if parents.is_empty() && children.is_empty() {
            without += 1;
        } else {
            with_relatives += 1;
        }
    }
    assert!(
        with_relatives >= floor && without >= floor,
        "{with_relatives} identifiers with relatives, {without} without"
    );
    eprintln!(
        "hostile input to the hierarchy on {}: {with_relatives} with relatives, {without} without",
        S::NAME
    );
}

/// The walk `Disk` defines, to radius `k` about `centre`, composed from the lists `list`
/// gives: each ring from the previous ring's zones in their own order, each zone's list in
/// its own order, first arrivals only. Returns the walk and, for each radius, the length of
/// the walk to it.
fn disk_walk(
    centre: u64,
    k: usize,
    mut list: impl FnMut(u64) -> Vec<u64>,
) -> (Vec<u64>, Vec<usize>) {
    let mut walk = vec![centre];
    let mut seen = HashSet::from([centre]);
    let mut ring = vec![centre];
    let mut ends = vec![1usize];
    for _ in 0..k {
        let mut next = Vec::new();
        for &z in &ring {
            for n in list(z) {
                if seen.insert(n) {
                    next.push(n);
                }
            }
        }
        walk.extend(&next);
        ends.push(walk.len());
        ring = next;
    }
    (walk, ends)
}

/// The disk of radius 0, 1 and 2 about every tenth zone of the main sample.
///
/// The engine has no disk, but it has the neighbour lists a disk is composed from. Each
/// walk must be the one `Disk` defines, composed from the engine's own lists as
/// [`engine_list`] gives them, which on aperture 7 leave out the four kinds of entry this
/// crate drops: the same zones, in the same order. It is also checked against the same
/// walk composed from this crate's own lists; against itself, on a second call; and
/// against its own rings. It must yield no zone twice, and the walk of each radius must be
/// the beginning of the walk of radius 2, since it is the same walk stopped at an earlier
/// ring.
pub fn disk_is_composed_neighbours<S: Subject>() {
    let mut floor = Floor::new("disk", 500);
    let g = S::grid();
    let mut null_centres = 0usize;
    for id in zones::<S>().into_iter().step_by(10) {
        let text = g.text_id(id);
        if !engine_may_be_asked::<S>(id.0) {
            // Only this crate's null zone may stand here (see `Neighbours::compare_list`):
            // its disk is itself alone, and the engine is not asked for its neighbours.
            let walk: Vec<ZoneId> = g.zone(id).disk(2).map(|z| z.id()).collect();
            assert_eq!(
                (id, walk),
                (ZoneId::NULL, vec![ZoneId::NULL]),
                "the engine cannot read {text} back (seed {SEED:#x})"
            );
            null_centres += 1;
            floor.hit();
            continue;
        }

        let (engine, _) = disk_walk(id.0, 2, |z| engine_list::<S>(z).0);
        let (rule, ring_ends) = disk_walk(id.0, 2, |z| {
            g.neighbors(ZoneId(z)).into_iter().map(|n| n.0).collect()
        });
        let whole: Vec<u64> = g.zone(id).disk(2).map(|z| z.id().0).collect();
        let texts =
            |w: &[u64]| -> Vec<String> { w.iter().map(|&z| g.text_id(ZoneId(z))).collect() };
        assert_eq!(
            whole,
            engine,
            "disk(2) of {text}: ours {:?}, composed from the engine's lists {:?} (seed {SEED:#x})",
            texts(&whole),
            texts(&engine)
        );
        assert_eq!(
            whole, rule,
            "disk(2) of {text} is not in the order defined (seed {SEED:#x})"
        );
        let again: Vec<u64> = g.zone(id).disk(2).map(|z| z.id().0).collect();
        assert_eq!(
            again, whole,
            "two calls of disk(2) of {text} differ (seed {SEED:#x})"
        );
        let rings: Vec<Vec<ZoneId>> = g.disk(id, 2).rings().collect();
        assert_eq!(
            rings.concat().iter().map(|z| z.0).collect::<Vec<_>>(),
            whole,
            "the rings of disk(2) of {text} are not its walk (seed {SEED:#x})"
        );
        for k in 0..=2usize {
            let walk: Vec<u64> = g.zone(id).disk(k as u32).map(|z| z.id().0).collect();
            assert_eq!(
                walk[..],
                whole[..ring_ends[k]],
                "disk({k}) of {text} is not the walk of disk(2) to ring {k} (seed {SEED:#x})"
            );
            let got: HashSet<u64> = walk.iter().copied().collect();
            assert_eq!(
                got.len(),
                walk.len(),
                "disk({k}) of {text} yields a zone twice (seed {SEED:#x})"
            );
        }
        floor.hit();
    }
    if null_centres > 0 {
        eprintln!(
            "disk: the null zone {null_centres} time(s), where the engine cannot read its own"
        );
    }
}

/// `AnyGrid`, the handle of a caller that chooses its grid at run time, reached through
/// `get_grid` by the grid's name in lower case, as a caller might spell it, and compared
/// both with the typed grid it names and with the engine.
///
/// Over the main sample at every resolution, its quantisation must equal the typed grid's,
/// identifier for identifier, and meets the engine's through the same accounting as
/// `zone_from_geo_text_id_and_resolution`, so that a disagreement is admitted only as that
/// test admits it. Over every fifth zone of `zones`, every other operation it offers must
/// equal the typed grid's answer exactly, coordinates bit for bit; and where the engine
/// can judge, the engine's too: the text identifier, its reading back and the resolution;
/// the centroid and vertices bit for bit, and the vertex count as the pentagon flag;
/// the neighbour list, as `neighbours_as_the_engine_lists_them` compares it; and the parent and
/// children, against the engine's own lists (see `engine_parent_and_children`). The
/// sub-zone methods, the hierarchy predicates and the disk, which the engine is not asked
/// to judge here, must equal the typed grid's answers.
pub fn any_grid_answers_as_the_typed_grid_and_the_engine<S: Subject>() {
    let g = S::grid();
    let any = rs4dggs::get_grid(&S::NAME.to_lowercase())
        .unwrap_or_else(|e| panic!("get_grid({}) in lower case: {e}", S::NAME));
    assert_eq!(any.name(), S::NAME, "the grid get_grid names");
    assert!(
        rs4dggs::GRID_NAMES.contains(&S::NAME),
        "GRID_NAMES lacks {}",
        S::NAME
    );
    assert_eq!(any.max_resolution(), g.max_resolution());
    assert_eq!(i32::from(any.max_resolution()), o::max_level(S::ORACLE));

    let mut floor = Floor::new("AnyGrid zone_from_geo", main_sample_floor::<S>());
    let mut acc = Accounted::default();
    for (lat, lon) in inputs::<S>() {
        for res in 0..=g.max_resolution() {
            let id = any.zone_from_geo(lat, lon, res).unwrap_or_else(|e| {
                panic!("AnyGrid zone_from_geo({lat}, {lon}, {res}): {e} (seed {SEED:#x})")
            });
            assert_eq!(
                id,
                g.zone_from_geo(lat, lon, res).unwrap().id(),
                "AnyGrid zone_from_geo({lat}, {lon}, {res}) (seed {SEED:#x})"
            );
            let theirs = o::zone_from_geo(S::ORACLE, lat, lon, i32::from(res));
            acc.quantisation::<S>("AnyGrid", lat, lon, res, id, theirs);
            floor.hit();
        }
    }
    drop(floor);
    eprintln!(
        "AnyGrid zone_from_geo: {} nulls agreed, {} listed divergences met, class-rule \
         divergences by rule and resolution {:?}",
        acc.nulls_agreed,
        acc.listed.values().sum::<usize>(),
        acc.by_class
    );
    assert_main_accounting::<S>(&acc);

    let bits = |ps: &[rs4dggs::GeoPoint]| -> Vec<(u64, u64)> {
        ps.iter()
            .map(|p| (p.lat.to_bits(), p.lon.to_bits()))
            .collect()
    };
    // Every fifth zone of `zones`: measured, 1,122 zones on IGEO7, and on aperture 3, where
    // each grid's projection makes its own zones of the main sample, 1,176 on ISEA3H and
    // IVEA3H and 1,175 on RTEA3H.
    let min = match (S::T::APERTURE, S::NAME) {
        (7, _) => 1100,
        (3, "ISEA3H") => 1_176,
        (3, "IVEA3H") => 1_176,
        (3, "RTEA3H") => 1_175,
        (a, name) => panic!("no AnyGrid floor is recorded for {name}, of aperture {a}"),
    };
    let mut floor = Floor::new("AnyGrid zones", min);
    let mut u = Neighbours::default();
    let mut worst = 0.0f64;
    let mut windings = Windings::default();
    for id in zones::<S>().into_iter().step_by(5) {
        let z = g.zone(id);
        let text = g.text_id(id);
        let at = format!("{text} (seed {SEED:#x})");

        // The identifier.
        assert_eq!(any.text_id(id), text, "text of {at}");
        assert_eq!(
            o::text_id(S::ORACLE, id.0),
            text,
            "the engine's text of {at}"
        );
        assert_eq!(any.zone_from_text(&text), Ok(id), "reading back {at}");
        assert_eq!(
            o::zone_from_text(S::ORACLE, &text),
            id.0,
            "the engine reading back {at}"
        );
        assert_eq!(any.resolution(id), g.resolution(id), "resolution of {at}");
        if id == ZoneId::NULL {
            // The one deliberate difference, documented at `Grid::resolution`.
            assert_eq!(
                (any.resolution(id), o::level(S::ORACLE, id.0)),
                (0, engine_null_level::<S>()),
                "the null zone's resolution"
            );
        } else {
            assert_eq!(
                i32::from(any.resolution(id)),
                o::level(S::ORACLE, id.0),
                "level of {at}"
            );
        }

        // The geometry. The null zone's departs from the engine's by decision on aperture
        // 3, and is compared with the typed grid's and re-checked on the engine's side.
        let c = any.centroid(id);
        let tc = g.centroid(id);
        assert_eq!(bits(&[c]), bits(&[tc]), "centroid of {at}");
        let vs = any.vertices(id);
        assert_eq!(bits(&vs), bits(&g.vertices(id)), "vertices of {at}");
        assert_eq!(
            any.is_pentagon(id),
            g.is_pentagon(id),
            "pentagon flag of {at}"
        );
        if null_geometry_departs::<S>(id) {
            div::engine_null_zone_geometry::<S>().unwrap_or_else(|why| panic!("{why}"));
        } else {
            assert_bit_identical(
                &format!("AnyGrid centroid of {at}"),
                (c.lat, c.lon),
                o::centroid(S::ORACLE, id.0),
                &mut worst,
            );
            let theirs = windings.engine_ring::<S>(id);
            assert_eq!(vs.len(), theirs.len(), "vertex count of {at}");
            for (a, b) in vs.iter().zip(&theirs) {
                assert_bit_identical(
                    &format!("AnyGrid vertex of {at}"),
                    (a.lat, a.lon),
                    *b,
                    &mut worst,
                );
            }
            assert_eq!(
                any.is_pentagon(id),
                theirs.len() == 5,
                "pentagon flag of {at} against the engine's vertex count"
            );
        }

        // The neighbours.
        let ns = any.neighbors(id);
        assert_eq!(ns, g.neighbors(id), "neighbours of {at}");
        u.compare_list::<S>(id, &ns);

        // The hierarchy.
        let parent = any.parent(id);
        assert_eq!(parent, z.parent().map(|p| p.id()), "parent of {at}");
        assert_eq!(any.parents(id), g.parents(id), "parents of {at}");
        assert_eq!(
            any.centroid_parent(id),
            g.centroid_parent(id),
            "centroid parent of {at}"
        );
        assert_eq!(
            any.is_centroid_child(id),
            g.is_centroid_child(id),
            "centroid child {at}"
        );
        if let Some(p) = parent {
            let pz = g.zone(p);
            assert_eq!(
                any.is_immediate_child_of(id, p),
                z.is_immediate_child_of(&pz),
                "{at} as an immediate child of its parent"
            );
            assert_eq!(
                any.is_ancestor_of(p, id),
                pz.is_ancestor_of(&z),
                "the parent of {at} as its ancestor"
            );
        }
        let children = any.children(id);
        assert_eq!(children, g.children(id), "children of {at}");
        engine_parent_and_children::<S>(id, parent, &children, &at);
        if let Some(&n) = ns.first() {
            assert_eq!(
                any.is_sibling_of(id, n),
                z.is_sibling_of(&g.zone(n)),
                "{at} and a neighbour as siblings"
            );
            assert_eq!(
                any.sub_zone_index(id, n),
                g.sub_zone_index(id, n),
                "sub-zone index of a neighbour of {at}"
            );
        }

        // The sub-zones, and the disk.
        for depth in [0u8, 1] {
            assert_eq!(
                any.count_sub_zones(id, depth),
                g.count_sub_zones(id, depth),
                "count_sub_zones({depth}) of {at}"
            );
            assert_eq!(
                any.first_sub_zone(id, depth),
                g.first_sub_zone(id, depth),
                "first_sub_zone({depth}) of {at}"
            );
            assert_eq!(
                any.sub_zones(id, depth),
                g.sub_zones(id, depth),
                "sub_zones({depth}) of {at}"
            );
            assert_eq!(
                any.sub_zone_at_index(id, depth, 0),
                g.sub_zone_at_index(id, depth, 0),
                "sub_zone_at_index({depth}, 0) of {at}"
            );
        }
        assert_eq!(
            any.sub_zone_index(id, id),
            g.sub_zone_index(id, id),
            "sub-zone index of {at} in itself"
        );
        let disk: Vec<ZoneId> = any.disk(id, 1).collect();
        assert_eq!(disk, g.disk(id, 1).collect::<Vec<_>>(), "disk(1) of {at}");
        let rings: Vec<Vec<ZoneId>> = any.disk(id, 1).rings().collect();
        assert_eq!(
            rings,
            g.disk(id, 1).rings().collect::<Vec<_>>(),
            "rings of disk(1) of {at}"
        );
        floor.hit();
    }
    u.print::<S>("AnyGrid neighbours");
    eprintln!("AnyGrid: centroids and vertices, worst gap {worst:e} degrees");
    // The polar pentagons among every fifth zone: measured, 8 on ISEA3H, 4 on IVEA3H and 6 on
    // RTEA3H.
    let polar_min = match S::NAME {
        "ISEA3H" => 8,
        "IVEA3H" => 4,
        "RTEA3H" => 6,
        _ => 0,
    };
    let no_area = match S::NAME {
        "IVEA7H" => 1,
        "RTEA7H" => 2,
        _ => 0,
    };
    windings.finish::<S>("AnyGrid vertices", polar_min, no_area, 0);
}

/// The engine's side of the parent and children of `id` that `AnyGrid` answered, for
/// `any_grid_answers_as_the_typed_grid_and_the_engine`. They are the engine's own lists,
/// asked where [`engine_may_be_asked`] allows: the parent is the first of its parents, and
/// the children are its children in its order. On aperture 7 the lists are taken less the
/// three kinds of entry this crate drops, as in
/// `parents_and_children_as_the_engine_lists_them`. On aperture 3 they are taken whole, but
/// at the finest resolution, where this crate answers no children.
fn engine_parent_and_children<S: Subject>(
    id: ZoneId,
    parent: Option<ZoneId>,
    children: &[ZoneId],
    at: &str,
) {
    match S::T::APERTURE {
        7 => {
            let (their_parents, their_children, _) = engine_parents_and_children::<S>(id);
            assert_eq!(
                parent.map(|p| p.0),
                their_parents.first().copied(),
                "parent of {at}, the first of the engine's parents less the three kinds"
            );
            let ours: Vec<u64> = children.iter().map(|c| c.0).collect();
            assert_eq!(
                ours, their_children,
                "children of {at}, against the engine's less the three kinds"
            );
        }
        3 => {
            if !engine_may_be_asked::<S>(id.0) {
                // Only this crate's null zone may stand here, which has no relations.
                assert!(
                    id == ZoneId::NULL && parent.is_none() && children.is_empty(),
                    "the engine cannot read {at} back"
                );
                return;
            }
            let theirs = o::parents(S::ORACLE, id.0);
            assert_eq!(
                parent.map(|p| p.0),
                theirs.first().copied(),
                "parent of {at}, the first of the engine's parents"
            );
            let theirs = o::children(S::ORACLE, id.0);
            if S::grid().resolution(id) < S::grid().max_resolution() {
                let ours: Vec<u64> = children.iter().map(|c| c.0).collect();
                assert_eq!(ours, theirs, "children of {at}, against the engine's");
            } else {
                // At the finest resolution this crate answers none, where the engine's are
                // identifiers of level 34 that it cannot read back (see
                // `children_at_resolution_33_are_none`).
                assert!(children.is_empty(), "{at} has children");
                for c in theirs {
                    assert_eq!(
                        o::level(S::ORACLE, c),
                        34,
                        "the engine's child {c:#018x} of {at}"
                    );
                    let ct = o::text_id(S::ORACLE, c);
                    assert_eq!(
                        o::zone_from_text(S::ORACLE, &ct),
                        o::NULL_ZONE,
                        "the engine's own parser reads its child {ct} of {at} back"
                    );
                }
            }
        }
        a => panic!("no engine hierarchy is recorded for aperture {a}"),
    }
}

/// The sub-zone order of the aperture-7 grids at depth 1, against the engine's: the count,
/// and the order entry by entry, in which an entry that the engine gives as its null zone,
/// or as an identifier that it cannot read back, is the null zone here, at the same place.
/// The first sub-zone and the zone at each index are the entries of that order; the engine
/// is not asked for the zone at an index, which it answers wrongly at a pentagon at an odd
/// depth. At depth 0 every answer is the zone itself. A zone of the finest resolution has
/// no sub-zones below it, where the engine answers nothing, and the request is refused.
/// The null zone is refused as an invalid zone, before its depth is looked at, by
/// `sub_zone_index` as by the others.
///
/// `sub_zone_index` of the first child of every zone is the engine's `getSubZoneIndex`,
/// asked live, and the order holds the child at it. It is the child's place in the order
/// at every zone whose children lie coarser than level 15. From level 15 the engine's walk
/// may find no index for a child that its order names, and then there is none here, or
/// answer an index at which its order holds another zone, where there is none here either:
/// both are counted exactly, per grid, and each such zone is held to the band about the
/// broken seams.
///
/// `sub_zone_index` of two different zones at one level is a characterised divergence,
/// re-checked live at every zone of the sample that has a neighbour: the engine answers
/// 0, having compared only the two levels (`RI7H.ec:229-230`), and this crate answers
/// `None`, since at depth 0 a zone's only sub-zone is itself.
pub fn sub_zones_at_depth_1_as_the_engine_lists_them<S: Subject<T = HexA7, I = Z7>>() {
    the_limit_on_a_list_of_sub_zones_is_its_ceiling();
    let mut floor = Floor::new("sub-zones", 1000);
    let mut null_zone_checked = 0usize;
    // The two branches that not every zone reaches, counted so that neither can fall
    // silent: a first child to ask about, and a neighbour at the zone's own level.
    let mut children_checked = 0usize;
    let mut same_level_checked = 0usize;
    // The orders compared with the engine's and their entries; the entries that are the
    // null zone; and the zones of the finest resolution, which have no order below them.
    let mut orders = Floor::new("sub-zone orders at depth 1", 900);
    let (mut entries, mut null_entries, mut finest) = (0usize, 0usize, 0usize);
    // The first children that the order names and that have no index, as in the engine;
    // and those for which the engine answers the place of another zone, with none here.
    let (mut no_index, mut index_of_another) = (0usize, 0usize);
    let g = S::grid();
    for id in zones::<S>()
        .into_iter()
        .step_by(5)
        .chain(std::iter::once(ZoneId::NULL))
    {
        let z = g.zone(id);
        let text = z.text_id();
        let invalid = |r: std::result::Result<(), Error>| matches!(r, Err(Error::InvalidZone(_)));
        // The order at depth 1, where there is one.
        let mut order: Vec<ZoneId> = Vec::new();
        if id == ZoneId::NULL {
            for depth in [0u8, 1, 2, 19] {
                assert!(
                    invalid(z.count_sub_zones(depth).map(drop)),
                    "count_sub_zones({depth}) of {text}"
                );
                assert!(
                    invalid(z.first_sub_zone(depth).map(drop)),
                    "first_sub_zone({depth}) of {text}"
                );
                assert!(
                    invalid(z.sub_zones(depth).map(drop)),
                    "sub_zones({depth}) of {text}"
                );
                assert!(
                    invalid(z.sub_zone_at_index(depth, 0).map(drop)),
                    "sub_zone_at_index({depth}) of {text}"
                );
            }
        } else {
            // Depth 0 is the zone itself.
            assert_eq!(g.count_sub_zones(id, 0), Ok(1), "{text}");
            assert_eq!(g.sub_zones(id, 0), Ok(vec![id]), "{text}");
            assert_eq!(g.first_sub_zone(id, 0), Ok(id), "{text}");
            assert_eq!(g.sub_zone_at_index(id, 0, 0), Ok(id), "{text}");
            if g.resolution(id) == g.max_resolution() {
                assert!(
                    invalid(g.count_sub_zones(id, 1).map(drop)),
                    "count_sub_zones(1) of {text}"
                );
                assert!(
                    invalid(g.first_sub_zone(id, 1).map(drop)),
                    "first_sub_zone(1) of {text}"
                );
                assert!(
                    invalid(g.sub_zones(id, 1).map(drop)),
                    "sub_zones(1) of {text}"
                );
                assert!(
                    invalid(g.sub_zone_at_index(id, 1, 0).map(drop)),
                    "sub_zone_at_index(1) of {text}"
                );
                finest += 1;
            } else {
                assert!(
                    engine_may_be_asked::<S>(id.0),
                    "the engine cannot read {text} back (seed {SEED:#x})"
                );
                let count = g.count_sub_zones(id, 1).unwrap();
                assert_eq!(
                    count,
                    o::count_sub_zones(S::ORACLE, id.0, 1),
                    "count_sub_zones(1) of {text}, against the engine's"
                );
                order = g.sub_zones(id, 1).unwrap();
                let theirs: Vec<ZoneId> = o::sub_zones(S::ORACLE, id.0, 1)
                    .into_iter()
                    .map(|e| {
                        if e == o::NULL_ZONE || !o::engine_can_read(S::ORACLE, e) {
                            ZoneId::NULL
                        } else {
                            ZoneId(e)
                        }
                    })
                    .collect();
                assert_eq!(
                    order, theirs,
                    "sub_zones(1) of {text}, against the engine's"
                );
                assert_eq!(
                    order.len() as u64,
                    count,
                    "the length of the order of {text}"
                );
                assert_eq!(
                    g.first_sub_zone(id, 1),
                    Ok(order[0]),
                    "first_sub_zone(1) of {text}"
                );
                for (index, &sub) in order.iter().enumerate() {
                    assert_eq!(
                        g.sub_zone_at_index(id, 1, index as u64),
                        Ok(sub),
                        "sub_zone_at_index(1, {index}) of {text}"
                    );
                }
                assert!(
                    matches!(
                        g.sub_zone_at_index(id, 1, count),
                        Err(Error::IndexOutOfRange { .. })
                    ),
                    "sub_zone_at_index(1, {count}) of {text}"
                );
                entries += order.len();
                null_entries += order.iter().filter(|&&s| s == ZoneId::NULL).count();
                orders.hit();
            }
        }
        if let Some(c) = z.children().first() {
            // The index of a child is the engine's, and the order holds the child at it:
            // its place in the order, the first where the order names it twice. Where
            // there is none, the order does not hold the child and the engine finds none
            // either; or, in the broken seams, the engine finds none for a child that its
            // order names, or answers the place of another zone.
            let place = order.iter().position(|&s| s == c.id()).map(|i| i as u64);
            let engine = o::sub_zone_index(S::ORACLE, id.0, c.id().0);
            let ct = c.text_id();
            match z.sub_zone_index(c) {
                Ok(Some(i)) => assert_eq!(
                    (Some(i), i as i64),
                    (place, engine),
                    "sub_zone_index of {ct}, the first child of {text}, against its place in \
                     the order and against the engine's"
                ),
                Ok(None) => {
                    if engine == -1 {
                        no_index += usize::from(place.is_some());
                    } else {
                        assert_ne!(
                            usize::try_from(engine).ok().and_then(|i| order.get(i)),
                            Some(&c.id()),
                            "sub_zone_index of {ct}, the first child of {text}, is none, and \
                             the engine's is {engine}, where the order holds the child"
                        );
                        index_of_another += 1;
                    }
                    if place.is_some() || engine != -1 {
                        let centre = g.centroid(id);
                        assert!(
                            c.resolution() >= 15
                                && div::near_broken_seam::<S>(centre.lat, centre.lon),
                            "{ct}, the first child of {text}, has no index: the engine's is \
                             {engine}, the order holds it at {place:?}, and the zone lies \
                             {:e} degrees from the broken seams (seed {SEED:#x})",
                            broken_seam_distance_deg::<S>(centre.lat, centre.lon)
                        );
                    }
                }
                Err(e) => panic!("sub_zone_index of {ct}, the first child of {text}: {e}"),
            }
            children_checked += 1;
        }
        if let Some(&w) = S::grid().neighbors(id).first() {
            // A characterised divergence, re-checked live: a different zone at the same
            // level is no sub-zone of this one, where the engine answers 0 for it.
            let wt = S::grid().text_id(w);
            assert_eq!(
                z.sub_zone_index(&S::grid().zone(w)),
                Ok(None),
                "sub_zone_index of {wt}, a neighbour, in {text}"
            );
            assert_eq!(
                o::sub_zone_index(S::ORACLE, id.0, w.0),
                0,
                "the engine's sub-zone index of {wt}, a neighbour, in {text}"
            );
            same_level_checked += 1;
        }
        if id == ZoneId::NULL {
            // A characterised divergence: the engine answers 0 for its null zone within
            // itself, comparing the two levels before anything else, where this crate
            // validates both identifiers first and refuses. Re-checked live.
            assert!(
                matches!(z.sub_zone_index(&z), Err(Error::InvalidZone(_))),
                "sub_zone_index of {text} in itself"
            );
            assert_eq!(
                o::sub_zone_index(S::ORACLE, id.0, id.0),
                0,
                "the engine's sub-zone index of {text} in itself"
            );
            null_zone_checked += 1;
        } else {
            // Pinned as it stands, and reported: at the zone's own resolution, and
            // above it, `sub_zone_index` answers without asking for an order, as
            // py4dggs's grid.py does, so a zone is its own sub-zone at index 0 even
            // here.
            assert_eq!(
                z.sub_zone_index(&z),
                Ok(Some(0)),
                "sub_zone_index of {text} in itself"
            );
        }
        floor.hit();
    }
    assert!(
        null_zone_checked > 0,
        "the null-zone branch never ran (seed {SEED:#x})"
    );
    assert!(
        children_checked > 0,
        "no zone had a child to ask about (seed {SEED:#x})"
    );
    assert!(
        same_level_checked > 0,
        "no zone had a neighbour to ask about (seed {SEED:#x})"
    );
    eprintln!("sub-zones: null zone checked {null_zone_checked} time(s)");
    eprintln!(
        "sub-zones: a first child asked about at {children_checked} zones, a neighbour at \
         {same_level_checked}"
    );
    assert!(
        finest > 0,
        "no zone of the finest resolution was asked (seed {SEED:#x})"
    );
    eprintln!(
        "sub-zones: {} orders at depth 1 as the engine's, {entries} entries, {null_entries} of \
         them the null zone; {finest} zones of the finest resolution refused; {no_index} first \
         children that the order names have no index, as in the engine, and {index_of_another} \
         have none where the engine answers the place of another zone",
        orders.n
    );
    // Held exactly, per grid, so that the reading of the engine's null zone and of an
    // identifier it cannot read back as the null zone cannot fall silent: the orders, their
    // entries, the entries that are the null zone, and the zones of the finest resolution;
    // and, of the first children, those that have no index as in the engine, and those
    // that have none where the engine answers the place of another zone.
    let recorded = match S::NAME {
        "IGEO7" => [1_046, 13_530, 31, 76, 0, 0],
        "IVEA7H" => [1_021, 13_203, 27, 95, 0, 1],
        "RTEA7H" => [1_021, 13_207, 27, 94, 0, 1],
        other => panic!(
            "nothing is recorded of the sub-zone orders of {other} at depth 1: the counts \
             printed above are what this run met"
        ),
    };
    assert_eq!(
        [
            orders.n,
            entries,
            null_entries,
            finest,
            no_index,
            index_of_another
        ],
        recorded,
        "the orders at depth 1, their entries, the entries that are the null zone, the zones \
         of the finest resolution, the first children without an index as in the engine, and \
         those without one where the engine answers the place of another zone, against what \
         is recorded of the grid (seed {SEED:#x})"
    );
}

/// Input the engine answers meaninglessly and this crate refuses, each refusal re-checked
/// against the engine's recorded answer on every run (see the account of the third kind
/// in `expected_divergences`):
/// - every non-finite coordinate at every resolution, which the engine answers with the
///   pentagon of base cell 1, whose centroid is the icosahedron's vertex 1, and this
///   crate refuses, naming the coordinate at fault;
/// - 168 texts of twenty digits, which the engine reads as zones without geometry, or as
///   its null zone on a pentagon's deleted child, and this crate refuses by their length;
///   and one of twenty-one digits, which both refuse;
/// - `sub_zone_index` with an identifier that has no cell behind it, which the engine
///   answers, 0 within itself and -1 against a real zone or the reverse, and this crate
///   refuses.
///
/// Every count is fixed by construction and asserted exactly.
pub fn refusals_depart_from_the_engine_as_recorded<S: Subject<T = HexA7, I = Z7>>() {
    let mut floor = Floor::new("refusals", 140 + 169 + 12);

    let pentagon = icosahedron_vertices()[1];
    let mut non_finite = 0usize;
    for (lat, lon, name) in div::NON_FINITE_COORDINATES {
        for res in 0..=S::grid().max_resolution() {
            let at = format!("({lat}, {lon}) at {res}");
            assert_eq!(
                S::grid().zone_from_geo(lat, lon, res).err(),
                Some(Error::NonFinite { name }),
                "our answer at {at}"
            );
            let e = o::zone_from_geo(S::ORACLE, lat, lon, i32::from(res));
            assert_eq!(
                o::text_id(S::ORACLE, e),
                div::engine_answer_to_a_non_finite_coordinate(res),
                "the engine no longer answers as recorded at {at}"
            );
            let c = o::centroid(S::ORACLE, e);
            assert!(
                arc_deg(c, pentagon) < TOL_DEG,
                "the engine's answer at {at} is centred at {c:?}, not at vertex 1 {pentagon:?}"
            );
            non_finite += 1;
            floor.hit();
        }
    }
    assert_eq!(non_finite, 7 * 20);

    let (mut read, mut deleted) = (0usize, 0usize);
    for text in div::level_20_texts() {
        match div::level_20_text_refused::<S>(&text) {
            Ok(true) => read += 1,
            Ok(false) => deleted += 1,
            Err(why) => panic!("{why}"),
        }
        floor.hit();
    }
    assert_eq!(
        (read, deleted),
        (div::LEVEL_20_TEXTS_READ, div::LEVEL_20_TEXTS_DELETED),
        "level-20 texts the engine reads, and reads as its null zone"
    );
    // One digit more is past the packing, and the engine refuses it too: the boundary
    // of the departure is the twentieth digit and no other.
    let past = format!("01{}", "1".repeat(21));
    assert_eq!(o::zone_from_text(S::ORACLE, &past), o::NULL_ZONE, "{past}");
    assert!(
        matches!(S::grid().zone_from_text(&past), Err(Error::InvalidZone(_))),
        "{past}"
    );
    floor.hit();

    let real = S::grid().zone_from_text("0064156").unwrap().id();
    for bad in div::FABRICATED {
        for (id, sub, engine) in [(bad, bad, 0), (real.0, bad, -1), (bad, real.0, -1)] {
            let at = format!("{id:#018x} and {sub:#018x}");
            assert!(
                matches!(
                    S::grid().sub_zone_index(ZoneId(id), ZoneId(sub)),
                    Err(Error::InvalidZone(_))
                ),
                "our sub-zone index for {at}"
            );
            assert_eq!(
                o::sub_zone_index(S::ORACLE, id, sub),
                engine,
                "the engine's sub-zone index for {at}"
            );
            floor.hit();
        }
    }
    eprintln!(
        "refusals: {non_finite} non-finite coordinates, {read} level-20 texts read and \
         {deleted} read as null by the engine, {} fabricated sub-zone indices",
        3 * div::FABRICATED.len()
    );
}

/// The null zone agrees with DGGAL's in everything but its resolution, which this port
/// reports as 0 where DGGAL reports -1 (see `Grid::resolution`).
pub fn the_null_zone_matches_dggal_but_for_its_resolution<S: Subject<T = HexA7, I = Z7>>() {
    let n = ZoneId::NULL;
    assert_eq!(n.0, o::NULL_ZONE);
    assert_eq!(S::grid().text_id(n), o::text_id(S::ORACLE, o::NULL_ZONE));
    assert_eq!(
        o::zone_from_text(S::ORACLE, rs4dggs::NULL_TEXT),
        o::NULL_ZONE
    );
    assert_eq!(
        S::grid().zone_from_text(rs4dggs::NULL_TEXT).unwrap().id(),
        n
    );
    let c = S::grid().centroid(n);
    assert_eq!((c.lat, c.lon), o::centroid(S::ORACLE, o::NULL_ZONE));
    let ours: Vec<(f64, f64)> = S::grid()
        .vertices(n)
        .iter()
        .map(|p| (p.lat, p.lon))
        .collect();
    assert_eq!(ours, o::vertices(S::ORACLE, o::NULL_ZONE));
    assert!(S::grid().neighbors(n).is_empty() && o::neighbors(S::ORACLE, o::NULL_ZONE).is_empty());
    assert_eq!(
        (S::grid().resolution(n), o::level(S::ORACLE, o::NULL_ZONE)),
        (0, -1)
    );
}

pub fn level_20_null_geometry_matches_dggal<S: Subject<T = HexA7, I = Z7>>() {
    use rs4dggs::Address;
    let id = Z7::encode(&Address::new(1, &[1; 20]));
    let c = S::grid().centroid(id);
    assert_eq!((c.lat, c.lon), o::centroid(S::ORACLE, id.0));
    assert_eq!((c.lat, c.lon), (0.0, 0.0));
    let ours: Vec<(f64, f64)> = S::grid()
        .vertices(id)
        .iter()
        .map(|p| (p.lat, p.lon))
        .collect();
    assert_eq!(ours, o::vertices(S::ORACLE, id.0));
    assert!(o::neighbors(S::ORACLE, id.0).is_empty() && S::grid().neighbors(id).is_empty());
}

/// Re-checks each listed divergence against the live engine: our answer is deterministic
/// and the listed one, the engine still gives its listed answer and gives it every time,
/// and the evidence the entry rests on still holds.
pub fn listed_divergences_still_hold_with_their_evidence<S: Subject>() {
    // An entry whose grid is misspelt would belong to no grid and never be re-checked.
    for d in DIVERGENCES {
        assert!(
            super::subject::NAMES.contains(&d.grid),
            "the listed divergence at ({}, {}, {}) names the grid {:?}, which is no subject's",
            d.lat,
            d.lon,
            d.res,
            d.grid
        );
    }
    // Every entry listed for this grid; the table holds other grids' entries too. The
    // count the grid is held to is its own `Measured` value and not the length of the
    // table, which this test would otherwise take its floor from and check at once:
    // deleting every entry for a grid would then leave a floor of nothing, met by a
    // loop that never runs.
    let listed = DIVERGENCES.iter().filter(|d| d.grid == S::NAME).count();
    print_observed::<S>(&[("listed_divergences", &listed)], listed);
    let expected = measured::<S>().listed_divergences;
    assert_eq!(
        listed,
        expected,
        "the table lists {listed} divergences for {}, where {expected} are recorded for it",
        S::NAME
    );
    let mut floor = Floor::new("listed divergences", expected);
    for d in DIVERGENCES.iter().filter(|d| d.grid == S::NAME) {
        let at = format!("({}, {}, {})", d.lat, d.lon, d.res);
        let a = S::grid().zone_from_geo(d.lat, d.lon, d.res).unwrap().id();
        let b = S::grid().zone_from_geo(d.lat, d.lon, d.res).unwrap().id();
        assert_eq!(a, b, "not deterministic at {at}");
        assert_eq!(S::grid().text_id(a), d.ours, "our answer at {at}");
        let e = o::zone_from_geo(S::ORACLE, d.lat, d.lon, i32::from(d.res));
        assert_eq!(
            e,
            o::zone_from_geo(S::ORACLE, d.lat, d.lon, i32::from(d.res)),
            "DGGAL unstable at {at}"
        );
        assert_eq!(
            o::text_id(S::ORACLE, e),
            d.dggal,
            "DGGAL no longer answers as listed at {at}: stale entry"
        );
        match d.evidence {
            Evidence::EngineCannotReadItsOwnAnswer => {
                assert_eq!(
                    o::zone_from_text(S::ORACLE, d.dggal),
                    o::NULL_ZONE,
                    "DGGAL reads {} at {at}",
                    d.dggal
                );
                assert!(
                    S::grid().zone_from_text(d.dggal).is_err(),
                    "we read {} at {at}",
                    d.dggal
                );
                // On aperture 7 the unreadable answers lie inside the strip in which the
                // odd resolutions either side answer the null zone; on aperture 3 they lie
                // along two seams, where the resolutions either side need not.
                let either_side: &[u8] = if S::T::APERTURE == 7 {
                    &[d.res - 1, d.res + 1]
                } else {
                    &[]
                };
                for &r in either_side {
                    assert_eq!(
                        o::zone_from_geo(S::ORACLE, d.lat, d.lon, i32::from(r)),
                        o::NULL_ZONE,
                        "DGGAL at {r}, {at}"
                    );
                    assert_eq!(
                        S::grid().zone_from_geo(d.lat, d.lon, r).unwrap().id(),
                        ZoneId::NULL,
                        "ours at {r}, {at}"
                    );
                }
                let far = arc_deg((d.lat, d.lon), o::centroid(S::ORACLE, e));
                assert!(
                    far > 1.0,
                    "DGGAL's answer lies only {far} degrees from {at}"
                );
                let rule = class_rules::<S>()
                    .into_iter()
                    .find(|&(name, _)| name == UNREADABLE)
                    .map(|(_, rule)| rule)
                    .expect("a class rule for unreadable answers");
                rule(d.lat, d.lon, d.res, a.0, e).unwrap_or_else(|why| panic!("{at}: {why}"));
            }
            Evidence::BoundaryTie => {
                // Two real cells that share an edge, two vertices in common, and the
                // nudge evidence of the class rule.
                let ours = o::vertices(S::ORACLE, a.0);
                let theirs = o::vertices(S::ORACLE, e);
                let shared = ours
                    .iter()
                    .filter(|p| theirs.iter().any(|q| arc_deg(**p, *q) < TOL_DEG))
                    .count();
                assert!(shared >= 2, "the two cells at {at} share {shared} vertices");
                div::boundary_tie_on_a_mirror_plane::<S>(d.lat, d.lon, d.res, a.0, e)
                    .unwrap_or_else(|why| panic!("{at}: {why}"));
            }
        }
        floor.hit();
    }
}

// ---------------------------------------------------------------------------------
// Folded in from the former `tests/quantize_vs_dggal.rs`, which guarded the odd-level
// quantisation fix. Its sample, its null probe and both of its guard assertions are kept
// as they were; its table of carve-outs is now the first five entries of
// `common::expected_divergences::DIVERGENCES`, where their evidence is re-checked.
// ---------------------------------------------------------------------------------

/// A point 0.0009 degrees from the north pole, where the eC's own candidate search finds
/// no containing cell at odd resolutions of 15 and above, and the engine answers an
/// identifier it cannot read back at 18. It is the first listed divergence.
const NULL_PROBE: (f64, f64) = (89.999_129_147_255_16, 9.333_761_199_597_859);

/// A small linear congruential generator, so that this sample stays exactly as it was.
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
    }
}

/// A point 0.01 degrees from the north pole, at which the aperture-3 engine answers its
/// null zone at resolutions 0 and 1, as this crate does: the point py4dggs takes for it.
const APERTURE_3_NULL_PROBE: (f64, f64) = (89.99, -169.0);

/// A uniform spread, a band around each pole, a band around the orientation axis, and
/// the null probe; on aperture 3, its own null probe too, since the aperture-7 probe lies in
/// no null answer there. The random bands do not reach the strip on their own, which is why
/// the probe is a fixed point.
fn quantize_sample<S: Subject>() -> Vec<(f64, f64)> {
    let mut pts = vec![
        (38.7223, -9.1393),
        (0.0, 0.0),
        (-33.9, 151.2),
        (89.9, 45.0),
        (-89.9, -45.0),
        (45.0, -90.0),
        NULL_PROBE,
    ];
    let mut rng = Lcg(0x5EED);
    for _ in 0..80 {
        pts.push((rng.next() * 180.0 - 90.0, rng.next() * 360.0 - 180.0));
    }
    for _ in 0..80 {
        let d = rng.next() * 3.0;
        let lat = if rng.next() < 0.5 {
            90.0 - d
        } else {
            -90.0 + d
        };
        pts.push((lat, rng.next() * 360.0 - 180.0));
    }
    for _ in 0..80 {
        pts.push((
            31.717_474_411_461_1 + rng.next() * 2.0 - 1.0,
            -11.2 + rng.next() * 2.0 - 1.0,
        ));
    }
    if S::T::APERTURE == 3 {
        pts.push(APERTURE_3_NULL_PROBE);
    }
    pts
}

/// Quantisation at the topology layer, below `Grid`, at every resolution the grid admits.
///
/// py4dggs quantised wrongly at odd resolutions from 13, beyond the reach of its own
/// fuzzing and golden tables, so a comparison against it would have agreed with the
/// defect; this compares against the engine.
pub fn quantize_matches_dggal_at_every_resolution<S: Subject>() {
    use rs4dggs::{GridConfig, Projection};

    let g = S::P::build_geometry(&GridConfig::default());
    let pts = quantize_sample::<S>();
    let mut failures: Vec<String> = Vec::new();
    let mut nulls_agreed = 0usize;
    let mut carve_outs_seen = 0usize;
    // Measured: on aperture 7, 247 points at 20 resolutions, 4,940 quantisations; on
    // aperture 3, 248 points at 34 resolutions, 8,432, the whole sample.
    let min = match S::T::APERTURE {
        7 => 4900,
        3 => 8_432,
        a => panic!("no topology-layer floor is recorded for aperture {a}"),
    };
    let mut floor = Floor::new("topology-layer quantisation", min);

    for res in 0..=S::I::MAX_RESOLUTION {
        for &(lat, lon) in &pts {
            floor.hit();
            // A null answer maps to the sentinel, whose text is "(null)", exactly what
            // the oracle prints for it.
            let ours = S::I::to_text(match S::T::quantize(S::P::forward(&g, lat, lon), res) {
                Some(a) => S::I::encode(&a),
                None => ZoneId::NULL,
            });
            let theirs = o::text_id(
                S::ORACLE,
                o::zone_from_geo(S::ORACLE, lat, lon, i32::from(res)),
            );
            if ours == theirs {
                if ours == rs4dggs::NULL_TEXT {
                    nulls_agreed += 1;
                }
            } else if let Some(d) = div::lookup(S::NAME, lat, lon, res) {
                // A listed point admits only the answers recorded for it, on both sides,
                // as `Accounted::quantisation` requires of the same entries.
                if (ours.as_str(), theirs.as_str()) == (d.ours, d.dggal) {
                    carve_outs_seen += 1;
                } else {
                    failures.push(format!(
                        "resolution {res} at ({lat}, {lon}): ours {ours}, DGGAL {theirs}, where \
                         the listed divergence records ours {}, DGGAL {}",
                        d.ours, d.dggal
                    ));
                }
            } else {
                failures.push(format!(
                    "resolution {res} at ({lat}, {lon}): ours {ours}, DGGAL {theirs}"
                ));
            }
        }
    }

    // The floor is checked and printed here, before anything below can fail, so that its
    // count is on record in any case and its line comes first, as it always has.
    let cases = floor.n;
    drop(floor);
    print_observed::<S>(
        &[
            ("topology_min_nulls_agreed", &nulls_agreed),
            ("topology_carve_outs", &carve_outs_seen),
        ],
        cases,
    );
    assert!(
        failures.is_empty(),
        "{} of {} quantisations disagree with DGGAL and are not accounted for by the \
         characterised divergences in common/expected_divergences.rs: either no entry lists \
         the point, or its entry records other answers.\n\
         A new entry belongs there only once it has been characterised the way those \
         were: check whether DGGAL returns the same identifier for several distinct \
         points, whether its centroid is anywhere near the input, whether repeated calls \
         agree, and above all whether `getZoneFromTextID` parses back the text \
         `getZoneTextID` prints. If the engine's answer survives those checks, this port \
         is wrong and the fix belongs in its topology, not there.\nFirst few:\n{}",
        failures.len(),
        pts.len() * usize::from(S::I::MAX_RESOLUTION + 1),
        failures
            .iter()
            .take(10)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
    let m = measured::<S>();
    // The null path has to be exercised, not merely described.
    assert!(
        nulls_agreed >= m.topology_min_nulls_agreed,
        "the sample no longer reaches the engine's null zone (agreed on {nulls_agreed}); \
         for IGEO7, NULL_PROBE should give it at resolutions 15, 17 and 19"
    );
    // And the carve-outs have to stay live.
    assert_eq!(
        carve_outs_seen, m.topology_carve_outs,
        "expected exactly the measured carve-outs, for IGEO7 the resolution-18 divergence \
         at NULL_PROBE; saw {carve_outs_seen}"
    );
}
