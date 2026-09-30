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
//! zones at one level, in `sub_zone_methods_refuse_on_aperture_7`); or it lies inside one
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
//! On aperture 7, DGGAL's `getZoneParents` and `getZoneChildren` are its geometric
//! hierarchy, a different relation from the congruent Z7 digit path this crate
//! implements, and its `countSubZones` and relatives answer an order that the aperture-7
//! grids here refuse to define. None of them is used as an oracle there. The hierarchy is
//! checked through text identifiers instead, which are the digit path itself: see
//! `the_text_route_reads_digit_paths_as_the_hierarchy`. DGGAL offers `isZoneAncestorOf`
//! and `areZonesSiblings` (`dggrs.ec:352`, `:360`), but they answer over its geometric
//! hierarchy, not the congruent one this crate implements, so they cannot serve as the
//! oracle for `is_ancestor_of` and `is_sibling_of`, which are checked by the same text-id
//! route instead, in `ancestry_and_siblings_via_text_ids`.

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

/// The route by which the hierarchy is checked, verified before it is relied upon.
///
/// The claim is that DGGAL's `getZoneFromTextID` reads a Z7 text identifier as the digit
/// path it spells, so that a truncated text names the congruent parent and an extended
/// one a child. Three things are checked. The engine rejects what is not a path: a digit
/// 7, a base cell 12, and a deleted pentagon child. The engine prints back exactly the
/// path it read, here for each truncated text and in `check_identifier` for every zone
/// of the main and seam samples. And geometrically,
/// the zone the engine reads from a truncated text contains the child's centroid, in
/// the sense that the child's centroid lies nearer that zone's centroid than any of that
/// zone's neighbours': so a truncated text names the cell the child lies in, which is
/// what a parent is, and not merely some identifier.
pub fn the_text_route_reads_digit_paths_as_the_hierarchy<S: Subject<T = HexA7, I = Z7>>() {
    assert_eq!(
        o::zone_from_text(S::ORACLE, "00641567"),
        o::NULL_ZONE,
        "digit 7"
    );
    assert_eq!(
        o::zone_from_text(S::ORACLE, "12"),
        o::NULL_ZONE,
        "base cell 12"
    );
    assert_eq!(
        o::zone_from_text(S::ORACLE, "002"),
        o::NULL_ZONE,
        "north pentagon, digit 2"
    );
    assert_eq!(
        o::zone_from_text(S::ORACLE, "065"),
        o::NULL_ZONE,
        "south pentagon, digit 5"
    );
    assert_ne!(
        o::zone_from_text(S::ORACLE, "003"),
        o::NULL_ZONE,
        "north pentagon, digit 3"
    );

    let mut floor = Floor::new("text route", 2000);
    // The neighbours of a parent compared with the child's distance to it: counted, so
    // that an engine list that came back empty could not pass the comparison unseen.
    let mut neighbours_compared = 0usize;
    for (lat, lon) in sphere_points(SEED, 400) {
        for res in [1u8, 4, 7, 10, 12] {
            let child = o::zone_from_geo(S::ORACLE, lat, lon, i32::from(res));
            let text = o::text_id(S::ORACLE, child);
            let parent_text = &text[..text.len() - 1];
            let parent = o::zone_from_text(S::ORACLE, parent_text);
            assert_eq!(
                o::text_id(S::ORACLE, parent),
                parent_text,
                "DGGAL prints back {parent_text}"
            );
            let cc = o::centroid(S::ORACLE, child);
            let pc = o::centroid(S::ORACLE, parent);
            let to_parent = arc_deg(cc, pc);
            for n in o::neighbors(S::ORACLE, parent) {
                let to_other = arc_deg(cc, o::centroid(S::ORACLE, n));
                assert!(
                    to_parent < to_other,
                    "{text}: its centroid is nearer {} than {parent_text} (seed {SEED:#x})",
                    o::text_id(S::ORACLE, n)
                );
                neighbours_compared += 1;
            }
            floor.hit();
        }
    }
    eprintln!("text route: {neighbours_compared} neighbours of a parent compared");
    assert!(
        neighbours_compared > 0,
        "no neighbour of any parent was compared (seed {SEED:#x})"
    );
}

pub fn congruent_hierarchy_via_text_ids<S: Subject<T = HexA7, I = Z7>>() {
    let mut floor = Floor::new("hierarchy", 5000);
    for id in zones::<S>() {
        let z = S::grid().zone(id);
        let text = z.text_id();
        let res = z.resolution();
        match z.parent() {
            Some(p) => {
                assert_eq!(
                    p.id().0,
                    o::zone_from_text(S::ORACLE, &text[..text.len() - 1]),
                    "parent of {text} (seed {SEED:#x})"
                );
                assert_eq!(z.parents(), vec![p], "parents of {text} (seed {SEED:#x})");
                assert_eq!(
                    z.centroid_parent(),
                    Some(p),
                    "centroid parent of {text} (seed {SEED:#x})"
                );
                assert!(
                    z.is_centroid_child() && z.is_immediate_child_of(&p),
                    "{text} is not its parent's centroid child (seed {SEED:#x})"
                );
            }
            None => assert_eq!(res, 0, "{text} has no parent (seed {SEED:#x})"),
        }
        let children: Vec<u64> = z.children().iter().map(|c| c.id().0).collect();
        if id == ZoneId::NULL || res == S::grid().max_resolution() {
            assert!(children.is_empty(), "{text} has children (seed {SEED:#x})");
        } else {
            let omit = if text[..2].parse::<u8>().unwrap() <= 5 {
                2
            } else {
                5
            };
            let digits: Vec<u8> = (0..7).filter(|&d| !z.is_pentagon() || d != omit).collect();
            let want: Vec<u64> = digits
                .iter()
                .map(|d| o::zone_from_text(S::ORACLE, &format!("{text}{d}")))
                .collect();
            assert_eq!(children, want, "children of {text} (seed {SEED:#x})");
            assert!(
                !want.contains(&o::NULL_ZONE),
                "a child of {text} is null (seed {SEED:#x})"
            );
            if z.is_pentagon() {
                assert_eq!(
                    o::zone_from_text(S::ORACLE, &format!("{text}{omit}")),
                    o::NULL_ZONE,
                    "the deleted pentagon child of {text} (seed {SEED:#x})"
                );
            }
        }
        floor.hit();
    }
}

/// `is_ancestor_of` and `is_sibling_of` have no DGGAL counterpart for the congruent
/// hierarchy, so they are checked against a composition of the engine's answers, read
/// through its `getZoneFromTextID`, which takes a Z7 text as the digit path it spells
/// (see `the_text_route_reads_digit_paths_as_the_hierarchy`).
///
/// Ancestry is prefix order on the text. Every proper prefix of a zone's text names, in
/// the engine, a zone of which ours must say it is an ancestor, and the zone the engine
/// reads from the text extended by one digit is a descendant. The negatives come from
/// the engine as well: no neighbour of the zone's parent, in the engine's own list, is
/// an ancestor, nor is the zone itself.
///
/// Siblings are texts of the same length that agree but for their last digit. Every
/// child the engine reads under the zone's parent is a sibling but the zone itself,
/// which is not its own sibling, and under a pentagon the engine must refuse the
/// deleted child; each of the engine's neighbours of the zone is a sibling exactly
/// when its text shares the parent's prefix. A zone at resolution 0 has no parent and
/// so no sibling, itself included.
///
/// The null zone, which the sample holds at odd resolutions of 15 and above, has no
/// parents, so it is ancestor and sibling of nothing, and nothing is either to it.
pub fn ancestry_and_siblings_via_text_ids<S: Subject<T = HexA7, I = Z7>>() {
    let g = S::grid();
    let mut floor = Floor::new("ancestry and siblings", 1000);
    let root = g.zone_from_text("00").unwrap();
    let mut null_zone_checked = 0usize;
    for id in zones::<S>()
        .into_iter()
        .step_by(5)
        .chain(std::iter::once(ZoneId::NULL))
    {
        let z = g.zone(id);
        let text = z.text_id();
        if id == ZoneId::NULL {
            assert!(
                !z.is_ancestor_of(&root)
                    && !root.is_ancestor_of(&z)
                    && !z.is_sibling_of(&root)
                    && !root.is_sibling_of(&z)
                    && !z.is_sibling_of(&z),
                "the null zone is related to something (seed {SEED:#x})"
            );
            null_zone_checked += 1;
            floor.hit();
            continue;
        }
        let engine = |t: &str| g.zone(ZoneId(o::zone_from_text(S::ORACLE, t)));

        // Ancestors: every proper prefix, and nothing the other way round.
        for len in 2..text.len() {
            let prefix = &text[..len];
            let a = engine(prefix);
            assert_ne!(
                a.id().0,
                o::NULL_ZONE,
                "DGGAL reads no zone from {prefix} (seed {SEED:#x})"
            );
            assert!(
                a.is_ancestor_of(&z) && !z.is_ancestor_of(&a),
                "{prefix} is not an ancestor of {text} (seed {SEED:#x})"
            );
        }
        assert!(
            !z.is_ancestor_of(&z),
            "{text} is its own ancestor (seed {SEED:#x})"
        );
        if z.resolution() < g.max_resolution() {
            let child_text = format!("{text}0");
            let c = engine(&child_text);
            assert!(
                z.is_ancestor_of(&c) && !c.is_ancestor_of(&z),
                "{text} is not an ancestor of {child_text} (seed {SEED:#x})"
            );
        }

        let Some(prefix) = text.get(..text.len() - 1).filter(|p| p.len() >= 2) else {
            // Resolution 0: no parent, so no sibling, not even itself.
            assert!(
                !z.is_sibling_of(&z),
                "{text} at resolution 0 is its own sibling (seed {SEED:#x})"
            );
            for n in o::neighbors(S::ORACLE, id.0) {
                assert!(
                    !z.is_sibling_of(&g.zone(ZoneId(n))),
                    "{text} at resolution 0 has the sibling {} (seed {SEED:#x})",
                    o::text_id(S::ORACLE, n)
                );
            }
            floor.hit();
            continue;
        };

        // Not an ancestor: the engine's neighbours of the parent.
        for n in o::neighbors(S::ORACLE, o::zone_from_text(S::ORACLE, prefix)) {
            let nt = o::text_id(S::ORACLE, n);
            let expected = nt.len() < text.len() && text.starts_with(&nt);
            assert_eq!(
                g.zone(ZoneId(n)).is_ancestor_of(&z),
                expected,
                "is {nt}, a neighbour of {prefix}, an ancestor of {text} (seed {SEED:#x})"
            );
        }

        // Siblings: the parent's children, as the engine reads them, but the zone's own
        // digit, which names the zone itself and so is no sibling of it.
        let omit = if text[..2].parse::<u8>().unwrap() <= 5 {
            '2'
        } else {
            '5'
        };
        let parent_is_pentagon = prefix[2..].bytes().all(|b| b == b'0');
        for d in '0'..='6' {
            let st = format!("{prefix}{d}");
            let s = engine(&st);
            if parent_is_pentagon && d == omit {
                assert_eq!(
                    s.id().0,
                    o::NULL_ZONE,
                    "DGGAL reads the deleted pentagon child {st} (seed {SEED:#x})"
                );
                assert!(
                    !z.is_sibling_of(&s),
                    "{text} has the null zone as a sibling (seed {SEED:#x})"
                );
            } else if st == text {
                assert_eq!(
                    s.id(),
                    z.id(),
                    "the engine reads {st} as something other than {text} itself (seed {SEED:#x})"
                );
                assert!(
                    !z.is_sibling_of(&s),
                    "{text} is its own sibling (seed {SEED:#x})"
                );
            } else {
                assert_ne!(
                    s.id().0,
                    o::NULL_ZONE,
                    "DGGAL reads no zone from {st} (seed {SEED:#x})"
                );
                assert!(
                    z.is_sibling_of(&s) && s.is_sibling_of(&z),
                    "{st} is not a sibling of {text} (seed {SEED:#x})"
                );
            }
        }

        // Siblings or not: the engine's neighbours of the zone, by their text. The engine
        // sometimes lists the zone itself among its own neighbours (one of the kinds of
        // entry in `DROPPED`), which is excluded here too, since a zone is not its own
        // sibling.
        for n in o::neighbors(S::ORACLE, id.0) {
            let nt = o::text_id(S::ORACLE, n);
            let expected = n != id.0 && nt.len() == text.len() && nt.starts_with(prefix);
            assert_eq!(
                z.is_sibling_of(&g.zone(ZoneId(n))),
                expected,
                "is {nt}, a neighbour of {text}, a sibling (seed {SEED:#x})"
            );
        }
        floor.hit();
    }
    assert!(
        null_zone_checked > 0,
        "the null-zone branch never ran (seed {SEED:#x})"
    );
    eprintln!("ancestry and siblings: null zone checked {null_zone_checked} time(s)");
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
/// children, read through the engine's text route as in
/// `congruent_hierarchy_via_text_ids`. The sub-zone methods, the hierarchy predicates and
/// the disk, which the engine cannot judge, must equal the typed grid's answers.
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

        // The sub-zones, which the aperture-7 grids refuse, and the disk.
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
/// `any_grid_answers_as_the_typed_grid_and_the_engine`. On aperture 7 they are read through
/// the text route, as in `congruent_hierarchy_via_text_ids`: the parent is the zone the
/// engine reads from the text shortened by its last digit, and each child is a digit path
/// one digit longer that the engine reads as that very child. On aperture 3 they are the
/// engine's own lists, asked where [`engine_may_be_asked`] allows: the parent is the first
/// of its parents, and the children are its children in its order, but at the finest
/// resolution, where this crate answers none.
fn engine_parent_and_children<S: Subject>(
    id: ZoneId,
    parent: Option<ZoneId>,
    children: &[ZoneId],
    at: &str,
) {
    let text = S::grid().text_id(id);
    match S::T::APERTURE {
        7 => {
            if let Some(p) = parent {
                assert_eq!(
                    p.0,
                    o::zone_from_text(S::ORACLE, &text[..text.len() - 1]),
                    "parent of {at}, read by the engine from the truncated text"
                );
            }
            for &ch in children {
                let ct = S::grid().text_id(ch);
                assert!(
                    ct.len() == text.len() + 1 && ct.starts_with(&text),
                    "{ct} is not a digit path below {at}"
                );
                assert_eq!(
                    o::zone_from_text(S::ORACLE, &ct),
                    ch.0,
                    "child {ct} of {at}, read by the engine"
                );
            }
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

/// The aperture-7 grids define no sub-zone order, so the sub-zone methods refuse rather
/// than answer. The null zone is refused as an invalid zone, before the question of
/// order arises, by `sub_zone_index` as by the others.
///
/// `sub_zone_index` of two different zones at one level is a characterised divergence,
/// re-checked live at every zone of the sample that has a neighbour: the engine answers
/// 0, having compared only the two levels (`RI7H.ec:229-230`), and this crate answers
/// `None`, since at depth 0 a zone's only sub-zone is itself.
pub fn sub_zone_methods_refuse_on_aperture_7<S: Subject<T = HexA7, I = Z7>>() {
    let mut floor = Floor::new("sub-zones", 1000);
    let mut null_zone_checked = 0usize;
    // The two branches that not every zone reaches, counted so that neither can fall
    // silent: a first child to ask about, and a neighbour at the zone's own level.
    let mut children_checked = 0usize;
    let mut same_level_checked = 0usize;
    for id in zones::<S>()
        .into_iter()
        .step_by(5)
        .chain(std::iter::once(ZoneId::NULL))
    {
        let z = S::grid().zone(id);
        let text = z.text_id();
        let refused = |r: std::result::Result<(), Error>| match r {
            Err(Error::NoSubZoneOrder) => id != ZoneId::NULL,
            Err(Error::InvalidZone(_)) => id == ZoneId::NULL,
            _ => false,
        };
        for depth in [0u8, 1, 2, 19] {
            assert!(
                refused(z.count_sub_zones(depth).map(drop)),
                "count_sub_zones({depth}) of {text}"
            );
            assert!(
                refused(z.first_sub_zone(depth).map(drop)),
                "first_sub_zone({depth}) of {text}"
            );
            assert!(
                refused(z.sub_zones(depth).map(drop)),
                "sub_zones({depth}) of {text}"
            );
            assert!(
                refused(z.sub_zone_at_index(depth, 0).map(drop)),
                "sub_zone_at_index({depth}) of {text}"
            );
        }
        if let Some(c) = z.children().first() {
            assert_eq!(
                z.sub_zone_index(c),
                Err(Error::NoSubZoneOrder),
                "sub_zone_index of {text}"
            );
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
