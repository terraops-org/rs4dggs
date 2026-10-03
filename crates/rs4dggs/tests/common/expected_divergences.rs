//! Where this port answers differently from DGGAL and its answer is the defensible one.
//!
//! An entry is admitted only with evidence of what the engine does and of why our answer
//! is better, and the evidence is re-checked against the live engine on every run by
//! `listed_divergences_still_hold_with_their_evidence`, so that an entry cannot outlive
//! the behaviour it describes and quietly turn into documentation of a defect.
//! Disagreements for which neither side has been shown right do not belong here; they
//! are in `unresolved`.
//!
//! Each listed entry names its grid, and each class rule checks its evidence against the
//! grid under test, a `Subject`. Everything described below was measured on IGEO7, and
//! measured again on IVEA7H and on RTEA7H, which show the same phenomena at the same
//! places, with the same identifiers; the figures given are IGEO7's unless another grid's
//! are named. The strip that the rule for unreadable answers requires is one of the grid's
//! unresolved regions, so that rule admits nothing on a grid for which no such region is
//! recorded.
//!
//! There are two kinds.
//!
//! **The engine cannot read back its own answer.** Along two icosahedron edges, the one
//! from base cell 0 to base cell 1 over the north pole and the one from base cell 1 to
//! base cell 6, DGGAL's odd resolutions from 15 have a strip in which its own candidate
//! search finds no cell, and it answers its null zone; this port does the same, and the
//! two agree there. At the even resolutions 16 and 18 inside that strip the engine does
//! not answer null. The eC's `computeParents` has no early exit (RI7H_Z7.ec:289-294):
//! where `parent0` finds no parent it stores the null zone and keeps walking, and
//! `from7H` packs whatever those null entries yield. The result is one fixed identifier
//! per resolution, the same for points tens of degrees apart:
//! `012222222222222220` (`0x1492492492490fff`) at 16 and `01222222222222222220`
//! (`0x149249249249243f`) at 18, that is base cell 1 followed by the digit 2 repeated
//! and a final 0. Under a northern pentagon the digit 2 is the deleted child, so this
//! path names no zone: DGGAL's own `getZoneFromTextID` returns its null zone for the
//! very text `getZoneTextID` printed, this crate's parser rejects it as a non-canonical
//! pentagon path, and the engine places its centroid near (85.44, -88.60) at resolution
//! 18, several degrees from any of the inputs. Repeated calls return the same value, so
//! it is not a tie-break that could fall either way between runs. This port answers the
//! null zone instead, which is what the engine answers at the odd resolutions either
//! side, and which unlike the engine's value can be printed and read back.
//!
//! Five points were characterised this way when the null zone was added to the
//! topology protocol, all at resolution 18; they are listed below with a witness at
//! resolution 16, which that census did not reach. The oracle suite meets the same
//! phenomenon at many more points along the two edges, and those are admitted by
//! [`engine_cannot_read_its_own_answer`], which checks the same evidence live at each
//! point rather than trusting a list. IVEA7H and RTEA7H answer the same two identifiers
//! in the same strip; two of the points are listed for each, one at each resolution.
//!
//! At the outer margin of the strip, some 5e-5 to 7e-5 degrees from the edge, the engine
//! at resolution 19 answers, at some points, an identifier it cannot read back as well,
//! `012222222222222222201` or `012222222222222222204`, where this port answers null.
//! This was found on all three grids by a sweep across the two edges, and at none of the
//! points this suite samples. The rule does not admit it, having no witness here, so a
//! sample that met it would fail and call for its characterisation.
//!
//! **Boundary ties by symmetry.** At exact singularities, the poles and the vertex-0
//! meridian among them, an independent port and the engine can break a tie differently,
//! and the engine's own choice there has been seen to change when reference tables were
//! regenerated from it, so such a point is admitted by the evidence of the tie rather
//! than by a fixed answer. The grid has two kinds of symmetry, and they differ by
//! parity.
//!
//! - At even resolutions the grid is mirror-symmetric about three planes: the vertex-0
//!   meridian (11.2 and -168.8), which holds the edge from base cell 0 to 1 and so
//!   passes over the north pole; the meridian at right angles to it (101.2 and -78.8);
//!   and the equator. Where a cell boundary lies in one of these planes, the cells
//!   either side are mirror images, and every point of the boundary is an exact tie.
//! - At odd resolutions the grid is chiral, and no reflection holds. A half-turn about
//!   the axis where two of the planes meet does hold, at every resolution. Those axes
//!   meet the sphere at six points: the two poles, (0, 11.2), (0, -168.8), (0, 101.2)
//!   and (0, -78.8). At those points alone, an odd resolution has ties, and they are
//!   ties by rotation, not reflection.
//!
//! Both were measured while this suite was written, on centroids of cells around each
//! of the six points, carried by the symmetry and quantised again by the engine:
//! - under the half-turns every centroid returned to a cell's centroid within 2.4e-13
//!   degrees, at every resolution from 0 to 11;
//! - under the reflections the same held at the even resolutions, within 2.6e-13;
//! - at the odd resolutions the reflections missed by 10.6, 1.46, 0.21, 0.030, 4.3e-3
//!   and 6.1e-4 degrees at resolutions 1 to 11.
//!
//! The same measurement on IVEA7H gives the same picture: the half-turns hold within
//! 2.3e-13 degrees at every resolution from 0 to 11, the reflections hold within 2.4e-13
//! at the even ones, and at the odd ones they miss by 10.5, 1.47, 0.21, 0.030, 4.3e-3
//! and 6.1e-4 degrees. So does it on RTEA7H: the half-turns hold within 2.1e-13 degrees,
//! the reflections within 2.4e-13 at the even resolutions, and at the odd ones they miss
//! by 10.4, 1.48, 0.21, 0.030, 4.3e-3 and 6.1e-4 degrees.
//!
//! The engine and this port broke some of these ties differently, from resolution 0
//! upwards, and at a few one side's candidate search accepted neither cell and answered
//! null. **They no longer do, in any sample the suite takes, on any of the three
//! grids**: the ties were decided by the last bits of the planar point, and once the
//! projection was put into the order the compiled engine performs it, this port's planar
//! points became the engine's own and both sides name the same cell at every one of them.
//! The four witnesses once listed for the class were removed then, and the suite now
//! holds its admissions at nothing; the rule is kept, and so is the account of it below,
//! so that a tie that reappears, under another build of the engine for instance, is read
//! against its history rather than admitted. The engine also places a boundary on a
//! mirror plane only to within about 1e-9 degrees, so an input a hair to one side of the
//! plane can land in the other cell there. The rule is
//! [`boundary_tie_on_a_mirror_plane`], and its evidence is decisive without any appeal to
//! geometry:
//! - a small step towards our cell, the engine answers our cell;
//! - the same step towards the engine's cell, we answer the engine's;
//! - where one side answered null, the same step away from the real answer, both sides
//!   answer the same other cell instead.
//!
//! The real answer's boundary therefore passes within that step of the input, by the
//! reckoning of both implementations.
//!
//! **Input this crate refuses.** A third kind is not a disagreement over an answer but
//! the refusal of one. Where the engine answers input that has no meaning with an
//! answer that looks like one, this crate returns an error, deliberately, and documents
//! the departure where a caller meets it. There are three such inputs, and each is
//! re-checked against the live engine on every run, by
//! `refusals_depart_from_the_engine_as_recorded` and, for the null zone's sub-zone
//! index, by `sub_zones_at_depth_1_as_the_engine_lists_them`:
//! - a latitude or longitude that is NaN or infinite, which the engine answers with the
//!   pentagon of base cell 1 at the requested resolution, whatever the input
//!   ([`engine_answer_to_a_non_finite_coordinate`]);
//! - the text of a level-20 identifier, twenty digits after the base cell, which the
//!   engine reads as a zone with no geometry, centred at (0, 0), and which this crate
//!   refuses as longer than the finest resolution allows ([`level_20_texts`]);
//! - `sub_zone_index` given an identifier with no cell behind it, which the engine
//!   answers without looking at the identifier, 0 within itself and -1 against a real
//!   zone, and which this crate refuses as an invalid zone ([`FABRICATED`]).
//!
//! The fourth refusal, of a `GridConfig` with a non-finite angle, has no counterpart in
//! the engine, which fixes each grid's orientation in its own source (`ri5x6.ec:106`)
//! and takes none from a caller. There is no answer of the engine's to re-check, and it
//! is covered by the crate's unit tests alone.
//!
//! **An answer the engine gives without looking.** Asked for the sub-zone index of one
//! zone within another of the same level, the engine compares the two levels and answers
//! 0 before it looks at either zone (`getSubZoneIndex`, `RI7H.ec:229-230`), so that it
//! reports a zone's neighbour as its sub-zone at index 0. This crate answers `None`,
//! since at depth 0 a zone's only sub-zone is the zone itself. The divergence is
//! re-checked against the live engine by `sub_zones_at_depth_1_as_the_engine_lists_them`, for a
//! neighbour of every zone of its sample that has one.
//!
//! **On aperture 3 the engine cannot read back its own answer along two seams.** A point
//! within about `1e-6 / 3^(res / 2)` of either seam of the planar layout is re-anchored
//! twice by the engine's quantiser (`RI3H.ec:1374-1413`) and leaves it with the root of a
//! polar pentagon, 10 or 11, without the modulo that the quantiser's later branches apply,
//! and with a non-zero index, which no zone has there. Root 10 comes along the icosahedron
//! edge from base cell 0 to base cell 1, over the north pole, and root 11 along the edge
//! from base cell 6 to base cell 11, which ends at the southern polar vertex, at every
//! resolution from 2 to 21. The engine's own `getZoneFromTextID` refuses the text it
//! prints for such an identifier, and it draws it as the polar pentagon of its root, up to
//! some 63 degrees from the input. This crate answers the null zone there. The rule is
//! [`engine_answers_an_unreadable_polar_root`]; four points are listed as its witnesses on
//! each aperture-3 grid, one on each seam at resolution 2 and one on each at resolution
//! 21. The phenomenon was measured on ISEA3H, IVEA3H and RTEA3H, each over its own samples,
//! and is the same on the three: the same two seams, the same resolutions, and a reach from
//! the edge within the same bound. No tie has been met on aperture 3, where the tie rule is
//! not offered.
//!
//! The null zone departs on aperture 3 as well: the engine computes a centroid and six
//! vertices for it from the bits of its identifier, where this crate answers the point
//! (0, 0) and six zero vertices, as on aperture 7, whose engine does the same. The
//! engine's side is re-checked by [`engine_null_zone_geometry`], wherever the null zone
//! is met.
//!
//! So do the aperture-3 engine's sub-zone orders at some zones on the edges of the rhombi,
//! where a unit in the last place carries a point across an interruption of the layout that it
//! should only reach, and the engine's order names zones twice, names zones that are not
//! descendants of the zone and omits descendants. This crate's order there is the zone's
//! descendants, each once, in the generator's own scanline order. At every order in which the
//! two differ, [`engine_sub_zone_order_is_faulty`] re-checks both sides by the engine's own
//! tests: the engine's order fails them, and this crate's passes them.
//!
//! **The sub-zone order of the aperture-7 grids** is the engine's own, entry for entry, and
//! departs from it in five ways, each this crate's stated answer and each re-checked against
//! the live engine by `common::sub_zone_order`:
//! - an entry that the engine gives as an identifier which it cannot read back is the null
//!   zone here, at the same place, as an entry that is the engine's null zone is; this is
//!   met at the broken seams alone, from sub-zones of level 15, and counted by level;
//! - the first sub-zone at depth 0 is the zone itself, where the engine's `getFirstSubZone`
//!   names a neighbour of the zone, or its null zone, at all but a few zones;
//! - the zone at an index is the entry of the order at that index. The engine's
//!   `getSubZoneAtIndex` answers another zone at a pentagon at an odd depth, from the
//!   scanline after the pentagon's own; it is not asked here, since the vendored binding
//!   declares the call without its depth, and the comparison is with the engine's order;
//! - the first sub-zone of the pentagon of the south pole at an odd level, at depth 2, is
//!   the first entry of the engine's order, where the engine's `getFirstSubZone` ends the
//!   process ([`engine_first_sub_zone_ends_the_process`]): the call is not made there, and
//!   the orders it is kept from are counted;
//! - the index of a sub-zone is the engine's `getSubZoneIndex` wherever the order holds the
//!   zone at that index, and none otherwise. At the broken seams the engine's walk, which
//!   compares centroids, answers for some zones an index at which its own order names
//!   another zone; such a zone has no index here, and each is counted by level. Where the
//!   engine answers -1 for a zone that its order names, there is none here either, and
//!   those are counted too. The index is not sought in an order longer than a list may be,
//!   where the engine walks an order of any length; the suites ask none so long.
//!
//! No order longer than [`LONGEST_ORDER_ASKED`] is asked of the engine, which ends the
//! process at `2^28` sub-zones on these grids.
//!
//! **Identifiers this crate has no geometry for.** A zone without geometry (the null zone, a Z7
//! identifier of twenty digits, an identifier that no cell of the grid has behind it, and on
//! aperture 3 the sub-hexes C and D of a polar root, which the engine reads and this crate takes
//! for no cell) answers the empty ring, no extent and no area. The engine answers
//! each in its own way, and [`NoGeometry`] names the four: the null zone (the cleared extent, the
//! area `+inf`); an identifier it does not draw (the cleared extent, no ring, a finite area);
//! the same with the area `+inf`, where its count of zones leaves a divisor of zero; and the
//! polar sub-hexes, which it draws, with a ring, a finite extent and a finite area, so that all
//! three answers depart. Each is classified by [`engine_answers_without_geometry`], which
//! re-checks the engine's side live at every instance, and the instances are counted exactly, per
//! grid, by `geometry::the_refined_boundary_and_the_extent_are_the_engines`, which asserts that
//! each is met and that the three answers of this crate are the empty ones. On aperture 3 the
//! ring of an identifier the engine cannot read is not asked of it: its answer is a null array,
//! on which the vendored binding crashes.
//!
//! **The zones of a box.** `zones_in_box` yields every zone of the level whose extent meets the
//! box, by the engine's own test of two extents, where the engine's `listZones` leaves some of
//! them out. On the aperture-3 grids the engine looks for zones about a sample of points of the
//! box, and misses a zone that belongs to another rhombus than the one its samples fall in: the
//! one zone of level 2 that holds the whole box about Lisbon, and every zone of a small cap
//! about a pole, for which it returns nothing at all. On the aperture-7 grids it descends from
//! the twelve base cells, and was measured to miss one zone in six boxes of some forty
//! thousand; the suites' boxes meet none. Wherever the two answers differ,
//! [`engine_omits_zones_of_a_box`] re-checks the evidence against the live engine: its answer
//! is this crate's less some zones, in the same order, and each zone it leaves out meets the
//! box by the engine's own extent of it, so that by the engine's own rule the zone belongs
//! there. The instances are counted exactly, per grid and by kind of box, by
//! `zones::the_zones_of_a_box_are_the_engines_and_those_it_leaves_out`, and the boxes known by
//! name by `zones::the_boxes_known_by_name_depart_from_the_engine_as_recorded`.

use dggal_oracle as o;
use rs4dggs::ZoneId;

use super::boxes;
use super::subject::Subject;
use super::unresolved::regions;
use super::{arc_deg, dist_to_arc_deg, icosahedron_vertices, nudge};

/// The kind of evidence an entry rests on; the live check differs by kind.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Evidence {
    /// DGGAL answers an identifier that its own `getZoneFromTextID` cannot read back, and
    /// this port answers the null zone.
    EngineCannotReadItsOwnAnswer,
    /// The input lies on a cell boundary to within rounding, on a mirror plane at an even
    /// resolution or at a half-turn axis point at an odd one; the engine names one cell
    /// and this port the other. No entry rests on it at present: the four that did ceased
    /// to diverge when the projection was ported, and were removed. Its live check is kept
    /// for a new one.
    BoundaryTie,
}

/// One characterised point at which this port and DGGAL disagree.
#[derive(Debug)]
pub struct Divergence {
    pub grid: &'static str,
    pub lat: f64,
    pub lon: f64,
    pub res: u8,
    /// Our text identifier at this input.
    pub ours: &'static str,
    /// DGGAL's text identifier at this input.
    pub dggal: &'static str,
    pub evidence: Evidence,
    pub reason: &'static str,
}

/// The engine's answer at resolution 18 inside the strip, which it cannot read back.
pub const UNREADABLE_18: &str = "01222222222222222220";
/// The engine's answer at resolution 16 inside the strip, which it cannot read back.
pub const UNREADABLE_16: &str = "012222222222222220";

/// The characterised points. The first five are those the null-zone work found; each
/// is 0.001 to 0.23 degrees from the north pole, on the polar edge.
pub static DIVERGENCES: &[Divergence] = &[
    Divergence {
        grid: "IGEO7",
        lat: 89.999_129_147_255_16,
        lon: 9.333_761_199_597_859,
        res: 18,
        ours: rs4dggs::NULL_TEXT,
        dggal: UNREADABLE_18,
        evidence: Evidence::EngineCannotReadItsOwnAnswer,
        reason: "the null-zone probe: null at 15, 17 and 19 on both sides, the unreadable \
                 identifier from the engine at 18",
    },
    Divergence {
        grid: "IGEO7",
        lat: 89.990_330_249_427_36,
        lon: 11.086_742_952_786_494,
        res: 18,
        ours: rs4dggs::NULL_TEXT,
        dggal: UNREADABLE_18,
        evidence: Evidence::EngineCannotReadItsOwnAnswer,
        reason: "same null-zone parent walk, same identifier",
    },
    Divergence {
        grid: "IGEO7",
        lat: 89.896_274_363_623_95,
        lon: -168.805_578_920_963_94,
        res: 18,
        ours: rs4dggs::NULL_TEXT,
        dggal: UNREADABLE_18,
        evidence: Evidence::EngineCannotReadItsOwnAnswer,
        reason: "same null-zone parent walk, same identifier",
    },
    Divergence {
        grid: "IGEO7",
        lat: 89.949_875_972_429_51,
        lon: 11.218_898_183_735_66,
        res: 18,
        ours: rs4dggs::NULL_TEXT,
        dggal: UNREADABLE_18,
        evidence: Evidence::EngineCannotReadItsOwnAnswer,
        reason: "same null-zone parent walk, same identifier",
    },
    Divergence {
        grid: "IGEO7",
        lat: 89.778_410_001_068_01,
        lon: -168.796_000_793_232_4,
        res: 18,
        ours: rs4dggs::NULL_TEXT,
        dggal: UNREADABLE_18,
        evidence: Evidence::EngineCannotReadItsOwnAnswer,
        reason: "same null-zone parent walk, same identifier",
    },
    Divergence {
        grid: "IGEO7",
        lat: 89.9975,
        lon: -170.0,
        res: 16,
        ours: rs4dggs::NULL_TEXT,
        dggal: UNREADABLE_16,
        evidence: Evidence::EngineCannotReadItsOwnAnswer,
        reason: "the same parent walk one even level coarser, found by a polar sweep: \
                 null at 15 and 17 on both sides, the unreadable identifier from the \
                 engine at 16",
    },
    // IVEA7H. Of the IGEO7 witnesses above, four of the five at resolution 18 and the one
    // at 16 are witnesses on this grid too, while at (89.778, -168.796) IVEA7H's engine
    // answers a real cell at 18 and agrees with this port. One of each resolution is
    // listed here; the class rule admits the rest.
    Divergence {
        grid: "IVEA7H",
        lat: 89.999_129_147_255_16,
        lon: 9.333_761_199_597_859,
        res: 18,
        ours: rs4dggs::NULL_TEXT,
        dggal: UNREADABLE_18,
        evidence: Evidence::EngineCannotReadItsOwnAnswer,
        reason: "the null-zone probe on IVEA7H: null at 15, 17 and 19 on both sides, the \
                 unreadable identifier from the engine at 18, as on IGEO7",
    },
    Divergence {
        grid: "IVEA7H",
        lat: 89.9975,
        lon: -170.0,
        res: 16,
        ours: rs4dggs::NULL_TEXT,
        dggal: UNREADABLE_16,
        evidence: Evidence::EngineCannotReadItsOwnAnswer,
        reason: "the resolution-16 witness on IVEA7H: null at 15 and 17 on both sides, the \
                 unreadable identifier from the engine at 16, as on IGEO7",
    },
    // RTEA7H. Of the IGEO7 witnesses above, the same four of the five at resolution 18 and
    // the one at 16 are witnesses on this grid too, while at (89.778, -168.796) RTEA7H's
    // engine answers a real cell at 18 and agrees with this port, as IVEA7H's does. One of
    // each resolution is listed here; the class rule admits the rest.
    Divergence {
        grid: "RTEA7H",
        lat: 89.999_129_147_255_16,
        lon: 9.333_761_199_597_859,
        res: 18,
        ours: rs4dggs::NULL_TEXT,
        dggal: UNREADABLE_18,
        evidence: Evidence::EngineCannotReadItsOwnAnswer,
        reason: "the null-zone probe on RTEA7H: null at 15, 17 and 19 on both sides, the \
                 unreadable identifier from the engine at 18, as on IGEO7",
    },
    Divergence {
        grid: "RTEA7H",
        lat: 89.9975,
        lon: -170.0,
        res: 16,
        ours: rs4dggs::NULL_TEXT,
        dggal: UNREADABLE_16,
        evidence: Evidence::EngineCannotReadItsOwnAnswer,
        reason: "the resolution-16 witness on RTEA7H: null at 15 and 17 on both sides, the \
                 unreadable identifier from the engine at 16, as on IGEO7",
    },
    // ISEA3H: one witness on each seam at the coarsest resolution the rule admits and one at
    // the finest (see `engine_answers_an_unreadable_polar_root`); the rule admits the rest.
    Divergence {
        grid: "ISEA3H",
        lat: 89.999_979_814_075_29,
        lon: 101.199_999_744_738_58,
        res: 2,
        ours: rs4dggs::NULL_TEXT,
        dggal: "BA-2-A",
        evidence: Evidence::EngineCannotReadItsOwnAnswer,
        reason: "2e-5 degrees from the north pole, on the edge from base cell 0 to 1: root \
                 10 with index 2, drawn as the pentagon at 58.4 N, 31.6 degrees away",
    },
    Divergence {
        grid: "ISEA3H",
        lat: -27.777_308_938_516_747,
        lon: -146.892_522_341_526_25,
        res: 2,
        ours: rs4dggs::NULL_TEXT,
        dggal: "BB-6-A",
        evidence: Evidence::EngineCannotReadItsOwnAnswer,
        reason: "on the edge from base cell 6 to 11: root 11 with index 6, drawn as the \
                 pentagon at 58.4 S, 34.2 degrees away",
    },
    Divergence {
        grid: "ISEA3H",
        lat: 80.308_344_729_849_85,
        lon: 11.200_000_006_396_612,
        res: 21,
        ours: rs4dggs::NULL_TEXT,
        dggal: "KA-95EE-B",
        evidence: Evidence::EngineCannotReadItsOwnAnswer,
        reason: "on the edge from base cell 0 to 1 at the finest resolution the rule \
                 admits, 21.9 degrees from the pentagon the engine draws",
    },
    Divergence {
        grid: "ISEA3H",
        lat: -30.111_251_718_264_086,
        lon: -147.894_842_553_269_26,
        res: 21,
        ours: rs4dggs::NULL_TEXT,
        dggal: "KB-67EA811D-B",
        evidence: Evidence::EngineCannotReadItsOwnAnswer,
        reason: "on the edge from base cell 6 to 11 at the finest resolution the rule \
                 admits, 31.7 degrees from the pentagon the engine draws",
    },
    // IVEA3H and RTEA3H: on each, the same point a hair from the north pole as on ISEA3H at
    // resolution 2, on the seam over the pole, and three points of the sample of
    // `unreadable_answers_along_two_seams_are_the_null_zone` lifted by the grid's own
    // projection, the farthest there from the pentagon the engine draws: one on the southern
    // seam at resolution 2, and one on each seam at 21, the finest resolution the rule
    // admits.
    Divergence {
        grid: "IVEA3H",
        lat: 89.999_979_814_075_29,
        lon: 101.199_999_744_738_58,
        res: 2,
        ours: rs4dggs::NULL_TEXT,
        dggal: "BA-2-A",
        evidence: Evidence::EngineCannotReadItsOwnAnswer,
        reason: "2e-5 degrees from the north pole, on the edge from base cell 0 to 1: root \
                 10 with index 2, drawn as the pentagon at 58.4 N, 31.6 degrees away",
    },
    Divergence {
        grid: "IVEA3H",
        lat: -2.991_702_715_904_28,
        lon: -138.051_193_378_182_03,
        res: 2,
        ours: rs4dggs::NULL_TEXT,
        dggal: "BB-3-A",
        evidence: Evidence::EngineCannotReadItsOwnAnswer,
        reason: "on the edge from base cell 6 to 11: root 11 with index 3, drawn as the \
                 pentagon at 58.4 S, 60.4 degrees away",
    },
    Divergence {
        grid: "IVEA3H",
        lat: 61.521_976_349_644_866,
        lon: -168.800_000_002_330_53,
        res: 21,
        ours: rs4dggs::NULL_TEXT,
        dggal: "KA-B89-B",
        evidence: Evidence::EngineCannotReadItsOwnAnswer,
        reason: "on the edge from base cell 0 to 1 at the finest resolution the rule \
                 admits, 60.1 degrees from the pentagon the engine draws",
    },
    Divergence {
        grid: "IVEA3H",
        lat: -2.991_719_510_762_800_7,
        lon: -138.051_175_715_144_58,
        res: 21,
        ours: rs4dggs::NULL_TEXT,
        dggal: "KB-A64B371-B",
        evidence: Evidence::EngineCannotReadItsOwnAnswer,
        reason: "on the edge from base cell 6 to 11 at the finest resolution the rule \
                 admits, 60.4 degrees from the pentagon the engine draws",
    },
    Divergence {
        grid: "RTEA3H",
        lat: 89.999_979_814_075_29,
        lon: 101.199_999_744_738_58,
        res: 2,
        ours: rs4dggs::NULL_TEXT,
        dggal: "BA-2-A",
        evidence: Evidence::EngineCannotReadItsOwnAnswer,
        reason: "2e-5 degrees from the north pole, on the edge from base cell 0 to 1: root \
                 10 with index 2, drawn as the pentagon at 58.4 N, 31.6 degrees away",
    },
    Divergence {
        grid: "RTEA3H",
        lat: -3.097_648_492_855_740_6,
        lon: -138.085_553_524_148_62,
        res: 2,
        ours: rs4dggs::NULL_TEXT,
        dggal: "BB-3-A",
        evidence: Evidence::EngineCannotReadItsOwnAnswer,
        reason: "on the edge from base cell 6 to 11: root 11 with index 3, drawn as the \
                 pentagon at 58.4 S, 60.3 degrees away",
    },
    Divergence {
        grid: "RTEA3H",
        lat: 74.360_572_371_260_89,
        lon: -168.800_000_004_082_05,
        res: 21,
        ours: rs4dggs::NULL_TEXT,
        dggal: "KA-39AB-B",
        evidence: Evidence::EngineCannotReadItsOwnAnswer,
        reason: "on the edge from base cell 0 to 1 at the finest resolution the rule \
                 admits, 47.2 degrees from the pentagon the engine draws",
    },
    Divergence {
        grid: "RTEA3H",
        lat: -15.273_306_569_463_816,
        lon: -142.149_991_871_530_8,
        res: 21,
        ours: rs4dggs::NULL_TEXT,
        dggal: "KB-33F5B3E3-B",
        evidence: Evidence::EngineCannotReadItsOwnAnswer,
        reason: "on the edge from base cell 6 to 11 at the finest resolution the rule \
                 admits, 47.5 degrees from the pentagon the engine draws",
    },
];

/// The listed entry for this grid, input and resolution, if there is one.
pub fn lookup(grid: &str, lat: f64, lon: f64, res: u8) -> Option<&'static Divergence> {
    DIVERGENCES
        .iter()
        .find(|d| d.grid == grid && d.lat == lat && d.lon == lon && d.res == res)
}

/// How far from the two affected edges the strip reaches, in degrees.
///
/// Measured again on each grid on 2026-09-23, over a sweep of both edges at resolutions
/// 13 to 19: 24 stations along each edge, and at each station every offset across the
/// edge out to 1e-3 degrees either side in steps of 1e-6, that is 672,336 points a grid.
/// The engine's own null answers lie within 7.4e-5 degrees of the edge on IGEO7 and
/// within 7.5e-5 on IVEA7H and RTEA7H, and the inputs at which the two sides answer
/// differently within 6.7e-5, 6.5e-5 and 6.5e-5. Resolutions 13 and 14 show neither.
///
/// The zones of the oracle samples whose neighbours differ reach a little further, their
/// centroids lying within 7.7e-5 degrees of an edge on IGEO7 and on IVEA7H and within
/// 9.1e-5 on RTEA7H.
///
/// Measured again on 2026-09-26, after the projection was put into the order the compiled
/// engine performs it, on every disagreement the band admits in the three suites: the
/// zones whose neighbour sets differ lie within 7.67e-5 degrees of an edge on IGEO7,
/// 7.66e-5 on IVEA7H and 9.09e-5 on RTEA7H (293, 283 and 294 admissions, over the three
/// neighbour tests); the inputs admitted as answers the engine cannot read back lie within
/// 5.73e-5 on each grid (55, 46 and 42, over the main and seam samples and the listed
/// witnesses); every disk that differs, at radius 1 or 2, holds a zone within 3e-15 of an
/// edge; and no quantisation is admitted as one that neither answer contains. The reach
/// is the one measured before the port, unchanged. The bound, which was 2e-4, is 1.5e-4:
/// 1.65 times the widest of these reaches, RTEA7H's 9.09e-5, and twice the strip's own
/// width of 7.5e-5.
///
/// Measured again on 2026-09-28, on the deep digit paths under the polar pentagons that
/// `deep_polar_digit_paths_answer_as_the_engine` samples: the 460 zones under base cell
/// 0 whose neighbour sets differ, the same zones on every grid, lie within 9.45e-5
/// degrees of an edge on IGEO7, 9.08e-5 on IVEA7H and 8.77e-5 on RTEA7H. For 23 of them
/// on IGEO7 and 21 on the other two the nearest point is not on the edge from base cell
/// 0 to 1 but at its end, the pentagon of base cell 0, within 5.84e-5 degrees of it, so
/// that the band is there a disc about the pentagon. The widest reach is now IGEO7's
/// 9.45e-5, and the bound of 1.5e-4 is 1.59 times it and still covers it.
///
/// Since 2026-09-29 the band admits no neighbour list and no disk: the aperture-7 grids
/// list neighbours by the engine's own algorithm, and the suites compare every list with
/// the engine's less four kinds of entry, wherever the zone lies. The reaches of the
/// differing neighbour sets above are therefore a record of what the band once had to
/// cover. It now bounds the two quantisation rules alone, whose inputs lie within the
/// strip's own width of 7.5e-5 degrees and, for the answers the engine cannot read back,
/// within 5.73e-5; the bound is left as it stands, twice that width.
///
/// Since 2026-10-02 the band bounds a third thing: every zone that the comparison of the
/// aperture-7 hierarchy counts in a class of the broken seams (an entry left out of its
/// parents or of its children, a child that names other parents, no parent, and the rest)
/// must have its centroid within it. Measured then: among the zones that positions give,
/// within 1.11e-5 degrees of an edge on IGEO7, 1.10e-5 on IVEA7H and 1.12e-5 on RTEA7H
/// (49, 39 and 38 zones); among the identifiers built by text under base cell 0, at levels
/// 14 to 19, within 8.28e-5 on IGEO7, 7.98e-5 on IVEA7H and 7.67e-5 on RTEA7H (329 zones
/// on each), and the same reaches over every such identifier with a tail of five digits.
pub const BROKEN_SEAM_BAND_DEG: f64 = 1.5e-4;

/// The two icosahedron edges along which the engine's odd resolutions from 15 have a
/// strip of no cells, as pairs of base cells: 0 to 1 runs over the north pole along the
/// vertex-0 meridian and its antimeridian, 1 to 6 runs south-east from base cell 1. The
/// other twenty-eight edges were swept at resolutions 13 to 19 and show neither null
/// answers nor disagreements.
pub const BROKEN_SEAMS: [(usize, usize); 2] = [(0, 1), (1, 6)];

/// Whether a point lies in the strip of one of the grid's unresolved regions, at
/// whatever resolution: for IGEO7, within [`BROKEN_SEAM_BAND_DEG`] of one of
/// [`BROKEN_SEAMS`]. On a grid with no region recorded, nowhere.
pub fn near_broken_seam<S: Subject>(lat: f64, lon: f64) -> bool {
    regions::<S>().iter().any(|r| r.near(lat, lon))
}

/// The class rule of which the listed [`Evidence::EngineCannotReadItsOwnAnswer`] entries
/// are witnesses. It admits a disagreement only when all of the evidence holds at this
/// very input, checked against the live engine now: our answer is the null zone; the
/// engine's is not; the engine cannot read back the text it prints for its own answer,
/// and neither can this crate; the resolution is 16 or 18; and the input lies in the
/// strip along one of the two affected edges. Anything else is returned as an error for
/// the caller to fail on, since it is either a regression or a new phenomenon that needs
/// its own investigation.
///
/// Makes its oracle calls one after another, never one inside another.
pub fn engine_cannot_read_its_own_answer<S: Subject>(
    lat: f64,
    lon: f64,
    res: u8,
    ours: u64,
    engine: u64,
) -> Result<(), String> {
    let text = o::text_id(S::ORACLE, engine);
    let mut missing = Vec::new();
    if ours != o::NULL_ZONE {
        missing.push("our answer is not the null zone");
    }
    if engine == o::NULL_ZONE {
        missing.push("the engine's answer is the null zone");
    }
    if o::zone_from_text(S::ORACLE, &text) != o::NULL_ZONE {
        missing.push("the engine reads its own answer back");
    }
    if S::grid().zone_from_text(&text).is_ok() {
        missing.push("this crate reads the engine's answer");
    }
    if res != 16 && res != 18 {
        missing.push("the resolution is neither 16 nor 18");
    }
    if !near_broken_seam::<S>(lat, lon) {
        missing.push("the input is not in the strip along either affected edge");
    }
    if missing.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "engine {text} ({engine:#018x}), ours {ours:#018x}: {}",
            missing.join("; ")
        ))
    }
}

/// The unit normals of the three planes about which this orientation of the icosahedron
/// is mirror-symmetric: the vertex-0 meridian (11.2 east, continuing as -168.8), which
/// holds the edge from base cell 0 to 1 and so passes over the north pole; the meridian
/// at right angles to it (101.2, continuing as -78.8); and the equator. The twelve
/// vertices pair off across each of them: across the equator, for instance, base cell 0
/// with 8, 1 with 11, 2 with 7 and 5 with 9. 11.2 is `-180 - orientation_lon`, the same
/// `lon1` from which the eC's `fixPoles` takes its longitudes.
pub fn mirror_normals() -> [[f64; 3]; 3] {
    let (a, b) = (11.2f64.to_radians(), 101.2f64.to_radians());
    [
        [-a.sin(), a.cos(), 0.0],
        [-b.sin(), b.cos(), 0.0],
        [0.0, 0.0, 1.0],
    ]
}

/// How near a mirror plane an input must lie for [`boundary_tie_on_a_mirror_plane`] to
/// consider it, in degrees: about a millimetre.
pub const MIRROR_BAND_DEG: f64 = 1e-8;

/// Whether a point lies within [`MIRROR_BAND_DEG`] of one of the three mirror planes.
pub fn on_a_mirror_plane(lat: f64, lon: f64) -> bool {
    let p = super::unit_vector(lat, lon);
    mirror_normals().iter().any(|n| {
        let off = (p[0] * n[0] + p[1] * n[1] + p[2] * n[2])
            .abs()
            .min(1.0)
            .asin();
        off.to_degrees() <= MIRROR_BAND_DEG
    })
}

/// The six points where the axes of the grid's half-turns meet the sphere: the lines in
/// which the three mirror planes meet one another. Each is the midpoint of an
/// icosahedron edge: the north pole of the edge from base cell 0 to 1, the south pole of
/// 8 to 11, (0, 11.2) of 3 to 4, (0, -168.8) of 6 to 10, (0, 101.2) of 5 to 9 and
/// (0, -78.8) of 2 to 7. A half-turn about any of them carries the grid onto itself at
/// every resolution, odd ones included; see the module note.
pub const AXIS_POINTS: [(f64, f64); 6] = [
    (90.0, 0.0),
    (-90.0, 0.0),
    (0.0, 11.2),
    (0.0, -168.8),
    (0.0, 101.2),
    (0.0, -78.8),
];

/// Whether a point lies within [`MIRROR_BAND_DEG`] of one of the [`AXIS_POINTS`].
pub fn at_an_axis_point(lat: f64, lon: f64) -> bool {
    AXIS_POINTS
        .iter()
        .any(|&p| arc_deg((lat, lon), p) <= MIRROR_BAND_DEG)
}

/// The fraction of our cell's circumradius by which [`boundary_tie_on_a_mirror_plane`]
/// moves the input to each side, unless that is less than [`TIE_FLOOR_DEG`].
pub const TIE_NUDGE: f64 = 1e-6;

/// The least distance, in degrees, by which [`boundary_tie_on_a_mirror_plane`] moves
/// the input. Measured: at (-10, 11.2 + 1e-12) the engine's answer changes side more
/// than once between 1e-12 and 1e-9 degrees east of the vertex-0 meridian, where this
/// port's changes once, at the meridian; so the engine places a boundary on a mirror
/// plane only to within about 1e-9 degrees. The floor allows ten times that.
pub const TIE_FLOOR_DEG: f64 = 1e-8;

/// How near a pole, in degrees, an input must lie for [`boundary_tie_on_a_mirror_plane`]
/// to treat it as the pole itself: rounding, and nothing more.
pub const POLE_TIE_DEG: f64 = 1e-12;

/// The class rule for boundary ties, of which any listed [`Evidence::BoundaryTie`] entry
/// is a witness. It admits nothing at present: see the last paragraphs below.
///
/// Where the ties lie depends on parity (see the module note). At an even resolution
/// the grid is mirror-symmetric about the three planes of [`mirror_normals`], so a cell
/// boundary lying in one of them is an exact tie along its whole length. At an odd
/// resolution the grid is chiral, and only the half-turns about the six points of
/// [`AXIS_POINTS`] hold, so there are ties at those points and nowhere else. The engine
/// and this port break some of these ties differently, and at a few one side's candidate
/// search accepts neither cell and answers null. Keeping the rule to the planes at even
/// resolutions and to the axis points at odd ones stops it from absorbing, elsewhere, a
/// genuine containment defect small enough to hide inside the nudge.
///
/// It admits a disagreement only when, at this very input: every real answer is a zone
/// the engine reads back; the input lies within [`MIRROR_BAND_DEG`] of a mirror plane
/// and, at an odd resolution, of an axis point too; and the nudges below agree. The step
/// is [`TIE_NUDGE`] of the circumradius of the cell stepped towards, or
/// [`TIE_FLOOR_DEG`], whichever is larger.
///
/// - The input is a pole, to within [`POLE_TIE_DEG`], both answers are real, and by
///   the engine's own vertices each of the two cells has a vertex exactly on that pole,
///   which the pole snap puts there. The pole then belongs to both cells, and no nudge
///   is needed; none is possible either, since at resolution 19 the snap collapses the
///   polar cells' centroids onto the pole and leaves no direction to step in. Where the
///   pole is not a vertex of both, it lies on an edge they share, and the nudges decide.
/// - Both answers real: a step towards our cell's centroid, the engine answers our
///   cell, and a step towards the engine's cell's centroid, we answer the engine's.
/// - One answer null, which happens where one side's candidate search accepts neither
///   of two mirror cells at an exact point: a step towards the real answer's centroid,
///   the side that answered null answers that cell, and a step away from it both sides
///   answer the same real cell, a different one.
///
/// Either way the real answer's boundary passes within one step of the input, by the
/// reckoning of both implementations. Anything else is returned as an error.
///
/// **Where our answer is null, the tie is the one kind in which our answer is the
/// inferior one**: "no cell" for a point that lies on the boundary of a real cell is
/// worse than naming the adjacent cell, which is what the engine does. It is admitted
/// under this rule, but counted apart, so that a new instance shows. **That class, and then the
/// whole of the rule's, is now empty on all three grids**, in both the main sample and
/// the seam sweep, and the ceilings have followed the observation down to nothing; the
/// paragraphs below record what the rule held and why it closed, so that a new instance
/// is read against its history.
///
/// The inferior ties closed in two steps, and the rest of the class in a third. The eC's
/// side test in its written source order,
/// `A * x + B * y + C` (`RI7H.ec:1360-1365`), gives a value a hair below zero on a point
/// lying on the edge two candidate cells share, -5.4e-20 at (4.5, 5.5) at resolution 9,
/// and so rejects both cells. DGGAL is built with `-O2 -ffast-math`
/// (`Makefile.dggal:133`) and `containsPoint` is not one of the functions the eC guards
/// against unsafe maths, so gcc re-associated that sum; the order it emitted was read out
/// of the shipped library and ported, and at (4.5, 5.5) it returns the exact zero that
/// accepts the cell. That took the instances at the midpoint of the edge from base cell
/// 6 to 10, (0, -168.8), at resolutions 9 and 19, which all three grids showed.
///
/// What survived that step was decided not by the arithmetic but by the polygon the
/// arithmetic was given: this port's vertices in the 5x6 CRS differed from
/// `getZoneCRSVertices` by one to a few units in the last place at those stations, so the
/// side test, performed instruction for instruction as the engine performs it, was asked
/// about a slightly different edge. The vertex generation is unguarded fast-math code in
/// DGGAL as well; the order gcc emitted for it was read out of the shipped library too
/// and ported, and every vertex this port computes in the 5x6 CRS is now bit-identical to
/// the engine's over a census of 16,163 zones. That took the remaining instances: the
/// south pole itself at resolution 15 on all three grids, the north pole at 11 on RTEA7H,
/// and (0, 11.2) at resolution 1 on IGEO7.
///
/// The engine did the same to us in the other direction at a few stations, answering null
/// where this port named a cell, and those were counted as ordinary ties; the last of
/// them, at the north pole itself at resolution 1 on IVEA7H, went in the third step, this
/// port now answering the null zone there as the engine does.
///
/// Whether the engine reaches its own answer from a given point depends on the forward
/// projection, and the two are not the same question. Measured against the live engine,
/// with the geographic point set in radians: at IGEO7's and IVEA7H's station on the edge
/// from 6 to 10, at the south-pole station at 15 on all three grids, and at the midpoint
/// of the edge from 3 to 4 at 17 where the engine answers null, the two sides' planar
/// points are identical to the last bit, so the difference was the quantiser's alone. At
/// RTEA7H's station on the edge from 6 to 10 they are one unit in the last place apart in
/// each coordinate, at (-90, 0) at 15 three units and one, and at 11 nine units in the
/// abscissa at (90, 0) and one at the seam's north-pole station; there the engine's own
/// answer belonged to its own point, and the ordinary ties that remained at those stations
/// were differences of the forward projection, not of the quantiser.
///
/// The third step was the projection. Its guarded kernels call unguarded helpers, which
/// gcc compiled under `-ffast-math` as well, and the orders it emitted for them were read
/// out of the shipped library and ported (see `rs4dggs`'s `projections::icovertex`),
/// together with the departures of this port's own from the eC. The forward projection,
/// which gives the planar point every quantisation starts from, was then the engine's own
/// at every station of a census of 3,184 and at 300,000 random points on each grid, and
/// every ordinary tie of both samples went with it: 29 in the main sample and 84 along the
/// seams on IGEO7, 15 and 58 on IVEA7H, 13 and 92 on RTEA7H.
///
/// The centroids are this port's, which the suite asserts to be the engine's to the last
/// bit; they only choose the direction of each nudge, and the decisive answers are the
/// quantisations. Makes its oracle calls one after another.
pub fn boundary_tie_on_a_mirror_plane<S: Subject>(
    lat: f64,
    lon: f64,
    res: u8,
    ours: u64,
    engine: u64,
) -> Result<(), String> {
    let g = S::grid();
    let (ot, et) = (g.text_id(ZoneId(ours)), g.text_id(ZoneId(engine)));
    let why = |what: &str| Err(format!("ours {ot}, engine {et}: {what}"));
    let null = o::NULL_ZONE;
    if ours != null && o::zone_from_text(S::ORACLE, &ot) != ours {
        return why("the engine does not read our answer back");
    }
    if engine != null && o::zone_from_text(S::ORACLE, &et) != engine {
        return why("the engine does not read its own answer back");
    }
    if !on_a_mirror_plane(lat, lon) {
        return why("the input is not on a mirror plane");
    }
    if res % 2 == 1 && !at_an_axis_point(lat, lon) {
        return why(
            "at an odd resolution the grid is chiral, so a tie holds only at one of the six \
             axis points of its half-turns, and the input is at none of them",
        );
    }
    if ours != null && engine != null && 90.0 - lat.abs() <= POLE_TIE_DEG {
        let pole = 90.0f64.copysign(lat);
        let on_pole = |z: u64| o::vertices(S::ORACLE, z).iter().any(|v| v.0 == pole);
        if on_pole(ours) && on_pole(engine) {
            return Ok(());
        }
        // Otherwise the pole lies on an edge the two cells share, as it does at the
        // coarse resolutions, and the nudges below decide.
    }
    let engine_at = |p: (f64, f64)| o::zone_from_geo(S::ORACLE, p.0, p.1, i32::from(res));
    let ours_at = |p: (f64, f64)| g.zone_from_geo(p.0, p.1, res).unwrap().id().0;
    // The centroid of a real answer, and the step to take towards it.
    let aim = |z: u64| {
        let c = g.centroid(ZoneId(z));
        let c = (c.lat, c.lon);
        let r = g
            .vertices(ZoneId(z))
            .iter()
            .map(|v| arc_deg(c, (v.lat, v.lon)))
            .fold(0.0f64, f64::max);
        (c, (TIE_NUDGE * r).max(TIE_FLOOR_DEG))
    };
    let p = (lat, lon);
    if ours != null && engine != null {
        let (co, eps) = aim(ours);
        let there = engine_at(nudge(p, co, eps));
        if there != ours {
            return why(&format!(
                "{eps:e} degrees towards our cell, the engine answers {}",
                o::text_id(S::ORACLE, there)
            ));
        }
        let (ce, eps) = aim(engine);
        let there = ours_at(nudge(p, ce, eps));
        if there != engine {
            return why(&format!(
                "{eps:e} degrees towards the engine's cell, we answer {}",
                g.text_id(ZoneId(there))
            ));
        }
    } else {
        let (real, nulls_side): (u64, &dyn Fn((f64, f64)) -> u64) = if ours != null {
            (ours, &engine_at)
        } else {
            // The inferior tie: we answer no cell where the engine names one. See the
            // note above; the suite counts it apart.
            (engine, &ours_at)
        };
        let (c, eps) = aim(real);
        let there = nulls_side(nudge(p, c, eps));
        if there != real {
            return why(&format!(
                "{eps:e} degrees towards the real answer, the side that answered null answers {}",
                g.text_id(ZoneId(there))
            ));
        }
        let away = nudge(p, c, -eps);
        let (a, b) = (ours_at(away), engine_at(away));
        if a != b || a == real || a == null {
            return why(&format!(
                "{eps:e} degrees away from the real answer, we answer {} and the engine {}",
                g.text_id(ZoneId(a)),
                o::text_id(S::ORACLE, b)
            ));
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------------
// Input this crate refuses. See the module's account of the third kind.
// ---------------------------------------------------------------------------------

/// Non-finite coordinates, each with the parameter this crate names in refusing it: the
/// latitude where it is at fault, and the longitude only where the latitude is finite.
pub const NON_FINITE_COORDINATES: [(f64, f64, &str); 7] = [
    (f64::NAN, 0.0, "latitude"),
    (f64::INFINITY, 0.0, "latitude"),
    (f64::NEG_INFINITY, 0.0, "latitude"),
    (0.0, f64::NAN, "longitude"),
    (0.0, f64::INFINITY, "longitude"),
    (0.0, f64::NEG_INFINITY, "longitude"),
    (f64::NAN, f64::NAN, "latitude"),
];

/// The engine's recorded answer to every non-finite coordinate at resolution `res`, on
/// each of the three aperture-7 grids: the pentagon of base cell 1, the base cell's
/// digit path of zeros, which is one fixed, real cell in the northern Pacific.
pub fn engine_answer_to_a_non_finite_coordinate(res: u8) -> String {
    format!("01{}", "0".repeat(usize::from(res)))
}

/// Texts of twenty digits after the base cell, built for every base cell on two paths:
/// nineteen ones and then each digit, and nineteen zeros, the pentagon's own path, and
/// then each digit. The engine reads each as a level-20 zone, but for the one digit on
/// the pentagon path of each base cell that is the pentagon's deleted child, which it
/// reads as its null zone; that makes [`LEVEL_20_TEXTS_READ`] and
/// [`LEVEL_20_TEXTS_DELETED`] of the 168. This crate refuses all of them, by their length.
pub fn level_20_texts() -> Vec<String> {
    let mut v = Vec::new();
    for base in 0..12 {
        for path in ["1", "0"] {
            for d in 0..7 {
                v.push(format!("{base:02}{}{d}", path.repeat(19)));
            }
        }
    }
    v
}

/// Of [`level_20_texts`], those the engine reads as a level-20 zone.
pub const LEVEL_20_TEXTS_READ: usize = 156;
/// Of [`level_20_texts`], those the engine reads as its null zone: one per base cell,
/// the pentagon's deleted child at the twentieth level.
pub const LEVEL_20_TEXTS_DELETED: usize = 12;

/// Identifiers with no cell behind them: the null zone, and the three others whose base
/// field, 12 to 14, names no base cell, each with every digit slot unused.
pub const FABRICATED: [u64; 4] = [
    0xCFFF_FFFF_FFFF_FFFF,
    0xDFFF_FFFF_FFFF_FFFF,
    0xEFFF_FFFF_FFFF_FFFF,
    o::NULL_ZONE,
];

/// The rule for a level-20 text: the engine reads it, live, as a zone at level 20 with no
/// geometry, the zero point as its centroid, six zero vertices and no neighbours, or, on
/// the pentagon's deleted child, as its null zone; and this crate refuses it. Returns
/// whether the engine read a zone, or an error naming what failed.
///
/// Makes its oracle calls one after another, never one inside another.
pub fn level_20_text_refused<S: Subject>(text: &str) -> Result<bool, String> {
    let ours = S::grid().zone_from_text(text);
    if !matches!(&ours, Err(rs4dggs::Error::InvalidZone(m)) if m.contains("longer than")) {
        return Err(format!(
            "{text}: we answer {ours:?}, not a refusal by length"
        ));
    }
    let e = o::zone_from_text(S::ORACLE, text);
    if e == o::NULL_ZONE {
        // Only the pentagon's deleted child may be read as the null zone: the path of
        // zeros ending in 2 under a northern base cell, 0 to 5, or in 5 under a
        // southern one.
        let base: u32 = text[..2]
            .parse()
            .map_err(|_| format!("{text}: no base cell"))?;
        let deleted = if base <= 5 { '2' } else { '5' };
        let path = &text[2..];
        if path[..19].bytes().all(|c| c == b'0') && path.ends_with(deleted) {
            return Ok(false);
        }
        return Err(format!(
            "{text}: the engine reads its null zone, and the text is not a deleted child"
        ));
    }
    let level = o::level(S::ORACLE, e);
    let centroid = o::centroid(S::ORACLE, e);
    let vertices = o::vertices(S::ORACLE, e);
    let neighbours = o::neighbors(S::ORACLE, e);
    if level != 20
        || centroid != (0.0, 0.0)
        || vertices != vec![(0.0, 0.0); 6]
        || !neighbours.is_empty()
        || o::text_id(S::ORACLE, e) != text
    {
        return Err(format!(
            "{text}: the engine reads level {level}, centroid {centroid:?}, {} vertices and \
             {} neighbours, not a level-20 zone without geometry",
            vertices.len(),
            neighbours.len()
        ));
    }
    Ok(true)
}

// ---------------------------------------------------------------------------------
// The departures of aperture 3. See the module's account of them.
// ---------------------------------------------------------------------------------

/// The icosahedron edges along which the aperture-3 engine answers an identifier it cannot
/// read back, by the polar root it names: root 10 along the edge from base cell 0 to 1,
/// over the north pole, and root 11 along the edge from base cell 6 to 11, which ends at
/// the southern polar vertex, as pairs of base cells. The other twenty-eight edges show
/// neither in any sample of the suite.
pub const UNREADABLE_POLAR_SEAMS: [(u64, (usize, usize)); 2] = [(10, (0, 1)), (11, (6, 11))];

/// How far from its edge an input the rule for unreadable polar roots admits may lie, in
/// degrees. Measured on each aperture-3 grid over the samples of its suite, the widest in
/// `unreadable_answers_along_two_seams_are_the_null_zone`: the inputs it admits lie within
/// 2.24e-5 degrees of the edge from base cell 0 to 1, the farthest a hair from the north
/// pole, and within 2.05e-5 of the edge from 6 to 11 on ISEA3H; within 2.18e-5 and 2.00e-5
/// on IVEA3H; and within 2.19e-5 and 1.74e-5 on RTEA3H. The bound is 4e-5, 1.8 times the
/// widest reach, ISEA3H's, and 1.83 times the widest of each of the other two.
pub const UNREADABLE_POLAR_BAND_DEG: f64 = 4e-5;

/// The resolutions at which the aperture-3 engine answers an identifier it cannot read back
/// along the two seams; at 0 and 1 and from 22 on it answers none.
pub const UNREADABLE_POLAR_RESOLUTIONS: std::ops::RangeInclusive<u8> = 2..=21;

/// The class rule of which the aperture-3 entries of [`DIVERGENCES`] are witnesses. It
/// admits a disagreement only when all of the evidence holds at this very input, checked
/// against the live engine now: our answer is the null zone; the engine's has the root of a
/// polar pentagon, 10 or 11, and a non-zero index, and the level asked for; the engine
/// cannot read it back (`engine_can_read` is false) and neither can this crate; the
/// resolution is one of [`UNREADABLE_POLAR_RESOLUTIONS`]; the input lies within
/// [`UNREADABLE_POLAR_BAND_DEG`] of the edge of [`UNREADABLE_POLAR_SEAMS`] that goes with
/// the root; and the engine draws its answer as that polar pentagon, centred on the very
/// double it gives the pentagon's own centroid, with five vertices. Anything else is
/// returned as an error for the caller to fail on.
///
/// The engine is asked nothing relational about its answer: only its text, its reading of
/// that text, its level, its centroid and its vertices, none of which harms the process.
/// Makes its oracle calls one after another, never one inside another.
pub fn engine_answers_an_unreadable_polar_root<S: Subject>(
    lat: f64,
    lon: f64,
    res: u8,
    ours: u64,
    engine: u64,
) -> Result<(), String> {
    let text = o::text_id(S::ORACLE, engine);
    let why = |what: String| {
        Err(format!(
            "engine {text} ({engine:#018x}), ours {ours:#018x}: {what}"
        ))
    };
    if engine == o::NULL_ZONE {
        return why("the engine's answer is the null zone".into());
    }
    let (root, index) = ((engine >> 53) & 0xF, (engine >> 2) & ((1 << 51) - 1));
    let Some(&(_, (i, j))) = UNREADABLE_POLAR_SEAMS.iter().find(|&&(r, _)| r == root) else {
        return why(format!(
            "the engine's answer has root {root}, not a polar one"
        ));
    };
    if index == 0 {
        return why("the engine's answer is a polar pentagon, at index 0".into());
    }
    let mut missing = Vec::new();
    if ours != o::NULL_ZONE {
        missing.push("our answer is not the null zone".to_string());
    }
    if o::level(S::ORACLE, engine) != i32::from(res) {
        missing.push("the engine's answer is not at the level asked for".into());
    }
    if o::engine_can_read(S::ORACLE, engine) {
        missing.push("the engine reads its own answer back".into());
    }
    if S::grid().zone_from_text(&text).is_ok() {
        missing.push("this crate reads the engine's answer".into());
    }
    if !UNREADABLE_POLAR_RESOLUTIONS.contains(&res) {
        missing.push("the resolution is not one of 2 to 21".into());
    }
    let v = icosahedron_vertices();
    let off = dist_to_arc_deg((lat, lon), v[i], v[j]);
    if off > UNREADABLE_POLAR_BAND_DEG {
        missing.push(format!(
            "the input lies {off:e} degrees from the edge from base cell {i} to {j}"
        ));
    }
    let pentagon = o::zone_from_text(S::ORACLE, &format!("A{root:X}-0-A"));
    let (c, pc) = (
        o::centroid(S::ORACLE, engine),
        o::centroid(S::ORACLE, pentagon),
    );
    if (c.0.to_bits(), c.1.to_bits()) != (pc.0.to_bits(), pc.1.to_bits()) {
        missing.push(format!(
            "the engine centres its answer at {c:?}, not on the polar pentagon's {pc:?}"
        ));
    }
    let n = o::vertices(S::ORACLE, engine).len();
    if n != 5 {
        missing.push(format!("the engine draws its answer with {n} vertices"));
    }
    if missing.is_empty() {
        Ok(())
    } else {
        why(missing.join("; "))
    }
}

/// The centroid the aperture-3 engine gives its null zone on the grid under test, read from
/// the engine through this suite's oracle: the planar point that the arithmetic of the
/// `centroid` property makes of the null zone's bits, a root of 15 and an index of
/// `2^51 - 1`, carried to the sphere by the grid's own projection, so that it differs by
/// grid.
pub fn aperture_3_engine_null_centroid<S: Subject>() -> (f64, f64) {
    match S::NAME {
        "ISEA3H" => (-21.209_249_684_002_74, 35.706_408_821_777_72),
        "IVEA3H" => (-21.274_705_264_329_572, 35.681_841_845_149_46),
        "RTEA3H" => (-21.513_812_911_472_883, 35.591_901_296_748_87),
        other => panic!("no null-zone centroid of the engine is recorded for {other}"),
    }
}

/// The engine's side of the null zone's geometry on aperture 3, re-checked live: its
/// identifier is the one `(null)` names, the engine cannot read it at a level of 33 or
/// less (it reads it at 63), and it gives it the centroid
/// [`aperture_3_engine_null_centroid`], to the last bit, and six vertices within a
/// millionth of a degree of it. This crate answers the point (0, 0) and six zero vertices.
///
/// Makes its oracle calls one after another, and none that is relational.
pub fn engine_null_zone_geometry<S: Subject>() -> Result<(), String> {
    let null = o::NULL_ZONE;
    let mut missing = Vec::new();
    if o::zone_from_text(S::ORACLE, rs4dggs::NULL_TEXT) != null {
        missing.push("the engine reads (null) as another identifier".to_string());
    }
    if o::engine_can_read(S::ORACLE, null) || o::level(S::ORACLE, null) != 63 {
        missing.push(format!(
            "the engine reads its null zone at level {}",
            o::level(S::ORACLE, null)
        ));
    }
    let c = o::centroid(S::ORACLE, null);
    let want = aperture_3_engine_null_centroid::<S>();
    if (c.0.to_bits(), c.1.to_bits()) != (want.0.to_bits(), want.1.to_bits()) {
        missing.push(format!("the engine centres its null zone at {c:?}"));
    }
    let vs = o::vertices(S::ORACLE, null);
    if vs.len() != 6 || vs.iter().any(|&p| arc_deg(p, want) > 1e-6) {
        missing.push(format!("the engine draws its null zone as {vs:?}"));
    }
    if missing.is_empty() {
        Ok(())
    } else {
        Err(format!("the engine's null zone: {}", missing.join("; ")))
    }
}

/// The largest sub-zone order that the engine is asked for whole, by
/// [`engine_sub_zone_order_is_faulty`] on aperture 3 and by the comparison of the orders on
/// aperture 7 (`common::sub_zone_order`): well below the `2^28` sub-zones at which its
/// `getSubZones` answers nothing on aperture 3 and ends the process on aperture 7, and quick
/// to build.
pub const LONGEST_ORDER_ASKED: u64 = 1 << 20;

/// Whether the aperture-7 engine's `getFirstSubZone` of `id` at `depth` ends the process: at
/// the pentagon of the south pole (base cell 11, every digit 0) of an odd level, at depth 2,
/// where it subtracts through a null pointer (`RI7H.ec:3246-3255`). No Rust mechanism can
/// catch that, so the call is not made there, and the comparison of the orders counts the
/// pairs it leaves out. The engine's `getSubZones` answers at the same pairs, and its first
/// entry is what this crate's `first_sub_zone` is held to.
pub fn engine_first_sub_zone_ends_the_process<S: Subject>(id: ZoneId, depth: u8) -> bool {
    let text = S::grid().text_id(id);
    let (base, digits) = text.split_at(2);
    depth == 2 && base == "11" && digits.len() % 2 == 1 && digits.bytes().all(|d| d == b'0')
}

/// What is wrong with the aperture-3 engine's sub-zone order of a zone, by its own tests: the
/// places that repeat a zone already named; the distinct zones that its `zoneHasSubZone`
/// rejects; and its `countSubZones` less the distinct zones that its `zoneHasSubZone` accepts,
/// the descendants it omits.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FaultyOrder {
    pub repeats: usize,
    pub rejected: usize,
    pub omitted: u64,
}

/// The rule for a sub-zone order of this crate that differs from the aperture-3 engine's. It
/// holds only when all of the evidence holds, checked against the live engine now:
/// - the engine's order of `parent` at `depth`, asked for only where it holds at most
///   [`LONGEST_ORDER_ASKED`] sub-zones, is not the zone's descendants: it names a zone twice,
///   names one that its own `zoneHasSubZone` rejects, or holds fewer distinct zones than its
///   own `countSubZones`;
/// - this crate's order, `ours`, is: it names each zone once, each accepted by the engine's
///   `zoneHasSubZone`, as many as the engine's count, and at depth 1 the engine's
///   `getZoneChildren` as a set.
///
/// The engine is then at fault, by its own tests, and this crate's order is the consistent one.
/// What is wrong with the engine's order is returned, so that the caller can report it.
///
/// Asks the engine about no zone it cannot read back, and about no sub-zone at level 33, where
/// its `zoneHasSubZone` fails on the needle's centroid child, of level 34: such an order is
/// refused as one that cannot be checked. Makes its oracle calls one after another.
pub fn engine_sub_zone_order_is_faulty<S: Subject>(
    parent: u64,
    depth: u8,
    ours: &[u64],
) -> Result<FaultyOrder, String> {
    use std::collections::{BTreeMap, BTreeSet};
    let at = format!("{} at depth {depth}", o::text_id(S::ORACLE, parent));
    if !o::engine_can_read(S::ORACLE, parent) {
        return Err(format!("{at}: the engine cannot read the zone back"));
    }
    if o::level(S::ORACLE, parent) + i32::from(depth) >= 33 {
        return Err(format!(
            "{at}: the sub-zones lie at level 33, where the engine's zoneHasSubZone fails"
        ));
    }
    let count = o::count_sub_zones(S::ORACLE, parent, i32::from(depth));
    if count > LONGEST_ORDER_ASKED {
        return Err(format!(
            "{at}: the engine's order holds {count} sub-zones, too many to be asked for"
        ));
    }
    // The engine's `zoneHasSubZone` of each distinct zone, asked once.
    let mut accepted: BTreeMap<u64, bool> = BTreeMap::new();
    let mut has = |s: u64| -> Result<bool, String> {
        if let Some(&a) = accepted.get(&s) {
            return Ok(a);
        }
        if !o::engine_can_read(S::ORACLE, s) {
            return Err(format!(
                "{at}: the engine cannot read {} back",
                o::text_id(S::ORACLE, s)
            ));
        }
        let a = o::zone_has_sub_zone(S::ORACLE, parent, s);
        accepted.insert(s, a);
        Ok(a)
    };

    let theirs = o::sub_zones(S::ORACLE, parent, i32::from(depth));
    let distinct: BTreeSet<u64> = theirs.iter().copied().collect();
    let (mut rejected, mut kept) = (0, 0u64);
    for &s in &distinct {
        if has(s)? {
            kept += 1;
        } else {
            rejected += 1;
        }
    }
    let fault = FaultyOrder {
        repeats: theirs.len() - distinct.len(),
        rejected,
        omitted: count.saturating_sub(kept),
    };
    if fault.repeats == 0 && fault.rejected == 0 && fault.omitted == 0 {
        return Err(format!(
            "{at}: the engine's order is the zone's descendants, each once, and this crate's \
             departs from it"
        ));
    }

    let mut missing = Vec::new();
    let ours_distinct: BTreeSet<u64> = ours.iter().copied().collect();
    if ours_distinct.len() != ours.len() {
        missing.push(format!(
            "this crate's order names {} places twice",
            ours.len() - ours_distinct.len()
        ));
    }
    if ours.len() as u64 != count {
        missing.push(format!(
            "this crate's order holds {} sub-zones, the engine's count {count}",
            ours.len()
        ));
    }
    for &s in &ours_distinct {
        if !has(s)? {
            missing.push(format!(
                "the engine's zoneHasSubZone rejects {} of this crate's order",
                o::text_id(S::ORACLE, s)
            ));
        }
    }
    if depth == 1 {
        let children: BTreeSet<u64> = o::children(S::ORACLE, parent).into_iter().collect();
        if children != ours_distinct {
            missing.push("this crate's order is not the engine's children".into());
        }
    }
    if missing.is_empty() {
        Ok(fault)
    } else {
        Err(format!("{at}: {}", missing.join("; ")))
    }
}

/// The aperture-3 engine's recorded answer to every non-finite coordinate at resolution
/// `res`, on each of the three aperture-3 grids: the pentagon of root rhombus 0 at that
/// level, `A0-0-A` at resolution 0, `A0-0-B` at 1, `C0-0-B` at 5 and `Q0-0-B` at 33, one
/// fixed, real cell centred on the icosahedron's vertex of base cell 1.
pub fn aperture_3_engine_answer_to_a_non_finite_coordinate(res: u8) -> String {
    let sub_hex = if res % 2 == 0 { 'A' } else { 'B' };
    format!("{}0-0-{sub_hex}", char::from(b'A' + res / 2))
}

/// How the engine answers an identifier for which this crate has no geometry, whose ring is
/// the empty one, whose extent and whose area are `None`. Each is a departure of
/// this crate from the engine in at least one of the three answers, and the engine's side of
/// each is re-checked live by [`engine_answers_without_geometry`] at every instance.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum NoGeometry {
    /// The null zone: the engine gives it no ring and the cleared extent, as this crate does,
    /// and the area `+inf`, where this crate gives none.
    NullZone,
    /// An identifier the engine does not draw: no ring (empty on aperture 7, null and so not
    /// asked on aperture 3), the cleared extent, and the area of the formula of its level, a
    /// finite number. On aperture 7, a Z7 identifier of twenty digits; on aperture 3, an
    /// identifier it cannot read back from its own text. The ring and the extent agree with
    /// this crate's; the area does not.
    Undrawn,
    /// As [`NoGeometry::Undrawn`], with the area `+inf` where the engine's count of the zones
    /// of the level leaves a divisor of zero: an identifier with no cell behind it.
    UndrawnWithoutArea,
    /// A zone the engine reads and draws, ring, extent and area, and for which this crate has
    /// no cell: the sub-hexes C and D of a polar root of aperture 3. All three answers depart.
    Drawn,
}

/// What the engine answers for the identifier `id` of the grid under test, asked of the live
/// engine: whether its ring was asked (on aperture 3, only where the engine can read the
/// identifier), the points of that ring, whether its extent is the cleared one, and its area.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EngineAnswers {
    pub ring_asked: bool,
    pub ring_points: usize,
    pub extent_cleared: bool,
    pub extent_is_finite: bool,
    pub area: f64,
}

/// The engine's three answers for `id`, one oracle call for each, none inside another.
pub fn engine_answers<S: Subject>(id: u64) -> EngineAnswers {
    use rs4dggs::Topology;
    let x = f64::from_bits(0x7f91_df46_a252_9d38);
    let ring_asked = S::T::APERTURE != 3 || o::engine_can_read(S::ORACLE, id);
    let ring_points = if ring_asked {
        o::refined_vertices(S::ORACLE, id, 0).len()
    } else {
        0
    };
    let extent = o::extent(S::ORACLE, id);
    let cleared = [x, x, -x, -x].map(f64::to_bits);
    EngineAnswers {
        ring_asked,
        ring_points,
        extent_cleared: extent.map(f64::to_bits) == cleared,
        extent_is_finite: extent.iter().all(|m| m.is_finite() && m.abs() < 4.0),
        area: o::area(S::ORACLE, id),
    }
}

/// The sub-hexes C and D of the two polar roots of aperture 3 at the odd levels 1 to 33, as the
/// engine reads them from their text: this crate has no cell for any. 68 of them. The only
/// identifiers that [`NoGeometry::Drawn`] admits.
pub fn polar_sub_hexes<S: Subject>() -> Vec<u64> {
    use rs4dggs::Topology;
    use std::collections::HashMap;
    use std::sync::{Mutex, OnceLock};
    // There are none on aperture 7, where the text `AA-0-C` is not an identifier. On aperture 3
    // the list is built once per grid, the first time it is asked for.
    static BUILT: OnceLock<Mutex<HashMap<&'static str, Vec<u64>>>> = OnceLock::new();
    if S::T::APERTURE != 3 {
        return Vec::new();
    }
    let mut built = BUILT.get_or_init(Default::default).lock().unwrap();
    built
        .entry(S::ORACLE)
        .or_insert_with(|| {
            let mut ids = Vec::new();
            for level in 0..=16u8 {
                for root in ['A', 'B'] {
                    for sub_hex in ['C', 'D'] {
                        let text = format!("{}{root}-0-{sub_hex}", char::from(b'A' + level));
                        let id = o::zone_from_text(S::ORACLE, &text);
                        assert_ne!(id, o::NULL_ZONE, "the engine no longer reads {text}");
                        ids.push(id);
                    }
                }
            }
            ids
        })
        .clone()
}

/// Whether `id` belongs to the class of identifiers the engine does not draw, which
/// [`NoGeometry::Undrawn`] and [`NoGeometry::UndrawnWithoutArea`] admit: the fabricated
/// identifiers; on aperture 7 a Z7 identifier of level 20; on aperture 3 an identifier of a polar
/// root (10 or 11) with a non-zero index that the engine cannot read, the class its quantiser
/// hands out along two seams.
fn in_the_undrawn_class<S: Subject>(id: u64) -> bool {
    use rs4dggs::Topology;
    if id != o::NULL_ZONE && FABRICATED.contains(&id) {
        return true;
    }
    match S::T::APERTURE {
        7 => o::level(S::ORACLE, id) == 20,
        _ => (id >> 53) & 0xF >= 10 && !o::engine_can_read(S::ORACLE, id),
    }
}

/// The kind of the departure at `id`, an identifier this crate has no geometry for, with the
/// engine's side of it re-checked live; or what failed. The evidence of each kind:
/// - the null zone: the cleared extent, the area `+inf`, no ring at any refinement from 0 to 3
///   where it is asked;
/// - an identifier the engine does not draw: the cleared extent, no ring where it is asked, a
///   finite and positive area, or `+inf` for one with no cell behind it; and the identifier is
///   of the class ([`in_the_undrawn_class`]), so that an ordinary zone left empty fails here;
/// - one the engine draws: a ring of at least five points, a finite extent, a finite and
///   positive area, and the identifier is one of [`polar_sub_hexes`], so that on aperture 7
///   nothing is `Drawn` and an ordinary zone left empty fails here, naming the zone.
pub fn engine_answers_without_geometry<S: Subject>(id: u64) -> Result<NoGeometry, String> {
    let a = engine_answers::<S>(id);
    let at = format!("{id:#018x} on {}: {a:?}", S::NAME);
    if !a.ring_asked || a.ring_points == 0 {
        // No ring, or one not asked: the engine does not draw it.
        if a.ring_asked {
            for r in 1..=3 {
                if !o::refined_vertices(S::ORACLE, id, r).is_empty() {
                    return Err(format!("{at}: no ring at 0 and one at a refinement of {r}"));
                }
            }
        }
        if !a.extent_cleared {
            return Err(format!(
                "{at}: no ring and an extent that is not the cleared one"
            ));
        }
        if id == o::NULL_ZONE {
            return if a.area == f64::INFINITY {
                Ok(NoGeometry::NullZone)
            } else {
                Err(format!("{at}: the null zone's area is not +inf"))
            };
        }
        if !in_the_undrawn_class::<S>(id) {
            return Err(format!(
                "{at}: not of the class of identifiers the engine does not draw"
            ));
        }
        return if a.area == f64::INFINITY {
            Ok(NoGeometry::UndrawnWithoutArea)
        } else if a.area.is_finite() && a.area > 0.0 {
            Ok(NoGeometry::Undrawn)
        } else {
            Err(format!("{at}: an area that is neither finite nor +inf"))
        };
    }
    if a.ring_points >= 5
        && a.extent_is_finite
        && !a.extent_cleared
        && a.area.is_finite()
        && a.area > 0.0
        && id != o::NULL_ZONE
        && polar_sub_hexes::<S>().contains(&id)
    {
        Ok(NoGeometry::Drawn)
    } else {
        Err(format!("{at}: not an answer of any recorded kind"))
    }
}

/// The rule for an answer of this crate to a box that differs from the engine's `listZones`.
/// `ours` and `theirs` are the two answers for the box `bbox`, in degrees, at `level`. It holds
/// only when all of the evidence holds, checked against the live engine now:
/// - the two answers differ, and the engine's is this crate's less some zones: every zone of
///   the engine's is found in this crate's, in the same order;
/// - each zone that the engine leaves out is one it reads back, of the level asked, and its
///   extent, the engine's own, meets the box by the engine's own test of two extents
///   ([`boxes::meets`]).
///
/// The engine then leaves out zones that its own rule puts in the box, and this crate's answer
/// is the consistent one. How many it leaves out is returned. A zone of the engine's that this
/// crate lacks is never admitted.
pub fn engine_omits_zones_of_a_box<S: Subject>(
    level: u8,
    bbox: &[f64; 4],
    ours: &[u64],
    theirs: &[u64],
) -> Result<usize, String> {
    if ours == theirs {
        return Err("the two answers are one, and nothing departs".into());
    }
    let mut rest = ours.iter();
    if let Some(&lacking) = theirs.iter().find(|&&t| !rest.any(|&z| z == t)) {
        return Err(format!(
            "{} is in the engine's answer of {} zones and not, at its place, in this crate's of {}",
            o::text_id(S::ORACLE, lacking),
            theirs.len(),
            ours.len()
        ));
    }
    let theirs: std::collections::HashSet<u64> = theirs.iter().copied().collect();
    let r = boxes::radians(bbox);
    let mut left_out = 0;
    for &zone in ours.iter().filter(|z| !theirs.contains(z)) {
        let text = o::text_id(S::ORACLE, zone);
        if !o::engine_can_read(S::ORACLE, zone) {
            return Err(format!("the engine cannot read {text} back"));
        }
        if o::level(S::ORACLE, zone) != i32::from(level) {
            return Err(format!("{text} is not of level {level}"));
        }
        let extent = o::extent(S::ORACLE, zone);
        if !boxes::meets(&extent, &r) {
            return Err(format!(
                "{text} is not in the engine's answer, and its extent {extent:?}, the engine's \
                 own, does not meet the box {r:?}"
            ));
        }
        left_out += 1;
    }
    Ok(left_out)
}
