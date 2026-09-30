//! Disagreements with the engine that are measured and bounded but not resolved.
//!
//! This is not a divergence: this port's answer has not been shown to be the right
//! one, so it must not be read as a ruling. The regions are findings about one grid
//! each, and [`regions`] gives the set for the grid under test; no grid is assumed to
//! share another's. IGEO7 has one region, [`BROKEN_SEAMS`], where the engine is not
//! consistent with itself, and IVEA7H and RTEA7H have the same one, measured afresh on
//! each grid: [`IVEA7H_BROKEN_SEAMS`] and [`RTEA7H_BROKEN_SEAMS`]. Inside it one
//! quantisation may differ, the case in which neither answer contains the input, and
//! the rule for the answers the engine cannot read back applies only there. The tests
//! count every such case, print the counts and hold them under ceilings, so that the
//! region cannot grow unnoticed. Outside it the answers must match exactly.
//!
//! No region admits a neighbour list or a disk any longer. The aperture-7 grids list a
//! zone's neighbours by the engine's own algorithm, and the suites compare every list with
//! the engine's, entry for entry, less four kinds of entry this crate drops, each
//! re-checked against the engine; the disks are composed from the same lists. The region
//! admitted a list that differed only on the evidence of the engine's own list for the
//! zone, while the lists came from a geometric search of this crate's own.
//!
//! Two earlier regions recorded defects of this port and were removed once the defects
//! were fixed. One was the pole snap, the eC's `fixPoles`, now ported; the polar
//! coordinates match the engine exactly, and `the_poles_snap_as_the_engine_does`
//! asserts it. The other was the degenerate poles, where that geometric search
//! under-reported near the poles.

use dggal_oracle as o;

use super::expected_divergences::{BROKEN_SEAM_BAND_DEG, BROKEN_SEAMS as BROKEN_EDGES};
use super::subject::Subject;
use super::{arc_deg, dist_to_arc_deg, icosahedron_vertices};

/// One unresolved region: a strip along some of the icosahedron's edges, from some
/// resolution on.
#[derive(Debug)]
pub struct Region {
    pub name: &'static str,
    pub account: &'static str,
    /// The icosahedron edges the strip runs along, as pairs of base cells.
    pub edges: &'static [(usize, usize)],
    /// How far from those edges the strip reaches, in degrees.
    pub band_deg: f64,
    /// The coarsest resolution at which the region applies.
    pub min_res: u8,
}

impl Region {
    /// Whether a point lies in the strip, at whatever resolution.
    pub fn near(&self, lat: f64, lon: f64) -> bool {
        let v = icosahedron_vertices();
        self.edges
            .iter()
            .any(|&(i, j)| dist_to_arc_deg((lat, lon), v[i], v[j]) <= self.band_deg)
    }

    /// Whether a point lies in the region at this resolution.
    pub fn contains(&self, lat: f64, lon: f64, res: u8) -> bool {
        res >= self.min_res && self.near(lat, lon)
    }
}

/// The unresolved regions of the grid under test. The broken seams were measured on
/// IGEO7 and are not assumed of any other grid: a grid has a region here only once the
/// evidence for it has been measured against the engine, and until then every
/// disagreement on that grid must be a characterised divergence or fail.
pub fn regions<S: Subject>() -> &'static [Region] {
    match S::NAME {
        "IGEO7" => &IGEO7_REGIONS,
        "IVEA7H" => &IVEA7H_REGIONS,
        "RTEA7H" => &RTEA7H_REGIONS,
        "ISEA3H" => &ISEA3H_REGIONS,
        "IVEA3H" => &IVEA3H_REGIONS,
        "RTEA3H" => &RTEA3H_REGIONS,
        other => panic!("no set of unresolved regions is recorded for {other}"),
    }
}

/// The regions of ISEA3H: none. Every disagreement with its engine that the samples meet is
/// either a characterised divergence or a failure; its neighbour lists, its hierarchy and its
/// sub-zones agree with the engine's wherever the engine can be asked.
static ISEA3H_REGIONS: [Region; 0] = [];

/// The regions of IVEA3H: none, as on ISEA3H, and measured so on this grid.
static IVEA3H_REGIONS: [Region; 0] = [];

/// The regions of RTEA3H: none, as on ISEA3H, and measured so on this grid.
static RTEA3H_REGIONS: [Region; 0] = [];

/// The regions of IGEO7.
static IGEO7_REGIONS: [Region; 1] = [BROKEN_SEAMS];

/// The regions of IVEA7H.
static IVEA7H_REGIONS: [Region; 1] = [IVEA7H_BROKEN_SEAMS];

/// The regions of RTEA7H.
static RTEA7H_REGIONS: [Region; 1] = [RTEA7H_BROKEN_SEAMS];

/// Along the two icosahedron edges where the engine's odd resolutions have a strip of
/// no cells (see `expected_divergences`), from resolution 15.
pub const BROKEN_SEAMS: Region = Region {
    name: "broken seams",
    account: "Within 1.5e-4 degrees of the edge from base cell 0 to 1 (over the north pole) \
              or of the edge from 1 to 6, at resolutions 15 to 19. The engine is not \
              consistent with itself there: its quantisation answers its null zone at the \
              odd resolutions and identifiers it cannot read back at 16 and 18, and its \
              neighbour lists contain its null zone, identifiers it cannot read back, the \
              zone itself, repeats and cells some 36 degrees away, and they are not \
              reciprocal. Two distinct identifiers can name the same polygon. Where the two \
              sides name two different cells of which neither contains the input, neither \
              side is an oracle for the other, and the case is counted. The neighbour lists \
              are not admitted here: they are the engine's own, less four kinds of entry, \
              everywhere. Resolutions 13 and 14 match exactly along the same edges, and the \
              other twenty-eight edges match at every resolution.",
    edges: &BROKEN_EDGES,
    band_deg: BROKEN_SEAM_BAND_DEG,
    min_res: BROKEN_SEAMS_MIN_RES,
};

/// The broken seams of IVEA7H, measured on that grid and found along IGEO7's: the same
/// two edges and the same resolutions, and a strip no wider, so the region takes the same
/// bounds (see `BROKEN_SEAM_BAND_DEG`). The engine's IVEA7H class is the
/// same aperture-7 topology as its IGEO7 class with another projection in it
/// (`IVEA7H_Z7.ec:12-16` against `ISEA7H_Z7.ec:12-16`), and the two place the base cells
/// on the same twelve points.
pub const IVEA7H_BROKEN_SEAMS: Region = Region {
    name: "broken seams",
    account: "Within 1.5e-4 degrees of the edge from base cell 0 to 1 (over the north pole) \
              or of the edge from 1 to 6, at resolutions 15 to 19, measured on IVEA7H. A \
              sweep across all thirty edges at resolutions 13 to 19, out to 1e-3 degrees \
              either side, finds the strip along these two edges only, from resolution 15 \
              and within 7.5e-5 degrees of the edge: there the engine answers its null zone \
              at the odd resolutions and identifiers it cannot read back at 18, and at 16 \
              along the edge from 0 to 1, and at 19 some of those too, at the strip's outer \
              margin. Elsewhere the sweep finds nothing but boundary ties on the mirror \
              planes. Off the samples, between 1.8e-5 and 2.5e-5 degrees from the north \
              pole at resolution 16, the two sides can also name two different cells of \
              which neither contains the input, as on IGEO7 and over much the same range; \
              neither side is an oracle for the other there. The neighbour lists are not \
              admitted here: they are the engine's own, less four kinds of entry, \
              everywhere.",
    edges: &BROKEN_EDGES,
    band_deg: BROKEN_SEAM_BAND_DEG,
    min_res: BROKEN_SEAMS_MIN_RES,
};

/// The broken seams of RTEA7H, measured on that grid and found along IGEO7's: the same
/// two edges, the same resolutions, and a strip no wider, so the region takes the same
/// bounds (see `BROKEN_SEAM_BAND_DEG`). The engine's RTEA7H class is the same aperture-7
/// topology as its IGEO7 class with another projection in it (`RTEA7H_Z7.ec:12-16`
/// against `ISEA7H_Z7.ec:12-16`), and the two place the base cells on the same twelve
/// points.
pub const RTEA7H_BROKEN_SEAMS: Region = Region {
    name: "broken seams",
    account: "Within 1.5e-4 degrees of the edge from base cell 0 to 1 (over the north pole) \
              or of the edge from 1 to 6, at resolutions 15 to 19, measured on RTEA7H. A \
              sweep across all thirty edges at resolutions 13 to 19, out to 1e-3 degrees \
              either side, finds the strip along these two edges only, from resolution 15 \
              and within 7.5e-5 degrees of the edge: there the engine answers its null zone \
              at the odd resolutions and identifiers it cannot read back at 18, and at 16 \
              along the edge from 0 to 1, and at 19 some of those too, at the strip's outer \
              margin. Elsewhere the sweep finds nothing but boundary ties on the mirror \
              planes. Where the two sides of IGEO7 and IVEA7H name two different cells of \
              which neither contains the input, between 1.8e-5 and 2.5e-5 degrees from the \
              north pole at resolution 16, the two sides of RTEA7H name the same cell, which \
              by the engine's own vertices does not contain the input either, the input \
              lying 2.6 of that cell's circumradii outside it at IGEO7's own point; a polar \
              scan at resolutions 15 to 19 finds no two differing real answers there of \
              which neither contains the input. The neighbour lists are not admitted here: \
              they are the engine's own, less four kinds of entry, everywhere.",
    edges: &BROKEN_EDGES,
    band_deg: BROKEN_SEAM_BAND_DEG,
    min_res: BROKEN_SEAMS_MIN_RES,
};

/// How far inside a polygon a point lies, in degrees: positive inside, negative outside,
/// the distance to the nearest edge either way. The polygon is taken as convex and
/// projected gnomonically onto the plane tangent at the point, where its edges are
/// straight, so the measure holds for cells of any size at any latitude, the poles
/// included. Checked while this suite was written: every one of 3,027 engine answers
/// over a seeded sample, polar points included, contains its input by this measure, bar
/// the single case [`neither_answer_contains_the_input`] records.
pub fn depth_inside_deg(p: (f64, f64), polygon: &[(f64, f64)]) -> f64 {
    let u = super::unit_vector(p.0, p.1);
    let helper = if u[2].abs() < 0.9 {
        [0.0, 0.0, 1.0]
    } else {
        [1.0, 0.0, 0.0]
    };
    let cross = |a: [f64; 3], b: [f64; 3]| {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    };
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let e1 = cross(helper, u);
    let n1 = dot(e1, e1).sqrt();
    let e1 = [e1[0] / n1, e1[1] / n1, e1[2] / n1];
    let e2 = cross(u, e1);
    let flat: Vec<(f64, f64)> = polygon
        .iter()
        .map(|v| {
            let w = super::unit_vector(v.0, v.1);
            let s = dot(w, u);
            (dot(w, e1) / s, dot(w, e2) / s)
        })
        .collect();
    let n = flat.len();
    let area: f64 = (0..n)
        .map(|i| {
            let (a, b) = (flat[i], flat[(i + 1) % n]);
            a.0 * b.1 - a.1 * b.0
        })
        .sum();
    let sign = area.signum();
    (0..n)
        .filter_map(|i| {
            let (a, b) = (flat[i], flat[(i + 1) % n]);
            let len = (b.0 - a.0).hypot(b.1 - a.1);
            // A degenerate edge, two vertices at the same place, bounds nothing.
            (len > 0.0)
                .then(|| sign * ((b.0 - a.0) * (0.0 - a.1) - (b.1 - a.1) * (0.0 - a.0)) / len)
        })
        .fold(f64::INFINITY, f64::min)
        .to_degrees()
}

/// The margin, as a fraction of a cell's circumradius, by which the input must lie
/// outside each answer for [`neither_answer_contains_the_input`] to hold.
pub const EXCLUSION_MARGIN: f64 = 1e-3;

/// Whether one of the grid's regions accounts for a quantisation disagreement: inside
/// the region, for IGEO7 the strip of [`BROKEN_SEAMS`] at resolution 15 or finer, both
/// answers real, and the input outside both of them by the engine's own vertices, each
/// by more than [`EXCLUSION_MARGIN`] of the cell's circumradius. Neither side is right
/// there, so this is counted, not admitted as a divergence. Measured again on 2026-09-23,
/// after the projection work that moved polar planar points: one case in the main sample
/// of IGEO7, 2.02e-5 degrees from the north pole at resolution 16, where the input lies
/// 1.90 of our cell's circumradii outside our answer and 2.74 of the engine's outside its
/// own.
///
/// Makes its oracle calls one after another.
pub fn neither_answer_contains_the_input<S: Subject>(
    lat: f64,
    lon: f64,
    res: u8,
    ours: u64,
    engine: u64,
) -> Result<(), String> {
    if ours == o::NULL_ZONE || engine == o::NULL_ZONE {
        return Err("one answer is the null zone".into());
    }
    if region_containing::<S>((lat, lon), res).is_none() {
        return Err("the input is not in one of the grid's regions at this resolution".into());
    }
    for (who, z) in [("ours", ours), ("the engine's", engine)] {
        let poly = o::vertices(S::ORACLE, z);
        let c = o::centroid(S::ORACLE, z);
        let radius = poly.iter().map(|v| arc_deg(c, *v)).fold(0.0f64, f64::max);
        let depth = depth_inside_deg((lat, lon), &poly);
        if depth > -EXCLUSION_MARGIN * radius {
            return Err(format!("the input lies {depth:e} degrees inside {who}"));
        }
    }
    Ok(())
}

/// The resolution from which [`BROKEN_SEAMS`] applies. Measured: 14 and below match.
pub const BROKEN_SEAMS_MIN_RES: u8 = 15;

/// The region of the grid under test that a point lies in at this resolution, if any.
pub fn region_containing<S: Subject>(p: (f64, f64), res: u8) -> Option<&'static Region> {
    regions::<S>().iter().find(|r| r.contains(p.0, p.1, res))
}
