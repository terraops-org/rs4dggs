//! Shared helpers for the oracle tests: seeded sampling, the adversarial inputs, the
//! icosahedron's edges, the exact comparison of coordinates with the engine's, the
//! tolerance kept for the engine compared with itself, and the compared-case floors;
//! the grid under test as a [`subject::Subject`]; and the suite itself, which
//! [`oracle_suite!`] runs for one grid.
//!
//! Every helper here that asks the engine something makes exactly one oracle call per
//! question and never calls the oracle from inside another oracle call: the oracle runs
//! every request on a single worker thread, so a nested call would wait for itself.
#![allow(dead_code)]
pub mod aperture3;
pub mod boxes;
pub mod expected_divergences;
pub mod facts;
pub mod geometry;
pub mod sub_zone_order;
pub mod subject;
pub mod suite;
pub mod unresolved;
pub mod zones;

/// The whole oracle suite for one grid: one `#[test]` for each test of [`suite`] that
/// applies to the grid's aperture, each a call to it for the subject given, a type
/// implementing [`subject::Subject`]. The aperture is named at the call, as
/// `oracle_suite!(Subject, aperture 7)`, and each arm lists its tests; a test written for one
/// aperture only binds its subject's topology and indexing, so that listing it under the
/// other arm does not compile.
macro_rules! oracle_suite {
    ($subject:ty, aperture 7) => {
        crate::common::oracle_tests!($subject; suite:
            max_resolution_matches,
            the_three_engines_agree_on_the_icosahedron_vertices,
            zone_from_geo_text_id_and_resolution,
            seams_between_faces_quantise_like_dggal,
            centroids,
            vertices_and_pentagons,
            the_poles_snap_as_the_engine_does,
            deep_polar_digit_paths_answer_as_the_engine,
            neighbours_as_the_engine_lists_them,
            neighbours_near_the_poles_and_along_the_seams,
            parents_and_children_as_the_engine_lists_them,
            the_hierarchy_at_seam_identifiers_built_by_text,
            the_hierarchy_does_not_panic_on_hostile_input,
            disk_is_composed_neighbours,
            any_grid_answers_as_the_typed_grid_and_the_engine,
            sub_zones_at_depth_1_as_the_engine_lists_them,
            the_null_zone_matches_dggal_but_for_its_resolution,
            level_20_null_geometry_matches_dggal,
            refusals_depart_from_the_engine_as_recorded,
            listed_divergences_still_hold_with_their_evidence,
            quantize_matches_dggal_at_every_resolution,
        );
        crate::common::oracle_tests!($subject; sub_zone_order:
            sub_zone_orders_of_the_hexagons_as_the_engine_lists_them,
            sub_zone_orders_of_the_pentagons_as_the_engine_lists_them,
            deep_sub_zone_orders_as_the_engine_lists_them,
            sub_zone_orders_at_seam_identifiers_built_by_text,
            first_sub_zone_at_depth_0_departs_from_the_engine_as_recorded,
            the_sub_zone_methods_do_not_panic_on_hostile_input,
        );
        crate::common::oracle_tests!($subject; facts:
            the_area_of_every_sampled_zone_is_the_engines,
            an_identifier_with_no_cell_has_no_area,
            identifiers_the_engine_cannot_read_or_cannot_draw_have_the_areas_decided,
            the_grid_facts_are_the_engines,
            the_levels_from_an_area_are_the_engines,
        );
        crate::common::oracle_tests!($subject; geometry:
            the_refined_boundary_and_the_extent_are_the_engines,
            the_finest_refinements_are_the_engines,
            the_new_methods_do_not_panic_on_hostile_input,
        );
        crate::common::oracle_tests!($subject; zones:
            the_zones_of_a_box_are_the_engines_and_those_it_leaves_out,
            the_boxes_known_by_name_depart_from_the_engine_as_recorded,
            at_a_level_taken_whole_the_zones_of_a_box_are_those_of_the_rule,
            the_margin_of_the_test_of_two_extents_is_the_librarys,
            the_estimate_of_a_large_box_is_not_under_its_count,
            the_compaction_on_the_sample_is_the_engines,
            the_enumerations_the_estimate_and_the_compaction_withstand_hostile_input,
        );
    };
    ($subject:ty, aperture 3) => {
        crate::common::oracle_tests!($subject; suite:
            max_resolution_matches,
            zone_from_geo_text_id_and_resolution,
            seams_between_faces_quantise_like_dggal,
            centroids,
            vertices_and_pentagons,
            the_poles_snap_as_the_engine_does,
            neighbours_as_the_engine_lists_them,
            neighbours_near_the_poles_and_along_the_seams,
            disk_is_composed_neighbours,
            any_grid_answers_as_the_typed_grid_and_the_engine,
            listed_divergences_still_hold_with_their_evidence,
            quantize_matches_dggal_at_every_resolution,
        );
        crate::common::oracle_tests!($subject; aperture3:
            the_engine_places_its_pentagons_on_the_icosahedron_vertices,
            hierarchy_and_predicates_match_the_engine,
            sub_zones_and_paging_match_the_engine,
            the_null_zone_departs_from_the_engine_as_recorded,
            unreadable_answers_along_two_seams_are_the_null_zone,
            children_at_resolution_33_are_none,
            phantom_polar_texts_are_refused,
            sub_zone_answers_depart_from_the_engine_as_recorded,
            non_finite_coordinates_are_refused,
        );
        crate::common::oracle_tests!($subject; facts:
            the_area_of_every_sampled_zone_is_the_engines,
            an_identifier_with_no_cell_has_no_area,
            the_grid_facts_are_the_engines,
            the_levels_from_an_area_are_the_engines,
        );
        crate::common::oracle_tests!($subject; geometry:
            the_refined_boundary_and_the_extent_are_the_engines,
            the_finest_refinements_are_the_engines,
            the_new_methods_do_not_panic_on_hostile_input,
        );
        crate::common::oracle_tests!($subject; zones:
            the_zones_of_a_box_are_the_engines_and_those_it_leaves_out,
            the_boxes_known_by_name_depart_from_the_engine_as_recorded,
            at_a_level_taken_whole_the_zones_of_a_box_are_those_of_the_rule,
            the_margin_of_the_test_of_two_extents_is_the_librarys,
            the_estimate_of_a_large_box_is_not_under_its_count,
            the_compaction_on_the_sample_is_the_engines,
            the_enumerations_the_estimate_and_the_compaction_withstand_hostile_input,
        );
    };
}
pub(crate) use oracle_suite;

/// One `#[test]` for each test named, from the module named ([`suite`] or [`aperture3`]),
/// calling it for the subject given.
macro_rules! oracle_tests {
    ($subject:ty; $module:ident: $($name:ident),* $(,)?) => {
        $(
            #[test]
            fn $name() {
                crate::common::$module::$name::<$subject>()
            }
        )*
    };
}
pub(crate) use oracle_tests;

use std::sync::OnceLock;

use dggal_oracle as o;

/// A coordinate tolerance in degrees, for the comparisons that do not set this crate
/// against the engine: the engine against itself, as when two of its cells are asked
/// which vertices they share, or against a point derived from its own answers, as the
/// icosahedron's vertices are. Wherever this crate's coordinates meet the engine's, they
/// must be the same doubles: see [`assert_bit_identical`].
pub const TOL_DEG: f64 = 1e-9;

/// splitmix64: dependency-free, seeded and reproducible.
pub struct Rng(pub u64);

impl Rng {
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// A uniform value in `[0, 1)`.
    pub fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// Holds a test that builds lists of sub-zones to the limit under its ceiling: with the limit
/// lowered by the environment, a list that the suites need would be refused.
pub fn the_limit_on_a_list_of_sub_zones_is_its_ceiling() {
    assert_eq!(
        rs4dggs::grid::max_materialised_sub_zones(),
        rs4dggs::grid::MAX_MATERIALISED_SUB_ZONES,
        "the limit on a list of sub-zones is lowered: unset \
         RS4DGGS_MAX_MATERIALISED_SUB_ZONES to run the suites"
    );
}

/// Uniform points on the sphere, as (latitude, longitude) in degrees. Trigonometry in a
/// test is fine: the rule that confines it to `math.rs` covers the crate's `src/` only.
pub fn sphere_points(seed: u64, n: usize) -> Vec<(f64, f64)> {
    let mut r = Rng(seed);
    (0..n)
        .map(|_| {
            let lat = (2.0 * r.unit() - 1.0).asin().to_degrees();
            let lon = 360.0 * r.unit() - 180.0;
            (lat, lon)
        })
        .collect()
}

/// The inputs where porting defects hide. Always included, never skipped.
///
/// Both poles, from every side; the vertex-0 meridian at 11.2 degrees east and its
/// antimeridian, which together are the icosahedron edge that runs over the north pole;
/// the antimeridian itself; the twelve pentagon centres, four of which (bases 8 to 11)
/// set bit 63 of their identifiers; and two points a hair from the north pole, across
/// the polar edge, that lie on a cell boundary to within rounding (see
/// `expected_divergences`). The seams between the other faces come from
/// [`seam_points`].
pub fn adversarial_points() -> Vec<(f64, f64)> {
    let mut v = vec![
        (89.9999999, 0.0),
        (-89.9999999, 0.0),
        (89.9999999, 11.2),
        (-89.9999999, -168.8),
        (90.0, 0.0),
        (-90.0, 0.0),
        (0.0, 180.0),
        (0.0, -180.0),
        (45.0, 180.0),
        (-45.0, -180.0),
        (89.999_979_814_075_29, 101.199_999_744_738_58),
        (89.999_993_277_379_59, 101.199_999_234_215_7),
    ];
    for lat in [-80.0, -45.0, -10.0, 0.0, 10.0, 45.0, 80.0] {
        v.push((lat, 11.2)); // the vertex-0 meridian
        v.push((lat, -168.8)); // its antimeridian
        v.push((lat, 11.2 + 1e-12));
    }
    // The icosahedron vertices, that is the twelve pentagon centres: DGGAL's own
    // resolution-0 centroids.
    v.extend(icosahedron_vertices().iter().copied());
    v
}

/// Unit vector of a geographic point in degrees.
pub fn unit_vector(lat: f64, lon: f64) -> [f64; 3] {
    let (la, lo) = (lat.to_radians(), lon.to_radians());
    [la.cos() * lo.cos(), la.cos() * lo.sin(), la.sin()]
}

/// Geographic point in degrees of a (not necessarily unit) vector. The latitude comes
/// from an arc tangent rather than an arc sine, which near a pole would round
/// colatitudes below about 1e-6 degrees to nothing.
pub fn geographic(v: [f64; 3]) -> (f64, f64) {
    (
        v[2].atan2(v[0].hypot(v[1])).to_degrees(),
        v[1].atan2(v[0]).to_degrees(),
    )
}

/// The point `eps` degrees from `p` along the great circle towards `c`.
pub fn nudge(p: (f64, f64), c: (f64, f64), eps: f64) -> (f64, f64) {
    let (u, w) = (unit_vector(p.0, p.1), unit_vector(c.0, c.1));
    let t = eps / arc_deg(p, c);
    geographic([
        u[0] + t * (w[0] - u[0]),
        u[1] + t * (w[1] - u[1]),
        u[2] + t * (w[2] - u[2]),
    ])
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn normalised(a: [f64; 3]) -> [f64; 3] {
    let n = dot(a, a).sqrt();
    [a[0] / n, a[1] / n, a[2] / n]
}

/// Great-circle distance between two geographic points, in degrees.
pub fn arc_deg(a: (f64, f64), b: (f64, f64)) -> f64 {
    let (u, v) = (unit_vector(a.0, a.1), unit_vector(b.0, b.1));
    cross(u, v)
        .iter()
        .map(|c| c * c)
        .sum::<f64>()
        .sqrt()
        .atan2(dot(u, v))
        .to_degrees()
}

/// Distance in degrees from `p` to the great-circle arc between `a` and `b`: the distance
/// to the great circle where the foot of the perpendicular falls within the arc, and to
/// the nearer end otherwise.
pub fn dist_to_arc_deg(p: (f64, f64), a: (f64, f64), b: (f64, f64)) -> f64 {
    let (u, v, w) = (
        unit_vector(a.0, a.1),
        unit_vector(b.0, b.1),
        unit_vector(p.0, p.1),
    );
    let n = normalised(cross(u, v));
    let within = dot(cross(u, w), n) >= 0.0 && dot(cross(w, v), n) >= 0.0;
    if within {
        dot(w, n).abs().clamp(0.0, 1.0).asin().to_degrees()
    } else {
        arc_deg(p, a).min(arc_deg(p, b))
    }
}

/// The twelve icosahedron vertices, indexed by base cell: DGGAL's own resolution-0
/// centroids, asked of the engine once, one call per vertex.
///
/// They are IGEO7's, and serve every aperture-7 grid: the three share the orientation of
/// the icosahedron, and every grid's suite checks that the three engines place their base
/// cells on the same twelve points, bit for bit, rather than assuming it (see
/// `suite::the_three_engines_agree_on_the_icosahedron_vertices`).
pub fn icosahedron_vertices() -> &'static [(f64, f64); 12] {
    static V: OnceLock<[(f64, f64); 12]> = OnceLock::new();
    V.get_or_init(|| {
        let mut v = [(0.0, 0.0); 12];
        for (base, slot) in v.iter_mut().enumerate() {
            let z = o::zone_from_text(o::IGEO7, &format!("{base:02}"));
            *slot = o::centroid(o::IGEO7, z);
        }
        v
    })
}

/// The thirty icosahedron edges, as pairs of base cells whose vertices are adjacent.
pub fn icosahedron_edges() -> Vec<(usize, usize)> {
    let v = icosahedron_vertices();
    let mut e = Vec::new();
    for i in 0..12 {
        for j in (i + 1)..12 {
            // Adjacent vertices are 63.4 degrees apart; the next nearest are 116.6.
            if arc_deg(v[i], v[j]) < 90.0 {
                e.push((i, j));
            }
        }
    }
    assert_eq!(e.len(), 30, "an icosahedron has thirty edges");
    e
}

/// Points along every one of the thirty seams between projection faces: eleven stations
/// along each edge, each taken on the edge itself and at four small offsets across it
/// (1e-9 and 1e-6 radians either side), so that a point sits on the seam at the coarse
/// resolutions and just beside it at the fine ones. The station halfway along the edge
/// between bases 0 and 1 is the north pole.
pub fn seam_points() -> Vec<(f64, f64)> {
    let v = icosahedron_vertices();
    let mut pts = Vec::new();
    for (i, j) in icosahedron_edges() {
        let a = unit_vector(v[i].0, v[i].1);
        let b = unit_vector(v[j].0, v[j].1);
        let n = normalised(cross(a, b));
        for t in 1..12 {
            let f = f64::from(t) / 12.0;
            let p = normalised([
                a[0] * (1.0 - f) + b[0] * f,
                a[1] * (1.0 - f) + b[1] * f,
                a[2] * (1.0 - f) + b[2] * f,
            ]);
            for e in [-1e-6, -1e-9, 0.0, 1e-9, 1e-6] {
                pts.push(geographic([
                    p[0] + e * n[0],
                    p[1] + e * n[1],
                    p[2] + e * n[2],
                ]));
            }
        }
    }
    pts
}

/// Asserts that this crate's geographic point is the engine's to the last bit, latitude
/// and longitude alike, and records in `worst` the gap between the two.
///
/// The claim is exact equality on x86-64 Linux with glibc, where the engine and this
/// crate call the same `libm` and where the suites are run; nothing is claimed of any
/// other target. The gap is measured all the same, so that a failure says how far apart
/// the points lie and each test can print the largest it met, which is zero while the
/// claim holds. It is the great-circle distance between the points, not the larger of the
/// latitude and longitude differences: near a pole a longitude difference means almost
/// nothing, and at the pole itself nothing at all.
pub fn assert_bit_identical(what: &str, ours: (f64, f64), theirs: (f64, f64), worst: &mut f64) {
    let d = arc_deg(ours, theirs);
    *worst = worst.max(d);
    let bits = |p: (f64, f64)| (p.0.to_bits(), p.1.to_bits());
    let (o, t) = (bits(ours), bits(theirs));
    assert!(
        o == t,
        "{what}: ours {ours:?} ({:#018x}, {:#018x}), DGGAL {theirs:?} ({:#018x}, {:#018x}), \
         gap {d:e} degrees; the suite asserts the engine's own doubles, a claim made for \
         x86-64 Linux with glibc only",
        o.0,
        o.1,
        t.0,
        t.1
    );
}

/// Which way a ring turns, as seen from outside the sphere.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Winding {
    Anticlockwise,
    Clockwise,
    /// The ring encloses no area: its vertices lie on one point, or on one great circle.
    NoArea,
}

/// The winding of a ring of `(lat, lon)` vertices: the sign of the sum of
/// `(d_i x d_(i+1)) . m`, where `m` is the mean of the vertices as vectors and `d_i` is
/// vertex `i` less `m`. The sum equals that of `(v_i x v_(i+1)) . m`, but keeps its precision
/// in the smallest cells, some 1e-8 radians across, where products of unit vectors would
/// lose it to rounding. A vertex on a pole is the pole itself, whatever its longitude, so
/// that a ring whose every vertex the engine has snapped onto a pole is a single point. A
/// ring whose sum is not larger than a millionth of the largest squared `|d_i|` encloses no
/// area. Every ring with an area in a probe of some 490,000 zones on the six grids, from
/// resolution 0 to the finest, gave at least 0.98 of that square: six orders of magnitude
/// clear of the threshold.
pub fn winding(ring: &[(f64, f64)]) -> Winding {
    let at = |p: &(f64, f64)| {
        if p.0.abs() == 90.0 {
            [0.0, 0.0, p.0.signum()]
        } else {
            unit_vector(p.0, p.1)
        }
    };
    let v: Vec<[f64; 3]> = ring.iter().map(at).collect();
    let n = v.len() as f64;
    let m = [0, 1, 2].map(|k| v.iter().map(|p| p[k]).sum::<f64>() / n);
    let d: Vec<[f64; 3]> = v.iter().map(|p| [0, 1, 2].map(|k| p[k] - m[k])).collect();
    let s: f64 = (0..d.len())
        .map(|i| dot(cross(d[i], d[(i + 1) % d.len()]), m))
        .sum();
    let r2 = d.iter().map(|p| dot(*p, *p)).fold(0.0, f64::max);
    if s.abs() <= 1e-6 * r2 {
        Winding::NoArea
    } else if s > 0.0 {
        Winding::Anticlockwise
    } else {
        Winding::Clockwise
    }
}

/// Whether `id` is an aperture-3 polar pentagon: root 10 or 11, in the fields of the I3H
/// packing.
pub fn is_polar_pentagon(id: rs4dggs::ZoneId) -> bool {
    (id.0 >> 53) & 0xF >= 10
}

/// The engine's ring of zone `id` in this crate's order. On aperture 3 the engine's
/// `[v0, v1, ..., v(n-1)]` is `[v0, v(n-1), ..., v1]`, turned about its first vertex, for
/// every zone but the two polar pentagons, whose ring keeps the engine's order, as does
/// every ring on aperture 7.
pub fn in_crate_order<S: subject::Subject>(
    id: rs4dggs::ZoneId,
    mut ring: Vec<(f64, f64)>,
) -> Vec<(f64, f64)> {
    use rs4dggs::Topology;
    if S::T::APERTURE == 3 && !is_polar_pentagon(id) {
        ring[1..].reverse();
    }
    ring
}

/// The rings a test compares with the engine's, counted by the way they turn.
#[derive(Default)]
pub struct Windings {
    /// Rings compared.
    pub rings: usize,
    /// Aperture-7 rings whose every vertex lies on a pole, which enclose no area.
    pub no_area: usize,
    /// Aperture-7 rings among whose vertices the engine puts the point (0, 0).
    pub zero_vertex: usize,
    /// Aperture-3 polar pentagons, at which alone the engine's ring runs anticlockwise.
    pub polar: usize,
}

impl Windings {
    /// The engine's ring of `id` in this crate's order ([`in_crate_order`]), for the caller to
    /// compare with this crate's ring vertex for vertex. One oracle call.
    ///
    /// Asserts that it runs anticlockwise, with two exceptions on aperture 7 alone, each
    /// counted: a ring whose every vertex the engine snaps onto a pole, which encloses no
    /// area; and a ring among whose vertices the engine puts the point (0, 0), under the
    /// pentagon of base cell 0 at resolutions 15, 17 and 19, where it answers a planar vertex
    /// on the diagonal of the 5x6 layout with the zero point. The rest of such a ring must
    /// run anticlockwise. On aperture 3 it asserts too that the engine's own ring runs
    /// anticlockwise at a polar pentagon and clockwise at every other zone, the exception
    /// the crate's rule rests on, so that a change in the engine's order fails here. The null
    /// zone has no ring, and its six zero points are given as the engine gives them.
    pub fn engine_ring<S: subject::Subject>(&mut self, id: rs4dggs::ZoneId) -> Vec<(f64, f64)> {
        use rs4dggs::Topology;
        let engine = o::vertices(S::ORACLE, id.0);
        if id == rs4dggs::ZoneId::NULL {
            return engine;
        }
        let text = S::grid().text_id(id);
        if S::T::APERTURE == 3 {
            let polar = is_polar_pentagon(id);
            let expected = if polar {
                Winding::Anticlockwise
            } else {
                Winding::Clockwise
            };
            assert_eq!(
                winding(&engine),
                expected,
                "the engine's own ring of {text}, {}a polar pentagon: {engine:?}",
                if polar { "" } else { "not " }
            );
            self.polar += usize::from(polar);
        }
        let ring = in_crate_order::<S>(id, engine);
        let mut rest = ring.clone();
        if S::T::APERTURE == 7 && ring.contains(&(0.0, 0.0)) {
            assert!(
                text.starts_with("00"),
                "the engine puts the point (0, 0) among the vertices of {text}, which is not \
                 under base cell 0: {ring:?}"
            );
            let resolution = S::grid().resolution(id);
            let zeros = ring.iter().filter(|v| **v == (0.0, 0.0)).count();
            assert!(
                matches!(resolution, 15 | 17 | 19) && matches!(zeros, 1 | 3),
                "the point (0, 0) among the vertices of {text}: resolution {resolution}, \
                 {zeros} such vertices, where 15, 17 or 19 and 1 or 3 were measured: {ring:?}"
            );
            rest.retain(|v| *v != (0.0, 0.0));
            self.zero_vertex += 1;
        }
        match winding(&rest) {
            Winding::Anticlockwise => {}
            Winding::NoArea if S::T::APERTURE == 7 && rest.iter().all(|v| v.0.abs() == 90.0) => {
                self.no_area += 1;
            }
            w => panic!("the ring of {text} in this crate's order runs {w:?}: {ring:?}"),
        }
        self.rings += 1;
        ring
    }

    /// Prints the counts under `what`, and asserts that the test met at least `polar_min`
    /// polar pentagons on aperture 3, so that both sides of the engine's exception are seen,
    /// and, on aperture 7, exactly `no_area` rings with every vertex on a pole and exactly
    /// `zero_vertex` rings with the point (0, 0), so that neither class can grow nor vanish
    /// unnoticed.
    pub fn finish<S: subject::Subject>(
        &self,
        what: &str,
        polar_min: usize,
        no_area: usize,
        zero_vertex: usize,
    ) {
        use rs4dggs::Topology;
        if S::T::APERTURE == 3 {
            assert!(
                self.polar >= polar_min,
                "{what}: {} polar pentagons met, floor is {polar_min}",
                self.polar
            );
            eprintln!(
                "{what}: {} rings anticlockwise in this crate's order; the engine's own ring \
                 anticlockwise at {} polar pentagons and clockwise at the {} other zones",
                self.rings,
                self.polar,
                self.rings - self.polar
            );
        } else {
            assert_eq!(
                self.no_area, no_area,
                "{what}: rings with every vertex on a pole"
            );
            assert_eq!(
                self.zero_vertex, zero_vertex,
                "{what}: rings with the point (0, 0) among their vertices"
            );
            eprintln!(
                "{what}: {} rings compared, {} of them with every \
                 vertex on a pole and no area, and {} with the point (0, 0) among their \
                 vertices, anticlockwise without it",
                self.rings, self.no_area, self.zero_vertex
            );
        }
    }
}

/// Count floor: every oracle test proves how much it compared, so that coverage cannot
/// quietly shrink.
pub struct Floor {
    name: &'static str,
    min: usize,
    pub n: usize,
}

impl Floor {
    pub fn new(name: &'static str, min: usize) -> Self {
        Floor { name, min, n: 0 }
    }

    pub fn hit(&mut self) {
        self.n += 1;
    }
}

impl Drop for Floor {
    fn drop(&mut self) {
        if !std::thread::panicking() {
            assert!(
                self.n >= self.min,
                "{}: compared {} cases, floor is {}",
                self.name,
                self.n,
                self.min
            );
            eprintln!("{}: {} cases", self.name, self.n);
        }
    }
}
