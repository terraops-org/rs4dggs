use dggal_oracle::*;

#[test]
fn lisbon_res1_is_006() {
    let z = zone_from_geo(IGEO7, 38.7223, -9.1393, 1);
    assert_eq!(text_id(IGEO7, z), "006");
}

#[test]
fn max_level_is_19() {
    assert_eq!(max_level(IGEO7), 19);
}

#[test]
fn bit63_zone_neighbours_survive() {
    // Base cells 8-11 set bit 63: the pydggal int64 bug. u64 all the way here.
    let z = zone_from_text(IGEO7, "1001234");
    assert!(z >> 63 == 1, "expected a bit-63 zone, got {z:#x}");
    let n = neighbors(IGEO7, z);
    assert_eq!(n.len(), 6);
    assert!(n.iter().all(|&v| v != NULL_ZONE));
}

#[test]
fn centroid_round_trips_through_zone_from_geo() {
    let z = zone_from_text(IGEO7, "0064156");
    let (lat, lon) = centroid(IGEO7, z);
    assert!((lat - 38.61687542349096).abs() < 1e-9 && (lon - -8.874508932649084).abs() < 1e-9);
    assert_eq!(zone_from_geo(IGEO7, lat, lon, 5), z);
    assert_eq!(vertices(IGEO7, z).len(), 6);
}

#[test]
fn the_three_aperture3_grids_quantise_lisbon() {
    for grid in [ISEA3H, IVEA3H, RTEA3H] {
        let z = zone_from_geo(grid, 38.7223, -9.1393, 1);
        assert_ne!(
            z, NULL_ZONE,
            "{grid}: Lisbon at level 1 quantised to NULL_ZONE"
        );
        assert!(
            !text_id(grid, z).is_empty(),
            "{grid}: empty text id for zone {z:#x}"
        );
    }
}

#[test]
fn zone_from_5x6_reads_as_b6_5_a() {
    let z = zone_from_5x6(ISEA3H, 2, 3.6666666666666665, 3.333333333333333);
    assert_eq!(text_id(ISEA3H, z), "B6-5-A");
}

#[test]
fn crs_centroid_5x6_round_trips_bit_for_bit() {
    let z = zone_from_text(ISEA3H, "B6-5-A");
    let (x, y) = crs_centroid_5x6(ISEA3H, z);
    assert_eq!(
        (x.to_bits(), y.to_bits()),
        (
            3.6666666666666665f64.to_bits(),
            3.333333333333333f64.to_bits()
        )
    );
}

#[test]
fn sub_zone_crs_centroids_5x6_are_what_get_sub_zones_quantises() {
    // An even and an odd parent, a hexagon and a pentagon; one depth of each parity.
    for (text, depth) in [("B6-5-A", 3), ("A6-0-C", 2), ("A5-0-B", 3)] {
        let z = zone_from_text(ISEA3H, text);
        let centroids = sub_zone_crs_centroids_5x6(ISEA3H, z, depth);
        let subs = sub_zones(ISEA3H, z, depth);
        assert_eq!(centroids.len() as u64, count_sub_zones(ISEA3H, z, depth));
        assert_eq!(centroids.len(), subs.len(), "{text} at {depth}");
        let sz_level = level(ISEA3H, z) + depth;
        for (i, (&(x, y), &sub)) in centroids.iter().zip(&subs).enumerate() {
            assert_eq!(
                zone_from_5x6(ISEA3H, sz_level, x, y),
                sub,
                "{text} at {depth}, [{i}]"
            );
        }
    }
    // At depth 0, the zone's own centroid, bit for bit.
    let z = zone_from_text(ISEA3H, "B6-5-A");
    let (x, y) = crs_centroid_5x6(ISEA3H, z);
    let at0: Vec<(u64, u64)> = sub_zone_crs_centroids_5x6(ISEA3H, z, 0)
        .iter()
        .map(|p| (p.0.to_bits(), p.1.to_bits()))
        .collect();
    assert_eq!(at0, vec![(x.to_bits(), y.to_bits())]);
}

#[test]
fn engine_can_read_guards_the_null_zone_and_the_unreadable_north_pole_answer() {
    assert!(!engine_can_read(ISEA3H, NULL_ZONE));

    // The unreadable north-pole class's resolution-2 identifier: the engine's own text of
    // it is rejected by getZoneFromTextID (root above 9 requires index 0), so the guard
    // must refuse it even though its level, taken alone, is well within bounds.
    let north_pole = 0x0340_0000_0000_0008u64;
    assert_eq!(text_id(ISEA3H, north_pole), "BA-2-A");
    assert!(!engine_can_read(ISEA3H, north_pole));

    let z = zone_from_text(ISEA3H, "B6-5-A");
    assert!(engine_can_read(ISEA3H, z));
}

#[test]
fn aperture3_relational_predicates_are_self_consistent() {
    let parent = zone_from_text(ISEA3H, "B6-5-A");
    assert!(engine_can_read(ISEA3H, parent));

    let kids = children(ISEA3H, parent);
    assert!(
        kids.len() >= 2,
        "expected at least two children of {parent:#x} to compare as siblings, got {kids:?}"
    );
    for &k in &kids {
        assert!(engine_can_read(ISEA3H, k), "child {k:#x} must be readable");
    }

    let cc = centroid_child(ISEA3H, parent);
    assert!(
        kids.contains(&cc),
        "the centroid child must be among the engine's children"
    );
    assert!(zone_has_sub_zone(ISEA3H, parent, cc));
    assert!(is_immediate_child_of(ISEA3H, cc, parent));
    assert!(is_ancestor_of(ISEA3H, parent, cc));

    let sibling_a = kids[0];
    let sibling_b = kids[1];
    assert!(are_siblings(ISEA3H, sibling_a, sibling_b));
    assert!(!are_siblings(ISEA3H, parent, sibling_a));

    let vs = crs_vertices_5x6(ISEA3H, parent);
    assert!(
        vs.len() == 5 || vs.len() == 6,
        "expected a pentagon or hexagon, got {} vertices",
        vs.len()
    );
}

#[test]
fn worker_survives_a_panicking_job() {
    // An embedded NUL byte panics inside the vendored binding's own
    // `CString::new(zoneID).unwrap()` (vendor/dggal/src/lib.rs, getZoneFromTextID): a panic
    // raised in Rust, exactly what with()'s catch_unwind/resume_unwind is meant to survive.
    // If a future DGGAL binding version replaces that unwrap with a Result, this trigger
    // stops panicking and the assertion below fails loudly, which is the right outcome:
    // a silently-vanished trigger would otherwise stop proving anything.
    let result = std::panic::catch_unwind(|| zone_from_text(IGEO7, "\0"));
    assert!(result.is_err(), "expected the embedded NUL byte to panic");

    // The worker thread, and every later oracle call in this binary, must still work.
    assert_eq!(max_level(IGEO7), 19);
    assert_eq!(
        zone_from_geo(IGEO7, 38.7223, -9.1393, 1),
        zone_from_text(IGEO7, "006")
    );
}

// The literals below are read from the engine, measured
// live against libdggal.so, BuildID e75d6ab1...: the wrappers are thin, and these pin that they
// hand back what the engine gives, in radians and square metres, bit for bit.

#[test]
fn the_area_of_a_zone_is_the_engines_own_bits() {
    // (grid, zone, level, area bits): a hexagon and a pentagon at level 5 on each aperture.
    for (grid, text, lv, bits) in [
        (IGEO7, "0313522", 5, 0x41e69c7fa6cc8c7fu64),
        (IGEO7, "0100000", 5, 0x41e2d7bfb5aa7515),
        (ISEA3H, "C4-49-B", 5, 0x42486f9cfb414d5e),
        (ISEA3H, "C0-0-B", 5, 0x42445d02d16115cf),
    ] {
        let z = zone_from_text(grid, text);
        assert_eq!(level(grid, z), lv, "{grid} {text}");
        assert_eq!(area(grid, z).to_bits(), bits, "{grid} {text}");
    }
    // Without geometry: the null zone's area is infinite on every grid; an unreadable north-pole
    // identifier of aperture 3 has the ordinary area of its level.
    for grid in [IGEO7, IVEA7H, RTEA7H, ISEA3H, IVEA3H, RTEA3H] {
        assert_eq!(area(grid, NULL_ZONE), f64::INFINITY, "{grid}");
    }
    assert_eq!(
        area(IGEO7, 0x329cb814e5c0a72e).to_bits(),
        0x3f44f2599691f084
    );
    // The ordinary hexagon area of its level, 2, though the engine cannot read it.
    assert_eq!(
        area(ISEA3H, 0x0340_0000_0000_0008).to_bits(),
        0x42949e2c73ff1947
    );
}

#[test]
fn the_refined_ring_has_the_engines_number_of_points() {
    // One zone of each kind at level 5; columns are the refinements
    // 0, 1, 2, 3, 4, 12 and 20.
    let refinements = [0, 1, 2, 3, 4, 12, 20];
    for (grid, text, counts) in [
        (IGEO7, "0313522", [72, 6, 12, 18, 24, 72, 120]),
        (IGEO7, "0100000", [61, 5, 11, 16, 21, 61, 101]),
        (ISEA3H, "C4-49-B", [60, 6, 12, 18, 24, 72, 120]),
        (ISEA3H, "C0-0-B", [50, 5, 10, 15, 20, 60, 100]),
    ] {
        let z = zone_from_text(grid, text);
        for (n, want) in refinements.into_iter().zip(counts) {
            let ring = refined_vertices(grid, z, n);
            assert_eq!(ring.len(), want, "{grid} {text} at refinement {n}");
            // Radians, as the engine gives them.
            for &(lat, lon) in &ring {
                assert!(lat.abs() <= std::f64::consts::FRAC_PI_2 + 1e-9, "{lat}");
                assert!(lon.abs() <= 2.0 * std::f64::consts::PI, "{lon}");
            }
        }
    }
}

#[test]
fn the_ring_of_an_aperture7_identifier_the_engine_cannot_read_is_still_asked() {
    // `01222222222222222220`: the engine writes it for a neighbour and cannot read it back, yet
    // gives it an ordinary hexagon; a Z7 identifier of twenty digits has an empty ring.
    let unreadable = 0x149249249249243fu64;
    assert!(!engine_can_read(IGEO7, unreadable));
    assert_eq!(refined_vertices(IGEO7, unreadable, 0).len(), 72);
    assert!(refined_vertices(IGEO7, 0x329cb814e5c0a72e, 0).is_empty());
    assert!(refined_vertices(IGEO7, NULL_ZONE, 0).is_empty());
}

#[test]
#[should_panic(expected = "BA-2-A")]
fn the_ring_of_an_aperture3_identifier_the_engine_cannot_read_is_refused() {
    // The engine answers a null array, and the vendored binding dereferences it.
    refined_vertices(ISEA3H, 0x0340_0000_0000_0008, 0);
}

#[test]
#[should_panic(expected = "negative")]
fn a_negative_refinement_is_refused() {
    // The engine answers 0, 1 or 2 points from a null backing array, which the vendored
    // binding's `slice::from_raw_parts` aborts the debug build on: no panic, no recovery.
    refined_vertices(IGEO7, zone_from_text(IGEO7, "0313522"), -1);
}

#[test]
fn the_extent_of_a_zone_without_geometry_is_the_cleared_extent() {
    // ll = (X, X), ur = (-X, -X), X = 0x7f91df46a2529d38 radians.
    let x = f64::from_bits(0x7f91df46a2529d38);
    let cleared = [x, x, -x, -x].map(f64::to_bits);
    for (grid, zone) in [
        (IGEO7, NULL_ZONE),
        (IGEO7, 0x0000_0000_0000_0000),
        (IGEO7, 0x329cb814e5c0a72e),
        (ISEA3H, NULL_ZONE),
        // Asked of the engine on aperture 3, where only the ring is unsafe.
        (ISEA3H, 0x0340_0000_0000_0008),
    ] {
        assert_eq!(
            extent(grid, zone).map(f64::to_bits),
            cleared,
            "{grid} {zone:#x}"
        );
    }
}

#[test]
fn the_extent_of_a_zone_is_the_engines_own_in_radians() {
    // The unreadable aperture-7 neighbour has an ordinary finite extent, read from the engine:
    // 85.438249661, -88.603247248 to 85.438251704, -88.603225043 degrees.
    let [a, b, c, d] = extent(IGEO7, 0x149249249249243f);
    for (got, want) in [
        (a, 85.438249661),
        (b, -88.603247248),
        (c, 85.438251704),
        (d, -88.603225043),
    ] {
        assert!(
            (got.to_degrees() - want).abs() < 1e-8,
            "{got} against {want}"
        );
    }
    // The extent is the loop over the ring at refinement 0: the four bounds.
    let z = zone_from_text(IGEO7, "0313522");
    let ring = refined_vertices(IGEO7, z, 0);
    let fold = |pick: fn(&(f64, f64)) -> f64, op: fn(f64, f64) -> f64, start| {
        ring.iter().map(pick).fold(start, op)
    };
    let want = [
        fold(|p| p.0, f64::min, f64::INFINITY),
        fold(|p| p.1, f64::min, f64::INFINITY),
        fold(|p| p.0, f64::max, f64::NEG_INFINITY),
        fold(|p| p.1, f64::max, f64::NEG_INFINITY),
    ];
    assert_eq!(extent(IGEO7, z).map(f64::to_bits), want.map(f64::to_bits));
}

#[test]
fn the_grid_facts_are_the_engines_own() {
    // Note 4.4: ratio, max parents, max children, max neighbours, 64K depth, max depth,
    // countZones(max), getRefZoneArea(0) and (max); the three grids of an aperture agree.
    let a7 = (
        7,
        2,
        13,
        6,
        6,
        9,
        19,
        113988951853731432u64,
        0x3f72540e63bfb274u64,
    );
    let a3 = (
        3,
        3,
        7,
        6,
        10,
        17,
        33,
        55590605665555232u64,
        0x3f82ca8c52e6248au64,
    );
    for (grids, (ratio, parents, children, neighbours, d64k, depth, max, count, area_max)) in [
        ([IGEO7, IVEA7H, RTEA7H], a7),
        ([ISEA3H, IVEA3H, RTEA3H], a3),
    ] {
        for grid in grids {
            assert_eq!(refinement_ratio(grid), ratio, "{grid}");
            assert_eq!(max_parents(grid), parents, "{grid}");
            assert_eq!(max_children(grid), children, "{grid}");
            assert_eq!(max_neighbors(grid), neighbours, "{grid}");
            assert_eq!(depth_64k(grid), d64k, "{grid}");
            assert_eq!(max_depth(grid), depth, "{grid}");
            assert_eq!(max_level(grid), max, "{grid}");
            assert_eq!(count_zones(grid, 0), 12, "{grid}");
            assert_eq!(count_zones(grid, max), count, "{grid}");
            assert_eq!(count_zones(grid, -1), 2, "{grid}");
            assert_eq!(
                ref_zone_area(grid, 0).to_bits(),
                0x42c35449acbf27b3,
                "{grid}"
            );
            assert_eq!(ref_zone_area(grid, max).to_bits(), area_max, "{grid}");
            assert_eq!(
                ref_zone_area(grid, -1).to_bits(),
                0x42ecfe6e831ebb8c,
                "{grid}"
            );
        }
    }
    // Levels 1 and 9 of aperture 7, as the engine counts them.
    assert_eq!(count_zones(IGEO7, 1), 72);
    assert_eq!(ref_zone_area(IGEO7, 1).to_bits(), 0x4299c5b790fedf99);
    assert_eq!(count_zones(IGEO7, 9), 403536072);
}

#[test]
fn the_levels_from_an_area_and_from_a_size_are_the_engines_own() {
    // Note 4.4: not the inverse of getRefZoneArea; on ISEA3H the area of level 9 answers 10.
    assert_eq!(level_from_ref_zone_area(IGEO7, ref_zone_area(IGEO7, 1)), 1);
    assert_eq!(level_from_ref_zone_area(IGEO7, ref_zone_area(IGEO7, 9)), 9);
    assert_eq!(
        level_from_ref_zone_area(ISEA3H, ref_zone_area(ISEA3H, 9)),
        10
    );
    // 25 or 43 for an area of 0 or NaN; 0 for a negative or an infinite one.
    for (grid, top) in [(IGEO7, 25), (ISEA3H, 43)] {
        assert_eq!(level_from_ref_zone_area(grid, 0.0), top, "{grid}");
        assert_eq!(level_from_ref_zone_area(grid, f64::NAN), top, "{grid}");
        assert_eq!(level_from_ref_zone_area(grid, -1.0), 0, "{grid}");
        assert_eq!(level_from_ref_zone_area(grid, f64::INFINITY), 0, "{grid}");
    }
    // Metres per sub-zone: sqrt of the reference area of the sub-zone level, and
    // the level is the one of the squared size, less the depth, and never below 0.
    for grid in [IGEO7, ISEA3H] {
        for (lv, depth) in [(0, 0), (3, 2), (7, 5), (max_level(grid), 0)] {
            let m = meters_per_sub_zone(grid, lv, depth);
            assert_eq!(
                m.to_bits(),
                ref_zone_area(grid, lv + depth).sqrt().to_bits(),
                "{grid} {lv} {depth}"
            );
            assert_eq!(
                level_from_meters_per_sub_zone(grid, m, depth),
                (level_from_ref_zone_area(grid, m * m) - depth).max(0),
                "{grid} {lv} {depth}"
            );
        }
    }
    assert_eq!(meters_per_sub_zone(IGEO7, 19, 0), 6.689314985989693e-2);
    assert_eq!(meters_per_sub_zone(ISEA3H, 33, 0), 9.578826868562763e-2);
}

#[test]
fn the_planar_refined_ring_takes_its_refinement() {
    let z = zone_from_text(IGEO7, "0313522");
    // The edge refinement of one, which every earlier caller passes, is the hexagon itself.
    assert_eq!(refined_crs_vertices_5x6(IGEO7, z, 1).len(), 6);
    assert_eq!(refined_crs_vertices_5x6(IGEO7, z, 3).len(), 18);
}

/// A box made in degrees, handed over as the crate converts: the product with `Pi / 180`.
fn deg(b: [f64; 4]) -> [f64; 4] {
    b.map(|x| x * (std::f64::consts::PI / 180.0))
}

fn ids(grid: &str, texts: &[&str]) -> Vec<u64> {
    texts.iter().map(|t| zone_from_text(grid, t)).collect()
}

fn texts(grid: &str, zones: &[u64]) -> Vec<String> {
    zones.iter().map(|&z| text_id(grid, z)).collect()
}

#[test]
fn lisbon_box_on_isea3h() {
    let b = Some(deg([38.0, -10.0, 39.0, -9.0]));
    // The engine omits the zone that holds the box at level 2, and one of two at level 3.
    assert!(list_zones(ISEA3H, 2, b).is_empty());
    assert_eq!(texts(ISEA3H, &list_zones(ISEA3H, 3, b)), ["B2-5-C"]);
    assert_eq!(list_zones(ISEA3H, 6, b).len(), 3);
}

#[test]
fn level_zero_of_igeo7_in_the_engines_order() {
    let z = list_zones(IGEO7, 0, None);
    assert_eq!(
        texts(IGEO7, &z),
        ["01", "06", "02", "07", "03", "08", "04", "09", "05", "10", "00", "11"]
    );
}

/// `None` is the whole world: every zone of the level, as many as the engine counts.
#[test]
fn the_whole_world_is_every_zone_of_the_level() {
    for (grid, level, count) in [(ISEA3H, 4, 812), (IGEO7, 2, 492)] {
        let a = list_zones(grid, level, None);
        assert_eq!(a.len(), count);
        assert_eq!(a.len() as u64, count_zones(grid, level));
    }
}

#[test]
fn compaction_of_the_sub_zones_of_one_zone() {
    let sub = sub_zones(ISEA3H, zone_from_text(ISEA3H, "B4-2-A"), 2);
    assert_eq!(sub.len(), 13);
    let c = compact_zones(ISEA3H, &sub);
    assert_eq!(c.len(), 7);
    let mut want = ids(
        ISEA3H,
        &[
            "B4-2-A", "C2-11-A", "C2-19-A", "C2-2C-A", "C4-E-A", "C4-11-A", "C4-19-A",
        ],
    );
    want.sort_unstable();
    assert_eq!(c, want);

    let sub7 = sub_zones(IGEO7, zone_from_text(IGEO7, "0064"), 1);
    let c7 = compact_zones(IGEO7, &sub7);
    assert_eq!((sub7.len(), c7.len()), (13, 13));
    assert!(c7.contains(&zone_from_text(IGEO7, "0064")));
    assert!(!c7.contains(&zone_from_text(IGEO7, "00640")));
}

#[test]
fn compaction_of_mixed_levels_loses_the_coarser_zone() {
    let c = compact_zones(IGEO7, &ids(IGEO7, &["00", "0064"]));
    assert_eq!(texts(IGEO7, &c), ["0064"]);
}

#[test]
#[should_panic(expected = "a coordinate that is not finite")]
fn list_zones_refuses_a_box_that_is_not_finite() {
    list_zones(ISEA3H, 4, Some([0.0, 0.0, f64::INFINITY, 0.1]));
}

#[test]
#[should_panic(expected = "a longitude beyond Pi")]
fn list_zones_refuses_a_longitude_beyond_pi() {
    list_zones(IGEO7, 2, Some(deg([10.0, 170.0, 20.0, -181.0])));
}

#[test]
#[should_panic(expected = "a sliver")]
fn list_zones_refuses_a_sliver() {
    list_zones(
        ISEA3H,
        10,
        Some([
            f64::from_bits(0x3fe6666666666666),
            f64::from_bits(0x4008000000000000),
            f64::from_bits(0x3fe6666666ab1dd0),
            f64::from_bits(0x4008000000000063),
        ]),
    );
}

#[test]
#[should_panic(expected = "a box of one point on a pole")]
fn list_zones_refuses_a_point_on_a_pole() {
    list_zones(
        IGEO7,
        3,
        Some([
            std::f64::consts::FRAC_PI_2,
            0.0,
            std::f64::consts::FRAC_PI_2,
            0.0,
        ]),
    );
}

#[test]
#[should_panic(expected = "a box of no height on a pole or of no width on the antimeridian")]
fn list_zones_refuses_a_box_of_no_height_on_a_pole_on_aperture_7() {
    list_zones(IGEO7, 2, Some(deg([90.0, 0.0, 90.0, 10.0])));
}

#[test]
#[should_panic(expected = "a box of no height on a pole or of no width on the antimeridian")]
fn list_zones_refuses_a_box_of_no_width_on_the_antimeridian_on_aperture_7() {
    list_zones(IGEO7, 5, Some(deg([10.0, -180.0, 20.0, -180.0])));
}

#[test]
#[should_panic(expected = "a box of no height on a pole or of no width on the antimeridian")]
fn list_zones_refuses_a_box_of_no_width_over_the_antimeridian_on_aperture_7() {
    list_zones(IGEO7, 3, Some(deg([10.0, 180.0, 20.0, -180.0])));
}

/// What aperture 7 is refused, aperture 3 answers, with an empty array or a zone; and beside the
/// refused boxes aperture 7 answers too: a box of no width a millionth of a degree within the
/// antimeridian, one that reaches it from either side, a box of no height below the pole.
#[test]
fn the_flat_boxes_that_the_engine_answers() {
    assert!(list_zones(ISEA3H, 2, Some(deg([90.0, 0.0, 90.0, 10.0]))).is_empty());
    assert!(list_zones(ISEA3H, 5, Some(deg([10.0, -180.0, 20.0, -180.0]))).is_empty());
    assert_eq!(
        list_zones(ISEA3H, 2, Some(deg([10.0, 180.0, 10.0, 180.0]))).len(),
        1
    );
    for (level, b, count) in [
        (3, [10.0, 179.999999, 20.0, 179.999999], 5),
        (3, [10.0, 179.0, 20.0, 180.0], 5),
        (3, [10.0, 180.0, 20.0, -179.0], 8),
        (3, [10.0, 179.0, 20.0, -180.0], 5),
        (4, [89.999999, 0.0, 89.999999, 10.0], 1),
        (3, [10.0, 0.0, 20.0, 0.0], 5),
    ] {
        assert_eq!(list_zones(IGEO7, level, Some(deg(b))).len(), count, "{b:?}");
    }
}

#[test]
#[should_panic(expected = "latitude inverted")]
fn list_zones_refuses_inverted_latitudes() {
    list_zones(IGEO7, 3, Some(deg([40.0, 10.0, 30.0, 20.0])));
}

#[test]
#[should_panic(expected = "outside the grid's")]
fn list_zones_refuses_a_level_beyond_the_grids() {
    list_zones(IGEO7, 20, None);
}
