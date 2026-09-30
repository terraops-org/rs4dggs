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
