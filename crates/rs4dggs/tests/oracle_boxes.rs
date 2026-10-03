//! The boxes the suites put to the engine's `listZones`, and the test of which of them it may be
//! asked, shown to do their work against the live engine.
// The suites that follow use what this one does not.
#[allow(dead_code)]
#[path = "common/boxes.rs"]
mod boxes;

use boxes::{
    GRIDS, KINDS, Rng, engine_survives, is_aperture7, make_box, radians, roots, seam, sites,
    zone_width_deg,
};
use dggal_oracle::{IGEO7, ISEA3H, IVEA3H, IVEA7H, RTEA3H, RTEA7H, list_zones};

/// A box in degrees about `(lat, lon)`, `half` zone widths of `level` to each side on the ground.
fn about(grid: &str, level: i32, (lat, lon): (f64, f64), half: f64) -> [f64; 4] {
    let h = half * zone_width_deg(grid, level);
    let w = h / lat.to_radians().cos();
    [lat - h, lon - w, lat + h, lon + w]
}

/// The box of the polar cap that takes the engine without bound on aperture 3 at a deep level,
/// and a box of one degree on aperture 7 at level 14, which holds thousands of millions of
/// zones, are not put to the engine; the box about Lisbon is.
#[test]
fn the_engine_is_not_asked_what_it_does_not_survive() {
    let cap = [89.997928885, 11.131869680, 90.0, 108.081826856];
    assert!(!engine_survives(ISEA3H, 20, &cap));
    assert!(engine_survives(ISEA3H, 8, &cap));
    assert!(!engine_survives(IGEO7, 14, &[38.0, -10.0, 39.0, -9.0]));
    assert!(engine_survives(ISEA3H, 6, &[38.0, -10.0, 39.0, -9.0]));
    assert!(engine_survives(IGEO7, 5, &[38.0, -10.0, 39.0, -9.0]));
    // What the wrapper refuses with a panic is never let through.
    for (level, b) in [
        (4, [0.0, 0.0, f64::INFINITY, 5.73]),
        (4, [10.0, 170.0, 20.0, -181.0]),
        (4, [40.0, 10.0, 30.0, 20.0]),
        (4, [95.0, 10.0, 96.0, 20.0]),
        (4, [90.0, 0.0, 90.0, 0.0]),
        (4, [40.0, 10.0, 40.0, 10.0 + 1e-10]),
        (20, [0.0, 0.0, 1.0, 1.0]),
        (-1, [0.0, 0.0, 1.0, 1.0]),
    ] {
        assert!(!engine_survives(IGEO7, level, &b), "{b:?} at level {level}");
    }
    // A long thin box holds few zones and still costs the engine the rectangle it spans.
    assert!(!engine_survives(ISEA3H, 14, &[0.0, 0.0, 0.1, 100.0]));
    // On aperture 7 the engine is asked to the level at which the crate's box query ends.
    assert!(engine_survives(IGEO7, 14, &[89.9999997, -32.3, 90.0, 36.5]));
    assert!(!engine_survives(
        IGEO7,
        15,
        &[89.9999997, -32.3, 90.0, 36.5]
    ));
    assert!(engine_survives(IGEO7, 14, &[0.0, 0.0, 0.01, 0.01]));
}

/// On aperture 7 the engine answers a box that no extent can meet with a null array, at every
/// level, and the asking process aborts: a box of no height on a pole, a box of no width, or a
/// point, on the antimeridian, and the same within the engine's margin of either. They are kept
/// from it; their neighbours are not, and the engine answers those. On aperture 3 it answers
/// them all.
#[test]
fn a_box_that_no_extent_can_meet_is_not_put_to_aperture_7() {
    let below = |x: f64, units: u32| (0..units).fold(x, |x, _| x.next_down());
    let kept = [
        [90.0, 0.0, 90.0, 10.0],
        [-90.0, -30.0, -90.0, 40.0],
        [90.0, -180.0, 90.0, 180.0],
        [10.0, 180.0, 20.0, 180.0],
        [10.0, -180.0, 20.0, -180.0],
        [10.0, 180.0, 10.0, 180.0],
        // Over the antimeridian from +180 to -180 degrees: both pieces have no width.
        [10.0, 180.0, 20.0, -180.0],
        // Two units in the last place within the pole and the antimeridian: within the margin.
        [below(90.0, 2), 0.0, below(90.0, 2), 10.0],
        [10.0, below(180.0, 2), 20.0, below(180.0, 2)],
    ];
    let asked = [
        // A hundred units in the last place within: beyond the margin.
        [below(90.0, 100), 0.0, below(90.0, 100), 10.0],
        [10.0, below(180.0, 100), 20.0, below(180.0, 100)],
        // Boxes with a width that reach the antimeridian, from either side and over it.
        [10.0, 179.0, 20.0, 180.0],
        [10.0, 180.0, 20.0, -179.0],
        [10.0, 179.0, 20.0, -180.0],
        // A box with a height that reaches the pole; no width elsewhere; a point elsewhere.
        [89.9, 0.0, 90.0, 10.0],
        [10.0, 11.2, 80.0, 11.2],
        [0.0, 0.0, 0.0, 0.0],
    ];
    for grid in [IGEO7, IVEA7H, RTEA7H] {
        for level in [0, 2, 5, 14] {
            for b in &kept {
                assert!(!engine_survives(grid, level, b), "{grid} {level} {b:?}");
            }
        }
        for b in &asked {
            assert!(engine_survives(grid, 4, b), "{grid} {b:?}");
            assert!(
                !list_zones(grid, 4, Some(radians(b))).is_empty(),
                "{grid} {b:?}"
            );
        }
    }
    // Aperture 3 answers them all. A point on a pole is another matter: no grid is asked it.
    for b in kept.iter().chain(&asked) {
        assert!(engine_survives(ISEA3H, 4, b), "{b:?}");
        list_zones(ISEA3H, 4, Some(radians(b)));
    }
    assert!(!engine_survives(ISEA3H, 4, &[90.0, 0.0, 90.0, 0.0]));
}

/// The seam is the engine's own: its four vertices are where the engine's cost was measured,
/// on each grid of aperture 3.
#[test]
fn the_seam_runs_where_it_was_measured() {
    let measured = [
        (58.397, 11.2),
        (58.397, -168.8),
        (0.0, -137.083),
        (-58.397, -168.8),
    ];
    for grid in [ISEA3H, IVEA3H, RTEA3H] {
        for ((lat, lon), (mlat, mlon)) in seam(grid).into_iter().zip(measured) {
            assert!(
                (lat - mlat).abs() < 1e-3 && (lon - mlon).abs() < 1e-3,
                "{grid}: ({lat}, {lon}) against ({mlat}, {mlon})"
            );
        }
    }
}

/// On aperture 3 a box across the seam of the engine's layout is kept from it beyond level 8,
/// and so is a box across one of the eight other cut edges far from the end at which the layout
/// is whole; a box at one of the eight vertices that are not of the seam is put to it at the
/// finest levels, and so is a box on the south pole as far as level 16, and the engine answers
/// them.
#[test]
fn the_cuts_of_the_engines_layout_on_aperture_3() {
    // The middles of the two edges of the seam that do not pass over the pole: a box of four
    // zone widths there is not answered at level 16.
    for mid in [(30.0574, -147.8712), (-30.0574, -147.8712)] {
        for grid in [ISEA3H, IVEA3H, RTEA3H] {
            assert!(!engine_survives(grid, 16, &about(grid, 16, mid, 2.0)));
            assert!(!engine_survives(grid, 9, &about(grid, 9, mid, 2.0)));
            assert!(engine_survives(grid, 8, &about(grid, 8, mid, 2.0)));
        }
    }
    // The third edge of the seam passes over the north pole.
    assert!(!engine_survives(ISEA3H, 16, &[89.99, -180.0, 90.0, 180.0]));
    assert!(!engine_survives(
        ISEA3H,
        16,
        &about(ISEA3H, 16, (75.0, 11.2), 2.0)
    ));
    assert!(!engine_survives(
        ISEA3H,
        16,
        &about(ISEA3H, 16, (75.0, -168.8), 2.0)
    ));
    // The four vertices of the seam.
    for v in seam(ISEA3H) {
        assert!(
            !engine_survives(ISEA3H, 20, &about(ISEA3H, 20, v, 2.0)),
            "{v:?}"
        );
    }
    // Sixty zone widths to the side of the seam the engine is asked, and answers.
    let beside = [30.0574, -147.8712 + 60.0 * zone_width_deg(ISEA3H, 16)];
    let b = about(ISEA3H, 16, (beside[0], beside[1]), 2.0);
    assert!(engine_survives(ISEA3H, 16, &b));
    assert!(!list_zones(ISEA3H, 16, Some(radians(&b))).is_empty());

    // The eight vertices that are not of the seam, at the finest levels: asked and answered.
    for grid in [ISEA3H, IVEA3H, RTEA3H] {
        let on_seam = seam(grid);
        let others: Vec<(f64, f64)> = roots(grid)
            .into_iter()
            .filter(|v| !on_seam.contains(v))
            .collect();
        assert_eq!(others.len(), 8, "{grid}");
        for v in others {
            for level in [12, 16, 20, 26, 33] {
                let b = about(grid, level, v, 2.0);
                assert!(engine_survives(grid, level, &b), "{grid} {level} {v:?}");
                assert!(
                    !list_zones(grid, level, Some(radians(&b))).is_empty(),
                    "{grid} {level} {v:?}"
                );
            }
        }
    }

    // The south pole lies in the middle of a cut edge that is not of the seam: 3,240 zone widths
    // from its whole end at level 16, where the engine is asked and answers, and 29,000 at
    // level 20, where it is not asked.
    let south = |level: i32| {
        [
            -90.0,
            -180.0,
            -90.0 + 2.0 * zone_width_deg(ISEA3H, level),
            180.0,
        ]
    };
    for level in [9, 12, 16] {
        assert!(engine_survives(ISEA3H, level, &south(level)), "{level}");
        list_zones(ISEA3H, level, Some(radians(&south(level))));
    }
    for level in [17, 20, 26, 33] {
        assert!(!engine_survives(ISEA3H, level, &south(level)), "{level}");
    }
    // Along such an edge, the one that runs south from (-58.397, 11.2) to the pole: at level 26
    // a box 1,000 zone widths from the whole end is asked and answered, one 100,000 is not asked.
    let (lat, lon) = roots(ISEA3H)[5];
    let at = |widths: f64| (lat - widths * zone_width_deg(ISEA3H, 26), lon);
    let near = about(ISEA3H, 26, at(1_000.0), 2.0);
    assert!(engine_survives(ISEA3H, 26, &near));
    assert!(!list_zones(ISEA3H, 26, Some(radians(&near))).is_empty());
    assert!(!engine_survives(
        ISEA3H,
        26,
        &about(ISEA3H, 26, at(100_000.0), 2.0)
    ));
}

/// Every kind of box on one grid of each aperture, at levels of every class: each box that the
/// engine may be asked is put to it, and it answers in the ascending order of its identifiers
/// without a duplicate. The box is handed to it in radians as the crate converts.
#[test]
fn every_kind_of_box_at_every_class_of_level_is_answered() {
    for (grid, levels) in [
        (ISEA3H, [1, 4, 8, 12, 16, 24, 33].as_slice()),
        (IGEO7, [0, 2, 5, 10, 14].as_slice()),
    ] {
        for &level in levels {
            let size = zone_width_deg(grid, level);
            let mut rng = Rng::new(1_000 + level as u64);
            let sites = sites(grid, level, &mut rng);
            let (mut asked, mut kept) = (0, 0);
            let mut slowest = (0.0, String::new());
            for kind in KINDS {
                for _ in 0..20 {
                    let b = make_box(kind, &mut rng, size, &sites);
                    if !engine_survives(grid, level, &b) {
                        kept += 1;
                        continue;
                    }
                    let t = std::time::Instant::now();
                    let z = list_zones(grid, level, Some(radians(&b)));
                    let took = t.elapsed().as_secs_f64();
                    if took > slowest.0 {
                        slowest = (took, format!("{kind:?} {b:?} ({} zones)", z.len()));
                    }
                    assert!(
                        z.windows(2).all(|w| w[0] < w[1]) || is_aperture7(grid),
                        "{grid} level {level} {kind:?} {b:?}"
                    );
                    asked += 1;
                }
            }
            eprintln!(
                "{grid} level {level}: {asked} asked, {kept} kept; slowest {:.3} s: {}",
                slowest.0, slowest.1
            );
            assert!(
                asked > 0,
                "{grid} level {level}: nothing asked ({kept} kept)"
            );
        }
    }
    assert_eq!(GRIDS.len(), 6);
}
