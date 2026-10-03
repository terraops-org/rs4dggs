//! An upper bound, by arithmetic, on the number of zones of a level whose extent meets a
//! bounding box: held here against the count itself, which is the length of the answer of
//! `zones_in_box`, on boxes of every size at every level, beside the poles, over the
//! antimeridian and about the vertices of the icosahedron, and on boxes of hundreds of
//! thousands of zones at high latitudes. None of it needs the engine.
use rs4dggs::indexings::{I3h, Z7};
use rs4dggs::projections::Isea;
use rs4dggs::topologies::{HexA3, HexA7};
use rs4dggs::{AnyGrid, Error, Extent, GeoPoint, Grid, GridConfig};

const APERTURE_3: [&str; 3] = ["ISEA3H", "IVEA3H", "RTEA3H"];
const APERTURE_7: [&str; 3] = ["IGEO7", "IVEA7H", "RTEA7H"];

/// What the estimate adds to the product of the areas, so that it is never nought.
const ALLOWANCE: u64 = 3;

/// The box about Lisbon, `[south, west, north, east]` in degrees, as every box here.
const LISBON: [f64; 4] = [38.0, -10.0, 39.0, -9.0];

/// The whole world.
const WORLD: [f64; 4] = [-90.0, -180.0, 90.0, 180.0];

fn grid(name: &str) -> AnyGrid {
    rs4dggs::get_grid(name).unwrap()
}

/// A box `[south, west, north, east]`, in degrees, as the extent a caller hands over.
fn degrees(b: &[f64; 4]) -> Extent {
    Extent {
        ll: GeoPoint {
            lat: b[0],
            lon: b[1],
        },
        ur: GeoPoint {
            lat: b[2],
            lon: b[3],
        },
    }
}

/// How many zones of `level` have an extent that meets the box: the answer, counted.
fn count(grid: AnyGrid, level: u8, b: &[f64; 4]) -> u64 {
    grid.zones_in_box(level, &degrees(b)).unwrap().count() as u64
}

fn estimate(grid: AnyGrid, level: u8, b: &[f64; 4]) -> u64 {
    grid.estimate_zones_in_box(level, &degrees(b)).unwrap()
}

/// The width of a zone of `level`, in degrees: the square root of the sphere's area over the
/// count of zones.
fn zone_width(grid: AnyGrid, level: u8) -> f64 {
    let zones = grid.count_zones(level).unwrap() as f64;
    (4.0 * std::f64::consts::PI / zones).sqrt().to_degrees()
}

/// The box about the point `(lat, lon)` that reaches `half` degrees of arc from it on every
/// side, cut at the poles, brought within -180 to 180 degrees of longitude (so that it may run
/// over the antimeridian), and the whole circle where it is as wide as that. With `half`
/// nought it is the point itself.
fn about(lat: f64, lon: f64, half: f64) -> [f64; 4] {
    let (south, north) = ((lat - half).max(-90.0), (lat + half).min(90.0));
    let cos = lat.to_radians().cos();
    if half >= 180.0 * cos {
        return [south, -180.0, north, 180.0];
    }
    let within = |x: f64| {
        if x < -180.0 {
            x + 360.0
        } else if x > 180.0 {
            x - 360.0
        } else {
            x
        }
    };
    [
        south,
        within(lon - half / cos),
        north,
        within(lon + half / cos),
    ]
}

/// The estimate over the count, by the size of the answer: 1 to 9 zones, 10 to 99, 100 to 999,
/// a thousand and more.
#[derive(Debug)]
struct Ratios {
    /// The least ratio met, over every answer of at least one zone.
    least: f64,
    /// The greatest ratio met in each class of size.
    greatest: [f64; 4],
    /// How many answers of each class were counted.
    answers: [usize; 4],
    /// How many boxes held no zone.
    empty: usize,
}

impl Ratios {
    fn new() -> Self {
        Ratios {
            least: f64::INFINITY,
            greatest: [0.0; 4],
            answers: [0; 4],
            empty: 0,
        }
    }

    /// Holds the estimate of a box against its count, and records their ratio.
    fn hold(&mut self, grid: AnyGrid, level: u8, b: &[f64; 4]) {
        let (counted, estimated) = (count(grid, level, b), estimate(grid, level, b));
        assert!(
            estimated >= counted,
            "{} level {level}, the box {b:?}: {counted} zones, estimated at {estimated}",
            grid.name()
        );
        assert!(estimated >= ALLOWANCE);
        assert!(estimated <= grid.count_zones(level).unwrap() + ALLOWANCE);
        if counted == 0 {
            self.empty += 1;
            return;
        }
        let class = match counted {
            1..=9 => 0,
            10..=99 => 1,
            100..=999 => 2,
            _ => 3,
        };
        let ratio = estimated as f64 / counted as f64;
        self.least = self.least.min(ratio);
        self.greatest[class] = self.greatest[class].max(ratio);
        self.answers[class] += 1;
    }
}

/// The half-sides of the boxes made about each place, in widths of a zone of the level: a
/// point, a box smaller than a zone, and boxes of some nine, some hundred and fifty and some
/// thirteen hundred zones.
const HALF_SIDES: [f64; 5] = [0.0, 0.3, 1.5, 6.0, 18.0];

/// The estimate against the count at every level of a grid, on boxes of each of
/// [`HALF_SIDES`] about each of nine places: Lisbon; the two poles; two points of the
/// antimeridian; two vertices of the icosahedron, where a pentagon lies; a point of a side of
/// a rhombus; and a point of the equator.
fn never_under(name: &str) -> Ratios {
    let g = grid(name);
    let mut vertices = g.zones(0).unwrap().map(|z| g.centroid(z));
    let (first, eighth) = (vertices.next().unwrap(), vertices.nth(6).unwrap());
    let places = [
        (38.7223, -9.1393),
        (90.0, 0.0),
        (-90.0, 0.0),
        (0.0, 180.0),
        (-45.0, -179.999),
        (first.lat, first.lon),
        (eighth.lat, eighth.lon),
        (70.0, 11.2),
        (0.0, 0.0),
    ];
    let mut ratios = Ratios::new();
    for level in 0..=g.max_box_level() {
        let width = zone_width(g, level);
        for (lat, lon) in places {
            for half in HALF_SIDES {
                ratios.hold(g, level, &about(lat, lon, half * width));
            }
        }
    }
    eprintln!("{name}: {ratios:?}");
    ratios
}

/// What every grid must show: no estimate under its count, answers of every class of size in
/// number, and no estimate farther above its count than the greatest that was measured, with
/// a little room. On these boxes, which are as high as they are wide, an answer of under ten
/// zones is estimated at up to twenty times its size, and one of a thousand zones or more at
/// no more than a fifth above it.
fn is_within_what_was_measured(ratios: &Ratios, levels: usize) {
    assert!(ratios.least >= 1.0, "{ratios:?}");
    for (class, floor) in [10, 5, 5, 4].into_iter().enumerate() {
        assert!(
            ratios.answers[class] >= floor * levels,
            "class {class} of {ratios:?}"
        );
    }
    for (class, ceiling) in [20.0, 3.5, 1.6, 1.2].into_iter().enumerate() {
        assert!(
            ratios.greatest[class] <= ceiling,
            "class {class} of {ratios:?}"
        );
    }
}

#[test]
fn the_estimate_is_never_under_the_count_on_isea3h() {
    is_within_what_was_measured(&never_under("ISEA3H"), 34);
}

#[test]
fn the_estimate_is_never_under_the_count_on_ivea3h() {
    is_within_what_was_measured(&never_under("IVEA3H"), 34);
}

#[test]
fn the_estimate_is_never_under_the_count_on_rtea3h() {
    is_within_what_was_measured(&never_under("RTEA3H"), 34);
}

#[test]
fn the_estimate_is_never_under_the_count_on_igeo7() {
    is_within_what_was_measured(&never_under("IGEO7"), 15);
}

#[test]
fn the_estimate_is_never_under_the_count_on_ivea7h() {
    is_within_what_was_measured(&never_under("IVEA7H"), 15);
}

#[test]
fn the_estimate_is_never_under_the_count_on_rtea7h() {
    is_within_what_was_measured(&never_under("RTEA7H"), 15);
}

/// The boxes whose zones are known by name elsewhere: the box about Lisbon, the small caps
/// about the poles, the box of 355 degrees over the antimeridian, a segment of a meridian and
/// one of a parallel, at every level at which the answer is small enough to count.
#[test]
fn the_estimate_is_never_under_the_count_of_the_boxes_whose_zones_are_known() {
    let boxes = [
        LISBON,
        [89.5, -180.0, 90.0, 180.0],
        [-90.0, -180.0, -89.5, 180.0],
        [89.9, 0.0, 90.0, 90.0],
        [89.997928885, 11.131869680, 90.0, 108.081826856],
        [30.0, 10.0, 40.0, 5.0],
        [38.0, -9.0, 39.0, -9.0],
        [38.5, -10.0, 38.5, -9.0],
    ];
    let mut counted = 0;
    for name in APERTURE_3.into_iter().chain(APERTURE_7) {
        let g = grid(name);
        let mut ratios = Ratios::new();
        for level in 0..=g.max_box_level() {
            for b in &boxes {
                // An answer of more than some five thousand zones is left to the tests of
                // the large boxes: the estimate itself says which.
                if estimate(g, level, b) <= 5_000 {
                    ratios.hold(g, level, b);
                    counted += 1;
                }
            }
        }
        assert!(ratios.least >= 1.0, "{name}: {ratios:?}");
    }
    assert!(counted >= 650, "{counted} boxes counted");
}

/// The estimate of a box on the twin of a grid that takes latitudes as they come, as
/// latitudes of the sphere, and converts none.
fn estimate_without_conversion(name: &str, level: u8, b: &[f64; 4]) -> u64 {
    let mut plain = GridConfig::default();
    plain.authalic = false;
    let b = degrees(b);
    let estimated = match name {
        "IGEO7" => Grid::<Isea, HexA7, Z7>::new(plain, "IGEO7")
            .unwrap()
            .estimate_zones_in_box(level, &b),
        "ISEA3H" => Grid::<Isea, HexA3, I3h>::new(plain, "ISEA3H")
            .unwrap()
            .estimate_zones_in_box(level, &b),
        other => panic!("no twin of {other} is built here"),
    };
    estimated.unwrap()
}

/// A large box at a high latitude, where the estimate must take the latitudes of the box on
/// the sphere of equal area: a cap above 70 degrees is some 0.8 per cent larger there than
/// its geodetic latitudes say, which from some two hundred thousand zones on exceeds what the
/// estimate adds about the rim.
///
/// `zones` is the count, which is asserted. `within` is how far above the count the estimate
/// may lie, as a factor: on boxes of this size the rim weighs little, and the estimate is a
/// few parts in a hundred, or in a thousand, above the count. `gain` is the part of the
/// estimate that the conversion of the latitudes brings, measured against the twin of the
/// grid that converts none, and held to two parts in ten thousand: without the conversion
/// the estimate of the two caps of hundreds of thousands of zones is under their count, and
/// that of the bands and of the smaller caps is not, the rim being enough there.
fn a_large_box_is_not_under_estimated(
    name: &str,
    level: u8,
    b: &[f64; 4],
    zones: u64,
    within: f64,
    gain: f64,
) {
    let g = grid(name);
    let (counted, estimated) = (count(g, level, b), estimate(g, level, b));
    let unconverted = estimate_without_conversion(name, level, b);
    let gained = estimated as f64 / unconverted as f64 - 1.0;
    eprintln!(
        "{name} level {level}, {b:?}: {counted} zones, estimated at {estimated}, {:.5} of it; \
         {unconverted} without the conversion, which brings {gained:.5}",
        estimated as f64 / counted as f64
    );
    assert_eq!(counted, zones, "{name} level {level}, {b:?}");
    assert!(
        estimated >= counted,
        "{name} level {level}, {b:?}: {counted} zones, estimated at {estimated}"
    );
    assert!(
        estimated as f64 <= within * counted as f64,
        "{name} level {level}, {b:?}: {counted} zones, estimated at {estimated}"
    );
    assert!(
        (gained - gain).abs() < 2e-4,
        "{name} level {level}, {b:?}: the conversion brings {gained}"
    );
}

const NORTHERN_CAP: [f64; 4] = [70.0, -180.0, 90.0, 180.0];
const SOUTHERN_CAP: [f64; 4] = [-90.0, -180.0, -70.0, 180.0];
const HIGH_BAND: [f64; 4] = [60.0, -180.0, 80.0, 180.0];
const MIDDLE_BAND: [f64; 4] = [40.0, -180.0, 50.0, 180.0];

#[test]
fn the_cap_above_70_degrees_north_at_level_7_of_igeo7_is_not_under_estimated() {
    a_large_box_is_not_under_estimated("IGEO7", 7, &NORTHERN_CAP, 251_422, 1.009, 0.0081);
}

// The band from 60 to 80 degrees at level 7 of IGEO7 holds 494,686 zones and is estimated at
// 498,586, the conversion of the latitudes bringing 0.0070 of that. It is not run here: to
// count it takes longer than the two caps together, and it holds nothing that the same band
// at level 6, below, does not. As a call:
// `a_large_box_is_not_under_estimated("IGEO7", 7, &HIGH_BAND, 494_686, 1.009, 0.0070)`.

#[test]
fn the_cap_above_70_degrees_north_at_level_14_of_isea3h_is_not_under_estimated() {
    a_large_box_is_not_under_estimated("ISEA3H", 14, &NORTHERN_CAP, 1_456_540, 1.004, 0.0082);
}

#[test]
fn the_caps_and_the_bands_at_level_6_of_igeo7_are_not_under_estimated() {
    a_large_box_is_not_under_estimated("IGEO7", 6, &NORTHERN_CAP, 36_180, 1.022, 0.0080);
    a_large_box_is_not_under_estimated("IGEO7", 6, &SOUTHERN_CAP, 36_180, 1.022, 0.0080);
    a_large_box_is_not_under_estimated("IGEO7", 6, &HIGH_BAND, 71_178, 1.022, 0.0070);
    a_large_box_is_not_under_estimated("IGEO7", 6, &MIDDLE_BAND, 74_278, 1.044, 0.0022);
}

/// The conversion of the latitudes to the sphere of equal area, seen without counting a zone:
/// on a grid that takes latitudes as they come, the same cap is estimated at fewer zones, by
/// the 0.8 per cent that the cap gains on the sphere of equal area.
#[test]
fn the_latitudes_of_the_box_are_taken_on_the_sphere_of_equal_area() {
    let mut plain = GridConfig::default();
    plain.authalic = false;
    let plain: Grid<Isea, HexA7, Z7> = Grid::new(plain, "IGEO7").unwrap();
    let cap = degrees(&NORTHERN_CAP);
    let on_the_sphere = plain.estimate_zones_in_box(7, &cap).unwrap() as f64;
    let converted = rs4dggs::igeo7().estimate_zones_in_box(7, &cap).unwrap() as f64;
    let gain = converted / on_the_sphere - 1.0;
    assert!((0.0075..0.0090).contains(&gain), "{gain}");
    // The whole world is the whole world on either.
    let world = degrees(&WORLD);
    assert_eq!(
        plain.estimate_zones_in_box(7, &world).unwrap(),
        rs4dggs::igeo7().estimate_zones_in_box(7, &world).unwrap()
    );
}

/// Beside a pole at the finest levels the sines of the two latitudes of a small cap are the
/// same number, or differ by one unit in the last place, and the cap holds zones all the
/// same: the estimate is not the allowance alone, and is not under the count.
#[test]
fn a_small_cap_about_a_pole_at_the_finest_level_is_not_estimated_at_nought() {
    for name in APERTURE_3.into_iter().chain(APERTURE_7) {
        let g = grid(name);
        let level = g.max_box_level();
        let width = zone_width(g, level);
        for (south, north) in [(90.0 - 6.0 * width, 90.0), (-90.0, -90.0 + 6.0 * width)] {
            let cap = [south, -180.0, north, 180.0];
            let (counted, estimated) = (count(g, level, &cap), estimate(g, level, &cap));
            assert!(counted >= 100, "{name} level {level}, {cap:?}: {counted}");
            assert!(
                estimated >= counted,
                "{name} level {level}, {cap:?}: {counted} zones, estimated at {estimated}"
            );
            assert!(estimated <= 3 * counted);
        }
    }
}

/// Where the estimate stands nearest its count: a narrow wedge upon a pole at the finest level
/// of the aperture-3 grids, a millionth of a zone width high and one degree wide. Ten zones
/// have an extent that reaches the north pole there, from centroids up to 1.7 zone widths
/// away, and eight of them meet the wedge, which lies across the meridian of 11.2 degrees
/// east; the grown box, which near a pole is grown less in longitude than the zones reach by
/// arc, has the area of some five zones. The estimate is 9, and holds by its allowance alone,
/// with one zone to spare. On the south pole the same wedge has 6 zones.
#[test]
fn a_narrow_wedge_upon_a_pole_at_the_finest_level_is_estimated_one_zone_above_its_count() {
    for name in APERTURE_3 {
        let g = grid(name);
        let height = 1e-6 * zone_width(g, 33);
        let north = [90.0 - height, 11.0, 90.0, 12.0];
        let (counted, estimated) = (count(g, 33, &north), estimate(g, 33, &north));
        assert_eq!((counted, estimated), (8, 9), "{name}, {north:?}");
        let south = [-90.0, 11.0, -90.0 + height, 12.0];
        let (counted, estimated) = (count(g, 33, &south), estimate(g, 33, &south));
        assert_eq!(counted, 6, "{name}, {south:?}");
        assert!(
            estimated >= counted,
            "{name}, {south:?}: {counted} zones, estimated at {estimated}"
        );
    }
}

/// The whole world is estimated at the count of the level, and at most the allowance above
/// it, at every level: at the finest levels of aperture 3 the count is beyond what a `f64`
/// holds exactly, and the estimate is still not under it.
#[test]
fn the_whole_world_is_estimated_at_the_count_of_the_level() {
    for name in APERTURE_3.into_iter().chain(APERTURE_7) {
        let g = grid(name);
        for level in 0..=g.max_box_level() {
            let zones = g.count_zones(level).unwrap();
            for world in [WORLD, [-90.0, 0.0, 90.0, -1e-300]] {
                let estimated = estimate(g, level, &world);
                assert!(
                    (zones..=zones + ALLOWANCE).contains(&estimated),
                    "{name} level {level}, {world:?}: {zones} zones, estimated at {estimated}"
                );
            }
        }
        // At a level whose zones can be counted, the count of the box is that of the level.
        assert_eq!(count(g, 3, &WORLD), g.count_zones(3).unwrap());
    }
}

/// A point is estimated at a few zones, and never at fewer than hold it: one, two or three.
#[test]
fn a_point_is_estimated_at_a_few_zones() {
    for name in APERTURE_3.into_iter().chain(APERTURE_7) {
        let g = grid(name);
        for level in 0..=g.max_box_level() {
            for (lat, lon) in [(38.7223, -9.1393), (-33.9, 151.2), (12.3, 179.5)] {
                let point = [lat, lon, lat, lon];
                let (counted, estimated) = (count(g, level, &point), estimate(g, level, &point));
                assert!(
                    (1..=3).contains(&counted),
                    "{name} level {level}: {counted}"
                );
                assert!(
                    (counted..=20).contains(&estimated),
                    "{name} level {level}, {point:?}: {counted} zones, estimated at {estimated}"
                );
            }
        }
    }
}

/// The estimate refuses what `zones_in_box` refuses, with the same error, the level before the
/// box, and computes nothing first.
#[test]
fn the_estimate_refuses_what_the_zones_of_a_box_refuse() {
    let refused = [
        ("ISEA3H", 34, LISBON),
        ("ISEA3H", 255, LISBON),
        ("IGEO7", 15, LISBON),
        ("IGEO7", 19, LISBON),
        ("IGEO7", 20, LISBON),
        ("ISEA3H", 5, [f64::NAN, -10.0, 39.0, -9.0]),
        ("ISEA3H", 5, [38.0, f64::INFINITY, 39.0, -9.0]),
        ("ISEA3H", 5, [38.0, -10.0, f64::NEG_INFINITY, -9.0]),
        ("ISEA3H", 5, [38.0, -10.0, 39.0, f64::NAN]),
        ("IGEO7", 5, [-90.5, -10.0, 39.0, -9.0]),
        ("IGEO7", 5, [38.0, -10.0, 90.000001, -9.0]),
        ("IGEO7", 5, [38.0, -180.5, 39.0, -9.0]),
        ("IGEO7", 5, [38.0, -10.0, 39.0, 181.0]),
        ("IGEO7", 5, [38.0, 190.0, 39.0, -9.0]),
        ("RTEA3H", 5, [39.0, -10.0, 38.0, -9.0]),
        // The level is refused before the box.
        ("RTEA3H", 34, [39.0, -10.0, 38.0, f64::NAN]),
        ("RTEA7H", 15, [f64::NAN, -10.0, 38.0, -9.0]),
    ];
    for (name, level, b) in refused {
        let g = grid(name);
        let by_the_estimate = g.estimate_zones_in_box(level, &degrees(&b)).unwrap_err();
        let by_the_zones = g.zones_in_box(level, &degrees(&b)).unwrap_err();
        assert_eq!(by_the_estimate, by_the_zones, "{name} level {level}, {b:?}");
    }
    assert_eq!(
        grid("IGEO7").estimate_zones_in_box(15, &degrees(&LISBON)),
        Err(Error::ResolutionOutOfRange { res: 15, max: 14 })
    );
    // The limits themselves are accepted.
    assert!(
        grid("IGEO7")
            .estimate_zones_in_box(14, &degrees(&WORLD))
            .is_ok()
    );
    assert!(
        grid("ISEA3H")
            .estimate_zones_in_box(33, &degrees(&WORLD))
            .is_ok()
    );
}

/// The grid taken by its type and the handle on it give one estimate.
#[test]
fn the_grid_and_its_handle_give_one_estimate() {
    let iberia = degrees(&[36.0, -10.0, 44.0, 3.0]);
    for level in [0, 7, 14] {
        assert_eq!(
            rs4dggs::isea3h().estimate_zones_in_box(level, &iberia),
            grid("ISEA3H").estimate_zones_in_box(level, &iberia)
        );
        assert_eq!(
            rs4dggs::rtea7h().estimate_zones_in_box(level, &iberia),
            grid("RTEA7H").estimate_zones_in_box(level, &iberia)
        );
    }
    // Sixty zones of level 7 of ISEA3H meet the box.
    assert_eq!(count(grid("ISEA3H"), 7, &[36.0, -10.0, 44.0, 3.0]), 60);
    let estimated = grid("ISEA3H").estimate_zones_in_box(7, &iberia).unwrap();
    assert!((60..=120).contains(&estimated), "{estimated}");
}
