//! The documents the encodings write, compared byte for byte with the files under
//! `tests/goldens/`, each of which was read by eye once when it was made.
//!
//! When a document differs from its file, the test names the first byte that differs and leaves
//! what was written under the tests' temporary directory, to be looked at and, if it is right,
//! copied over the file by hand. A missing file fails the test.

use std::path::Path;

use rs4dggs::{AnyGrid, ZoneId};
use rs4dggs_ogc::{
    Coordinates, Link, LinkTemplate, ZoneGeometry, ZoneListGeoJson, ZoneListJson, rel,
    write_zone_feature,
};

/// Holds `written` to the file `tests/goldens/<name>`.
fn assert_golden(name: &str, written: &[u8]) {
    assert_golden_keeping(name, name, written);
}

/// Holds `written` to the file `tests/goldens/<name>`, and keeps it as `kept` when it differs.
fn assert_golden_keeping(name: &str, kept: &str, written: &[u8]) {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/goldens")
        .join(name);
    let golden = std::fs::read(&path)
        .unwrap_or_else(|e| panic!("the golden {} cannot be read: {e}", path.display()));
    if golden == written {
        return;
    }
    let kept = Path::new(env!("CARGO_TARGET_TMPDIR")).join(kept);
    std::fs::write(&kept, written).expect("what was written can be kept");
    let first = golden
        .iter()
        .zip(written)
        .position(|(a, b)| a != b)
        .unwrap_or(golden.len().min(written.len()));
    panic!(
        "{name}: the document differs from its golden at byte {first} (the golden has {} bytes, \
         the document {}); what was written is in {}",
        golden.len(),
        written.len(),
        kept.display()
    );
}

fn isea3h() -> AnyGrid {
    rs4dggs::get_grid("ISEA3H").unwrap()
}

fn ids(grid: AnyGrid, texts: &[&str]) -> Vec<ZoneId> {
    texts
        .iter()
        .map(|t| grid.zone_from_text(t).unwrap())
        .collect()
}

/// The sum of the zones' areas, in the order given.
fn area(grid: AnyGrid, zones: &[ZoneId]) -> f64 {
    zones.iter().map(|&z| grid.area(z).unwrap()).sum()
}

/// The 32 zones of level 1 of ISEA3H, in the order of the standard's example of a zone query
/// (its Annex C), which is also DGGAL's.
const ISEA3H_LEVEL_1: [&str; 32] = [
    "A0-0-B", "A0-0-C", "A0-0-D", "A1-0-B", "A1-0-C", "A1-0-D", "A2-0-B", "A2-0-C", "A2-0-D",
    "A3-0-B", "A3-0-C", "A3-0-D", "A4-0-B", "A4-0-C", "A4-0-D", "A5-0-B", "A5-0-C", "A5-0-D",
    "A6-0-B", "A6-0-C", "A6-0-D", "A7-0-B", "A7-0-C", "A7-0-D", "A8-0-B", "A8-0-C", "A8-0-D",
    "A9-0-B", "A9-0-C", "A9-0-D", "AA-0-B", "AB-0-B",
];

/// The links and the template the golden of level 1 carries.
const LEVEL_1_LINKS: [Link<'static>; 3] = [
    Link {
        rel: rel::DGGRS,
        href: "/collections/gebco/dggs/ISEA3H",
        title: Some("ISEA3H DGGS for GEBCO"),
        media_type: None,
    },
    Link {
        rel: rel::DGGRS_DEFINITION,
        href: "/dggrs/ISEA3H",
        title: Some("ISEA3H DGGRS definition"),
        media_type: Some("application/json"),
    },
    Link {
        rel: rel::GEODATA,
        href: "/collections/gebco",
        title: None,
        media_type: None,
    },
];
const LEVEL_1_TEMPLATES: [LinkTemplate<'static>; 1] = [LinkTemplate {
    rel: rel::DGGRS_ZONE_DATA,
    uri_template: "/collections/gebco/dggs/ISEA3H/zones/{zoneId}/data",
    title: Some("ISEA3H data for GEBCO"),
}];

fn list(grid: AnyGrid, zones: &[ZoneId], area: f64, links: &[Link], t: &[LinkTemplate]) -> Vec<u8> {
    let mut out = Vec::new();
    let mut list = ZoneListJson::new(&mut out, grid);
    for &z in zones {
        list.zone(z).unwrap();
    }
    assert_eq!(
        list.finish(Some(area), links, t).unwrap(),
        zones.len() as u64
    );
    out
}

#[test]
fn two_zones_of_isea3h_with_their_area_and_the_two_links_required() {
    let grid = isea3h();
    let zones = ids(grid, &["A6-0-C", "AA-0-B"]);
    let area = area(grid, &zones);
    assert_eq!(area, 31170676883138.758);
    let links = [
        Link {
            rel: rel::DGGRS,
            href: "/dggs/ISEA3H",
            title: None,
            media_type: None,
        },
        Link {
            rel: rel::DGGRS_DEFINITION,
            href: "/dggrs/ISEA3H",
            title: None,
            media_type: None,
        },
    ];
    assert_golden(
        "zones-ISEA3H-two.json",
        &list(grid, &zones, area, &links, &[]),
    );
}

#[test]
fn the_level_1_of_isea3h_with_its_area_links_and_template() {
    let grid = isea3h();
    let zones = ids(grid, &ISEA3H_LEVEL_1);
    let area = area(grid, &zones);
    // The standard prints 510065621724088.5; the last digits depend on the order of the sum.
    assert!((area - 510065621724088.5).abs() <= 1.0, "{area}");
    assert_golden(
        "zones-ISEA3H-level1.json",
        &list(grid, &zones, area, &LEVEL_1_LINKS, &LEVEL_1_TEMPLATES),
    );
}

#[test]
fn the_librarys_zones_of_level_1_give_the_golden_of_level_1() {
    let grid = isea3h();
    let zones: Vec<ZoneId> = grid.zones(1).unwrap().collect();
    let area = area(grid, &zones);
    assert_golden(
        "zones-ISEA3H-level1.json",
        &list(grid, &zones, area, &LEVEL_1_LINKS, &LEVEL_1_TEMPLATES),
    );
}

#[test]
fn a_compacted_page_of_a_box_carries_the_area_of_the_page_before_compaction() {
    let grid = isea3h();
    let bbox = rs4dggs::Extent {
        ll: rs4dggs::GeoPoint {
            lat: 30.0,
            lon: -10.0,
        },
        ur: rs4dggs::GeoPoint {
            lat: 50.0,
            lon: 10.0,
        },
    };
    let page: Vec<ZoneId> = grid.zones_in_box(8, &bbox).unwrap().take(200).collect();
    assert_eq!(page.len(), 200);
    let compacted = grid.compact_zones(&page).unwrap();
    assert!(compacted.len() < page.len(), "{}", compacted.len());
    let page_area = area(grid, &page);
    // On aperture 3 a zone has about the area of three of its children, so that the seven that
    // compact into it sum to more than it: the areas differ.
    assert_ne!(page_area, area(grid, &compacted));
    let texts: Vec<String> = compacted
        .iter()
        .map(|&z| format!("\"{}\"", grid.text_id(z)))
        .collect();
    let expected = format!(
        "{{\"zones\":[{}],\"returnedAreaMetersSquare\":{page_area},\"links\":[]}}",
        texts.join(",")
    );
    assert!(page_area > 1e12 && page_area < 1e17);
    assert_eq!(
        String::from_utf8(list(grid, &compacted, page_area, &[], &[])).unwrap(),
        expected
    );
}

const REGION: ZoneGeometry = ZoneGeometry::Region { edge_refinement: 0 };

/// The coordinates in their shortest exact form, so that the documents hold the library's own
/// doubles.
const SHORTEST: Coordinates = Coordinates { decimals: None };

/// The zones of `grid` named by `texts` as a GeoJSON list, the features numbered from 1, with
/// the coordinates in their shortest exact form.
fn geojson(grid: AnyGrid, texts: &[&str], geometry: ZoneGeometry) -> Vec<u8> {
    let mut out = Vec::new();
    let mut list = ZoneListGeoJson::new(&mut out, grid, geometry, SHORTEST, 1);
    for t in texts {
        list.zone(grid.zone_from_text(t).unwrap()).unwrap();
    }
    let written = list.finish().unwrap();
    assert_eq!(written.features, texts.len() as u64);
    assert_eq!(written.without_geometry, 0);
    out
}

/// The text of `doc` from the first `"geometry":` to the `,"properties"` after it.
fn first_geometry(doc: &str) -> &str {
    let start = doc.find(r#""geometry":"#).unwrap() + r#""geometry":"#.len();
    &doc[start..start + doc[start..].find(r#","properties""#).unwrap()]
}

#[test]
fn the_zone_list_of_e6_317_a_as_geojson() {
    let doc = geojson(isea3h(), &["E6-317-A"], REGION);
    assert_golden("zones-ISEA3H-E6-317-A.geojson", &doc);
    let doc = String::from_utf8(doc).unwrap();
    let geometry = first_geometry(&doc);
    let first = "[34.254939992054695,45.71280373883526]";
    assert!(
        geometry.starts_with(&format!(
            r#"{{"type":"Polygon","coordinates":[[{first},[34.230483446088726,45.65500194237048],"#
        )),
        "{geometry}"
    );
    assert!(geometry.ends_with(&format!(",{first}]]}}")), "{geometry}");
    assert_eq!(geometry.matches("],[").count() + 1, 49, "positions");
}

#[test]
fn the_zone_list_of_e6_317_a_as_centroids() {
    let doc = geojson(isea3h(), &["E6-317-A"], ZoneGeometry::Centroid);
    assert_golden("zones-ISEA3H-E6-317-A.centroid.geojson", &doc);
    assert_eq!(
        first_geometry(&String::from_utf8(doc).unwrap()),
        r#"{"type":"Point","coordinates":[34.78016915103061,45.42937741847805]}"#
    );
}

#[test]
fn a_missing_golden_fails_and_a_differing_document_fails_and_is_kept() {
    let missing = std::panic::catch_unwind(|| assert_golden("no-such-golden.json", b"{}"));
    assert!(missing.is_err());
    // Kept under a name of its own, so that it never stands where a real failure leaves its
    // document.
    let kept = Path::new(env!("CARGO_TARGET_TMPDIR")).join("a-differing-document.json");
    let _ = std::fs::remove_file(&kept);
    let differing = std::panic::catch_unwind(|| {
        assert_golden_keeping(
            "zones-ISEA3H-two.json",
            "a-differing-document.json",
            b"{\"zones\":[]}",
        )
    });
    assert!(differing.is_err());
    assert_eq!(std::fs::read(&kept).unwrap(), b"{\"zones\":[]}");
}

// Zone data.

use rs4dggs_ogc::{Dggrs, Property, Values, ZoneData, write_dggs_json};

#[path = "common/data_cases.rs"]
mod data_cases;

#[test]
fn zone_data_as_dggs_json_is_each_golden() {
    let mut seen = 0;
    data_cases::each(|name, data| {
        let mut out = Vec::new();
        write_dggs_json(&mut out, data).unwrap();
        assert_golden(name, &out);
        seen += 1;
    });
    assert_eq!(seen, 8);
}

use rs4dggs_ogc::write_dggs_ubjson;

/// The data of `C2-23-C` at depth 1, the values 0 to 6 as `f64`, the third missing if asked.
fn f64_band(missing_third: bool) -> Vec<u8> {
    let dggrs = Dggrs::from_id("ISEA3H").unwrap();
    let mut band: Vec<Option<f64>> = (0..7).map(|i| Some(f64::from(i))).collect();
    if missing_third {
        band[2] = None;
    }
    let mut out = Vec::new();
    write_dggs_ubjson(
        &mut out,
        &ZoneData {
            dggrs,
            zone: dggrs.grid().zone_from_text("C2-23-C").unwrap(),
            depths: &[1],
            properties: &[Property {
                name: "i",
                depths: &[Values::F64(&band)],
            }],
        },
    )
    .unwrap();
    out
}

#[test]
fn zone_data_as_dggs_ubjson_is_typed_when_nothing_is_missing() {
    let out = f64_band(false);
    assert_eq!(out.len(), 216);
    // `{`, the key `dggrs`, a string of 48 bytes.
    assert!(out.starts_with(&[
        0x7b, 0x69, 0x05, b'd', b'g', b'g', b'r', b's', 0x53, 0x69, 0x30
    ]));
    // `depths`: a counted array of one `U` 1.
    assert!(out.windows(12).any(|w| w == b"depths[#i\x01U\x01"));
    // `data`: float64, typed, seven of them; the second is 1.
    let at = out.windows(6).position(|w| w == b"[$D#i\x07").unwrap();
    assert_eq!(out[at + 14..at + 22], [0x3f, 0xf0, 0, 0, 0, 0, 0, 0]);
    assert_golden("data-ISEA3H-C2-23-C-depth1.f64.ubj", &out);
}

#[test]
fn zone_data_as_dggs_ubjson_is_counted_with_a_null_where_a_value_is_missing() {
    let out = f64_band(true);
    assert_eq!(out.len(), 213);
    let at = out.windows(4).position(|w| w == b"[#i\x07").unwrap();
    // `D` and eight bytes a value, and `Z` for the third.
    assert_eq!(out[at + 4], b'D');
    assert_eq!(out[at + 22], b'Z');
    assert_eq!(out[at + 23], b'D');
    assert_golden("data-ISEA3H-C2-23-C-depth1.f64.null-at-2.ubj", &out);
}

#[test]
fn the_zone_e6_317_a_as_a_feature() {
    let grid = isea3h();
    let zone: ZoneId = grid.zone_from_text("E6-317-A").unwrap();
    let mut out = Vec::new();
    assert!(write_zone_feature(&mut out, grid, zone, 1, REGION, SHORTEST).unwrap());
    assert_golden("zone-ISEA3H-E6-317-A.feature.geojson", &out);
}

/// `01012` lies across the antimeridian and is cut into three polygons; `005` has the north pole
/// as a vertex and is one polygon that touches it.
#[test]
fn the_zone_list_of_01012_and_005_on_isea7h_z7() {
    let grid = rs4dggs::get_grid("IGEO7").unwrap();
    let doc = geojson(grid, &["01012", "005"], REGION);
    assert_golden("zones-ISEA7H_Z7-01012-005.geojson", &doc);
    let doc = String::from_utf8(doc).unwrap();
    let (first, second) = doc.split_at(doc.find(r#"},{"type":"Feature","id":2,"#).unwrap());
    let multi = first_geometry(first);
    assert!(multi.starts_with(r#"{"type":"MultiPolygon","#), "{multi}");
    assert_eq!(multi.matches("]]],[[[").count() + 1, 3, "rings");
    let single = first_geometry(second);
    assert!(single.starts_with(r#"{"type":"Polygon","#), "{single}");
    assert!(single.contains(",90]"), "the pole: {single}");
}

// Zone data as GeoJSON.

use rs4dggs_ogc::{Written, write_zone_data_geojson};

/// The features numbered in `doc`, in order.
fn feature_ids(doc: &str) -> Vec<u64> {
    doc.split(r#"{"type":"Feature","id":"#)
        .skip(1)
        .map(|f| f[..f.find(',').unwrap()].parse().unwrap())
        .collect()
}

#[test]
fn zone_data_of_c2_23_c_as_geojson_with_centroids() {
    let dggrs = Dggrs::from_id("ISEA3H").unwrap();
    let elevation = [
        Some(-14.209265189621),
        Some(57.2747655762701),
        None,
        Some(120.5),
        Some(0.25),
        Some(-3.0),
        Some(1000.0),
    ];
    let data = ZoneData {
        dggrs,
        zone: dggrs.grid().zone_from_text("C2-23-C").unwrap(),
        depths: &[1],
        properties: &[Property {
            name: "Elevation",
            depths: &[Values::F64(&elevation)],
        }],
    };
    let mut out = Vec::new();
    let written =
        write_zone_data_geojson(&mut out, &data, ZoneGeometry::Centroid, SHORTEST).unwrap();
    assert_eq!(
        written,
        Written {
            features: 7,
            without_geometry: 0,
            null_zones_skipped: 0
        }
    );
    assert_golden("data-ISEA3H-C2-23-C-depth1.centroids.geojson", &out);
    let doc = String::from_utf8(out).unwrap();
    assert_eq!(doc.matches(r#""geometry":{"type":"Point","#).count(), 7);
    assert!(doc.contains(r#","properties":{"Elevation":-14.209265189621,"zoneID":"#));
    assert!(doc.contains(r#","properties":{"Elevation":null,"zoneID":"#));
}

/// The cases of zone data of several depths and at a broken seam, written without geometry.
#[test]
fn zone_data_as_geojson_without_geometry_is_each_golden() {
    let mut seen = 0;
    data_cases::each(|name, data| {
        let (features, skipped) = match name {
            "data-ISEA3H-C2-23-C-depths1-2.json" => (20, 0),
            "data-ISEA7H_Z7-00055353260226021-depth1.seam.json" => (11, 2),
            _ => return,
        };
        let mut out = Vec::new();
        let written =
            write_zone_data_geojson(&mut out, data, ZoneGeometry::None, Coordinates::default())
                .unwrap();
        assert_eq!(
            written,
            Written {
                features,
                without_geometry: features,
                null_zones_skipped: skipped
            },
            "{name}"
        );
        assert_golden(&name.replace(".json", ".none.geojson"), &out);
        let doc = String::from_utf8(out).unwrap();
        assert_eq!(feature_ids(&doc), (1..=features).collect::<Vec<_>>());
        seen += 1;
    });
    assert_eq!(seen, 2);
}
