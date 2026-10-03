# rs4dggs-ogc

The encodings of OGC API - Discrete Global Grid Systems, Part 1: Core (OGC 21-038r1), over the six
grids of the `rs4dggs` library. It writes the zone list as JSON, as GeoJSON and as 64-bit
integers; zone data as DGGS-JSON, DGGS-UBJSON and GeoJSON; a zone's geometry as a valid GeoJSON
geometry; and it holds the six grids as a service names them, with their URIs and their definition
documents.

Every encoding is a pure function that writes one document to a sink, anything that implements
`std::io::Write`: it reads no environment, opens no connection, knows nothing of HTTP and builds no
document in memory. A service passes the sink of its response; a test passes a `Vec<u8>`. The
crate depends on `rs4dggs` alone, which supplies every value written (the identifiers, the
sub-zone order, the refined ring, the centroid, the area); what this crate adds is the shape of
each document, the table of the grids, the two byte orders, and the making of a valid GeoJSON
geometry from a zone's ring. It has no `unsafe` code.

It does not write the documents that belong to a service: the zone information document, the
description of a grid, the list of grids, the landing page, the conformance declaration. It
offers the pieces they share (a zone's geometry, a zone as a feature, the relation types, the
grids' URIs and definitions). It reads nothing: there is no reader of any encoding, and no parser
of a request's parameters. It does not write DGGS-JSON-FG, JSON-FG, HTML, GeoTIFF or
CoverageJSON. It describes the standard's encodings; it claims no conformance to the standard,
which only the OGC can certify.

## What it writes

| Document | Media type | Writer |
|---|---|---|
| Zone list, JSON | `application/json` | `ZoneListJson` |
| Zone list, GeoJSON | `application/geo+json` | `ZoneListGeoJson` |
| Zone list, 64-bit integers | `application/x-binary` | `write_zone_list_uint64` |
| A zone's geometry, alone or as a feature | `application/geo+json` | `write_zone_geometry`, `write_zone_feature` |
| Zone data, DGGS-JSON | `application/json` | `write_dggs_json` |
| Zone data, DGGS-UBJSON | `application/ubjson` | `write_dggs_ubjson` |
| Zone data, GeoJSON | `application/geo+json` | `write_zone_data_geojson` |
| A grid's definition | `application/json` | `Dggrs::definition` |

The media types are constants in `media_type`, and the relation types of the standard's links,
in their `https` form, are constants in `rel`. Every writer returns a `Result`; no input of a
caller makes it panic.

## A list of zones in three forms

The writers of lists take the library's own grid, `rs4dggs::AnyGrid`, and the zones one at a time
or as a slice the caller has bounded. The area of the list, `returnedAreaMetersSquare`, is the
caller's to give (see "What a caller must know"). The features of a GeoJSON list are numbered from
the number the caller gives, so that a service that pages a list continues the count of the page
before; each carries the zone's text identifier as its property `zoneID`.

```rust
use rs4dggs_ogc::{
    Coordinates, Link, ZoneGeometry, ZoneListGeoJson, ZoneListJson, rel, write_zone_list_uint64,
};

let grid = rs4dggs::get_grid("ISEA3H")?;
let zones = [grid.zone_from_text("A6-0-C")?, grid.zone_from_text("AA-0-B")?];
let area: f64 = zones.iter().filter_map(|&zone| grid.area(zone)).sum();

// JSON: the identifiers as text, the area, and the two links the standard requires.
let links = [
    Link::new(rel::DGGRS, "/dggs/ISEA3H"),
    Link::new(rel::DGGRS_DEFINITION, "/dggrs/ISEA3H"),
];
let mut json = Vec::new();
let mut list = ZoneListJson::new(&mut json, grid);
for &zone in &zones {
    list.zone(zone)?;
}
assert_eq!(list.finish(Some(area), &links, &[])?, 2);
assert!(json.starts_with(br#"{"zones":["A6-0-C","AA-0-B"],"returnedAreaMetersSquare":31170676883138.758,"#));

// GeoJSON: one feature a zone, here with its centroid, numbered from 1; the coordinates with
// eight decimals, the default.
let mut geojson = Vec::new();
let centroid = ZoneGeometry::Centroid;
let mut list = ZoneListGeoJson::new(&mut geojson, grid, centroid, Coordinates::default(), 1);
for &zone in &zones {
    list.zone(zone)?;
}
let written = list.finish()?;
assert_eq!((written.features, written.without_geometry), (2, 0));
let geojson = String::from_utf8(geojson)?;
assert!(geojson.starts_with(r#"{"type":"FeatureCollection","features":[{"type":"Feature","id":1,"#));
assert!(geojson.contains(r#""properties":{"zoneID":"AA-0-B"}"#));

// 64-bit integers: the count, then each identifier, eight bytes each, little endian.
let mut binary = Vec::new();
write_zone_list_uint64(&mut binary, &zones)?;
assert_eq!(binary.len(), 8 + 2 * 8);
assert_eq!(binary[..8], 2u64.to_le_bytes());
assert_eq!(binary[8..16], zones[0].0.to_le_bytes());
# Ok::<(), Box<dyn std::error::Error>>(())
```

## Zone data in three forms

The writers of zone data take a `ZoneData`: the grid as the service names it (a `&'static Dggrs`,
whose URI the documents carry), the zone, the relative depths asked, and one or more properties,
each with one slice of values for each depth, in the order of the library's `sub_zones`. A
missing value is `None`.

```rust
use rs4dggs_ogc::{
    Coordinates, Dggrs, Property, Values, ZoneData, ZoneGeometry, write_dggs_json,
    write_dggs_ubjson, write_zone_data_geojson,
};

let dggrs = Dggrs::from_id("ISEA3H").unwrap();
let zone = dggrs.grid().zone_from_text("C2-23-C")?;
assert_eq!(dggrs.grid().count_sub_zones(zone, 1)?, 7);
let height = [Some(-12.5), Some(0.0), None, Some(1234.5678), Some(7.0), Some(f64::NAN), Some(1e17)];
let values = [Values::F64(&height)];
let properties = [Property::new("height", &values)];
let data = ZoneData::new(dggrs, zone, &[1], &properties);

// DGGS-JSON: a missing value, and one that is not finite, are written null.
let mut json = Vec::new();
write_dggs_json(&mut json, &data)?;
assert_eq!(
    String::from_utf8(json)?,
    concat!(
        r#"{"dggrs":"https://www.opengis.net/def/dggrs/OGC/1.0/ISEA3H","zoneId":"C2-23-C","#,
        r#""depths":[1],"values":{"height":[{"depth":1,"shape":{"count":7,"subZones":7},"#,
        r#""data":[-12.5,0,null,1234.5678,7,null,1e17]}]}}"#,
    )
);

// DGGS-UBJSON: the same document, every number big endian. With a value missing, the data are
// a counted container, a marker before each value and `Z` for the missing ones.
let mut ubjson = Vec::new();
write_dggs_ubjson(&mut ubjson, &data)?;
assert!(ubjson.windows(4).any(|w| w == b"[#i\x07"));

// GeoJSON: one feature a sub-zone, its value then its identifier, here without geometry.
let mut geojson = Vec::new();
let written = write_zone_data_geojson(&mut geojson, &data, ZoneGeometry::None, Coordinates::default())?;
assert_eq!((written.features, written.null_zones_skipped), (7, 0));
let geojson = String::from_utf8(geojson)?;
assert!(geojson.contains(r#""id":3,"geometry":null,"properties":{"height":null,"zoneID":"#));
# Ok::<(), Box<dyn std::error::Error>>(())
```

Several depths give one entry a depth in DGGS-JSON and DGGS-UBJSON, as the standard requires, and
in GeoJSON one collection with the features of every depth in the order of the depths, numbered
from 1 throughout.

## The six grids

`Dggrs::all()` lists them in this order, and `Dggrs::from_id` finds one by the identifier a
request names, exactly as written. `grid()` gives the library's grid, built directly: no text of a
client ever reaches `rs4dggs::get_grid`.

| Identifier | The library's name | Title | Registered with the OGC | `crs()` |
|---|---|---|---|---|
| `ISEA3H` | `ISEA3H` | ISEA Aperture 3 Hexagonal | yes | `https://www.opengis.net/def/crs/OGC/0/1534` |
| `IVEA3H` | `IVEA3H` | IVEA Aperture 3 Hexagonal | yes | none |
| `RTEA3H` | `RTEA3H` | RTEA Aperture 3 Hexagonal | no | none |
| `ISEA7H_Z7` | `IGEO7` | ISEA Aperture 7 Hexagonal / Z7 Indexing | no | `https://www.opengis.net/def/crs/OGC/0/1534` |
| `IVEA7H_Z7` | `IVEA7H` | IVEA Aperture 7 Hexagonal / Z7 Indexing | no | none |
| `RTEA7H_Z7` | `RTEA7H` | RTEA Aperture 7 Hexagonal / Z7 Indexing | no | none |

**Three names are refused**: `ISEA7H`, `IVEA7H` and `RTEA7H`. In DGGAL they denote grids of
another indexing (the rhombus identifiers, with the letters A to H), and so do the first two in
the OGC register; yet `rs4dggs::get_grid` accepts `IVEA7H` and `RTEA7H` for the library's Z7
grids. A service that passed a client's text to the library would therefore answer for the wrong
grid. `IGEO7` is accepted as a second name of `ISEA7H_Z7`; the entry's `id()` is then not the text
asked, which tells a service that a second name was used. The match is exact: `isea3h` is
refused, and a service that wants other casings puts the text in upper case before it asks.

**Four of the six are not registered with the OGC** (of the hexagonal grids of this family, the
register, read on 2026-10-02, holds `ISEA3H`, `IVEA3H`, `ISEA7H` and `IVEA7H`): `RTEA3H`,
`ISEA7H_Z7`, `IVEA7H_Z7` and `RTEA7H_Z7`. No definition of them was found, at the OGC or elsewhere,
so their definitions in `definitions/` were written for this crate, from the registered definitions
of IVEA3H, ISEA7H and IVEA7H and from DGGAL; each says so in its own `description`, and that it
holds modifications that the OGC has neither approved nor adopted. The definitions of `ISEA3H` and
`IVEA3H` are the register's files, byte for byte. The URI of each of the six, `uri()`, is of the
registered form, `https://www.opengis.net/def/dggrs/OGC/1.0/{id}`, and is the one a data document
carries: the last segment is the name DGGAL's own tools recognise. For the four it anticipates a
registration that has not been made. `registered_uri()` is `Some` for `ISEA3H` and `IVEA3H` alone,
and is what a grid's description shows, since the standard asks for that member only of a grid
registered with an authority.

## Byte orders

The two binary encodings have opposite byte orders, each as its own specification requires: the
list of 64-bit integers is **little endian** (a count, then each identifier, eight bytes each),
and DGGS-UBJSON, as every UBJSON document, is **big endian**.

## The types of DGGS-UBJSON

Neither the standard nor DGGAL says which UBJSON types the values take, and DGGAL 0.0.6 reads no
UBJSON. The choice made here:

- **The container of the data.** Where every value of a slice is present and finite, a typed
  container (`[$d#n` for `f32`, `[$D#n` for `f64`, and so on): the type once, then the bare
  values. Otherwise the counted container (`[#n`), a type marker before each value and `Z` for a
  value that is missing or not finite, since a typed container cannot hold a null. A reader of
  UBJSON must accept both forms. The difference is not small: 5,536 values of `f32` take 22,144
  bytes of payload in a typed container, and 27,668 in a counted one with three of them missing.
- **A band of floats** keeps its own type: `d` for `f32`, `D` for `f64`.
- **A band of integers** is written in the smallest signed type that holds every value of the slice:
  `i` (8 bits), `I` (16), `l` (32) or `L` (64). UBJSON has no unsigned type beyond 8 bits, and its
  unsigned 8-bit type is never used for data, since `py-ubjson` reads a typed container of it as a
  string of bytes. The cost is one pass over the slice for its least and greatest values, and a type
  that may differ between two responses of one band (a band of `u16` is written `i` or `I` while its
  values stay below 32,768, `l` otherwise); the values are the same. The other course, one fixed
  type for each kind of band (`u8` as `I`, `u16` as `l`, `u32` as `L`), would keep one type always,
  at the price of four bytes a value for a band of `u16` and eight for a band of `u32`, twice their
  size.
- **The other members**: `depths` a counted array of `U` (unsigned 8-bit), each depth a `U`, the
  counts and lengths in the smallest type that holds them.

These types were read back by `py-ubjson` 0.16.1, which returns every document as the DGGS-JSON of
the same data (each `f32` widened to a double, equal once rounded back); no other reader has been
tried.

## What a caller must know

- **The null zone** (`ZoneId::NULL`) is refused with `Error::NullZone` by every writer that takes
  a zone. A list that may hold it is filtered by the caller. In zone data the crate walks the
  sub-zone order itself: there, an entry of the order that is the null zone (on an aperture-7 grid,
  among the sub-zones of level 15 and finer at the two broken seams of the engine) keeps its
  position in DGGS-JSON and DGGS-UBJSON with whatever value the caller gives it (`None` writes
  `null`), and in GeoJSON is passed over, without a number, and counted in
  `Written::null_zones_skipped`. A client that links DGGAL reads those two positions as a zone
  the engine names there, far away (at one seam tested, a zone at 85.4 N 88.6 W).
- **A zone the grid does not define** is refused by the writers of zone data, with
  `Error::Grid`, before anything is written. The writers of lists, and of a zone's geometry or
  feature, check nothing but the null zone: an identifier the grid does not define is written as
  the library's `text_id` renders it, which may be the text of another zone, and in the binary
  list as the integer it is; its geometry is what the library gives for it (in a trial, no region
  and the point (0, 0) as its centroid). Pass these writers only zones the library gave: from
  `zone_from_text`, `zones`, `zones_in_box`, `compact_zones` or `sub_zones`.
- **A zone that cannot be drawn** (a ring that is empty, lies wholly on a pole, has fewer than
  three distinct points, or a coordinate that is not finite; a centroid that is not finite) is
  written in a list or in zone data as a feature with `"geometry":null`, and counted in
  `Written::without_geometry`, so that a list holds a feature for every zone given.
  `write_zone_geometry` writes nothing for it and answers `false`, so that the service leaves the
  member out. In the sample described under "How it is verified", three zones of level 19 on each
  aperture-7 grid have a ring that is one point repeated.
- **Coordinates are written with eight decimals** by default, about a millimetre on the ground,
  and never with fewer: `Coordinates::decimals(n)` with `n` below eight writes eight.
  `Coordinates::shortest()` writes the shortest form that reads back to the same double;
  `places()` answers the decimals applied.
  The floor is measured: rounded more coarsely, the rings of fine zones were written invalid in
  great numbers, the rounding making them cross or touch themselves (with six decimals, 17 to 27
  per cent of the aperture-7 zones of levels 17 to 19, and 16 to 23 per cent of the aperture-3
  zones of levels 31 to 33, in a sample at middle latitudes; with seven, 19 per cent of the
  aperture-7 zones of level 19). With decimals the ring is rounded before it is made a polygon,
  and a vertex at which the rounded ring turns back on itself to within half a rounding step is
  removed, taking a sliver narrower than the precision asked. In the sample described under
  "Rings the engine draws crossing", the only shapes a strict reader refuses at eight decimals are
  those of rings the engine itself draws crossing, and twelve zones of the finest aperture-3
  levels at the poles are written without geometry, their rounded ring enclosing nothing.
- **A coordinate that is not finite** cannot be written in GeoJSON. A geometry in which the
  cutting at the antimeridian made one from finite coordinates is refused with
  `Error::NotFinite`; no ring of the library's grids does so.
- **Across the antimeridian**, a zone's region is cut into as many polygons as it has pieces on
  either side, a `MultiPolygon` where there are two or more, every longitude within -180 to 180;
  a ring about a pole is closed along the pole. This holds for every document that carries a
  geometry, the zone information's `geometry` member among them, where the standard's schema
  names a `Polygon` alone.
- **Rings the engine draws crossing.** Some refined rings, as the engine gives them, cross or
  touch themselves; they are cut and closed as any other, not mended, and a strict reader of
  GeoJSON may refuse the shape. In the sample of 14,238 zones described below, 32 such rings were
  found (4, 4, 7, 4, 4 and 9 on ISEA3H, IVEA3H, RTEA3H, ISEA7H_Z7, IVEA7H_Z7 and RTEA7H_Z7): 16
  that cross themselves away from a pole, and 16 beside a pole, most of them rings nearly all of
  whose points the engine puts on the pole itself. All lie at levels 15 to 19 on aperture 7, and
  at levels 26, 32 and 33 on aperture 3. A wider sample of 70,562 zones, tried once, held 53, of
  the same kind. These are measurements of samples that seek the hard places, not of whole
  levels; every other ring of both samples gave a valid shape. These figures are of the shortest
  form. With eight decimals, the removal of a ring's spikes makes valid eight of these rings on
  the RTEA grids, and leaves the four of each aperture-3 grid that go out and back along the pole
  with nothing to draw, so that they are written without geometry: 0, 0, 0, 4, 4 and 4 shapes
  are then refused.
- **Bit 63.** The 64-bit identifiers are unsigned. Those of the zones under the zones `08` to `11`
  of level 0 of the three aperture-7 grids have their most significant bit set, so that a client
  that reads the eight bytes as a signed integer sees a negative number for them.
- **The area of a list** (`returnedAreaMetersSquare`) is the caller's to compute: the sum of the
  library's `area` over the zones of the page **before** compaction, which counts once an area
  that compacted zones would count twice. For `A6-0-C` and `AA-0-B` on ISEA3H it is
  31170676883138.758 (the standard's Annex C prints 31170676883138.8); for the 32 zones of level
  1 summed in Annex C's order, 510065621724088.75 (Annex C: 510065621724088.5). A value that is
  not finite is left out, the member being optional.
- **Refused before the first byte, or left unfinished.** Data are held first to what the
  document requires and refused before anything is written: values whose number is not that of
  the sub-zones at their depth (`Error::Shape`), which a client would otherwise read in silence
  with each value on another zone; a property without one slice for each depth
  (`Error::Depths`); a property named `zoneID` (`Error::ReservedName`), which zone data as
  GeoJSON cannot carry beside the member of that name, in all three encodings alike; two
  properties of one name (`Error::RepeatedName`); a zone or a depth the library refuses
  (`Error::Grid`), and, in GeoJSON, a list of sub-zones or a refinement beyond the library's
  limits. A failure of the sink leaves the document truncated where it arose. Each value is a
  small write: a service passes a buffered sink.
- **DGGAL's own reader reads one depth.** `dgg togeo` reads only the first depth of a DGGS-JSON
  document that holds several; the other depths are not wrong for that.
- **The public types may grow in a later release.** The structs a caller gives the writers
  (`Coordinates`, `Link`, `LinkTemplate`, `Property`, `ZoneData`) are built by their
  constructors, not from their fields: `Link::new(rel, href)` with `with_title` and
  `with_media_type`, `LinkTemplate::new(rel, uri_template)` with `with_title`,
  `Property::new(name, depths)` and `ZoneData::new(dggrs, zone, depths, properties)`. The counts
  of `Written` are read by name. A `match` on `ZoneGeometry`, `Values` or `geometry::Shape`
  needs an arm for the variants not yet known. Since the constructors borrow their slices, the
  values and properties of data used in more than one statement are bound first, as in the
  examples above.

## How it is verified

Without DGGAL, in every run of the gate (and on GitHub):

- Every document is compared byte for byte with a golden file under `tests/goldens/`, read once by
  eye when it was made; the goldens of DGGS-JSON were composed independently of the writer. Every
  JSON and GeoJSON golden and definition is parsed as JSON by Python's own reader, which refuses
  `NaN`, `Infinity`, a repeated member and invalid UTF-8 here.
- The UBJSON documents are read back by a reader in the crate's own tests into the JSON text of the
  same data, byte for byte equal to the DGGS-JSON goldens.
- The making of a geometry is held by a census: on each of the six grids, at every level, the zone
  at each of fifteen points (both poles, four points beside them, seven on the antimeridian, the
  vertex of the icosahedron at 58.4 N 11.2 E, and Lisbon) and its neighbours, 14,238 zones in all.
  A judge in the tests, which needs no geometry library, holds every shape: each ring closed, of
  four positions or more, within range, anticlockwise and simple. The counts above are asserted
  exactly, zone by zone. The judge was compared with GEOS (Shapely 2.1.1, GEOS 3.13.1) over
  198,016 shapes, again over 28,458 after a later correction, and again over the census with
  eight decimals and over a further 12,192 zones written three ways: it agreed on every one.
- The registered definitions and the standard's two files kept as fixtures (its DGGS-JSON schema
  and its first example) are held, in the gate, to the checksums of the copies fetched from the
  OGC.

With DGGAL v0.0.6 configured (`--features oracle`; the build is the one the library's own suites
use), the crate is compared with DGGAL's own tool, `dgg`:

- the binary list: for 100 zones of each grid (600 in all, among them the zones of level 0,
  pentagons, and zones with bit 63 set), the identifier written is `dgg info`'s 64-bit integer;
- the JSON list: `dgg list 1` on ISEA3H gives the golden's 32 zones in its order;
- the six grids: `dgg` accepts the last segment of each URI as the grid's name, and agrees on its
  refinement ratio, finest level and default depth;
- a zone's region: for 50 zones of each grid away from the antimeridian and the poles (300 in
  all), the polygon written has as many positions as `dgg geom` draws, each within 1e-9 degrees
  (the largest difference found was 9.15e-11);
- DGGS-JSON: `dgg togeo` reads 84 documents (28 zones of the six grids, hexagons and pentagons, at
  depths 1, 2 and 3, more than 5,000 values) and puts each value on the zone the library has at
  that position, which verifies the writer and the sub-zone order together; at a broken seam it
  names a zone at the two positions where the library has the null zone, as said above;
- zone data as GeoJSON: the features of the same 84 documents equal `dgg togeo`'s, identifier and
  value, feature for feature;
- the DGGS-JSON schema: every golden of zone data, and the standard's example, are validated
  against the standard's own schema by Python's `jsonschema`. The schema as published refuses the
  `null` the standard requires for a missing value (it says `nullable`, a word of OpenAPI 3.0
  that a validator of JSON Schema 2020-12 ignores), so it is corrected in memory in that one place
  and the file stays as published. **This step runs only with DGGAL**: a gate without DGGAL, and
  GitHub, never validate against the schema; the member names are held there by the goldens.

DGGAL writes no DGGS-UBJSON and no binary list of its own: for UBJSON there is no engine to
compare with. The standard's printed responses (its Annex C) serve for the names of members,
identifiers and their order alone: its coordinates are those of another build of the engine,
1.3e-9 degrees from these.

## Rust version

`rust-version = "1.85"`, as the library's.

## Licence

MIT (see `LICENSE`). `NOTICE` gives the OGC's notice and licence for the two definitions and the
two files of the standard that the crate carries, and for the four definitions written from them,
and DGGAL's notice, which applies through the library.
