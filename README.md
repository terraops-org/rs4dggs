# rs4dggs

A pure Rust port of [py4dggs](https://github.com/terraops-org/py4dggs), which is itself a Python
port of DGGAL's discrete global grid systems: IGEO7, IVEA7H, RTEA7H, ISEA3H, IVEA3H and RTEA3H.
DGGAL, not py4dggs, is the reference for both ports (see Why it agrees, below). The crate has no
dependencies and no `unsafe` code.

## Getting started

Add the crate with `cargo add rs4dggs`, or depend on the repository by `git` or by `path`.

A grid is chosen once: at compile time with `igeo7()`, `ivea7h()`, `rtea7h()`, `isea3h()`,
`ivea3h()` or `rtea3h()`, or by name at run time with `get_grid("igeo7")`, in any casing.
Everything else starts from a zone, either a coordinate quantised at a resolution, from 0 to 19 on
the aperture-7 grids and from 0 to 33 on the aperture-3 grids, or a text identifier read back.
From a zone come its centroid and boundary, the boundary again with its edges refined, its extent
and its area, its neighbours and the disk around it, its place in the hierarchy and its
sub-zones. From a grid come the facts of its levels: how many zones a level
holds, what area a zone of it is given, and which level suits a pixel of a given size; and from a
grid come as well the zones themselves: every zone of a level, the zones of a level that meet a
bounding box, each in DGGAL's own order and resumable after any zone, an upper bound on how many
zones a box holds, and DGGAL's compaction of a set of zones. The first example below shows the
zone and the facts in turn, and the second the zones of a box.

```rust
use rs4dggs::{GRID_NAMES, get_grid, igeo7, isea3h, ivea3h, ivea7h, rtea3h, rtea7h};

// A grid named at compile time.
let grid = igeo7();
let zone = grid.zone_from_geo(38.7223, -9.1393, 5)?;
assert_eq!(zone.text_id(), "0064156");

let centroid = zone.centroid();          // WGS84 degrees
let boundary = zone.vertices();          // six points for a hexagon, five for a pentagon
assert_eq!(boundary.len(), 6);

let neighbours = zone.neighbors();       // the six (or five) edge neighbours
assert_eq!(neighbours.len(), 6);

// The hierarchy is DGGAL's, and it is not congruent: a zone need not be contained in the zone
// that its identifier names with one digit fewer. This one lies across the boundary of two zones
// of the resolution above and has both as parents, the one its digits name first; and it has
// thirteen children, the seven its digits name and six that it shares with its neighbours.
let parents: Vec<String> = zone.parents().iter().map(|p| p.text_id()).collect();
assert_eq!(parents, ["006415", "006414"]);
let parent = zone.parent().unwrap();     // the first of them
assert_eq!(parent.text_id(), "006415");
assert!(parent.is_ancestor_of(&zone));
assert!(!zone.is_sibling_of(&zone));     // a zone is not its own sibling, as DGGAL has it
let children = zone.children();
assert_eq!(children.len(), 13);
assert_eq!(children[0].text_id(), "00641560");   // the centroid child, with the one parent
assert!(children[0].is_centroid_child() && !zone.is_centroid_child());
assert_eq!(children[0].parents(), [zone]);

// The disk of radius 2: a lazy walk from the centre outwards, ring by ring, in a
// defined order. Only what is taken is computed, so any radius may be asked for.
let disk: Vec<_> = zone.disk(2).collect();
assert_eq!(disk.len(), 19);
assert_eq!(disk[0], zone);
let ring_sizes: Vec<usize> = zone.disk(2).rings().map(|ring| ring.len()).collect();
assert_eq!(ring_sizes, [1, 6, 12]);
let first_seven: Vec<_> = zone.disk(u32::MAX).take(7).collect();
assert_eq!(first_seven[..], disk[..7]);

// The boundary with its edges divided into parts, as the engine draws a zone: 0 asks for the
// engine's own choice, twelve parts to an edge at this resolution, and so 72 points. The ring
// runs anticlockwise, is not closed, and its longitudes may pass 180 degrees (see What a
// caller must know).
let ring = zone.refined_vertices(0)?;
assert_eq!(ring.len(), 72);
assert_eq!(ring[0], boundary[0]);

// The extent: the least and the greatest latitude and the western and eastern ends, in degrees.
let extent = zone.extent().unwrap();
assert!(extent.ll.lat < 38.7223 && 38.7223 < extent.ur.lat);
assert!(extent.ll.lon < -9.1393 && -9.1393 < extent.ur.lon);

// The area in square metres is nominal, the same for every hexagon of a resolution.
assert!((zone.area().unwrap() - 3.03484e9).abs() < 1e4);

// The facts of a grid: the zones of a level, and the level for a pixel of 100 metres.
assert_eq!(grid.count_zones(5)?, 168_072);
assert_eq!(grid.level_from_meters_per_sub_zone(100.0, 0), 12);
assert_eq!(grid.max_neighbors(), 6);

// The same point under the other two projections. At this resolution the three cells
// happen to coincide; four resolutions further down they no longer do, since the three
// projections place the same coordinates differently within the icosahedron's faces.
assert_eq!(ivea7h().zone_from_geo(38.7223, -9.1393, 5)?.text_id(), "0064156");
assert_eq!(rtea7h().zone_from_geo(38.7223, -9.1393, 5)?.text_id(), "0064156");
assert_eq!(igeo7().zone_from_geo(38.7223, -9.1393, 10)?.text_id(), "006415654630");
assert_eq!(ivea7h().zone_from_geo(38.7223, -9.1393, 10)?.text_id(), "006415650050");
assert_eq!(rtea7h().zone_from_geo(38.7223, -9.1393, 10)?.text_id(), "006415650342");

// The sub-zones of the aperture-7 zone, in DGGAL's order of scanlines, which has no relation
// to the digits of the identifiers: the centroid child stands in the middle of the thirteen,
// and two levels down there are 55 sub-zones, where the digits name 49 zones.
let zone = igeo7().zone_from_geo(38.7223, -9.1393, 5)?;
let subs = zone.sub_zones(1)?;
assert_eq!(subs.len(), 13);
assert_eq!(subs[6].text_id(), "00641560");
assert_eq!(zone.sub_zone_index(&subs[6])?, Some(6));
assert_eq!(zone.count_sub_zones(2)?, 55);
assert_eq!(zone.sub_zone_at_index(14, 1)?.resolution(), 19);

// An aperture-3 grid. Its hierarchy is DGGAL's as well: this zone has three parents and
// seven children.
let zone = isea3h().zone_from_geo(38.7223, -9.1393, 5)?;
assert_eq!(zone.text_id(), "C2-23-C");
let parents: Vec<String> = zone.parents().iter().map(|p| p.text_id()).collect();
assert_eq!(parents, ["C2-23-A", "C4-5-A", "C4-6-A"]);
assert_eq!(zone.children().len(), 7);
assert_eq!(ivea3h().zone_from_geo(38.7223, -9.1393, 15)?.text_id(), "H2-19C2D2-B");
assert_eq!(rtea3h().zone_from_geo(38.7223, -9.1393, 15)?.text_id(), "H2-19C2D3-D");

// Its sub-zones, in the grid's own order: listed whole at a modest depth, or paged by index at
// any depth down to the finest resolution without being listed.
let subs = zone.sub_zones(2)?;
assert_eq!(subs.len() as u64, zone.count_sub_zones(2)?);
assert_eq!(zone.sub_zone_at_index(2, 5)?, subs[5]);
assert_eq!(zone.sub_zone_index(&subs[5])?, Some(5));
assert_eq!(zone.sub_zone_at_index(28, 1_000_000)?.resolution(), 33);

// Or a grid named at run time, in any case, which answers in identifiers rather than
// zones.
assert_eq!(
    GRID_NAMES,
    ["IGEO7", "IVEA7H", "RTEA7H", "ISEA3H", "IVEA3H", "RTEA3H"]
);
let any = get_grid("rtea7h")?;
assert_eq!(any.name(), "RTEA7H");
assert_eq!(any.max_resolution(), 19);
let id = any.zone_from_geo(38.7223, -9.1393, 10)?;
assert_eq!(any.text_id(id), "006415650342");
assert_eq!(get_grid("ISEA3H")?.max_resolution(), 33);
# Ok::<(), rs4dggs::Error>(())
```

The zones of a box are an iterator in DGGAL's order, taken in pages: a page ends at a zone, and
the next page is entered after it, at a cost that does not grow with its place in the answer.

```rust
use rs4dggs::{Extent, GeoPoint, igeo7, isea3h};

// The box about Lisbon: from 38 to 39 degrees north, and from 10 to 9 degrees west.
let lisbon = Extent {
    ll: GeoPoint { lat: 38.0, lon: -10.0 },
    ur: GeoPoint { lat: 39.0, lon: -9.0 },
};
let grid = isea3h();

// An upper bound on how many zones of level 12 meet the box, by arithmetic alone, before any
// zone is visited: what a service refuses a request by.
assert_eq!(grid.estimate_zones_in_box(12, &lisbon)?, 187);

// A first page of fifty zones; then the next fifty, entered after the last zone of the first.
let first: Vec<_> = grid.zones_in_box(12, &lisbon)?.take(50).collect();
let last = *first.last().unwrap();
let second: Vec<_> = grid.zones_in_box(12, &lisbon)?.after(last)?.take(50).collect();

// The pages are the answer in order, each zone once.
let all: Vec<_> = grid.zones_in_box(12, &lisbon)?.collect();
assert_eq!(all.len(), 128);
assert_eq!(first[..], all[..50]);
assert_eq!(second[..], all[50..100]);

// DGGAL's compaction of the answer: coarser zones in the place of finer ones.
assert_eq!(grid.compact_zones(&all)?.len(), 88);

// The finest level of the two enumerations: the last of the aperture-3 grids, and 14 on the
// aperture-7 grids, whose last level is 19.
assert_eq!(grid.max_box_level(), 33);
let grid = igeo7();
assert_eq!(grid.max_box_level(), 14);

// Every zone of a level, in the engine's order, which on the aperture-7 grids is not that of
// the identifiers.
let level_0: Vec<String> = grid.zones(0)?.map(|z| grid.text_id(z)).collect();
assert_eq!(
    level_0,
    ["01", "06", "02", "07", "03", "08", "04", "09", "05", "10", "00", "11"]
);
# Ok::<(), rs4dggs::Error>(())
```

Four small programs show the same in context, each run with `cargo run --example <name>` from a
clone of the repository, which needs nothing but Rust; for instance
`cargo run --example quantise -- 5 38.7223 -9.1393`:

- `quantise`: points on all six grids at one resolution, from 0 to 19, the range they share;
- `hierarchy`: a zone's ancestors and children on IGEO7;
- `geojson`: a zone of IGEO7 and its neighbours as GeoJSON, with no dependency;
- `service`: a grid chosen by name, and a disk served ring by ring.

The crate reads two environment variables, each of which lowers a limit and never raises it:
`RS4DGGS_MAX_EDGE_REFINEMENT`, the greatest edge refinement that `refined_vertices` accepts,
described under Departures from the engine; and `RS4DGGS_MAX_MATERIALISED_SUB_ZONES`, the
greatest number of sub-zones that `sub_zones` lists, described under What a caller must know.

`cargo doc --open` renders the full reference. The entry points a newcomer meets first, quantising,
reading an identifier, the boundary, the neighbours, the hierarchy, the disk and choosing a grid by
name, each carry an example of their own.

## The workspace's other crates

The workspace holds three crates: this library and two over it.

`rs4dggs-cli`, in `crates/rs4dggs-cli`, is a command-line tool, `rs4dggs`; its README describes
every command, the formats and batch use. Three examples:

```text
rs4dggs igeo7 zone 38.72,-9.14 10
rs4dggs igeo7 neighbours 006415654636
rs4dggs isea3h geom C2-23-C > cell.geojson
```

`rs4dggs-ogc`, in `crates/rs4dggs-ogc`, writes the encodings of OGC API - DGGS over this library's
six grids: the zone list as JSON, GeoJSON and 64-bit integers, zone data as DGGS-JSON, DGGS-UBJSON
and GeoJSON, and a zone's geometry as a valid GeoJSON geometry, cut at the antimeridian; it holds
the six grids under the identifiers a service gives them, with their URIs and definition
documents. Its README describes each encoding and what a caller must know. The tool takes from it
the making of its GeoJSON polygons.

## Status

All six grids are implemented. The three aperture-7 grids, **IGEO7**, **IVEA7H** and **RTEA7H**,
share one aperture-7 hexagonal topology and one Z7 indexing; the three aperture-3 grids,
**ISEA3H**, **IVEA3H** and **RTEA3H**, share one aperture-3 hexagonal topology and one I3H
indexing. Within each family the grids differ only in the projection laid over the icosahedron,
ISEA, IVEA and RTEA respectively. Each grid offers quantisation from geographic coordinates,
centroids, cell boundaries, edge neighbours, its hierarchy and the relational predicates over it,
the disk of radius k, and text identifiers, at every resolution its indexing can hold: 0 to 19 for
Z7 and 0 to 33 for I3H. Each also gives the boundary of a zone with its edges refined, its
geographic extent and its area, as DGGAL's own accessors give them, and the facts that DGGAL
states of a grid: the number of zones at a level, the reference area of a zone and the size of its
sub-zones at each level, the levels that suit a given area or size, and the greatest numbers of
parents, children and neighbours a zone has.

Each grid also enumerates its zones: every zone of a level (`zones`), and the zones of a level
whose extent meets a bounding box (`zones_in_box`), each as a lazy iterator in the order of
DGGAL's `listZones`, which can be entered after any zone of the level at a cost that does not
grow with the zone's place in the sequence, so that an answer is taken in pages and a cursor is
one zone identifier. With them come an upper bound on the number of zones of a box, by
arithmetic and before any work (`estimate_zones_in_box`), and DGGAL's compaction of a set of
zones of one level (`compact_zones`), the engine's answer zone for zone and in its order. The two
enumerations answer at every level of the aperture-3 grids and to level 14 of the aperture-7
grids (`max_box_level`; see Known limitations).

The hierarchy is DGGAL's own on every grid, and on none of them is it congruent: a zone need not
be contained in one parent. On the aperture-7 grids a zone has one parent or two, and thirteen
children, or eleven under a pentagon; on the aperture-3 grids one parent or three, and seven
children, or six under a pentagon. Every grid orders its sub-zones at any depth, in DGGAL's own
order, and pages through that order by index without building it. On the aperture-3 grids the
order is corrected where the engine's own is at fault; on the aperture-7 grids it is the
engine's own, entry for entry, with the null zone where the engine's entry is an identifier it
cannot read back, and what departs from the engine there is the answer of single calls that the
engine itself gets wrong (see Departures from the engine).

The identifiers `neighbors` and `quantize` follow py4dggs's spelling.

## What changes from 0.1.0

Version 0.1.0 answered the hierarchy of the three aperture-7 grids from the digits of the Z7
identifier: a zone had one parent, the identifier with its last digit dropped, and seven
children, one for each digit, or six under a pentagon. That is a congruent hierarchy, and it is
not the grid's. A hexagon cannot be cut into seven smaller hexagons, so a child is not contained
in the zone that its digits name: measured on the engine, 91.4 per cent of the area of the IGEO7
zone `006415654636` lies in `00641565463` and 8.6 per cent, one corner, in `00641565462`. DGGAL
names as a zone's parents the coarser zones that hold a part of it, and from this version the
aperture-7 grids answer as DGGAL does, under the names the crate already had, as the aperture-3
grids have done from the start. No signature changes. On IGEO7, IVEA7H and RTEA7H these answers
change:

- **`parents`** lists one zone or two, where it listed one. `006415654636` has `00641565463` and
  `00641565462`. The first is the zone that the identifier names with its last digit dropped, so
  that `Zone::parent` answers as before, save at some identifiers of the two broken edges (see
  Known limitations). A centroid child, whose last digit is 0, lies wholly within that zone and
  has it alone.
- **`children`** lists thirteen zones, or eleven under a pentagon, where it listed seven, or six,
  and lists them in the engine's order and no longer in the order of the digits. Those of
  `006415654636` are `0064156546360`, `...61`, `...65`, `...64`, `...66`, `...62` and `...63`, the
  seven that its digits name, the centroid child first, and then `0064156546342`,
  `0064156546033`, `0064156546251`, `0064156546215`, `0064156546324` and `0064156546306`, which
  lie across its boundary and are children of a neighbour as well.
- **`is_centroid_child`** is true where the last digit of the identifier is 0, and false of a
  zone of resolution 0; it was true of every zone. `006415654636` is none.
- **`centroid_parent`** is the first of the parents that is itself a centroid child, and `None`
  where neither is one; it was the one parent. `006415654636` has none, and `00641565463601` has
  `0064156546360`.
- **`is_ancestor_of`, `is_immediate_child_of` and `is_sibling_of`** follow the parents: a zone is
  the immediate child of each of its parents, an ancestor is reached through either, and two
  zones are siblings where they share either. `0064156` has `00422` and `00645` among its
  ancestors beside `00641`, and `0064141` is its sibling, which under the digits it was not.
- **A Z7 identifier of twenty digits**, bound directly with `Grid::zone`, has no parents, as in
  the engine; it had the parent that its digits name. It has no children, as before.
- **`count_sub_zones`, `first_sub_zone`, `sub_zones` and `sub_zone_at_index`** answer, where they
  returned `Error::NoSubZoneOrder`, and `sub_zone_index` answers for a finer zone, where it
  returned the same error. The sub-zones of `006415654636` at depth 1 are its thirteen children,
  in the order of the engine's scanlines. No grid of the crate returns `Error::NoSubZoneOrder`
  any longer; the variant stays, for a topology that defined no order.

On the aperture-3 grids no answer changes. On all six, the limit on a list of sub-zones can now
be lowered by an environment variable (see What a caller must know), which changes nothing while
it is unset. The methods added since 0.1.0, which change no earlier answer, are the refined
boundary, the extent, the area and the facts of a grid, and the enumerations of the zones of a
level and of a box, the estimate and the compaction, all named under Status. With the last come
two iterators, `Zones` and `ZonesInBox`, and four variants of `Error`, which stays
`#[non_exhaustive]`: `LatitudeOutOfRange` and `LongitudeOutOfRange` refuse a bounding box with a
coordinate out of range and name that coordinate, `"south"`, `"west"`, `"north"` or `"east"`, as
`NonFinite` does, and `InvertedBox`, which carries nothing, refuses a box whose south lies above
its north; `MixedLevels` refuses a set of more than one level given to `compact_zones`.
`NonFinite` now also names the coordinate of a box at fault, `"south"`, `"west"`, `"north"` or `"east"`, and `ResolutionOutOfRange` answers a
level beyond `max_box_level` with that limit as its `max`. The command-line tool follows the
library, and its own README lists what that changes in what it prints.

## Verification

Every operation this crate offers is checked against a live DGGAL v0.0.6, directly or through
compositions of the engine's answers, by way of a test-only crate, `dggal-oracle`, that links
DGGAL's own Rust bindings and is never a dependency of the published crate. One suite, generic
over the grid, runs six times, each time against that grid's own DGGAL class, `ISEA7H_Z7`,
`IVEA7H_Z7`, `RTEA7H_Z7`, `ISEA3H`, `IVEA3H` and `RTEA3H`, with tests of its own for what only
the aperture-3 grids have, and it compares `AnyGrid`, the grid chosen at run time, with both the
typed grid and the engine. No ceiling and no admitted disagreement is borrowed from one grid for
another: each was measured on its grid against its own engine.

- **Quantisation** from geographic coordinates is compared at every resolution the grid holds,
  over a main sample of adversarial inputs (both poles, the meridian of the icosahedron's first
  vertex at 11.2 degrees east and its antimeridian, the antimeridian itself, the twelve pentagon
  centres, points a hair from the north pole) and 400 seeded uniform points, and over a sweep of
  1,650 points on the icosahedron's thirty edges and a hair either side of them: 33,000
  quantisations of the sweep on each aperture-7 grid and 56,100 on each aperture-3 grid. On the
  aperture-3 grids a further 352 points along the two seams where the engine answers identifiers
  it cannot read back are quantised at all 34 resolutions (see Departures from the engine).
  Identifiers and their text compare exactly, and both sides read that text back to the same
  identifier.
- **Centroids, vertices and neighbours** are compared over the zones the main sample names at
  fifteen resolutions spread from 0 to 19 on the aperture-7 grids, every one from 13 upwards
  among them, and at seventeen spread from 0 to 33 on the aperture-3 grids, the finest three
  among them, and the neighbours again near both poles and along every edge. On the aperture-7
  grids a further 2,860 digit paths deep under the two poles, 1,430 under each, zeros followed by
  at most three other digits at resolutions 15 to 19, are compared the same way. Identifiers
  compare exactly, and so do neighbour lists, in the engine's own order, on every grid: on the
  aperture-7 grids the engine's own list less the four kinds this crate never hands out (see
  Departures from the engine), each dropped entry re-checked live, and on the aperture-3 grids
  the engine's list entire. So do coordinates: every latitude and longitude must be the engine's
  own double, bit for bit. Vertices are compared in order, on the aperture-3 grids with the
  engine's ring reordered as Departures from the engine describes, and every ring compared must
  run anticlockwise as seen from outside the sphere, save on the aperture-7 grids the rings whose
  every vertex the engine snaps onto a pole, which enclose no area, and the rings among whose
  vertices it puts the point (0, 0) (see What a caller must know), whose other vertices must run
  anticlockwise; both are counted, the second exactly. On the aperture-3 grids the engine's own
  ring must run clockwise at every zone compared but the two polar pentagons, and anticlockwise at
  those. Each test prints the largest gap it measured all the same.
- **The refined boundary, the extent, the area and the grid facts.** At every zone of the main
  sample, widened by a fixed list of identifiers that meets each site where the engine's compiled
  order of operations differs from its source, the ring at edge refinements of 0, 1, 2 and 3 is
  compared as a sequence, in length and in every double, and so are the extent's four doubles
  and the area, through `Grid`, `AnyGrid` and `Zone`. That is 5,624, 5,596 and 5,591 zones with
  geometry on IGEO7, IVEA7H and RTEA7H and 5,937, 5,936 and 5,932 on ISEA3H, IVEA3H and RTEA3H,
  and 1,311,496, 1,302,454 and 1,301,632 coordinates on the first three and 967,582, 964,550 and
  967,096 on the others, with no difference. At one zone of each of eleven kinds on every grid
  the ring is compared again at 2,000, 20,000 and 100,000 parts to an edge, 8,150,414,
  8,148,308 and 8,149,677 points on the aperture-7 grids and 18,130,433, 13,252,704 and
  18,141,833 on the aperture-3 grids, and the refusals above that bound are checked. The grid
  facts are compared at every level, 256 comparisons on each aperture-7 grid and 669 on each
  aperture-3 grid, and `level_from_ref_zone_area` and `level_from_meters_per_sub_zone`, which
  the crate mirrors for every input, over a census of areas spread logarithmically over
  thirty-two decades, each level's own area and its neighbours to the bit, and the special
  values: 801,836 comparisons on each aperture-7 grid and 802,844 on each aperture-3 grid. The
  zones without geometry are classified against the live engine and counted exactly, and what
  each grid's sample meets of the degenerate rings (see What a caller must know) is held
  to exact counts too, so that a change in any of them fails. The suite never asks the engine
  for the ring of an identifier it cannot read on the aperture-3 grids, 34 in each sample, where
  it answers a null array that the binding the suite links does not survive; the extent and the
  area of those are asked and compared.
- **The hierarchy** is the engine's own on every grid, and parents, children, the centroid
  parent and `is_centroid_child` are compared with its answers directly, in its order. On the
  aperture-7 grids the parents and the children are compared as sequences with the engine's
  `getZoneParents` and `getZoneChildren` less the three kinds this crate never hands out (see
  Departures from the engine), each dropped entry re-checked live, over two samples: the zones
  that the main sample's positions give, 5,608, 5,580 and 5,575 on IGEO7, IVEA7H and RTEA7H, and
  2,276 identifiers built by text beside a broken edge, at resolutions 14 to 19, the same on the
  three grids. `is_centroid_child` is compared exactly; the centroid parent with the first of
  the engine's own parents that the engine calls a centroid child (at 91 of the identifiers
  built by text that parent is one the engine cannot read back, and this crate's answer is
  checked there against the parents that remain), while the engine's `getZoneCentroidParent`
  is asserted to answer its null zone at every zone, so that a corrected engine would be
  noticed; and `is_immediate_child_of`, `is_sibling_of` and `is_ancestor_of`
  with the engine's `isZoneImmediateChildOf`, `areZonesSiblings` and `isZoneAncestorOf`, on
  pairs drawn in rotation from each zone's own lists, 50,081, 5,584 and 5,905 pairs on IGEO7,
  49,883, 5,572 and 5,874 on IVEA7H and 49,854, 5,560 and 5,893 on RTEA7H among the zones from
  positions, and 28,968, 2,172 and 2,298 on each grid among the identifiers built by text, each
  held to a floor. What the two broken edges do to these lists is counted by resolution and held
  to an exact table for each grid, and every zone so counted must lie within the band of those
  edges (see Known limitations). A further 3,010 hostile identifiers on each grid are put to the
  same methods without the engine, for the agreement of `Grid`, `AnyGrid` and `Zone` and for
  the absence of a panic. On the aperture-3 grids `is_ancestor_of`, `is_immediate_child_of` and
  `is_sibling_of` are compared over some 37,000 distinct pairs of zones on each grid, each pair
  asked once and the engine's answers counted exactly. `disk`, which has no counterpart in
  DGGAL, is compared in order at radii 0, 1 and 2 with a walk composed of the engine's own
  neighbour lists, less, on the aperture-7 grids, the four kinds dropped from those lists, on
  every grid.
- **The sub-zone methods on the aperture-7 grids.** The count is compared with the engine's
  `countSubZones` and the order with its `getSubZones`, entry by entry: 8,905, 8,866 and 8,890
  orders holding 812,303, 809,714 and 813,512 entries on IGEO7, IVEA7H and RTEA7H. They are the
  orders below the hexagons of the main sample at depth 1, below one in six of them at depth 2
  and one in forty at depth 3 (5,845, 5,806 and 5,830 orders); below the twelve pentagons of
  every resolution from 0 to 18 at depths 1 to 3, and of resolutions 0 and 1 at depths 4 and 5
  as well (696 orders, 466,956 entries); ten orders at depth 4 and six at depth 5 below
  hexagons, one of them beside a broken edge (126,808 entries, the longest order 17,053); and
  below the 2,276 identifiers built by text beside a
  broken edge, at depths 1 to 3 (2,348 orders, 72,426 entries). `sub_zone_at_index` is asked at
  every index of every one of those orders and compared with the entry of the engine's order,
  never with the engine's own search by index, to which the binding the suite links passes no
  depth and which is at fault at the pentagons (see Departures from the engine). `first_sub_zone`
  is compared with the first entry, and the engine's `getFirstSubZone` with the first entry of
  its own order, save at nine pairs on each grid where that call would end the process, which
  are counted and at which this crate's answer is compared with the first entry of the engine's
  order. `sub_zone_index` is compared with the engine's `getSubZoneIndex` for 97,553, 97,137
  and 97,912 entries of those orders and for 42,760, 42,555 and 42,742 zones of a neighbouring
  order, most of them no sub-zones of the zone asked. Below sub-zone resolution 15 every order
  compared is as long as its count, names no zone twice, holds no null zone and, at depth 1, is
  the set of the children, and every entry asked has its place as its index, here and in the
  engine; from resolution 15, at the two broken edges, each class is counted by resolution and
  held to an exact table for each grid. Depth 0 is compared at every zone of the main sample,
  and the engine's answers there classified and counted exactly. A further 1,046, 1,021 and
  1,021 orders at depth 1 are compared beside the other sub-zone checks of the common suite,
  and 3,010 hostile identifiers on each grid are put to every sub-zone method without the
  engine.
- **The sub-zone methods on the aperture-3 grids.** The counts, the first sub-zones, the orders
  and the indices are compared with the engine's, over some 1,500 orders holding some 27,000
  sub-zones on each grid and over sub-zones at every depth down to resolution 33.
  `sub_zone_at_index` is compared with the order it indexes, for the same reason as on aperture
  7; a separate probe found the engine's own `getSubZoneAtIndex` the same at all 4,800 pairs it
  compared. The input this crate refuses, and every answer in which it departs from the
  engine, is checked against the engine's own answer to it, save the aperture-3 sub-zone order
  corrected at the rhombus edges: there the order is checked against the engine's own descendant
  test rather than its answer, and the engine's index, which refers to its own order, is not asked
  (see Departures from the engine). On the aperture-7 grids three of the engine's answers are
  not asked by the suites either: its search by index below a pentagon, which the tests of the
  command-line tool ask of DGGAL's own tool instead; its first sub-zone where the call would end
  the process; and the walk that does not end, from which it does not return.
- **The zones of a box, every zone of a level, the estimate and the compaction.** Each suite
  draws boxes of nine kinds (anywhere, smaller than a zone, large, over the antimeridian, about a
  pole, about a vertex of the icosahedron, a point, a segment, and about the corners of zones),
  in degrees as a caller gives them, at every level that `zones_in_box` answers: 612 boxes on
  each aperture-3 grid and 135 on each aperture-7 grid. Each is put to the engine's `listZones`
  unless the engine is known not to return from it or would spend more on it than the suite
  affords (see Departures from the engine): 494 and 125 are put to it, and the 118 and 10 kept from
  it are counted by kind, exactly, and checked by this crate's own invariants. Every zone of the
  engine's answer is in this crate's, in the same order. Where this crate's answer is the larger
  (see Departures from the engine), each zone the engine leaves out is re-checked live to be a
  zone of the level whose extent, as the engine gives it, meets the box by the engine's own test,
  and the boxes so met are counted by kind, exactly: 67, 66 and 59 boxes on ISEA3H, IVEA3H and
  RTEA3H, in 13 of which the engine answers no zone, and 101, 101 and 90 zones left out; none on
  the aperture-7 grids. At the levels next beyond those that the tests of the library take whole,
  7 and 8 on the aperture-3 grids and 5 on the aperture-7 grids, the answers of 100 and 50 boxes
  (from 30,113 to 31,091 zones on a grid) are compared with the rule taken over the engine's
  whole list of the level and its own extents, and `zones` with that list. Boxes laid to overlap
  the extent of a zone by less than 2e-15 radian, by between that and the 2.2e-15 of the
  engine's source, and by more find the zone in this crate's answer exactly where the margin
  that the library holds puts it; where the two margins part, the engine's own answer holds the
  zone, save at 2, 2 and 1 boxes on the aperture-3 grids, where it is a zone the engine leaves
  out, as above. The
  estimate is held not under the count of 559 whole answers on each aperture-3 grid and 125 on
  each aperture-7 grid, and of the cap above 70 degrees north, 54,362, 54,356 and 54,330 zones of
  level 11 on ISEA3H, IVEA3H and RTEA3H, estimated at 55,284, and 36,180, 36,186 and 36,170 of
  level 6 on IGEO7, IVEA7H and RTEA7H, estimated at 36,934. `compact_zones` is compared with the
  engine's `compactZones`, zone for zone and in its order, on the answers of the boxes (224, 225
  and 215 sets on the aperture-3 grids and 97, 94 and 92 on the aperture-7 grids, and 53, 54 and
  54, and 23, 25 and 24, of them again with one zone in twenty taken out) and on the children or
  the sub-zones of zones of the suite's own sample, which reaches every level, the poles, the
  seams and the pentagons (417 sets on each aperture-3 grid and 433, 439 and 427 on the others).
  Every count above is held exactly. A sweep of hostile input, 2,000 boxes, 9,030 cursors and
  from 9,430 to 9,451 slices on each grid, finds an error or an answer and never a panic.

A run of 2026-10-03, of 42 tests on each aperture-7 grid and 35 on each aperture-3 grid,
compared 909,189 cases on IGEO7, 908,818 on IVEA7H, 908,821 on RTEA7H, 1,015,533 on ISEA3H,
1,015,432 on IVEA3H and 1,015,183 on RTEA3H, 5,772,976 in all, each test asserting a floor on
the number it compared; of these, the levels from an area account for 801,836 on each aperture-7
grid and 802,844 on each aperture-3 grid. A case of that sum is one zone, point, order or input
compared. The sum leaves out the tests of the aperture-7 sub-zone order, of the aperture-7
hierarchy at identifiers built by text, of the zones of a box, the estimate and the compaction,
and of hostile input, which count orders, entries, pairs, boxes, sets and identifiers, and whose
figures are given above. In that run no
tie remained in any sample, that is, no input on a cell boundary at which the two sides name
different cells, and the suite holds its admission of ties at none, so that a tie which
reappears fails. The suite asserts that every geographic
point it compares with the engine's is the engine's own to the last bit, latitude and longitude
alike. When they were last counted, in the run of 2026-10-02, all 355,162 such points were
(67,466 on IGEO7, 67,227 on IVEA7H, 67,179 on RTEA7H, 51,117 on ISEA3H, 51,106 on IVEA3H and
51,067 on RTEA3H); the tests that compare them have not changed since, and in the run of
2026-10-03 the largest coordinate gap measured by each geometric test was zero again. The
boundary and extent test, which
compares raw bits of its own, adds the coordinates counted above. That assertion is made for
x86-64 Linux with glibc, where the suites run, and for no other platform.

A disagreement is admitted in two ways only. Either it is a characterised divergence, in which
this crate's answer has been shown to be the defensible one and the engine's behaviour is
re-checked live at every instance on every run; or it lies inside the region along two
icosahedron edges where the aperture-7 engine is not consistent with itself (see Known
limitations), where each case is counted and held under a ceiling. On the aperture-3 grids no
such region was met, and every disagreement admitted there is a characterised divergence.
Anything else fails.

Beside the six suites, four unit tests compare what no public method returns. The first compares
the refined cell boundary over which the aperture-7 containment test is performed with the
nearest accessor DGGAL exposes, at 360 zones fixed by construction on each aperture-7 grid. The
2,666 vertices of 343 of those zones are compared bit for bit, with one exception pinned by zone,
index and bit pattern; the other 17 zones, at which DGGAL's accessor subdivides a side twentyfold
and so traces a different polygon, are counted exactly, so that none can leave the comparison
unnoticed. The second compares the points from which the aperture-3 sub-zone orders are
generated with the engine's own, bit for bit: 4,105,210 of them, over 17,710 pairs of a zone and
a depth on ISEA3H, save at 6,884 points in 98 of those pairs where the engine's own order is at
fault and this crate's is corrected (see Departures from the engine), and there checked against
the definition of a sub-zone instead. The third does the same for the aperture-7 orders, with no
exception: the centroids that the generator gives in the plane of the layout are the engine's
`getSubZoneCRSCentroids`, bit for bit, 12,567,960 of them over 3,240 pairs of a zone and a depth
on IGEO7 (every zone of resolutions 0 to 2 and the twelve pentagons of each resolution from 3 to
8, at depths 1 to 5), and each is found again by its index. The fourth compares the two parts of
the search for an index on aperture 7, the test of whether a zone is a sub-zone of another and
the walk that finds its place, with the engine's `zoneHasSubZone` and `getSubZoneIndex`, answer
for answer, at 157,016 pairs of zones on IGEO7.

A further unit test on each grid compares the refined boundary with the engine's, as the suites
do but over another census: the ring at edge refinements of 0 to 3 and the extent, bit for bit,
over every zone to level 3 on the aperture-7 grids and to level 4 on the aperture-3 grids, and at
every level up to the finest over the zones that contain both poles, 136 points beside them, 87 on
and beside the antimeridian, points along each of the thirty edges of the icosahedron and
pseudo-random points. It met 10,107, 10,090 and 10,103 zones and 2,533,434, 2,529,576 and
2,533,310 coordinates on IGEO7, IVEA7H and RTEA7H, and 7,613, 7,616 and 7,626 zones and
1,340,292, 1,337,242 and 1,344,734 coordinates on ISEA3H, IVEA3H and RTEA3H, and no difference.

Tests of their own compare the enumerations, the estimate and the compaction more widely than
the suites, against the engine where it is configured and by this crate's own invariants
everywhere. `zones` gives the engine's whole list of a level, entry for entry and never sorted,
at levels 0 to 8 on the aperture-3 grids and 0 to 5 on the aperture-7 grids, 883,560 entries,
whole and in pages entered after the last zone of the one before. At levels 0 to 6 on the
aperture-3 grids and 0 to 4 on the aperture-7 grids, the answers of 31,980 boxes of the nine
kinds are the rule
taken over the engine's whole list and its own extents, set and order; at the levels beyond,
to the last that each grid enumerates, every zone of the answers of 4,995 boxes meets the box
by the engine's own extent, the zones beside the answer that meet it are in it, and the pages
are the answer; and some 3,000 of those answers, the ones taken whole of the boxes with a height
and a width (at least 730 on each aperture-3 grid and 280 on each aperture-7 grid), hold the
zone at the centre of the box. Of all those boxes, 9,817 are put to
the engine's `listZones`, and its answer lies within this crate's, in its order, in every one.
The estimate is not under the count of any of the 36,431 answers of those boxes, nor of 1,530
boxes on each aperture-3 grid and 675 on each aperture-7 grid drawn at every level, nor of caps
about the north pole of 251,422 zones at level 7 of IGEO7 and 1,456,540 at level 14 of ISEA3H,
estimated at 253,432 and 1,461,434. `compact_zones` gives the engine's answer, zone for zone
and in its order, on 6,536 sets on the aperture-3 grids, at levels 2 to 8, 12, 20 and 33, and
1,505 on the aperture-7 grids, at levels 2 to 6: the zones of boxes, the same with one zone in
twenty taken out, the sub-zones of zones and whole levels; 402 of the aperture-3 sets, all with
a gap, and none of the aperture-7 sets are covered beyond themselves. On the aperture-7 grids
it gives the engine's answer, too, on 118 sets of level 16 on each grid beside a broken edge,
and on 370 identifiers of levels 17 and 19 built there from digits, each compacted alone, less,
at 35 of them, the identifier of the engine's answer that the engine cannot read back.

**All of this holds on x86-64 Linux with glibc**, where the engine and this crate call the same
`libm`. Nothing is known about wasm32, or about any target whose floating-point functions come
from another library.

The suites ship in the published package, under `tests/`, but run only inside the source
repository, with DGGAL installed and the `oracle` feature enabled. That feature is meaningful only
there: it compiles the suites against `dggal-oracle`, a test-only crate that the published
package does not carry, so in an unpacked copy of the package `cargo test` runs the unit tests and
the doctests alone, and needs no DGGAL. To run the verification, clone
<https://github.com/terraops-org/rs4dggs>, point `DGGAL_SITE_PACKAGES` at a Python environment
holding `dggal==0.0.6` (`pip install dggal==0.0.6`), and from the clone run:

```sh
DGGAL_SITE_PACKAGES=/path/to/site-packages ./scripts/check.sh
```

## Why it agrees

DGGAL is compiled with `-O2 -ffast-math`, and the eC marks the functions that must escape it with
`__attribute__ ((optimize("-fno-unsafe-math-optimizations")))`. In the functions that carry no
such mark, and in the helpers that marked functions call, which gcc compiles out of line under the
file's own options, the compiler was free to re-associate, and in places the shipped library
computes in an order its source does not describe. Where that order was found to differ from the
source, this crate as a rule performs the compiled order, even where no answer on the registered
grids turns on it, as in the wrap of the longitude that ends the inverse projection, and each such
site cites its address in the library, whose BuildID is `e75d6ab18f8b460713ebcd21df6f07922fb909c3`
(DGGAL v0.0.6 as distributed in the `dggal` Python wheel). The few sites that keep another form,
the source's own or py4dggs's, say so, each with the evidence that it changes no answer on the
registered grids. Everywhere else the eC's source text decides, and py4dggs is a guide rather than
an authority. Probes that read the engine's own members and entry points, on x86-64 Linux with
glibc, found the projection bit-identical to the engine's, on all three projections, once those
orders were ported: its cached geometry, the forward projection at 3,184 geographic stations and at
300,000 random points, and the inverse at the 654 of 1,831 planar stations and the 129,523 of
300,000 random points that lie on a face of the layout.

The aperture-3 topology is in the same position. None of the members it ports carries the mark,
and at eight sites where the library computes in another order the difference can change an
answer, a zone or a coordinate; those, and the sites beside them where it cannot, are performed
in the compiled order. A probe of the planar quantiser at 4,707,626 points, built on cell
boundaries, seams, both poles and far outside the layout at every resolution from 0 to 33, found
its answer the engine's at every point but eight, where this crate departs from the engine on
purpose (see Departures from the engine), and the centroids and vertices of the 297,986 zones it
named the engine's to the last bit.

So are the aperture-7 parents, children and sub-zone order: none of the functions they port
carries the mark. At one site the compiled order changes an answer. The offset from the vertex of
a zone that stands uppermost on the icosahedral net to its first sub-zone centroid is written in
the source as `k / (3.0 * p)`; the library holds `k / 3.0` folded as the double `fl(k * fl(1/3))`,
which for two of its six values of `k` is not the double nearest to `k/3`, and divides it by `p`.
The first centroid, from which the scanlines of the order set out, then differs from the source's
by a unit in the last place at some zones, and unit tests on the engine's own centroids fail when
the site is put back to the source's order. At every other site the compiled order is followed
all the same, each site with its address, although it was found to change no answer. The sites of
the walk across the interruptions of the layout, and those of the scanlines of the order, change
none on the arguments of real orders, whose origins lie inside the layout. Those of the searches
for the parents and the children changed no list in an ablation, the sites put back to the
source's order, over the parents and children of 62,876 zones of IGEO7. The three of the walk
that finds the index of a sub-zone are equal to the source by algebra, and part from it at no
number, at a NaN alone.

The local gate reads the BuildID of the `libdggal.so` the oracle links and fails if it is another.
A different build of DGGAL, a rebuilt 0.0.6 or a later release, may compile these functions
differently, and the orders would then have to be read again from that library.

## Departures from the engine

Where the engine answers input that has no meaning with an answer that looks like one, answers a
question without examining it, answers in contradiction of its own answers, or leaves out of an
answer what its own test admits, this crate departs from it on purpose; and it departs from the
order of the engine's vertices on the aperture-3 grids, so that a ring runs the same way on every
grid. Each departure is documented on the method that makes it. Where the suite re-checks the
engine's side of a departure on every run, so that the departure stays as described, the item
says so.

- **A NaN or infinite coordinate.** `zone_from_geo` returns `Error::NonFinite`, naming the
  latitude or the longitude. The engine answers one fixed pentagon at whatever resolution was
  asked and whatever the input, centred at about (58.4, -168.8) in the northern Pacific, on all
  six grids: that of base cell 1 on the aperture-7 grids, and that of root rhombus 0 on the
  aperture-3 grids. Re-checked by the suite. A finite coordinate out of range is not refused: a
  latitude beyond 90 degrees either way, or a longitude beyond 180, is projected as it stands, as
  the engine projects it, and answered with some zone, or the null zone, rather than an error, so
  a caller that needs such input refused must check it first. At 18,120 quantisations of such
  coordinates on the three aperture-7 grids, a sample the suite does not repeat, every answer was
  the engine's own.
- **A `GridConfig` with a NaN or infinite angle.** `Grid::new` returns `Error::NonFinite`,
  naming the field. The engine takes no orientation from a caller, since it fixes each grid's in
  its own source, so there is no answer of its to compare.
- **The text of a Z7 identifier finer than resolution 19.** `zone_from_text` refuses text longer
  than an identifier at resolution 19, the finest. The engine reads twenty digits, the whole
  packing, as a level-20 zone and gives it no geometry: the zero point as its centroid, six zero
  vertices and no neighbours, which a caller could not tell from a genuine cell at (0, 0).
  Re-checked by the suite. This crate still prints such an identifier when one is bound directly
  with `Grid::zone`, and that text does not read back, so a caller that binds raw identifiers
  should compare their resolution with `max_resolution`.
- **The text of a polar pentagon's sub-hexagon C or D**, such as `AA-0-C`, on the aperture-3
  grids. `zone_from_text` refuses it. The engine reads it, although no neighbour, parent, child or
  sub-zone of any zone is such a cell. Re-checked by the suite, which also checks that no relation
  of a polar pentagon reaches one.
- **The null zone on the aperture-3 grids** answers as it does on the aperture-7 grids: the zero
  point as its centroid, six zero vertices, no neighbours, parents or children, and resolution 0;
  and the sub-zone methods refuse it with `Error::InvalidZone`. The engine gives it a centroid
  that belongs to no cell, near (-21.2, 35.7) on ISEA3H and a little way off it on the other two,
  six vertices about that point, a level of 63, and a parent chain on which the engine divides by
  zero. Re-checked by the suite, which never asks the engine for the relations that crash it.
- **`sub_zone_index` with an identifier that has no cell behind it**, the null zone among them.
  This crate refuses it with `Error::InvalidZone`, as the other sub-zone methods do. On the
  aperture-7 grids the engine compares the two zones' levels before anything else, and answers 0
  for such an identifier within itself and -1 against a real zone; re-checked by the suite. On
  the aperture-3 grids, asked for the index of the null zone within a real zone, the engine never
  returns, and the suite does not ask it.
- **`sub_zone_index` for two different zones at one level.** This crate answers `None`, since at
  depth 0 a zone's only sub-zone is itself. The engine answers 0 whenever the two share a level,
  and so reports a zone's neighbour as its sub-zone at index 0. The suite re-checks that answer,
  on every grid, on neighbours of the zones it samples.
- **The first sub-zone at depth 0.** `first_sub_zone(z, 0)` and `sub_zone_at_index(z, 0, 0)`
  answer the zone itself on every grid, as `sub_zones(z, 0)` does. The engine's `getFirstSubZone`
  has no case for depth 0: it quantises a vertex of the zone at the zone's own level and so names
  a neighbour. On the aperture-3 grids it does so for about two thirds of zones (763 of the 1,176
  sampled on ISEA3H, 760 of 1,175 on IVEA3H and 757 of 1,175 on RTEA3H), and the suite finds its
  answer among the zone's neighbours. On the aperture-7 grids it does so at nearly every zone: of
  the 5,607, 5,579 and 5,574 zones compared on IGEO7, IVEA7H and RTEA7H it names a neighbour at
  5,594, 5,570 and 5,566, its null zone at 6, 5 and 4, the zone itself at 3, 1 and 1 and, at the
  two broken edges, another identifier that is no neighbour at 4, 3 and 3. At those edges the
  engine's own `getSubZones` at depth 0 holds another identifier than the zone at 10, 5 and 5 of
  the zones compared. Re-checked by the suite, which classifies the engine's answer at every
  zone and holds each class to its exact count.
- **The sub-zone at an index below an aperture-7 pentagon.** `sub_zone_at_index(z, d, i)` is entry
  `i` of `sub_zones(z, d)`, at every zone. The engine's `getSubZoneAtIndex` is not, below a
  pentagon at an odd depth: from the scanline after the pentagon's own it answers another zone
  than its `getSubZones` lists at that index, often one that is no sub-zone of the pentagon at
  all. For `00` at depth 1 it answers `030` and `034` at the indices 9 and 10, where its list
  holds `006` and `043`. Measured through DGGAL's own command-line tool `dgg`, against which
  this crate's tool is tested: of the 10,292 indices of the orders below the twelve pentagons of
  resolutions 0 and 1 and below three hexagons, at depths 1 to 3, the engine's answer by index is
  not the entry of its own list at 2,038 on each of the three grids, every one of them below a
  pentagon at depth 1 or 3. The fault is not spread evenly. A pentagon of resolution 0 has two
  such indices at depth 1 and 115 at depth 3, or one and 54 if its base cell is one of 6 to 10; a
  pentagon of resolution 1 has none at depth 1 and 97 at depth 3, or 53 for the same five base
  cells; and no index is at fault at depth 2 or below a hexagon. Re-checked on every run by the
  tool's tests, which hold that table exactly; the library's suites compare `sub_zone_at_index`
  with the engine's order and not with its search, at every index of 457 orders below pentagons
  at odd depths, among the others (see Verification).
- **The first sub-zone where the engine's call ends the process.** Below the pentagon of the
  south pole (base cell 11, every digit 0) at an odd resolution, at depth 2, the engine's
  `getFirstSubZone` ends the process that calls it, while its `getSubZones` of the same zone and
  depth answers. `first_sub_zone` answers the first entry of that list there, `11022` for `110`.
  The suite does not make the engine's call at such a pair: nine on each grid, at the
  resolutions 1 to 17, counted exactly, at each of which this crate's answer is compared with
  the first entry of the engine's list.
- **The index of a sub-zone on the aperture-7 grids is confirmed, found where the engine's test
  accepts the zone, and bounded.** `sub_zone_index` finds an index as the engine's `getSubZoneIndex`
  does, by the engine's test of whether the zone is a sub-zone and its walk of the order's
  centroids, and answers it where `sub_zone_at_index` at that index gives back the zone asked about.
  Where the test accepts the zone and the walk gives no such index, the zone is sought in the order
  itself, generated entry by entry and never built, and its first place there is the index, or
  `None` where the order does not hold it. Where the test refuses the zone the answer is `None` at
  once, as the engine's -1. The engine's walk fails at the two broken edges, among sub-zones of
  resolution 15 and finer. There it answers, for some entries of its own order that its test
  accepts, -1, and for some an index at which its order holds another zone: in the suites' samples,
  59, 67 and 63 entries of the first kind on IGEO7, IVEA7H and RTEA7H, and 403, 379 and 411 of the
  second, each of which has its place here. Entries that the test refuses have no index, here as in
  the engine: 64, 108 and 70 in the main samples and 2,988 among the identifiers built by text, all
  of resolutions 15 to 19. And 18, 16 and 22 times the engine answers an index for a zone that the
  order does not hold, which has none here. Where an order names a zone twice, the index is the
  place that the walk finds, if it finds one, and the first place otherwise. Re-checked by the
  suite, which at every such entry compares the engine's answer and its test, and finds the zone at
  the place answered. The walk costs as the index does, and the search in the order more, so the
  request is refused with `Error::TooManySubZones` where the order is longer than the limit on a
  list, which the engine does not do: it walks an order of any length.
- **A walk of the aperture-7 sub-zone generator that does not end is given up.** The engine
  steps from one centroid of an order to the next across the interruptions of its layout by a
  function whose loop has no bound, and from some arguments it does not return. This crate gives
  such a walk up after 64 turns of that loop, and the entry it was to reach is then the null
  zone. No call made for a real order is known to reach the bound: over 1,022,074,348 calls made
  on IGEO7 for the orders of a census (every zone of resolutions 0 to 3 and a sample of each
  resolution to 18, at depths 1 to 3 and, for the pentagons and one zone in forty, to 5), and
  over a thousand million on each of the other two grids, no walk took more than four turns.
  Arguments that reach the bound exist all the same: of 186,104 invented ones the engine did not
  return from 1,593, which are exactly those at which this crate gives up, and one such call the
  engine's own source sets aside as endless. A unit test gives that call to this crate alone,
  since the engine does not return from it.
- **The centroid parent on the aperture-7 grids.** `centroid_parent` answers the first of the
  zone's parents that is itself a centroid child, and `None` where none is: the meaning that the
  engine's source gives it, and the one it has on the aperture-3 grids. The engine's own
  `getZoneCentroidParent` answers its null zone for every zone of these three grids, because it
  asks for the zone's parents through the class of the grid, which reads the identifier it is
  handed in another indexing than the one it is written in. Re-checked by the suite, which finds
  the engine's null zone at every zone it compares, and compares this crate's answer with the
  first of the engine's own parents that the engine calls a centroid child.
- **Three kinds of entry dropped from the aperture-7 parents and children.** `parents` and
  `children` list the engine's own `getZoneParents` and `getZoneChildren`, entry for entry and in
  its order, less its null zone, an identifier the engine cannot read back, and a repeat of an
  entry already listed. All three occur only at the two broken edges (see Known limitations),
  and each dropped entry is re-checked live as one of the three kinds. Among the zones that the
  main sample's positions give, 104, 76 and 76 null zones and 24, 8 and 8 repeats are dropped
  from lists of children on IGEO7, IVEA7H and RTEA7H, and nothing from any list of parents.
  Among the 2,276 identifiers built by text beside a broken edge, the same on the three grids,
  459 null zones, 1,002 unreadable identifiers and 207 repeats are dropped from lists of
  children, and 91 unreadable identifiers and one repeat from lists of parents. A hexagon can
  therefore list fewer than thirteen children, and a pentagon fewer than eleven, although the
  engine's own list has thirteen entries, or eleven, at every zone compared below the finest
  resolution. Where the first of the engine's parents that it calls a centroid child is an
  identifier it cannot read back, that entry is dropped, and `centroid_parent` is the first
  centroid child among the parents that remain, or `None`: 91 of those 2,276 identifiers, and
  none of the zones from positions.
- **The index of a sub-zone at resolution 33**, on the aperture-3 grids. `sub_zone_index` answers
  it, where the engine's `getSubZoneIndex` answers -1 for most such sub-zones: for 7,300 of the
  8,725 compared on ISEA3H, 7,277 of 8,725 on IVEA3H and 7,270 of 8,711 on RTEA3H, some 83 to 84
  per cent, the sub-zones at resolution 33 of every zone of the suite's samples at resolutions 31
  and 32. The engine first asks its own `zoneHasSubZone`, which fails there on the needle's
  level-34 centroid child, whose index overflows its field. Re-checked by the suite: wherever the
  engine answers -1 its `zoneHasSubZone` is false, and wherever it answers an index, the index is
  this crate's.
- **The aperture-3 sub-zone order at some zones on the edges of the rhombi.** There the engine's
  order is at fault: a rounding of one unit in the last place decides whether a point of the
  generator's scan is carried across an interruption of the layout that it should only reach, and
  where it is carried wrongly the engine's order names a zone more than once, names a zone that
  is not a descendant of the zone, or omits a descendant. `sub_zones` there answers each
  descendant once, as many as the engine's own count, equal to `children` at depth 1, in the
  generator's own scanline order; `first_sub_zone`, `sub_zone_at_index` and `sub_zone_index`
  follow the same corrected order, so a zone has one, unique index. Every other order is the
  engine's, down to the bits of the generator's centroids. Measured over the edges of every
  rhombus, the fault falls on them alone: at resolutions 13, 19 and 21, on some sub-hexes along
  the top edges of the northern rhombi 2, 4 and 8 and the left edges of the southern rhombi 3, 5
  and 9, at every depth measured (1 to 8); on the left-edge hexagons of rhombus 5 at resolutions
  16, 20, 22, 24 and 26, and of rhombi 1 and 9 at 24, 26 and 28, at some even depths; and on the
  southern polar pentagon at resolutions 11 and 25, at depths 3, 5 and 7. The correction follows
  the mechanism, not this list. Re-checked by the suite, which finds the engine's own order
  failing its own tests and this crate's passing them, at eight corrected orders on each grid, the
  same eight, and counts them exactly.
- **The order of the vertices on the aperture-3 grids.** Every ring runs anticlockwise as seen
  from outside the sphere, on every grid (bar the aperture-7 exceptions stated under "What a
  caller must know": rings with every vertex on a pole, and rings that hold the point (0, 0)). The aperture-7 engine gives its rings so. The aperture-3
  engine gives them clockwise, save at the two polar pentagons, roots 10 and 11 (the zones whose
  text has `A` or `B` for its second character, such as `AA-0-A` and `CB-0-B`), whose rings run
  anticlockwise at every resolution: so it was found at every zone of resolutions 0 to 7, at the
  twelve pentagons of every resolution from 0 to 33, and at some 102,000 random zones of those
  resolutions on each of the three grids. For every other zone this crate gives the engine's ring
  `[v0, v1, ..., v(n-1)]` as `[v0, v(n-1), ..., v1]`: it starts at the engine's own first vertex
  and turns the other way. The polar pentagons keep the engine's order. Every coordinate is the
  engine's own, to the last bit, and the sub-zones are generated from the engine's own order, as
  the engine generates them. Re-checked by the suite, which finds the engine's ring clockwise at
  every zone it compares but the polar pentagons, and anticlockwise at those.
- **Identifiers the engine cannot read back.** On the aperture-7 grids, at resolutions 16 and 18,
  in a narrow strip along two icosahedron edges (see Known limitations), the engine answers with
  an identifier that its own text parser cannot read back, one fixed identifier per resolution
  for points tens of degrees apart; this crate answers the null zone there, which is what the
  engine answers at the odd resolutions either side. At a few points of the strip's margin at
  resolution 19 the engine does the same and this crate again answers the null zone; those points
  were found by a separate sweep, and no sample of the suite meets them. On the aperture-3 grids,
  at every resolution from 2 to 21, the engine does the same along two seams of its layout: along
  the icosahedron edge that runs over the north pole, on the meridians 11.2 E and 168.8 W from the
  pole down to 58.4 N, and in the southern hemisphere from the equator to 58.4 S, between the
  meridians 137 W and 169 W, along the edge that ends at the southern polar vertex. There it
  answers an identifier with the root of a polar pentagon and a non-zero index, which its parser
  refuses and to which it gives the geometry of a polar pentagon tens of degrees away; this crate
  answers the null zone. The bands reach some 2e-5 degrees from the edges on each grid. On every
  grid the suite re-checks, at every such point it meets, that the engine cannot read its answer
  back. An aperture-7 sub-zone order holds the null zone likewise, at the engine's own place,
  where the engine's order holds such an identifier: 275, 260 and 260 entries of the orders the
  suite compares on IGEO7, IVEA7H and RTEA7H, among sub-zones of resolutions 16 and 18 and,
  rarely, 19, each re-checked live.
- **Four kinds of entry dropped from an aperture-7 neighbour list.** `neighbors` lists the
  engine's own `getZoneNeighbors`, entry for entry and in its order, less its null zone, since the
  null zone is no neighbour of anything; the zone itself; a repeat of an entry already listed; and
  an identifier the engine cannot read back, which this crate never hands out, as elsewhere in
  this section. All four occur only inside the broken seams (see Known limitations), where the
  engine's own lists are not consistent with themselves, and each dropped entry is re-checked live
  as one of the four kinds. Over a census of 267,435, 267,537 and 267,612 engine-readable zones on
  IGEO7, IVEA7H and RTEA7H, 848, 838 and 851 zones per grid have a drop,
  and the entries dropped are 806, 792 and 806 null zones, 106, 105 and 110 entries naming the zone
  itself, 748, 747 and 745 repeats of another entry, and 25, 23 and 21 identifiers the engine
  cannot read back. A hexagon can therefore list fewer than six neighbours, and a pentagon fewer
  than five.
- **The children of a zone at resolution 33**, the finest on the aperture-3 grids. This crate
  answers none, as it answers none at resolution 19 on the aperture-7 grids. The engine answers
  seven identifiers of level 34, or six under a pentagon, which its own parser cannot read back
  and on which its parent chain divides by zero. Re-checked by the suite.
- **A sub-zone depth that would reach past resolution 33**, the finest on the aperture-3 grids.
  `count_sub_zones`, `first_sub_zone`, `sub_zones` and `sub_zone_at_index` refuse it with
  `Error::InvalidZone`. The engine answers regardless: its `getSubZonesCount` gives its ordinary
  formula for a level past its own maximum, its `getFirstSubZone` quantises a vertex at that
  level and names one of the level-34 identifiers of the previous bullet, which its own parser
  cannot read back, and only its `getSubZones` answers nothing, since it checks the level before
  generating anything. Not re-checked live: the suite never asks a depth that would reach past
  resolution 33. On the aperture-7 grids a depth that would reach past resolution 19 is refused
  in the same way, and the suite does not ask the engine there either: its `getSubZones` then
  answers a null array, which the binding the suite links does not survive.
- **The null zone's resolution** is reported as 0 on every grid, so that the null zone reads
  alike on all six, where DGGAL reports -1 on the aperture-7 grids and, reading the sentinel's
  fields as they stand, 63 on the aperture-3 grids. A caller who needs to tell a resolution-0 zone
  from the null zone compares the identifier with `ZoneId::NULL`, or its text with `NULL_TEXT`.
  The suite pins both values, and checks the null zone against the engine in everything else.
- **A zone without geometry has no ring, no extent and no area.** `refined_vertices` answers the
  empty ring, `extent` and `area` answer `None`, for the null zone, for a Z7 identifier of twenty
  digits, and for an identifier that no cell of the grid has behind it: the very identifiers for
  which `vertices` gives six zero points, and on the aperture-3 grids they include those of a polar
  root with a non-zero index, such as `BA-2-A`. An aperture-7 identifier that has an address but
  does not read back as text keeps its ordinary ring, extent and area. The engine answers its null
  zone with an empty ring, an extent cleared to `1.7976931348623155e308` and its negation, and an
  area of `+inf`; an identifier it does not draw it gives the area of its level's formula, a finite
  number, or `+inf` where the divisor of that formula is zero. Re-checked by the suite at every
  such identifier it meets, which classifies the engine's answer and counts the classes exactly:
  in the main sample and the fixed list, the null zone once on every grid, undrawn identifiers
  with a finite area 5 on each aperture-7 grid and 32 on each aperture-3 grid, and one on each
  aperture-3 grid whose area is `+inf`. A further class, which the engine does draw, is the
  subject of the first of Known limitations.
- **The edge refinement is bounded.** The engine accepts any `int` as the number of parts into
  which an edge is divided; `refined_vertices` takes a `u32` and refuses, with
  `Error::EdgeRefinementOutOfRange`, one above `rs4dggs::grid::max_edge_refinement()`. That limit
  is `rs4dggs::grid::MAX_EDGE_REFINEMENT`, 100,000, the greatest refinement at which the boundary
  was compared with the engine's: a hexagon then has 600,000 points, and one with an edge beside
  a pole, which the engine divides twenty times as finely, 2,500,000. Beyond some 107 million
  the engine's own count of divisions wraps. The environment variable
  `RS4DGGS_MAX_EDGE_REFINEMENT` lowers the limit and never raises it, for an operator who must
  bound what a caller can ask. It is read once, at the first use, and kept for the life of the
  process; its value is read as a `u32`, once surrounding whitespace is trimmed; a number above
  the ceiling gives the ceiling; and a value that is empty, not a number, negative or too large
  for a `u32` is ignored in silence, since a library prints nothing. On a target without an
  environment, `wasm32-unknown-unknown` for one, the variable is absent. A refinement of 0, the
  engine's own choice, is always accepted whatever the limit, so a limit below 20 does not bound
  the ring it gives: the engine's choice divides an ordinary edge into at most twenty parts, and
  an edge beside a pole into twenty times as many. The argument is examined before the zone, so
  the refusal does not depend on the identifier. The refusals at 100,001, 2^20 and `u32::MAX` are
  checked on `Grid`, `AnyGrid` and `Zone`. A contributor who runs the suites must leave the
  variable unset, and they refuse to run otherwise.
- **The levels asked of the grid's facts are bounded.** `count_zones`, `ref_zone_area` and
  `meters_per_sub_zone` take a level from 0 to `max_resolution`, and the sum of level and depth
  for the last, and refuse another with `Error::ResolutionOutOfRange`; the engine goes on past
  its finest level from tables that run out, with counts that wrap and then settle at 2.
  `level_from_ref_zone_area` and `level_from_meters_per_sub_zone`, which answer a level for an
  area or a size, are not bounded: they mirror the engine for every input, and so answer up to
  `max_resolution + depth_64k`, 25 on the aperture-7 grids and 43 on the aperture-3 grids, which
  is beyond the finest resolution; an area of 0 or NaN answers that greatest level, a negative or
  infinite one answers 0, and the answer is not always the inverse of `ref_zone_area`: on ISEA3H
  the reference area of level 9 answers 10. A caller that needs a level the grid has takes the
  lesser of the answer and `max_resolution`. Both are checked against the engine over a census
  of 801,836 cases on each aperture-7 grid and 802,844 on each aperture-3 grid (see Verification), and the engine's own answers at
  levels this crate refuses are asked too, so that nothing in the arithmetic behind the bound
  goes untested.
- **Refusals name their cause, at bounded length.** `zone_from_text` returns
  `Error::InvalidZone` for text that names no zone, where the engine answers its null zone for
  any text it cannot parse (`RI7H_Z7.ec:418`, `RI3H.ec:908-931`). A rejected text identifier is
  quoted in the message to its first thirty-two characters, followed by its length, and an
  unknown grid name is kept to its first thirty-two characters, so that no error carries caller
  input at unbounded length. The suite does not compare these.
- **The zones of a box are those whose extent meets it, which on the aperture-3 grids is more
  than the engine's `listZones` returns.** `zones_in_box` yields a zone of the level if and only
  if its extent meets the box by the engine's own test of two extents. On the aperture-3 grids
  the engine looks for zones about points sampled within the box, in the plane of its layout,
  and leaves out a zone that belongs to another rhombus than the one in which those points fall,
  among them the zone that holds the whole of a small box, and every zone of a small cap about a
  pole. The box about Lisbon, from 38 to 39 degrees north and from 10 to 9 degrees west, is one:
  at level 2 the engine returns no zone where `B4-2-A` holds the whole box, and at levels 3 and 4
  one zone, `B2-5-C` and `C2-23-A`, where the box meets two (`B2-5-C` and `B4-2-B`, `C2-23-A`
  and `C4-6-A`); at levels 0 and 1, and from level 5 to level 10, the two answers are one, on
  ISEA3H, IVEA3H and RTEA3H alike. The cap from 89.5 degrees north to the pole is answered with
  no zone at level 6, where two zones cover it. The engine returns no zone outside the rule. Its
  zones come in its order, and this crate's answer holds them in that order with those it
  leaves out among them, so that it departs from the engine only by addition. Re-checked by the
  suite at every box where the two answers differ, as Verification describes: 67, 66 and 59 of
  494 boxes on each aperture-3 grid. On the aperture-7 grids the engine descends from the twelve
  base cells through the children of each zone and returns the zones of the rule, save, rarely,
  a zone whose extent reaches beyond the extents of both its parents, which it does not reach.
  That was measured at six boxes in 38,880, of the same nine kinds, at levels 1 to 5; no box of
  the suites meets one, so that on the aperture-7 grids this departure is stated here and not
  re-checked by them.
- **A bounding box that is none is refused, before any work.** `zones_in_box` and
  `estimate_zones_in_box` return `Error::NonFinite` for a coordinate that is not finite,
  `Error::LatitudeOutOfRange` for a latitude beyond a pole, `Error::LongitudeOutOfRange` for a
  longitude outside -180 to 180 degrees, each naming the coordinate at fault, and
  `Error::InvertedBox` for a south above the north.
  The engine answers such a box with something, with an empty list or with a null array; on the
  aperture-3 grids it does not return from a box with an infinite latitude; and on both
  apertures it overflows its stack on a box over the antimeridian with a longitude beyond 180
  degrees. A west above the east is no error, here as in the engine: the box runs eastwards over
  the antimeridian. A box one of whose sides is shorter than 1e-11 radian without being nought,
  the other being longer, is answered by the rule, a box from which the engine's aperture-3
  search does not return. The suites never hand the engine such boxes, nor those on which its
  aperture-3 search spends seconds or gigabytes (among them a box across one of the edges of the
  icosahedron at which the plane of its layout is cut, and a box of more than 20,000 zones), nor,
  on the aperture-7 grids, a box in which it finds no zone, which it answers with a null array
  that the binding the suites link does not survive; they count those they keep from it.
- **The compaction takes the zones of one level, and only zones.** `compact_zones` refuses a
  slice that holds zones of more than one level with `Error::MixedLevels`, which carries the
  level of the first zone and that of the first zone of another level. The engine takes such a
  set and loses zones of it: a zone that has no coarser zone to give way to is dropped as soon as
  the set holds a finer one, so that `00 0064` compacts to `0064` on IGEO7, and an answer of its
  own, which is of several levels, is not a set that it compacts to itself. The null zone, and on
  the aperture-7 grids a Z7 identifier of twenty digits, are refused with `Error::InvalidZone`,
  where the engine leaves both out of the set, as is every identifier that names no zone of the
  grid. And on the aperture-7 grids, where the engine's answer holds an identifier that the
  engine cannot itself read back, this crate's answer lacks it, as its lists of parents and
  children do: that was met only for identifiers of resolutions 17 and 19 built from digits
  beside a broken edge (see Known limitations), which are not the zone found at their own
  centroid, and for no set made of zones found from positions. Each is checked against the
  engine's own answer by the tests of the compaction.

## What a caller must know

The measurements quoted in this section were made on x86-64 Linux with glibc, as were those
above.

- **A disk is an ordered walk, not a set.** `disk(k)` returns a lazy iterator, `Disk`. The centre
  comes first, as ring 0; each further ring is built from the zones of the ring before, in the
  order they were yielded, each zone's neighbours being taken in the order `neighbors` returns
  them and each zone yielded on its first arrival only. The order is the same on every call.
  Nothing is computed until the disk is consumed, so a caller pays for what it takes and for
  nothing more, and `disk(u32::MAX)` is safe to construct. The walk remembers every zone it has
  yielded, for the reason given in the next point, so the memory it holds grows with what is
  taken: of the order of 3k² identifiers for a disk of radius k consumed whole. `Disk::rings`
  gives the same walk a ring at a time.
- **The neighbour relation is not symmetric everywhere.** On the aperture-7 grids, along the two
  icosahedron edges where the engine is not consistent with itself (see Known limitations), at
  resolutions 15 to 19, a zone may list a neighbour that does not list it back. This crate's lists
  there are the engine's own, less the four kinds documented under Departures from the engine, so
  these one-way pairs are the engine's own one-way pairs, not a separate phenomenon to weigh
  against it. Measured over some 949,000 to 951,000 zones sampled on each aperture-7 grid (every
  zone of resolutions 0 to 4, the pentagons and both poles, the two broken edges and every other
  icosahedron edge, the two-step neighbourhoods of those zones at resolutions 13 and up, and
  some 66,500 zones from uniform random points at every resolution), there are 118,920 such
  one-way pairs on IGEO7, 117,659 on IVEA7H and 119,448 on RTEA7H, most of them nearer the edge
  from base cell 0 to 1 than the one from 1 to 6, every one within 1.3e-5 degrees of one of the
  two edges, and none beyond 2e-4 degrees or anywhere else in that sample.
  The 2,860 digit paths deep under the two poles that Verification names reach further: there,
  among the 1,673 one-way pairs per grid, a pair may lie up to 3.57e-5 degrees from the nearer
  edge. On the aperture-3 grids the neighbour lists are the engine's own, in its order.
- **An identifier carries no grid.** The three aperture-7 grids share the Z7 packing, so an
  identifier, or its text, taken from one of them is accepted by the others without complaint and
  names another cell there. Measured over 60,000 identifiers drawn from all three, the centroids
  that one identifier names on two grids lie some 7 to 14 km apart on average, according to the
  pair, and up to 56 km; the pentagons coincide on all three. The three aperture-3 grids likewise
  share the I3H packing: the zone that contains Lisbon at resolution 5 is `C2-23-C` on each, but
  at resolution 15 each names a different zone for it. Nothing in the identifier reveals the
  mistake, so a caller that keeps identifiers of more than one grid must keep each grid's name
  beside them.
- **Vertices run anticlockwise** as seen from outside the sphere, on every grid, with no closing
  vertex appended, save the few aperture-7 rings that hold the point (0, 0), described below. On the
  aperture-7 grids that is the engine's order: measured on some 35,700 zones of each, every ring
  that encloses any area runs anticlockwise, and every ring stands in the engine's order. On the
  aperture-3 grids the ring is the engine's reversed from its first vertex, save at the two polar
  pentagons, where it is the engine's in its order (see Departures from the engine), so a caller
  that compares rings with DGGAL's there must reorder one of them. Each longitude is given on its
  own, between -180 and 180, so a ring that straddles the antimeridian is not unwrapped, and a ring
  about a pole runs through every longitude; a caller that draws or intersects rings on a plane must
  split or shift them itself. A ring can also name one point twice. Near the poles at resolutions 18
  and 19 of the aperture-7 grids, the engine's pole snap, which this crate reproduces, puts some
  consecutive vertices on the pole, and at 19 some cells have every vertex there. Beside the
  pentagons of base cells 0 and 1 at resolutions 16 and 18, a ring of six can end on the point it
  began at, and a test of the first and last vertices for equality would take it for a closed
  pentagon. Under the pentagon of base cell 0 at resolutions 15, 17 and 19 the engine puts the point
  (0, 0) among the vertices of some cells, one or three of them, where a planar vertex lies on the
  diagonal of its layout and the inverse projection answers the zero point; the rest of such a ring
  runs anticlockwise. The suite meets sixteen such cells on each aperture-7 grid. Each of these is
  the engine's own answer.
- **The refined boundary and the extent.** `refined_vertices` gives the ring the engine draws, its
  own sequence and every coordinate its own to the last bit. It runs anticlockwise as seen from
  outside the sphere and it is never closed: its last point does not repeat its first, and none
  of the rings measured does. On the aperture-7 grids it begins at the first vertex of
  `vertices`; on the aperture-3 grids at the second, the first coming last, and the two polar
  pentagons alone begin elsewhere, the northern at its fifth vertex and the southern at its
  first. **Away from the poles its longitudes are continuous**, and so lie beyond 180 degrees
  either way where a zone passes the antimeridian, as the engine gives them, so that the ring does
  not jump by a turn from one point to the next: 94, 95 and 95 of the zones met on IGEO7, IVEA7H and RTEA7H, and
  119, 117 and 117 on ISEA3H, IVEA3H and RTEA3H. A caller that needs each longitude between -180
  and 180 wraps them itself. The extent, `Extent { ll, ur }` in degrees, has both its
  longitudes in that interval, and **over the antimeridian the western end is the greater**,
  `ll.lon > ur.lon`: the same zones, one for one, so that it reads as an extent which runs east
  from `ll.lon` through 180 degrees to `ur.lon` and is not an interval. A zone that touches a
  pole has `ur.lat` of exactly 90, or `ll.lat` of exactly -90, and an extent about 180 degrees
  wide, which the engine does not widen to the whole circle: 60, 61 and 60 zones on the
  aperture-7 grids and 67 on each aperture-3 grid. The extent is that of the ring at the
  engine's own refinement, and not of the vertices alone.
  The count of points is not always the count of edges times the refinement: an edge across an
  interruption of the layout gains a point there from a refinement of 2, an edge beside a pole
  has twenty times as many parts, two consecutive points within 1e-11 radians of one another
  are given once, and a point whose inverse projection fails is left out, so that a refinement
  of 1 need not give the vertices of `vertices`. Beside the poles the engine's rings degenerate,
  and these are its rings, mirrored: at resolution 19 of the aperture-7 grids a ring of a zone
  centred on a pole is one point, the pole, 3 such rings on each grid in the sample; on the
  aperture-3 grids a ring holds two or more points at a pole, each with a latitude of exactly 90
  degrees (38, 36 and 38 rings in the sample), and some name one point at two places (14, 12 and
  11). A vertex on a pole carries no meaningful longitude, and the ring may jump by half a turn
  across it: on ISEA3H, at the engine's own refinement, `A8-0-C` runs from (89.946, 191.2) to
  (90, 11.2) and on to (89.946, 11.2), and `F0-79-A` has five points on the pole, at -168.8,
  between sides at -78.8 and -258.8; on IGEO7 the zone `013` has (90, 191.2) between longitudes of
  243.6 and 63.6. A caller that reads the ring as planar coordinates replaces each run of pole
  vertices by two points at that pole, with the longitudes of the neighbouring vertices before
  and after the run. A ring of fewer than three points, as at resolution 19 above, is not a
  polygon. On the rings examined near the poles the pole was always a vertex of the ring, never
  enclosed by it. Four rings on each aperture-3 grid, at
  resolution 33 and centred on or beside a pole, enclose no area. Some rings cross themselves.
  On the sphere, with the pole one point, three on each aperture-7 grid do and none on the
  aperture-3 grids; read as planar longitudes and latitudes, as a drawing or a
  GeoJSON reader reads them, 3, 3 and 6 rings on the aperture-7 grids do, beside the poles at
  resolutions 15 to 19, and 17, 12 and 19 on the aperture-3 grids, from resolution 10. A caller
  that draws or intersects rings on a plane must expect them. All these counts are of the main
  sample and the fixed list of the suites, held exactly there, and are not counts of every zone.
- **The area is nominal.** `area` is the engine's `getZoneArea`: every hexagon of a resolution is
  the same share of the area of the whole ellipsoid, and a pentagon five sixths of one, whatever
  the projection and however the zone is drawn. It is not an area measured from the boundary.
  `ref_zone_area` is the whole ellipsoid divided by the count of zones, which is neither the
  hexagon's nor the pentagon's.
- **The null zone is an answer, not an error.** `zone_from_geo` succeeds with the null zone,
  whose identifier is `ZoneId::NULL`, where a point falls in no cell. On the aperture-7 grids the
  engine answers it in a narrow strip along the same two edges, at the odd resolutions from 15. On
  the aperture-3 grids it answers it near both poles at resolutions 0 and 1, and at the coarsest
  resolutions in a narrow band beside each of the ten icosahedron edges that end at the vertices of
  the two polar pentagons, 58.4 N 11.2 E and 58.4 S 168.8 W: a probe met it there at resolutions
  0 to 7, reaching between 1e-6 and 2e-6 radians from the edge at resolutions 0 and 1, and less
  than 1e-7 at 6 and 7. This crate answers as the engine does there, at every such point the
  suites sample. On every grid it answers the null zone too at the points described under
  Departures from the engine, where the engine answers an identifier it cannot read back. `Err`
  is kept for what the caller got wrong.
- **The aperture-7 hierarchy is not the digits of the identifier.** The digits name a path, and
  not a containment. The first of a zone's parents is the zone that its identifier names with
  the last digit dropped, and the zones that the identifier names with one digit more are among
  its children: the first seven of the thirteen, or the first six of the eleven, taken as a set,
  for the engine lists them in another order than the digits'. A zone whose last digit is not 0
  has a second parent, and every zone has six children, or five, that the digits make children
  of its neighbours. `Zone::parent` followed upwards therefore drops one digit at a time and
  gives one of the lines by which a coarser zone is an ancestor; a caller that wants every zone
  of which a zone is a part follows `parents`. That is so away from the two broken edges; what
  the engine's lists do at them is said, once, under Known limitations.
- **The aperture-7 sub-zone order has no relation to the digits.** The engine generates the
  sub-zones as centroids, in scanlines across the zone that set out beside the vertex of it that
  stands uppermost on the icosahedral net: at depth 1 a hexagon has scanlines of 1, 4, 3, 4 and 1
  zones. The centroid descendant, which the identifier names with as many zeros more as the
  depth, stands in the middle of the order and never at its head: at index 6 of 13, 27 of 55 and
  189 of 379. `count_sub_zones` exceeds `7^d`, since the order holds the zones that lie across
  the zone's boundary as well as those within it: 13, 55, 379 and 2,449 at depths 1 to 4 below a
  hexagon, and 11, 46, 316 and 2,041 below a pentagon. At depth 1 the sub-zones are the
  children, in another order. From depth 3 a zone that the digits name below a zone need not be
  a sub-zone of it at all, and the order holds zones that the digits name below its neighbours;
  a caller that wants the digit descendants builds them from the text, and one that wants the
  zones of a finer resolution that cover a zone asks for its sub-zones. An entry of an order may
  be `ZoneId::NULL`, at the two broken edges alone (see Known limitations): it is counted by
  `count_sub_zones`, it has no geometry, every method that validates an identifier refuses it,
  and a caller who walks an order skips it.
- **Sub-zones cost what they list.** `count_sub_zones` answers at any depth. `sub_zones` builds
  the whole list and refuses one longer than `MAX_MATERIALISED_SUB_ZONES`, four million, with
  `Error::TooManySubZones`: below a hexagon that is any depth beyond 13 on the aperture-3 grids,
  and any depth beyond 7 on the aperture-7 grids, where an order has 825,259 entries at depth 7
  and 5,767,201 at depth 8, or 687,716 and 4,806,001 below a pentagon. `sub_zone_at_index` pages
  through an order without building it, at any depth down to the finest resolution, and that
  limit does not bind it. On the aperture-3 grids it takes a number of steps of the order of
  3^(d/2) at depth d; a single call at depth 33 below a zone of resolution 0 takes under a tenth
  of a second in a release build. On the aperture-7 grids **it costs by the scanlines it passes,
  and so not the same at every index**: the search reckons the width of each scanline before the
  one that holds the index, and then generates one sub-zone. Measured in a release build, below
  a hexagon at depth 14, whose order has 678,223,896,391 entries, it takes some 6 microseconds
  at index 1, a millisecond at the middle and two milliseconds at the last index; and below a
  zone of resolution 0 at depth 19, the deepest order there is, between 0.3 and 0.4 seconds at
  the last index, as the engine's own search does. The cost is bounded by the number of
  scanlines, some 134 million at most, and never by the number of sub-zones. `sub_zone_index` is
  bound by the limit on a list, on every grid, and above it is refused with that error rather
  than searched for. On the aperture-3 grids it searches the list that `sub_zones` builds. On
  the aperture-7 grids it builds no list: it walks the centroids of the order from its head, so
  that its cost grows with the index, some 12 milliseconds to the middle of the 825,259
  sub-zones of a hexagon at depth 7 in a release build and twice that to the last, and
  microseconds for a zone that is no sub-zone. Where the engine's test accepts a zone that the
  walk does not find, the order itself is searched, at some 3.7 microseconds an entry in a release
  build, since each entry is quantised: up to 3 seconds through an order at depth 7, the longest
  within the limit. The limit in force is
  `rs4dggs::grid::max_materialised_sub_zones()`, which the error carries: the environment
  variable `RS4DGGS_MAX_MATERIALISED_SUB_ZONES` lowers it from four million and never raises it,
  for an operator who must bound what a caller can ask, on every grid. It is read as
  `RS4DGGS_MAX_EDGE_REFINEMENT` is: once, at the first use, and kept for the life of the
  process; as a `u64`, once surrounding whitespace is trimmed; a number above the ceiling gives
  the ceiling; and a value that is empty, not a number, negative or too large for a `u64` is
  ignored in silence. `count_sub_zones`, `first_sub_zone` and `sub_zone_at_index` are not bound
  by it, and a depth of 0, whose list is the zone itself, is answered whatever the limit. A
  contributor who runs the tests must leave the variable unset: the suites refuse to run
  under a lowered limit, and the unit tests that build lists fail under one (19 of 253 at a limit
  of 20, and 7 at 1,000, when it was tried). The local gate, `scripts/check.sh`, refuses to run
  while this variable or `RS4DGGS_MAX_EDGE_REFINEMENT` is set, and says which.
- **The zones of a box are those whose extent meets it, and the extent is a bounding box.**
  `zones_in_box` yields a zone where its extent, the least and greatest latitude and the
  westernmost and easternmost longitude of its boundary, meets the box by more than 2e-15 radian
  in each direction; two extents that only touch do not meet. So the answer covers a little more
  than the box: a zone whose extent meets the box at a corner is yielded although its hexagon
  does not touch the box. A west above the east is a box that runs eastwards over the
  antimeridian, so that a west of 10 and an east of 5 make a box 355 degrees wide (663 zones at
  level 6 of ISEA3H). A box may have no area: a segment yields the zones whose extents it
  strictly crosses, and a point those whose extents strictly hold it, from one to three away
  from the lines below, or none. **A box of no height exactly on a pole, a point on a pole among
  them, and a box of no width exactly on the antimeridian (a west and an east both of 180
  degrees, or both of -180) yield no zone at all**, on every grid and at every level, since no
  extent reaches beyond a pole and the pieces of an extent cut at the antimeridian end there: a
  segment along the antimeridian crosses no extent strictly, whereas one a thousandth of a degree
  from it yields zones. **A point exactly on the equator, or on one of the four meridians named
  below at an odd level of the aperture-3 grids, may yield none**, and a segment along such a
  line may lack the zone that holds one of its points, the extents of the zones on either side
  ending at the line or short of it: on the equator a tenth to a
  quarter of the points tried yielded none, on the aperture-3 grids at every level tried and on
  the aperture-7 grids at the even levels; on those meridians a third of the points tried at the
  odd levels of the aperture-3 grids, where the engine's extent of the zone that holds the point
  ends up to 2.6e-10 radian before the meridian. `zone_from_geo` is the call for the zone at a
  point. The whole world as a box gives every zone of the level, which `zones` gives more
  cheaply.
- **The order, and what a cursor is.** `zones` and `zones_in_box` yield the zones in the order
  of DGGAL's `listZones`: the zones of a lattice of ten rhombi and two polar roots, rhombus after
  rhombus, row after row, cell after cell, and the two polar zones last. It is the same on every
  call. On the aperture-3 grids it is the ascending order of the identifiers. On the aperture-7
  grids it is the order of the engine's internal identifier, and not that of the Z7
  identifiers, which no sorting of them gives: level 0 comes as `01 06 02 07 03 08 04 09 05 10
  00 11`. A cursor is one zone identifier, the last zone of the page before: `after(zone)` gives
  the zones of the sequence that follow it, and the pages so taken are the whole answer, each
  zone once. The zone need not be in the answer, only of the level; an iterator already
  consumed is entered anew, and after the last zone comes nothing. `after` refuses with
  `Error::InvalidZone` an identifier that the enumeration of that level does not itself yield:
  a zone of another level, the null zone, a Z7 identifier with a digit written after its last,
  or an identifier that names no cell of the lattice, such as the child that a pentagon does not
  have. A level beyond `max_box_level` is refused by both enumerations with
  `Error::ResolutionOutOfRange`, whose `max` is that limit.
- **What a step of the zones of a box costs.** Nothing is computed until the iterator is
  consumed, and no list is held. A step never costs by its place in the answer, so that a page
  far into an answer costs what the first page costs, nor by the area of the box. It costs by the
  length, in cells, of the stretch of the rim of the box that runs along a row of a rhombus of
  the lattice, which is a few cells where the rim crosses the rows and grows as the square root
  of the side of the lattice where the rim turns back. Measured on an optimised build, a cell of
  such a stretch costs from 7 to 47 microseconds: 7 to 9 on the aperture-3 grids to level 26,
  some 18 at level 32 and 35 at level 33, and 11 to 21 on the aperture-7 grids at the even
  levels and 40 to 47 at the odd ones; where the rim turns back, one step was measured at up to
  seven tenths of a second at the finest level of the aperture-3 grids (0.10, 0.23 and 0.68
  seconds at levels 31, 32 and 33). **The worst case
  is a side of the box that lies on one of the meridians along which the rhombi are joined**, to
  within some two cells, for then the stretch is the whole of that side. On the six grids of the
  registry there are four: 11.2 degrees east and 168.8 degrees west, beyond 58.4 degrees of
  latitude in either hemisphere, and 78.8 degrees west and 101.2 degrees east, between 31.7
  degrees south and north; values of one decimal, which a box typed by hand, or a tiling by
  fifths of a degree, lands on exactly. For a segment of one degree of the meridian of 11.2
  degrees east the dearest step took 7 milliseconds at level 20 of ISEA3H, 228 at level 26, and
  286 at level 14 of IGEO7; a third of a degree away from the meridian, 0.3 millisecond. A side
  of `n` cells puts at least `n` zones in the answer, so **a caller that refuses a box whose
  estimate exceeds a limit `N` of its own bounds a step by `N` cells**; one that sets no such
  limit can be made to spend seconds on one step. An ordinary step takes a fraction of a
  millisecond: a page of a thousand zones of a box was measured at 4.6 to 5.6 milliseconds at
  level 12 of ISEA3H and 15.9 to 17.7 at level 13 of IGEO7, the first page and one far into the
  answer alike. `zones`, which tests
  no extent, costs some 12 to 27 nanoseconds a zone on the aperture-3 grids and 2 to 5
  microseconds on the aperture-7 grids, where each zone is converted to its Z7 identifier.
- **What the estimate promises.** `estimate_zones_in_box` is the area of the box grown by 1.75
  zone widths on every side, taken on the sphere on which the zones are equal in area, over the
  area of a zone, rounded up, plus 3, and never more than the count of the level plus 3. It is an
  upper bound by measurement, not by proof: it was never under the count on any box on which the
  two were compared, at every level that `zones_in_box` answers, beside the poles, over the
  antimeridian and about the vertices of the icosahedron, from a single point to caps and bands
  of half a million and of a million and a half zones. Its thinnest margin was measured upon a
  pole at the finest level of the aperture-3 grids: a wedge on the north pole at level 33, a
  millionth of a zone width high and one degree wide, holds 8 zones and is estimated at 9, on
  ISEA3H, IVEA3H and RTEA3H alike, so that there the bound holds by the 3 that it adds. At
  levels 31 and 32 the thinnest margin measured is 6 zones, and on the aperture-7 grids 7. How
  far above the count it lies weighs by the size of the answer: an answer of under ten zones may
  be estimated at up to some twenty-one times its size (a point away from the poles, which one
  to three zones hold, is estimated at 16 zones, and at 12 to 16 at the two coarsest levels);
  one of ten to a hundred zones at up to some five and a half times; one of a hundred to a
  thousand at up to some twice; one of a thousand to ten thousand at no more than a third above;
  one of ten thousand or more at no more than an eighth; and a band about the world of half a
  million zones at eight parts in a thousand. A caller that must not be surprised keeps a limit
  of its own on what it takes from the iterator as well.
- **What the compaction gives.** `compact_zones` is DGGAL's `compactZones` for the zones of one
  level, the engine's answer zone for zone and in its order, and its cost is that of the slice.
  On the aperture-3 grids it goes two levels at a time: a zone two levels up takes the place of
  the seven sub-zones that lie wholly within it, so that the zones of the answer are of the
  level of the set or an even number of levels above, and every zone of a level becomes the
  twelve of level 0. On the aperture-7 grids it goes one level at a time: a zone takes the place
  of its children where the set has all thirteen, or eleven, and a child gives way only when
  every one of its parents was taken, so that a lone zone with its children replaces the child
  at its centre alone. **The zones of the answer overlap**: the thirteen sub-zones two levels
  below `B4-2-A` on ISEA3H compact to that zone and the six on its vertices, whose areas sum to
  fifteen thirteenths of the set's; the thirteen children of `0064` on IGEO7 compact to thirteen
  zones, `0064` in the place of `00640`, whose areas sum to nineteen thirteenths. **The area that
  a compacted list covers is the area of the set, not the sum of the areas of the list.** The
  order of the answer is the engine's: on the aperture-3 grids the ascending order of the
  identifiers, coarser zones first; on the aperture-7 grids the order of the engine's internal
  identifier, in which a level and the odd level above it share their cells, so that their
  zones alternate (`0064`, of level 2, stands between `03323` and `00641`, of level 3); a caller
  that wants the coarser zones first sorts the answer by level, stably. **On the aperture-3
  grids a set with a gap may be covered beyond itself**, by two routes: the engine takes a
  coarser zone when the six neighbours of the sub-zone at its centre are in the set, and does
  not look for that sub-zone itself; and from the second pass on it counts a sub-zone on a
  vertex of the coarser zone as present where the answer already holds the one zone, two levels
  finer, at its centre. A set without a gap, the zones of a box, the sub-zones of a zone or a
  whole level, is covered exactly. A page of a longer answer is a set with a gap where it ends
  before a zone whose neighbours it holds, the zone on a pole, which comes last, among them: if
  each page is compacted by itself, the pages together still cover exactly the whole answer and
  nothing outside it, but a zone may be covered by two pages, as was measured in one sample of
  five boxes at levels 5 to 8 of ISEA3H, so that a client that sums over the expanded pages removes the repeats. On the
  aperture-7 grids, away from the two broken edges (see Known limitations), no set was found
  covered beyond itself, a coarser zone being taken only where the set has every one of its
  children. There an identifier that the engine reads as
  another zone comes back as that zone's: an identifier given by its bits that names the child a
  pentagon does not have comes back as the engine reads it, so that `0x05FF_FFFF_FFFF_FFFF`,
  whose bits spell `002` (a text that `zone_from_text` refuses, the pentagon `00` having no such
  child), comes back as `003`. A compacted list cannot be entered at a position: compaction is a
  function of the whole set, and the cursor of `zones_in_box` pages the answer before it.
  Measured on an optimised build, a thousand zones of level 14 of IGEO7 are compacted in some
  60 milliseconds, as the engine compacts them, the cost of a zone growing with the level.
- **Names.** `GRID_NAMES` lists the names `get_grid` accepts, in the order it tries them.
  `get_grid` matches a name without regard to case, and a grid's `name()` answers in the
  canonical form whatever spelling reached it.
- **The contract.** This crate's stated contract is that caller input to `Grid`, `Zone` and
  `AnyGrid` returns an error rather than panicking. The traits `Projection`, `Topology` and
  `Indexing`, which are sealed, and the implementations of them this crate ships are its own parts
  made visible so that a caller can compose a grid from them; the contract does not extend to
  them, and their documentation states what each expects. The methods that `Grid` asks of its
  parts for itself alone, such as the compaction of a topology, the walk of the lattice in which
  it enumerates the zones of a level, or the inverse of a projection in radians, are not reachable
  from outside the crate. Where one of them answers differently
  from a DGGAL accessor of similar name, as `Topology::planar_centroid` does from
  `getZoneCRSCentroid`, its documentation says so.

## Known limitations

- **The engine draws 68 polar sub-hexagons on each aperture-3 grid that this crate gives no
  geometry.** They are the sub-hexagons C and D of the two polar roots, such as `CA-0-C`, at each
  of seventeen odd resolutions: two roots, two sub-hexagons, seventeen resolutions. The engine
  reads their text and draws them, with a ring, an extent and an area. No zone of the crate has
  ever reached them, since no neighbour, parent, child or sub-zone of any zone compared is such a cell,
  and `zone_from_text` refuses them (see Departures from the engine), so `refined_vertices`,
  `extent` and `area` answer the empty ring, `None` and `None` for them where the engine has an
  answer to each. This is the largest departure of the methods that give geometry. It is
  characterised live on each grid, all 68 of them counted, and whether to give these zones a cell
  and geometry has not been decided.
- **The aperture-7 engine contradicts itself along two icosahedron edges.** Along the edge from
  base cell 0 to base cell 1 over the north pole, and the edge from base cell 1 to base cell 6,
  the engine's answers stop being consistent with one another from resolution 15 upwards: its
  odd resolutions leave a strip with no cells, and its own neighbour lists there include its null
  zone, identifiers it cannot read back, the zone itself and cells far away, and are not
  reciprocal. Its lists of parents and of children and its sub-zone orders are in the same case
  there, the children from resolution 14, and a later item of this list describes them. This
  crate's neighbour lists there are now the engine's own, less the four kinds
  documented under Departures from the engine, and Verification compares them exactly rather than
  admitting a disagreement; there is nothing left there for a ceiling to hold. The suite still
  admits a disagreement inside the band on quantisation alone: a quantisation at which the two
  sides name different cells, neither of which contains the input, measured at none in the
  suites' own samples on 2026-09-29; and an identifier the engine cannot read back, measured the
  same day at up to 3 points of the main sample and up to 46, 41 and 37 of the seam sample on
  IGEO7, IVEA7H and RTEA7H. The engine's own null answers lie within 7.5e-5 degrees of the edges
  (measured 2026-09-23), and the unreadable answers admitted lay within 5.73e-5 degrees when
  last measured (2026-09-26); neither reach has been re-measured since. The band itself stays
  1.5e-4 degrees either side of each edge, which about each end of an edge becomes a disc of
  that radius about the pentagon there. The other twenty-eight edges match
  at every resolution, and so do resolutions 13 and 14 along these two. On the aperture-3 grids
  the suites met no such region in any of their samples.
- **The engine's own far neighbours are passed through, and some hexagons list fewer than six.**
  Because this crate's aperture-7 lists are now the engine's own, an entry it keeps can lie far
  from the zone: on IGEO7 the list of `00515151515151510530`, at resolution 18, holds
  `00454545454545405216`, `...212` and `...213`, each 35.88 degrees away, three of its six
  neighbours. That is the engine's own answer, not a construction of this crate's, and it is not
  corrected. Measured on 2026-09-29 over the suites' own samples: none of the 2,860 deep digit
  paths under the two poles, nor the 2,012 to 2,020 zones sampled near the south pole, nor the
  6,435 to 6,480 zones sampled along the thirty icosahedron edges, lists an entry beyond one
  degree away; near the north pole, 1,788 zones sampled on each grid, 31 or 32 zones per grid
  do, 69 or 70 such entries in all, none nearer than 35.88 degrees. Dropping the four kinds also
  shortens some lists: 461 hexagons and 3 pentagons of the deep digit paths, 184 to 186 hexagons
  near the north pole and 72 to 81 of those sampled along the edges, all of them on the two
  broken edges, list fewer neighbours than they have sides. No hexagon lists fewer than two, and
  no list near the south pole is short.
  Neither side is an oracle in that band, and a caller that needs neighbours there should expect
  both the short lists and the far entries.
- **The aperture-7 parents, children and sub-zone orders at the two broken edges are the
  engine's, and are not corrected.** This crate drops from a list of parents or children what is
  no zone (see Departures from the engine) and corrects nothing else there: a hierarchy that
  held at those edges would need a quantiser that works on them, which neither the engine nor
  this crate has. What follows was measured against the live engine, and is the same on IGEO7,
  IVEA7H and RTEA7H wherever one figure is given, and in that order wherever three are.
  - **Parents, at identifiers that a text can name.** From resolution 16 some identifiers
    beside the edges have no parent, although they lie above resolution 0, and some have a
    first parent that is not the zone their identifier names with its last digit dropped:
    `00000000000000001311` has none, and `0000000000000001644` has `000000000000000005` and
    `000000000000000000`. Of the 84,036 identifiers that `00`, zeros and five further digits
    name at resolutions 14 to 19, 3,941 have no parent and 3,556 another first parent; the suite
    holds 2,276 identifiers, one in 37 of those and eight chosen ones, among which they are 107
    and 98, to exact counts at each resolution. **Not one of those identifiers is the zone found
    at its own centroid**, by the engine's quantisation or by this crate's, which the suite
    asserts of every one it meets;
    and among the zones that the main sample's positions give, some 5,600 on each grid at
    fifteen resolutions between 0 and 19, there is none of either kind. As measured, then, a
    caller who obtains zones from positions does not meet them, and one who writes identifiers
    as text can. It is a measurement of one region, under base cell 0, and no proof; the edge
    from base cell 1 to base cell 6 was not sampled by text.
  - **Parents, at zones that positions give.** One anomaly of the parents is met there, and is
    not to be grouped with the two above: a zone that is no centroid child and has one parent
    all the same. It occurs at 4, 3 and 2 zones of the main sample, at resolutions 15 to 19, and
    at 100 of the 2,276 identifiers built by text.
  - **Children.** A list of children shortened by the three kinds dropped does occur at zones
    that are the zone found at their own centroid, from resolution 14: at 33, 26 and 26 zones of
    the main sample, at resolutions 14 to 18, and at 882 of the 84,036 identifiers. At those
    edges a zone need not be named among the parents of a zone it lists as a child (87, 66 and
    66 such pairs in the main sample), and the zones that an identifier names with one digit
    more need not be the first of its children (150 of the 2,276 identifiers built by text, and
    none of the zones from positions).
  - **Sub-zone orders.** Among sub-zones of resolution 15 and finer, an order at those edges may
    hold the null zone, where the engine's own holds its null zone (3,078, 2,982 and 3,072
    entries of the orders the suite compares) or an identifier it cannot read back (275, 260
    and 260); it may name a zone twice (35 entries on each grid); and it may name a zone that is
    no descendant of the zone, or leave a descendant out, which the suite compares with the
    engine's order and does not count apart. All of that is the engine's and is kept: the length
    of an order is `count_sub_zones` and the place of every entry is the engine's, so that a
    client that holds the engine's order reads the same zone at the same place. Of the 17,053
    sub-zones of `010004000400` at depth 5, 99 are the null zone. Of the orders the suite
    compares, 265, 249 and 252 hold the null zone or a repeat.
  - **The sub-zones at depth 1 are not always the children there.** The engine's `getSubZones`
    and its `getZoneChildren` part at some zones of those edges, among sub-zones of resolutions
    16 and 18: at 5, 3 and 3 hexagons of the main sample and at 57 of the identifiers built by
    text. The order of `0000000000000001644` holds thirteen zones, and its children are four.
  - **The index of a sub-zone.** The engine's `getSubZoneIndex` answers -1 for some entries of an
    order at those edges. Where its own test of a sub-zone refuses the entry, the entry has no index
    here either: at 25, 69 and 31 of the entries asked below hexagons of the main sample, 39 below
    the pentagons and 2,988 below the identifiers built by text. Where the test accepts it, this
    crate finds it in the order itself, at its first place (see Departures from the engine): at 3,
    11 and 7 below hexagons, 43 below the pentagons and 13 below the identifiers built by text. Of
    the order of `0000000000000005` at depth 1, the entry at 4 has no index, and the entry at 9 has
    its place as its index, where the engine has none. Where an order names a zone twice, its index
    is one of its places.
  - **The compaction.** `compact_zones` gives the engine's answer there as elsewhere, from the
    engine's lists of parents and children as they are. A zone for which the engine finds no
    parent leaves the answer and no zone stands for it (`00000000000000001311`, at resolution
    18), so that a caller can lose a zone there in silence; a zone is taken although the engine
    lists its null zone among its children; the sub-zones of the answer need not be the set; and
    an identifier that the engine reads as another zone comes back as that zone's
    (`0000000000000001644` as `0000000000000000056`).
  None of this was met at a zone coarser than resolution 14, nor in an order whose sub-zones are
  coarser than resolution 15: below those resolutions the suites assert that every order
  compared is sound and that every class above is absent.
- **The zones of a level and of a box end at level 14 on the aperture-7 grids.** `zones` and
  `zones_in_box` enumerate the lattice in which the engine orders the zones of a level, and on
  the aperture-7 grids that lattice is not whole from level 16: beside the same two edges there
  are cells that host no zone under an identifier of its own, some twelve million of the 3.3e14
  cells of levels 16 and 17 and some four thousand million of the 1.6e16 of levels 18 and 19,
  and addresses there that are no zone of the lattice. None was found to level 15: the lattice
  is the engine's own list of every zone to level 5, a census of every cell to level 9 found
  each hosting its zones, and samples along the borders of every rhombus beyond found none
  empty.
  The search for the zones of a box also takes a zone to lie within one side of a cell of the
  corner of the cell that hosts it: to level 14 none was found beyond 0.44 of a side, and from
  level 15 there are zones along those edges that lie farther, which the search would not find.
  So both enumerations, and the estimate with them, refuse a level beyond `max_box_level`, 14,
  with `Error::ResolutionOutOfRange`: zones of some 75 square metres, where the grids reach
  level 19. On the aperture-3 grids they answer at every level, to 33. The zones of one zone at
  a finer level, which a service lists for a parent, are its sub-zones, and `sub_zones` and
  `sub_zone_at_index` answer them at every level.
- **A sub-zone of an identifier that is not the zone found at its own centroid has no index.**
  The sub-zone methods accept an identifier that has an address although no position quantises
  to it, such as the child that a pentagon lacks, as the engine does: `sub_zones` lists entries
  below it, and `sub_zone_index` answers `None` for each of them, at any resolution, because the
  engine's test of a sub-zone finds the real zone at the coarser level and not the identifier
  asked about. In the suites' sweep of hostile identifiers that is 26 of the 3,307 entries asked
  for their index, on each grid.
- **Only the registered configurations are verified.** `Grid::new` accepts any finite
  `GridConfig`, but the engine takes no orientation from a caller, so a grid built on another
  configuration has no engine to be compared with. What is verified is the six grids that
  `igeo7`, `ivea7h`, `rtea7h`, `isea3h`, `ivea3h`, `rtea3h` and `get_grid` return.

## Rust version

`rust-version = "1.85"` is the minimum Rust toolchain that builds and uses `rs4dggs`, the
published crate, on its own; the test-only oracle that exercises it against a live DGGAL needs
Rust 1.87, which concerns contributors running the full local gate, never a consumer of the
published crate.

## Licence

MIT (see `LICENSE`). `NOTICE` reproduces the BSD 3-Clause licence of DGGAL, whose algorithms this
crate ports.
