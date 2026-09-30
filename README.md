# rs4dggs

A pure Rust port of [py4dggs](https://github.com/terraops-org/py4dggs), which is itself a Python
port of DGGAL's discrete global grid systems: IGEO7, IVEA7H, RTEA7H, ISEA3H, IVEA3H and RTEA3H.
DGGAL, not py4dggs, is the reference for both ports (see Why it agrees, below). The crate has no
dependencies and no `unsafe` code.

## Getting started

Add the crate with `cargo add rs4dggs` once it is published; until then, depend on the repository
by `git` or by `path`.

A grid is chosen once: at compile time with `igeo7()`, `ivea7h()`, `rtea7h()`, `isea3h()`,
`ivea3h()` or `rtea3h()`, or by name at run time with `get_grid("igeo7")`, in any casing.
Everything else starts from a zone, either a coordinate quantised at a resolution, from 0 to 19 on
the aperture-7 grids and from 0 to 33 on the aperture-3 grids, or a text identifier read back.
From a zone come its centroid and boundary, its neighbours and the disk around it, its place in the
hierarchy and, on the aperture-3 grids, its sub-zones. The example below shows each of these in
turn.

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

let parent = zone.parent().unwrap();     // the congruent Z7 digit path, one digit shorter
assert_eq!(parent.text_id(), "006415");
assert!(parent.is_ancestor_of(&zone));
assert!(!zone.is_sibling_of(&zone));     // a zone is not its own sibling, as DGGAL has it

// The disk of radius 2: a lazy walk from the centre outwards, ring by ring, in a
// defined order. Only what is taken is computed, so any radius may be asked for.
let disk: Vec<_> = zone.disk(2).collect();
assert_eq!(disk.len(), 19);
assert_eq!(disk[0], zone);
let ring_sizes: Vec<usize> = zone.disk(2).rings().map(|ring| ring.len()).collect();
assert_eq!(ring_sizes, [1, 6, 12]);
let first_seven: Vec<_> = zone.disk(u32::MAX).take(7).collect();
assert_eq!(first_seven[..], disk[..7]);

// The same point under the other two projections. At this resolution the three cells
// happen to coincide; four resolutions further down they no longer do, since the three
// projections place the same coordinates differently within the icosahedron's faces.
assert_eq!(ivea7h().zone_from_geo(38.7223, -9.1393, 5)?.text_id(), "0064156");
assert_eq!(rtea7h().zone_from_geo(38.7223, -9.1393, 5)?.text_id(), "0064156");
assert_eq!(igeo7().zone_from_geo(38.7223, -9.1393, 10)?.text_id(), "006415654630");
assert_eq!(ivea7h().zone_from_geo(38.7223, -9.1393, 10)?.text_id(), "006415650050");
assert_eq!(rtea7h().zone_from_geo(38.7223, -9.1393, 10)?.text_id(), "006415650342");

// An aperture-3 grid. Its hierarchy is DGGAL's, and not congruent: this zone has three
// parents and seven children.
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

Four small programs show the same in context, each run with `cargo run --example <name>` from a
clone of the repository, which needs nothing but Rust; for instance
`cargo run --example quantise -- 5 38.7223 -9.1393`:

- `quantise`: points on all six grids at one resolution, from 0 to 19, the range they share;
- `hierarchy`: a zone's ancestors and children on IGEO7;
- `geojson`: a zone of IGEO7 and its neighbours as GeoJSON, with no dependency;
- `service`: a grid chosen by name, and a disk served ring by ring.

`cargo doc --open` renders the full reference. The entry points a newcomer meets first, quantising,
reading an identifier, the boundary, the neighbours, the hierarchy, the disk and choosing a grid by
name, each carry an example of their own.

## Command-line tool

The workspace also holds `rs4dggs`, a command-line tool over this library, in `crates/rs4dggs-cli`;
its README describes every command, the formats and batch use. Three examples:

```text
rs4dggs igeo7 zone 38.72,-9.14 10
rs4dggs igeo7 neighbours 006415654636
rs4dggs isea3h geom C2-23-C > cell.geojson
```

## Status

All six grids are implemented. The three aperture-7 grids, **IGEO7**, **IVEA7H** and **RTEA7H**,
share one aperture-7 hexagonal topology and one Z7 indexing; the three aperture-3 grids,
**ISEA3H**, **IVEA3H** and **RTEA3H**, share one aperture-3 hexagonal topology and one I3H
indexing. Within each family the grids differ only in the projection laid over the icosahedron,
ISEA, IVEA and RTEA respectively. Each grid offers quantisation from geographic coordinates,
centroids, cell boundaries, edge neighbours, its hierarchy and the relational predicates over it,
the disk of radius k, and text identifiers, at every resolution its indexing can hold: 0 to 19 for
Z7 and 0 to 33 for I3H.

The two hierarchies differ. On the aperture-7 grids it is the congruent Z7 digit path, in which a
zone has one parent. On the aperture-3 grids it is DGGAL's own, which is not congruent: a zone has
one parent or three, and seven children, or six under a pentagon. The aperture-3 grids also order
their sub-zones at any depth, in DGGAL's order, corrected where the engine's own is at fault (see
Departures from the engine), and page through that order by index without building it.

The identifiers `neighbors` and `quantize` follow py4dggs's spelling.

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
- **The hierarchy.** On the aperture-7 grids the congruent hierarchy, which DGGAL's geometric
  hierarchy does not match, is checked through the engine's reading of text identifiers, which
  spell the digit path itself. DGGAL's `isZoneAncestorOf` and `areZonesSiblings` answer over its
  geometric hierarchy, so `is_ancestor_of` and `is_sibling_of` are checked by the same route
  instead. On the aperture-3 grids the hierarchy is the engine's own, and parents, children, the
  centroid parent and `is_centroid_child` are compared with its answers directly, in its order,
  as are `is_ancestor_of`, `is_immediate_child_of` and `is_sibling_of` over some 37,000 distinct
  pairs of zones on each grid, each pair asked once and the engine's answers counted exactly.
  `disk`, which has no counterpart in DGGAL, is compared in order at radii 0, 1 and 2 with a walk
  composed of the engine's own neighbour lists, less, on the aperture-7 grids, the same four
  kinds, on every grid.
- **The sub-zone methods** are checked on the aperture-7 grids only as refusals (see Known
  limitations). On the aperture-3 grids the counts, the first sub-zones, the orders and the
  indices are compared with the engine's, over some 1,500 orders holding some 27,000 sub-zones on
  each grid and over sub-zones at every depth down to resolution 33. `sub_zone_at_index` is
  compared with the order it indexes, since the binding the suite links does not pass a depth to
  the engine's own `getSubZoneAtIndex`; a separate probe found the two the same at all 4,800
  pairs it compared. The input this crate refuses, and every answer in which it departs from the
  engine, is checked against the engine's own answer to it, save the aperture-3 sub-zone order
  corrected at the rhombus edges: there the order is checked against the engine's own descendant
  test rather than its answer, and the engine's index, which refers to its own order, is not asked
  (see Departures from the engine).

A run of 2026-09-29 compared 97,877 cases on IGEO7, 97,581 on IVEA7H, 97,593 on RTEA7H,
199,905 on ISEA3H, 199,806 on IVEA3H and 199,565 on RTEA3H, 892,327 in all, each test asserting
a floor on the number it compared. In that run no tie remained in any sample, that is, no input
on a cell boundary at which the two sides name different cells, and the suite holds its admission
of ties at none, so that a tie which reappears fails. The suite asserts that every geographic
point it compares with the engine's is the engine's own to the last bit, latitude and longitude
alike, and in that run all 355,162 such points were (67,466 on IGEO7, 67,227 on IVEA7H, 67,179
on RTEA7H, 51,117 on ISEA3H, 51,106 on IVEA3H and 51,067 on RTEA3H), so that the largest
coordinate gap measured by each geometric test was zero. That assertion is made for x86-64 Linux
with glibc, where the suites run, and for no other platform.

A disagreement is admitted in two ways only. Either it is a characterised divergence, in which
this crate's answer has been shown to be the defensible one and the engine's behaviour is
re-checked live at every instance on every run; or it lies inside the region along two
icosahedron edges where the aperture-7 engine is not consistent with itself (see Known
limitations), where each case is counted and held under a ceiling. On the aperture-3 grids no
such region was met, and every disagreement admitted there is a characterised divergence.
Anything else fails.

Beside the six suites, two unit tests compare what no public method returns. The first compares
the refined cell boundary over which the aperture-7 containment test is performed with the
nearest accessor DGGAL exposes, at 360 zones fixed by construction on each aperture-7 grid. The
2,666 vertices of 343 of those zones are compared bit for bit, with one exception pinned by zone,
index and bit pattern; the other 17 zones, at which DGGAL's accessor subdivides a side twentyfold
and so traces a different polygon, are counted exactly, so that none can leave the comparison
unnoticed. The second compares the points from which the aperture-3 sub-zone orders are
generated with the engine's own, bit for bit: 4,105,210 of them, over 17,710 pairs of a zone and
a depth on ISEA3H, save at 6,884 points in 98 of those pairs where the engine's own order is at
fault and this crate's is corrected (see Departures from the engine), and there checked against
the definition of a sub-zone instead.

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

The local gate reads the BuildID of the `libdggal.so` the oracle links and fails if it is another.
A different build of DGGAL, a rebuilt 0.0.6 or a later release, may compile these functions
differently, and the orders would then have to be read again from that library.

## Departures from the engine

Where the engine answers input that has no meaning with an answer that looks like one, answers a
question without examining it, or answers in contradiction of its own answers, this crate departs
from it on purpose; and it departs from the order of the engine's vertices on the aperture-3 grids,
so that a ring runs the same way on every grid. Each departure is documented on the method that
makes it. Where the suite re-checks the engine's side of a departure on every run, so that the
departure stays as described, the item says so.

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
- **The first sub-zone at depth 0**, on the aperture-3 grids. `first_sub_zone(z, 0)` and
  `sub_zone_at_index(z, 0, 0)` answer the zone itself, as `sub_zones(z, 0)` does on both sides.
  The engine's `getFirstSubZone` has no case for depth 0: it quantises a vertex of the zone at the
  zone's own level and so names a neighbour, for about two thirds of zones (763 of the 1,176
  sampled on ISEA3H, 760 of 1,175 on IVEA3H and 757 of 1,175 on RTEA3H). Re-checked by the suite,
  which finds the engine's answer among the zone's neighbours.
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
  back.
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
  resolution 33.
- **The null zone's resolution** is reported as 0 on every grid, so that the null zone reads
  alike on all six, where DGGAL reports -1 on the aperture-7 grids and, reading the sentinel's
  fields as they stand, 63 on the aperture-3 grids. A caller who needs to tell a resolution-0 zone
  from the null zone compares the identifier with `ZoneId::NULL`, or its text with `NULL_TEXT`.
  The suite pins both values, and checks the null zone against the engine in everything else.
- **Refusals name their cause, at bounded length.** `zone_from_text` returns
  `Error::InvalidZone` for text that names no zone, where the engine answers its null zone for
  any text it cannot parse (`RI7H_Z7.ec:418`, `RI3H.ec:908-931`). A rejected text identifier is
  quoted in the message to its first thirty-two characters, followed by its length, and an
  unknown grid name is kept to its first thirty-two characters, so that no error carries caller
  input at unbounded length. The suite does not compare these.

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
- **Sub-zones cost what they list.** On the aperture-3 grids `sub_zone_at_index` pages through an
  order without building it, at any depth down to resolution 33, in a number of steps of the order
  of 3^(d/2) at depth d; a single call at depth 33 below a zone of resolution 0 takes under a
  tenth of a second in a release build. `count_sub_zones` answers at any depth. `sub_zones` builds the
  whole list and refuses one longer than `MAX_MATERIALISED_SUB_ZONES`, four million, with
  `Error::TooManySubZones`, which below a hexagon means any depth beyond 13. `sub_zone_index`
  searches the list that `sub_zones` builds, and is bound by the same limit: the index of a
  sub-zone more than 13 levels below a hexagon is refused with that error rather than searched
  for.
- **Names.** `GRID_NAMES` lists the names `get_grid` accepts, in the order it tries them.
  `get_grid` matches a name without regard to case, and a grid's `name()` answers in the
  canonical form whatever spelling reached it.
- **The contract.** This crate's stated contract is that caller input to `Grid`, `Zone` and
  `AnyGrid` returns an error rather than panicking. The traits `Projection`, `Topology` and
  `Indexing`, which are sealed, and the implementations of them this crate ships are its own parts
  made visible so that a caller can compose a grid from them; the contract does not extend to
  them, and their documentation states what each expects. Where one of them answers differently
  from a DGGAL accessor of similar name, as `Topology::planar_centroid` does from
  `getZoneCRSCentroid`, its documentation says so.

## Known limitations

- **The aperture-7 engine contradicts itself along two icosahedron edges.** Along the edge from
  base cell 0 to base cell 1 over the north pole, and the edge from base cell 1 to base cell 6,
  the engine's answers stop being consistent with one another from resolution 15 upwards: its
  odd resolutions leave a strip with no cells, and its own neighbour lists there include its null
  zone, identifiers it cannot read back, the zone itself and cells far away, and are not
  reciprocal. This crate's neighbour lists there are now the engine's own, less the four kinds
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
- **No sub-zone order on the aperture-7 grids.** The aperture-7 topology, like py4dggs's, defines
  no sub-zone order over the congruent Z7 hierarchy it implements, so `sub_zones`,
  `count_sub_zones`, `first_sub_zone` and `sub_zone_at_index` return `Error::NoSubZoneOrder`.
  DGGAL's own IGEO7 does define one, over its geometric hierarchy (`RI7H_Z7.ec:492`, `:576`); a Z7
  sub-zone order is future work. `sub_zone_index` answers without an order only where none is
  needed: a zone is its own sub-zone at index 0, and a different zone at the same resolution, or
  a coarser one, is none; a finer one reaches the missing order, and the error.
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
