# rs4dggs-cli

`rs4dggs` is a command-line tool over the `rs4dggs` library, a pure-Rust port of the discrete
global grid systems of DGGAL (IGEO7, IVEA7H, RTEA7H, ISEA3H, IVEA3H and RTEA3H). It finds the zone
that contains a point, prints a zone's card, boundary, neighbours and disk, relates two zones,
lists a zone's sub-zones in order, and answers for many points or zones read from standard input,
as text, CSV or GeoJSON. Its dependencies are the library and its companion crate of encodings,
`rs4dggs-ogc`.

Its invocation follows DGGAL's own tool, `dgg`, and for every command that both have, the
identifiers and coordinates it prints are those `dgg` prints, save where the library departs from
the engine on purpose (see "Differences from `dgg`").

## Installing

From a clone of the repository:

```
cargo install --path crates/rs4dggs-cli
```

and, once the crate is published, `cargo install rs4dggs-cli`. The binary is named `rs4dggs`. The
minimum Rust version is 1.85.

## Invocation

```
rs4dggs <grid> <command> [options] [arguments]
rs4dggs grids
rs4dggs help
rs4dggs <grid> help [command]
```

- **Grids:** `igeo7`, `ivea7h`, `rtea7h` (aperture 7), `isea3h`, `ivea3h`, `rtea3h` (aperture 3),
  in any case. DGGAL's names `isea7h_z7`, `ivea7h_z7` and `rtea7h_z7` are accepted for the
  aperture-7 grids.
- **A point** is latitude then longitude in degrees (WGS84), as one argument, `38.72,-9.14`, or as
  two, `38.72 -9.14`. The latitude must lie within [-90, 90] and the longitude within
  [-180, 180]; anything else is refused.
- **A zone** is its text identifier (`006415654636`, `C2-23-C`), its 64-bit identifier in
  decimal, or in hexadecimal after `0x`. A decimal identifier that also reads as a Z7 text
  identifier is read as text; `0x` is unambiguous.
- **Options** may stand anywhere after the command, in `dgg`'s single-dash form or with two
  dashes:

| Option | Meaning | Default |
|---|---|---|
| `-f`, `--format text\|csv\|geojson` | the output format | text for one input, CSV in batch, GeoJSON for `geom` |
| `-o`, `--output <file>` | write to a file | standard output |
| `-precision <n>` | `n` decimal places for every coordinate, in every format (0 to 17) | the shortest exact form |
| `-depth <n>` | the relative depth for `sub` | 1 |
| `-centroids` | GeoJSON points at the centroids instead of polygons | polygons |

`-h`, `--help` and `help` print help anywhere, and `-V` and `--version` the version. The command
`rs4dggs <grid> help <command>` (or `rs4dggs help <grid> <command>`) explains one command with
examples on that grid. The version printed is the tool's, which is also the library's.

`-precision` and `-centroids` apply where a coordinate is printed. They are refused, with exit
code 2, for `grids`, the grid's `info`, `rel`, `sub` and `index`, which write text with no
coordinate.

`-o` creates the file only when there is something to write: a command refused before it has an
answer leaves an existing file untouched. A batch run writes its header, or an empty collection,
over the file even when every line is bad, and then exits 1.

## Commands

### `grids` and `info`

```
$ rs4dggs grids
grid    aperture  resolutions  projection  indexing
IGEO7   7         0 to 19      ISEA        Z7
IVEA7H  7         0 to 19      IVEA        Z7
RTEA7H  7         0 to 19      RTEA        Z7
ISEA3H  3         0 to 33      ISEA        I3H
IVEA3H  3         0 to 33      IVEA        I3H
RTEA3H  3         0 to 33      RTEA        I3H
```

```
$ rs4dggs igeo7 info
IGEO7
  aperture     7 (each zone has about 7 times the area of its children)
  resolutions  0 to 19
  projection   ISEA (icosahedral Snyder equal-area)
  indexing     Z7
```

### `zone`: the zone containing a point

With a resolution, the answer is the zone's card:

```
$ rs4dggs igeo7 zone 38.72,-9.14 10
IGEO7 zone 006415654636  (integer 940643638281502719, hex 0D0DD667BFFFFFFF)
  resolution  10, hexagon
  centroid    38.71963792895191, -9.140391954153882   (lat, lon)
  parents     2: 00641565463, 00641565462
  children    13: 0064156546360 (centre), 0064156546361, 0064156546365, 0064156546364, 0064156546366, 0064156546362, 0064156546363, 0064156546342, 0064156546033, 0064156546251, 0064156546215, 0064156546324, 0064156546306
  neighbours  6: 006415654634, 006415654603, 006415654625, 006415654621, 006415654632, 006415654630
  vertices    6:
              38.72191423979572, -9.141012890880889
              38.720580521455155, -9.143387999135763
              38.71830421654527, -9.142767009043583
              38.71736161646542, -9.139771072963079
              38.71869526727886, -9.137395980908451
              38.720971585697846, -9.138016808726869
```

Without a resolution, the answer is that point's zone at every resolution of the grid:

```
$ rs4dggs igeo7 zone 38.72,-9.14 | head -8
IGEO7 zones containing 38.72, -9.14
  res  zone
    0  00
    1  006
    2  0064
    3  00641
    4  006415
    5  0064156
```

Where the grid assigns no cell to a point (a few strips of the aperture-7 grids at fine
resolutions, described in the library's documentation), `zone` says so in words, the zone field of
the CSV is empty, and GeoJSON writes no feature for it.

### `info <zone>`: the zone's card

The card is the one shown above, and `rs4dggs igeo7 info 006415654636` prints it. The parents and
the children are DGGAL's own on every grid, in its order. On the aperture-3 grids a zone has one
parent or three, and seven children, or six under a pentagon. The parent that is itself a centroid
child is marked, as `dgg` marks it, and so is the child at the zone's centre:

```
$ rs4dggs isea3h info C2-23-C
ISEA3H zone C2-23-C  (integer 306244774661193870, hex 044000000000008E)
  resolution  5, hexagon
  centroid    37.731557553208035, -9.22932570934645   (lat, lon)
  parents     3: C2-23-A, C4-5-A, C4-6-A (centroid parent)
  children    7: D2-128-A (centre), D2-10C-A, D2-10D-A, D4-11-A, D4-10-A, D2-143-A, D2-127-A
  neighbours  6: C4-6-B, C2-23-D, C2-1A-D, C4-5-C, C2-23-B, C4-5-B
  vertices    6:
              39.714594062128235, -11.51357162615407
              37.43488060579155, -12.509252745189645
              35.47408521066121, -10.284053776669873
              35.71681487704544, -7.06851413841143
              37.93396678448876, -5.914425532835252
              39.97279511302393, -8.120762410070837
```

On the aperture-7 grids a zone has one parent or two, and thirteen children, or eleven under a
pentagon: first those whose identifiers continue the zone's own (seven, or six under a pentagon),
then those that lie under a neighbouring zone as well. The marks are the same:

```
$ rs4dggs igeo7 info 00641565463601 | head -5
IGEO7 zone 00641565463601  (integer 940643637241315327, hex 0D0DD66781FFFFFF)
  resolution  12, hexagon
  centroid    38.72015364223723, -9.140141353964504   (lat, lon)
  parents     2: 0064156546360 (centroid parent), 0064156546361
  children    13: 006415654636010 (centre), 006415654636011, 006415654636015, 006415654636014, 006415654636016, 006415654636012, 006415654636013, 006415654636162, 006415654636053, 006415654636001, 006415654636035, 006415654636344, 006415654636126
```

Where no parent is a centroid child, as for `006415654636` above, none is marked. On every grid a
zone of resolution 0 has no parent, and a zone of the finest resolution no children.

### `geom <zone>`: the boundary

The boundary is written as GeoJSON (section "GeoJSON" below); the format option `csv` is refused
with exit code 2. For instance, at three decimal places:

```
$ rs4dggs igeo7 geom 10656 -precision 3
{"type":"FeatureCollection","features":[{"type":"Feature","properties":{"zone":"10656","grid":"IGEO7","resolution":3},"geometry":{"type":"MultiPolygon","coordinates":[[[[178.612,1.754],[179.288,-0.438],[180.000,-0.595],[180.000,3.476],[179.895,3.500],[178.612,1.754]]],[[[-180.000,-0.595],[-178.744,-0.873],[-177.424,0.871],[-178.117,3.051],[-180.000,3.476],[-180.000,-0.595]]]]}}]}
```

### `neighbours <zone>`

The zones sharing an edge with the zone, in DGGAL's order (`neighbors` is also accepted):

```
$ rs4dggs igeo7 neighbours 006415654636
IGEO7 zone 006415654636: 6 neighbours
  006415654634
  006415654603
  006415654625
  006415654621
  006415654632
  006415654630
```

### `disk <zone> <k>`

Every zone within `k` neighbour steps, ring by ring, the zone itself first:

```
$ rs4dggs igeo7 disk 006415654636 1
IGEO7 disk of radius 1 about 006415654636
  ring 0 (1): 006415654636
  ring 1 (6): 006415654634, 006415654603, 006415654625, 006415654621, 006415654632, 006415654630
```

Answers are written as they are found. A radius that covers the whole sphere, at a fine
resolution, is a very long walk and runs until interrupted; the radius is not capped.

### `rel <zone> <zone>`

```
$ rs4dggs igeo7 rel 00641565463 006415654634
00641565463 is coarser than 006415654634 by 1 resolution
00641565463 and 006415654634 are not neighbours
00641565463 is an immediate parent of 006415654634
00641565463 is an ancestor of 006415654634
006415654634 is sub-zone 9 of 00641565463, at depth 1
```

The relations stated are: identical or not, which is coarser and by how many resolutions,
neighbours or not, parent and child, ancestor and descendant, siblings, and the depth and index
of a zone that is a sub-zone of the other. The index is left unsaid where the library gives
none, and where the order is longer than the library lists (see `sub` below).

### `sub` and `index`

The sub-zones of a zone at a relative depth, in DGGAL's order, or the one at a given index; and
the index of a sub-zone within a zone, at the depth their resolutions set:

```
$ rs4dggs isea3h sub A4-0-A
ISEA3H zone A4-0-A: 6 sub-zones at depth 1
      0  A2-0-D
      1  A2-0-C
      2  A3-0-C
      3  A4-0-B
      4  A4-0-C
      5  A4-0-D
$ rs4dggs isea3h sub A4-0-A 8 -depth 3
sub-zone 8 of A4-0-A at depth 3: B2-5-C
$ rs4dggs isea3h index A4-0-A B2-5-C
B2-5-C is sub-zone 8 of A4-0-A, at depth 3
```

On the aperture-7 grids the order bears no relation to the digits of the identifiers: it runs
across the zone line by line, as the engine's does:

```
$ rs4dggs igeo7 sub 00
IGEO7 zone 00: 11 sub-zones at depth 1
      0  013
      1  023
      2  005
      3  001
      4  053
      5  004
      6  000
      7  003
      8  033
      9  006
     10  043
$ rs4dggs igeo7 sub 0064156 27 -depth 2
sub-zone 27 of 0064156 at depth 2: 006415600
$ rs4dggs igeo7 index 0064156 006415600
006415600 is sub-zone 27 of 0064156, at depth 2
```

A position of an order may hold no zone. Along the two broken seams of the aperture-7 grids, at
fine resolutions, DGGAL's own order has the null zone at some positions, and the library keeps
every position as the engine has it. The tool prints `(no cell)` there, as `zone` does for a
resolution at which a point lies in no cell:

```
$ rs4dggs igeo7 sub 010004000400 -depth 5 | sed -n '205,207p'
    203  01000400040032330
    204  (no cell)
    205  01000400040032323
```

Of a zone that has no index in the order, `index` says so, with exit code 0:

```
$ rs4dggs isea3h index A4-0-A B6-5-C
B6-5-C has no index among the sub-zones of A4-0-A
```

A list of more than 4,000,000 sub-zones is refused, with exit code 1, and so is the index of a
sub-zone within an order so long; the environment variable `RS4DGGS_MAX_MATERIALISED_SUB_ZONES`
lowers that limit, and never raises it. One sub-zone by its index is answered at any depth. Both
commands write text alone.

## Formats

**Text** is written for a person. **CSV** has a header line and these columns:

| Command | Columns |
|---|---|
| `zone <point> <res>` | `lat,lon,zone` |
| `zone <point>` | `lat,lon,resolution,zone`, one row per resolution |
| `info <zone>` | `zone,resolution,shape,centroid_lat,centroid_lon,parent` (the parents separated by spaces) |
| `neighbours <zone>` | `zone,neighbour`, one row per pair |
| `disk <zone> <k>` | `zone,ring,member`, one row per member |

```
$ rs4dggs igeo7 zone 38.72,-9.14 10 -f csv
lat,lon,zone
38.72,-9.14,006415654636
$ rs4dggs igeo7 disk 006415654636 1 -f csv | head -4
zone,ring,member
006415654636,0,006415654636
006415654636,1,006415654634
006415654636,1,006415654603
```

**Numbers.** Every coordinate is written in the shortest form that reads back to the very same
double, so a person sees exactly what the library computed, and every coordinate the tool prints
is that double exactly. `-precision n` writes `n` decimal places instead, in every format (six
places are about 11 cm); a small negative number that rounds to zero is written as `0.000`, not
`-0.000`. A coarse `-precision` can collapse the finest cells about a pole into invalid polygons.

## Batch

`zone`, `info`, `geom`, `neighbours` and `disk` accept `-` in place of their point or zone, and
then read one per line from standard input:

```
$ printf '38.72,-9.14\n0,0\n' | rs4dggs igeo7 zone - 4 -f geojson -centroids -precision 3
{"type":"FeatureCollection","features":[{"type":"Feature","properties":{"zone":"006415","grid":"IGEO7","resolution":4,"lat":38.720,"lon":-9.140},"geometry":{"type":"Point","coordinates":[-8.878,39.128]}},{"type":"Feature","properties":{"zone":"036565","grid":"IGEO7","resolution":4,"lat":0.000,"lon":0.000},"geometry":{"type":"Point","coordinates":[-0.049,0.000]}}]}
```

- Blank lines and lines beginning with `#` are skipped. Line endings may be `\n` or `\r\n`, and
  the last line needs no newline. Bytes that are not UTF-8 are read as replacement characters and
  the line is reported as bad.
- A point line holds two numbers separated by a comma, spaces or a tab.
- A first line that does not parse and contains a letter is taken as a header and skipped. A
  zone identifier such as `C2-23-C` contains letters but parses, and so is read as input.
- A bad line is reported on standard error with its line number and the reason, and the run goes
  on; the exit code is 1 at the end if any line was bad:

```
$ printf 'lat,lon\n# c\n38.72,-9.14\r\nx y\n40 -8' | rs4dggs igeo7 zone - 8 -f csv
line 4: "x" is not a latitude
lat,lon,zone
38.72,-9.14,0064156546
40,-8,0065462446
```

- Answers are written as each line is read, so memory stays flat however long the input.
- `rel`, `sub` and `index` do not read standard input.
- When the reader closes the pipe early (`rs4dggs ... | head`), the tool stops quietly with exit
  code 0.

## GeoJSON

GeoJSON is written by `geom`, `info`, `zone`, `neighbours` and `disk`, for one input or in batch:
a FeatureCollection with one feature per zone, and the properties `zone`, `grid` and
`resolution`. With `-centroids` each feature is a point.

Each feature also names what it answers, so that a batch can be joined back to its input:

| Command | Added properties |
|---|---|
| `zone` | `lat` and `lon`, the point as the tool read it (`-precision` applies, as to every coordinate) |
| `neighbours` | `neighbour_of`, the zone asked about |
| `disk` | `ring`, and `centre`, the zone asked about |
| `geom`, `info` | none: the feature is the zone asked about |

The polygons are made clean for QGIS, GDAL, GEOS and web maps, as RFC 7946 asks. GeoJSON draws
edges as straight lines in latitude and longitude, so the tool cuts each ring where it meets the
antimeridian or a pole:

- A cell that crosses the antimeridian is written as a MultiPolygon of two polygons that meet at
  longitudes 180 and -180. The zone `10656` above is an example.
- A cell that reaches a pole reaches it along the pole's latitude, so that its outline runs along
  latitude 90 (or -90) between two longitudes, as the first polygon of zone `013` shows:

```
$ rs4dggs igeo7 zone 90,0 1 -f geojson -precision 3
{"type":"FeatureCollection","features":[{"type":"Feature","properties":{"zone":"013","grid":"IGEO7","resolution":1,"lat":90.000,"lon":0.000},"geometry":{"type":"MultiPolygon","coordinates":[[[[101.200,69.180],[134.624,61.696],[164.907,63.112],[180.000,66.571],[180.000,90.000],[63.367,90.000],[63.367,82.392],[101.200,69.180]]],[[[-180.000,66.571],[-159.068,71.368],[-116.633,82.392],[-116.633,90.000],[-180.000,90.000],[-180.000,66.571]]]]}}]}
```

  At resolutions 0 to 3 a pole lies on an edge shared by two cells, never inside one.
- Every ring runs anticlockwise, closed, with its first point repeated at the end, as the library
  gives its rings on every grid.
- A zone whose every vertex lies on a pole (the engine snaps some vertices to a pole; IGEO7 at
  resolution 19 has such zones) has no area. It is left out of GeoJSON, as the null zone is, and
  the tool says so on standard error; text and CSV still answer for it:

```
$ rs4dggs igeo7 zone 90,0 19 -f geojson
1 zone has no area (every vertex on a pole) and was left out
{"type":"FeatureCollection","features":[]}
```

Text output shows a zone's vertices as the library gives them, anticlockwise, and CSV carries no
geometry.

## Exit codes

| Code | Meaning |
|---|---|
| 0 | success, including a pipe closed early by the reader, and help |
| 1 | the input named something that is not so (a zone that does not exist, a point out of range, a resolution beyond the grid's finest, a list of sub-zones longer than the limit), a file that cannot be written, or any bad line in batch |
| 2 | a malformed invocation: an unknown command, grid or option, a missing argument |

Every error message says what was wrong and, for a malformed invocation, gives a correct example.

## Differences from `dgg`

For the commands both tools have, the identifiers and coordinates equal `dgg`'s as `dgg` prints
them wherever the library agrees with the engine, checked on a sample of zones and points on all
six grids. Where the library departs from the engine on purpose, the tool follows the library; the
library's README lists every such departure under "Departures from the engine". These concern a
user of `dgg` directly:

- the centroid parent on the aperture-7 grids: `dgg` marks no parent there, since DGGAL answers
  no centroid parent for any zone of those grids, and the tool marks the first parent that is
  itself a centroid child, as both tools do on aperture 3;
- the parents and children on the aperture-7 grids, which omit the null zone, repeated entries
  and identifiers the engine cannot read back, all three found only along two broken seams at
  fine resolutions; everywhere else the lists are `dgg`'s own, in its order;
- the neighbours on the aperture-7 grids, which omit the null zone, the zone itself, repeated
  entries and identifiers the engine cannot read back, all four found only along the same two
  seams;
- the sub-zone at an index on the aperture-7 grids. At a pentagon, at an odd depth,
  `dgg <grid> sub <zone> <index>` answers at some indices a zone that is not the entry of its own
  list: sub-zone 9 of `00` at depth 1 is `006` in its list, and `030` by index. At depth 0, where
  the order is the zone alone, `dgg` by index answers another zone. At the southern polar
  pentagon of an odd resolution (`110`, for instance) at depth 2, `dgg` asked for sub-zone 0 ends
  with a segmentation fault. In all three cases the tool answers the entry of the list;
- the index of a sub-zone on the aperture-7 grids, along the two seams: there the engine's index
  may be a position at which its own list holds another zone (`dgg` gives 6 for
  `000000000000000001` within `00000000000000000`, whose list holds it at position 3), or none for
  a zone that its list holds and its own test of a sub-zone accepts, and the tool, as the
  library, states the position at which the order holds the sub-zone, the first where it holds
  it twice; a zone that the engine's test refuses has none, as in `dgg`;
- the index of a zone among its own sub-zones, at depth 0, on every grid:
  `dgg <grid> index <zone> <zone>` says `sub-zone <zone> not found within parent <zone>`, with
  exit code 1, and the tool answers index 0, since at depth 0 the order is the zone alone;
- the sub-zone order on the aperture-3 grids at some zones on the edges of the rhombi, where the
  engine's own order names a zone twice, names one that is no descendant, or omits one, and the
  tool gives each descendant once;
- the null zone, which the tool answers in the narrow strips where the engine answers an
  identifier that it cannot itself read back;
- the order of a zone's vertices on the aperture-3 grids: `dgg` prints them clockwise, save at the
  two polar pentagons (such as `AA-0-A`), and the tool, as the library, starts from the same first
  vertex and takes the others in reverse, so that every ring runs anticlockwise; the coordinates
  are the same.

What `dgg` has and this tool has not yet:

- the area of a zone and its geodesic extent, and the `level` command;
- `list` and `grid`, the zones within a bounding box;
- `compact` and `decompact`;
- `togeo` and DGGS-JSON;
- other output coordinate systems (`-crs`: the 5x6 plane, the icosahedron net);
- neighbour directions (`dgg` prints `(direction n)` beside each neighbour);
- densified edges in GeoJSON: straight lines in latitude and longitude depart visibly from the
  true edges of cells at resolutions 0 to 2.

## Licence

MIT, as the library's.
