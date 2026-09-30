//! Test-only live oracle over DGGAL v0.0.6's own Rust bindings.
//!
//! Degrees in, degrees out; zone ids are `u64` throughout (DGGAL's `DGGRSZone`), which is
//! why this path does not have pydggal's signed-int64 overflow on Z7 zones with bit 63 set.
//!
//! DGGAL's thread safety is unverified, so every DGGAL object (the `Application`, the
//! `DGGAL` handle, and each grid) is created by, and only ever touched from, one dedicated
//! worker thread that lives for the process's lifetime. Callers on any thread hand the
//! worker a job over a channel and block for the reply; `with()` below is the single choke
//! point through which every public function in this crate runs. This confines DGGAL
//! objects to their creating thread without requiring any `unsafe impl Send`.
//!
//! A job that panics in Rust (an `.unwrap()` in the vendored binding layer, or an
//! assertion of ours) is caught on the worker and the panic is resumed on the caller's own
//! thread, so the failing call still fails with its original message and location, while
//! the worker and every later call in the same test binary survive. This guards against
//! panics raised in Rust that unwind only through Rust frames; it is not a guard against a
//! crash inside eC or DGGAL's own C code (an `abort()`, a failed eC assertion, a signal),
//! which no Rust mechanism can catch, nor against a Rust panic that would have to unwind
//! across an `extern "C"` boundary, which the vendored bindings declare their C call sites
//! as: since Rust 1.81, such a panic aborts the process at that boundary instead of
//! propagating. Either way the process goes down; catch_unwind only ever sees panics that
//! stay within Rust frames, such as the ones this crate's own code and the vendored Rust
//! binding layer can raise.
extern crate dggal;
extern crate dggal_sys;
extern crate ecrt;

use std::any::Any;
use std::collections::HashMap;
use std::panic::{self, AssertUnwindSafe};
use std::sync::mpsc::{self, Sender};
use std::sync::{Mutex, OnceLock};
use std::thread;

use dggal::{CRS, DGGAL, DGGRS};
use ecrt::Application;

pub const IGEO7: &str = "ISEA7H_Z7";
pub const IVEA7H: &str = "IVEA7H_Z7";
pub const RTEA7H: &str = "RTEA7H_Z7";
pub const ISEA3H: &str = "ISEA3H";
pub const IVEA3H: &str = "IVEA3H";
pub const RTEA3H: &str = "RTEA3H";
pub const NULL_ZONE: u64 = u64::MAX;

/// Owned solely by the worker thread spawned in `worker()`: never sent, shared or
/// touched from any other thread.
struct State {
    _app: Application,
    dggal: DGGAL,
    grids: HashMap<String, DGGRS>,
}

/// A unit of work run on the worker thread with exclusive access to `State`.
type Job = Box<dyn FnOnce(&mut State) + Send>;

/// The channel to the worker thread, guarded by a mutex so that concurrent callers (the
/// default test harness runs tests on several threads) queue their jobs one at a time
/// rather than racing to send on a shared `Sender`.
fn worker() -> &'static Mutex<Sender<Job>> {
    static SENDER: OnceLock<Mutex<Sender<Job>>> = OnceLock::new();
    SENDER.get_or_init(|| {
        let (tx, rx) = mpsc::channel::<Job>();
        thread::Builder::new()
            .name("dggal-oracle".to_string())
            .spawn(move || {
                let app = Application::new(&vec!["dggal-oracle".to_string()]);
                let dggal = DGGAL::new(&app);
                let mut state = State {
                    _app: app,
                    dggal,
                    grids: HashMap::new(),
                };
                for job in rx {
                    job(&mut state);
                }
            })
            .expect("dggal-oracle: failed to spawn the DGGAL worker thread");
        Mutex::new(tx)
    })
}

/// Runs `f` against the named grid's `DGGRS` handle on the dedicated worker thread,
/// creating and caching the grid there on first use, and blocks the calling thread for
/// the result. The single chokepoint every public function in this crate goes through.
///
/// A panic raised in Rust while running `f`, or while creating the grid, is caught on the
/// worker (see the module doc) and resumed here, on the caller's own thread, so the
/// failing call fails with its original message and location without taking the worker,
/// or any later call in this process, down with it. The closure passed to `catch_unwind`
/// is asserted unwind-safe: `state.grids` is only ever mutated by the `insert` below,
/// which runs after `DGGRS::new` has already succeeded, so a panic during grid creation or
/// during `f` (which only borrows the grid) never leaves it partially mutated; `f` and
/// `grid` are moved into the closure and consumed there, so neither is observable again
/// after a panic either way.
///
/// The worker runs one job at a time (see `worker`'s doc), so `f` must never call `with`
/// itself: that second call would wait on the reply from the very thread that is running
/// `f`, and neither side would ever proceed. A function that composes several public,
/// `with`-wrapped calls (`engine_can_read` is one) must do so at the top level, one call
/// dispatched and answered before the next begins, never from inside another call's `f`.
fn with<R: Send + 'static>(grid: &str, f: impl FnOnce(&DGGRS) -> R + Send + 'static) -> R {
    let grid = grid.to_string();
    let (reply_tx, reply_rx) = mpsc::channel::<Result<R, Box<dyn Any + Send>>>();
    let job: Job = Box::new(move |state: &mut State| {
        let outcome = panic::catch_unwind(AssertUnwindSafe(|| {
            if !state.grids.contains_key(&grid) {
                let g = DGGRS::new(&state.dggal, &grid)
                    .unwrap_or_else(|e| panic!("DGGAL grid {grid}: {e}"));
                state.grids.insert(grid.clone(), g);
            }
            f(&state.grids[&grid])
        }));
        // If the caller's receiver is already gone (e.g. it panicked while waiting),
        // there is nothing to propagate the result to; move on.
        let _ = reply_tx.send(outcome);
    });
    worker()
        .lock()
        .unwrap_or_else(|e| e.into_inner())
        .send(job)
        .unwrap_or_else(|_| {
            panic!("dggal-oracle: the DGGAL worker thread has stopped, most likely after a panic")
        });
    match reply_rx.recv() {
        Ok(Ok(value)) => value,
        Ok(Err(payload)) => panic::resume_unwind(payload),
        Err(_) => panic!("dggal-oracle: the DGGAL worker thread stopped before replying, most likely after a panic"),
    }
}

fn geo(lat: f64, lon: f64) -> dggal_sys::GeoPoint {
    dggal_sys::GeoPoint {
        lat: lat.to_radians(),
        lon: lon.to_radians(),
    }
}

pub fn zone_from_geo(grid: &str, lat: f64, lon: f64, level: i32) -> u64 {
    with(grid, move |g| {
        g.getZoneFromWGS84Centroid(level, &geo(lat, lon))
    })
}
pub fn zone_from_text(grid: &str, text: &str) -> u64 {
    let text = text.to_string();
    with(grid, move |g| g.getZoneFromTextID(&text))
}
pub fn text_id(grid: &str, zone: u64) -> String {
    with(grid, move |g| g.getZoneTextID(zone))
}
pub fn level(grid: &str, zone: u64) -> i32 {
    with(grid, move |g| g.getZoneLevel(zone))
}
pub fn max_level(grid: &str) -> i32 {
    with(grid, |g| g.getMaxDGGRSZoneLevel())
}
pub fn centroid(grid: &str, zone: u64) -> (f64, f64) {
    with(grid, move |g| {
        let p = g.getZoneWGS84Centroid(zone);
        (p.lat.to_degrees(), p.lon.to_degrees())
    })
}
pub fn vertices(grid: &str, zone: u64) -> Vec<(f64, f64)> {
    with(grid, move |g| {
        g.getZoneWGS84Vertices(zone)
            .iter()
            .map(|p| (p.lat.to_degrees(), p.lon.to_degrees()))
            .collect()
    })
}
/// DGGAL's planar 5x6 CRS, `CRS { ogc, 153456 }`, shared by every accessor below that takes
/// or returns coordinates in it rather than in WGS84.
fn crs_5x6() -> dggal::CRS {
    dggal::CRS!(dggal::ogc, 153_456)
}

/// The zone at `level` whose planar cell, in the 5x6 CRS, contains the point `(x, y)`: the
/// planar census, `getZoneFromCRSCentroid`. Exact on a point that is itself a cell's
/// centroid; a quantiser on any other point.
pub fn zone_from_5x6(grid: &str, level: i32, x: f64, y: f64) -> u64 {
    with(grid, move |g| {
        g.getZoneFromCRSCentroid(level, crs_5x6(), &ecrt::Pointd { x, y })
    })
}
/// The zone's centroid in the planar 5x6 CRS, before any projection to WGS84:
/// `getZoneCRSCentroid`.
pub fn crs_centroid_5x6(grid: &str, zone: u64) -> (f64, f64) {
    with(grid, move |g| {
        let p = g.getZoneCRSCentroid(zone, crs_5x6());
        (p.x, p.y)
    })
}
/// The centroids of the sub-zones `depth` levels below `zone`, in the planar 5x6 CRS and in the
/// engine's sub-zone order, before any quantisation: `getSubZoneCRSCentroids`, which in this CRS
/// hands back the aperture-3 generators' output untouched (`RI3H.ec:676-705`, the
/// `case CRS { ogc, 153456 }: break;` of `:684`); at depth 0, the zone's own centroid
/// (`I3HSubZones.ec:1789-1790`).
///
/// Call it only where the engine answers, below `2^28` sub-zones and at a sub-zone level of 33
/// or less: where it answers null, the binding builds its list from a null pointer, which the
/// debug build of Rust's standard library stops with an abort rather than a panic.
pub fn sub_zone_crs_centroids_5x6(grid: &str, zone: u64, depth: i32) -> Vec<(f64, f64)> {
    with(grid, move |g| {
        g.getSubZoneCRSCentroids(zone, crs_5x6(), depth)
            .iter()
            .map(|p| (p.x, p.y))
            .collect()
    })
}
/// The zone's boundary in the planar 5x6 CRS, before any projection to WGS84:
/// `getZoneCRSVertices`.
pub fn crs_vertices_5x6(grid: &str, zone: u64) -> Vec<(f64, f64)> {
    with(grid, move |g| {
        g.getZoneCRSVertices(zone, crs_5x6())
            .iter()
            .map(|p| (p.x, p.y))
            .collect()
    })
}
/// The zone's centroid child, `getZoneCentroidChild`: on the aperture-3 grids, the one
/// child that shares the parent's own centroid. The null zone where there is none.
pub fn centroid_child(grid: &str, zone: u64) -> u64 {
    with(grid, move |g| g.getZoneCentroidChild(zone))
}
/// Whether `sub` lies anywhere under `zone` in the hierarchy, `zoneHasSubZone`. Unlike the
/// predicates below, this does not reach `getZoneParents`, so `engine_can_read` need not be
/// checked first (`dggrs.ec:151-189`, not overridden on RI3H): it compares levels, then
/// walks down from `sub` through its centroid child, or its own vertices, and planar
/// containment, never up through parents.
pub fn zone_has_sub_zone(grid: &str, zone: u64, sub: u64) -> bool {
    with(grid, move |g| g.zoneHasSubZone(zone, sub))
}
/// Whether `a` is an ancestor of `b` at any depth, `isZoneAncestorOf` with `maxDepth` `0`,
/// the eC's own spelling of "unbounded" (`dggrs.ec:331`, `!maxDepth`).
///
/// Recurses upward from `b` through `getZoneParents` (`dggrs.ec:328-349`), so `b`, and every
/// zone the recursion passes through on the way to `a`, must satisfy `engine_can_read`; call
/// this only once `engine_can_read(grid, b)` holds.
pub fn is_ancestor_of(grid: &str, a: u64, b: u64) -> bool {
    with(grid, move |g| g.isZoneAncestorOf(a, b, 0))
}
/// Whether `a` and `b` share a parent, `areZonesSiblings`.
///
/// Calls `getZoneParents` on both `a` and `b` (`dggrs.ec:360-378`), so both must satisfy
/// `engine_can_read` before this is called.
pub fn are_siblings(grid: &str, a: u64, b: u64) -> bool {
    with(grid, move |g| g.areZonesSiblings(a, b))
}
/// Whether `child` is an immediate child of `parent`, `isZoneImmediateChildOf`, itself
/// `isZoneAncestorOf(parent, child, 1)` (`dggrs.ec:405-408`).
///
/// Reaches `getZoneParents` on `child` exactly as `is_ancestor_of` does on its own second
/// argument, so `child` must satisfy `engine_can_read` before this is called.
pub fn is_immediate_child_of(grid: &str, child: u64, parent: u64) -> bool {
    with(grid, move |g| g.isZoneImmediateChildOf(child, parent))
}
/// True only when the engine's own text of `zone`, read back through `zone_from_text`,
/// returns `zone` again, and `zone`'s level does not exceed 33.
///
/// Guards the calls known to crash or hang the engine outright, a failure no Rust mechanism
/// can catch (see the module documentation): anything that computes a parent, whether through
/// `getZoneParents` directly (`parents`, `centroid_parent`, whose own `parent0`/`getParents`
/// arithmetic reaches it just the same, `RI3H.ec:1203-1223`) or by recursing through it
/// (`is_ancestor_of`, `are_siblings`, `is_immediate_child_of`), dies with `SIGFPE` on the null
/// zone's own parents, on the unreadable north-pole identifiers (root 10 with a non-zero
/// index), and on level-34 identifiers; `neighbors` (`getZoneNeighbors`) on the unreadable
/// north-pole identifiers writes seven entries into the binding's six-slot stack arrays,
/// undefined behaviour in the test process; and `sub_zone_index` (`getSubZoneIndex`) does not
/// return for a sub-zone above resolution 33, the null zone included.
///
/// This function itself never makes any of those calls: it composes `text_id`,
/// `zone_from_text` and `level`, none of which is known to crash or hang it, each a separate
/// top-level dispatch to the worker thread rather than a nested one (see `with`'s
/// documentation on why a call from inside another `with` closure would deadlock).
///
/// Two failure shapes both answer `false`, for different reasons: the null zone reads its
/// own text, `"(null)"`, back to itself, but at level 63, so it fails on level alone; the
/// unreadable north-pole identifiers fail the round trip itself, because `getZoneFromTextID`
/// rejects a root above 9 unless its index is 0 (`RI3H:919-920`), which their text is not.
pub fn engine_can_read(grid: &str, zone: u64) -> bool {
    zone_from_text(grid, &text_id(grid, zone)) == zone && level(grid, zone) <= 33
}

/// The zone's refined boundary in DGGAL's planar 5x6 CRS, `CRS { ogc, 153456 }`, at
/// `edgeRefinement` one.
///
/// This is the one accessor that reaches the routine the eC's own containment test is
/// given. `containsPoint` calls `getBaseRefinedVerticesNoAlloc(false, 1, ...)`
/// (RI7H.ec:956), for which the binding exposes nothing; `getZoneRefinedCRSVertices` in
/// this CRS routes to `getIcoNetRefinedVertices` and thence to the allocating twin of that
/// routine, so the polygon it returns is built by the same boundary tracing. The two are
/// not interchangeable at every zone: on a zone whose boundary crosses a rhombus
/// interruption the allocating side intersperses its intermediate points differently and
/// returns a different number of them. The caller must therefore compare counts before
/// coordinates, and never discard a zone merely because the counts differ.
///
/// Coordinates come back exactly as DGGAL computed them, in 5x6 units, with no rounding
/// or conversion on this side.
pub fn refined_crs_vertices_5x6(grid: &str, zone: u64) -> Vec<(f64, f64)> {
    with(grid, move |g| {
        g.getZoneRefinedCRSVertices(zone, crs_5x6(), 1)
            .iter()
            .map(|p| (p.x, p.y))
            .collect()
    })
}
pub fn neighbors(grid: &str, zone: u64) -> Vec<u64> {
    with(grid, move |g| {
        let mut nb_types = [0i32; 6];
        g.getZoneNeighbors(zone, &mut nb_types)
    })
}
pub fn parents(grid: &str, zone: u64) -> Vec<u64> {
    with(grid, move |g| g.getZoneParents(zone))
}
pub fn children(grid: &str, zone: u64) -> Vec<u64> {
    with(grid, move |g| g.getZoneChildren(zone))
}
pub fn centroid_parent(grid: &str, zone: u64) -> u64 {
    with(grid, move |g| g.getZoneCentroidParent(zone))
}
pub fn is_centroid_child(grid: &str, zone: u64) -> bool {
    with(grid, move |g| g.isZoneCentroidChild(zone))
}
pub fn count_sub_zones(grid: &str, zone: u64, depth: i32) -> u64 {
    with(grid, move |g| g.countSubZones(zone, depth))
}
pub fn first_sub_zone(grid: &str, zone: u64, depth: i32) -> u64 {
    with(grid, move |g| g.getFirstSubZone(zone, depth))
}
pub fn sub_zones(grid: &str, zone: u64, depth: i32) -> Vec<u64> {
    with(grid, move |g| g.getSubZones(zone, depth))
}
/// The engine's `getSubZoneIndex`: the position of `sub` in `parent`'s sub-zone order,
/// or -1 where it finds none. On the aperture-7 grids it answers 0 whenever the two
/// zones share a level, before it consults either (`RI7H.ec:229-230`).
pub fn sub_zone_index(grid: &str, parent: u64, sub: u64) -> i64 {
    with(grid, move |g| g.getSubZoneIndex(parent, sub))
}
