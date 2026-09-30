//! A zone and its neighbours as a GeoJSON FeatureCollection, written with no dependency.
//!
//! ```text
//! cargo run --example geojson -- [<lat> <lon> <resolution>]   (default: Lisbon at resolution 5)
//! ```
//!
//! Each ring is closed by repeating its first vertex, and keeps the crate's anticlockwise
//! order, which is what GeoJSON asks of an exterior ring. A cell is written as the crate
//! returns it, every longitude between -180 and 180: so consecutive vertices of a cell that
//! crosses the antimeridian can differ by almost 360 degrees, and the ring of a cell about a
//! pole runs through every longitude. That is valid JSON, but a renderer must split or shift
//! such a ring before drawing it.
//!
//! A point in no cell is refused, rather than written: the crate answers the null zone there,
//! whose geometry is the point (0, 0) and which no map should draw.

use std::process::ExitCode;

use rs4dggs::{GeoPoint, ZoneId, igeo7};

const USAGE: &str = "usage: geojson [<lat> <lon> <resolution>]";

fn ring(vertices: &[GeoPoint]) -> String {
    let mut points: Vec<String> = vertices
        .iter()
        .map(|p| format!("[{},{}]", p.lon, p.lat))
        .collect();
    points.push(points[0].clone());
    format!("[[{}]]", points.join(","))
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let (lat, lon, res): (f64, f64, u8) = match args.as_slice() {
        [] => (38.7223, -9.1393, 5),
        [lat, lon, res] => match (lat.parse(), lon.parse(), res.parse()) {
            (Ok(lat), Ok(lon), Ok(res)) => (lat, lon, res),
            _ => {
                eprintln!("{USAGE}");
                return ExitCode::FAILURE;
            }
        },
        _ => {
            eprintln!("{USAGE}");
            return ExitCode::FAILURE;
        }
    };
    let zone = match igeo7().zone_from_geo(lat, lon, res) {
        Ok(zone) => zone,
        Err(e) => {
            eprintln!("error: {e}");
            return ExitCode::FAILURE;
        }
    };
    if zone.id() == ZoneId::NULL {
        eprintln!("error: ({lat}, {lon}) lies in no cell at resolution {res}");
        return ExitCode::FAILURE;
    }
    let features: Vec<String> = std::iter::once(zone)
        .chain(zone.neighbors())
        .map(|z| {
            format!(
                r#"{{"type":"Feature","properties":{{"id":"{}","resolution":{}}},"geometry":{{"type":"Polygon","coordinates":{}}}}}"#,
                z.text_id(),
                z.resolution(),
                ring(&z.vertices())
            )
        })
        .collect();
    println!(
        r#"{{"type":"FeatureCollection","features":[{}]}}"#,
        features.join(",")
    );
    ExitCode::SUCCESS
}
