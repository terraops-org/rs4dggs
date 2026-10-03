//! Walk the hierarchy of a zone of IGEO7: its primary parents up to the base cell, where the
//! engine names them, each an ancestor, and its children. A zone that is no centroid child
//! has a second parent, which this walk does not follow.
//!
//! ```text
//! cargo run --example hierarchy -- [<zone id>]    (default: 0064156, Lisbon at resolution 5)
//! ```

use std::process::ExitCode;

use rs4dggs::{ZoneId, igeo7};

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let text = match args.as_slice() {
        [] => "0064156".to_string(),
        [text] => text.clone(),
        _ => {
            eprintln!("usage: hierarchy [<zone id>]");
            return ExitCode::FAILURE;
        }
    };
    let zone = match igeo7().zone_from_text(&text) {
        Ok(zone) => zone,
        Err(e) => {
            eprintln!("error: {e}");
            return ExitCode::FAILURE;
        }
    };
    // The null zone's own text reads back as the null zone, which has no place in the
    // hierarchy: refuse it rather than describe it.
    if zone.id() == ZoneId::NULL {
        eprintln!("error: {text:?} names no zone");
        return ExitCode::FAILURE;
    }
    let kind = if zone.is_pentagon() {
        "a pentagon"
    } else {
        "a hexagon"
    };
    println!(
        "{} at resolution {}, {kind}",
        zone.text_id(),
        zone.resolution()
    );
    let mut ancestor = zone.parent();
    while let Some(a) = ancestor {
        println!(
            "  ancestor {} at resolution {}",
            a.text_id(),
            a.resolution()
        );
        ancestor = a.parent();
    }
    let children = zone.children();
    if children.is_empty() {
        println!(
            "  no children: resolution {} is the finest",
            zone.resolution()
        );
    } else {
        println!("  {} children:", children.len());
        for child in &children {
            println!("    {}", child.text_id());
        }
    }
    ExitCode::SUCCESS
}
