//! The text layouts: two-space indent, labels padded to twelve columns.

use std::io::{self, Write};

use rs4dggs::GeoPoint;

use crate::commands::{Card, Named, coord};

fn point(p: &GeoPoint, precision: Option<usize>) -> String {
    format!("{}, {}", coord(p.lat, precision), coord(p.lon, precision))
}

fn named(list: &[Named], mark: &str) -> String {
    let items: Vec<String> = list
        .iter()
        .map(|n| {
            if n.mark {
                format!("{} ({mark})", n.text)
            } else {
                n.text.clone()
            }
        })
        .collect();
    format!("{}: {}", list.len(), items.join(", "))
}

/// The card of one zone.
pub fn card(out: &mut dyn Write, c: &Card, precision: Option<usize>) -> io::Result<()> {
    writeln!(
        out,
        "{} zone {}  (integer {}, hex {:016X})",
        c.grid, c.text, c.id.0, c.id.0
    )?;
    let shape = if c.pentagon { "pentagon" } else { "hexagon" };
    writeln!(out, "  resolution  {}, {shape}", c.resolution)?;
    writeln!(
        out,
        "  centroid    {}   (lat, lon)",
        point(&c.centroid, precision)
    )?;
    if c.parents.is_empty() {
        writeln!(out, "  parents     none")?;
    } else {
        writeln!(
            out,
            "  parents     {}",
            named(&c.parents, "centroid parent")
        )?;
    }
    if c.children.is_empty() {
        writeln!(out, "  children    none")?;
    } else {
        writeln!(out, "  children    {}", named(&c.children, "centre"))?;
    }
    writeln!(
        out,
        "  neighbours  {}: {}",
        c.neighbours.len(),
        c.neighbours.join(", ")
    )?;
    writeln!(out, "  vertices    {}:", c.vertices.len())?;
    for v in &c.vertices {
        writeln!(out, "              {}", point(v, precision))?;
    }
    Ok(())
}

/// The zones containing a point, one per resolution.
pub fn zones_at(
    out: &mut dyn Write,
    grid: &str,
    (lat, lon): (f64, f64),
    rows: &[(u8, String)],
    precision: Option<usize>,
) -> io::Result<()> {
    writeln!(
        out,
        "{grid} zones containing {}, {}",
        coord(lat, precision),
        coord(lon, precision)
    )?;
    writeln!(out, "  res  zone")?;
    for (r, z) in rows {
        writeln!(out, "{r:>5}  {z}")?;
    }
    Ok(())
}

pub fn neighbours(out: &mut dyn Write, grid: &str, zone: &str, list: &[String]) -> io::Result<()> {
    writeln!(out, "{grid} zone {zone}: {} neighbours", list.len())?;
    for n in list {
        writeln!(out, "  {n}")?;
    }
    Ok(())
}

pub fn disk_head(out: &mut dyn Write, grid: &str, k: u32, zone: &str) -> io::Result<()> {
    writeln!(out, "{grid} disk of radius {k} about {zone}")
}

pub fn ring(out: &mut dyn Write, i: usize, zones: &[String]) -> io::Result<()> {
    writeln!(out, "  ring {i} ({}): {}", zones.len(), zones.join(", "))
}
