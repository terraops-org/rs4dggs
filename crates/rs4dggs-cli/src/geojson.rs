//! A GeoJSON FeatureCollection, written feature by feature.

use std::io::{self, Write};

use rs4dggs::{AnyGrid, ZoneId};
use rs4dggs_ogc::geometry::{Shape, clean};

use crate::commands::coord;

pub struct Collection<'a> {
    out: &'a mut dyn Write,
    precision: Option<usize>,
    centroids: bool,
    started: bool,
    features: u64,
    /// Zones with no cell to draw: the null zone.
    in_no_cell: u64,
    /// Zones with a coordinate that is not finite, which GeoJSON cannot hold.
    not_finite: u64,
    /// Zones whose every vertex lies on a pole: they have no area to draw.
    no_area: u64,
}

impl<'a> Collection<'a> {
    pub fn new(out: &'a mut dyn Write, precision: Option<usize>, centroids: bool) -> Self {
        Collection {
            out,
            precision,
            centroids,
            started: false,
            features: 0,
            in_no_cell: 0,
            not_finite: 0,
            no_area: 0,
        }
    }

    /// Opens the collection, once. It is opened with the first feature, or at the end, so that an
    /// input refused outright writes nothing.
    pub fn begin(&mut self) -> io::Result<()> {
        if !self.started {
            self.started = true;
            write!(self.out, r#"{{"type":"FeatureCollection","features":["#)?;
        }
        Ok(())
    }

    fn position(&self, p: [f64; 2]) -> String {
        format!(
            "[{},{}]",
            coord(p[0], self.precision),
            coord(p[1], self.precision)
        )
    }

    fn ring(&self, ring: &[[f64; 2]]) -> String {
        let points: Vec<String> = ring.iter().map(|p| self.position(*p)).collect();
        format!("[{}]", points.join(","))
    }

    /// The geometry of `id`: its polygon, or with `-centroids` its centroid; `None` if a
    /// coordinate of it is not finite.
    fn geometry(&self, grid: &AnyGrid, id: ZoneId) -> Option<String> {
        if self.centroids {
            let c = grid.centroid(id);
            let p = [c.lon, c.lat];
            return p
                .iter()
                .all(|x| x.is_finite())
                .then(|| format!(r#"{{"type":"Point","coordinates":{}}}"#, self.position(p)));
        }
        let shape = clean(&grid.vertices(id));
        let all: Vec<&[f64; 2]> = match &shape {
            Shape::Polygon(r) => r.iter().collect(),
            Shape::MultiPolygon(rs) => rs.iter().flatten().collect(),
        };
        if all.iter().any(|p| !(p[0].is_finite() && p[1].is_finite())) {
            return None;
        }
        Some(match &shape {
            Shape::Polygon(r) => {
                format!(r#"{{"type":"Polygon","coordinates":[{}]}}"#, self.ring(r))
            }
            Shape::MultiPolygon(rs) => {
                let polygons: Vec<String> =
                    rs.iter().map(|r| format!("[{}]", self.ring(r))).collect();
                format!(
                    r#"{{"type":"MultiPolygon","coordinates":[{}]}}"#,
                    polygons.join(",")
                )
            }
        })
    }

    /// Writes the feature of zone `id`, with `extra` after its properties: members already
    /// written, each led by a comma (`,"ring":1`), which say what the feature answers. The null
    /// zone, and a zone with a coordinate that is not finite, are counted and left out.
    pub fn zone(&mut self, grid: &AnyGrid, id: ZoneId, extra: &str) -> io::Result<()> {
        if id == ZoneId::NULL {
            self.in_no_cell += 1;
            return Ok(());
        }
        if !self.centroids && grid.vertices(id).iter().all(|v| v.lat.abs() == 90.0) {
            self.no_area += 1;
            return Ok(());
        }
        let Some(geometry) = self.geometry(grid, id) else {
            self.not_finite += 1;
            return Ok(());
        };
        self.begin()?;
        if self.features > 0 {
            write!(self.out, ",")?;
        }
        self.features += 1;
        write!(
            self.out,
            r#"{{"type":"Feature","properties":{{"zone":"{}","grid":"{}","#,
            grid.text_id(id),
            grid.name(),
        )?;
        write!(
            self.out,
            r#""resolution":{}{extra}}},"geometry":{geometry}}}"#,
            grid.resolution(id),
        )
    }

    /// Closes the collection, and says on `err` what was left out.
    pub fn end(&mut self, err: &mut dyn Write) -> io::Result<()> {
        self.begin()?;
        writeln!(self.out, "]}}")?;
        match self.in_no_cell {
            0 => {}
            1 => writeln!(err, "1 input lies in no cell and was left out")?,
            n => writeln!(err, "{n} inputs lie in no cell and were left out")?,
        }
        match self.no_area {
            0 => {}
            1 => writeln!(
                err,
                "1 zone has no area (every vertex on a pole) and was left out"
            )?,
            n => writeln!(
                err,
                "{n} zones have no area (every vertex on a pole) and were left out"
            )?,
        }
        match self.not_finite {
            0 => {}
            1 => writeln!(
                err,
                "1 zone with a coordinate that is not finite was left out"
            )?,
            n => writeln!(
                err,
                "{n} zones with a coordinate that is not finite were left out"
            )?,
        }
        Ok(())
    }
}
