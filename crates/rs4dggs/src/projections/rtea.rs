//! The RTEA projection: the shared icosahedral equal-area kernel bound to the RTEA vertex
//! assignment (the eC `rtea` case of `VGCRadialVertex`, `icoVertexGreatCircle.ec:147`), in
//! which the sub-triangle's edge-midpoint vertex becomes the radial vertex A, the three
//! sides AB, AC and BC are cyclically permuted with respect to ISEA, and the right angle
//! moves from A to B.
use crate::interfaces::sealed::{Private, ProjectionPlumbing};
use crate::projections::icovertex::{self, IcoGeometry, RadialVertex};
use crate::{GeoPoint, GridConfig, PlanarPoint, Projection};

/// The rhombic-triacontahedral equal-area projection, as RTEA7H uses it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct Rtea;

impl crate::interfaces::sealed::Sealed for Rtea {}

impl Projection for Rtea {
    type Geometry = IcoGeometry;
    /// The eC's `RTEAProjection` (`icoVertexGreatCircle.ec:111-114`), which sets its radial vertex.
    fn build_geometry(config: &GridConfig) -> IcoGeometry {
        icovertex::build_geometry(config, RadialVertex::Rtea)
    }
    /// The eC's `RI5x6Projection::forward` (`ri5x6.ec:394-459`), which `RTEAProjection` inherits.
    fn forward(geom: &IcoGeometry, lat: f64, lon: f64) -> PlanarPoint {
        icovertex::forward(geom, lat, lon)
    }
    /// The eC's `RI5x6Projection::inverse` (`ri5x6.ec:511-556`), which `RTEAProjection` inherits.
    fn inverse(geom: &IcoGeometry, p: PlanarPoint, odd_grid: bool) -> GeoPoint {
        icovertex::inverse(geom, p, odd_grid)
    }
}

impl ProjectionPlumbing for Rtea {
    /// The eC's `RI5x6Projection::inverse` (`ri5x6.ec:511-556`) as it answers: in radians, and
    /// `None` where it returns `false`.
    fn inverse_radians(
        _: Private,
        geom: &IcoGeometry,
        p: PlanarPoint,
        odd_grid: bool,
    ) -> Option<(f64, f64)> {
        icovertex::inverse_radians(geom, p, odd_grid)
    }
    /// The authalic latitude where the geometry converts latitudes, as `forward` does
    /// (`ri5x6.ec:252-255`, `:397`), and the latitude itself where it does not.
    fn sphere_latitude(_: Private, geom: &IcoGeometry, lat: f64) -> f64 {
        icovertex::sphere_latitude(geom, lat)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{GridConfig, Projection};

    /// The planar points of DGGAL v0.0.6's own `RTEAProjection`, read from
    /// `RI5x6Projection_forward` with each point set in radians as `deg * (PI / 180)`, the
    /// product [`Rtea::forward`] forms. They were taken from py4dggs until the projection
    /// was put into the order the compiled engine performs it; against the engine, one
    /// had moved, the abscissa at (-33.9, 151.2), by one unit in the last place. The face
    /// is this port's; the engine's point lies on it.
    #[test]
    fn forward_anchors_match_dggal() {
        let g = Rtea::build_geometry(&GridConfig::default());
        let cases = [
            (
                (38.7223, -9.1393),
                (1, 1.9595370656698252, 1.3528965187978825),
            ),
            ((0.0, 0.0), (2, 2.321451057055118, 2.321451057055117)),
            ((-33.9, 151.2), (18, 3.3801191698750452, 4.735633200562001)),
            ((89.9, 45.0), (4, 4.99914032858217, 4.498290844165996)),
        ];
        for ((lat, lon), (face, x, y)) in cases {
            let p = Rtea::forward(&g, lat, lon);
            assert_eq!((p.face, p.x, p.y), (face, x, y), "({lat},{lon})");
        }
    }
}
