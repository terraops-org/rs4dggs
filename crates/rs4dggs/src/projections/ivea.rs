//! The IVEA projection: the shared icosahedral equal-area kernel bound to the IVEA vertex
//! assignment (the eC `ivea` case of `VGCRadialVertex`, `icoVertexGreatCircle.ec:137`), in
//! which the triangle's B and C vertices, its angles beta and gamma and its sides AB and AC
//! are exchanged with respect to ISEA.
use crate::projections::icovertex::{self, IcoGeometry, RadialVertex};
use crate::{GeoPoint, GridConfig, PlanarPoint, Projection};

/// The icosahedral vertex-oriented great-circle equal-area projection, as IVEA7H uses
/// it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct Ivea;

impl crate::interfaces::sealed::Sealed for Ivea {}

impl Projection for Ivea {
    type Geometry = IcoGeometry;
    /// The eC's `IVEAProjection` (`icoVertexGreatCircle.ec:102-104`), which keeps the default
    /// radial vertex, `ivea` (`icoVertexGreatCircle.ec:118`).
    fn build_geometry(config: &GridConfig) -> IcoGeometry {
        icovertex::build_geometry(config, RadialVertex::Ivea)
    }
    /// The eC's `RI5x6Projection::forward` (`ri5x6.ec:394-459`), which `IVEAProjection` inherits.
    fn forward(geom: &IcoGeometry, lat: f64, lon: f64) -> PlanarPoint {
        icovertex::forward(geom, lat, lon)
    }
    /// The eC's `RI5x6Projection::inverse` (`ri5x6.ec:511-556`), which `IVEAProjection` inherits.
    fn inverse(geom: &IcoGeometry, p: PlanarPoint, odd_grid: bool) -> GeoPoint {
        icovertex::inverse(geom, p, odd_grid)
    }
    /// The eC's `RI5x6Projection::inverse` (`ri5x6.ec:511-556`) as it answers: in radians, and
    /// `None` where it returns `false`.
    fn inverse_radians(geom: &IcoGeometry, p: PlanarPoint, odd_grid: bool) -> Option<(f64, f64)> {
        icovertex::inverse_radians(geom, p, odd_grid)
    }
    /// The authalic latitude where the geometry converts latitudes, as `forward` does
    /// (`ri5x6.ec:252-255`, `:397`), and the latitude itself where it does not.
    fn sphere_latitude(geom: &IcoGeometry, lat: f64) -> f64 {
        icovertex::sphere_latitude(geom, lat)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{GridConfig, Projection};

    /// The planar points of DGGAL v0.0.6's own `IVEAProjection`, read from
    /// `RI5x6Projection_forward` with each point set in radians as `deg * (PI / 180)`, the
    /// product [`Ivea::forward`] forms. They were taken from py4dggs until the projection
    /// was put into the order the compiled engine performs it; against the engine, one
    /// had moved, the ordinate at (89.9, 45), by minus one unit in the last place. The
    /// face is this port's; the engine's point lies on it.
    #[test]
    fn forward_anchors_match_dggal() {
        let g = Ivea::build_geometry(&GridConfig::default());
        let cases = [
            (
                (38.7223, -9.1393),
                (1, 1.9589788258393153, 1.3526515920175908),
            ),
            ((0.0, 0.0), (2, 2.3258606230304615, 2.32586062303046)),
            ((-33.9, 151.2), (18, 3.380488933513674, 4.7367089221166685)),
            ((89.9, 45.0), (4, 4.999118287978329, 4.498276931171835)),
        ];
        for ((lat, lon), (face, x, y)) in cases {
            let p = Ivea::forward(&g, lat, lon);
            assert_eq!((p.face, p.x, p.y), (face, x, y), "({lat},{lon})");
        }
    }
}
