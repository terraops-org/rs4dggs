//! ISEA projection: the shared icosahedral kernel bound to the ISEA vertex assignment
//! (eC `isea` case of `VGCRadialVertex`, `icoVertexGreatCircle.ec:127`). Port of py4dggs
//! `projections/isea.py`.
use crate::interfaces::sealed::{Private, ProjectionPlumbing};
use crate::projections::icovertex::{self, IcoGeometry, RadialVertex};
use crate::{GeoPoint, GridConfig, PlanarPoint, Projection};

/// The icosahedral Snyder equal-area projection, as IGEO7 uses it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct Isea;

impl crate::interfaces::sealed::Sealed for Isea {}

impl Projection for Isea {
    type Geometry = IcoGeometry;
    /// The eC's `ISEAProjection` (`icoVertexGreatCircle.ec:106-109`), which sets its radial vertex.
    fn build_geometry(config: &GridConfig) -> IcoGeometry {
        icovertex::build_geometry(config, RadialVertex::Isea)
    }
    /// The eC's `RI5x6Projection::forward` (`ri5x6.ec:394-459`), which `ISEAProjection` inherits.
    fn forward(geom: &IcoGeometry, lat: f64, lon: f64) -> PlanarPoint {
        icovertex::forward(geom, lat, lon)
    }
    /// The eC's `RI5x6Projection::inverse` (`ri5x6.ec:511-556`), which `ISEAProjection` inherits.
    fn inverse(geom: &IcoGeometry, p: PlanarPoint, odd_grid: bool) -> GeoPoint {
        icovertex::inverse(geom, p, odd_grid)
    }
}

impl ProjectionPlumbing for Isea {
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
    use crate::{GridConfig, PlanarPoint, Projection};

    fn geom() -> <Isea as Projection>::Geometry {
        Isea::build_geometry(&GridConfig::default())
    }

    /// The planar points of DGGAL v0.0.6's own `ISEAProjection`, read from
    /// `RI5x6Projection_forward` with each point set in radians as `deg * (PI / 180)`, the
    /// product [`Isea::forward`] forms. They were taken from py4dggs until the projection
    /// was put into the order the compiled engine performs it; against the engine, two
    /// had moved: Lisbon's abscissa by two units in the last place and its ordinate by
    /// minus two, and the ordinate at (89.9, 45) by one. The face is this port's; the
    /// engine's point lies on it.
    #[test]
    fn forward_anchors_match_dggal() {
        let g = geom();
        let cases = [
            (
                (38.7223, -9.1393),
                (1, 1.9584673695128463, 1.3535585389790992),
            ),
            ((0.0, 0.0), (2, 2.3272842366389384, 2.327284236638937)),
            ((-33.9, 151.2), (18, 3.380067423047124, 4.73575777367487)),
            ((89.9, 45.0), (4, 4.9991188727490625, 4.498276297953374)),
        ];
        for ((lat, lon), (face, x, y)) in cases {
            let p = Isea::forward(&g, lat, lon);
            assert_eq!((p.face, p.x, p.y), (face, x, y), "({lat},{lon})");
        }
    }

    /// The centroid of zone `0064156` from its planar centroid, as DGGAL v0.0.6's
    /// `RI5x6Projection_inverse` gives it on an odd grid, converted to degrees with
    /// `180 / PI`; py4dggs gives the same two doubles.
    #[test]
    fn inverse_anchors_match_dggal() {
        let g = geom();
        let q = Isea::inverse(
            &g,
            PlanarPoint {
                face: -1,
                x: 1.9620991253644315,
                y: 1.3556851311953353,
            },
            // Resolution 5 is odd; the flag matters only at the poles.
            true,
        );
        assert_eq!((q.lat, q.lon), (38.61687542349096, -8.874508932649084));
    }

    #[test]
    fn round_trip_over_seeded_sample() {
        let g = geom();
        let mut s: u64 = 0x9E37_79B9_7F4A_7C15;
        for _ in 0..10_000 {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            let u = (s >> 11) as f64 / (1u64 << 53) as f64;
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            let v = (s >> 11) as f64 / (1u64 << 53) as f64;
            let lat = -89.0 + 178.0 * u;
            let lon = -179.0 + 358.0 * v;
            // The sample stays a degree clear of the poles, where the flag has no effect.
            let q = Isea::inverse(&g, Isea::forward(&g, lat, lon), false);
            assert!(
                (q.lat - lat).abs() < 1e-9 && (q.lon - lon).abs() < 1e-9,
                "({lat},{lon}) -> {q:?}"
            );
        }
    }

    #[test]
    fn get_face_outside_layout_is_minus_one() {
        assert_eq!(crate::projections::icovertex::get_face(-1.0, 3.0), -1);
    }
}
