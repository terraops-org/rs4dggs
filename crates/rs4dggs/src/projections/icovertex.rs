//! The shared icosahedral great-circle equal-area projection kernel.
//!
//! This is the projection half of the ISEA7H / IVEA7H family: it maps a geographic
//! latitude and longitude onto a point `(face, x, y)` in DGGAL's oblique 5x6 grid, the
//! planar layout that the hexagonal topology then quantises, and back again. It is
//! equal-area by construction, following Snyder's 1992 great-circle vector method as
//! refined in DGGAL's `icoVertexGreatCircle.ec`.
//!
//! DGGAL's `SliceAndDiceGreatCircleIcosahedralProjection` is a single kernel with a
//! `radialVertex` enumeration (`isea`, `ivea`, `rtea`) that selects the per-variant
//! spherical-triangle edge angles and the A-versus-BC sub-triangle assignment. This
//! module mirrors that arrangement: every function here is shared, and the variant is
//! carried by [`RadialVertex`] on the geometry. The variant-specific state is two
//! things, exactly as in the eC property setter (`icoVertexGreatCircle.ec:120-164`): the
//! triangle constants computed by `variant_consts`, and the `b_is_a` flag, which XORs
//! `radial == Ivea` into the sub-triangle test.
//!
//! Provenance, cited inline block by block:
//!   - `projections/ri5x6.ec`, for the orientation, the twelve icosahedron vertices, the
//!     twenty-face 5x6 layout and `getFace`.
//!   - `projections/icoVertexGreatCircle.ec`, for the equal-area forward and inverse
//!     vector projection (`sqrtOneMinusDotOver2`, the midpoint-cross `sphericalTriArea`,
//!     the barycentric sub-triangle decomposition and the `radialVertex` selection).
//!   - `projections/authalic.ec`, for the geodetic to authalic latitude conversion
//!     (Karney 2022, a Clenshaw series over the WGS84 ellipsoid).
//!   - `projections/barycentric5x6.ec` is a different projection class (Goldberg); the
//!     barycentric, spherical-area and slerp helpers ported here live in `ri5x6.ec`,
//!     which is what the citations below name.
//!
//! Ported from py4dggs's `projections/icovertex.py`. A handful of py4dggs expressions
//! substitute an analytically equal form for a runtime computation of the eC; this port
//! restores the eC form in each of those places and marks it with a `// faithful:`
//! comment giving the eC file and line. The same comments mark every place where the
//! two sources genuinely disagree.
//!
//! One class of departure goes further, and against the eC's own text, as in the
//! hexagonal aperture-7 topology, [`crate::topologies::HexA7`]. DGGAL v0.0.6 is
//! compiled with `-O2 -ffast-math`
//! (`Makefile.dggal:133`). The four kernels of `icoVertexGreatCircle.ec` carry the
//! `-fno-unsafe-math-optimizations` attribute, but gcc will not inline a function across
//! a mismatch of floating-point options, so the helpers they call, `Vector3D::Normalize`,
//! the `length` property, `cartesianToBary`, `sphericalTriArea` and
//! `sqrtOneMinusDotOver2` among them, stay out of line and were compiled with the rest of
//! the file's options; so were the constructor, the authalic series, `getFace` and
//! `fixPoles`, which carry no attribute at all. Where gcc re-associated one of those
//! operations in a way that can change an answer, this module performs the order read
//! out of the shipped `libdggal.so`, BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`,
//! and each such site carries a `// faithful:` comment quoting its address and the eC
//! line it replaces. The longitude wrap that ends the inverse, the eC's `wrapLon`, which
//! also carries no attribute, is performed in its compiled form too, although neither of
//! its branches is entered on the registered grids and no answer there turns on it (see
//! `wrap_lon`); and so is its fellow `wrapLonAt`, with which the engine unwraps the refined
//! boundary of a zone about its centroid (see `wrap_lon_at`). Elsewhere the eC's source text
//! decides; the few sites that keep another
//! form, py4dggs's among them, say so, each with the evidence that it changes no answer
//! on the registered grids. With these orders every cached geometry field is
//! bit-identical to the engine's own, and so are the forward projection at each of the
//! 3,184 geographic stations of the census the comments below cite, and the inverse at
//! each of the 654 of its 1,831 planar stations that lie on a face of the layout. All of
//! it was measured on x86-64 Linux with glibc, where the engine and this crate share one
//! `libm`; that says nothing of wasm32 or of any target whose floating-point functions
//! come from another library.

use crate::math;
use crate::{GeoPoint, GridConfig, PlanarPoint};

/// The golden ratio, eC `phi` (`ri5x6.ec:79`). A square root, so a function rather
/// than a constant: every transcendental in this crate goes through [`crate::math`].
/// Only the tests read it, since the three arcs below that the eC derives from it are
/// pinned as the library holds them, and so it is compiled for the tests alone.
#[cfg(test)]
pub(crate) fn phi() -> f64 {
    (1.0 + math::sqrt(5.0)) / 2.0
}

// --- Authalic latitude conversion (Karney 2022) over WGS84 (authalic.ec) ---
/// WGS84 semi-major axis, eC `wgs84Major`.
pub(crate) const WGS84_A: f64 = 6378137.0;
/// WGS84 semi-minor axis, eC `wgs84Minor`.
pub(crate) const WGS84_B: f64 = 6356752.314245179;

/// Geodetic to authalic series coefficients, eC `Cxiphi` (`authalic.ec:19`).
#[rustfmt::skip]
pub(crate) const CXIPHI: [f64; 21] = [
    -4.0 / 3.0,    -4.0 / 45.0,    88.0 / 315.0,       538.0 / 4725.0,     20824.0 / 467775.0,      -44732.0 / 2837835.0,
                   34.0 / 45.0,     8.0 / 105.0,     -2482.0 / 14175.0,   -37192.0 / 467775.0,   -12467764.0 / 212837625.0,
                                -1532.0 / 2835.0,     -898.0 / 14175.0,    54968.0 / 467775.0,   100320856.0 / 1915538625.0,
                                                      6007.0 / 14175.0,    24496.0 / 467775.0,    -5884124.0 / 70945875.0,
                                                                          -23356.0 / 66825.0,      -839792.0 / 19348875.0,
                                                                                                  570284222.0 / 1915538625.0,
];

/// Authalic to geodetic series coefficients, eC `Cphixi` (`authalic.ec:30`).
#[rustfmt::skip]
pub(crate) const CPHIXI: [f64; 21] = [
    4.0 / 3.0,  4.0 / 45.0,   -16.0 / 35.0,  -2582.0 / 14175.0,  60136.0 / 467775.0,    28112932.0 / 212837625.0,
               46.0 / 45.0,  152.0 / 945.0, -11966.0 / 14175.0, -21016.0 / 51975.0,    251310128.0 / 638512875.0,
                            3044.0 / 2835.0,  3802.0 / 14175.0, -94388.0 / 66825.0,     -8797648.0 / 10945935.0,
                                              6059.0 / 4725.0,  41072.0 / 93555.0,   -1472637812.0 / 638512875.0,
                                                               768272.0 / 467775.0,   -455935736.0 / 638512875.0,
                                                                                      4210684958.0 / 1915538625.0,
];

/// Face vertex indices, eC `icoIndices[20][3]` (`ri5x6.ec:17`).
#[rustfmt::skip]
pub(crate) const FACE_VERTICES: [[usize; 3]; 20] = [
    [0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 5], [0, 5, 1],
    [6, 2, 1], [7, 3, 2], [8, 4, 3], [9, 5, 4], [10, 1, 5],
    [2, 6, 7], [3, 7, 8], [4, 8, 9], [5, 9, 10], [1, 10, 6],
    [11, 7, 6], [11, 8, 7], [11, 9, 8], [11, 10, 9], [11, 6, 10],
];

/// The 5x6 planar coordinates of each face's three vertices, eC `vertices5x6[20][3]`
/// (`ri5x6.ec:48`).
#[rustfmt::skip]
pub(crate) const FACE_5X6: [[[f64; 2]; 3]; 20] = [
    [[1.0, 0.0], [0.0, 0.0], [1.0, 1.0]], [[2.0, 1.0], [1.0, 1.0], [2.0, 2.0]], [[3.0, 2.0], [2.0, 2.0], [3.0, 3.0]], [[4.0, 3.0], [3.0, 3.0], [4.0, 4.0]], [[5.0, 4.0], [4.0, 4.0], [5.0, 5.0]],
    [[0.0, 1.0], [1.0, 1.0], [0.0, 0.0]], [[1.0, 2.0], [2.0, 2.0], [1.0, 1.0]], [[2.0, 3.0], [3.0, 3.0], [2.0, 2.0]], [[3.0, 4.0], [4.0, 4.0], [3.0, 3.0]], [[4.0, 5.0], [5.0, 5.0], [4.0, 4.0]],
    [[1.0, 1.0], [0.0, 1.0], [1.0, 2.0]], [[2.0, 2.0], [1.0, 2.0], [2.0, 3.0]], [[3.0, 3.0], [2.0, 3.0], [3.0, 4.0]], [[4.0, 4.0], [3.0, 4.0], [4.0, 5.0]], [[5.0, 5.0], [4.0, 5.0], [5.0, 6.0]],
    [[0.0, 2.0], [1.0, 2.0], [0.0, 1.0]], [[1.0, 3.0], [2.0, 3.0], [1.0, 2.0]], [[2.0, 4.0], [3.0, 4.0], [2.0, 3.0]], [[3.0, 5.0], [4.0, 5.0], [3.0, 4.0]], [[4.0, 6.0], [5.0, 6.0], [4.0, 5.0]],
];

/// The scalar triple product `A . (B x C)` of any sub-triangle frame, which is the
/// same for all 120 of them: the volume of the parallelepiped they span.
///
/// The eC computes it at run time (`icoVertexGreatCircle.ec:203`), and so does this
/// port; py4dggs substitutes this closed form instead. It is kept here because the
/// equality of the two is worth pinning, and `parallelepiped_v_matches_every_frame`
/// does precisely that. Nothing outside that test reads it, so it is compiled for the
/// tests alone.
#[cfg(test)]
pub(crate) fn parallelepiped_v() -> f64 {
    math::sqrt((5.0 - 2.0 * math::sqrt(5.0)) / 15.0)
}

/// The area of one icosahedral sub-triangle: 6 degrees, that is pi/30 radians. The eC
/// hard-codes it as `Degrees { 6 }` rather than calling `sphericalTriArea`
/// (`icoVertexGreatCircle.ec:179`, `:243` and `:365`).
// faithful: DGGAL v0.0.6 as built with `-O2 -ffast-math` (`Makefile.dggal:133`), BuildID
// `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, loads this area from `.rodata` at `0x4a450`
// at all four of its uses (`0x3f6e4` in the forward, `0x3eb10`, `0x3f084` and `0x3f098`
// in the inverse), and the value there is `0x3fbacee9f37bebd7`, one unit in the last
// place above `6.0 * (PI / 180.0)`, which is the product every conversion at run time
// gives (ecrt's `Degrees` to `Radians` property answers `0x3fbacee9f37bebd6` for 6). The
// value is the double nearest to `0.1047197551196598`, the angle written out to sixteen
// decimal places, which is how the eC compiler emits a `Degrees` literal into the C it
// hands to gcc; the text is the inference, the bit pattern is the fact, and so the bit
// pattern is what is pinned here.
pub(crate) const SDT_AREA: f64 = f64::from_bits(0x3fba_cee9_f37b_ebd7);

// --------------------------------------------------------------------------- //
// Per-variant spherical-triangle constants (icoVertexGreatCircle.ec:120-164)
// --------------------------------------------------------------------------- //

/// Which vertex of the base sub-triangle the great circles radiate from: the eC's
/// `VGCRadialVertex` (`icoVertexGreatCircle.ec:100`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RadialVertex {
    /// Snyder's 1992 projection on the icosahedron; the centre becomes the radial
    /// vertex A.
    Isea,
    /// The vertex-oriented projection of the paper; the corner becomes A.
    Ivea,
    /// Snyder extended to the rhombic triacontahedron.
    Rtea,
}

impl RadialVertex {
    /// The `(va, vb, vc)` indices into the sub-triangle vertex array
    /// `[mid(0), corner(1), centre(2)]`. The eC then assigns `a = vb`,
    /// `b = bIsA ? va : vc` and `c = bIsA ? vc : va`
    /// (`icoVertexGreatCircle.ec:128`, `:138`, `:148`).
    pub(crate) fn vertex_order(self) -> (usize, usize, usize) {
        match self {
            RadialVertex::Isea => (0, 2, 1),
            RadialVertex::Ivea => (0, 1, 2),
            RadialVertex::Rtea => (1, 0, 2),
        }
    }
}

/// The spherical-triangle sides and angles of one variant, with the trigonometry the
/// kernel reuses. Mirrors the fields the eC `radialVertex` property setter assigns
/// (`icoVertexGreatCircle.ec:120-158`); only the fields this kernel actually reads are
/// carried.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VariantConsts {
    /// Side AB.
    pub(crate) ab: f64,
    /// Side AC.
    pub(crate) ac: f64,
    /// Side BC.
    pub(crate) bc: f64,
    /// Angle at A.
    pub(crate) alpha: f64,
    /// Angle at B.
    pub(crate) beta: f64,
    /// Angle at C.
    pub(crate) gamma: f64,
    pub(crate) cos_ab: f64,
    pub(crate) cos_ac: f64,
    pub(crate) sin_ac: f64,
    pub(crate) cos_bc: f64,
}

// The three arcs of the icosahedron that every sub-triangle is built from, the sides
// the eC's `radialVertex` setter assigns to `AB`, `AC` and `BC` in an order that depends
// on the variant (`icoVertexGreatCircle.ec:127-163`).
// faithful: DGGAL v0.0.6 as built with `-O2 -ffast-math` (`Makefile.dggal:133`), BuildID
// `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, the `radialVertex` setter at `0x3e410`, which
// replaces the `acos`, `sqrt` and `atan` the eC writes at `:132-134`, `:142-144` and
// `:152-154`. gcc folded each of the three cases whole at compile time, `phi` being a
// `define` (`ri5x6.ec:79`), and the setter stores the sides from `.rodata`: ISEA's from
// `0x49e28` to `0x49e38`, IVEA's from `0x49e58` to `0x49e68`, RTEA's from `0x49e88` to
// `0x49e98`. The engine therefore calls no function of `libm` for them, and neither does
// this crate: each side is the bit pattern the library holds. On x86-64 Linux with glibc
// the expressions the eC writes give the same three doubles at run time, so no answer
// changes there; the pinned pattern keeps them so on any target. The same setter folds the
// angles, and the sines, cosines and half-angle tangent of the sides, into constants as
// well (`0x49e10` to `0x49ea8`, `0x4a3f0` to `0x4a428`). The angles are formed below as
// products with `DEG2RAD`, which call no `libm` and equal the library's constants; of the
// rest this crate reads three cosines and one sine, which it computes at run time through
// [`crate::math`] and which equal the library's constants on glibc.
/// `acos(sqrt((phi + 1) / 3))`, the arc from the centre of a face to the midpoint of one
/// of its edges.
pub(crate) const CENTRE_TO_MIDPOINT: f64 = f64::from_bits(0x3fd7_59ed_d04f_68da);
/// `atan(1 / phi)`, the arc from the midpoint of an edge to either of its ends: half an
/// edge.
pub(crate) const MIDPOINT_TO_CORNER: f64 = f64::from_bits(0x3fe1_b6e1_92eb_be44);
/// `atan(2 / (phi * phi))`, the arc from the centre of a face to one of its corners.
pub(crate) const CENTRE_TO_CORNER: f64 = f64::from_bits(0x3fe4_e01e_2d74_e77d);

/// Select the spherical-triangle sides and angles for `radial`, the eC `radialVertex`
/// cases (`icoVertexGreatCircle.ec:127-163`). `Isea` reproduces the historical module
/// constants exactly, `Ivea` swaps AB with AC, and `Rtea` permutes all three sides and
/// moves the right angle from A to B.
pub(crate) fn variant_consts(radial: RadialVertex) -> VariantConsts {
    let (ab, ac, bc, alpha, beta, gamma) = match radial {
        RadialVertex::Isea => (
            CENTRE_TO_MIDPOINT,
            MIDPOINT_TO_CORNER,
            CENTRE_TO_CORNER,
            90.0 * math::DEG2RAD,
            60.0 * math::DEG2RAD,
            36.0 * math::DEG2RAD,
        ),
        RadialVertex::Ivea => (
            MIDPOINT_TO_CORNER,
            CENTRE_TO_MIDPOINT,
            CENTRE_TO_CORNER,
            90.0 * math::DEG2RAD,
            36.0 * math::DEG2RAD,
            60.0 * math::DEG2RAD,
        ),
        RadialVertex::Rtea => (
            MIDPOINT_TO_CORNER,
            CENTRE_TO_CORNER,
            CENTRE_TO_MIDPOINT,
            36.0 * math::DEG2RAD,
            90.0 * math::DEG2RAD,
            60.0 * math::DEG2RAD,
        ),
    };
    VariantConsts {
        ab,
        ac,
        bc,
        alpha,
        beta,
        gamma,
        cos_ab: math::cos(ab),
        cos_ac: math::cos(ac),
        sin_ac: math::sin(ac),
        cos_bc: math::cos(bc),
    }
}

/// Precomputed per-configuration geometry, built once by `build_geometry`.
///
/// This is the ISEA projection's `Geometry`, so the type itself is public, and a
/// caller obtains one from [`crate::Projection::build_geometry`]. Its fields are
/// `pub(crate)`, since the functions that read them are this crate's own. A
/// caller reaches the type itself through its re-export at [`crate::projections`],
/// never through this module, which is private.
///
/// Mirrors the eC's own cached state (the rotated icosahedron vertices, the face
/// centroid tables and the authalic Clenshaw coefficient pair) and additionally carries
/// `azimuth_deg` and `authalic` from the [`GridConfig`], because `forward` and
/// `inverse` receive only the geometry and so must read those two off it. The
/// `radial` and `consts` fields likewise thread the variant through, so that the
/// kernel needs no globals.
#[derive(Debug, Clone, PartialEq)]
pub struct IcoGeometry {
    /// `(geodetic -> authalic, authalic -> geodetic)` Clenshaw coefficients.
    pub(crate) authalic_cp: ([f64; 6], [f64; 6]),
    /// The twelve icosahedron vertices, in DGGAL's Y-up convention.
    pub(crate) ico_vertices: [[f64; 3]; 12],
    /// The twenty face centroids.
    pub(crate) face_centroids: [[f64; 3]; 20],
    /// Per face, the unit normals of the three planes through the origin and each of the
    /// face's edges, eC `icoFacePlanes` (`ri5x6.ec:111-119`), which decide the face a
    /// point is projected on.
    pub(crate) face_planes: [[[f64; 3]; 3]; 20],
    /// Per face, the three edge midpoints, normalised onto the sphere.
    pub(crate) ico3rd_mids: [[[f64; 3]; 3]; 20],
    /// Per face, the planar centre of its 5x6 triangle.
    pub(crate) ico56_center: [[f64; 2]; 20],
    /// Per face, the three planar edge midpoints of its 5x6 triangle.
    pub(crate) ico56_mids: [[[f64; 2]; 3]; 20],
    /// The vertex-2 azimuth offset, from the [`GridConfig`].
    pub(crate) azimuth_deg: f64,
    /// The orientation's longitude in degrees, from the [`GridConfig`]. Only
    /// `fix_poles` reads it, for the eC's `lon1 = -180 - orientation.lon`
    /// (`ri5x6.ec:466`).
    pub(crate) orientation_lon_deg: f64,
    /// Whether the geodetic to authalic conversion is applied.
    pub(crate) authalic: bool,
    /// The projection variant, the eC `radialVertex`.
    pub(crate) radial: RadialVertex,
    /// The variant's spherical-triangle constants.
    pub(crate) consts: VariantConsts,
}

// --------------------------------------------------------------------------- //
// Authalic latitude (Karney 2022 Clenshaw series, authalic.ec)
// --------------------------------------------------------------------------- //

/// Fold the twenty-one flat series coefficients into six Clenshaw coefficients for one
/// ellipsoid, through its third flattening `n`. Done once per geometry and direction
/// (`authalic.ec:40`).
pub(crate) fn precompute_coefficients(a: f64, b: f64, c: &[f64; 21]) -> [f64; 6] {
    let n = (a - b) / (a + b);
    let mut d = n;
    let mut cp = [0.0f64; 6];
    cp[0] = (((((c[5] * n + c[4]) * n + c[3]) * n + c[2]) * n + c[1]) * n + c[0]) * d;
    d *= n;
    cp[1] = ((((c[10] * n + c[9]) * n + c[8]) * n + c[7]) * n + c[6]) * d;
    d *= n;
    cp[2] = (((c[14] * n + c[13]) * n + c[12]) * n + c[11]) * d;
    d *= n;
    cp[3] = ((c[17] * n + c[16]) * n + c[15]) * d;
    d *= n;
    cp[4] = (c[19] * n + c[18]) * d;
    // faithful: DGGAL v0.0.6 as built with `-O2 -ffast-math` (`Makefile.dggal:133`),
    // BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, `precomputeCoefficients` at
    // `0x3d930`: the eC's last step, `d *= n` and then `cp[5] = C[20] * d`
    // (`authalic.ec:50-51`), is compiled as `(n * C[20]) * d` on the `d` of the step
    // before, `n` to the fifth (`0x3da2f`, `0x3da37`), which puts `cp[5]` one unit in the
    // last place from the source's in both directions of the conversion. The third
    // flattening itself stays a division (`0x3d93c`), and the other five coefficients are
    // compiled in the source's order.
    cp[5] = (n * c[20]) * d;
    cp
}

/// Evaluate the Clenshaw series at latitude `phi`, giving the converted latitude
/// (`authalic.ec:54`).
pub(crate) fn apply_coefficients(cp: &[f64; 6], phi: f64) -> f64 {
    let szeta = math::sin(phi);
    let czeta = math::cos(phi);
    // 2 * cos(2 * phi), kept in the eC's factored form. The library forms it as
    // `(2 * (s + c)) * (c - s)` (`0x3dad1` to `0x3dae1`), which is the same double: both
    // round the one exact product of `2`, `fl(c + s)` and `fl(c - s)`.
    let x = 2.0 * (czeta - szeta) * (czeta + szeta);
    let mut u0 = x * cp[5] + cp[4];
    let mut u1 = x * u0 + cp[3];
    // faithful: DGGAL v0.0.6 as built with `-O2 -ffast-math` (`Makefile.dggal:133`),
    // BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, `applyCoefficients` at
    // `0x3da90`, the one routine behind both `latGeodeticToAuthalic` and
    // `latAuthalicToGeodetic`: the eC's three recurrences `X * u1 - u0 + cp[2]`,
    // `X * u0 - u1 + cp[1]` and `X * u1 - u0 + cp[0]` (`authalic.ec:63-65`) are each
    // compiled as `X * u + (cp[k] - u')`, the difference of the coefficient and the older
    // term being formed first (`0x3daf6`, `0x3db10` and `0x3db24` to `0x3db31`). Measured
    // against the engine's own function over 200,000 latitudes drawn uniformly from
    // [-pi/2, pi/2], the source's order differs from the engine's at 170 in the geodetic
    // to authalic direction and at 162 in the other, by a unit in the last place; this
    // order agrees at all of them. The closing `phi + 2 * s * c * u0` is compiled as
    // `(2 * u0) * (s * c) + phi`, which rounds the same exact product and is kept as
    // written.
    u0 = x * u1 + (cp[2] - u0);
    u1 = x * u0 + (cp[1] - u1);
    u0 = x * u1 + (cp[0] - u0);
    phi + 2.0 * szeta * czeta * u0
}

// --------------------------------------------------------------------------- //
// Vector math (three components, DGGAL's Y-up convention)
// --------------------------------------------------------------------------- //

/// The dot product of two vectors, eC `Vector3D::DotProduct` (`src/ecere3D/Vector3D.ec:17-20`).
pub(crate) fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// The cross product of two vectors, eC `Vector3D::CrossProduct`
/// (`src/ecere3D/Vector3D.ec:22-27`).
pub(crate) fn cross(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// The unit vector along `v`, or the zero vector when `v` has no length, eC
/// `Vector3D::Normalize` (`src/ecere3D/Vector3D.ec:29-40`) as the engine was compiled.
///
/// Every normalisation of the projection goes through this one function in the library,
/// out of line, from the face planes of the constructor to the foot of the great circle
/// in `forwardVector`, so this port has one function too.
pub(crate) fn normalize(v: &[f64; 3]) -> [f64; 3] {
    // faithful: DGGAL v0.0.6 as built with `-O2 -ffast-math` (`Makefile.dggal:133`),
    // BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, `Vector3D::Normalize` at
    // `0x42460`, which replaces `Vector3D.ec:31-36`. The squares are summed as
    // `(z*z + y*y) + x*x` (`0x42472` to `0x4248a`), where the eC writes
    // `x*x + y*y + z*z`, and the three divisions by the length become one reciprocal and
    // three multiplications (`0x424a0` to `0x424ba`). The reciprocal alone moves 44 of
    // the 180 components of the face planes; the order of the sum moves no cached field,
    // but without it 528 of the census's 3,184 ISEA forward stations leave the engine's
    // answer.
    let m2 = (v[2] * v[2] + v[1] * v[1]) + v[0] * v[0];
    // faithful: the eC's test is `if(m)` on the length, a length different from zero
    // exactly (`Vector3D.ec:32`); py4dggs, and this port until the engine was read, tested
    // `length < 1e-15` instead, which the 108 stations nearest the poles of the RTEA
    // census reach. The library tests the sum of squares, which is zero exactly when the
    // length is, with `comisd` and a `je` that no parity test guards (`0x4248e`,
    // `0x42492`): under `-ffinite-math-only` an unordered comparison is taken for an
    // equal one, and a NaN is answered with the zero vector, as here.
    if m2 == 0.0 || m2.is_nan() {
        return [0.0, 0.0, 0.0];
    }
    let r = 1.0 / math::sqrt(m2);
    [v[0] * r, v[1] * r, v[2] * r]
}

/// The Euclidean length of `v`, eC `Vector3D`'s `length` property
/// (`src/ecere3D/Vector3D.ec:42`) as the engine was compiled.
pub(crate) fn vec_length(v: &[f64; 3]) -> f64 {
    // faithful: DGGAL v0.0.6 as built with `-O2 -ffast-math` (`Makefile.dggal:133`),
    // BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, the `length` getter at
    // `0x42430`, which replaces `Vector3D.ec:42`: `x*x + (z*z + y*y)` (`0x42439` to
    // `0x42452`), where the eC writes `x*x + y*y + z*z`.
    math::sqrt(v[0] * v[0] + (v[2] * v[2] + v[1] * v[1]))
}

// --------------------------------------------------------------------------- //
// Quaternion math, used to rotate the base icosahedron into the configured
// orientation (yaw is the orientation longitude, pitch the orientation latitude).
// --------------------------------------------------------------------------- //

/// The quaternion `[w, x, y, z]` of a yaw and a pitch, eC `Quaternion::YawPitch`
/// (`src/ecere3D/Quaternion.ec:10-21`), in the eC's own member order.
///
/// The eC calls it with the orientation negated, `q.YawPitch(-orientation.lon,
/// -orientation.lat)` (`ri5x6.ec:271`), and [`compute_vertices`] does the same. The
/// library's copy at `0x42670` multiplies each angle by `0.5` where the eC divides it by
/// 2, the same double, and forms the four products in the eC's pairs; nothing here needs
/// the compiled order.
pub(crate) fn quaternion_yaw_pitch(yaw: f64, pitch: f64) -> [f64; 4] {
    let s_yaw = math::sin(yaw / 2.0);
    let c_yaw = math::cos(yaw / 2.0);
    let s_pitch = math::sin(pitch / 2.0);
    let c_pitch = math::cos(pitch / 2.0);
    [
        c_pitch * c_yaw,
        s_pitch * c_yaw,
        c_pitch * s_yaw,
        s_pitch * s_yaw,
    ]
}

/// The vector `s` rotated by the quaternion `quat = [w, x, y, z]`, eC
/// `Vector3D::MultQuaternion` (`src/ecere3D/Vector3D.ec:44-58`), which is the eC's own
/// formula rather than the textbook rotation py4dggs writes: the same rotation, but not
/// the same doubles.
pub(crate) fn mult_quaternion(s: &[f64; 3], quat: &[f64; 4]) -> [f64; 3] {
    let [w, x, y, z] = *quat;
    // faithful: DGGAL v0.0.6 as built with `-O2 -ffast-math` (`Makefile.dggal:133`),
    // BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, `Vector3D::MultQuaternion` at
    // `0x422d0`, which replaces `a = w*w - (v.x*v.x+v.y*v.y+v.z*v.z)` at `Vector3D.ec:47`:
    // the squares of the vector part are subtracted in two steps, `(w*w - (y*y + x*x))`
    // and then `z*z` (`0x42325`, `0x42342` and `0x4236e`). The dot product, the cross
    // product and the three sums of the result are compiled in the eC's order; the
    // library forms `(2 * x) * dotVS` where the eC writes `2 * dotVS * x`, and the two
    // round the same exact product. Of the 36 components of the rotated vertices, the
    // textbook rotation left 12 off the engine's, by at most 1.1e-16, and the eC's formula
    // with the source's `a` 20.
    let a = (w * w - (y * y + x * x)) - z * z;
    let dot_vs = x * s[0] + y * s[1] + z * s[2];
    let cross = [
        s[1] * z - s[2] * y,
        s[2] * x - s[0] * z,
        s[0] * y - s[1] * x,
    ];
    [
        2.0 * dot_vs * x + a * s[0] + 2.0 * w * cross[0],
        2.0 * dot_vs * y + a * s[1] + 2.0 * w * cross[1],
        2.0 * dot_vs * z + a * s[2] + 2.0 * w * cross[2],
    ]
}

// --------------------------------------------------------------------------- //
// Coordinate conversions (DGGAL Y-up cartesian <-> geographic radians)
// --------------------------------------------------------------------------- //

/// Cartesian to geographic radians, eC `cartesianToGeo` (`ri5x6.ec:145`). Returns
/// `[lat, lon]`.
pub(crate) fn dggal_to_geo(x: f64, y: f64, z: f64) -> [f64; 2] {
    let p = math::sqrt(x * x + z * z);
    [math::atan2(-y, p), math::atan2(x, -z)]
}

/// Geographic radians to cartesian, eC `geoToCartesian` (`ri5x6.ec:137`).
pub(crate) fn geo_to_dggal_cart(lat: f64, lon: f64) -> [f64; 3] {
    let clat = math::cos(lat);
    [
        math::sin(lon) * clat,
        -math::sin(lat),
        -math::cos(lon) * clat,
    ]
}

/// Build the twelve icosahedron vertices, two poles and two staggered pentagons, then
/// rotate them into the configured orientation (eC `getVertices`, `ri5x6.ec:257`).
///
/// The eC writes the southern band's trigonometry as `-sin(-t)` and `cos(-t)`, which
/// are the same doubles as `sin(t)` and `cos(t)` used here.
///
/// `ta` keeps py4dggs's grouping, `(-270 + 72 * i) * DEG2RAD`, rather than the eC's
/// `Degrees { -180 - 36/2 - 72 } + s * i` (`ri5x6.ec:275`), and `ba` adds `step / 2` to
/// it as the eC does. The library, where `getVertices` is inlined into the constructor,
/// computes neither of those: it folds the first pentagon's two vertices to constants at
/// compile time, and for the other four forms `ta` as `i * s - fl(3 pi / 2)` and `ba` as
/// `i * s - fl(3 pi / 2 - s / 2)` (`0x3a5b3` to `0x3a5e2`, with the constants at
/// `0x4a358`, `0x4a360` and `0x4a368`), re-associating the half step into the constant.
/// The angles are the five values of `i` and no others, so the three forms were compared
/// on all of them: they give the same ten doubles, and the vertices are bit-identical to
/// the engine's, as measured on x86-64 Linux with glibc, where the engine and this crate
/// share one `libm`. The rewriting cannot change an answer, and this form stands.
pub(crate) fn compute_vertices(
    orientation_lat_rad: f64,
    orientation_lon_rad: f64,
) -> [[f64; 3]; 12] {
    let t = math::atan(0.5);
    let sin_t = math::sin(t);
    let cos_t = math::cos(t);
    let step = math::TWO_PI / 5.0;
    let mut verts = [[0.0f64; 3]; 12];

    verts[0] = [0.0, -1.0, 0.0];
    verts[11] = [0.0, 1.0, 0.0];

    for i in 0..5 {
        let ta = (-270.0 + 72.0 * i as f64) * math::DEG2RAD;
        let ba = ta + step / 2.0;
        verts[1 + i] = [math::cos(ta) * cos_t, -sin_t, math::sin(ta) * cos_t];
        verts[6 + i] = [math::cos(ba) * cos_t, sin_t, math::sin(ba) * cos_t];
    }

    let q = quaternion_yaw_pitch(-orientation_lon_rad, -orientation_lat_rad);
    for v in verts.iter_mut() {
        *v = mult_quaternion(v, &q);
    }

    verts
}

/// `sin(angle / 2)` between the unit vectors `a` and `b`, computed stably through the
/// normalised-midpoint cross product (eC `sqrtOneMinusDotOver2`,
/// `icoVertexGreatCircle.ec:340`), which avoids cancelling 1 against a dot product
/// close to 1.
pub(crate) fn sqrt_one_minus_dot_over2(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let mid = normalize(&[
        (a[0] + b[0]) / 2.0,
        (a[1] + b[1]) / 2.0,
        (a[2] + b[2]) / 2.0,
    ]);
    let c = cross(a, &mid);
    let d = vec_length(&c);
    if d < 1e-8 {
        return vec_length(&[a[0] - b[0], a[1] - b[1], a[2] - b[2]]) / 2.0;
    }
    d
}

/// The signed spherical excess of triangle ABC through the midpoint-cross form (eC
/// `sphericalTriArea`, `ri5x6.ec:224`), which stays accurate for the very small
/// sub-triangles this kernel works with.
// The chain is the eC's `Max(-1.0, Min(1.0, x))` written out; `clamp` would differ from
// it on a NaN input, which the chain, like the eC macro and py4dggs, resolves to 1.0.
#[allow(clippy::manual_clamp)]
pub(crate) fn spherical_tri_area(a: &[f64; 3], b: &[f64; 3], c: &[f64; 3]) -> f64 {
    let mid_ab = normalize(&[
        (a[0] + b[0]) / 2.0,
        (a[1] + b[1]) / 2.0,
        (a[2] + b[2]) / 2.0,
    ]);
    let mid_bc = normalize(&[
        (b[0] + c[0]) / 2.0,
        (b[1] + c[1]) / 2.0,
        (b[2] + c[2]) / 2.0,
    ]);
    let mid_ca = normalize(&[
        (c[0] + a[0]) / 2.0,
        (c[1] + a[1]) / 2.0,
        (c[2] + a[2]) / 2.0,
    ]);
    let cr = cross(&mid_bc, &mid_ca);
    math::asin(dot(&mid_ab, &cr).min(1.0).max(-1.0)) * 2.0
}

// --------------------------------------------------------------------------- //
// Planar barycentric helpers (eC `cartesianToBary` at `ri5x6.ec:211` and
// `baryToCartesian` at `ri5x6.ec:237`; `barycentric5x6.ec` is a separate Goldberg
// projection class that merely calls them, with a different `knownOneOverDet`)
// --------------------------------------------------------------------------- //

/// The barycentric coordinates of the planar point `p` in triangle `(p1, p2, p3)`.
/// A non-zero `known_one_over_det` short-circuits the determinant when the caller
/// already knows it, which is -6 or -1 here.
pub(crate) fn cartesian_to_bary(
    p: &[f64; 2],
    p1: &[f64; 2],
    p2: &[f64; 2],
    p3: &[f64; 2],
    known_one_over_det: f64,
) -> [f64; 3] {
    let d31x = p1[0] - p3[0];
    let d31y = p1[1] - p3[1];
    let d23x = p3[0] - p2[0];
    let d23y = p3[1] - p2[1];
    let d3px = p[0] - p3[0];
    let d3py = p[1] - p3[1];
    let o_det = if known_one_over_det != 0.0 {
        known_one_over_det
    } else {
        1.0 / (d23x * d31y - d23y * d31x)
    };
    let b0 = (d23x * d3py - d23y * d3px) * o_det;
    let b1 = (d31x * d3py - d31y * d3px) * o_det;
    // faithful: DGGAL v0.0.6 as built with `-O2 -ffast-math` (`Makefile.dggal:133`),
    // BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, `cartesianToBary` at `0x3c050`,
    // which replaces `b[2] = 1 - b[0] - b[1]` at `ri5x6.ec:221`: the two coordinates are
    // added first and their sum subtracted from 1 (`0x3c0d6` to `0x3c0de`). Both the
    // face's own barycentrics and those of `inverseVector` come through this one
    // function; without the port, 70 of the 654 census stations that lie on an ISEA face
    // leave the engine's answer.
    [b0, b1, 1.0 - (b0 + b1)]
}

/// The planar point of the barycentric coordinates `b` in triangle `(p1, p2, p3)`, eC
/// `baryToCartesian` (`ri5x6.ec:237-244`).
pub(crate) fn bary_to_cartesian(
    b: &[f64; 3],
    p1: &[f64; 2],
    p2: &[f64; 2],
    p3: &[f64; 2],
) -> [f64; 2] {
    [
        b[0] * p1[0] + b[1] * p2[0] + b[2] * p3[0],
        b[0] * p1[1] + b[1] * p2[1] + b[2] * p3[1],
    ]
}

/// Spherical interpolation by `movement` radians along the `distance`-long great-circle
/// arc running from `p0` towards `p1` (eC `slerpAngle`, `ri5x6.ec:176`).
pub(crate) fn slerp_angle(p0: &[f64; 3], p1: &[f64; 3], distance: f64, movement: f64) -> [f64; 3] {
    let s_distance = math::sin(distance);
    // faithful: the eC divides by `sin(distance)` unconditionally (`ri5x6.ec:180-181`),
    // and so does DGGAL v0.0.6 (BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`),
    // where this guarded function is compiled as written: `sin` at `0x3c16d` and the
    // `divsd` at `0x3c186`, with no test between. py4dggs returns `p0` when
    // `|sin(distance)| < 1e-15`, a branch the eC does not have.
    let o_o_sin = 1.0 / s_distance;
    let l0 = math::sin(distance - movement);
    let l1 = math::sin(movement);
    [
        (l0 * p0[0] + l1 * p1[0]) * o_o_sin,
        (l0 * p0[1] + l1 * p1[1]) * o_o_sin,
        (l0 * p0[2] + l1 * p1[2]) * o_o_sin,
    ]
}

// --------------------------------------------------------------------------- //
// Face and sub-triangle finding.
//
// The face is found as the vendored DGGAL v0.06 finds it: the faces in order from 0,
// the first whose three edge planes the point does not straddle (`ri5x6.ec:408-411`).
// py4dggs instead takes the face with the nearest centroid, with an early-accept
// threshold, from a newer DGGAL lineage. The two agree inside a face, but not on an edge
// that two faces share, where the eC's order picks the lower-numbered face and the
// nearest centroid can pick either. On an edge the 5x6 layout cuts, such as the edge
// from base cell 0 to 1 over the north pole, the point has two planar images, one on
// each face: the engine's own forward projection puts the north pole on face 0, at
// (0.5, 0), where the nearest centroid had put it on face 4, at (5, 4.5), and the two
// images quantise to different cells. On an edge the layout does not cut, the two faces
// compute the same planar point through different vertices, and so differ in its last
// bits, which on a cell boundary can decide the cell, as they did at (0, -168.8), on the
// edge from base cell 6 to 10.
//
// The sub-triangle within the face is found in the same manner and for the same reason:
// the eC's three thirds in order, the first that contains the point, and then the half
// within it (`icoVertexGreatCircle.ec:417-446`). py4dggs takes the nearest centroid here
// too, with two further early-accept thresholds absent from v0.06. Unlike a face edge, a
// sub-triangle boundary is not cut in the 5x6 layout, so the two choices give the same
// planar point to within rounding, and the point never moves to another place; but the
// shift is not nothing, and at a cell boundary it decides the cell. Measured over a
// 32,213-point sample rich in face medians, seams and poles, on the three projections:
// the choice differs at 1,176 of the 96,639 forward calls, moves the planar point at 950
// of them, and by at most 1.80e-12.
// --------------------------------------------------------------------------- //

/// eC `WITHIN_TRI_THRESHOLD` (`ri5x6.ec:330`): a point nearer a face's edge plane than
/// this is taken to lie on it, and that plane does not decide.
const WITHIN_TRI_THRESHOLD: f64 = 1e-12;

/// The angle between two unit vectors, in radians, eC `angleBetweenUnitVectors`
/// (`ri5x6.ec:122-135`), which takes the arc sine of half a chord, from whichever of `v`
/// and `-v` is nearer `u`.
// The bound is a min/max chain rather than `clamp`, which would differ from it on a NaN;
// see the second citation below.
#[allow(clippy::manual_clamp)]
fn angle_between_unit_vectors(u: &[f64; 3], v: &[f64; 3]) -> f64 {
    // faithful: DGGAL v0.0.6 as built with `-O2 -ffast-math` (`Makefile.dggal:133`),
    // BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, `angleBetweenUnitVectors` at
    // `0x3c350`, which replaces `u.x * v.x + u.y * v.y + u.z * v.z < 0` at `ri5x6.ec:125`:
    // the dot product is summed as `(u.x*v.x + u.z*v.z) + u.y*v.y` (`0x3c362` to
    // `0x3c392`). Its sign chooses the chord, and it can differ from the source's only for
    // vectors at a right angle to within rounding, where both chords give a right angle
    // and the face and sub-triangle tests that consume it reject the point on their
    // planes anyway; the order is ported as the library performs it all the same.
    let d = (u[0] * v[0] + u[2] * v[2]) + u[1] * v[1];
    // faithful: the eC's `s < -1 ? -1 : s > 1 ? 1 : s` (`ri5x6.ec:128`, `:133`) is
    // compiled as `minsd` with 1 and then `maxsd` with -1 (`0x3c3c2`, `0x3c3ca`), which
    // answers a NaN with 1, as `f64::min` and `f64::max` do; the conditional chain would
    // let it through.
    if d < 0.0 {
        let s = vec_length(&[-v[0] - u[0], -v[1] - u[1], -v[2] - u[2]]) / 2.0;
        math::PI - 2.0 * math::asin(s.min(1.0).max(-1.0))
    } else {
        let s = vec_length(&[v[0] - u[0], v[1] - u[1], v[2] - u[2]]) / 2.0;
        2.0 * math::asin(s.min(1.0).max(-1.0))
    }
}

/// Whether the unit vector `v` lies within the spherical triangle whose edge planes are
/// `planes` and whose first vertex is `v1`, eC `vertexWithinSphericalTriPlanes`
/// (`ri5x6.ec:362-383`): no more than a right angle from `v1`, which rules out the
/// antipodal triangle, and on the same side of every plane it is not within
/// [`WITHIN_TRI_THRESHOLD`] of. A point on an edge therefore lies within both faces that
/// share it.
fn vertex_within_spherical_tri_planes(v: &[f64; 3], planes: &[[f64; 3]; 3], v1: &[f64; 3]) -> bool {
    // faithful: the eC declares the angle's function as returning `Degrees` and casts its
    // result to `Radians` here (`ri5x6.ec:122`, `:368`), and this port used to send the
    // angle through degrees and back on that reading. The engine performs no conversion
    // at all: in DGGAL v0.0.6 as built with `-O2 -ffast-math` (`Makefile.dggal:133`),
    // BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, the function at `0x3c350`
    // returns the arc sine doubled, and both callers compare that value with
    // `fl(PI / 2)` from `0x4a048` directly (`0x3c74d` and `0x3c876`). The eC compiler
    // emits no conversion for the returned unit, which is the reading the machine code
    // bears out, so the angle is compared as it comes.
    if angle_between_unit_vectors(v, v1) > math::PI / 2.0 {
        return false;
    }
    let mut sgn = 0i32;
    for plane in planes {
        let d = plane[0] * v[0] + plane[1] * v[1] + plane[2] * v[2];
        if d.abs() > WITHIN_TRI_THRESHOLD {
            // `Sgn` of a value that is not zero.
            let s = if d > 0.0 { 1 } else { -1 };
            if sgn != 0 && s != sgn {
                return false;
            }
            sgn = s;
        }
    }
    true
}

/// Whether the vector `v` lies within the spherical triangle `(v1, v2, v3)`, eC
/// `vertexWithinSphericalTri` (`ri5x6.ec:333-360`): the same test as
/// [`vertex_within_spherical_tri_planes`], with the three edge planes built on the spot
/// from the triangle's own vertices instead of taken precomputed. The eC's
/// `Plane::FromPoints({ 0, 0, 0 }, a, b)` gives `normalize(b x a)`
/// (`src/ecere3D/Plane.ec:14-24`), as in [`build_geometry`]'s `face_planes`.
fn vertex_within_spherical_tri(v: &[f64; 3], v1: &[f64; 3], v2: &[f64; 3], v3: &[f64; 3]) -> bool {
    let planes = [
        normalize(&cross(v2, v1)),
        normalize(&cross(v3, v2)),
        normalize(&cross(v1, v3)),
    ];
    vertex_within_spherical_tri_planes(v, &planes, v1)
}

/// The `(top, bottom)` face pair of a rhombus of the 5x6 layout, or `None` when the
/// rhombus index names no rhombus (eC `getFace`'s switch, `ri5x6.ec:309-322`).
pub(crate) fn rhombus_faces(rhombus: i64) -> Option<(i32, i32)> {
    match rhombus {
        0 => Some((0, 5)),
        2 => Some((1, 6)),
        4 => Some((2, 7)),
        6 => Some((3, 8)),
        8 => Some((4, 9)),
        1 => Some((10, 15)),
        3 => Some((11, 16)),
        5 => Some((12, 17)),
        7 => Some((13, 18)),
        9 => Some((14, 19)),
        _ => None,
    }
}

/// The icosahedron face the unit vector `p` is projected on: the first face, from 0,
/// within whose edge planes it lies (eC `forward`, `ri5x6.ec:408-411`), or `None` if
/// there is none, which for a unit vector does not happen, the faces covering the sphere.
// faithful: the eC's loop order, which on an edge shared by two faces, and so on the
// edges the 5x6 layout cuts, decides which of the point's two planar images it gets. See
// the section note above.
pub(crate) fn find_face(geom: &IcoGeometry, p: &[f64; 3]) -> Option<usize> {
    (0..20).find(|&face| {
        vertex_within_spherical_tri_planes(
            p,
            &geom.face_planes[face],
            &geom.ico_vertices[FACE_VERTICES[face][0]],
        )
    })
}

/// Locate the sub-triangle 0 to 5 of `face` containing the direction `v` (eC
/// `forwardIcoFace`, `icoVertexGreatCircle.ec:417-446`): the three thirds about the face
/// centre in turn, the first that contains `v`, and then the half of that third. See the
/// section note above on what the choice decides.
// faithful: the eC runs both tests on the `vCenter` and `vMid` it has just built and has
// not yet normalised (`icoVertexGreatCircle.ec:405-409`, `:420`, `:430`, `:440`; the two
// normalisations come afterwards, at `:448-449`, for the frame alone). The selection
// therefore works on a centre of length about 0.795 and a midpoint of about 0.951, and it
// is ported on those same quantities. The frame downstream takes the normalised centre
// and midpoints that [`build_geometry`] caches, which are the same doubles.
pub(crate) fn find_sub_tri(geom: &IcoGeometry, face: usize, v: &[f64; 3]) -> usize {
    let [ii0, ii1, ii2] = FACE_VERTICES[face];
    let v1 = geom.ico_vertices[ii0];
    let v2 = geom.ico_vertices[ii1];
    let v3 = geom.ico_vertices[ii2];
    let v_center = [
        (v1[0] + v2[0] + v3[0]) / 3.0,
        (v1[1] + v2[1] + v3[1]) / 3.0,
        (v1[2] + v2[2] + v3[2]) / 3.0,
    ];
    let mid = |a: &[f64; 3], b: &[f64; 3]| {
        [
            (a[0] + b[0]) / 2.0,
            (a[1] + b[1]) / 2.0,
            (a[2] + b[2]) / 2.0,
        ]
    };

    if vertex_within_spherical_tri(v, &v_center, &v2, &v3) {
        let v_mid = mid(&v2, &v3);
        if vertex_within_spherical_tri(v, &v_center, &v_mid, &v3) {
            0
        } else {
            1
        }
    } else if vertex_within_spherical_tri(v, &v_center, &v3, &v1) {
        let v_mid = mid(&v3, &v1);
        // faithful: the eC's second test in this branch closes on `v3`, not on `v1`
        // (`icoVertexGreatCircle.ec:432`), unlike the other two, which close on the
        // vertex of the third they did not reach by the midpoint.
        if vertex_within_spherical_tri(v, &v_center, &v_mid, &v3) {
            2
        } else {
            3
        }
    } else if vertex_within_spherical_tri(v, &v_center, &v1, &v2) {
        let v_mid = mid(&v1, &v2);
        if vertex_within_spherical_tri(v, &v_center, &v_mid, &v2) {
            4
        } else {
            5
        }
    } else {
        // The eC has no fourth branch: `subTri` keeps the 0 it was declared with while
        // `vMid` stays uninitialised and the frame's second vertex a null pointer
        // (`icoVertexGreatCircle.ec:410-414`), so what it would compute here is not
        // defined. No point can reach it: the three thirds cover the face, their outer
        // sides are the face's own edges, and `find_face` admits only a point within
        // those; a point on a shared side lies within both thirds, the test allowing a
        // plane distance of up to [`WITHIN_TRI_THRESHOLD`]. Measured with this branch
        // made to panic: never entered by the 32,213-point forward sample on the three
        // projections, nor by any case of the six oracle suites or the library tests.
        0
    }
}

// --------------------------------------------------------------------------- //
// The equal-area vector projection itself (icoVertexGreatCircle.ec)
// --------------------------------------------------------------------------- //

/// The forward Snyder map of the unit vector `v` in the spherical triangle `(a, b, c)`
/// to planar barycentrics in the triangle `(pa, pb, pc)` (eC `forwardVector`,
/// `icoVertexGreatCircle.ec:357`). It is equal-area: the planar area fraction equals the
/// spherical one, `area_abp / SDT_AREA`.
pub(crate) fn forward_vector(
    v: &[f64; 3],
    a: &[f64; 3],
    b: &[f64; 3],
    c: &[f64; 3],
    pa: &[f64; 2],
    pb: &[f64; 2],
    pc: &[f64; 2],
) -> [f64; 2] {
    let c1 = cross(a, v);
    let c2 = cross(b, c);
    // faithful: the eC goes straight from the normalisation to `sphericalTriArea`
    // (`icoVertexGreatCircle.ec:383-385`), and so does the library, where this function is
    // inlined into `forwardIcoFace` and the call to `Normalize` at `0x3f672` is followed
    // by the call to `sphericalTriArea` at `0x3f680` with nothing between. py4dggs turns
    // `p` round when `p . B < 0`, and this port did too, on the argument that the test is
    // inert for a point inside its own sub-triangle, where it is. It is not inert for a
    // point that coincides with A to within rounding: `A x v` is then rounding residue,
    // `p` a direction of no meaning, and its sign whatever the residue makes it. The
    // census of 3,184 stations meets one such point, (0, -168.8) on RTEA, the midpoint of
    // the edge from base cell 6 to 10 and so the vertex A of its sub-triangle; there the
    // turn gave (4.5, 5.499999999999998) where the engine answers (4.5, 5.5) exactly, and
    // without it this port does too.
    let p = normalize(&cross(&c1, &c2));

    let area_abp = spherical_tri_area(a, b, &p).max(0.0);
    let h = sqrt_one_minus_dot_over2(a, v) / sqrt_one_minus_dot_over2(a, &p);

    let b0 = 1.0 - h;
    let b2 = h.min(h * area_abp / SDT_AREA);
    let b1 = h - b2;

    bary_to_cartesian(&[b0, b1, b2], pa, pb, pc)
}

/// The inverse Snyder map: the planar point `pi`, in triangle `(pa, pb, pc)`, back to a
/// unit vector in the spherical triangle `(a, b, c)` (eC `inverseVector`,
/// `icoVertexGreatCircle.ec:174`). The closed form solves for the foot `p` on edge BC
/// and then slerps from A; a fallback of two slerps handles the near-edge and
/// near-singular cases. `consts` supplies the per-variant sides and angles.
#[allow(clippy::too_many_arguments)]
// See the note on `spherical_tri_area` for why the fallback's bound is a min/max chain
// rather than `clamp`.
#[allow(clippy::manual_clamp)]
pub(crate) fn inverse_vector(
    pi: &[f64; 2],
    pa: &[f64; 2],
    pb: &[f64; 2],
    pc: &[f64; 2],
    a: &[f64; 3],
    b: &[f64; 3],
    c: &[f64; 3],
    b_is_a: bool,
    consts: &VariantConsts,
) -> [f64; 3] {
    let bary = cartesian_to_bary(pi, pa, pb, pc, -6.0);

    if bary[0] > 1.0 - 1e-15 {
        return *a;
    }
    if bary[1] > 1.0 - 1e-15 {
        return *b;
    }
    if bary[2] > 1.0 - 1e-15 {
        return *c;
    }

    let c1 = cross(b, c);

    let h = 1.0 - bary[0];
    let b2oh = bary[2] / h;
    let b2oh_abc = b2oh * SDT_AREA;
    // faithful: eC icoVertexGreatCircle.ec:198 takes sin(b2ohABC) directly; py4dggs
    // substitutes the analytically equal double-angle form 2 * halfC * sqrt(1 - halfC^2).
    let s = math::sin(b2oh_abc);
    let c01 = if b_is_a { consts.cos_ab } else { consts.cos_bc };
    let c12 = consts.cos_ac;
    let c20 = if b_is_a { consts.cos_bc } else { consts.cos_ab };
    let s12 = consts.sin_ac;
    // faithful: eC icoVertexGreatCircle.ec:203 computes the scalar triple product of
    // A, B and C at run time; py4dggs substitutes the analytically equal closed form
    // sqrt((5 - 2 * sqrt(5)) / 15), which this port keeps as `parallelepiped_v` and
    // pins against every frame in the tests below.
    let v = a[0] * c1[0] + a[1] * c1[1] + a[2] * c1[2];
    let half_c = math::sin(b2oh_abc / 2.0);
    // Half-angle identity; 1 - sqrt(1 - S * S) would be far less precise.
    let cc = 2.0 * half_c * half_c;
    let f = s * v + cc * (c01 * c12 - c20);
    let g = cc * s12 * (1.0 + c01);
    let f2 = f * f;
    let g2 = g * g;
    let gf = g * f;
    let numerator = s12 * (f2 - g2) - 2.0 * gf * c12;
    let divisor = s12 * (f2 + g2);

    if numerator.abs() > 1e-9 && divisor.abs() > 1e-9 {
        // The trigonometry-free branch, equivalent to two slerps, for the
        // non-degenerate cases.
        let o_o_divisor = 1.0 / divisor;
        let ap = (numerator * o_o_divisor).max(0.0);
        let bp = (2.0 * gf * o_o_divisor).min(1.0);
        let px = ap * b[0] + bp * c[0];
        let py = ap * b[1] + bp * c[1];
        let pz = ap * b[2] + bp * c[2];

        let av = a[0] * px + a[1] * py + a[2] * pz;
        let bv = 1.0 + h * h * (av - 1.0);
        let bvp = h * math::sqrt((1.0 + bv) / (1.0 + av));
        let avp = bv - av * bvp;

        [
            avp * a[0] + bvp * px,
            avp * a[1] + bvp * py,
            avp * a[2] + bvp * pz,
        ]
    } else {
        // The two-slerp fallback for the degenerate cases where the linear-algebra
        // form breaks down (eC `icoVertexGreatCircle.ec:238-271`).
        //
        // faithful: eC :246 reads `beta` and `gamma` as the triangle's angles, not its
        // sides: 60 and 36 degrees for ISEA, with `alpha` 90 degrees for ISEA and IVEA
        // but 36 for RTEA; the comment beside `areaABC = Degrees { 6 }` at :243, which
        // reads `beta + gamma + alpha - Pi`, settles it, since 60 + 36 + 90 - 180 = 6.
        // py4dggs instead passes the sides AB and AC as beta and gamma and hard-codes
        // alpha as pi/2, which does not satisfy that identity; this port follows the eC.
        let beta = consts.beta;
        let gamma = consts.gamma;
        let alpha = consts.alpha;
        let b1pb2 = bary[1] + bary[2];
        let up_over_up_pvp = if b1pb2 < 1e-11 {
            0.0
        } else if b_is_a {
            bary[1] / b1pb2
        } else {
            bary[2] / b1pb2
        };
        let rho_plus_delta = beta + gamma - up_over_up_pvp * SDT_AREA;
        // T - U = rho + delta + alpha - Pi
        let area_abd = rho_plus_delta + alpha - math::PI;

        let (d, bd) = if area_abd.abs() < 1e-11 {
            // The B or C vertex at the angle alpha, a right angle for ISEA and IVEA.
            (if b_is_a { *b } else { *c }, consts.ab)
        } else if (area_abd - SDT_AREA).abs() < 1e-13 {
            // The B or C vertex at the angle gamma.
            (if b_is_a { *c } else { *b }, consts.bc)
        } else {
            let ad = 2.0 * math::atan2(g, f);
            let d = slerp_angle(b, c, consts.ac, ad);
            // The eC's `Max(-1.0, Min(1.0, x))`; see the note on `spherical_tri_area`
            // for why this is not `clamp`.
            let bd = math::acos(dot(a, &d).min(1.0).max(-1.0));
            (d, bd)
        };

        // A is the vertex the great circles radiate from, at the angle beta; the
        // identity is x' / (x' + y') = 1 - b_0.
        // faithful: the eC calls `asin` on the product unclamped
        // (`icoVertexGreatCircle.ec:270`), and so does DGGAL v0.0.6 (BuildID
        // `e75d6ab18f8b460713ebcd21df6f07922fb909c3`): `sin` at `0x3ef5a`, the product,
        // and `asin` at `0x3ef64`. py4dggs bounds the argument by 1 first, a step the eC
        // does not take.
        let x = 2.0 * math::asin((1.0 - bary[0]) * math::sin(bd / 2.0));
        slerp_angle(a, &d, bd, x)
    }
}

/// Resolve the sub-triangle `sub_tri` of `face` to its projection frame: the spherical
/// vertex triple and the 5x6 vertex triple, both reordered into the eC's `(a, b, c)`
/// roles for this variant, plus the `b_is_a` flag.
///
/// [`forward_ico_face`] and [`inverse_ico_face`] differ only in how they find
/// `sub_tri`, by the spherical containment tests or by planar barycentrics; everything
/// downstream is this shared frame selection. It is pure construction and integer index
/// selection, with no floating-point arithmetic, so the output stays bit-identical. The
/// eC assigns `a = vb`, `b = bIsA ? va : vc` and `c = bIsA ? vc : va`, with
/// `bIsA = (radialVertex == ivea) ^ (subTri == 0 || subTri == 3 || subTri == 4)`
/// (`icoVertexGreatCircle.ec:334` and `:451`).
pub(crate) fn sub_tri_frame(
    geom: &IcoGeometry,
    face: usize,
    sub_tri: usize,
) -> ([[f64; 3]; 3], [[f64; 2]; 3], bool) {
    let [ii0, ii1, ii2] = FACE_VERTICES[face];
    let v1 = geom.ico_vertices[ii0];
    let v2 = geom.ico_vertices[ii1];
    let v3 = geom.ico_vertices[ii2];
    let p1 = FACE_5X6[face][0];
    let p2 = FACE_5X6[face][1];
    let p3 = FACE_5X6[face][2];
    let tri3rd = sub_tri >> 1;

    let p5x6 = [
        geom.ico56_mids[face][tri3rd],
        match sub_tri {
            0 | 2 => p3,
            1 | 4 => p2,
            _ => p1,
        },
        geom.ico56_center[face],
    ];
    let v3d = [
        geom.ico3rd_mids[face][tri3rd],
        match sub_tri {
            0 | 2 => v3,
            1 | 4 => v2,
            _ => v1,
        },
        geom.face_centroids[face],
    ];

    let b_is_a = (geom.radial == RadialVertex::Ivea) ^ matches!(sub_tri, 0 | 3 | 4);
    let (va, vb, vc) = geom.radial.vertex_order();
    let a = vb;
    let b_idx = if b_is_a { va } else { vc };
    let c_idx = if b_is_a { vc } else { va };
    (
        [v3d[a], v3d[b_idx], v3d[c_idx]],
        [p5x6[a], p5x6[b_idx], p5x6[c_idx]],
        b_is_a,
    )
}

/// Project the direction `v` onto `face`: find its sub-triangle, then map within the
/// centre-anchored sub-triangle to 5x6 planar coordinates (eC `forwardIcoFace`,
/// `icoVertexGreatCircle.ec:396`).
pub(crate) fn forward_ico_face(geom: &IcoGeometry, face: usize, v: &[f64; 3]) -> [f64; 2] {
    let sub_tri = find_sub_tri(geom, face, v);
    let (v3d, p5x6, _) = sub_tri_frame(geom, face, sub_tri);
    forward_vector(v, &v3d[0], &v3d[1], &v3d[2], &p5x6[0], &p5x6[1], &p5x6[2])
}

/// The inverse of [`forward_ico_face`]: from a 5x6 planar point on `face`, pick the
/// sub-triangle by planar barycentrics, then run the inverse vector map (eC
/// `inverseIcoFace`, `icoVertexGreatCircle.ec:277`).
pub(crate) fn inverse_ico_face(geom: &IcoGeometry, face: usize, v56: &[f64; 2]) -> [f64; 3] {
    let p1 = FACE_5X6[face][0];
    let p2 = FACE_5X6[face][1];
    let p3 = FACE_5X6[face][2];

    let b = cartesian_to_bary(v56, &p1, &p2, &p3, -1.0);
    let sub_tri = if b[0] <= b[1] && b[0] <= b[2] {
        if b[1] < b[2] { 0 } else { 1 }
    } else if b[1] <= b[0] && b[1] <= b[2] {
        if b[0] < b[2] { 2 } else { 3 }
    } else if b[0] < b[1] {
        4
    } else {
        5
    };

    let (v3d, p5x6, b_is_a) = sub_tri_frame(geom, face, sub_tri);
    inverse_vector(
        v56,
        &p5x6[0],
        &p5x6[1],
        &p5x6[2],
        &v3d[0],
        &v3d[1],
        &v3d[2],
        b_is_a,
        &geom.consts,
    )
}

/// Compute the icosahedron and the authalic data for `config` and `radial`.
///
/// The authalic coefficients are fixed to WGS84; the `authalic` flag only decides
/// whether the conversion is applied, in [`forward`] and [`inverse`]. `azimuth_deg` and
/// `authalic` are carried on the geometry so that the pure forward and inverse need no
/// configuration argument (eC `RI5x6Projection` constructor, `ri5x6.ec:97`).
pub(crate) fn build_geometry(config: &GridConfig, radial: RadialVertex) -> IcoGeometry {
    let authalic_cp = (
        precompute_coefficients(WGS84_A, WGS84_B, &CXIPHI),
        precompute_coefficients(WGS84_A, WGS84_B, &CPHIXI),
    );

    let ico_vertices = compute_vertices(
        config.orientation_lat_deg * math::DEG2RAD,
        config.orientation_lon_deg * math::DEG2RAD,
    );

    // faithful: the eC builds `vCenter` by dividing the vertex sum by three and only then
    // normalising (`icoVertexGreatCircle.ec:288-292` and `:405-409`, then `:330` and
    // `:448`), in both `forwardIcoFace` and `inverseIcoFace`, at every call; this caches
    // the same quantity once. py4dggs normalises the bare sum instead, and this port did
    // too, on the ground that the division is a scaling the normalisation undoes; it is,
    // analytically, but not to the last bit, and without the division 1,951 of the
    // census's 3,184 ISEA forward stations leave the engine's answer. `ico3rd_mids`
    // below halves the pair sum before normalising,
    // as the eC does at `:420` and `:449`, and [`find_sub_tri`] builds the same divided
    // sum for its selection, which runs before either normalisation.
    let mut face_centroids = [[0.0f64; 3]; 20];
    for f in 0..20 {
        let [i0, i1, i2] = FACE_VERTICES[f];
        let v0 = ico_vertices[i0];
        let v1 = ico_vertices[i1];
        let v2 = ico_vertices[i2];
        face_centroids[f] = normalize(&[
            (v0[0] + v1[0] + v2[0]) / 3.0,
            (v0[1] + v1[1] + v2[1]) / 3.0,
            (v0[2] + v1[2] + v2[2]) / 3.0,
        ]);
    }

    // The eC's `icoFacePlanes` (`ri5x6.ec:111-119`): for each face, the planes through the
    // origin and the edges v1 v2, v2 v3 and v3 v1, each built by
    // `Plane::FromPoints({ 0, 0, 0 }, a, b)`, whose normal is the normalised cross product
    // of `b - 0` with `a - 0` (`src/ecere3D/Plane.ec:14-24`). Subtracting the origin
    // changes no bit, so the normal is `normalize(b x a)`.
    let mut face_planes = [[[0.0f64; 3]; 3]; 20];
    for (f, planes) in face_planes.iter_mut().enumerate() {
        let [i0, i1, i2] = FACE_VERTICES[f];
        let (v1, v2, v3) = (ico_vertices[i0], ico_vertices[i1], ico_vertices[i2]);
        *planes = [
            normalize(&cross(&v2, &v1)),
            normalize(&cross(&v3, &v2)),
            normalize(&cross(&v1, &v3)),
        ];
    }

    let mut ico3rd_mids = [[[0.0f64; 3]; 3]; 20];
    let mut ico56_center = [[0.0f64; 2]; 20];
    let mut ico56_mids = [[[0.0f64; 2]; 3]; 20];

    for f in 0..20 {
        let [ii0, ii1, ii2] = FACE_VERTICES[f];
        let v1 = ico_vertices[ii0];
        let v2 = ico_vertices[ii1];
        let v3 = ico_vertices[ii2];
        let p1 = FACE_5X6[f][0];
        let p2 = FACE_5X6[f][1];
        let p3 = FACE_5X6[f][2];

        ico56_center[f] = [(p1[0] + p2[0] + p3[0]) / 3.0, (p1[1] + p2[1] + p3[1]) / 3.0];

        ico3rd_mids[f] = [
            normalize(&[
                (v2[0] + v3[0]) / 2.0,
                (v2[1] + v3[1]) / 2.0,
                (v2[2] + v3[2]) / 2.0,
            ]),
            normalize(&[
                (v3[0] + v1[0]) / 2.0,
                (v3[1] + v1[1]) / 2.0,
                (v3[2] + v1[2]) / 2.0,
            ]),
            normalize(&[
                (v1[0] + v2[0]) / 2.0,
                (v1[1] + v2[1]) / 2.0,
                (v1[2] + v2[2]) / 2.0,
            ]),
        ];
        ico56_mids[f] = [
            [(p2[0] + p3[0]) / 2.0, (p2[1] + p3[1]) / 2.0],
            [(p3[0] + p1[0]) / 2.0, (p3[1] + p1[1]) / 2.0],
            [(p1[0] + p2[0]) / 2.0, (p1[1] + p2[1]) / 2.0],
        ];
    }

    IcoGeometry {
        authalic_cp,
        ico_vertices,
        face_centroids,
        face_planes,
        ico3rd_mids,
        ico56_center,
        ico56_mids,
        azimuth_deg: config.azimuth_deg,
        orientation_lon_deg: config.orientation_lon_deg,
        authalic: config.authalic,
        radial,
        consts: variant_consts(radial),
    }
}

/// Geodetic latitude to authalic latitude (eC `latGeodeticToAuthalic`,
/// `authalic.ec:9`). The Karney coefficients live on the geometry.
pub(crate) fn geodetic_to_authalic(geom: &IcoGeometry, lat_rad: f64) -> f64 {
    apply_coefficients(&geom.authalic_cp.0, lat_rad)
}

/// The latitude on the sphere the projection works on, for a geodetic latitude in radians:
/// [`geodetic_to_authalic`] where the geometry converts latitudes, as [`forward`] does (the
/// eC's own `latGeodeticToAuthalic` of the projection, `ri5x6.ec:252-255`, called at
/// `ri5x6.ec:397`), and the latitude itself where it takes them as latitudes on the sphere.
pub(crate) fn sphere_latitude(geom: &IcoGeometry, lat_rad: f64) -> f64 {
    if geom.authalic {
        geodetic_to_authalic(geom, lat_rad)
    } else {
        lat_rad
    }
}

/// Authalic latitude to geodetic latitude (eC `latAuthalicToGeodetic`,
/// `authalic.ec:14`). The Karney coefficients live on the geometry.
pub(crate) fn authalic_to_geodetic(geom: &IcoGeometry, lat_rad: f64) -> f64 {
    apply_coefficients(&geom.authalic_cp.1, lat_rad)
}

/// Derive the face identifier from a 5x6 planar point alone (eC `getFace`,
/// `ri5x6.ec:292`). Edge and corner points are nudged inwards by `epsilon` so that the
/// rhombus lookup is unambiguous; a point outside the 5x6 layout gives -1.
///
/// The engine calls `getFace` only from `inverse`, after the wrap, and so does this port,
/// through [`face_after_wrap`]; this is the same function on a point no wrap has touched,
/// for the tests.
#[cfg(test)]
pub(crate) fn get_face(x: f64, y: f64) -> i32 {
    face_after_wrap(x, y, Wrap::None, x)
}

/// Which of the edge wraps of the eC's `inverse` (`ri5x6.ec:516-531`) a planar point
/// went through, as far as it matters to `getFace`, which the library inlines into
/// `inverse` and compiles together with the wraps: see [`face_after_wrap`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Wrap {
    /// No wrap, or one whose abscissa the library nudges as the eC writes it: the first
    /// branch, `x < 0 && y < 0`, and the fourth, fifth and eighth.
    None,
    /// The second and third branches, which add 5 to an abscissa below 0 when the
    /// ordinate is below 1, or at most `1e-10` above it.
    AddFive,
    /// The sixth branch, which takes 5 from both coordinates above 5.
    SubtractFiveFromBoth,
    /// The seventh branch, which takes 5 from an abscissa above 5 and sets the ordinate
    /// to 0.
    SubtractFiveOrdinateZero,
}

/// eC `getFace` (`ri5x6.ec:292-328`) on the point `(x, y)` that `inverse`'s wrap has
/// produced from an abscissa `x_before`, as the library performs the two together.
///
/// The nudge of the abscissa follows the eC's text except at the sites below, and the
/// nudge of the ordinate and the rhombus lookup follow it everywhere; a transcription of
/// the library's code from `0x3cd4c` to `0x3d557` was compared with the eC's text over
/// 921,899 planar points built on every threshold of both, and the two choose different
/// faces at nine of them, all at the two thresholds that move. The engine's own inverse
/// agrees with this function at all 921,899 on the three projections: at the 280,945 that
/// land on a face, to the last bit of the geographic point, and at the other 640,954 in
/// rejecting the point. That was measured on x86-64 Linux with glibc, where the engine and
/// this crate share one `libm`, and says nothing of another target.
fn face_after_wrap(x: f64, y: f64, wrap: Wrap, x_before: f64) -> i32 {
    let epsilon = 1e-11;
    let mut x = x;
    let mut y = y;
    // faithful: DGGAL v0.0.6 as built with `-O2 -ffast-math` (`Makefile.dggal:133`),
    // BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, `inverse` at `0x3cd20`, with
    // `getFace` inlined: on four of the wrap branches gcc folded the wrap's `+ 5` or
    // `- 5` into the nudge's `+ 1e-11` or `- 1e-11`, or into the test that guards it,
    // re-associating the two into one constant, `fl(5 + 1e-11)` at `0x4a188` or
    // `fl(5 + fl(5 - 1e-11))` at `0x4a3b8`, where the eC writes them apart
    // (`ri5x6.ec:297-298` on the point of `:516-531`).
    x = match wrap {
        // `(x + 5) + 1e-11` as `x + fl(5 + 1e-11)` at `0x3d37d` and `0x3ce7f`, on the
        // abscissa below -5 that the nudge's `x < 0` becomes (`0x3d368`, `0x3ce6a`, an
        // exact rewriting). The first branch keeps the eC's order (`0x3ced5`). Such a
        // point is outside the layout, since its ordinate is 5 or more, and gets no face
        // either way.
        Wrap::AddFive if x_before < -5.0 => x_before + (5.0 + epsilon),
        // `(x - 5) < fl(5 - 1e-11)` as `x < fl(5 + fl(5 - 1e-11))` at `0x3cd97`, the one
        // threshold that moves: at the abscissa `fl(5 + fl(5 - 1e-11))` itself, the eC's
        // test holds and the library's does not. And `(x - 5) - 1e-11` as
        // `x - fl(5 + 1e-11)` above 10 at `0x3d50e`. The two give the same double up to
        // 13.000000000010001, where they first differ, and differ at about three doubles
        // in five of [13, 16) and at every double drawn from [16, 20); at every one of
        // them the nudged abscissa exceeds 5 on either reading, the point leaves the
        // layout and gets no face, so no face depends on which is taken.
        Wrap::SubtractFiveFromBoth => {
            if y > x && x_before < 5.0 + (5.0 - epsilon) {
                x + epsilon
            } else if x_before > 10.0 {
                x_before - (5.0 + epsilon)
            } else if y < x && x > epsilon {
                x - epsilon
            } else {
                x
            }
        }
        // `(x - 5) > 1e-11` as `x > fl(5 + 1e-11)` at `0x3d528`, and `(x - 5) - 1e-11` as
        // `x - fl(5 + 1e-11)` at `0x3d52e` and, above 10, at `0x3d4de`. The threshold
        // moves by one double: at the abscissa `fl(5 + 1e-11)` the eC nudges and the
        // library does not, and the two name different faces there.
        Wrap::SubtractFiveOrdinateZero => {
            if x_before > 5.0 + epsilon {
                x_before - (5.0 + epsilon)
            } else {
                x
            }
        }
        Wrap::None | Wrap::AddFive => {
            if x < 0.0 || (y > x && x < 5.0 - epsilon) {
                x + epsilon
            } else if x > 5.0 || (y < x && x > 0.0 + epsilon) {
                x - epsilon
            } else {
                x
            }
        }
    };
    if y < 0.0 || (x > y && y < 6.0 - epsilon) {
        y += epsilon;
    } else if y > 6.0 || (x < y && y > 0.0 + epsilon) {
        y -= epsilon;
    }

    if (0.0..=5.0).contains(&x) && (0.0..=6.0).contains(&y) {
        let ix = (x.floor() as i64).clamp(0, 4);
        let iy = (y.floor() as i64).clamp(0, 5);
        if iy == ix || iy == ix + 1 {
            let rhombus = ix + iy;
            let top = x - ix as f64 > y - iy as f64;
            if let Some((t, bt)) = rhombus_faces(rhombus) {
                return if top { t } else { bt };
            }
        }
    }
    -1
}

/// Geographic degrees to a planar [`PlanarPoint`] in 5x6 space (eC `forward`,
/// `ri5x6.ec:394`). It applies the vertex-2 azimuth offset, the geodetic to authalic
/// conversion when enabled, finds the face, and then runs the equal-area face
/// projection.
pub(crate) fn forward(geom: &IcoGeometry, lat: f64, lon: f64) -> PlanarPoint {
    // vertex2Azimuth (ri5x6.ec:397)
    let lon_rad = (lon - geom.azimuth_deg) * math::DEG2RAD;
    let lat_rad = lat * math::DEG2RAD;
    let lat_conv = if geom.authalic {
        geodetic_to_authalic(geom, lat_rad)
    } else {
        lat_rad
    };
    let v3d = geo_to_dggal_cart(lat_conv, lon_rad);
    let Some(face) = find_face(geom, &v3d) else {
        // The eC clears the point to (0, 0) and reports failure, which its caller ignores
        // (`ri5x6.ec:405-406` and `RI7H.ec:315`); the face is left to be derived.
        return PlanarPoint {
            face: -1,
            x: 0.0,
            y: 0.0,
        };
    };
    let [x, y] = forward_ico_face(geom, face, &v3d);
    PlanarPoint {
        face: face as i32,
        x,
        y,
    }
}

/// The eC's `epsilon5x6` (`ri5x6.ec:465`), verbatim: how near a pole's position in the
/// 5x6 layout, on each axis, a planar point must lie for [`fix_poles`] to move it onto
/// the pole. The eC's own comment records that 1e-5 had caused imprecisions of about
/// 111 m near the geographic poles.
const EPSILON_5X6: f64 = 1e-8;

/// eC `fixPoles` (`ri5x6.ec:463-509`): the geographic point in radians that a planar
/// point within [`EPSILON_5X6`] of a pole's 5x6 position is moved to, or `None` for any
/// other point.
///
/// The pole appears at four places in the 5x6 layout, two for each pole: (1.5, 3) and
/// (2, 3.5) for the south, (5, 4.5) and (0.5, 0) for the north. A point near one of them
/// becomes exactly that pole, at longitude `qOffset + lon1 + 180 * add180`, where `lon1`
/// is `-180 - orientation.lon`, `qOffset` is 0 on an odd grid and 90 on an even one, and
/// `add180` depends on the side of the position the point lies and, at two of the four,
/// on the grid's parity. On the canonical orientation that gives -168.8 or 11.2 on an odd
/// grid and -78.8 or 101.2 on an even one.
///
/// The eC keeps, commented out, an older way of deciding `oddGrid` from the longitude
/// (`ri5x6.ec:468-496`); it is not ported, since the live code takes the flag from the
/// caller.
fn fix_poles(geom: &IcoGeometry, vx: f64, vy: f64, odd_grid: bool) -> Option<(f64, f64)> {
    // faithful: the eC's `Degrees` in this function hold radians once compiled, and the
    // arithmetic is done on them. In DGGAL v0.0.6 as built with `-O2 -ffast-math`
    // (`Makefile.dggal:133`), BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`,
    // `fixPoles` at `0x3cb20` forms `lon1` as `fl(-PI)` less the orientation's longitude
    // in radians (`0x3cb3c` to `0x3cb50`), takes `qOffset` as 0 or `fl(PI / 2)` and the
    // latitude as `fl(PI / 2)` or its negation from `.rodata`, and converts `180 * add180`
    // at run time through ecrt's `Degrees` to `Radians` property (the call at `0x3cbdd`),
    // which multiplies by `fl(PI / 180)`. It then adds, re-associated, `(that + qOffset) +
    // lon1` (`0x3cbea`, `0x3cbee`), where the eC writes `qOffset + lon1 + (add180 * 180)`
    // (`ri5x6.ec:508`). This port formed the longitude in degrees and converted the sum,
    // 1.776e-15 degrees from the engine's at both pole images. The orientation's
    // longitude in radians is the product [`build_geometry`] gives [`compute_vertices`],
    // which is the engine's own on the canonical orientation.
    let lon1 = -math::PI - geom.orientation_lon_deg * math::DEG2RAD;
    let q_offset = if odd_grid { 0.0 } else { math::PI / 2.0 };
    let near = |x: f64, y: f64| (vx - x).abs() < EPSILON_5X6 && (vy - y).abs() < EPSILON_5X6;
    let (north, add180) = if near(1.5, 3.0) {
        (false, vx > 1.5)
    } else if near(2.0, 3.5) {
        (false, (vy > 3.5) ^ odd_grid)
    } else if near(5.0, 4.5) {
        (true, vy < 4.5)
    } else if near(0.5, 0.0) {
        (true, (vx < 0.5) ^ odd_grid)
    } else {
        return None;
    };
    let lat = if north {
        math::PI / 2.0
    } else {
        -(math::PI / 2.0)
    };
    let turn = math::radians(if add180 { 180.0 } else { 0.0 });
    Some((lat, (turn + q_offset) + lon1))
}

/// A planar [`PlanarPoint`] back to a geographic [`GeoPoint`] in degrees: the eC's `inverse`
/// (`ri5x6.ec:511-556`), which is [`inverse_radians`], with each member of its answer
/// multiplied by `fl(180 / Pi)`, as ecrt converts an angle for the engine's callers.
///
/// A point for which the eC's function fails, because no face can be derived for it, comes
/// back as the zero point, to which the eC clears its result before it returns `false`
/// (`ri5x6.ec:554-555`); so does a point that names a face above 19, which the icosahedron
/// does not have.
pub(crate) fn inverse(geom: &IcoGeometry, p: PlanarPoint, odd_grid: bool) -> GeoPoint {
    match inverse_radians(geom, p, odd_grid) {
        Some((lat, lon)) => GeoPoint {
            lat: lat * math::RAD2DEG,
            lon: lon * math::RAD2DEG,
        },
        None => GeoPoint { lat: 0.0, lon: 0.0 },
    }
}

/// A planar [`PlanarPoint`] back to a geographic point in radians, latitude then longitude
/// (eC `inverse`, `ri5x6.ec:511-556`), or `None` where the eC's function returns `false`.
/// A negative `face` means "derive it from x and y". The edge-wrap
/// fixups come first, then the inverse face projection, then the snap onto the poles,
/// the azimuth offset, the authalic to geodetic conversion and the wrap of the
/// longitude, in the eC's order (`ri5x6.ec:544-550`).
///
/// The eC's function answers in radians and says whether it succeeded, and this is that
/// function as it stands. [`inverse`] converts it for the callers that want degrees and
/// take the zero point for a failure. The refined boundary of a zone wants neither: the
/// eC's `getRefinedVertices` (`RI7H.ec:536-612`, `RI3H.ec:451-673`) leaves out a point whose
/// inverse fails, and compares neighbouring points in radians.
///
/// It fails for a point whose face cannot be derived, because it lies outside the 5x6
/// layout, as in the eC (`ri5x6.ec:533-534`); and for a point that names a face above 19,
/// which the icosahedron does not have.
///
/// `odd_grid` is the eC's own `oddGrid` argument. It says whether the point belongs to a
/// grid of odd resolution, and it matters only to [`fix_poles`], which uses it to choose
/// the longitude of a point it moves onto a pole.
///
/// The snap onto the poles is ported, although py4dggs leaves it out. py4dggs's premise
/// was that the topology never feeds an exact pole point to the projection, but
/// `fixPoles` does not wait for an exact one: it moves any point within 1e-8 of a pole's
/// 5x6 position, about 6e-7 degrees. On IGEO7 the oracle suite measured it catching a
/// single coordinate at resolution 17, one vertex of each polar cell at 18, and at 19
/// the centroid of each polar cell and, for some cells, every vertex; the live engine
/// reports all of them exactly on the pole. Without the snap this port placed those at
/// 18 and 19 from 3e-7 to 7e-7 degrees away. The omission was also the mechanism
/// behind the IVEA and RTEA pole and vertex-0-meridian conformance exceptions recorded
/// in py4dggs.
///
/// The azimuth is added in radians before the longitude is wrapped, as the eC adds
/// `vertex2Azimuth`, and not in degrees after the wrap, as py4dggs does; so a grid built
/// with any azimuth answers a longitude within half a turn of zero, where py4dggs's
/// order could answer one beyond 180 degrees. The bound is the eC's own: `wrapLon`
/// leaves a longitude no further than `Pi + 10 * DBL_EPSILON` from zero, so one that
/// lands within five units in the last place beyond `fl(Pi)` stands, and is answered as
/// at most 180.0000000000001 degrees. At the azimuth of 0 that the six registered
/// grids use, the two orders give the same doubles.
// The first two edge-wrap branches below do share a body; they are kept apart because
// the eC writes them apart, and collapsing them would hide that the third branch, which
// py4dggs drops, sits between the second and the rest of the chain. The fourth branch tests
// the first's condition again (`vy < 0 && vx < 0`), so it is never taken; it stays because the
// eC writes it there (ri5x6.ec:516-531), and Clippy's `ifs_same_cond` is allowed for that reason.
#[allow(clippy::if_same_then_else, clippy::ifs_same_cond)]
pub(crate) fn inverse_radians(
    geom: &IcoGeometry,
    p: PlanarPoint,
    odd_grid: bool,
) -> Option<(f64, f64)> {
    let mut vx = p.x;
    let mut vy = p.y;
    // The eC's eight edge-wrap branches, in its own order (ri5x6.ec:516-531).
    // faithful: py4dggs collapses the first two into `x < 0 and y < 1`, which is
    // equivalent, but omits the third branch entirely, so a seam point with x below 0
    // and y just above 1 falls through to "outside the layout" there.
    let wrap = if vx < 0.0 && vy < 0.0 {
        vx += 5.0;
        vy += 5.0;
        Wrap::None
    } else if vx < 0.0 && vy < 1.0 {
        vx += 5.0;
        vy += 5.0;
        Wrap::AddFive
    } else if vx < 0.0 && vy < 1.0 + 1e-10 {
        vx += 5.0;
        vy = 6.0;
        Wrap::AddFive
    } else if vy < 0.0 && vx < 0.0 {
        vx += 5.0;
        vy += 5.0;
        Wrap::None
    } else if vy < 0.0 && vx < 1e-10 {
        vx = 5.0;
        vy += 5.0;
        Wrap::None
    } else if vx > 5.0 && vy > 5.0 {
        vx -= 5.0;
        vy -= 5.0;
        Wrap::SubtractFiveFromBoth
    } else if vx > 5.0 && vy > 5.0 - 1e-10 {
        vx -= 5.0;
        vy = 0.0;
        Wrap::SubtractFiveOrdinateZero
    } else if vy > 6.0 && vx > 5.0 - 1e-10 {
        vy -= 5.0;
        vx = 0.0;
        Wrap::None
    } else {
        Wrap::None
    };

    let mut face = p.face;
    if face < 0 {
        face = face_after_wrap(vx, vy, wrap, p.x);
    }
    // A face above 19 names no face of the icosahedron, and it is answered as a face
    // that cannot be derived is, with a failure. `face_after_wrap` never returns such a
    // face, so no faithful answer changes; the bound only keeps a caller's out-of-range
    // face from indexing past the twenty-entry face tables below.
    if !(0..=19).contains(&face) {
        return None;
    }

    let p3d = inverse_ico_face(geom, face as usize, &[vx, vy]);
    let [mut lat, mut lon] = dggal_to_geo(p3d[0], p3d[1], p3d[2]);
    // faithful: fixPoles runs on the wrapped 5x6 point, after cartesianToGeo and before
    // the azimuth and the authalic conversion (ri5x6.ec:544-549).
    if let Some((la, lo)) = fix_poles(geom, vx, vy, odd_grid) {
        lat = la;
        lon = lo;
    }

    // faithful: the eC adds `vertex2Azimuth` to the longitude in radians, after
    // `fixPoles` and before the latitude is converted and the longitude wrapped
    // (`ri5x6.ec:547-550`), and so does DGGAL v0.0.6 (BuildID
    // `e75d6ab18f8b460713ebcd21df6f07922fb909c3`), where this unguarded function performs
    // it as written: one `addsd` of the member at offset `0x1a0` (`0x3d0d7` to `0x3d0e7`),
    // which the constructor sets to zero (`0x3a4fb`, `ri5x6.ec:105`) and the forward
    // subtracts in the same unit (`0x3c98b`), then the authalic conversion (`0x3d10e`) and
    // `wrapLon` (`0x3d120`). The configuration holds the azimuth in degrees, and it is
    // converted as every angle of this crate is, by the product with `DEG2RAD`. py4dggs
    // wraps first and adds the azimuth in degrees afterwards, which at any azimuth but zero
    // can answer a longitude beyond 180 degrees.
    lon += math::radians(geom.azimuth_deg);

    let lat_geo = if geom.authalic {
        authalic_to_geodetic(geom, lat)
    } else {
        lat
    };

    Some((lat_geo, wrap_lon(lon)))
}

/// A longitude in radians brought back within half a turn of zero, the eC's `wrapLon`
/// (`GeoExtent.ec:26-33`), which `inverse` applies last (`ri5x6.ec:550`). A longitude no
/// further from zero than `Pi + 10 * DBL_EPSILON`, or NaN, is returned untouched.
///
/// It changes no answer on the six registered grids, because neither branch is ever
/// entered there: the longitude that reaches it is an `atan2` result, within
/// `[-Pi, Pi]`, or one of the four longitudes `fix_poles` assigns, all within the same
/// range, and the azimuth the eC adds before the wrap is zero on all six. With both
/// branches made to panic, none did over the six oracle suites and every library test
/// but the two that drive the wrap on purpose, one directly and one through a grid
/// turned by an azimuth; nor over the centroids and vertices of 1,819,645 zones on the
/// three aperture-7 grids and the neighbours of the zones among them that were
/// quantised from points on the edges of the icosahedron and near the poles, or over
/// 12,002,940 inverse projections of random planar points and of points at and beside
/// the poles' planar images and the layout's corners.
///
/// The port is the library's function on every input, not only on those: called
/// directly at its address in the loaded library, the library's `wrapLon` and this
/// function answered the same doubles, to the last bit, at 356,333 inputs, 283,975 of
/// which enter a branch: the thresholds and `Pi` itself, unit by unit in the last place,
/// every multiple of `Pi` up to `4000 * Pi` and three units either side, uniform draws
/// over ranges from 4 to 1e300 wide, and random bit patterns. The source's own form,
/// with its division and `floor`, answers otherwise at 6,639 of them.
fn wrap_lon(x: f64) -> f64 {
    // faithful: DGGAL v0.0.6 as built with `-O2 -ffast-math` (`Makefile.dggal:133`),
    // BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, `wrapLon` at `0x43050`, called
    // out of line from `inverse` at `0x3d120`, which replaces `GeoExtent.ec:26-33`. The
    // function carries no attribute. The thresholds `-Pi - radEpsilon` and
    // `Pi + radEpsilon` are folded into the constants `0xc00921fb54442d1d` at `0x4a540`
    // and `0x400921fb54442d1d` at `0x4a550`, which are exact, `10 * DBL_EPSILON` being
    // five units in the last place of `Pi`; the rest gcc rewrote in two ways. The division
    // by `2*Pi` becomes a multiplication by `fl(1 / (2*Pi))`, `0x3fc45f306dc9c883` at
    // `0x4a548` (`0x4307a` and `0x430c6`). And `floor` becomes a truncation towards zero,
    // `cvttsd2si` then `cvtsi2sd`, taken only when the quotient's magnitude is below 2^52
    // (`0x4a090`), above which a double is already whole and is used as it stands
    // (`0x43082` to `0x4308a`, `0x430ce` to `0x430d6`). The truncation is the floor
    // because the quotient is positive in both branches. The count of turns multiplies
    // `2*Pi` at `0x49f90`, and the product is added to or subtracted from the longitude.
    const LIMIT: f64 = f64::from_bits(0x4009_21fb_5444_2d1d);
    if x < -LIMIT {
        x + truncated_turns((math::PI - x) * ONE_OVER_TWO_PI) * math::TWO_PI
    } else if x > LIMIT {
        x - truncated_turns((x + math::PI) * ONE_OVER_TWO_PI) * math::TWO_PI
    } else {
        x
    }
}

/// `fl(1 / (2*Pi))`, by which the library multiplies where the eC's `wrapLon` and `wrapLonAt`
/// divide by `2*Pi` (`GeoExtent.ec:29`, `:31`, `:39`, `:41`); `.rodata` at `0x4a548`.
const ONE_OVER_TWO_PI: f64 = f64::from_bits(0x3fc4_5f30_6dc9_c883);

/// 2^52, at and above which every double is whole: there the library's inlined `floor` of
/// `wrapLon` and `wrapLonAt` (`GeoExtent.ec:29`, `:31`, `:39`, `:41`) takes the count of turns
/// as it stands; `.rodata` at `0x4a090`.
const WHOLE_ABOVE: f64 = 4_503_599_627_370_496.0;

/// `fl(Pi + 1e-7)`, the eC's `Pi + Radians { epsilon }` of `wrapLonAt` (`GeoExtent.ec:38`,
/// `:40`, with `epsilon` at `:13`), folded into one constant; `.rodata` at `0x4a560`, and its
/// negation, the folded `-Pi - Radians { epsilon }`, at `0x4a558`.
const WRAP_LON_AT_LIMIT: f64 = math::PI + 1e-7;

/// A count of turns truncated towards zero, as the library truncates it in `wrapLon` and
/// `wrapLonAt` where the eC writes `floor` (`GeoExtent.ec:29`, `:31`, `:39`, `:41`): a
/// `cvttsd2si` to 64 bits and a `cvtsi2sd` back, taken only where the magnitude is below
/// 2^52, above which the count, or a NaN, is used as it stands.
fn truncated_turns(t: f64) -> f64 {
    if t.abs() < WHOLE_ABOVE {
        t as i64 as f64
    } else {
        t
    }
}

/// The eC's `wrapLonAt` with `q == -1` (`GeoExtent.ec:35-53`): how far the longitude `lon`
/// lies from the longitude `c_lon`, both in radians, brought within half a turn of zero. It
/// answers that offset and does not add `c_lon` back, which the eC's own callers remark and
/// do themselves (`RI7H.ec:573`, `RI3H.ec:517`). An offset no further from zero than
/// `Pi + 1e-7`, or NaN, is returned untouched.
///
/// The engine builds the refined boundary of a zone with it: each point's longitude is
/// unwrapped about the zone's centroid, so that the ring of a zone over the antimeridian is
/// continuous. Both of those callers pass `q == -1`, and the arm for another `q`
/// (`GeoExtent.ec:42-51`, `0x43191` to `0x432cf`), which chooses a turn by the quadrant of
/// `c_lon`, is not ported.
///
/// The port is the library's function on every input: called directly at its address in the
/// loaded library with `q == -1`, the library's `wrapLonAt` and this function answered the
/// same doubles, to the last bit, at 370,208 pairs of a longitude and a centre, of which
/// 123,293 enter the first branch and 123,457 the second: the two thresholds and `Pi`,
/// `3 * Pi` and `5 * Pi` of either sign, each about twelve centres and forty units in the last
/// place either side; every multiple of `Pi` up to `400 * Pi` about four centres, and three
/// units either side; uniform draws within 4 to 1e300 of zero; pairs drawn within half a turn
/// of zero; and random bit patterns. The source's own form answers otherwise at 5,004 of
/// them. None of those is a pair the boundary of a zone presents, a longitude and a centroid
/// each within half a turn of zero, of which the sample holds 150,066, with 33,820 that enter
/// a branch: there one turn is counted on either reading. So neither of the two sites below
/// moves a ring: with either put back in the source's order alone, none of the 413,414 rings
/// sampled on the six grids left the engine's.
pub(crate) fn wrap_lon_at(lon: f64, c_lon: f64) -> f64 {
    // faithful: DGGAL v0.0.6 as built with `-O2 -ffast-math` (`Makefile.dggal:133`),
    // BuildID `e75d6ab18f8b460713ebcd21df6f07922fb909c3`, `wrapLonAt` at `0x43130`, which
    // carries no attribute and is called out of line (on the six grids' path, from
    // `getRefinedVertices` at `0x2aab1` and `0xf789`, among the library's five call sites). The offset `lon - cLon` is formed as written (`0x43140`) and compared
    // with the two folded thresholds (`0x43144`, `0x431b8`). Both branches multiply by
    // `fl(1 / (2*Pi))` where the eC divides by `2*Pi`, and take their `floor` inline, only
    // where the magnitude is below 2^52 (`0x43176`, `0x431ea`). The count of turns
    // multiplies `2*Pi` at `0x49f90`, and the product is added to the offset or subtracted
    // from it (`0x43180` to `0x43188`, `0x431f4` to `0x431fc`, `0x4329e` to `0x432a6`).
    let d = lon - c_lon;
    if d < -WRAP_LON_AT_LIMIT {
        // faithful: site WL-1, `GeoExtent.ec:39`, at `0x4314a` to `0x43188`. The eC's
        // `Pi - lon`, on the offset, is taken on the arguments, re-associated, as
        // `(Pi + cLon) - lon` (`0x4315a`, `0x4315e`), before the product with the
        // reciprocal (`0x43162`). The floor is a true one here: the truncation, less one
        // where it exceeds the count (`0x43208` to `0x4322b`); the count is positive on this
        // branch, so that the truncation never exceeds it and the correction is inert.
        let t = ((math::PI + c_lon) - lon) * ONE_OVER_TWO_PI;
        let truncated = truncated_turns(t);
        let turns = if truncated > t {
            truncated - 1.0
        } else {
            truncated
        };
        d + turns * math::TWO_PI
    } else if d > WRAP_LON_AT_LIMIT {
        // faithful: site WL-2, `GeoExtent.ec:41`, at `0x431c2` to `0x431fc`: `Pi` and the
        // offset are added as written (`0x431da`), the sum multiplies the reciprocal
        // (`0x431de`), and the floor is the bare truncation (`0x43290` to `0x43299`), which
        // is the floor, the count being positive there.
        d - truncated_turns((math::PI + d) * ONE_OVER_TWO_PI) * math::TWO_PI
    } else {
        d
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The eC computes the scalar triple product `A . (B x C)` of the frame at run
    /// time, which is what [`inverse_vector`] does; py4dggs substitutes the closed form
    /// [`parallelepiped_v`]. Pin the equality of the two over the 120 frames of each of
    /// the three variants, 360 in all. It shows at the same time that every frame is
    /// right-handed.
    #[test]
    fn parallelepiped_v_matches_every_frame() {
        let expected = parallelepiped_v();
        for radial in [RadialVertex::Isea, RadialVertex::Ivea, RadialVertex::Rtea] {
            let geom = build_geometry(&GridConfig::default(), radial);
            for face in 0..20 {
                for sub_tri in 0..6 {
                    let (v3d, _, _) = sub_tri_frame(&geom, face, sub_tri);
                    let v = dot(&v3d[0], &cross(&v3d[1], &v3d[2]));
                    assert!(
                        v > 0.0 && (v - expected).abs() < 1e-15,
                        "{radial:?} face {face} sub_tri {sub_tri}: {v}"
                    );
                }
            }
        }
    }

    /// The variant constants must satisfy the eC's own identity for the sub-triangle
    /// area, `beta + gamma + alpha - Pi = 6 degrees`
    /// (`icoVertexGreatCircle.ec:243`), for every variant.
    #[test]
    fn variant_angles_sum_to_the_sub_triangle_area() {
        for radial in [RadialVertex::Isea, RadialVertex::Ivea, RadialVertex::Rtea] {
            let c = variant_consts(radial);
            assert!(
                (c.beta + c.gamma + c.alpha - math::PI - SDT_AREA).abs() < 1e-15,
                "{radial:?}"
            );
        }
    }

    /// With an azimuth other than zero, every longitude the inverse answers stays within
    /// half a turn of zero, because the azimuth is added before the wrap, as the eC adds
    /// it. The sample runs over the three projections, seven azimuths from -179 to 179
    /// degrees, and points on the antimeridian, 1e-7 and 0.1 degrees from it on either
    /// side, and at three other longitudes, at five latitudes, both through
    /// `inverse(forward(..))`, for both parities of the pole snap, and through the
    /// centroids and vertices of the zones of a `Grid` built on that azimuth: 27,396
    /// longitudes, of which 6,080 lie within 0.01 degrees of the antimeridian, and at
    /// least 1,000 must. With the azimuth added after the wrap, in degrees, 10,948 of them
    /// lay beyond the bound, as far as 183.1 degrees.
    ///
    /// The bound is the eC's, not 180: `wrapLon` leaves standing a longitude up to five
    /// units in the last place beyond `fl(Pi)`, `180.0000000000001` degrees, and 140 of
    /// the sample land there, from inputs on the antimeridian itself, whose round trip
    /// comes back a few units in the last place past it. The engine, given the same
    /// azimuth, would answer the same: its inverse before the azimuth is this crate's to
    /// the last bit, the addition is one `addsd`, and its `wrapLon` leaves those values
    /// as `the_longitude_wrap_is_the_engines` pins.
    #[test]
    fn a_turned_grid_answers_longitudes_within_half_a_turn() {
        use crate::grid::Grid;
        use crate::indexings::Z7;
        use crate::topologies::HexA7;
        use crate::{Indexing, Projection, Topology};

        fn zone_longitudes<P: Projection, T: Topology, I: Indexing>(
            grid: &Grid<P, T, I>,
            lat: f64,
            lon: f64,
            out: &mut Vec<f64>,
        ) {
            for res in [0u8, 3, 9, 15] {
                let zone = grid.zone_from_geo(lat, lon, res).unwrap();
                out.push(zone.centroid().lon);
                out.extend(zone.vertices().iter().map(|v| v.lon));
            }
        }

        let azimuths = [30.0, -30.0, 90.0, -90.0, 179.0, -179.0, 45.5];
        let lats = [-60.0, -1.0, 0.0, 33.3, 75.0];
        let lons = [
            -180.0,
            -179.9999999,
            -179.9,
            -150.0,
            0.0,
            100.0,
            179.9,
            179.9999999,
            180.0,
        ];
        // The furthest from zero that `wrapLon` leaves a longitude, in degrees.
        let bound = f64::from_bits(0x4009_21fb_5444_2d1d) * math::RAD2DEG;
        let (mut seen, mut near, mut beyond) = (0usize, 0usize, 0usize);
        for az in azimuths {
            let config = GridConfig {
                azimuth_deg: az,
                ..GridConfig::default()
            };
            let isea = Grid::<crate::projections::Isea, HexA7, Z7>::new(config, "IGEO7").unwrap();
            let ivea = Grid::<crate::projections::Ivea, HexA7, Z7>::new(config, "IVEA7H").unwrap();
            let rtea = Grid::<crate::projections::Rtea, HexA7, Z7>::new(config, "RTEA7H").unwrap();
            for lat in lats {
                for lon in lons {
                    let mut out = Vec::new();
                    for radial in [RadialVertex::Isea, RadialVertex::Ivea, RadialVertex::Rtea] {
                        let g = build_geometry(&config, radial);
                        for odd_grid in [false, true] {
                            out.push(inverse(&g, forward(&g, lat, lon), odd_grid).lon);
                        }
                    }
                    zone_longitudes(&isea, lat, lon, &mut out);
                    zone_longitudes(&ivea, lat, lon, &mut out);
                    zone_longitudes(&rtea, lat, lon, &mut out);
                    for l in out {
                        assert!(
                            l.abs() <= bound,
                            "azimuth {az}, ({lat}, {lon}): a longitude of {l}, beyond the wrap"
                        );
                        seen += 1;
                        near += usize::from(l.abs() > 179.99);
                        beyond += usize::from(l.abs() > 180.0);
                    }
                }
            }
        }
        assert!(
            near >= 1000,
            "only {near} of {seen} longitudes near the antimeridian"
        );
        assert_eq!(
            (seen, beyond),
            (27_396, 140),
            "the sample, and the longitudes past 180"
        );
    }

    /// The sides of each variant are the constants the library's `radialVertex` setter
    /// at `0x3e410` stores from `.rodata`, as read there by address, in the roles it
    /// stores them in: ISEA's `AB`, `AC` and `BC` at `0x49e28`, `0x49e30` and `0x49e38`,
    /// IVEA's at `0x49e58`, `0x49e60` and `0x49e68`, RTEA's at `0x49e88`, `0x49e90` and
    /// `0x49e98`. On x86-64 Linux with glibc they are also what the eC's expressions give
    /// at run time, which is why pinning them changes no answer there.
    #[test]
    fn the_sides_are_the_engines_constants() {
        #[rustfmt::skip]
        let rodata: [(RadialVertex, [u64; 3]); 3] = [
            (RadialVertex::Isea, [0x3fd7_59ed_d04f_68da, 0x3fe1_b6e1_92eb_be44, 0x3fe4_e01e_2d74_e77d]),
            (RadialVertex::Ivea, [0x3fe1_b6e1_92eb_be44, 0x3fd7_59ed_d04f_68da, 0x3fe4_e01e_2d74_e77d]),
            (RadialVertex::Rtea, [0x3fe1_b6e1_92eb_be44, 0x3fe4_e01e_2d74_e77d, 0x3fd7_59ed_d04f_68da]),
        ];
        for (radial, sides) in rodata {
            let c = variant_consts(radial);
            assert_eq!(
                [c.ab.to_bits(), c.ac.to_bits(), c.bc.to_bits()],
                sides,
                "{radial:?}"
            );
        }
        #[cfg(all(target_arch = "x86_64", target_os = "linux", target_env = "gnu"))]
        {
            let phi = phi();
            let at_run_time = [
                math::acos(math::sqrt((phi + 1.0) / 3.0)),
                math::atan(1.0 / phi),
                math::atan(2.0 / (phi * phi)),
            ];
            assert_eq!(
                at_run_time.map(f64::to_bits),
                [CENTRE_TO_MIDPOINT, MIDPOINT_TO_CORNER, CENTRE_TO_CORNER].map(f64::to_bits)
            );
        }
    }

    /// The north pole lies on the edge from base cell 0 to 1, which the 5x6 layout cuts,
    /// so it has two planar images, one near (0.5, 0) on face 0 and one near (5, 4.5) on
    /// face 4. The eC takes the first face in its order, and so face 0; so does the
    /// engine, whose forward projections (`ISEAProjection`, `IVEAProjection` and
    /// `RTEAProjection` of DGGAL 0.0.6) give the points recorded here. Taking the face
    /// with the nearest centroid had put the pole on face 4.
    ///
    /// The engine's points were read with the geographic point set in radians, as
    /// [`crate::Grid::zone_from_geo`] and the oracle set it. They must not be read
    /// through the binding's degree setter, whose constant is `0.0174532925199433`
    /// rather than `PI / 180` (`site-packages/ecrt/ecrt.py:3400-3401`): that setter puts
    /// 90 degrees 4.44e-16 radians past the pole and moves every one of these points.
    ///
    /// All three are the engine's to the last bit. Before the projection was put into the
    /// order the compiled engine performs it, IVEA's ordinate was zero where the engine's
    /// is 3.700743415417188e-17, and RTEA's abscissa nine units in the last place from the
    /// engine's.
    #[test]
    fn the_north_pole_is_projected_on_the_face_the_engine_takes() {
        let cases: [(RadialVertex, (f64, f64)); 3] = [
            (
                RadialVertex::Isea,
                (0.49999999999999867, -1.4802973661668753e-16),
            ),
            (
                RadialVertex::Ivea,
                (0.4999999999999989, 3.700743415417188e-17),
            ),
            (RadialVertex::Rtea, (0.4999999999999988, 0.0)),
        ];
        for (radial, engine) in cases {
            let g = build_geometry(&GridConfig::default(), radial);
            let p = forward(&g, 90.0, 0.0);
            assert_eq!(p.face, 0, "{radial:?}");
            assert_eq!(
                (p.x.to_bits(), p.y.to_bits()),
                (engine.0.to_bits(), engine.1.to_bits()),
                "{radial:?}: {p:?} against the engine's {engine:?}"
            );
        }
    }

    /// The zone at the north pole at every resolution, as the engine answers it on each
    /// of the three grids (DGGAL 0.0.6, `getZoneFromWGS84Centroid`, the same on all
    /// three): base cell 1, then the digits 3 and 2 in turn. With the face taken as the
    /// nearest centroid, RTEA7H answered base cell 0 instead at every resolution but 1.
    ///
    /// Two answers of RTEA7H differed from the engine's after that, and neither does now.
    /// Resolution 11 went when the vertices of the 5x6 CRS were ported from the library
    /// (see `topologies::hex_a7::add_non_polar_vertices_refined` and its fellows), and
    /// resolution 19, where each side's quantiser, given the other side's planar point,
    /// answered as the other side does, went when the projection was: the pole's planar
    /// point lay nine units in the last place from the engine's, on a cell boundary, and
    /// is now the engine's own (see the test above).
    #[test]
    fn the_north_pole_quantises_as_the_engine_does() {
        let digits = "3232323232323232323";
        let engine = |res: usize| format!("01{}", &digits[..res]);
        for res in 0..=19usize {
            let r = res as u8;
            let want = engine(res);
            let igeo7 = crate::igeo7()
                .zone_from_geo(90.0, 0.0, r)
                .unwrap()
                .text_id();
            let ivea7h = crate::ivea7h()
                .zone_from_geo(90.0, 0.0, r)
                .unwrap()
                .text_id();
            let rtea7h = crate::rtea7h()
                .zone_from_geo(90.0, 0.0, r)
                .unwrap()
                .text_id();
            assert_eq!(igeo7, want, "IGEO7 at {res}");
            assert_eq!(ivea7h, want, "IVEA7H at {res}");
            assert_eq!(rtea7h, want, "RTEA7H at {res}");
        }
    }

    /// The longitude wrap answers as the library's compiled `wrapLon` at `0x43050` does,
    /// called directly at that address, at points where the source's text or py4dggs's
    /// wrap would answer otherwise. At 65.97344572538564 the source's division and
    /// `floor` give 3.1415926535897825, and the library, multiplying by `fl(1 / (2*Pi))`
    /// and truncating, counts one turn more. One unit in the last place above `fl(Pi)`,
    /// py4dggs subtracts a turn and the library leaves the value alone, as it does up to
    /// its threshold and not a unit beyond.
    #[test]
    fn the_longitude_wrap_is_the_engines() {
        let cases: [(u64, u64); 7] = [
            (0x4050_7e4c_ef4c_bd97, 0xc009_21fb_5444_2d20),
            (0xc050_7e4c_ef4c_bd97, 0x4009_21fb_5444_2d20),
            (0x4009_21fb_5444_2d19, 0x4009_21fb_5444_2d19),
            (0x4009_21fb_5444_2d1d, 0x4009_21fb_5444_2d1d),
            (0x4009_21fb_5444_2d1e, 0xc009_21fb_5444_2d12),
            (0xc009_21fb_5444_2d1e, 0x4009_21fb_5444_2d12),
            (0x4010_0000_0000_0000, 0xc002_43f6_a888_5a30),
        ];
        for (x, engine) in cases {
            let x = f64::from_bits(x);
            let ours = wrap_lon(x);
            assert_eq!(
                ours.to_bits(),
                engine,
                "{x:e}: {ours:e} against the engine's {:e}",
                f64::from_bits(engine)
            );
        }
        assert!(wrap_lon(f64::NAN).is_nan());
    }

    /// The constants of the two wraps are the library's own, read from its `.rodata`: the
    /// threshold of `wrapLonAt` and its negation, which are `fl(Pi + 1e-7)` and its mirror,
    /// the reciprocal by which both wraps multiply, and the bound above which a double is
    /// whole. They pin WL-1 and WL-2, the two compiled sites of `wrapLonAt`, which move no
    /// ring of the six grids and so cannot be pinned by one.
    #[test]
    fn the_constants_of_the_longitude_wraps_match_the_engines_rodata() {
        assert_eq!(WRAP_LON_AT_LIMIT.to_bits(), 0x4009_21fb_61b0_2665); // 0x4a560
        assert_eq!((-WRAP_LON_AT_LIMIT).to_bits(), 0xc009_21fb_61b0_2665); // 0x4a558
        assert_eq!(ONE_OVER_TWO_PI.to_bits(), 0x3fc4_5f30_6dc9_c883); // 0x4a548
        assert_eq!(ONE_OVER_TWO_PI, 1.0 / math::TWO_PI);
        assert_eq!(WHOLE_ABOVE.to_bits(), 0x4330_0000_0000_0000); // 0x4a090
        assert_eq!(math::PI.to_bits(), 0x4009_21fb_5444_2d18); // 0x4a040
        assert_eq!(math::TWO_PI.to_bits(), 0x4019_21fb_5444_2d18); // 0x49f90
    }

    /// The wrap about a centre answers as the library's compiled `wrapLonAt` at `0x43130`
    /// does, called directly at that address with `q == -1`: each row is a longitude, a
    /// centre and the engine's answer, as bits. The first four are the threshold and one
    /// unit in the last place beyond it, on either side. The next four are points of WL-1
    /// where the source's text answers the other end of the half turn: two where its grouping
    /// `Pi - (lon - cLon)` and the library's `(Pi + cLon) - lon` give different doubles, and two
    /// where its division by `2*Pi` and the library's product with `fl(1 / (2*Pi))` count
    /// different turns. The next two are the same division at WL-2. The last three are
    /// ordinary: a longitude and a centre on either side of the antimeridian, both ways, and a
    /// pair that needs no wrap.
    #[test]
    fn the_longitude_wrap_about_a_centre_is_the_engines() {
        let cases: [(u64, u64, u64); 13] = [
            (0x4009_21fb_61b0_2665, 0, 0x4009_21fb_61b0_2665),
            (0x4009_21fb_61b0_2666, 0, 0xc009_21fb_46d8_33ca),
            (0xc009_21fb_61b0_2665, 0, 0xc009_21fb_61b0_2665),
            (0xc009_21fb_61b0_2666, 0, 0x4009_21fb_46d8_33ca),
            (
                0xc022_d9ee_e0f7_30af,
                0xbf4c_9871_03b7_633a,
                0x4009_21fb_5444_2d1c,
            ),
            (
                0xc019_239c_d570_8113,
                0x4009_1eb8_51eb_851f,
                0x4009_21fb_5444_2d1c,
            ),
            (
                0xc04f_6aae_597a_c2dc,
                0x4009_1eb8_51eb_851f,
                0x4009_21fb_5444_2d20,
            ),
            (0xc050_7e4c_ef4c_bd97, 0, 0x4009_21fb_5444_2d20),
            (
                0x404f_8966_ab66_47fa,
                0xc007_3333_3333_3333,
                0xc009_21fb_5444_2d20,
            ),
            (0x4050_7e4c_ef4c_bd97, 0, 0xc009_21fb_5444_2d20),
            (
                0xc008_cccc_cccc_cccd,
                0x4008_cccc_cccc_cccd,
                0x3fb5_4ba1_ddd8_12c0,
            ),
            (
                0x4008_cccc_cccc_cccd,
                0xc008_cccc_cccc_cccd,
                0xbfb5_4ba1_ddd8_12c0,
            ),
            (
                0x3ff0_0000_0000_0000,
                0x3fe0_0000_0000_0000,
                0x3fe0_0000_0000_0000,
            ),
        ];
        for (lon, c_lon, engine) in cases {
            let (lon, c_lon) = (f64::from_bits(lon), f64::from_bits(c_lon));
            let ours = wrap_lon_at(lon, c_lon);
            assert_eq!(
                ours.to_bits(),
                engine,
                "{lon:e} about {c_lon:e}: {ours:e} against the engine's {:e}",
                f64::from_bits(engine)
            );
        }
        assert!(wrap_lon_at(f64::NAN, 0.0).is_nan());
        assert!(wrap_lon_at(0.0, f64::NAN).is_nan());
    }

    /// The sub-triangle area is the constant the library loads from `0x4a450`, one unit
    /// in the last place above the product a conversion at run time gives, and the double
    /// nearest to the sixteen decimal places the eC compiler writes for `Degrees { 6 }`.
    #[test]
    fn the_sub_triangle_area_is_the_engines_constant() {
        assert_eq!(SDT_AREA.to_bits(), 0x3fba_cee9_f37b_ebd7);
        assert_eq!(SDT_AREA, "0.1047197551196598".parse::<f64>().unwrap());
        assert_eq!(SDT_AREA.to_bits(), (6.0 * math::DEG2RAD).to_bits() + 1);
    }

    /// The icosahedron and the authalic coefficients as the engine caches them, read from
    /// a live `ISEAProjection` of DGGAL v0.0.6 (BuildID
    /// `e75d6ab18f8b460713ebcd21df6f07922fb909c3`) through `RI5x6Projection`'s members:
    /// the twelve vertices, component by component, and the two coefficient sets. The
    /// face planes, the other cached field, were compared in the same way, all 180
    /// components on each of the three projections, and are bit-identical too.
    #[test]
    fn the_geometry_is_the_engines() {
        #[rustfmt::skip]
        const VERTICES: [u64; 36] = [
            0x3fba_2436_1b94_90b2, 0xbfeb_3888_0b46_03dd, 0xbfe0_80c4_6325_e9ab,
            0xbfba_2436_1b94_909a, 0xbfeb_3888_0b46_03ec, 0x3fe0_80c4_6325_e994,
            0xbfea_b3d1_1b5a_b26b, 0xbfe0_d2ca_0da1_530c, 0xbfc5_261c_db24_7664,
            0xbfd6_6e7a_58b9_981d, 0x3cd8_0000_0000_0000, 0xbfed_f857_decd_447f,
            0x3fe5_ca4b_99ef_0730, 0x3cd7_0000_0000_0000, 0xbfe7_6f4a_57e8_2058,
            0x3fea_b3d1_1b5a_b26d, 0xbfe0_d2ca_0da1_530d, 0x3fc5_261c_db24_7622,
            0xbfe5_ca4b_99ef_0730, 0xbcd7_0000_0000_0000, 0x3fe7_6f4a_57e8_2057,
            0xbfea_b3d1_1b5a_b26c, 0x3fe0_d2ca_0da1_530d, 0xbfc5_261c_db24_7626,
            0x3fba_2436_1b94_9092, 0x3feb_3888_0b46_03ec, 0xbfe0_80c4_6325_e994,
            0x3fea_b3d1_1b5a_b26a, 0x3fe0_d2ca_0da1_530d, 0x3fc5_261c_db24_765f,
            0x3fd6_6e7a_58b9_981f, 0xbcd8_0000_0000_0000, 0x3fed_f857_decd_4480,
            0xbfba_2436_1b94_90b2, 0x3feb_3888_0b46_03dd, 0x3fe0_80c4_6325_e9ab,
        ];
        #[rustfmt::skip]
        const COEFFICIENTS: [u64; 12] = [
            0xbf62_57f6_a0d5_0b21, 0x3ec1_dffd_081a_49bd, 0xbe25_fbdd_fed0_2d60,
            0x3d8d_a502_0676_031a, 0xbcf5_0550_38aa_3227, 0x3c5e_c85d_69dd_e8f9,
            0x3f62_57f6_2d10_62c3, 0x3ec8_2f9e_c17c_da47, 0x3e35_d85b_dd44_ac07,
            0x3da6_6f1c_1c0b_0cf1, 0x3d18_abf4_96e6_48da, 0x3c8c_6911_c632_c726,
        ];
        let g = build_geometry(&GridConfig::default(), RadialVertex::Isea);
        for (i, v) in g.ico_vertices.iter().enumerate() {
            for k in 0..3 {
                assert_eq!(
                    v[k].to_bits(),
                    VERTICES[i * 3 + k],
                    "vertex {i}, component {k}"
                );
            }
        }
        for k in 0..6 {
            assert_eq!(
                g.authalic_cp.0[k].to_bits(),
                COEFFICIENTS[k],
                "geodetic to authalic {k}"
            );
            assert_eq!(
                g.authalic_cp.1[k].to_bits(),
                COEFFICIENTS[6 + k],
                "authalic to geodetic {k}"
            );
        }
    }

    /// Latitudes at which the eC's written order of the Clenshaw recurrences and the
    /// library's disagree by a unit in the last place; the expected values are the
    /// engine's own `latGeodeticToAuthalic` and `latAuthalicToGeodetic`.
    #[test]
    fn the_authalic_series_is_summed_as_the_engine_sums_it() {
        let g = build_geometry(&GridConfig::default(), RadialVertex::Isea);
        let g2a = [
            (0x3fec_9404_038e_da78_u64, 0x3fec_8216_a5fb_c91f_u64),
            (0xbfd0_1148_96d7_8720, 0xbfcf_ff50_055f_cc20),
        ];
        let a2g = [
            (0x3fe5_16a7_d1f5_e9a0_u64, 0x3fe5_286d_7c71_de62_u64),
            (0xbfea_9d90_1b1f_8766, 0xbfea_afd2_cb0b_b67e),
        ];
        for (phi, engine) in g2a {
            let got = geodetic_to_authalic(&g, f64::from_bits(phi));
            assert_eq!(got.to_bits(), engine, "geodetic to authalic at {phi:#018x}");
        }
        for (phi, engine) in a2g {
            let got = authalic_to_geodetic(&g, f64::from_bits(phi));
            assert_eq!(got.to_bits(), engine, "authalic to geodetic at {phi:#018x}");
        }
    }

    /// The engine's `Normalize` answers the zero vector only for a vector whose squares
    /// sum to zero, a NaN included; a vector shorter than 1e-15 is normalised. The
    /// expected components are the library's own, called at `0x42460`.
    #[test]
    fn normalize_answers_zero_only_where_the_engine_does() {
        let tiny = normalize(&[1e-160, -1e-161, 3e-162]);
        let engine = [
            0x3fef_d3e2_2075_8ddb,
            0xbfb9_764e_805e_0b16,
            0x3f9e_8df7_cd3d_a6e7,
        ];
        assert_eq!(tiny.map(f64::to_bits), engine);
        let short = normalize(&[1e-16, 2e-16, -3e-17]);
        let engine = [
            0x3fdc_5e14_9002_4f56,
            0x3fec_5e14_9002_4f56,
            0xbfc1_053f_899a_fc67,
        ];
        assert_eq!(short.map(f64::to_bits), engine);
        // The smallest subnormal squares to zero, and so answers zero.
        assert_eq!(normalize(&[f64::from_bits(1), 0.0, 0.0]), [0.0; 3]);
        assert_eq!(normalize(&[0.0; 3]), [0.0; 3]);
        assert_eq!(normalize(&[f64::NAN, 1.0, 0.0]), [0.0; 3]);
    }

    /// The two constants gcc folded into `inverse`'s wrap, as Rust folds the same
    /// expressions: `fl(5 + 1e-11)` at `0x4a188` and `fl(5 + fl(5 - 1e-11))` at `0x4a3b8`.
    #[test]
    fn the_folded_wrap_constants_are_the_engines() {
        let epsilon: f64 = 1e-11;
        assert_eq!((5.0 + epsilon).to_bits(), 0x4014_0000_0000_2bfb);
        assert_eq!((5.0 + (5.0 - epsilon)).to_bits(), 0x4023_ffff_ffff_ea02);
    }

    /// The two abscissae at which the library's folded wrap names another face than the
    /// eC's text: `fl(5 + 1e-11)` on the seventh wrap branch, where the eC nudges the
    /// wrapped abscissa and the library does not, and `fl(5 + fl(5 - 1e-11))` on the
    /// sixth, where the eC nudges it one way and the library the other. The expected
    /// points are the engine's own `RI5x6Projection_inverse` on `ISEAProjection` for an
    /// odd grid, in degrees. At the three points the engine takes faces 0, 9 and 19,
    /// where the eC's text would take faces 5, 4 and 14.
    #[test]
    fn the_wrap_names_the_face_the_engine_names() {
        let g = build_geometry(&GridConfig::default(), RadialVertex::Isea);
        assert_eq!(
            face_after_wrap(
                5.0 + 1e-11 - 5.0,
                0.0,
                Wrap::SubtractFiveOrdinateZero,
                5.0 + 1e-11
            ),
            0
        );
        let ka = 5.0 + (5.0 - 1e-11);
        assert_eq!(
            face_after_wrap(ka - 5.0, 6.0, Wrap::SubtractFiveFromBoth, ka),
            19
        );
        let cases = [
            (
                (5.0 + 1e-11, 5.0),
                (0x404d_32d5_ad57_5f51_u64, 0xc065_1999_9999_9999_u64),
            ),
            ((ka, 10.0), (0x404d_32d5_ad55_0842, 0xc065_1999_9999_f659)),
            ((ka, 11.0), (0xbe03_cf7c_20a7_acc2, 0xc061_22a4_0cb4_1518)),
        ];
        for ((x, y), engine) in cases {
            let q = inverse(&g, PlanarPoint { face: -1, x, y }, true);
            assert_eq!((q.lat.to_bits(), q.lon.to_bits()), engine, "({x}, {y})");
        }
    }

    /// The pole images as the engine answers them on `ISEAProjection`, whose `fixPoles`
    /// forms the longitude in radians: before this port formed it in degrees, each of
    /// these longitudes was a unit or two in the last place from the engine's.
    #[test]
    fn the_poles_are_placed_where_the_engine_places_them() {
        let g = build_geometry(&GridConfig::default(), RadialVertex::Isea);
        let cases = [
            ((0.5, 0.0), true, (90.0, 11.19999999999999)),
            ((0.5, 0.0), false, (90.0, -78.80000000000001)),
            ((5.0, 4.5), true, (90.0, -168.8)),
            ((1.5, 3.0), true, (-90.0, -168.8)),
            ((2.0, 3.5), true, (-90.0, 11.19999999999999)),
            ((0.5 - 1e-9, 0.0), false, (90.0, 101.19999999999999)),
        ];
        for ((x, y), odd, engine) in cases {
            let q = inverse(&g, PlanarPoint { face: -1, x, y }, odd);
            assert_eq!((q.lat, q.lon), engine, "({x}, {y}), odd {odd}");
        }
    }

    /// The inverse in radians is the eC's `inverse` with its own return value
    /// (`ri5x6.ec:511-556`): `None` where the eC answers `false`, which is where no face can be
    /// derived, and the point in radians elsewhere. The inverse in degrees is that one and no
    /// more: the zero point for a failure, and each member's product with `fl(180 / Pi)` for a
    /// success. The sample is a lattice of quarter steps that overruns the 5x6 layout on every
    /// side, with the four pole images, on each of the three projections and under both
    /// parities.
    #[test]
    fn the_inverse_in_degrees_is_the_inverse_in_radians_converted() {
        let mut sample: Vec<(f64, f64)> = vec![(0.5, 0.0), (5.0, 4.5), (1.5, 3.0), (2.0, 3.5)];
        for i in -2..=22 {
            for j in -2..=26 {
                sample.push((f64::from(i) * 0.25, f64::from(j) * 0.25));
            }
        }
        let (mut answered, mut failed) = (0, 0);
        for radial in [RadialVertex::Isea, RadialVertex::Ivea, RadialVertex::Rtea] {
            let geom = build_geometry(&GridConfig::default(), radial);
            for &(x, y) in &sample {
                for odd in [false, true] {
                    let p = PlanarPoint { face: -1, x, y };
                    let degrees = inverse(&geom, p, odd);
                    let got = (degrees.lat.to_bits(), degrees.lon.to_bits());
                    match inverse_radians(&geom, p, odd) {
                        Some((lat, lon)) => {
                            let want = (lat * math::RAD2DEG, lon * math::RAD2DEG);
                            assert_eq!(got, (want.0.to_bits(), want.1.to_bits()), "({x}, {y})");
                            answered += 1;
                        }
                        None => {
                            assert_eq!(got, (0.0_f64.to_bits(), 0.0_f64.to_bits()), "({x}, {y})");
                            failed += 1;
                        }
                    }
                }
            }
        }
        assert_eq!((answered, failed), (1_446, 2_928));

        let g = build_geometry(&GridConfig::default(), RadialVertex::Isea);
        // Outside the layout: no face, so a failure, where the degrees answer the zero point.
        let outside = PlanarPoint {
            face: -1,
            x: 3.0,
            y: 0.5,
        };
        assert_eq!(inverse_radians(&g, outside, true), None);
        // A face the icosahedron does not have is a failure too.
        let no_such_face = PlanarPoint {
            face: 20,
            x: 0.5,
            y: 0.5,
        };
        assert_eq!(inverse_radians(&g, no_such_face, true), None);
        // The north pole's image is the pole itself, in radians.
        let pole = PlanarPoint {
            face: -1,
            x: 0.5,
            y: 0.0,
        };
        let (lat, _) = inverse_radians(&g, pole, true).unwrap();
        assert_eq!(lat.to_bits(), (math::PI / 2.0).to_bits());
    }
}
