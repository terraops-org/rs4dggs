//! Shared value types: geographic and planar points, grid configuration, the
//! packed zone identifier and the digit-path address an `Indexing` encodes it
//! from. Ported from py4dggs's `types.py`.
//!
//! Source: py4dggs's `types.py`; `Address`, which py4dggs passes as a base cell and a list of
//! digits, is this crate's own.

/// A geographic coordinate in decimal degrees (WGS84 latitude/longitude).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GeoPoint {
    /// Latitude in decimal degrees, positive to the north.
    pub lat: f64,
    /// Longitude in decimal degrees, positive to the east.
    pub lon: f64,
}

/// A point in a projection's planar space.
///
/// `face` carries the icosahedron face id for an icosahedral projection, or a
/// projection-specific sentinel: `-1` means "derive the face from `x` and
/// `y`". `x` and `y` are that projection's own planar coordinates (the 5x6
/// oblique grid used by ISEA, IVEA and RTEA).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarPoint {
    /// The icosahedron face, or the sentinel `-1`, "derive the face from `x` and
    /// `y`".
    pub face: i32,
    /// The first planar coordinate, in the projection's own units.
    pub x: f64,
    /// The second planar coordinate, in the projection's own units.
    pub y: f64,
}

/// Orientation and datum handling for an icosahedral grid, mirroring DGGAL.
///
/// The `Default` implementation gives the canonical IGEO7 orientation: the
/// same icosahedron pole placement and azimuth that DGGAL uses for its own
/// IGEO7 grid, on the authalic (equal-area) sphere.
///
/// Marked `#[non_exhaustive]`, so a caller outside this crate builds a custom
/// orientation from [`GridConfig::default`] and assignment to the public fields,
/// rather than a struct literal; a field added later, for datum handling or the
/// like, must not break that caller.
///
/// Its angles must be finite: [`crate::Grid::new`] refuses a configuration in which
/// any of them is NaN or infinite, naming the field, and the check therefore lives
/// there rather than in a constructor of this type, which a caller outside the
/// crate does not use.
#[derive(Clone, Copy, Debug, PartialEq)]
#[non_exhaustive]
pub struct GridConfig {
    /// The pitch, in degrees, of the rotation that orients the icosahedron, DGGAL's
    /// `orientation.lat` (`ri5x6.ec:106`, `:271`). The canonical 31.7174744114611 is
    /// 90 less the authalic latitude of the icosahedron's first vertex,
    /// 58.2825255885389.
    pub orientation_lat_deg: f64,
    /// The yaw, in degrees, of the same rotation, DGGAL's `orientation.lon`. The
    /// canonical -11.2 puts the first vertex on the meridian 11.2 degrees east.
    pub orientation_lon_deg: f64,
    /// DGGAL's `vertex2Azimuth`, in degrees: subtracted from every longitude before
    /// the forward projection and added back by the inverse before it wraps the
    /// longitude (`ri5x6.ec:397`, `:547`), and so a turn of the whole grid about the
    /// polar axis. Zero on the canonical orientation.
    pub azimuth_deg: f64,
    /// Whether latitudes are converted between geodetic, on the WGS84 ellipsoid, and
    /// authalic, on the sphere of the same area, as DGGAL converts them. When false,
    /// a latitude is taken as a latitude on the sphere.
    pub authalic: bool,
}

impl Default for GridConfig {
    fn default() -> Self {
        GridConfig {
            orientation_lat_deg: 31.7174744114611,
            orientation_lon_deg: -11.20,
            azimuth_deg: 0.0,
            authalic: true,
        }
    }
}

/// A DGGRS zone identifier: the packed integer an `Indexing` encodes an
/// `Address` into, and the type every public API that names a zone passes
/// around. Always `u64`, never a signed type: the Z7 packing uses bit 63 at
/// the coarsest resolutions, which does not fit in an `i64`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct ZoneId(pub u64);

impl ZoneId {
    /// DGGAL's `nullZone` sentinel (`dggrs.ec:9`): "no such zone".
    ///
    /// A quantisation that finds no cell for the input point returns this, and
    /// it is a real, reachable outcome rather than only a caller error. On the
    /// three aperture-7 grids it arises in a strip along two of the
    /// icosahedron's edges, the edge from base cell 0 to base cell 1 that runs
    /// over the north pole and the edge from base cell 1 to base cell 6: at the
    /// odd resolutions from 15 upwards the engine's own candidate search finds
    /// no containing cell there, and this crate answers as the engine does. The
    /// strip is narrow but it is not polar. Measured on each grid over a sweep
    /// across both edges, the engine's null answers lie within 7.5e-5 degrees
    /// of the edge and are met at stations along its whole length, tens of
    /// degrees from either pole. It arises also at isolated stations on the
    /// grid's mirror planes, where the input lies exactly on a cell boundary and
    /// the engine finds no cell for it. At resolutions 16 and 18, and at a few
    /// points of the strip's outer margin at 19, this crate answers the null
    /// zone where the engine answers an identifier its own parser cannot read
    /// back; that is a deliberate divergence, and the null zone is what the
    /// engine itself answers at the neighbouring odd resolutions (see "Known
    /// limitations" in the crate README).
    ///
    /// On the three aperture-3 grids the engine answers it near both poles at
    /// resolutions 0 and 1, and at the coarsest resolutions in a narrow band
    /// beside each of the ten icosahedron edges that end at the vertices of the
    /// two polar pentagons: a probe met it there at resolutions 0 to 7, never
    /// further than 2e-6 radians from the edge. This crate answers as the engine
    /// does. Along two seams of the layout, at every resolution from 2 to 21, it
    /// also answers the null zone where the engine answers an identifier its own
    /// parser cannot read back, a deliberate divergence of the same kind as on
    /// the aperture-7 grids, which [`crate::topologies::HexA3`] describes at its
    /// `quantize`.
    ///
    /// In every sample taken there is no point left at which this crate answers
    /// the null zone where the engine names a cell: that class, which held three
    /// points on the aperture-7 grids, is empty in each of the suites' samples,
    /// whose counts of it are pinned at zero on every grid. One case of that
    /// shape does remain below this level, on the aperture-7 planar quantiser at
    /// the corner of the planar layout, which no sample of the grids reaches;
    /// [`crate::Topology::quantize`] records it.
    ///
    /// A caller can recognise it in several ways: by comparing the identifier
    /// with this constant, by its text, which is [`NULL_TEXT`], and through
    /// `is_pentagon`, which is false for it. Its `resolution` reads as 0 on every
    /// grid, where DGGAL reports -1 on the aperture-7 grids and 63 on the
    /// aperture-3 grids, so the resolution alone cannot tell it from a genuine
    /// resolution-0 zone. On every grid it has the null geometry DGGAL gives it
    /// on the aperture-7 grids: the zero point as its centroid, six zero
    /// vertices, and no neighbours, parents or children. On the aperture-3 grids
    /// that is a deliberate departure: the engine gives its null zone there a
    /// centroid that belongs to no cell, vertices about that point, and a parent
    /// chain on which it divides by zero.
    pub const NULL: ZoneId = ZoneId(u64::MAX);
}

/// What DGGAL's `getZoneTextID` prints for [`ZoneId::NULL`].
///
/// This is DGGAL's only signal that a zone is the null sentinel, so the
/// text-id serialisation mirrors it exactly rather than stringifying the
/// sentinel's bit fields into something that would look like an ordinary,
/// decodable id.
pub const NULL_TEXT: &str = "(null)";

/// The largest number of digits an [`Address`] can hold.
///
/// Z7 needs twenty: nineteen for its finest resolution, and a twentieth for
/// the packing level that is representable but has no geometry. The constant
/// was once set higher, on the assumption that the I3H indexing of the
/// aperture-3 grids would need a digit for each of its 33 levels. It does
/// not: that indexing carries the whole cell in `base` and leaves the digits
/// empty, as py4dggs's `indexings/i3h.py` does, so twenty suffices for every
/// grid this crate ports.
pub const ADDRESS_CAPACITY: usize = 20;

/// A zone's address in a topology's own terms: a base cell plus the
/// direction digits that refine it, one per resolution level below the base.
///
/// The digits live in a fixed-size array rather than a `Vec` so that
/// `Address` is `Copy` and needs no allocation; unused array slots beyond
/// `len` are always kept zeroed (`new` and `push` are the only ways to set
/// them), which is what makes the derived `PartialEq`/`Hash` below correct:
/// two addresses with the same `base` and the same number of digits either
/// agree on every digit up to `len` and on every zeroed slot after it, or
/// they already differ within the significant digits.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Address {
    /// The base cell, or for a topology that packs a whole cell into it, that cell.
    pub base: u64,
    len: u8,
    digits: [u8; ADDRESS_CAPACITY],
}

impl Address {
    /// Builds an address from a base cell and its digits, most significant
    /// (coarsest) first.
    ///
    /// # Panics
    ///
    /// Panics if `digits.len()` exceeds [`ADDRESS_CAPACITY`]. This is not a
    /// caller-input bound to validate defensively: every topology's own
    /// aperture and every indexing's own `MAX_RESOLUTION` already bound how
    /// many digits a real address can hold, so a caller that overruns
    /// `ADDRESS_CAPACITY` has an internal bug (a wrong constant, or a caller
    /// building digits that do not correspond to any real resolution), and a
    /// panic reports it at the point of construction rather than letting a
    /// silently truncated address masquerade as a valid one further down the
    /// pipeline.
    pub fn new(base: u64, digits: &[u8]) -> Address {
        assert!(
            digits.len() <= ADDRESS_CAPACITY,
            "Address::new: {} digits exceeds capacity {ADDRESS_CAPACITY}",
            digits.len()
        );
        let mut stored = [0u8; ADDRESS_CAPACITY];
        stored[..digits.len()].copy_from_slice(digits);
        Address {
            base,
            len: digits.len() as u8,
            digits: stored,
        }
    }

    /// The significant digits, most significant first. Never includes the
    /// unused, zeroed tail of the backing array.
    pub fn digits(&self) -> &[u8] {
        &self.digits[..self.len as usize]
    }

    /// The number of significant digits.
    pub fn len(&self) -> usize {
        self.len as usize
    }

    /// Whether the address has no digits at all, i.e. is the base cell
    /// itself at resolution 0.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Appends one more (least significant) digit.
    ///
    /// # Panics
    ///
    /// Panics if the address is already at [`ADDRESS_CAPACITY`] digits, for
    /// the same reason as `new`: this is an internal invariant on how deep a
    /// real address can go, not a bound on caller-supplied input.
    pub fn push(&mut self, d: u8) {
        let len = self.len as usize;
        assert!(
            len < ADDRESS_CAPACITY,
            "Address::push: already at capacity {ADDRESS_CAPACITY}"
        );
        self.digits[len] = d;
        self.len += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The registry builds its grids from these values and relies on their being
    /// finite, which a literal cannot fail to be; the last assertion makes that
    /// reliance explicit.
    #[test]
    fn grid_config_defaults_are_igeo7() {
        let c = GridConfig::default();
        assert_eq!(c.orientation_lat_deg, 31.7174744114611);
        assert_eq!(c.orientation_lon_deg, -11.20);
        assert_eq!(c.azimuth_deg, 0.0);
        assert!(c.authalic);
        assert!(
            [c.orientation_lat_deg, c.orientation_lon_deg, c.azimuth_deg]
                .iter()
                .all(|v| v.is_finite())
        );
    }

    #[test]
    fn address_push_and_digits() {
        let mut a = Address::new(3, &[6, 5]);
        a.push(3);
        assert_eq!(a.base, 3);
        assert_eq!(a.digits(), &[6, 5, 3]);
        assert_eq!(a.len(), 3);
        assert!(Address::new(0, &[]).is_empty());
    }

    #[test]
    fn null_zone_is_all_ones() {
        assert_eq!(ZoneId::NULL.0, 0xFFFF_FFFF_FFFF_FFFF);
        assert_eq!(NULL_TEXT, "(null)");
    }

    /// The unused, zeroed tail of the backing array must never leak into
    /// equality: an address assembled digit by digit via `push` must equal
    /// the same address built in one call to `new`.
    #[test]
    fn address_equality_ignores_unused_slots() {
        let built_at_once = Address::new(1, &[2]);
        let mut built_incrementally = Address::new(1, &[]);
        built_incrementally.push(2);
        assert_eq!(built_at_once, built_incrementally);
    }
}
