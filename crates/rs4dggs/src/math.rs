//! The single entry point for transcendental math and square roots.
//!
//! Every trig and sqrt call in the crate goes through here so that swapping `std` for
//! `libm`, should a WebAssembly build call for it, is a one-file change.
//! `scripts/check.sh` fails if a call appears anywhere else.
pub const PI: f64 = std::f64::consts::PI;
pub const TWO_PI: f64 = 2.0 * PI;
/// Python `math.radians` multiplies by `pi / 180`; kept as that exact product.
pub const DEG2RAD: f64 = PI / 180.0;
/// Python `math.degrees` multiplies by `180 / pi`; kept as that exact product.
pub const RAD2DEG: f64 = 180.0 / PI;

/// No eC counterpart: a wrapper over `f64::sin`, kept here so that every such call can
/// move to `libm` in one place.
#[inline]
pub fn sin(x: f64) -> f64 {
    x.sin()
}
/// No eC counterpart: a wrapper over `f64::cos`, kept here so that every such call can
/// move to `libm` in one place.
#[inline]
pub fn cos(x: f64) -> f64 {
    x.cos()
}
/// No eC counterpart: a wrapper over `f64::asin`, kept here so that every such call can
/// move to `libm` in one place.
#[inline]
pub fn asin(x: f64) -> f64 {
    x.asin()
}
/// No eC counterpart: a wrapper over `f64::acos`, kept here so that every such call can
/// move to `libm` in one place.
#[inline]
pub fn acos(x: f64) -> f64 {
    x.acos()
}
/// No eC counterpart: a wrapper over `f64::atan`, kept here so that every such call can
/// move to `libm` in one place.
#[inline]
pub fn atan(x: f64) -> f64 {
    x.atan()
}
/// No eC counterpart: a wrapper over `f64::atan2`, kept here so that every such call can
/// move to `libm` in one place.
#[inline]
pub fn atan2(y: f64, x: f64) -> f64 {
    y.atan2(x)
}
/// No eC counterpart: a wrapper over `f64::sqrt`, kept here so that every such call can
/// move to `libm` in one place.
#[inline]
pub fn sqrt(x: f64) -> f64 {
    x.sqrt()
}
/// No eC counterpart: a wrapper over `f64::ln`, kept here so that every such call can
/// move to `libm` in one place.
#[inline]
pub fn ln(x: f64) -> f64 {
    x.ln()
}
/// `x` raised to an integer power, as Python's `7 ** n` in the out-of-table
/// fallback of `hex_a7.pow7`. It lives here rather than at the call site
/// because `f64::powi` is a libm-backed operation like the trigonometry above,
/// so it has to move together with them if `std` is ever swapped for `libm`.
///
/// No eC counterpart: a wrapper over `f64::powi`, kept here with the trigonometry.
#[inline]
pub fn powi(x: f64, n: i32) -> f64 {
    x.powi(n)
}

/// `x * (pi / 180)`, as Python's `math.radians`. Not `f64::to_radians`, which
/// takes a different constant path and so is not guaranteed to round to the
/// same bit pattern; the port must match Python's own product exactly.
///
/// No eC counterpart: it reproduces Python's `math.radians`, as py4dggs computes it.
#[inline]
pub fn radians(deg: f64) -> f64 {
    deg * DEG2RAD
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn radians_matches_python_math() {
        // Python: math.radians(x) = x * (pi / 180)
        assert_eq!(radians(38.7223), 38.7223 * (PI / 180.0));
    }
}
