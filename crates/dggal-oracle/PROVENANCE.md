# Provenance

Test-only. Never a dependency of a published crate.

| Vendored file | Upstream source (unmodified, tracked in the upstream repository) |
|---|---|
| `vendor/ecrt-sys/src/lib.rs` | `github.com/ecere/eC`, commit `8a633973133bb007e446bead00f30f5492d5a23b`, `bindings/rust/ecrt_cffi.rs` |
| `vendor/ecrt/src/lib.rs` | `github.com/ecere/eC`, same commit, `bindings/rust/ecrt.rs` |
| `vendor/ecrt-sys/c/ecrt.{c,h}` | `github.com/ecere/eC`, same commit, `bindings/c/` |
| `vendor/dggal-sys/src/lib.rs` | `github.com/ecere/dggal`, tag `v0.0.6` (commit `c323c4c444522f16fd6eba0c56ec65714a147c8c`), `bindings/rust/dggal_cffi.rs` |
| `vendor/dggal/src/lib.rs` | `github.com/ecere/dggal`, tag `v0.0.6`, `bindings/rust/dggal.rs` |
| `vendor/dggal-sys/c/dggal.{c,h}` | `github.com/ecere/dggal`, tag `v0.0.6`, `bindings/c/` |

DGGAL version 0.0.6, BSD-3-Clause (`LICENSE`, Ecere Corporation). eC is built from the `main`
branch, ahead of its own `0.0.5` tag, at the commit named above; DGGAL is built at its own tagged
release, `v0.0.6`, following the build instructions in DGGAL's `BUILDING.md`: clone both
repositories into a `dgbuild` directory, then build eC before DGGAL. Checksums, verified against
both repositories at the commits named above: `vendor/SHA256SUMS`.

Runtime libraries are not vendored: `libdggal.so` and `libecrt.so` come from the `dggal==0.0.6`
pip wheel via `DGGAL_SITE_PACKAGES`. A `dgbuild` built from source instead, as `BUILDING.md`
describes, places 32-bit libraries under `dgbuild/eC/obj/linux.x32` and
`dgbuild/dggal/obj/release.linux.x32`, which do not link on x86-64.

DGGAL class names: IGEO7 = `ISEA7H_Z7`, IVEA7H = `IVEA7H_Z7`, RTEA7H = `RTEA7H_Z7`, and `ISEA3H`,
`IVEA3H`, `RTEA3H`. Angles in the raw bindings are radians.
