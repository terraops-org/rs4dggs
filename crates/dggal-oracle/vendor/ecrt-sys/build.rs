// Compiles the eC runtime's C binding layer and links the wheel's libecrt.so.
fn main() {
    // Re-run whenever the variable changes, including when it is set after a build without it.
    println!("cargo:rerun-if-env-changed=DGGAL_SITE_PACKAGES");
    // Without the engine there is nothing to build against. Only the oracle's own tests and the
    // crate's `oracle` feature link this crate, and scripts/check.sh refuses to run without the
    // variable, so a clone without DGGAL can still build, test and document the published crate.
    let Ok(site) = std::env::var("DGGAL_SITE_PACKAGES") else {
        println!("cargo:warning=DGGAL_SITE_PACKAGES is not set: the DGGAL oracle is not built (see scripts/check.sh)");
        return;
    };
    let dir = format!("{site}/ecrt/lib");
    assert!(std::path::Path::new(&format!("{dir}/libecrt.so")).exists(), "missing {dir}/libecrt.so");
    cc::Build::new().file("c/ecrt.c").include("c").flag("-w").opt_level(2).compile("ecrt_c");
    println!("cargo:rustc-link-search=native={dir}");
    println!("cargo:rustc-link-lib=dylib=ecrt");
}
