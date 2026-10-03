# Contributing to rs4dggs

## The local gate

Every change is verified on the contributor's own machine, by `scripts/check.sh`: formatting,
lints, the tests (among them the oracle suites, which compare the library with a live DGGAL v0.0.6
bit for bit, and the tool and the encodings with DGGAL's own `dgg`), the documentation, the
packages, the examples and the build on Rust 1.85, the
declared minimum. GitHub runs the gate on every push to `main` and every pull request
(`.github/workflows/ci.yml`), and once more on the tagged commit of a release, but without the
oracle suites (`scripts/check.sh --without-dggal`): PyPI's DGGAL wheel is another build than the one
the crate's compiled operation orders were read from. The oracle therefore runs locally alone, and a
commit that has not passed the full local gate has not been checked against DGGAL at all.

Without DGGAL, `scripts/check.sh --without-dggal` runs everything else. The oracle suites need the
very DGGAL build whose compiled operation orders the crate follows: the gate compares the BuildID of
the `libdggal.so` it links with the one recorded in the sources, and refuses any other. That build
is the maintainer's, compiled from `github.com/ecere/pydggal` at `v0.0.6`; a DGGAL 0.0.6 from PyPI,
or compiled elsewhere, is another build. The gate with DGGAL also needs `python3` with the module
`jsonschema` (`pip install jsonschema`), with which it checks the encodings' documents of zone data
against the standard's DGGS-JSON schema; the gate without DGGAL needs no more than Python's own
`json`. With the recorded build installed:

```sh
export DGGAL_SITE_PACKAGES=/path/to/its/site-packages
scripts/check.sh
```

The gate refuses to run while `RS4DGGS_MAX_MATERIALISED_SUB_ZONES` or
`RS4DGGS_MAX_EDGE_REFINEMENT` is set in the environment, with or without DGGAL: each lowers a
limit of the crate, and the tests assume both limits at their ceilings. The tests that compare the
command-line tool and the encodings with DGGAL's own `dgg` run it under the `timeout` of the
coreutils, which must therefore be on the path; they run only where DGGAL is configured.

## The pre-commit hook

The repository carries a hook that runs the gate before every commit that changes the code or a
published file. Switch it on once in each clone:

```sh
git config core.hooksPath .githooks
git config rs4dggs.dggalSitePackages "$DGGAL_SITE_PACKAGES"
```

The hook checks the working tree, so changes left unstaged are checked with the staged ones. A
commit that touches neither the code nor a published file passes without the gate. Where no DGGAL
is configured, the hook runs `scripts/check.sh --without-dggal` and says so. `git commit
--no-verify` skips the hook when there is good reason.

## Releases

A release is made by pushing a tag that names the workspace version in `Cargo.toml`, for example
`v0.1.0` for version `0.1.0`. The workflow `.github/workflows/release.yml` then runs the gate, save
the oracle suites, on a clean machine, packages the three crates, builds the command-line tool for
Linux, macOS and Windows, and publishes a GitHub Release carrying every file and its SHA-256
checksum. The crates are not published to crates.io by the workflow.
