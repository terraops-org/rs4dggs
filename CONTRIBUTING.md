# Contributing to rs4dggs

## The local gate

Every change is verified on the contributor's own machine, by `scripts/check.sh`: formatting,
lints, the tests (among them the oracle suites, which compare the crate with a live DGGAL v0.0.6
bit for bit), the documentation, the packages, the examples and the build on Rust 1.85, the
declared minimum. GitHub runs the same gate again on every push to `main` and every pull request
(`.github/workflows/ci.yml`), and once more on the tagged commit of a release; the local gate is
the quicker of the three, and the one to pass before pushing.

The oracle needs DGGAL's own libraries, from its Python wheel:

```sh
python3 -m venv ~/.venvs/dggal && ~/.venvs/dggal/bin/pip install dggal==0.0.6
export DGGAL_SITE_PACKAGES=$(~/.venvs/dggal/bin/python -c 'import site; print(site.getsitepackages()[0])')
scripts/check.sh
```

## The pre-commit hook

The repository carries a hook that runs the gate before every commit that changes the code or a
published file. Switch it on once in each clone:

```sh
git config core.hooksPath .githooks
git config rs4dggs.dggalSitePackages "$DGGAL_SITE_PACKAGES"
```

The hook checks the working tree, so changes left unstaged are checked with the staged ones. A
commit that touches neither the code nor a published file passes without the gate, and
`git commit --no-verify` skips it when there is good reason.

## Releases

A release is made by pushing a tag that names the workspace version in `Cargo.toml`, for example
`v0.1.0` for version `0.1.0`. The workflow `.github/workflows/release.yml` then runs the gate on a
clean machine, packages both crates, builds the command-line tool for Linux, macOS and Windows, and
publishes a GitHub Release carrying every file and its SHA-256 checksum. The crates are not
published to crates.io by the workflow.
