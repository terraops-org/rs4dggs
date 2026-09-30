#!/usr/bin/env bash
# The local gate: a milestone is done only when this passes. Actions minutes cost money,
# so this, not CI, is where verification happens.
set -euo pipefail
cd "$(dirname "$0")/.."

: "${DGGAL_SITE_PACKAGES:?set DGGAL_SITE_PACKAGES to a site-packages with dggal==0.0.6 installed (pip install dggal==0.0.6)}"
for lib in "$DGGAL_SITE_PACKAGES/dggal/lib/libdggal.so" "$DGGAL_SITE_PACKAGES/ecrt/lib/libecrt.so"; do
  [[ -f "$lib" ]] || { echo "missing $lib" >&2; exit 1; }
done
export LD_LIBRARY_PATH="$DGGAL_SITE_PACKAGES/dggal/lib:$DGGAL_SITE_PACKAGES/ecrt/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

# The operation orders that crates/rs4dggs takes from the compiled engine, in
# src/topologies/hex_a7.rs, src/topologies/hex_a3.rs, src/topologies/hex_a3_subzones.rs,
# src/fivebysix.rs and src/projections/icovertex.rs, were read out of one build of DGGAL, and
# every citation there names an address in it. Another build may
# compile those same functions differently, and the symptom would be disagreements with no
# stated cause, so the library the oracle links is compared with the BuildID the citations
# themselves record. The recorded value is read from every file under crates/rs4dggs that
# cites one, the README among them, and those files must agree, so that the value is kept in
# the citations and nowhere else.
# A citation is a forty-digit hexadecimal string in backticks, which in those files is only ever
# a BuildID.
echo "== DGGAL build identity"
CITED_IN=()
while IFS= read -r f; do CITED_IN+=("$f"); done < <(
  grep -rlE '`[0-9a-f]{40}`' crates/rs4dggs/src crates/rs4dggs/README.md | sort || true)
if (( ${#CITED_IN[@]} == 0 )); then
  echo "no BuildID is recorded under crates/rs4dggs, so the compiled orders name no library" >&2; exit 1
fi
CITING="${CITED_IN[*]}"
RECORDED=$(grep -ohE '`[0-9a-f]{40}`' "${CITED_IN[@]}" | tr -d '`' | sort -u)
if [[ $(wc -l <<<"$RECORDED") -ne 1 ]]; then
  echo "the files that record the DGGAL build the compiled orders were read from disagree:" >&2
  for f in "${CITED_IN[@]}"; do
    echo "  $f: $(grep -ohE '`[0-9a-f]{40}`' "$f" | tr -d '`' | sort -u | paste -sd' ')" >&2
  done
  echo "One build was read, so one BuildID is right; correct the others." >&2
  exit 1
elif ! command -v readelf >/dev/null 2>&1; then
  echo "DGGAL build identity not checked: readelf is not on the PATH (install binutils)"
else
  LINKED=$(readelf -n "$DGGAL_SITE_PACKAGES/dggal/lib/libdggal.so" 2>/dev/null \
    | grep -oE 'Build ID: [0-9a-f]{40}' | grep -oE '[0-9a-f]{40}' || true)
  if [[ -z "$LINKED" ]]; then
    echo "$DGGAL_SITE_PACKAGES/dggal/lib/libdggal.so carries no BuildID note, so it cannot be" >&2
    echo "the build the compiled operation orders were read from, which does carry one; that" >&2
    echo "build is recorded in $CITING. Read the orders out of this library and correct the" >&2
    echo "citations, or point DGGAL_SITE_PACKAGES at the recorded build, before trusting the suites." >&2
    exit 1
  elif [[ "$LINKED" != "$RECORDED" ]]; then
    echo "the oracle links libdggal.so with BuildID $LINKED, where the compiled operation orders" >&2
    echo "were read from BuildID $RECORDED, as recorded in $CITING." >&2
    echo "Those orders, and the addresses cited beside them, belong to that build; this one may" >&2
    echo "compile the same functions differently. Read the orders out of this library and correct" >&2
    echo "the citations before trusting the suites, rather than changing the recorded value alone." >&2
    exit 1
  else
    echo "libdggal.so BuildID $LINKED, the build the compiled orders were read from, as recorded in"
    echo "  $CITING"
  fi
fi

# Our crates only: the vendored DGGAL bindings under crates/dggal-oracle/vendor are unmodified
# upstream code and must not be reformatted, linted or doctested. -p already scopes fmt and
# test to our own crates, but Cargo runs clippy's RUSTC_WORKSPACE_WRAPPER over every
# workspace member regardless of -p, so without --no-deps the vendored path dependencies
# would be linted (and fail) too; --no-deps keeps clippy's lint pass on the selected
# packages only, while the vendored crates still build normally for linking.
CRATES=(-p rs4dggs -p rs4dggs-cli -p dggal-oracle)
# `cargo fmt` is scoped to the owned crates deliberately. The four crates under
# crates/dggal-oracle/vendor are verbatim upstream DGGAL and eCere binding sources, kept
# byte-identical so that they can be diffed against a new release; reformatting them would
# destroy that property and create a diff nobody can review.
echo "== fmt";     cargo fmt "${CRATES[@]}" --check
# --features rs4dggs/oracle, not a bare --features oracle: with several packages selected by
# -p, an unqualified feature name is only accepted when every selected package declares
# it, and the others do not. The oracle feature is rs4dggs's own; the qualified form
# says so and reaches the two inline oracle-dependent tests, in topologies/hex_a7.rs and
# topologies/hex_a3_subzones.rs, and the six integration suites under tests/, one per grid,
# so that clippy and test cover them here, in the one place the feature is meaningful,
# exactly as they did before it existed.
echo "== clippy";  cargo clippy "${CRATES[@]}" --all-targets --no-deps --features rs4dggs/oracle,rs4dggs-cli/oracle -- -D warnings
echo "== test";    cargo test "${CRATES[@]}" --features rs4dggs/oracle,rs4dggs-cli/oracle
echo "== doc";      RUSTDOCFLAGS="-D warnings" cargo doc -p rs4dggs --no-deps
echo "== trig only in math.rs"
# Both call forms are caught: the method form `x.sin()` and the path form `f64::sin(x)`.
if grep -rnE '(\.|f64::)(sin|cos|tan|asin|acos|atan|atan2|sqrt|hypot|cbrt|exp|ln|powf|powi|sinh|cosh|tanh|to_radians|to_degrees)\(' \
     crates/rs4dggs/src --include='*.rs' | grep -v '^crates/rs4dggs/src/math.rs:'; then
  echo "trig/sqrt/powi call outside math.rs (see above)" >&2; exit 1
fi
echo "== no algebraic floating point"
# Rust's algebraic methods (`algebraic_add`, `algebraic_mul` and their siblings) allow the
# compiler to reassociate and contract floating-point arithmetic as it sees fit, which is what
# C's -ffast-math does. The crate agrees with one gcc build of DGGAL bit for bit because it
# performs that build's own order explicitly, in strict IEEE arithmetic; with these methods the
# order would be LLVM's, and could change with the compiler, the optimisation level or the
# target. Not even math.rs may use them.
if grep -rnE '(\.|::)algebraic_[a-z]+[[:space:]]*\(' crates/rs4dggs/src --include='*.rs'; then
  echo "algebraic floating-point call (see above): the compiled order must be written out" >&2; exit 1
fi
echo "== package"; cargo package -p rs4dggs --allow-dirty --quiet
# Cargo cannot verify the tool's package until rs4dggs is on crates.io, because it
# resolves the path dependency's version against the registry; so the file set is checked now
# and the full check waits for publication.
cargo package -p rs4dggs-cli --allow-dirty --list >/dev/null
echo "rs4dggs-cli: file set checked; its full package check waits for rs4dggs on crates.io"
echo "== README"
# The repository's front page and the crate's are one document kept in two places: GitHub
# renders the root copy, crates.io and docs.rs the crate's. The crate's is the one to edit.
if ! cmp -s README.md crates/rs4dggs/README.md; then
  echo "README.md differs from crates/rs4dggs/README.md; edit the crate's copy and copy it to the root" >&2
  exit 1
fi
echo "== without DGGAL"
# A clone without the engine must still build, test, document and run the crate itself:
# only the test-only oracle needs DGGAL. A target directory of its own, so that the
# oracle's build scripts, which read the variable, are not re-run on every gate.
nodggal() { env -u DGGAL_SITE_PACKAGES -u LD_LIBRARY_PATH CARGO_TARGET_DIR=target/no-dggal "$@"; }
nodggal cargo test -q >/dev/null
nodggal cargo doc -q --no-deps
nodggal cargo run -q --example quantise -- 5 38.7223 -9.1393 >/dev/null
echo "== examples"
# Each program runs, answers what it should, and refuses bad input without panicking.
# Output is captured before it is searched: piping into `grep -q` could fail the pipeline
# under `set -o pipefail` if grep exits before the program has finished writing.
ex() { cargo run -q -p rs4dggs --example "$@"; }
out=$(ex quantise -- 5 38.7223 -9.1393); grep -q $'\t0064156\t' <<<"$out"
ex hierarchy >/dev/null
out=$(ex hierarchy -- 00); grep -q '6 children' <<<"$out"
out=$(ex hierarchy -- 006415600000000000000); grep -q 'no children: resolution 19' <<<"$out"
# Each polygon's ring must be closed, its last point repeating its first, as GeoJSON requires.
rings_closed='import json,sys
for f in json.load(sys.stdin)["features"]:
    for r in f["geometry"]["coordinates"]:
        assert len(r) >= 4 and r[0] == r[-1], f["properties"]'
ex geojson | python3 -c "$rings_closed"
ex geojson -- 0 179.99 3 | python3 -c "$rings_closed"
out=$(ex service -- rtea7h 006415650342 2); grep -q '^ring 2: ' <<<"$out"
# A refusal is exit status 1, which `cargo run` passes through; a panic or a failed build is
# 101, and success is 0, so anything but 1 fails the gate. The last two cases are a point in
# no cell and the null zone's own text, which must be refused rather than drawn at (0, 0).
for bad in "quantise" "quantise -- x 1 2" "service -- H3 00 1" "service -- igeo7 0064158 1" \
           "geojson -- a b c" "hierarchy -- 0064158" "hierarchy -- 00 00" \
           "geojson -- 89.99912914725516 9.333761199597859 17" "hierarchy -- (null)"; do
  rc=0
  # shellcheck disable=SC2086  # the words of $bad are the program's arguments
  ex $bad >/dev/null 2>&1 || rc=$?
  [[ $rc -eq 1 ]] || { echo "example did not refuse cleanly (exit $rc): $bad" >&2; exit 1; }
done
echo "== the tool"
# The tool runs, and refuses a malformed invocation with exit status 2 rather than a panic.
cargo run -q -p rs4dggs-cli -- igeo7 info >/dev/null
for bad in "h3 info" "igeo7 frobnicate"; do
  rc=0
  # shellcheck disable=SC2086  # the words of $bad are the tool's arguments
  cargo run -q -p rs4dggs-cli -- $bad >/dev/null 2>&1 || rc=$?
  [[ $rc -eq 2 ]] || { echo "the tool did not refuse cleanly (exit $rc): $bad" >&2; exit 1; }
done
echo "== declared minimum (rust-version = 1.85, rs4dggs and rs4dggs-cli; see Cargo.toml)"
# Not the whole gate above: the test-only oracle (dggal-oracle, its vendored ecrt) needs
# rust 1.87 for integer_sign_cast, so this checks the published crates on their own, and
# only when 1.85 is actually installed. Never silent: absent, it says so rather than
# skipping without a trace.
if ! command -v rustup >/dev/null 2>&1; then
  echo "declared minimum not checked: rustup is not on the PATH, so no toolchain can be selected by name"
elif ! rustup toolchain list 2>/dev/null | grep -q '^1\.85'; then
  echo "declared minimum not checked: rust 1.85 is not installed (rustup toolchain install 1.85 --profile minimal)"
else
  cargo +1.85 check -p rs4dggs -p rs4dggs-cli
fi
# Local additions this checkout keeps outside the published set, run here, where present.
if [[ -x check.local.sh ]]; then ./check.local.sh; fi
echo "check.sh: all green"
